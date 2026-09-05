#include "yoga_c_interface.h"
#include <Tracer.h>
#include <parfait/UgridReader.h>
#include <stdio.h>
#include <t-infinity/PluginLoader.h>
#include <t-infinity/VectorFieldAdapter.h>
#include <array>
#include <parfait/SyncField.h>
#include "CartesianLoadBalancer.h"
#include "YogaConfiguration.h"
#include "DistanceFieldAdapter.h"
#include "MeshSystemInfo.h"
#include "PartitionInfo.h"
#include "ReceptorUpdate.h"
#include "WeightBasedInterpolator.h"
#include "YogaInstance.h"
#include "YogaMesh.h"
#include "AssemblyViaExchange.h"
#include "OversetData.h"
#include "Connectivity.h"
#include "GraphColoring.h"
#include "ParallelColorCombinator.h"
#include "GhostSyncPatternBuilder.h"
#include "ColorSyncer.h"
#include "ComponentGridIdentifier.h"
#include <parfait/FloatingPointExceptions.h>
#include "Diagnostics.h"
#include <DomainConnectivityInfo.h>
#include <FUN3DAdjointData.h>
#include <FUN3DComplexChecker.h>
#include <C_InterfaceHelpers.h>

#ifdef YOGA_WITH_ZMQ
    #include "ZMQMessager.h"
    #include "AssemblyViaZMQPostMan.h"
#endif

using namespace YOGA;
using namespace Parfait;

void* yoga_create_instance(int comm,
                           int (*numberOfNodes)(),
                           int (*numberOfCells)(),
                           int (*numberOfNodesInCell)(int),
                           void (*getNodesInCell)(int, int*),
                           long (*getGlobalNodeId)(int),
                           int (*getComponentId)(int),
                           void (*getPoint)(int, double*),
                           int (*getRankOfNodeOwner)(int),
                           int (*getNumberOfBoundaryFaces)(),
                           int (*numberOfNodesInBoundaryFace)(int),
                           void (*getNodesInFace)(int, int*),
                           int (*getBoundaryCondition)(int),
                           int numberOfSolutionVars,
                           void (*getSolutionInCellAt)(int, double, double, double, double*),
                           void (*getSolutionAtNode)(int, double*)) {
    auto yoga_comm = MPI_Comm_f2c(comm);
    MessagePasser mp(yoga_comm);
    auto config = YogaConfiguration(mp);
    Tracer::setDebug();
    if(mp.Rank() == 0) {
        printf("Yoga, registering functions from solver\n");
        fflush(stdout);
    }
    bool is_complex = FUN3DComplexChecker::isComplex(getPoint);
    if(mp.Rank() == 0){
        if(is_complex)
            printf("Yoga: detected FUN3D running in complex-mode\n");
    }

    auto complexInterceptor = [&](int node_id,double* xyz){
        Parfait::Point<std::complex<double>> p;
        getPoint(node_id, (double*)p.data());
        xyz[0] = p[0].real();
        xyz[1] = p[1].real();
        xyz[2] = p[2].real();
    };

    Tracer::begin("Set up YogaMesh");
    auto mesh = YogaMesh();
    mesh.setNodeCount(numberOfNodes());
    if(is_complex) {
        mesh.setXyzForNodes(complexInterceptor);
        mesh.setComplexGetter(getPoint);
    } else {
        mesh.setXyzForNodes(getPoint);
    }
    mesh.setGlobalNodeIds(getGlobalNodeId);
    mesh.setOwningRankForNodes(getRankOfNodeOwner);
    mesh.setComponentIdsForNodes(getComponentId);
    mesh.fixInvalidFun3DComponentIds(mp);
    mesh.setCellCount(numberOfCells());
    mesh.setCells(numberOfNodesInCell, getNodesInCell);
    mesh.setFaceCount(getNumberOfBoundaryFaces());
    mesh.setFaces(numberOfNodesInBoundaryFace, getNodesInFace);
    mesh.setBoundaryConditions(getBoundaryCondition,numberOfNodesInBoundaryFace);
    Tracer::end("Set up YogaMesh");

    if(config.shouldDumpPartFile()){
        writeFUN3DPartitionFile(mp,"fun3d_part_file.data",mesh,true);
    }
    return new YogaInstance(mp,
                            std::move(mesh),
                            numberOfSolutionVars,
                            getSolutionInCellAt,
                            getSolutionAtNode,
                            is_complex);
}

YogaInstance& extractReference(void* yoga_instance){
    return *(YogaInstance*)yoga_instance;
}

void yoga_perform_domain_assembly(void* yoga_instance) {
    auto& instance = extractReference(yoga_instance);
    auto mp = instance.mp;
    if (mp.Rank() == 0) printf("Yoga: in perform_domain_assembly\n");
    auto config = YOGA::YogaConfiguration(mp);
    auto& mesh = instance.mesh;
    int load_balancer = config.selectedLoadBalancer();
    int target_voxel_size = config.selectedTargetVoxelSize();
    if (instance.mp.Rank() == 0) {
        printf("load balancer: %i target voxel size: %i\n", load_balancer, target_voxel_size);
    }
    int extra_layers = config.numberOfExtraLayersForInterpBcs();
    int rcb_agglom_ncells = config.rcbAgglomerationSize();
    bool should_add_max_receptors = config.shouldAddExtraReceptors();
    auto component_grid_importance = config.getComponentGridImportance();
    Parfait::disableFloatingPointExceptions();
    if(config.shouldUseZMQPath()) {
#ifdef YOGA_WITH_ZMQ
        instance.oversetData = YOGA::assemblyViaZMQPostMan(mp, mesh, load_balancer, target_voxel_size,
                                                           extra_layers,
                                                           rcb_agglom_ncells,
                                                           should_add_max_receptors,
                                                           instance.doesCellContainNode);
#else
        throw std::logic_error("Yoga not compiled with ZMQ support");
#endif
    }
    else{
        instance.oversetData = YOGA::assemblyViaExchange(mp,
                                                         mesh,
                                                         load_balancer,
                                                         target_voxel_size,
                                                         extra_layers,
                                                         rcb_agglom_ncells,
                                                         should_add_max_receptors,
                                                         component_grid_importance,
                                                         instance.doesCellContainNode);
    }
    Tracer::begin("Generate inverse receptors");
    if(instance.is_complex) {
        instance.inverse_receptors_complex =
            YOGA::generateInverseReceptors<std::complex<double>>(mp, instance.oversetData->receptors, mesh);
        instance.fun3d_adjoint_data_complex = std::make_shared<FUN3DAdjointData<std::complex<double>>>(mp,mesh,
                                                                                                       instance.oversetData->statuses,
                                                                                                       instance.oversetData->receptors);
        instance.interpolator_complex = std::make_shared<WeightBasedInterpolator<std::complex<double>>>(
            instance.nSolutionVariables, instance.inverse_receptors_complex, mesh, instance.calcWeightsForReceptor);
    } else {
        instance.inverse_receptors = YOGA::generateInverseReceptors<double>(mp, instance.oversetData->receptors, mesh);
        instance.fun3d_adjoint_data = std::make_shared<FUN3DAdjointData<double>>(mp,mesh,
                                                                                 instance.oversetData->statuses,
                                                                                 instance.oversetData->receptors);
        instance.interpolator = std::make_shared<WeightBasedInterpolator<double>>(
            instance.nSolutionVariables, instance.inverse_receptors, mesh, instance.calcWeightsForReceptor);
    }
    Tracer::end("Generate inverse receptors");

    // just leave them off...
    //Parfait::enableFloatingPointExceptions();
}

void yoga_update_receptor_solutions(void* yoga_instance) {
    auto& instance = extractReference(yoga_instance);
    Parfait::disableFloatingPointExceptions();
    Tracer::begin("Update receptor solutions");
    int nvar = instance.nSolutionVariables;
    getUpdatedSolutionsFromInterpolator(instance.mp,
                                        instance.oversetData->receptors,
                                        instance.inverse_receptors,
                                        nvar,
                                        *instance.interpolator.get(),
                                        instance.getSolutionFromSolverAtNode,
                                        instance.solution_at_nodes);
    Tracer::end("Update receptor solutions");
    // just leave them off
    //Parfait::enableFloatingPointExceptions();
}

int yoga_get_node_status(void* yoga_instance, int node_id) {
    auto& instance = extractReference(yoga_instance);
    int n = instance.mesh.nodeCount();
    if (node_id >= n) throw std::logic_error("Yoga: get_node_status: invalid node id");
    switch (instance.oversetData->statuses[node_id]) {
        case FringeNode:
            return -1;
        case InNode:
            return 1;
        case OutNode:
            return 0;
        case Orphan:
            return -2;
        default:
            throw std::logic_error(std::string(__FILE__) + ": unknown status");
    }
}

void yoga_get_solution_for_receptor(void* yoga_instance, int node_id, double* Q) {
    auto& instance = extractReference(yoga_instance);
    int nvar = instance.nSolutionVariables;
    double* q = &instance.solution_at_nodes.at(nvar * node_id);
    std::copy(q,q+nvar,Q);
}

void yoga_destroy_instance(void* yoga_instance) {
    YogaInstance* p = (YogaInstance*)yoga_instance;
    if (p->mp.Rank() == 0) printf("Yoga: cleanup\n");
    delete p;
}


int yoga_get_donor_count(void* yoga_instance, int node_id) {
    auto& yoga = extractReference(yoga_instance);
    if(yoga.is_complex)
        return yoga.fun3d_adjoint_data_complex->donorCount(node_id);
    else
        return yoga.fun3d_adjoint_data->donorCount(node_id);
}

void yoga_fill_donor_local_ids(void* yoga_instance, int node_id, int* donor_ids) {
    auto& yoga = extractReference(yoga_instance);
    if(yoga.is_complex){
        auto donors = yoga.fun3d_adjoint_data_complex->donorIndices(node_id);
        std::copy(donors.begin(), donors.end(), donor_ids);
    }
    else {
        auto donors = yoga.fun3d_adjoint_data->donorIndices(node_id);
        std::copy(donors.begin(), donors.end(), donor_ids);
    }
}

void yoga_fill_donor_owning_ranks(void* yoga_instance, int node_id, int* owning_ranks) {
    auto& yoga = extractReference(yoga_instance);
    if(yoga.is_complex){
        auto owners = yoga.fun3d_adjoint_data_complex->donorOwners(node_id);
        std::copy(owners.begin(),owners.end(),owning_ranks);
    }
    else{
        auto owners = yoga.fun3d_adjoint_data->donorOwners(node_id);
        std::copy(owners.begin(),owners.end(),owning_ranks);
    }
}

void yoga_fill_donor_weights(void* yoga_instance, int node_id, double* weights) {
    auto& yoga = extractReference(yoga_instance);
    if(yoga.is_complex){
        auto w = yoga.fun3d_adjoint_data_complex->donorWeights(node_id);
        auto nasty_ptr = (std::complex<double>*)weights;
        std::copy(w.begin(), w.end(), nasty_ptr);
    }
    else {
        auto w = yoga.fun3d_adjoint_data->donorWeights(node_id);
        std::copy(w.begin(), w.end(), weights);
    }
}

double yoga_partial_receptor_x(void* stencil, int weight_index) {
    auto ptr = (StencilWrapper*)stencil;
    int n_donors = ptr->n_donors;
    switch (n_donors) {
        case 4: return receptorPartial<4>(ptr,weight_index,Partial::DX);
        case 5: return receptorPartial<5>(ptr,weight_index,Partial::DX);
        case 6: return receptorPartial<6>(ptr,weight_index,Partial::DX);
        case 8: return receptorPartial<8>(ptr,weight_index,Partial::DX);
    }
    return 0;
}
double yoga_partial_receptor_y(void* stencil, int weight_index) {
    auto ptr = (StencilWrapper*)stencil;
    int n_donors = ptr->n_donors;
    switch (n_donors) {
        case 4: return receptorPartial<4>(ptr,weight_index,Partial::DY);
        case 5: return receptorPartial<5>(ptr,weight_index,Partial::DY);
        case 6: return receptorPartial<6>(ptr,weight_index,Partial::DY);
        case 8: return receptorPartial<8>(ptr,weight_index,Partial::DY);
    }
    return 0;
}
double yoga_partial_receptor_z(void* stencil, int weight_index) {
    auto ptr = (StencilWrapper*)stencil;
    int n_donors = ptr->n_donors;
    switch (n_donors) {
        case 4: return receptorPartial<4>(ptr,weight_index,Partial::DZ);
        case 5: return receptorPartial<5>(ptr,weight_index,Partial::DZ);
        case 6: return receptorPartial<6>(ptr,weight_index,Partial::DZ);
        case 8: return receptorPartial<8>(ptr,weight_index,Partial::DZ);
    }
    return 0;
}

double yoga_partial_donor_x(void* stencil, int weight_index, int donor_index) {
    auto ptr = (StencilWrapper*)stencil;
    int n_donors = ptr->n_donors;
    switch (n_donors) {
        case 4: return donorPartial<4>(ptr, weight_index, donor_index, Partial::DX);
        case 5: return donorPartial<5>(ptr, weight_index, donor_index, Partial::DX);
        case 6: return donorPartial<6>(ptr, weight_index, donor_index, Partial::DX);
        case 8: return donorPartial<8>(ptr, weight_index, donor_index, Partial::DX);
    }
    return 0.0;
}

double yoga_partial_donor_y(void* stencil, int weight_index, int donor_index){
    auto ptr = (StencilWrapper*)stencil;
    int n_donors = ptr->n_donors;
    switch (n_donors) {
        case 4: return donorPartial<4>(ptr, weight_index, donor_index, Partial::DY);
        case 5: return donorPartial<5>(ptr, weight_index, donor_index, Partial::DY);
        case 6: return donorPartial<6>(ptr, weight_index, donor_index, Partial::DY);
        case 8: return donorPartial<8>(ptr, weight_index, donor_index, Partial::DY);
    }
    return 0.0;
}
double yoga_partial_donor_z(void* stencil, int weight_index, int donor_index){
    auto ptr = (StencilWrapper*)stencil;
    int n_donors = ptr->n_donors;
    switch (n_donors) {
        case 4: return donorPartial<4>(ptr, weight_index, donor_index, Partial::DZ);
        case 5: return donorPartial<5>(ptr, weight_index, donor_index, Partial::DZ);
        case 6: return donorPartial<6>(ptr, weight_index, donor_index, Partial::DZ);
        case 8: return donorPartial<8>(ptr, weight_index, donor_index, Partial::DZ);
    }
    return 0.0;
}

void* yoga_create_linearized_donor_stencil(int n_donors, double* donor_xyzs, double* query_xyz){
    switch(n_donors) {
        case 4:
            return (void*)makeWrappedStencil<4>(donor_xyzs,query_xyz);
        case 5:
            return (void*)makeWrappedStencil<5>(donor_xyzs,query_xyz);
        case 6:
            return (void*)makeWrappedStencil<6>(donor_xyzs,query_xyz);
        case 8:
            return (void*)makeWrappedStencil<8>(donor_xyzs,query_xyz);
        default:
            throw std::logic_error("Unrecognized donor stencil of size "+std::to_string(n_donors));
    }
    return nullptr;
}

void yoga_destroy_linearized_donor_stencil(void* stencil){
    auto ptr = (StencilWrapper*)stencil;
    int n_donors = ptr->n_donors;
    switch (n_donors) {
        case 4:
            delete (LinearizedStencil<4>*) ptr->stencil;
            break;
        case 5:
            delete (LinearizedStencil<5>*) ptr->stencil;
            break;
        case 6:
            delete (LinearizedStencil<6>*) ptr->stencil;
            break;
        case 8:
            delete (LinearizedStencil<8>*) ptr->stencil;
            break;
    }
    delete ptr;
    stencil = nullptr;
}
