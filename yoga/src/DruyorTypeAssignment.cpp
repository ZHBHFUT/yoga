#include "yoga_win_msvc.h"
#include "DruyorTypeAssignment.h"
#include <parfait/SyncField.h>
#include "CartesianLoadBalancer.h"
#include "YogaConfiguration.h"
#include "RootPrinter.h"
#include "ScalableHoleMap.h"
#include "YogaStatuses.h"
#include "FloodFill.h"
#include "YogaToTInfinityAdapter.h"
#include "CartBlockGenerator.h"
#include "YogaMemoryProbe.h"
#include <t-infinity/Shortcuts.h>
#include <t-infinity/MeshSanityChecker.h>
#include <t-infinity/CartMesh.h>
using namespace inf;
namespace YOGA {

void visualizeHoleMaps(MessagePasser mp,const std::vector<ScalableHoleMap>& hole_maps){
    if(0 == mp.Rank()) {
        MessagePasser root_only(MPI_COMM_SELF);
        int i = 0;
        for (auto hm : hole_maps) {
            int nx = hm.block.numberOfCells_X();
            int ny = hm.block.numberOfCells_Y();
            int nz = hm.block.numberOfCells_Z();

            auto cart_mesh = CartMesh::create(nx,ny,nz,hm.block);
            auto filter = FilterFactory::volumeFilter(root_only.getCommunicator(),cart_mesh);
            auto volume_only_mesh = filter->getMesh();
            auto cell_statuses = std::make_shared<VectorFieldAdapter>("cell-statuses",FieldAttributes::Cell(),hm.cell_statuses);
            shortcut::visualize("hole-map"+std::to_string(i++),root_only,volume_only_mesh,{cell_statuses});
        }
    }
}

DruyorTypeAssignment::DruyorTypeAssignment(const YogaMesh& mesh,
                                           std::vector<Receptor>& receptors,
                                           const std::map<long, int>& g2l,
                                           Parfait::SyncPattern& syncPattern,
                                           const PartitionInfo& partitionInfo,
                                           const MeshSystemInfo& mesh_system_info,
                                           const std::vector<Parfait::Extent<double>>& component_grid_extents,
                                           const std::vector<ScalableHoleMap>& hole_maps,
                                           int extra_layers,
                                           bool should_add_max_receptors,
                                           int max_cart_image_cells,
                                           MessagePasser mp)
    : mesh(mesh),
      extra_layers_for_bcs(extra_layers),
      receptors(receptors),
      globalToLocal(g2l),
      sync_pattern(syncPattern),
      partition_info(partitionInfo),
      mesh_system_info(mesh_system_info),
      max_cart_image_cells(std::max(1, max_cart_image_cells)),
      mp(mp),
      node_to_node(Connectivity::nodeToNode(mesh))
{
    yogaMemProbe("Druyor ctor enter (after n2n)", mp.Rank());
    auto candidate_hole_points = getIdsOfHoleNodes(mesh,hole_maps);
    std::map<long,int> global_to_receptor_index;
    for(int i=0;i<long(receptors.size());i++){
        auto& r = receptors[i];
        global_to_receptor_index[r.globalId] = i;
    }
    std::vector<bool> remove_hole_point(candidate_hole_points.size(),false);
    for(size_t i=0;i<candidate_hole_points.size();i++){
        auto& hole_point = candidate_hole_points[i];
        int local_node_id = hole_point.first;
        long gid = mesh.globalNodeId(local_node_id);
        int cutting_component = hole_point.second;
        if(global_to_receptor_index.count(gid) != 0){
            int receptor_index = global_to_receptor_index.at(gid);
            auto& receptor = receptors[receptor_index];
            for (auto& d : receptor.candidateDonors) {
                if (d.component == cutting_component) {
                    remove_hole_point[i] = true;
                }
            }
        }
    }
    std::vector<long> actual_hole_ids;
    for(size_t i=0;i<candidate_hole_points.size();i++){
        if(not remove_hole_point[i]){
            int local_id = candidate_hole_points[i].first;
            actual_hole_ids.push_back(mesh.globalNodeId(local_id));
        }
    }

    std::vector<int> is_hole(mesh.nodeCount(),0);
    for(long gid:actual_hole_ids){
        auto it = globalToLocal.find(gid);
        if (it != globalToLocal.end()) {
            is_hole[it->second] = 1;
        }
    }
    syncVector(mp, is_hole, globalToLocal, sync_pattern);
    for(int i=0;i<mesh.nodeCount();i++){
        if(is_hole[i]){
            hole_nodes.push_back(mesh.globalNodeId(i));
        }
    }
    yogaMemProbe("Druyor ctor done", mp.Rank());
}

void DruyorTypeAssignment::turnOutIntoReceptorIfValidDonorsExist(std::vector<StatusKeeper>& node_statuses) {
    for (auto& r : receptors) {
        int local_id = globalToLocal.at(r.globalId);
        if (node_statuses[local_id].value() == OutNode and hasValidDonor(r))
            node_statuses[local_id].transition(FringeNode);
    }
    syncStatuses(node_statuses);
}

std::vector<StatusKeeper> DruyorTypeAssignment::getNodeStatuses(const YogaMesh& mesh,
                                                              std::vector<Receptor>& receptors,
                                                              const std::map<long, int>& g2l,
                                                              Parfait::SyncPattern& syncPattern,
                                                              const PartitionInfo& partitionInfo,
                                                              const MeshSystemInfo& system_info,
                                                              const std::vector<Parfait::Extent<double>>& component_grid_extents,
                                                              const std::vector<ScalableHoleMap>& hole_maps,
                                                              int extra_layers,
                                                              bool should_add_max_receptors,
                                                              int max_cart_image_cells,
                                                              MessagePasser mp) {
    Tracer::begin("Construct type assignment");
    DruyorTypeAssignment typeAssignment(
        mesh, receptors, g2l,
        syncPattern, partitionInfo,system_info,component_grid_extents,
        hole_maps, extra_layers,should_add_max_receptors,max_cart_image_cells,mp);
    Tracer::end("Construct type assignment");

    Tracer::begin("Determine statuses");
    auto node_statuses = typeAssignment.determineNodeStatuses2();
    Tracer::end("Determine statuses");

    if (should_add_max_receptors) typeAssignment.turnOutIntoReceptorIfValidDonorsExist(node_statuses);

    return node_statuses;
}

void DruyorTypeAssignment::performSanityChecks(const std::vector<StatusKeeper>& statuses) {
    std::vector<int> cell;
    int count = 0;
    for (int i = 0; i < mesh.numberOfCells(); i++) {
        cell.resize(mesh.numberOfNodesInCell(i));
        mesh.getNodesInCell(i, cell.data());
        if (doSelectedNodesContain(cell, InNode, statuses) and doSelectedNodesContain(cell, OutNode, statuses)) {
            count++;
            //throw std::logic_error("Found cell with both In and Out nodes");
        }
    }
    int n = 0;
    for(auto s:statuses){
        if(s.value() == NodeStatus::MandatoryReceptor) n++;
    }
    for(auto status:statuses) {
        auto s = status.value();
        if (s != NodeStatus::InNode and s != NodeStatus::OutNode and s != FringeNode and s != Orphan){
            printf("Found %i MandatoryReceptors\n",n);
            throw std::logic_error("Found a node with invalid status: "+std::to_string(s));
        }
    }
    if(count > 0){
        printf("Found %i cells with in and out nodes\n",count);
    }
}

void DruyorTypeAssignment::convertUnknownToOut(std::vector<StatusKeeper>& statuses) const {
    for (auto& s : statuses) {
        if (Unknown == s.value())
            s.transition(OutNode);
    }
    syncStatuses(statuses);
}

void DruyorTypeAssignment::convertRemainingCandidates(std::vector<StatusKeeper>& statuses) const {
    enum DistanceClass { kDistOut, kDistAdjacent, kDistNear, kDistFar };
    std::vector<int> distance_to_out(statuses.size(), kDistFar);
    for(int id=0;id<mesh.nodeCount();id++){
        if(statuses[id].value() == NodeStatus::OutNode){
            distance_to_out[id] = kDistOut;
        }
    }
    for(int id=0;id<mesh.nodeCount();id++){
        auto& nbrs = node_to_node[id];
        if(distance_to_out[id] == kDistFar and doSelectedNodesContain(nbrs,NodeStatus::OutNode,statuses)){
            distance_to_out[id] = kDistAdjacent;
        }
    }
    syncVector(mp, distance_to_out, globalToLocal, sync_pattern);
    for(int id=0;id<mesh.nodeCount();id++){
        auto& nbrs = node_to_node[id];
        bool has_adj_nbr = false;
        for(int nbr:nbrs){
            if(distance_to_out[nbr] == kDistAdjacent) has_adj_nbr = true;
        }
        if(distance_to_out[id] == kDistFar and has_adj_nbr){
            distance_to_out[id] = kDistNear;
        }
    }
    syncVector(mp, distance_to_out, globalToLocal, sync_pattern);

    int n_changed = 0;
    do {
        for (int id = 0; id < mesh.nodeCount(); id++) {
            auto& s = statuses[id];
            if (ReceptorCandidate == s.value()) {
                auto& nbrs = node_to_node[id];
                if (doSelectedNodesContain(nbrs, NodeStatus::InNode, statuses) and distance_to_out[id] == kDistFar) {
                    s.transition(InNode);
                    n_changed++;
                }
            }
        }
        syncStatuses(statuses);
        n_changed = mp.ParallelSum(n_changed);
    }while(n_changed>0);

    for(auto& s:statuses){
        if(NodeStatus::ReceptorCandidate == s.value()){
            s.transition(Orphan);
        }
    }
    syncStatuses(statuses);
}
void DruyorTypeAssignment::convertMandatoryReceptors(std::vector<StatusKeeper>& statuses,
    const std::vector<int>& receptor_indices,
    const std::vector<bool>& is_mine) {
    for (size_t i = 0; i < statuses.size(); i++) {
        if(not is_mine[i]) continue;
        if (statuses[i].value() == MandatoryReceptor) {
            int receptor_index = receptor_indices[i];
            if (receptor_index != -1) {
                auto& r = receptors[receptor_index];
                if (hasValidDonor(r)) {
                    statuses[i].transition(FringeNode);
                } else {
                    statuses[i].transition(Orphan);
                }
            } else {
                statuses[i].transition(Orphan);
            }
        }
    }
    syncStatuses(statuses);
}
std::vector<int> DruyorTypeAssignment::getReceptorIndices() const {
    std::vector<int> receptor_indices(mesh.nodeCount(), -1);
    for (int i = 0; i < long(receptors.size()); i++) {
        auto& r = receptors[i];
        int local_node_id = globalToLocal.at(r.globalId);
        receptor_indices[local_node_id] = i;
    }
    return receptor_indices;
}

void DruyorTypeAssignment::updateCounts(const std::vector<StatusKeeper>& node_statuses,
    DruyorTypeAssignment::StatusCounts& counts,
    const std::vector<bool>& is_mine) {
    Tracer::begin("update counts");
    counts.candidate = countStatus(node_statuses, ReceptorCandidate,is_mine);
    counts.out = countStatus(node_statuses, OutNode,is_mine);
    counts.in = countStatus(node_statuses, InNode,is_mine);
    counts.unknown = countStatus(node_statuses, Unknown,is_mine);
    counts.receptor = countStatus(node_statuses, FringeNode,is_mine);
    counts.orphan = countStatus(node_statuses, Orphan,is_mine);
    counts.mandatory_receptor = countStatus(node_statuses,MandatoryReceptor,is_mine);

    counts.candidate = mp.ParallelSum(counts.candidate);
    counts.out = mp.ParallelSum(counts.out);
    counts.in = mp.ParallelSum(counts.in);
    counts.unknown = mp.ParallelSum(counts.unknown);
    counts.orphan = mp.ParallelSum(counts.orphan);
    counts.receptor = mp.ParallelSum(counts.receptor);
    counts.mandatory_receptor = mp.ParallelSum(counts.mandatory_receptor);
    Tracer::end("update counts");
}
void DruyorTypeAssignment::printCounts(const std::string& msg,const DruyorTypeAssignment::StatusCounts& counts) const {
    if (mp.Rank() == 0) {
        printf("%s\n",msg.c_str());
        printf("in: %li out: %li receptor: %li orphan: %li cand: %li unk: %li, mand: %li\n",
               counts.in,
               counts.out,
               counts.receptor,
               counts.orphan,
               counts.candidate,
               counts.unknown,
               counts.mandatory_receptor);
    }
}

void DruyorTypeAssignment::markMandatoryReceptors(std::vector<StatusKeeper>& statuses,
                const std::vector<bool>& is_mine,
                const YogaMesh& m){
    std::set<int> remote_ids_to_switch;
    for(int i=0;i<m.numberOfBoundaryFaces();i++){
        if(m.getBoundaryCondition(i) == BoundaryConditions::Interpolation){
            auto face = m.getNodesInBoundaryFace(i);
            for(int node_id:face){
                if(is_mine[node_id])
                    statuses[node_id].transition(MandatoryReceptor);
                else
                    remote_ids_to_switch.insert(node_id);
            }
        }
    }
    std::vector<std::vector<long>> gids_for_ranks(mp.NumberOfProcesses());
    for(int id:remote_ids_to_switch){
        long gid = mesh.globalNodeId(id);
        int owning_rank = mesh.nodeOwner(id);
        gids_for_ranks[owning_rank].push_back(gid);
    }
    auto gids_from_ranks = mp.Exchange(gids_for_ranks);
    for(auto& gids:gids_from_ranks){
        for(long gid:gids){
            int local_id = globalToLocal.at(gid);
            statuses[local_id].transition(MandatoryReceptor);
        }
    }
    syncStatuses(statuses);
}

void DruyorTypeAssignment::markNeighborsOfMandatoryReceptors(std::vector<StatusKeeper>& statuses,
    const std::vector<bool>& is_mine){
    std::set<int> ids_to_switch;
    for(size_t i=0;i<statuses.size();i++){
        if(statuses[i].value() == MandatoryReceptor){
            for(int id:node_to_node[i]){
                if(statuses[id].value() == Unknown){
                    ids_to_switch.insert(id);
                }
            }
        }
    }
    for(int id:ids_to_switch)
        statuses[id].transition(MandatoryReceptor);
    syncStatuses(statuses);
}

void DruyorTypeAssignment::markHolePointsOut(std::vector<StatusKeeper>& node_statuses,
                                             const std::vector<bool>& is_node_mine,
                                             const std::vector<long>& hole_gids) {
    for (long global_id : hole_gids) {
        int local_id = globalToLocal.at(global_id);
        if (is_node_mine[local_id]){
            node_statuses[local_id].transition(OutNode);
        }
    }
    syncStatuses(node_statuses);
}

void DruyorTypeAssignment::markSurfaceNodesIn(std::vector<StatusKeeper>& node_statuses,
    const std::vector<bool>& is_node_mine){
    for(size_t i=0;i<node_statuses.size();i++){
        if(BoundaryConditions::Solid == partition_info.getBoundaryConditionForNode(i)){
            node_statuses[i].transition(InNode);
        }
    }
    syncStatuses(node_statuses);
}

void DruyorTypeAssignment::markDefiniteInPoints(std::vector<StatusKeeper>& node_statuses,
                                                const std::vector<bool>& is_node_mine) {
    std::set<int> ids_to_mark;
    auto is_receptor = getIsReceptor(is_node_mine);
    for (size_t i = 0; i < node_statuses.size(); i++) {
        if (not is_node_mine[i]) continue;
        if (node_statuses[i].value() == MandatoryReceptor) continue;
        const auto& nbrs = node_to_node[i];
        if (doSelectedNodesContain(nbrs, OutNode, node_statuses)){
            continue;
        }
        if (not is_receptor[i] and Unknown == node_statuses[i].value()) {
            ids_to_mark.insert(i);
        }
    }
    for(int id:ids_to_mark){
        node_statuses[id].transition(InNode);
    }
    syncStatuses(node_statuses);
}

std::vector<bool> DruyorTypeAssignment::getIsReceptor(const std::vector<bool>& is_node_mine) {
    std::vector<bool> is_receptor(is_node_mine.size(), false);
    for (auto& r : receptors) {
        int local_id = globalToLocal.at(r.globalId);
        if (is_node_mine[local_id]) is_receptor[local_id] = true;
    }
    return is_receptor;
}

void DruyorTypeAssignment::markNodesClosestToTheirGeometry(std::vector<StatusKeeper>& node_statuses,
                                                           const std::vector<bool>& is_node_mine,
                                                           const PartitionInfo& partitionInfo) {
    for (auto & r : receptors) {
        int localId = globalToLocal.at(r.globalId);
        if (not is_node_mine[localId]) continue;
        double nodeDistance = r.distance;
        if (node_statuses[localId].value() == MandatoryReceptor) continue;
        auto& nbrs = node_to_node[localId];
        if (doSelectedNodesContain(nbrs, OutNode, node_statuses)) continue;
        if (getMinDonorDistance(localId, r.candidateDonors) > nodeDistance) {
            node_statuses[localId].transition(InNode);
        }
    }
    syncStatuses(node_statuses);
}

void DruyorTypeAssignment::markCandidateReceptors(std::vector<StatusKeeper>& node_statuses,
                                                  const std::vector<bool>& is_node_mine) {
    std::vector<int> cell;
    std::set<int> ids_to_mark;
    for (int i = 0; i < mesh.numberOfCells(); i++) {
        cell.resize(mesh.numberOfNodesInCell(i));
        mesh.getNodesInCell(i, cell.data());
        bool has_in = doSelectedNodesContain(cell, InNode, node_statuses);
        bool has_unknown = doSelectedNodesContain(cell, Unknown, node_statuses);
        if (has_in and has_unknown) {
            for (auto& id : cell) {
                if (not is_node_mine[id]) continue;
                if (Unknown == node_statuses[id].value()) {
                    ids_to_mark.insert(id);
                }
            }
        }
    }
    for(int id:ids_to_mark){
        node_statuses[id].transition(ReceptorCandidate);
    }
    syncStatuses(node_statuses);
}

bool DruyorTypeAssignment::doSelectedNodesContain(const std::vector<int>& selected_nodes,
                                                  const NodeStatus s,
                                                  const std::vector<StatusKeeper>& node_statuses) const {
    for (int id : selected_nodes)
        if (s == node_statuses[id].value()) return true;
    return false;
}

void DruyorTypeAssignment::convertCandidatesToReceptorsIfHaveValidDonor(std::vector<StatusKeeper>& statuses, const std::vector<bool>& is_node_mine) {
    for (auto& r : receptors) {
        int node_id = globalToLocal.at(r.globalId);
        if (not is_node_mine[node_id]) continue;
        if (ReceptorCandidate == statuses[node_id].value()) {
            if (hasValidDonor(r)) {
                statuses[node_id].transition(FringeNode);
            }
        }
    }
    syncStatuses(statuses);
}

bool DruyorTypeAssignment::hasValidDonor(const Receptor& r) {
    for (auto& d : r.candidateDonors)
        if (d.isValid) return true;
    return false;
}

void DruyorTypeAssignment::updateDonorValidity(const YogaMesh& mesh,
                                               std::vector<Receptor>& receptors,
                                               std::vector<StatusKeeper>& node_statuses,
                                               const std::map<long, int>& g2l
                                               ) {
    std::vector<std::vector<int>> cell_ids_to_check(mp.NumberOfProcesses());
    std::vector<std::vector<int>> receptor_indices(mp.NumberOfProcesses());
    std::vector<std::vector<Parfait::Point<double>>> receptor_xyz(mp.NumberOfProcesses());
    std::vector<std::vector<int>> donor_indices(mp.NumberOfProcesses());
    for (size_t i=0;i<receptors.size();i++) {
        auto& r = receptors[i];
        for(size_t j=0;j<r.candidateDonors.size();j++){
            auto& d = r.candidateDonors[j];
            int target_rank = d.cellOwner;
            cell_ids_to_check[target_rank].push_back(d.cellId);
            receptor_indices[target_rank].push_back(i);
            receptor_xyz[target_rank].push_back(mesh.getNode<double>(g2l.at(r.globalId)));
            donor_indices[target_rank].push_back(j);
        }
    }

    auto cell_ids_to_check_from_ranks = mp.Exchange(cell_ids_to_check);
    auto xyzs_to_check_from_ranks = mp.Exchange(receptor_xyz);

    auto donor_validity = verifyDonorValidityLocally(cell_ids_to_check_from_ranks,node_statuses);
    applyDonorValidityUpdate(donor_validity, cell_ids_to_check, receptor_indices, donor_indices);
}

void DruyorTypeAssignment::applyDonorValidityUpdate(const std::vector<std::vector<int>>& donor_validity,
                                                    const std::vector<std::vector<int>>& cell_ids_to_check,
                                                    const std::vector<std::vector<int>>& receptor_indices,
                                                    const std::vector<std::vector<int>>& donor_indices){
    auto reply = mp.Exchange(donor_validity);
    for (int rank = 0; rank < mp.NumberOfProcesses(); rank++) {
        for (size_t i = 0; i < cell_ids_to_check[rank].size(); i++) {
            int receptor_index = receptor_indices[rank][i];
            int donor_index = donor_indices[rank][i];
            if (0 == reply[rank][i]) {
                receptors[receptor_index].candidateDonors[donor_index].isValid = false;
            } else {
                receptors[receptor_index].candidateDonors[donor_index].isValid = true;
            }
        }
    }
}

std::vector<std::vector<int>> DruyorTypeAssignment::verifyDonorValidityLocally(const std::vector<std::vector<int>>& cell_ids_to_check_from_ranks,
    const std::vector<StatusKeeper>& node_statuses){

    std::vector<std::vector<int>> donor_validity(mp.NumberOfProcesses());
    for (int rank = 0; rank < mp.NumberOfProcesses(); rank++) {
        std::vector<int> cell;
        auto& cell_ids = cell_ids_to_check_from_ranks[rank];
        for (int cell_id : cell_ids) {
            cell.resize(mesh.numberOfNodesInCell(cell_id));
            mesh.getNodesInCell(cell_id, cell.data());
            if (has_at_least_one_valid_node(node_statuses, cell)) {
                donor_validity[rank].push_back(1);
            }
            else {
                donor_validity[rank].push_back(0);
            }
        }
    }

    return donor_validity;
}

int DruyorTypeAssignment::countStatus(const std::vector<StatusKeeper>& statuses, const NodeStatus s,
    const std::vector<bool>& is_mine) {
    int n = 0;
    for (size_t i=0;i<statuses.size();i++)
        if (s == statuses[i].value() and is_mine[i]) n++;
    return n;
}

double DruyorTypeAssignment::getMinDonorDistance(int local_node_id,const std::vector<CandidateDonor>& candidates) {
    double d = candidates[0].distance;
    for(auto& donor:candidates){
        d = std::min(d,donor.distance);
    }
    return d;
}

bool DruyorTypeAssignment::has_at_least_one_valid_node(const std::vector<StatusKeeper>& node_statuses,
                                                       const std::vector<int>& cell) {
    for (int id : cell)
        if (InNode == node_statuses[id].value()) return true;
    return false;
}

std::vector<bool> DruyorTypeAssignment::getIsNodeMine() {
    std::vector<bool> is_mine(mesh.nodeCount(), true);
    for (int i = 0; i < mesh.nodeCount(); i++) {
        if (mp.Rank() != mesh.nodeOwner(i)) is_mine[i] = false;
    }
    return is_mine;
}


std::vector<std::pair<int, int>> DruyorTypeAssignment::getIdsOfHoleNodes(const YogaMesh& mesh,
                                                                     const std::vector<ScalableHoleMap>& h) {
    std::set<std::pair<int, int>> holes;
    for (int i = 0; i < mesh.nodeCount(); i++) {
        auto xyz = mesh.getNode<double>(i);
        int associated_component = mesh.getAssociatedComponentId(i);
        Parfait::Extent<double> e{xyz, xyz};
        for (auto& holeMap : h) {
            if (holeMap.getAssociatedComponentId() != associated_component)
                if (holeMap.doesOverlapHole(e)) holes.insert({i, holeMap.getAssociatedComponentId()});
        }
    }
    return std::vector<std::pair<int, int>>(holes.begin(), holes.end());
}

void DruyorTypeAssignment::markNodesOfStraddlingCellsIn(std::vector<StatusKeeper>& statuses,
    const std::vector<bool>& is_mine){
    std::vector<bool> is_in_overlap(mesh.nodeCount(),false);
    for(int i=0;i<mesh.nodeCount();i++){
        int component = mesh.getAssociatedComponentId(i);
        auto p = mesh.getNode<double>(i);
        for(int j=0;j<mesh_system_info.numberOfComponents();j++){
            if(j == component) continue;
            if(mesh_system_info.getComponentExtent(j).intersects(p)){
                is_in_overlap[i] = true;
            }
        }
    }

    std::set<int> ids_to_flip;
    for(int i=0;i<mesh.nodeCount();i++){
        if(not is_in_overlap[i]){
            for(int nbr:node_to_node[i]){
                if(is_in_overlap[nbr]){
                    ids_to_flip.insert(nbr);
                }
            }
        }
    }

    int layers = 0;
    for(int i=0;i<layers;i++) {
        auto original_status = statuses;
        for (int id : ids_to_flip) {
            statuses[id].transition(InNode);
        }
        syncStatuses(statuses);
        for(size_t j=0;j<statuses.size();j++){
            if(statuses[j].value() != original_status[j].value()){
                ids_to_flip.insert(j);
            }
        }
        std::set<int> nbrs_to_flip;
        for (int id : ids_to_flip) {
            for (int nbr : node_to_node[id]) {
                if (is_in_overlap[nbr]) {
                    nbrs_to_flip.insert(nbr);
                }
            }
        }
        ids_to_flip.swap(nbrs_to_flip);
    }
}

std::vector<Parfait::Extent<double>> DruyorTypeAssignment::generateNodeExtents(){
    Tracer::begin("create node extents");
    std::vector<Parfait::Extent<double>> extents(mesh.nodeCount());
    for(int i=0;i<mesh.nodeCount();i++){
        auto& e = extents[i];
        auto p = mesh.getNode<double>(i);
        e = {p,p};
        for(int nbr:node_to_node[i]){
            Parfait::ExtentBuilder::add(e,mesh.getNode<double>(nbr));
        }
    }
    Tracer::end("create node extents");
    return extents;
}

std::vector<bool> DruyorTypeAssignment::createMandatoryReceptorMask(const Parfait::CartBlock& block,
                                                                   int component_id,
                                                                   const std::vector<StatusKeeper>& statuses,
                                                                   const std::vector<Parfait::Extent<double>>& node_extents){
    Tracer::begin("mask");
    std::vector<char> tally(block.numberOfCells(),0);
    for(size_t id=0;id<statuses.size();id++){
        if(component_id != mesh.getAssociatedComponentId(id)) continue;
        if(statuses[id].value() == MandatoryReceptor){
            auto& e = node_extents[id];
            auto range = block.getRangeOfOverlappingCells(e);
            for(int i=range.lo[0];i<range.hi[0];i++){
                for(int j=range.lo[1];j<range.hi[1];j++){
                    for(int k=range.lo[2];k<range.hi[2];k++){
                        int cell = block.convert_ijk_ToCellId(i,j,k);
                        tally[cell] = 1;
                    }
                }
            }
        }
    }
    Tracer::begin("max");
    tally = mp.ElementalMax(tally);
    Tracer::end("max");
    Tracer::end("mask");
    return std::vector<bool>(tally.begin(),tally.end());
}

std::vector<int> toIntVector(const std::vector<StatusKeeper>& s){
    std::vector<int> v;
    v.reserve(s.size());
    for(auto& keeper:s)
        v.push_back(keeper.value());
    return v;
}

std::shared_ptr<FieldInterface> getFilteredStatusField(
    const std::vector<StatusKeeper>& statuses,
    std::shared_ptr<CellIdFilter> filter){
    auto status_vector = toIntVector(statuses);
    auto status_field = std::make_shared<VectorFieldAdapter>("node-statuses",
                                                                       FieldAttributes::Node(),
                                                                       status_vector);
    return filter->apply(status_field);
}

std::vector<StatusKeeper> DruyorTypeAssignment::determineNodeStatuses2() {
    RootPrinter root_printer(mp.Rank());
    root_printer.print("Yoga: assigning node statuses\n");
    yogaMemProbe("determineNodeStatuses2 enter", mp.Rank());

    auto is_mine = getIsNodeMine();
    yogaMemProbe("after getIsNodeMine", mp.Rank());
    std::vector<StatusKeeper> statuses(mesh.nodeCount());
    StatusCounts counts;
    updateCounts(statuses,counts,is_mine);
    printCounts("Initialize statuses",counts);

    //auto tinf_mesh = std::make_shared<YogaToTInfinityAdapter>(mesh,globalToLocal,mp);
    //auto selector = std::make_shared<CompositeSelector>();
    //selector->add(std::make_shared<NodeValueSelector>(0,0,component_ids));
    //selector->add(std::make_shared<DimensionalitySelector>(3));
    //Parfait::Point<double> origin {0,0,0};
    //Parfait::Point<double> y_normal {0,1,0};
    //selector->add(std::make_shared<PlaneSelector>(origin,y_normal));
    //auto filter = std::make_shared<CellIdFilter>(mp.getCommunicator(),tinf_mesh,selector);
    //auto component_0_mesh = filter->getMesh();

    //auto filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_0.vtk",mp,component_0_mesh,{filtered_status});

    auto receptor_indices = getReceptorIndices();
    markMandatoryReceptors(statuses, is_mine, mesh);
    yogaMemProbe("after markMandatoryReceptors", mp.Rank());
    updateCounts(statuses,counts,is_mine);
    printCounts("Mark mandatory receptors",counts);

    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_1.vtk",mp,component_0_mesh,{filtered_status});

    for(int i=0;i<extra_layers_for_bcs;i++) {
        markNeighborsOfMandatoryReceptors(statuses,is_mine);
    }
    updateCounts(statuses,counts,is_mine);
    printCounts("Add receptor layers",counts);
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_2.vtk",mp,component_0_mesh,{filtered_status});

    markHolePointsOut(statuses, is_mine, hole_nodes);
    updateCounts(statuses,counts,is_mine);
    printCounts("Mark hole points",counts);
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_3.vtk",mp,component_0_mesh,{filtered_status});


    tryToMakeNodesInIfTheyOverlapMandatoryReceptors(statuses);
    yogaMemProbe("after tryToMakeNodesIn", mp.Rank());
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_4.vtk",mp,component_0_mesh,{filtered_status});

    updateCounts(statuses,counts,is_mine);
    printCounts("Improve multi-overlap regions",counts);
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_5.vtk",mp,component_0_mesh,{filtered_status});



    markNodesOfStraddlingCellsIn(statuses, is_mine);
    updateCounts(statuses,counts,is_mine);
    printCounts("Mark nodes in straddling cells",counts);
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_6.vtk",mp,component_0_mesh,{filtered_status});

    markSurfaceNodesIn(statuses, is_mine);
    updateCounts(statuses,counts,is_mine);
    printCounts("Mark surface nodes",counts);
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_7.vtk",mp,component_0_mesh,{filtered_status});

    markNodesClosestToTheirGeometry(statuses, is_mine,partition_info);
    updateCounts(statuses,counts,is_mine);
    printCounts("Apply distance criteria",counts);
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_8.vtk",mp,component_0_mesh,{filtered_status});

    markDefiniteInPoints(statuses, is_mine);
    updateCounts(statuses,counts,is_mine);
    printCounts("Mark definite in points",counts);
    updateDonorValidity(mesh,receptors,statuses,globalToLocal);
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_9.vtk",mp,component_0_mesh,{filtered_status});

    markCandidateReceptors(statuses, is_mine);
    updateCounts(statuses,counts,is_mine);
    printCounts("Mark Candidate receptors",counts);
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_10.vtk",mp,component_0_mesh,{filtered_status});


    updateDonorValidity(mesh,receptors,statuses,globalToLocal);
    convertCandidatesToReceptorsIfHaveValidDonor(statuses,is_mine);
    updateCounts(statuses,counts,is_mine);
    printCounts("Convert candidates",counts);
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_11.vtk",mp,component_0_mesh,{filtered_status});

    updateDonorValidity(mesh,receptors,statuses,globalToLocal);
    convertMandatoryReceptorsToOutIfDonorHasBetterDistance(statuses,receptor_indices,is_mine);
    updateCounts(statuses,counts,is_mine);
    printCounts("Discard mandatory receptors based on distance criteria",counts);
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_12.vtk",mp,component_0_mesh,{filtered_status});

    convertMandatoryReceptors(statuses,receptor_indices,is_mine);
    updateCounts(statuses,counts,is_mine);
    printCounts("Convert mandatory receptors",counts);
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_13.vtk",mp,component_0_mesh,{filtered_status});

    convertUnknownToOut(statuses);
    convertRemainingCandidates(statuses);
    yogaMemProbe("after convertRemainingCandidates", mp.Rank());
    updateCounts(statuses,counts,is_mine);
    printCounts("Convert remaining",counts);
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_14.vtk",mp,component_0_mesh,{filtered_status});

    filterOrphansThatAreOutsideComputationalDomain(statuses);
    yogaMemProbe("after filterOrphans", mp.Rank());
    updateCounts(statuses,counts,is_mine);
    printCounts("Filter orphans outside computational domain",counts);
    //filtered_status = getFilteredStatusField(statuses,filter);
    //shortcut::visualize("step_15.vtk",mp,component_0_mesh,{filtered_status});
    
    performSanityChecks(statuses);
    yogaMemProbe("determineNodeStatuses2 done", mp.Rank());

    return statuses;
}
void DruyorTypeAssignment::tryToMakeNodesInIfTheyOverlapMandatoryReceptors(std::vector<StatusKeeper>& statuses) {
    auto node_extents = generateNodeExtents();
    yogaMemProbe("after generateNodeExtents", mp.Rank());
    const int max_cells = max_cart_image_cells;
    std::vector<Parfait::CartBlock> mandatory_receptor_images;
    std::vector<std::vector<bool>> image_counts;
    for (int component = 0; component < mesh_system_info.numberOfComponents(); component++) {
        int has_mandatory = 0;
        for (size_t i = 0; i < statuses.size(); i++) {
            if (component != mesh.getAssociatedComponentId(i)) continue;
            if (statuses[i].value() == MandatoryReceptor) {
                has_mandatory = 1;
                break;
            }
        }
        has_mandatory = mp.ParallelMax(has_mandatory);
        auto e = mesh_system_info.getComponentExtent(component);
        mandatory_receptor_images.push_back(generateCartBlock(e, max_cells));
        if (has_mandatory) {
            auto& block = mandatory_receptor_images.back();
            image_counts.emplace_back(createMandatoryReceptorMask(block, component, statuses, node_extents));
        } else {
            image_counts.emplace_back(mandatory_receptor_images.back().numberOfCells(), false);
        }
    }
    yogaMemProbe("after mandatory receptor images", mp.Rank());

    std::vector<int> ids_to_make_in;
    for (size_t i = 0; i < statuses.size(); i++) {
        auto& s = statuses[i];
        if (s.value() == Unknown) {
            int component = mesh.getAssociatedComponentId(i);
            auto p = mesh.getNode<double>(i);
            for (int j = 0; j < mesh_system_info.numberOfComponents(); j++) {
                if (j == component) continue;
                if (not mesh_system_info.getComponentExtent(j).intersects(p)) continue;
                auto& block = mandatory_receptor_images[j];
                auto& mask = image_counts[j];
                if (not block.intersects(p)) continue;
                int overlapping_image_cell = block.getIdOfContainingCell(p.data());
                if (overlapping_image_cell >= 0 &&
                    static_cast<size_t>(overlapping_image_cell) < mask.size() &&
                    mask[overlapping_image_cell]) {
                    ids_to_make_in.push_back(i);
                    break;
                }
            }
        }
    }

    for (int id : ids_to_make_in) {
        bool has_out_nbr = false;
        for (int nbr : node_to_node[id]) {
            if (statuses[nbr].value() == OutNode) has_out_nbr = true;
        }
        if (not has_out_nbr) statuses[id].transition(InNode);
    }
    syncStatuses(statuses);
}

void DruyorTypeAssignment::filterOrphansThatAreOutsideComputationalDomain(std::vector<StatusKeeper>& statuses){
    FloodFill floodFill(node_to_node);

    std::vector<int> status_vector;
    for(auto s:statuses)
        status_vector.push_back(s.value());
    std::set<int> seeds;
    std::set<int> receptor_ids;
    for(auto& r:receptors){
        long gid = r.globalId;
        int lid = globalToLocal.at(gid);
        receptor_ids.insert(lid);
    }
    for(int i=0;i<mesh.numberOfBoundaryFaces();i++){
        if(BoundaryConditions::Interpolation == mesh.getBoundaryCondition(i)) {
            auto face = mesh.getNodesInBoundaryFace(i);
            for(int node:face)
                if(isOrphanAndHasNoDonorCandidates(status_vector, receptor_ids, node))
                    seeds.insert(node);
        }
    }
    int n = int(seeds.size());
    n = mp.ParallelSum(n);
    if(mp.Rank() == 0)
        printf("There were %i orphans on interpolation boundaries who might be inside an unclosed body\n",n);
    yogaMemProbe("filterOrphans before flood fill", mp.Rank());
    auto sync = [&](){
        Tracer::begin("syncVector");
        syncVector(mp, status_vector, globalToLocal, sync_pattern);
        Tracer::end("syncVector");
    };
    auto parallel_max = [&](int n){return mp.ParallelMax(n);};
    std::map<int,int> allowed_transitions;
    allowed_transitions[NodeStatus::Orphan] = NodeStatus::OutNode;
    allowed_transitions[NodeStatus::InNode] = NodeStatus::OutNode;
    floodFill.fill(status_vector,seeds,allowed_transitions,sync,parallel_max);
    yogaMemProbe("filterOrphans after flood fill", mp.Rank());

    for(size_t i=0;i<statuses.size();i++)
        statuses[i].current = static_cast<NodeStatus>(status_vector[i]);
}
bool DruyorTypeAssignment::isOrphanAndHasNoDonorCandidates(const std::vector<int>& status_vector,
                                                           const std::set<int>& receptor_ids,
                                                           int node) const {
    return Orphan == status_vector[node] and receptor_ids.count(node)==0;
}

void DruyorTypeAssignment::convertMandatoryReceptorsToOutIfDonorHasBetterDistance(std::vector<StatusKeeper>& statuses,
    const std::vector<int>& receptor_indices,
    const std::vector<bool>& is_mine){
    std::set<int> ids_to_mark_out;
    for(size_t i=0;i<statuses.size();i++){
        if(not is_mine[i]) continue;
        const auto& s = statuses[i];
        if(MandatoryReceptor == s.value()) {
            int receptor_index = receptor_indices[i];
            if(receptor_index < 0){
                statuses[i].transition(Orphan);
                continue;
            }
            auto& r = receptors.at(receptor_index);
            double distance = r.distance;
            bool has_valid_donor_with_better_distance = false;
            for(auto& d:r.candidateDonors){
                if(d.distance < distance and d.isValid){
                    has_valid_donor_with_better_distance = true;
                }
            }
            if(has_valid_donor_with_better_distance){
                bool has_no_in_nbrs = true;
                for(auto& nbr:node_to_node[i]){
                    if(statuses[nbr].value() == InNode)
                        has_no_in_nbrs = false;
                }
                if(has_no_in_nbrs){
                    ids_to_mark_out.insert(i);
                }
            }
        }
    }
    for(int id:ids_to_mark_out){
        statuses[id].transition(OutNode);
    }
    syncStatuses(statuses);

}

void DruyorTypeAssignment::syncStatuses(std::vector<StatusKeeper>& statuses) const {
    std::vector<NodeStatus> status_vector(statuses.size());
    for(int i=0;i<int(statuses.size());i++)
        status_vector[i] = statuses[i].value();
    syncVector(mp, status_vector, globalToLocal, sync_pattern);
    for(int i=0;i<int(statuses.size());i++)
        if(statuses[i].value() != status_vector[i])
            statuses[i].transition(status_vector[i]);
}
}
