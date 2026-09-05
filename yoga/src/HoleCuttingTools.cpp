#include "HoleCuttingTools.h"
#include "YogaMemoryProbe.h"
#include <string>

namespace YOGA {

std::vector<ScalableHoleMap> createHoleMaps(MessagePasser mp,
                                            const YogaMesh& view,
                                            const PartitionInfo& partitionInfo,
                                            const MeshSystemInfo& meshSystemInfo,
                                            int max_cells_per_hole_map) {
    std::vector<ScalableHoleMap> hole_maps;
    for (int i = 0; i < meshSystemInfo.numberOfBodies(); i++) {
        hole_maps.emplace_back(ScalableHoleMap(mp, view, partitionInfo, meshSystemInfo, i, max_cells_per_hole_map));
        yogaMemProbe(("after hole map body " + std::to_string(i)).c_str(), mp.Rank());
    }
    return hole_maps;
}

std::vector<Parfait::Extent<double>> getComponentExtents(const MeshSystemInfo& info) {
    std::vector<Parfait::Extent<double>> extents;
    for (int i = 0; i < info.numberOfComponents(); i++) extents.push_back(info.getComponentExtent(i));
    return extents;
}

bool doesOverlapOtherComponent(const std::vector<Parfait::Extent<double>>& component_extents,
                               const Parfait::Point<double>& p,
                               int component_of_point) {
    for (int i = 0; i < long(component_extents.size()); i++) {
        if (i != component_of_point and component_extents[i].intersects(p)) return true;
    }
    return false;
}

    bool doesOverlapOtherComponent(const std::vector<Parfait::Extent<double>>& component_extents,
                                   const Parfait::Extent<double>& e,
                                   int component_of_point) {
        for (int i = 0; i < long(component_extents.size()); i++) {
            if (i != component_of_point and component_extents[i].intersects(e)) return true;
        }
        return false;
    }

std::vector<Parfait::Point<double>> extractPointsForComponentGrid(
    const YogaMesh& mesh, const std::vector<Parfait::Extent<double>>& component_extents, int id) {
    std::vector<Parfait::Point<double>> points;
    for (int i = 0; i < mesh.nodeCount(); i++) {
        if (mesh.getAssociatedComponentId(i) == id) {
            auto p = mesh.getNode<double>(i);
            if (doesOverlapOtherComponent(component_extents, p, id)) points.push_back(p);
        }
    }
    return points;
}

std::vector<bool> getInitialWorkUnitMask(MessagePasser mp, std::shared_ptr<LoadBalancer> load_balancer) {
    std::vector<int> mask(mp.NumberOfProcesses(),1);
    int n = load_balancer->getRemainingVoxelCount();
    if(mp.NumberOfProcesses() > 100) {
        mask[0] = 0;
    }
    for(int i=n+1;i<mp.NumberOfProcesses();i++){
        mask[i] = 0;
    }
    mp.Broadcast(mask, 0);
    return std::vector<bool>(mask.begin(),mask.end());
}

std::vector<Parfait::Extent<double>> getInitialWorkUnits(MessagePasser mp, std::shared_ptr<LoadBalancer> load_balancer) {
    std::vector<Parfait::Extent<double>> initial_units(mp.NumberOfProcesses());
    if (mp.Rank() == 0) {
        int first_worker = 0;
        if(mp.NumberOfProcesses() > 100) first_worker = 1;
        for (int i = first_worker; i < mp.NumberOfProcesses(); i++){
            if(load_balancer->getRemainingVoxelCount()==0) break;
            initial_units[i] = load_balancer->getWorkVoxel();
        }
    }
    mp.Broadcast(initial_units, 0);
    return initial_units;
}

std::vector<Receptor> removeNonReceptors(const std::vector<Receptor>& receptors,
                                         const std::vector<NodeStatus>& node_statuses,
                                         const std::map<long, int>& global_to_local_node_id) {
    std::vector<Receptor> actual_receptors;
    for (auto& r : receptors) {
        int node_id = global_to_local_node_id.at(r.globalId);
        if (node_statuses[node_id] == FringeNode) actual_receptors.push_back(r);
    }
    return actual_receptors;
}
}
