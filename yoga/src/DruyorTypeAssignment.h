
#pragma once

#ifdef _WIN32
#pragma push_macro("interface")
#pragma push_macro("IN")
#pragma push_macro("OUT")
#undef interface
#undef IN
#undef OUT
#endif

#include <parfait/SyncPattern.h>
#include "Connectivity.h"
#include "DonorCollector.h"
#include "PartitionInfo.h"
#include "ScalableHoleMap.h"
#include "YogaMesh.h"
#include "YogaStatuses.h"

namespace YOGA {
class DruyorTypeAssignment {
  public:
    DruyorTypeAssignment(const YogaMesh& mesh,
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
                         MessagePasser mp);

    void turnOutIntoReceptorIfValidDonorsExist(std::vector<StatusKeeper>& node_statuses);

    static std::vector<StatusKeeper> getNodeStatuses(const YogaMesh& mesh,
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
                                                   MessagePasser mp);

    void performSanityChecks(const std::vector<StatusKeeper>& statuses);

    std::vector<StatusKeeper> determineNodeStatuses2();

    std::vector<Parfait::Extent<double>> generateNodeExtents();

    void markNodesOfStraddlingCellsIn(std::vector<StatusKeeper>& statuses,
        const std::vector<bool>& is_mine);

    void markMandatoryReceptors(std::vector<StatusKeeper>& node_statuses,
            const std::vector<bool>& is_node_mine,
            const YogaMesh& m);
    void markHolePointsOut(std::vector<StatusKeeper>& node_statuses,
                           const std::vector<bool>& is_node_mine,
                           const std::vector<long>& hole_nodes);

    void markSurfaceNodesIn(std::vector<StatusKeeper>& node_statuses,
                    const std::vector<bool>& is_node_mine);
    void markDefiniteInPoints(std::vector<StatusKeeper>& statuses,
                              const std::vector<bool>& is_node_mine);
    std::vector<bool> getIsReceptor(const std::vector<bool>& is_node_mine);

    void markNodesClosestToTheirGeometry(std::vector<StatusKeeper>& node_statuses,
                                         const std::vector<bool>& is_node_mine,
                                         const PartitionInfo& partitionInfo
                                         );

    std::vector<bool> createMandatoryReceptorMask(const Parfait::CartBlock& block,
                                                 int component_id,
                                                 const std::vector<StatusKeeper>& statuses,
                                                 const std::vector<Parfait::Extent<double>>& node_extents);
    void markCandidateReceptors(std::vector<StatusKeeper>& node_statuses, const std::vector<bool>& is_node_mine);

    bool doSelectedNodesContain(const std::vector<int>& selected_nodes,
                                const NodeStatus s,
                                const std::vector<StatusKeeper>& node_statuses) const;

    bool hasValidDonor(const Receptor& r);

    void updateDonorValidity(const YogaMesh& mesh,
                             std::vector<Receptor>& receptors,
                             std::vector<StatusKeeper>& node_statuses,
                             const std::map<long, int>& g2l);

    void markNeighborsOfMandatoryReceptors(std::vector<StatusKeeper>& statuses,
        const std::vector<bool>& is_mine);

    struct StatusCounts{
        long in = 0;
        long out = 0;
        long candidate = 0;
        long receptor = 0;
        long orphan = 0;
        long mandatory_receptor = 0;
        long unknown = 0;
    };
  private:
    const YogaMesh& mesh;
    const int extra_layers_for_bcs;
    std::vector<long> hole_nodes;
    std::vector<Receptor>& receptors;
    const std::map<long, int>& globalToLocal;
    const Parfait::SyncPattern& sync_pattern;
    const PartitionInfo& partition_info;
    const MeshSystemInfo& mesh_system_info;
    const int max_cart_image_cells;
    MessagePasser mp;
    std::vector<std::vector<int>> node_to_node;


    int countStatus(const std::vector<StatusKeeper>& statuses, const NodeStatus s,const std::vector<bool>& is_mine);
    double getMinDonorDistance(int local_node_id,const std::vector<CandidateDonor>& candidates);
    bool has_at_least_one_valid_node(const std::vector<YOGA::StatusKeeper>& node_statuses, const std::vector<int>& cell);
    std::vector<bool> getIsNodeMine();
    std::vector<std::vector<int>> verifyDonorValidityLocally(const std::vector<std::vector<int>>& cell_ids_to_check_from_ranks,
                                                                                   const std::vector<StatusKeeper>& node_statuses);
    void applyDonorValidityUpdate(const std::vector<std::vector<int>>& donor_validity,
                                  const std::vector<std::vector<int>>& cell_ids_to_check,
                                  const std::vector<std::vector<int>>& receptor_indices,
                                  const std::vector<std::vector<int>>& donor_indices);

    std::vector<std::pair<int, int>> getIdsOfHoleNodes(const YogaMesh& mesh,
                                                                             const std::vector<ScalableHoleMap>& h);
    void printCounts(const std::string& msg,const StatusCounts& counts) const;
    std::vector<std::vector<int>> buildNodeToNode(const YogaMesh& m) const;
    void updateCounts(const std::vector<StatusKeeper>& node_statuses, StatusCounts& counts, const std::vector<bool>& is_mine);
    std::vector<int> getReceptorIndices() const;
    void convertMandatoryReceptors(std::vector<StatusKeeper>& statuses,
        const std::vector<int>& receptor_indices,
        const std::vector<bool>& is_mine);
    void convertUnknownToOut(std::vector<StatusKeeper>& statuses) const;
    void convertRemainingCandidates(std::vector<StatusKeeper>& statuses) const;
    void convertCandidatesToReceptorsIfHaveValidDonor(std::vector<StatusKeeper>& statuses, const std::vector<bool>& is_node_mine);
    void filterOrphansThatAreOutsideComputationalDomain(std::vector<StatusKeeper>& statuses);
    void convertMandatoryReceptorsToOutIfDonorHasBetterDistance(std::vector<StatusKeeper>& statuses,
            const std::vector<int>& receptor_indices,
            const std::vector<bool>& is_mine);
void tryToMakeNodesInIfTheyOverlapMandatoryReceptors(std::vector<StatusKeeper>& statuses);

    void syncStatuses(std::vector<StatusKeeper>& statuses) const;
bool isOrphanAndHasNoDonorCandidates(const std::vector<int>& status_vector,
                                     const std::set<int>& receptor_ids,
                                     int node) const;
};
#ifdef _WIN32
#pragma pop_macro("interface")
#pragma pop_macro("IN")
#pragma pop_macro("OUT")
#endif
}
