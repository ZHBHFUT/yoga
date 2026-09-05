#include "FragmentBalancer.h"
#include <map>
#include <stdexcept>
#include <string>
#include <MessagePasser/MessagePasser.h>
#include <parfait/RecursiveBisection.h>
#include "YogaMesh.h"
#include "PartitionInfo.h"
#include "MeshSystemInfo.h"
#include "parfait/Inspector.h"
#include "VoxelFragment.h"
#include "OverDecomposer.h"

namespace YOGA {

#ifdef _WIN32
// MS-MPI + sparse Balance partitions can trigger findOwner failures in parallel RCB.
static std::vector<int> windowsGatherSerialRcb(MessagePasser mp,
                                               const Agglomeration& agglomeration,
                                               int target_partitions,
                                               double tol) {
    const int local_n = static_cast<int>(agglomeration.points.size());
    auto counts = mp.Gather(local_n);
    int global_n = 0;
    for (int c : counts) {
        global_n += c;
    }
    int begin = 0;
    for (int r = 0; r < mp.Rank(); ++r) {
        begin += counts[r];
    }

    std::vector<Parfait::Point<double>> global_points;
    std::vector<int> global_part;
    auto points_by_rank = mp.Gather(agglomeration.points);
    if (mp.Rank() == 0) {
        global_points.reserve(static_cast<size_t>(global_n));
        for (int r = 0; r < mp.NumberOfProcesses(); ++r) {
            for (const auto& p : points_by_rank[r]) {
                global_points.push_back(p);
            }
        }
        if (global_n > 0) {
            global_part = Parfait::recursiveBisection(global_points, target_partitions, tol);
        }
    }
    mp.Broadcast(global_part, global_n, 0);

    std::vector<int> local_part;
    if (local_n > 0) {
        local_part.assign(global_part.begin() + begin, global_part.begin() + begin + local_n);
    }
    return local_part;
}
#endif

std::pair<FragmentMap, AffinityMap> createAndBalanceFragments(MessagePasser mp,
                                                              const YogaMesh& mesh,
                                                              const PartitionInfo& partition_info,
                                                              const MeshSystemInfo& mesh_system_info,
                                                              const std::map<long, int>& g2l,
                                                              int rcb_agglom_ncells,
                                                              Parfait::Inspector& inspector) {
    auto agglomeration = agglomerateCells(mesh, inspector, partition_info, mesh_system_info, rcb_agglom_ncells);
    int target_partitions = mp.NumberOfProcesses();
    Tracer::begin("parallel rcb");
    std::vector<int> part;
#if defined(_WIN32)
    if (mp.NumberOfProcesses() > 1) {
        part = windowsGatherSerialRcb(mp, agglomeration, target_partitions, 1.0e-4);
    } else {
        part = Parfait::recursiveBisection(agglomeration.points, target_partitions, 1.0e-4);
    }
#else
    try {
        part = Parfait::recursiveBisection(mp, agglomeration.points, target_partitions, 1.0e-4);
    } catch (const std::exception&) {
        throw;
    }
#endif
    Tracer::end("parallel rcb");
    Tracer::begin("map to ranks");
    auto cell_ids_for_partitions = mapCellIdsToRanks(agglomeration.ids, part);
    Tracer::end("map to ranks");
    Tracer::begin("create exchange fragments");
    auto fragments_for_ranks = extractFragmentsForRanks(mp, mesh, cell_ids_for_partitions);
    auto node_fragment_affinity = buildNodeAffinities(mesh, g2l, fragments_for_ranks);
    Tracer::end("create exchange fragments");
    Tracer::begin("exchange");
    auto frags_from_ranks = mp.Exchange(fragments_for_ranks, VoxelFragment::pack, VoxelFragment::unpack);
    auto affinities = exchangeNodeAffinities(mp, node_fragment_affinity);
    Tracer::end("exchange");
    return {frags_from_ranks, affinities};
}

std::map<int, std::vector<bool>> buildNodeAffinities(const YogaMesh& view,
                                                     const std::map<long, int>& g2l,
                                                     const std::map<int, VoxelFragment>& fragments_for_ranks) {
    std::vector<bool> is_node_claimed_by_fragment(view.nodeCount(), false);
    std::map<int, std::vector<bool>> node_fragment_affinity;
    for (auto& pair : fragments_for_ranks) {
        int rank = pair.first;
        auto& frag = pair.second;
        node_fragment_affinity[rank].resize(frag.transferNodes.size(), false);
        for (size_t i = 0; i < frag.transferNodes.size(); i++) {
            auto& node = frag.transferNodes[i];
            int local_id = g2l.at(node.globalId);
            if (not is_node_claimed_by_fragment[local_id]) {
                is_node_claimed_by_fragment[local_id] = true;
                node_fragment_affinity[rank][i] = true;
            }
        }
    }
    return node_fragment_affinity;
}

Agglomeration agglomerateCells(const YogaMesh& view,
                               Parfait::Inspector& inspector,
                               const PartitionInfo& partition_info,
                               const MeshSystemInfo& mesh_system_info,
                                                               int rcb_agglom_ncells) {
    inspector.begin("overdecompose");
    std::vector<int> owned_cell_ids = OverDecomposer::identifyOwnedCells(view, partition_info, mesh_system_info);
    int n_sub_partitions = OverDecomposer::calcNumberOfPartitions(owned_cell_ids.size(), rcb_agglom_ncells);
    auto cell_centers = OverDecomposer::generateCellCenters(view, owned_cell_ids);
    std::vector<int> sub_partitions;
    if (n_sub_partitions > 0) {
        sub_partitions = Parfait::recursiveBisection(cell_centers, n_sub_partitions, 1.0e-4);
    }
    inspector.end("overdecompose");

    Tracer::begin("create agglomerated nodes");
    Agglomeration agglomeration;
    agglomeration.points.resize(n_sub_partitions, {0, 0, 0});
    agglomeration.ids.resize(n_sub_partitions);
    std::vector<double> denominator(n_sub_partitions, 0.0);
    for (size_t i = 0; i < cell_centers.size(); i++) {
        int part_id = sub_partitions[i];
        agglomeration.points[part_id] += cell_centers[i];
        agglomeration.ids[part_id].push_back(owned_cell_ids[i]);
        denominator[part_id] += 1.0;
    }
    for (int i = 0; i < n_sub_partitions; i++) {
        agglomeration.points[i] *= 1.0 / std::max(1.0,denominator[i]);
    }

    Tracer::end("create agglomerated nodes");
    return agglomeration;
}

std::map<int, std::vector<bool>> exchangeNodeAffinities(
    const MessagePasser& mp, const std::map<int, std::vector<bool>>& node_fragment_affinity) {
    auto pack_bools = [](MessagePasser::Message& msg, bool v) { msg.pack(v); };
    auto unpack_bools = [](MessagePasser::Message& msg, bool& v) { msg.unpack(v); };
    std::map<int, std::vector<bool>> affinities = mp.Exchange(node_fragment_affinity, pack_bools, unpack_bools);
    return affinities;
}

std::map<int, std::vector<int>> mapCellIdsToRanks(const std::vector<std::vector<int>>& agglomerated_ids,
                                                  const std::vector<int>& part) {
    std::map<int, std::vector<int>> cell_ids_for_partitions;
    for (size_t i = 0; i < part.size(); i++) {
        int part_id = part[i];
        auto& ids_for_partition = cell_ids_for_partitions[part_id];
        auto& ids_in_blob = agglomerated_ids[i];
        ids_for_partition.insert(ids_for_partition.end(), ids_in_blob.begin(), ids_in_blob.end());
    }
    return cell_ids_for_partitions;
}

std::map<int, VoxelFragment> extractFragmentsForRanks(const MessagePasser& mp,
                                                      const YogaMesh& view,
                                                      const std::map<int, std::vector<int>>& cell_ids_for_partitions) {
    std::set<int> nodes_in_overlap_cells;
    std::map<int, VoxelFragment> fragments_for_ranks;
    auto node_bcs = PartitionInfo::createNodeBcs(view);
    for (auto& pair : cell_ids_for_partitions) {
        int target_rank = pair.first;
        auto& ids = pair.second;
        for (int id : ids) {
            std::vector<int> cell(view.numberOfNodesInCell(id));
            view.getNodesInCell(id, cell.data());
            for (int node : cell) {
                if (view.nodeOwner(node) == mp.Rank()) nodes_in_overlap_cells.insert(node);
            }
        }
        fragments_for_ranks[target_rank] = VoxelFragment(view, node_bcs, ids, mp.Rank());
    }
    return fragments_for_ranks;
}

}
