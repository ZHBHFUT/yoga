#include <Tracer.h>
#include <parfait/ExtentWriter.h>
#include <parfait/RecursiveBisection.h>
#include <parfait/Timing.h>
#include <memory>
#include "AlternateMapBuilder.h"
#include "CartesianLoadBalancer.h"
#include "YogaConfiguration.h"
#include "Connectivity.h"
#include "DistanceFieldAdapter.h"
#include "DistributedLoadBalancer.h"
#include "DruyorTypeAssignment.h"
#include "ExchangeBasedAssembly.h"
#include "GhostSyncPatternBuilder.h"
#include "GlobalToLocal.h"
#include "GridFetcher.h"
#include "HoleCutStatPrinter.h"
#include "HoleCuttingTools.h"
#include <parfait/Inspector.h>
#include "LoadBalancer.h"
#include "MeshSystemInfo.h"
#include "NanoFlannDistanceCalculator.h"
#include "OverDecomposer.h"
#include "OverlapMask.h"
#include "OversetData.h"
#include "ParallelSurface.h"
#include "PartitionInfo.h"
#include "RootPrinter.h"
#include "WorkVoxelBuilder.h"
#include "AssemblyViaExchange.h"
#include "InspectorPrinter.h"
#include "OverlapDetector.h"
#include "FragmentBalancer.h"
#include "ChunkedPointGatherer.h"
#include <parfait/LinearPartitioner.h>
#include "YogaMemoryProbe.h"

namespace YOGA{


std::vector<NodeStatus> generateNodeStatuses(MessagePasser mp,
                                             const PartitionInfo& partition_info,
                                             const MeshSystemInfo& mesh_system_info,
                                             const YogaMesh& mesh,
                                             std::vector<Receptor>& receptors,
                                             const std::map<long,int>& g2l,
                                             int extra_layers,
                                             bool should_add_max_receptors);

std::vector<NodeStatus> extractPlainStatuses(const std::vector<StatusKeeper>& status_keepers){
    std::vector<NodeStatus> plain_statuses(status_keepers.size());
    for(int i=0;i<int(status_keepers.size());i++)
        plain_statuses[i] = status_keepers[i].value();
    return plain_statuses;
}

void addNodeNeighborsToReceptors(std::vector<Receptor>& receptors,const YogaMesh& mesh, const std::map<long,int>& g2l){
    Tracer::begin("n2n");
    auto n2n = Connectivity::nodeToNode(mesh);
    Tracer::traceMemory();
    for(auto& r:receptors) {
        int local_id = g2l.at(r.globalId);
        for (int nbr : n2n[local_id]) r.nbr_gids.push_back(mesh.globalNodeId(nbr));
    }
    Tracer::traceMemory();
    Tracer::end("n2n");
}

std::map<int,std::vector<std::pair<int,int>>> mapNodeIdsToRanks(int rank,
                                                                const std::map<int,VoxelFragment>& fragments,
                                                                const OverlapDetector& overlap_detector,
                                                                const std::map<int,std::vector<bool>>& affinities);

void addWallDistanceToTransferNodes(MessagePasser mp,
                                    const YogaMesh& mesh,
                                    FragmentMap& frags_from_ranks);

void addWallDistanceToTransferNodesWithChunkedSurfaces(MessagePasser mp,
                                                       const YogaMesh& mesh,
                                                       FragmentMap& frags_from_ranks);

std::vector<std::vector<TransferNode>> buildAndExchangeQueryPoints(
    const MessagePasser& mp,
    Parfait::Inspector& inspector,
    const FragmentMap& frags_from_ranks,
    std::map<int, std::vector<std::pair<int, int>>>& node_keys_for_ranks);

std::map<int, ReceptorCollection> buildReceptorCollectionsForRanks(
    const MessagePasser& mp,
    FragmentDonorFinder& donor_finder,
    const std::vector<std::vector<TransferNode>>& query_pts_from_ranks);

std::map<int, std::vector<std::pair<int, int>>> buildNodeKeysForRanks(const MessagePasser& mp,
                                                                      const FragmentMap& frags_from_ranks,
                                                                      const AffinityMap& affinities,
                                                                      const FragmentDonorFinder& donor_finder);
std::vector<Receptor> exchangeAndUnpackDonors(const MessagePasser& mp,
                                            std::map<long, int>& g2l,
                                            const std::map<int, ReceptorCollection>& receptor_collections_for_ranks);


std::vector<Receptor> performDonorSearchViaSingleExchange(MessagePasser mp,
                                                          Parfait::Inspector& inspector,
                                                          const std::map<int,VoxelFragment>& frags_from_ranks,
                                                          std::map<int, std::vector<std::pair<int, int>>>&  node_keys_for_ranks,
                                                          FragmentDonorFinder& donor_finder,
                                                          std::map<long,int>& g2l);

std::vector<Receptor> performDonorSearchInChunksToReducePeakMemory(MessagePasser mp,
                                                                   Parfait::Inspector& inspector,
                                                                   const std::map<int,VoxelFragment>& frags_from_ranks,
                                                                   std::map<int, std::vector<std::pair<int, int>>>&  node_keys_for_ranks,
                                                                   FragmentDonorFinder& donor_finder,
                                                                   std::map<long,int>& g2l
                                                                   );

void modifyDistanceBasedOnComponentImportance(FragmentMap& map,const std::vector<int>& component_grid_importance);

std::shared_ptr<OversetData> assemblyViaExchange(MessagePasser mp,
                                                 YogaMesh& view,
                                                 int load_balancer_algorithm,
                                                 int target_voxel_size,
                                                 int extra_layers,
                                                 int rcb_agglom_ncells,
                                                 bool should_add_max_receptors,
                                                 const std::vector<int>& component_grid_importance,
                                                 std::function<bool(double*, int, double*)> is_in_cell) {
    yogaMemProbe("assemblyViaExchange enter", mp.Rank());
    mp.Barrier();
    auto before_assembly = Parfait::Now();
    Tracer::begin("Domain Assembly");
    Tracer::traceMemory();
    RootPrinter rootPrinter(mp.Rank());
    rootPrinter.print("Yoga: starting domain assembly\n");
    auto beginTrace = [](const std::string& s){
        Tracer::begin(s);
    };
    auto endTrace = [](const std::string& s){
      Tracer::end(s);
    };
    Parfait::Inspector inspector(mp,0,beginTrace,endTrace);
    Tracer::begin("partition info");
    PartitionInfo partition_info(view, mp.Rank());
    Tracer::end("partition info");
    yogaMemProbe("after PartitionInfo", mp.Rank());
    Tracer::traceMemory();

    Tracer::begin("build mesh system info");
    MeshSystemInfo mesh_system_info(mp, partition_info);
    Tracer::end("build mesh system info");
    yogaMemProbe("after MeshSystemInfo", mp.Rank());
    Tracer::traceMemory();


    auto g2l = GlobalToLocal::buildMap(view);
    auto fragments_and_affinities = createAndBalanceFragments(mp,
                                                              view,
                                                              partition_info,
                                                              mesh_system_info,
                                                              g2l,
                                                              rcb_agglom_ncells,
                                                              inspector);
    auto& frags_from_ranks = fragments_and_affinities.first;
    auto& affinities = fragments_and_affinities.second;
    yogaMemProbe("after createAndBalanceFragments", mp.Rank());

    FragmentDonorFinder donor_finder(frags_from_ranks,is_in_cell);
    Tracer::traceMemory();

    //addWallDistanceToTransferNodes(mp, view, frags_from_ranks);
    addWallDistanceToTransferNodesWithChunkedSurfaces(mp,view,frags_from_ranks);
    if(mesh_system_info.numberOfComponents() == int(component_grid_importance.size())) {
        modifyDistanceBasedOnComponentImportance(frags_from_ranks, component_grid_importance);
    }
    auto node_keys_for_ranks =
        buildNodeKeysForRanks(mp, frags_from_ranks, affinities, donor_finder);
    auto receptors = performDonorSearchViaSingleExchange(mp,
                                                         inspector,
                                                         frags_from_ranks,
                                                         node_keys_for_ranks,
                                                         donor_finder,
                                                         g2l);
    yogaMemProbe("after donor search", mp.Rank());
    addNodeNeighborsToReceptors(receptors,view,g2l);
    yogaMemProbe("after addNodeNeighbors", mp.Rank());
    node_keys_for_ranks.clear();
    frags_from_ranks.clear();
    donor_finder.clear();
    auto statuses = generateNodeStatuses(mp,partition_info,mesh_system_info,
                                         view,receptors,g2l,extra_layers,should_add_max_receptors);
    yogaMemProbe("after generateNodeStatuses", mp.Rank());
    printStats(view, statuses, rootPrinter, mp);



    Tracer::begin("receptor redux");
    auto reduced_receptors = removeNonReceptors(receptors, statuses, g2l);
    Tracer::end("receptor redux");
    Tracer::traceMemory();

    Tracer::begin("create overset-data");
    auto ptr = std::make_shared<OversetData>(std::move(statuses), std::move(reduced_receptors), std::move(g2l));
    Tracer::end("create overset-data");

    Tracer::end("Domain Assembly");

    auto after_assembly = Parfait::Now();
    std::string timing_msg = "Total assembly time: ";
    timing_msg += std::to_string(Parfait::elapsedTimeInSeconds(before_assembly,after_assembly));
    timing_msg += " seconds\n";
    rootPrinter.print(timing_msg);
    Tracer::traceMemory();

    inspector.collect();
    if(0 == mp.Rank()) {
        printInspectorResults(inspector);
    }

    return ptr;
}
void modifyDistanceBasedOnComponentImportance(FragmentMap& fragments, const std::vector<int>& component_grid_importance) {
    for(auto& pair:fragments){
        auto& fragment = pair.second;
        for(auto& node:fragment.transferNodes){
            auto level = component_grid_importance[node.associatedComponentId];
            if(level >= 1)
                node.distanceToWall /= (level*1.1);
        }
    }
}

std::vector<long> countQueryPointsSentToEachRank(MessagePasser mp,
                                                 std::map<int, std::vector<std::pair<int,int>>>& node_keys_for_ranks){
    std::vector<long> query_point_counts(mp.NumberOfProcesses(),0);
    for(auto& key:node_keys_for_ranks) {
        int target_rank = key.first;
        long number_of_query_points_i_will_send = key.second.size();
        query_point_counts[target_rank] += number_of_query_points_i_will_send;
    }
    mp.ElementalSum(query_point_counts,0);
    mp.Broadcast(query_point_counts,0);
    return query_point_counts;
}

std::map<int, std::vector<std::pair<int,int>>> extractNodeKeysForCommRound(int current_round,
                                                                            int total_rounds,
                                                                            std::map<int, std::vector<std::pair<int,int>>>& all_node_keys_for_ranks){
    std::map<int, std::vector<std::pair<int,int>>> node_key_subset;
    for(auto& pair:all_node_keys_for_ranks){
        int target_rank = pair.first;
        auto& keys = pair.second;
        int query_point_count = keys.size();
        auto range = Parfait::LinearPartitioner::getRangeForWorker(current_round,query_point_count,total_rounds);
        node_key_subset[target_rank] = std::vector<std::pair<int,int>>(keys.begin()+range.start,keys.begin()+range.end);
    }
    return node_key_subset;
}

std::vector<Receptor> performDonorSearchViaSingleExchange(MessagePasser mp,
                                                          Parfait::Inspector& inspector,
                                                          const std::map<int,VoxelFragment>& frags_from_ranks,
                                                          std::map<int, std::vector<std::pair<int, int>>>&  node_keys_for_ranks,
                                                          FragmentDonorFinder& donor_finder,
                                                          std::map<long,int>& g2l
){
    auto query_pts_from_ranks = buildAndExchangeQueryPoints(mp, inspector, frags_from_ranks, node_keys_for_ranks);
    inspector.begin("buildReceptorCollections");
    auto receptor_collections_for_ranks = buildReceptorCollectionsForRanks(mp, donor_finder, query_pts_from_ranks);
    inspector.end("buildReceptorCollections");
    query_pts_from_ranks.clear();
    query_pts_from_ranks.shrink_to_fit();
    return exchangeAndUnpackDonors(mp, g2l, receptor_collections_for_ranks);
}

std::vector<Receptor> performDonorSearchInChunksToReducePeakMemory(MessagePasser mp,
                                                                   Parfait::Inspector& inspector,
                                                                   const std::map<int,VoxelFragment>& frags_from_ranks,
                                                                   std::map<int, std::vector<std::pair<int, int>>>&  node_keys_for_ranks,
                                                                   FragmentDonorFinder& donor_finder,
                                                                   std::map<long,int>& g2l
                                                                   ){
    auto query_point_counts = countQueryPointsSentToEachRank(mp,node_keys_for_ranks);
    inspector.addExternalField("query-points-processed",query_point_counts);
    long max_query_points_per_round = 50000;
    int n_comm_rounds = countRequiredCommRounds(max_query_points_per_round,query_point_counts);

    std::vector<Receptor> all_receptors;
    RootPrinter rp(mp.Rank());
    for(int i=0;i<n_comm_rounds;i++) {
        std::string s = "round "+std::to_string(i)+" of "+std::to_string(n_comm_rounds);
        Tracer::begin(s);
        rp.print("-query point exchange "+s+"\n");
        inspector.begin("extractNodeKeys");
        auto node_keys = extractNodeKeysForCommRound(i,n_comm_rounds,node_keys_for_ranks);
        inspector.end("extractNodeKeys");
        auto query_pts_from_ranks = buildAndExchangeQueryPoints(mp, inspector, frags_from_ranks, node_keys);
        inspector.begin("buildReceptorCollections");
        auto receptor_collections_for_ranks = buildReceptorCollectionsForRanks(mp, donor_finder, query_pts_from_ranks);
        inspector.end("buildReceptorCollections");
        auto receptors = exchangeAndUnpackDonors(mp, g2l, receptor_collections_for_ranks);
        for(auto& receptor:receptors)
            all_receptors.push_back(receptor);
        Tracer::end(s);
    }
    return all_receptors;
}

std::vector<Receptor> exchangeAndUnpackDonors(const MessagePasser& mp,
                                            std::map<long, int>& g2l,
                                            const std::map<int, ReceptorCollection>& receptor_collections_for_ranks) {
    Tracer::begin("exchange donorGids");
    auto collection_packer = [](MessagePasser::Message& msg,ReceptorCollection& collection){
      collection.pack(msg);
    };
    auto collection_unpacker = [](MessagePasser::Message& msg,ReceptorCollection& collection){
      collection.unpack(msg);
    };
    auto receptor_collections_from_ranks = mp.Exchange(receptor_collections_for_ranks,
        collection_packer,collection_unpacker);
    Tracer::end("exchange donorGids");
    Tracer::traceMemory();

    Tracer::begin("unpack donor info");
    std::map<int,Receptor> receptor_map;
    for(auto& pair:receptor_collections_from_ranks){
        auto& collection = pair.second;
        for(size_t i=0;i<collection.size();i++){
            auto r = collection.get(i);
            int local_id = g2l.at(r.globalId);
            if(receptor_map.count(local_id) == 0){
                receptor_map[local_id] = r;
            }
            else{
                auto& donor_list = receptor_map[local_id].candidateDonors;
                donor_list.insert(donor_list.end(),
                    r.candidateDonors.begin(),r.candidateDonors.end());
            }
        }
    }
    Tracer::traceMemory();
    std::vector<Receptor> receptors;
    for(auto& pair:receptor_map){
        receptors.push_back(pair.second);
    }
    Tracer::end("unpack donor info");
    Tracer::traceMemory();
    return receptors;
}
std::map<int, std::vector<std::pair<int, int>>> buildNodeKeysForRanks(const MessagePasser& mp,
                                                                      const FragmentMap& frags_from_ranks,
                                                                      const AffinityMap& affinities,
                                                                      const FragmentDonorFinder& donor_finder) {
    Tracer::begin("gather donor finder extents");
    auto donor_finder_extents = mp.Gather(donor_finder.getExtent());
    Tracer::end("gather donor finder extents");

    Tracer::begin("build overlap detector");
    OverlapDetector overlap_detector(donor_finder_extents);
    Tracer::end("build overlap detector");

    auto node_keys_for_ranks = mapNodeIdsToRanks(mp.Rank(),
        frags_from_ranks,overlap_detector,affinities);
    return node_keys_for_ranks;
}
std::vector<NodeStatus> generateNodeStatuses(MessagePasser mp,
                                             const PartitionInfo& partition_info,
                                             const MeshSystemInfo& mesh_system_info,
                                             const YogaMesh& mesh,
                                             std::vector<Receptor>& receptors,
                                             const std::map<long,int>& g2l,
                                             int extra_layers,
                                                 bool should_add_max_receptors) {
    Tracer::begin("type assignment");
    std::vector<Parfait::Extent<double>> component_grid_extents;
    for(int i=0;i<mesh_system_info.numberOfComponents();i++)
        component_grid_extents.push_back(mesh_system_info.getComponentExtent(i));

    Tracer::begin("build sync pattern");
    auto sync_pattern = GhostSyncPatternBuilder::build(mesh, mp);
    Tracer::end("build sync pattern");
    Tracer::traceMemory();
    mp.Barrier();

    auto config = YogaConfiguration(mp);
    auto hole_maps = createHoleMaps(mp,
                                    mesh,
                                    partition_info,
                                    mesh_system_info,
                                    config.maxHoleMapCells());
    yogaMemProbe("after hole maps", mp.Rank());

    yogaMemProbe("before DruyorTypeAssignment", mp.Rank());

    auto statuses = DruyorTypeAssignment::getNodeStatuses(mesh,
                                                          receptors,
                                                          g2l,
                                                          sync_pattern,
                                                          partition_info,
                                                          mesh_system_info,
                                                          component_grid_extents,
                                                          hole_maps,
                                                          extra_layers,
                                                          should_add_max_receptors,
                                                          config.maxHoleMapCells(),
                                                          mp);
    Tracer::end("type assignment");
    Tracer::traceMemory();
    return extractPlainStatuses(statuses);
}
std::map<int, ReceptorCollection> buildReceptorCollectionsForRanks(
    const MessagePasser& mp,
    FragmentDonorFinder& donor_finder,
    const std::vector<std::vector<TransferNode>>& query_pts_from_ranks) {
    std::map<int,ReceptorCollection> receptor_collections_for_ranks;
    for(int rank=0;rank<mp.NumberOfProcesses();rank++){
        auto& query_pts = query_pts_from_ranks[rank];
            auto&& receptors = donor_finder.generateCandidateReceptors(query_pts);
            for (auto& r : receptors) {
                int owner = r.owner;
                auto& collection = receptor_collections_for_ranks[owner];
                collection.insert(r);
            }
    }
    Tracer::traceMemory();
    return receptor_collections_for_ranks;
}
std::vector<std::vector<TransferNode>> buildAndExchangeQueryPoints(
    const MessagePasser& mp,
    Parfait::Inspector& inspector,
    const FragmentMap& frags_from_ranks,
    std::map<int, std::vector<std::pair<int, int>>>& node_keys_for_ranks) {
    std::vector<std::vector<TransferNode>> query_pts_for_ranks(mp.NumberOfProcesses());
    std::set<long> unique_query_points;
    inspector.begin("pack");
    for(auto& pair:node_keys_for_ranks){
        int target_rank = pair.first;
        for(auto node_key:pair.second) {
            auto& tn = frags_from_ranks.at(node_key.first).transferNodes[node_key.second];
            query_pts_for_ranks[target_rank].push_back(tn);
            unique_query_points.insert(tn.globalId);
        }
    }
    inspector.end("pack");
    Tracer::end("pack query points");
    Tracer::traceMemory();

    Tracer::begin("exchange");
    auto query_pts_from_ranks = mp.Exchange(query_pts_for_ranks);
    Tracer::end("exchange");
    Tracer::traceMemory();
    return query_pts_from_ranks;
}

void addWallDistanceToTransferNodes(MessagePasser mp,
                                    const YogaMesh& mesh,
                                    FragmentMap& frags_from_ranks) {
    Tracer::begin("build surfaces");
    Tracer::traceMemory();
    auto surfaces = ParallelSurface::buildSurfaces(mp, mesh);
    Tracer::end("build surfaces");
    Tracer::traceMemory();
    Tracer::begin("distance");
    NanoFlannDistanceCalculator distanceCalculator(surfaces);
    Tracer::traceMemory();
    for (auto& pair : frags_from_ranks) {
        for (auto& node : pair.second.transferNodes) {
            node.distanceToWall = distanceCalculator.calcDistance(
                node.xyz, node.associatedComponentId);
        }
    }
    Tracer::end("distance");
}

void initializeWallDistanceForTransferNodes(FragmentMap& frags_from_ranks){
    double massive = std::sqrt(std::numeric_limits<double>::max());
    for(auto& pair: frags_from_ranks)
        for(auto& node: pair.second.transferNodes)
            node.distanceToWall = massive;
}

void overrideDistanceIfSmallerIsFound(FragmentMap& frags_from_ranks,
                                      std::vector<Parfait::Point<double>>&& chunk,
                                      int n_components,
                                      int current_component){
    std::vector<std::vector<Parfait::Point<double>>> surface_points_for_components(n_components);
    surface_points_for_components[current_component] = chunk;
    NanoFlannDistanceCalculator distance_calculator(surface_points_for_components);
    for(auto& pair: frags_from_ranks){
        for(auto& node: pair.second.transferNodes){
            if(node.associatedComponentId == current_component){
                double d = distance_calculator.calcDistance(node.xyz,current_component);
                if(d < node.distanceToWall)
                    node.distanceToWall = d;
            }
        }
    }

}

void addWallDistanceToTransferNodesWithChunkedSurfaces(MessagePasser mp,
                                    const YogaMesh& mesh,
                                    FragmentMap& frags_from_ranks) {
    initializeWallDistanceForTransferNodes(frags_from_ranks);
    int n = ParallelSurface::countComponents(mp,mesh);
    for(int component=0;component<n;component++){
        std::string s = "component_" + std::to_string(component);
        Tracer::begin(s);
        auto surface_points = ParallelSurface::getLocalSurfacePointsInComponent(mesh,component);
        long total_points = mp.ParallelSum(long(surface_points.size()));
        long max_points_per_chunk = 50000;
        max_points_per_chunk = std::max(max_points_per_chunk, total_points / 5);
        ChunkedPointGatherer gatherer(mp,surface_points,max_points_per_chunk);
        for(int i=0;i<gatherer.chunkCount();i++){
            Tracer::begin("surface-chunk");
            auto chunk = gatherer.extractChunkInParallel(i);
            overrideDistanceIfSmallerIsFound(frags_from_ranks,std::move(chunk),n,component);
            Tracer::end("surface-chunk");
        }
        Tracer::end(s);
    }
}

std::map<int,std::vector<std::pair<int,int>>> mapNodeIdsToRanks(int rank,
                                                                const std::map<int,VoxelFragment>& fragments,
                                                                const OverlapDetector& overlap_detector,
                                                                const std::map<int,std::vector<bool>>& affinities) {

    std::map<int, std::vector<std::pair<int,int>>> node_ids_to_ranks;
    std::vector<int> ranks;
    Tracer::begin("identify");
    for (auto& pair:fragments) {
        int frag_rank = pair.first;
        auto& frag = pair.second;
        for(size_t i=0;i<frag.transferNodes.size();i++) {
            bool is_uniquely_claimed_by_fragment = affinities.at(frag_rank)[i];
            if(not is_uniquely_claimed_by_fragment) continue;
            auto& node = frag.transferNodes[i];
            ranks.clear();
            overlap_detector.getOverlappingRanks(node.xyz, ranks);
            for (int r : ranks) {
                node_ids_to_ranks[r].push_back({frag_rank,i});
            }
        }
    }
    Tracer::end("identify");
    return node_ids_to_ranks;
}


int countRequiredCommRounds(long max_per_round,const std::vector<long>& query_point_counts){
    int max_rounds = 0;
    for(long n:query_point_counts){
        int rounds_needed_for_rank = n/max_per_round + 1;
        max_rounds = std::max(rounds_needed_for_rank, max_rounds);
    }
    return max_rounds;
}

}
