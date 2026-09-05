#include "AssemblyViaZMQPostMan.h"
#include <Tracer.h>
#include <parfait/RecursiveBisection.h>
#include <memory>
#include "AsyncWorker.h"
#include "CartesianLoadBalancer.h"
#include "YogaConfiguration.h"
#include "Connectivity.h"
#include "DistanceFieldAdapter.h"
#include "DistributedLoadBalancer.h"
#include "DonorDistributor.h"
#include "DruyorTypeAssignment.h"
#include "GhostSyncPatternBuilder.h"
#include "GlobalToLocal.h"
#include "GridFetcher.h"
#include "HoleCutStatPrinter.h"
#include "HoleCuttingTools.h"
#include <parfait/Inspector.h>
#include "LoadBalancer.h"
#include "MeshSystemInfo.h"
#include "OverDecomposer.h"
#include "OverlapMask.h"
#include "ParallelSurface.h"
#include "PartitionInfo.h"
#include "RootPrinter.h"
#include "WorkVoxelBuilder.h"
#include "InspectorPrinter.h"
#include "OverlapDetector.h"
#include "FragmentBalancer.h"

namespace YOGA{

void registerCallbacks(ZMQPostMan& post_man, ZMQPostMan::MailBox& mailbox, DistributedLoadBalancer& distributed_load_balancer);
Parfait::Extent<double> getExtent(std::vector<VoxelFragment> vector);
std::shared_ptr<OversetData> assemblyViaZMQPostMan(MessagePasser mp,
                                                   YogaMesh& view,
                                                   int load_balancer_algorithm,
                                                   int target_voxel_size,
                                                   int extra_layers,
                                                   int rcb_agglom_ncells,
                                                   bool should_add_max_receptors,
                                                   std::function<bool(double*, int, double*)> is_in_cell) {
    mp.Barrier();
    auto before_assembly = Parfait::Now();
    Tracer::begin("Domain Assembly");
    Tracer::traceMemory();
    RootPrinter rootPrinter(mp.Rank());
    rootPrinter.print("\nYoga: starting domain assembly\n");
    Parfait::Inspector inspector(mp,0);
    Tracer::begin("build surfaces");
    Tracer::traceMemory();
    auto surfaces = ParallelSurface::buildSurfaces(mp, view);
    Tracer::end("build surfaces");
    Tracer::traceMemory();

    Tracer::begin("partition info");
    PartitionInfo partition_info(view, mp.Rank());
    Tracer::end("partition info");
    Tracer::traceMemory();

    Tracer::begin("build mesh system info");
    MeshSystemInfo mesh_system_info(mp, partition_info);
    Tracer::end("build mesh system info");
    Tracer::traceMemory();

    auto g2l = GlobalToLocal::buildMap(view);

    auto fragments_and_affinities = createAndBalanceFragments(mp,view,
                                                              partition_info,
                                                              mesh_system_info,
                                                              g2l,
                                                              rcb_agglom_ncells,
                                                              inspector);
    auto& frags_from_ranks = fragments_and_affinities.first;
    std::vector<VoxelFragment> fragments;
    for(auto& frags : frags_from_ranks)
        fragments.emplace_back(frags.second);

    auto extents_for_ranks = mp.Gather(getExtent(fragments));
    OverlapDetector overlap_detector(extents_for_ranks);

    Tracer::traceMemory();
    auto config = YOGA::YogaConfiguration(mp);
    std::shared_ptr<LoadBalancer> loadBalancer = nullptr;
    Tracer::begin("build load balancer (old)");
    loadBalancer = std::make_shared<CartesianLoadBalancer>(mp, view, mesh_system_info, config.maxHoleMapCells());
    Tracer::end("build load balancer (old)");
    rootPrinter.print("Total work units: "+std::to_string(loadBalancer->getRemainingVoxelCount())+"\n");

    auto hole_maps = createHoleMaps(mp, view, partition_info, mesh_system_info, config.maxHoleMapCells());

    auto initial_work_unit_mask = getInitialWorkUnitMask(mp, loadBalancer);
    auto initial_work_units = getInitialWorkUnits(mp, loadBalancer);


    Tracer::traceMemory();

    std::set<int> storage_message_types {MessageTypes::DciUpdate,
                                         MessageTypes::DciReceipt,
                                         MessageTypes::GridFragment,
                                         MessageTypes::WorkExtentBox,
                                         MessageTypes::WorkUnitsComplete};
    std::set<int> callback_message_types {MessageTypes::GridRequest,MessageTypes::WorkRequest};
    callback_message_types.insert(storage_message_types.begin(),storage_message_types.end());
    ZMQPostMan post_man(mp,callback_message_types);
    ZMQPostMan::MailBox mailbox(storage_message_types);

    std::vector<int> host_ranks = {0};
    for(int i=1;i<mp.NumberOfProcesses()/50;i++) {
        int host_id = i*50;
        host_ranks.push_back(host_id);
    }
    DistributedLoadBalancer distributed_load_balancer(mp,*loadBalancer,host_ranks);
    GridFetcher grid_fetcher(fragments);
    long request_count = 0;

    registerCallbacks(post_man, mailbox, distributed_load_balancer);

    long total_fragments = mp.ParallelSum(fragments.size(),0);
    long ranks_with_no_fragments = 0;
    if(fragments.size()==0) ranks_with_no_fragments = 1;
    ranks_with_no_fragments = mp.ParallelSum(ranks_with_no_fragments,0);
    rootPrinter.print("Total mesh fragments: "+std::to_string(total_fragments)+"\n");
    rootPrinter.print("Ranks with no fragments: "+std::to_string(ranks_with_no_fragments)+"\n");
    auto grid_request_callback = [&](MessagePasser::Message&& msg){
        request_count++;
        if(request_count%1000 == 0) printf("Rank %i: %li grid requests so far...\n",mp.Rank(),request_count);
        int sending_rank = 0;
        msg.unpack(sending_rank);
        auto reply = grid_fetcher.doWork(msg);
        post_man.push(sending_rank,MessageTypes::GridFragment,std::move(reply));
    };

    post_man.registerCallBack(MessageTypes::GridRequest,grid_request_callback);

    post_man.start();

    std::shared_ptr<AsyncWorker> worker = nullptr;

    worker = std::make_shared<AsyncWorker>(mp.Rank(),mp.NumberOfProcesses(),mp,post_man, mailbox,
            surfaces, is_in_cell,
            overlap_detector,
            initial_work_unit_mask,
            initial_work_units,
            host_ranks,
            grid_fetcher);

    Tracer::traceMemory();
    Tracer::begin("donor search");
    worker->workUntilFinished();
    mp.Barrier();
    Tracer::end("donor search");
    Tracer::traceMemory();

    post_man.stop();
    Tracer::traceMemory();
    mp.Barrier();
    Tracer::begin("build sync pattern");
    auto sync_pattern = GhostSyncPatternBuilder::build(view,mp);
    Tracer::end("build sync pattern");
    Tracer::traceMemory();
    mp.Barrier();

    auto vecs = mp.Gather(worker->getVoxelSizes(), 0);
    if(0 == mp.Rank())
        printBinStats(vecs);

    Tracer::begin("type assignment");


    std::vector<Receptor> receptors = worker->receptors;
    std::vector<Parfait::Extent<double>> component_grid_extents;
    for(int i=0;i<mesh_system_info.numberOfComponents();i++)
        component_grid_extents.push_back(mesh_system_info.getComponentExtent(i));
    auto nodeStatuses = DruyorTypeAssignment::getNodeStatuses(view,
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

    Tracer::counter("Memory (MB)", mp.Rank());
    std::vector<NodeStatus> plain_statuses(nodeStatuses.size());
    for(int i=0;i<int(nodeStatuses.size());i++)
        plain_statuses[i] = nodeStatuses[i].value();
    printStats(view, plain_statuses, rootPrinter, mp);

    auto reduced_receptors = removeNonReceptors(receptors, plain_statuses, g2l);

    Tracer::end("Domain Assembly");
    inspector.collect();
    auto fragment_counts = mp.Gather(double(fragments.size()),0);
    auto request_counts = mp.Gather(double(request_count),0);
    inspector.addExternalField("fragments", fragment_counts);
    inspector.addExternalField("grid-requests", request_counts);

    if(mp.Rank() == 0){
        printInspectorResults(inspector);
    }
    if(config.shouldDumpStats() and 0 == mp.Rank()) {
        std::string timing_file = "inspector";
        printf("Dumping timing statistics to: %s.csv\n",timing_file.c_str());
        auto csv = inspector.asCSV();
        FILE *f = fopen(timing_file.c_str(),"w");
        fwrite(csv.data(),1,csv.size(),f);
        fclose(f);
    }

    auto ptr = std::make_shared<OversetData>(std::move(plain_statuses), std::move(receptors), std::move(g2l));
    auto after_assembly = Parfait::Now();
    std::string timing_msg = "Total assembly time: ";
    timing_msg += std::to_string(Parfait::elapsedTimeInSeconds(before_assembly,after_assembly));
    timing_msg += " seconds\n";
    rootPrinter.print(timing_msg);
    return ptr;
}

Parfait::Extent<double> getExtent(std::vector<VoxelFragment> fragments) {
    auto e = Parfait::ExtentBuilder::createEmptyBuildableExtent<double>();
    for(auto&frag:fragments)
        for(auto& node : frag.transferNodes)
            Parfait::ExtentBuilder::add(e,node.xyz);
    return e;
}

void registerCallbacks(ZMQPostMan& post_man, ZMQPostMan::MailBox& mailbox, DistributedLoadBalancer& distributed_load_balancer) {
    auto load_balance_callback = [&](MessagePasser::Message&& msg) {
        MessagePasser::Message reply;
        int sending_rank = 0;
        msg.unpack(sending_rank);
        int n_remaining = distributed_load_balancer.getRemainingVoxelCount();
        reply.pack(n_remaining);
        reply.pack(distributed_load_balancer.getWorkVoxel());
        post_man.push(sending_rank, WorkExtentBox, std::move(reply));
    };

    post_man.registerCallBack(WorkRequest, load_balance_callback);
    post_man.registerCallBack(DciUpdate, [&](MessagePasser::Message&& m) {
        int sending_rank = 0;
        m.unpack(sending_rank);
        mailbox.store(DciUpdate, std::move(m));
        post_man.push(sending_rank, DciReceipt, MessagePasser::Message());
    });
    post_man.registerCallBack(DciReceipt, [&](MessagePasser::Message&& m) { mailbox.store(DciReceipt, std::move(m)); });
    post_man.registerCallBack(GridFragment,
                              [&](MessagePasser::Message&& m) { mailbox.store(GridFragment, std::move(m)); });
    post_man.registerCallBack(WorkExtentBox,
                              [&](MessagePasser::Message&& m) { mailbox.store(WorkExtentBox, std::move(m)); });
    post_man.registerCallBack(WorkUnitsComplete,
                              [&](MessagePasser::Message&& m) { mailbox.store(WorkUnitsComplete, std::move(m)); });
}

void printBinStats(const std::vector<std::vector<int>>& vecs) {
    std::vector<int> bins(7, 0);
    int min = 1e6, max = 0, avg = 0;
    int total_units = 0;
    for (auto& v : vecs) {
        for (int n : v) {
            if (n < 100) {
                bins[0]++;
            } else if (n < 1000) {
                bins[1]++;
            } else if (n < 5000) {
                bins[2]++;
            } else if (n < 10000) {
                bins[3]++;
            } else if (n < 20000) {
                bins[4]++;
            } else if (n < 100000) {
                bins[5]++;
            } else {
                bins[6]++;
            }
            total_units++;
            avg += n;
            min = std::min(n, min);
            max = std::max(n, max);
        }
    }
    avg /= total_units;
    printf("Total work units processed: %i\n", total_units);
    printf("Average nodes/voxel:        %i\n", avg);
    printf("Minimum:                    %i\n", min);
    printf("Maximum:                    %i\n", max);
    printf("0-100   100-1k   1k-5k     5k-10k      10k-20k      20k-100k    100k+\n");
    for (auto x : bins) {
        printf("%i   |   ", x);
    }
    printf("\n");
}
}
