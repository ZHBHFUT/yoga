#include "CartesianLoadBalancer.h"
#include <Tracer.h>
#include <algorithm>
#include <parfait/CartBlockVisualize.h>
#include <parfait/LinearPartitioner.h>
#include <parfait/RecursiveBisection.h>
#include "DensityEstimator.h"
#include "Sleep.h"

namespace YOGA {

CartesianLoadBalancer::CartesianLoadBalancer(MessagePasser mp, const YogaMesh& mesh, MeshSystemInfo& info,
                                             int max_cart_image_cells)
    : targetNodesPerVoxel(calcTargetNodesPerVoxel(mesh.nodeCount())) {
    const int max_cells = std::max(1, max_cart_image_cells);
    auto overlap_extent = getExtentOfSystem(info);
    Parfait::CartBlock block(overlap_extent, 1, 1, 1);
    std::vector<Parfait::Extent<double>> component_extents;
    for (int i = 0; i < info.numberOfComponents(); i++) {
        component_extents.emplace_back(info.getComponentExtent(i));
    }
    auto nodeCountPerCell = MeshDensityEstimator::tallyNodesContainedByCartCells(mp, mesh, component_extents, block);
    std::vector<std::pair<int, Parfait::Extent<double>>> voxels;
    for (int i = 0; i < block.numberOfCells(); ++i) {
        voxels.push_back(std::make_pair(nodeCountPerCell[i], block.createExtentFromCell(i)));
    }
    auto mesh_image = generateCartBlock(overlap_extent, max_cells);
    auto image_node_counts =
        MeshDensityEstimator::tallyNodesContainedByCartCells(mp, mesh, component_extents, mesh_image);
    // refine(mp, mesh, mesh_image, image_node_counts, info, voxels);
    // refine2(mp, mesh, mesh_image, image_node_counts, info, voxels);
    refine3(mp, mesh, mesh_image, image_node_counts, info, voxels);
}

inline int CartesianLoadBalancer::calcTargetNodesPerVoxel(int localNodeCount) {
    // int totalNodesInSystem = MessagePasser::ParallelSum(localNodeCount);
    // int target = totalNodesInSystem/(MessagePasser::NumberOfProcesses()*2);
    // return std::min(50000,target);
    return 10000;
}

inline Parfait::Extent<double> CartesianLoadBalancer::getWorkVoxel() {
    auto e = workVoxels.top().second;
    workVoxels.pop();
    return e;
}

inline Parfait::Extent<double> CartesianLoadBalancer::getExtentOfSystem(MeshSystemInfo& info) {
    auto e = Parfait::ExtentBuilder::createEmptyBuildableExtent<double>();
    int biggest_grid = 0;
    double biggest_volume = Parfait::Extent<double>::volume(info.getComponentExtent(0));

    for (int i = 0; i < info.numberOfComponents(); i++) {
        auto e2 = info.getComponentExtent(i);
        double v = Parfait::Extent<double>::volume(e2);
        if (v > biggest_volume) {
            biggest_volume = v;
            biggest_grid = i;
        }
    }

    bool has_background_grid = true;
    auto background_extent = info.getComponentExtent(biggest_grid);
    for (int i = 0; i < info.numberOfComponents(); i++) {
        if (biggest_grid == i) continue;
        auto e2 = info.getComponentExtent(i);
        if (not background_extent.intersects(e2.lo)) has_background_grid = false;
        if (not background_extent.intersects(e2.hi)) has_background_grid = false;
    }
    if (not has_background_grid) {
        throw std::logic_error(
            "Yoga: could not identify background grid (the largest component does not contain the others)");
    }
    for (int i = 0; i < info.numberOfComponents(); i++) {
        if (has_background_grid and biggest_grid == i) continue;  // skip background grid
        auto e2 = info.getComponentExtent(i);
        Parfait::ExtentBuilder::add(e, e2);
    }
    return e;
}

void CartesianLoadBalancer::refine(MessagePasser mp,
                                   const YogaMesh& mesh,
                                   const Parfait::CartBlock& mesh_image,
                                   const std::vector<int>& image_counts,
                                   MeshSystemInfo& info,
                                   const std::vector<std::pair<int, Parfait::Extent<double>>>& initial_voxels) {
    std::vector<std::pair<int, Parfait::Extent<double>>> voxel_candidates = initial_voxels;
    if (mp.Rank() == 0) {
        int nbins = 1000;
        std::vector<int> bins(nbins);
        while (voxel_candidates.size() > 0) {
            Tracer::begin("refine");
            Parfait::Extent<double> refineCell;
            auto& candidate = voxel_candidates.back();
            int density = candidate.first;
            int targetChunkCount = std::max(2, density / targetNodesPerVoxel);  // at least cut the voxel in half
            auto refineDimensions = calcRefineDimensions(targetChunkCount, candidate.second);
            refineCell = candidate.second;
            voxel_candidates.pop_back();
            Parfait::CartBlock block(refineCell, refineDimensions[0], refineDimensions[1], refineDimensions[2]);
            Tracer::begin("tally");
            estimateNodesPerCell(mesh, block, mesh_image, image_counts, bins);
            Tracer::end("tally");
            for (int i = 0; i < block.numberOfCells(); ++i) {
                if (bins[i] < targetNodesPerVoxel) {
                    workVoxels.push(std::make_pair(bins[i], block.createExtentFromCell(i)));
                } else {
                    auto e = block.createExtentFromCell(i);
                    voxel_candidates.push_back({bins[i], e});
                }
            }
            Tracer::end("refine");
        }
    }
}

void CartesianLoadBalancer::estimateNodesPerCell(const YogaMesh& mesh,
                                                 const Parfait::CartBlock& block,
                                                 const Parfait::CartBlock& mesh_image,
                                                 const std::vector<int>& image_counts,
                                                 std::vector<int>& bins) {
    std::fill(bins.begin(), bins.end(), 0);
    auto thread_chunk = [&](int start, int end) {
        for (int cell = start; cell < end; cell++) {
            auto cell_extent = block.createExtentFromCell(cell);
            auto crossing = mesh_image.getCellIdsInExtent(cell_extent);
            auto slice = mesh_image.getRangeOfOverlappingCells(cell_extent);
            for (int i = slice.lo[0]; i < slice.hi[0]; i++) {
                for (int j = slice.lo[1]; j < slice.hi[1]; j++) {
                    for (int k = slice.lo[2]; k < slice.hi[2]; k++) {
                        int id = mesh_image.convert_ijk_ToCellId(i, j, k);
                        // auto crossing_extent = mesh_image.createExtentFromCell(id);
                        // auto intersection = cell_extent.intersection(crossing_extent);
                        // double factor = intersection.volume(intersection) / crossing_extent.volume(crossing_extent);
                        bins[cell] += image_counts[id];
                    }
                }
            }
        }
    };
    // int nthreads = 1;
    // std::vector<std::shared_ptr<std::thread>> threads;
    // for(int i=0;i<nthreads;i++){
    //    auto range = Parfait::LinearPartitioner::getRangeForWorker(i,block.numberOfCells(),nthreads);
    //    threads.push_back(std::make_shared<std::thread>([&]() { thread_chunk(range.start, range.end); }));
    //}
    auto range = Parfait::LinearPartitioner::getRangeForWorker(0, block.numberOfCells(), 1);
    thread_chunk(range.start, range.end);
    // for(auto& t:threads){
    //    t->join();
    //}
}

inline std::array<int, 3> CartesianLoadBalancer::calcRefineDimensions(int targetChunkCount,
                                                                      const Parfait::Extent<double>& voxel) {
    int n = targetChunkCount;
    if (n < 8) {
        double dx = voxel.getLength_X();
        double dy = voxel.getLength_Y();
        double dz = voxel.getLength_Z();
        if (dx > dy and dx > dz) return {n, 1, 1};
        if (dy > dx and dy > dz) return {1, n, 1};
        return {1, 1, n};
    } else if (n < 27) {
        return {2, 2, 2};
    } else if (n < 64) {
        return {3, 3, 3};
    } else if (n < 125) {
        return {4, 4, 4};
    } else {
        return {5, 5, 5};
    }
}


int countOwnedNodes(const YogaMesh& mesh, int rank) {
    int sum = 0;
    for (int i = 0; i < mesh.nodeCount(); i++) {
        if (mesh.nodeOwner(i) == rank) {
            sum++;
        }
    }
    return sum;
}

void CartesianLoadBalancer::refine2(MessagePasser mp,
                                    const YogaMesh& mesh,
                                    const Parfait::CartBlock& mesh_image,
                                    const std::vector<int>& image_counts,
                                    MeshSystemInfo& info,
                                    const std::vector<std::pair<int, Parfait::Extent<double>>>& initial_voxels) {
    std::vector<std::pair<int, Parfait::Extent<double>>> voxel_candidates = initial_voxels;
    long total_nodes = mp.ParallelSum(countOwnedNodes(mesh, mp.Rank()), 0);
    int npart = total_nodes / calcTargetNodesPerVoxel(mesh.nodeCount());
    std::vector<int> target_partitions = {npart};
    const int nthread = 1;

    //std::mutex candidate_mtx;
    //std::mutex work_voxel_mtx;
    //std::array<std::mutex,nthread> worker_mutex;
    //std::array<std::condition_variable,nthread> condition_var;
    std::array<int,nthread> thread_part_sizes;
    std::array<Parfait::Extent<double>,nthread> thread_extents;
    //std::array<bool,nthread> ready,processed;

    //std::array<std::thread,nthread> threads;

    if (mp.Rank() == 0) {
        Tracer::begin("refine2");
        while (voxel_candidates.size() > 0) {
            int n_units = std::min(nthread,int(voxel_candidates.size()));
            for(int i=0;i<n_units;i++) {
                auto& candidate = voxel_candidates.back();
                thread_part_sizes[i] = target_partitions.back();
                thread_extents[i] = candidate.second;
                voxel_candidates.pop_back();
                target_partitions.pop_back();
            }
            for(int i=0;i<n_units;i++) {
                // set cond var

                processRefineCell(
                    mesh, mesh_image, image_counts, thread_extents[i],
                    thread_part_sizes[i], voxel_candidates, target_partitions);
            }
        }
        Tracer::end("refine2");
    }
}

void CartesianLoadBalancer::processRefineCell(const YogaMesh& mesh,
    const Parfait::CartBlock& mesh_image,
    const std::vector<int>& image_counts,
    Parfait::Extent<double> refine_cell,
    const int num_partitions,
    std::vector<std::pair<int, Parfait::Extent<double>>>& voxel_candidates,
    std::vector<int>& target_partitions
    ) {
    // get longest dim
    int longest_dim = refine_cell.longestDimension();
    std::array<double, 3> dx;
    dx[0] = refine_cell.getLength_X();
    dx[1] = refine_cell.getLength_Y();
    dx[2] = refine_cell.getLength_Z();
    std::array<int, 3> refine_dimensions{1, 1, 1};
    int nbins = 1000;
    std::vector<int> bins(nbins);
    refine_dimensions[longest_dim] = nbins;
    // create # bins along that axis
    Parfait::CartBlock block(refine_cell, refine_dimensions[0], refine_dimensions[1], refine_dimensions[2]);
    // sum in parallel
    Tracer::begin("count bins-");
    estimateNodesPerCell(mesh, block, mesh_image, image_counts, bins);
    Tracer::end("count bins-");

    // use to find center
    long total = std::accumulate(bins.begin(), bins.end(), 0);
    double cut_ratio = Parfait::calcCutRatio(num_partitions);

    long half = cut_ratio * total;
    long so_far = 0;
    double split = refine_cell.center()[longest_dim];
    for (int i = 0; i < nbins; i++) {
        so_far += bins[i];
        if (so_far > half) {
            split = refine_cell.lo[longest_dim] + i * dx[longest_dim] / double(nbins);
            break;
        }
    }
    // split left/right
    auto left = refine_cell;
    left.hi[longest_dim] = split;
    auto right = refine_cell;
    right.lo[longest_dim] = split;
    // if < target, add to work-voxels, else push back onto candidates
    int nleft = Parfait::numberOfLeftPartitions(num_partitions);
    int nright = Parfait::numberOfRightPartitions(num_partitions);
    if (nleft == 1) {
        workVoxels.push(std::make_pair(half, left));
    } else if (nleft > 1) {
        voxel_candidates.push_back({half, left});
        target_partitions.push_back(nleft);
    }
    if (nright == 1) {
        workVoxels.push({half, right});
    } else if (nright > 1) {
        voxel_candidates.push_back({half, right});
        target_partitions.push_back(nright);
    }
}


void CartesianLoadBalancer::refine3(MessagePasser mp,
                                    const YogaMesh& mesh,
                                    const Parfait::CartBlock& mesh_image,
                                        const std::vector<double>& image_counts,
                                        MeshSystemInfo& info,
                                        const std::vector<std::pair<int, Parfait::Extent<double>>>& initial_voxels) {
        std::vector<std::pair<int, Parfait::Extent<double>>> voxel_candidates = initial_voxels;
        long total_nodes = mp.ParallelSum(countOwnedNodes(mesh,mp.Rank()),0);
        if(mp.Rank() == 0) {
            int npart = total_nodes / calcTargetNodesPerVoxel(mesh.nodeCount());
            std::vector<int> target_partitions = {npart};
            Tracer::begin("build estimator");
            DensityEstimator estimator(mesh_image,npart,image_counts);
            Tracer::end("build estimator");

            for(size_t i=0;i<estimator.extents.size();i++){
                auto e = estimator.extents[i];
                int n = estimator.counts[i];
                workVoxels.push({n,e});
            }
#if 0
            int nbins = 1000;
            std::vector<int> bins(nbins,0);
            Tracer::begin("refine3");
            std::vector<Parfait::Extent<double>> plot_extents;
            std::vector<double> density;
            while (voxel_candidates.size() > 0) {
                auto& candidate = voxel_candidates.back();
                int num_partitions = target_partitions.back();
                auto refineCell = candidate.second;
                voxel_candidates.pop_back();
                target_partitions.pop_back();
                // get longest dim
                int longest_dim = Parfait::longestDimension(refineCell);
                std::array<double,3> dx;
                dx[0] = refineCell.getLength_X();
                dx[1] = refineCell.getLength_Y();
                dx[2] = refineCell.getLength_Z();
                std::array<int,3> refine_dimensions {1,1,1};
                refine_dimensions[longest_dim] = nbins;
                // create # bins along that axis
                Parfait::CartBlock block(refineCell, refine_dimensions[0],refine_dimensions[1],refine_dimensions[2]);
                // sum in parallel
                //Tracer::begin("count bins");
                estimator.estimateContainedNodeCounts(block,bins);
                //for(int i=0;i<nbins;i++){
                //    auto e = block.createExtentFromCell(i);
                //    bins[i] = estimator.estimateContainedNodeCount(e);
                //}
                //Tracer::end("count bins");



                // use to find center
                long total = std::accumulate(bins.begin(),bins.end(),0);
                //printf("Binned total: %li, actual: %i, max_error: %i\n",total,actual_sum,max_error);
                double cut_ratio = Parfait::calcCutRatio(num_partitions);

                long half = cut_ratio*total;
                long so_far = 0;
                double split = refineCell.center()[longest_dim];
                for(int i=0;i<nbins;i++) {
                    so_far += bins[i];
                    if(so_far > half){
                        split = refineCell.lo[longest_dim] + i*dx[longest_dim]/double(nbins);
                        break;
                    }
                }
                // split left/right
                auto left = refineCell;
                left.hi[longest_dim] = split;
                auto right = refineCell;
                right.lo[longest_dim] = split;
                // if < target, add to work-voxels, else push back onto candidates
                int nleft = Parfait::numberOfLeftPartitions(num_partitions);
                int nright = Parfait::numberOfRightPartitions(num_partitions);
                if(nleft == 1){
                    workVoxels.push(std::make_pair(half,left));
                    plot_extents.push_back(left);
                    density.push_back(half);
                }
                else{
                    voxel_candidates.push_back({half,left});
                    target_partitions.push_back(nleft);
                }
                if(nright == 1){
                    workVoxels.push({half, right});
                    plot_extents.push_back(right);
                    density.push_back(half);
                }
                else{
                    voxel_candidates.push_back({half,right});
                    target_partitions.push_back(nright);
                }
            }
            Tracer::end("refine3");
            Parfait::ExtentWriter::write("work-units",plot_extents,density);
#endif
        }
}

}
