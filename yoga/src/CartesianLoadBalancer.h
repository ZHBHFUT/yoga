#pragma once
#include <parfait/CartBlock.h>
#include <queue>
#include "LoadBalancer.h"
#include "MeshDensityEstimator.h"
#include "MeshSystemInfo.h"
#include "YogaMesh.h"
#include "CartBlockGenerator.h"

namespace YOGA {
class CartesianLoadBalancer : public LoadBalancer {
  public:
    CartesianLoadBalancer(
        MessagePasser mp, const YogaMesh& mesh, MeshSystemInfo& info, int max_cart_image_cells);
    Parfait::Extent<double> getWorkVoxel();
    int getRemainingVoxelCount() { return workVoxels.size(); }

  private:
    typedef Parfait::Extent<double> Extent;
    typedef std::pair<int, Extent> Pair;
    class Compare {
      public:
        bool operator()(Pair& A, Pair& B) { return A.first < B.first; }
    };
    std::priority_queue<Pair, std::vector<Pair>, Compare> workVoxels;
    int targetNodesPerVoxel;
    int calcTargetNodesPerVoxel(int localNodeCount);
    Parfait::Extent<double> getExtentOfSystem(MeshSystemInfo& info);
    void refine(MessagePasser mp, const YogaMesh& mesh,
            const Parfait::CartBlock& mesh_image,
            const std::vector<int>& image_counts,
            MeshSystemInfo& info,
            const std::vector<std::pair<int,Parfait::Extent<double>>>& initial_voxels);
    void refine2(MessagePasser mp, const YogaMesh& mesh,
                const Parfait::CartBlock& mesh_image,
                const std::vector<int>& image_counts,
                MeshSystemInfo& info,
                const std::vector<std::pair<int,Parfait::Extent<double>>>& initial_voxels);
    void refine3(MessagePasser mp, const YogaMesh& mesh,
                 const Parfait::CartBlock& mesh_image,
                 const std::vector<double>& image_counts,
                 MeshSystemInfo& info,
                 const std::vector<std::pair<int,Parfait::Extent<double>>>& initial_voxels);
    std::array<int,3> calcRefineDimensions(int targetChunkCount, const Parfait::Extent<double>& voxel);
    void estimateNodesPerCell(const YogaMesh& mesh,
        const Parfait::CartBlock& block, const Parfait::CartBlock& mesh_image, const std::vector<int>& image_counts,
        std::vector<int>& bins);

    void processRefineCell(const YogaMesh& mesh,
                           const Parfait::CartBlock& mesh_image,
                           const std::vector<int>& image_counts,
                           Parfait::Extent<double> refine_cell,
                           const int num_partitions,
                           std::vector<std::pair<int, Parfait::Extent<double>>>& voxel_candidates,
                           std::vector<int>& target_partitions);
};
}
