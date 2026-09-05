#include "DistanceCalculator.h"
#include "BetterDistanceCalculator.h"

#include <Tracer.h>
#include <parfait/DistanceTree.h>
#include <parfait/ExtentWriter.h>
#include <parfait/ParallelExtent.h>
#include <parfait/PointWriter.h>
#include <parfait/STL.h>
#include <parfait/Timing.h>

namespace inf {

void traceMaxMemory(MessagePasser mp) {
    auto current_mem_usage = mp.ParallelMax(Tracer::usedMemoryMB(), 0);
    Tracer::counter("Max Mem (MB)", current_mem_usage);
    Tracer::traceMemory();
}

inline void printMaxMemory(MessagePasser mp, std::string message) {
    auto mem = Tracer::usedMemoryMB();
    auto max = mp.ParallelMax(mem);
    mp_rootprint("%s. Root: %ld MB, Max in Parallel: %ld\n", message.c_str(), mem, max);
}

class DistanceCalculator {
  public:
    template <typename FillPoint, typename FillCell, typename CellType, typename GetMetaData>
    DistanceCalculator(MessagePasser mp,
                       const Parfait::Extent<double>& query_bounds,
                       FillPoint fillPoint,
                       FillCell fillCell,
                       CellType cellType,
                       GetMetaData getMetaData,
                       int cell_count,
                       int max_depth = 10,
                       int max_objects_per_voxel = 20)
        : mp(mp), tree({{1, 1, 1}, {1, 1, 1}}) {
        addFacets(query_bounds, fillPoint, fillCell, cellType, getMetaData, cell_count);

        Tracer::begin("construct tree");
        traceMaxMemory(mp);
        tree = Parfait::DistanceTree(tree_domain);
        traceMaxMemory(mp);
        Tracer::end("construct tree");
        tree.setMaxDepth(max_depth);
        tree.setMaxObjectsPerVoxel(max_objects_per_voxel);

        printMaxMemory(mp, "Create empty tree");
        int c = 0;
        auto a = Parfait::Now();
        for (const auto& f : facets_p1) {
            tree.insert(&f);
            if (c++ == 1000) {
                c = 0;
            }
        }
        printMaxMemory(mp, "Insert done");

        Tracer::begin("finalize tree");
        traceMaxMemory(mp);
        tree.finalize();
        printMaxMemory(mp, "Finalize done");
        traceMaxMemory(mp);
        Tracer::end("finalize tree");
        auto b = Parfait::Now();
        auto s = Parfait::readableElapsedTimeAsString(a, b);
        if (mp.Rank() == 0) {
            printf("building time %s\n", s.c_str());
            tree.printVoxelStatistics();
        }
    }

    inline Parfait::Point<double> closest(const Parfait::Point<double>& p) const {
        return tree.closestPoint(p);
    }

    double distance(const Parfait::Point<double>& p) const { return (closest(p) - p).magnitude(); }

    std::tuple<Parfait::Point<double>, Distance::MetaData> closestProjectionWithMetaData(
        const Parfait::Point<double>& p) const {
        Parfait::Point<double> projected_point;
        int local_facet_index;
        std::tie(projected_point, local_facet_index) = tree.closestPointAndIndex(p);

        return {projected_point, meta_data[local_facet_index]};
    }

    std::vector<Parfait::Facet> getFacetsInsideVoxel(size_t id) const {
        std::vector<Parfait::Facet> facets;
        auto& inside = tree.voxels[id].inside_objects;
        for (int facet_index : inside) {
            facets.push_back(facets_p1[facet_index]);
        }
        return facets;
    }

  public:
    MessagePasser mp;
    Parfait::DistanceTree tree;
    Parfait::Extent<double> tree_domain;
    std::vector<Parfait::FacetSegment> facets_p1;
    std::vector<Distance::MetaData> meta_data;
    std::vector<Parfait::TriP2Segment> facets_p2;

    template <typename FillPoint, typename FillCell, typename CellType, typename GetMetaData>
    void addFacets(const Parfait::Extent<double>& query_bounds,
                   FillPoint,
                   FillCell,
                   CellType,
                   GetMetaData,
                   int);
};

template <typename Facet, typename FillPoint, typename FillCell, typename GetMetaData>
void appendTriFacet(std::vector<Facet>& facets,
                    std::vector<Distance::MetaData>& metadata,
                    FillPoint fillPoint,
                    FillCell fillCell,
                    GetMetaData getMetaData,
                    int cell_id) {
    std::array<int, 3> cell;
    fillCell(cell_id, cell.data());
    Parfait::Facet f;
    int i = 0;
    for (int n : cell) {
        fillPoint(n, f[i++].data());
    }

    facets.push_back(f);
    auto m = getMetaData(cell_id);
    metadata.push_back(m);
}
template <typename Facet, typename FillPoint, typename FillCell, typename GetMetaData>
void appendQuadFacets(std::vector<Facet>& facets,
                      std::vector<Distance::MetaData>& metadata,
                      FillPoint fillPoint,
                      FillCell fillCell,
                      GetMetaData getMetaData,
                      int cell_id) {
    std::array<int, 4> cell;
    fillCell(cell_id, cell.data());

    std::array<Parfait::Point<double>, 4> p;
    int i = 0;
    for (int n : cell) {
        fillPoint(n, p[i++].data());
    }

    facets.emplace_back(p[0], p[1], p[2]);
    facets.emplace_back(p[2], p[3], p[0]);
    auto m = getMetaData(cell_id);
    metadata.push_back(m);
    metadata.push_back(m);
}

template <typename Facet, typename FillPoint, typename FillCell, typename GetMetaData>
void appendTri6Facets(std::vector<Facet>& facets,
                      std::vector<Distance::MetaData>& metadata,
                      FillPoint fillPoint,
                      FillCell fillCell,
                      GetMetaData getMetaData,
                      int cell_id) {
    std::array<int, 6> cell;
    fillCell(cell_id, cell.data());

    std::array<Parfait::Point<double>, 6> p;
    int i = 0;
    for (int n : cell) {
        fillPoint(n, p[i++].data());
    }

    facets.emplace_back(p[0], p[3], p[5]);
    facets.emplace_back(p[3], p[1], p[4]);
    facets.emplace_back(p[3], p[4], p[5]);
    facets.emplace_back(p[4], p[2], p[5]);

    auto m = getMetaData(cell_id);
    metadata.push_back(m);
    metadata.push_back(m);
    metadata.push_back(m);
    metadata.push_back(m);
}

template <typename FillPoint, typename FillCell, typename CellType, typename GetMetaData>
void DistanceCalculator::addFacets(const Parfait::Extent<double>& query_bounds,
                                   FillPoint fillPoint,
                                   FillCell fillCell,
                                   CellType cellType,
                                   GetMetaData getMetaData,
                                   int cell_count) {
    std::vector<Parfait::Facet> send_facets_p1;
    std::vector<Parfait::Facet> recv_facets_p1;
    std::vector<Distance::MetaData> send_metadata;
    std::vector<Distance::MetaData> recv_metadata;

    Tracer::begin("add Facets");
    traceMaxMemory(mp);
    printMaxMemory(mp, "begin original add facets");

    bool used_p2_tri = false;
    for (int cell_id = 0; cell_id < cell_count; cell_id++) {
        auto type = cellType(cell_id);

        if (type == inf::MeshInterface::TRI_3)
            appendTriFacet(
                send_facets_p1, send_metadata, fillPoint, fillCell, getMetaData, cell_id);
        if (type == inf::MeshInterface::QUAD_4)
            appendQuadFacets(
                send_facets_p1, send_metadata, fillPoint, fillCell, getMetaData, cell_id);
        if (type == inf::MeshInterface::TRI_6) {
            appendTri6Facets(
                send_facets_p1, send_metadata, fillPoint, fillCell, getMetaData, cell_id);
            used_p2_tri = true;
        }
    }
    used_p2_tri = bool(mp.ParallelSum(int(used_p2_tri)));
    if (used_p2_tri) {
        if (mp.Rank() == 0) {
            PARFAIT_WARNING(
                "Approximated P2 Triangle with 4 linear triangles for distance search.");
        }
    }
    Tracer::begin("Gathering p1 facets");
    traceMaxMemory(mp);
    mp.Gather(send_facets_p1, recv_facets_p1);
    mp.Gather(send_metadata, recv_metadata);
    traceMaxMemory(mp);
    printMaxMemory(mp, "Gathered facets");
    send_facets_p1.clear();
    send_facets_p1.shrink_to_fit();
    printMaxMemory(mp, "Remove duplicate facets");
    traceMaxMemory(mp);
    Tracer::end("Gathering p1 facets");

    if (recv_facets_p1.empty()) {
        PARFAIT_WARNING(
            "Rank %d gathered no surface facets.  This may indicate there were no surface elements "
            "with matching tags");
        return;
    }

    Tracer::begin("culling far facets");
    printMaxMemory(mp, "Culling 0");
    facets_p1.reserve(recv_facets_p1.size());
    tree_domain = recv_facets_p1[0].getExtent();
    facets_p1.emplace_back(recv_facets_p1[0]);
    double closest_possible_furthest_distance = std::numeric_limits<double>::max();
    for (const auto& f : recv_facets_p1) {
        auto fe = f.getExtent();
        double d = furthestDistanceBetween(query_bounds, fe);
        closest_possible_furthest_distance = std::min(d, closest_possible_furthest_distance);
    }
    traceMaxMemory(mp);
    printMaxMemory(mp, "Culling 1");
    auto furthest_possible = closest_possible_furthest_distance;
    for (size_t i = 0; i < recv_facets_p1.size(); i++) {
        const auto& f = recv_facets_p1[i];
        auto fe = f.getExtent();
        auto closest = closestDistanceBetween(query_bounds, fe);
        if (closest < furthest_possible) {
            Parfait::ExtentBuilder::add(tree_domain, fe);
            facets_p1.emplace_back(f);
            meta_data.emplace_back(recv_metadata[i]);
        }
    }
    traceMaxMemory(mp);
    printMaxMemory(mp, "Culling 2");
    Tracer::end("culling far facets");
    Tracer::end("add Facets");
}

DistanceCalculator constructFromInfinityMesh(MessagePasser mp,
                                             const inf::MeshInterface& mesh,
                                             std::set<int> tags,
                                             int max_depth,
                                             int max_objects_per_voxel) {
    tags = mp.ParallelUnion(tags);

    std::vector<int> compact_to_cell_id;
    for (int cell_id = 0; cell_id < mesh.cellCount(); cell_id++) {
        if (mesh.cellOwner(cell_id) == mp.Rank()) {
            auto type = mesh.cellType(cell_id);
            if (inf::MeshInterface::is2DCellType(type)) {
                int tag = mesh.cellTag(cell_id);
                if (tags.count(tag)) {
                    compact_to_cell_id.push_back(cell_id);
                }
            }
        }
    }

    auto fillPoint = [&](int id, double* p) { mesh.nodeCoordinate(id, p); };
    auto getPoint = [&](int id) -> Parfait::Point<double> { return mesh.node(id); };

    auto fillCell = [&](int id, int* face) {
        int cell_id = compact_to_cell_id[id];
        mesh.cell(cell_id, face);
    };

    auto getCellType = [&](int id) {
        int cell_id = compact_to_cell_id[id];
        auto type = mesh.cellType(cell_id);
        switch (type) {
            case inf::MeshInterface::TRI_3:
            case inf::MeshInterface::QUAD_4:
            case inf::MeshInterface::TRI_6:
                return type;
            default:
                PARFAIT_THROW("Distance calculator can only support TRI_3, TRI_6, and QUAD_4");
        }
    };

    auto getMetaData = [&](int id) {
        Distance::MetaData m{mesh.globalCellId(id), mesh.cellOwner(id)};
        return m;
    };

    auto query_bounds = Parfait::ExtentBuilder::build(getPoint, mesh.nodeCount());

    DistanceCalculator distanceCalculator(mp,
                                          query_bounds,
                                          fillPoint,
                                          fillCell,
                                          getCellType,
                                          getMetaData,
                                          compact_to_cell_id.size(),
                                          max_depth,
                                          max_objects_per_voxel);

    return distanceCalculator;
}

std::vector<std::tuple<Parfait::Point<double>, Distance::MetaData>> calculateDistanceAndMetaData(
    MessagePasser mp,
    const inf::MeshInterface& mesh,
    std::set<int> tags,
    int max_depth,
    int max_objects_per_voxel) {
    auto distanceCalculator =
        constructFromInfinityMesh(mp, mesh, tags, max_depth, max_objects_per_voxel);

    std::vector<std::tuple<Parfait::Point<double>, Distance::MetaData>> out(mesh.nodeCount());
    for (int n = 0; n < mesh.nodeCount(); n++) {
        auto p = mesh.node(n);
        out[n] = distanceCalculator.closestProjectionWithMetaData(p);
    }

    Tracer::traceMemory();

    mp.Barrier();
    return out;
}

std::vector<double> calculateDistance(MessagePasser mp,
                                      const inf::MeshInterface& mesh,
                                      std::set<int> tags,
                                      int max_depth,
                                      int max_objects_per_voxel) {
    printMaxMemory(mp, "starting original distance");
    auto distanceCalculator =
        constructFromInfinityMesh(mp, mesh, tags, max_depth, max_objects_per_voxel);
    printMaxMemory(mp, "Original calculator constructed");

    std::vector<double> distance(mesh.nodeCount(), std::numeric_limits<double>::max());
    for (int n = 0; n < mesh.nodeCount(); n++) {
        auto p = mesh.node(n);
        distance[n] = distanceCalculator.distance(p);
    }

    Tracer::traceMemory();

    mp.Barrier();
    return distance;
}

std::vector<double> calcDistance(MessagePasser mp,
                                 const inf::MeshInterface& mesh,
                                 std::set<int> tags,
                                 const std::vector<Parfait::Point<double>>& points,
                                 int max_depth,
                                 int max_objects_per_voxel) {
    auto closest =
        calcNearestSurfacePoint(mp, mesh, tags, points, max_depth, max_objects_per_voxel);
    std::vector<double> distance(points.size());
    for (size_t i = 0; i < points.size(); i++) distance[i] = (points[i] - closest[i]).magnitude();
    return distance;
}

std::vector<Parfait::Point<double>> calcNearestSurfacePoint(
    MessagePasser mp,
    const MeshInterface& mesh,
    std::set<int> tags,
    const std::vector<Parfait::Point<double>>& points,
    int max_depth,
    int max_objects_per_voxel) {
    auto begin = Parfait::Now();
    auto max_memory = Tracer::usedMemoryMB();
    auto distanceCalculator =
        constructFromInfinityMesh(mp, mesh, tags, max_depth, max_objects_per_voxel);
    max_memory = std::max(Tracer::usedMemoryMB(), max_memory);

    std::vector<Parfait::Point<double>> closest(points.size());
    for (int n = 0; n < int(points.size()); n++) {
        auto p = points[n];
        closest[n] = distanceCalculator.closest(p);
    }

    mp.Barrier();
    max_memory = mp.ParallelMax(max_memory);
    auto end = Parfait::Now();
    mp_rootprint("Finished searching distance.\n");
    mp_rootprint("                 Search Time: %s\n",
                 Parfait::readableElapsedTimeAsString(begin, end).c_str());
    mp_rootprint("    Maximum Core Memory (MB): %s\n",
                 Parfait::bigNumberToStringWithCommas(max_memory).c_str());
    return closest;
}
}
