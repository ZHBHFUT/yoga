#pragma once
#include <MessagePasser/MessagePasser.h>
#include <t-infinity/Extract.h>
#include <parfait/Facet.h>
#include <parfait/ToString.h>
#include <parfait/DistanceTree.h>
#include <parfait/Flatten.h>
#include <parfait/LinearPartitioner.h>
#include <parfait/PointWriter.h>
#include <parfait/ExtentWriter.h>
#include <parfait/Timing.h>
#include <Tracer.h>
#include <t-infinity/MeasureLoadBalance.h>

namespace inf {

template <typename Container>
Parfait::Extent<double> findSurroundingDomain(const Container& facets) {
    Parfait::Extent<double> domain = Parfait::ExtentBuilder::createEmptyBuildableExtent<double>();
    for (auto& f : facets) {
        Parfait::ExtentBuilder::add(domain, f.getExtent());
    }
    return domain;
}

template <typename Container, typename MetaData>
std::tuple<std::vector<Parfait::Facet>, std::vector<MetaData>> extractFacetsForThisTree(
    const Container& facets,
    const std::vector<MetaData>& metadata,
    int tree_number,
    int number_of_trees) {
    std::vector<Parfait::Facet> out;
    std::vector<MetaData> metadata_out;

    for (long i = 0; i < long(facets.size()); i++) {
        if (i % number_of_trees == tree_number) {
            out.push_back(Parfait::Facet(facets[i]));
            metadata_out.push_back(metadata[i]);
        }
    }
    return {out, metadata_out};
}

Parfait::Point<double> furthestPossiblePoint();

double furthestDistanceBetween(const Parfait::Extent<double>& e1,
                               const Parfait::Extent<double>& e2);

double closestDistanceBetween(const Parfait::Extent<double>& e1, const Parfait::Extent<double>& e2);

template <typename MetaData>
std::tuple<std::vector<Parfait::FacetSegment>, std::vector<MetaData>>
extractNeededFacetsAsSearchable(const std::vector<Parfait::Facet>& facets,
                                const std::vector<MetaData>& metadata,
                                const Parfait::Extent<double>& query_bounds) {
    double closest_possible_furthest_distance = std::numeric_limits<double>::max();
    for (const auto& f : facets) {
        auto fe = f.getExtent();
        double d = furthestDistanceBetween(query_bounds, fe);
        closest_possible_furthest_distance = std::min(d, closest_possible_furthest_distance);
    }
    auto furthest_possible = closest_possible_furthest_distance;

    std::vector<Parfait::FacetSegment> searchableFacets;
    std::vector<MetaData> metadata_out;
    searchableFacets.reserve(facets.size() / 2);
    metadata_out.reserve(facets.size() / 2);
    for (size_t i = 0; i < facets.size(); i++) {
        const auto& f = facets[i];
        auto fe = f.getExtent();
        auto closest = closestDistanceBetween(query_bounds, fe);
        if (closest < furthest_possible) {
            searchableFacets.push_back(Parfait::FacetSegment(f));
            metadata_out.push_back(metadata[i]);
        }
    }
    return {searchableFacets, metadata_out};
}

namespace dist {
    void traceMemoryParallel(MessagePasser mp);
}
inline std::vector<double> calcBetterDistance2D(
    MessagePasser mp,
    const std::vector<Parfait::LineSegment>& segments_partitioned,
    const std::vector<Parfait::Point<double>>& points,
    int number_of_trees = 5,
    int max_depth = 10,
    int max_facets_per_voxel = 20) {
    std::vector<std::array<Parfait::Point<double>, 2>> segments_as_arrays(
        segments_partitioned.size());
    for (int i = 0; i < int(segments_partitioned.size()); i++) {
        segments_as_arrays[i][0] = segments_partitioned[i].a;
        segments_as_arrays[i][1] = segments_partitioned[i].b;
    }

    segments_as_arrays = Parfait::flatten(mp.Gather(segments_as_arrays));
    std::vector<Parfait::LineSegment> segments(segments_as_arrays.size());
    for (int i = 0; i < int(segments_as_arrays.size()); i++) {
        segments[i].a = segments_as_arrays[i][0];
        segments[i].b = segments_as_arrays[i][1];
    }

    auto far_point = furthestPossiblePoint();
    std::vector<Parfait::Point<double>> projected_point_locations(points.size(), far_point);

    auto domain = Parfait::ExtentBuilder::createEmptyBuildableExtent<double>();
    for (auto s : segments) {
        Parfait::ExtentBuilder::add(domain, s.a);
        Parfait::ExtentBuilder::add(domain, s.b);
    }
    auto query_bounds = Parfait::ExtentBuilder::build(points);
    Parfait::DistanceTree tree(query_bounds);
    tree.setMaxDepth(max_depth);
    tree.setMaxObjectsPerVoxel(max_facets_per_voxel);
    for (size_t i = 0; i < segments.size(); i++) {
        tree.insert(&segments[i]);
    }
    tree.finalize();

    for (size_t i = 0; i < points.size(); i++) {
        int index;
        std::tie(projected_point_locations[i], index) =
            tree.closestPointAndIndex(points[i], projected_point_locations[i]);
    }
    std::vector<double> distance(points.size());
    for (size_t i = 0; i < points.size(); i++) {
        distance[i] = (points[i] - projected_point_locations[i]).magnitude();
    }
    return distance;
}

template <typename MetaData>
std::tuple<std::vector<double>, std::vector<double>, std::vector<MetaData>>
calcBetterDistanceWithMetaData(MessagePasser mp,
                               const std::vector<Parfait::Facet>& facets,
                               const std::vector<Parfait::Point<double>>& points,
                               const std::vector<MetaData>& metadata,
                               int number_of_trees = 5,
                               int max_depth = 10,
                               int max_facets_per_voxel = 20) {
    static_assert(std::is_trivially_copyable<MetaData>::value,
                  "Distance Calculator MetaData must be trivially copyable.");

    Tracer::begin("wall distance search");
    dist::traceMemoryParallel(mp);
    auto begin = Parfait::Now();
    auto starting_mem = Tracer::usedMemoryMB();
    auto max_memory = starting_mem;

    auto far_point = furthestPossiblePoint();
    std::vector<Parfait::Point<double>> projected_point_locations(points.size(), far_point);
    std::vector<MetaData> projected_metadata(points.size(), -99);

    std::vector<double> search_cost(points.size(), 0.0);
    auto query_bounds = Parfait::ExtentBuilder::build(points);
    double my_search_time = 0.0;
    for (int tree_number = 0; tree_number < number_of_trees; tree_number++) {
        Tracer::begin("tree " + std::to_string(tree_number));
        dist::traceMemoryParallel(mp);
        auto begin_tree = Parfait::Now();
        std::vector<Parfait::Facet> facets_send;
        std::vector<MetaData> metadata_send;
        std::tie(facets_send, metadata_send) =
            extractFacetsForThisTree(facets, metadata, tree_number, number_of_trees);
        max_memory = std::max(max_memory, Tracer::usedMemoryMB());

        std::vector<Parfait::Facet> facets_recv;
        std::vector<MetaData> metadata_recv;

        mp.Gather(facets_send, facets_recv);
        mp.Gather(metadata_send, metadata_recv);
        if (facets_recv.size() != metadata_recv.size()) {
            PARFAIT_THROW("Require metadata for each facet.");
        }
        dist::traceMemoryParallel(mp);
        max_memory = std::max(max_memory, Tracer::usedMemoryMB());
        facets_send.clear();
        facets_send.shrink_to_fit();
        std::vector<Parfait::FacetSegment> facets_for_this_tree;
        std::vector<MetaData> metadata_for_this_tree;
        std::tie(facets_for_this_tree, metadata_for_this_tree) =
            extractNeededFacetsAsSearchable(facets_recv, metadata_recv, query_bounds);

        if (facets_for_this_tree.size() != metadata_for_this_tree.size()) {
            PARFAIT_THROW("Require metadata for each facet in tree.");
        }
        max_memory = std::max(max_memory, Tracer::usedMemoryMB());
        facets_recv.clear(), facets_recv.shrink_to_fit();
        metadata_recv.clear(), metadata_for_this_tree.shrink_to_fit();
        max_memory = std::max(max_memory, Tracer::usedMemoryMB());

        auto domain = findSurroundingDomain(facets_for_this_tree);
        Parfait::DistanceTree tree(domain);
        tree.setMaxDepth(max_depth);
        tree.setMaxObjectsPerVoxel(max_facets_per_voxel);
        for (size_t i = 0; i < facets_for_this_tree.size(); i++) {
            tree.insert(&facets_for_this_tree[i]);
        }
        max_memory = std::max(max_memory, Tracer::usedMemoryMB());
        tree.finalize();
        max_memory = std::max(max_memory, Tracer::usedMemoryMB());

        for (size_t i = 0; i < points.size(); i++) {
            auto b = Parfait::Now();
            int index;
            std::tie(projected_point_locations[i], index) =
                tree.closestPointAndIndex(points[i], projected_point_locations[i]);
            if (index >= 0) {
                projected_metadata[i] = metadata_for_this_tree[index];
            }
            auto e = Parfait::Now();
            search_cost[i] += Parfait::elapsedTimeInSeconds(b, e);
        }
        auto end_tree = Parfait::Now();
        my_search_time += Parfait::elapsedTimeInSeconds(begin_tree, end_tree);
        Tracer::end("tree " + std::to_string(tree_number));
    }

    std::vector<double> distance(points.size());
    max_memory = std::max(max_memory, Tracer::usedMemoryMB());
    dist::traceMemoryParallel(mp);
    for (size_t i = 0; i < points.size(); i++) {
        distance[i] = (points[i] - projected_point_locations[i]).magnitude();
    }
    max_memory = mp.ParallelMax(max_memory, 0);
    auto end = Parfait::Now();
    Tracer::end("wall distance search");

    inf::LoadBalance::printStatistics(mp, "Distance Search (s)", my_search_time);

    mp_rootprint("Finished searching distance.\n");
    mp_rootprint("                 Search Time: %s\n",
                 Parfait::readableElapsedTimeAsString(begin, end).c_str());
    mp_rootprint("    Maximum Core Memory (MB): %s\n",
                 Parfait::bigNumberToStringWithCommas(max_memory).c_str());
    return {distance, search_cost, projected_metadata};
}

std::tuple<std::vector<double>, std::vector<double>> calcBetterDistance(
    MessagePasser mp,
    const std::vector<Parfait::Facet>& facets,
    const std::vector<Parfait::Point<double>>& points,
    int number_of_trees = 5,
    int max_depth = 10,
    int max_facets_per_voxel = 20);

std::tuple<std::vector<double>, std::vector<double>, std::vector<int>> calcBetterDistance(
    MessagePasser mp,
    const inf::MeshInterface& mesh,
    std::set<int> tags,
    const std::vector<Parfait::Point<double>>& points,
    int number_of_trees = 5,
    int max_depth = 10,
    int max_facets_per_voxel = 20);

std::tuple<std::vector<double>, std::vector<double>, std::vector<int>> calcDistanceToNodes(
    MessagePasser mp,
    const inf::MeshInterface& mesh,
    const std::set<int>& tags,
    int number_of_trees = 5,
    int max_depth = 10,
    int max_facets_per_voxel = 20);
}
