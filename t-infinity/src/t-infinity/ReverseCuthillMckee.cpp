#include "ReverseCutthillMckee.h"
#include <numeric>

using namespace inf;

std::vector<int> ReverseCutthillMckee::calcNewIds(const std::vector<std::vector<int>>& adjacency) {
    std::vector<bool> do_own(adjacency.size(), true);
    return calcNewIds(adjacency, do_own);
}
int ReverseCutthillMckee::computeBandwidth(const std::vector<int>& row) {
    if (row.size() == 0)  // this can occur when a cell of one type is not connected to any other
                          // cells of that type
        return 0;
    auto min = *std::min_element(row.begin(), row.end());
    auto max = *std::max_element(row.begin(), row.end());
    return max - min;
}
std::vector<int> ReverseCutthillMckee::calcNewIdsQ(const std::vector<std::vector<int>>& adjacency,
                                                   double p) {
    int num_rows = adjacency.size();

    std::mt19937 gen;
    gen.seed(42);
    std::vector<int> new_ordering(num_rows);
    std::iota(new_ordering.begin(), new_ordering.end(), 0);

    for (int row = 0; row < num_rows;) {
        int num_remaining = num_rows - row;
        int bandwidth = computeBandwidth(adjacency[row]);
        int num_rows_to_shuffle = int(bandwidth / p);
        num_rows_to_shuffle = std::max(num_rows_to_shuffle, 1);
        num_rows_to_shuffle = std::min(num_rows_to_shuffle, num_remaining);

        std::shuffle(
            new_ordering.begin() + row, new_ordering.begin() + row + num_rows_to_shuffle, gen);
        row += num_rows_to_shuffle;
    }
    return new_ordering;
}
std::vector<int> ReverseCutthillMckee::calcNewIds(const std::vector<std::vector<int>>& adjacency,
                                                  const std::vector<bool>& do_own) {
    size_t num_rows = adjacency.size();
    std::set<int> already_assigned_nodes;
    std::set<int> remaining_nodes;
    for (size_t i = 0; i < num_rows; i++) {
        remaining_nodes.insert(i);
    }

    std::vector<int> new_ordering;
    int starting_node = Internal::findHighestValence(adjacency);
    new_ordering.push_back(starting_node);
    already_assigned_nodes.insert(starting_node);
    remaining_nodes.erase(starting_node);

    size_t current_index = 0;

    while (already_assigned_nodes.size() < num_rows) {
        int seed;
        if (current_index == new_ordering.size()) {
            seed = findHighestValence(adjacency, remaining_nodes);
        } else {
            seed = new_ordering[current_index++];
        }
        if (already_assigned_nodes.count(seed) != 1) {
            already_assigned_nodes.insert(seed);
            remaining_nodes.erase(seed);
            new_ordering.push_back(seed);
        }
        auto adj = Internal::getSortedAdjacency(adjacency[seed], adjacency);
        for (int id : adj) {
            if (already_assigned_nodes.count(id) != 1) {
                already_assigned_nodes.insert(id);
                remaining_nodes.erase(id);
                new_ordering.push_back(id);
            }
        }
    }

    if (already_assigned_nodes.size() != adjacency.size()) {
        PARFAIT_THROW("CRM reordering is not the correct size.  Sould be: " +
                      std::to_string(adjacency.size()) + " but is " +
                      std::to_string(already_assigned_nodes.size()));
    }
    new_ordering = putOwnedNodesFirst(new_ordering, do_own);
    std::vector<int> new_ids(num_rows);
    for (size_t i = 0; i < num_rows; i++) new_ids[new_ordering.at(i)] = i;
    return new_ids;
}
std::vector<int> ReverseCutthillMckee::putOwnedNodesFirst(
    const std::vector<int>& old_ids_in_new_ordering, const std::vector<bool>& do_own) {
    std::vector<int> owned_nodes_in_order;
    std::vector<int> ghost_nodes_in_order;
    for (size_t i = 0; i < old_ids_in_new_ordering.size(); i++) {
        int node_id = old_ids_in_new_ordering[i];
        if (do_own[i])
            owned_nodes_in_order.push_back(node_id);
        else
            ghost_nodes_in_order.push_back(node_id);
    }
    owned_nodes_in_order.insert(
        owned_nodes_in_order.end(), ghost_nodes_in_order.begin(), ghost_nodes_in_order.end());
    return owned_nodes_in_order;
}
int ReverseCutthillMckee::findHighestValence(const std::vector<std::vector<int>>& adjacency,
                                             const std::set<int>& possibles) {
    int max_connected_id = -1;
    int max_count = -1;
    for (auto node : possibles) {
        if (int(adjacency[node].size()) > max_count) {
            max_count = adjacency[node].size();
            max_connected_id = node;
        }
    }

    return max_connected_id;
}
