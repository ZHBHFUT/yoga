#pragma once
#include <numeric>
#include <vector>

namespace MessagePasserImpl {

inline size_t findOwner(const std::vector<size_t>& counts, long gid) {
    if (gid < 0) {
        throw std::logic_error("Could not find owner");
    }
    for (size_t i = 0; i < counts.size(); i++) {
        const auto bucket = counts[i];
        if (bucket == 0) {
            continue;
        }
        if (gid < static_cast<long>(bucket)) {
            return i;
        }
        gid -= static_cast<long>(bucket);
    }
    throw std::logic_error("Could not find owner");
}

template <typename T>
std::vector<T> flattenAndDelete(std::vector<std::vector<T>>& things) {
    std::vector<T> out;
    size_t count = 0;
    for (const auto& row : things) {
        count += row.size();
    }
    out.reserve(count);
    for (auto& row : things) {
        for (const auto& t : row) {
            out.push_back(t);
        }
        row.clear();
        row.shrink_to_fit();
    }
    return out;
}

inline std::vector<size_t> roughlyEqualCounts(size_t num_items, size_t num_buckets) {
    size_t equal = num_items / num_buckets;
    size_t remainder = num_items % num_buckets;
    std::vector<size_t> out(num_buckets, equal);
    for (size_t i = 0; i < remainder; i++) {
        out[i]++;
    }
    return out;
}

inline void shiftToRequestedRankRange(std::vector<size_t>& counts, int range_range_start) {
    size_t last = counts.size() - 1;
    size_t last_to_copy = last - static_cast<size_t>(range_range_start);
    size_t num_to_copy = counts.size() - static_cast<size_t>(range_range_start);
    for (size_t i = 0; i < num_to_copy; i++) {
        auto index = last - i;
        counts[index] = counts[last_to_copy - i];
    }
    for (size_t i = 0; i < static_cast<size_t>(range_range_start); i++) counts[i] = 0;
}

}

template <typename T>
std::vector<T> MessagePasser::Balance(std::vector<T>&& input, int rank_range_start, int rank_range_end) const {
    std::vector<T> values = std::move(input);
    int num_target_ranks = rank_range_end - rank_range_start;
    std::vector<size_t> proc_counts = Gather(values.size());
    size_t offset = 0;
    for (int i = 0; i < Rank(); i++) {
        offset += proc_counts[i];
    }
    auto total_count = std::accumulate(proc_counts.begin(), proc_counts.end(), size_t(0));
    proc_counts.clear();
    proc_counts.shrink_to_fit();
    auto equal_counts = MessagePasserImpl::roughlyEqualCounts(total_count, static_cast<size_t>(num_target_ranks));

    std::vector<std::vector<T>> send_to_ranks(NumberOfProcesses());
    for (size_t i = 0; i < values.size(); i++) {
        auto gid = offset + i;
        size_t owner = MessagePasserImpl::findOwner(equal_counts, static_cast<long>(gid));
        owner += static_cast<size_t>(rank_range_start);
        send_to_ranks[owner].push_back(values[i]);
    }
    equal_counts.clear();
    equal_counts.shrink_to_fit();
    values.clear();
    values.shrink_to_fit();

    auto incoming = Exchange(send_to_ranks);
    send_to_ranks.clear();
    return MessagePasserImpl::flattenAndDelete(incoming);
}

template <typename T>
std::vector<T> MessagePasser::Balance(const std::vector<T>& values, int rank_range_start, int rank_range_end) const {
    std::vector<T> values_copy(values);
    return Balance(std::move(values_copy), rank_range_start, rank_range_end);
}

template <typename T>
std::vector<T> MessagePasser::Balance(const std::vector<T>& values) const {
    return Balance(values, 0, NumberOfProcesses());
}
