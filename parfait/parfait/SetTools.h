#pragma once
#include <set>
#include <algorithm>
#include <iterator>

namespace Parfait {
namespace SetTools {
    template <typename T>
    std::set<T> Union(const std::set<T>& a, std::set<T> b) {
        for (const auto& aa : a) {
            b.insert(aa);
        }
        return b;
    }
    template <typename T>
    std::set<T> Difference(const std::set<T>& a, const std::set<T>& b) {
        std::set<T> result;
        std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::inserter(result, result.end()));
        return result;
    }
    template <typename T>
    std::set<T> Intersection(const std::set<T>& a, const std::set<T>& b) {
        std::set<T> result;
        std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::inserter(result, result.end()));
        return result;
    }

    template <typename T>
    std::set<T> range(T start, T end) {
        std::set<T> my_set;
        for (T t = start; t < end; t++) {
            my_set.insert(t);
        }
        return my_set;
    }

    template <typename T, typename Container>
    std::set<T> toSet(const Container& c) {
        std::set<T> s(c.begin(), c.end());
        return s;
    }
}
}
