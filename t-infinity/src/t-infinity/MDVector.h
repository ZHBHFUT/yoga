#pragma once
#include <stddef.h>
#include <vector>
#include <array>
#include <numeric>
#include <functional>

namespace inf {

struct LayoutLeft {};
struct LayoutRight {};

template <size_t Dim, typename Layout>
class MDIndex {};

template <size_t Dim>
class MDIndex<Dim, LayoutLeft> {
  public:
    explicit MDIndex(const std::array<int, Dim>& d) : d(d) {}
    template <typename... Args>
    size_t operator()(int i, Args... indicies) const {
        return calcIndex(1, i, indicies...);
    }
    template <typename... Args>
    size_t calcIndex(int rank, size_t index, int i, Args... args) const {
        index += static_cast<size_t>(i) *
                 std::accumulate(d.begin() + rank, d.end(), size_t{1}, std::multiplies<size_t>());
        return calcIndex(++rank, index, args...);
    }
    size_t calcIndex(int rank, size_t index) const { return index; }
    std::array<int, Dim> d;
};

template <size_t Dim>
class MDIndex<Dim, LayoutRight> {
  public:
    explicit MDIndex(const std::array<int, Dim>& d) : d(d) {}
    template <typename... Args>
    size_t operator()(int i, Args... indicies) const {
        size_t index = 0;
        calcIndex(1, index, i, indicies...);
        return index;
    }
    template <typename... Args>
    void calcIndex(int rank, size_t& index, int i, Args... args) const {
        index += static_cast<size_t>(i) *
                 std::accumulate(d.begin() + rank, d.end(), size_t{1}, std::multiplies<size_t>());
        rank += 1;
        calcIndex(rank, index, args...);
    }
    void calcIndex(int rank, size_t& index, int i) const { index += i; }
    std::array<int, Dim> d;
};

template <typename T, size_t Dim, typename Layout = LayoutRight>
class MDVector {
  public:
    explicit MDVector(const std::array<int, Dim>& dimensions)
        : indexing(dimensions), vec(getVecSize()) {}
    MDVector(const std::array<int, Dim>& dimensions, const T& initial_value)
        : indexing(dimensions), vec(getVecSize(), initial_value) {
        static_assert(not std::is_same<T, bool>::value, "Do not use bool with MDVector.");
    }
    MDVector() : MDVector(std::array<int, Dim>{}) {}

    template <typename... Args>
    T& operator()(Args... indicies) {
        return vec[indexing(indicies...)];
    }

    template <typename... Args>
    const T& operator()(Args... indicies) const {
        return vec[indexing(indicies...)];
    }

    template <typename... Args>
    const T& at(Args... indicies) const {
        return vec.at(indexing(indicies...));
    }

    template <typename... Args>
    T& at(Args... indicies) {
        return vec.at(indexing(indicies...));
    }

    int extent(int dimension) const { return indexing.d[dimension]; }
    std::array<int, Dim> dimensions() const { return indexing.d; }

    T* data() { return vec.data(); }
    const T* data() const { return vec.data(); }

    void fill(const T& value) { std::fill(vec.begin(), vec.end(), value); }

    MDIndex<Dim, Layout> indexing;
    std::vector<T> vec;

  private:
    int getVecSize() {
        return std::accumulate(indexing.d.begin(), indexing.d.end(), 1, std::multiplies<>());
    }
};

template <typename T>
using Vector3D = MDVector<T, 3>;

template <typename T>
using Vector4D = MDVector<T, 4>;

template <typename Functor>
void loopMDRange(const std::array<int, 2>& begin, const std::array<int, 2>& end, const Functor& f) {
    for (int i = begin[0]; i < end[0]; ++i) {
        for (int j = begin[1]; j < end[1]; ++j) {
            f(i, j);
        }
    }
}
template <typename Functor>
void loopMDRange(const std::array<int, 3>& begin, const std::array<int, 3>& end, const Functor& f) {
    for (int i = begin[0]; i < end[0]; ++i) {
        for (int j = begin[1]; j < end[1]; ++j) {
            for (int k = begin[2]; k < end[2]; ++k) {
                f(i, j, k);
            }
        }
    }
}
template <typename Functor>
void loopMDRange(const std::array<int, 4>& begin, const std::array<int, 4>& end, const Functor& f) {
    for (int i = begin[0]; i < end[0]; ++i) {
        for (int j = begin[1]; j < end[1]; ++j) {
            for (int k = begin[2]; k < end[2]; ++k) {
                for (int l = begin[3]; l < end[3]; ++l) {
                    f(i, j, k, l);
                }
            }
        }
    }
}

template <typename T, size_t rank>
auto padMDVector(const MDVector<T, rank>& v, int padding, T initial_value = 0.0) {
    auto dimensions = v.dimensions();
    for (auto& d : dimensions) d += 2 * padding;
    MDVector<T, rank> padded(dimensions, initial_value);
    loopMDRange(std::array<int, rank>{}, v.dimensions(), [&](int i, int j, int k) {
        padded(i + padding, j + padding, k + padding) = v(i, j, k);
    });
    return padded;
}

}
