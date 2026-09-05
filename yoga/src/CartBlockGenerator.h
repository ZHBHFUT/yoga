#pragma once
#include <algorithm>
#include <cmath>
#include <parfait/CartBlock.h>

namespace YOGA {
inline Parfait::CartBlock generateCartBlock(const Parfait::Extent<double>& e, int max_cells) {
    if (max_cells < 1) {
        max_cells = 1;
    }
    constexpr double min_len = 1.0e-12;
    const double lx = std::max(e.getLength_X(), min_len);
    const double ly = std::max(e.getLength_Y(), min_len);
    const double lz = std::max(e.getLength_Z(), min_len);
    const double dxdy = lx / ly;
    const double dxdz = lx / lz;
    int nx, ny, nz;
    nx = std::max(1, static_cast<int>(std::cbrt(static_cast<double>(max_cells))));
    ny = std::max(1, static_cast<int>(nx / dxdy));
    nz = std::max(1, static_cast<int>(nx / dxdz));
    while (static_cast<long long>(nx) * ny * nz <= max_cells) {
        nx++;
        ny = std::max(1, static_cast<int>(nx / dxdy));
        nz = std::max(1, static_cast<int>(nx / dxdz));
    }
    nx--;
    ny = std::max(1, static_cast<int>(nx / dxdy));
    nz = std::max(1, static_cast<int>(nx / dxdz));
    // Guard against near-degenerate extents producing billions of cells.
    ny = std::min(ny, max_cells);
    nz = std::min(nz, max_cells);
    nx = std::min(nx, max_cells);
    while (static_cast<long long>(nx) * ny * nz > max_cells) {
        if (nz > 1) {
            nz--;
        } else if (ny > 1) {
            ny--;
        } else if (nx > 1) {
            nx--;
        } else {
            break;
        }
    }
    return {e, nx, ny, nz};
}
}
