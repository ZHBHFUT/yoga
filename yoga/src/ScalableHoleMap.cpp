#include "ScalableHoleMap.h"

namespace YOGA {
bool ScalableHoleMap::doesOverlapHole(Parfait::Extent<double>& e) const {
    auto slice = block.getRangeOfOverlappingCells(e);
    for (int i = slice.lo[0]; i < slice.hi[0]; i++) {
        for (int j = slice.lo[1]; j < slice.hi[1]; j++) {
            for (int k = slice.lo[2]; k < slice.hi[2]; k++) {
                int id = block.convert_ijk_ToCellId(i, j, k);
                if (id < 0 || static_cast<size_t>(id) >= cell_statuses.size()) {
                    continue;
                }
                if (CartBlockFloodFill::Crossing == cell_statuses[id] or
                    CartBlockFloodFill::InHole == cell_statuses[id]) return true;
            }
        }
    }
    return false;
}

void ScalableHoleMap::blankLocally() {
    for (int i = 0; i < mesh.numberOfBoundaryFaces(); i++)
        if (partition_info.getAssociatedComponentIdForFace(i) == associatedComponentId)
            if (Solid == mesh.getBoundaryCondition(i)) blankWithExtent(partition_info.getExtentForFace(i));
}

void ScalableHoleMap::blankWithExtent(const Parfait::Extent<double>& e) {
    auto ids = block.getCellIdsInExtent(e);
    for (auto& id : ids) cell_statuses[id] = CartBlockFloodFill::Crossing;
}

void ScalableHoleMap::sync() {
    cell_statuses = mp.ElementalMax(cell_statuses);
}

void ScalableHoleMap::floodFill() {
    auto seeds = CartBlockFloodFill::identifyAndMarkOuterCells(block, cell_statuses);
    SymmetryFinder symmetryFinder(mp, mesh);
    auto planes = symmetryFinder.findSymmetryPlanes();
    for (auto plane : planes) {
        if (plane.componentId() == associatedComponentId) {
            CartBlockFloodFill::removeSeedsOnPlane(plane, block, seeds, cell_statuses);
            break;
        }
    }
    CartBlockFloodFill::stackFill(seeds, block, cell_statuses);
}

}
