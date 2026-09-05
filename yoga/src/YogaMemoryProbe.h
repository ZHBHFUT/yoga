#pragma once
#include <Tracer.h>
#include <fstream>
#include <string>

namespace YOGA {

// Append process working-set (MB) to yoga_mem.log for Windows memory triage.
inline void yogaMemProbe(const char* label, int rank = -1) {
#ifdef _WIN32
    const size_t mb = Tracer::usedMemoryMB();
    std::ofstream f("yoga_mem.log", std::ios::app);
    if (!f) {
        return;
    }
    if (rank >= 0) {
        f << "[rank " << rank << "] ";
    }
    f << label << ": " << mb << " MB";
    if (mb > 2048) {
        f << "  *** OVER 2GB ***";
    }
    f << '\n';
#else
    (void)label;
    (void)rank;
#endif
}

}  // namespace YOGA
