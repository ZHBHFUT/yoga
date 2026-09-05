#pragma once
#include <parfait/FileTools.h>
#include <unistd.h>

namespace Parfait {

class Console {
  public:
    int new_file_descriptor;
    fpos_t pos;

    inline void redirect(std::string filename, std::string write_mode = "w") {
        fflush(stdout);
        fgetpos(stdout, &pos);
        new_file_descriptor = dup(fileno(stdout));
        freopen(filename.c_str(), write_mode.c_str(), stdout);
    }
    inline void reset() {
        fflush(stdout);
        dup2(new_file_descriptor, fileno(stdout));
#ifdef _WIN32
        _close(new_file_descriptor);
#else
        close(new_file_descriptor);
#endif
        clearerr(stdout);
        fsetpos(stdout, &pos);
    }
};

}