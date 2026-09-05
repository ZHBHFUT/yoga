#pragma once
#include <string>
#ifdef _WIN32
namespace Parfait {
inline std::string stackTrace() { return {}; }
}
#else
#include <sstream>
#include <iterator>
#include <cxxabi.h>
#include <execinfo.h>

namespace Parfait {
namespace stack_trace_impl {
    inline std::string demangle(const std::string& s) {
        int status = 0;
        auto buffer = abi::__cxa_demangle(s.data(), NULL, NULL, &status);
        std::string out = s;
        if (0 == status) {
            out = buffer;
        }
        free(buffer);
        return out;
    }

    inline std::vector<std::string> breakSymbolLineIntoWords(const std::string& line) {
        std::vector<std::string> words;
        std::istringstream iss(line);
        std::vector<std::string> results((std::istream_iterator<std::string>(iss)),
                                         std::istream_iterator<std::string>());
        return results;
    }
}
inline std::string stackTrace() {
    std::string trace;
    void* trace_elems[150];
    int trace_elem_count(backtrace(trace_elems, 150));
    char** stack_syms(backtrace_symbols(trace_elems, trace_elem_count));
    for (int i = 0; i < trace_elem_count; ++i) {
        auto words = stack_trace_impl::breakSymbolLineIntoWords(stack_syms[i]);
        for (auto& word : words) {
            trace += stack_trace_impl::demangle(word) + " ";
        }
        trace += "\n";
    }
    free(stack_syms);
    return trace;
}
}
#endif