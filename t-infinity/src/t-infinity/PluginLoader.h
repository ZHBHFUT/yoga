#pragma once
#include <dlfcn.h>
#include <sys/stat.h>
#include <unistd.h>
#include <memory>
#include <stdexcept>
#include <functional>
#include <string>
#include <vector>
#include <parfait/FileTools.h>
#include "Tracer.h"

namespace inf {
namespace plugins {

    inline std::vector<std::string> possibleExtensions() {
#if __APPLE__
        return {".dylib", ".so"};
#elif defined(_WIN32)
        return {".dll"};
#else
        return {".so"};
#endif
    }

    inline std::vector<std::string> getPossibleLibraryNames(const std::string& base) {
        std::vector<std::string> names;
        for (const auto& e : possibleExtensions()) {
            std::string name = base + e;
            names.push_back(name);
            names.push_back("lib" + name);
        }
        return names;
    }

    inline std::string findLibraryWithBasename(const std::string& path, const std::string base) {
        for (auto& name : getPossibleLibraryNames(base))
            if (Parfait::FileTools::doesFileExist(path + "/" + name)) return name;
        return "";
    }

    inline void throwOnFailure(char* err) {
        if (err) throw std::logic_error("Error extracting creator from plugin " + std::string(err));
    }

    inline void throwIfNotFound(std::string plugin) {
        if (not Parfait::FileTools::doesFileExist(plugin))
            throw std::logic_error("Could not find library exists: " + plugin);
    }

    inline bool doesSymbolExistInSharedObject(void* shared_object_handle, std::string symbol) {
        void* const result = dlsym(shared_object_handle, symbol.c_str());
        return (not dlerror()) and (result != NULL);
    }

    inline bool doesSymbolExist(std::string plugin, std::string symbol) {
        if (not Parfait::FileTools::doesFileExist(plugin)) return false;
        auto handle = dlopen(plugin.c_str(), RTLD_LAZY);
        return doesSymbolExistInSharedObject(handle, symbol);
    }

    inline void forwardTracer(void* shared_object_handle) {
        if (doesSymbolExistInSharedObject(shared_object_handle, "tracer_set_handle")) {
            void* tracer_handle = nullptr;
            tracer_get_handle(&tracer_handle);
            auto setter = (void (*)(void*))dlsym(shared_object_handle, "tracer_set_handle");
            (*setter)(tracer_handle);
        }
    }

    template <typename T>
    std::function<T> loadSymbol(std::string plugin, std::string symbol) {
        throwIfNotFound(plugin);
        auto handle = dlopen(plugin.c_str(), RTLD_LOCAL | RTLD_LAZY);
        throwOnFailure(dlerror());
        void* const result = dlsym(handle, symbol.c_str());
        throwOnFailure(dlerror());
        forwardTracer(handle);
        return reinterpret_cast<T*>(result);
    }

    template <class T>
    auto loadFromPlugin(std::string path, std::string plugin_name, std::string creator_string) {
        auto plugin = findLibraryWithBasename(path, plugin_name);
        std::string full_name_with_path = path + "/" + plugin;
        std::function<T> creator = loadSymbol<T>(full_name_with_path, creator_string);
        return creator();
    }
}

template <typename T>
struct PluginLoader {
    static auto loadPlugin(std::string path, std::string plugin_name, std::string creator_string) {
        return plugins::loadFromPlugin<std::shared_ptr<T>()>(path, plugin_name, creator_string);
    }
    template <typename... Args>
    static auto loadPluginAndCreate(std::string path,
                                    std::string plugin_name,
                                    std::string creator_string,
                                    Args... args) {
        return plugins::loadFromPlugin<std::shared_ptr<T>()>(path, plugin_name, creator_string)
            ->create(args...);
    }
};
}