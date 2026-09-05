#pragma once
#include <parfait/JsonParser.h>
#include <unistd.h>
#include <pwd.h>
#include <parfait/SetTools.h>

namespace inf {

std::string getInfDir() {
    char* settings_dir = getenv("INF_SETTINGS_DIR");
    if (not settings_dir) settings_dir = getenv("HOME");
#if defined(_WIN32)
    if (not settings_dir) settings_dir = getenv("USERPROFILE");
#endif
    if (not settings_dir) settings_dir = getpwuid(getuid())->pw_dir;
    return std::string(settings_dir) + "/.infinity";
}

std::string getDefaultDir() { return getPluginDir() + "/../"; }

std::string getExtensionsSettingsFilename() {
    std::string prefix = DRIVER_PREFIX;
    return prefix + "_extensions.json";
}

std::vector<std::string> getAvailableExtensions(const std::string& prefix) {
    using namespace Parfait::StringTools;
    std::set<std::string> extensions;
    std::string path_to_plugins = getPathToTinfLib();
    for (auto& lib : Parfait::FileTools::listFilesInPath(path_to_plugins)) {
        auto words = Parfait::StringTools::split(lib, "_");
        if (words.size() > 2 and contains(words.front(), prefix) and
            contains(words[1], "SubCommand"))
            extensions.insert(words[2]);
    }
    return {extensions.begin(), extensions.end()};
}

std::set<std::string> getLoadedExtensionsFromFile(std::string file) {
    if (Parfait::FileTools::doesFileExist(file)) {
        auto contents = Parfait::FileTools::loadFileToString(file);
        auto dict = Parfait::JsonParser::parse(contents);
        if (dict.has("loaded")) {
            std::vector<std::string> loaded = dict.at("loaded").asStrings();
            return {loaded.begin(), loaded.end()};
        }
    }
    return {};
}

Parfait::Dictionary getExtensionSettings() {
    Parfait::Dictionary dict;
    std::string user_file = getInfDir() + "/" + getExtensionsSettingsFilename();
    std::string default_file = getDefaultDir() + "/" + getExtensionsSettingsFilename();

    std::set<std::string> loaded = {"core"};
    loaded = Parfait::SetTools::Union(loaded, getLoadedExtensionsFromFile(default_file));
    loaded = Parfait::SetTools::Union(loaded, getLoadedExtensionsFromFile(user_file));
    dict["loaded"] = std::vector<std::string>{loaded.begin(), loaded.end()};
    return dict;
}

std::vector<std::string> getLoadedExtensionsFromFile() {
    auto dict = getExtensionSettings();
    return dict["loaded"].asStrings();
}

}
