#include "FileTools.h"
#include <fstream>
#include <string.h>
#include <cstring>
#ifdef _WIN32
#include "getline.h"
#include <sys/stat.h>
#include <io.h>
#else
#include <sys/stat.h>
#ifndef S_IFDIR
#define S_IFDIR _S_IFDIR
#endif
#endif
#include <glob.h>
#include "StringTools.h"
#include "Throw.h"
#include <thread>
std::string Parfait::FileTools::readNextLine(FILE* fp) {
    char* line = NULL;
    size_t len = 0;

    std::string string_line;
    int read = getline(&line, &len, fp);
    bool good = read != -1;
    if (good) {
        string_line = std::string(line);
    }
    free(line);
    return string_line;
}
bool Parfait::FileTools::readNextLine(FILE* fp, std::string& string_line) {
    char* line = NULL;
    size_t len = 0;

    int read = getline(&line, &len, fp);
    bool good = read != -1;
    if (good) {
        string_line = Parfait::StringTools::strip(std::string(line));
    }
    free(line);
    return good;
}
std::string Parfait::FileTools::readLinesInFile(FILE* fp, int n) {
    std::string lines;
    for (int i = 0; i < n; i++) {
        std::string l;
        readNextLine(fp, l);
        lines.append(l);
        lines.append("\n");
    }
    return lines;
}
long Parfait::FileTools::countLinesInFile(std::string filename) {
    long n = 0;
    FILE* f = fopen(filename.c_str(), "r");
    if (f == nullptr) return 0;
    std::string line;
    while (readNextLine(f, line)) {
        n++;
    }
    fclose(f);
    return n;
}
std::string Parfait::FileTools::loadFileToString(std::string filename) {
    std::ifstream f(filename.c_str());
    if (not f.is_open()) {
        PARFAIT_THROW("Could not open file for reading: " + filename);
    }
    auto contents = std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
    f.close();
    return contents;
}
void Parfait::FileTools::writeStringToFile(std::string filename, std::string contents) {
    FILE* fp = fopen(filename.c_str(), "w");
    PARFAIT_ASSERT(fp != nullptr, "Could not open file " + filename);
    fwrite(contents.data(), sizeof(char), contents.size(), fp);
    fclose(fp);
}
void Parfait::FileTools::appendToFile(std::string filename, std::string message) {
    FILE* fp = fopen(filename.c_str(), "a+");
    if (fp == nullptr) {
        PARFAIT_THROW("Could not open file to append: " + filename);
    }

    fprintf(fp, "%s", message.c_str());

    fclose(fp);
}
bool Parfait::FileTools::doesDirectoryExist(std::string directory) {
#ifdef _WIN32
    struct stat info;
    std::memset(&info, 0, sizeof(info));
    if (stat(directory.c_str(), &info) != 0)
        return false;
    else if ((info.st_mode & S_IFDIR) != 0)
#else
    struct stat info;
    if (stat(directory.c_str(), &info) != 0)
        return false;
    else if (info.st_mode & S_IFDIR)
#endif
        return true;
    else
        return false;
}
bool Parfait::FileTools::doesFileExist(std::string filename) {
    bool good = false;
    FILE* fp = fopen(filename.c_str(), "r");
    if (fp == nullptr) {
        good = false;
    } else {
        good = true;
        fclose(fp);
    }
    return good;
}
std::vector<std::string> Parfait::FileTools::listFilesInPath(const std::string& path) {
    glob_t result;
    memset(&result, 0, sizeof(result));
    auto pattern = path + "*";
    int flag = glob(pattern.c_str(), GLOB_TILDE, nullptr, &result);
    std::vector<std::string> filenames;
    if (flag == 0)
        for (size_t i = 0; i < result.gl_pathc; i++) filenames.emplace_back(StringTools::stripPath(result.gl_pathv[i]));
    else
        throw std::logic_error("Parfait::FileTools: failed to glob +" + path);
    globfree(&result);
    return filenames;
}
void Parfait::FileTools::waitForFile(std::string filename, int max_seconds) {
    double slept_for = 0.0;
    while (not doesFileExist(filename) and slept_for < max_seconds) {
        slept_for += 0.1;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}
