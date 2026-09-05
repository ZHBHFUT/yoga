#pragma once
#include <string>
#include <mutex>
#include <map>
#include <fstream>
#include <chrono>

class Event {
  public:
    std::string category = "DEFAULT";
    int process_id = 0;
    int thread_id = -1;
    std::string phase = "i";
    std::string name = "DEFAULT";
    uint64_t timestamp;
    std::map<std::string, long> args;
    static uint64_t getThreadId();
    static std::string getTimeStamp();
};

class TracerImpl {
  public:
    TracerImpl(std::string file_name, int processId, bool enabled = true);
    TracerImpl();
    ~TracerImpl();
    void log(Event&& e);

    void begin(std::string name, std::string category);
    void end(std::string name, std::string category);
    void setThreadName(const std::string& name, int processId, int threadId);
    void beginFlowEvent(int eventId, int processId, int threadId, const std::string& timestamp);
    void endFlowEvent(int eventId, int processId);
    void setFast();
    void setDebug();
    void flush();
    void skipClosingBrace();
    bool isEnabled();

  public:
    std::string file_name;
    bool is_first = true;
    bool flush_always = false;
    bool skip_closing_brace = false;
    int process_id = 0;
    std::mutex lock;
    std::ofstream file;
    std::chrono::high_resolution_clock::duration clock_start;
    bool is_enabled = true;

    uint64_t timestamp();
    void printFileHeader();
    void printFileFooter();
    void addThreadId(Event& e) const;
    void addTimeStamp(Event& event);
    void writeEvent(const Event& e);
    void beginItem(std::string& output);
    void initialize();

};

