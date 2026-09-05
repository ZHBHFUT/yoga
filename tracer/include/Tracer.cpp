#include "Tracer.h"
#include "TracerImpl.h"
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <vector>
#include <memory>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#endif

#define INSTANCE                                 \
    if (global_tracer_handle == nullptr) return; \
    if (not global_tracer_handle->isEnabled())   \
        return;                                  \
    else                                         \
        (*global_tracer_handle)

std::shared_ptr<TracerImpl> global_tracer_handle = nullptr;

std::shared_ptr<TracerImpl> castHandle(void* handle) {
    auto ptr = *(std::shared_ptr<TracerImpl>*)handle;
    return ptr;
}

int tracer_set_handle(void* handle) {
    auto ptr = castHandle(handle);
    global_tracer_handle = ptr;
    return 0;
}
int tracer_get_handle(void** handle_ptr) {
    *handle_ptr = (void*)&global_tracer_handle;
    return 0;
}
void tracer_print_handle() {
    if (global_tracer_handle == nullptr) {
        printf("Tracer: handle is nullptr in print\n");
    } else {
        printf("Tracer: my tracer instance set.  Writing to: %s\n", global_tracer_handle->file_name.c_str());
    }
}

class JsonItem {
  public:
    void addPair(const char* key, const char* value);
    void addPair(const char* key, const std::string& value);
    void addPair(const std::string& key, const char* value);
    void addPair(const std::string& key, const std::string& value);
    void addPair(const std::string& key, const JsonItem& item);
    template <typename T>
    void addPair(std::string key, T value);
    std::string getString() const { return json + "}"; }

  private:
    bool empty = true;
    std::string json = "{";
    void prependCommaIfNeeded();
};

template <typename T>
void JsonItem::addPair(std::string key, T value) {
    prependCommaIfNeeded();
    json += "\"" + key + "\": " + std::to_string(value);
}

class TimeStampedEvent {
    TimeStampedEvent() : thread_id(Event::getThreadId()), timestamp(Event::getTimeStamp()) {}

    std::string category = "DEFAULT";
    int thread_id = -1;
    std::string phase = "i";
    std::string name = "DEFAULT";
    std::string timestamp = "";
};

uint64_t Event::getThreadId() {
    auto tid = std::this_thread::get_id();
    uint64_t* p = (uint64_t*)&tid;
    return *p;
}

std::string Event::getTimeStamp() {
    auto m = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now().time_since_epoch());
    return std::to_string(m.count());
}

void Tracer::setFileName(std::string name) { fileName() = name; }

std::string& Tracer::fileName() {
    static std::string file_name;
    return file_name;
}
void Tracer::setFast() { INSTANCE.setFast(); }
void Tracer::setDebug() { INSTANCE.setDebug(); }
void Tracer::flush() { INSTANCE.flush(); }
void Tracer::skipClosingBrace() { INSTANCE.skipClosingBrace(); }

int& Tracer::processId() {
    static int p_id = -1;
    return p_id;
}
void Tracer::log(std::string name, std::string category) {
    Event e;
    e.process_id = processId();
    e.category = std::move(category);
    e.name = std::move(name);
    e.phase = "i";
    INSTANCE.log(std::move(e));
}
void Tracer::begin(std::string name, std::string category) { INSTANCE.begin(name, category); }
void Tracer::end(std::string name, std::string category) { INSTANCE.end(name, category); }
void Tracer::setThreadName(std::string name) { INSTANCE.setThreadName(std::move(name), processId(), 0); }

void Tracer::endFlowEvent(int id) { INSTANCE.endFlowEvent(id, processId()); }

void Tracer::beginFlowEvent(int id) { INSTANCE.beginFlowEvent(id, processId(), 0, "derp"); }

void Tracer::counter(std::string name, long n) { counter(name, {{name, n}}); }

void Tracer::counter(std::string name, const std::map<std::string, long>& things) {
    Event e;
    e.process_id = processId();
    e.name = std::move(name);
    e.category = "category";
    e.phase = "C";
    e.args = things;
    INSTANCE.log(std::move(e));
}
size_t Tracer::totalMemoryMB() {
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * (page_size / 1024) / (1024);
}

#ifdef __APPLE__
size_t Tracer::availableMemoryMB() {
    mach_msg_type_number_t count = HOST_VM_INFO_COUNT;
    vm_statistics_data_t vmstat;
    if (KERN_SUCCESS != host_statistics(mach_host_self(), HOST_VM_INFO, (host_info_t)&vmstat, &count)) return 0;
    auto free_pages = vmstat.free_count;
    return free_pages * 4096 / 1024 / 1024;
}
#else
size_t Tracer::availableMemoryMB() {
    long pages = sysconf(_SC_AVPHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * (page_size / 1024) / (1024);
}
#endif

#ifdef _WIN32
size_t Tracer::usedMemoryMB() {
    PROCESS_MEMORY_COUNTERS pmc = {};
    pmc.cb = sizeof(pmc);
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        return static_cast<size_t>(pmc.WorkingSetSize / (1024U * 1024U));
    }
    return 0;
}
#elif defined(__APPLE__)
size_t Tracer::usedMemoryMB() {
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    if (KERN_SUCCESS != task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count)) {
        return -1;
    }
    return t_info.resident_size / (1024 * 1024);
}
#else
size_t parseLine(char* line) {
    int i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i - 3] = '\0';
    i = atoi(line);
    return i;
}

size_t Tracer::usedMemoryMB() {  // Note: this value is in MB!
    FILE* file = fopen("/proc/self/status", "r");
    size_t result = 0;
    if (!file) {
        return 0;
    }
    char line[128];

    while (fgets(line, 128, file) != NULL) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result / 1024;
}
#endif
double Tracer::usedMemoryPercent() {
    auto t = totalMemoryMB();
    auto u = usedMemoryMB();
    auto p = 100 * double(u) / double(t);
    return p;
}
double Tracer::availableMemoryPercent() {
    auto t = totalMemoryMB();
    auto a = availableMemoryMB();
    return 100 * double(a) / double(t);
}
void Tracer::mark(std::string m) {
    auto u = usedMemoryMB();
    log(m);
    counter("Memory (MB)", u);
}
void Tracer::traceMemory() {
    auto t = totalMemoryMB();
    auto u = usedMemoryMB();
    auto a = availableMemoryMB();
    auto p = size_t(100 * double(a) / double(t));
    counter("Rank Used Memory (MB)", u);
    counter("Node Remaining Memory (MB)", a);
    counter("Node Remaining Memory %", p);
}

TracerScopedEvent::TracerScopedEvent(const std::string& event_name) { Tracer::begin(event_name); }
TracerScopedEvent::~TracerScopedEvent() {
    if (not std::uncaught_exception()) Tracer::end(event_name);
}

std::string getExtension(std::string filename) { return filename.substr(filename.find_last_of(".") + 1); }

void reportMultipleInitializations() {
    std::string message = "[Tracer] WARNING: Muliple initialize calls detected.  Repeated calls ignored.";
    printf("%s\n", message.c_str());
    Tracer::log(message, "TRACER");
}

void Tracer::initialize(std::string filename) {
    if (isInitialized()) {
        return reportMultipleInitializations();
    }
    processId() = 0;
    auto ext = getExtension(filename);
    if (ext != "trace") filename += ".trace";
    global_tracer_handle = std::make_shared<TracerImpl>(filename, processId());
}

void Tracer::initialize(int process_id) {
    if (isInitialized()) {
        return reportMultipleInitializations();
    }
    processId() = process_id;
    std::string filename = "trace_" + std::to_string(process_id) + ".trace";
    global_tracer_handle = std::make_shared<TracerImpl>(filename, processId());
}
void Tracer::initialize(std::string filename, int process_id, bool enabled) {
    if (isInitialized()) {
        return reportMultipleInitializations();
    }
    processId() = process_id;
    filename += "_" + std::to_string(process_id) + ".trace";
    global_tracer_handle = std::make_shared<TracerImpl>(filename, processId(), enabled);
}
void Tracer::finalize() {
    global_tracer_handle.reset();
    global_tracer_handle = nullptr;
}
bool Tracer::isInitialized() { return global_tracer_handle != nullptr; }

bool Tracer::isEnabled() {
    if (not isInitialized()) return false;
    return global_tracer_handle->isEnabled();
}

TracerImpl::TracerImpl() : file_name("trace.trace") { initialize(); }

TracerImpl::TracerImpl(std::string file_name, int processId, bool enabled)
    : file_name(file_name), process_id(processId), is_enabled(enabled) {
    if (not is_enabled) return;
    initialize();
}

uint64_t TracerImpl::timestamp() {
    auto m = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now().time_since_epoch() - clock_start);
    return m.count();
}
TracerImpl::~TracerImpl() {
    if (not is_enabled) return;
    lock.lock();
    if (not skip_closing_brace) printFileFooter();
    file.close();
    lock.unlock();
}
void TracerImpl::setFast() {
    std::lock_guard<std::mutex> guard(lock);
    flush_always = false;
}
void TracerImpl::setDebug() {
    std::lock_guard<std::mutex> guard(lock);
    flush_always = true;
}
void TracerImpl::log(Event&& e) {
    addThreadId(e);
    addTimeStamp(e);
    writeEvent(e);
}

void TracerImpl::beginItem(std::string& output) {
    if (is_first)
        is_first = false;
    else
        output += ",";
    output += "\n";
}

void TracerImpl::writeEvent(const Event& e) {
    JsonItem item;
    item.addPair("cat", e.category);
    item.addPair("pid", e.process_id);
    item.addPair("tid", e.thread_id);
    item.addPair("ts", e.timestamp);
    item.addPair("ph", e.phase);
    item.addPair("name", e.name);
    if (!e.args.empty()) {
        JsonItem args;
        for (auto& arg : e.args) {
            args.addPair(arg.first, arg.second);
        }
        item.addPair("args", args);
    }
    std::string output = "";
    auto s = item.getString();
    std::lock_guard<std::mutex> guard(lock);
    beginItem(output);
    file << output + s;
    if (flush_always) file.flush();
}
void TracerImpl::flush() { file.flush(); }
void TracerImpl::skipClosingBrace() { skip_closing_brace = true; }
void TracerImpl::addThreadId(Event& e) const {
    if (e.thread_id == -1) {
        auto tid = std::this_thread::get_id();
        uint64_t* p = (uint64_t*)&tid;
        e.thread_id = *p;
    }
}
void TracerImpl::addTimeStamp(Event& e) { e.timestamp = timestamp(); }
void TracerImpl::printFileHeader() { file << "["; }
void TracerImpl::printFileFooter() { file << std::endl << "]"; }
void TracerImpl::initialize() {
    lock.lock();
    file.open(file_name);
    printFileHeader();
    lock.unlock();
    clock_start = std::chrono::high_resolution_clock::now().time_since_epoch();
}

void TracerImpl::setThreadName(const std::string& name, int processId, int threadId) {
    std::lock_guard<std::mutex> guard(lock);
    std::string output = "";
    beginItem(output);
    auto tid = std::this_thread::get_id();
    uint64_t* p = (uint64_t*)&tid;
    threadId = *p;

    JsonItem item;
    item.addPair("name", "thread_name");
    item.addPair("ph", "M");
    item.addPair("pid", processId);
    item.addPair("tid", threadId);
    JsonItem threadNamePair;
    threadNamePair.addPair("name", name);
    item.addPair("args", threadNamePair);
    file << (output + item.getString());
}

void TracerImpl::beginFlowEvent(int eventId, int processId, int threadId, const std::string& timestamp) {
    std::lock_guard<std::mutex> guard(lock);
    std::string output = "";
    beginItem(output);
    Event e;
    addThreadId(e);
    addTimeStamp(e);

    JsonItem item;
    item.addPair("cat", "foo");
    item.addPair("name", "flow_event");
    item.addPair("pid", processId);
    item.addPair("tid", e.thread_id);
    item.addPair("ts", e.timestamp);
    item.addPair("id", eventId);
    item.addPair("ph", "s");
    file << (output + item.getString());
}

void TracerImpl::endFlowEvent(int eventId, int processId) {
    std::lock_guard<std::mutex> guard(lock);
    std::string output = "";
    beginItem(output);
    Event e;
    addThreadId(e);
    addTimeStamp(e);
    JsonItem item;
    item.addPair("cat", "foo");
    item.addPair("name", "flow_event");
    item.addPair("pid", processId);
    item.addPair("tid", e.thread_id);
    item.addPair("ts", e.timestamp);
    item.addPair("id", eventId);
    item.addPair("ph", "f");
    item.addPair("bp", "e");
    file << item.getString();
}
void TracerImpl::begin(std::string name, std::string category) {
    Event e;
    e.process_id = process_id;
    e.name = std::move(name);
    e.category = std::move(category);
    e.phase = "B";
    log(std::move(e));
}
void TracerImpl::end(std::string name, std::string category) {
    Event e;
    e.process_id = process_id;
    e.name = std::move(name);
    e.category = std::move(category);
    e.phase = "E";
    log(std::move(e));
}
bool TracerImpl::isEnabled() { return is_enabled; }

void JsonItem::addPair(const std::string& key, const char* value) {
    prependCommaIfNeeded();
    json += "\"" + key + "\": \"" + value + "\"";
}

void JsonItem::prependCommaIfNeeded() {
    if (empty)
        empty = false;
    else
        json += ", ";
}

void JsonItem::addPair(const std::string& key, const std::string& value) {
    prependCommaIfNeeded();
    json += "\"" + key + "\": \"" + value + "\"";
}

void JsonItem::addPair(const std::string& key, const JsonItem& item) {
    prependCommaIfNeeded();
    json += "\"" + key + "\": " + item.getString();
}
void JsonItem::addPair(const char* key, const char* value) {
    prependCommaIfNeeded();
    json += "\"";
    json += key;
    json += "\": \"";
    json += value;
    json += "\"";
}
void JsonItem::addPair(const char* key, const std::string& value) {
    prependCommaIfNeeded();
    json += "\"";
    json += key;
    json += "\": \"";
    json += value;
    json += "\"";
}
