#ifdef _WIN32
#include <cstdint>
using u_int64_t = std::uint64_t;
#endif

#include "Snap.h"
#include "VectorFieldAdapter.h"
#include <MessagePasser/MessagePasser.h>
#include <parfait/StringTools.h>
#include <parfait/FileTools.h>
#include <parfait/Throw.h>
#include <parfait/OmlParser.h>
#include <parfait/LinearPartitioner.h>
#include <t-infinity/FileStreamer.h>
#include <cmath>
#include <string>
#include <Tracer.h>

using namespace inf;
Snap::Snap(MPI_Comm comm) : mp(comm) {
    if (mp.Rank() == 0) {
        if (Parfait::FileTools::doesFileExist("backdoor.poml")) {
            auto settings =
                Parfait::OmlParser::parse(Parfait::FileTools::loadFileToString("backdoor.poml"));
            if (settings.has("snap_chunk_size_mb")) {
                int target_chunk_size = settings.at("snap_chunk_size_mb").asInt();
                chunk_max_size_in_MB = std::max(target_chunk_size, 1);
            }
        }
    }
    mp.Broadcast(chunk_max_size_in_MB, 0);
}

Snap::Snap(MessagePasser m) : Snap(m.getCommunicator()) {}

void addCellTopology(inf::Snap& snap, const inf::MeshInterface& mesh) {
    std::vector<long> local_to_global_cell(mesh.cellCount());
    for (int c = 0; c < mesh.cellCount(); c++) {
        local_to_global_cell[c] = mesh.globalCellId(c);
    }
    std::vector<bool> do_own(mesh.cellCount(), true);
    for (int c = 0; c < mesh.cellCount(); c++) {
        do_own[c] = mesh.ownedCell(c);
    }
    snap.setTopology(FieldAttributes::Cell(), local_to_global_cell, do_own);
}

void addNodeTopology(inf::Snap& snap, const inf::MeshInterface& mesh) {
    std::vector<long> local_to_global(mesh.nodeCount());
    for (int c = 0; c < mesh.nodeCount(); c++) {
        local_to_global[c] = mesh.globalNodeId(c);
    }
    std::vector<bool> do_own(mesh.nodeCount(), true);
    for (int n = 0; n < mesh.nodeCount(); n++) {
        do_own[n] = mesh.ownedNode(n);
    }
    snap.setTopology(FieldAttributes::Node(), local_to_global, do_own);
}

void Snap::addMeshTopologies(const inf::MeshInterface& mesh) {
    addCellTopology(*this, mesh);
    addNodeTopology(*this, mesh);
}

void Snap::setTopology(std::string association,
                       const std::vector<long>& global_ids,
                       const std::vector<bool>& do_own) {
    topologies.emplace(association, Parfait::Topology(global_ids, do_own));
}
void Snap::load(std::string filename) {
    Tracer::begin(__FUNCTION__);
    auto extension = Parfait::StringTools::getExtension(filename);
    if (extension != "snap") {
        filename += ".snap";
    }
    fields = readFields(filename);
    Tracer::end(__FUNCTION__);
}

void Snap::add(std::shared_ptr<inf::FieldInterface> field) { fields[field->name()] = field; }

void Snap::writeFieldHeader(FileStreamer& f,
                            size_t n_global_items,
                            size_t entry_length,
                            std::unordered_map<std::string, std::string> attributes) const {
    Tracer::begin("write field header");
    u_int64_t field_length = sizeof(double) * n_global_items;
    u_int64_t n_global = n_global_items;
    uint64_t entry_length_uint64 = entry_length;
    u_int64_t length_remaining = field_length + sizeof(n_global) + sizeof(entry_length_uint64) +
                                 SnapIOHelpers::calcMapSize(attributes);
    f.write(&length_remaining, sizeof(u_int64_t), 1);
    f.write(&n_global, sizeof(u_int64_t), 1);
    f.write(&entry_length_uint64, sizeof(u_int64_t), 1);
    SnapIOHelpers::writeMapToFile(f, attributes);
    Tracer::end("write field header");
}

std::unordered_map<std::string, std::string> Snap::readFieldHeader(FileStreamer& f,
                                                                   size_t& n_global_items,
                                                                   size_t& entry_length) {
    std::unordered_map<std::string, std::string> attributes;
    if (version_read == 2) {
        u_int64_t n = 0;
        f.read(&n, sizeof(u_int64_t), 1);
        std::string name;
        name.resize(n);
        f.read(&name[0], sizeof(char), name.size());

        u_int64_t length_remaining;
        u_int64_t n_global;
        f.read(&length_remaining, sizeof(u_int64_t), 1);
        f.read(&n_global, sizeof(u_int64_t), 1);
        n_global_items = n_global;
        int association_int;
        f.read(&association_int, sizeof(int), 1);
        std::string association = getAssociationFromInt(association_int);
        attributes[FieldAttributes::Association()] = association;
        attributes[FieldAttributes::DataType()] = FieldAttributes::Float64();
        attributes[FieldAttributes::name()] = name;
        entry_length = 1;  // previous to Snap version 3 all fields must be scalars.
    } else if (version_read == 3) {
        u_int64_t length_remaining;
        u_int64_t n_global;
        u_int64_t entry_length_uint64;
        f.read(&length_remaining, sizeof(u_int64_t), 1);
        f.read(&n_global, sizeof(u_int64_t), 1);
        f.read(&entry_length_uint64, sizeof(u_int64_t), 1);

        entry_length = entry_length_uint64;
        n_global_items = n_global;
        attributes = SnapIOHelpers::readMapFromFile(f);
    }

    return attributes;
}

int getAssocationAsInt(std::string assocaition_string) {
    if (assocaition_string == FieldAttributes::Node()) {
        return -1;
    } else if (assocaition_string == FieldAttributes::Cell()) {
        return -2;
    } else {
        return -9;
    }
}

void Snap::extractFieldDataInRangeByStride(const FieldInterface& field,
                                           long start,
                                           long end,
                                           std::vector<double>& data,
                                           std::vector<long>& gids) const {
    auto association = field.attribute(FieldAttributes::Association());
    auto& top = getTopology(association);
    gids.clear();
    data.clear();

    int stride = field.blockSize();

    std::vector<double> d_temp(stride);
    for (int local = 0; local < int(top.global_ids.size()); local++) {
        auto global = top.global_ids[local];
        if (top.do_own[local] and (global >= start and global < end)) {
            field.value(local, d_temp.data());
            gids.push_back(global);
            for (auto d : d_temp) {
                data.push_back(d);
            }
        }
    }
}

std::map<long, double> Snap::extractScalarFieldAsMapInRange(const FieldInterface& field,
                                                            long start,
                                                            long end) const {
    auto association = field.attribute(FieldAttributes::Association());
    auto& top = getTopology(association);
    std::map<long, double> out;
    for (auto pair : top.global_to_local) {
        long global = pair.first;
        auto local = *pair.second.begin();
        if (top.do_own[local] and (global >= start and global < end)) {
            double d;
            field.value(local, &d);
            out[global] = d;
        }
    }
    return out;
}

void Snap::writeFile(std::string filename) const {
    std::shared_ptr<FileStreamer> f = nullptr;
    auto extension = Parfait::StringTools::getExtension(filename);
    if (extension != "snap") {
        filename += ".snap";
    }

    if (0 == mp.Rank()) {
        f = FileStreamer::create("default");
        if (not f->openWrite(filename)) {
            PARFAIT_THROW("Could not open file for snap reading: " + filename);
        }
        writeHeader(filename, *f);
    }

    for (auto& field : fields) {
        std::string name = field.first;
        writeField(*f, *field.second);
    }

    if (0 == mp.Rank()) f->close();
}

void Snap::writeHeader(const std::string& filename, FileStreamer& f) const {
    Tracer::begin("write header");
    f.write(&latest_version, sizeof(u_int64_t), 1);
    u_int64_t nfields = fields.size();
    f.write(&nfields, sizeof(u_int64_t), 1);
    Tracer::end("write header");
}

std::map<std::string, std::shared_ptr<inf::FieldInterface>> Snap::readFields(
    const std::string& filename) {
    Tracer::begin(__FUNCTION__);
    if (topologies.empty())
        PARFAIT_THROW(
            "Cannot read snap file, you haven't registered any topology (NODE, CELL, FACE, "
            "etc...)");
    std::map<std::string, std::shared_ptr<inf::FieldInterface>> fields_from_file;
    int nfields = 0;
    std::shared_ptr<FileStreamer> f = nullptr;
    if (0 == mp.Rank()) {
        f = FileStreamer::create("default");
        openFileAndReadHeader(filename, nfields, *f);
    }
    mp.Broadcast(nfields, 0);
    for (int i = 0; i < nfields; i++) {
        auto field = readField(*f);
        fields_from_file[field->name()] = field;
    }

    if (mp.Rank() == 0) f->close();
    Tracer::end(__FUNCTION__);
    return fields_from_file;
}

int64_t SnapIOHelpers::calcMapSize(const std::unordered_map<std::string, std::string>& map) {
    uint64_t buffer_size = 0;
    buffer_size += sizeof(int64_t);  // num pairs integer
    for (const auto& pair : map) {
        uint64_t size = pair.first.size() * sizeof(char) + pair.second.size() * sizeof(char) +
                        sizeof(uint64_t) * 2;
        buffer_size += size;
    }
    return buffer_size;
}
std::unordered_map<std::string, std::string> SnapIOHelpers::readMapFromFile(FileStreamer& f) {
    uint64_t count;
    f.read(&count, sizeof(uint64_t), 1);

    std::unordered_map<std::string, std::string> my_map;

    for (uint64_t i = 0; i < count; i++) {
        uint64_t length;
        f.read(&length, sizeof(uint64_t), 1);
        std::vector<char> key(length);
        f.read(key.data(), sizeof(char), length);
        f.read(&length, sizeof(uint64_t), 1);
        std::vector<char> value(length);
        f.read(value.data(), sizeof(char), length);

        std::string s_key(key.begin(), key.end());
        std::string s_value(value.begin(), value.end());

        my_map[s_key] = s_value;
    }
    return my_map;
}

void SnapIOHelpers::writeMapToFile(FileStreamer& f,
                                   const std::unordered_map<std::string, std::string>& map) {
    uint64_t num_pairs = map.size();
    f.write(&num_pairs, sizeof(num_pairs), 1);

    for (const auto& pair : map) {
        std::string key = pair.first;
        std::string value = pair.second;
        uint64_t size = key.size();
        f.write(&size, sizeof(size), 1);
        f.write(key.data(), sizeof(char), size);
        size = value.size();
        f.write(&size, sizeof(size), 1);
        f.write(value.data(), sizeof(char), size);
    }
}

void Snap::openFileAndReadHeader(const std::string& filename, int& nfields, FileStreamer& f) {
    Tracer::begin(__FUNCTION__);
    if (not f.openRead(filename)) {
        PARFAIT_THROW("Couldn't open file: " + filename);
    }
    f.read(&version_read, sizeof(u_int64_t), 1);
    u_int64_t num_fields;
    f.read(&num_fields, sizeof(u_int64_t), 1);
    nfields = num_fields;

    throwIfVersionUnknown();
    Tracer::end(__FUNCTION__);
}

std::map<long, std::set<int>> Snap::buildGlobalToLocals(const std::vector<long>& gids) const {
    std::map<long, std::set<int>> g2l;
    for (size_t i = 0; i < gids.size(); i++) {
        auto global = gids[i];
        auto local = i;
        g2l[global].insert(local);
    }
    return g2l;
}

std::vector<std::string> Snap::availableFields() const {
    std::vector<std::string> available;
    for (const auto& field : fields) available.emplace_back(field.first);
    return available;
}

long Snap::countGlobalItems(const std::vector<long>& global_ids) const {
    long biggest = 0;
    for (long id : global_ids) biggest = std::max(biggest, id);
    biggest = mp.ParallelMax(biggest);
    return biggest + 1;
}

std::shared_ptr<FieldInterface> Snap::retrieve(std::string name) const {
    if (fields.count(name) == 0) {
        PARFAIT_WARNING("Available fields are: " +
                        Parfait::StringTools::join(availableFields(), ", "));
        PARFAIT_THROW("Snap has no field: " + name);
    }
    return fields.at(name);
}
std::shared_ptr<inf::FieldInterface> Snap::readField(FileStreamer& f) {
    Tracer::begin(__FUNCTION__);
    size_t n_global_items;
    size_t entry_length;
    std::unordered_map<std::string, std::string> attributes;
    if (mp.Rank() == 0) {
        attributes = readFieldHeader(f, n_global_items, entry_length);
    }
    mp.Broadcast(n_global_items, 0);
    mp.Broadcast(entry_length, 0);
    mp.Broadcast(attributes, 0);

    auto& top = getTopology(attributes.at(FieldAttributes::Association()));
    std::vector<double> field_data(top.global_ids.size() * entry_length);

    size_t megabyte = 1024 * 1024;
    size_t chunk_size = chunk_max_size_in_MB * megabyte / sizeof(double);
    size_t num_read = 0;
    size_t chunk_number = 1;
    std::vector<double> chunk_buffer;
    while (num_read < n_global_items) {
        Tracer::begin("chunk " + std::to_string(chunk_number));
        long num_to_read = std::min(chunk_size, n_global_items - num_read);
        long start = num_read;
        long end = start + num_to_read;
        if (mp.Rank() == 0) {
            chunk_buffer.resize(num_to_read * entry_length);
            f.read(chunk_buffer.data(), sizeof(double), num_to_read * entry_length);
        }
        mp.Broadcast(chunk_buffer, 0);
        for (auto global : top.global_ids) {
            if (global >= start and global < end) {
                int i = global - start;
                for (auto local : top.global_to_local.at(global)) {
                    for (size_t e = 0; e < entry_length; e++) {
                        field_data.at(local * entry_length + e) =
                            chunk_buffer[i * entry_length + e];
                    }
                }
            }
        }
        num_read += num_to_read;
        chunk_number++;
        Tracer::end("chunk " + std::to_string(chunk_number));
    }

    auto output_field =
        std::make_shared<inf::VectorFieldAdapter>(attributes.at(FieldAttributes::name()),
                                                  attributes.at(FieldAttributes::Association()),
                                                  entry_length,
                                                  field_data);
    output_field->setAdapterAttributes(attributes);
    Tracer::end(__FUNCTION__);
    return output_field;
}

std::vector<long> countGIDsInRange(const std::map<long, std::set<int>>& g2l, long start, long end) {
    std::vector<long> gids;
    for (auto pair : g2l) {
        long global = pair.first;
        if (global >= start and global < end) gids.push_back(global);
    }
    return gids;
}

template <typename T, typename I>
std::vector<T> GatherAndSort(const MessagePasser& mp,
                             const std::vector<T>& vec,
                             int stride,
                             const std::vector<I>& ordinal_ids,
                             int root) {
    if (vec.size() / stride != ordinal_ids.size())
        throw std::domain_error(
            "MessagePasser::GatherAndSort requires equal number of objects and ids");
    Tracer::begin("count items");
    long max_ordinal = 0;
    for (auto& gid : ordinal_ids) {
        max_ordinal = std::max(max_ordinal, gid);
    }
    max_ordinal = mp.ParallelMax(max_ordinal);
    size_t total_objects = max_ordinal + 1;
    Tracer::end("count items");

    auto objects_from_ranks = mp.Gather(vec, root);
    auto gids_from_all_ranks = mp.Gather(ordinal_ids, root);
    if (mp.Rank() == root) {
        Tracer::begin("process");
        std::vector<T> output(total_objects * stride);
        for (size_t r = 0; r < objects_from_ranks.size(); r++) {
            auto& gids_from_rank = gids_from_all_ranks[r];
            auto& objects_from_rank = objects_from_ranks[r];
            for (size_t object_index = 0; object_index < gids_from_rank.size(); object_index++) {
                long gid = gids_from_rank[object_index];
                for (int i = 0; i < stride; i++) {
                    output[stride * gid + i] = objects_from_rank[stride * object_index + i];
                }
            }
        }
        Tracer::end("process");
        return output;
    } else {
        return {};
    }
}

void Snap::writeField(FileStreamer& f, const FieldInterface& field) const {
    Tracer::begin("write field " + field.name());
    std::string association = field.attribute(FieldAttributes::Association());
    std::string name = field.attribute(FieldAttributes::name());

    auto& top = getTopology(association);
    size_t n_global_items = countGlobalItems(top.global_ids);
    if (mp.Rank() == 0) {
        writeFieldHeader(f, n_global_items, field.blockSize(), field.getAllAttributes());
    }

    size_t using_chunk_size = chunk_max_size_in_MB;
    auto memory_remaining = Tracer::availableMemoryMB();
    if (memory_remaining < using_chunk_size) {
        using_chunk_size = memory_remaining / 2;
    }
    using_chunk_size = mp.ParallelMin(using_chunk_size);
    if (using_chunk_size <= 0) {
        PARFAIT_WARNING("Free Memory available (" + std::to_string(using_chunk_size * 2) +
                        "MB) may not be enough to write snap file\n");
    }
    size_t megabyte = 1024 * 1024;
    size_t chunk_size = using_chunk_size * megabyte / (sizeof(double) * field.blockSize());
    size_t num_written = 0;
    size_t chunk_number = 1;
    std::vector<double> chunk_data;
    long start = 0;
    while (num_written < n_global_items) {
        Tracer::begin("chunk " + std::to_string(chunk_number));
        long num_to_write = std::min(chunk_size, n_global_items - start);
        long end = start + num_to_write;
        std::vector<double> chunk_data_i_send;
        std::vector<long> gids_in_range;
        Tracer::begin("extract field data in range");
        extractFieldDataInRangeByStride(field, start, end, chunk_data_i_send, gids_in_range);
        Tracer::end("extract field data in range");
        for (auto& gid : gids_in_range) gid -= start;
        Tracer::begin("GatherAndSort");
        chunk_data = GatherAndSort(mp, chunk_data_i_send, field.blockSize(), gids_in_range, 0);
        Tracer::end("GatherAndSort");
        if (mp.Rank() == 0) {
            Tracer::begin("write chunk");
            f.write(chunk_data.data(), sizeof(double), chunk_data.size());
            Tracer::end("write chunk");
        }
        num_written += num_to_write;
        start = end;
        chunk_number++;
        Tracer::end("chunk " + std::to_string(chunk_number));
    }
    Tracer::end("write field " + field.name());
}

void Snap::setMaxChunkSizeInMB(size_t chunk_size) { chunk_max_size_in_MB = chunk_size; }

const Parfait::Topology& Snap::getTopology(std::string association) const {
    if (topologies.count(association) == 0) {
        PARFAIT_THROW("Attempting to access Topology for Unknown Association: " + association);
    }
    return topologies.at(association);
}

Parfait::Topology& Snap::getTopology(std::string association) {
    if (topologies.count(association) == 0) {
        PARFAIT_THROW("Attempting to access Topology for Unknown Association: " + association);
    }
    return topologies.at(association);
}
void Snap::throwIfVersionUnknown() {
    if (version_read != 2 and version_read != 3) {
        PARFAIT_THROW(
            "Encountered snap file with unsupported version.\nSupported versions 2, 3.  "
            "Encountered version " +
            std::to_string(version_read));
    }
}
Parfait::Topology Snap::copyTopology(std::string association) const {
    if (topologies.count(association) == 0) {
        PARFAIT_THROW("Attempting to access Topology for Unknown Association: " + association);
    }
    return topologies.at(association);
}
bool Snap::has(std::string field_name) const { return fields.count(field_name) == 1; }
std::vector<std::unordered_map<std::string, std::string>> Snap::readSummary(std::string filename) {
    std::vector<std::unordered_map<std::string, std::string>> summary;
    int num_fields;
    if (mp.Rank() == 0) {
        auto f = FileStreamer::create("default");
        openFileAndReadHeader(filename, num_fields, *f);

        size_t n_global_items;
        size_t entry_length;
        for (int i = 0; i < num_fields; i++) {
            auto attributes = readFieldHeader(*f, n_global_items, entry_length);
            attributes["block-size"] = std::to_string(entry_length);
            attributes["size"] = std::to_string(n_global_items);
            summary.push_back(attributes);
            size_t skip_bytes = n_global_items * entry_length * sizeof(double);
            f->skip(skip_bytes);
        }
        f->close();
    }
    mp.Broadcast(num_fields, 0);
    summary.resize(num_fields);
    for (int i = 0; i < num_fields; i++) {
        mp.Broadcast(summary[i], 0);
    }

    return summary;
}
void Snap::setDefaultTopology(std::string filename) {
    long cell_count = -1;
    long node_count = -1;
    int num_fields;
    if (mp.Rank() == 0) {
        auto f = FileStreamer::create("default");
        openFileAndReadHeader(filename, num_fields, *f);

        size_t n_global_items;
        size_t entry_length;
        for (int i = 0; i < num_fields; i++) {
            auto attributes = readFieldHeader(*f, n_global_items, entry_length);
            auto association = attributes.at(FieldAttributes::Association());
            if (association == FieldAttributes::Cell()) cell_count = n_global_items;
            if (association == FieldAttributes::Node()) node_count = n_global_items;

            size_t skip_bytes = n_global_items * entry_length * sizeof(double);
            f->skip(skip_bytes);
        }
        f->close();
    }
    mp.Broadcast(cell_count, 0);
    mp.Broadcast(node_count, 0);
    mp.Broadcast(num_fields, 0);

    if (node_count >= 0) {
        auto range = Parfait::LinearPartitioner::getRangeForWorker(
            mp.Rank(), node_count, mp.NumberOfProcesses());
        std::vector<long> global_ids(range.end - range.start);
        std::iota(global_ids.begin(), global_ids.end(), range.start);
        std::vector<bool> do_own(global_ids.size(), true);
        setTopology(FieldAttributes::Node(), global_ids, do_own);
    }
    if (cell_count >= 0) {
        auto range = Parfait::LinearPartitioner::getRangeForWorker(
            mp.Rank(), cell_count, mp.NumberOfProcesses());
        std::vector<long> global_ids(range.end - range.start);
        std::iota(global_ids.begin(), global_ids.end(), range.start);
        std::vector<bool> do_own(global_ids.size(), true);
        setTopology(FieldAttributes::Cell(), global_ids, do_own);
    }
}
