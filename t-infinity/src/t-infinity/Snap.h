#pragma once
#include <MessagePasser/MessagePasser.h>
#include <memory>
#include <cstdint>
#include "Communicator.h"
#include "MeshInterface.h"
#include <parfait/Throw.h>
#include <parfait/Topology.h>
#include <t-infinity/FieldInterface.h>
#include "FileStreamer.h"

namespace inf {

class Snap {
  public:
    explicit Snap(MPI_Comm comm);
    explicit Snap(MessagePasser mp);
    void addMeshTopologies(const inf::MeshInterface& mesh);
    void setTopology(std::string association,
                     const std::vector<long>& global_ids,
                     const std::vector<bool>& do_own);
    void load(std::string filename);
    void add(std::shared_ptr<inf::FieldInterface> field);
    std::shared_ptr<FieldInterface> retrieve(std::string) const;
    void writeFile(std::string filename) const;
    std::vector<std::string> availableFields() const;
    void setMaxChunkSizeInMB(size_t chunk_size);
    bool has(std::string field_name) const;
    inline void clear() { fields.clear(); }

    Parfait::Topology copyTopology(std::string association) const;
    void setDefaultTopology(std::string filename);

    inline static std::string getAssociationFromInt(int a) {
        if (a == -1) return inf::FieldAttributes::Node();
        if (a == -2) return inf::FieldAttributes::Cell();
        if (a == -9) return inf::FieldAttributes::Unassigned();
        PARFAIT_THROW("Encountered unsupported association int");
    }

    std::vector<std::unordered_map<std::string, std::string>> readSummary(std::string filename);

  private:
    using Association = std::string;
    MessagePasser mp;
    std::map<Association, Parfait::Topology> topologies;
    std::map<std::string, std::shared_ptr<inf::FieldInterface>> fields;
    size_t chunk_max_size_in_MB = 10;
#ifdef _WIN32
    std::uint64_t version_read = 0;
    std::uint64_t latest_version = 3;
#else
    u_int64_t version_read = 0;
    u_int64_t latest_version = 3;
#endif

    std::map<long, std::set<int>> buildGlobalToLocals(const std::vector<long>& gids) const;
    void writeField(FileStreamer& f, const FieldInterface& field) const;
    std::map<std::string, std::shared_ptr<inf::FieldInterface>> readFields(
        const std::string& filename);
    long countGlobalItems(const std::vector<long>& global_ids) const;
    void writeHeader(const std::string& filename, FileStreamer& f) const;
    void writeFieldHeader(FileStreamer& f,
                          size_t n_global_items,
                          size_t entry_length,
                          std::unordered_map<std::string, std::string> attributes) const;
    std::map<long, double> extractScalarFieldAsMapInRange(const FieldInterface& field,
                                                          long start,
                                                          long end) const;
    std::shared_ptr<inf::FieldInterface> readField(FileStreamer& f);
    void openFileAndReadHeader(const std::string& filename, int& nfields, FileStreamer& f);
    std::unordered_map<std::string, std::string> readFieldHeader(FileStreamer& f,
                                                                 size_t& n_global_items,
                                                                 size_t& entry_length);
    const Parfait::Topology& getTopology(std::string association) const;
    Parfait::Topology& getTopology(std::string association);
    void throwIfVersionUnknown();
    void extractFieldDataInRangeByStride(const FieldInterface& field,
                                         long start,
                                         long end,
                                         std::vector<double>& data,
                                         std::vector<long>& gids) const;
};

namespace SnapIOHelpers {
    void writeMapToFile(FileStreamer& f, const std::unordered_map<std::string, std::string>& map);
    std::unordered_map<std::string, std::string> readMapFromFile(FileStreamer& f);
    int64_t calcMapSize(const std::unordered_map<std::string, std::string>& map);
}
}
