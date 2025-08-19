#pragma once

#include <memory>
#include <stdint.h>

#include "../util/require.hpp"

namespace camspork
{

constexpr uint32_t sync_bit = 0x8000'0000;

struct assignment_record_id
{
    uint32_t node_id = 0;
    operator bool() const
    {
        return node_id != 0;
    }
};

struct barrier_id
{
    uint32_t data = 0;
};

struct syncv_init_t
{
    // TODO remove unused stuff (or use it).
    const char* filename;
    uint64_t memory_budget;
    uint64_t debug_assignment_id;
    uint64_t debug_operation_id;
};

struct SyncvTable;

struct SyncvTableDeleter
{
    void operator() (SyncvTable* victim) const;
};

using SyncvTable_unique_ptr = std::unique_ptr<SyncvTable, SyncvTableDeleter>;

struct SigthreadInterval;

struct ThreadCuboid
{
    uint32_t task_index = 0;
    uint32_t dim = 0;
    static constexpr uint32_t max_dim = 8;
    uint32_t domain[max_dim];
    uint32_t offset[max_dim];
    uint32_t box[max_dim];

    // TODO add TlSigCuboid and replace SigthreadInterval.
    SigthreadInterval with_timeline(uint32_t bits) const;
};

template <typename Stream>
Stream&& operator<<(Stream&& s, const ThreadCuboid& cuboid)
{
    CAMSPORK_REQUIRE_CMP(cuboid.dim, <=, cuboid.max_dim, "tried to print ThreadCuboid with too many dimensions");
    const uint32_t dim = cuboid.dim;
    auto print_list = [&s, dim] (const uint32_t* p)
    {
        s << "[";
        if (dim > 0) {
            s << p[0];
            for (uint32_t i = 1; i < dim; ++i) {
                s << ", " << p[i];
            }
        }
        s << "]";
    };
    s << "{\"task_index\": " << cuboid.task_index;
    s << ", \"domain\": ";
    print_list(cuboid.domain);
    s << ", \"offset\": ";
    print_list(cuboid.offset);
    s << ", \"box\": ";
    print_list(cuboid.box);
    s << "}";
    return static_cast<Stream&&>(s);
}

}
