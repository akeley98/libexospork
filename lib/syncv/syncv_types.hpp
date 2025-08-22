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

    operator bool() const
    {
        return data != 0;
    }
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
    // Note: might replace these with std::vector or something.
    // Hence it's best to access this data with the accessor functions,
    // and don't assume uint32_t* for the future iterator type.
    uint32_t dim_data = 0;
    static constexpr uint32_t max_dim = 8;
    uint32_t domain_data[max_dim];
    uint32_t offset_data[max_dim];
    uint32_t box_data[max_dim];

    // State access.
    uint32_t dim() const { return dim_data; }
    uint32_t* domain() { return domain_data; }
    const uint32_t* domain() const { return domain_data; }
    uint32_t* offset() { return offset_data; }
    const uint32_t* offset() const { return offset_data; }
    uint32_t* box() { return box_data; }
    const uint32_t* box() const { return box_data; }

    // domain[dim_idx: dim_idx + 1] = [domain_0, domain_1]
    // offset[dim_idx: dim_idx + 1] = [offset_0, offset_1]
    // box[dim_idx: dim_idx + 1] = [box_0, box_1]
    void split_replace(
        uint32_t dim_idx,
        uint32_t domain_0, uint32_t domain_1,
        uint32_t offset_0, uint32_t offset_1,
        uint32_t box_0, uint32_t box_1)
    {
        const uint32_t new_dim = dim() + 1;
        CAMSPORK_REQUIRE_CMP(new_dim, <=, max_dim, "implementation limit: ThreadCuboid::max_dim exceeded");
        dim_data = new_dim;
        for (uint32_t dst = new_dim - 1; dst > dim_idx; --dst) {
            domain_data[dst] = domain_data[dst - 1];
            offset_data[dst] = offset_data[dst - 1];
            box_data[dst] = box_data[dst - 1];
        }
        domain_data[dim_idx] = domain_0;
        domain_data[dim_idx + 1] = domain_1;
        offset_data[dim_idx] = offset_0;
        offset_data[dim_idx + 1] = offset_1;
        box_data[dim_idx] = box_0;
        box_data[dim_idx + 1] = box_1;
    }

    uint32_t domain_num_threads() const
    {
        uint32_t prod = 1;
        for (uint32_t i = 0; i < dim(); ++i) {
            prod *= domain()[i];
        }
        return prod;
    };

    // TODO add TlSigInput and replace SigthreadInterval.
    SigthreadInterval with_timeline(uint32_t bitfield) const;

    // Initialize (end-begin)-dimensional domain with all threads active
    // i.e. offset = 0, box = domain.
    template <typename Iterator>
    static ThreadCuboid full(Iterator begin, Iterator end)
    {
        ThreadCuboid cuboid;
        cuboid.task_index = 0;
        const ptrdiff_t dim = end - begin;
        CAMSPORK_REQUIRE_CMP(dim, >=, 0, "iterators in wrong order?");
        CAMSPORK_REQUIRE_CMP(dim, <=, max_dim, "implementation limit: ThreadCuboid::max_dim exceeded");
        cuboid.dim_data = uint32_t(dim);
        for (Iterator it = begin; it != end; ++it) {
            const auto i = it - begin;
            const uint32_t c = uint32_t(*it);
            cuboid.domain_data[i] = c;
            cuboid.offset_data[i] = 0;
            cuboid.box_data[i] = c;
        }
        return cuboid;
    }
};

template <typename Stream>
Stream&& operator<<(Stream&& s, const ThreadCuboid& cuboid)
{
    const uint32_t dim = cuboid.dim();
    CAMSPORK_REQUIRE_CMP(dim, <=, cuboid.max_dim, "tried to print ThreadCuboid with too many dimensions");
    auto print_list = [dim, &s] (auto p)
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
    print_list(cuboid.domain());
    s << ", \"offset\": ";
    print_list(cuboid.offset());
    s << ", \"box\": ";
    print_list(cuboid.box());
    s << "}";
    return static_cast<Stream&&>(s);
}

}
