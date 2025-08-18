#pragma once

#include <memory>
#include <stdint.h>

namespace camspork
{

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

}
