#pragma once

#include <cassert>
#include <memory>
#include <stddef.h>
#include <stdint.h>
#include <vector>

#include "syncv_types.hpp"
#include "sigthread.hpp"

namespace camspork
{

struct SyncvTable;

// We record pending barrier awaits as (barrier ID, counter) pairs.
// barrier_id_bits many bits are used for the barrier ID.
// (32 - barrier_id_bits) bits are used for the counter.
//
// This limits the number of live barriers and the number of times
// a barrier can be used. The latter is capped in a real CUDA program
// by the number of times an mbarrier can be used: pow(2, 20).
constexpr uint32_t barrier_id_bits = 12;
constexpr uint32_t max_live_barriers = 1u << barrier_id_bits;
using pending_await_t = uint32_t;

inline uint32_t pending_await_barrier_id(pending_await_t id)
{
    return id & ((1u << barrier_id_bits) - 1u);
}

inline uint32_t pending_await_counter(pending_await_t id)
{
    return id >> barrier_id_bits;
}

inline pending_await_t pack_pending_await(uint32_t barrier_id, uint32_t counter)
{
    const uint32_t id = barrier_id | counter << barrier_id_bits;
    assert(pending_await_barrier_id(id) == barrier_id);
    assert(pending_await_counter(id) == counter);
    return id;
}

struct VisRecordDebugData
{
    uint8_t original_actor_signature;
    std::vector<SigthreadInterval> visibility_set;
    std::vector<pending_await_t> pending_await_list;
};



// *** Primary Implemented Interface ***
SyncvTable* new_syncv_table(const syncv_init_t& init);
SyncvTable* copy_syncv_table(const SyncvTable* table);
void delete_syncv_table(SyncvTable* table);
void on_r(SyncvTable* table, size_t N, assignment_record_id* array, SigthreadInterval accessor_set);
void on_rw(SyncvTable* table, size_t N, assignment_record_id* array, SigthreadInterval accessor_set);
void clear_visibility(SyncvTable* table, size_t N, assignment_record_id* array);
void alloc_barriers(SyncvTable* table, size_t N, barrier_id* bar);
void free_barriers(SyncvTable* table, size_t N, barrier_id* bar);
void on_fence(SyncvTable* table, SigthreadInterval V1, SigthreadInterval V2_full, SigthreadInterval V2_temporal, bool transitive);
void on_arrive(SyncvTable* table, barrier_id* bar, SigthreadInterval V1, bool transitive);
void on_await(SyncvTable* table, barrier_id* bar, SigthreadInterval V2_full, SigthreadInterval V2_temporal);
void begin_no_checking(SyncvTable* table);
void end_no_checking(SyncvTable* table);



// *** Debug Inspection Interface ***
void debug_register_records(SyncvTable* table, size_t N, assignment_record_id* records);
void debug_unregister_records(SyncvTable* table, size_t N, assignment_record_id* records);
void debug_get_read_vis_record_data(const SyncvTable* table, uint32_t id, VisRecordDebugData* out);
void debug_get_mutate_vis_record_data(const SyncvTable* table, uint32_t id, VisRecordDebugData* out);
void debug_validate_state(SyncvTable* table);
void debug_pre_delete_check(SyncvTable* table);

}
