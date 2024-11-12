#pragma once

#include <cassert>
#include <memory>
#include <stddef.h>
#include <stdint.h>
#include <vector>

#include "sigthread.hpp"

#include "../../include/exospork/syncv.h"

namespace exospork
{

struct SyncEnv;

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
}

struct VisRecordDebugData
{
    uint8_t original_actor_signature;
    std::vector<SigthreadInterval> visibility_set;
    std::vector<pending_await_t> pending_await_list;
};



// *** Primary Implemented Interface ***
SyncEnv* new_sync_env(const exospork_syncv_init_t& init);
void delete_sync_env(SyncEnv* p_env);
void on_r(SyncEnv* p_env, size_t N, exospork_syncv_value_t* values, SigthreadInterval accessor_set);
void on_rw(SyncEnv* p_env, size_t N, exospork_syncv_value_t* values, SigthreadInterval accessor_set);
void clear_values(SyncEnv* p_env, size_t N, exospork_syncv_value_t* values);
void on_fence(SyncEnv* p_env, SigthreadInterval V1, SigthreadInterval V2, bool transitive);



// *** Debug Inspection Interface ***
void debug_register_values(SyncEnv* p_env, size_t N, exospork_syncv_value_t* values);
void debug_unregister_values(SyncEnv* p_env, size_t N, exospork_syncv_value_t* values);
uint32_t debug_get_write_vis_record_id(const exospork_syncv_value_t* p_assignment_record);


struct SyncEnvDeleter
{
    void operator() (SyncEnv* p) const
    {
        delete_sync_env(p);
    }
};

using SyncEnv_unique_ptr = std::unique_ptr<SyncEnv, SyncEnvDeleter>;

}
