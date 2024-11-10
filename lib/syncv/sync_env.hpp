#pragma once

#include <cassert>
#include <stddef.h>
#include <stdint.h>

namespace exospork
{

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

inline uint32_t pack_pending_await(uint32_t barrier_id, uint32_t counter)
{
    const uint32_t id = barrier_id | counter << barrier_id_bits;
    assert(pending_await_barrier_id(id) == barrier_id);
    assert(pending_await_counter(id) == counter);
}

}
