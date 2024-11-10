#pragma once

#include <cassert>
#include <stdint.h>

namespace exospork
{

// A single sigthread consists (conceptually) of a pair of
// (thread ID, actor signature). We never store this directly.
// Instead, we work with sets of sigthreads.
//
// A sigthread interval is the cross product
//     [tid_lo, tid_hi) \times A
// where A is a set of actor signatures defined by the bits set in
//     bitfield & sigbits_mask
//
// The async bit is used elsewhere, to help store V_A and V_S compactly.
// (The async visibility set and sync visibility set).
struct SigthreadInterval
{
    // Thread index range [tid_lo, tid_hi)
    uint32_t tid_lo, tid_hi;

    // sigbits | (async_only ? 0u : sync_bit)
    uint32_t bitfield;

    // Top bit reserved, other bits used for sigbits.
    static constexpr uint32_t sync_bit = 0x8000'0000;

    bool async_only() const
    {
        return (bitfield & sync_bit) == 0;
    }

    void assert_valid() const
    {
        assert(sigbits() != 0);
        assert(tid_lo <= tid_hi);
    }

    uint32_t sigbits() const
    {
        return bitfield & ~sync_bit;
    }

    bool operator==(SigthreadInterval other) const
    {
        uint32_t diff = tid_lo ^ other.tid_lo;
        diff |= tid_hi ^ other.tid_hi;
        diff |= bitfield ^ other.bitfield;
        return diff == 0;
    }

    bool operator!=(SigthreadInterval other) const
    {
        return !(*this == other);
    }

    bool intersects(const SigthreadInterval& other) const
    {
        // <= due to tid_hi being an exclusive bound.
        const bool tid_disjoint = tid_hi <= other.tid_lo || other.tid_hi <= tid_lo;
        const uint32_t this_sigbits = sigbits();
        const uint32_t other_sigbits = other.sigbits();
        assert(this_sigbits != 0);
        assert(other_sigbits != 0);
        return !tid_disjoint && 0 != (this_sigbits & other_sigbits);
    }
};


}  // end namespace
