#include "sync_env.hpp"

#include <algorithm>
#include <cassert>
#include <functional>
#include <memory>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <tuple>

#include "sigthread.hpp"
#include "../util/node_pool.hpp"

#include "../../include/exospork/syncv.h"

// Maybe replace later
#include <unordered_map>
template <typename K, typename V> using Map = std::unordered_map<K, V>;

namespace exospork
{

namespace
{

struct PendingAwaitListNode
{
    pending_await_t await_id;
    nodepool::id<PendingAwaitListNode> exospork_next_id;
};

// We encode a visibility set as a list of sorted, minimal
// sigthread intervals. The intervals are sorted in that
// a.tid_hi <= b.tid_lo for a before b in the list, and the list
// is minimal in that no more intervals are used than needed
// (mostly by merging adjacent intervals with the same bitfield).
//
// We encode both V_A (async visibility set) and V_S (sync visibility
// set) in compressed form, taking advantage of how V_S \subseteq V_A.
// V_A is the union of all intervals.
// V_S is the union of all intervals with sync_bit set.
struct SigthreadIntervalListNode
{
    SigthreadInterval data;
    nodepool::id<SigthreadIntervalListNode> exospork_next_id;
};

static_assert(sizeof(SigthreadIntervalListNode) == 16, "Check that you meant to change this perf-critical struct");

struct VisRecord
{
    // Owning reference to singly-linked list.
    nodepool::id<SigthreadIntervalListNode> visibility_set;

    // Owning reference to singly-linked list.
    nodepool::id<PendingAwaitListNode> pending_await_list;

    uint8_t original_actor_signature;

    // This has nothing to do with the main purpose of the struct; only needed for assignment_record_remove_duplicates.
    // This should be in VisRecordListNode conceptually, but that would waste 4 bytes.
    uint8_t tmp_is_duplicate;
};

struct VisRecordListNode
{
    // Count of owning references.
    // AssignmentRecord references (and ReadVisRecordListNode) are owning.
    // Forwarding references are owning.
    // Memoization table references are non-owning.
    uint32_t refcnt;

    // If in base state, this is the next node in the memoization bucket
    // If in the forwarding state, this is an owning reference to the
    // forwarded-to visibility record.
    nodepool::id<VisRecordListNode> exospork_next_id;

    // If the visibility record is in the base state, this is the valid data.
    // If the visibility record is in the forwarding state, the data is that of the record at get(exospork_next_id).
    VisRecord base_data;

    // Empty visibility set is never valid.
    // We will use that to indicate forwarding.
    bool is_forwarded() const
    {
        return base_data.visibility_set._1_index == 0;
    }
};

static_assert(sizeof(VisRecordListNode) == 20, "Check that you meant to change this perf-critical struct");

struct ReadVisRecordListNode
{
    nodepool::id<VisRecordListNode> vis_record_id;
    nodepool::id<ReadVisRecordListNode> exospork_next_id;
};

// Assignment record: 1 write visibility record + collection of read visibility records.
// This is associated for each variable of the program undergoing synchronization validation.
// i.e. the user should have parallel arrays of AssignmentRecord and program variable values.
struct AssignmentRecord
{
    // Single write visibility record (always exists)
    nodepool::id<VisRecordListNode> write_vis_record_id;

    // Zero or more read visibility records.
    nodepool::id<ReadVisRecordListNode> read_vis_records_head_id;

    // Unique ID for this assignment
    uint64_t assignment_id : 52;

    // See assignment_record_remove_duplicates.
    uint64_t last_sync_counter_bits : 12;
};

static_assert(sizeof(AssignmentRecord) == sizeof(exospork_syncv_value_t));

inline AssignmentRecord* from_api(exospork_syncv_value_t* p)
{
    return reinterpret_cast<AssignmentRecord*>(p);
}

inline const AssignmentRecord* from_api(const exospork_syncv_value_t* p)
{
    return reinterpret_cast<const AssignmentRecord*>(p);
}

inline exospork_syncv_value_t* to_api(AssignmentRecord* p)
{
    return reinterpret_cast<exospork_syncv_value_t*>(p);
}

struct BarrierState
{
    uint32_t max_pipelining = 0;
    uint32_t arrive_count = 0;
    uint32_t await_count = 0;
    SigthreadInterval arrive_sigthreads = {};
    SigthreadInterval await_sigthreads = {};
};

template <uint32_t Level> constexpr uint64_t bucket_level_size = 0;
template<> constexpr uint64_t bucket_level_size<0>  = 1;
template<> constexpr uint64_t bucket_level_size<1>  = 32;
template<> constexpr uint64_t bucket_level_size<2>  = 128;
template<> constexpr uint64_t bucket_level_size<3>  = 256;
template<> constexpr uint64_t bucket_level_size<4>  = 512;
template<> constexpr uint64_t bucket_level_size<5>  = 1024;
template<> constexpr uint64_t bucket_level_size<6>  = 0x1000;
template<> constexpr uint64_t bucket_level_size<7>  = 0x1'0000;
template<> constexpr uint64_t bucket_level_size<8>  = 0x10'0000;
template<> constexpr uint64_t bucket_level_size<9>  = 0x100'0000;
template<> constexpr uint64_t bucket_level_size<10> = 0x1000'0000;
template<> constexpr uint64_t bucket_level_size<11> = 0x1'0000'0000;
constexpr uint32_t bucket_level_count         = 12;

template <uint32_t BucketLevel> struct IntervalBucket;

template <uint32_t BucketLevel>
struct IntervalBucketParentPointer
{
    static_assert(BucketLevel < bucket_level_count);
    uint32_t child_index_in_parent = 0;
    IntervalBucket<BucketLevel + 1>* p_parent = nullptr;
};

template <>
struct IntervalBucketParentPointer<bucket_level_count - 1>
{
};

// Let Sz = bucket_level_size<BucketLevel>; an interval bucket for a given
// BucketLevel encompasses the interval of thread IDs [I * Sz, (I+1) * Sz - 1]
// for some index I. The idea is to store visibility records in the smallest
// possible (i.e. most specific) bucket that is still a superset.
//
// Unless BucketLevel == 0 (single thread buckets),
// the bucket has child interval buckets of one lower level,
// with the original interval evenly subdivided into N-many child buckets
// owned by the parent bucket.
//
// All buckets except the top-level bucket store a back pointer to
// their parent; see IntervalBucketParentPointer (This is a non-owning
// back ptr; child can't outlive parent).
//
// Buckets should be removed from the tree when empty, see
// delete_interval_bucket_if_empty.
template <uint32_t BucketLevel>
struct IntervalBucket : IntervalBucketParentPointer<BucketLevel>
{
    static_assert(BucketLevel != 0, "BucketLevel = 0 for illustration only");
    static_assert(BucketLevel != 1, "Should be template specialization");
    static_assert(BucketLevel < bucket_level_count);

    static constexpr uint32_t bucket_level = BucketLevel;
    static constexpr uint32_t child_count = bucket_level_size<BucketLevel> / bucket_level_size<BucketLevel - 1>;
    std::unique_ptr<IntervalBucket<BucketLevel - 1>> child_interval_buckets[child_count];

    // Count of non-null pointers in child_interval_buckets.
    uint32_t nonempty_child_count = 0;

    // Bucket for this interval.
    nodepool::id<VisRecordListNode> bucket = {0};

    // For making code re-entrant (prevent de-allocation while being visited).
    uint32_t visitor_count = 0;
};

template <>
struct IntervalBucket<1> : IntervalBucketParentPointer<1>
{
    static constexpr uint32_t bucket_level = 1;
    static constexpr uint32_t child_count = bucket_level_size<1>;
    nodepool::id<VisRecordListNode> child_interval_buckets[child_count] = {};

    // Count of non-null single_thread_buckets.
    uint32_t nonempty_child_count = 0;

    // Bucket for this interval.
    nodepool::id<VisRecordListNode> bucket = {0};

    // For making code re-entrant (prevent de-allocation while being visited).
    uint32_t visitor_count = 0;
};


// TODO consider further sub-bucketing, e.g. by sigbits.
// If we do this, we have to be careful not to double-consider items when
// moving between buckets. e.g. bucket by lowest to highest sigbit count.


template <uint32_t BucketLevel>
bool interval_bucket_is_empty(const IntervalBucket<BucketLevel>& bucket) noexcept
{
    return bucket.nonempty_child_count == 0 && !bucket.bucket && !bucket.visitor_count;
}

// De-allocate the given bucket if it's empty and not the top-level bucket.
// We presume that the bucket is owned by its parent (unique_ptr tree).
//
// We do not make any modifications to the parent except for nulling out the pointer.
// In particular, we don't change nonempty_child_count, or handle deleting the parent
// if it too is now empty.
template <uint32_t BucketLevel>
void delete_interval_bucket_if_empty(IntervalBucket<BucketLevel>* p) noexcept
{
    if (interval_bucket_is_empty(*p)) {
        for (const auto& child : p->child_interval_buckets) {
            assert(!child);  // nonempty_child_count was wrong.
        }

        if constexpr (BucketLevel < bucket_level_count - 1) {
            // Parent pointer should be correct.
            IntervalBucket<BucketLevel + 1>* p_parent = p->p_parent;
            assert(p_parent);
            const uint32_t child_index = p->child_index_in_parent;
            assert(child_index < p_parent->child_count);
            assert(p_parent->nonempty_child_count > 0);

            // p is invalidated after this (unique_ptr reset).
            assert(p_parent->child_interval_buckets[child_index].get() == p);
            p_parent->child_interval_buckets[child_index].reset();
        }
    }
}

}  // end namespace



// "Everything" struct that implements the synchronization environment.
struct SyncEnv
{
    // Failure flag
    // This is to be set upon an exception being thrown through the non-private
    // interface. If this happens, the internal state may be inconsistent, but memory
    // shouldn't be formally leaked since we'll delete the memory pools later anyway.
    // Ignore all further SyncEnv commands if failed = true.
    bool failed = false;

    // Counters for operations
    uint64_t assignment_counter = 0;  // Write counter
    uint64_t sync_counter = 0;        // Sync operation counter
    uint64_t operation_counter = 0;   // All operations counter

    // Memory pool state.
    uintptr_t original_memory_budget = 0;
    uintptr_t current_memory_budget = 0;
    std::tuple<
        nodepool::Pool<SigthreadIntervalListNode>,
        nodepool::Pool<PendingAwaitListNode>,
        nodepool::Pool<VisRecordListNode>,
        nodepool::Pool<ReadVisRecordListNode>> pool_tuple;

    // Barrier state. TODO
    // The Nth bit is set if N is allocated as a barrier ID.
    uint64_t live_barrier_bits[max_live_barriers / 64] = {0};
    BarrierState barrier_states[max_live_barriers];

    // Memoization table state.
    IntervalBucket<bucket_level_count - 1> top_level_bucket;

    // Debugging/Testing
    // Record pointers to registered AssignmentRecord arrays (plus the array sizes).
    // Ordinarily we rely on the user alone to remember these arrays, but we need to cache them if we are
    // doing consistency checks.
    Map<AssignmentRecord*, size_t> debug_registered_assignment_records;

    uint64_t debug_assignment_id = 0;
    uint64_t debug_operation_id = 0;



    // *** Memory Pool Allocators; Linked List Manipulation ***
    // See nodepool for more info.



    template <typename ListNode>
    ListNode& alloc_default_node(nodepool::id<ListNode>* out_id) noexcept
    {
        using TypedPool = nodepool::Pool<ListNode>;
        TypedPool& pool = std::get<TypedPool>(pool_tuple);
        return pool.alloc_default_node(&current_memory_budget, out_id);
    }

    template <typename ListNode>
    void extend_free_list(nodepool::id<ListNode> head_id) noexcept
    {
        using TypedPool = nodepool::Pool<ListNode>;
        TypedPool& pool = std::get<TypedPool>(pool_tuple);
        pool.extend_free_list(head_id);
    }

    template <typename ListNode>
    void insert_next_node(nodepool::id<ListNode>* p_insert_after, nodepool::id<ListNode> insert_me) noexcept
    {
        using TypedPool = nodepool::Pool<ListNode>;
        TypedPool& pool = std::get<TypedPool>(pool_tuple);
        return pool.insert_next_node(p_insert_after, insert_me);
    }

    // Given a pointer to the exospork_next_id member of a node in a list,
    // but don't add it to the free chain: the node is returned to the caller.
    template <typename ListNode>
    [[nodiscard]] nodepool::id<ListNode> remove_next_node(nodepool::id<ListNode>* p_id) noexcept
    {
        using TypedPool = nodepool::Pool<ListNode>;
        TypedPool& pool = std::get<TypedPool>(pool_tuple);
        return pool.remove_next_node(p_id);
    }

    // Given a pointer to the exospork_next_id member of a node in a list,
    // remove the next node of the list and add its memory to the free chain.
    // This shouldn't be used if the ListNode itself owns stuff.
    template <typename ListNode>
    void remove_and_free_next_node(nodepool::id<ListNode>* p_id) noexcept
    {
        nodepool::id<ListNode> victim_id = remove_next_node(p_id);
        assert(!get(victim_id).exospork_next_id);  // Should have been removed from list
        extend_free_list(victim_id);               // so this free only adds 1 node to free chain.
    }

    template <typename ListNode>
    ListNode& get(nodepool::id<ListNode> id) noexcept
    {
        using TypedPool = nodepool::Pool<ListNode>;
        return std::get<TypedPool>(pool_tuple).get(id);
    }

    template <typename ListNode>
    const ListNode& get(nodepool::id<ListNode> id) const noexcept
    {
        using TypedPool = nodepool::Pool<ListNode>;
        return std::get<TypedPool>(pool_tuple).get(id);
    }

    template <typename ListNode>
    uint32_t debug_node_pool_size() const noexcept
    {
        using TypedPool = nodepool::Pool<ListNode>;
        const uint32_t sz{std::get<TypedPool>(pool_tuple).size()};
        return sz;
    }

    // Increment reference count of visibility record.
    void incref(nodepool::id<VisRecordListNode> id) noexcept
    {
        VisRecordListNode& node = get(id);
        assert(node.refcnt != 0);
        node.refcnt++;
        assert(node.refcnt != 0);  // Overflow check
    }

    // Decrement reference count of visibility record,
    // and handle necessary free-ing in case of 0 refcnt.
    void decref(nodepool::id<VisRecordListNode> id) noexcept
    {
        assert(id);
        VisRecordListNode& node = get(id);
        if (0 == --node.refcnt) {
            if (node.is_forwarded()) {
                // Add physical storage of victim visibility record to free chain,
                // then decref owning reference to forwarded-to visibility record,
                auto fwd_id = node.exospork_next_id;
                assert(fwd_id);
                free_single_vis_record(id);
                decref(fwd_id);  // Hope for tail call.
            }
            else {
                // Non-forwarded (base) visibility record must be removed from memoization first.
                auto memoized_id = remove_memoized(&node);
                assert(id == memoized_id);
                assert(get(memoized_id).exospork_next_id == 0);  // Should have been removed from bucket's list.
                free_single_vis_record(memoized_id);
            }
        }
    }

    void reset_vis_record_data(VisRecord* p_data) noexcept
    {
        static_assert(sizeof(*p_data) == 12, "update me");
        p_data->original_actor_signature = ~0;
        extend_free_list(p_data->visibility_set);
        p_data->visibility_set = {};
        extend_free_list(p_data->pending_await_list);
        p_data->pending_await_list = {};
    }

    void free_single_vis_record(nodepool::id<VisRecordListNode> id) noexcept
    {
        assert(id);
        VisRecordListNode& node = get(id);
        assert(node.refcnt == 0);
        reset_vis_record_data(&node.base_data);
        node.exospork_next_id = {};  // Avoid freeing entire list.
        extend_free_list(id);
    }

    void reset_assignment_record(AssignmentRecord* p_record) noexcept
    {
        // Decref write visibility record.
        if (p_record->write_vis_record_id) {
            decref(p_record->write_vis_record_id);
            p_record->write_vis_record_id = {};
        }

        // Decref read visibility records.
        auto read_id = p_record->read_vis_records_head_id;
        while (read_id) {
            ReadVisRecordListNode& node = get(read_id);
            decref(node.vis_record_id);
            read_id = node.exospork_next_id;
        }

        // Free physical storage of read vis record list.
        extend_free_list(p_record->read_vis_records_head_id);
        p_record->read_vis_records_head_id = {};

        p_record->assignment_id = 0;
    }



    // *** Operations on Visibility Records ***
    // If the visibility record is modified, you need to be careful to update the memoization table.



    // Insert to pending await list
    void add_pending_await(VisRecord* p, pending_await_t await_id)
    {
        nodepool::id<PendingAwaitListNode> new_node_id;
        auto& node = alloc_default_node(&new_node_id);
        assert(!node.exospork_next_id);
        node.await_id = await_id;
        insert_next_node(&p->pending_await_list, new_node_id);
    }

    // Remove from pending await if found (return found flag).
    bool remove_pending_await(VisRecord* p, pending_await_t await_id) noexcept
    {
        nodepool::id<PendingAwaitListNode>* p_node_id = &p->pending_await_list;
        nodepool::id<PendingAwaitListNode> node_id;

        while ((node_id = *p_node_id)) {
            PendingAwaitListNode& node = get(node_id);
            if (node.await_id == await_id) {
                remove_and_free_next_node(p_node_id);
                return true;
            }
            p_node_id = &node.exospork_next_id;
        }
        return false;
    }

    // Allocate a new visibility record whose visibility set consists of one interval.
    // This will later need to be added to the memoization table.
    nodepool::id<VisRecordListNode> alloc_visibility_record(SigthreadInterval init_interval)
    {
        assert(init_interval.tid_hi > init_interval.tid_lo);

        nodepool::id<SigthreadIntervalListNode> sigthreads_id;
        SigthreadIntervalListNode& sigthreads_node = alloc_default_node(&sigthreads_id);
        sigthreads_node.data = init_interval;

        nodepool::id<VisRecordListNode> vis_record_id;
        VisRecordListNode& vis_record = alloc_default_node(&vis_record_id);
        vis_record.refcnt = 1;
        vis_record.base_data.original_actor_signature = init_interval.get_unique_actor_signature();
        vis_record.base_data.visibility_set = sigthreads_id;

        assert(equal(vis_record.base_data, init_interval));

        return vis_record_id;
    }

    // Union sigthread interval into the visibility set(s).
    // Caller will have to make changes to the memoization table afterwards.
    // Recall V_A (async visibility set) and V_S (sync visibility set) has V_S \subseteq V_A
    // and sync_bit on an interval means to include it in both V_S and V_A, and not just V_A.
    void union_sigthread_interval(VisRecord* p, SigthreadInterval input)
    {
        // The current model doesn't allow for augmenting only V_A without V_S, remove this if that changes.
        assert(!input.async_only());

        // Non-empty input check (cross product of non-empty thread interval and non-empty actor signature set).
        assert(input.tid_hi > input.tid_lo);
        assert(0 != (input.bitfield & ~input.sync_bit));
        static_assert(input.sync_bit == 0x8000'0000, "review this code for needed changes");
        using node_id = nodepool::id<SigthreadIntervalListNode>;

        // Visibility set must not be created empty (or VisRecord is in forwarding state). See alloc_visibility_record.
        assert(p->visibility_set);

        // Modify and/or add intervals.
        // We can view each pointer-to-node_id as a "gap" between intervals (imagine an arrow between nodes = a gap).
        // We abuse language to consider there to be a "gap" before the leftmost and after the rightmost intervals.
        // This runs for N+1 iterations where N is the current interval count.
        uint32_t gap_tid_lo = 0;
        for (node_id* p_id = &p->visibility_set; 1;) {
            // Must remember the next iteration's node now, so we don't get confused by insertion.
            const node_id original_next_node_id = *p_id;

            uint32_t gap_tid_hi;
            if (original_next_node_id) {
                gap_tid_hi = get(original_next_node_id).data.tid_lo;
                assert(gap_tid_lo <= gap_tid_hi);  // intervals were out of order.
            }
            else {
                gap_tid_hi = input.tid_hi;  // For "gap" after the rightmost interval.
            }

            // If the gap is non-empty and overlaps the input interval, we need to insert a new interval.
            SigthreadInterval new_interval{};
            new_interval.tid_lo = std::max(gap_tid_lo, input.tid_lo);
            new_interval.tid_hi = std::min(gap_tid_hi, input.tid_hi);
            new_interval.bitfield = input.bitfield;
            if (new_interval.tid_hi > new_interval.tid_lo) {
                node_id new_node_id{};
                SigthreadIntervalListNode& new_node = alloc_default_node(&new_node_id);
                new_node.data = new_interval;
                insert_next_node(p_id, new_node_id);
            }

            if (!original_next_node_id) {
                break;
            }

            // This runs N times (leftmost interval is the "next node" of the imaginary before-left gap)
            // We will modify each original interval, which is the "next node" of the gap just processed.
            //
            // The interval may be subdivided into up to 3 intervals depending on overlap with input interval,
            // since only the overlapped portion should have its bitfield modified.
            //
            // 1st interval: keeps original bits, left of intersection.
            // 2nd interval: bitfield augmented, footprint of intersection.
            // 3rd interval: keeps original bits, right of intersection.
            SigthreadIntervalListNode& next_node = get(original_next_node_id);
            const SigthreadInterval original_data = next_node.data;
            const uint32_t intersect_tid_lo = std::max(original_data.tid_lo, input.tid_lo);
            const uint32_t intersect_tid_hi = std::min(original_data.tid_hi, input.tid_hi);
            const uint32_t added_bits = input.bitfield & ~original_data.bitfield;
            const bool change_needed = added_bits != 0 && (intersect_tid_lo < intersect_tid_hi);

            // Possibly add 1st interval.
            if (change_needed && original_data.tid_lo < intersect_tid_lo) {
                node_id new_node_id{};
                SigthreadIntervalListNode& new_node = alloc_default_node(&new_node_id);
                new_node.data.tid_lo = original_data.tid_lo;
                new_node.data.tid_hi = intersect_tid_lo;
                new_node.data.bitfield = original_data.bitfield;
                assert(*p_id == original_next_node_id);
                insert_next_node(p_id, new_node_id);
            }

            // Now update iteration state (we do this now so the following insertions work)
            p_id = &next_node.exospork_next_id;
            gap_tid_lo = next_node.data.tid_hi;

            if (change_needed) {
                // 2nd interval; we recycle the existing node since the 2nd interval is guaranteed non-empty
                next_node.data.tid_lo = intersect_tid_lo;
                next_node.data.tid_hi = intersect_tid_hi;
                next_node.data.bitfield |= added_bits;

                // Possibly add 3rd interval, insert after 2nd interval.
                if (intersect_tid_hi < original_data.tid_hi) {
                    node_id new_node_id{};
                    SigthreadIntervalListNode& new_node = alloc_default_node(&new_node_id);
                    new_node.data.tid_lo = intersect_tid_hi;
                    new_node.data.tid_hi = original_data.tid_hi;
                    new_node.data.bitfield = original_data.bitfield;
                    insert_next_node(p_id, new_node_id);
                    p_id = &new_node.exospork_next_id;  // Need to point to the real gap (after original_data.tid_hi)
                }
            }
        }

        // Merge redundant intervals
        // For each "current node", we try to merge it with the next node, if it exists.
        assert(p->visibility_set);
        SigthreadIntervalListNode* p_current_node = &get(p->visibility_set);
        for (node_id next_id; (next_id = p_current_node->exospork_next_id); ) {
            SigthreadIntervalListNode* p_next_node = &get(next_id);

            SigthreadInterval& current = p_current_node->data;
            const SigthreadInterval& next = p_next_node->data;
            assert(current.tid_lo < current.tid_hi);
            assert(current.tid_hi <= next.tid_lo);
            assert(next.tid_lo < next.tid_hi);

            if (current.tid_hi == next.tid_lo && current.bitfield == next.bitfield) {
                // Merge next node into current node, and remove next node from list.
                current.tid_hi = next.tid_hi;
                assert(p_current_node->exospork_next_id == next_id);
                remove_and_free_next_node(&p_current_node->exospork_next_id);
            }
            else {
                // If not merged, we need to process the next node in the next iteration.
                // If we did merge, we keep the same current node, since we may need to merge with the new next node.
                p_current_node = p_next_node;
            }
        }
    }

    bool valid_adjacent(SigthreadInterval first, SigthreadInterval second) const
    {
        // Check that the two sigthread intervals can be adjacent in a valid visibility set encoding.
        return first.tid_lo < first.tid_hi && first.tid_hi <= second.tid_lo && second.tid_lo < second.tid_hi
          && (first.tid_hi < second.tid_lo || first.bitfield != second.bitfield);
    }

    // Check if visibility records are equal.
    bool equal(const VisRecord& a, const VisRecord& b)
    {
        static_assert(sizeof(a) == 12, "Update me");

        // Must not have empty visibility set (forwarding state passed?)
        assert(a.visibility_set);
        assert(b.visibility_set);

        // Check equal original actor signature.
        if (a.original_actor_signature != b.original_actor_signature) {
            return false;
        }

        // Check equal intervals. We rely on (and enforce) the non-redundant encoding requirement.
        {
            using node_id = nodepool::id<SigthreadIntervalListNode>;
            node_id id_a = a.visibility_set;
            node_id id_b = b.visibility_set;

            while (id_a && id_b) {
                const SigthreadIntervalListNode& current_a = get(id_a);
                const SigthreadIntervalListNode& current_b = get(id_b);
                id_a = current_a.exospork_next_id;
                id_b = current_b.exospork_next_id;
                assert(!id_a || valid_adjacent(current_a.data, get(id_a).data));
                assert(!id_b || valid_adjacent(current_b.data, get(id_b).data));

                if (current_a.data != current_b.data) {
                    return false;
                }
            }

            if (id_a != id_b) {
                return false;  // Lists have different lengths.
            }
        }

        // Check equal pending awaits.
        // TODO explain why it's OK not to check for re-ordering. Should be OK because there's a total ordering
        // of how pending awaits are simulated and added to lists.
        using node_id = nodepool::id<PendingAwaitListNode>;
        node_id id_a = a.pending_await_list;
        node_id id_b = b.pending_await_list;

        while (id_a && id_b) {
            const PendingAwaitListNode& current_a = get(id_a);
            const PendingAwaitListNode& current_b = get(id_b);
            id_a = current_a.exospork_next_id;
            id_b = current_b.exospork_next_id;

            if (current_a.await_id != current_b.await_id) {
                return false;
            }
        }

        return id_a == id_b;  // Check lists had the same length.
    }

    // Check if a visibility record matches what would have been constructed
    // from the given sigthread interval set by alloc_visibility_record.
    bool equal(const VisRecord& a, SigthreadInterval interval)
    {
        static_assert(sizeof(a) == 12, "Update me");

        // Must not have empty visibility set (forwarding state passed?)
        assert(a.visibility_set);

        if (a.pending_await_list) {
            return false;
        }

        // Check if only one interval, and it equals the input interval, with correct original actor signature.
        const SigthreadIntervalListNode& node = get(a.visibility_set);
        return !node.exospork_next_id
                 && node.data == interval
                 && (interval.sigbits() >> a.original_actor_signature == 1u);
    }

    template <bool SyncOnly, bool Transitive>
    bool visible_to_impl(const VisRecord& vis_record, SigthreadInterval access_set)
    {
        // Must not have empty visibility set (forwarding state passed?)
        assert(vis_record.visibility_set);

        const uint32_t sigbits_mask = Transitive ? ~uint32_t(0) : uint32_t(1) << vis_record.original_actor_signature;

        nodepool::id<SigthreadIntervalListNode> id = vis_record.visibility_set;
        while (id) {
            const SigthreadIntervalListNode& current_node = get(id);
            id = current_node.exospork_next_id;
            assert(!id || valid_adjacent(current_node.data, get(id).data));

            if (current_node.data.intersects(access_set, sigbits_mask)) {
                if (!SyncOnly || !current_node.data.async_only()) {
                    return true;
                }
            }
        }
        return false;
    }

    // Check if the visibility record is visible-to an access with the given sigthread access set.
    bool visible_to(const VisRecord& vis_record, SigthreadInterval accessor_set)
    {
        return visible_to_impl<true, true>(vis_record, accessor_set);
    }

    // Check if the visibility record synchronizes-with a synchronization statement with the given first visibility set.
    template <bool Transitive>
    bool synchronizes_with(const VisRecord& vis_record, SigthreadInterval V1)
    {
        return visible_to_impl<false, Transitive>(vis_record, V1);
    }



    // *** Memoization ***



    // Get the smallest possible sigthread interval that is a superset of
    // the given visibility set (ignore sync_bit); assumes non-empty input set.
    // This is needed to index into the correct bucket (the smallest one possible containing the visibility set).
    // Note, at time of writing the sigbits aren't used for bucketing, but maybe they should be.
    SigthreadInterval minimal_superset_interval(nodepool::id<SigthreadIntervalListNode> id) const
    {
        assert(id);
        const SigthreadIntervalListNode* p_node = &get(id);
        p_node->data.assert_valid();
        SigthreadInterval ret = p_node->data;

        while (1) {
            id = p_node->exospork_next_id;
            if (!id) {
                ret.bitfield &= ~ret.sync_bit;
                ret.assert_valid();
                return ret;
            }

            p_node = &get(id);
            assert(p_node->data.tid_lo >= ret.tid_hi);  // Not sorted?
            p_node->data.assert_valid();
            ret.tid_hi = p_node->data.tid_hi;
            ret.bitfield |= p_node->data.bitfield;
        }
    }

    // "Remove forwarding"; replace ID of forwarded visibility record
    // with ID of base visibility record that the original record forwarded to.
    // Assumes the given ID is intended as an owning ID.
    // Return record data.
    VisRecord& remove_forwarding(nodepool::id<VisRecordListNode>* p_id) noexcept
    {
        const nodepool::id<VisRecordListNode> old_id = *p_id;
        nodepool::id<VisRecordListNode> id = old_id;
        assert(id);
        VisRecordListNode* p_node = &get(id);
        assert(p_node->refcnt != 0);

        if (!p_node->is_forwarded()) {
            return p_node->base_data;  // No ID change
        }

        // Resolve the forwarding.
        do {
            id = p_node->exospork_next_id;
            assert(id);
            p_node = &get(id);
            assert(p_node->refcnt != 0);
        } while (p_node->is_forwarded());

        assert(id != old_id);
        incref(id);
        decref(old_id);  // Will take care of deallocating chain of forwarding if needed.
        assert(*p_id == old_id);
        *p_id = id;
        return p_node->base_data;
    }

    // Like remove_forwarding but non-destructive, i.e., don't actually replace the ID of a forwarding visibility
    // record with that of the forwarded-to base visibility record.
    VisRecord const_resolve_forwarding(nodepool::id<VisRecordListNode> id) const noexcept
    {
        assert(id);
        const VisRecordListNode* p_node = &get(id);
        assert(p_node->refcnt != 0);

        while (p_node->is_forwarded()) {
            id = p_node->exospork_next_id;
            assert(id);
            p_node = &get(id);
            assert(p_node->refcnt != 0);
        }

        return p_node->base_data;
    }


    enum class BucketProcessType
    {
        Find = 0,
        Insert = 1,
        MapAll = 2,
    };

    // Skeleton code for modifying the memoization table, while maintaining
    // internal consistency.
    // This is based on this function being available
    //
    //     this->process_bucket(nodepool::id<VisRecordListNode>*, Command)
    //
    // which may modify or delete the bucket (linked list) that has been passed.
    //
    // Only buckets that intersect the minimal_superset are processed.
    // Furthermore, if the operation is "exact", only the smallest bucket
    // containing the minimal_superset is processed.
    //
    // The operation details depend on BucketProcessType:
    //
    // Find: exact; callback skipped if bucket empty.
    //   Returns ID given by process_bucket.
    //
    // Insert: exact; create and process new empty child bucket if needed.
    //   Returns ID given by process_bucket.
    //
    // MapAll: not exact; returns 0 ID.
    //   We process smaller buckets before larger buckets, on the assumption
    //   that process_bucket may move items from smaller to larger buckets
    //   (so we need to avoid double-processing). This is a subtle thing to
    //   account for if we modify the bucketing scheme.
    template <BucketProcessType Type, typename Command>
    nodepool::id<VisRecordListNode> for_buckets(SigthreadInterval minimal_superset, const Command& command)
    {
        return this->for_buckets_impl<Type>(&top_level_bucket,
                                            minimal_superset.tid_lo,
                                            minimal_superset.tid_hi, command);
    }

    template <BucketProcessType Type, uint32_t BucketLevel, typename Command>
    nodepool::id<VisRecordListNode> for_buckets_impl(IntervalBucket<BucketLevel>* p_bucket,
                                                     int64_t relative_tid_lo, int64_t relative_tid_hi,
                                                     const Command& command)
    {
        if constexpr (Type != BucketProcessType::Insert && BucketLevel < bucket_level_count - 1) {
            // Left behind empty bucket that should have been de-allocated.
            assert(!interval_bucket_is_empty(*p_bucket));
        }

        assert(relative_tid_lo < relative_tid_hi);  // Input interval needs to be non-empty
        constexpr bool ExactType = Type != BucketProcessType::MapAll;
        nodepool::id<VisRecordListNode> result_id{};

        // Calculate range of child buckets that intersect the input interval.
        constexpr uint32_t child_size{bucket_level_size<BucketLevel - 1>};
        const uint32_t child_min_index = relative_tid_lo < 0 ? 0u : uint32_t(relative_tid_lo) / child_size;
        const uint32_t child_max_index = std::min(uint32_t(relative_tid_hi - 1) / child_size,
                                                  uint32_t(p_bucket->child_count - 1));

        auto visit_child = [this, p_bucket, relative_tid_lo, relative_tid_hi, &command] (uint32_t child_index)
        {
            nodepool::id<VisRecordListNode> lambda_result_id = {};
            assert(child_index < p_bucket->child_count);
            auto& child_ref = p_bucket->child_interval_buckets[child_index];

            if (!child_ref && Type != BucketProcessType::Insert) {
                // Skip empty bucket if not inserting.
                return lambda_result_id;
            }

            try {
                if (!child_ref && Type == BucketProcessType::Insert) {
                    // Speculate that the child bucket will be filled.
                    // We will undo this later if wrong.
                    p_bucket->nonempty_child_count++;

                    if constexpr (BucketLevel != 1) {
                        // Create child interval bucket (unique_ptr).
                        child_ref.reset(new IntervalBucket<BucketLevel - 1>);
                        child_ref->p_parent = p_bucket;
                        child_ref->child_index_in_parent = child_index;
                    }
                }

                // Process the child bucket.
                // This may result in child_ref being nulled out.
                if constexpr (BucketLevel == 1) {
                    if constexpr (ExactType) {
                        lambda_result_id = this->process_bucket(&child_ref, command);
                    }
                    else {
                        this->process_bucket(&child_ref, command);
                    }
                }
                else {
                    const uint64_t offset = child_index * child_size;
                    lambda_result_id = this->for_buckets_impl<Type>(child_ref.get(),
                                                                    relative_tid_lo - offset, relative_tid_hi - offset,
                                                                    command);
                }
            }
            catch (...) {
                if (!child_ref) {
                    assert(p_bucket->nonempty_child_count > 0);
                    p_bucket->nonempty_child_count--;
                }
                throw;
            }

            // Child bucket may have been deallocated for being empty.
            if (!child_ref) {
                assert(p_bucket->nonempty_child_count);
                p_bucket->nonempty_child_count--;
            }
            return lambda_result_id;
        };

        try {
            p_bucket->visitor_count++;

            if constexpr (ExactType) {
                if (child_min_index == child_max_index) {
                    // tid interval fits in child bucket; visit it.
                    result_id = visit_child(child_min_index);
                }
                else if (Type == BucketProcessType::Insert || p_bucket->bucket) {
                    // tid interval doesn't fit in child, so the current bucket is the correct (exact) one.
                    result_id = this->process_bucket(&p_bucket->bucket, command);
                }
            }
            else {
                // Non-exact; we process smaller (child) buckets first before larger (this level's) buckets.
                for (uint32_t child_index = child_min_index; child_index <= child_max_index; ++child_index) {
                    visit_child(child_index);
                }
                if (p_bucket->bucket) {
                    this->process_bucket(&p_bucket->bucket, command);
                }
            }
        }
        catch (...) {
            // For the most part, I don't care for exception safety in this code
            // but I make a defensive exception here.
            p_bucket->visitor_count--;
            delete_interval_bucket_if_empty(p_bucket);
            throw;
        }

        p_bucket->visitor_count--;
        delete_interval_bucket_if_empty(p_bucket);
        return result_id;
    }

    // Find visibility record in memoization bucket for which lambda(const VisRecord&) returns true.
    // Returns pointer to ID of record found (non-owning), or null if not found.
    template <typename Lambda>
    nodepool::id<VisRecordListNode>* bucket_search(nodepool::id<VisRecordListNode>* p_bucket_head, Lambda&& lambda)
    {
        using node_id = nodepool::id<VisRecordListNode>;
        node_id* p_id = p_bucket_head;

        for (node_id id; (id = *p_id); ) {
            VisRecordListNode& node = get(id);
            assert(!node.is_forwarded());  // Should not be in memoization table.

            if (lambda(node.base_data)) {
                return p_id;
            }

            p_id = &node.exospork_next_id;
        }
        return nullptr;
    }

    struct NewVisRecordCommand
    {
        SigthreadInterval init_interval;
        uint32_t added_refcnt;
    };

    // Add a new visibility record, or return existing memoized one, constructed from the given sigthread interval set.
    // The returned ID is an owning reference (ownership count given by added_refcnt).
    [[nodiscard]] nodepool::id<VisRecordListNode> memoize_new_vis_record(SigthreadInterval init_interval,
                                                                         uint32_t added_refcnt)
    {
        assert(init_interval.tid_hi > init_interval.tid_lo);
        assert(init_interval.sigbits() != 0);
        assert(added_refcnt != 0);

        NewVisRecordCommand command{init_interval, added_refcnt};
        nodepool::id<VisRecordListNode> id = for_buckets<BucketProcessType::Insert>(init_interval, command);
        assert(id);
        assert(equal(const_resolve_forwarding(id), init_interval));
        return id;
    }

    nodepool::id<VisRecordListNode> process_bucket(nodepool::id<VisRecordListNode>* p_bucket_head,
                                                   NewVisRecordCommand command)
    {
        auto lambda = [this, command] (const VisRecord& record) {
            return equal(record, command.init_interval);
        };
        nodepool::id<VisRecordListNode>* p_found_id = bucket_search(p_bucket_head, lambda);

        nodepool::id<VisRecordListNode> new_id;

        if (p_found_id) {
            // Existing memoized entry found.
            new_id = *p_found_id;
            assert(new_id);
            get(new_id).refcnt += command.added_refcnt;
        }
        else {
            // Add memoized base visibility set entry to bucket of memoization table.
            new_id = alloc_visibility_record(command.init_interval);
            assert(!get(new_id).exospork_next_id);
            get(new_id).refcnt = command.added_refcnt;
            assert(!get(new_id).is_forwarded());
            insert_next_node(p_bucket_head, new_id);
        }
        assert(new_id);
        return new_id;
    }

    struct RemoveMemoizedCommand
    {
        const VisRecordListNode* p_node;
    };

    // This removes the given node from the memoization table, but does not decrement the reference count or free it.
    // Recall that the memoization table does not own (reference count) the VisRecords contained.
    [[nodiscard]] nodepool::id<VisRecordListNode> remove_memoized(const VisRecordListNode* p_node) noexcept
    {
        assert(p_node);
        assert(p_node->refcnt == 0);
        assert(!p_node->is_forwarded());

        RemoveMemoizedCommand command{p_node};
        auto bucket_key = minimal_superset_interval(p_node->base_data.visibility_set);
        return for_buckets<BucketProcessType::Find>(bucket_key, command);
    }

    // Find and remove node in bucket.
    nodepool::id<VisRecordListNode> process_bucket(nodepool::id<VisRecordListNode>* p_bucket_head,
                                                   RemoveMemoizedCommand command)
    {
        auto lambda = [this, command] (const VisRecord& record) {
            return equal(record, command.p_node->base_data);
        };
        nodepool::id<VisRecordListNode>* p_id = bucket_search(p_bucket_head, lambda);
        if (p_id) {
            return remove_next_node(p_id);
        }
        else {
            return {};
        }
    }

    struct MemoizeOrForwardCommand
    {
        nodepool::id<VisRecordListNode> input_id;
    };

    // Given an existing visibility record in the base state that's not in the memoization table, either
    //   * Add it to the memoization table, if it's unique.
    //   * Put it in the forwarding state (discard existing state) and forward to equal already-memoized record.
    void memoize_or_forward(nodepool::id<VisRecordListNode> id)
    {
        assert(id);
        VisRecordListNode& node = get(id);
        assert(node.refcnt != 0);
        assert(!node.is_forwarded());
        assert(!node.exospork_next_id);  // shouldn't be in any linked list (memoization bucket or forwarded?)

        MemoizeOrForwardCommand command{id};
        auto bucket_key = minimal_superset_interval(node.base_data.visibility_set);
        for_buckets<BucketProcessType::Insert>(bucket_key, command);
    }

    nodepool::id<VisRecordListNode> process_bucket(nodepool::id<VisRecordListNode>* p_bucket_head,
                                                   MemoizeOrForwardCommand command)
    {
        VisRecordListNode& input_node = get(command.input_id);
        VisRecord input_vis_record = input_node.base_data;

        auto lambda = [this, input_vis_record] (const VisRecord& record) {
            return equal(record, input_vis_record);
        };

        nodepool::id<VisRecordListNode>* p_id = bucket_search(p_bucket_head, lambda);
        if (p_id) {
            // If equivalent memoized node found in bucket, forward input node to it.
            const nodepool::id fwd_id = *p_id;
            assert(fwd_id);
            assert(fwd_id != command.input_id);  // Trying to memoize something already in the memoization table.

            reset_vis_record_data(&input_node.base_data);  // Clear data to put visibility record into forwarding state.
            input_node.exospork_next_id = fwd_id;
            assert(input_node.is_forwarded());
            incref(fwd_id);  // Forwarding reference is owning.
        }
        else {
            // Insert input node to memoization bucket. No refcnt changes needed for memoization.
            // IMPORTANT: this memoization is at the start of the bucket. This means if the caller of this function
            // is processing this bucket, the caller probably won't encounter this node. See process_buckets_for_sync.
            insert_next_node(p_bucket_head, command.input_id);
        }

        return command.input_id;
    }

    template <bool Transitive>
    struct FenceUpdateCommand
    {
        SigthreadInterval V1;
        SigthreadInterval V2;

        void update_for_sync(SyncEnv& env, VisRecord* p_record) const
        {
            if (env.synchronizes_with<Transitive>(*p_record, V1)) {
                env.union_sigthread_interval(p_record, V2);
            }
        };
    };

    template <bool Transitive>
    struct ArriveUpdateCommand
    {
        SigthreadInterval V1;
        pending_await_t await_id;

        void update_for_sync(SyncEnv& env, VisRecord* p_record) const
        {
            if (env.synchronizes_with<Transitive>(*p_record, V1)) {
                env.add_pending_await(p_record, await_id);
            }
        }
    };

    struct AwaitUpdateCommand
    {
        SigthreadInterval V1;
        SigthreadInterval V2;
        pending_await_t await_id;

        void update_for_sync(SyncEnv& env, VisRecord* p_record) const
        {
            if (env.remove_pending_await(p_record, await_id)) {
                assert(env.synchronizes_with<true>(*p_record, V1));
                env.union_sigthread_interval(p_record, V2);
            }
        }
    };

    // Big payoff for all this code: function that performs the effects of a synchronization statement with the given
    // sync type and given first/second visibility sets. This affects all visibility records whose visibility set
    // intersects with the first visibility set of the synchronization statement.
    //
    // The real entrypoints are the ones specialized for fence, arrive, await.
    template <typename Command>
    void update_vis_records_for_sync_impl(const Command& command)
    {
        ++sync_counter;

        // Only visibility sets that intersect the first visibility set (V1) can be updated by this sync.
        // This is even the case for Await, assuming the V1 for the corresponding Arrive was correctly given.
        const SigthreadInterval minimal_superset = command.V1;
        for_buckets<BucketProcessType::MapAll>(minimal_superset, command);
    }

    template <typename Command>
    void process_bucket_for_sync_impl(nodepool::id<VisRecordListNode>* p_bucket_head, const Command& command)
    {
        // The bucket update process for handling the effects of synchronization on visibility records is quite
        // risky actually. When we modify a visibility record, we temporarily remove it from the memoization bucket,
        // modify it, then attempt to re-insert it. This is fundamentally needed since the modification may
        // cause a duplicate to be created, or a node to be in the wrong bucket.
        //
        // However, this re-insertion may cause the memoization table to be modified unexpectedly.
        // We have to be very careful when traversing the bucket's linked list, and this also explains
        // the IntervalBucket::visitor_count value, if it still exists; (see for_buckets_impl).
        //
        // Note, everything will be left in an inconsistent state in case an exception is thrown.

        using node_id = nodepool::id<VisRecordListNode>;
        node_id* p_id = p_bucket_head;

        // This might be really confusing. p_id is a pointer to a node ID.
        // It could be a pointer to the bucket (itself the ID of the head of the bucket node list) or it
        // could be the pointer to the exospork_next_id member of the PREVIOUS node (relative to current_node).
        while (node_id current_node_id = *p_id) {
            // Now temporarily remove the current node from the bucket linked list.
            // ("next_node" reflects the "pointer to previous node" viewpoint explained above).
            // *p_id will now be the ID of the node that formerly was after current_node, which (if not ID = 0)
            // is the node that we should process on the next iteration.
            VisRecordListNode& current_node = get(remove_next_node(p_id));
            assert(!current_node.exospork_next_id);// Should have been removed from list.
            assert(!current_node.is_forwarded());  // Invalid empty visibility set (forwarding state memoized?)

            // Update the visibility record stored in the node.
            command.update_for_sync(*this, &current_node.base_data);
            assert(p_id != &current_node.exospork_next_id);

            // This is where the node might get re-inserted to the memoization table.
            // *p_id might change value here again, but it's guaranteed p_id doesn't point inside &current_node.
            memoize_or_forward(current_node_id);

            // This part is dicey. We removed the node from the memoization table, then possibly re-inserted it,
            // either into another bucket, or at the head of this bucket. See the weird assert below.
            // Also see process_bucket(, MemoizeOrForwardCommand).
            // It's possible we re-inserted exactly into its old place, so we need to do some special logic
            // to avoid getting stuck in an infinite loop.
            if (*p_id == current_node_id) {
                // I'm fairly sure this is the only reason this branch should happen, due to how we insert nodes
                // only at the head of buckets. If this assert goes off, the code may still be correct;
                // this is just a warning-to-self to check that my mental model is correct.
                assert(p_id == p_bucket_head);
                p_id = &get(current_node_id).exospork_next_id;
            }
        }
    }

    // Augment all visibility records that synchronize with the first visibility set of the fence.
    void update_vis_records_for_fence(SigthreadInterval V1, SigthreadInterval V2, bool transitive)
    {
        V2.bitfield |= SigthreadInterval::sync_bit;  // Augment both V_A and V_S.
        if (transitive) {
            FenceUpdateCommand<true> command{V1, V2};
            update_vis_records_for_sync_impl(command);
        }
        else {
            FenceUpdateCommand<false> command{V1, V2};
            update_vis_records_for_sync_impl(command);
        }
    }

    template <bool Transitive>
    void process_bucket(nodepool::id<VisRecordListNode>* p_bucket_head, const FenceUpdateCommand<Transitive>& command)
    {
        process_bucket_for_sync_impl(p_bucket_head, command);
    }

    // TODO arrive, await



    // *** Access Safety Checking (read/write safety) ***



    template <bool IsWrite>
    void on_access_impl(size_t N, AssignmentRecord* p_assignment_records, SigthreadInterval accessor_set)
    {
        // We will memoize the new visibility record once
        if (N == 0) {
            return;
        }
        nodepool::id<VisRecordListNode> vis_record_id = memoize_new_vis_record(accessor_set, uint32_t(N));

        for (size_t i = 0; i < N; ++i) {
            AssignmentRecord& assignment_record = p_assignment_records[i];

            // Check against previous write visibility record
            if (assignment_record.write_vis_record_id) {
                const VisRecord& write_record = remove_forwarding(&assignment_record.write_vis_record_id);
                assert(visible_to(write_record, accessor_set));  // TODO shouldn't be assert: RAW or WAW
            }

            // If the access is a write, also check against the list of previous read visibility records.
            if constexpr (IsWrite) {
                nodepool::id<ReadVisRecordListNode> read_id = assignment_record.read_vis_records_head_id;
                while (read_id) {
                    ReadVisRecordListNode& node = get(read_id);
                    const VisRecord& read_record = remove_forwarding(&node.vis_record_id);
                    assert(visible_to(read_record, accessor_set));  // TODO shouldn't be assert: WAR hazard
                    read_id = get(read_id).exospork_next_id;
                }
            }

            // Add new visibility record (either as new write visibility record, or appended read visibility record).

            if constexpr (IsWrite) {
                // Clear out assignment record on write and add the single write visibility record.
                reset_assignment_record(&assignment_record);
                assignment_record.write_vis_record_id = vis_record_id;
                assignment_record.assignment_id = ++assignment_counter;
                assert(assignment_record.assignment_id != 0);
                assignment_record.last_sync_counter_bits = sync_counter;
            }
            else {
                // Add the new visibility record to the list of read visibility records.
                nodepool::id<ReadVisRecordListNode> read_id;
                ReadVisRecordListNode& read_node = alloc_default_node(&read_id);
                read_node.vis_record_id = vis_record_id;
                insert_next_node(&assignment_record.read_vis_records_head_id, read_id);

                // If we leave things as-is, read vis records may build up indefinitely for variables that are written
                // once and read many times. We fix this by removing duplicates; however, this is really expensive,
                // so we only do it once after each synchronization event (synchronization is when memoization kicks
                // in to potentially allow us to recognize duplicates due to duplicated IDs).
                const auto old_bits = assignment_record.last_sync_counter_bits;
                assignment_record.last_sync_counter_bits = sync_counter;
                if (old_bits != assignment_record.last_sync_counter_bits) {
                    // This could fail if the bits of sync_counter overflow exactly.
                    // However, this is unlikely, and is only a performance issue if so (we fail to remove duplicates).
                    assignment_record_remove_duplicates(&assignment_record);
                }
            }
        }
    }

    void on_r(size_t N, AssignmentRecord* p_assignment_records, SigthreadInterval accessor_set)
    {
        on_access_impl<false>(N, p_assignment_records, accessor_set);
    }

    void on_rw(size_t N, AssignmentRecord* p_assignment_records, SigthreadInterval accessor_set)
    {
        on_access_impl<true>(N, p_assignment_records, accessor_set);
    }

    void clear_values(size_t N, AssignmentRecord* p_assignment_records)
    {
        for (size_t i = 0; i < N; ++i) {
            reset_assignment_record(p_assignment_records + i);
        }
    }

    // Resolve forwarding and remove duplicate read visibility records.
    // Removing forwarding causes two equivalent read visibility records to have identical IDs
    // (both referring to the shared entry in the memoization table).
    void assignment_record_remove_duplicates(AssignmentRecord* p_assignment_record)
    {
        using node_id = nodepool::id<ReadVisRecordListNode>;

        // Clear out tmp_is_duplicate to 0.
        for (node_id id = p_assignment_record->read_vis_records_head_id; id; ) {
            const ReadVisRecordListNode node = get(id);
            get(node.vis_record_id).base_data.tmp_is_duplicate = 0;
            id = node.exospork_next_id;
        }

        // Remove duplicates, using remove_forwarding (unique ID iff unique record)
        // and tmp_is_duplicate to recognize duplicates.
        node_id* p_read_id = &p_assignment_record->read_vis_records_head_id;
        while (node_id next_id = *p_read_id) {
            ReadVisRecordListNode& next_node = get(next_id);
            uint8_t& is_duplicate = remove_forwarding(&next_node.vis_record_id).tmp_is_duplicate;

            if (is_duplicate) {
                // Duplicate, remove next node from list (decrements refcount for duplicated vis record).
                // This causes (next_id = *p_read_id) to change, so we don't have to update p_read_id.
                // i.e. since we removed the next node, we're ready to process a new next node next iteration.
                node_id victim_id = remove_next_node(p_read_id);
                assert(victim_id == next_id);
                assert(next_node.exospork_next_id == 0);
                decref(next_node.vis_record_id);
                extend_free_list(victim_id);
            }
            else {
                // If next node survives, remember the visibility set ID and move on.
                is_duplicate = 1;
                p_read_id = &get(next_id).exospork_next_id;
            }
        }
    }



    // *** Debugging / Testing ***



    // Get IDs of read visibility records of assignment record.
    void debug_get_read_vis_record_ids(const AssignmentRecord& record, std::vector<uint32_t>* out) const
    {
        out->clear();
        nodepool::id<ReadVisRecordListNode> id = record.read_vis_records_head_id;
        while (id) {
            const ReadVisRecordListNode& node = get(id);
            out->push_back(node.vis_record_id._1_index);
            id = node.exospork_next_id;
        }
    }

    // Get info for a given visibility record.
    void debug_get_vis_record_data(uint32_t id, VisRecordDebugData* out) const
    {
        assert(id);
        const VisRecord record = const_resolve_forwarding(nodepool::id<VisRecordListNode>{id});

        out->original_actor_signature = record.original_actor_signature;

        out->visibility_set.clear();
        for (nodepool::id<SigthreadIntervalListNode> node_id = record.visibility_set; node_id;) {
            const SigthreadIntervalListNode& node = get(node_id);
            out->visibility_set.push_back(node.data);
            node_id = node.exospork_next_id;
        }

        out->pending_await_list.clear();
        for (nodepool::id<PendingAwaitListNode> node_id = record.pending_await_list; node_id; ) {
            const PendingAwaitListNode& node = get(node_id);
            out->pending_await_list.push_back(node.await_id);
            node_id = node.exospork_next_id;
        }
    }

    void debug_register_values(size_t N, AssignmentRecord* values)
    {
        [[maybe_unused]] const bool unique = debug_registered_assignment_records.insert({values, N}).second;
        assert(unique);
    }

    void debug_unregister_values(size_t N, AssignmentRecord* values)
    {
        auto it = debug_registered_assignment_records.find(values);
        assert(it != debug_registered_assignment_records.end());
        assert(N == it->second);
        debug_registered_assignment_records.erase(it);
    }

    void validate_state() const
    {
        // TODO
    }
};



// *** Primary Implemented Interface ***



#define INTERFACE_PROLOGUE(p_env) \
try { \
    if (p_env->failed) { \
        return; \
    } \
    p_env->operation_counter++;

#define INTERFACE_EPILOGUE(p_env) \
} \
catch (...) { \
    p_env->failed = true; \
    throw; \
} \
if (p_env->operation_counter == p_env->debug_operation_id) { \
    p_env->validate_state(); \
}

SyncEnv* new_sync_env(const exospork_syncv_init_t& init)
{
    SyncEnv* p_env = new SyncEnv;
    p_env->original_memory_budget = init.memory_budget;
    p_env->current_memory_budget = init.memory_budget;
    p_env->debug_assignment_id = init.debug_assignment_id;
    p_env->debug_operation_id = init.debug_operation_id;
    return p_env;
}

void delete_sync_env(SyncEnv* p_env)
{
    // XXX temporary assert. This could go off due to the user not free-ing their own stuff.
    assert(p_env->failed || interval_bucket_is_empty(p_env->top_level_bucket));

    delete p_env;
}

void on_r(SyncEnv* p_env, size_t N, exospork_syncv_value_t* values, SigthreadInterval accessor_set)
{
    INTERFACE_PROLOGUE(p_env)
    p_env->on_r(N, from_api(values), accessor_set);
    INTERFACE_EPILOGUE(p_env)
}

void on_rw(SyncEnv* p_env, size_t N, exospork_syncv_value_t* values, SigthreadInterval accessor_set)
{
    INTERFACE_PROLOGUE(p_env)
    p_env->on_rw(N, from_api(values), accessor_set);
    INTERFACE_EPILOGUE(p_env)
}

void clear_values(SyncEnv* p_env, size_t N, exospork_syncv_value_t* values)
{
    INTERFACE_PROLOGUE(p_env)
    p_env->clear_values(N, from_api(values));
    INTERFACE_EPILOGUE(p_env)
}

void on_fence(SyncEnv* p_env, SigthreadInterval V1, SigthreadInterval V2, bool transitive)
{
    INTERFACE_PROLOGUE(p_env)
    p_env->update_vis_records_for_fence(V1, V2, transitive);
    INTERFACE_EPILOGUE(p_env)
}

// TODO arrive/await/barrier



// *** Debug Inspection Interface ***



void debug_register_values(SyncEnv* p_env, size_t N, exospork_syncv_value_t* values)
{
    p_env->debug_register_values(N, from_api(values));
}

void debug_unregister_values(SyncEnv* p_env, size_t N, exospork_syncv_value_t* values)
{
    p_env->debug_unregister_values(N, from_api(values));
}

uint32_t debug_get_write_vis_record_id(const exospork_syncv_value_t* p_assignment_record)
{
    return from_api(p_assignment_record)->write_vis_record_id._1_index;
}

void debug_get_read_vis_record_ids(const SyncEnv* p_env, const exospork_syncv_value_t* p_assignment_record,
                                   std::vector<uint32_t>* out)
{
    p_env->debug_get_read_vis_record_ids(*from_api(p_assignment_record), out);
}

void debug_get_vis_record_data(const SyncEnv* p_env, uint32_t id, VisRecordDebugData* out)
{
    p_env->debug_get_vis_record_data(id, out);
}

}  // end namespace
