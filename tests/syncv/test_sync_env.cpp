#include "test_syncv_util.hpp"

#include <stdio.h>
#include <stdint.h>

namespace {
constexpr uint32_t sig_generic = 1u;
constexpr uint32_t sig_async = 2u;

camspork::ThreadCuboid simple_cuboid(uint32_t dim, uint32_t tid_lo, uint32_t tid_hi)
{
    auto cuboid = camspork::ThreadCuboid::full(&dim, 1 + &dim);
    cuboid.box()[0] = tid_hi - tid_lo;
    cuboid.offset()[0] = tid_lo;
    return cuboid;
}

}

using namespace camspork;

void test_sync_env()
{
    const bool do_validate = false;
    constexpr uint32_t warp_count = 32;

    syncv_init_t init{};
    init.filename = __FILE__;
    init.memory_budget = 1u << 30;
    SyncvTable_unique_ptr table_unique_ptr(new_syncv_table(init));

    SyncvTable* table = table_unique_ptr.get();

    assignment_record_id foobar[32 * warp_count] = {};
    debug_register_records(table, 32 * warp_count, foobar);

    auto maybe_validate = [do_validate] (SyncvTable* table)
    {
        if (do_validate) {
            debug_validate_state(table);
        }
    };

    for (uint32_t w = 0; w < warp_count; ++w) {
        for (uint32_t tid = w * 32; tid < w * 32 + 32; ++tid) {
            ThreadCuboid cuboid = simple_cuboid(32 * warp_count, tid, tid+1);
            on_rw(table, 1, &foobar[tid], cuboid, sig_generic | TlSigInterval::ordered_bits);
        }

        for (uint32_t tid = w * 32; tid < w * 32 + 32; ++tid) {
            ThreadCuboid cuboid = simple_cuboid(32 * warp_count, tid, tid+1);
            on_r(table, 1, &foobar[tid], cuboid, sig_generic | TlSigInterval::ordered_bits);
        }
    }

    for (int i = 0; i < 3; ++i) {
        // Simple barrier
        if (true) {
            ThreadCuboid cuboid = simple_cuboid(32 * warp_count, 0, 32 * warp_count);
            on_fence(table, true, cuboid, sig_generic, sig_generic, sig_generic);
        }

        maybe_validate(table);

        for (uint32_t tid = 0; tid < 32 * warp_count; ++tid) {
            ThreadCuboid cuboid = simple_cuboid(32 * warp_count, tid, tid+1);
            for (uint32_t i = 0; i < 32 * warp_count; ++i) {
                on_r(table, 1, &foobar[i], cuboid, sig_generic | TlSigInterval::ordered_bits);
            }

            if (tid > 3 && tid < 33) {
                maybe_validate(table);
            }
        }

        if (true) {
            ThreadCuboid cuboid = simple_cuboid(32 * warp_count, 0, 32 * warp_count);
            const auto bits = sig_generic | sig_async;
            on_fence(table, true, cuboid, bits, bits, bits);
        }

        maybe_validate(table);

        // This should fail if the above barrier is skipped (WAR)
        for (uint32_t tid = 0; tid < 32 * warp_count; ++tid) {
            ThreadCuboid cuboid = simple_cuboid(32 * warp_count, tid, tid+1);
            on_rw(table, 1, &foobar[tid], cuboid, sig_generic | TlSigInterval::ordered_bits);
        }

        maybe_validate(table);
    }

    // Deep copy (unique_ptr deletes original)
    table_unique_ptr.reset(copy_syncv_table(table));
    table = table_unique_ptr.get();

    if (true) {
        ThreadCuboid cuboid = simple_cuboid(32 * warp_count, 0, warp_count);  // ???
        on_fence(table, false, cuboid, sig_generic, sig_async, sig_async);
    }

    barrier_id bar{};
    alloc_barriers(table, 1, &bar);

    maybe_validate(table);

    // // Producer
    // {
    //     ThreadCuboid cuboid = simple_cuboid(32 * warp_count, 0, 32);
    //     // on_rw(table, 32, &foobar[0], V1);
    //     on_rw(table, 32, &foobar[0], cuboid, sig_async);
    //     on_arrive(table, &bar, V1, false);
    //     maybe_validate(table);
    // }

    // // Consumer
    // {
    //     maybe_validate(table);
    //     TlSigInterval V2{32, 64, sig_generic | TlSigInterval::ordered_bits};
    //     on_await(table, &bar, V2, V2);
    //     maybe_validate(table);
    //     on_r(table, 32, &foobar[0], V2);
    //     on_rw(table, 32, &foobar[0], V2);
    //     maybe_validate(table);
    // }


    free_barriers(table, 1, &bar);
    clear_visibility(table, 32 * warp_count, foobar);
    maybe_validate(table);
    debug_unregister_records(table, 32 * warp_count, foobar);
    debug_pre_delete_check(table);
}

