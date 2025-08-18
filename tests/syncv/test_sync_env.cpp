#include "test_syncv_util.hpp"

#include <stdio.h>

namespace {
constexpr uint32_t sig_generic = 1u;
constexpr uint32_t sig_async = 2u;
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
            SigthreadInterval accessor_set{tid, tid+1, sig_generic | SigthreadInterval::sync_bit};
            on_rw(table, 1, &foobar[tid], accessor_set);
        }

        for (uint32_t tid = w * 32; tid < w * 32 + 32; ++tid) {
            SigthreadInterval accessor_set{tid, tid+1, sig_generic | SigthreadInterval::sync_bit};
            on_r(table, 1, &foobar[tid], accessor_set);
        }
    }

    for (int i = 0; i < 3; ++i) {
        // Simple barrier
        if (true) {
            SigthreadInterval V1{0, warp_count * 32, sig_generic};
            SigthreadInterval V2{0, warp_count * 32, sig_generic};
            on_fence(table, V1, V2, V2, true);
        }

        maybe_validate(table);

        for (uint32_t tid = 0; tid < 32 * warp_count; ++tid) {
            SigthreadInterval accessor_set{tid, tid+1, sig_generic | SigthreadInterval::sync_bit};
            for (uint32_t i = 0; i < 32 * warp_count; ++i) {
                on_r(table, 1, &foobar[i], accessor_set);
            }

            if (tid > 3 && tid < 33) {
                maybe_validate(table);
            }
        }

        if (true) {
            SigthreadInterval V{0, warp_count * 32, sig_generic | sig_async};
            on_fence(table, V, V, V, true);
        }

        maybe_validate(table);

        // This should fail if the above barrier is skipped (WAR)
        for (uint32_t tid = 0; tid < 32 * warp_count; ++tid) {
            SigthreadInterval accessor_set{tid, tid+1, sig_generic | SigthreadInterval::sync_bit};
            on_rw(table, 1, &foobar[tid], accessor_set);
        }

        maybe_validate(table);
    }

    // Deep copy (unique_ptr deletes original)
    table_unique_ptr.reset(copy_syncv_table(table));
    table = table_unique_ptr.get();

    if (true) {
        SigthreadInterval V1{0, warp_count, sig_generic};
        SigthreadInterval V2{0, warp_count, sig_async};
        on_fence(table, V1, V2, V2, false);
    }

    barrier_id bar{};
    alloc_barriers(table, 1, &bar);

    maybe_validate(table);

    // Producer
    {
        SigthreadInterval V1{0, 32, sig_async};
        // on_rw(table, 32, &foobar[0], V1);
        on_rw(table, 32, &foobar[0], V1);
        on_arrive(table, &bar, V1, false);
        maybe_validate(table);
    }

    // Consumer
    {
        maybe_validate(table);
        SigthreadInterval V2{32, 64, sig_generic | SigthreadInterval::sync_bit};
        on_await(table, &bar, V2, V2);
        maybe_validate(table);
        on_r(table, 32, &foobar[0], V2);
        on_rw(table, 32, &foobar[0], V2);
        maybe_validate(table);
    }


    free_barriers(table, 1, &bar);
    clear_visibility(table, 32 * warp_count, foobar);
    maybe_validate(table);
    debug_unregister_records(table, 32 * warp_count, foobar);
    debug_pre_delete_check(table);
}

