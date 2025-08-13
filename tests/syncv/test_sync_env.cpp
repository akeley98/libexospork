#include "test_syncv_util.hpp"

#include <stdio.h>

namespace {
constexpr uint32_t sig_generic = 1u;
constexpr uint32_t sig_async = 2u;
}

void test_sync_env()
{
    constexpr uint32_t warp_count = 32;

    exospork_syncv_init_t init{};
    init.filename = __FILE__;
    init.memory_budget = 1u << 30;
    exospork::SyncEnv_unique_ptr env_unique_ptr(exospork::new_sync_env(init));

    exospork::SyncEnv* p_env = env_unique_ptr.get();

    exospork_syncv_value_t values[32 * warp_count] = {};
    debug_register_values(p_env, 32 * warp_count, values);

    for (uint32_t w = 0; w < warp_count; ++w) {
        for (uint32_t tid = w * 32; tid < w * 32 + 32; ++tid) {
            exospork::SigthreadInterval accessor_set{tid, tid+1, sig_generic | EXOSPORK_SYNC_ACCESS_BIT};
            on_rw(p_env, 1, &values[tid], accessor_set);
        }

        for (uint32_t tid = w * 32; tid < w * 32 + 32; ++tid) {
            exospork::SigthreadInterval accessor_set{tid, tid+1, sig_generic | EXOSPORK_SYNC_ACCESS_BIT};
            on_r(p_env, 1, &values[tid], accessor_set);
        }
    }

    for (int i = 0; i < 3; ++i) {
        // Simple barrier
        if (true) {
            exospork::SigthreadInterval V1{0, warp_count * 32, sig_generic};
            exospork::SigthreadInterval V2{0, warp_count * 32, sig_generic};
            on_fence(p_env, V1, V2, V2, true);
        }

        debug_validate_state(p_env);

        for (uint32_t tid = 0; tid < 32 * warp_count; ++tid) {
            exospork::SigthreadInterval accessor_set{tid, tid+1, sig_generic | EXOSPORK_SYNC_ACCESS_BIT};
            for (uint32_t i = 0; i < 32 * warp_count; ++i) {
                on_r(p_env, 1, &values[i], accessor_set);
            }

            if (tid > 3 && tid < 33) {
                debug_validate_state(p_env);
            }
        }

        if (true) {
            exospork::SigthreadInterval V{0, warp_count * 32, sig_generic | sig_async};
            on_fence(p_env, V, V, V, true);
        }

        debug_validate_state(p_env);

        // This should fail if the above barrier is skipped (WAR)
        for (uint32_t tid = 0; tid < 32 * warp_count; ++tid) {
            exospork::SigthreadInterval accessor_set{tid, tid+1, sig_generic | EXOSPORK_SYNC_ACCESS_BIT};
            on_rw(p_env, 1, &values[tid], accessor_set);
        }

        debug_validate_state(p_env);
    }

    // Deep copy (unique_ptr deletes original)
    env_unique_ptr.reset(copy_sync_env(p_env));
    p_env = env_unique_ptr.get();

    if (true) {
        exospork::SigthreadInterval V1{0, warp_count, sig_generic};
        exospork::SigthreadInterval V2{0, warp_count, sig_async};
        on_fence(p_env, V1, V2, V2, false);
    }

    exospork_syncv_barrier_t bar{};
    alloc_barrier(p_env, &bar);

    debug_validate_state(p_env);

    // Producer
    {
        exospork::SigthreadInterval V1{0, 32, sig_async};
        // on_rw(p_env, 32, &values[0], V1);
        on_rw(p_env, 32, &values[0], V1);
        on_arrive(p_env, &bar, V1, false);
        debug_validate_state(p_env);
    }

    // Consumer
    {
        debug_validate_state(p_env);
        exospork::SigthreadInterval V2{32, 64, sig_generic | EXOSPORK_SYNC_ACCESS_BIT};
        on_await(p_env, &bar, V2, V2);
        debug_validate_state(p_env);
        on_r(p_env, 32, &values[0], V2);
        on_rw(p_env, 32, &values[0], V2);
        debug_validate_state(p_env);
    }


    free_barrier(p_env, &bar);
    clear_values(p_env, 32 * warp_count, values);
    debug_validate_state(p_env);
    debug_unregister_values(p_env, 32 * warp_count, values);
    debug_pre_delete_check(p_env);
}
