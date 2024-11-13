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
        {
            exospork::SigthreadInterval V1{0, warp_count * 32, sig_generic};
            exospork::SigthreadInterval V2{0, warp_count * 32, sig_generic};
            on_fence(p_env, V1, V2, true);
        }

        for (uint32_t tid = 0; tid < 32 * warp_count; ++tid) {
            exospork::SigthreadInterval accessor_set{tid, tid+1, sig_generic | EXOSPORK_SYNC_ACCESS_BIT};
            for (uint32_t i = 0; i < 32 * warp_count; ++i) {
                on_r(p_env, 1, &values[i], accessor_set);
            }
        }

        if (true) {
            exospork::SigthreadInterval V1{0, warp_count * 32, sig_generic};
            exospork::SigthreadInterval V2{0, warp_count * 32, sig_generic};
            on_fence(p_env, V1, V2, true);
        }

        // This should fail if the above barrier is skipped (WAR)
        for (uint32_t tid = 0; tid < 32 * warp_count; ++tid) {
            exospork::SigthreadInterval accessor_set{tid, tid+1, sig_generic | EXOSPORK_SYNC_ACCESS_BIT};
            on_rw(p_env, 1, &values[tid], accessor_set);
        }
    }

    clear_values(p_env, 32 * warp_count, values);
    debug_unregister_values(p_env, 32 * warp_count, values);
}
