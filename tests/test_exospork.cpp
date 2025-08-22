#include "syncv/test_sync_env.hpp"
#include "test_cuboid_util.hpp"

#include "../lib/program/builder.hpp"
#include "../lib/program/exec.hpp"

namespace camspork
{

void fib_program_test()
{
    const uint32_t c_fib_size = 20;
    ProgramBuilder builder;
    ExprRef _0 = builder.add_Const(0);
    ExprRef _1 = builder.add_Const(1);
    ExprRef fib_size = builder.add_Const(c_fib_size);
    Varname fib = builder.add_variable("fib");
    Varname iter = builder.add_variable("iter");

    builder.add_ValueEnvAlloc(fib, 1, &fib_size);
    builder.add_MutateValue(fib, 1, &_0, binop::Assign, _0);
    builder.add_MutateValue(fib, 1, &_1, binop::Assign, _1);
    ExprRef read_iter = builder.add_ReadValue(iter, 0, nullptr);

    builder.push_SeqFor(iter, builder.add_Const(2), fib_size);
    {
        ExprRef iter_1 = builder.add_BinOp(binop::Sub, read_iter, _1);
        ExprRef iter_2 = builder.add_BinOp(binop::Sub, read_iter, builder.add_Const(2));
        ExprRef e1 = builder.add_ReadValue(fib, 1, &iter_1);
        ExprRef e2 = builder.add_ReadValue(fib, 1, &iter_2);
        builder.add_MutateValue(fib, 1, &read_iter, binop::Assign, builder.add_BinOp(binop::Add, e1, e2));
    }
    builder.pop_body();

    Varname dst = builder.add_variable("dst");
    builder.add_ValueEnvAlloc(dst, 1, &fib_size);
    builder.push_SeqFor(iter, _0, fib_size);
    {
        ExprRef read_fib = builder.add_ReadValue(fib, 1, &read_iter);
        ExprRef cond = builder.add_BinOp(binop::Mod, read_fib, builder.add_Const(5));
        builder.push_If(cond);
        {
            // Not divisible by 5
            builder.add_MutateValue(fib, 1, &read_iter, binop::Assign, builder.add_USub(read_fib));
            builder.add_MutateValue(fib, 1, &read_iter, binop::Mul, builder.add_Const(10000));
            // Is divisible by 5
            builder.begin_orelse();
            builder.add_MutateValue(fib, 1, &read_iter, binop::Div, builder.add_Const(5));
        }
        builder.pop_body();
    }
    builder.pop_body();

    builder.finish();

    ProgramEnv env(builder.size(), builder.data());
    env.exec();

    // FILE* f = fopen("test_exospork.camspork", "wb");
    // fwrite(builder.data(), 1, builder.size(), f);
    // fclose(f);

    auto fib_value_slot = env.value_slot(fib);
    for (uint32_t i = 0; i < c_fib_size; ++i) {
        printf("[%u] %i\n", i, fib_value_slot.idx(&i, 1 + &i));
    }
}

void threads_program_test()
{
    ProgramBuilder builder;
    const ExprRef _0 = builder.add_Const(0);

    const uint32_t initial_domain[2] = {8, 384};
    builder.push_ParallelBlock(2, initial_domain);
    {
        const Varname task_v = builder.add_variable("task");
        builder.push_TasksFor(task_v, _0, builder.add_Const(2));
        {
            // Loop over warpgroups (all CTAs active)
            const Varname warpgroup_v = builder.add_variable("warpgroup");
            builder.push_ThreadsFor(warpgroup_v, _0, builder.add_Const(3), 1, 0, 128);
            {
                // Split cluster into 4 x 2 grid of CTAs
                builder.push_DomainSplit(0, 2);
                {
                    // Loop over fast CTA dimension.
                    const Varname cta_slow_v = builder.add_variable("cta_slow");
                    const Varname cta_fast_v = builder.add_variable("cta_fast");
                    builder.push_ThreadsFor(cta_fast_v, _0, builder.add_Const(2), 1, 0, 1);
                    {
                        // Loop over slow CTA dimension.
                        builder.push_ThreadsFor(cta_slow_v, _0, builder.add_Const(4), 0, 0, 1);
                        {
                            // Loop over warps in warpgroup
                            const Varname warp_v = builder.add_variable("warp");
                            builder.push_ThreadsFor(warp_v, _0, builder.add_Const(4), 2, 0, 32);
                            {
                            }
                            builder.pop_body();
                        }
                        builder.pop_body();
                    }
                    builder.pop_body();
                }
                builder.pop_body();
            }
            builder.pop_body();

            // Loop over CTAs linearly.
            const Varname cta_v = builder.add_variable("cta");
            const Varname warp_lo_v = builder.add_variable("warp_lo");
            const Varname warpgroup_hi_v = builder.add_variable("warpgroup_hi");
            builder.push_ThreadsFor(cta_v, _0, builder.add_Const(8), 0, 0, 1);
            {
                // Loop over warps [0, 3]
                builder.push_ThreadsFor(warp_lo_v, _0, builder.add_Const(4), 1, 0, 32);
                {
                }
                builder.pop_body();
                // Loop over warpgroups [4, 11]
                builder.push_ThreadsFor(warpgroup_hi_v, _0, builder.add_Const(2), 1, 128, 128);
                {
                }
                builder.pop_body();
            }
            builder.pop_body();
        }
        // End tasks loop
        builder.pop_body();
    }
    // End ParallelBlock
    builder.pop_body();

    builder.finish();
    ProgramEnv env(builder.size(), builder.data());
    env.exec();
}

void simple_fence_test(uint32_t n_tasks, bool have_fence)
{
    ProgramBuilder builder;
    const ExprRef _0 = builder.add_Const(0);
    const ExprRef _1 = builder.add_Const(1);
    const ExprRef _256 = builder.add_Const(256);

    const Varname thread_v = builder.add_variable("thread");
    const Varname gmem_v = builder.add_variable("gmem");
    builder.add_SyncEnvAlloc(gmem_v, 1, &_256);

    const uint32_t initial_domain[1] = {256};
    builder.push_ParallelBlock(1, initial_domain);
    {
        const Varname task_v = builder.add_variable("task");
        builder.push_TasksFor(task_v, _0, builder.add_Const(n_tasks));
        {
            const Varname smem_v = builder.add_variable("smem");
            builder.add_SyncEnvAlloc(smem_v, 1, &_256);
            builder.push_ThreadsFor(thread_v, _0, _256, 0, 0, 1);
            {
                OffsetExtentExpr idx;
                idx.offset_e = builder.add_ReadValue(thread_v, 0, nullptr);
                idx.extent_e = _1;
                builder.add_SyncEnvAccess(smem_v, 1, &idx, true, false, 1, 1);
            }
            builder.pop_body();
            if (have_fence) {
                builder.add_Fence(true, 1, 1, 1);
            }
            builder.push_ThreadsFor(thread_v, _0, _256, 0, 0, 1);
            {
                const Varname i_v = builder.add_variable("i");
                builder.push_SeqFor(i_v, _0, _256);
                {
                    OffsetExtentExpr idx;
                    idx.offset_e = builder.add_ReadValue(i_v, 0, nullptr);
                    idx.extent_e = _1;
                    builder.add_SyncEnvAccess(smem_v, 1, &idx, false, false, 1, 1);
                }
                builder.pop_body();
                OffsetExtentExpr idx;
                idx.offset_e = builder.add_ReadValue(thread_v, 0, nullptr);
                idx.extent_e = _1;
                builder.add_SyncEnvAccess(gmem_v, 1, &idx, true, false, 1, 1);
            }
            builder.pop_body();
            if (have_fence) {
                builder.add_Fence(true, 1, 1, 1);
            }
        }
        builder.pop_body();
    }
    // End ParallelBlock
    builder.pop_body();

    builder.finish();

    const bool expected_error = n_tasks > 1 || (n_tasks >= 1 && !have_fence);
    bool had_error = false;
    ProgramEnv env(builder.size(), builder.data());
    try {
        env.exec();
    }
    catch (SyncvCheckFail& fail) {
        had_error = true;
    }
    CAMSPORK_REQUIRE_CMP(had_error, ==, expected_error, "test failed");
}

}

int main()
{
    camspork::simple_fence_test(1, true);
    camspork::simple_fence_test(2, true);
    camspork::simple_fence_test(1, false);
    camspork::threads_program_test();
    camspork::test_cuboid_util();
    test_sync_env();
    camspork::fib_program_test();
    return 0;
}
