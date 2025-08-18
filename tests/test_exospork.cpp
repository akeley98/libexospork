#include "syncv/test_sync_env.hpp"

#include "../lib/program/builder.hpp"
#include "../lib/program/exec.hpp"

namespace camspork
{

void Main()
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

}

int main()
{
    test_sync_env();
    camspork::Main();
    return 0;
}
