#include "exec.hpp"

#include <sstream>
#include <stdio.h>
#include <utility>

#include "grammar.hpp"

namespace camspork
{

class SwapThreadCuboid
{
    ThreadCuboid saved;
    ThreadCuboid* p_restore;
  public:
    // Sets *p_cuboid = new_cuboid, and restores old value upon destruction.
    [[nodiscard]] SwapThreadCuboid(ThreadCuboid* p_cuboid, ThreadCuboid new_cuboid)
    {
        saved = *p_cuboid;
        p_restore = p_cuboid;
        *p_cuboid = new_cuboid;
    };

    [[nodiscard]] SwapThreadCuboid(ThreadCuboid* p_cuboid)
    {
        saved = *p_cuboid;
        p_restore = p_cuboid;
    };

    SwapThreadCuboid(SwapThreadCuboid&&) = delete;

    ~SwapThreadCuboid()
    {
        *p_restore = saved;
    }
};

// Borrowed reference wrapper around ProgramEnv, to implement actual per-node-type execution.
class ProgramExec
{
    size_t buffer_size;
    const char* p_buffer;
    ProgramEnv& env;

  public:
    ProgramExec(ProgramEnv* p_self)
      : buffer_size(p_self->program_buffer_size)
      , p_buffer(p_self->p_program_buffer.get())
      , env(*p_self)
    {
    }

    // ******************************************************************************************
    // Many nodes define array indices as a VLA of ExprRef.
    // We provide a stripped-down iterator over these exprs, evaluated as values.
    // ******************************************************************************************
    struct ExprIterator
    {
        const ExprRef* p_node_ref;
        const ProgramExec* p_exec;

        intptr_t operator-(ExprIterator other) const
        {
            return p_node_ref - other.p_node_ref;
        }

        ExprIterator operator+(intptr_t i) const
        {
            return ExprIterator{p_node_ref + i, p_exec};
        }

        ExprIterator& operator++ ()
        {
            p_node_ref++;
            return *this;
        }

        value_t operator* () const
        {
            return p_exec->eval(*p_node_ref);
        }

        bool operator==(ExprIterator other) const
        {
            return p_node_ref == other.p_node_ref;
        }

        bool operator!=(ExprIterator other) const
        {
            return p_node_ref != other.p_node_ref;
        }
    };

    template <typename Node>
    ExprIterator expr_vla_begin(const Node* node) const
    {
        return ExprIterator{&node_vla_get_unsafe(node, 0), this};
    };

    template <typename Node>
    ExprIterator expr_vla_end(const Node* node) const
    {
        return expr_vla_begin(node) + node->camspork_vla_size;
    };

    template <typename Node>
    std::vector<extent_t> eval_extent(const Node* node) const
    {
        ExprIterator e_begin = expr_vla_begin(node);
        ExprIterator e_end = expr_vla_end(node);
        std::vector<extent_t> extent;
        extent.reserve(e_end - e_begin);
        for (ExprIterator e = e_begin; e != e_end; ++e) {
            extent.push_back(extent_t(*e));
        }
        return extent;
    }

    // ******************************************************************************************
    // EXECUTE STATEMENT
    // ******************************************************************************************
    __attribute__((always_inline))
    void exec(StmtRef s)
    {
        if (s) {
            s.dispatch(*this, buffer_size, p_buffer);
        }
    }

    template <bool IsOOO, bool IsMutate>
    void operator() (const SyncEnvAccessNode<IsOOO, IsMutate>*)
    {
        CAMSPORK_REQUIRE(0, "TODO: implement SyncEnvAccessNode");
    }

    void operator() (const MutateValue* node)
    {
        value_t& lhs = env.value_slot(node->name).idx(expr_vla_begin(node), expr_vla_end(node));
        const value_t rhs = eval(node->rhs);
        lhs = eval_binop(node->op, lhs, rhs);
    }

    void operator() (const Fence*)
    {
        CAMSPORK_REQUIRE(0, "TODO: implement Fence");
    }

    void operator() (const Arrive*)
    {
        CAMSPORK_REQUIRE(0, "TODO: implement Arrive");
    }

    void operator() (const Await*)
    {
        CAMSPORK_REQUIRE(0, "TODO: implement Await");
    }

    void operator() (const ValueEnvAlloc* node)
    {
        env.value_slot(node->name) = VarSlotEntry<value_t>(eval_extent(node));
    }

    void operator() (const SyncEnvAlloc* node)
    {
        env.sync_slot(node->name) = VarSlotEntry<assignment_record_id>(eval_extent(node));
    }

    void operator() (const SyncEnvFreeShard*)
    {
        CAMSPORK_REQUIRE(0, "TODO: implement SyncEnvFreeShard");
    }

    void operator() (const BarrierEnvAlloc*)
    {
        CAMSPORK_REQUIRE(0, "TODO: implement BarrierEnvAlloc");
    }

    void operator() (const BarrierEnvFree*)
    {
        CAMSPORK_REQUIRE(0, "TODO: implement BarrierEnvFree");
    }

    struct BodyExecImpl
    {
        ProgramExec& exec;
        uint32_t stmts_left;
        const StmtRef* p_stmts;

        template <typename Stmt>
        void operator() (const Stmt* node)
        {
            exec(node);
            stmts_left--;
            p_stmts++;
            if (stmts_left <= 0) {
                return;
            }
            // Tail call to execute next stmt.
            // Since this is per-Stmt type, the branch predictor may learn correlations
            // of what the next statement type is with respect to this statement type.
            p_stmts->dispatch(*this, exec.buffer_size, exec.p_buffer);
        }
    };

    void operator() (const StmtBody* node)
    {
        uint32_t num_stmts = node->camspork_vla_size;
        if (num_stmts > 0) {
            const StmtRef* p_stmts = &node_vla_get(node, 0);
            BodyExecImpl impl{*this, num_stmts, p_stmts};
            p_stmts->dispatch(impl, buffer_size, p_buffer);
        }
    }

    void operator() (const If* node)
    {
        StmtRef s = eval(node->cond) ? node->body : node->orelse;
        exec(s);  // Inlined only once
    }

    template <typename Step>
    void exec_for_impl(const BaseForStmt* node, Step&& step)
    {
        const auto lo = eval(node->lo);
        const auto hi = eval(node->hi);
        env.alloc_scalar_value(node->iter, lo);
        for (value_t i = lo; i < hi; ++i) {
            // Look up Varslot each time in case the loop body did something bad!
            env.value_slot(node->iter).scalar() = i;
            exec(node->body);
            step();
        }
    }

    void operator() (const SeqFor* node)
    {
        exec_for_impl(node, [] {});
    }

    void operator() (const TasksFor* node)
    {
        auto* p_task_index = &env.thread_cuboid.task_index;
        auto step = [p_task_index] {++*p_task_index;};
        exec_for_impl(node, step);
    }

    void operator() (const ThreadsFor* node)
    {
        const uint32_t dim_idx = node->dim_idx;

        const uint32_t offset_c = node->offset;
        const uint32_t box_c = node->box;
        const auto lo = eval(node->lo);
        const auto hi = eval(node->hi);
        env.alloc_scalar_value(node->iter, lo);

        // This shouldn't hard to change, but just test it quickly if you change this.
        CAMSPORK_REQUIRE_CMP(lo, ==, 0, "Expected ThreadsFor loop to start from 0 for now");

        // Restores thread cuboid before returning.
        SwapThreadCuboid swap(&env.thread_cuboid);

        CAMSPORK_REQUIRE_CMP(dim_idx, <, env.thread_cuboid.dim, "ThreadsFor::dim_idx out of range");
        CAMSPORK_REQUIRE_CMP(offset_c + (hi - lo) * box_c, <=, env.thread_cuboid.box[dim_idx],
                             "ThreadsFor consumes more threads than exists in the current thread box");

        env.thread_cuboid.offset[dim_idx] += offset_c;
        env.thread_cuboid.box[dim_idx] = box_c;

        for (value_t i = lo; i < hi; ++i) {
            if (true) {
                const char* var_c_name = env.var_slots[node->iter.slot].name.c_str();
                printf("%s = %i, %s\n", var_c_name, i, (std::stringstream() << env.thread_cuboid).str().c_str());
            }

            // Look up Varslot each time in case the loop body did something bad!
            env.value_slot(node->iter).scalar() = i;
            exec(node->body);

            // Slide thread box over for the next iteration.
            env.thread_cuboid.offset[dim_idx] += box_c;
        }
    }

    void operator() (const ParallelBlock* node)
    {
        ThreadCuboid new_cuboid;
        const uint32_t dim = node->camspork_vla_size;
        CAMSPORK_REQUIRE_CMP(dim, <= , ThreadCuboid::max_dim, "implementation limit: domain too big");
        new_cuboid.dim = dim;
        for (uint32_t i = 0; i < dim; ++i) {
            const auto dim_extent = node_vla_get(node, i);
            new_cuboid.domain[i] = dim_extent;
            new_cuboid.offset[i] = 0;
            new_cuboid.box[i] = dim_extent;
        }

        // Execute body with new thread cuboid, and restore before returning (~SwapThreadCuboid).
        SwapThreadCuboid swap(&env.thread_cuboid, new_cuboid);
        exec(node->body);
    }

    void operator() (const DomainSplit* node)
    {
        ThreadCuboid new_cuboid = env.thread_cuboid;
        const uint32_t split_idx = node->dim_idx;
        const uint32_t split_factor = node->split_factor;
        CAMSPORK_REQUIRE_CMP(split_idx, <, new_cuboid.dim, "out-of-range DomainSplit::dim_idx");
        CAMSPORK_REQUIRE_CMP(split_factor, >=, 1, "invalid DomainSplit::split_factor");

        const uint32_t domain_c = new_cuboid.domain[split_idx];
        if (domain_c == split_factor || split_factor == 1) {
            // Unchanged.
        }
        else {
            const uint32_t new_dim = new_cuboid.dim + 1;
            CAMSPORK_REQUIRE_CMP(new_dim, <=, ThreadCuboid::max_dim, "implementation limit: domain too big");
            const uint32_t offset_c = new_cuboid.offset[split_idx];
            const uint32_t box_c = new_cuboid.box[split_idx];
            CAMSPORK_REQUIRE_CMP(domain_c % split_factor, ==, 0, "Invalid DomainSplit::split_factor for current env");
            CAMSPORK_REQUIRE_CMP(offset_c % split_factor, ==, 0, "Invalid DomainSplit::split_factor for current env");

            const uint32_t offset_0 = offset_c / split_factor;
            const uint32_t offset_1 = 0;
            const uint32_t domain_0 = domain_c / split_factor;
            const uint32_t domain_1 = split_factor;
            uint32_t box_0, box_1;
            if (box_c < domain_c) {
                box_0 = 1;
                box_1 = box_c;
            }
            else {
                CAMSPORK_REQUIRE_CMP(box_c % split_factor, ==, 0, "Invalid DomainSplit::split_factor for current env");
                box_0 = box_c / split_factor;
                box_1 = split_factor;
            }

            // Insert new domain/offset/box coordinates in place of the old ones.
            new_cuboid.dim = new_dim;
            for (uint32_t dst = new_dim - 1; dst > split_idx; --dst) {
                new_cuboid.domain[dst] = new_cuboid.domain[dst - 1];
                new_cuboid.offset[dst] = new_cuboid.offset[dst - 1];
                new_cuboid.box[dst] = new_cuboid.box[dst - 1];
            }
            new_cuboid.domain[split_idx] = domain_0;
            new_cuboid.domain[split_idx + 1] = domain_1;
            new_cuboid.offset[split_idx] = offset_0;
            new_cuboid.offset[split_idx + 1] = offset_1;
            new_cuboid.box[split_idx] = box_0;
            new_cuboid.box[split_idx + 1] = box_1;
        }

        // Execute body with new thread cuboid, and restore before returning (~SwapThreadCuboid).
        SwapThreadCuboid swap(&env.thread_cuboid, new_cuboid);
        exec(node->body);
    }

    // ******************************************************************************************
    // EVALUATE EXPR
    // ******************************************************************************************
    __attribute__((always_inline))
    value_t eval(ExprRef e) const
    {
        return e.dispatch(*this, buffer_size, p_buffer);
    }

    __attribute__((always_inline))
    value_t operator() (const ReadValue* node) const
    {
        return env.value_slot(node->name).idx(expr_vla_begin(node), expr_vla_end(node));
    }

    __attribute__((always_inline))
    value_t operator() (const Const* node) const
    {
        return node->value;
    }

    value_t operator() (const USub* node) const
    {
        return -eval(node->arg);
    }

    value_t operator() (const BinOp* node) const
    {
        CAMSPORK_REQUIRE_CMP(int(node->op), !=, int(binop::Assign), "binop::Assign is only allowed in MutateValue");
        return eval_binop(node->op, eval(node->lhs), eval(node->rhs));
    }

    static value_t eval_binop(binop op, value_t lhs, value_t rhs)
    {
        switch (op) {
          case binop::Assign:
            return rhs;
          case binop::Add:
            return lhs + rhs;
          case binop::Sub:
            return lhs - rhs;
          case binop::Mul:
            return lhs * rhs;
          case binop::Div:
            // Python-style division
            CAMSPORK_REQUIRE_CMP(rhs, >, 0, "Can only divide by positive numbers");
            return (lhs < 0 ? lhs + rhs - 1 : lhs) / rhs;
          case binop::Mod:
            CAMSPORK_REQUIRE_CMP(rhs, >, 0, "Can only modulo by positive numbers");
            return (lhs < 0 ? lhs + rhs - 1 : lhs) % rhs;
          case binop::Less:
            return lhs < rhs;
          case binop::Leq:
            return lhs <= rhs;
          case binop::Greater:
            return lhs > rhs;
          case binop::Geq:
            return lhs >= rhs;
          case binop::Eq:
            return lhs == rhs;
          case binop::Neq:
            return lhs != rhs;
        }
        return 0;  // XXX should do something
    }

    // ******************************************************************************************
    // First time startup, initialize variable table
    // ******************************************************************************************
    void init_vars(VarConfigTableRef table)
    {
        table.dispatch(*this, buffer_size, p_buffer);
    }

    void operator() (const VarConfigTable* table)
    {
        // Initialize variable tables, then iterate over the variable length array
        // to initialize all the variable slots.
        const auto num_slots = table->camspork_vla_size;
        for (uint32_t i = 0; i < num_slots; ++i) {
            VarConfigRef config = node_vla_get(table, i);
            config.dispatch(*this, buffer_size, p_buffer);
        }
    }

    void operator() (const VarConfig* config)
    {
        const auto num_bytes = config->camspork_vla_size;
        const char* p_str = &node_vla_get(config, 0);
        env.var_slots.push_back({std::string(p_str, p_str + num_bytes), {}, {}, {}});
    }
};



ProgramEnv::ProgramEnv(size_t buffer_size, const char* buffer)
  : program_buffer_size(buffer_size)
  , p_program_buffer(make_shared_program_buffer(buffer_size, buffer))
  , header(ProgramHeader::validate(buffer_size, p_program_buffer.get()))
{
    ProgramExec(this).init_vars(header.var_config_table);
};

void ProgramEnv::exec(StmtRef stmt)
{
    ProgramExec(this).exec(stmt);
}

}
