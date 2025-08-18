#include "exec.hpp"

#include <utility>

#include "grammar.hpp"

namespace camspork
{

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
    }

    void operator() (const MutateValue* node)
    {
        value_t& lhs = env.value_slot(node->name).idx(expr_vla_begin(node), expr_vla_end(node));
        const value_t rhs = eval(node->rhs);
        lhs = eval_binop(node->op, lhs, rhs);
    }

    void operator() (const Fence*)
    {
    }

    void operator() (const Arrive*)
    {
    }

    void operator() (const Await*)
    {
    }

    void operator() (const ValueEnvAlloc* node)
    {
        env.value_slot(node->name) = VarSlotEntry<value_t>(eval_extent(node));
    }

    void operator() (const SyncEnvAlloc*)
    {
    }

    void operator() (const SyncEnvFreeShard*)
    {
    }

    void operator() (const BarrierEnvAlloc*)
    {
    }

    void operator() (const BarrierEnvFree*)
    {
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

    template <typename Init, typename Step>
    void exec_for_impl(const BaseForStmt* node, Init&& init, Step&& step)
    {
        const auto lo = eval(node->lo);
        const auto hi = eval(node->hi);
        env.alloc_scalar_value(node->iter, lo);
        init();
        for (value_t i = lo; i < hi; ++i) {
            // Look up Varslot each time in case the loop body did something bad!
            env.value_slot(node->iter).scalar() = i;
            exec(node->body);
            step();
        }
    }

    void operator() (const SeqFor* node)
    {
        exec_for_impl(node, [] {}, [] {});
    }

    void operator() (const TasksFor* node)
    {
        value_t* p_task_index = &env.task_index;
        exec_for_impl(node, [] {}, [p_task_index] {++*p_task_index;});
    }

    void operator() (const ThreadsFor*)
    {
    }

    void operator() (const DomainDefine*)
    {
    }

    void operator() (const DomainSplit*)
    {
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
        // to initialize all the variable names.
        const auto num_slots = table->camspork_vla_size;
        // TODO use just one vector instead of parallel vectors?
        env.value_env_slots.resize(num_slots);
        env.sync_env_slots.resize(num_slots);
        env.barrier_env_slots.resize(num_slots);
        env.variable_names.reserve(num_slots);  // Will be emplaced_back below.
        for (uint32_t i = 0; i < num_slots; ++i) {
            VarConfigRef config = node_vla_get(table, i);
            config.dispatch(*this, buffer_size, p_buffer);
        }
    }

    void operator() (const VarConfig* config)
    {
        const auto num_bytes = config->camspork_vla_size;
        const char* p_str = &node_vla_get(config, 0);
        env.variable_names.emplace_back(p_str, p_str + num_bytes);
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
