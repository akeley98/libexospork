#pragma once

#include <string>
#include <vector>

#include "grammar.hpp"
#include "../util/require.hpp"

namespace camspork
{

class ProgramBuilder;

struct BodyBuilder
{
    // When dispatched, this sets the vector of stmts to be the body of the body_of statement,
    // unless is_orelse is true, then we set the vector to be the orelse of the body_of If statement.
    ProgramBuilder* p_program_builder;
    StmtRef body_of;
    std::vector<StmtRef> stmts;
    bool is_orelse = false;

    template <uint32_t StmtType>
    void operator() (stmt<StmtType>*) const
    {
        // Fallback for node types not specifically targetted below.
        CAMSPORK_REQUIRE_CMP(StmtType, ==, -1, "Internal error: invalid node type for BodyBuilder");
    }

    void operator() (If* node) const
    {
        StmtRef s = body_to_nursery();
        if (is_orelse) {
            node->orelse = s;
        }
        else {
            node->body = s;
        }
    }

    template <typename Node>
    void set_body_common(Node* node) const
    {
        StmtRef s = body_to_nursery();
        CAMSPORK_REQUIRE(!is_orelse, "Only If statements may have an orelse");
        node->body = s;
    }

    void operator() (SeqFor* node) const { set_body_common(node); }
    void operator() (TasksFor* node) const { set_body_common(node); }
    void operator() (ThreadsFor* node) const { set_body_common(node); }
    void operator() (ParallelBlock* node) const { set_body_common(node); }
    void operator() (DomainSplit* node) const { set_body_common(node); }

    void begin_orelse();
    StmtRef body_to_nursery() const;
};

class ProgramBuilder
{
    std::vector<std::string> variable_slot_names;
    NodeNursery nursery;
    // 0th entry in body_stack corresponds to the top-level program
    // Further levels are used while building If, For, etc.
    std::vector<BodyBuilder> body_stack;
    bool is_finished = false;

  public:
    ProgramBuilder();

    // Finalize the program. All push_* must have been paired with pop_body().
    void finish();

    size_t size() const
    {
        CAMSPORK_REQUIRE(is_finished, "call ProgramBuilder::finish()");
        return nursery.size();
    }

    const char* data() const
    {
        CAMSPORK_REQUIRE(is_finished, "call ProgramBuilder::finish()");
        return nursery.data();
    }

    // ******************************************************************************************
    // Add variables to the program.
    // ******************************************************************************************
    Varname add_variable(const char* name)
    {
        variable_slot_names.push_back(name);
        return Varname{uint32_t(variable_slot_names.size()) - 1};
    }

    // ******************************************************************************************
    // Add expressions to the program.
    // ******************************************************************************************
    ExprRef add_ReadValue(Varname name, size_t num_idx, const ExprRef* idx);
    ExprRef add_Const(value_t value);
    ExprRef add_USub(ExprRef arg);
    ExprRef add_BinOp(binop op, ExprRef lhs, ExprRef rhs);

    // ******************************************************************************************
    // Add statements that don't have a body to the program.
    // ******************************************************************************************
    StmtRef add_SyncEnvAccess(
        Varname name, size_t num_idx, const OffsetExtentExpr* idx,
        bool is_mutate, bool is_ooo, qual_bits_t initial_qual_bit, qual_bits_t extended_qual_bits);
    StmtRef add_MutateValue(Varname name, size_t num_idx, const ExprRef* idx, binop op, ExprRef rhs);
    StmtRef add_Fence(
        uint32_t V1_transitive, qual_bits_t L1_qual_bits,
        qual_bits_t L2_full_qual_bits, qual_bits_t L2_temporal_qual_bits);
    StmtRef add_ValueEnvAlloc(Varname name, size_t num_dims, const ExprRef* extent);
    StmtRef add_SyncEnvAlloc(Varname name, size_t num_dims, const ExprRef* extent);

    // ******************************************************************************************
    // Add statements with a body to the program.
    // Further statements go into the body of the new statement, until you call pop_body().
    // Use begin_orelse() to switch to adding statements to the orelse of an If statement.
    // ******************************************************************************************
    StmtRef push_If(ExprRef cond);
    void begin_orelse();
    StmtRef push_SeqFor(Varname iter, ExprRef lo, ExprRef hi);
    StmtRef push_TasksFor(Varname iter, ExprRef lo, ExprRef hi);
    StmtRef push_ThreadsFor(Varname iter, ExprRef lo, ExprRef hi, uint32_t dim_idx, uint32_t offset, uint32_t box);
    StmtRef push_ParallelBlock(size_t dim, const uint32_t* domain);
    StmtRef push_DomainSplit(uint32_t dim_idx, uint32_t split_factor);

  private:
    void check_not_finished() const;

    template <typename...Args>
    StmtRef append_impl(Args... a);

    template <typename...Args>
    StmtRef push_impl(Args... a);

  public:
    void pop_body();

    // For use by BodyBuilder.
    void end_body_builder(const BodyBuilder& body_builder);
    StmtRef body_to_nursery(const std::vector<StmtRef>& stmts);
};

}
