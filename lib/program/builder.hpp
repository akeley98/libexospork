#pragma once

#include <string>
#include <vector>

#include "grammar.hpp"
#include "../util/api_util.hpp"
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
        return Varname{uint32_t(variable_slot_names.size())};
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
        uint32_t is_mutate, uint32_t is_ooo, qual_bits_t initial_qual_bit, qual_bits_t extended_qual_bits);
    StmtRef add_MutateValue(Varname name, size_t num_idx, const ExprRef* idx, binop op, ExprRef rhs);
    StmtRef add_Fence(
        uint32_t V1_transitive, qual_bits_t L1_qual_bits,
        qual_bits_t L2_full_qual_bits, qual_bits_t L2_temporal_qual_bits);
    StmtRef add_ValueEnvAlloc(Varname name, size_t num_dims, const ExprRef* extent);
    StmtRef add_SyncEnvAlloc(Varname name, size_t num_dims, const ExprRef* extent);
    StmtRef add_BarrierEnvAlloc(Varname name, size_t num_dims, const ExprRef* extent);

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

struct BinOpTableEntry
{
    int second_char = -1;
    binop op = static_cast<binop>(0);
};

class BinOpTable
{
    BinOpTableEntry entries_by_char[256][2];
  public:
    BinOpTable();
    binop get(const char* p_str) const;
};

extern const BinOpTable binop_table;

}  // end namespace

// 0 or null returns signal an error.

CAMSPORK_EXPORT camspork::ProgramBuilder* camspork_new_ProgramBuilder();
CAMSPORK_EXPORT void camspork_delete_ProgramBuilder(camspork::ProgramBuilder* p_builder);
CAMSPORK_EXPORT int camspork_finish_ProgramBuilder(camspork::ProgramBuilder* p_builder);
CAMSPORK_EXPORT camspork::Varname camspork_add_variable(camspork::ProgramBuilder* p_builder, const char* p_name);

CAMSPORK_EXPORT camspork::ExprRef camspork_add_ReadValue(camspork::ProgramBuilder* p_builder,
    camspork::Varname name, uint32_t num_idx, const camspork::ExprRef* idx);
CAMSPORK_EXPORT camspork::ExprRef camspork_add_Const(camspork::ProgramBuilder* p_builder,
    camspork::value_t value);
CAMSPORK_EXPORT camspork::ExprRef camspork_add_USub(camspork::ProgramBuilder* p_builder,
    camspork::ExprRef arg);
CAMSPORK_EXPORT camspork::ExprRef camspork_add_BinOp(camspork::ProgramBuilder* p_builder,
    camspork::binop op, camspork::ExprRef lhs, camspork::ExprRef rhs);

CAMSPORK_EXPORT camspork::StmtRef camspork_add_SyncEnvAccess(camspork::ProgramBuilder* p_builder,
    camspork::Varname name, uint32_t num_idx, const camspork::OffsetExtentExpr* idx, uint32_t is_mutate, uint32_t is_ooo,
    camspork::qual_bits_t initial_qual_bit, camspork::qual_bits_t extended_qual_bits);
CAMSPORK_EXPORT camspork::StmtRef camspork_add_MutateValue(camspork::ProgramBuilder* p_builder,
    camspork::Varname name, uint32_t num_idx, const camspork::ExprRef* idx, camspork::binop op, camspork::ExprRef rhs);
CAMSPORK_EXPORT camspork::StmtRef camspork_add_Fence(camspork::ProgramBuilder* p_builder,
    uint32_t V1_transitive, camspork::qual_bits_t L1_qual_bits,
    camspork::qual_bits_t L2_full_qual_bits, camspork::qual_bits_t L2_temporal_qual_bits);
CAMSPORK_EXPORT camspork::StmtRef camspork_add_ValueEnvAlloc(camspork::ProgramBuilder* p_builder,
    camspork::Varname name, uint32_t num_dims, const camspork::ExprRef* extent);
CAMSPORK_EXPORT camspork::StmtRef camspork_add_SyncEnvAlloc(camspork::ProgramBuilder* p_builder,
    camspork::Varname name, uint32_t num_dims, const camspork::ExprRef* extent);
CAMSPORK_EXPORT camspork::StmtRef camspork_add_BarrierEnvAlloc(camspork::ProgramBuilder* p_builder,
    camspork::Varname name, uint32_t num_dims, const camspork::ExprRef* extent);

CAMSPORK_EXPORT camspork::StmtRef camspork_push_If(camspork::ProgramBuilder* p_builder,
    camspork::ExprRef cond);
CAMSPORK_EXPORT int camspork_begin_orelse(camspork::ProgramBuilder* p_builder);
CAMSPORK_EXPORT camspork::StmtRef camspork_push_SeqFor(camspork::ProgramBuilder* p_builder,
    camspork::Varname iter, camspork::ExprRef lo, camspork::ExprRef hi);
CAMSPORK_EXPORT camspork::StmtRef camspork_push_TasksFor(camspork::ProgramBuilder* p_builder,
    camspork::Varname iter, camspork::ExprRef lo, camspork::ExprRef hi);
CAMSPORK_EXPORT camspork::StmtRef camspork_push_ThreadsFor(camspork::ProgramBuilder* p_builder,
    camspork::Varname iter, camspork::ExprRef lo, camspork::ExprRef hi, uint32_t dim_idx, uint32_t offset, uint32_t box);
CAMSPORK_EXPORT camspork::StmtRef camspork_push_ParallelBlock(camspork::ProgramBuilder* p_builder,
    uint32_t dim, const uint32_t* domain);
CAMSPORK_EXPORT camspork::StmtRef camspork_push_DomainSplit(camspork::ProgramBuilder* p_builder,
    uint32_t dim_idx, uint32_t split_factor);
CAMSPORK_EXPORT int camspork_pop_body(camspork::ProgramBuilder* p_builder);


CAMSPORK_EXPORT camspork::binop camspork_binop_from_str(const char* p_str);
