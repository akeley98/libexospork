#include "builder.hpp"

namespace camspork
{

void BodyBuilder::begin_orelse()
{
    CAMSPORK_REQUIRE(!is_orelse, "Already called begin_orelse.");
    CAMSPORK_REQUIRE(body_of.type_id() == If_ID, "Must be defining If-statement body to call begin_orelse");
    p_program_builder->end_body_builder(*this);
    stmts.clear();
    is_orelse = true;
}

StmtRef BodyBuilder::body_to_nursery() const
{
    return p_program_builder->body_to_nursery(stmts);
}

ProgramBuilder::ProgramBuilder()
{
    // Start out the generated binary with the ProgramHeader.
    ProgramHeader header_template;
    for (const uint32_t& magic : ProgramHeader::expected_magic_numbers) {
        header_template.magic_numbers[&magic - &ProgramHeader::expected_magic_numbers[0]] = magic;
    }
    nursery.add_blob(sizeof header_template, &header_template);

    // Top-level BodyBuilder is special; we never pop it, and it is for building the top-level stmt body.
    body_stack.push_back({this, {}, {}});
}

void ProgramBuilder::finish()
{
    if (is_finished) {
        return;
    }

    // Package up top-level statements ... requires all If/For were finished.
    CAMSPORK_REQUIRE_CMP(body_stack.size(), ==, 1, "missing pop_body() before finish()");
    StmtRef top_level_stmt = body_to_nursery(body_stack[0].stmts);

    // Put together variable table.
    std::vector<VarConfigRef> var_configs;
    for (const std::string& name : variable_slot_names) {
        var_configs.push_back(nursery.add_node<VarConfigRef>(VarConfig{}, name.size(), name.data()));
    }
    auto var_config_table = nursery.add_node<VarConfigTableRef>(
        VarConfigTable{}, var_configs.size(), var_configs.data());

    // Fill header; don't access it until we are sure we never call add_node again!
    auto& header = reinterpret_cast<ProgramHeader&>(nursery.data()[0]);
    header.top_level_stmt = top_level_stmt;
    header.var_config_table = var_config_table;

    is_finished = true;
}

ExprRef ProgramBuilder::add_ReadValue(Varname name, size_t num_idx, const ExprRef* idx)
{
    check_not_finished();
    return nursery.add_node<ExprRef>(ReadValue{name}, num_idx, idx);
}

ExprRef ProgramBuilder::add_Const(value_t value)
{
    check_not_finished();
    return nursery.add_node<ExprRef>(Const{value});
}

ExprRef ProgramBuilder::add_USub(ExprRef arg)
{
    check_not_finished();
    return nursery.add_node<ExprRef>(USub{arg});
}

ExprRef ProgramBuilder::add_BinOp(binop op, ExprRef lhs, ExprRef rhs)
{
    check_not_finished();
    return nursery.add_node<ExprRef>(BinOp{op, lhs, rhs});
}

StmtRef ProgramBuilder::add_MutateValue(Varname name, size_t num_idx, const ExprRef* idx, binop op, ExprRef rhs)
{
    return append_impl(MutateValue{name, op, rhs}, num_idx, idx);
}

StmtRef ProgramBuilder::add_ValueEnvAlloc(Varname name, size_t num_dims, const ExprRef* extent)
{
    return append_impl(ValueEnvAlloc{name}, num_dims, extent);
}

StmtRef ProgramBuilder::push_If(ExprRef cond)
{
    return push_impl(If{cond, {}, {}});
}

void ProgramBuilder::begin_orelse()
{
    body_stack.back().begin_orelse();
}

StmtRef ProgramBuilder::push_SeqFor(Varname iter, ExprRef lo, ExprRef hi)
{
    SeqFor node;
    node.iter = iter;
    node.lo = lo;
    node.hi = hi;
    return push_impl(node);
}

StmtRef ProgramBuilder::push_TasksFor(Varname iter, ExprRef lo, ExprRef hi)
{
    TasksFor node;
    node.iter = iter;
    node.lo = lo;
    node.hi = hi;
    return push_impl(node);
}

void ProgramBuilder::check_not_finished() const
{
    CAMSPORK_REQUIRE(!is_finished, "Cannot modify program after finish()");
}

template <typename...Args>
StmtRef ProgramBuilder::append_impl(Args... a)
{
    check_not_finished();
    StmtRef ref = nursery.add_node<StmtRef>(a...);
    body_stack.back().stmts.push_back(ref);
    return ref;
}

template <typename Stmt>
StmtRef ProgramBuilder::push_impl(Stmt s)
{
    check_not_finished();
    StmtRef body_of = nursery.add_node<StmtRef>(s);
    body_stack.back().stmts.push_back(body_of);
    body_stack.push_back({this, body_of, {}});
    return body_of;
}

void ProgramBuilder::pop_body()
{
    CAMSPORK_REQUIRE_CMP(body_stack.size(), >=, 2, "Cannot pop last statement body; use finish()");
    end_body_builder(body_stack.back());
    body_stack.pop_back();
}

void ProgramBuilder::end_body_builder(const BodyBuilder& builder)
{
    CAMSPORK_REQUIRE(builder.body_of, "Internal error, null builder.body_of");
    builder.body_of.dispatch(builder, nursery.size(), nursery.data());
};

StmtRef ProgramBuilder::body_to_nursery(const std::vector<StmtRef>& stmts)
{
    check_not_finished();
    if (stmts.size() == 1) {
        return stmts[0];
    }
    else if (stmts.size() == 0) {
        return StmtRef{};
    }
    else {
        return nursery.add_node<StmtRef>(StmtBody{}, stmts);
    }
}

}
