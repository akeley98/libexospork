#pragma once

#include <cassert>
#include <new>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#define CAMSPORK_NODE_VLA_MEMBER(T) \
    static constexpr bool camspork_vla_member = true; \
    uint32_t camspork_vla_size = 0; \
    using camspork_vla_type = T; \
    static_assert(alignof(T) == 4); \
    /* Bytes needed for variable-length-array, rounded up to 4 */ \
    uint32_t camspork_vla_bytes() const { return (sizeof(T) * camspork_vla_size + 3) & ~uint32_t(3); } \
    uint32_t camspork_total_bytes() const { return sizeof(*this) + camspork_vla_bytes(); }

#define CAMSPORK_NODE_NO_VLA() \
    static constexpr bool camspork_vla_member = false; \
    uint32_t camspork_vla_bytes() const { return 0; } \
    uint32_t camspork_total_bytes() const { return sizeof(*this) + camspork_vla_bytes(); }

namespace camspork
{

template <typename Node>
const typename Node::camspork_vla_type& node_vla_get(const Node* p_node, uint32_t i)
{
    static_assert(alignof(Node) == 4);
    static_assert(Node::camspork_vla_member);
    const char* p_vla = static_cast<const char*>(p_node) + sizeof(Node);
    assert(i < p_node->camspork_vla_size);
    return reinterpret_cast<const typename Node::camspork_vla_type*>(p_vla)[i];
}

template <template<uint32_t> typename NodeType, uint32_t NumTypes>
struct NodeRef
{
    static_assert(NumTypes <= 32);

    uint32_t raw_data;

    // Bottom 5 bits holds the ID of the node type.
    uint32_t type_id() const
    {
        return raw_data & 31;
    }

    // Top 27 bits holds the address of the node in the file, multiplied by 4 bytes.
    uint32_t byte_offset() const
    {
        return (raw_data >> 3) & ~uint32_t(3);
    }

    void set_type_byte_offset(uint32_t type, size_t byte_offset)
    {
        assert(type < 32);
        raw_data = type | byte_offset << 3;
        assert(this->byte_offset() == byte_offset);
    }

    template <typename Callable>
    void dispatch(Callable&& callable, size_t buffer_size, char* buffer)
    {
        const size_t byte_offset = this->byte_offset();

        #define CAMSPORK_DISPATCH_CASE(N) \
          case N: \
            if constexpr (N < NumTypes) { \
                using Node = NodeType<N>; \
                assert(byte_offset + sizeof(Node) <= buffer_size); \
                auto p_node = reinterpret_cast<const Node*>(buffer + byte_offset); \
                assert(byte_offset + p_node->camspork_total_bytes() <= buffer_size); \
                callable(p_node); \
            }

        switch (type_id()) {
            CAMSPORK_DISPATCH_CASE(0)
            CAMSPORK_DISPATCH_CASE(1)
            CAMSPORK_DISPATCH_CASE(2)
            CAMSPORK_DISPATCH_CASE(3)
            CAMSPORK_DISPATCH_CASE(4)
            CAMSPORK_DISPATCH_CASE(5)
            CAMSPORK_DISPATCH_CASE(6)
            CAMSPORK_DISPATCH_CASE(7)
            CAMSPORK_DISPATCH_CASE(8)
            CAMSPORK_DISPATCH_CASE(9)
            CAMSPORK_DISPATCH_CASE(10)
            CAMSPORK_DISPATCH_CASE(11)
            CAMSPORK_DISPATCH_CASE(12)
            CAMSPORK_DISPATCH_CASE(13)
            CAMSPORK_DISPATCH_CASE(14)
            CAMSPORK_DISPATCH_CASE(15)
            CAMSPORK_DISPATCH_CASE(16)
            CAMSPORK_DISPATCH_CASE(17)
            CAMSPORK_DISPATCH_CASE(18)
            CAMSPORK_DISPATCH_CASE(19)
            CAMSPORK_DISPATCH_CASE(20)
            CAMSPORK_DISPATCH_CASE(21)
            CAMSPORK_DISPATCH_CASE(22)
            CAMSPORK_DISPATCH_CASE(23)
            CAMSPORK_DISPATCH_CASE(24)
            CAMSPORK_DISPATCH_CASE(25)
            CAMSPORK_DISPATCH_CASE(26)
            CAMSPORK_DISPATCH_CASE(27)
            CAMSPORK_DISPATCH_CASE(28)
            CAMSPORK_DISPATCH_CASE(29)
            CAMSPORK_DISPATCH_CASE(30)
            CAMSPORK_DISPATCH_CASE(31)
        }
    }
};

class NodeNursery
{
    // p_nursery_data is an allocation of size nursery_capacity bytes.
    // Of that, p_nursery_data[0:nursery_size] holds real data.
    uint32_t nursery_size = 0;
    uint32_t nursery_capacity = 0;
    char* p_nursery_data = 0;

    NodeNursery() = default;
    NodeNursery(NodeNursery&&) = delete;
    ~NodeNursery()
    {
        free(p_nursery_data);
    }

    template <template<uint32_t> typename NodeType, uint32_t TypeID>
    NodeRef<NodeType, TypeID>
    add_node(
        NodeType<TypeID> node,
        const std::vector<typename NodeType<TypeID>::camspork_vla_type>& vla)
    {
        node.camspork_vla_size = uint32_t(vla.size());
        const uint32_t offset = add_blob(sizeof(node), &node);
        add_blob(node.camspork_vla_bytes(), vla.data());
        NodeRef<NodeType, TypeID> node_ref;
        node_ref.set_type_byte_offset(TypeID, offset);
        return node_ref;
    }

    template <template<uint32_t> typename NodeType, uint32_t TypeID>
    NodeRef<NodeType, TypeID>
    add_node(NodeType<TypeID> node)
    {
        static_assert(!node.camspork_vla_member);
        const uint32_t offset = add_blob(sizeof(node), &node);
        NodeRef<NodeType, TypeID> node_ref;
        node_ref.set_type_byte_offset(TypeID, offset);
        return node_ref;
    }

    uint32_t add_blob(size_t bytes, const char* p_blob)
    {
        const uint32_t offset = nursery_size;
        assert(bytes % 4 == 0);
        reserve_bytes(nursery_size + bytes);
        memcpy(p_nursery_data + nursery_size, p_blob, bytes);
        size_t new_size = nursery_size + bytes;
        assert(new_size <= UINT32_MAX);
        nursery_size = uint32_t(new_size);
        return offset;
    }

  private:
    void reserve_bytes(size_t bytes)
    {
        if (bytes < nursery_capacity) {
            bytes = (bytes + 4095) & ~size_t(4095);
            char* p_new = static_cast<char*>(realloc(p_nursery_data, bytes));
            if (p_new == nullptr) {
                throw std::bad_alloc();
            }
            p_nursery_data = p_new;
            assert(bytes <= UINT32_MAX);
            nursery_capacity = uint32_t(bytes);
        }
    }
};


struct Varname
{
    uint32_t slot;
};

enum class binop
{
    assign,
    add,
    sub,
    mul,
    div,
    mod,
};

template <uint32_t TypeID>
struct expr
{
};

static constexpr uint32_t NumExprTypes = 4;

using ExprRef = NodeRef<expr, NumExprTypes>;

// ReadValue(Varname name, expr* idx)
using ReadValue = expr<0>;
template <>
struct expr<0>
{
    Varname name;
    CAMSPORK_NODE_VLA_MEMBER(ExprRef)

    uint32_t dim() const
    {
        return camspork_vla_size;
    }
};


// Const(int value)
using Const = expr<1>;
template <>
struct expr<1>
{
    int32_t value;
    CAMSPORK_NODE_NO_VLA()
};


// USub(expr arg)
using USub = expr<2>;
template<>
struct expr<2>
{
    ExprRef arg;
    CAMSPORK_NODE_NO_VLA()
};


// BinOp(binop op, expr lhs, expr rhs)
using BinOp = expr<3>;
template<>
struct expr<3>
{
    binop op;
    ExprRef lhs;
    ExprRef rhs;
    CAMSPORK_NODE_NO_VLA()
};

// Update this if you add more expr node types.
static_assert(NumExprTypes == 4);


struct OffsetExtentExpr
{
    ExprRef offset_e;
    ExprRef extent_e;
};



using qual_bits_t = uint32_t;

template <uint32_t TypeID>
struct stmt
{
};

static constexpr uint32_t NumStmtTypes = 19;

using StmtRef = NodeRef<stmt, NumStmtTypes>;

template <bool IsOOO, bool IsMutate>
struct SyncEnvAccessNode
{
    static constexpr bool is_ooo = IsOOO;
    static constexpr bool is_mutate = IsMutate;
    Varname name;
    uint32_t initial_qual_bit;
    uint32_t extended_qual_bits;
    CAMSPORK_NODE_VLA_MEMBER(OffsetExtentExpr)
};

// SyncEnvRead(Varname name, qual_tl initial_qual_bit, qual_tl* extended_qual_bits, expr* offset, expr* extent)
using SyncEnvRead = stmt<0>;
template <>
struct stmt<0> : SyncEnvAccessNode<false, false>
{
};

// SyncEnvReadOOO(Varname name, qual_tl initial_qual_bit, qual_tl* extended_qual_bits, expr* offset, expr* extent)
using SyncEnvReadOOO = stmt<1>;
template <>
struct stmt<1> : SyncEnvAccessNode<true, false>
{
};

// SyncEnvMutate(Varname name, qual_tl initial_qual_bit, qual_tl* extended_qual_bits, expr* offset, expr* extent)
using SyncEnvMutate = stmt<2>;
template <>
struct stmt<2> : SyncEnvAccessNode<false, true>
{
};

// SyncEnvMutateOOO(Varname name, qual_tl initial_qual_bit, qual_tl* extended_qual_bits, expr* offset, expr* extent)
using SyncEnvMutateOOO = stmt<3>;
template <>
struct stmt<3> : SyncEnvAccessNode<true, true>
{
};

// MutateValue(Varname name, expr e, binop op, expr* idx)
using MutateValue = stmt<4>;
template<>
struct stmt<4>
{
    Varname name;
    ExprRef e;
    binop op;
    CAMSPORK_NODE_VLA_MEMBER(ExprRef)
};

// Fence(qual_tl* L1_qual_bits, qual_tl* L2_full_qual_bits, qual_tl* L2_temporal_qual_bits)
using Fence = stmt<5>;
template<>
struct stmt<5>
{
    qual_bits_t L1_qual_bits;
    qual_bits_t L2_full_qual_bits;
    qual_bits_t L2_temporal_qual_bits;
    CAMSPORK_NODE_NO_VLA()
};

struct ArriveIdx
{
    Varname idx;
    uint32_t multicast_per_expr;
};

// Arrive(Varname name, bool V1_transitive, qual_tl* L1_qual_bits, Varname* idx, multicast_flag* multicast_flags)
using Arrive = stmt<6>;
template<>
struct stmt<6>
{
    Varname name;
    uint32_t V1_transitive;
    uint32_t L1_qual_bits;
    uint32_t num_barrier_exprs;
    CAMSPORK_NODE_VLA_MEMBER(ArriveIdx)
};

// Await(Varname name, qual_tl* L2_full_qual_bits, qual_tl* L2_temporal_qual_bits, Varname* idx)
using Await = stmt<7>;
template<>
struct stmt<7>
{
    Varname name;
    uint32_t L2_full_qual_bits;
    uint32_t L2_temporal_qual_bits;
    CAMSPORK_NODE_VLA_MEMBER(Varname)
};

// ValueEnvAlloc(Varname name, expr* extent)
using ValueEnvAlloc = stmt<8>;
template<>
struct stmt<8>
{
    Varname name;
    CAMSPORK_NODE_VLA_MEMBER(ExprRef)
};

// SyncEnvAlloc(Varname name, expr* extent)
using SyncEnvAlloc = stmt<9>;
template<>
struct stmt<9>
{
    Varname name;
    CAMSPORK_NODE_VLA_MEMBER(ExprRef)
};

// SyncEnvFreeShard(Varname name, Varname* distributed_iters)
using SyncEnvFreeShard = stmt<10>;
template<>
struct stmt<10>
{
    Varname name;
    CAMSPORK_NODE_VLA_MEMBER(Varname)
};

// BarrierEnvAlloc(Varname name, expr* extent)
using BarrierEnvAlloc = stmt<11>;
template<>
struct stmt<11>
{
    Varname name;
    CAMSPORK_NODE_VLA_MEMBER(ExprRef)
};

// BarrierEnvFree(Varname name)
using BarrierEnvFree = stmt<12>;
template<>
struct stmt<12>
{
    Varname name;
    CAMSPORK_NODE_NO_VLA()
};

// If(expr cond, stmt* body, stmt* orelse)
using If = stmt<13>;
template<>
struct stmt<13>
{
    // body = statements [0 : camspork_vla_size - num_orelse_stmts]
    // orelse = statements [camspork_vla_size - num_orelse_stmts : camspork_vla_size]
    ExprRef cond;
    uint32_t num_orelse_stmts;
    CAMSPORK_NODE_VLA_MEMBER(StmtRef)
};

// SeqFor(Varname iter, expr lo, expr hi, stmt* body)
using SeqFor = stmt<14>;
template<>
struct stmt<14>
{
    Varname iter;
    ExprRef lo;
    ExprRef hi;
    CAMSPORK_NODE_VLA_MEMBER(StmtRef)
};

// TasksFor(Varname iter, expr lo, expr hi, stmt* body)
using TasksFor = stmt<15>;
template<>
struct stmt<15>
{
    Varname iter;
    ExprRef lo;
    ExprRef hi;
    CAMSPORK_NODE_VLA_MEMBER(StmtRef)
};

// ThreadsFor(Varname iter, expr lo, expr hi, int dim_idx, int offset, int box_size, stmt* body)
using ThreadsFor = stmt<16>;
template<>
struct stmt<16>
{
    Varname iter;
    ExprRef lo;
    ExprRef hi;
    uint32_t dim_idx;
    uint32_t offset;
    uint32_t box_size;
    CAMSPORK_NODE_VLA_MEMBER(StmtRef)
};

// DomainDefine(stmt child_stmt, int* extent)
using DomainDefine = stmt<17>;
template<>
struct stmt<17>
{
    StmtRef child_stmt;
    CAMSPORK_NODE_VLA_MEMBER(uint32_t)
};

// DomainSplit(stmt child_stmt, int dim_idx, int split_factor)
using DomainSplit = stmt<18>;
template<>
struct stmt<18>
{
    StmtRef child_stmt;
    uint32_t dim_idx;
    uint32_t split_factor;
};

// Update this if you add more stmt node types.
static_assert(NumStmtTypes == 19);

}
