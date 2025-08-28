import os
from ctypes import *
from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple

lib = cdll.LoadLibrary(os.path.join(os.path.split(__name__)[0], "bin/libexospork.so"))


class BuilderExpr:
    __slots__ = []

    @staticmethod
    def typecheck(item):
        if isinstance(item, int):
            return BuilderConst(item)
        elif isinstance(item, ExprRef):
            return item
        else:
            assert isinstance(item, BuilderExpr), "expected int, ExprRef, or BuilderExpr"
            return item

    def __add__(self, other):
        return BuilderBinOp(binop_Add, self, self.typecheck(other))

    def __radd__(self, other):
        return BuilderBinOp(binop_Add, self.typecheck(other), self)

    def __sub__(self, other):
        return BuilderBinOp(binop_Sub, self, self.typecheck(other))

    def __rsub__(self, other):
        return BuilderBinOp(binop_Sub, self.typecheck(other), self)

    def __mul__(self, other):
        return BuilderBinOp(binop_Mul, self, self.typecheck(other))

    def __rmul__(self, other):
        return BuilderBinOp(binop_Mul, self, self.typecheck(other))

    def __truediv__(self, other):
        return BuilderBinOp(binop_Div, self, self.typecheck(other))

    def __floordiv__(self, other):
        return BuilderBinOp(binop_Div, self, self.typecheck(other))

    def __mod__(self, other):
        return BuilderBinOp(binop_Mod, self, self.typecheck(other))

    def __lt__(self, other):
        return BuilderBinOp(binop_Less, self, self.typecheck(other))

    def __le__(self, other):
        return BuilderBinOp(binop_Leq, self, self.typecheck(other))

    def __gt__(self, other):
        return BuilderBinOp(binop_Greater, self, self.typecheck(other))

    def __ge__(self, other):
        return BuilderBinOp(binop_Geq, self, self.typecheck(other))

    def __eq__(self, other):
        return BuilderBinOp(binop_Eq, self, self.typecheck(other))

    def __neq__(self, other):
        return BuilderBinOp(binop_Neq, self, self.typecheck(other))

    def __neg__(self):
        return BuilderUSub(self)


def check_return(code):
    if not code:
        raise ValueError(str(_thread_local_message_c_str(), "utf-8"))
    return code

class VoidPtr(c_void_p):
    pass


# These can only be passed to the builder that produced them!
class ExprRef(Structure, BuilderExpr):
    _fields_ = [("raw_data", c_uint32)]
    def __bool__(self):
        return self.raw_data != 0  # 0 used to signal error (use check_return)

    def build_expr(self, builder):
        return self

class StmtRef(Structure):
    _fields_ = [("raw_data", c_uint32)]
    def __bool__(self):
        return self.raw_data != 0  # 0 used to signal error (use check_return)

class Varname(Structure, BuilderExpr):
    _fields_ = [("slot_1_index", c_uint32)]
    def __bool__(self):
        return self.slot_1_index != 0  # 0 used to signal error (use check_return)

    def as_index_expr(self):
        return BuilderIndexExpr(self, ())

    def build_expr(self, builder):
        return BuilderIndexExpr(self, ()).build_expr(builder)

    def __getitem__(self, idxs):
        return BuilderIndexExpr(self, ())[idxs]

class OffsetExtentExpr(Structure):
    _fields_ = [("offset", ExprRef), ("extent", ExprRef)]

# binop enum values are always the same for a given operator (_binop_from_str)
class binop(Structure):
    _fields_ = [("enum_value", c_uint32)]
    def __bool__(self):
        return self.enum_value != 0  # 0 used to signal error (use check_return)


ptr_uint32 = POINTER(c_uint32)
ptr_ExprRef = POINTER(ExprRef)
ptr_OffsetExtentExpr = POINTER(OffsetExtentExpr)

_thread_local_message_c_str = lib.camspork_thread_local_message_c_str
_thread_local_message_c_str.restype = c_char_p
_thread_local_message_c_str.argtypes = ()

_thread_local_print_program = lib.camspork_thread_local_print_program
_thread_local_print_program.restype = c_int
_thread_local_print_program.argtypes = (c_size_t, c_void_p)

_new_ProgramBuilder = lib.camspork_new_ProgramBuilder
_new_ProgramBuilder.restype = VoidPtr
_new_ProgramBuilder.argtypes = ()

_delete_ProgramBuilder = lib.camspork_delete_ProgramBuilder
_delete_ProgramBuilder.restype = None
_delete_ProgramBuilder.argtypes = (c_void_p, )

_finish_ProgramBuilder = lib.camspork_finish_ProgramBuilder
_finish_ProgramBuilder.restype = c_int
_finish_ProgramBuilder.argtypes = (c_void_p, )

_ProgramBuilder_size = lib.camspork_ProgramBuilder_size
_ProgramBuilder_size.restype = c_size_t
_ProgramBuilder_size.argtypes = (c_void_p, )

_ProgramBuilder_data = lib.camspork_ProgramBuilder_data
_ProgramBuilder_data.restype = VoidPtr
_ProgramBuilder_data.argtypes = (c_void_p, )

_add_variable = lib.camspork_add_variable
_add_variable.restype = Varname
_add_variable.argtypes = (c_void_p, c_char_p)

_add_ReadValue = lib.camspork_add_ReadValue
_add_ReadValue.restype = ExprRef
_add_ReadValue.argtypes = (c_void_p, Varname, c_uint32, ptr_ExprRef)

_add_Const = lib.camspork_add_Const
_add_Const.restype = ExprRef
_add_Const.argtypes = (c_void_p, c_int32)

_add_USub = lib.camspork_add_USub
_add_USub.restype = ExprRef
_add_USub.argtypes = (c_void_p, ExprRef)

_add_BinOp = lib.camspork_add_BinOp
_add_BinOp.restype = ExprRef
_add_BinOp.argtypes = (c_void_p, binop, ExprRef, ExprRef)

_add_SyncEnvAccess = lib.camspork_add_SyncEnvAccess
_add_SyncEnvAccess.restype = StmtRef
_add_SyncEnvAccess.argtypes = (c_void_p, Varname, c_uint32, ptr_OffsetExtentExpr, c_uint32, c_uint32, c_uint32, c_uint32)

_add_MutateValue = lib.camspork_add_MutateValue
_add_MutateValue.restype = StmtRef
_add_MutateValue.argtypes = (c_void_p, Varname, c_uint32, ptr_ExprRef, binop, ExprRef)

_add_Fence = lib.camspork_add_Fence
_add_Fence.restype = StmtRef
_add_Fence.argtypes = (c_void_p, c_uint32, c_uint32, c_uint32, c_uint32)

_add_ValueEnvAlloc = lib.camspork_add_ValueEnvAlloc
_add_ValueEnvAlloc.restype = StmtRef
_add_ValueEnvAlloc.argtypes = (c_void_p, Varname, c_uint32, ptr_ExprRef)

_add_SyncEnvAlloc = lib.camspork_add_SyncEnvAlloc
_add_SyncEnvAlloc.restype = StmtRef
_add_SyncEnvAlloc.argtypes = (c_void_p, Varname, c_uint32, ptr_ExprRef)

_add_BarrierEnvAlloc = lib.camspork_add_BarrierEnvAlloc
_add_BarrierEnvAlloc.restype = StmtRef
_add_BarrierEnvAlloc.argtypes = (c_void_p, Varname, c_uint32, ptr_ExprRef)

_push_If = lib.camspork_push_If
_push_If.restype = StmtRef
_push_If.argtypes = (c_void_p, ExprRef)

_begin_orelse = lib.camspork_begin_orelse
_begin_orelse.restype = c_int
_begin_orelse.argtypes = (c_void_p,)

_push_SeqFor = lib.camspork_push_SeqFor
_push_SeqFor.restype = StmtRef
_push_SeqFor.argtypes = (c_void_p, Varname, ExprRef, ExprRef)

_push_TasksFor = lib.camspork_push_TasksFor
_push_TasksFor.restype = StmtRef
_push_TasksFor.argtypes = (c_void_p, Varname, ExprRef, ExprRef)

_push_ThreadsFor = lib.camspork_push_ThreadsFor
_push_ThreadsFor.restype = StmtRef
_push_ThreadsFor.argtypes = (c_void_p, Varname, ExprRef, ExprRef, c_uint32, c_uint32, c_uint32)

_push_ParallelBlock = lib.camspork_push_ParallelBlock
_push_ParallelBlock.restype = StmtRef
_push_ParallelBlock.argtypes = (c_void_p, c_uint32, ptr_uint32)

_push_DomainSplit = lib.camspork_push_DomainSplit
_push_DomainSplit.restype = StmtRef
_push_DomainSplit.argtypes = (c_void_p, c_uint32, c_uint32)

_pop_body = lib.camspork_pop_body
_pop_body.restype = c_int
_pop_body.argtypes = (c_void_p,)

_binop_from_str = lib.camspork_binop_from_str
_binop_from_str.restype = binop
_binop_from_str.argtypes = (c_char_p,)

_binop_to_str = lib.camspork_binop_to_str
_binop_to_str.restype = c_char_p
_binop_to_str.arg_types = (binop,)


binop_Assign = check_return(_binop_from_str(b"="))
binop_Add = check_return(_binop_from_str(b"+"))
binop_Sub = check_return(_binop_from_str(b"-"))
binop_Mul = check_return(_binop_from_str(b"*"))
binop_Div = check_return(_binop_from_str(b"/"))
binop_Mod = check_return(_binop_from_str(b"%"))
binop_Less = check_return(_binop_from_str(b"<"))
binop_Leq = check_return(_binop_from_str(b"<="))
binop_Greater = check_return(_binop_from_str(b">"))
binop_Geq = check_return(_binop_from_str(b">="))
binop_Eq = check_return(_binop_from_str(b"=="))
binop_Neq = check_return(_binop_from_str(b"!="))


def to_binop(op):
    if isinstance(op, binop):
        return op
    elif isinstance(op, str):
        return check_return(_binop_from_str(bytes(op, "utf8")))
    else:
        assert isinstance(op, bytes)
        return check_return(_binop_from_str(op))


class BodyCtx:
    __slots__ = ["_builder", "_on_enter", "stmt", "body", "orelse"]
    _builder: VoidPtr
    _on_enter: Callable[[VoidPtr], None]
    stmt: StmtRef
    body: StmtRef
    orelse: StmtRef

    def __init__(self, builder, on_enter):
        self._builder = builder
        self._on_enter = on_enter

    def __enter__(self, *a):
        stmt = check_return(self._on_enter(self._builder))
        assert isinstance(stmt, StmtRef)
        self.stmt = stmt

    def __exit__(self, *a):
        _pop_body(self._builder)


@dataclass(slots=True)
class BuilderIndexExpr(BuilderExpr):
    _varname: Varname
    _idx: Tuple[BuilderExpr | ExprRef]

    def c_dim_idxs(self, builder):
        dim = len(self._idx)
        if dim == 0:
            return 0, ptr_ExprRef()
        else:
            e = (ExprRef * dim)()
            for i, tmp in enumerate(self._idx):
                e[i] = tmp.build_expr(builder)
            return dim, e

    def build_expr(self, builder) -> ExprRef:
        # When interpreted as an expression, generate ReadValue"""
        dim, e = self.c_dim_idxs(builder)
        return check_return(_add_ReadValue(builder, self._varname, dim, e))

    def as_index_expr(self):
        return self

    def __getitem__(self, a):
        if isinstance(a, tuple):
            a = tuple(self.typecheck(v) for v in a)
        else:
            a = (self.typecheck(a),)
        return BuilderIndexExpr(self._varname, self._idx + a)

@dataclass(slots=True)
class BuilderConst(BuilderExpr):
    _value: int

    def build_expr(self, builder) -> ExprRef:
        return check_return(_add_Const(builder, self._value))


@dataclass(slots=True)
class BuilderUSub(BuilderExpr):
    _arg: BuilderExpr | ExprRef

    def build_expr(self, builder) -> ExprRef:
        return check_return(_add_USub(builder, self._arg.build_expr(builder)))


@dataclass(slots=True)
class BuilderBinOp(BuilderExpr):
    _binop: binop
    _lhs: BuilderExpr | ExprRef
    _rhs: BuilderExpr | ExprRef

    def build_expr(self, builder) -> ExprRef:
        return check_return(_add_BinOp(builder, self._binop, self._lhs.build_expr(builder), self._rhs.build_expr(builder)))


class ProgramBuilder:
    __slots__ = ["_builder", "_varname_dict"]

    _builder: VoidPtr
    _varname_dict: Dict[object, Varname]

    def __init__(self):
        self._builder = check_return(_new_ProgramBuilder())
        self._varname_dict = {}

    def __del__(self):
        _delete_ProgramBuilder(self._builder)
        self._builder = 0

    def finish(self):
        check_return(_finish_ProgramBuilder(self._builder))

    def add_variable(self, name, to_ascii=lambda name: bytes(str(name), "utf8")) -> Varname:
        assert name not in self._varname_dict, f"Duplicate variable name {name!r}"
        varname = check_return(_add_variable(self._builder, to_ascii(name)))
        self._varname_dict[name] = varname
        return varname

    def get_varname(self, var):
        if isinstance(var, Varname):
            return var
        else:
            return self._varname_dict[var]

    def __getitem__(self, var):
        return self.get_varname(var)

    def build_expr(self, e) -> ExprRef:
        return BuilderExpr.typecheck(e).build_expr(self._builder)

    def MutateValue(self, dst: BuilderIndexExpr, op, rhs) -> StmtRef:
        dim, idxs = dst.c_dim_idxs(self._builder)
        check_return(_add_MutateValue(self._builder, dst._varname, dim, idxs, to_binop(op), self.build_expr(rhs)))

    def Fence(self, V1_transitive: bool, L1_qual_bits: int, L2_full_qual_bits: int, L2_temporal_qual_bits: int):
        check_return(_add_Fence(self._builder, V1_transitive, L1_qual_bits, L2_full_qual_bits, L2_temporal_qual_bits))

    def ValueEnvAlloc(self, e: Varname | BuilderIndexExpr):
        self._add_alloc(_add_ValueEnvAlloc, e)

    def SyncEnvAlloc(self, e: Varname | BuilderIndexExpr):
        self._add_alloc(_add_SyncEnvAlloc, e)

    def BarrierEnvAlloc(self, e: Varname | BuilderIndexExpr):
        self._add_alloc(_add_SyncEnvAlloc, e)

    def _add_alloc(self, c_adder, e):
        e = e.as_index_expr()
        dim, idxs = e.c_dim_idxs(self._builder)
        c_adder(self._builder, e._varname, dim, idxs)

    def If(self, cond) -> BodyCtx:
        cond = self.build_expr(cond)
        return BodyCtx(self._builder, lambda builder: _push_If(builder, cond))

    def begin_orelse(self):
        check_return(_begin_orelse(self._builder))

    def SeqFor(self, var, lo, hi) -> BodyCtx:
        var = self.get_varname(var)
        lo = self.build_expr(lo)
        hi = self.build_expr(hi)
        return BodyCtx(self._builder, lambda builder: _push_SeqFor(builder, var, lo, hi))

    def TasksFor(self, var, lo, hi) -> BodyCtx:
        var = self.get_varname(var)
        lo = self.build_expr(lo)
        hi = self.build_expr(hi)
        return BodyCtx(self._builder, lambda builder: _push_TasksFor(builder, var, lo, hi))

    def ThreadsFor(self, var, lo, hi, dim_idx: int, offset: int, box: int) -> BodyCtx:
        var = self.get_varname(var)
        lo = self.build_expr(lo)
        hi = self.build_expr(hi)
        assert isinstance(dim_idx, int)
        assert isinstance(offset, int)
        assert isinstance(box, int)
        return BodyCtx(self._builder, lambda builder: _push_ThreadsFor(builder, var, lo, hi, dim_idx, offset, box))

    def ParallelBlock(self, *coords):
        dim = len(coords)
        array = (c_uint32 * dim)(coords)
        return BodyCtx(self._builder, lambda builder: _push_ParallelBlock(builder, dim, array))

    def DomainSplit(self, dim_idx: int, split_factor: int):
        assert isinstance(dim_idx, int)
        assert isinstance(split_factor, int)
        return BodyCtx(self._builder, lambda builder: _push_DomainSplit(builder, dim_idx, split_factor))


if __name__ == "__main__":
    b = ProgramBuilder()
    fib_size = 20
    _fib = b.add_variable("fib")
    _iter = b.add_variable("iter")
    b.ValueEnvAlloc(_fib[fib_size])
    b.MutateValue(_fib[0], "=", 0)
    b.MutateValue(_fib[1], "=", 1)
    with b.SeqFor(_iter, 2, fib_size):
        b.MutateValue(_fib[_iter], "=", _fib[_iter-1] + _fib[_iter-2])

    _dst = b.add_variable("dst")
    b.ValueEnvAlloc(_dst[fib_size])
    with b.SeqFor(_iter, 0, fib_size):
      with b.If(_fib % 5):
        b.MutateValue(_fib[_iter], "=", -_fib)
        b.MutateValue(_fib[_iter], "*", 10000)
        b.begin_orelse()
        b.MutateValue(_fib[_iter], "/", 5)

    b.finish()
    check_return(_thread_local_print_program(_ProgramBuilder_size(b._builder), _ProgramBuilder_data(b._builder)))
    print(str(_thread_local_message_c_str(), "utf8"))
