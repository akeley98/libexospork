from ctypes import *

from typing import Dict

lib = cdll.LoadLibrary("bin/libexospork.so")

def check_return(code):
    if not code:
        raise ValueError(str(_thread_local_message_c_str(), "utf-8"))
    return code

class VoidPtr(c_void_p):
    pass

ptr_uint32 = POINTER(c_uint32)

_thread_local_message_c_str = lib.camspork_thread_local_message_c_str
_thread_local_message_c_str.restype = c_char_p
_thread_local_message_c_str.argtypes = ()

_new_ProgramBuilder = lib.camspork_new_ProgramBuilder
_new_ProgramBuilder.restype = VoidPtr
_new_ProgramBuilder.argtypes = ()

_delete_ProgramBuilder = lib.camspork_delete_ProgramBuilder
_delete_ProgramBuilder.restype = None
_delete_ProgramBuilder.argtypes = (c_void_p, )

_finish_ProgramBuilder = lib.camspork_finish_ProgramBuilder
_finish_ProgramBuilder.restype = c_int
_finish_ProgramBuilder.argtypes = (c_void_p, )

_add_variable = lib.camspork_add_variable
_add_variable.restype = c_uint32
_add_variable.argtypes = (c_void_p, c_char_p)

_add_ReadValue = lib.camspork_add_ReadValue
_add_ReadValue.restype = c_uint32
_add_ReadValue.argtypes = (c_void_p, c_uint32, c_uint32, ptr_uint32)

_add_Const = lib.camspork_add_Const
_add_Const.restype = c_uint32
_add_Const.argtypes = (c_void_p, c_int32)

_add_USub = lib.camspork_add_USub
_add_USub.restype = c_uint32
_add_USub.argtypes = (c_void_p, c_uint32)

_add_BinOp = lib.camspork_add_BinOp
_add_BinOp.restype = c_uint32
_add_BinOp.argtypes = (c_void_p, c_uint32, c_uint32, c_uint32)

_add_SyncEnvAccess = lib.camspork_add_SyncEnvAccess
_add_SyncEnvAccess.restype = c_uint32
_add_SyncEnvAccess.argtypes = (c_void_p, c_uint32, c_uint32, ptr_uint32, c_uint32, c_uint32, c_uint32, c_uint32)

_add_MutateValue = lib.camspork_add_MutateValue
_add_MutateValue.restype = c_uint32
_add_MutateValue.argtypes = (c_void_p, c_uint32, c_uint32, ptr_uint32, c_uint32, c_uint32)

_add_Fence = lib.camspork_add_Fence
_add_Fence.restype = c_uint32
_add_Fence.argtypes = (c_void_p, c_uint32, c_uint32, c_uint32, c_uint32)

_add_ValueEnvAlloc = lib.camspork_add_ValueEnvAlloc
_add_ValueEnvAlloc.restype = c_uint32
_add_ValueEnvAlloc.argtypes = (c_void_p, c_uint32, c_uint32, ptr_uint32)

_add_SyncEnvAlloc = lib.camspork_add_SyncEnvAlloc
_add_SyncEnvAlloc.restype = c_uint32
_add_SyncEnvAlloc.argtypes = (c_void_p, c_uint32, c_uint32, ptr_uint32)

_add_BarrierEnvAlloc = lib.camspork_add_BarrierEnvAlloc
_add_BarrierEnvAlloc.restype = c_uint32
_add_BarrierEnvAlloc.argtypes = (c_void_p, c_uint32, c_uint32, ptr_uint32)

_push_If = lib.camspork_push_If
_push_If.restype = c_uint32
_push_If.argtypes = (c_void_p, c_uint32)

_begin_orelse = lib.camspork_begin_orelse
_begin_orelse.restype = c_uint32
_begin_orelse.argtypes = (c_void_p,)

_push_SeqFor = lib.camspork_push_SeqFor
_push_SeqFor.restype = c_uint32
_push_SeqFor.argtypes = (c_void_p, c_uint32, c_uint32, c_uint32)

_push_TasksFor = lib.camspork_push_TasksFor
_push_TasksFor.restype = c_uint32
_push_TasksFor.argtypes = (c_void_p, c_uint32, c_uint32, c_uint32)

_push_ThreadsFor = lib.camspork_push_ThreadsFor
_push_ThreadsFor.restype = c_uint32
_push_ThreadsFor.argtypes = (c_void_p, c_uint32, c_uint32, c_uint32, c_uint32, c_uint32, c_uint32)

_push_ParallelBlock = lib.camspork_push_ParallelBlock
_push_ParallelBlock.restype = c_uint32
_push_ParallelBlock.argtypes = (c_void_p, c_uint32, ptr_uint32)

_push_DomainSplit = lib.camspork_push_DomainSplit
_push_DomainSplit.restype = c_uint32
_push_DomainSplit.argtypes = (c_void_p, c_uint32, c_uint32)

_pop_body = lib.camspork_pop_body
_pop_body.restype = c_uint32
_pop_body.argtypes = (c_void_p,)

_binop_from_str = lib.camspork_binop_from_str
_binop_from_str.restype = c_uint32
_binop_from_str.argtypes = (c_char_p,)
