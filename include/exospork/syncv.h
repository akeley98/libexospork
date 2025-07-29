// Synchronization validation
#pragma once
#ifndef EXOSPORK_SYNCV_H_
#define EXOSPORK_SYNCV_H_

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EXOSPORK_SYNC_ACCESS_BIT 0x80000000

struct exospork_syncv_init_t
{
    const char* filename;
    uint64_t memory_budget;
    uint64_t debug_assignment_id;
    uint64_t debug_operation_id;
};

// Opaque value associated with each unique memory location of the program
// undergoing synchronization validation.
// i.e. there should be a parallel array of exospork_syncv_value_t for each
// array of values in the program.
struct exospork_syncv_value_t
{
    uint32_t node_id;
#ifdef __cplusplus
    operator bool() const
    {
        return node_id != 0;
    }
#endif
};

// Opaque value for each barrier.
struct exospork_syncv_barrier_t
{
    uint32_t data;
};

void exospork_syncv_init(size_t sizeof_struct,
                         const exospork_syncv_init_t*);

void exospork_syncv_deinit();

void exospork_syncv_r(int lineno, size_t N, exospork_syncv_value_t* values,
                      uint32_t tid_lo, uint32_t tid_hi, uint32_t bitfield);

void exospork_syncv_rw(int lineno, size_t N, exospork_syncv_value_t* values,
                       uint32_t tid_lo, uint32_t tid_hi, uint32_t bitfield);

void exospork_syncv_clear(int lineno, size_t N, exospork_syncv_value_t* values);

void exospork_syncv_alloc_barrier(int lineno, exospork_syncv_barrier_t* bar);

void exospork_syncv_free_barrier(int lineno, exospork_syncv_barrier_t* bar);

void exospork_syncv_fence(
        int lineno, uint32_t tid_lo, uint32_t tid_hi,
        uint32_t v1_sigbits, uint32_t v2_sigbits, int transitive);

void exospork_syncv_arrive(
        int lineno, exospork_syncv_barrier_t* bar,
        uint32_t v1_tid_lo, uint32_t v1_tid_hi, uint32_t v1_sigbits,
        int transitive);

void exospork_syncv_await(
        int lineno, exospork_syncv_barrier_t* bar,
        uint32_t v2_tid_lo, uint32_t v2_tid_hi, uint32_t v2_sigbits);



// *** Debugging / Testing Backdoors ***

void exospork_syncv_debug_register_values(
        uint32_t N, exospork_syncv_value_t* values);

void exospork_syncv_debug_unregister_values(
        uint32_t N, exospork_syncv_value_t* values);



#ifdef __cplusplus
}  // extern "C"
#endif

#endif
