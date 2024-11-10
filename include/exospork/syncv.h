#pragma once
#ifndef EXOSPORK_SYNCV_H_
#define EXOSPORK_SYNCV_H_

#include <stddef.h>
#include <stdint.h>

struct exospork_syncv_init_t
{
    const char* filename;
    uint64_t log_assignment_index;
};

struct exospork_syncv_value_t
{
    uint64_t data_0, data_1;
};

struct exospork_syncv_bar_t
{
    uint64_t data;
};

void exospork_syncv_init(size_t sizeof_struct,
                         const exospork_syncv_init_t*);

void exospork_syncv_deinit();

void exospork_syncv_r(int lineno, size_t N, exospork_syncv_value_t* values,
                      uint32_t tid_lo, uint32_t tid_hi, uint32_t bitfield);

void exospork_syncv_rw(int lineno, size_t N, exospork_syncv_value_t* values,
                       uint32_t tid_lo, uint32_t tid_hi, uint32_t bitfield);

void exospork_syncv_alloc_bar(int lineno, exospork_syncv_bar_t* bar,
                              int max_pipeline);

void exospork_syncv_free_bar(int lineno, exospork_syncv_bar_t* bar);

void exospork_syncv_fence(
        int lineno,
        uint32_t v1_tid_lo, uint32_t v1_tid_hi, uint32_t v1_sigbits,
        uint32_t v2_tid_lo, uint32_t v2_tid_hi, uint32_t v2_sigbits);

void exospork_syncv_arrive(
        int lineno, exospork_syncv_bar_t* bar,
        uint32_t v1_tid_lo, uint32_t v1_tid_hi, uint32_t v1_sigbits);

void exospork_syncv_await(
        int lineno, exospork_syncv_bar_t* bar,
        uint32_t v2_tid_lo, uint32_t v2_tid_hi, uint32_t v2_sigbits);

#endif
