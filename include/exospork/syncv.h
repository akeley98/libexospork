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


#ifdef __cplusplus
}  // extern "C"
#endif

#endif
