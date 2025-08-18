#pragma once

#include <memory>
#include <stdint.h>
#include <string.h>
#include <utility>
#include <vector>

#include "../util/require.hpp"
#include "grammar.hpp"

namespace camspork
{

template <typename T>
class VarSlotEntry
{
    T* p_data;
    T self_data;  // p_data = &self_data if of size 1; otherwise heap allocated.
    std::vector<extent_t> extent;

  public:
    VarSlotEntry(std::vector<extent_t> extent_arg)
    {
        const size_t alloc_size = get_alloc_size(extent_arg);
        if (alloc_size <= 1) {
            p_data = &self_data;
        }
        else {
            p_data = new T[alloc_size];
        }
        extent = std::move(extent_arg);
    };

    VarSlotEntry(T scalar_init = 0)
    {
        p_data = &self_data;
        self_data = scalar_init;
    }

    ~VarSlotEntry()
    {
        free_if_allocated();
    }

    VarSlotEntry(const VarSlotEntry& other)
    {
        const size_t alloc_size = get_alloc_size(other.extent);
        p_data = other.p_data == &other.self_data ? &self_data : new T[alloc_size];
        memcpy(p_data, other.p_data, alloc_size * sizeof(p_data[0]));
        extent = other.extent;
    }

    VarSlotEntry(VarSlotEntry&& other) noexcept
    {
        move_from(std::move(other));
    };

    VarSlotEntry& operator=(VarSlotEntry other) noexcept
    {
        free_if_allocated();
        move_from(std::move(other));
        return *this;
    }

    // C-style multidimensional array indexing.
    // Provide the indices as a [begin, end) iterator pair.
    template <typename IdxIterator>
    T& idx(const IdxIterator begin, const IdxIterator end)
    {
        CAMSPORK_REQUIRE_CMP(extent.size(), ==, size_t(end - begin), "wrong index count used to read VarSlotEntry");
        size_t linear_idx = 0;
        for (IdxIterator iter = begin ; iter != end; ++iter) {
            const auto dim = end - begin;
            const size_t idx = *iter;
            CAMSPORK_REQUIRE_CMP(idx, <, extent[dim], "out-of-bounds access in abstract machine program");
            linear_idx = linear_idx * extent[dim] + idx;
        }
        return p_data[linear_idx];
    }

    template <typename IdxIterator>
    const T& idx(const IdxIterator begin, const IdxIterator end) const
    {
        return const_cast<VarSlotEntry*>(this)->idx(begin, end);
    }

    T& scalar()
    {
        CAMSPORK_REQUIRE_CMP(extent.size(), ==, 0, "tried to read tensor VarSlotEntry as scalar");
        return p_data[0];
    }

    const T& scalar() const
    {
        CAMSPORK_REQUIRE_CMP(extent.size(), ==, 0, "tried to read tensor VarSlotEntry as scalar");
        return p_data[0];
    }

  private:
    static size_t get_alloc_size(const std::vector<extent_t>& extent_arg)
    {
        size_t prod = 1;
        for (size_t n : extent_arg) {
            prod *= n;
        }
        return prod;
    }

    void free_if_allocated()
    {
        if (p_data && p_data != &self_data) {
            delete[] p_data;
        }
    }

    void move_from(VarSlotEntry&& other)
    {
        p_data = other.p_data == &other.self_data ? &self_data : other.p_data;
        self_data = other.self_data;
        extent = std::move(other.extent);

        other.extent.clear();
        other.p_data = &other.self_data;
        other.self_data = 0;
    }
};

class ProgramEnv
{
    size_t program_buffer_size;
    std::shared_ptr<char[]> p_program_buffer;
    const ProgramHeader& header;  // Validated from p_program_buffer
    std::vector<std::string> variable_names;
    std::vector<VarSlotEntry<value_t>> value_env_slots;
    std::vector<VarSlotEntry<sync_id_t>> sync_env_slots;
    value_t task_index = 0;

  public:
    friend class ProgramExec;

    ProgramEnv(size_t buffer_size, const char* buffer);

    // Currently moves are the same as copies.
    ProgramEnv(const ProgramEnv&) = default;
    ProgramEnv& operator=(const ProgramEnv&) = default;
    ~ProgramEnv() = default;

    void exec()
    {
        exec(header.top_level_stmt);
    }

    void exec(StmtRef stmt);

    void alloc_values(Varname name, std::vector<extent_t> extent)
    {
        CAMSPORK_C_BOUNDSCHECK(name.slot, value_env_slots.size());
        value_env_slots[name.slot] = VarSlotEntry<int32_t>(std::move(extent));
    }

    void alloc_scalar_value(Varname name, int32_t value)
    {
        CAMSPORK_C_BOUNDSCHECK(name.slot, value_env_slots.size());
        value_env_slots[name.slot] = VarSlotEntry<int32_t>(value);
    }

    VarSlotEntry<int32_t>& value_slot(Varname name)
    {
        CAMSPORK_C_BOUNDSCHECK(name.slot, value_env_slots.size());
        return value_env_slots[name.slot];
    }

    const VarSlotEntry<int32_t>& value_slot(Varname name) const
    {
        CAMSPORK_C_BOUNDSCHECK(name.slot, value_env_slots.size());
        return value_env_slots[name.slot];
    }

  private:
    std::shared_ptr<char[]> make_shared_program_buffer(size_t buffer_size, const char* buffer)
    {
        std::shared_ptr<char[]> p_result(new char[buffer_size]);
        memcpy(p_result.get(), buffer, buffer_size);
        return p_result;
    }
};

}
