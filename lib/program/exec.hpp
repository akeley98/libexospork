#pragma once

#include <memory>
#include <stdint.h>
#include <string.h>
#include <utility>
#include <vector>

#include "grammar.hpp"
#include "../syncv/syncv_table.hpp"
#include "../syncv/syncv_types.hpp"
#include "../util/require.hpp"

namespace camspork
{

template <typename T>
class VarSlotEntry
{
    T* p_data;
    T self_data;  // p_data = &self_data if of size 1; otherwise heap allocated.
    std::vector<extent_t> _extent;

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
        _extent = std::move(extent_arg);
    };

    VarSlotEntry(T scalar_init = T{})
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
        const size_t alloc_size = get_alloc_size(other._extent);
        p_data = other.p_data == &other.self_data ? &self_data : new T[alloc_size];
        memcpy(p_data, other.p_data, alloc_size * sizeof(p_data[0]));
        _extent = other._extent;
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

    void reset()
    {
        free_if_allocated();
        p_data = &self_data;
        _extent.clear();
        self_data = T{};
    }

    const std::vector<extent_t>& extent() const
    {
        return _extent;
    }

    // C-style multidimensional array indexing.
    // Provide the indices as a [begin, end) iterator pair.
    template <typename IdxIterator>
    T& idx(const IdxIterator begin, const IdxIterator end)
    {
        CAMSPORK_REQUIRE_CMP(_extent.size(), ==, size_t(end - begin), "wrong index count used to read VarSlotEntry");
        size_t linear_idx = 0;
        for (IdxIterator iter = begin ; iter != end; ++iter) {
            const auto dim = iter - begin;
            const size_t idx = *iter;
            CAMSPORK_REQUIRE_CMP(idx, <, _extent[dim], "out-of-bounds access in abstract machine program");
            linear_idx = linear_idx * _extent[dim] + idx;
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
        CAMSPORK_REQUIRE_CMP(_extent.size(), ==, 0, "tried to read tensor VarSlotEntry as scalar");
        return p_data[0];
    }

    const T& scalar() const
    {
        CAMSPORK_REQUIRE_CMP(_extent.size(), ==, 0, "tried to read tensor VarSlotEntry as scalar");
        return p_data[0];
    }

    T* data()
    {
        return p_data;
    }

    const T* data() const
    {
        return p_data;
    }

    size_t size() const
    {
        return get_alloc_size(_extent);
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
        _extent = std::move(other._extent);

        other._extent.clear();
        other.p_data = &other.self_data;
        other.self_data = T{};
    }
};

struct VarSlotEnvs
{
    std::string name;
    // For this variable name: Value env, Sync env, Barrier env.
    // Conceptually, this is mapping from indices to data.
    VarSlotEntry<value_t> value;
    VarSlotEntry<assignment_record_id> sync;
    VarSlotEntry<barrier_id> barrier;
};

class ProgramEnvSyncvTable
{
    SyncvTable* raw_ptr;
  public:
    ProgramEnvSyncvTable(SyncvTable* arg = nullptr) : raw_ptr(arg)
    {
    }
    ProgramEnvSyncvTable(const ProgramEnvSyncvTable& other)
    {
        raw_ptr = other.raw_ptr ? copy_syncv_table(other.raw_ptr) : nullptr;
    }
    ProgramEnvSyncvTable& operator=(ProgramEnvSyncvTable other)
    {
        std::swap(raw_ptr, other.raw_ptr);
        return *this;
    }
    ~ProgramEnvSyncvTable()
    {
        if (raw_ptr) {
            delete_syncv_table(raw_ptr);
        }
    }

    SyncvTable* get() const
    {
        return raw_ptr;
    }
};

class ProgramEnv
{
    size_t program_buffer_size;
    std::shared_ptr<char[]> p_program_buffer;
    const ProgramHeader& header;  // Validated from p_program_buffer
    ProgramEnvSyncvTable p_syncv_table;
    ThreadCuboid thread_cuboid;
    std::vector<VarSlotEnvs> var_slots;

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
        CAMSPORK_C_BOUNDSCHECK(name.slot(), var_slots.size());
        var_slots[name.slot()].value = VarSlotEntry<value_t>(std::move(extent));
    }

    void alloc_scalar_value(Varname name, value_t value)
    {
        CAMSPORK_C_BOUNDSCHECK(name.slot(), var_slots.size());
        var_slots[name.slot()].value = VarSlotEntry<value_t>(value);
    }

    const std::string& get_varname(Varname name)
    {
        CAMSPORK_C_BOUNDSCHECK(name.slot(), var_slots.size());
        return var_slots[name.slot()].name;
    }

    VarSlotEntry<value_t>& value_slot(Varname name)
    {
        CAMSPORK_C_BOUNDSCHECK(name.slot(), var_slots.size());
        return var_slots[name.slot()].value;
    }

    const VarSlotEntry<value_t>& value_slot(Varname name) const
    {
        CAMSPORK_C_BOUNDSCHECK(name.slot(), var_slots.size());
        return var_slots[name.slot()].value;
    }

    void alloc_sync(Varname name, std::vector<extent_t> extent)
    {
        CAMSPORK_C_BOUNDSCHECK(name.slot(), var_slots.size());
        var_slots[name.slot()].sync = VarSlotEntry<assignment_record_id>(std::move(extent));
    }

    VarSlotEntry<assignment_record_id>& sync_slot(Varname name)
    {
        CAMSPORK_C_BOUNDSCHECK(name.slot(), var_slots.size());
        return var_slots[name.slot()].sync;
    }

    const VarSlotEntry<assignment_record_id>& sync_slot(Varname name) const
    {
        CAMSPORK_C_BOUNDSCHECK(name.slot(), var_slots.size());
        return var_slots[name.slot()].sync;
    }

    VarSlotEntry<barrier_id>& barrier_slot(Varname name)
    {
        CAMSPORK_C_BOUNDSCHECK(name.slot(), var_slots.size());
        return var_slots[name.slot()].barrier;
    }

    const VarSlotEntry<barrier_id>& barrier_slot(Varname name) const
    {
        CAMSPORK_C_BOUNDSCHECK(name.slot(), var_slots.size());
        return var_slots[name.slot()].barrier;
    }

  private:
    std::shared_ptr<char[]> make_shared_program_buffer(size_t buffer_size, const char* buffer)
    {
        std::shared_ptr<char[]> p_result(new char[buffer_size]);
        memcpy(p_result.get(), buffer, buffer_size);
        return p_result;
    }
};

}  // end namespace camspork

// Return 0
