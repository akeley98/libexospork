#pragma once

#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string>

namespace camspork
{

std::string& thread_local_message_ref();

template <typename LHS, typename RHS>
[[noreturn]] void camspork_require_cmp_fail(
    const char* file, int line,
    LHS lhs, RHS rhs,
    const char* lhs_name, const char* op, const char* rhs_name, const char* msg)
{
    try {
        std::stringstream s;
        s << "Failed requirement: " << lhs_name << " " << op << " " << rhs_name;
        s << " (with " << lhs_name << " = " << lhs << "; ";
        s << rhs_name << " = " << rhs << "; (";
        s << msg << " @ " << file << ":" << line << ")";
        thread_local_message_ref() = s.str();
    }
    catch (...) {
        fprintf(stderr, "Exception while formatting error: %s:%i\n", __FILE__, __LINE__);
    }
    throw std::runtime_error("camspork_require_cmp_fail");
}

}

#define CAMSPORK_REQUIRE_CMP(lhs, op, rhs, msg) do { \
    if (!(lhs op rhs)) { \
        camspork_require_cmp_fail(__FILE__, __LINE__, lhs, rhs, #lhs, #op, #rhs, msg); \
    } \
} while (0)

#define CAMSPORK_C_BOUNDSCHECK(i, bounds) \
    CAMSPORK_REQUIRE_CMP(size_t(i), <, size_t(bounds), "this may be a bug in the C library, or a corrupt input file")
