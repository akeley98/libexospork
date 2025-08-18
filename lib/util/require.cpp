#include "require.hpp"

namespace camspork
{

void camspork_require_fail(
    const char* file, int line, const char* expr_str, const char* msg)
{
    try {
        std::stringstream s;
        s << "Failed requirement: " << expr_str;
        s << " (" << msg << " @ " << file << ":" << line << ")";
        thread_local_message_ref() = s.str();
    }
    catch (...) {
        fprintf(stderr, "Exception while formatting error: %s:%i\n", __FILE__, __LINE__);
    }
    throw std::runtime_error("camspork_require_fail");
}

std::string& thread_local_message_ref()
{
    thread_local std::string s;
    return s;
}

}
