#include "require.hpp"

namespace camspork
{

std::string& thread_local_message_ref()
{
    thread_local std::string s;
    return s;
}

}
