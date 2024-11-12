#include "../../include/exospork/syncv.h"

#include "sigthread.hpp"

static_assert(EXOSPORK_SYNC_ACCESS_BIT == exospork::SigthreadInterval::sync_bit);

extern "C" void exospork_syncv_init(size_t sizeof_struct, const exospork_syncv_init_t*)
{

}
