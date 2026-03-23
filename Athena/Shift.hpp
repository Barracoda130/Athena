#pragma once

#include "mpfr-impl.h"

mp_limb_t
atn_rshift( mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt );