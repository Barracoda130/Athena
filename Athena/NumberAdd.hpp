#pragma once

#include "MpfrInclude.hpp"

int atn_add1( mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, Athena::round_t rnd_mode );
int atn_add1sp( mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_rnd_t rnd_mode );