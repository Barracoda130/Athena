#include "mpfr-impl.h"

MPFR_THREAD_VAR( mpfr_flags_t, __gmpfr_flags, 0 )
MPFR_THREAD_VAR( mpfr_exp_t, __gmpfr_emin, MPFR_EMIN_DEFAULT )
MPFR_THREAD_VAR( mpfr_exp_t, __gmpfr_emax, MPFR_EMAX_DEFAULT )