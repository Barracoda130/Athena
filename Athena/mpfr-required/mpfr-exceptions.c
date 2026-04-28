// ALL CODE IN THE FOLLOWING FILE WAS COPIED FROM MPFR 4.2.0, WHICH IS LICENSED UNDER THE GNU LGPL, VERSION 3 OR LATER.
// Below is the original copyright notice from the MPFR library, which applies to this file.

/* Utilities for MPFR developers, not exported.

Copyright 1999-2025 Free Software Foundation, Inc.
Contributed by the Pascaline and Caramba projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.
If not, see <https://www.gnu.org/licenses/>. */

#include "mpfr-impl.h"

MPFR_THREAD_VAR( mpfr_flags_t, __gmpfr_flags, 0 )
MPFR_THREAD_VAR( mpfr_exp_t, __gmpfr_emin, MPFR_EMIN_DEFAULT )
MPFR_THREAD_VAR( mpfr_exp_t, __gmpfr_emax, MPFR_EMAX_DEFAULT )