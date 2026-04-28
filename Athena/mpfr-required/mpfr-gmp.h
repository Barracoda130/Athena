// ALL CODE IN THE FOLLOWING FILE WAS COPIED FROM GMP, WHICH IS LICENSED UNDER THE GNU LGPL, VERSION 3 OR LATER.
// Below is the original copyright notice from the GMP library, which applies to this file.

/*
Copyright 1991, 1993, 1994, 1996, 2000 - 2002 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and /or modify
it under the terms of either :

*the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    or

    *the GNU General Public License as published by the Free Software
    Foundation; either version 2 of the License, or (at your option) any
    later version.

    or both in parallel, as here.

    The GNU MP Library is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License
    for more details.

        You should have received copies of the GNU General Public License and the
        GNU Lesser General Public License along with the GNU MP Library.If not,
        see https ://www.gnu.org/licenses/.  */

#ifndef __ATN_GMPFR_GMP_H__
#define __ATN_GMPFR_GMP_H__
#include "mpfr-impl.h"

#if defined(__cplusplus)
extern "C" {
#endif

void
mpfr_assert_fail( const char* filename, int linenum,
                  const char* expr );

struct tmp_marker
{
    void* ptr;
    size_t size;
    struct tmp_marker* next;
};

# define _MPFR_PROTO(x) ()
#define __GMP_DECLSPEC_EXPORT  __declspec(dllexport)
//# define __MPFR_DECLSPEC __GMP_DECLSPEC_EXPORT

__MPFR_DECLSPEC void* atn_tmp_allocate (struct tmp_marker**,
    size_t);
__MPFR_DECLSPEC void atn_tmp_free (struct tmp_marker*);

/* Do not define TMP_SALLOC (see the test in mpfr-impl.h)! */
#define TMP_ALLOC(n) (MPFR_LIKELY ((n) < 16384) ?       \
                      alloca (n) : atn_tmp_allocate (&tmp_marker, (n)))
#define TMP_DECL(m) struct tmp_marker *tmp_marker
#define TMP_MARK(m) (tmp_marker = 0)
#define TMP_FREE(m) atn_tmp_free (tmp_marker)

#if defined(__cplusplus)
}
#endif

#define ASSERT_FAIL(expr)  mpfr_assert_fail (__FILE__, __LINE__, #expr)

#define MPN_ZERO(dst, n) memset((dst), 0, (n)*MPFR_BYTES_PER_MP_LIMB)
#define MPN_COPY(dst,src,n) \
  do                                                                  \
    {                                                                 \
      if ((dst) != (src))                                             \
        {                                                             \
          MPFR_ASSERTD ((char *) (dst) >= (char *) (src) +            \
                                     (n) * MPFR_BYTES_PER_MP_LIMB ||  \
                        (char *) (src) >= (char *) (dst) +            \
                                     (n) * MPFR_BYTES_PER_MP_LIMB);   \
          memcpy ((dst), (src), (n) * MPFR_BYTES_PER_MP_LIMB);        \
        }                                                             \
    }                                                                 \
  while (0)

#define MP_LIMB_T_MAX (~(mp_limb_t)0)
typedef mp_limb_t UWtype;
#define W_TYPE_SIZE GMP_NUMB_BITS
#define MPN_COPY_DECR(dst,src,n) memmove((dst),(src),(n)*MPFR_BYTES_PER_MP_LIMB)

#endif