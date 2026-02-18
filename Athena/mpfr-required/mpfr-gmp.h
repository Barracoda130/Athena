#ifndef __GMPFR_GMP_H__
#define __GMPFR_GMP_H__

#if defined(__cplusplus)
extern "C" {
#endif

void
mpfr_assert_fail( const char* filename, int linenum,
                  const char* expr );

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

#define TMP_ALLOC(n) alloca (n)

#endif