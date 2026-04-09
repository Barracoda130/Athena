#include "Number.hpp"

#include "Athena.hpp"

#include <cassert>
#include <iostream>
#include <chrono>
#include <random>
#include <iomanip>
#include <algorithm>
#include <bit>
#include <limits>
#include <type_traits>

#include "mpfr-impl.h"
#include "Shift.hpp"

#include "IntrinsicMath.hpp"
#include "NumberAdd.hpp"

#if defined(_MSC_VER)
#include <intrin.h>
#include <immintrin.h>
#endif

using namespace Athena;


#define MPN_OVERLAP_P(xp, xsize, yp, ysize)				\
  ((xp) + (xsize) > (yp) && (yp) + (ysize) > (xp))
#define MPN_SAME_OR_INCR2_P(dst, dsize, src, ssize)			\
  ((dst) <= (src) || ! MPN_OVERLAP_P (dst, dsize, src, ssize))
#define MPN_SAME_OR_INCR_P(dst, src, size)				\
  MPN_SAME_OR_INCR2_P(dst, size, src, size)

# define DEBUG(x)

namespace
{
	inline mpfr_prec_t countLeadingZeros( mp_limb_t value )
	{
		using UnsignedLimb = std::make_unsigned_t<mp_limb_t>;
		const auto limbValue = static_cast<UnsignedLimb>( value );
		return static_cast<mpfr_prec_t>( std::countl_zero( limbValue ) );
	}
}

int
atn_sub1sp( mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_rnd_t rnd_mode )
{
	mpfr_exp_t bx, cx;
	mpfr_uexp_t d;
	mpfr_prec_t p, sh, cnt;
	mp_size_t n;
	mp_limb_t* ap, * bp, * cp;
	mp_limb_t limb;
	int inexact;
	mp_limb_t bcp, bcp1; /* Cp and C'p+1 */
	mp_limb_t bbcp = (mp_limb_t)-1, bbcp1 = (mp_limb_t)-1; /* Cp+1 and C'p+2,
	  gcc claims that they might be used uninitialized. We fill them with invalid
	  values, which should produce a failure if so. See README.dev file. */

	MPFR_TMP_DECL( marker );

	MPFR_TMP_MARK( marker );

	MPFR_ASSERTD( MPFR_PREC( a ) == MPFR_PREC( b ) && MPFR_PREC( b ) == MPFR_PREC( c ) );
	MPFR_ASSERTD( MPFR_IS_PURE_FP( b ) );
	MPFR_ASSERTD( MPFR_IS_PURE_FP( c ) );

	/* Read prec and num of limbs */
	p = MPFR_PREC( b );
	n = MPFR_PREC2LIMBS( p );

	/* Fast cmp of |b| and |c|*/
	bx = MPFR_GET_EXP( b );
	cx = MPFR_GET_EXP( c );
	if ( MPFR_UNLIKELY( bx == cx ) )
	{
		mp_size_t k = n - 1;
		/* Check mantissa since exponent are equals */
		bp = MPFR_MANT( b );
		cp = MPFR_MANT( c );
		while ( k >= 0 && MPFR_UNLIKELY( bp[k] == cp[k] ) )
			k--;
		if ( MPFR_UNLIKELY( k < 0 ) )
			/* b == c ! */
		{
			/* Return exact number 0 */
			if ( rnd_mode == MPFR_RNDD )
				MPFR_SET_NEG( a );
			else
				MPFR_SET_POS( a );
			MPFR_SET_ZERO( a );
			MPFR_RET( 0 );
		}
		else if ( bp[k] > cp[k] )
			goto BGreater;
		else
		{
			MPFR_ASSERTD( bp[k] < cp[k] );
			goto CGreater;
		}
	}
	else if ( MPFR_UNLIKELY( bx < cx ) )
	{
		/* Swap b and c and set sign */
		mpfr_srcptr t;
		mpfr_exp_t tx;
	CGreater:
		MPFR_SET_OPPOSITE_SIGN( a, b );
		t = b;  b = c;  c = t;
		tx = bx; bx = cx; cx = tx;
	}
	else
	{
		/* b > c */
	BGreater:
		MPFR_SET_SAME_SIGN( a, b );
	}

	/* Now b > c */
	MPFR_ASSERTD( bx >= cx );
	d = (mpfr_uexp_t)bx - cx;
	DEBUG( printf( "New with diff=%lu\n", (unsigned long)d ) );

	if ( MPFR_UNLIKELY( d <= 1 ) )
	{
		if ( MPFR_LIKELY( d < 1 ) )
		{
			/* <-- b -->
			   <-- c --> : exact sub */
			ap = MPFR_MANT( a );
			atn_sub_n( ap, MPFR_MANT( b ), MPFR_MANT( c ), n );
			/* Normalize */
		ExactNormalize:
			limb = ap[n - 1];
			if ( MPFR_LIKELY( limb ) )
			{
				/* First limb is not zero. */
               cnt = countLeadingZeros( limb );
				/* cnt could be == 0 <= SubD1Lose */
				if ( MPFR_LIKELY( cnt ) )
				{
					mpn_lshift( ap, ap, n, cnt ); /* Normalize number */
					bx -= cnt; /* Update final expo */
				}
				/* Last limb should be ok */
				MPFR_ASSERTD( !(ap[0] & MPFR_LIMB_MASK( (unsigned int)(-p)
					% GMP_NUMB_BITS )) );
			}
			else
			{
				/* First limb is zero */
				mp_size_t k = n - 1, len;
				/* Find the first limb not equal to zero.
				   FIXME:It is assume it exists (since |b| > |c| and same prec)*/
				do
				{
					MPFR_ASSERTD( k > 0 );
					limb = ap[--k];
				} while ( limb == 0 );
				MPFR_ASSERTD( limb != 0 );
               cnt = countLeadingZeros( limb );
				k++;
				len = n - k; /* Number of last limb */
				MPFR_ASSERTD( k >= 0 );
				if ( MPFR_LIKELY( cnt ) )
					mpn_lshift( ap + len, ap, k, cnt ); /* Normalize the High Limb*/
				else
				{
					/* Must use DECR since src and dest may overlap & dest>=src*/
					MPN_COPY_DECR( ap + len, ap, k );
				}
				MPN_ZERO( ap, len ); /* Zeroing the last limbs */
				bx -= cnt + len * GMP_NUMB_BITS; /* Update Expo */
				/* Last limb should be ok */
				MPFR_ASSERTD( !(ap[len] & MPFR_LIMB_MASK( (unsigned int)(-p)
					% GMP_NUMB_BITS )) );
			}
			/* Check expo underflow */
			if ( MPFR_UNLIKELY( bx < __gmpfr_emin ) )
			{
				MPFR_TMP_FREE( marker );
				/* inexact=0 */
				DEBUG( printf( "(D==0 Underflow)\n" ) );
				if ( rnd_mode == MPFR_RNDN &&
					(bx < __gmpfr_emin - 1 ||
						(/*inexact >= 0 &&*/ mpfr_powerof2_raw( a ))) )
					rnd_mode = MPFR_RNDZ;
				return mpfr_underflow( a, rnd_mode, MPFR_SIGN( a ) );
			}
			MPFR_SET_EXP( a, bx );
			/* No rounding is necessary since the result is exact */
			MPFR_ASSERTD( ap[n - 1] > ~ap[n - 1] );
			MPFR_TMP_FREE( marker );
			return 0;
		}
		else /* if (d == 1) */
		{
			/* | <-- b -->
			   |  <-- c --> */
			mp_limb_t c0, mask;
			mp_size_t k;
			MPFR_UNSIGNED_MINUS_MODULO( sh, p );
			/* If we lose at least one bit, compute 2*b-c (Exact)
			 * else compute b-c/2 */
			bp = MPFR_MANT( b );
			cp = MPFR_MANT( c );
			k = n - 1;
			limb = bp[k] - cp[k] / 2;
			if ( limb > MPFR_LIMB_HIGHBIT )
			{
				/* We can't lose precision: compute b-c/2 */
				/* Shift c in the allocated temporary block */
			SubD1NoLose:
				c0 = cp[0] & (MPFR_LIMB_ONE << sh);
				cp = MPFR_TMP_LIMBS_ALLOC( n );
				mpn_rshift( cp, MPFR_MANT( c ), n, 1 );
				if ( MPFR_LIKELY( c0 == 0 ) )
				{
					/* Result is exact: no need of rounding! */
					ap = MPFR_MANT( a );
					atn_sub_n( ap, bp, cp, n );
					MPFR_SET_EXP( a, bx ); /* No expo overflow! */
					/* No truncate or normalize is needed */
					MPFR_ASSERTD( ap[n - 1] > ~ap[n - 1] );
					/* No rounding is necessary since the result is exact */
					MPFR_TMP_FREE( marker );
					return 0;
				}
				ap = MPFR_MANT( a );
				mask = ~MPFR_LIMB_MASK( sh );
				cp[0] &= mask; /* Delete last bit of c */
				atn_sub_n( ap, bp, cp, n );
				MPFR_SET_EXP( a, bx );                 /* No expo overflow! */
				MPFR_ASSERTD( !(ap[0] & ~mask) );    /* Check last bits */
				/* No normalize is needed */
				MPFR_ASSERTD( ap[n - 1] > ~ap[n - 1] );
				/* Rounding is necessary since c0 = 1*/
				/* Cp =-1 and C'p+1=0 */
				bcp = 1; bcp1 = 0;
				if ( MPFR_LIKELY( rnd_mode == MPFR_RNDN ) )
				{
					/* Even Rule apply: Check Ap-1 */
					if ( MPFR_LIKELY( (ap[0] & (MPFR_LIMB_ONE << sh)) == 0 ) )
						goto truncate;
					else
						goto sub_one_ulp;
				}
				MPFR_UPDATE_RND_MODE( rnd_mode, MPFR_IS_NEG( a ) );
				if ( rnd_mode == MPFR_RNDZ )
					goto sub_one_ulp;
				else
					goto truncate;
			}
			else if ( MPFR_LIKELY( limb < MPFR_LIMB_HIGHBIT ) )
			{
				/* We lose at least one bit of prec */
				/* Calcul of 2*b-c (Exact) */
				/* Shift b in the allocated temporary block */
			SubD1Lose:
				bp = MPFR_TMP_LIMBS_ALLOC( n );
				mpn_lshift( bp, MPFR_MANT( b ), n, 1 );
				ap = MPFR_MANT( a );
				atn_sub_n( ap, bp, cp, n );
				bx--;
				goto ExactNormalize;
			}
			else
			{
				/* Case: limb = 100000000000 */
				/* Check while b[k] == c'[k] (C' is C shifted by 1) */
				/* If b[k]<c'[k] => We lose at least one bit*/
				/* If b[k]>c'[k] => We don't lose any bit */
				/* If k==-1 => We don't lose any bit
				   AND the result is 100000000000 0000000000 00000000000 */
				mp_limb_t carry;
				do
				{
					carry = cp[k] & MPFR_LIMB_ONE;
					k--;
				} while ( k >= 0 &&
					bp[k] == (carry = cp[k] / 2 + (carry << (GMP_NUMB_BITS - 1))) );
				if ( MPFR_UNLIKELY( k < 0 ) )
				{
					/*If carry then (sh==0 and Virtual c'[-1] > Virtual b[-1]) */
					if ( MPFR_UNLIKELY( carry ) ) /* carry = cp[0]&MPFR_LIMB_ONE */
					{
						/* FIXME: Can be faster? */
						MPFR_ASSERTD( sh == 0 );
						goto SubD1Lose;
					}
					/* Result is a power of 2 */
					ap = MPFR_MANT( a );
					MPN_ZERO( ap, n );
					ap[n - 1] = MPFR_LIMB_HIGHBIT;
					MPFR_SET_EXP( a, bx ); /* No expo overflow! */
					/* No Normalize is needed*/
					/* No Rounding is needed */
					MPFR_TMP_FREE( marker );
					return 0;
				}
				/* carry = cp[k]/2+(cp[k-1]&1)<<(GMP_NUMB_BITS-1) = c'[k]*/
				else if ( bp[k] > carry )
					goto SubD1NoLose;
				else
				{
					MPFR_ASSERTD( bp[k] < carry );
					goto SubD1Lose;
				}
			}
		}
	}
	else if ( MPFR_UNLIKELY( d >= p ) )
	{
		ap = MPFR_MANT( a );
		MPFR_UNSIGNED_MINUS_MODULO( sh, p );
		/* We can't set A before since we use cp for rounding... */
		/* Perform rounding: check if a=b or a=b-ulp(b) */
		if ( MPFR_UNLIKELY( d == p ) )
		{
			/* cp == -1 and c'p+1 = ? */
			bcp = 1;
			/* We need Cp+1 later for a very improbable case. */
			bbcp = (MPFR_MANT( c )[n - 1] & (MPFR_LIMB_ONE << (GMP_NUMB_BITS - 2)));
			/* We need also C'p+1 for an even more unprobable case... */
			if ( MPFR_LIKELY( bbcp ) )
				bcp1 = 1;
			else
			{
				cp = MPFR_MANT( c );
				if ( MPFR_UNLIKELY( cp[n - 1] == MPFR_LIMB_HIGHBIT ) )
				{
					mp_size_t k = n - 1;
					do
					{
						k--;
					} while ( k >= 0 && cp[k] == 0 );
					bcp1 = (k >= 0);
				}
				else
					bcp1 = 1;
			}
			DEBUG( printf( "(D=P) Cp=-1 Cp+1=%d C'p+1=%d \n", bbcp != 0, bcp1 != 0 ) );
			bp = MPFR_MANT( b );

			/* Even if src and dest overlap, it is ok using MPN_COPY */
			if ( MPFR_LIKELY( rnd_mode == MPFR_RNDN ) )
			{
				if ( MPFR_UNLIKELY( bcp && bcp1 == 0 ) )
					/* Cp=-1 and C'p+1=0: Even rule Apply! */
					/* Check Ap-1 = Bp-1 */
					if ( (bp[0] & (MPFR_LIMB_ONE << sh)) == 0 )
					{
						MPN_COPY( ap, bp, n );
						goto truncate;
					}
				MPN_COPY( ap, bp, n );
				goto sub_one_ulp;
			}
			MPFR_UPDATE_RND_MODE( rnd_mode, MPFR_IS_NEG( a ) );
			if ( rnd_mode == MPFR_RNDZ )
			{
				MPN_COPY( ap, bp, n );
				goto sub_one_ulp;
			}
			else
			{
				MPN_COPY( ap, bp, n );
				goto truncate;
			}
		}
		else
		{
			/* Cp=0, Cp+1=-1 if d==p+1, C'p+1=-1 */
			bcp = 0; bbcp = (d == p + 1); bcp1 = 1;
			DEBUG( printf( "(D>P) Cp=%d Cp+1=%d C'p+1=%d\n", bcp != 0, bbcp != 0, bcp1 != 0 ) );
			/* Need to compute C'p+2 if d==p+1 and if rnd_mode=NEAREST
			   (Because of a very improbable case) */
			if ( MPFR_UNLIKELY( d == p + 1 && rnd_mode == MPFR_RNDN ) )
			{
				cp = MPFR_MANT( c );
				if ( MPFR_UNLIKELY( cp[n - 1] == MPFR_LIMB_HIGHBIT ) )
				{
					mp_size_t k = n - 1;
					do
					{
						k--;
					} while ( k >= 0 && cp[k] == 0 );
					bbcp1 = (k >= 0);
				}
				else
					bbcp1 = 1;
				DEBUG( printf( "(D>P) C'p+2=%d\n", bbcp1 != 0 ) );
			}
			/* Copy mantissa B in A */
			MPN_COPY( ap, MPFR_MANT( b ), n );
			/* Round */
			if ( MPFR_LIKELY( rnd_mode == MPFR_RNDN ) )
				goto truncate;
			MPFR_UPDATE_RND_MODE( rnd_mode, MPFR_IS_NEG( a ) );
			if ( rnd_mode == MPFR_RNDZ )
				goto sub_one_ulp;
			else /* rnd_mode = AWAY */
				goto truncate;
		}
	}
	else
	{
		mpfr_uexp_t dm;
		mp_size_t m;
		mp_limb_t mask;

		/* General case: 2 <= d < p */
		MPFR_UNSIGNED_MINUS_MODULO( sh, p );
		cp = MPFR_TMP_LIMBS_ALLOC( n );

		/* Shift c in temporary allocated place */
		dm = d % GMP_NUMB_BITS;
		m = d / GMP_NUMB_BITS;
		if ( MPFR_UNLIKELY( dm == 0 ) )
		{
			/* dm = 0 and m > 0: Just copy */
			MPFR_ASSERTD( m != 0 );
			MPN_COPY( cp, MPFR_MANT( c ) + m, n - m );
			MPN_ZERO( cp + n - m, m );
		}
		else if ( MPFR_LIKELY( m == 0 ) )
		{
			/* dm >=2 and m == 0: just shift */
			MPFR_ASSERTD( dm >= 2 );
			mpn_rshift( cp, MPFR_MANT( c ), n, dm );
		}
		else
		{
			/* dm > 0 and m > 0: shift and zero  */
			mpn_rshift( cp, MPFR_MANT( c ) + m, n - m, dm );
			MPN_ZERO( cp + n - m, m );
		}

		DEBUG( mpfr_print_mant_binary( "Before", MPFR_MANT( c ), p ) );
		DEBUG( mpfr_print_mant_binary( "B=    ", MPFR_MANT( b ), p ) );
		DEBUG( mpfr_print_mant_binary( "After ", cp, p ) );

		/* Compute bcp=Cp and bcp1=C'p+1 */
		if ( MPFR_LIKELY( sh ) )
		{
			/* Try to compute them from C' rather than C (FIXME: Faster?) */
			bcp = (cp[0] & (MPFR_LIMB_ONE << (sh - 1)));
			if ( MPFR_LIKELY( cp[0] & MPFR_LIMB_MASK( sh - 1 ) ) )
				bcp1 = 1;
			else
			{
				/* We can't compute C'p+1 from C'. Compute it from C */
				/* Start from bit x=p-d+sh in mantissa C
				   (+sh since we have already looked sh bits in C'!) */
				mpfr_prec_t x = p - d + sh - 1;
				if ( MPFR_LIKELY( x > p ) )
					/* We are already looked at all the bits of c, so C'p+1 = 0*/
					bcp1 = 0;
				else
				{
					mp_limb_t* tp = MPFR_MANT( c );
					mp_size_t kx = n - 1 - (x / GMP_NUMB_BITS);
					mpfr_prec_t sx = GMP_NUMB_BITS - 1 - (x % GMP_NUMB_BITS);
					DEBUG( printf( "(First) x=%lu Kx=%ld Sx=%lu\n",
						(unsigned long)x, (long)kx,
						(unsigned long)sx ) );
					/* Looks at the last bits of limb kx (if sx=0 does nothing)*/
					if ( tp[kx] & MPFR_LIMB_MASK( sx ) )
						bcp1 = 1;
					else
					{
						/*kx += (sx==0);*/
						/*If sx==0, tp[kx] hasn't been checked*/
						do
						{
							kx--;
						} while ( kx >= 0 && tp[kx] == 0 );
						bcp1 = (kx >= 0);
					}
				}
			}
		}
		else
		{
			/* Compute Cp and C'p+1 from C with sh=0 */
			mp_limb_t* tp = MPFR_MANT( c );
			/* Start from bit x=p-d in mantissa C */
			mpfr_prec_t  x = p - d;
			mp_size_t   kx = n - 1 - (x / GMP_NUMB_BITS);
			mpfr_prec_t sx = GMP_NUMB_BITS - 1 - (x % GMP_NUMB_BITS);
			MPFR_ASSERTD( p >= d );
			bcp = (tp[kx] & (MPFR_LIMB_ONE << sx));
			/* Looks at the last bits of limb kx (If sx=0, does nothing)*/
			if ( tp[kx] & MPFR_LIMB_MASK( sx ) )
				bcp1 = 1;
			else
			{
				/*kx += (sx==0);*/ /*If sx==0, tp[kx] hasn't been checked*/
				do
				{
					kx--;
				} while ( kx >= 0 && tp[kx] == 0 );
				bcp1 = (kx >= 0);
			}
		}
		DEBUG( printf( "sh=%lu Cp=%d C'p+1=%d\n", sh, bcp != 0, bcp1 != 0 ) );

		/* Check if we can lose a bit, and if so compute Cp+1 and C'p+2 */
		bp = MPFR_MANT( b );
		if ( MPFR_UNLIKELY( (bp[n - 1] - cp[n - 1]) <= MPFR_LIMB_HIGHBIT ) )
		{
			/* We can lose a bit so we precompute Cp+1 and C'p+2 */
			/* Test for trivial case: since C'p+1=0, Cp+1=0 and C'p+2 =0 */
			if ( MPFR_LIKELY( bcp1 == 0 ) )
			{
				bbcp = 0;
				bbcp1 = 0;
			}
			else /* bcp1 != 0 */
			{
				/* We can lose a bit:
				   compute Cp+1 and C'p+2 from mantissa C */
				mp_limb_t* tp = MPFR_MANT( c );
				/* Start from bit x=(p+1)-d in mantissa C */
				mpfr_prec_t x = p + 1 - d;
				mp_size_t kx = n - 1 - (x / GMP_NUMB_BITS);
				mpfr_prec_t sx = GMP_NUMB_BITS - 1 - (x % GMP_NUMB_BITS);
				MPFR_ASSERTD( p > d );
				DEBUG( printf( "(pre) x=%lu Kx=%ld Sx=%lu\n",
					(unsigned long)x, (long)kx,
					(unsigned long)sx ) );
				bbcp = (tp[kx] & (MPFR_LIMB_ONE << sx));
				/* Looks at the last bits of limb kx (If sx=0, does nothing)*/
				/* If Cp+1=0, since C'p+1!=0, C'p+2=1 ! */
				if ( MPFR_LIKELY( bbcp == 0 || (tp[kx] & MPFR_LIMB_MASK( sx )) ) )
					bbcp1 = 1;
				else
				{
					/*kx += (sx==0);*/ /*If sx==0, tp[kx] hasn't been checked*/
					do
					{
						kx--;
					} while ( kx >= 0 && tp[kx] == 0 );
					bbcp1 = (kx >= 0);
					DEBUG( printf( "(Pre) Scan done for %ld\n", (long)kx ) );
				}
			} /*End of Bcp1 != 0*/
			DEBUG( printf( "(Pre) Cp+1=%d C'p+2=%d\n", bbcp != 0, bbcp1 != 0 ) );
		} /* End of "can lose a bit" */

	  /* Clean shifted C' */
		mask = ~MPFR_LIMB_MASK( sh );
		cp[0] &= mask;

		/* Subtract the mantissa c from b in a */
		ap = MPFR_MANT( a );
		atn_sub_n( ap, bp, cp, n );
		DEBUG( mpfr_print_mant_binary( "Sub=  ", ap, p ) );

		/* Normalize: we lose at max one bit*/
		if ( MPFR_UNLIKELY( MPFR_LIMB_MSB( ap[n - 1] ) == 0 ) )
		{
			/* High bit is not set and we have to fix it! */
			/* Ap >= 010000xxx001 */
			mpn_lshift( ap, ap, n, 1 );
			/* Ap >= 100000xxx010 */
			if ( MPFR_UNLIKELY( bcp != 0 ) ) /* Check if Cp = -1 */
				/* Since Cp == -1, we have to substract one more */
			{
				mpn_sub_1( ap, ap, n, MPFR_LIMB_ONE << sh );
				MPFR_ASSERTD( MPFR_LIMB_MSB( ap[n - 1] ) != 0 );
			}
			/* Ap >= 10000xxx001 */
			/* Final exponent -1 since we have shifted the mantissa */
			bx--;
			/* Update bcp and bcp1 */
			MPFR_ASSERTN( bbcp != (mp_limb_t)-1 );
			MPFR_ASSERTN( bbcp1 != (mp_limb_t)-1 );
			bcp = bbcp;
			bcp1 = bbcp1;
			/* We dont't have anymore a valid Cp+1!
			   But since Ap >= 100000xxx001, the final sub can't unnormalize!*/
		}
		MPFR_ASSERTD( !(ap[0] & ~mask) );

		/* Rounding */
		if ( MPFR_LIKELY( rnd_mode == MPFR_RNDN ) )
		{
			if ( MPFR_LIKELY( bcp == 0 ) )
				goto truncate;
			else if ( (bcp1) || ((ap[0] & (MPFR_LIMB_ONE << sh)) != 0) )
				goto sub_one_ulp;
			else
				goto truncate;
		}

		/* Update rounding mode */
		MPFR_UPDATE_RND_MODE( rnd_mode, MPFR_IS_NEG( a ) );
		if ( rnd_mode == MPFR_RNDZ && (MPFR_LIKELY( bcp || bcp1 )) )
			goto sub_one_ulp;
		goto truncate;
	}
	MPFR_RET_NEVER_GO_HERE();

	/* Sub one ulp to the result */
sub_one_ulp:
	mpn_sub_1( ap, ap, n, MPFR_LIMB_ONE << sh );
	/* Result should be smaller than exact value: inexact=-1 */
	inexact = -1;
	/* Check normalisation */
	if ( MPFR_UNLIKELY( MPFR_LIMB_MSB( ap[n - 1] ) == 0 ) )
	{
		/* ap was a power of 2, and we lose a bit */
		/* Now it is 0111111111111111111[00000 */
		mpn_lshift( ap, ap, n, 1 );
		bx--;
		/* And the lost bit x depends on Cp+1, and Cp */
		/* Compute Cp+1 if it isn't already compute (ie d==1) */
		/* FIXME: Is this case possible? */
		if ( MPFR_UNLIKELY( d == 1 ) )
			bbcp = 0;
		DEBUG( printf( "(SubOneUlp)Cp=%d, Cp+1=%d C'p+1=%d\n", bcp != 0, bbcp != 0, bcp1 != 0 ) );
		/* Compute the last bit (Since we have shifted the mantissa)
		   we need one more bit!*/
		MPFR_ASSERTN( bbcp != (mp_limb_t)-1 );
		if ( (rnd_mode == MPFR_RNDZ && bcp == 0)
			|| (rnd_mode == MPFR_RNDN && bbcp == 0)
			|| (bcp && bcp1 == 0) ) /*Exact result*/
		{
			ap[0] |= MPFR_LIMB_ONE << sh;
			if ( rnd_mode == MPFR_RNDN )
				inexact = 1;
			DEBUG( printf( "(SubOneUlp) Last bit set\n" ) );
		}
		/* Result could be exact if C'p+1 = 0 and rnd == Zero
		   since we have had one more bit to the result */
		   /* Fixme: rnd_mode == MPFR_RNDZ needed ? */
		if ( bcp1 == 0 && rnd_mode == MPFR_RNDZ )
		{
			DEBUG( printf( "(SubOneUlp) Exact result\n" ) );
			inexact = 0;
		}
	}

	goto end_of_sub;

truncate:
	/* Check if the result is an exact power of 2: 100000000000
	   in which cases, we could have to do sub_one_ulp due to some nasty reasons:
	   If Result is a Power of 2:
		+ If rnd = AWAY,
		|  If Cp=-1 and C'p+1 = 0, SubOneUlp and the result is EXACT.
		   If Cp=-1 and C'p+1 =-1, SubOneUlp and the result is above.
		   Otherwise truncate
		+ If rnd = NEAREST,
		   If Cp= 0 and Cp+1  =-1 and C'p+2=-1, SubOneUlp and the result is above
		   If cp=-1 and C'p+1 = 0, SubOneUlp and the result is exact.
		   Otherwise truncate.
		X bit should always be set if SubOneUlp*/
	if ( MPFR_UNLIKELY( ap[n - 1] == MPFR_LIMB_HIGHBIT ) )
	{
		mp_size_t k = n - 1;
		do
		{
			k--;
		} while ( k >= 0 && ap[k] == 0 );
		if ( MPFR_UNLIKELY( k < 0 ) )
		{
			/* It is a power of 2! */
			/* Compute Cp+1 if it isn't already compute (ie d==1) */
			/* FIXME: Is this case possible? */
			if ( d == 1 )
				bbcp = 0;
			DEBUG( printf( "(Truncate) Cp=%d, Cp+1=%d C'p+1=%d C'p+2=%d\n", \
				bcp != 0, bbcp != 0, bcp1 != 0, bbcp1 != 0 ) );
			MPFR_ASSERTN( bbcp != (mp_limb_t)-1 );
			MPFR_ASSERTN( (rnd_mode != MPFR_RNDN) || (bcp != 0) || (bbcp == 0) || (bbcp1 != (mp_limb_t)-1) );
			if ( ((rnd_mode != MPFR_RNDZ) && bcp)
				||
				((rnd_mode == MPFR_RNDN) && (bcp == 0) && (bbcp) && (bbcp1)) )
			{
				DEBUG( printf( "(Truncate) Do sub\n" ) );
				mpn_sub_1( ap, ap, n, MPFR_LIMB_ONE << sh );
				mpn_lshift( ap, ap, n, 1 );
				ap[0] |= MPFR_LIMB_ONE << sh;
				bx--;
				/* FIXME: Explain why it works (or why not)... */
				inexact = (bcp1 == 0) ? 0 : (rnd_mode == MPFR_RNDN) ? -1 : 1;
				goto end_of_sub;
			}
		}
	}

	/* Calcul of Inexact flag.*/
	inexact = MPFR_LIKELY( bcp || bcp1 ) ? 1 : 0;

end_of_sub:
	/* Update Expo */
	/* FIXME: Is this test really useful?
		If d==0      : Exact case. This is never called.
		if 1 < d < p : bx=MPFR_EXP(b) or MPFR_EXP(b)-1 > MPFR_EXP(c) > emin
		if d == 1    : bx=MPFR_EXP(b). If we could lose any bits, the exact
					   normalisation is called.
		if d >=  p   : bx=MPFR_EXP(b) >= MPFR_EXP(c) + p > emin
	   After SubOneUlp, we could have one bit less.
		if 1 < d < p : bx >= MPFR_EXP(b)-2 >= MPFR_EXP(c) > emin
		if d == 1    : bx >= MPFR_EXP(b)-1 = MPFR_EXP(c) > emin.
		if d >=  p   : bx >= MPFR_EXP(b)-1 > emin since p>=2.
	*/
	MPFR_ASSERTD( bx >= __gmpfr_emin );
	/*
	  if (MPFR_UNLIKELY(bx < __gmpfr_emin))
	  {
		DEBUG( printf("(Final Underflow)\n") );
		if (rnd_mode == MPFR_RNDN &&
			(bx < __gmpfr_emin - 1 ||
			 (inexact >= 0 && mpfr_powerof2_raw (a))))
		  rnd_mode = MPFR_RNDZ;
		MPFR_TMP_FREE(marker);
		return mpfr_underflow (a, rnd_mode, MPFR_SIGN(a));
	  }
	*/
	MPFR_SET_EXP( a, bx );

	MPFR_TMP_FREE( marker );
	MPFR_RET( inexact * MPFR_INT_SIGN( a ) );
}


int
atn_sub1( mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_rnd_t rnd_mode )
{
	int sign;
	mpfr_uexp_t diff_exp;
	mpfr_prec_t cancel, cancel1;
	mp_size_t cancel2, an, bn, cn, cn0;
	mp_limb_t* ap, * bp, * cp;
	mp_limb_t carry, bb, cc;
	int inexact, shift_b, shift_c, add_exp = 0;
	int cmp_low = 0; /* used for rounding to nearest: 0 if low(b) = low(c),
						negative if low(b) < low(c), positive if low(b)>low(c) */
	int sh, k;
	MPFR_TMP_DECL( marker );

	MPFR_TMP_MARK( marker );
	ap = MPFR_MANT( a );
	an = MPFR_LIMB_SIZE( a );

	sign = mpfr_cmp2( b, c, &cancel );
	if ( MPFR_UNLIKELY( sign == 0 ) )
	{
		if ( rnd_mode == MPFR_RNDD )
			MPFR_SET_NEG( a );
		else
			MPFR_SET_POS( a );
		MPFR_SET_ZERO( a );
		MPFR_RET( 0 );
	}

	/*
	 * If subtraction: sign(a) = sign * sign(b)
	 * If addition: sign(a) = sign of the larger argument in absolute value.
	 *
	 * Both cases can be simplidied in:
	 * if (sign>0)
	 *    if addition: sign(a) = sign * sign(b) = sign(b)
	 *    if subtraction, b is greater, so sign(a) = sign(b)
	 * else
	 *    if subtraction, sign(a) = - sign(b)
	 *    if addition, sign(a) = sign(c) (since c is greater)
	 *      But if it is an addition, sign(b) and sign(c) are opposed!
	 *      So sign(a) = - sign(b)
	 */

	if ( sign < 0 ) /* swap b and c so that |b| > |c| */
	{
		mpfr_srcptr t;
		MPFR_SET_OPPOSITE_SIGN( a, b );
		t = b; b = c; c = t;
	}
	else
		MPFR_SET_SAME_SIGN( a, b );

	/* Check if c is too small.
	   A more precise test is to replace 2 by
		(rnd == MPFR_RNDN) + mpfr_power2_raw (b)
		but it is more expensive and not very useful */
	if ( MPFR_UNLIKELY( MPFR_GET_EXP( c ) <= MPFR_GET_EXP( b )
		- (mpfr_exp_t)MAX( MPFR_PREC( a ), MPFR_PREC( b ) ) - 2 ) )
	{
		/* Remember, we can't have an exact result! */
		/*   A.AAAAAAAAAAAAAAAAA
		   = B.BBBBBBBBBBBBBBB
			-                     C.CCCCCCCCCCCCC */
			/* A = S*ABS(B) +/- ulp(a) */
		MPFR_SET_EXP( a, MPFR_GET_EXP( b ) );
		MPFR_RNDRAW_EVEN( inexact, a, MPFR_MANT( b ), MPFR_PREC( b ),
			rnd_mode, MPFR_SIGN( a ),
			if ( MPFR_UNLIKELY( ++MPFR_EXP( a ) > __gmpfr_emax ) )
				inexact = mpfr_overflow( a, rnd_mode, MPFR_SIGN( a ) ) );
		/* inexact = mpfr_set4 (a, b, rnd_mode, MPFR_SIGN (a));  */
		if ( inexact == 0 )
		{
			/* a = b (Exact)
			   But we know it isn't (Since we have to remove `c')
			   So if we round to Zero, we have to remove one ulp.
			   Otherwise the result is correctly rounded. */
			if ( MPFR_IS_LIKE_RNDZ( rnd_mode, MPFR_IS_NEG( a ) ) )
			{
				mpfr_nexttozero( a );
				MPFR_RET( -MPFR_INT_SIGN( a ) );
			}
			MPFR_RET( MPFR_INT_SIGN( a ) );
		}
		else
		{
			/*   A.AAAAAAAAAAAAAA
			   = B.BBBBBBBBBBBBBBB
				-                   C.CCCCCCCCCCCCC */
				/* It isn't exact so Prec(b) > Prec(a) and the last
				   Prec(b)-Prec(a) bits of `b' are not zeros.
				   Which means that removing c from b can't generate a carry
				   execpt in case of even rounding.
				   In all other case the result and the inexact flag should be
				   correct (We can't have an exact result).
				   In case of EVEN rounding:
					 1.BBBBBBBBBBBBBx10
				   -                     1.CCCCCCCCCCCC
				   = 1.BBBBBBBBBBBBBx01  Rounded to Prec(b)
				   = 1.BBBBBBBBBBBBBx    Nearest / Rounded to Prec(a)
				   Set gives:
					 1.BBBBBBBBBBBBB0   if inexact == EVEN_INEX  (x == 0)
					 1.BBBBBBBBBBBBB1+1 if inexact == -EVEN_INEX (x == 1)
				   which means we get a wrong rounded result if x==1,
				   i.e. inexact= MPFR_EVEN_INEX */
			if ( MPFR_UNLIKELY( inexact == MPFR_EVEN_INEX * MPFR_INT_SIGN( a ) ) )
			{
				mpfr_nexttozero( a );
				inexact = -MPFR_INT_SIGN( a );
			}
			MPFR_RET( inexact );
		}
	}

	diff_exp = (mpfr_uexp_t)MPFR_GET_EXP( b ) - MPFR_GET_EXP( c );

	/* reserve a space to store b aligned with the result, i.e. shifted by
	   (-cancel) % GMP_NUMB_BITS to the right */
	bn = MPFR_LIMB_SIZE( b );
	MPFR_UNSIGNED_MINUS_MODULO( shift_b, cancel );
	cancel1 = (cancel + shift_b) / GMP_NUMB_BITS;

	/* the high cancel1 limbs from b should not be taken into account */
	if ( MPFR_UNLIKELY( shift_b == 0 ) )
	{
		bp = MPFR_MANT( b ); /* no need of an extra space */
		/* Ensure ap != bp */
		if ( MPFR_UNLIKELY( ap == bp ) )
		{
			bp = MPFR_TMP_LIMBS_ALLOC( bn );
			MPN_COPY( bp, ap, bn );
		}
	}
	else
	{
		bp = MPFR_TMP_LIMBS_ALLOC( bn + 1 );
		bp[0] = mpn_rshift( bp + 1, MPFR_MANT( b ), bn++, shift_b );
	}

	/* reserve a space to store c aligned with the result, i.e. shifted by
		(diff_exp-cancel) % GMP_NUMB_BITS to the right */
	cn = MPFR_LIMB_SIZE( c );
	if ( (UINT_MAX % GMP_NUMB_BITS) == (GMP_NUMB_BITS - 1)
		&& ((-(unsigned)1) % GMP_NUMB_BITS > 0) )
		shift_c = ((mpfr_uexp_t)diff_exp - cancel) % GMP_NUMB_BITS;
	else
	{
		shift_c = diff_exp - (cancel % GMP_NUMB_BITS);
		shift_c = (shift_c + GMP_NUMB_BITS) % GMP_NUMB_BITS;
	}
	MPFR_ASSERTD( shift_c >= 0 && shift_c < GMP_NUMB_BITS );

	if ( MPFR_UNLIKELY( shift_c == 0 ) )
	{
		cp = MPFR_MANT( c );
		/* Ensure ap != cp */
		if ( ap == cp )
		{
			cp = MPFR_TMP_LIMBS_ALLOC( cn );
			MPN_COPY( cp, ap, cn );
		}
	}
	else
	{
		cp = MPFR_TMP_LIMBS_ALLOC( cn + 1 );
		cp[0] = mpn_rshift( cp + 1, MPFR_MANT( c ), cn++, shift_c );
	}

#ifdef DEBUG
	printf( "rnd=%s shift_b=%d shift_c=%d diffexp=%lu\n",
		mpfr_print_rnd_mode( rnd_mode ), shift_b, shift_c,
		(unsigned long)diff_exp );
#endif

	MPFR_ASSERTD( ap != cp );
	MPFR_ASSERTD( bp != cp );

	/* here we have shift_c = (diff_exp - cancel) % GMP_NUMB_BITS,
		  0 <= shift_c < GMP_NUMB_BITS
	   thus we want cancel2 = ceil((cancel - diff_exp) / GMP_NUMB_BITS) */

	   /* Possible optimization with a C99 compiler (i.e. well-defined
		  integer division): if MPFR_PREC_MAX is reduced to
		  ((mpfr_prec_t)((mpfr_uprec_t)(~(mpfr_uprec_t)0)>>1) - GMP_NUMB_BITS + 1)
		  and diff_exp is of type mpfr_exp_t (no need for mpfr_uexp_t, since
		  the sum or difference of 2 exponents must be representable, as used
		  by the multiplication code), then the computation of cancel2 could
		  be simplified to
			cancel2 = (cancel - (diff_exp - shift_c)) / GMP_NUMB_BITS;
		  because cancel, diff_exp and shift_c are all non-negative and
		  these variables are signed. */

	MPFR_ASSERTD( cancel >= 0 );
	if ( cancel >= diff_exp )
		/* Note that cancel is signed and will be converted to mpfr_uexp_t
		   (type of diff_exp) in the expression below, so that this will
		   work even if cancel is very large and diff_exp = 0. */
		cancel2 = (cancel - diff_exp + (GMP_NUMB_BITS - 1)) / GMP_NUMB_BITS;
	else
		cancel2 = -(mp_size_t)((diff_exp - cancel) / GMP_NUMB_BITS);
	/* the high cancel2 limbs from b should not be taken into account */
#ifdef DEBUG
	printf( "cancel=%lu cancel1=%lu cancel2=%ld\n",
		(unsigned long)cancel, (unsigned long)cancel1, (long)cancel2 );
#endif

	/*               ap[an-1]        ap[0]
			   <----------------+-----------|---->
			   <----------PREC(a)----------><-sh->
   cancel1
   limbs        bp[bn-cancel1-1]
   <--...-----><----------------+-----------+----------->
	cancel2
	limbs       cp[cn-cancel2-1]                                    cancel2 >= 0
	  <--...--><----------------+----------------+---------------->
				  (-cancel2)                                        cancel2 < 0
					 limbs      <----------------+---------------->
	*/

	/* first part: put in ap[0..an-1] the value of high(b) - high(c),
	   where high(b) consists of the high an+cancel1 limbs of b,
	   and high(c) consists of the high an+cancel2 limbs of c.
	 */

	 /* copy high(b) into a */
	if ( MPFR_LIKELY( an + (mp_size_t)cancel1 <= bn ) )
		/* a: <----------------+-----------|---->
		   b: <-----------------------------------------> */
		MPN_COPY( ap, bp + bn - (an + cancel1), an );
	else
		/* a: <----------------+-----------|---->
		   b: <-------------------------> */
		if ( (mp_size_t)cancel1 < bn ) /* otherwise b does not overlap with a */
		{
			MPN_ZERO( ap, an + cancel1 - bn );
			MPN_COPY( ap + (an + cancel1 - bn), bp, bn - cancel1 );
		}
		else
			MPN_ZERO( ap, an );

	/* subtract high(c) */
	if ( MPFR_LIKELY( an + cancel2 > 0 ) ) /* otherwise c does not overlap with a */
	{
		mp_limb_t* ap2;

		if ( cancel2 >= 0 )
		{
			if ( an + cancel2 <= cn )
				/* a: <----------------------------->
				   c: <-----------------------------------------> */
				atn_sub_n( ap, ap, cp + cn - (an + cancel2), an );
			else
				/* a: <---------------------------->
				   c: <-------------------------> */
			{
				ap2 = ap + an + (cancel2 - cn);
				if ( cn > cancel2 )
					atn_sub_n( ap2, ap2, cp, cn - cancel2 );
			}
		}
		else /* cancel2 < 0 */
		{
			mp_limb_t borrow;

			if ( an + cancel2 <= cn )
				/* a: <----------------------------->
				   c: <-----------------------------> */
				borrow = atn_sub_n( ap, ap, cp + cn - (an + cancel2),
					an + cancel2 );
			else
				/* a: <---------------------------->
				   c: <----------------> */
			{
				ap2 = ap + an + cancel2 - cn;
				borrow = atn_sub_n( ap2, ap2, cp, cn );
			}
			ap2 = ap + an + cancel2;
			mpn_sub_1( ap2, ap2, -cancel2, borrow );
		}
	}

	/* now perform rounding */
	sh = (mpfr_prec_t)an * GMP_NUMB_BITS - MPFR_PREC( a );
	/* last unused bits from a */
	carry = ap[0] & MPFR_LIMB_MASK( sh );
	ap[0] -= carry;

	if ( MPFR_LIKELY( rnd_mode == MPFR_RNDN ) )
	{
		if ( MPFR_LIKELY( sh ) )
		{
			/* can decide except when carry = 2^(sh-1) [middle]
			   or carry = 0 [truncate, but cannot decide inexact flag] */
			if ( carry > (MPFR_LIMB_ONE << (sh - 1)) )
				goto add_one_ulp;
			else if ( (0 < carry) && (carry < (MPFR_LIMB_ONE << (sh - 1))) )
			{
				inexact = -1; /* result if smaller than exact value */
				goto truncate;
			}
			/* now carry = 2^(sh-1), in which case cmp_low=2,
			   or carry = 0, in which case cmp_low=0 */
			cmp_low = (carry == 0) ? 0 : 2;
		}
	}
	else /* directed rounding: set rnd_mode to RNDZ iff toward zero */
	{
		if ( MPFR_IS_RNDUTEST_OR_RNDDNOTTEST( rnd_mode, MPFR_IS_NEG( a ) ) )
			rnd_mode = MPFR_RNDZ;

		if ( carry )
		{
			if ( rnd_mode == MPFR_RNDZ )
			{
				inexact = -1;
				goto truncate;
			}
			else /* round away */
				goto add_one_ulp;
		}
	}

	/* we have to consider the low (bn - (an+cancel1)) limbs from b,
	   and the (cn - (an+cancel2)) limbs from c. */
	bn -= an + cancel1;
	cn0 = cn;
	cn -= an + cancel2;

#ifdef DEBUG
	printf( "last sh=%d bits from a are %lu, bn=%ld, cn=%ld\n",
		sh, (unsigned long)carry, (long)bn, (long)cn );
#endif

	/* for rounding to nearest, we couldn't conclude up to here in the following
	   cases:
	   1. sh = 0, then cmp_low=0: we can either truncate, subtract one ulp
		  or add one ulp: -1 ulp < low(b)-low(c) < 1 ulp
	   2. sh > 0 but the low sh bits from high(b)-high(c) equal 2^(sh-1):
		  -0.5 ulp <= -1/2^sh < low(b)-low(c)-0.5 < 1/2^sh <= 0.5 ulp
		  we can't decide the rounding, in that case cmp_low=2:
		  either we truncate and flag=-1, or we add one ulp and flag=1
	   3. the low sh>0 bits from high(b)-high(c) equal 0: we know we have to
		  truncate but we can't decide the ternary value, here cmp_low=0:
		  -0.5 ulp <= -1/2^sh < low(b)-low(c) < 1/2^sh <= 0.5 ulp
		  we always truncate and inexact can be any of -1,0,1
	*/

	/* note: here cn might exceed cn0, in which case we consider a zero limb */
	for ( k = 0; (bn > 0) || (cn > 0); k = 1 )
	{
		/* if cmp_low < 0, we know low(b) - low(c) < 0
		   if cmp_low > 0, we know low(b) - low(c) > 0
			  (more precisely if cmp_low = 2, low(b) - low(c) = 0.5 ulp so far)
		   if cmp_low = 0, so far low(b) - low(c) = 0 */

		   /* get next limbs */
		bb = (bn > 0) ? bp[--bn] : 0;
		if ( (cn > 0) && (cn-- <= cn0) )
			cc = cp[cn];
		else
			cc = 0;

		/* cmp_low compares low(b) and low(c) */
		if ( cmp_low == 0 ) /* case 1 or 3 */
			cmp_low = (bb < cc) ? -2 + k : (bb > cc) ? 1 : 0;

		/* Case 1 for k=0 splits into 7 subcases:
		   1a: bb > cc + half
		   1b: bb = cc + half
		   1c: 0 < bb - cc < half
		   1d: bb = cc
		   1e: -half < bb - cc < 0
		   1f: bb - cc = -half
		   1g: bb - cc < -half

		   Case 2 splits into 3 subcases:
		   2a: bb > cc
		   2b: bb = cc
		   2c: bb < cc

		   Case 3 splits into 3 subcases:
		   3a: bb > cc
		   3b: bb = cc
		   3c: bb < cc
		*/

		/* the case rounding to nearest with sh=0 is special since one couldn't
		   subtract above 1/2 ulp in the trailing limb of the result */
		if ( rnd_mode == MPFR_RNDN && sh == 0 && k == 0 ) /* case 1 for k=0 */
		{
			mp_limb_t half = MPFR_LIMB_HIGHBIT;

			/* add one ulp if bb > cc + half
			   truncate if cc - half < bb < cc + half
			   sub one ulp if bb < cc - half
			*/

			if ( cmp_low < 0 ) /* bb < cc: -1 ulp < low(b) - low(c) < 0,
								cases 1e, 1f and 1g */
			{
				if ( cc >= half )
					cc -= half;
				else /* since bb < cc < half, bb+half < 2*half */
					bb += half;
				/* now we have bb < cc + half:
				   we have to subtract one ulp if bb < cc,
				   and truncate if bb > cc */
			}
			else if ( cmp_low >= 0 ) /* bb >= cc, cases 1a to 1d */
			{
				if ( cc < half )
					cc += half;
				else /* since bb >= cc >= half, bb - half >= 0 */
					bb -= half;
				/* now we have bb > cc - half: we have to add one ulp if bb > cc,
				   and truncate if bb < cc */
				if ( cmp_low > 0 )
					cmp_low = 2;
			}
		}

#ifdef DEBUG
		printf( "k=%u bb=%lu cc=%lu cmp_low=%d\n", k,
			(unsigned long)bb, (unsigned long)cc, cmp_low );
#endif
		if ( cmp_low < 0 ) /* low(b) - low(c) < 0: either truncate or subtract
							one ulp */
		{
			if ( rnd_mode == MPFR_RNDZ )
				goto sub_one_ulp; /* set inexact=-1 */
			else if ( rnd_mode != MPFR_RNDN ) /* round away */
			{
				inexact = 1;
				goto truncate;
			}
			else /* round to nearest */
			{
				/* If cmp_low < 0 and bb > cc, then -0.5 ulp < low(b)-low(c) < 0,
				   whatever the value of sh.
				   If sh>0, then cmp_low < 0 implies that the initial neglected
				   sh bits were 0 (otherwise cmp_low=2 initially), thus the
				   weight of the new bits is less than 0.5 ulp too.
				   If k > 0 (and sh=0) this means that either the first neglected
				   limbs bb and cc were equal (thus cmp_low was 0 for k=0),
				   or we had bb - cc = -0.5 ulp or 0.5 ulp.
				   The last case is not possible here since we would have
				   cmp_low > 0 which is sticky.
				   In the first case (where we have cmp_low = -1), we truncate,
				   whereas in the 2nd case we have cmp_low = -2 and we subtract
				   one ulp.
				*/
				if ( bb > cc || sh > 0 || cmp_low == -1 )
				{  /* -0.5 ulp < low(b)-low(c) < 0,
					  bb > cc corresponds to cases 1e and 1f1
					  sh > 0 corresponds to cases 3c and 3b3
					  cmp_low = -1 corresponds to case 1d3 (also 3b3) */
					inexact = 1;
					goto truncate;
				}
				else if ( bb < cc ) /* here sh = 0 and low(b)-low(c) < -0.5 ulp,
									 this corresponds to cases 1g and 1f3 */
					goto sub_one_ulp;
				/* the only case where we can't conclude is sh=0 and bb=cc,
				   i.e., we have low(b) - low(c) = -0.5 ulp (up to now), thus
				   we don't know if we must truncate or subtract one ulp.
				   Note: for sh=0 we can't have low(b) - low(c) = -0.5 ulp up to
				   now, since low(b) - low(c) > 1/2^sh */
			}
		}
		else if ( cmp_low > 0 ) /* 0 < low(b) - low(c): either truncate or
								 add one ulp */
		{
			if ( rnd_mode == MPFR_RNDZ )
			{
				inexact = -1;
				goto truncate;
			}
			else if ( rnd_mode != MPFR_RNDN ) /* round away */
				goto add_one_ulp;
			else /* round to nearest */
			{
				if ( bb > cc )
				{
					/* if sh=0, then bb>cc means that low(b)-low(c) > 0.5 ulp,
					   and similarly when cmp_low=2 */
					if ( cmp_low == 2 ) /* cases 1a, 1b1, 2a and 2b1 */
						goto add_one_ulp;
					/* sh > 0 and cmp_low > 0: this implies that the sh initial
					   neglected bits were 0, and the remaining low(b)-low(c)>0,
					   but its weight is less than 0.5 ulp */
					else /* 0 < low(b) - low(c) < 0.5 ulp, this corresponds to
							cases 3a, 1d1 and 3b1 */
					{
						inexact = -1;
						goto truncate;
					}
				}
				else if ( bb < cc ) /* 0 < low(b) - low(c) < 0.5 ulp, cases 1c,
									 1b3, 2b3 and 2c */
				{
					inexact = -1;
					goto truncate;
				}
				/* the only case where we can't conclude is bb=cc, i.e.,
				   low(b) - low(c) = 0.5 ulp (up to now), thus we don't know
				   if we must truncate or add one ulp. */
			}
		}
		/* after k=0, we cannot conclude in the following cases, we split them
		   according to the values of bb and cc for k=1:
		   1b. sh=0 and cmp_low = 1 and bb-cc = half [around 0.5 ulp]
			   1b1. bb > cc: add one ulp, inex = 1
			   1b2: bb = cc: cannot conclude
			   1b3: bb < cc: truncate, inex = -1
		   1d. sh=0 and cmp_low = 0 and bb-cc = 0 [around 0]
			   1d1: bb > cc: truncate, inex = -1
			   1d2: bb = cc: cannot conclude
			   1d3: bb < cc: truncate, inex = +1
		   1f. sh=0 and cmp_low = -1 and bb-cc = -half [around -0.5 ulp]
			   1f1: bb > cc: truncate, inex = +1
			   1f2: bb = cc: cannot conclude
			   1f3: bb < cc: sub one ulp, inex = -1
		   2b. sh > 0 and cmp_low = 2 and bb=cc [around 0.5 ulp]
			   2b1. bb > cc: add one ulp, inex = 1
			   2b2: bb = cc: cannot conclude
			   2b3: bb < cc: truncate, inex = -1
		   3b. sh > 0 and cmp_low = 0 [around 0]
			   3b1. bb > cc: truncate, inex = -1
			   3b2: bb = cc: cannot conclude
			   3b3: bb < cc: truncate, inex = +1
		*/
	}

	if ( (rnd_mode == MPFR_RNDN) && cmp_low != 0 )
	{
		/* even rounding rule */
		if ( (ap[0] >> sh) & 1 )
		{
			if ( cmp_low < 0 )
				goto sub_one_ulp;
			else
				goto add_one_ulp;
		}
		else
			inexact = (cmp_low > 0) ? -1 : 1;
	}
	else
		inexact = 0;
	goto truncate;

sub_one_ulp: /* sub one unit in last place to a */
	mpn_sub_1( ap, ap, an, MPFR_LIMB_ONE << sh );
	inexact = -1;
	goto end_of_sub;

add_one_ulp: /* add one unit in last place to a */
	if ( MPFR_UNLIKELY( mpn_add_1( ap, ap, an, MPFR_LIMB_ONE << sh ) ) )
		/* result is a power of 2: 11111111111111 + 1 = 1000000000000000 */
	{
		ap[an - 1] = MPFR_LIMB_HIGHBIT;
		add_exp = 1;
	}
	inexact = 1; /* result larger than exact value */

truncate:
	if ( MPFR_UNLIKELY( (ap[an - 1] >> (GMP_NUMB_BITS - 1)) == 0 ) )
		/* case 1 - epsilon */
	{
		ap[an - 1] = MPFR_LIMB_HIGHBIT;
		add_exp = 1;
	}

end_of_sub:
	/* we have to set MPFR_EXP(a) to MPFR_EXP(b) - cancel + add_exp, taking
	   care of underflows/overflows in that computation, and of the allowed
	   exponent range */
	if ( MPFR_LIKELY( cancel ) )
	{
		mpfr_exp_t exp_a;

		cancel -= add_exp; /* OK: add_exp is an int equal to 0 or 1 */
		exp_a = MPFR_GET_EXP( b ) - cancel;
		if ( MPFR_UNLIKELY( exp_a < __gmpfr_emin ) )
		{
			MPFR_TMP_FREE( marker );
			if ( rnd_mode == MPFR_RNDN &&
				(exp_a < __gmpfr_emin - 1 ||
					(inexact >= 0 && mpfr_powerof2_raw( a ))) )
				rnd_mode = MPFR_RNDZ;
			return mpfr_underflow( a, rnd_mode, MPFR_SIGN( a ) );
		}
		MPFR_SET_EXP( a, exp_a );
	}
	else /* cancel = 0: MPFR_EXP(a) <- MPFR_EXP(b) + add_exp */
	{
		/* in case cancel = 0, add_exp can still be 1, in case b is just
		   below a power of two, c is very small, prec(a) < prec(b),
		   and rnd=away or nearest */
		mpfr_exp_t exp_b;

		exp_b = MPFR_GET_EXP( b );
		if ( MPFR_UNLIKELY( add_exp && exp_b == __gmpfr_emax ) )
		{
			MPFR_TMP_FREE( marker );
			return mpfr_overflow( a, rnd_mode, MPFR_SIGN( a ) );
		}
		MPFR_SET_EXP( a, exp_b + add_exp );
	}
	MPFR_TMP_FREE( marker );
	/* check that result is msb-normalized */
	MPFR_ASSERTD( ap[an - 1] > ~ap[an - 1] );
	MPFR_RET( inexact * MPFR_INT_SIGN( a ) );
}


int
atn_sub( mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_rnd_t rnd_mode )
{
	MPFR_LOG_FUNC
	( ("b[%Pu]=%.*Rg c[%Pu]=%.*Rg rnd=%d",
		mpfr_get_prec( b ), mpfr_log_prec, b,
		mpfr_get_prec( c ), mpfr_log_prec, c, rnd_mode),
		("a[%Pu]=%.*Rg", mpfr_get_prec( a ), mpfr_log_prec, a) );

	if ( MPFR_ARE_SINGULAR( b, c ) )
	{
		if ( MPFR_IS_NAN( b ) || MPFR_IS_NAN( c ) )
		{
			MPFR_SET_NAN( a );
			MPFR_RET_NAN;
		}
		else if ( MPFR_IS_INF( b ) )
		{
			if ( !MPFR_IS_INF( c ) || MPFR_SIGN( b ) != MPFR_SIGN( c ) )
			{
				MPFR_SET_INF( a );
				MPFR_SET_SAME_SIGN( a, b );
				MPFR_RET( 0 ); /* exact */
			}
			else
			{
				MPFR_SET_NAN( a ); /* Inf - Inf */
				MPFR_RET_NAN;
			}
		}
		else if ( MPFR_IS_INF( c ) )
		{
			MPFR_SET_INF( a );
			MPFR_SET_OPPOSITE_SIGN( a, c );
			MPFR_RET( 0 ); /* exact */
		}
		else if ( MPFR_IS_ZERO( b ) )
		{
			if ( MPFR_IS_ZERO( c ) )
			{
				int sign = rnd_mode != MPFR_RNDD
					? ((MPFR_IS_NEG( b ) && MPFR_IS_POS( c )) ? -1 : 1)
					: ((MPFR_IS_POS( b ) && MPFR_IS_NEG( c )) ? 1 : -1);
				MPFR_SET_SIGN( a, sign );
				MPFR_SET_ZERO( a );
				MPFR_RET( 0 ); /* 0 - 0 is exact */
			}
			else
				return mpfr_neg( a, c, rnd_mode );
		}
		else
		{
			MPFR_ASSERTD( MPFR_IS_ZERO( c ) );
			return mpfr_set( a, b, rnd_mode );
		}
	}

	MPFR_ASSERTD( MPFR_IS_PURE_FP( b ) );
	MPFR_ASSERTD( MPFR_IS_PURE_FP( c ) );

	if ( MPFR_LIKELY( MPFR_SIGN( b ) == MPFR_SIGN( c ) ) )
	{ /* signs are equal, it's a real subtraction */
		if ( MPFR_LIKELY( MPFR_PREC( a ) == MPFR_PREC( b )
			&& MPFR_PREC( b ) == MPFR_PREC( c ) ) )
			return atn_sub1sp( a, b, c, rnd_mode );
		else
			return atn_sub1( a, b, c, rnd_mode );
	}
	else
	{ /* signs differ, it's an addition */
		if ( MPFR_GET_EXP( b ) < MPFR_GET_EXP( c ) )
		{ /* exchange rounding modes toward +/- infinity */
			int inexact;
			rnd_mode = MPFR_INVERT_RND( rnd_mode );
			if ( MPFR_LIKELY( MPFR_PREC( a ) == MPFR_PREC( b )
				&& MPFR_PREC( b ) == MPFR_PREC( c ) ) )
				inexact = atn_add1sp( a, c, b, rnd_mode );
			else
				inexact = atn_add1( a, c, b, rnd_mode );
			MPFR_CHANGE_SIGN( a );
			return -inexact;
		}
		else
		{
			if ( MPFR_LIKELY( MPFR_PREC( a ) == MPFR_PREC( b )
				&& MPFR_PREC( b ) == MPFR_PREC( c ) ) )
				return atn_add1sp( a, b, c, rnd_mode );
			else
				return atn_add1( a, b, c, rnd_mode );
		}
	}
}


namespace Athena
{
    void sub( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round )
    {
        atn_sub( a_Result.m_Value, a_Num1.m_Value, a_Num2.m_Value, a_Round );
    }
}