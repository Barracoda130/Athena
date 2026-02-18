#include "Number.hpp"

#include <cassert>
#include <iostream>

#include "mpfr-impl.h"


using namespace Athena;

namespace
{
    // TODO: Understand these macros
	// Macros from gmp-6.2.1/mpn/generic/add_n.c
#define MPN_OVERLAP_P(xp, xsize, yp, ysize)				\
  ((xp) + (xsize) > (yp) && (yp) + (ysize) > (xp))
#define MPN_SAME_OR_INCR2_P(dst, dsize, src, ssize)			\
  ((dst) <= (src) || ! MPN_OVERLAP_P (dst, dsize, src, ssize))
#define MPN_SAME_OR_INCR_P(dst, src, size)				\
  MPN_SAME_OR_INCR2_P(dst, size, src, size)

# define DEBUG(x)

    int atn_add1( mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, round_t rnd_mode );

    mp_limb_t atn_add_n( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {
        mp_limb_t ul, vl, sl, rl, cy, cy1, cy2;

        assert ( n >= 1 );
        assert ( MPN_SAME_OR_INCR_P( rp, up, n ) );
        assert ( MPN_SAME_OR_INCR_P( rp, vp, n ) );

        cy = 0;
        do
        {
            ul = *up++;
            vl = *vp++;
            sl = ul + vl;
            cy1 = sl < ul;
            rl = sl + cy;
            cy2 = rl < sl;
            cy = cy1 | cy2;
            *rp++ = rl;
            
        } while ( --n != 0 );

        return cy;
    }

    

    MPFR_HOT_FUNCTION_ATTR int
        atn_add1( mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, round_t rnd_mode )
    {
        mp_limb_t* aMant, * bMant, * cMant;
        precision_t aPrec, bPrec, cPrec, aMaxPrec;
        mp_size_t aNumLimbs, bNumLimbs, cNumLimbs;
        exponent_t difw, exp, diff_exp;
        int nonSigBitsRemaining, rb, fb, inex;
        // MPFR_TMP_DECL( marker );

		// Ensure operands b and c are unbounded floating-point numbers (UBF) and that the exponent of b is not less than that of c.
        // Only run in debug mode
        MPFR_ASSERTD( MPFR_IS_PURE_UBF( b ) );
        MPFR_ASSERTD( MPFR_IS_PURE_UBF( c ) );
        MPFR_ASSERTD( !MPFR_UBF_EXP_LESS_P( b, c ) );

		// Check if b is an unbounded floating-point number (UBF) and 
        // handle potential overflow if its exponent exceeds the maximum allowed exponent.
        if ( MPFR_UNLIKELY( MPFR_IS_UBF( b ) ) )
        {
			// If b is an unbounded floating-point number (UBF), retrieve its exponent using MPFR_UBF_GET_EXP and 
            // check if it exceeds the maximum allowed exponent (__gmpfr_emax). If it does, call mpfr_overflow to handle the overflow condition and return the appropriate value based on the sign of b.
            exp = MPFR_UBF_GET_EXP( b );
            if ( exp > __gmpfr_emax )
                return mpfr_overflow( a, rnd_mode, MPFR_SIGN( b ) );;
        }
        else
			// If b is not an unbounded floating-point number (UBF), retrieve its exponent using MPFR_GET_EXP
            exp = MPFR_GET_EXP( b );

		// In debug mode check that the exponent of b does not exceed the maximum allowed exponent (__gmpfr_emax).
        MPFR_ASSERTD( exp <= __gmpfr_emax );


        aPrec = MPFR_GET_PREC( a );
        bPrec = MPFR_GET_PREC( b );
        cPrec = MPFR_GET_PREC( c );

        aNumLimbs = MPFR_PREC2LIMBS( aPrec ); /* number of limbs of a */
        /**
         * Calculates the precision in bits of a multiprecision number by multiplying
         * the number of limbs (an) by the number of bits per limb (GMP_NUMB_BITS),
         * and stores the result in aq2 as an mpfr_prec_t type.
         */
         /**
          * @brief Calculates the precision in bits of operand 'a' multiplied by GMP_NUMB_BITS
          *
          * @details
          * aq2 represents the total number of bits available in the limb array of operand 'a'.
          * It is computed by multiplying the number of limbs (an) by GMP_NUMB_BITS (typically 32 or 64).
          *
          * @note
          * The difference between aq and aq2:
          * - aq: The actual precision (significant bits) of operand 'a' that are being used
          * - aq2: The total capacity (bits) of the limb array storing operand 'a'
          *
          * aq2 is typically >= aq since the allocated limbs may have more storage capacity
          * than the precision actually requires.
          */
        aMaxPrec = (mpfr_prec_t) aNumLimbs * GMP_NUMB_BITS;
        nonSigBitsRemaining = aMaxPrec - aPrec;                  /* non-significant bits in low limb */

        bNumLimbs = MPFR_PREC2LIMBS( bPrec ); /* number of limbs of b */
        cNumLimbs = MPFR_PREC2LIMBS( cPrec ); /* number of limbs of c */

        aMant = MPFR_MANT( a );
        bMant = MPFR_MANT( b );
        cMant = MPFR_MANT( c );

        // Optimisation, we ignore this for now
        /*
        if ( MPFR_UNLIKELY( aMant == bMant ) )
        {
            // If the two mantissas are equal, just copy b mantissa into a
            bMant = MPFR_TMP_LIMBS_ALLOC( bNumLimbs );
			MPN_COPY( bMant, aMant, bNumLimbs );    // Basically uses memcpy - copies bNumLimbs limbs from aMant to bMant

			// If the mantissa of a is the same as that of c, we need to update c's mantissa to point to 
            // the newly allocated memory for b's mantissa. This is necessary because we will be modifying a's mantissa 
            // (which is the same as b's), and we want to ensure that c's mantissa does not get 
            // inadvertently modified when we perform operations on a.
            if ( aMant == cMant )
            {
                cMant = bMant;
            }
        }
        else if ( aMant == cMant )
        {
            cMant = MPFR_TMP_LIMBS_ALLOC( cNumLimbs );
            MPN_COPY( cMant, aMant, cNumLimbs );
        }

        */

		// Set sign and rounding mode
        MPFR_SET_SAME_SIGN( a, b );
        MPFR_UPDATE2_RND_MODE( rnd_mode, MPFR_SIGN( b ) );

        /* now rnd_mode is either MPFR_RNDN, MPFR_RNDZ, MPFR_RNDA or MPFR_RNDF. */
        if ( MPFR_UNLIKELY( MPFR_IS_UBF( c ) ) )
        {
            MPFR_STAT_STATIC_ASSERT( MPFR_EXP_MAX > MPFR_PREC_MAX );
            diff_exp = mpfr_ubf_diff_exp( b, c );
        }
        else
            diff_exp = exp - MPFR_GET_EXP( c );

        MPFR_ASSERTD( diff_exp >= 0 );

        /*
         * 1. Compute the significant part A', the non-significant bits of A
         * are taken into account.
         *
         * 2. Perform the rounding. At each iteration, we remember:
         *     _ r = rounding bit
         *     _ f = following bits (same value)
         * where the result has the form: [number A]rfff...fff + a remaining
         * value in the interval [0,2) ulp. We consider the most significant
         * bits of the remaining value to update the result; a possible carry
         * is immediately taken into account and A is updated accordingly. As
         * soon as the bits f don't have the same value, A can be rounded.
         * Variables:
         *     _ rb = rounding bit (0 or 1).
         *     _ fb = following bits (0 or 1), then sticky bit.
         * If fb == 0, the only thing that can change is the sticky bit.
         */

        rb = fb = -1; /* means: not initialized */

        if ( MPFR_UNLIKELY( MPFR_UEXP( aMaxPrec ) <= diff_exp ) )
        { /* c does not overlap with a' */
            if ( MPFR_UNLIKELY( aNumLimbs > bNumLimbs ) )
            { /* a has more limbs than b */
              /* copy b to the most significant limbs of a */
                MPN_COPY( aMant + ( aNumLimbs - bNumLimbs ), bMant, bNumLimbs );
                /* zero the least significant limbs of a */
                MPN_ZERO( aMant, aNumLimbs - bNumLimbs );
            }
            else /* an <= bn */
            {
                /* copy the most significant limbs of b to a */
                MPN_COPY( aMant, bMant + ( bNumLimbs - aNumLimbs ), aNumLimbs );
            }
        }
        else /* aq2 > diff_exp */
        { /* c overlaps with a' */
            mp_limb_t* a2p;
            mp_limb_t cc;
            mpfr_prec_t dif;
            mp_size_t difn, k;
            int shift;

            /* copy c (shifted) into a */

            dif = aMaxPrec - diff_exp;
            /* dif is the number of bits of c which overlap with a' */

            difn = MPFR_PREC2LIMBS( dif );
            /* only the highest difn limbs from c have to be considered */
            if ( MPFR_UNLIKELY( difn > cNumLimbs ) )
            {
                /* c doesn't have enough limbs; take into account the virtual
                   zero limbs now by zeroing the least significant limbs of a' */
                MPFR_ASSERTD( difn - cNumLimbs <= aNumLimbs );
                MPN_ZERO( aMant, difn - cNumLimbs );
                difn = cNumLimbs;
            }
            k = diff_exp / GMP_NUMB_BITS;

            /* zero the most significant k limbs of a */
            a2p = aMant + ( aNumLimbs - k );
            MPN_ZERO( a2p, k );

            shift = diff_exp % GMP_NUMB_BITS;

            if ( MPFR_LIKELY( shift ) )
            {
                MPFR_ASSERTD( a2p - difn >= aMant );
                cc = mpn_rshift( a2p - difn, cMant + ( cNumLimbs - difn ), difn, shift );
                if ( MPFR_UNLIKELY( a2p - difn > aMant ) )
                    *( a2p - difn - 1 ) = cc;
            }
            else
                MPN_COPY( a2p - difn, cMant + ( cNumLimbs - difn ), difn );

            /* add b to a */
            cc = aNumLimbs > bNumLimbs
                ? atn_add_n( aMant + ( aNumLimbs - bNumLimbs ), aMant + ( aNumLimbs - bNumLimbs ), bMant, bNumLimbs )
                : atn_add_n( aMant, aMant, bMant + ( bNumLimbs - aNumLimbs ), aNumLimbs );

            if ( MPFR_UNLIKELY( cc ) ) /* carry */
            {
                if ( MPFR_UNLIKELY( exp == __gmpfr_emax ) )
                {
                    inex = mpfr_overflow( a, rnd_mode, MPFR_SIGN( a ) );
                    goto end_of_add;
                }
                exp++;
                rb = ( aMant[0] >> nonSigBitsRemaining ) & 1; /* LSB(a) --> rounding bit after the shift */
                if ( MPFR_LIKELY( nonSigBitsRemaining ) )
                {
                    mp_limb_t mask, bb;

                    mask = MPFR_LIMB_MASK( nonSigBitsRemaining );
                    bb = aMant[0] & mask;
                    aMant[0] &= MPFR_LIMB_LSHIFT( ~mask, 1 );
                    if ( bb == 0 )
                        fb = 0;
                    else if ( bb == mask )
                        fb = 1;
                }
                mpn_rshift( aMant, aMant, aNumLimbs, 1 );
                aMant[aNumLimbs - 1] += MPFR_LIMB_HIGHBIT;
                if ( nonSigBitsRemaining && fb < 0 )
                    goto rounding;
            } /* cc */
        } /* aq2 > diff_exp */

      /* zero the non-significant bits of a */
        if ( MPFR_LIKELY( rb < 0 && nonSigBitsRemaining ) )
        {
            mp_limb_t mask, bb;

            mask = MPFR_LIMB_MASK( nonSigBitsRemaining );
            bb = aMant[0] & mask;
            aMant[0] &= ~mask;
            rb = bb >> ( nonSigBitsRemaining - 1 );
            if ( MPFR_LIKELY( nonSigBitsRemaining > 1 ) )
            {
                mask >>= 1;
                bb &= mask;
                if ( bb == 0 )
                    fb = 0;
                else if ( bb == mask )
                    fb = 1;
                else
                    goto rounding;
            }
        }

        /* Determine rounding and sticky bits (and possible carry).
           In faithful rounding, we may stop two bits after ulp(a):
           the approximation is regarded as the number formed by a,
           the rounding bit rb and an additional bit fb; and the
           corresponding error is < 1/2 ulp of the unrounded result. */

        difw = (mpfr_exp_t) aNumLimbs - (mpfr_exp_t) ( diff_exp / GMP_NUMB_BITS );
        /* difw is the number of limbs from b (regarded as having an infinite
           precision) that have already been combined with c; -n if the next
           n limbs from b won't be combined with c. */

        if ( MPFR_UNLIKELY( bNumLimbs > aNumLimbs ) )
        { /* there are still limbs from b that haven't been taken into account */
            mp_size_t bk;

            if ( fb == 0 && difw <= 0 )
            {
                fb = 1; /* c hasn't been taken into account ==> sticky bit != 0 */
                goto rounding;
            }

            bk = bNumLimbs - aNumLimbs; /* index of lowest considered limb from b, > 0 */
            while ( difw < 0 )
            { /* ulp(next limb from b) > msb(c) */
                mp_limb_t bb;

                bb = bMant[--bk];

                MPFR_ASSERTD( fb != 0 );
                if ( fb > 0 )
                {
                    /* Note: Here, we can round to nearest, but the loop may still
                       be necessary to determine whether there is a carry from c,
                       which will have an effect on the ternary value. However, in
                       faithful rounding, we do not have to determine the ternary
                       value, so that we can end the loop here. */
                    if ( bb != MPFR_LIMB_MAX || rnd_mode == MPFR_RNDF )
                        goto rounding;
                }
                else /* fb not initialized yet */
                {
                    if ( rb < 0 ) /* rb not initialized yet */
                    {
                        rb = bb >> ( GMP_NUMB_BITS - 1 );
                        bb |= MPFR_LIMB_HIGHBIT;
                    }
                    fb = 1;
                    if ( bb != MPFR_LIMB_MAX )
                        goto rounding;
                }

                if ( bk == 0 )
                { /* b has entirely been read */
                    fb = 1; /* c hasn't been taken into account
                               ==> sticky bit != 0 */
                    goto rounding;
                }

                difw++;
            } /* while */
            MPFR_ASSERTD( bk > 0 && difw >= 0 );

            if ( difw <= cNumLimbs )
            {
                mp_size_t ck;
                mp_limb_t cprev;
                int difs;

                ck = cNumLimbs - difw;
                difs = diff_exp % GMP_NUMB_BITS;

                if ( difs == 0 && ck == 0 )
                    goto c_read;

                cprev = ck == cNumLimbs ? 0 : cMant[ck];

                if ( fb < 0 )
                {
                    mp_limb_t bb, cc;

                    if ( difs )
                    {
                        cc = cprev << ( GMP_NUMB_BITS - difs );
                        if ( --ck >= 0 )
                        {
                            cprev = cMant[ck];
                            cc += cprev >> difs;
                        }
                    }
                    else
                        cc = cMant[--ck];

                    bb = bMant[--bk] + cc;

                    if ( bb < cc /* carry */
                         && ( rb < 0 || ( rb ^= 1 ) == 0 )
                         && mpn_add_1( aMant, aMant, aNumLimbs, MPFR_LIMB_ONE << nonSigBitsRemaining ) )
                    {
                        if ( exp == __gmpfr_emax )
                        {
                            inex = mpfr_overflow( a, rnd_mode, MPFR_SIGN( a ) );
                            goto end_of_add;
                        }
                        exp++;
                        aMant[aNumLimbs - 1] = MPFR_LIMB_HIGHBIT;
                        rb = 0;
                    }

                    if ( rb < 0 ) /* rb not initialized yet */
                    {
                        rb = bb >> ( GMP_NUMB_BITS - 1 );
                        bb <<= 1;
                        bb |= bb >> ( GMP_NUMB_BITS - 1 );
                    }

                    fb = bb != 0;
                    if ( fb && bb != MPFR_LIMB_MAX )
                        goto rounding;
                } /* fb < 0 */

              /* At least two bits after ulp(a) have been read, which is
                 sufficient for faithful rounding, as we do not need to
                 determine on which side of a breakpoint the result is. */
                if ( rnd_mode == MPFR_RNDF )
                    goto rounding;

                while ( bk > 0 )
                {
                    mp_limb_t bb, cc;

                    if ( difs )
                    {
                        if ( ck < 0 )
                            goto c_read;
                        cc = cprev << ( GMP_NUMB_BITS - difs );
                        if ( --ck >= 0 )
                        {
                            cprev = cMant[ck];
                            cc += cprev >> difs;
                        }
                    }
                    else
                    {
                        if ( ck == 0 )
                            goto c_read;
                        cc = cMant[--ck];
                    }

                    bb = bMant[--bk] + cc;
                    if ( bb < cc ) /* carry */
                    {
                        fb ^= 1;
                        if ( fb )
                            goto rounding;
                        rb ^= 1;
                        if ( rb == 0 && mpn_add_1( aMant, aMant, aNumLimbs, MPFR_LIMB_ONE << nonSigBitsRemaining ) )
                        {
                            if ( MPFR_UNLIKELY( exp == __gmpfr_emax ) )
                            {
                                inex = mpfr_overflow( a, rnd_mode, MPFR_SIGN( a ) );
                                goto end_of_add;
                            }
                            exp++;
                            aMant[aNumLimbs - 1] = MPFR_LIMB_HIGHBIT;
                        }
                    } /* bb < cc */

                    if ( !fb && bb != 0 )
                    {
                        fb = 1;
                        goto rounding;
                    }
                    if ( fb && bb != MPFR_LIMB_MAX )
                        goto rounding;
                } /* while */

              /* b has entirely been read */

                if ( fb || ck < 0 )
                    goto rounding;
                if ( difs && MPFR_LIMB_LSHIFT( cprev, GMP_NUMB_BITS - difs ) != 0 )
                {
                    fb = 1;
                    goto rounding;
                }
                while ( ck )
                {
                    if ( cMant[--ck] )
                    {
                        fb = 1;
                        goto rounding;
                    }
                } /* while */
            } /* difw <= cn */
            else
            { /* c has entirely been read */
            c_read:
                if ( fb < 0 ) /* fb not initialized yet */
                {
                    mp_limb_t bb;

                    MPFR_ASSERTD( bk > 0 );
                    bb = bMant[--bk];
                    if ( rb < 0 ) /* rb not initialized yet */
                    {
                        rb = bb >> ( GMP_NUMB_BITS - 1 );
                        bb &= ~MPFR_LIMB_HIGHBIT;
                    }
                    fb = bb != 0;
                } /* fb < 0 */
                if ( fb || rnd_mode == MPFR_RNDF )
                    goto rounding;
                while ( bk )
                {
                    if ( bMant[--bk] )
                    {
                        fb = 1;
                        goto rounding;
                    }
                } /* while */
            } /* difw > cn */
        } /* bn > an */
        else if ( fb != 1 ) /* if fb == 1, the sticky bit is 1 (no possible carry) */
        { /* b has entirely been read */
            if ( difw > cNumLimbs )
            { /* c has entirely been read */
                if ( rb < 0 )
                    rb = 0;
                fb = 0;
            }
            else if ( diff_exp > MPFR_UEXP( aMaxPrec ) )
            { /* b is followed by at least a zero bit, then by c */
                if ( rb < 0 )
                    rb = 0;
                fb = 1;
            }
            else
            {
                mp_size_t ck;
                int difs;

                MPFR_ASSERTD( difw >= 0 && cNumLimbs >= difw );
                ck = cNumLimbs - difw;
                difs = diff_exp % GMP_NUMB_BITS;

                if ( difs == 0 && ck == 0 )
                { /* c has entirely been read */
                    if ( rb < 0 )
                        rb = 0;
                    fb = 0;
                }
                else
                {
                    mp_limb_t cc;

                    cc = difs ? ( MPFR_ASSERTD( ck < cNumLimbs ),
                                  cMant[ck] << ( GMP_NUMB_BITS - difs ) ) : cMant[--ck];
                    if ( rb < 0 )
                    {
                        rb = cc >> ( GMP_NUMB_BITS - 1 );
                        cc &= ~MPFR_LIMB_HIGHBIT;
                    }
                    if ( cc == 0 && rnd_mode == MPFR_RNDF )
                    {
                        fb = 0;
                        goto rounding;
                    }
                    while ( cc == 0 )
                    {
                        if ( ck == 0 )
                        {
                            fb = 0;
                            goto rounding;
                        }
                        cc = cMant[--ck];
                    } /* while */
                    fb = 1;
                }
            }
        } /* fb != 1 */

    rounding:
        /* rnd_mode should be one of MPFR_RNDN, MPFR_RNDF, MPFR_RNDZ or MPFR_RNDA */
        if ( MPFR_LIKELY( rnd_mode == MPFR_RNDN || rnd_mode == MPFR_RNDF ) )
        {
            if ( fb == 0 )
            {
                if ( rb == 0 )
                {
                    inex = 0;
                    goto set_exponent;
                }
                /* round to even */
                if ( aMant[0] & ( MPFR_LIMB_ONE << nonSigBitsRemaining ) )
                    goto rndn_away;
                else
                    goto rndn_zero;
            }
            if ( rb == 0 )
            {
            rndn_zero:
                inex = MPFR_IS_NEG( a ) ? 1 : -1;
                goto set_exponent;
            }
            else
            {
            rndn_away:
                inex = MPFR_IS_POS( a ) ? 1 : -1;
                goto add_one_ulp;
            }
        }
        else if ( rnd_mode == MPFR_RNDZ )
        {
            inex = rb || fb ? ( MPFR_IS_NEG( a ) ? 1 : -1 ) : 0;
            goto set_exponent;
        }
        else
        {
            MPFR_ASSERTN( rnd_mode == MPFR_RNDA );
            inex = rb || fb ? ( MPFR_IS_POS( a ) ? 1 : -1 ) : 0;
            if ( inex )
                goto add_one_ulp;
            else
                goto set_exponent;
        }

    add_one_ulp: /* add one unit in last place to a */
        if ( MPFR_UNLIKELY( mpn_add_1( aMant, aMant, aNumLimbs, MPFR_LIMB_ONE << nonSigBitsRemaining ) ) )
        {
            if ( MPFR_UNLIKELY( exp == __gmpfr_emax ) )
            {
                inex = mpfr_overflow( a, rnd_mode, MPFR_SIGN( a ) );
                goto end_of_add;
            }
            exp++;
            aMant[aNumLimbs - 1] = MPFR_LIMB_HIGHBIT;
        }

    set_exponent:
        if ( MPFR_UNLIKELY( exp < __gmpfr_emin ) )  /* possible if b and c are UBF's */
        {
            if ( rnd_mode == MPFR_RNDN &&
                 ( exp < __gmpfr_emin - 1 ||
                   ( inex >= 0 && mpfr_powerof2_raw( a ) ) ) )
                rnd_mode = MPFR_RNDZ;
            inex = mpfr_underflow( a, rnd_mode, MPFR_SIGN( a ) );
            goto end_of_add;
        }
        MPFR_SET_EXP( a, exp );

    end_of_add:
        MPFR_RET( inex );
    }

    int
        atn_add1sp( mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_rnd_t rnd_mode )
    {
        mpfr_uexp_t d;
        mpfr_prec_t p;
        unsigned int sh;
        mp_size_t n;
        mp_limb_t* ap, * cp;
        mpfr_exp_t bx;
        mp_limb_t limb;
        int inexact;

        MPFR_ASSERTD( MPFR_PREC( a ) == MPFR_PREC( b ) && MPFR_PREC( b ) == MPFR_PREC( c ) );
        MPFR_ASSERTD( MPFR_IS_PURE_FP( b ) );
        MPFR_ASSERTD( MPFR_IS_PURE_FP( c ) );
        MPFR_ASSERTD( MPFR_GET_EXP( b ) >= MPFR_GET_EXP( c ) );

        /* Read prec and num of limbs */
        p = MPFR_PREC( b );
        n = MPFR_PREC2LIMBS( p );
        MPFR_UNSIGNED_MINUS_MODULO( sh, p );
        bx = MPFR_GET_EXP( b );
        d = (mpfr_uexp_t)(bx - MPFR_GET_EXP( c ));

        DEBUG( printf( "New add1sp with diff=%lu\n", (unsigned long)d ) );

        if ( MPFR_UNLIKELY( d == 0 ) )
        {
            /* d==0 */
            DEBUG( mpfr_print_mant_binary( "C= ", MPFR_MANT( c ), p ) );
            DEBUG( mpfr_print_mant_binary( "B= ", MPFR_MANT( b ), p ) );
            bx++;                                /* exp + 1 */
            ap = MPFR_MANT( a );
            limb = mpn_add_n( ap, MPFR_MANT( b ), MPFR_MANT( c ), n );
            DEBUG( mpfr_print_mant_binary( "A= ", ap, p ) );
            MPFR_ASSERTD( limb != 0 );             /* There must be a carry */
            limb = ap[0];                        /* Get LSB (In fact, LSW) */
            mpn_rshift( ap, ap, n, 1 );            /* Shift mantissa A */
            ap[n - 1] |= MPFR_LIMB_HIGHBIT;        /* Set MSB */
            ap[0] &= ~MPFR_LIMB_MASK( sh );      /* Clear LSB bit */
            if ( MPFR_LIKELY( (limb & (MPFR_LIMB_ONE << sh)) == 0 ) ) /* Check exact case */
            {
                inexact = 0; goto set_exponent;
            }
            /* Zero: Truncate
               Nearest: Even Rule => truncate or add 1
               Away: Add 1 */
            if ( MPFR_LIKELY( rnd_mode == MPFR_RNDN ) )
            {
                if ( MPFR_LIKELY( (ap[0] & (MPFR_LIMB_ONE << sh)) == 0 ) )
                {
                    inexact = -1; goto set_exponent;
                }
                else
                    goto add_one_ulp;
            }
            MPFR_UPDATE_RND_MODE( rnd_mode, MPFR_IS_NEG( b ) );
            if ( rnd_mode == MPFR_RNDZ )
            {
                inexact = -1; goto set_exponent;
            }
            else
                goto add_one_ulp;
        }
        else if ( MPFR_UNLIKELY( d >= p ) )
        {
            if ( MPFR_LIKELY( d > p ) )
            {
                /* d > p : Copy B in A */
                /* Away:    Add 1
                   Nearest: Trunc
                   Zero:    Trunc */
                if ( MPFR_LIKELY( rnd_mode == MPFR_RNDN
                    || MPFR_IS_LIKE_RNDZ( rnd_mode, MPFR_IS_NEG( b ) ) ) )
                {
                copy_set_exponent:
                    ap = MPFR_MANT( a );
                    MPN_COPY( ap, MPFR_MANT( b ), n );
                    inexact = -1;
                    goto set_exponent;
                }
                else
                {
                copy_add_one_ulp:
                    ap = MPFR_MANT( a );
                    MPN_COPY( ap, MPFR_MANT( b ), n );
                    goto add_one_ulp;
                }
            }
            else
            {
                /* d==p : Copy B in A */
                /* Away:    Add 1
                   Nearest: Even Rule if C is a power of 2, else Add 1
                   Zero:    Trunc */
                if ( MPFR_LIKELY( rnd_mode == MPFR_RNDN ) )
                {
                    /* Check if C was a power of 2 */
                    cp = MPFR_MANT( c );
                    if ( MPFR_UNLIKELY( cp[n - 1] == MPFR_LIMB_HIGHBIT ) )
                    {
                        mp_size_t k = n - 1;
                        do
                        {
                            k--;
                        } while ( k >= 0 && cp[k] == 0 );
                        if ( MPFR_UNLIKELY( k < 0 ) )
                            /* Power of 2: Even rule */
                            if ( (MPFR_MANT( b )[0] & (MPFR_LIMB_ONE << sh)) == 0 )
                                goto copy_set_exponent;
                    }
                    /* Not a Power of 2 */
                    goto copy_add_one_ulp;
                }
                else if ( MPFR_IS_LIKE_RNDZ( rnd_mode, MPFR_IS_NEG( b ) ) )
                    goto copy_set_exponent;
                else
                    goto copy_add_one_ulp;
            }
        }
        else
        {
            mp_limb_t mask;
            mp_limb_t bcp, bcp1; /* Cp and C'p+1 */

            /* General case: 1 <= d < p */
            cp = MPFR_TMP_LIMBS_ALLOC( n );

            /* Shift c in temporary allocated place */
            {
                mpfr_uexp_t dm;
                mp_size_t m;

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
                    /* dm >=1 and m == 0: just shift */
                    MPFR_ASSERTD( dm >= 1 );
                    mpn_rshift( cp, MPFR_MANT( c ), n, dm );
                }
                else
                {
                    /* dm > 0 and m > 0: shift and zero  */
                    mpn_rshift( cp, MPFR_MANT( c ) + m, n - m, dm );
                    MPN_ZERO( cp + n - m, m );
                }
            }

            DEBUG( mpfr_print_mant_binary( "Before", MPFR_MANT( c ), p ) );
            DEBUG( mpfr_print_mant_binary( "B=    ", MPFR_MANT( b ), p ) );
            DEBUG( mpfr_print_mant_binary( "After ", cp, p ) );

            /* Compute bcp=Cp and bcp1=C'p+1 */
            if ( MPFR_LIKELY( sh > 0 ) )
            {
                /* Try to compute them from C' rather than C */
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
            else /* sh == 0 */
            {
                /* Compute Cp and C'p+1 from C with sh=0 */
                mp_limb_t* tp = MPFR_MANT( c );
                /* Start from bit x=p-d in mantissa C */
                mpfr_prec_t  x = p - d;
                mp_size_t   kx = n - 1 - (x / GMP_NUMB_BITS);
                mpfr_prec_t sx = GMP_NUMB_BITS - 1 - (x % GMP_NUMB_BITS);
                MPFR_ASSERTD( p >= d );
                bcp = tp[kx] & (MPFR_LIMB_ONE << sx);
                /* Looks at the last bits of limb kx (If sx=0, does nothing)*/
                if ( tp[kx] & MPFR_LIMB_MASK( sx ) )
                    bcp1 = 1;
                else
                {
                    do
                    {
                        kx--;
                    } while ( kx >= 0 && tp[kx] == 0 );
                    bcp1 = (kx >= 0);
                }
            }
            DEBUG( printf( "sh=%u Cp=%lu C'p+1=%lu\n", sh,
                (unsigned long)bcp, (unsigned long)bcp1 ) );

            /* Clean shifted C' */
            mask = ~MPFR_LIMB_MASK( sh );
            cp[0] &= mask;

            /* Add the mantissa c from b in a */
            ap = MPFR_MANT( a );
            limb = mpn_add_n( ap, MPFR_MANT( b ), cp, n );
            DEBUG( mpfr_print_mant_binary( "Add=  ", ap, p ) );

            /* Check for overflow */
            if ( MPFR_UNLIKELY( limb ) )
            {
                limb = ap[0] & (MPFR_LIMB_ONE << sh); /* Get LSB */
                mpn_rshift( ap, ap, n, 1 );          /* Shift mantissa*/
                bx++;                               /* Fix exponent */
                ap[n - 1] |= MPFR_LIMB_HIGHBIT;       /* Set MSB */
                ap[0] &= mask;                    /* Clear LSB bit */
                bcp1 |= bcp;                     /* Recompute C'p+1 */
                bcp = limb;                    /* Recompute Cp */
                DEBUG( printf( "(Overflow) Cp=%lu C'p+1=%lu\n",
                    (unsigned long)bcp, (unsigned long)bcp1 ) );
                DEBUG( mpfr_print_mant_binary( "Add=  ", ap, p ) );
            }

            /* Round:
                Zero: Truncate but could be exact.
                Away: Add 1 if Cp or C'p+1 !=0
                Nearest: Truncate but could be exact if Cp==0
                         Add 1 if C'p+1 !=0,
                         Even rule else */
            if ( MPFR_LIKELY( rnd_mode == MPFR_RNDN ) )
            {
                if ( MPFR_LIKELY( bcp == 0 ) )
                {
                    inexact = MPFR_LIKELY( bcp1 ) ? -1 : 0; goto set_exponent;
                }
                else if ( MPFR_UNLIKELY( bcp1 == 0 ) && (ap[0] & (MPFR_LIMB_ONE << sh)) == 0 )
                {
                    inexact = -1; goto set_exponent;
                }
                else
                    goto add_one_ulp;
            }
            MPFR_UPDATE_RND_MODE( rnd_mode, MPFR_IS_NEG( b ) );
            if ( rnd_mode == MPFR_RNDZ )
            {
                inexact = MPFR_LIKELY( bcp || bcp1 ) ? -1 : 0;
                goto set_exponent;
            }
            else
            {
                if ( MPFR_UNLIKELY( bcp == 0 && bcp1 == 0 ) )
                {
                    inexact = 0; goto set_exponent;
                }
                else
                    goto add_one_ulp;
            }
        }
        MPFR_ASSERTN( 0 );

    add_one_ulp:
        /* add one unit in last place to a */
        DEBUG( printf( "AddOneUlp\n" ) );
        if ( MPFR_UNLIKELY( mpn_add_1( ap, ap, n, MPFR_LIMB_ONE << sh ) ) )
        {
            /* Case 100000x0 = 0x1111x1 + 1*/
            DEBUG( printf( "Pow of 2\n" ) );
            bx++;
            ap[n - 1] = MPFR_LIMB_HIGHBIT;
        }
        inexact = 1;

    set_exponent:
        if ( MPFR_UNLIKELY( bx > __gmpfr_emax ) ) /* Check for overflow */
        {
            DEBUG( printf( "Overflow\n" ) );
            MPFR_SET_SAME_SIGN( a, b );
            return mpfr_overflow( a, rnd_mode, MPFR_SIGN( a ) );
        }
        MPFR_SET_EXP( a, bx );
        MPFR_SET_SAME_SIGN( a, b );

        MPFR_RET( inexact * MPFR_INT_SIGN( a ) );
    }


}


namespace Athena
{
    MPFR_HOT_FUNCTION_ATTR int
        atn_add( mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_rnd_t rnd_mode )
    {
        MPFR_LOG_FUNC
        ( ( "b[%Pd]=%.*Rg c[%Pd]=%.*Rg rnd=%d",
            mpfr_get_prec( b ), mpfr_log_prec, b,
            mpfr_get_prec( c ), mpfr_log_prec, c, rnd_mode ),
          ( "a[%Pd]=%.*Rg", mpfr_get_prec( a ), mpfr_log_prec, a ) );

        if ( MPFR_ARE_SINGULAR_OR_UBF( b, c ) )
        {
            if ( MPFR_IS_NAN( b ) || MPFR_IS_NAN( c ) )
            {
                MPFR_SET_NAN( a );
                MPFR_RET_NAN;
            }
            /* neither b nor c is NaN here */
            else if ( MPFR_IS_INF( b ) )
            {
                if ( !MPFR_IS_INF( c ) || MPFR_SIGN( b ) == MPFR_SIGN( c ) )
                {
                    MPFR_SET_INF( a );
                    MPFR_SET_SAME_SIGN( a, b );
                    MPFR_RET( 0 ); /* exact */
                }
                else
                {
                    MPFR_SET_NAN( a );
                    MPFR_RET_NAN;
                }
            }
            else if ( MPFR_IS_INF( c ) )
            {
                MPFR_SET_INF( a );
                MPFR_SET_SAME_SIGN( a, c );
                MPFR_RET( 0 ); /* exact */
            }
            /* now both b and c are finite numbers */
            else if ( MPFR_IS_ZERO( b ) )
            {
                if ( MPFR_IS_ZERO( c ) )
                {
                    /* for round away, we take the same convention for 0 + 0
                       as for round to zero or to nearest: it always gives +0,
                       except (-0) + (-0) = -0. */
                    MPFR_SET_SIGN( a,
                                   ( rnd_mode != MPFR_RNDD ?
                                     ( MPFR_IS_NEG( b ) && MPFR_IS_NEG( c ) ?
                                       MPFR_SIGN_NEG : MPFR_SIGN_POS ) :
                                     ( MPFR_IS_POS( b ) && MPFR_IS_POS( c ) ?
                                       MPFR_SIGN_POS : MPFR_SIGN_NEG ) ) );
                    MPFR_SET_ZERO( a );
                    MPFR_RET( 0 ); /* 0 + 0 is exact */
                }
                return mpfr_set( a, c, rnd_mode );
            }
            else if ( MPFR_IS_ZERO( c ) )
            {
                return mpfr_set( a, b, rnd_mode );
            }
            else
            {
                MPFR_ASSERTD( MPFR_IS_PURE_UBF( b ) );
                MPFR_ASSERTD( MPFR_IS_PURE_UBF( c ) );
                /* mpfr_sub1sp and mpfr_add1sp are not intended to support UBF,
                   for which optimization is less important. */
                if ( MPFR_SIGN( b ) != MPFR_SIGN( c ) )
                    return mpfr_sub1( a, b, c, rnd_mode );
                else if ( MPFR_UBF_EXP_LESS_P( b, c ) )
                    return atn_add1( a, c, b, rnd_mode );
                else
                    return atn_add1( a, b, c, rnd_mode );
            }
        }

        MPFR_ASSERTD( MPFR_IS_PURE_FP( b ) );
        MPFR_ASSERTD( MPFR_IS_PURE_FP( c ) );

        if ( MPFR_UNLIKELY( MPFR_SIGN( b ) != MPFR_SIGN( c ) ) )
        { /* signs differ, it is a subtraction */
            if ( MPFR_LIKELY( MPFR_PREC( a ) == MPFR_PREC( b )
                              && MPFR_PREC( b ) == MPFR_PREC( c ) ) )
                return mpfr_sub1sp( a, b, c, rnd_mode );
            else
                return mpfr_sub1( a, b, c, rnd_mode );
        }
        else
        { /* signs are equal, it's an addition */
            if ( MPFR_LIKELY( MPFR_PREC( a ) == MPFR_PREC( b )
                              && MPFR_PREC( b ) == MPFR_PREC( c ) ) )
                return atn_add1sp( a, b, c, rnd_mode );
            else
                if ( MPFR_GET_EXP( b ) < MPFR_GET_EXP( c ) )
                    return atn_add1( a, c, b, rnd_mode );
                else
                    return atn_add1( a, b, c, rnd_mode );
        }
    }

    void add( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round )
    {
        atn_add( a_Result.m_Value, a_Num1.m_Value, a_Num2.m_Value, a_Round );
    }
} // namespace Athena