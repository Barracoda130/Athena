#include "Number.hpp"

#include <cassert>
#include <iostream>
#include <chrono>
#include <random>
#include <iomanip>
#include <algorithm>

#include "mpfr-impl.h"

#if defined(_MSC_VER)
#include <intrin.h>
#endif


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
    mp_limb_t atn_add_n_intrinsic( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );

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

    mp_limb_t atn_add_n_intrinsic( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {
        assert( n >= 1 );
        assert( MPN_SAME_OR_INCR_P( rp, up, n ) );
        assert( MPN_SAME_OR_INCR_P( rp, vp, n ) );

#if defined(_MSC_VER)
        unsigned char carry = 0;

#if GMP_NUMB_BITS == 64
        auto* r = reinterpret_cast<unsigned __int64*>(rp);
        auto* u = reinterpret_cast<const unsigned __int64*>(up);
        auto* v = reinterpret_cast<const unsigned __int64*>(vp);

        while ( n >= 4 )
        {
            carry = _addcarry_u64( carry, u[0], v[0], &r[0] );
            carry = _addcarry_u64( carry, u[1], v[1], &r[1] );
            carry = _addcarry_u64( carry, u[2], v[2], &r[2] );
            carry = _addcarry_u64( carry, u[3], v[3], &r[3] );
            u += 4; v += 4; r += 4; n -= 4;
        }
        while ( n-- > 0 )
            carry = _addcarry_u64( carry, *u++, *v++, r++ );

#elif GMP_NUMB_BITS == 32
        auto* r = reinterpret_cast<unsigned int*>(rp);
        auto* u = reinterpret_cast<const unsigned int*>(up);
        auto* v = reinterpret_cast<const unsigned int*>(vp);

        while ( n >= 4 )
        {
            carry = _addcarry_u32( carry, u[0], v[0], &r[0] );
            carry = _addcarry_u32( carry, u[1], v[1], &r[1] );
            carry = _addcarry_u32( carry, u[2], v[2], &r[2] );
            carry = _addcarry_u32( carry, u[3], v[3], &r[3] );
            u += 4; v += 4; r += 4; n -= 4;
        }
        while ( n-- > 0 )
            carry = _addcarry_u32( carry, *u++, *v++, r++ );

#else
        mp_limb_t ul, vl, sl, rl, cy, cy1, cy2;
        cy = 0;
        mp_size_t i = 0;
        do
        {
            ul = up[i];
            vl = vp[i];
            sl = ul + vl;
            cy1 = sl < ul;
            rl = sl + cy;
            cy2 = rl < sl;
            cy = cy1 | cy2;
            rp[i] = rl;
            i++;
        } while ( i < n );
        carry = static_cast<unsigned char>( cy );
#endif

        return static_cast<mp_limb_t>( carry );
#elif defined(__GNUC__) || defined(__clang__)
        mp_limb_t carry = 0;
        while ( n >= 4 )
        {
            mp_limb_t sum;
            unsigned char c1 = __builtin_add_overflow( up[0], vp[0], &sum );
            unsigned char c2 = __builtin_add_overflow( sum, carry, &rp[0] );
            carry = c1 | c2;

            c1 = __builtin_add_overflow( up[1], vp[1], &sum );
            c2 = __builtin_add_overflow( sum, carry, &rp[1] );
            carry = c1 | c2;

            c1 = __builtin_add_overflow( up[2], vp[2], &sum );
            c2 = __builtin_add_overflow( sum, carry, &rp[2] );
            carry = c1 | c2;

            c1 = __builtin_add_overflow( up[3], vp[3], &sum );
            c2 = __builtin_add_overflow( sum, carry, &rp[3] );
            carry = c1 | c2;

            up += 4; vp += 4; rp += 4; n -= 4;
        }
        while ( n-- > 0 )
        {
            mp_limb_t sum;
            unsigned char c1 = __builtin_add_overflow( *up++, *vp++, &sum );
            unsigned char c2 = __builtin_add_overflow( sum, carry, rp );
            carry = c1 | c2;
            ++rp;
        }
        return carry;
#else
        mp_limb_t ul, vl, sl, rl, cy, cy1, cy2;
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
#endif
    }

    

    int
        atn_add1( mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, round_t rnd_mode )
    {
        mp_limb_t* ap, * bp, * cp;
        mpfr_prec_t aq, bq, cq, aq2;
        mp_size_t an, bn, cn;
        mpfr_exp_t difw, exp;
        int sh, rb, fb, inex;
        mpfr_uexp_t diff_exp;
        MPFR_TMP_DECL( marker );

        MPFR_ASSERTD( MPFR_IS_PURE_FP( b ) );
        MPFR_ASSERTD( MPFR_IS_PURE_FP( c ) );

        MPFR_TMP_MARK( marker );

        aq = MPFR_PREC( a );
        bq = MPFR_PREC( b );
        cq = MPFR_PREC( c );

        an = MPFR_PREC2LIMBS( aq ); /* number of limbs of a */
        aq2 = (mpfr_prec_t)an * GMP_NUMB_BITS;
        sh = aq2 - aq;                  /* non-significant bits in low limb */

        bn = MPFR_PREC2LIMBS( bq ); /* number of limbs of b */
        cn = MPFR_PREC2LIMBS( cq ); /* number of limbs of c */

        ap = MPFR_MANT( a );
        bp = MPFR_MANT( b );
        cp = MPFR_MANT( c );

        if ( MPFR_UNLIKELY( ap == bp ) )
        {
            bp = MPFR_TMP_LIMBS_ALLOC( bn );
            MPN_COPY( bp, ap, bn );
            if ( ap == cp )
            {
                cp = bp;
            }
        }
        else if ( MPFR_UNLIKELY( ap == cp ) )
        {
            cp = MPFR_TMP_LIMBS_ALLOC( cn );
            MPN_COPY( cp, ap, cn );
        }

        exp = MPFR_GET_EXP( b );
        MPFR_SET_SAME_SIGN( a, b );
        MPFR_UPDATE2_RND_MODE( rnd_mode, MPFR_SIGN( b ) );
        /* now rnd_mode is either MPFR_RNDN, MPFR_RNDZ or MPFR_RNDA */
        /* Note: exponents can be negative, but the unsigned subtraction is
           a modular subtraction, so that one gets the correct result. */
        diff_exp = (mpfr_uexp_t)exp - MPFR_GET_EXP( c );

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

        if ( MPFR_UNLIKELY( MPFR_UEXP( aq2 ) <= diff_exp ) )
        { /* c does not overlap with a' */
            if ( MPFR_UNLIKELY( an > bn ) )
            { /* a has more limbs than b */
              /* copy b to the most significant limbs of a */
                MPN_COPY( ap + (an - bn), bp, bn );
                /* zero the least significant limbs of a */
                MPN_ZERO( ap, an - bn );
            }
            else /* an <= bn */
            {
                /* copy the most significant limbs of b to a */
                MPN_COPY( ap, bp + (bn - an), an );
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

            dif = aq2 - diff_exp;
            /* dif is the number of bits of c which overlap with a' */

            difn = MPFR_PREC2LIMBS( dif );
            /* only the highest difn limbs from c have to be considered */
            if ( MPFR_UNLIKELY( difn > cn ) )
            {
                /* c doesn't have enough limbs; take into account the virtual
                   zero limbs now by zeroing the least significant limbs of a' */
                MPFR_ASSERTD( difn - cn <= an );
                MPN_ZERO( ap, difn - cn );
                difn = cn;
            }
            k = diff_exp / GMP_NUMB_BITS;

            /* zero the most significant k limbs of a */
            a2p = ap + (an - k);
            MPN_ZERO( a2p, k );

            shift = diff_exp % GMP_NUMB_BITS;

            if ( MPFR_LIKELY( shift ) )
            {
                MPFR_ASSERTD( a2p - difn >= ap );
                cc = mpn_rshift( a2p - difn, cp + (cn - difn), difn, shift );
                if ( MPFR_UNLIKELY( a2p - difn > ap ) )
                    *(a2p - difn - 1) = cc;
            }
            else
                MPN_COPY( a2p - difn, cp + (cn - difn), difn );

            /* add b to a */
            cc = MPFR_UNLIKELY( an > bn )
                ? atn_add_n( ap + (an - bn), ap + (an - bn), bp, bn )
                : atn_add_n( ap, ap, bp + (bn - an), an );

            if ( MPFR_UNLIKELY( cc ) ) /* carry */
            {
                if ( MPFR_UNLIKELY( exp == __gmpfr_emax ) )
                {
                    inex = mpfr_overflow( a, rnd_mode, MPFR_SIGN( a ) );
                    goto end_of_add;
                }
                exp++;
                rb = (ap[0] >> sh) & 1; /* LSB(a) --> rounding bit after the shift */
                if ( MPFR_LIKELY( sh ) )
                {
                    mp_limb_t mask, bb;

                    mask = MPFR_LIMB_MASK( sh );
                    bb = ap[0] & mask;
                    ap[0] &= (~mask) << 1;
                    if ( bb == 0 )
                        fb = 0;
                    else if ( bb == mask )
                        fb = 1;
                }
                mpn_rshift( ap, ap, an, 1 );
                ap[an - 1] += MPFR_LIMB_HIGHBIT;
                if ( sh && fb < 0 )
                    goto rounding;
            } /* cc */
        } /* aq2 > diff_exp */

      /* non-significant bits of a */
        if ( MPFR_LIKELY( rb < 0 && sh ) )
        {
            mp_limb_t mask, bb;

            mask = MPFR_LIMB_MASK( sh );
            bb = ap[0] & mask;
            ap[0] &= ~mask;
            rb = bb >> (sh - 1);
            if ( MPFR_LIKELY( sh > 1 ) )
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

        /* determine rounding and sticky bits (and possible carry) */

        difw = (mpfr_exp_t)an - (mpfr_exp_t)(diff_exp / GMP_NUMB_BITS);
        /* difw is the number of limbs from b (regarded as having an infinite
           precision) that have already been combined with c; -n if the next
           n limbs from b won't be combined with c. */

        if ( MPFR_UNLIKELY( bn > an ) )
        { /* there are still limbs from b that haven't been taken into account */
            mp_size_t bk;

            if ( fb == 0 && difw <= 0 )
            {
                fb = 1; /* c hasn't been taken into account ==> sticky bit != 0 */
                goto rounding;
            }

            bk = bn - an; /* index of lowest considered limb from b, > 0 */
            while ( difw < 0 )
            { /* ulp(next limb from b) > msb(c) */
                mp_limb_t bb;

                bb = bp[--bk];

                MPFR_ASSERTD( fb != 0 );
                if ( fb > 0 )
                {
                    if ( bb != MP_LIMB_T_MAX )
                    {
                        fb = 1; /* c hasn't been taken into account
                                   ==> sticky bit != 0 */
                        goto rounding;
                    }
                }
                else /* fb not initialized yet */
                {
                    if ( rb < 0 ) /* rb not initialized yet */
                    {
                        rb = bb >> (GMP_NUMB_BITS - 1);
                        bb |= MPFR_LIMB_HIGHBIT;
                    }
                    fb = 1;
                    if ( bb != MP_LIMB_T_MAX )
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

            if ( difw <= cn )
            {
                mp_size_t ck;
                mp_limb_t cprev;
                int difs;

                ck = cn - difw;
                difs = diff_exp % GMP_NUMB_BITS;

                if ( difs == 0 && ck == 0 )
                    goto c_read;

                cprev = ck == cn ? 0 : cp[ck];

                if ( fb < 0 )
                {
                    mp_limb_t bb, cc;

                    if ( difs )
                    {
                        cc = cprev << (GMP_NUMB_BITS - difs);
                        if ( --ck >= 0 )
                        {
                            cprev = cp[ck];
                            cc += cprev >> difs;
                        }
                    }
                    else
                        cc = cp[--ck];

                    bb = bp[--bk] + cc;

                    if ( bb < cc /* carry */
                        && (rb < 0 || (rb ^= 1) == 0)
                        && mpn_add_1( ap, ap, an, MPFR_LIMB_ONE << sh ) )
                    {
                        if ( exp == __gmpfr_emax )
                        {
                            inex = mpfr_overflow( a, rnd_mode, MPFR_SIGN( a ) );
                            goto end_of_add;
                        }
                        exp++;
                        ap[an - 1] = MPFR_LIMB_HIGHBIT;
                        rb = 0;
                    }

                    if ( rb < 0 ) /* rb not initialized yet */
                    {
                        rb = bb >> (GMP_NUMB_BITS - 1);
                        bb <<= 1;
                        bb |= bb >> (GMP_NUMB_BITS - 1);
                    }

                    fb = bb != 0;
                    if ( fb && bb != MP_LIMB_T_MAX )
                        goto rounding;
                } /* fb < 0 */

                while ( bk > 0 )
                {
                    mp_limb_t bb, cc;

                    if ( difs )
                    {
                        if ( ck < 0 )
                            goto c_read;
                        cc = cprev << (GMP_NUMB_BITS - difs);
                        if ( --ck >= 0 )
                        {
                            cprev = cp[ck];
                            cc += cprev >> difs;
                        }
                    }
                    else
                    {
                        if ( ck == 0 )
                            goto c_read;
                        cc = cp[--ck];
                    }

                    bb = bp[--bk] + cc;
                    if ( bb < cc ) /* carry */
                    {
                        fb ^= 1;
                        if ( fb )
                            goto rounding;
                        rb ^= 1;
                        if ( rb == 0 && mpn_add_1( ap, ap, an, MPFR_LIMB_ONE << sh ) )
                        {
                            if ( MPFR_UNLIKELY( exp == __gmpfr_emax ) )
                            {
                                inex = mpfr_overflow( a, rnd_mode, MPFR_SIGN( a ) );
                                goto end_of_add;
                            }
                            exp++;
                            ap[an - 1] = MPFR_LIMB_HIGHBIT;
                        }
                    } /* bb < cc */

                    if ( !fb && bb != 0 )
                    {
                        fb = 1;
                        goto rounding;
                    }
                    if ( fb && bb != MP_LIMB_T_MAX )
                        goto rounding;
                } /* while */

              /* b has entirely been read */

                if ( fb || ck < 0 )
                    goto rounding;
                if ( difs && cprev << (GMP_NUMB_BITS - difs) )
                {
                    fb = 1;
                    goto rounding;
                }
                while ( ck )
                {
                    if ( cp[--ck] )
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
                    bb = bp[--bk];
                    if ( rb < 0 ) /* rb not initialized yet */
                    {
                        rb = bb >> (GMP_NUMB_BITS - 1);
                        bb &= ~MPFR_LIMB_HIGHBIT;
                    }
                    fb = bb != 0;
                } /* fb < 0 */
                if ( fb )
                    goto rounding;
                while ( bk )
                {
                    if ( bp[--bk] )
                    {
                        fb = 1;
                        goto rounding;
                    }
                } /* while */
            } /* difw > cn */
        } /* bn > an */
        else if ( fb != 1 ) /* if fb == 1, the sticky bit is 1 (no possible carry) */
        { /* b has entirely been read */
            if ( difw > cn )
            { /* c has entirely been read */
                if ( rb < 0 )
                    rb = 0;
                fb = 0;
            }
            else if ( diff_exp > MPFR_UEXP( aq2 ) )
            { /* b is followed by at least a zero bit, then by c */
                if ( rb < 0 )
                    rb = 0;
                fb = 1;
            }
            else
            {
                mp_size_t ck;
                int difs;

                MPFR_ASSERTD( difw >= 0 && cn >= difw );
                ck = cn - difw;
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

                    cc = difs ? (MPFR_ASSERTD( ck < cn ),
                        cp[ck] << (GMP_NUMB_BITS - difs)) : cp[--ck];
                    if ( rb < 0 )
                    {
                        rb = cc >> (GMP_NUMB_BITS - 1);
                        cc &= ~MPFR_LIMB_HIGHBIT;
                    }
                    while ( cc == 0 )
                    {
                        if ( ck == 0 )
                        {
                            fb = 0;
                            goto rounding;
                        }
                        cc = cp[--ck];
                    } /* while */
                    fb = 1;
                }
            }
        } /* fb != 1 */

    rounding:
        /* rnd_mode should be one of MPFR_RNDN, MPFR_RNDZ or MPFR_RNDA */
        if ( MPFR_LIKELY( rnd_mode == MPFR_RNDN ) )
        {
            if ( fb == 0 )
            {
                if ( rb == 0 )
                {
                    inex = 0;
                    goto set_exponent;
                }
                /* round to even */
                if ( ap[0] & (MPFR_LIMB_ONE << sh) )
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
            inex = rb || fb ? (MPFR_IS_NEG( a ) ? 1 : -1) : 0;
            goto set_exponent;
        }
        else
        {
            MPFR_ASSERTN( rnd_mode == MPFR_RNDA );
            inex = rb || fb ? (MPFR_IS_POS( a ) ? 1 : -1) : 0;
            if ( inex )
                goto add_one_ulp;
            else
                goto set_exponent;
        }

    add_one_ulp: /* add one unit in last place to a */
        if ( MPFR_UNLIKELY( mpn_add_1( ap, ap, an, MPFR_LIMB_ONE << sh ) ) )
        {
            if ( MPFR_UNLIKELY( exp == __gmpfr_emax ) )
            {
                inex = mpfr_overflow( a, rnd_mode, MPFR_SIGN( a ) );
                goto end_of_add;
            }
            exp++;
            ap[an - 1] = MPFR_LIMB_HIGHBIT;
        }

    set_exponent:
        MPFR_SET_EXP( a, exp );

    end_of_add:
        MPFR_TMP_FREE( marker );
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
        MPFR_TMP_DECL( marker );

        MPFR_TMP_MARK( marker );

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
            limb = atn_add_n_intrinsic( ap, MPFR_MANT( b ), MPFR_MANT( c ), n );
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
            limb = atn_add_n_intrinsic( ap, MPFR_MANT( b ), cp, n );
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
            MPFR_TMP_FREE( marker );
            MPFR_SET_SAME_SIGN( a, b );
            return mpfr_overflow( a, rnd_mode, MPFR_SIGN( a ) );
        }
        MPFR_SET_EXP( a, bx );
        MPFR_SET_SAME_SIGN( a, b );

        MPFR_TMP_FREE( marker );
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
                if ( MPFR_GET_EXP( b ) < MPFR_GET_EXP( c ) )
                    return atn_add1sp( a, c, b, rnd_mode );
                else
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

    struct BenchmarkResult
    {
        double medianTimeNanoseconds;
        double minTimeNanoseconds;
        double maxTimeNanoseconds;
        size_t iterations;
    };

    namespace
    {
        template<typename Func>
        BenchmarkResult benchmark_function( Func func, size_t iterations )
        {
            BenchmarkResult result{ 0.0, std::numeric_limits<double>::max(), 0.0, iterations };
            std::vector<double> times;
            times.reserve( iterations );

            for ( size_t i = 0; i < iterations; ++i )
            {
                auto start = std::chrono::high_resolution_clock::now();
                func();
                auto end = std::chrono::high_resolution_clock::now();

                auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>( end - start ).count();
                double timeNs = static_cast<double>( duration );

                times.push_back( timeNs );
                result.minTimeNanoseconds = std::min( result.minTimeNanoseconds, timeNs );
                result.maxTimeNanoseconds = std::max( result.maxTimeNanoseconds, timeNs );
            }

            // Calculate median
            std::sort( times.begin(), times.end() );
            if ( iterations % 2 == 0 )
            {
                // Even number of elements: average of middle two
                result.medianTimeNanoseconds = (times[iterations / 2 - 1] + times[iterations / 2]) / 2.0;
            }
            else
            {
                // Odd number of elements: middle element
                result.medianTimeNanoseconds = times[iterations / 2];
            }

            return result;
        }
    
    }

    void benchmark_add_n_implementations( size_t limbCount, size_t iterations )
    {
        std::random_device rd;
        std::mt19937_64 gen( rd() );
        std::uniform_int_distribution<mp_limb_t> dist;

        std::cout << "=== Benchmarking atn_add_n implementations ===" << std::endl;
        std::cout << "Generating " << iterations << " sets of random test data..." << std::endl;

        // Pre-generate all test data to avoid timing contamination
        std::vector<std::vector<mp_limb_t>> test_data_a( iterations );
        std::vector<std::vector<mp_limb_t>> test_data_b( iterations );
        std::vector<std::vector<mp_limb_t>> results( iterations );

        for ( size_t iter = 0; iter < iterations; ++iter )
        {
            test_data_a[iter].resize( limbCount );
            test_data_b[iter].resize( limbCount );
            results[iter].resize( limbCount );

            for ( size_t i = 0; i < limbCount; ++i )
            {
                test_data_a[iter][i] = dist( gen );
                test_data_b[iter][i] = dist( gen );
            }
        }

        std::cout << "Data generation complete." << std::endl;
        std::cout << "Limb count: " << limbCount << std::endl;
        std::cout << "Iterations: " << iterations << std::endl;
        std::cout << "Limb size: " << GMP_NUMB_BITS << " bits" << std::endl;
        std::cout << "Total data size: " << (2 * iterations * limbCount * sizeof(mp_limb_t)) / (1024.0 * 1024.0) << " MB" << std::endl;

        // Show build configuration
#ifdef NDEBUG
        std::cout << "Build: RELEASE (optimized)" << std::endl;
#else
        std::cout << "Build: DEBUG (not optimized)" << std::endl;
#endif

#if defined(_MSC_VER)
        std::cout << "Compiler: MSVC " << _MSC_VER << std::endl;
#elif defined(__GNUC__)
        std::cout << "Compiler: GCC " << __GNUC__ << "." << __GNUC_MINOR__ << std::endl;
#elif defined(__clang__)
        std::cout << "Compiler: Clang " << __clang_major__ << "." << __clang_minor__ << std::endl;
#endif
        std::cout << std::endl;

        // Benchmark standard implementation with different data each iteration
        size_t iter_index = 0;
        auto result_standard = benchmark_function( [&]()
            {
                atn_add_n( results[iter_index].data(), 
                          test_data_a[iter_index].data(), 
                          test_data_b[iter_index].data(), 
                          limbCount );
                iter_index = (iter_index + 1) % iterations;
            }, iterations );

        std::cout << "atn_add_n (standard):" << std::endl;
        std::cout << "  Median:  " << result_standard.medianTimeNanoseconds << " ns" << std::endl;
        std::cout << "  Min:     " << result_standard.minTimeNanoseconds << " ns" << std::endl;
        std::cout << "  Max:     " << result_standard.maxTimeNanoseconds << " ns" << std::endl;
        std::cout << std::endl;

        // Benchmark intrinsic implementation with different data each iteration
        iter_index = 0;
        auto result_intrinsic = benchmark_function( [&]()
            {
                atn_add_n_intrinsic( results[iter_index].data(),
                    test_data_a[iter_index].data(),
                    test_data_b[iter_index].data(),
                    limbCount );
                iter_index = (iter_index + 1) % iterations;
            }, iterations );

        std::cout << "atn_add_n_intrinsic:" << std::endl;
        std::cout << "  Median:  " << result_intrinsic.medianTimeNanoseconds << " ns" << std::endl;
        std::cout << "  Min:     " << result_intrinsic.minTimeNanoseconds << " ns" << std::endl;
        std::cout << "  Max:     " << result_intrinsic.maxTimeNanoseconds << " ns" << std::endl;
        std::cout << std::endl;

        double speedup = result_standard.medianTimeNanoseconds / result_intrinsic.medianTimeNanoseconds;

        
        std::cout << "Performance comparison:" << std::endl;
        if ( speedup > 1.0 )
        {
            std::cout << "  Intrinsic version is " << std::fixed << std::setprecision( 2 ) << speedup << "x faster" << std::endl;
            std::cout << "  Speedup: " << std::fixed << std::setprecision( 1 ) << ((speedup - 1.0) * 100.0) << "%" << std::endl;
        }
        else
        {
            std::cout << "  Standard version is " << std::fixed << std::setprecision( 2 ) << (1.0 / speedup) << "x faster" << std::endl;
            std::cout << "  Slowdown: " << std::fixed << std::setprecision( 1 ) << ((1.0 - speedup) * 100.0) << "%" << std::endl;
        }
        std::cout << std::endl;

        // Verify correctness using first test set
        std::vector<mp_limb_t> verify_result1( limbCount );
        std::vector<mp_limb_t> verify_result2( limbCount );

        mp_limb_t carry1 = atn_add_n( verify_result1.data(), test_data_a[0].data(), test_data_b[0].data(), limbCount );
        mp_limb_t carry2 = atn_add_n_intrinsic( verify_result2.data(), test_data_a[0].data(), test_data_b[0].data(), limbCount );

        bool results_match = (carry1 == carry2) && 
                            std::equal( verify_result1.begin(), verify_result1.end(), verify_result2.begin() );

        if ( results_match )
        {
            std::cout << "Correctness: Both implementations produce identical results." << std::endl;
        }
        else
        {
            std::cout << "WARNING: Implementations produce different results!" << std::endl;
            std::cout << "  Standard carry: " << carry1 << std::endl;
            std::cout << "  Intrinsic carry: " << carry2 << std::endl;
            if ( carry1 == carry2 )
            {
                std::cout << "  Carry values match, but result limbs differ!" << std::endl;
            }
        }
    }

} // namespace Athena