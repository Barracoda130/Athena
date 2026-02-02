#include "Number.hpp"

#include <cassert>

namespace
{
	// Macros from gmp-6.2.1/mpn/generic/add_n.c
#define MPN_OVERLAP_P(xp, xsize, yp, ysize)				\
  ((xp) + (xsize) > (yp) && (yp) + (ysize) > (xp))
#define MPN_SAME_OR_INCR2_P(dst, dsize, src, ssize)			\
  ((dst) <= (src) || ! MPN_OVERLAP_P (dst, dsize, src, ssize))
#define MPN_SAME_OR_INCR_P(dst, src, size)				\
  MPN_SAME_OR_INCR2_P(dst, size, src, size)
    

    mp_limb_t mpn_add_n( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
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
}

namespace Athena
{
	void add( Number& a_Result, const Number& a_Num1, const Number& a_Num2, round_t a_Round )
	{
        // Assume both Num1 and Num2 are initialised

        mpfr_exp_t exp = a_Num1.m_Value->_mpfr_exp;

        // Assume exponent is less than max exponent


	}
} // namespace Athena