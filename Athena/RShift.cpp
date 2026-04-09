#include "mpfr-impl.h"

#include <cassert>

mp_limb_t
atn_rshift( mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt )
{
    mp_limb_t high_limb, low_limb;
    unsigned int tnc;
    mp_size_t i;
    mp_limb_t retval;

    assert( n >= 1 );
    assert( cnt >= 1 );
    assert( cnt < GMP_NUMB_BITS );
    assert( MPN_SAME_OR_INCR_P( rp, up, n ) );

    tnc = GMP_NUMB_BITS - cnt;
    high_limb = *up++;
    retval = (high_limb << tnc) & GMP_NUMB_MASK;
    low_limb = high_limb >> cnt;

    for ( i = n - 1; i != 0; i-- )
    {
        high_limb = *up++;
        *rp++ = low_limb | ((high_limb << tnc) & GMP_NUMB_MASK);
        low_limb = high_limb >> cnt;
    }
    *rp = low_limb;

    return retval;
}