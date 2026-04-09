#include "IntrinsicMath.hpp"

#include <cassert>
#include <iostream>

#if defined(_MSC_VER)
#include <intrin.h>
#include <immintrin.h>
#endif

// Each number refers to a different implementation of the add_n function
// 0 for the original implementation 
// 1 for the intrinsic implementation 
// 2 for the intrinsic implementation with ADX instructions
#define MATH_INTRINSIC_VERSION 2

namespace
{
    mp_limb_t atn_add_n_original( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );
    mp_limb_t atn_add_n_intrinsic( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );
    mp_limb_t atn_add_n_intrinsic_adx( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );
    mp_limb_t atn_sub_n_original( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );
    mp_limb_t atn_sub_n_intrinsic( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );
    mp_limb_t atn_sub_n_intrinsic_adx( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );

    mp_limb_t atn_add_n_original( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {
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
    }

    mp_limb_t atn_add_n_intrinsic( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {

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

        return static_cast<mp_limb_t>(carry);
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

    mp_limb_t atn_add_n_intrinsic_adx( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {

#if defined(_MSC_VER)
        unsigned char carry = 0;

#if GMP_NUMB_BITS == 64
        auto* r = reinterpret_cast<unsigned __int64*>(rp);
        auto* u = reinterpret_cast<const unsigned __int64*>(up);
        auto* v = reinterpret_cast<const unsigned __int64*>(vp);

        while ( n >= 4 )
        {
            carry = _addcarryx_u64( carry, u[0], v[0], &r[0] );
            carry = _addcarryx_u64( carry, u[1], v[1], &r[1] );
            carry = _addcarryx_u64( carry, u[2], v[2], &r[2] );
            carry = _addcarryx_u64( carry, u[3], v[3], &r[3] );
            u += 4; v += 4; r += 4; n -= 4;
        }
        while ( n-- > 0 )
            carry = _addcarryx_u64( carry, *u++, *v++, r++ );
#elif GMP_NUMB_BITS == 32
        auto* r = reinterpret_cast<unsigned int*>(rp);
        auto* u = reinterpret_cast<const unsigned int*>(up);
        auto* v = reinterpret_cast<const unsigned int*>(vp);

        while ( n >= 4 )
        {
            carry = _addcarryx_u32( carry, u[0], v[0], &r[0] );
            carry = _addcarryx_u32( carry, u[1], v[1], &r[1] );
            carry = _addcarryx_u32( carry, u[2], v[2], &r[2] );
            carry = _addcarryx_u32( carry, u[3], v[3], &r[3] );
            u += 4; v += 4; r += 4; n -= 4;
        }
        while ( n-- > 0 )
            carry = _addcarryx_u32( carry, *u++, *v++, r++ );
#else
        return atn_add_n_intrinsic( rp, up, vp, n );
#endif

        return static_cast<mp_limb_t>(carry);
#else
        return atn_add_n_intrinsic( rp, up, vp, n );
#endif
    }

    /****************************************************************************
     * Subtraction
     ****************************************************************************/

    mp_limb_t
        atn_sub_n_original( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {
        mp_limb_t ul, vl, sl, rl, cy, cy1, cy2;

        cy = 0;
        do
        {
            ul = *up++;
            vl = *vp++;
            sl = ul - vl;
            cy1 = sl > ul;
            rl = sl - cy;
            cy2 = rl > sl;
            cy = cy1 | cy2;
            *rp++ = rl;
        } while ( --n != 0 );

        return cy;
    }

    mp_limb_t atn_sub_n_intrinsic( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {
#if defined(_MSC_VER)
        unsigned char borrow = 0;

#if GMP_NUMB_BITS == 64
        auto* r = reinterpret_cast<unsigned __int64*>(rp);
        auto* u = reinterpret_cast<const unsigned __int64*>(up);
        auto* v = reinterpret_cast<const unsigned __int64*>(vp);

        while ( n >= 4 )
        {
            borrow = _subborrow_u64( borrow, u[0], v[0], &r[0] );
            borrow = _subborrow_u64( borrow, u[1], v[1], &r[1] );
            borrow = _subborrow_u64( borrow, u[2], v[2], &r[2] );
            borrow = _subborrow_u64( borrow, u[3], v[3], &r[3] );
            u += 4; v += 4; r += 4; n -= 4;
        }
        while ( n-- > 0 )
            borrow = _subborrow_u64( borrow, *u++, *v++, r++ );

#elif GMP_NUMB_BITS == 32
        auto* r = reinterpret_cast<unsigned int*>(rp);
        auto* u = reinterpret_cast<const unsigned int*>(up);
        auto* v = reinterpret_cast<const unsigned int*>(vp);

        while ( n >= 4 )
        {
            borrow = _subborrow_u32( borrow, u[0], v[0], &r[0] );
            borrow = _subborrow_u32( borrow, u[1], v[1], &r[1] );
            borrow = _subborrow_u32( borrow, u[2], v[2], &r[2] );
            borrow = _subborrow_u32( borrow, u[3], v[3], &r[3] );
            u += 4; v += 4; r += 4; n -= 4;
        }
        while ( n-- > 0 )
            borrow = _subborrow_u32( borrow, *u++, *v++, r++ );

#else
        mp_limb_t ul, vl, sl, rl, cy, cy1, cy2;
        cy = 0;
        mp_size_t i = 0;
        do
        {
            ul = up[i];
            vl = vp[i];
            sl = ul - vl;
            cy1 = sl > ul;
            rl = sl - cy;
            cy2 = rl > sl;
            cy = cy1 | cy2;
            rp[i] = rl;
            i++;
        } while ( i < n );
        borrow = static_cast<unsigned char>( cy );
#endif

        return static_cast<mp_limb_t>( borrow );
#elif defined(__GNUC__) || defined(__clang__)
        mp_limb_t borrow = 0;
        while ( n >= 4 )
        {
            mp_limb_t diff;
            unsigned char b1 = __builtin_sub_overflow( up[0], vp[0], &diff );
            unsigned char b2 = __builtin_sub_overflow( diff, borrow, &rp[0] );
            borrow = b1 | b2;

            b1 = __builtin_sub_overflow( up[1], vp[1], &diff );
            b2 = __builtin_sub_overflow( diff, borrow, &rp[1] );
            borrow = b1 | b2;

            b1 = __builtin_sub_overflow( up[2], vp[2], &diff );
            b2 = __builtin_sub_overflow( diff, borrow, &rp[2] );
            borrow = b1 | b2;

            b1 = __builtin_sub_overflow( up[3], vp[3], &diff );
            b2 = __builtin_sub_overflow( diff, borrow, &rp[3] );
            borrow = b1 | b2;

            up += 4; vp += 4; rp += 4; n -= 4;
        }
        while ( n-- > 0 )
        {
            mp_limb_t diff;
            unsigned char b1 = __builtin_sub_overflow( *up++, *vp++, &diff );
            unsigned char b2 = __builtin_sub_overflow( diff, borrow, rp );
            borrow = b1 | b2;
            ++rp;
        }
        return borrow;
#else
        mp_limb_t ul, vl, sl, rl, cy, cy1, cy2;
        cy = 0;
        do
        {
            ul = *up++;
            vl = *vp++;
            sl = ul - vl;
            cy1 = sl > ul;
            rl = sl - cy;
            cy2 = rl > sl;
            cy = cy1 | cy2;
            *rp++ = rl;
        } while ( --n != 0 );
        return cy;
#endif
    }

    mp_limb_t atn_sub_n_intrinsic_adx( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {
#if defined(_MSC_VER)
        unsigned char borrow = 0;

#if GMP_NUMB_BITS == 64
        auto* r = reinterpret_cast<unsigned __int64*>(rp);
        auto* u = reinterpret_cast<const unsigned __int64*>(up);
        auto* v = reinterpret_cast<const unsigned __int64*>(vp);

        while ( n >= 4 )
        {
            borrow = _subborrow_u64( borrow, u[0], v[0], &r[0] );
            borrow = _subborrow_u64( borrow, u[1], v[1], &r[1] );
            borrow = _subborrow_u64( borrow, u[2], v[2], &r[2] );
            borrow = _subborrow_u64( borrow, u[3], v[3], &r[3] );
            u += 4; v += 4; r += 4; n -= 4;
        }
        while ( n-- > 0 )
            borrow = _subborrow_u64( borrow, *u++, *v++, r++ );

#elif GMP_NUMB_BITS == 32
        auto* r = reinterpret_cast<unsigned int*>(rp);
        auto* u = reinterpret_cast<const unsigned int*>(up);
        auto* v = reinterpret_cast<const unsigned int*>(vp);

        while ( n >= 4 )
        {
            borrow = _subborrow_u32( borrow, u[0], v[0], &r[0] );
            borrow = _subborrow_u32( borrow, u[1], v[1], &r[1] );
            borrow = _subborrow_u32( borrow, u[2], v[2], &r[2] );
            borrow = _subborrow_u32( borrow, u[3], v[3], &r[3] );
            u += 4; v += 4; r += 4; n -= 4;
        }
        while ( n-- > 0 )
            borrow = _subborrow_u32( borrow, *u++, *v++, r++ );
#else
        return atn_sub_n_intrinsic( rp, up, vp, n );
#endif

        return static_cast<mp_limb_t>( borrow );
#else
        return atn_sub_n_intrinsic( rp, up, vp, n );
#endif
    }

}

namespace Athena
{
    void printAtnIntrinsicVersion()
    {
        std::cout <<
#if MATH_INTRINSIC_VERSION == 0
            "original"
#elif MATH_INTRINSIC_VERSION == 1
            "intrinsic"
#elif MATH_INTRINSIC_VERSION == 2
            "intrinsic_adx"
#endif
            << std::endl;
    }

    mp_limb_t atn_add_n( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {
#if MATH_INTRINSIC_VERSION == 0
        return atn_add_n_original( rp, up, vp, n );
#elif MATH_INTRINSIC_VERSION == 1
        return atn_add_n_intrinsic( rp, up, vp, n );
#elif MATH_INTRINSIC_VERSION == 2
        return atn_add_n_intrinsic_adx( rp, up, vp, n );
#else
        // Throw error
#error "Invalid INTRINSIC_VERSION. Expected 0, 1, or 2."
#endif
    }

    mp_limb_t atn_sub_n( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {
#if MATH_INTRINSIC_VERSION == 0
        return atn_sub_n_original( rp, up, vp, n );
#elif MATH_INTRINSIC_VERSION == 1
        return atn_sub_n_intrinsic( rp, up, vp, n );
#elif MATH_INTRINSIC_VERSION == 2
        return atn_sub_n_intrinsic_adx( rp, up, vp, n );
#else
#error "Invalid INTRINSIC_VERSION. Expected 0, 1, or 2."
#endif
    }
}