#include "IntrinsicMath.hpp"

#include <cassert>
#include <iostream>

#if defined(_MSC_VER)
#include <intrin.h>
#include <immintrin.h>
#endif

namespace
{
    mp_limb_t atn_add_n_original( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );
    mp_limb_t atn_add_n_intrinsic( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );
    mp_limb_t atn_add_n_intrinsic_adx( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );
    mp_limb_t atn_sub_n_original( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );
    mp_limb_t atn_sub_n_intrinsic( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );

    // The original implementation taken directly from the GMP library
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

	// The intrinsic-based implementation using _addcarry_u64 or _addcarry_u32 on MSVC, and __builtin_add_overflow on GCC/Clang.
    mp_limb_t atn_add_n_intrinsic( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {

		// Guard to detect compiler support for intrinsics. If not supported, fall back to the original implementation.
#if defined(_MSC_VER)
        unsigned char carry = 0;

#if GMP_NUMB_BITS == 64
        auto* r = reinterpret_cast<unsigned __int64*>(rp);
        auto* u = reinterpret_cast<const unsigned __int64*>(up);
        auto* v = reinterpret_cast<const unsigned __int64*>(vp);

		// Unroll the loop to process 4 limbs at a time, which can help reduce loop overhead and increase instruction-level parallelism.
        while ( n >= 4 )
        {
            carry = _addcarry_u64( carry, u[0], v[0], &r[0] );
            carry = _addcarry_u64( carry, u[1], v[1], &r[1] );
            carry = _addcarry_u64( carry, u[2], v[2], &r[2] );
            carry = _addcarry_u64( carry, u[3], v[3], &r[3] );
            u += 4; v += 4; r += 4; n -= 4;
        }

		// Handle any remaining limbs that didn't fit into the unrolled loop.
        while ( n-- > 0 )
            carry = _addcarry_u64( carry, *u++, *v++, r++ );

#elif GMP_NUMB_BITS == 32
		// If we're on a 32-bit platform, we can use the 32-bit intrinsics instead. The logic is the same, just with 32-bit types.
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
		// If we don't have access to the appropriate intrinsics, 
        // we can still use the original logic but with a more compact loop structure.
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

		// Unroll the loop to process 4 limbs at a time, which can help reduce loop overhead and increase instruction-level parallelism.
        while ( n >= 4 )
        {
            mp_limb_t sum;
            // Two operations required as doesn't include carry in the instruction
            unsigned char c1 = __builtin_add_overflow( up[0], vp[0], &sum );
            unsigned char c2 = __builtin_add_overflow( sum, carry, &rp[0] );
			// If either overflowed, we have a carry for the next limb
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

		// Handle any remaining limbs that didn't fit into the unrolled loop.
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
        // If we don't have access to the appropriate intrinsics, we can still use the original logic
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

	// An alternative implementation which uses the _addcarryx_u64 intrinsic on supported platforms
    mp_limb_t atn_add_n_intrinsic_adx( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {
		// Guard to detect compiler support for ADX intrinsics. If not supported, fall back to the non-ADX intrinsic implementation.
#if defined(_MSC_VER)
        unsigned char carry = 0;
        
#if GMP_NUMB_BITS == 64
        auto* r = reinterpret_cast<unsigned __int64*>(rp);
        auto* u = reinterpret_cast<const unsigned __int64*>(up);
        auto* v = reinterpret_cast<const unsigned __int64*>(vp);

		// Unroll the loop to process 4 limbs at a time, which can help reduce loop overhead and increase instruction-level parallelism.
        while ( n >= 4 )
        {
            carry = _addcarryx_u64( carry, u[0], v[0], &r[0] );
            carry = _addcarryx_u64( carry, u[1], v[1], &r[1] );
            carry = _addcarryx_u64( carry, u[2], v[2], &r[2] );
            carry = _addcarryx_u64( carry, u[3], v[3], &r[3] );
            u += 4; v += 4; r += 4; n -= 4;
        }
        
		// Handle any remaining limbs that didn't fit into the unrolled loop.
        while ( n-- > 0 )
            carry = _addcarryx_u64( carry, *u++, *v++, r++ );

#elif GMP_NUMB_BITS == 32
		// On 32-bit platforms, we can use the 32-bit version of the ADX intrinsics. The logic is the same, just with 32-bit types.
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
		// Fall back to the non-ADX intrinsic implementation if we're not on MSVC or don't have access to the ADX intrinsics.
		// This also includes the deafault GMP implementation for platforms that don't support intrinsics at all.
        return atn_add_n_intrinsic( rp, up, vp, n );
#endif
    }

    static inline __mmask8 unsigned_lt_epu64_mask( __m512i a, __m512i b )
    {
        // Unsigned compare a < b using a sign-bit bias transform
        const __m512i bias = _mm512_set1_epi64( 0x8000000000000000ULL );
        __m512i ax = _mm512_xor_si512( a, bias );
        __m512i bx = _mm512_xor_si512( b, bias );
        return _mm512_cmplt_epi64_mask( ax, bx );
    }

	// An AVX-512 implementation of add_n using carry-select. This processes 8 limbs at a time, and uses a small scalar prefix pass to resolve carries across the vector lanes.
    mp_limb_t
        mpn_add_n_avx512_carry_select( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {
        mp_limb_t carry = 0;

        auto* r = reinterpret_cast<std::uint64_t*>(rp);
        auto* u = reinterpret_cast<const std::uint64_t*>(up);
        auto* v = reinterpret_cast<const std::uint64_t*>(vp);

        const __m512i one = _mm512_set1_epi64( 1 );
        const __m512i all_ones = _mm512_set1_epi64( ~0ULL );

        while ( n >= 8 )
        {
            // Load 8 limbs from each operand.
            __m512i vu = _mm512_loadu_si512( reinterpret_cast<const void*>(u) );
            __m512i vv = _mm512_loadu_si512( reinterpret_cast<const void*>(v) );

            // Provisional sum assuming carry-in = 0
            __m512i s0 = _mm512_add_epi64( vu, vv );

            // Provisional sum assuming carry-in = 1
            __m512i s1 = _mm512_add_epi64( s0, one );

            // Generate mask:
            // g[i] = 1 iff u[i] + v[i] overflows 64 bits
            __mmask8 gmask = unsigned_lt_epu64_mask( s0, vu );

            // Propagate mask
            __mmask8 pmask = _mm512_cmpeq_epi64_mask( s0, all_ones );

            // Resolve the real per-lane carry chain.
            std::uint8_t g = static_cast<std::uint8_t>(gmask);
            std::uint8_t p = static_cast<std::uint8_t>(pmask);

            std::uint8_t cin_mask = 0;
            std::uint8_t c = static_cast<std::uint8_t>(carry);

            for ( int i = 0; i < 8; ++i )
            {
                // If current carry-in is 1, mark lane i to pick s1[i] later.
                if ( c )
                    cin_mask |= static_cast<std::uint8_t>( 1u << i );

                // Extract lane i generate/propagate bits from packed masks.
                std::uint8_t gi = (g >> i) & 1u;
                std::uint8_t pi = (p >> i) & 1u;

                // Advance carry to feed next lane: c becomes cout[i].
                c = static_cast<std::uint8_t>( gi | (pi & c) );
            }

            // Select s1 where the real carry-in was 1, else s0.
            __m512i res = _mm512_mask_blend_epi64( static_cast<__mmask8>(cin_mask), s0, s1 );
            _mm512_storeu_si512( reinterpret_cast<void*>(r), res );

            carry = c;
            u += 8;
            v += 8;
            r += 8;
            n -= 8;
        }

        // Scalar cleanup for any remaining limbs.
        while ( n-- > 0 )
        {
            unsigned __int64 out;
            carry = _addcarry_u64(
                static_cast<unsigned char>(carry),
                static_cast<unsigned __int64>(*u++),
                static_cast<unsigned __int64>(*v++),
                &out
            );
            *r++ = static_cast<std::uint64_t>(out);
        }

        return carry;
    }

    /****************************************************************************
     * Subtraction
     ****************************************************************************/

	 // Original implementation taken directly from the GMP library
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

	// The intrinsic-based implementation using _subborrow_u64 or _subborrow_u32 on MSVC, and __builtin_sub_overflow on GCC/Clang.
    mp_limb_t atn_sub_n_intrinsic( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {
		// Guard to detect compiler support for intrinsics. If not supported, fall back to the original implementation.
#if defined(_MSC_VER)
        unsigned char borrow = 0;

#if GMP_NUMB_BITS == 64
        auto* r = reinterpret_cast<unsigned __int64*>(rp);
        auto* u = reinterpret_cast<const unsigned __int64*>(up);
        auto* v = reinterpret_cast<const unsigned __int64*>(vp);

		// Unroll the loop to process 4 limbs at a time, which can help reduce loop overhead and increase instruction-level parallelism.
        while ( n >= 4 )
        {
            borrow = _subborrow_u64( borrow, u[0], v[0], &r[0] );
            borrow = _subborrow_u64( borrow, u[1], v[1], &r[1] );
            borrow = _subborrow_u64( borrow, u[2], v[2], &r[2] );
            borrow = _subborrow_u64( borrow, u[3], v[3], &r[3] );
            u += 4; v += 4; r += 4; n -= 4;
        }

		// Handle any remaining limbs that didn't fit into the unrolled loop.
        while ( n-- > 0 )
            borrow = _subborrow_u64( borrow, *u++, *v++, r++ );

#elif GMP_NUMB_BITS == 32
		// On 32-bit platforms, we can use the 32-bit version of the intrinsics instead. The logic is the same, just with 32-bit types.
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
		// If we don't have access to the appropriate intrinsics, we can still use the original logic
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
		// Unroll the loop to process 4 limbs at a time, which can help reduce loop overhead and increase instruction-level parallelism.
        mp_limb_t borrow = 0;
        while ( n >= 4 )
        {
            mp_limb_t diff;
			// Two operations required as doesn't include borrow in the instruction
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

		// Handle any remaining limbs that didn't fit into the unrolled loop.
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
		// If we don't have access to the appropriate intrinsics, we can still use the original logic
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
#elif MATH_INTRINSIC_VERSION == 3
			"intrinsic_avx512_carry_select"
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
#elif MATH_INTRINSIC_VERSION == 3
        return mpn_add_n_avx512_carry_select( rp, up, vp, n );
#else
        // Throw error
#error "Invalid INTRINSIC_VERSION. Expected 0, 1, or 2."
#endif
    }

    mp_limb_t atn_sub_n( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n )
    {
#if MATH_INTRINSIC_VERSION == 0
        return atn_sub_n_original( rp, up, vp, n );
#elif MATH_INTRINSIC_VERSION > 1
        return atn_sub_n_intrinsic( rp, up, vp, n );
#endif
    }
}