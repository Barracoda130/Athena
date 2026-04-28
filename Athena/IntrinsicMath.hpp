#pragma once
#include "MpfrInclude.hpp"
 
// Each number refers to a different implementation of the add_n function
// 0 for the original implementation 
// 1 for the intrinsic implementation 
// 2 for the intrinsic implementation with ADX instructions
// 3 for the carry select implementation
// Subtraction will always use the same intrinsic implementation unless set to 0 in which case it will use the original implementation
#define MATH_INTRINSIC_VERSION 3

namespace Athena
{
	mp_limb_t atn_add_n( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );
	void printAtnIntrinsicVersion();
	mp_limb_t atn_sub_n( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );
}