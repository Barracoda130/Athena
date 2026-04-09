#pragma once
#include "MpfrInclude.hpp"

namespace Athena
{
	mp_limb_t atn_add_n( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );
	void printAtnIntrinsicVersion();
	mp_limb_t atn_sub_n( mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n );
}