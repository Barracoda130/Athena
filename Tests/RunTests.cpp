#include "RunTests.hpp"

#include <iostream>
#include "CppUnitTest.h"

#include "TestDataFile.hpp"
#include "MpfrInclude.hpp"


void testPosInt()
{
	TestData::FileReader file( false, TestData::POSITIVE, false, TestData::POSITIVE );

	// Run through MPFR
	mpfr_t num1, num2;
	mpfr_init2( num1, 512 );
	mpfr_init2( num2, 512 );

	mpfr_set_str( num1, file.getNextA(), 10, MPFR_RNDD );
	mpfr_set_str( num2, file.getNextB(), 10, MPFR_RNDD );

	mpfr_add( num1, num1, num2, MPFR_RNDD );

	mpfr_exp_t exp;
	char* mpfrResult = mpfr_get_str( NULL, &exp, 10, 0, num1, MPFR_RNDD );

	std::cout << "str: " << mpfrResult << std::endl;
	std::cout << "exp: " << exp << std::endl;
}

void testNegInt()
{}

void testPosFloat()
{}

void testNegFloat()
{}

void testInt()
{}

void testFloat()
{}

void testAll()
{}
