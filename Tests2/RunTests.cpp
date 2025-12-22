#include "CppUnitTest.h"

#include <iostream>
#include <functional>

#include "TestDataFile.hpp"
#include "MpfrInclude.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;



static std::string mpfr_tToStr( mpfr_t a_Num )
{

	mpfr_exp_t exp;
	char* charArray = mpfr_get_str( NULL, &exp, 10, 0, a_Num, MPFR_RNDD );

	std::string result( charArray );

	// If number is negative add 1 to exp
	if ( result[0] == '-' )
		exp += 1;

	// If exp is negative we need to pad with 0s
	if ( exp < 0 )
	{
		// TODO
	}
	else if ( exp == 0 )
	{
		result.insert( 0, "0." );
	}
	else
	{
		result.insert( exp, "." );
	}


	// Look for trailing 0s after decimal point and remove them
	index_t lastZero = 0;
	index_t decimalPointIdx = 0;
	for ( index_t i = 0; i < result.size(); i++ )
	{
		if ( result[i] == '.' )
			decimalPointIdx = i;

		else if ( result[i] == '0' )
		{
			if ( lastZero == 0 && decimalPointIdx != 0 )	// Make sure we don't overrite a found 0
				lastZero = i;
		}
		else
			lastZero = 0;
	}

	if ( lastZero != 0 )
	{
		if ( decimalPointIdx + 1 == lastZero )
			result = result.substr( 0, lastZero - 1 );
		else
			result = result.substr( 0, lastZero );
	}

	return result;
}

void runTest( TestData::FileReader& a_File, 
			  std::function<void(mpfr_t, mpfr_t, mpfr_t, mpfr_rnd_t)> a_F1,
			  std::function<long long( long long, long long )> a_F2 )
{
	const char* n1 = a_File.getNextA();
	const char* n2 = a_File.getNextB();
	while ( !a_File.eofReached() )
	{
		// Run through MPFR
		mpfr_t mpfrResult;

		mpfr_t num1, num2;
		mpfr_init2( num1, 512 );
		mpfr_init2( num2, 512 );
		mpfr_init2( mpfrResult, 512 );

		mpfr_set_str( num1, n1, 10, MPFR_RNDD );
		mpfr_set_str( num2, n2, 10, MPFR_RNDD );

		a_F1( mpfrResult, num1, num2, MPFR_RNDD );

		std::string mpfrResultString = mpfr_tToStr( mpfrResult );

		long long n1Int = std::stoll( n1 );
		long long n2Int = std::stoll( n2 );

		std::string myResult = std::to_string( a_F2(n1Int, n2Int) );

		Assert::AreEqual( mpfrResultString, myResult );

		n1 = a_File.getNextA();
		n2 = a_File.getNextB();
	}
}

namespace Integers
{
	// Runs once before all tests to ensure that the file has been created
	TEST_MODULE_INITIALIZE( createTestData )
	{
		// Create a file reader, the read types don't matter
		TestData::FileReader file( false, TestData::POSITIVE, false, TestData::POSITIVE );

		// If eof is reached then the FileReader class has already closed the file
		// If the test data doesn't yet exist, generate it
		if ( file.eofReached() )
			TestData::generateData( 10 );
	}

	TEST_CLASS( Addition )
	{
		static long long basicAdd( long long a_N1, long long a_N2 )
		{
			return a_N1 + a_N2;
		}
		TEST_METHOD( positive )
		{
			TestData::FileReader file( false, TestData::POSITIVE, false, TestData::POSITIVE );

			runTest( file, mpfr_add, basicAdd );
		}

		TEST_METHOD( negative )
		{
			TestData::FileReader file( false, TestData::NEGATIVE, false, TestData::NEGATIVE );

			runTest( file, mpfr_add, basicAdd );
		}

		TEST_METHOD( positiveAndNegative )
		{
			TestData::FileReader file( false, TestData::POSITIVE, false, TestData::NEGATIVE );

			runTest( file, mpfr_add, basicAdd );

			// Now run it with negative + positive rather than positive + negative
			// to ensure both ways round work
			file.reset( false, TestData::NEGATIVE, false, TestData::POSITIVE );

			runTest( file, mpfr_add, basicAdd );
		}
	};

	TEST_CLASS( Subtraction )
	{
		static long long basicSub( long long a_N1, long long a_N2 )
		{
			return a_N1 - a_N2;
		}

		TEST_METHOD( positive )
		{
			TestData::FileReader file( false, TestData::POSITIVE, false, TestData::POSITIVE );

			runTest( file, mpfr_sub, basicSub );
		}

		TEST_METHOD( negative )
		{
			TestData::FileReader file( false, TestData::NEGATIVE, false, TestData::NEGATIVE );

			runTest( file, mpfr_sub, basicSub );
		}

		TEST_METHOD( positiveAndNegative )
		{
			TestData::FileReader file( false, TestData::POSITIVE, false, TestData::NEGATIVE );

			runTest( file, mpfr_sub, basicSub );

			// Now run it with negative + positive rather than positive + negative
			// to ensure both ways round work
			file.reset( false, TestData::NEGATIVE, false, TestData::POSITIVE );

			runTest( file, mpfr_sub, basicSub );
		}
	};

	TEST_CLASS( Multiplication )
	{
		static long long basicMult( long long a_N1, long long a_N2 )
		{
			return a_N1 * a_N2;
		}

		TEST_METHOD( positive )
		{
			TestData::FileReader file( false, TestData::POSITIVE, false, TestData::POSITIVE );

			runTest( file, mpfr_mul, basicMult );
		}

		TEST_METHOD( negative )
		{
			TestData::FileReader file( false, TestData::NEGATIVE, false, TestData::NEGATIVE );

			runTest( file, mpfr_mul, basicMult );
		}

		TEST_METHOD( positiveAndNegative )
		{
			TestData::FileReader file( false, TestData::POSITIVE, false, TestData::NEGATIVE );

			runTest( file, mpfr_mul, basicMult );

			// Now run it with negative + positive rather than positive + negative
			// to ensure both ways round work
			file.reset( false, TestData::NEGATIVE, false, TestData::POSITIVE );

			runTest( file, mpfr_mul, basicMult );
		}
	};

	TEST_CLASS( Division )
	{
		static long long basicSub( long long a_N1, long long a_N2 )
		{
			return a_N1 - a_N2;
		}

		TEST_METHOD( positive )
		{
			TestData::FileReader file( false, TestData::POSITIVE, false, TestData::POSITIVE );

			runTest( file, mpfr_div, basicSub );
		}

		TEST_METHOD( negative )
		{
			TestData::FileReader file( false, TestData::NEGATIVE, false, TestData::NEGATIVE );

			runTest( file, mpfr_div, basicSub );
		}

		TEST_METHOD( positiveAndNegative )
		{
			TestData::FileReader file( false, TestData::POSITIVE, false, TestData::NEGATIVE );

			runTest( file, mpfr_div, basicSub );

			// Now run it with negative + positive rather than positive + negative
			// to ensure both ways round work
			file.reset( false, TestData::NEGATIVE, false, TestData::POSITIVE );

			runTest( file, mpfr_div, basicSub );
		}
	};

}
