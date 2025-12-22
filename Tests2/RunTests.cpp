#include "CppUnitTest.h"

#include <iostream>

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

namespace Integers
{
	TEST_CLASS( Addition )
	{
		void simpleMpfrAdd( const char* a_Num1, const char* a_Num2, mpfr_t& a_Result )
		{
			mpfr_t num1, num2;
			mpfr_init2( num1, 512 );
			mpfr_init2( num2, 512 );
			mpfr_init2( a_Result, 512 );

			mpfr_set_str( num1, a_Num1, 10, MPFR_RNDD );
			mpfr_set_str( num2, a_Num2, 10, MPFR_RNDD );

			mpfr_add( a_Result, num1, num2, MPFR_RNDD );
		}

		// Runs once before all tests to ensure that the file has been created
		TEST_CLASS_INITIALIZE( createTestData )
		{
			// Create a file reader, the read types don't matter
			TestData::FileReader file( false, TestData::POSITIVE, false, TestData::POSITIVE );

			// If eof is reached then the FileReader class has already closed the file
			// If the test data doesn't yet exist, generate it
			if ( file.eofReached() )
				TestData::generateData( 10 );
		}
		TEST_METHOD( positive )
		{
			TestData::FileReader file( false, TestData::POSITIVE, false, TestData::POSITIVE );

			const char* n1 = file.getNextA();
			const char* n2 = file.getNextB();
			while ( !file.eofReached() )
			{
				// Run through MPFR
				mpfr_t mpfrResult;

				simpleMpfrAdd( n1, n2, mpfrResult );
				std::string mpfrResultString = mpfr_tToStr( mpfrResult );

				long long n1Int = std::stoll( n1 );
				long long n2Int = std::stoll( n2 );

				std::string myResult = std::to_string( n1Int + n2Int );

				Assert::AreEqual( mpfrResultString, myResult );

				n1 = file.getNextA();
				n2 = file.getNextB();
			}
		}

		TEST_METHOD( negative )
		{
			TestData::FileReader file( false, TestData::NEGATIVE, false, TestData::NEGATIVE );

			const char* n1 = file.getNextA();
			const char* n2 = file.getNextB();
			while ( !file.eofReached() )
			{
				// Run through MPFR
				mpfr_t mpfrResult;

				simpleMpfrAdd( n1, n2, mpfrResult );
				std::string mpfrResultString = mpfr_tToStr( mpfrResult );

				long long n1Int = std::stoll( n1 );
				long long n2Int = std::stoll( n2 );

				std::string myResult = std::to_string( n1Int + n2Int );

				Assert::AreEqual( mpfrResultString, myResult );

				n1 = file.getNextA();
				n2 = file.getNextB();
			}
		}

		TEST_METHOD( positiveAndNegative )
		{
			TestData::FileReader file( false, TestData::POSITIVE, false, TestData::NEGATIVE );

			const char* n1 = file.getNextA();
			const char* n2 = file.getNextB();
			while ( !file.eofReached() )
			{
				// Run through MPFR
				mpfr_t mpfrResult;

				simpleMpfrAdd( n1, n2, mpfrResult );
				std::string mpfrResultString = mpfr_tToStr( mpfrResult );

				long long n1Int = std::stoll( n1 );
				long long n2Int = std::stoll( n2 );

				std::string myResult = std::to_string( n1Int + n2Int );

				Assert::AreEqual( mpfrResultString, myResult );

				n1 = file.getNextA();
				n2 = file.getNextB();
			}

			// Now run it with negative + positive rather than positive + negative
			// to ensure both ways round work
			file.reset( false, TestData::NEGATIVE, false, TestData::POSITIVE );

			n1 = file.getNextA();
			n2 = file.getNextB();
			while ( !file.eofReached() )
			{
				// Run through MPFR
				mpfr_t mpfrResult;

				simpleMpfrAdd( n1, n2, mpfrResult );
				std::string mpfrResultString = mpfr_tToStr( mpfrResult );

				long long n1Int = std::stoll( n1 );
				long long n2Int = std::stoll( n2 );

				std::string myResult = std::to_string( n1Int + n2Int );

				Assert::AreEqual( mpfrResultString, myResult );

				n1 = file.getNextA();
				n2 = file.getNextB();
			}
		}
	};

}
