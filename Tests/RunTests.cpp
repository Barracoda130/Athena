#define CATCH_CONFIG_RUNNER

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_session.hpp>

#include <iostream>
#include <functional>
#include <format>

#include "TestDataFile.hpp"
#include "MpfrInclude.hpp"
#include "Athena.hpp"

#define TEST_PRECISION 512
#define NUM_TESTS 10000


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
			  std::function<void( mpfr_t, mpfr_t, mpfr_t, mpfr_rnd_t )> a_F1,
			  std::function<void( Athena::Number&, const Athena::Number&, const Athena::Number&, Athena::round_t )> a_F2 )
{
	const char* n1 = a_File.getNextA();
	const char* n2 = a_File.getNextB();
	while ( !a_File.eofReached() )
	{
		// Run through MPFR
		mpfr_t mpfrResult;

		mpfr_t num1, num2;
		mpfr_init2( num1, TEST_PRECISION );
		mpfr_init2( num2, TEST_PRECISION );
		mpfr_init2( mpfrResult, TEST_PRECISION );

		mpfr_set_str( num1, n1, 10, MPFR_RNDN );
		mpfr_set_str( num2, n2, 10, MPFR_RNDN );

		a_F1( mpfrResult, num1, num2, MPFR_RNDD );

		Athena::Number n1Atna( n1, TEST_PRECISION );
		Athena::Number n2Atna( n2, TEST_PRECISION );
		Athena::Number resultAtna( TEST_PRECISION );

		a_F2( resultAtna, n1Atna, n2Atna, MPFR_RNDD );

		if ( resultAtna == Athena::Number( mpfrResult, MPFR_RNDN ) )
			SUCCEED();
		else
		{
			INFO( std::format( "{}\n{}\nMPFR: {}\n ATNA: {}", 
				Athena::Number(n1, TEST_PRECISION ).str(),
				Athena::Number(n2, TEST_PRECISION ).str(),
				Athena::Number( mpfrResult, MPFR_RNDN ).str(),
				resultAtna.str() ) );
			FAIL();
		}

		n1 = a_File.getNextA();
		n2 = a_File.getNextB();
	}
}
	
void runTest( TestData::FileReader & a_File,
				std::function<void( mpfr_t, mpfr_t, mpfr_t, mpfr_rnd_t )> a_F1,
				std::function<long long( long long, long long )> a_F2 )
{
	const char* n1 = a_File.getNextA();
	const char* n2 = a_File.getNextB();
	while ( !a_File.eofReached() )
	{
		// Run through MPFR
		mpfr_t mpfrResult;

		mpfr_t num1, num2;
		mpfr_init2( num1, TEST_PRECISION );
		mpfr_init2( num2, TEST_PRECISION );
		mpfr_init2( mpfrResult, TEST_PRECISION );

		mpfr_set_str( num1, n1, 10, MPFR_RNDD );
		mpfr_set_str( num2, n2, 10, MPFR_RNDD );

		a_F1( mpfrResult, num1, num2, MPFR_RNDD );

		std::string mpfrResultString = mpfr_tToStr( mpfrResult );

		long long n1Int = std::stoll( n1 );
		long long n2Int = std::stoll( n2 );
		std::string longResult = std::to_string( a_F2(n1Int, n2Int) );

		REQUIRE( mpfrResultString == longResult );

		n1 = a_File.getNextA();
		n2 = a_File.getNextB();
	}
}

void initTest()
{
	// Create a file reader, the read types don't matter
	TestData::FileReader file( false, TestData::POSITIVE, false, TestData::POSITIVE );

	// If eof is reached then the FileReader class has already closed the file
	// If the test data doesn't yet exist, generate it
	if ( file.eofReached() )
		TestData::generateData( NUM_TESTS );
}

void cleanupTest() {}

TEST_CASE("Positive Addition")
{
	initTest();
	TestData::FileReader file( false, TestData::POSITIVE, false, TestData::POSITIVE );
	runTest( file, mpfr_add, Athena::add );
	cleanupTest();
}

TEST_CASE("Negative Addition")
{
	initTest();
	TestData::FileReader file( false, TestData::NEGATIVE, false, TestData::NEGATIVE );
	runTest( file, mpfr_add, Athena::add );
	cleanupTest();
}

TEST_CASE( "Positive and Negative Addition" )
{
	initTest();
	TestData::FileReader file( false, TestData::POSITIVE, false, TestData::NEGATIVE );
	runTest( file, mpfr_add, Athena::add );
	
	// Now run it with negative + positive rather than positive + negative
	// to ensure both ways round work
	file.reset( false, TestData::NEGATIVE, false, TestData::POSITIVE );
	
	runTest( file, mpfr_add, Athena::add );
	cleanupTest();
}

TEST_CASE( "Random Addition" )
{
	initTest();
	TestData::FileReader file( false, TestData::RANDOM, false, TestData::RANDOM );
	runTest( file, mpfr_add, Athena::add );
	cleanupTest();
}

// SUBTRACTION
TEST_CASE( "Positive Subtraction" )
{
	initTest();
	TestData::FileReader file( false, TestData::POSITIVE, false, TestData::POSITIVE );
	runTest( file, mpfr_sub, Athena::sub );
	cleanupTest();
}

TEST_CASE( "Negative Subtraction" )
{
	initTest();
	TestData::FileReader file( false, TestData::NEGATIVE, false, TestData::NEGATIVE );
	runTest( file, mpfr_sub, Athena::sub );
	cleanupTest();
}

TEST_CASE( "Positive and Negative Subtraction" )
{
	initTest();
	TestData::FileReader file( false, TestData::POSITIVE, false, TestData::NEGATIVE );
	runTest( file, mpfr_sub, Athena::sub );

	// Now run it with negative + positive rather than positive + negative
	// to ensure both ways round work
	file.reset( false, TestData::NEGATIVE, false, TestData::POSITIVE );

	runTest( file, mpfr_sub, Athena::sub );
	cleanupTest();
}

TEST_CASE( "Random Subtraction" )
{
	initTest();
	TestData::FileReader file( false, TestData::RANDOM, false, TestData::RANDOM );
	runTest( file, mpfr_sub, Athena::sub );
	cleanupTest();
}

// MULTIPLICATION
TEST_CASE( "Positive Multiplication" )
{
	initTest();
	TestData::FileReader file( false, TestData::POSITIVE, false, TestData::POSITIVE );
	runTest( file, mpfr_mul, Athena::mul );
	cleanupTest();
}

TEST_CASE( "Negative Multiplication" )
{
	initTest();
	TestData::FileReader file( false, TestData::NEGATIVE, false, TestData::NEGATIVE );
	runTest( file, mpfr_mul, Athena::mul );
	cleanupTest();
}

TEST_CASE( "Positive and Negative Multiplication" )
{
	initTest();
	TestData::FileReader file( false, TestData::POSITIVE, false, TestData::NEGATIVE );
	runTest( file, mpfr_mul, Athena::mul );

	// Now run it with negative + positive rather than positive + negative
	// to ensure both ways round work
	file.reset( false, TestData::NEGATIVE, false, TestData::POSITIVE );

	runTest( file, mpfr_mul, Athena::mul );
	cleanupTest();
}

TEST_CASE( "Random Multiplication" )
{
	initTest();
	TestData::FileReader file( false, TestData::RANDOM, false, TestData::RANDOM );
	runTest( file, mpfr_mul, Athena::mul );
	cleanupTest();
}

// DIVISION
TEST_CASE( "Positive Divison" )
{
	initTest();
	TestData::FileReader file( false, TestData::POSITIVE, false, TestData::POSITIVE );
	runTest( file, mpfr_div, Athena::div );
	cleanupTest();
}

TEST_CASE( "Negative Divison" )
{
	initTest();
	TestData::FileReader file( false, TestData::NEGATIVE, false, TestData::NEGATIVE );
	runTest( file, mpfr_div, Athena::div );
	cleanupTest();
}

TEST_CASE( "Positive and Negative Divison" )
{
	initTest();
	TestData::FileReader file( false, TestData::POSITIVE, false, TestData::NEGATIVE );
	runTest( file, mpfr_div, Athena::div );

	// Now run it with negative + positive rather than positive + negative
	// to ensure both ways round work
	file.reset( false, TestData::NEGATIVE, false, TestData::POSITIVE );

	runTest( file, mpfr_div, Athena::div );
	cleanupTest();
}

TEST_CASE( "Random Divison" )
{
	initTest();
	TestData::FileReader file( false, TestData::RANDOM, false, TestData::RANDOM );
	runTest( file, mpfr_div, Athena::div );
	cleanupTest();
}

int main( int argc, char* argv[] )
{
    Catch::Session session;
    
    // Override command line to run specific test
    //const char* customArgs[] = { argv[0], "Positive Addition" };
    int returnCode = session.applyCommandLine( 1, argv );
    if ( returnCode != 0 )
        return returnCode;

    returnCode = session.run();
    return returnCode;
}