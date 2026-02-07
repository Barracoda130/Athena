#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

int add( int a, int b )
{
	return a + b;
}

void t1()
{
	REQUIRE( 1 == 3 );
}

TEST_CASE( "Addition", "[math]" )
{
	t1();
	REQUIRE( add( 1, 2 ) == 2 );
	REQUIRE( add( -1, -2 ) == -3 );
	REQUIRE( add( -1, 2 ) == 1 );
	
}

TEST_CASE( "Multiplication", "[math]" )
{
	REQUIRE( 1 * 2 == 2 );
	REQUIRE( -1 * -2 == 2 );
	REQUIRE( -1 * 2 == -23 );
}