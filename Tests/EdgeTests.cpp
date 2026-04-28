#include <catch2/catch_test_macros.hpp>

#include "Athena.hpp"
#include "MpfrInclude.hpp"

namespace
{
    using AthenaBinaryFunction = void ( * )( Athena::Number&, const Athena::Number&, const Athena::Number&, Athena::round_t );
    using MpfrBinaryFunction = int ( * )( mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t );

    void requireMatchesMpfr(
        const char* lhs,
        const char* rhs,
        Athena::precision_t precision,
        MpfrBinaryFunction mpfrFunction,
        AthenaBinaryFunction athenaFunction )
    {
        mpfr_t mpfrLhs, mpfrRhs, mpfrExpected;
        mpfr_init2( mpfrLhs, precision );
        mpfr_init2( mpfrRhs, precision );
        mpfr_init2( mpfrExpected, precision );

        mpfr_set_str( mpfrLhs, lhs, 10, MPFR_RNDN );
        mpfr_set_str( mpfrRhs, rhs, 10, MPFR_RNDN );
        mpfrFunction( mpfrExpected, mpfrLhs, mpfrRhs, MPFR_RNDN );

        Athena::Number athenaLhs( lhs, precision );
        Athena::Number athenaRhs( rhs, precision );
        Athena::Number athenaResult( precision );

        athenaFunction( athenaResult, athenaLhs, athenaRhs, MPFR_RNDN );

        REQUIRE( athenaResult == mpfrExpected );

        mpfr_clear( mpfrExpected );
        mpfr_clear( mpfrRhs );
        mpfr_clear( mpfrLhs );
    }
}

TEST_CASE( "Edge cases for addition" )
{
    constexpr Athena::precision_t precision = 256;

    SECTION( "Adding exact opposites returns zero" )
    {
        requireMatchesMpfr( "123456789.123456789", "-123456789.123456789", precision, mpfr_add, Athena::add );
    }

    SECTION( "Adding zero keeps value" )
    {
        requireMatchesMpfr( "0", "-987654321.5", precision, mpfr_add, Athena::add );
        requireMatchesMpfr( "3141592653.589793", "0", precision, mpfr_add, Athena::add );
        requireMatchesMpfr( "0", "0", precision, mpfr_add, Athena::add );
    }

    SECTION( "Very different magnitudes" )
    {
        requireMatchesMpfr( "1e100", "1", precision, mpfr_add, Athena::add );
    }

    SECTION( "Tiny term is ignored at low precision" )
    {
        constexpr Athena::precision_t lowPrecision = 24;
        requireMatchesMpfr( "1e30", "1", lowPrecision, mpfr_add, Athena::add );
        requireMatchesMpfr( "-1e30", "1", lowPrecision, mpfr_add, Athena::add );
        requireMatchesMpfr( "1", "1e30", lowPrecision, mpfr_add, Athena::add );
    }

    SECTION( "Result larger than precision capacity rounds correctly" )
    {
        constexpr Athena::precision_t lowPrecision = 24;
        requireMatchesMpfr( "999999999999999999999999999999", "999999999999999999999999999999", lowPrecision, mpfr_add, Athena::add );
        requireMatchesMpfr( "1.111111111111111111111111111111e100", "8.888888888888888888888888888889e99", lowPrecision, mpfr_add, Athena::add );
    }
}

TEST_CASE( "Edge cases for subtraction" )
{
    constexpr Athena::precision_t precision = 256;

    SECTION( "Subtracting itself returns zero" )
    {
        requireMatchesMpfr( "-9999999.25", "-9999999.25", precision, mpfr_sub, Athena::sub );
    }

    SECTION( "Near cancellation" )
    {
        requireMatchesMpfr( "1.0000000000000001", "1.0", precision, mpfr_sub, Athena::sub );
    }

    SECTION( "Subtracting zero" )
    {
        requireMatchesMpfr( "123456789.5", "0", precision, mpfr_sub, Athena::sub );
        requireMatchesMpfr( "0", "123456789.5", precision, mpfr_sub, Athena::sub );
        requireMatchesMpfr( "0", "0", precision, mpfr_sub, Athena::sub );
    }

    SECTION( "Tiny term is ignored at low precision" )
    {
        constexpr Athena::precision_t lowPrecision = 24;
        requireMatchesMpfr( "1e30", "1", lowPrecision, mpfr_sub, Athena::sub );
        requireMatchesMpfr( "-1e30", "1", lowPrecision, mpfr_sub, Athena::sub );
        requireMatchesMpfr( "1", "1e30", lowPrecision, mpfr_sub, Athena::sub );
    }

    SECTION( "Result larger than precision capacity rounds correctly" )
    {
        constexpr Athena::precision_t lowPrecision = 24;
        requireMatchesMpfr( "-999999999999999999999999999999", "999999999999999999999999999999", lowPrecision, mpfr_sub, Athena::sub );
        requireMatchesMpfr( "1.000000000000000000000000000001e120", "-9.99999999999999999999999999999e119", lowPrecision, mpfr_sub, Athena::sub );
    }
}

