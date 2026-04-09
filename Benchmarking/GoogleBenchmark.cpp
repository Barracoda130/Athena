#include "GoogleBenchmark.hpp"

#include "MpfrInclude.hpp"

#include <array>
#include <cstring>
#include <vector>


#include <benchmark/benchmark.h>

#include "Athena.hpp"
#include "RandomNumber.hpp"

namespace
{
    using AthenaBinaryFunction = void ( * )( Athena::Number&, const Athena::Number&, const Athena::Number&, Athena::round_t );
    using MpfrBinaryFunction = int ( * )( mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t );

    void BM_Athena( benchmark::State& state, AthenaBinaryFunction a_Function )
    {
        const auto precision = static_cast<Athena::precision_t>( state.range( 0 ) );
        constexpr std::size_t poolSize = 16384;

        std::vector<Athena::Number> lhs;
        std::vector<Athena::Number> rhs;
        std::vector<Athena::Number> out;

        lhs.reserve( poolSize );
        rhs.reserve( poolSize );
        out.reserve( poolSize );

        for ( std::size_t i = 0; i < poolSize; ++i )
        {
            lhs.emplace_back( generateRandomNumber( static_cast<std::size_t>( precision ), true ) );
            rhs.emplace_back( generateRandomNumber( static_cast<std::size_t>( precision ), true ) );
            out.emplace_back( precision );
        }

        std::size_t idx = 0;
        for ( auto _ : state )
        {
            a_Function( out[idx], lhs[idx], rhs[idx], MPFR_RNDN );
            benchmark::DoNotOptimize( out[idx] );

            idx = ( idx + 1 ) % poolSize;
        }
    }

    void BM_mpfr( benchmark::State& state, MpfrBinaryFunction a_Function )
    {
        using MpfrStorage = std::array<__mpfr_struct, 1>;

        const auto precision = static_cast<mpfr_prec_t>( state.range( 0 ) );
        constexpr std::size_t poolSize = 16384;
        const std::size_t numDigits = ( static_cast<std::size_t>( precision ) * 3 ) / 10 + 1;

        std::vector<MpfrStorage> lhs( poolSize );
        std::vector<MpfrStorage> rhs( poolSize );
        std::vector<MpfrStorage> result( poolSize );

        std::vector<char> lhsBuffer( numDigits + 2 );
        std::vector<char> rhsBuffer( numDigits + 2 );

        for ( std::size_t i = 0; i < poolSize; ++i )
        {
            mpfr_init2( lhs[i].data(), precision );
            mpfr_init2( rhs[i].data(), precision );
            mpfr_init2( result[i].data(), precision );

            generateFloat( numDigits, lhsBuffer.data() );
            generateFloat( numDigits, rhsBuffer.data() );

            mpfr_set_str( lhs[i].data(), lhsBuffer.data(), 10, MPFR_RNDN );
            mpfr_set_str( rhs[i].data(), rhsBuffer.data(), 10, MPFR_RNDN );
        }

        std::size_t idx = 0;
        for ( auto _ : state )
        {
            a_Function( result[idx].data(), lhs[idx].data(), rhs[idx].data(), MPFR_RNDN );
            benchmark::DoNotOptimize( result[idx].data() );

            idx = ( idx + 1 ) % poolSize;
        }

        for ( std::size_t i = 0; i < poolSize; ++i )
        {
            mpfr_clear( lhs[i].data() );
            mpfr_clear( rhs[i].data() );
            mpfr_clear( result[i].data() );
        }
    }

    void BM_mpfrAdd( benchmark::State& state )
    {
        BM_mpfr( state, mpfr_add );
    }

    void BM_mpfrSub( benchmark::State& state )
    {
        BM_mpfr( state, mpfr_sub );
    }

    void BM_AthenaAdd( benchmark::State& state )
    {
        BM_Athena( state, Athena::add );
    }

    void BM_AthenaSub( benchmark::State& state )
    {
        BM_Athena( state, Athena::sub );
    }

    BENCHMARK( BM_AthenaAdd )
        ->Arg( 100 )
        ->Arg( 200 )
        ->Arg( 300 )
        ->Arg( 400 )
        ->Arg( 500 )
        ->Arg( 600 )
        ->Arg( 700 )
        ->Arg( 800 )
        ->Arg( 900 )
        ->Arg( 1000 );

    BENCHMARK( BM_mpfrAdd )
        ->Arg( 100 )
        ->Arg( 200 )
        ->Arg( 300 )
        ->Arg( 400 )
        ->Arg( 500 )
        ->Arg( 600 )
        ->Arg( 700 )
        ->Arg( 800 )
        ->Arg( 900 )
        ->Arg( 1000 );

    BENCHMARK( BM_AthenaSub )
        ->Arg( 100 )
        ->Arg( 200 )
        ->Arg( 300 )
        ->Arg( 400 )
        ->Arg( 500 )
        ->Arg( 600 )
        ->Arg( 700 )
        ->Arg( 800 )
        ->Arg( 900 )
        ->Arg( 1000 );

    BENCHMARK( BM_mpfrSub )
        ->Arg( 100 )
        ->Arg( 200 )
        ->Arg( 300 )
        ->Arg( 400 )
        ->Arg( 500 )
        ->Arg( 600 )
        ->Arg( 700 )
        ->Arg( 800 )
        ->Arg( 900 )
        ->Arg( 1000 );
}

namespace Benchmarking
{
    int runGoogleBenchmark( int argc, char** argv )
    {
        benchmark::Initialize( &argc, argv );
        if ( benchmark::ReportUnrecognizedArguments( argc, argv ) )
        {
            return 1;
        }

        benchmark::RunSpecifiedBenchmarks();
        return 0;
    }
}

