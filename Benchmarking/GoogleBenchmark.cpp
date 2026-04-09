#include "GoogleBenchmark.hpp"

#include "MpfrInclude.hpp"

#include <cstring>
#include <vector>


#include <benchmark/benchmark.h>

#include "Athena.hpp"

namespace
{
    void BM_AthenaAdd( benchmark::State& state )
    {
        const auto precision = static_cast<Athena::precision_t>( state.range( 0 ) );
        Athena::Number lhs( "12345.67890123456789", precision );
        Athena::Number rhs( "98765.43210987654321", precision );
        Athena::Number result( precision );

        for ( auto _ : state )
        {
            Athena::add( result, lhs, rhs, MPFR_RNDN );
            benchmark::DoNotOptimize( result );
        }
    }

    void BM_mpfrAdd( benchmark::State& state )
    {
        const auto precision = static_cast<mpfr_prec_t>(state.range( 0 ));
		mpfr_t lhs, rhs, result;
		mpfr_init2( lhs, precision );
		mpfr_init2( rhs, precision );
		mpfr_init2( result, precision );
		mpfr_set_str( lhs, "12345.67890123456789", 10, MPFR_RNDN );
		mpfr_set_str( rhs, "98765.43210987654321", 10, MPFR_RNDN );

        for ( auto _ : state )
        {
            mpfr_add( result, lhs, rhs, MPFR_RNDN );
            benchmark::DoNotOptimize( result );
        }
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
}

namespace Benchmarking
{
    int runGoogleBenchmark( int argc, char** argv )
    {
        std::vector<char*> filteredArgs;
        filteredArgs.reserve( static_cast<std::size_t>( argc ) );

        for ( int i = 0; i < argc; ++i )
        {
            if ( std::strcmp( argv[i], "--google-benchmark" ) != 0 )
            {
                filteredArgs.push_back( argv[i] );
            }
        }

        int filteredArgc = static_cast<int>( filteredArgs.size() );
        char** filteredArgv = filteredArgs.data();

        benchmark::Initialize( &filteredArgc, filteredArgv );
        if ( benchmark::ReportUnrecognizedArguments( filteredArgc, filteredArgv ) )
        {
            return 1;
        }

        benchmark::RunSpecifiedBenchmarks();
        return 0;
    }
}

