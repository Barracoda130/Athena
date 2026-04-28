#include "GoogleBenchmark.hpp"

#include "MpfrInclude.hpp"

#include <array>
#include <cstring>
#include <random>
#include <vector>

#include <benchmark/benchmark.h>

#include "Athena.hpp"
#include "RandomNumber.hpp"

#define POOL_SIZE 16384

namespace
{
    using AthenaBinaryFunction = void ( * )( Athena::Number&, const Athena::Number&, const Athena::Number&, Athena::round_t );
    using MpfrBinaryFunction = int ( * )( mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t );

	// Benchmark functions for Athena
    void BM_Athena( benchmark::State& state, AthenaBinaryFunction a_Function )
    {
		// Get the precision from the benchmark state
        const auto precision = static_cast<Athena::precision_t>( state.range( 0 ) );

        std::vector<Athena::Number> lhs;
        std::vector<Athena::Number> rhs;
        std::vector<Athena::Number> out;

        lhs.reserve( POOL_SIZE );
        rhs.reserve( POOL_SIZE );
        out.reserve( POOL_SIZE );

		// Generate random numbers for the benchmark
        for ( std::size_t i = 0; i < POOL_SIZE; ++i )
        {
            lhs.emplace_back( generateRandomNumber( static_cast<std::size_t>( precision ), true ) );
            rhs.emplace_back( generateRandomNumber( static_cast<std::size_t>( precision ), true ) );
            out.emplace_back( precision );
        }

		// Run the benchmark loop
        std::size_t idx = 0;
        for ( auto _ : state )
        {
            a_Function( out[idx], lhs[idx], rhs[idx], MPFR_RNDN );
            benchmark::DoNotOptimize( out[idx] );   // Make sure the compiler does not optimise away the result

            idx = ( idx + 1 ) % POOL_SIZE;
        }
    }

	// Benchmark functions for MPFR
    void BM_mpfr( benchmark::State& state, MpfrBinaryFunction a_Function )
    {
        using MpfrStorage = std::array<__mpfr_struct, 1>;

		// Get the precision from the benchmark state
        const auto precision = static_cast<mpfr_prec_t>( state.range( 0 ) );

		// Estimate the number of digits needed for the given precision (log10(2) ≈ 0.30103)
        const std::size_t numDigits = ( static_cast<std::size_t>( precision ) * 3 ) / 10 + 1;

        std::vector<MpfrStorage> lhs( POOL_SIZE );
        std::vector<MpfrStorage> rhs( POOL_SIZE );
        std::vector<MpfrStorage> result( POOL_SIZE );

        std::vector<char> lhsBuffer( numDigits + 2 );
        std::vector<char> rhsBuffer( numDigits + 2 );

		// Generate random numbers for the benchmark
        for ( std::size_t i = 0; i < POOL_SIZE; ++i )
        {
            mpfr_init2( lhs[i].data(), precision );
            mpfr_init2( rhs[i].data(), precision );
            mpfr_init2( result[i].data(), precision );

            generateFloat( numDigits, lhsBuffer.data() );
            generateFloat( numDigits, rhsBuffer.data() );

            mpfr_set_str( lhs[i].data(), lhsBuffer.data(), 10, MPFR_RNDN );
            mpfr_set_str( rhs[i].data(), rhsBuffer.data(), 10, MPFR_RNDN );
        }

		// Run the benchmark loop
        std::size_t idx = 0;
        for ( auto _ : state )
        {
            a_Function( result[idx].data(), lhs[idx].data(), rhs[idx].data(), MPFR_RNDN );
			benchmark::DoNotOptimize( result[idx].data() );     // Ensure the compiler does not optimise away the result

            idx = ( idx + 1 ) % POOL_SIZE;
        }

		// Clean up MPFR variables
        for ( std::size_t i = 0; i < POOL_SIZE; ++i )
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

	// Benchmark for Athena's atn_add_n function, which operates on raw limb arrays
    void BM_atnAddN( benchmark::State& state )
    {
        const auto limbCount = static_cast<mp_size_t>( state.range( 0 ) );

        std::vector<mp_limb_t> lhs( POOL_SIZE * limbCount );
        std::vector<mp_limb_t> rhs( POOL_SIZE * limbCount );
        std::vector<mp_limb_t> out( POOL_SIZE * limbCount );

        std::mt19937_64 rng( 0xC0FFEEULL );
        std::uniform_int_distribution<std::uint64_t> dist;

		// Generate random limb data for the benchmark
        for ( std::size_t i = 0; i < lhs.size(); ++i )
        {
            lhs[i] = static_cast<mp_limb_t>( dist( rng ) );
            rhs[i] = static_cast<mp_limb_t>( dist( rng ) );
        }

		// Run the benchmark loop
        std::size_t idx = 0;
        for ( auto _ : state )
        {
            const std::size_t offset = idx * static_cast<std::size_t>( limbCount );
            const auto carry = Athena::atn_add_n( out.data() + offset, lhs.data() + offset, rhs.data() + offset, limbCount );

            // Ensure result is not optimised away by compiler
            benchmark::DoNotOptimize( carry );
            benchmark::DoNotOptimize( out.data() + offset );

            idx = ( idx + 1 ) % POOL_SIZE;
        }

        state.SetItemsProcessed( state.iterations() );
    }

    BENCHMARK( BM_AthenaAdd )
        ->Arg( 128 )
        ->Arg( 256 )
        ->Arg( 384 )
        ->Arg( 512 )
        ->Arg( 640 )
        ->Arg( 768 )
        ->Arg( 896 )
        ->Arg( 1024 );

    BENCHMARK( BM_mpfrAdd )
        ->Arg( 128 )
        ->Arg( 256 )
        ->Arg( 384 )
        ->Arg( 512 )
        ->Arg( 640 )
        ->Arg( 768 )
        ->Arg( 896 )
        ->Arg( 1024 );

    BENCHMARK( BM_AthenaSub )
        ->Arg( 128 )
        ->Arg( 256 )
        ->Arg( 384 )
        ->Arg( 512 )
        ->Arg( 640 )
        ->Arg( 768 )
        ->Arg( 896 )
        ->Arg( 1024 );

    BENCHMARK( BM_mpfrSub )
        ->Arg( 128 )
        ->Arg( 256 )
        ->Arg( 384 )
        ->Arg( 512 )
        ->Arg( 640 )
        ->Arg( 768 )
        ->Arg( 896 )
        ->Arg( 1024 );

    BENCHMARK( BM_atnAddN )
        ->Arg( 128 )
        ->Arg( 256 )
        ->Arg( 384 )
        ->Arg( 512 )
        ->Arg( 640 )
        ->Arg( 768 )
        ->Arg( 896 )
        ->Arg( 1024 );
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

