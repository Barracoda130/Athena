#include "Benchmark.hpp"

#include <iostream>
#include <vector>

namespace Benchmarking
{
    struct BenchmarkResult
    {
        double meanTimeNanoseconds;
        double medianTimeNanoseconds;
        double minTimeNanoseconds;
        double maxTimeNanoseconds;
        size_t iterations;
    };

    template<typename Func>
    BenchmarkResult benchmark_function( Func func, size_t iterations )
    {
        BenchmarkResult result{ 0.0, 0.0, std::numeric_limits<double>::max(), 0.0, iterations };
        std::vector<double> times;
        times.reserve( iterations );
        double totalTime = 0.0;

        for ( size_t i = 0; i < iterations; ++i )
        {
            auto start = std::chrono::high_resolution_clock::now();
            func();
            auto end = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>( end - start ).count();
            double timeNs = static_cast<double>( duration );

            times.push_back( timeNs );
            totalTime += timeNs;
            result.minTimeNanoseconds = std::min( result.minTimeNanoseconds, timeNs );
            result.maxTimeNanoseconds = std::max( result.maxTimeNanoseconds, timeNs );
        }

        result.meanTimeNanoseconds = totalTime / static_cast<double>(iterations);

        // Calculate median
        std::sort( times.begin(), times.end() );
        if ( iterations % 2 == 0 )
        {
            // Even number of elements: average of middle two
            result.medianTimeNanoseconds = (times[iterations / 2 - 1] + times[iterations / 2]) / 2.0;
        }
        else
        {
            // Odd number of elements: middle element
            result.medianTimeNanoseconds = times[iterations / 2];
        }

        return result;
    }

    void benchmark_add_n_implementations( size_t limbCount, size_t iterations )
    {
        std::random_device rd;
        std::mt19937_64 gen( rd() );
        std::uniform_int_distribution<mp_limb_t> dist;

        std::cout << "=== Benchmarking atn_add_n implementations ===" << std::endl;
        std::cout << "Generating " << iterations << " sets of random test data..." << std::endl;

        // Pre-generate all test data to avoid timing contamination
        std::vector<std::vector<mp_limb_t>> test_data_a( iterations );
        std::vector<std::vector<mp_limb_t>> test_data_b( iterations );
        std::vector<std::vector<mp_limb_t>> results( iterations );

        for ( size_t iter = 0; iter < iterations; ++iter )
        {
            test_data_a[iter].resize( limbCount );
            test_data_b[iter].resize( limbCount );
            results[iter].resize( limbCount );

            for ( size_t i = 0; i < limbCount; ++i )
            {
                test_data_a[iter][i] = dist( gen );
                test_data_b[iter][i] = dist( gen );
            }
        }

        std::cout << "Data generation complete." << std::endl;
        std::cout << "Limb count: " << limbCount << std::endl;
        std::cout << "Iterations: " << iterations << std::endl;
        std::cout << "Limb size: " << GMP_NUMB_BITS << " bits" << std::endl;
        std::cout << "Total data size: " << (2 * iterations * limbCount * sizeof( mp_limb_t )) / (1024.0 * 1024.0) << " MB" << std::endl;

        // Show build configuration
#ifdef NDEBUG
        std::cout << "Build: RELEASE (optimized)" << std::endl;
#else
        std::cout << "Build: DEBUG (not optimized)" << std::endl;
#endif

#if defined(_MSC_VER)
        std::cout << "Compiler: MSVC " << _MSC_VER << std::endl;
#elif defined(__GNUC__)
        std::cout << "Compiler: GCC " << __GNUC__ << "." << __GNUC_MINOR__ << std::endl;
#elif defined(__clang__)
        std::cout << "Compiler: Clang " << __clang_major__ << "." << __clang_minor__ << std::endl;
#endif
        std::cout << std::endl;

        // Benchmark standard implementation with different data each iteration
        size_t iter_index = 0;
        auto result_standard = benchmark_function( [&]()
            {
                atn_add_n( results[iter_index].data(),
                    test_data_a[iter_index].data(),
                    test_data_b[iter_index].data(),
                    limbCount );
                iter_index = (iter_index + 1) % iterations;
            }, iterations );

        std::cout << "atn_add_n (standard):" << std::endl;
        std::cout << "  Mean:    " << result_standard.meanTimeNanoseconds << " ns" << std::endl;
        std::cout << "  Median:  " << result_standard.medianTimeNanoseconds << " ns" << std::endl;
        std::cout << "  Min:     " << result_standard.minTimeNanoseconds << " ns" << std::endl;
        std::cout << "  Max:     " << result_standard.maxTimeNanoseconds << " ns" << std::endl;
        std::cout << std::endl;

        // Benchmark intrinsic implementation with different data each iteration
        iter_index = 0;
        auto result_intrinsic = benchmark_function( [&]()
            {
                atn_add_n_intrinsic( results[iter_index].data(),
                    test_data_a[iter_index].data(),
                    test_data_b[iter_index].data(),
                    limbCount );
                iter_index = (iter_index + 1) % iterations;
            }, iterations );

        std::cout << "atn_add_n_intrinsic:" << std::endl;
        std::cout << "  Mean:    " << result_intrinsic.meanTimeNanoseconds << " ns" << std::endl;
        std::cout << "  Median:  " << result_intrinsic.medianTimeNanoseconds << " ns" << std::endl;
        std::cout << "  Min:     " << result_intrinsic.minTimeNanoseconds << " ns" << std::endl;
        std::cout << "  Max:     " << result_intrinsic.maxTimeNanoseconds << " ns" << std::endl;
        std::cout << std::endl;

        iter_index = 0;
        auto result_intrinsic_adx = benchmark_function( [&]()
            {
                atn_add_n_intrinsic_adx( results[iter_index].data(),
                    test_data_a[iter_index].data(),
                    test_data_b[iter_index].data(),
                    limbCount );
                iter_index = (iter_index + 1) % iterations;
            }, iterations );

        std::cout << "atn_add_n_intrinsic_adx:" << std::endl;
        std::cout << "  Mean:    " << result_intrinsic_adx.meanTimeNanoseconds << " ns" << std::endl;
        std::cout << "  Median:  " << result_intrinsic_adx.medianTimeNanoseconds << " ns" << std::endl;
        std::cout << "  Min:     " << result_intrinsic_adx.minTimeNanoseconds << " ns" << std::endl;
        std::cout << "  Max:     " << result_intrinsic_adx.maxTimeNanoseconds << " ns" << std::endl;
        std::cout << std::endl;

        double speedup = result_standard.medianTimeNanoseconds / result_intrinsic.medianTimeNanoseconds;


        std::cout << "Performance comparison:" << std::endl;
        if ( speedup > 1.0 )
        {
            std::cout << "  Intrinsic version is " << std::fixed << std::setprecision( 2 ) << speedup << "x faster" << std::endl;
            std::cout << "  Speedup: " << std::fixed << std::setprecision( 1 ) << ((speedup - 1.0) * 100.0) << "%" << std::endl;
        }
        else
        {
            std::cout << "  Standard version is " << std::fixed << std::setprecision( 2 ) << (1.0 / speedup) << "x faster" << std::endl;
            std::cout << "  Slowdown: " << std::fixed << std::setprecision( 1 ) << ((1.0 - speedup) * 100.0) << "%" << std::endl;
        }
        std::cout << std::endl;

        // Verify correctness using first test set
        std::vector<mp_limb_t> verify_result1( limbCount );
        std::vector<mp_limb_t> verify_result2( limbCount );
        std::vector<mp_limb_t> verify_result3( limbCount );

        mp_limb_t carry1 = atn_add_n( verify_result1.data(), test_data_a[0].data(), test_data_b[0].data(), limbCount );
        mp_limb_t carry2 = atn_add_n_intrinsic( verify_result2.data(), test_data_a[0].data(), test_data_b[0].data(), limbCount );
        mp_limb_t carry3 = atn_add_n_intrinsic_adx( verify_result3.data(), test_data_a[0].data(), test_data_b[0].data(), limbCount );

        bool results_match = (carry1 == carry2) &&
            (carry1 == carry3) &&
            std::equal( verify_result1.begin(), verify_result1.end(), verify_result2.begin() ) &&
            std::equal( verify_result1.begin(), verify_result1.end(), verify_result3.begin() );

        if ( results_match )
        {
            std::cout << "Correctness: Both implementations produce identical results." << std::endl;
        }
        else
        {
            std::cout << "WARNING: Implementations produce different results!" << std::endl;
            std::cout << "  Standard carry: " << carry1 << std::endl;
            std::cout << "  Intrinsic carry: " << carry2 << std::endl;
            std::cout << "  Intrinsic ADX carry: " << carry3 << std::endl;
            if ( carry1 == carry2 && carry1 == carry3 )
            {
                std::cout << "  Carry values match, but result limbs differ!" << std::endl;
            }
        }
    }

	void RunBenchmark()
	{
		
	}
}