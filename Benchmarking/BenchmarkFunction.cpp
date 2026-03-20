#include "BenchmarkFunction.hpp"

#include <algorithm>
#include <chrono>
#include <limits>
#include <utility>
#include <vector>

namespace Benchmarking
{
    BenchmarkFunction::BenchmarkFunction( BenchmarkedFunction function,
                                          DataGenerator generator,
                                          FunctionInvoker invoker )
        : m_Function( std::move( function ) )
        , m_Generator( std::move( generator ) )
        , m_Invoker( std::move( invoker ) )
    {
    }

    BenchmarkSummary BenchmarkFunction::run( std::size_t iterations )
    {
        if ( iterations == 0 )
        {
            return BenchmarkSummary{ 0.0, 0.0, 0.0, 0.0, 0 };
        }

        std::vector<DataItem> testData;
        testData.reserve( iterations );

        for ( std::size_t i = 0; i < iterations; ++i )
        {
            testData.push_back( m_Generator( i ) );
        }

        std::vector<double> times;
        times.reserve( iterations );

        BenchmarkSummary summary{
            0.0,
            0.0,
            std::numeric_limits<double>::max(),
            0.0,
            iterations
        };

        double totalTime = 0.0;

        for ( std::size_t i = 0; i < iterations; ++i )
        {
            auto start = std::chrono::high_resolution_clock::now();
            m_Invoker( m_Function, testData[i] );
            auto end = std::chrono::high_resolution_clock::now();

            double elapsedNs = static_cast<double>(
                std::chrono::duration_cast<std::chrono::nanoseconds>( end - start ).count() );

            times.push_back( elapsedNs );
            totalTime += elapsedNs;
            summary.minTimeNanoseconds = std::min( summary.minTimeNanoseconds, elapsedNs );
            summary.maxTimeNanoseconds = std::max( summary.maxTimeNanoseconds, elapsedNs );
        }

        summary.meanTimeNanoseconds = totalTime / static_cast<double>( iterations );

        std::sort( times.begin(), times.end() );
        if ( iterations % 2 == 0 )
        {
            summary.medianTimeNanoseconds =
                (times[iterations / 2 - 1] + times[iterations / 2]) / 2.0;
        }
        else
        {
            summary.medianTimeNanoseconds = times[iterations / 2];
        }

        return summary;
    }

    BenchmarkFunction::BenchmarkedFunction& BenchmarkFunction::function()
    {
        return m_Function;
    }
}
