#include "BenchmarkFunction.hpp"

#include <algorithm>
#include <chrono>
#include <limits>
#include <utility>
#include <vector>
#include <cassert>
#include <iostream>
#include <iomanip>

#include "RandomNumber.hpp"

namespace Benchmarking
{
    BenchmarkFunctionStdMath::BenchmarkFunctionStdMath( StdMathFunction a_Function,
                                                        Athena::precision_t a_Precision )
    {
        m_Function = a_Function;
        setPrecision( a_Precision );
        setIterations( 1000 );
    }

    void BenchmarkFunctionStdMath::setPrecisionRange( Athena::precision_t a_StartPrecision, Athena::precision_t a_EndPrecision, Athena::precision_t a_Step )
    {
		assert( a_Step > 0 && "Step must be greater than 0" );
		m_Precision = { a_StartPrecision, a_EndPrecision, a_Step };
    }

    void BenchmarkFunctionStdMath::runMedian()
    {
        m_UseJustMean = false;
        m_Summaries.clear();
        m_Summaries.reserve( (m_Precision.end - m_Precision.start) / m_Precision.step );
        for ( index_t prec = m_Precision.start; prec <= m_Precision.end; prec += m_Precision.step )
        {
            std::cout << "Running benchmark for precision " << prec << "..." << std::endl;
            generateTestData( m_Iterations + m_WarmupIterations, prec );
			BenchmarkSummary summary{ 0.0, 0.0, std::numeric_limits<double>::max(), 0.0, m_Iterations };

            std::vector<double> times;
            times.reserve( m_Iterations );
            double totalTime = 0.0;

            for ( index_t i = 0; i < m_WarmupIterations; i++ )
            {
                m_Function( m_TestResults[i], m_TestDataA[i], m_TestDataB[i], MPFR_RNDN );
			}
            
			for ( index_t i = 0; i < m_Iterations; i++ )
            {
				auto start = std::chrono::high_resolution_clock::now();
                m_Function( m_TestResults[i], m_TestDataA[i], m_TestDataB[i], MPFR_RNDN );
				auto end = std::chrono::high_resolution_clock::now();

                auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>( end - start ).count();
                double timeNs = static_cast<double>( duration );

                times.push_back( timeNs );
                totalTime += timeNs;
				summary.minTimeNanoseconds = std::min( summary.minTimeNanoseconds, timeNs );
                summary.maxTimeNanoseconds = std::max( summary.maxTimeNanoseconds, timeNs );
            }

			summary.meanTimeNanoseconds = totalTime / static_cast<double>(m_Iterations);
            
            // Calculate median
            std::sort( times.begin(), times.end() );
            if ( m_Iterations % 2 == 0 )
            {
                // Even number of elements: average of middle two
                summary.medianTimeNanoseconds = (times[m_Iterations / 2 - 1] + times[m_Iterations / 2]) / 2.0;
            }
            else
            {
                // Odd number of elements: middle element
                summary.medianTimeNanoseconds = times[m_Iterations / 2];
            }
			m_Summaries.emplace_back( summary );
        }

        printSummaries();
    }

    void BenchmarkFunctionStdMath::runMean()
    {
		m_UseJustMean = true;
        m_Summaries.clear();
        m_Summaries.reserve( (m_Precision.end - m_Precision.start) / m_Precision.step );
        for ( index_t prec = m_Precision.start; prec <= m_Precision.end; prec += m_Precision.step )
        {
            std::cout << "Running benchmark for precision " << prec << "..." << std::endl;
            generateTestData( m_Iterations + m_WarmupIterations, prec );
            BenchmarkSummary summary{ 0.0, 0.0, 0.0, 0.0, m_Iterations };

            double totalTime = 0.0;

			for ( index_t i = 0; i < m_WarmupIterations; i++ )
            {
                m_Function( m_TestResults[i], m_TestDataA[i], m_TestDataB[i], MPFR_RNDN );
            }

            auto start = std::chrono::high_resolution_clock::now();
            for ( index_t i = 0; i < m_Iterations; i++ )
            {
                m_Function( m_TestResults[i], m_TestDataA[i], m_TestDataB[i], MPFR_RNDN );
            }
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            double timeNs = static_cast<double>(duration);

            summary.meanTimeNanoseconds = timeNs / static_cast<double>(m_Iterations);
            m_Summaries.emplace_back( summary );
        }

        printSummaries();
    }

    void BenchmarkFunctionStdMath::runForProfiler()
    {
		generateTestData( 1, m_Precision.start );

        for ( index_t i = 0; i < m_Iterations; ++i )
        {
            m_Function( m_TestResults[0], m_TestDataA[0], m_TestDataB[0], MPFR_RNDN );
		}
    }

    void BenchmarkFunctionStdMath::printSummaries() const
    {
        if ( m_Summaries.empty() )
        {
            std::cout << "No benchmark summaries available. Run benchmark first." << std::endl;
            return;
        }

        constexpr int precisionWidth = 11;
        constexpr int iterationsWidth = 12;
        constexpr int timeWidth = 14;

        std::cout << std::endl;
        std::cout << "===== Benchmarking  results =====" << std::endl;
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

        if ( m_UseJustMean )
        {
            const auto printMeanSeparator = []()
            {
                std::cout << '+' << std::string( precisionWidth + 2, '-' )
                          << '+' << std::string( timeWidth + 2, '-' )
                          << '+' << std::endl;
            };

            printMeanSeparator();
            std::cout << "| " << std::setw( precisionWidth ) << std::left << "Precision"
                      << " | " << std::setw( timeWidth ) << std::left << "Mean (ns)"
                      << " |" << std::endl;
            printMeanSeparator();

            const std::size_t count = m_Summaries.size();
            for ( std::size_t i = 0; i < count; ++i )
            {
                const Athena::precision_t precision = m_Precision.start + static_cast<Athena::precision_t>( i ) * m_Precision.step;
                const BenchmarkSummary& summary = m_Summaries[i];

                std::cout << "| " << std::setw( precisionWidth ) << std::left << precision
                          << " | " << std::setw( timeWidth ) << std::left << std::fixed << std::setprecision( 2 ) << summary.meanTimeNanoseconds
                          << " |" << std::endl;
            }

            printMeanSeparator();
        }
        else
        {
            const auto printSeparator = []()
                {
                    std::cout << '+' << std::string( precisionWidth + 2, '-' )
                        << '+' << std::string( iterationsWidth + 2, '-' )
                        << '+' << std::string( timeWidth + 2, '-' )
                        << '+' << std::string( timeWidth + 2, '-' )
                        << '+' << std::string( timeWidth + 2, '-' )
                        << '+' << std::string( timeWidth + 2, '-' )
                        << '+' << std::endl;
                };

            printSeparator();
            std::cout << "| " << std::setw( precisionWidth ) << std::left << "Precision"
                << " | " << std::setw( iterationsWidth ) << std::left << "Iterations"
                << " | " << std::setw( timeWidth ) << std::left << "Mean (ns)"
                << " | " << std::setw( timeWidth ) << std::left << "Median (ns)"
                << " | " << std::setw( timeWidth ) << std::left << "Min (ns)"
                << " | " << std::setw( timeWidth ) << std::left << "Max (ns)"
                << " |" << std::endl;
            printSeparator();

            const std::size_t count = m_Summaries.size();
            for ( std::size_t i = 0; i < count; ++i )
            {
                const Athena::precision_t precision = m_Precision.start + static_cast<Athena::precision_t>( i ) * m_Precision.step;
                const BenchmarkSummary& summary = m_Summaries[i];

                std::cout << "| " << std::setw( precisionWidth ) << std::left << precision
                    << " | " << std::setw( iterationsWidth ) << std::left << summary.iterations
                    << " | " << std::setw( timeWidth ) << std::left << std::fixed << std::setprecision( 2 ) << summary.meanTimeNanoseconds
                    << " | " << std::setw( timeWidth ) << std::left << std::fixed << std::setprecision( 2 ) << summary.medianTimeNanoseconds
                    << " | " << std::setw( timeWidth ) << std::left << std::fixed << std::setprecision( 2 ) << summary.minTimeNanoseconds
                    << " | " << std::setw( timeWidth ) << std::left << std::fixed << std::setprecision( 2 ) << summary.maxTimeNanoseconds
                    << " |" << std::endl;
            }

            printSeparator();
        }
    }

	void BenchmarkFunctionStdMath::generateTestData( std::size_t a_Number, Athena::precision_t a_Precision )
    {
		m_TestDataA.resize( a_Number );
		m_TestDataB.resize( a_Number );
		m_TestResults.resize( a_Number );

		for ( index_t i = 0; i < a_Number; ++i )
        {
            m_TestDataA[i] = generateRandomNumber( a_Precision, true );
            m_TestDataB[i] = generateRandomNumber( a_Precision, true );
			m_TestResults[i] = Athena::Number( a_Precision );
        }
    }
}
