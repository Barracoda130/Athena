#pragma once

#include <any>
#include <cstddef>
#include <functional>

#include "Athena.hpp"

namespace Benchmarking
{
	struct BenchmarkSummary
	{
		double meanTimeNanoseconds;
		double medianTimeNanoseconds;
		double minTimeNanoseconds;
		double maxTimeNanoseconds;
		std::size_t iterations;
	};

	struct PrecisionRange
	{
		Athena::precision_t start;
		Athena::precision_t end;
		Athena::precision_t step;
	};

	using StdMathFunction = void ( * )( Athena::Number&, const Athena::Number&, const Athena::Number&, Athena::round_t );

	class BenchmarkFunctionStdMath
	{
	public:
		BenchmarkFunctionStdMath( StdMathFunction a_Function,
								  Athena::precision_t a_Precision = mpfr_get_default_prec() );

		void setPrecision( Athena::precision_t a_Precision ) { m_Precision = { a_Precision, a_Precision, 1 }; }
		void setPrecisionRange( Athena::precision_t a_StartPrecision, Athena::precision_t a_EndPrecision, Athena::precision_t a_Step = 1 );
		void setIterations( std::size_t a_Iterations ) { m_Iterations = a_Iterations; }
		void setWarmupIterations( std::size_t a_WarmupIterations ) { m_WarmupIterations = a_WarmupIterations; }

		void runMedian();
		void runMean();

		void runForProfiler();

		void printSummaries() const;

	private:
		void generateTestData( std::size_t a_Number, Athena::precision_t a_Precision );

		StdMathFunction m_Function;
		std::vector<Athena::Number> m_TestDataA;
		std::vector<Athena::Number> m_TestDataB;
		std::vector<Athena::Number> m_TestResults;
		PrecisionRange m_Precision;
		std::size_t m_Iterations;
		std::vector<BenchmarkSummary> m_Summaries;
		bool m_UseJustMean = false;
		std::size_t m_WarmupIterations = 1000;
	};
}


