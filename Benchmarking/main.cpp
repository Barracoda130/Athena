#include "BenchmarkFunction.hpp"
#include "Athena.hpp"

#include <iostream>

struct AddBenchmarkData
{
	Athena::Number lhs;
	Athena::Number rhs;
	Athena::Number result;
	Athena::round_t round;
};

int main()
{
	constexpr Athena::precision_t precision = 256;
	constexpr std::size_t iterations = 1000000;

	Benchmarking::BenchmarkFunctionStdMath benchmark( Athena::add );
	benchmark.setPrecisionRange( 100, 1000, 100 );
	benchmark.setIterations( iterations );

	//benchmark.runMean();
	//benchmark.runMedian();
	benchmark.runForProfiler();

	return 0;
}

