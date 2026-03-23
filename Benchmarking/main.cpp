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
	constexpr std::size_t iterations = 10000;

	Benchmarking::BenchmarkFunctionStdMath benchmark( Athena::add );
	benchmark.setInfoFunction( Athena::printAtn_add_nVersion );
	benchmark.setPrecisionRange( 100, 1000, 100 );
	benchmark.setIterations( iterations );

	benchmark.runMean();
	benchmark.runMedian();
	//benchmark.runForProfiler();

	return 0;
}

