#include "BenchmarkFunction.hpp"
#include "GoogleBenchmark.hpp"
#include "Athena.hpp"

#include <iostream>
#include <cstring>

struct AddBenchmarkData
{
	Athena::Number lhs;
	Athena::Number rhs;
	Athena::Number result;
	Athena::round_t round;
};

int main( int argc, char** argv )
{
   for ( int i = 1; i < argc; ++i )
	{
		if ( std::strcmp( argv[i], "--google-benchmark" ) == 0 )
		{
			
		}
	}
   return Benchmarking::runGoogleBenchmark( argc, argv );
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

