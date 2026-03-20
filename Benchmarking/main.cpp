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
	constexpr std::size_t iterations = 10000;

	Benchmarking::BenchmarkFunction benchmark(
		[]( Benchmarking::BenchmarkFunction::DataItem& dataItem )
		{
			auto& data = std::any_cast<AddBenchmarkData&>( dataItem );
			Athena::add( data.result, data.lhs, data.rhs, data.round );
		},
		[&]( std::size_t iterationIndex ) -> Benchmarking::BenchmarkFunction::DataItem
		{
			long long a = static_cast<long long>( iterationIndex % 100000 ) - 50000;
			long long b = static_cast<long long>( (iterationIndex * 37) % 100000 ) - 50000;

			return AddBenchmarkData{
				Athena::Number( a, precision ),
				Athena::Number( b, precision ),
				Athena::Number( precision ),
				MPFR_RNDN
			};
		},
		[]( Benchmarking::BenchmarkFunction::BenchmarkedFunction& function,
			Benchmarking::BenchmarkFunction::DataItem& dataItem )
		{
			function( dataItem );
		} );

	auto result = benchmark.run( iterations );

	std::cout << "Athena::add benchmark" << std::endl;
	std::cout << "Iterations: " << result.iterations << std::endl;
	std::cout << "Mean:       " << result.meanTimeNanoseconds << " ns" << std::endl;
	std::cout << "Median:     " << result.medianTimeNanoseconds << " ns" << std::endl;
	std::cout << "Min:        " << result.minTimeNanoseconds << " ns" << std::endl;
	std::cout << "Max:        " << result.maxTimeNanoseconds << " ns" << std::endl;

	return 0;
}

