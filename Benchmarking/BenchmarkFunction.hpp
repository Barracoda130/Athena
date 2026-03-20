#pragma once

#include <any>
#include <cstddef>
#include <functional>

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

	class BenchmarkFunction
	{
	public:
        using DataItem = std::any;
		using BenchmarkedFunction = std::function<void( DataItem& )>;
		using DataGenerator = std::function<DataItem( std::size_t iterationIndex )>;
		using FunctionInvoker = std::function<void( BenchmarkedFunction&, DataItem& )>;

		BenchmarkFunction( BenchmarkedFunction function,
						   DataGenerator generator,
						   FunctionInvoker invoker );

		BenchmarkSummary run( std::size_t iterations );

		BenchmarkedFunction& function();

	private:
        BenchmarkedFunction m_Function;
		DataGenerator m_Generator;
		FunctionInvoker m_Invoker;
	};
}
