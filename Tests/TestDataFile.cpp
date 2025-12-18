#include "TestDataFile.hpp"

#include <array>

static const char s_FileName[] = "TestData.csv";
static std::array<const char[7], 12> s_ColumnHeaders = {
	"posInt", "posInt",
	"negInt", "negInt",
	"rndInt", "rndInt",
	"posFlt", "posFlt",
	"negFlt", "negFlt",
	"rndFlt", "rndFlt"
};

namespace TestData
{
	const char* getFileName()
	{
		return s_FileName;
	}

	void generateData( std::size_t aNumRows )
	{
		// Generate aNumRows rows of test data
		
	}
}



