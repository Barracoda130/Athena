void generateInt( std::size_t a_Length, char* result )
{
	for ( index_t i = 0; i < a_Length; i++ )
	{
		result[i] = static_cast<char>( randInt( 0, 9 ) ) + '0';	// Get random number between 0 and 10
	}
	result[a_Length] = '\0';
}

// Generate a float between a_Min and a_Max
void generateFloat( std::size_t a_Length, char* result )
{
	int decimalPointIdx = randInt( 0, static_cast<int>( a_Length ) - 2 );
	index_t startIdx = 0;

	if ( decimalPointIdx == 0 )
	{
		result[0] = '0';
		result[1] = '.';
		startIdx = 2;
	}

	for ( index_t i = startIdx; i < a_Length + startIdx; i++ )
	{
		if ( decimalPointIdx == i )
			result[i] = '.';
		else
			result[i] = static_cast<char>( randInt( 0, 9 ) ) + '0';	// Get random number between 0 and 10
	}
	result[a_Length] = '\0';
}