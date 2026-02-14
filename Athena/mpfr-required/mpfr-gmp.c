#include "mpfr-gmp.h"

#include "mpfr-impl.h"

MPFR_COLD_FUNCTION_ATTR void
mpfr_assert_fail( const char* filename, int linenum,
                  const char* expr )
{
    if ( filename != NULL && filename[0] != '\0' )
    {
        fprintf( stderr, "%s:", filename );
        if ( linenum != -1 )
            fprintf( stderr, "%d: ", linenum );
    }
    fprintf( stderr, "MPFR assertion failed: %s\n", expr );
    abort();
}