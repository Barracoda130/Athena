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

void*
atn_default_allocate( size_t size )
{
    void* ret;
    ret = malloc( size );
    if ( ret == NULL )
    {
        fprintf( stderr, "MPFR: Can't allocate memory (size=%lu)\n",
            (unsigned long)size );
        abort();
    }
    return ret;
}

void*
atn_default_reallocate( void* oldptr, size_t old_size, size_t new_size )
{
    void* ret;
    ret = realloc( oldptr, new_size );
    if ( ret == NULL )
    {
        fprintf( stderr,
            "MPFR: Can't reallocate memory (old_size=%lu new_size=%lu)\n",
            (unsigned long)old_size, (unsigned long)new_size );
        abort();
    }
    return ret;
}

void
atn_default_free( void* blk_ptr, size_t blk_size )
{
    free( blk_ptr );
}

void*
atn_tmp_allocate( struct tmp_marker** tmp_marker, size_t size )
{
    struct tmp_marker* head;

    head = (struct tmp_marker*)
        atn_default_allocate( sizeof( struct tmp_marker ) );
    head->ptr = atn_default_allocate( size );
    head->size = size;
    head->next = *tmp_marker;
    *tmp_marker = head;
    return head->ptr;
}

void
atn_tmp_free( struct tmp_marker* tmp_marker )
{
    struct tmp_marker* t;

    while ( tmp_marker != NULL )
    {
        t = tmp_marker;
        atn_default_free( t->ptr, t->size );
        tmp_marker = t->next;
        atn_default_free( t, sizeof( struct tmp_marker ) );
    }
}