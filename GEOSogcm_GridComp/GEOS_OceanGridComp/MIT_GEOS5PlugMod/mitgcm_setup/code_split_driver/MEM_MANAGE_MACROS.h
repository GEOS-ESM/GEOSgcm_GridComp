
#define _DEALLOC(_a) IF ( ASSOCIATED( _a ) ) THEN; DEALLOCATE( _a ); NULLIFY( _a ); ENDIF
