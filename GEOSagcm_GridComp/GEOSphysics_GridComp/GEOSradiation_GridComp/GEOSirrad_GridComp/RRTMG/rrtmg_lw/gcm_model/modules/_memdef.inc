#ifdef _CUDA
#define dbreg(x) call dbal(x)
#define dbcop(x,y) call dbcp(x, cpointer);call c_f_pointer( cpointer, y, shape(x)) 
#define dbcopnp(x,y,t,u) call dbcp(x, cpointer);call c_f_pointer( cpointer, y, shape(x)) 
#define dreg(x,y,z) call ddbxeg(x,y,z,cpointer);call c_f_pointer( cpointer, x, [y,z] )
#define sreg(x,y,z) call ddbxeg(x,y,z,cpointer)
#define dbflushreg() call dbflushrg()
#define dbflushcop() call dbflushcp()
#else
#define dbreg(x) 
#define dbcop(x,y) y=>x
#define dbcopnp(x,y,u,v) if (allocated(y) .eqv. .true.) deallocate(y) ;allocate( y( u, v)); y=x
#define dbflushreg() 
#define dbflushcop() 
#define dreg(x,y,z) if (allocated(x) .eqv. .true.) deallocate(x) ;allocate( x( y , z))
#define sreg(x,y,z)
#endif


