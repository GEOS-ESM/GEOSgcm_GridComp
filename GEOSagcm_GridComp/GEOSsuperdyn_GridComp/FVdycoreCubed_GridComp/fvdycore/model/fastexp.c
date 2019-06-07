#include <math.h>
/*
 * Fast pow() reference implementation
 */

double fast_exp(double x, double y)
{
 float xx;
 float yy;
 float fret;
 double dret;
 xx=x;
 yy=y;
 fret = powf(xx,yy);
 dret = fret;
 return dret;
}

float fast_expf(float x, float y)
{
 return powf(x,y);
}

// Fortran bindings
double FASTEXP (double *x,double *y){return fast_exp(*x,*y);}
double FASTEXP_(double *x,double *y){return fast_exp(*x,*y);}
double fastexp (double *x,double *y){return fast_exp(*x,*y);}
double fastexp_(double *x,double *y){return fast_exp(*x,*y);}

float FASTEXPF (float *x,float *y){return fast_expf(*x,*y);}
float FASTEXPF_(float *x,float *y){return fast_expf(*x,*y);}
float fastexpf (float *x,float *y){return fast_expf(*x,*y);}
float fastexpf_(float *x,float *y){return fast_expf(*x,*y);}

