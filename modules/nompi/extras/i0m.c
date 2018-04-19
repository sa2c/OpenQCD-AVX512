
/*******************************************************************************
*
* File i0m.c
*
* Copyright (C) 2010 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the modified Bessel function I_0(x)
*
* The externally accessible function is
*
*   double i0m(double x)
*     This program returns exp(-x)*I_0(x) to machine precision
*     for x>=0. An error occurs if x is negative
*
* Notes:
*
* The Bessel function is calculated by evaluating the integral
* 
*   exp(-x)*I_0(x)=int_0^Pi (dt/Pi)*exp(-x*(1-cos(t))) 
*
* using Chebyshev polynomials
*
*******************************************************************************/

#define I0M_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "utils.h"
#include "extras.h"

static double pi,xs;


static double maxt(double x)
{
   double r;

   pi=4.0*atan(1.0);
   
   if (x<1.0)
      return pi;

   r=1.0-(0.5*log(2.0*pi*x)-log(DBL_EPSILON))/x;

   if (r>=1.0)
      return 0.0;
   else if (r<=-1.0)
      return pi;
   else 
      return acos(r);
}


static double f(double t)
{
   return exp(-xs*(1.0-cos(t)));
}


double i0m(double x)
{
   double a,b;
   
   if (x==0.0)
      return 1.0;

   error(x<0.0,1,"i0m [i0.c]","The argument x must be non-negative");

   a=0.0;
   b=maxt(x);
   xs=x;

   if (b==0.0)
      return (1.0/sqrt(2.0*pi*x))*(1.0+1.0/(8.0*x));

   return cheby_int(a,b,f,512,10.0*DBL_EPSILON)/pi;
}
