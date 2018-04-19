
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2008 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Zolotarev rational approximation to the function f(x)=1/|x|
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "random.h"
#include "utils.h"
#include "ratfcts.h"

static int ns=0;
static double As,*ars;


static void alloc_ars(int n)
{
   if (n<=ns)
      return;

   if (ns!=0)
      afree(ars);

   ars=amalloc(2*n*sizeof(double),3);

   error(ars==NULL,1,"alloc_ars [check3.c]",
         "Unable to allocate coefficient array");

   ns=n;
}


static double Zolo(int n,double y)
{
   int r;
   double p;

   p=1.0;

   for (r=0;r<n;r++)
      p*=((y+ars[2*r])/(y+ars[2*r+1]));

   return As*p;
}


int main(void)
{
   int n,r,ir;
   double a,b;
   double eps,delta;
   double y,dev,dmax;
   
   printf("\n");
   printf("Zolotarev rational approximation to the function f(x)=1/|x|\n");
   printf("-----------------------------------------------------------\n\n");

   rlxd_init(1,1234);   

   for (;;)
   {
      printf("Range [a,b] of |x| = ");

      ir=scanf("[");
      ir+=scanf("%lf",&a);
      ir+=scanf(",");
      ir+=scanf("%lf",&b);
      ir+=scanf("]");

      printf("Degree n = ");
      ir+=scanf("%d",&n);

      error(ir!=3,1,"main [check3.c]","Read error or invalid input data");
      error((a<=0.0)||(b<=a)||(n<1),1,"main [check3.c]",
            "Improper range or degree");

      eps=a/b;
      eps=eps*eps;

      alloc_ars(n);
      zolotarev(n,eps,&As,ars,&delta);

      dmax=0.0;

      for (r=0;r<100;r++)
      {
         ranlxd(&y,1);
         y=eps+(1.0-eps)*y;

         dev=fabs(1.0-sqrt(y)*Zolo(n,y));

         if (dev>dmax)
            dmax=dev;
      }

      printf("Relative error delta = %.1e (measured: %.1e)\n",delta,dmax);
      printf("Amplitude A: %.1e\n",As);
      printf("Coefficients a_r:     Numerator     Denominator\n");
      b=b*b;

      for (r=0;r<(2*n);r+=2)
      {
         printf("                      %.3e      %.3e\n",
                ars[r]*b,ars[r+1]*b);
      }
      
      printf("\n\n");
   }
   
   exit(0);
}
