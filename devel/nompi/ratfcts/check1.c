
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2008, 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the complete elliptic integral K(k)
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "random.h"
#include "utils.h"
#include "ratfcts.h"


static double Ksmall(double rk)
{
   double c0,c1,c2,c3;
   double k,p;

   c0=1.0;
   c1=1.0/4.0;
   c2=9.0/64.0;
   c3=25.0/256.0;
   
   k=(rk*rk)/(1.0+rk*rk);

   p=c2+k*c3;
   p=c1+k*p;
   p=c0+k*p;

   return 2.0*atan(1.0)*p;
}


int main(void)
{
   int n;
   double rk,k,kp,km,dev,dmax;
   
   printf("\n");
   printf("Computation of the complete elliptic integral K(k)\n");
   printf("--------------------------------------------------\n\n");

   rlxd_init(1,1234);   

   km=pow(DBL_EPSILON,0.125);
   dmax=fabs(1.0-Ksmall(0.0)/ellipticK(0.0));

   for (n=0;n<1000;n++)
   {
      ranlxd(&rk,1);
      rk*=km;

      dev=fabs(1.0-Ksmall(rk)/ellipticK(rk));

      if (dev>dmax)
         dmax=dev;
   }

   printf("Small k region: maximal relative error = %.1e\n",dmax);

   dmax=0.0;
   
   for (n=0;n<1000;n++)
   {
      ranlxd(&rk,1);
      rk=rk/(1.0-rk);

      k=rk/sqrt(1.0+rk*rk);
      kp=1.0/sqrt(1.0+rk*rk);
      
      dev=fabs(1.0-
               ellipticK(2.0*sqrt(k)*(1.0+k)/(kp*kp))/
               ((1.0+k)*ellipticK(rk)));

      if (dev>dmax)
         dmax=dev;
   }      
               
   printf("Gauss transformation: maximal relative error = %.1e\n\n",dmax);
   printf("Print values at specfied k/k'\n\n");
   
   for (;;)
   {
      printf("k/k' =  ");

      if (scanf("%lf",&rk)==1)
      {
         printf("k = %.8e,  k' = %.8e, K(k) = %.16e\n\n",
                rk/sqrt(1.0+rk*rk),1.0/sqrt(1.0+rk*rk),ellipticK(rk));
      }
      else
      {
         printf("Invalid input value, program stopped\n\n");
         break;
      }
   }
   
   exit(0);
}
