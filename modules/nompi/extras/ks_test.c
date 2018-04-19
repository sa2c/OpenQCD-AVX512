
/*******************************************************************************
*
* File ks_test.c
*
* Copyright (C) 2005, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Kolmogorov-Smirnov test
*
* The externally accessible functions are
*
*   void ks_test(int n,double f[],double *pkp,double *pkm)
*     For a given array f[0],f[1],...,f[n-1], the program calculates
*     the Kolmogorov-Smirnov statistics K_n^{+}=*pkp and K_n^{-}=*pkm
*
*   void ks_prob(int n,double kp,double km,double *pp,double *pm)
*     Computes the approximate probabilites *pp and *pm for the Kolmogorov-
*     Smirnov statistics K_n^{+} and K_n^{-} to be less than or equal to 
*     kp and km respectively (eq.(4) in the notes).
*
* Notes:
*
* See the notes
*
*   M. Luescher: Statistical tests
*   
* for a detailed description.
*
*******************************************************************************/

#define KS_TEST_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "extras.h"


void ks_test(int n, double f[],double *pkp,double *pkm)
{
   int *pn,k,i;
   double *pu,*pv,xn,sn,x,kp,km;

   if (n<=0)
   {
      printf("Error in ks_test: argument out of range\n");
      printf("Program aborted\n\n");
      exit(0);
   }

   pn=malloc((n+1)*sizeof(int));
   pu=malloc((n+1)*sizeof(double));
   pv=malloc((n+1)*sizeof(double));
   xn=(double)n;
 
   if (pn&&pu&&pv)
   {
      for (k=0;k<=n;k++)
      {
         pn[k]=0;
         pu[k]=xn;
         pv[k]=0.0;
      }
   }
   else
   {
      printf("Error in ks_test: could not allocate auxiliary arrays\n");
      printf("Program aborted\n\n");
      exit(0);
   }

   for (i=0;i<n;i++)
   {
      x=xn*f[i];
      if ((x<0)||(x>xn))
      {
         printf("Error in ks_test: argument out of range\n");
         printf("Program aborted\n\n");
         exit(0);
      }

      k=(int)x;
      pn[k]+=1;
      if (x<pu[k])
         pu[k]=x;
      if (x>pv[k])
         pv[k]=x;
   }

   sn=0.0;
   kp=0.0;
   km=0.0;

   for (k=0;k<=n;k++)
   {
      if (pn[k]>0)
      {
         x=pu[k]-sn;
         if (x>km)
            km=x;
         sn+=(double)pn[k];
         x=sn-pv[k];
         if (x>kp)
            kp=x;
      }
   }

   *pkp=kp/sqrt(xn);
   *pkm=km/sqrt(xn);

   free(pn);
   free(pu);
   free(pv);
}


void ks_prob(int n,double kp,double km,double *pp,double *pm)
{
   double xn;

   if (n<=0)
   {
      printf("Error in ks_prob: argument out of range\n");
      printf("Program aborted\n\n");
      exit(0);
   }

   xn=(double)n;

   if (kp<1e-8)
      *pp=0.0;
   else if (kp>3.5)
      *pp=1.0;
   else
      *pp=1.0-exp(-2.0*kp*kp)*(1.0-2.0*kp/(3.0*sqrt(xn)));

   if (km<1e-8)
      *pm=0.0;
   else if (km>3.5)
      *pm=1.0;
   else
      *pm=1.0-exp(-2.0*km*km)*(1.0-2.0*km/(3.0*sqrt(xn)));
}

