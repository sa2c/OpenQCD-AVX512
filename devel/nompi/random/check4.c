
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2005, 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Kolmogorov-Smirnov test of the random distribution produced by
* gauss and gauss_dble
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "utils.h"
#include "su3fcts.h"
#include "extras.h"


int main(void)
{
   int i,n;
   float *r;
   double *rd,*f,x;
   double kp,km,pp,pm;

   printf("\n");
   printf("Check of the distribution produced by gauss and gauss_dble\n");
   printf("----------------------------------------------------------\n");

   for (;;)
   {
      printf("\n");      
      printf("Specify number of trials (0 exits): ");

      if (scanf("%d",&n)==1)
      {
         printf("\n");

         if (n<=0)
            exit(0);

         f=amalloc(n*sizeof(double),3);
         r=amalloc(n*sizeof(float),3);

         gauss(r,n);
      
         for (i=0;i<n;i++)
         {
            x=(double)(r[i]);

            if (x>=0)
               f[i]=0.5+0.5*pchi_square(2.0*x*x,1);
            else
               f[i]=0.5-0.5*pchi_square(2.0*x*x,1);
         }
      
         ks_test(n,f,&kp,&km);
         ks_prob(n,kp,km,&pp,&pm);

         printf("Distribution produced by gauss\n");
         printf("Kolmogorov-Smirnov test: K+ = %4.2f, K- = %4.2f\n",kp,km);
         printf("This corresponds to Prob(K+) = %4.2f, Prob(K-) = %4.2f\n",
                pp,pm);
         printf("\n");

         afree(r);
         rd=amalloc(n*sizeof(double),3);      
         gauss_dble(rd,n);
      
         for (i=0;i<n;i++)
         {
            x=rd[i];

            if (x>=0)
               f[i]=0.5+0.5*pchi_square(2.0*x*x,1);
            else
               f[i]=0.5-0.5*pchi_square(2.0*x*x,1);
         }
      
         ks_test(n,f,&kp,&km);
         ks_prob(n,kp,km,&pp,&pm);

         printf("Distribution produced by gauss_dble\n");
         printf("Kolmogorov-Smirnov test: K+ = %4.2f, K- = %4.2f\n",kp,km);
         printf("This corresponds to Prob(K+) = %4.2f, Prob(K-) = %4.2f\n",
                pp,pm);
         printf("\n");      
      
         afree(f);
         afree(rd);
      }
      else
      {
         printf("Invalid input, program stopped\n\n");
         break;
      }
   }
   
   exit(0);
}
