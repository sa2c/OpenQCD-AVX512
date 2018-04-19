
/*******************************************************************************
*
* File check7.c
*
* Copyright (C) 2010 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Reweighting of gaussian distributions [statistical test of gauss_dble()]
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "utils.h"
#include "extras.h"

#define NSMALL 100
#define NLARGE 1000
#define NTEST  10000

#if (NLARGE<NSMALL)
#error : NLARGE must not be smaller than NSMALL!
#endif

static double avg[4],sig[4],dev[4],*rnd,*obs[4];


static void alloc_arrays(void)
{
   rnd=amalloc((NLARGE+4*NTEST)*sizeof(*rnd),4);

   error(rnd==NULL,1,"alloc_arrays [check7.c]",
         "Unable to allocate auxiliary arrays");

   obs[0]=rnd+NLARGE;
   obs[1]=obs[0]+NTEST;
   obs[2]=obs[1]+NTEST;
   obs[3]=obs[2]+NTEST;
}


static double set_rnd(int nms)
{
   int i;
   double sm;

   gauss_dble(rnd,nms);
   sm=0.0;

   for (i=0;i<nms;i++)
   {
      rnd[i]=2.0*rnd[i]*rnd[i];
      sm+=rnd[i];
   }

   return sm/(double)(nms);
}


static double get_obs(int nms,double sigsq)
{
   int i;
   double fact,sm;

   fact=-1.0/(2.0*sigsq);
   sm=0.0;

   for (i=0;i<nms;i++)
      sm+=exp(fact*rnd[i]);

   sm/=(double)(nms);

   return sigsq*(1.0/(sm*sm)-1.0);   
}


static void set_avgsig(void)
{
   int i,j;
   double d;

   for (i=0;i<4;i++)
   {
      avg[i]=0.0;
      sig[i]=0.0;
      dev[i]=0.0;
   }

   for (i=0;i<4;i++)
   {
      for (j=0;j<NTEST;j++)
         avg[i]+=obs[i][j];
   }

   for (i=0;i<4;i++)
      avg[i]/=(double)(NTEST);

   for (i=0;i<4;i++)
   {
      for (j=0;j<NTEST;j++)
      {
         d=fabs(obs[i][j]-avg[i]);
         sig[i]+=d*d;

         if (d>dev[i])
            dev[i]=d;
      }
   }

   for (i=0;i<4;i++)
   {
      sig[i]/=(double)(NTEST);
      sig[i]=sqrt(sig[i]);
   }
}
   

int main(void)
{
   int i;

   printf("\n");
   printf("Reweighting of gaussian distributions\n");
   printf("-------------------------------------\n\n");

   printf("Width of the distribution = 1.0\n");
   printf("Width of the observable = 2.0,3.0,4.0\n");
   printf("%d test simulations of size %d and %d\n\n",NTEST,NSMALL,NLARGE);

   alloc_arrays();
   
   printf("Sample size %d:\n\n",NSMALL);

   for (i=0;i<NTEST;i++)
   {
      obs[0][i]=set_rnd(NSMALL);
      obs[1][i]=get_obs(NSMALL,2.0);
      obs[2][i]=get_obs(NSMALL,3.0);
      obs[3][i]=get_obs(NSMALL,4.0);      
   }

   set_avgsig();

   printf("bias = % .4e, error = %.1e, maxerr = %.1e\n",
          avg[0]-1.0,sig[0],dev[0]);
   printf("bias = % .4e, error = %.1e, maxerr = %.1e (sigma=2.0)\n",
          avg[1]-1.0,sig[1],dev[1]);
   printf("bias = % .4e, error = %.1e, maxerr = %.1e (sigma=3.0)\n",
          avg[2]-1.0,sig[2],dev[2]);
   printf("bias = % .4e, error = %.1e, maxerr = %.1e (sigma=4.0)\n",
          avg[3]-1.0,sig[3],dev[3]);
   printf("\n");

   printf("Sample size %d:\n\n",NLARGE);

   for (i=0;i<NTEST;i++)
   {
      obs[0][i]=set_rnd(NLARGE);
      obs[1][i]=get_obs(NLARGE,2.0);
      obs[2][i]=get_obs(NLARGE,3.0);
      obs[3][i]=get_obs(NLARGE,4.0);      
   }

   set_avgsig();

   printf("bias = % .4e, error = %.1e, maxerr = %.1e\n",
          avg[0]-1.0,sig[0],dev[0]);
   printf("bias = % .4e, error = %.1e, maxerr = %.1e (sigma=2.0)\n",
          avg[1]-1.0,sig[1],dev[1]);
   printf("bias = % .4e, error = %.1e, maxerr = %.1e (sigma=3.0)\n",
          avg[2]-1.0,sig[2],dev[2]);
   printf("bias = % .4e, error = %.1e, maxerr = %.1e (sigma=4.0)\n",
          avg[3]-1.0,sig[3],dev[3]);
   printf("\n");

   exit(0);
}
