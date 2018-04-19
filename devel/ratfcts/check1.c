
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2009, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Initialization of rational functions
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "ratfcts.h"
#include "global.h"


static double eval_rat1(ratfct_t *rf,double x)
{
   int np,i;
   double *mu,*rmu,r;

   np=(*rf).np;
   mu=(*rf).mu;
   rmu=(*rf).rmu;
   r=0.0;

   for (i=0;i<np;i++)
      r+=rmu[i]/(x*x+mu[i]*mu[i]);

   return 1.0+r;
}


static double eval_rat2(ratfct_t *rf,double x)
{
   int np,i;
   double *nu,*rnu,r;
   complex_dble z;

   np=(*rf).np;
   nu=(*rf).nu;
   rnu=(*rf).rnu;
   z.re=1.0;
   z.im=0.0;

   for (i=0;i<np;i++)
   {
      r=rnu[i]/(x*x+nu[i]*nu[i]);
      z.re+=nu[i]*r;
      z.im+=x*r;
   }

   return z.re*z.re+z.im*z.im;
}


static double diff_rat1(double ra,double rb,ratfct_t *rf)
{
   int k;
   double r,x,d,dmax;

   dmax=0.0;
   
   for (k=0;k<1000;k++)
   {
      ranlxd(&r,1);
      x=ra+r*(rb-ra);

      d=fabs(1.0-(*rf).A*x*eval_rat1(rf,x));
      
      if (d>dmax)
         dmax=d;
   }

   d=dmax;
   MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   return dmax;
}


static double diff_rat2(double ra,double rb,ratfct_t *rf)
{
   int k;
   double r,x,d,dmax;

   dmax=0.0;
   
   for (k=0;k<1000;k++)
   {
      ranlxd(&r,1);
      x=ra+r*(rb-ra);

      d=fabs(1.0-eval_rat1(rf,x)*eval_rat2(rf,x));

      if (d>dmax)
         dmax=d;
   }

   d=dmax;
   MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   return dmax;
}


static double diff_rat3(double ra,double rb,ratfct_t *rf)
{
   int k;
   double r,x,d,dmax;

   dmax=0.0;
   
   for (k=0;k<1000;k++)
   {
      ranlxd(&r,1);
      x=ra+r*(rb-ra);

      d=fabs(1.0-(eval_rat1(rf+1,x)*eval_rat1(rf+2,x))/eval_rat1(rf,x));

      if (d>dmax)
         dmax=d;
   }

   d=dmax;
   MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   return dmax;
}


int main(int argc,char *argv[])
{
   int my_rank,irat[3];
   int np1,i,j;
   double dmax;
   rat_parms_t rp;
   ratfct_t rf[3];
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);
      fin=freopen("check1.in","r",stdin);
      
      printf("\n");
      printf("Initialization of rational functions\n");
      printf("------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   read_rat_parms(0);

   if (my_rank==0)
      fclose(fin);

   print_rat_parms();
   start_ranlux(0,123456);

   rp=rat_parms(0);
   irat[0]=0;
   irat[1]=0;
   irat[2]=rp.degree-1;
   rf[0]=ratfct(irat);

   if (my_rank==0)
   {
      printf("Complete rational function:\n");
      printf("np= %2d, A = %.2e, delta = %.2e\n",
             rf[0].np,rf[0].A,rf[0].delta);
   
      printf("  i      mu[i]       rmu[i]\n");

      for (i=0;i<rf[0].np;i++)
         printf("  %2d   %.4e  % .4e\n",
                i,rf[0].mu[i],rf[0].rmu[i]);

      printf("  i      nu[i]       rnu[i]\n");

      for (i=0;i<rf[0].np;i++)
         printf("  %2d   %.4e  % .4e\n",
                i,rf[0].nu[i],rf[0].rnu[i]);
   }
   
   dmax=diff_rat1(rp.range[0],rp.range[1],rf);

   if (my_rank==0)
      printf("Check approximation = %.1e (should be at most %.1e)\n",
             dmax,rf[0].delta);

   dmax=diff_rat2(rp.range[0],rp.range[1],rf);

   if (my_rank==0)
      printf("Check decomposition = %.1e\n\n",dmax);

   np1=0;

   for (i=0;i<rf[0].np;i++)
      np1+=(rf[0].mu[i]>0.1);

   if ((np1>0)&&(np1<rf[0].np))
   {
      irat[0]=0;
      irat[1]=0;
      irat[2]=np1-1;
      rf[1]=ratfct(irat);
   
      irat[0]=0;
      irat[1]=np1;
      irat[2]=rf[0].np-1;
      rf[2]=ratfct(irat);

      if (my_rank==0)
         printf("Factorization R0=R1*R2:\n\n");

      for (i=1;i<3;i++)
      {
         if (my_rank==0)
         {
            printf("Rational function no %d:\n",i);
            printf("np= %2d, A = %.2e, delta = %.2e\n",
                   rf[i].np,rf[i].A,rf[i].delta);

            printf("  i      mu[i]       rmu[i]\n");

            for (j=0;j<rf[i].np;j++)
               printf("  %2d   %.4e  % .4e\n",
                      j,rf[i].mu[j],rf[i].rmu[j]);

            printf("  i      nu[i]       rnu[i]\n");

            for (j=0;j<rf[i].np;j++)
               printf("  %2d   %.4e  % .4e\n",
                      j,rf[i].nu[j],rf[i].rnu[j]);

            printf("\n");
         }
      }

      dmax=diff_rat3(rp.range[0],rp.range[1],rf);

      if (my_rank==0)
         printf("Check factorization = %.1e\n\n",dmax);      
   }

   if (my_rank==0)
      fclose(flog);      
   
   MPI_Finalize();
   exit(0);
}
