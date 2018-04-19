
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2005, 2009, 2010, 2011 Martin Luescher, Filippo Palombi
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Checks of the programs in the module liealg
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "utils.h"
#include "lattice.h"
#include "linalg.h"
#include "global.h"

#define NMOM 100033

static double var[64],var_all[64];


int main(int argc,char *argv[])
{
   int my_rank,n,i,j;
   double dev,dmax,dmax_all;
   double nsq1,nsq2,sprod1,sprod2;
   double sm,r[8];
   double rn,cij,eij;
   su3_dble *M,*m,w;
   su3_alg_dble *X,*Y,*x;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);

      printf("\n");
      printf("Checks of the programs in the module liealg\n");
      printf("-------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      printf("Number of momenta: %d\n\n",NMOM);
   }

   start_ranlux(0,123456);
   geometry();

   X=amalloc(2*NMOM*sizeof(*X),4);
   M=amalloc(NMOM*sizeof(*M),4);
   error((X==NULL)||(M==NULL),1,
         "main [check1.c]","Unable to allocate field arrays");
   Y=X+NMOM;

   set_alg2zero(NMOM,X);
   dmax=0.0;

   for (n=0;n<NMOM;n++)
   {
      r[0]=X[n].c1;
      r[1]=X[n].c2;
      r[2]=X[n].c3;
      r[3]=X[n].c4;
      r[4]=X[n].c5;
      r[5]=X[n].c6;
      r[6]=X[n].c7;
      r[7]=X[n].c8;

      for (i=0;i<8;i++)
      {
         dev=fabs(r[i]);
         if (dev>dmax)
            dmax=dev;
      }

      X[n].c3=1.0;
   }

   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Check of set_alg2zero():\n\n");
      printf("max|X| = %.1e (should be 0.0)\n\n",dmax_all);
   }

   dmax=fabs(norm_square_alg(NMOM,1,X)-4.0*(double)(NMOM*NPROC));
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Check of norm_square_alg():\n\n");
      printf("Element count = %.1e (should be 0.0)\n\n",dmax_all);
   }

   sm=0.0;
   dmax=0.0;

   for (n=0;n<NMOM;n++)
   {
      x=X+n;
      m=M+n;

      random_alg(1,x);

      (*m).c11.re=0.0;
      (*m).c11.im=(*x).c1+(*x).c2;
      (*m).c22.re=0.0;
      (*m).c22.im=(*x).c2-2.0*(*x).c1;
      (*m).c33.re=0.0;
      (*m).c33.im=(*x).c1-2.0*(*x).c2;

      (*m).c12.re= (*x).c3;
      (*m).c12.im= (*x).c4;
      (*m).c21.re=-(*x).c3;
      (*m).c21.im= (*x).c4;

      (*m).c13.re= (*x).c5;
      (*m).c13.im= (*x).c6;
      (*m).c31.re=-(*x).c5;
      (*m).c31.im= (*x).c6;

      (*m).c23.re= (*x).c7;
      (*m).c23.im= (*x).c8;
      (*m).c32.re=-(*x).c7;
      (*m).c32.im= (*x).c8;

      su3xsu3(m,m,&w);
      nsq1=norm_square_alg(1,0,x);
      nsq2=-2.0*(w.c11.re+w.c22.re+w.c33.re);

      dev=fabs(1.0-nsq2/nsq1);
      if (dev>dmax)
         dmax=dev;

      sm+=nsq2;
   }

   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("|1.0+2*tr{X^2}/||X||^2| = %.1e (single elements)\n",dmax_all);
      printf("(should be less than %.1e or so)\n\n",DBL_EPSILON*sqrt(8.0));
   }

   dmax=fabs(1.0-sm/norm_square_alg(NMOM,0,X));
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("|1.0+2*tr{X^2}/||X||^2| = %.1e (whole vector)\n",dmax_all);
      printf("(should be less than %.1e or so)\n\n",
             DBL_EPSILON*sqrt(8.0*(double)(NMOM)));
   }

   random_alg(NMOM,X);
   random_alg(NMOM,Y);

   nsq1=norm_square_alg(NMOM,1,X);
   nsq2=norm_square_alg(NMOM,1,Y);
   sprod1=scalar_prod_alg(NMOM,1,X,Y);

   for (n=0;n<NMOM;n++)
   {
      X[n].c1+=Y[n].c1;
      X[n].c2+=Y[n].c2;
      X[n].c3+=Y[n].c3;
      X[n].c4+=Y[n].c4;
      X[n].c5+=Y[n].c5;
      X[n].c6+=Y[n].c6;
      X[n].c7+=Y[n].c7;
      X[n].c8+=Y[n].c8;
   }

   sprod2=0.5*(norm_square_alg(NMOM,1,X)-nsq1-nsq2);

   if (my_rank==0)
   {
      printf("Check of scalar_prod_alg(): %.1e\n",
             fabs(sprod1-sprod2)/sqrt(nsq1*nsq2));
      printf("(should be less than %.1e or so)\n\n",
             DBL_EPSILON*sqrt(8.0*(double)(NMOM)*(double)(NPROC)));
   }

   random_alg(NMOM,X);

   for (i=0;i<8;i++)
   {
      for (j=0;j<8;j++)
         var[8*i+j]=0.0;
   }

   for (n=0;n<NMOM;n++)
   {
      r[0]=X[n].c1;
      r[1]=X[n].c2;
      r[2]=X[n].c3;
      r[3]=X[n].c4;
      r[4]=X[n].c5;
      r[5]=X[n].c6;
      r[6]=X[n].c7;
      r[7]=X[n].c8;

      for (i=0;i<8;i++)
      {
         for (j=i;j<8;j++)
            var[8*i+j]+=r[i]*r[j];
      }
   }

   MPI_Reduce(var,var_all,64,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Check of random_alg():\n\n");
      dmax=0.0;
      rn=1.0/((double)(NMOM)*(double)(NPROC));

      for (i=0;i<8;i++)
      {
         for (j=i;j<8;j++)
         {
            if ((i==j)&&(i>1))
            {
               cij=1.0/4.0;
               eij=sqrt(2.0*rn)/4.0;
            }
            else if (i==j)
            {
               cij=1.0/9.0;
               eij=sqrt(2.0*rn)/9.0;
            }
            else if ((i==0)&&(j==1))
            {
               cij=1.0/18.0;
               eij=sqrt(5.0*rn)/18.0;
            }
            else if ((i<2)&&(j>1))
            {
               cij=0.0;
               eij=sqrt(rn)/6.0;
            }
            else
            {
               cij=0.0;
               eij=sqrt(rn)/4.0;
            }

            var_all[8*i+j]*=rn;

            if (cij!=0.0)
            {
               printf("<b[%d]*b[%d]> = % .4e, deviation = %.1e+-%.1e\n",
                      i,j,var_all[8*i+j],fabs(var_all[8*i+j]-cij),eij);
            }
            else
            {
               dev=fabs(var_all[8*i+j])/eij;

               if (dev>dmax)
                  dmax=dev;
            }
         }
      }

      eij=sqrt(rn)/4.0;
      printf("\n");
      printf("For all other i,j, ");
      printf("max|<b[i]*b[j]>| = %.1e (should be %.1e or so)\n\n",
             dmax*eij,2.0*eij);
   }

   rn=-1.2345;
   random_alg(NMOM,X);
   random_alg(NMOM,Y);

   nsq1=norm_square_alg(NMOM,1,X);
   nsq2=norm_square_alg(NMOM,1,Y);
   sprod1=scalar_prod_alg(NMOM,1,X,Y);

   muladd_assign_alg(NMOM,rn,X,Y);
   sm=norm_square_alg(NMOM,1,Y)-nsq2-rn*rn*nsq1-2.0*rn*sprod1;
   sm=fabs(sm)/nsq1;

   if (my_rank==0)
   {
      printf("Check of muladd_assign_alg(): %.1e\n",sm);
      printf("(should be less than 1.0e-15 or so)\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
