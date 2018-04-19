
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2005, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the programs in the module sflds.c.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "utils.h"
#include "lattice.h"
#include "linalg.h"
#include "sflds.h"
#include "global.h"

#define NFLDS 3

typedef union
{
   spinor s;
   float r[24];
} spin_t;

typedef union
{
   spinor_dble s;
   double r[24];
} spin_dble_t;

static float sig[NFLDS];
static double sigd[NFLDS];


int main(int argc,char *argv[])
{
   int my_rank,ie,k,i,ix;
   float *r;
   double *rd,var,var_all,d,dmax;
   spinor **ps;
   spinor_dble **psd;
   spin_t *sps;
   spin_dble_t *spsd;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);
      printf("\n");
      printf("Check of the programs in the module sflds.c\n");
      printf("-------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   start_ranlux(0,12345);
   geometry();
   alloc_ws(2*NFLDS);
   alloc_wsd(2*NFLDS);
   ps=reserve_ws(2*NFLDS);
   psd=reserve_wsd(2*NFLDS);
   ie=0;

   for (k=0;k<NFLDS;k++)
   {
      sps=(spin_t*)(ps[k]);

      for (ix=0;ix<NSPIN;ix++)
      {
         r=(*sps).r;

         for (i=0;i<24;i++)
            if (r[i]!=0.0f)
               ie=1;

         sps+=1;
      }
   }

   error(ie,1,"main [check1.c]",
         "Single-precision fields are not properly initialized");

   if (my_rank==0)
   {
      printf("Choose random single-precision spinor fields\n");
      ranlxs(sig,NFLDS);
   }

   MPI_Bcast(sig,NFLDS,MPI_FLOAT,0,MPI_COMM_WORLD);

   for (k=0;k<NFLDS;k++)
   {
      random_s(VOLUME,ps[k],sig[k]);
      var=0.0;
      sps=(spin_t*)(ps[k]);

      for (ix=0;ix<VOLUME;ix++)
      {
         r=(*sps).r;

         for (i=0;i<24;i++)
            var+=(double)(r[i]*r[i]);

         sps+=1;
      }

      MPI_Reduce(&var,&var_all,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      if (my_rank==0)
      {
         var_all/=((double)(12*NPROC)*(double)(VOLUME));
         printf("<s[%d]^2> = %.4e (sigma^2 = %.4e)\n",
                k,var_all,sig[k]*sig[k]);
      }
   }

   ie=0;

   for (k=0;k<NFLDS;k++)
   {
      set_s2zero(VOLUME,ps[k]);
      d=(double)(norm_square(VOLUME,1,ps[k]));
      ie|=(d!=0.0);
   }

   error(ie!=0,1,"main [check1.c]",
         "Fields are not properly set to zero by set_s2zero()");

   for (k=0;k<NFLDS;k++)
   {
      spsd=(spin_dble_t*)(psd[k]);

      for (ix=0;ix<NSPIN;ix++)
      {
         rd=(*spsd).r;

         for (i=0;i<24;i++)
            if (rd[i]!=0.0)
               ie=1;

         spsd+=1;
      }
   }

   error(ie,1,"main [check1.c]",
         "Double-precision fields are not properly initialized");

   if (my_rank==0)
   {
      printf("\n");
      printf("Choose random double-precision spinor fields\n");
      ranlxd(sigd,NFLDS);
   }

   MPI_Bcast(sigd,NFLDS,MPI_DOUBLE,0,MPI_COMM_WORLD);

   for (k=0;k<NFLDS;k++)
   {
      random_sd(VOLUME,psd[k],sigd[k]);
      var=0.0;
      spsd=(spin_dble_t*)(psd[k]);

      for (ix=0;ix<VOLUME;ix++)
      {
         rd=(*spsd).r;

         for (i=0;i<24;i++)
            var+=(rd[i]*rd[i]);

         spsd+=1;
      }

      MPI_Reduce(&var,&var_all,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      if (my_rank==0)
      {
         var_all/=((double)(12*NPROC)*(double)(VOLUME));
         printf("<sd[%d]^2> = %.4e (sigma^2 = %.4e)\n",
                k,var_all,sigd[k]*sigd[k]);
      }
   }

   ie=0;

   for (k=0;k<NFLDS;k++)
   {
      set_sd2zero(VOLUME,psd[k]);
      d=norm_square_dble(VOLUME,1,psd[k]);
      ie|=(d!=0.0);
   }

   error(ie!=0,1,"main [check1.c]",
         "Fields are not properly set to zero by set_sd2zero()");

   for (k=0;k<NFLDS;k++)
   {
      random_s(VOLUME,ps[k],1.0f);
      random_sd(VOLUME,psd[k],1.0);
      assign_s2s(VOLUME,ps[k],ps[k+NFLDS]);
      assign_sd2sd(VOLUME,psd[k],psd[k+NFLDS]);

      mulr_spinor_add(VOLUME,ps[k],ps[k+NFLDS],-1.0f);
      mulr_spinor_add_dble(VOLUME,psd[k],psd[k+NFLDS],-1.0);
      sps=(spin_t*)(ps[k]);
      spsd=(spin_dble_t*)(psd[k]);

      for (ix=0;ix<VOLUME;ix++)
      {
         r=(*sps).r;
         rd=(*spsd).r;

         for (i=0;i<24;i++)
         {
            if (r[i]!=0.0f)
               ie=1;
            if (rd[i]!=0.0)
               ie=2;
         }
      }
   }

   error(ie==1,1,"main [check1.c]","assign_s2s() is incorrect");
   error(ie==2,1,"main [check1.c]","assign_sd2sd() is incorrect");

   for (k=0;k<NFLDS;k++)
   {
      random_s(VOLUME,ps[k],1.0f);
      assign_s2sd(VOLUME,ps[k],psd[k]);
      assign_sd2s(VOLUME,psd[k],ps[k+NFLDS]);

      mulr_spinor_add(VOLUME,ps[k],ps[k+NFLDS],-1.0f);
      sps=(spin_t*)(ps[k]);

      for (ix=0;ix<VOLUME;ix++)
      {
         r=(*sps).r;

         for (i=0;i<24;i++)
            if (r[i]!=0.0f)
               ie=1;
      }
   }

   error(ie==1,1,"main [check1.c]",
         "assign_s2sd() or assign_sd2s() is incorrect");

   dmax=0.0;

   for (k=0;k<NFLDS;k++)
   {
      random_s(VOLUME,ps[k],1.0f);
      random_s(VOLUME,ps[k+NFLDS],1.0f);
      assign_s2sd(VOLUME,ps[k],psd[k]);
      assign_s2sd(VOLUME,ps[k+NFLDS],psd[k+NFLDS]);

      diff_s2s(VOLUME,ps[k],ps[k+NFLDS]);
      mulr_spinor_add_dble(VOLUME,psd[k+NFLDS],psd[k],-1.0);
      d=norm_square_dble(VOLUME,1,psd[k+NFLDS]);
      assign_s2sd(VOLUME,ps[k+NFLDS],psd[k]);
      mulr_spinor_add_dble(VOLUME,psd[k+NFLDS],psd[k],1.0);

      d=norm_square_dble(VOLUME,1,psd[k+NFLDS])/d;
      if (d>dmax)
         dmax=d;
   }

   if (my_rank==0)
   {
      printf("\n");
      printf("Relative deviations (should be less than 1.0e-7 or so):\n");
      printf("diff_s2s():  %.1e\n",sqrt(dmax));
   }

   dmax=0.0;

   for (k=0;k<NFLDS;k++)
   {
      random_s(VOLUME,ps[k],1.0f);
      random_sd(VOLUME,psd[k],1.0);
      assign_sd2sd(VOLUME,psd[k],psd[k+NFLDS]);

      add_s2sd(VOLUME,ps[k],psd[k]);
      d=norm_square_dble(VOLUME,1,psd[k]);
      mulr_spinor_add_dble(VOLUME,psd[k],psd[k+NFLDS],-1.0);
      assign_s2sd(VOLUME,ps[k],psd[k+NFLDS]);
      mulr_spinor_add_dble(VOLUME,psd[k],psd[k+NFLDS],-1.0);

      d=norm_square_dble(VOLUME,1,psd[k])/d;
      if (d>dmax)
         dmax=d;
   }

   if (my_rank==0)
      printf("add_s2sd():  %.1e\n",sqrt(dmax));

   dmax=0.0;

   for (k=0;k<NFLDS;k++)
   {
      random_sd(VOLUME,psd[k],1.0);
      random_sd(VOLUME,psd[k+NFLDS],1.0);

      diff_sd2s(VOLUME,psd[k],psd[k+NFLDS],ps[k]);
      mulr_spinor_add_dble(VOLUME,psd[k],psd[k+NFLDS],-1.0);
      d=norm_square_dble(VOLUME,1,psd[k]);
      assign_s2sd(VOLUME,ps[k],psd[k+NFLDS]);
      mulr_spinor_add_dble(VOLUME,psd[k],psd[k+NFLDS],-1.0);

      d=norm_square_dble(VOLUME,1,psd[k])/d;
      if (d>dmax)
         dmax=d;
   }

   if (my_rank==0)
   {
      printf("diff_sd2s(): %.1e\n\n",sqrt(dmax));
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
