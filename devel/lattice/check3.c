
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the programs bnd_s2zero() and bnd_sd2zero().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "random.h"
#include "su3fcts.h"
#include "utils.h"
#include "sflds.h"
#include "lattice.h"
#include "global.h"

#define N0 (NPROC0*L0)
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


static int is_zero(spinor *s)
{
   int i,ie;
   spin_t *sp;

   sp=(spin_t*)(s);
   ie=1;

   for (i=0;i<24;i++)
      ie&=((*sp).r[i]==0.0f);

   return ie;
}


static int is_zero_dble(spinor_dble *s)
{
   int i,ie;
   spin_dble_t *sp;

   sp=(spin_dble_t*)(s);
   ie=1;

   for (i=0;i<24;i++)
      ie&=((*sp).r[i]==0.0);

   return ie;
}


static int check_sbnd(ptset_t set,spinor *s)
{
   int bc,ix,t;
   int io,ie;

   bc=bc_type();
   ie=1;

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);
      io=((set==ALL_PTS)||
          ((set==EVEN_PTS)&&(ix<(VOLUME/2)))||
          ((set==ODD_PTS)&&(ix>=(VOLUME/2))));

      if ((io!=0)&&(((t==0)&&(bc!=3))||((t==(N0-1))&&(bc==0))))
         ie&=is_zero(s);
      else
         ie&=(is_zero(s)^0x1);

      s+=1;
   }

   return ie;
}


static int check_sbnd_dble(ptset_t set,spinor_dble *s)
{
   int bc,ix,t;
   int io,ie;

   bc=bc_type();
   ie=1;

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);
      io=((set==ALL_PTS)||
          ((set==EVEN_PTS)&&(ix<(VOLUME/2)))||
          ((set==ODD_PTS)&&(ix>=(VOLUME/2))));

      if ((io!=0)&&(((t==0)&&(bc!=3))||((t==(N0-1))&&(bc==0))))
         ie&=is_zero_dble(s);
      else
         ie&=(is_zero_dble(s)^0x1);

      s+=1;
   }

   return ie;
}


int main(int argc,char *argv[])
{
   int my_rank,bc,ie,is,k;
   double phi[2],phi_prime[2],theta[3];
   spinor **ps;
   spinor_dble **psd;
   ptset_t set;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);
      printf("\n");
      printf("Check of the programs bnd_s2zero() and bnd_sd2zero()\n");
      printf("----------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>]");
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime,theta);
   print_bc_parms(0);

   start_ranlux(0,12345);
   geometry();
   alloc_ws(NFLDS);
   alloc_wsd(NFLDS);

   ps=reserve_ws(NFLDS);
   psd=reserve_wsd(NFLDS);
   ie=1;

   for (is=0;is<4;is++)
   {
      if (is==0)
         set=EVEN_PTS;
      else if (is==1)
         set=ODD_PTS;
      else if (is==2)
         set=ALL_PTS;
      else
         set=NO_PTS;

      for (k=0;k<NFLDS;k++)
      {
         random_s(NSPIN,ps[k],1.0f);
         random_sd(NSPIN,psd[k],1.0);

         bnd_s2zero(set,ps[k]);
         bnd_sd2zero(set,psd[k]);

         ie&=check_sbnd(set,ps[k]);
         ie&=check_sbnd_dble(set,psd[k]);
      }
   }

   error(ie!=1,1,"main [check3.c]",
         "Spinor fields vanish on an incorrect set of points");

   if (my_rank==0)
   {
      printf("No errors detected --- all programs work correctly\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
