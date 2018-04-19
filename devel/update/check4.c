
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of add_chrono() and get_chrono().
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
#include "lattice.h"
#include "sflds.h"
#include "linalg.h"
#include "mdflds.h"
#include "update.h"
#include "global.h"


static void set_psi(spinor_dble **chi,spinor_dble *psi)
{
   int i;
   double t;
   complex_dble z;

   t=mdtime();
   assign_sd2sd(VOLUME,chi[0],psi);

   for (i=1;i<4;i++)
   {
      z.re=pow(t,(double)(i));
      z.im=0.0;
      mulc_spinor_add_dble(VOLUME,psi,chi[i],z);
   }
}


int main(int argc,char *argv[])
{
   int my_rank,i;
   int nop,iop,itu;
   int ncr,ifr,zero;
   double phi[2],phi_prime[2],theta[3];
   double kappa,mu,eps,dev;
   spinor_dble **chi,**wsd;
   mdstep_t *s,*sm;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);

      printf("\n");
      printf("Check of add_chrono() and get_chrono()\n");
      printf("--------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   mu=0.5;
   zero=0;
   ncr=4;

   kappa=0.1365;
   set_lat_parms(5.3,1.6667,1,&kappa,1.789);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;
   set_bc_parms(0,1.0,1.0,1.0,1.0,phi,phi_prime,theta);

   set_hmc_parms(0,NULL,1,1,&mu,2,2.0);
   ifr=0;
   set_mdint_parms(0,OMF4,0.0,1,1,&ifr);
   ifr=1;
   set_mdint_parms(1,OMF4,0.2,ncr,1,&ifr);

   set_force_parms(0,FRG,0,0,0,NULL,NULL,NULL);
   set_force_parms(1,FRF_TM1,0,0,0,&zero,&zero,&ncr);

   print_mdint_parms();
   print_force_parms();

   start_ranlux(0,1234);
   geometry();
   alloc_wsd(6);
   chi=reserve_wsd(4);
   wsd=reserve_wsd(2);

   setup_chrono();
   set_mdsteps();
   s=mdsteps(&nop,&itu);
   sm=s+nop;

   for (i=0;i<4;i++)
      random_sd(VOLUME,chi[i],1.0);

   for (;s<sm;s++)
   {
      iop=(*s).iop;
      eps=(*s).eps;

      if (iop==itu)
         step_mdtime(eps);
      else if (iop==1)
      {
         set_psi(chi,wsd[0]);

         if (get_chrono(1,wsd[1]))
         {
            mulr_spinor_add_dble(VOLUME,wsd[1],wsd[0],-1.0);
            dev=norm_square_dble(VOLUME,1,wsd[1])/
               norm_square_dble(VOLUME,1,wsd[0]);

            if (my_rank==0)
               printf("t = %.3f, dev = %.1e\n",mdtime(),sqrt(dev));
         }

         add_chrono(1,wsd[0]);
      }
   }

   if (my_rank==0)
   {
      printf("\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
