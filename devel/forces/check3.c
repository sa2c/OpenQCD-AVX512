
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2005, 2008-2013, 2016 Martin Luescher, Filippo Palombi
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the programs force0() and action0().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "mdflds.h"
#include "linalg.h"
#include "forces.h"
#include "global.h"

#define N0 (NPROC0*L0)


static void rot_ud(double eps)
{
   int bc,ix,t,ifc;
   su3_dble *u;
   su3_alg_dble *mom;
   mdflds_t *mdfs;

   bc=bc_type();
   mdfs=mdflds();
   mom=(*mdfs).mom;
   u=udfld();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         expXsu3(eps,mom,u);
         mom+=1;
         u+=1;

         if (bc!=0)
            expXsu3(eps,mom,u);
         mom+=1;
         u+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
               expXsu3(eps,mom,u);
            mom+=1;
            u+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
            expXsu3(eps,mom,u);
         mom+=1;
         u+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            expXsu3(eps,mom,u);
            mom+=1;
            u+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            expXsu3(eps,mom,u);
            mom+=1;
            u+=1;
         }
      }
   }

   set_flags(UPDATED_UD);
}


static int is_frc_zero(su3_alg_dble *f)
{
   int ie;

   ie=1;
   ie&=((*f).c1==0.0);
   ie&=((*f).c2==0.0);
   ie&=((*f).c3==0.0);
   ie&=((*f).c4==0.0);
   ie&=((*f).c5==0.0);
   ie&=((*f).c6==0.0);
   ie&=((*f).c7==0.0);
   ie&=((*f).c8==0.0);

   return ie;
}


static int check_bnd_frc(su3_alg_dble *frc)
{
   int bc,ix,t,ifc,ie;

   bc=bc_type();
   ie=0;

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((t==0)&&(bc==0))
      {
         ie|=is_frc_zero(frc);
         frc+=1;

         ie|=(is_frc_zero(frc)^0x1);
         frc+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            ie|=is_frc_zero(frc);
            frc+=1;
         }
      }
      else if ((t==0)&&(bc==1))
      {
         ie|=is_frc_zero(frc);
         frc+=1;

         ie|=is_frc_zero(frc);
         frc+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            ie|=(is_frc_zero(frc)^0x1);
            frc+=1;
         }
      }
      else if ((t==(N0-1))&&(bc==0))
      {
         ie|=(is_frc_zero(frc)^0x1);
         frc+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            ie|=is_frc_zero(frc);
            frc+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            ie|=is_frc_zero(frc);
            frc+=1;
         }
      }
   }

   return ie;
}


static double dSdt(double c)
{
   int ie;
   mdflds_t *mdfs;

   mdfs=mdflds();
   ie=check_bnd_frc((*mdfs).mom);
   error(ie!=0,1,"dSdt [check3.c]",
         "Momentum field vanishes on an incorrect set of links");

   force0(c);
   ie=check_bnd_frc((*mdfs).frc);
   error(ie!=0,1,"dSdt [check3.c]",
         "Force field vanishes on an incorrect set of links");

   return scalar_prod_alg(4*VOLUME,0,(*mdfs).mom,(*mdfs).frc);
}


static double chk_chs(double c)
{
   double dev;
   su3_alg_dble **wfd;
   mdflds_t *mdfs;

   wfd=reserve_wfd(1);
   mdfs=mdflds();

   random_ud();
   force0(c);
   assign_alg2alg(4*VOLUME,(*mdfs).frc,wfd[0]);

   set_ud_phase();
   force0(c);
   muladd_assign_alg(4*VOLUME,-1.0,(*mdfs).frc,wfd[0]);
   dev=norm_square_alg(4*VOLUME,0,wfd[0])/
      norm_square_alg(4*VOLUME,0,(*mdfs).frc);
   release_wfd();

   return sqrt(dev);
}


int main(int argc,char *argv[])
{
   int my_rank,k,ie,bc;
   double c,eps,act0,act1,dact,dsdt;
   double dev_frc,sig_loss,rdmy;
   double phi[2],phi_prime[2],theta[3];
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);

      printf("\n");
      printf("Check of the programs force0() and action0()\n");
      printf("--------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>]");
   }

   set_lat_parms(3.5,0.33,0,NULL,1.0);
   print_lat_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   theta[0]=0.38;
   theta[1]=-1.25;
   theta[2]=0.54;
   set_bc_parms(bc,0.9012,1.2034,1.0,1.0,phi,phi_prime,theta);
   print_bc_parms(3);

   start_ranlux(0,1234);
   geometry();
   alloc_wfd(1);
   c=0.789;
   dev_frc=chk_chs(c);

   if (my_rank==0)
      printf("Deviation of gauge force after a phase change = %.1e\n\n",
             dev_frc);

   for (k=0;k<4;k++)
   {
      random_ud();
      set_ud_phase();
      random_mom();
      dsdt=dSdt(c);

      eps=1.0e-4;
      rot_ud(eps);
      act0=2.0*action0(0)/3.0;
      rot_ud(-eps);

      rot_ud(-eps);
      act1=2.0*action0(0)/3.0;
      rot_ud(eps);

      rot_ud(2.0*eps);
      act0-=action0(0)/12.0;
      rot_ud(-2.0*eps);

      rot_ud(-2.0*eps);
      act1-=action0(0)/12.0;
      rot_ud(2.0*eps);

      act0*=c;
      act1*=c;

      dact=(act0-act1)/eps;
      dev_frc=dsdt-dact;
      sig_loss=-log10(fabs(1.0-act0/act1));

      rdmy=dsdt;
      MPI_Reduce(&rdmy,&dsdt,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&dsdt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      rdmy=dev_frc;
      MPI_Reduce(&rdmy,&dev_frc,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&dev_frc,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      rdmy=sig_loss;
      MPI_Reduce(&rdmy,&sig_loss,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&sig_loss,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      unset_ud_phase();
      ie=check_bc(0.0);
      error_root(ie!=1,1,"main [check3.c]",
                 "Operations did not preserve boundary conditions");

      if (my_rank==0)
      {
         printf("Relative deviation of dS/dt = %.2e ",fabs(dev_frc/dsdt));
         printf("[significance loss = %d digits]\n",(int)(sig_loss));
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
