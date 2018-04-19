
/*******************************************************************************
*
* File check8.c
*
* Copyright (C) 2012, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check and performance of the multi-shift CG solver.
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
#include "archive.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "forces.h"
#include "global.h"

static int my_rank,bc,first,last,step;
static int nmu,nmx;
static double kappa,csw,*mu,cF,cF_prime;
static double uphi[2],uphi_prime[2],theta[3],m0,*res;
static char cnfg_dir[NAME_SIZE],cnfg_file[NAME_SIZE],nbase[NAME_SIZE];


int main(int argc,char *argv[])
{
   int nsize,icnfg,status,k,ie;
   double nrm,del;
   double wt1,wt2,wdt;
   spinor_dble *eta,*chi,*phi,**psi,**wsd,**rsd;
   lat_parms_t lat;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check8.log","w",stdout);
      fin=freopen("check8.in","r",stdin);

      printf("\n");
      printf("Check and performance of the multi-shift CG solver\n");
      printf("--------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      find_section("Configurations");
      read_line("name","%s",nbase);
      read_line("cnfg_dir","%s",cnfg_dir);
      read_line("first","%d",&first);
      read_line("last","%d",&last);
      read_line("step","%d",&step);

      find_section("Lattice parameters");
      read_line("kappa","%lf",&kappa);
      read_line("csw","%lf",&csw);
      nmu=count_tokens("mu");

      find_section("Boundary conditions");
      read_line("type","%d",&bc);

      uphi[0]=0.0;
      uphi[1]=0.0;
      uphi_prime[0]=0.0;
      uphi_prime[1]=0.0;
      cF=1.0;
      cF_prime=1.0;

      if (bc==1)
         read_dprms("uphi",2,uphi);

      if ((bc==1)||(bc==2))
         read_dprms("uphi'",2,uphi_prime);

      if (bc!=3)
         read_line("cF","%lf",&cF);

      if (bc==2)
         read_line("cF'","%lf",&cF_prime);
      else
         cF_prime=cF;

      read_dprms("theta",3,theta);
   }

   MPI_Bcast(nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(&first,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&last,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&step,1,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&csw,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmu,1,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(uphi,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(uphi_prime,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cF,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cF_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(theta,3,MPI_DOUBLE,0,MPI_COMM_WORLD);

   mu=malloc(2*nmu*sizeof(*mu));
   error(mu==NULL,1,"main [check8.c]","Unable to allocate auxiliary arrays");
   res=mu+nmu;

   if (my_rank==0)
   {
      find_section("Lattice parameters");
      read_dprms("mu",nmu,mu);

      find_section("CG");
      read_line("nmx","%d",&nmx);
      error_root(nmu!=count_tokens("res"),1,"main [check8.c]",
                 "The numbers of twisted masses and residues do not match");
      read_dprms("res",nmu,res);

      fclose(fin);
   }

   MPI_Bcast(mu,nmu,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(res,nmu,MPI_DOUBLE,0,MPI_COMM_WORLD);

   lat=set_lat_parms(5.5,1.0,1,&kappa,csw);
   print_lat_parms();

   set_bc_parms(bc,1.0,1.0,cF,cF_prime,uphi,uphi_prime,theta);
   print_bc_parms(2);

   start_ranlux(0,1234);
   geometry();

   m0=lat.m0[0];
   set_sw_parms(m0);

   if (my_rank==0)
   {
      printf("mu = %.6f",mu[0]);
      for (k=1;k<nmu;k++)
         printf(", %.6f",mu[k]);
      printf("\n\n");

      printf("CG parameters:\n");
      printf("nmx = %d\n",nmx);
      printf("res = %.2e",res[0]);
      for (k=1;k<nmu;k++)
         printf(", %.2e",res[k]);
      printf("\n\n");

      printf("Configurations %sn%d -> %sn%d in steps of %d\n\n",
             nbase,first,nbase,last,step);
      fflush(flog);
   }

   if (nmu==1)
      alloc_wsd(8);
   else
      alloc_wsd(5+2*nmu);

   wsd=reserve_wsd(2);
   eta=wsd[0];
   chi=wsd[1];
   psi=reserve_wsd(nmu);

   error_root(((last-first)%step)!=0,1,"main [check8.c]",
              "last-first is not a multiple of step");
   check_dir_root(cnfg_dir);
   nsize=name_size("%s/%sn%d",cnfg_dir,nbase,last);
   error_root(nsize>=NAME_SIZE,1,"main [check8.c]",
              "configuration file name is too long");
   ie=0;

   for (icnfg=first;icnfg<=last;icnfg+=step)
   {
      sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,icnfg);
      import_cnfg(cnfg_file);

      if (my_rank==0)
      {
         printf("Configuration no %d\n\n",icnfg);
         fflush(flog);
      }

      set_ud_phase();
      random_sd(VOLUME,eta,1.0);
      bnd_sd2zero(ALL_PTS,eta);
      nrm=sqrt(norm_square_dble(VOLUME/2,1,eta));
      assign_sd2sd(VOLUME,eta,chi);

      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      tmcgm(nmx,res,nmu,mu,eta,psi,&status);

      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();
      wdt=wt2-wt1;

      mulr_spinor_add_dble(VOLUME,chi,eta,-1.0);
      del=norm_square_dble(VOLUME,1,chi);
      error_root(del!=0.0,1,"main [check8.c]",
                 "Source field is not preserved");

      if (my_rank==0)
      {
         printf("status = %d\n",status);
         printf("time = %.2e sec (total)\n",wdt);
         if (status>0)
            printf("     = %.2e usec (per point and CG iteration)\n",
                   (1.0e6*wdt)/((double)(status)*(double)(VOLUME)));
         fflush(flog);
         error_root(status<0,1,"main [check8.c]",
                    "Solver did not converge");
         printf("residues = ");
      }

      rsd=reserve_wsd(1);
      phi=rsd[0];
      status=0;

      for (k=0;k<nmu;k++)
      {
         Dwhat_dble(mu[k],psi[k],chi);
         mulg5_dble(VOLUME/2,chi);
         Dwhat_dble(-mu[k],chi,phi);
         mulg5_dble(VOLUME/2,phi);
         mulr_spinor_add_dble(VOLUME/2,phi,eta,-1.0);
         del=sqrt(norm_square_dble(VOLUME/2,1,phi))/nrm;

         if (del<res[k])
            status+=1;

         if (my_rank==0)
         {
            if (k==0)
               printf("%.2e",del);
            else
               printf(", %.2e",del);
         }
      }

      release_wsd();
      ie+=(status<nmu);

      if (my_rank==0)
      {
         printf("\n");

         if (status==nmu)
            printf("All residues are as required\n\n");
         else
            printf("ERROR: %d residues are too large\n\n",nmu-status);

         fflush(flog);
      }
   }

   if (my_rank==0)
   {
      if (ie==0)
         printf("No errors detected --- all seems fine!\n\n");
      else
         printf("ERROR: the residues are too large (%d configurations)\n\n",ie);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
