
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check and performance of the deflated SAP+GCR solver.
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
#include "dirac.h"
#include "dfl.h"
#include "global.h"

int my_rank,id,first,last,step;
int bs_sap[4],nmr_sap,ncy_sap,nkv_gcr,nmx_gcr;
int bs_dfl[4],Ns,nkv_dfl,nmx_dfl,nkv_dpr,nmx_dpr,eoflg,bc;
int ninv_dgn,nmr_dgn,ncy_dgn;
double kappa,csw,mu,cF,cF_prime;
double phi[2],phi_prime[2],theta[3],m0,res_gcr,res_dpr;
double kappa_dgn,mu_dgn;
char cnfg_dir[NAME_SIZE],cnfg_file[NAME_SIZE],nbase[NAME_SIZE];


int main(int argc,char *argv[])
{
   int nsize,icnfg,ncnfg;
   int status[2],avgstat[2];
   double rho,nrm,del,resm;
   double wt1,wt2,wdt,wta;
   spinor_dble **psd;
   lat_parms_t lat;
   tm_parms_t tm;
   dfl_pro_parms_t dpr;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);
      fin=freopen("check4.in","r",stdin);

      printf("\n");
      printf("Check and performance of the deflated SAP+GCR solver\n");
      printf("----------------------------------------------------\n\n");

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
      read_line("mu","%lf",&mu);
      read_line("eoflg","%d",&eoflg);

      find_section("Boundary conditions");
      read_line("type","%d",&bc);

      phi[0]=0.0;
      phi[1]=0.0;
      phi_prime[0]=0.0;
      phi_prime[1]=0.0;
      cF=1.0;
      cF_prime=1.0;

      if (bc==1)
         read_dprms("phi",2,phi);

      if ((bc==1)||(bc==2))
         read_dprms("phi'",2,phi_prime);

      if (bc!=3)
         read_line("cF","%lf",&cF);

      if (bc==2)
         read_line("cF'","%lf",&cF_prime);
      else
         cF_prime=cF;

      read_dprms("theta",3,theta);

      find_section("SAP");
      read_line("bs","%d %d %d %d",bs_sap,bs_sap+1,bs_sap+2,bs_sap+3);
      read_line("nmr","%d",&nmr_sap);
      read_line("ncy","%d",&ncy_sap);

      find_section("Deflation subspace");
      read_line("bs","%d %d %d %d",bs_dfl,bs_dfl+1,bs_dfl+2,bs_dfl+3);
      read_line("Ns","%d",&Ns);

      find_section("Deflation subspace generation");
      read_line("kappa","%lf",&kappa_dgn);
      read_line("mu","%lf",&mu_dgn);
      read_line("ninv","%d",&ninv_dgn);
      read_line("nmr","%d",&nmr_dgn);
      read_line("ncy","%d",&ncy_dgn);

      find_section("Deflation projection");
      read_line("nkv","%d",&nkv_dpr);
      read_line("nmx","%d",&nmx_dpr);
      read_line("res","%lf",&res_dpr);

      find_section("GCR");
      read_line("nkv","%d",&nkv_gcr);
      read_line("nmx","%d",&nmx_gcr);
      read_line("res","%lf",&res_gcr);

      fclose(fin);
   }

   MPI_Bcast(nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(&first,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&last,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&step,1,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&csw,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cF,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&mu,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&eoflg,1,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(phi,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(phi_prime,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cF,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cF_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(theta,3,MPI_DOUBLE,0,MPI_COMM_WORLD);

   MPI_Bcast(bs_sap,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmr_sap,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ncy_sap,1,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(bs_dfl,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Ns,1,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(&kappa_dgn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&mu_dgn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&ninv_dgn,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmr_dgn,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ncy_dgn,1,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(&nkv_dpr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmx_dpr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&res_dpr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   MPI_Bcast(&nkv_gcr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmx_gcr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&res_gcr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   lat=set_lat_parms(5.5,1.0,1,&kappa,csw);
   print_lat_parms();

   set_bc_parms(bc,1.0,1.0,cF,cF_prime,phi,phi_prime,theta);
   print_bc_parms(3);

   set_sap_parms(bs_sap,1,nmr_sap,ncy_sap);
   m0=lat.m0[0];
   set_sw_parms(m0);
   tm=set_tm_parms(eoflg);
   set_dfl_parms(bs_dfl,Ns);
   dpr=set_dfl_pro_parms(nkv_dpr,nmx_dpr,res_dpr);
   set_dfl_gen_parms(kappa_dgn,mu_dgn,ninv_dgn,nmr_dgn,ncy_dgn);

   if (my_rank==0)
   {
      printf("mu = %.6f\n",mu);
      printf("eoflg = %d\n\n",tm.eoflg);
   }

   print_sap_parms(1);
   print_dfl_parms(0);

   if (my_rank==0)
   {
      printf("GCR parameters:\n");
      printf("nkv = %d\n",nkv_gcr);
      printf("nmx = %d\n",nmx_gcr);
      printf("res = %.2e\n\n",res_gcr);

      printf("Configurations %sn%d -> %sn%d in steps of %d\n\n",
             nbase,first,nbase,last,step);
      fflush(flog);
   }

   start_ranlux(0,1234);
   geometry();

   if (Ns<=(2*nkv_gcr))
      alloc_ws(2*nkv_gcr+2);
   else
      alloc_ws(Ns+2);
   alloc_wsd(6);
   alloc_wv(2*dpr.nkv+2);
   alloc_wvd(4);
   psd=reserve_wsd(3);

   error_root(((last-first)%step)!=0,1,"main [check4.c]",
              "last-first is not a multiple of step");

   nsize=name_size("%s/%sn%d",cnfg_dir,nbase,last);
   error_root(nsize>=NAME_SIZE,1,"main [check4.c]",
              "cnfg_dir name is too long");

   ncnfg=0;
   avgstat[0]=0;
   avgstat[1]=0;
   resm=0.0;
   wta=0.0;

   for (icnfg=first;icnfg<=last;icnfg+=step)
   {
      sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,icnfg);
      import_cnfg(cnfg_file);
      set_ud_phase();

      if (my_rank==0)
      {
         printf("Configuration no %d\n",icnfg);
         fflush(flog);
      }

      dfl_modes(status);
      error_root(status[0]<0,1,"main [check4.c]",
                 "Subspace generation failed");
      random_sd(VOLUME,psd[0],1.0);
      bnd_sd2zero(ALL_PTS,psd[0]);
      nrm=sqrt(norm_square_dble(VOLUME,1,psd[0]));
      assign_sd2sd(VOLUME,psd[0],psd[2]);

      rho=dfl_sap_gcr(nkv_gcr,nmx_gcr,res_gcr,mu,psd[0],psd[1],status);

      mulr_spinor_add_dble(VOLUME,psd[2],psd[0],-1.0);
      del=norm_square_dble(VOLUME,1,psd[2]);
      error_root(del!=0.0,1,"main [check4.c]",
                 "Source field is not preserved");

      Dw_dble(mu,psd[1],psd[2]);
      mulr_spinor_add_dble(VOLUME,psd[2],psd[0],-1.0);
      del=sqrt(norm_square_dble(VOLUME,1,psd[2]));

      if (my_rank==0)
      {
         printf("status = %d,%d\n",status[0],status[1]);
         printf("rho   = %.2e, res   = %.2e\n",rho,res_gcr);
         printf("check = %.2e, check = %.2e\n",del,del/nrm);
         fflush(flog);
      }

      if ((status[0]>=0)&&(status[1]>=0))
      {
         ncnfg+=1;
         avgstat[0]+=status[0];
         avgstat[1]+=status[1];
         del/=nrm;

         if (del>resm)
            resm=del;

         MPI_Barrier(MPI_COMM_WORLD);
         wt1=MPI_Wtime();

         rho=dfl_sap_gcr(nkv_gcr,nmx_gcr,res_gcr,mu,psd[0],psd[0],status);

         MPI_Barrier(MPI_COMM_WORLD);
         wt2=MPI_Wtime();
         wdt=wt2-wt1;
         wta+=wdt;

         if (my_rank==0)
         {
            printf("time = %.2e sec (w/o preparatory steps)\n",wdt);
            if (status[0]>0)
               printf("     = %.2e usec (per point and GCR iteration)",
                      (1.0e6*wdt)/((double)(status[0])*(double)(VOLUME)));
            printf("\n\n");
            fflush(flog);
         }

         mulr_spinor_add_dble(VOLUME,psd[0],psd[1],-1.0);
         del=norm_square_dble(VOLUME,1,psd[0]);
         error_root(del!=0.0,1,"main [check4.c]","Incorrect result when "
                    "the input and output fields coincide");
      }
   }

   if (my_rank==0)
   {
      printf("Summary of results\n");
      printf("------------------\n\n");

      printf("Processed %d configurations\n",ncnfg);
      printf("Solver failed in %d cases\n",(last-first)/step+1-ncnfg);
      printf("Maximal relative residue = %.1e\n",resm);

      status[0]=(avgstat[0]+ncnfg/2)/ncnfg;
      status[1]=(avgstat[1]+ncnfg/2)/ncnfg;
      wta/=(double)(ncnfg);

      printf("Average status = %d,%d\n",status[0],status[1]);
      printf("Average time = %.2e sec (w/o preparatory steps)\n",wta);
      if (status[0]>0)
         printf("             = %.2e usec (per point and GCR iteration)",
                (1.0e6*wta)/((double)(status[0])*(double)(VOLUME)));
      printf("\n\n");

      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
