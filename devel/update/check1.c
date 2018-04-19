
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of set_mdsteps().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "mdflds.h"
#include "update.h"
#include "global.h"

static int my_rank;
static force_t force[]={FRG,FRF_TM1,FRF_TM1_EO,FRF_TM1_EO_SDET,
                        FRF_TM2,FRF_TM2_EO,FRF_RAT,FRF_RAT_SDET};


static void read_hmc_parms(void)
{
   int nlv;
   double tau;

   if (my_rank==0)
   {
      find_section("HMC parameters");
      read_line("nlv","%d",&nlv);
      read_line("tau","%lf",&tau);
   }

   MPI_Bcast(&nlv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&tau,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   set_hmc_parms(0,NULL,0,0,NULL,nlv,tau);
}


static void read_integrator(void)
{
   int nlv,i,j,k,idf;
   int irat[3],imu[4],isp[4],ncr[4];
   hmc_parms_t hmc;
   mdint_parms_t mdp;
   force_parms_t fp;
   char line[NAME_SIZE];

   for (i=0;i<3;i++)
      irat[i]=0;

   for (i=0;i<4;i++)
   {
      imu[i]=0;
      isp[i]=0;
      ncr[i]=0;
   }

   hmc=hmc_parms();
   nlv=hmc.nlv;

   for (i=0;i<nlv;i++)
   {
      read_mdint_parms(i);
      mdp=mdint_parms(i);

      for (j=0;j<mdp.nfr;j++)
      {
         k=mdp.ifr[j];
         fp=force_parms(k);

         if (fp.force==FORCES)
         {
            if (my_rank==0)
            {
               sprintf(line,"Force %d",k);
               find_section(line);
               read_line("force","%s",line);

               if (strcmp(line,"FRG")==0)
                  idf=0;
               else if (strcmp(line,"FRF_TM1")==0)
                  idf=1;
               else if (strcmp(line,"FRF_TM1_EO")==0)
                  idf=2;
               else if (strcmp(line,"FRF_TM1_EO_SDET")==0)
                  idf=3;
               else if (strcmp(line,"FRF_TM2")==0)
                  idf=4;
               else if (strcmp(line,"FRF_TM2_EO")==0)
                  idf=5;
               else if (strcmp(line,"FRF_RAT")==0)
                  idf=6;
               else if (strcmp(line,"FRF_RAT_SDET")==0)
                  idf=7;
               else
                  error_root(1,1,"read_integrator [check1.c]",
                             "Unknown force %s",line);
            }

            MPI_Bcast(&idf,1,MPI_INT,0,MPI_COMM_WORLD);
            set_force_parms(k,force[idf],0,0,irat,imu,isp,ncr);
         }
      }
   }
}


int main(int argc,char *argv[])
{
   int i,ie;
   int nop,itu,*iop;
   double phi[2],phi_prime[2],theta[3];
   double kappa,*eps;
   mdstep_t *mds;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);
      fin=freopen("check1.in","r",stdin);

      printf("\n");
      printf("Check of set_mdsteps()\n");
      printf("----------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

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

   read_hmc_parms();
   read_integrator();

   if (my_rank==0)
      fclose(fin);

   start_ranlux(0,1234);
   geometry();
   set_mdsteps();
   print_mdsteps(0x6);

   mds=mdsteps(&nop,&itu);
   i=0;

   while (mds[i].iop<=itu)
      i+=1;

   error((mds[i].iop!=(itu+1))||(i!=(nop-1)),1,"main [check1.c]",
         "Parameters nop or itu returned by mdsteps are incorrect");

   iop=malloc(nop*sizeof(*iop));
   eps=malloc(nop*sizeof(*eps));
   error((iop==NULL)||(eps==NULL),1,"main [check1.c]",
         "Unable to allocate auxiliary arrays");

   for (i=0;i<nop;i++)
   {
      iop[i]=mds[i].iop;
      eps[i]=mds[i].eps;
   }

   MPI_Bcast(iop,nop,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(eps,nop,MPI_DOUBLE,0,MPI_COMM_WORLD);

   ie=0;

   for (i=0;i<nop;i++)
   {
      ie|=(iop[i]!=mds[i].iop);
      ie|=(eps[i]!=mds[i].eps);
   }

   error(ie!=0,1,"main [check1.c]",
         "Integration steps are not globally the same");

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
