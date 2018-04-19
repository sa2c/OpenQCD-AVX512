
/*******************************************************************************
*
* File time1.c
*
* Copyright (C) 2007, 2008, 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of Awhat().
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
#include "uflds.h"
#include "sflds.h"
#include "vflds.h"
#include "linalg.h"
#include "dirac.h"
#include "dfl.h"
#include "little.h"
#include "global.h"


static void random_basis(int Ns)
{
   int i;
   spinor **ws;

   ws=reserve_ws(Ns);

   for (i=0;i<Ns;i++)
   {
      random_s(VOLUME,ws[i],1.0f);
      bnd_s2zero(ALL_PTS,ws[i]);
   }

   dfl_subspace(ws);
   release_ws();
}


int main(int argc,char *argv[])
{
   int my_rank,bc,count,nt;
   int i,nflds,ifail;
   int bs[4],Ns,nb,nbb,nv;
   double phi[2],phi_prime[2],theta[3];
   double mu,wt1,wt2,wdt;
   complex **wv;
   FILE *fin=NULL,*flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("time1.log","w",stdout);
      fin=freopen("check3.in","r",stdin);

      printf("\n");
      printf("Timing of Awhat()\n");
      printf("-----------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      read_line("Ns","%d",&Ns);
      fclose(fin);

      printf("bs = %d %d %d %d\n",bs[0],bs[1],bs[2],bs[3]);
      printf("Ns = %d\n\n",Ns);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [time1.c]",
                    "Syntax: time1 [-bc <type>]");
   }

   set_lat_parms(5.5,1.0,0,NULL,1.978);
   print_lat_parms();

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Ns,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   theta[0]=0.38;
   theta[1]=-1.25;
   theta[2]=0.54;
   set_bc_parms(bc,1.0,1.0,0.9012,1.2034,phi,phi_prime,theta);
   print_bc_parms(2);

   set_sw_parms(0.125);
   set_dfl_parms(bs,Ns);
   mu=0.0376;

   nb=VOLUME/(bs[0]*bs[1]*bs[2]*bs[3]);
   nbb=2*(FACE0/(bs[1]*bs[2]*bs[3])+
          FACE1/(bs[0]*bs[2]*bs[3])+
          FACE2/(bs[0]*bs[1]*bs[3])+
          FACE3/(bs[0]*bs[1]*bs[2]));
   nv=Ns*nb;

   start_ranlux(0,123456);
   geometry();

   alloc_ws(Ns);
   alloc_wvd(2);
   random_ud();
   set_ud_phase();
   random_basis(Ns);

   ifail=set_Awhat(mu);
   error(ifail!=0,1,"main [time1.c]","Inversion of Aee or Aoo failed");

   if (my_rank==0)
   {
      printf("Number of points = %d\n",VOLUME);
      printf("Number of blocks = %d\n",nb);
      printf("Number of points/block = %d\n",bs[0]*bs[1]*bs[2]*bs[3]);
      printf("Vector field size = %.2f KB\n",
             (double)(sizeof(complex)*nv)*1.0e-3);
      printf("Awhat array size = %.2f MB\n\n",
             (double)(sizeof(complex)*8*Ns*nv)*1.0e-6);
      fflush(flog);
   }

   nflds=(int)(1.0e6/(double)(sizeof(complex)*nv));
   if ((nflds%2)!=0)
      nflds+=1;
   if (nflds==0)
      nflds=2;

   alloc_wv(nflds);
   wv=reserve_wv(nflds);

   for (i=0;i<nflds;i++)
      random_v(nv,wv[i],1.0f);

   nt=(int)(1.0e3/(double)(8*Ns));
   if (nt<1)
      nt=1;
   wdt=0.0;

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<nt;count++)
      {
         for (i=0;i<nflds;i+=2)
            Awhat(wv[i],wv[i+1]);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      nt*=2;
   }

   wdt=4.0e6*wdt/(double)(nt*nflds);

   if (my_rank==0)
   {
      printf("Time per application of Awhat(), including communications:\n");
      printf("Total:     %4.3f msec\n",wdt*1.0e-3);
      printf("Per block: %4.3f usec",wdt/(double)(nb));
      printf(" (%d Mflops [%d bit arithmetic])\n",
             (int)(64.0*(double)(nb*Ns*Ns)/wdt),(int)(4*sizeof(complex)));
      printf("Per point: %4.3f usec\n\n",wdt/(double)(VOLUME));
      fflush(flog);
   }

   if (NPROC>1)
   {
      nt/=2;
      if (nt==0)
         nt=1;
      wdt=0.0;

      while (wdt<5.0)
      {
         for (i=0;i<nflds;i+=2)
            set_v2zero((nbb/2)*Ns,wv[i+1]+nv);

         MPI_Barrier(MPI_COMM_WORLD);
         wt1=MPI_Wtime();
         for (count=0;count<nt;count++)
         {
            for (i=0;i<nflds;i+=2)
            {
               cpv_int_bnd(wv[i]);
               cpv_ext_bnd(wv[i+1]);
            }
         }
         MPI_Barrier(MPI_COMM_WORLD);
         wt2=MPI_Wtime();

         wdt=wt2-wt1;
         nt*=2;
      }

      wdt=4.0e6*wdt/(double)(nt*nflds);

      if (my_rank==0)
      {
         printf("There are %d boundary blocks\n",nbb);
         printf("Time per application of Awhat() for the communications:\n");
         printf("Total:     %4.3f msec\n",wdt*1.0e-3);
         printf("Per block: %4.3f usec\n",wdt/(double)(nb));
         printf("Per point: %4.3f usec\n\n",wdt/(double)(VOLUME));
      }
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
