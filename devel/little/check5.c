
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2007, 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the program set_ltl_modes().
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


static double check_vd(int Ns,int nvh)
{
   int i,j;
   double d,dev;
   complex_dble **vd,z;

   vd=vdflds();
   dev=0.0;

   for (i=0;i<Ns;i++)
   {
      for (j=0;j<=i;j++)
      {
         z=vprod_dble(nvh,1,vd[i],vd[j]);

         if (i==j)
            z.re-=1.0;

         d=z.re*z.re+z.im*z.im;

         if (d>dev)
            dev=d;
      }
   }

   return sqrt(dev);
}


static double check_Awvd(int Ns,int nvh)
{
   int i;
   double d,dev;
   complex_dble **vd,**wvd,z;

   vd=vdflds();
   wvd=reserve_wvd(2);

   dev=0.0;
   z.re=-1.0;
   z.im=0.0;

   for (i=0;i<Ns;i++)
   {
      assign_vd2vd(nvh,vd[i],wvd[0]);
      Awhat_dble(wvd[0],wvd[1]);
      mulc_vadd_dble(nvh,wvd[1],vd[i]+nvh,z);
      d=vnorm_square_dble(nvh,1,wvd[1])/vnorm_square_dble(nvh,1,vd[i]+nvh);

      if (d>dev)
         dev=d;
   }

   release_wvd();

   return sqrt(dev);
}


static double check_ltl_matrix(int Ns,int nvh)
{
   int i,j,ie;
   double dev;
   complex_dble **vd,*amat,*bmat,*cmat,z;

   vd=vdflds();
   amat=ltl_matrix();
   bmat=amalloc(2*Ns*Ns*sizeof(*amat),ALIGN);
   error(bmat==NULL,1,"check_ltl_matrix [check5.c]",
         "Unable to allocate auxiliary arrays");
   cmat=bmat+Ns*Ns;

   for (i=0;i<Ns;i++)
   {
      for (j=0;j<Ns;j++)
         bmat[i*Ns+j]=vprod_dble(nvh,1,vd[i],vd[j]+nvh);
   }

   cmat_mul_dble(Ns,amat,bmat,cmat);
   dev=0.0;

   for (i=0;i<Ns;i++)
   {
      for (j=0;j<Ns;j++)
      {
         z.re=cmat[i*Ns+j].re;
         z.im=cmat[i*Ns+j].im;

         if (i==j)
            z.re-=1.0;

         dev+=(z.re*z.re+z.im*z.im);
      }
   }

   assign_vd2vd(Ns*Ns,amat,bmat);
   MPI_Bcast((double*)(amat),2*Ns*Ns,MPI_DOUBLE,0,MPI_COMM_WORLD);

   ie=0;

   for (i=0;i<(Ns*Ns);i++)
      if ((amat[i].re!=bmat[i].re)||(amat[i].im!=bmat[i].im))
         ie=1;

   error(ie!=0,1,"check_ltl_matrix [check5.c]",
         "Little matrix is not globally the same");

   return sqrt(dev)/(double)(Ns);
}


static double check_vflds(int Ns,int nvh)
{
   int i;
   double d,dev;
   complex **v;
   complex_dble **vd,**wvd,z;

   z.re=-1.0;
   z.im=0.0;
   dev=0.0;

   v=vflds();
   vd=vdflds();
   wvd=reserve_wvd(1);

   for (i=0;i<Ns;i++)
   {
      assign_v2vd(nvh,v[i],wvd[0]);
      mulc_vadd_dble(nvh,wvd[0],vd[i],z);
      d=vnorm_square_dble(nvh,1,wvd[0]);
      d/=vnorm_square_dble(nvh,1,vd[i]);
      if (d>dev)
         dev=d;
   }

   for (i=0;i<Ns;i++)
   {
      assign_v2vd(nvh,v[i]+nvh,wvd[0]);
      mulc_vadd_dble(nvh,wvd[0],vd[i]+nvh,z);
      d=vnorm_square_dble(nvh,1,wvd[0]);
      d/=vnorm_square_dble(nvh,1,vd[i]+nvh);
      if (d>dev)
         dev=d;
   }

   release_wvd();

   return sqrt(dev);
}


static double check_mds(int Ns,int nvh)
{
   int nv,k,l;
   double d,dev;
   complex **vs;
   complex_dble **vd,**wvd;
   spinor **mds,**ws;

   nv=2*nvh;
   mds=reserve_ws(Ns);
   ws=reserve_ws(1);
   vs=vflds();
   dev=0.0;

   for (k=0;k<Ns;k++)
   {
      dfl_v2s(vs[Ns+k],ws[0]);
      mulr_spinor_add(VOLUME,ws[0],mds[k],-1.0);

      d=(double)(norm_square(VOLUME,1,ws[0])/vnorm_square(nv,1,vs[Ns+k]));
      if (d>dev)
         dev=d;
   }

   release_ws();
   release_ws();

   vd=vdflds();
   wvd=reserve_wvd(1);

   for (k=0;k<Ns;k++)
   {
      assign_v2vd(nvh,vs[Ns+k],wvd[0]);

      for (l=0;l<Ns;l++)
         vproject_dble(nvh,1,wvd[0],vd[l]);

      d=vnorm_square_dble(nvh,1,wvd[0]);
      d/=(double)(vnorm_square(nvh,1,vs[Ns+k]));
      if (d>dev)
         dev=d;
   }

   release_wvd();

   return sqrt(dev);
}


int main(int argc,char *argv[])
{
   int my_rank,bc,ifail;
   int bs[4],Ns,nb,nvh;
   double phi[2],phi_prime[2],theta[3];
   double mu,dev;
   FILE *fin=NULL,*flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check5.log","w",stdout);
      fin=freopen("check3.in","r",stdin);

      printf("\n");
      printf("Check of the program set_ltl_modes()\n");
      printf("------------------------------------\n\n");

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
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check5.c]",
                    "Syntax: check5 [-bc <type>]");
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

   start_ranlux(0,123456);
   geometry();

   alloc_ws(Ns+1);
   alloc_wvd(3);

   nb=VOLUME/(bs[0]*bs[1]*bs[2]*bs[3]);
   nvh=Ns*(nb/2);

   random_ud();
   set_ud_phase();
   random_basis(Ns);
   ifail=set_Awhat(mu);
   error_root(ifail!=0,1,"main [check5.c]",
              "Computation of the little Dirac operator failed");

   if (my_rank==0)
      printf("Maximal relative deviations found:\n\n");

   dev=check_vd(Ns,nvh);

   if (my_rank==0)
      printf("Orthonormality of vdflds: %.2e\n",dev);

   dev=check_Awvd(Ns,nvh);

   if (my_rank==0)
      printf("Awhat*vdflds:             %.2e\n",dev);

   dev=check_ltl_matrix(Ns,nvh);

   if (my_rank==0)
      printf("Little-little matrix:     %.2e\n\n",dev);

   dev=check_vflds(Ns,nvh);

   if (my_rank==0)
      printf("Single-precision fields:  %.2e\n",dev);

   dev=check_mds(Ns,nvh);

   if (my_rank==0)
   {
      printf("Global deflation modes:   %.2e\n\n",dev);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
