
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2005, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the programs in the module vinit.c
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
#include "linalg.h"
#include "vflds.h"
#include "global.h"

#define NFLDS 5

static float sig[NFLDS];
static double sigd[NFLDS];


int main(int argc,char *argv[])
{
   int my_rank,ie,k,ix;
   int bs[4],Ns,nb,nv;
   double var,var_all,d,dmax;
   complex z;
   complex_dble zd;
   complex **wv;
   complex_dble **wvd;
   FILE *fin=NULL,*flog=NULL;  

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout); 
      fin=freopen("check1.in","r",stdin);
      
      printf("\n");
      printf("Check of the programs in the module vinit\n");
      printf("-----------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      read_line("Ns","%d",&Ns);
      fclose(fin);

      printf("bs = %d %d %d %d\n",bs[0],bs[1],bs[2],bs[3]);
      printf("Ns = %d\n\n",Ns);
      fflush(flog);      
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Ns,1,MPI_INT,0,MPI_COMM_WORLD);
   
   start_ranlux(0,12345);
   geometry();
   set_dfl_parms(bs,Ns);
   
   alloc_wv(2*NFLDS);
   alloc_wvd(2*NFLDS);
   wv=reserve_wv(2*NFLDS);
   wvd=reserve_wvd(2*NFLDS);   

   nb=VOLUME/(bs[0]*bs[1]*bs[2]*bs[3]);
   nv=Ns*nb;
   z.im=0.0f;
   zd.im=0.0;
   ie=0;

   if (my_rank==0)
   {
      printf("Choose random single-precision fields\n");   
      ranlxs(sig,NFLDS);
   }

   MPI_Bcast(sig,NFLDS,MPI_FLOAT,0,MPI_COMM_WORLD);

   for (k=0;k<NFLDS;k++)
   {
      random_v(nv,wv[k],sig[k]);
      var=0.0;

      for (ix=0;ix<nv;ix++)
         var+=(double)((wv[k][ix].re*wv[k][ix].re+
                        wv[k][ix].im*wv[k][ix].im));

      MPI_Reduce(&var,&var_all,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      
      if (my_rank==0)
      {
         var_all/=((double)(NPROC)*(double)(nv));               
         printf("<s[%d]^2> = %.4e (sigma^2 = %.4e)\n",
                k,var_all,sig[k]*sig[k]);
      }
   }

   if (my_rank==0)
   {
      printf("\n");
      printf("Choose random double-precision fields\n");   
      ranlxd(sigd,NFLDS);
   }

   MPI_Bcast(sigd,NFLDS,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   for (k=0;k<NFLDS;k++)
   {
      random_vd(nv,wvd[k],sigd[k]);
      var=0.0;

      for (ix=0;ix<nv;ix++)
         var+=(wvd[k][ix].re*wvd[k][ix].re+wvd[k][ix].im*wvd[k][ix].im);

      MPI_Reduce(&var,&var_all,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      
      if (my_rank==0)
      {
         var_all/=((double)(NPROC)*(double)(nv));
         printf("<sd[%d]^2> = %.4e (sigma^2 = %.4e)\n",
                k,var_all,sigd[k]*sigd[k]);
      }
   }

   for (k=0;k<NFLDS;k++)
   {
      random_v(nv,wv[k],1.0f);
      random_vd(nv,wvd[k],1.0);
      assign_v2v(nv,wv[k],wv[k+NFLDS]);
      assign_vd2vd(nv,wvd[k],wvd[k+NFLDS]);

      z.re=-1.0f;
      zd.re=-1.0;
      mulc_vadd(nv,wv[k],wv[k+NFLDS],z);
      mulc_vadd_dble(nv,wvd[k],wvd[k+NFLDS],zd);      
      
      for (ix=0;ix<nv;ix++)
      {
         if ((wv[k][ix].re!=0.0f)||(wv[k][ix].im!=0.0f))
            ie=1;
         if ((wvd[k][ix].re!=0.0)||(wvd[k][ix].im!=0.0))
            ie=2;
      }
   }

   error(ie==1,1,"main [check2.c]","assign_v2v() is incorrect");
   error(ie==2,1,"main [check2.c]","assign_vd2vd() is incorrect");

   for (k=0;k<NFLDS;k++)
   {
      random_v(nv,wv[k],1.0f);
      assign_v2vd(nv,wv[k],wvd[k]);
      assign_vd2v(nv,wvd[k],wv[k+NFLDS]);

      z.re=-1.0f;
      mulc_vadd(nv,wv[k],wv[k+NFLDS],z);

      for (ix=0;ix<nv;ix++)
      {
         if ((wv[k][ix].re!=0.0f)||(wv[k][ix].im!=0.0f))
            ie=1;
      }
   }

   error(ie==1,1,"main [check2.c]",
         "assign_v2vd() or assign_vd2v() is incorrect");

   dmax=0.0;
   
   for (k=0;k<NFLDS;k++)
   {
      random_v(nv,wv[k],1.0f);
      random_vd(nv,wvd[k],1.0);
      assign_vd2vd(nv,wvd[k],wvd[k+NFLDS]);

      add_v2vd(nv,wv[k],wvd[k]);
      d=vnorm_square_dble(nv,1,wvd[k]);
      zd.re=-1.0;
      mulc_vadd_dble(nv,wvd[k],wvd[k+NFLDS],zd);
      assign_v2vd(nv,wv[k],wvd[k+NFLDS]);
      mulc_vadd_dble(nv,wvd[k],wvd[k+NFLDS],zd);

      d=vnorm_square_dble(nv,1,wvd[k])/d;
      if (d>dmax)
         dmax=d;
   }

   if (my_rank==0)
   {
      printf("\n");
      printf("Relative deviations (should be less than 1.0e-7 or so):\n");
      printf("add_v2vd():  %.1e\n",sqrt(dmax));
   }
   
   dmax=0.0;
   
   for (k=0;k<NFLDS;k++)
   {
      random_vd(nv,wvd[k],1.0);
      random_vd(nv,wvd[k+NFLDS],1.0);

      diff_vd2v(nv,wvd[k],wvd[k+NFLDS],wv[k]);
      zd.re=-1.0;
      mulc_vadd_dble(nv,wvd[k],wvd[k+NFLDS],zd);
      d=vnorm_square_dble(nv,1,wvd[k]);
      assign_v2vd(nv,wv[k],wvd[k+NFLDS]);
      mulc_vadd_dble(nv,wvd[k],wvd[k+NFLDS],zd);

      d=vnorm_square_dble(nv,1,wvd[k])/d;
      if (d>dmax)
         dmax=d;
   }

   if (my_rank==0)
   {
      printf("diff_vd2v(): %.1e\n\n",sqrt(dmax));
      fclose(flog);
   }

   MPI_Finalize();   
   exit(0);
}
