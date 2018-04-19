
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2007, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Allocation and initialization of the global vector fields
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "vflds.h"
#include "global.h"
   
#define NFIELDS 7


int main(int argc,char *argv[])
{
   int my_rank,ie,k,ix;
   int bs[4],Ns;
   int nb,nbb,nv,nvec;
   complex **wv;
   complex_dble **wvd;
   dfl_parms_t dfl;
   FILE *fin=NULL,*flog=NULL;   

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);
      fin=freopen("check1.in","r",stdin);
      
      printf("\n");
      printf("Allocation and initialization of the global vector fields\n");
      printf("---------------------------------------------------------\n\n");

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
   
   start_ranlux(0,123456);   
   geometry();
   dfl=set_dfl_parms(bs,Ns);

   error((bs[0]!=dfl.bs[0])||(bs[1]!=dfl.bs[1])||
         (bs[2]!=dfl.bs[2])||(bs[3]!=dfl.bs[3])||(Ns!=dfl.Ns),1,
         "main [check1.c]","Parameter bs[4] or Ns are incorrectly set");

   alloc_wv(NFIELDS);
   alloc_wvd(NFIELDS);
   wv=reserve_wv(NFIELDS);
   wvd=reserve_wvd(NFIELDS);   

   nb=VOLUME/(bs[0]*bs[1]*bs[2]*bs[3]);
   nbb=(FACE0/(bs[1]*bs[2]*bs[3])+FACE1/(bs[0]*bs[2]*bs[3])+
        FACE2/(bs[0]*bs[1]*bs[3])+FACE3/(bs[0]*bs[1]*bs[2]));      
   nv=Ns*nb;
   nvec=Ns*(nb+nbb);
   ie=0;

   for (k=1;k<NFIELDS;k++)
   {
      if ((wv[k]-wv[0])!=(k*nvec))
         ie=1;
      if ((wvd[k]-wvd[0])!=(k*nvec))
         ie=2;      
   }

   error(ie==1,1,"main [check1.c]",
         "Field addresses reserved by wv_reserve() are incorrect");
   error(ie==2,1,"main [check1.c]",
         "Field addresses reserved by wvd_reserve() are incorrect");

   for (k=0;k<NFIELDS;k++)
   {
      for (ix=0;ix<nvec;ix++)
      {
         if ((wv[k][ix].re!=0.0f)||(wv[k][ix].im!=0.0f))
            ie=1;

         if ((wvd[k][ix].re!=0.0)||(wvd[k][ix].im!=0.0))
            ie=2;         
      }
   }

   error(ie==1,1,"main [check1.c]",
         "Fields reserved by wv_reserve() are not correctly initialized");
   error(ie==2,1,"main [check1.c]",
         "Fields reserved by wvd_reserve() are not correctly initialized");

   release_wv();
   release_wvd();

   wv=vflds();
   wvd=vdflds();

   for (k=1;k<(2*Ns);k++)
      if ((wv[k]-wv[0])!=(k*nv))
         ie=1;

   for (k=1;k<Ns;k++)
      if ((wvd[k]-wvd[0])!=(k*nv))
         ie=2;      
   
   error(ie==1,1,"main [check1.c]",
         "Field addresses returned by vflds() are incorrect");
   error(ie==2,1,"main [check1.c]",
         "Field addresses returned by vdflds() are incorrect");

   for (k=0;k<(2*Ns);k++)
   {
      for (ix=0;ix<nv;ix++)
         if ((wv[k][ix].re!=0.0f)||(wv[k][ix].im!=0.0f))
            ie=1;
   }

   for (k=0;k<Ns;k++)
   {
      for (ix=0;ix<nv;ix++)
         if ((wvd[k][ix].re!=0.0)||(wvd[k][ix].im!=0.0))
            ie=2;         
   }   

   error(ie==1,1,"main [check1.c]",
         "Fields allocated by vflds() are not correctly initialized");
   error(ie==2,1,"main [check1.c]",
         "Fields allocated by vdflds() are not correctly initialized");
   
   if (my_rank==0)
   {
      printf("No errors detected\n\n");
      fclose(flog);
   }
   
   MPI_Finalize();   
   exit(0);
}
