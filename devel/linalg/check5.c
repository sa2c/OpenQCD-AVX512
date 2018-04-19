
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2007, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Checks on the programs in the module valg_dble
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
#include "vflds.h"
#include "linalg.h"
#include "global.h"

static complex_dble v[25];
static complex_dble *ppk[5];


static complex_dble sp(int vol,complex_dble *pk,complex_dble *pl)
{
   int ix;
   double x,y;
   complex_dble z;

   x=0.0;
   y=0.0;

   for (ix=0;ix<vol;ix++)
   {
      x+=(double)((*pk).re*(*pl).re+(*pk).im*(*pl).im);
      y+=(double)((*pk).re*(*pl).im-(*pk).im*(*pl).re);
      pk+=1;
      pl+=1;
   }

   z.re=x;
   z.im=y;

   return z;
}


int main(int argc,char *argv[])
{
   int my_rank,i,j,vol,off;
   int bs[4],Ns,nb,nv;
   int icom,ieo;
   double r,zsq;
   double d,dmax,dall;
   complex_dble w,z;
   complex_dble **wvd,*pk,*pl;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check5.log","w",stdout);
      fin=freopen("check4.in","r",stdin);

      printf("\n");
      printf("Checks on the programs in the module valg_dble\n");
      printf("----------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      fclose(fin);

      printf("bs = %d %d %d %d\n\n",bs[0],bs[1],bs[2],bs[3]);
      fflush(flog);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);

   start_ranlux(0,12345);
   geometry();

   Ns=4;
   set_dfl_parms(bs,Ns);
   nb=VOLUME/(bs[0]*bs[1]*bs[2]*bs[3]);
   nv=Ns*nb;

   alloc_wvd(10);
   wvd=reserve_wvd(10);
   dall=0.0;

   for (icom=0;icom<2;icom++)
   {
      if ((icom==0)||(NPROC>1))
      {
         if (my_rank==0)
         {
            if (icom==1)
            {
               printf("Checks with global summation\n");
               printf("============================\n\n");
            }
            else
            {
               printf("Checks without global summation\n");
               printf("===============================\n\n");
            }
         }

         for (ieo=0;ieo<3;ieo++)
         {
            if (my_rank==0)
            {
               if (ieo==0)
                  printf("First case: full lattice\n\n");
               else if (ieo==1)
                  printf("Second case: even points\n\n");
               else
                  printf("Third case: odd points\n\n");
            }

            vol=nv/2;
            off=0;

            if (ieo==0)
               vol=nv;
            if (ieo==2)
               off=nv/2;

            for (i=0;i<10;i++)
               random_vd(vol,wvd[i]+off,1.0f);

            dmax=0.0;

            for (i=0;i<10;i++)
            {
               pk=wvd[i]+off;
               pl=wvd[9-i]+off;

               if (icom==1)
               {
                  z=sp(vol,pk,pl);
                  MPI_Reduce(&z.re,&w.re,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
                  MPI_Bcast(&w.re,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
               }
               else
                  w=sp(vol,pk,pl);

               z=vprod_dble(vol,icom,pk,pl);
               r=vnorm_square_dble(vol,icom,pk)*vnorm_square_dble(vol,icom,pl);
               d=(z.re-w.re)*(z.re-w.re)+(z.im-w.im)*(z.im-w.im);
               d=sqrt(d/r);
               if (d>dmax)
                  dmax=d;

               z=vprod_dble(vol,icom,pk,pk);
               r=vnorm_square_dble(vol,icom,pk);

               d=fabs(z.im/r);
               if (d>dmax)
                  dmax=d;

               d=fabs(z.re/r-1.0);
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Check of vprod_dble and vnorm_square_dble: %.2e\n\n",
                      dmax);
            }

            dmax=0.0;
            z.re= 0.345;
            z.im=-0.876;
            zsq=z.re*z.re+z.im*z.im;

            for (i=0;i<9;i++)
            {
               pk=wvd[i]+off;
               pl=wvd[i+1]+off;

               w=vprod_dble(vol,icom,pk,pl);
               r=vnorm_square_dble(vol,icom,pk)+
                  zsq*vnorm_square_dble(vol,icom,pl)
                  +2.0f*(z.re*w.re-z.im*w.im);
               mulc_vadd_dble(vol,pk,pl,z);

               d=fabs(r/vnorm_square_dble(vol,icom,pk)-1.0);
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Consistency of vprod_dble, vnorm_square_dble\n");
               printf("and mulc_vadd_dble: %.2e\n\n",dmax);
            }

            for (i=0;i<10;i++)
               random_vd(vol,wvd[i]+off,1.0f);

            dmax=0.0;

            for (i=0;i<10;i++)
            {
               pk=wvd[i]+off;

               if (i>0)
               {
                  pl=wvd[i-1]+off;
                  vproject_dble(vol,icom,pk,pl);
                  z=vprod_dble(vol,icom,pk,pl);

                  d=(fabs(z.re)+fabs(z.im))/
                     sqrt(vnorm_square_dble(vol,icom,pk));

                  if (d>dmax)
                     dmax=d;
               }

               vnormalize_dble(vol,icom,pk);
               r=vnorm_square_dble(vol,icom,pk);

               d=fabs(r-1.0);
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Consistency of vprod_dble, vnorm_square_dble,\n");
               printf("vnormalize_dble and vproject_dble: %.2e\n\n",dmax);
            }

            for (i=0;i<5;i++)
            {
               pk=wvd[i]+off;
               pl=wvd[i+5]+off;

               random_vd(vol,wvd[i]+off,1.0f);
               assign_vd2vd(vol,pk,pl);

               for (j=0;j<5;j++)
               {
                  v[5*i+j].re=0.1234*(double)(i^2)-0.8976*(double)(j);
                  v[5*i+j].im=0.2231*(double)(i)+0.9922*(double)(j^2);
               }

               ppk[i]=pl;
            }

            vrotate_dble(vol,5,ppk,v);
            dmax=0.0;

            for (i=5;i<10;i++)
            {
               pk=wvd[i]+off;

               for (j=0;j<5;j++)
               {
                  z.re=-v[5*j+(i-5)].re;
                  z.im=-v[5*j+(i-5)].im;

                  pl=wvd[j]+off;
                  mulc_vadd_dble(vol,pk,pl,z);
               }

               r=vnorm_square_dble(vol,icom,pk);

               d=fabs(r);
               if (d>dmax)
                  dmax=d;
            }

            dmax/=vnorm_square_dble(vol,icom,wvd[0]+off);
            dmax=sqrt(dmax);

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Consistency of mulc_vadd_dble\n");
               printf("and vrotate_dble: %.2e\n\n",dmax);
            }
         }
      }
   }

   if (my_rank==0)
   {
      printf("Maximal deviation in all tests: %.2e\n\n",dall);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
