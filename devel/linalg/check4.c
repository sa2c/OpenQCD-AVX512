
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2007, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency checks on the programs in the module valg
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


static complex v[25];
static complex *ppk[5];


static complex sp(int vol,complex *pk,complex *pl)
{
   int ix;
   double x,y;
   complex z;

   x=0.0;
   y=0.0;

   for (ix=0;ix<vol;ix++)
   {
      x+=(double)((*pk).re*(*pl).re+(*pk).im*(*pl).im);
      y+=(double)((*pk).re*(*pl).im-(*pk).im*(*pl).re);
      pk+=1;
      pl+=1;
   }

   z.re=(float)(x);
   z.im=(float)(y);

   return z;
}


int main(int argc,char *argv[])
{
   int my_rank,i,j,vol,off;
   int bs[4],Ns,nb,nv;
   int icom,ieo;
   float r,zsq;
   double d,dmax,dall;
   complex w,z;
   complex **wv,*pk,*pl;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);
      fin=freopen("check4.in","r",stdin);

      printf("\n");
      printf("Checks on the programs in the module valg\n");
      printf("-----------------------------------------\n\n");

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

   alloc_wv(10);
   wv=reserve_wv(10);
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
               random_v(vol,wv[i]+off,1.0f);

            dmax=0.0;

            for (i=0;i<10;i++)
            {
               pk=wv[i]+off;
               pl=wv[9-i]+off;

               if (icom==1)
               {
                  z=sp(vol,pk,pl);
                  MPI_Reduce(&z.re,&w.re,2,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
                  MPI_Bcast(&w.re,2,MPI_FLOAT,0,MPI_COMM_WORLD);
               }
               else
                  w=sp(vol,pk,pl);

               z=vprod(vol,icom,pk,pl);
               r=vnorm_square(vol,icom,pk)*vnorm_square(vol,icom,pl);
               d=(double)((z.re-w.re)*(z.re-w.re)+(z.im-w.im)*(z.im-w.im));
               d=sqrt(d/(double)(r));
               if (d>dmax)
                  dmax=d;

               z=vprod(vol,icom,pk,pk);
               r=vnorm_square(vol,icom,pk);

               d=fabs((double)(z.im/r));
               if (d>dmax)
                  dmax=d;

               d=fabs((double)(z.re/r-1.0f));
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Check of vprod and vnorm_square: %.2e\n\n",dmax);
            }

            dmax=0.0;
            z.re= 0.345f;
            z.im=-0.876f;
            zsq=z.re*z.re+z.im*z.im;

            for (i=0;i<9;i++)
            {
               pk=wv[i]+off;
               pl=wv[i+1]+off;

               w=vprod(vol,icom,pk,pl);
               r=vnorm_square(vol,icom,pk)+zsq*vnorm_square(vol,icom,pl)
                  +2.0f*(z.re*w.re-z.im*w.im);
               mulc_vadd(vol,pk,pl,z);

               d=fabs((double)(r/vnorm_square(vol,icom,pk)-1.0f));
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Consistency of vprod, vnorm_square\n");
               printf("and mulc_vadd: %.2e\n\n",dmax);
            }

            for (i=0;i<10;i++)
               random_v(vol,wv[i]+off,1.0f);

            dmax=0.0;

            for (i=0;i<10;i++)
            {
               pk=wv[i]+off;

               if (i>0)
               {
                  pl=wv[i-1]+off;
                  vproject(vol,icom,pk,pl);
                  z=vprod(vol,icom,pk,pl);

                  d=(fabs((double)(z.re))+
                     fabs((double)(z.im)))/
                     sqrt((double)(vnorm_square(vol,icom,pk)));

                  if (d>dmax)
                     dmax=d;
               }

               vnormalize(vol,icom,pk);
               r=vnorm_square(vol,icom,pk);

               d=fabs((double)(r-1.0f));
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Consistency of vprod, vnorm_square,\n");
               printf("vnormalize and vproject: %.2e\n\n",dmax);
            }

            for (i=0;i<5;i++)
            {
               pk=wv[i]+off;
               pl=wv[i+5]+off;

               random_v(vol,wv[i]+off,1.0f);
               assign_v2v(vol,pk,pl);

               for (j=0;j<5;j++)
               {
                  v[5*i+j].re=0.1234f*(float)(i^2)-0.8976f*(float)(j);
                  v[5*i+j].im=0.2231f*(float)(i)+0.9922f*(float)(j^2);
               }

               ppk[i]=pl;
            }

            vrotate(vol,5,ppk,v);
            dmax=0.0;

            for (i=5;i<10;i++)
            {
               pk=wv[i]+off;

               for (j=0;j<5;j++)
               {
                  z.re=-v[5*j+(i-5)].re;
                  z.im=-v[5*j+(i-5)].im;

                  pl=wv[j]+off;
                  mulc_vadd(vol,pk,pl,z);
               }

               r=vnorm_square(vol,icom,pk);

               d=fabs((double)(r));
               if (d>dmax)
                  dmax=d;
            }

            dmax/=(double)(vnorm_square(vol,icom,wv[0]+off));
            dmax=sqrt(dmax);

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Consistency of mulc_vadd\n");
               printf("and vrotate: %.2e\n\n",dmax);
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
