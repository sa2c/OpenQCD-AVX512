
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2005, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Checks on the programs in the module salg.c
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "utils.h"
#include "lattice.h"
#include "sflds.h"
#include "linalg.h"
#include "global.h"

#define _acc_sp(z,x,y) \
   (z).re+=(double)((x).re*(y).re+(x).im*(y).im); \
   (z).im+=(double)((x).re*(y).im-(x).im*(y).re)

static complex v[25];
static spinor *ppk[5];


static complex sp(int vol,spinor *pk,spinor *pl)
{
   complex w;
   complex_dble z;
   spinor *pm;

   z.re=0.0;
   z.im=0.0;
   pm=pk+vol;

   for (;pk<pm;pk++)
   {
      _acc_sp(z,(*pk).c1.c1,(*pl).c1.c1);
      _acc_sp(z,(*pk).c1.c2,(*pl).c1.c2);
      _acc_sp(z,(*pk).c1.c3,(*pl).c1.c3);

      _acc_sp(z,(*pk).c2.c1,(*pl).c2.c1);
      _acc_sp(z,(*pk).c2.c2,(*pl).c2.c2);
      _acc_sp(z,(*pk).c2.c3,(*pl).c2.c3);

      _acc_sp(z,(*pk).c3.c1,(*pl).c3.c1);
      _acc_sp(z,(*pk).c3.c2,(*pl).c3.c2);
      _acc_sp(z,(*pk).c3.c3,(*pl).c3.c3);

      _acc_sp(z,(*pk).c4.c1,(*pl).c4.c1);
      _acc_sp(z,(*pk).c4.c2,(*pl).c4.c2);
      _acc_sp(z,(*pk).c4.c3,(*pl).c4.c3);

      pl+=1;
   }

   w.re=(float)(z.re);
   w.im=(float)(z.im);

   return w;
}


int main(int argc,char *argv[])
{
   int my_rank,i,j,vol,off;
   int icom,ieo;
   float r,zsq;
   double d,dmax,dall;
   complex z,w;
   spinor **ps,*pk,*pl,*pj;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);

      printf("\n");
      printf("Consistency of the programs in the module salg\n");
      printf("----------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   start_ranlux(0,12345);
   geometry();
   alloc_ws(10);
   ps=reserve_ws(10);
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

            vol=VOLUME/2;
            off=0;

            if (ieo==0)
               vol=VOLUME;
            if (ieo==2)
               off=VOLUME/2;

            for (i=0;i<10;i++)
               random_s(vol,ps[i]+off,1.0f);

            dmax=0.0;

            for (i=0;i<10;i++)
            {
               pk=ps[i]+off;
               pl=ps[9-i]+off;

               if (icom==1)
               {
                  z=sp(vol,pk,pl);
                  MPI_Reduce(&z.re,&w.re,2,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
                  MPI_Bcast(&w.re,2,MPI_FLOAT,0,MPI_COMM_WORLD);
               }
               else
                  w=sp(vol,pk,pl);

               z=spinor_prod(vol,icom,pk,pl);
               r=norm_square(vol,icom,pk)*norm_square(vol,icom,pl);
               d=(double)((z.re-w.re)*(z.re-w.re)+(z.im-w.im)*(z.im-w.im));
               d=sqrt(d/(double)(r));
               if (d>dmax)
                  dmax=d;

               r=spinor_prod_re(vol,icom,pk,pl);
               d=fabs((double)(z.re/r-1.0f));
               if (d>dmax)
                  dmax=d;

               z=spinor_prod(vol,icom,pk,pk);
               r=norm_square(vol,icom,pk);

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
               printf("Check of spinor_prod, spinor_prod_re\n");
               printf("and norm_square: %.2e\n\n",dmax);
            }

            dmax=0.0;
            z.re= 0.345f;
            z.im=-0.876f;
            zsq=z.re*z.re+z.im*z.im;

            for (i=0;i<9;i++)
            {
               pk=ps[i]+off;
               pl=ps[i+1]+off;

               w=spinor_prod(vol,icom,pk,pl);
               r=norm_square(vol,icom,pk)+zsq*norm_square(vol,icom,pl)
                  +2.0f*(z.re*w.re-z.im*w.im);
               mulc_spinor_add(vol,pk,pl,z);

               d=fabs((double)(r/norm_square(vol,icom,pk)-1.0f));
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Consistency of spinor_prod, norm_square\n");
               printf("and mulc_spinor_add: %.2e\n\n",dmax);
            }

            for (i=0;i<10;i++)
               random_s(vol,ps[i]+off,1.0f);

            dmax=0.0;
            r=-1.234f;
            z.re=-r;
            z.im=0.0f;

            for (i=0;i<8;i+=3)
            {
               pk=ps[i]+off;
               pl=ps[i+1]+off;
               pj=ps[i+2]+off;

               assign_s2s(vol,pk,pj);
               mulr_spinor_add(vol,pk,pl,r);
               mulc_spinor_add(vol,pk,pl,z);
               mulr_spinor_add(vol,pk,pj,-1.0);

               d=(double)(norm_square(vol,icom,pk)/norm_square(vol,icom,pj));
               d=sqrt(d);
               if (d>dmax)
                  dmax=d;

               assign_s2s(vol,pl,pk);
               scale(vol,r,pk);
               mulc_spinor_add(vol,pk,pl,z);

               d=(double)(norm_square(vol,icom,pk)/norm_square(vol,icom,pl));
               d=sqrt(d);
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Consistency of mulr_spinor_add, scale\n");
               printf("and mulc_spinor_add: %.2e\n\n",dmax);
            }

            for (i=0;i<10;i++)
               random_s(vol,ps[i]+off,1.0f);

            dmax=0.0;

            for (i=0;i<10;i++)
            {
               pk=ps[i]+off;

               if (i>0)
               {
                  pl=ps[i-1]+off;
                  project(vol,icom,pk,pl);
                  z=spinor_prod(vol,icom,pk,pl);

                  d=(fabs((double)(z.re))+
                     fabs((double)(z.im)))/
                     sqrt((double)(norm_square(vol,icom,pk)));

                  if (d>dmax)
                     dmax=d;
               }

               normalize(vol,icom,pk);
               r=norm_square(vol,icom,pk);

               d=fabs((double)(r-1.0f));
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Consistency of spinor_prod, norm_square,\n");
               printf("normalize and project: %.2e\n\n",dmax);
            }

            for (i=0;i<5;i++)
            {
               pk=ps[i]+off;
               pl=ps[i+5]+off;

               random_s(vol,ps[i]+off,1.0f);
               assign_s2s(vol,pk,pl);

               for (j=0;j<5;j++)
               {
                  v[5*i+j].re=0.1234f*(float)(i^2)-0.8976f*(float)(j);
                  v[5*i+j].im=0.2231f*(float)(i)+0.9922f*(float)(j^2);
               }

               ppk[i]=pl;
            }

            rotate(vol,5,ppk,v);
            dmax=0.0;

            for (i=5;i<10;i++)
            {
               pk=ps[i]+off;

               for (j=0;j<5;j++)
               {
                  z.re=-v[5*j+(i-5)].re;
                  z.im=-v[5*j+(i-5)].im;

                  pl=ps[j]+off;
                  mulc_spinor_add(vol,pk,pl,z);
               }

               r=norm_square(vol,icom,pk);

               d=fabs((double)(r));
               if (d>dmax)
                  dmax=d;
            }

            dmax/=(double)(norm_square(vol,icom,ps[0]+off));
            dmax=sqrt(dmax);

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Consistency of mulc_spinor_add\n");
               printf("and rotate: %.2e\n\n",dmax);
            }

            dmax=0.0;

            for (i=0;i<5;i++)
            {
               pk=ps[i]+off;
               pl=ps[9-i]+off;
               random_s(vol,pk,1.0f);
               assign_s2s(vol,pk,pl);
               mulg5(vol,pk);
               mulg5(vol,pk);

               z.re=-1.0f;
               z.im=0.0f;

               mulc_spinor_add(vol,pl,pk,z);
               r=norm_square(vol,icom,pl)/norm_square(vol,icom,pk);
               d=sqrt((double)(r));
               if (d>dmax)
                  dmax=d;

               random_s(vol,pl,1.0f);
               z=spinor_prod(vol,icom,pk,pl);
               mulg5(vol,pk);
               mulg5(vol,pl);
               w=spinor_prod(vol,icom,pk,pl);

               d=(fabs((double)(z.re-w.re))+fabs((double)(z.im-w.im)))/
                  (fabs((double)(z.re))+fabs((double)(z.im)));
               if (d>dmax)
                  dmax=d;

               random_s(vol,pk,1.0f);
               assign_s2s(vol,pk,pl);
               mulg5(vol,pk);
               mulmg5(vol,pk);

               z.re=1.0f;
               z.im=0.0f;

               mulc_spinor_add(vol,pl,pk,z);
               r=norm_square(vol,icom,pl)/norm_square(vol,icom,pk);
               d=sqrt((double)(r));
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Check of mulg5 and mulmg5: %.2e\n\n",dmax);
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
