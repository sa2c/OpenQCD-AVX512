
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2005, 2011, 2012, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Checks on the programs in the module salg_dble.c
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

static complex_dble v[25];
static spinor_dble *ppk[5];


static complex_dble sp(int vol,spinor_dble *pk,spinor_dble *pl)
{
   complex_dble z;
   spinor_dble *pm;

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

   return z;
}


int main(int argc,char *argv[])
{
   int my_rank,i,j,vol,off;
   int icom,ieo;
   double r,cs,cr,zsq,d,dmax,dall;
   complex_dble z,w;
   spinor_dble **psd,*pk,*pl,*pj;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);

      printf("\n");
      printf("Consistency of the programs in the module salg_dble\n");
      printf("---------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   start_ranlux(0,12345);
   geometry();
   alloc_wsd(10);
   psd=reserve_wsd(10);
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
               random_sd(vol,psd[i]+off,1.0);

            dmax=0.0;

            for (i=0;i<10;i++)
            {
               pk=psd[i]+off;
               pl=psd[9-i]+off;

               if (icom==1)
               {
                  z=sp(vol,pk,pl);
                  MPI_Reduce(&z.re,&w.re,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
                  MPI_Bcast(&w.re,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
               }
               else
                  w=sp(vol,pk,pl);

               z=spinor_prod_dble(vol,icom,pk,pl);
               r=norm_square_dble(vol,icom,pk)*norm_square_dble(vol,icom,pl);
               d=(z.re-w.re)*(z.re-w.re)+(z.im-w.im)*(z.im-w.im);
               d=sqrt(d/r);
               if (d>dmax)
                  dmax=d;

               r=spinor_prod_re_dble(vol,icom,pk,pl);

               d=fabs(z.re/r-1.0);
               if (d>dmax)
                  dmax=d;

               z=spinor_prod_dble(vol,icom,pk,pk);
               r=norm_square_dble(vol,icom,pk);

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
               printf("Check of spinor_prod, spinor_prod_re\n");
               printf("and norm_square: %.2e\n\n",dmax);
            }

            dmax=0.0;

            for (i=0;i<10;i++)
            {
               pk=psd[i]+off;
               pl=psd[9-i]+off;

               z=spinor_prod5_dble(vol,icom,pk,pl);
               mulg5_dble(vol,pl);
               w=spinor_prod_dble(vol,icom,pk,pl);

               r=norm_square_dble(vol,icom,pk)*norm_square_dble(vol,icom,pl);
               d=(z.re-w.re)*(z.re-w.re)+(z.im-w.im)*(z.im-w.im);
               d=sqrt(d/r);
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Consistency check of spinor_prod5, mulg5\n");
               printf("and spinor_prod: %.2e\n\n",dmax);
            }

            dmax=0.0;
            z.re= 0.345;
            z.im=-0.876;
            zsq=z.re*z.re+z.im*z.im;

            for (i=0;i<9;i++)
            {
               pk=psd[i]+off;
               pl=psd[i+1]+off;

               w=spinor_prod_dble(vol,icom,pk,pl);
               r=norm_square_dble(vol,icom,pk)+zsq*norm_square_dble(vol,icom,pl)
                  +2.0*(z.re*w.re-z.im*w.im);
               mulc_spinor_add_dble(vol,pk,pl,z);

               d=fabs(r/norm_square_dble(vol,icom,pk)-1.0);
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
               random_sd(vol,psd[i]+off,1.0);

            dmax=0.0;
            r=-1.234;
            z.re=-r;
            z.im=0.0;

            for (i=0;i<8;i+=3)
            {
               pk=psd[i]+off;
               pl=psd[i+1]+off;
               pj=psd[i+2]+off;

               assign_sd2sd(vol,pk,pj);
               mulr_spinor_add_dble(vol,pk,pl,r);
               mulc_spinor_add_dble(vol,pk,pl,z);
               mulr_spinor_add_dble(vol,pk,pj,-1.0);

               d=norm_square_dble(vol,icom,pk)/norm_square_dble(vol,icom,pj);
               d=sqrt(d);
               if (d>dmax)
                  dmax=d;

               assign_sd2sd(vol,pl,pk);
               scale_dble(vol,r,pk);
               mulc_spinor_add_dble(vol,pk,pl,z);

               d=norm_square_dble(vol,icom,pk)/norm_square_dble(vol,icom,pl);
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
               random_sd(vol,psd[i]+off,1.0);

            dmax=0.0;
            cs=0.785;
            cr=-1.567;

            for (i=0;i<8;i+=3)
            {
               pk=psd[i]+off;
               pl=psd[i+1]+off;
               pj=psd[i+2]+off;

               assign_sd2sd(vol,pk,pj);
               combine_spinor_dble(vol,pk,pl,cs,cr);
               scale_dble(vol,cs,pj);
               mulr_spinor_add_dble(vol,pj,pl,cr);
               mulr_spinor_add_dble(vol,pk,pj,-1.0);

               d=norm_square_dble(vol,icom,pk)/norm_square_dble(vol,icom,pj);
               d=sqrt(d);
               if (d>dmax)
                  dmax=d;
            }

            if (my_rank==0)
            {
               if (dmax>dall)
                  dall=dmax;
               printf("Consistency of mulr_spinor_add, scale\n");
               printf("and combine_spinor: %.2e\n\n",dmax);
            }

            for (i=0;i<10;i++)
               random_sd(vol,psd[i]+off,1.0);

            dmax=0.0;

            for (i=0;i<10;i++)
            {
               pk=psd[i]+off;

               if (i>0)
               {
                  pl=psd[i-1]+off;
                  project_dble(vol,icom,pk,pl);
                  z=spinor_prod_dble(vol,icom,pk,pl);

                  d=(fabs(z.re)+fabs(z.im))/sqrt(norm_square_dble(vol,icom,pk));

                  if (d>dmax)
                     dmax=d;
               }

               normalize_dble(vol,icom,pk);
               r=norm_square_dble(vol,icom,pk);

               d=fabs(r-1.0);
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
               pk=psd[i]+off;
               pl=psd[i+5]+off;

               random_sd(vol,psd[i]+off,1.0);
               assign_sd2sd(vol,pk,pl);

               for (j=0;j<5;j++)
               {
                  v[5*i+j].re=0.1234*(double)(i^2)-0.8976*(double)(j);
                  v[5*i+j].im=0.2231*(double)(i)+0.9922*(double)(j^2);
               }

               ppk[i]=pl;
            }

            rotate_dble(vol,5,ppk,v);
            dmax=0.0;

            for (i=5;i<10;i++)
            {
               pk=psd[i]+off;

               for (j=0;j<5;j++)
               {
                  z.re=-v[5*j+(i-5)].re;
                  z.im=-v[5*j+(i-5)].im;

                  pl=psd[j]+off;
                  mulc_spinor_add_dble(vol,pk,pl,z);
               }

               r=norm_square_dble(vol,icom,pk);

               d=fabs(r);
               if (d>dmax)
                  dmax=d;
            }

            dmax/=norm_square_dble(vol,icom,psd[0]+off);
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
               pk=psd[i]+off;
               pl=psd[9-i]+off;
               random_sd(vol,pk,1.0);
               assign_sd2sd(vol,pk,pl);
               mulg5_dble(vol,pk);
               mulg5_dble(vol,pk);

               z.re=-1.0;
               z.im=0.0;

               mulc_spinor_add_dble(vol,pl,pk,z);
               r=norm_square_dble(vol,icom,pl)/norm_square_dble(vol,icom,pk);
               d=sqrt(r);
               if (d>dmax)
                  dmax=d;

               random_sd(vol,pl,1.0);
               z=spinor_prod_dble(vol,icom,pk,pl);
               mulg5_dble(vol,pk);
               mulg5_dble(vol,pl);
               w=spinor_prod_dble(vol,icom,pk,pl);

               d=(fabs(z.re-w.re)+fabs(z.im-w.im))/
                  (fabs(z.re)+fabs(z.im));
               if (d>dmax)
                  dmax=d;

               random_sd(vol,pk,1.0);
               assign_sd2sd(vol,pk,pl);
               mulg5_dble(vol,pk);
               mulmg5_dble(vol,pk);

               z.re=1.0;
               z.im=0.0;

               mulc_spinor_add_dble(vol,pl,pk,z);
               r=norm_square_dble(vol,icom,pl)/norm_square_dble(vol,icom,pk);
               d=sqrt(r);
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
