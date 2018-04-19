
/*******************************************************************************
*
* File time2.c
*
* Copyright (C) 2005, 2008, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of the salg_dble routines
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

static complex_dble *vmat,*wmat;
static spinor_dble **psd,*ppk[5];


static double wt_spinor_prod_dble(int nflds,int icom)
{
   int my_rank,nmax,n,i,ib;
   double wt1,wt2,wdt,wtav;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   for (i=0;i<nflds;i++)
      random_sd(VOLUME,psd[i],1.0);

   nmax=1;

   for (ib=0;ib<1;nmax*=2)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      for (n=0;n<nmax;n++)
      {
         for (i=0;i<nflds;i+=2)
            spinor_prod_dble(VOLUME,icom,psd[i],psd[i+1]);
      }

      wt2=MPI_Wtime();
      wdt=wt2-wt1;

      MPI_Reduce(&wdt,&wtav,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      if (my_rank==0)
      {
         wtav/=(double)(NPROC);

         if (wtav>2.0)
            ib=1;

         wtav/=(double)((nmax*nflds)/2);
      }

      MPI_Bcast(&ib,1,MPI_INT,0,MPI_COMM_WORLD);
   }

   MPI_Bcast(&wtav,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   return wtav;
}


static double wt_norm_square_dble(int nflds,int icom)
{
   int my_rank,nmax,n,i,ib;
   double wt1,wt2,wdt,wtav;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   nmax=1;

   for (ib=0;ib<1;nmax*=2)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      for (n=0;n<nmax;n++)
      {
         for (i=0;i<nflds;i++)
            (void)(norm_square_dble(VOLUME,icom,psd[i]));
      }

      wt2=MPI_Wtime();
      wdt=wt2-wt1;

      MPI_Reduce(&wdt,&wtav,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      if (my_rank==0)
      {
         wtav/=(double)(NPROC);

         if (wtav>2.0)
            ib=1;

         wtav/=(double)(nmax*nflds);
      }

      MPI_Bcast(&ib,1,MPI_INT,0,MPI_COMM_WORLD);
   }

   MPI_Bcast(&wtav,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   return wtav;
}


static double wt_normalize_dble(int nflds,int icom)
{
   int my_rank,nmax,n,i,ib;
   double wt1,wt2,wdt,wtav;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   nmax=1;

   for (ib=0;ib<1;nmax*=2)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      for (n=0;n<nmax;n++)
      {
         for (i=0;i<nflds;i++)
            (void)(normalize_dble(VOLUME,icom,psd[i]));
      }

      wt2=MPI_Wtime();
      wdt=wt2-wt1;

      MPI_Reduce(&wdt,&wtav,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      if (my_rank==0)
      {
         wtav/=(double)(NPROC);

         if (wtav>2.0)
            ib=1;

         wtav/=(double)(nmax*nflds);
      }

      MPI_Bcast(&ib,1,MPI_INT,0,MPI_COMM_WORLD);
   }

   MPI_Bcast(&wtav,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   return wtav;
}


static double wt_mulc_spinor_add_dble(int nflds)
{
   int my_rank,nmax,n,i,ib;
   complex_dble z;
   double wt1,wt2,wdt,wtav;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   z.re=0.123;
   z.im=0.456;
   nmax=1;

   for (ib=0;ib<1;nmax*=2)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      for (n=0;n<nmax;n++)
      {
         for (i=0;i<nflds;i+=2)
            mulc_spinor_add_dble(VOLUME,psd[i],psd[i+1],z);

         z.re-=z.re;
         z.im-=z.im;
      }

      wt2=MPI_Wtime();
      wdt=wt2-wt1;

      MPI_Reduce(&wdt,&wtav,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      if (my_rank==0)
      {
         wtav/=(double)(NPROC);

         if (wtav>2.0)
            ib=1;

         wtav/=(double)((nmax*nflds)/2);
      }

      MPI_Bcast(&ib,1,MPI_INT,0,MPI_COMM_WORLD);
   }

   MPI_Bcast(&wtav,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   return wtav;
}


static double wt_project_dble(int nflds,int icom)
{
   int my_rank,nmax,n,i,ib;
   double wt1,wt2,wdt,wtav;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   for (i=0;i<nflds;i++)
   {
      random_sd(VOLUME,psd[i],1.0);
      (void)(normalize_dble(VOLUME,icom,psd[i]));
   }

   nmax=1;

   for (ib=0;ib<1;nmax*=2)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      for (n=0;n<nmax;n++)
      {
         for (i=0;i<nflds;i+=2)
            project_dble(VOLUME,icom,psd[i],psd[i+1]);
      }

      wt2=MPI_Wtime();
      wdt=wt2-wt1;

      MPI_Reduce(&wdt,&wtav,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      if (my_rank==0)
      {
         wtav/=(double)(NPROC);

         if (wtav>2.0)
            ib=1;

         wtav/=(double)((nmax*nflds)/2);
      }

      MPI_Bcast(&ib,1,MPI_INT,0,MPI_COMM_WORLD);
   }

   MPI_Bcast(&wtav,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   return wtav;
}


static void gram_schmidt(int n,spinor_dble **s)
{
   int i,j,k;

   for (i=0;i<n;i++)
   {
      for (k=0;k<2;k++)
      {
         for (j=0;j<i;j++)
            project_dble(VOLUME,0,s[i],s[j]);
      }

      (void)(normalize_dble(VOLUME,0,s[i]));
   }
}


static void random_unitary(int n)
{
   int i,j,k;

   vmat=amalloc(2*n*n*sizeof(*vmat),4);
   wmat=vmat+n*n;
   error(vmat==NULL,1,"random_unitary [time2.c]",
         "Unable to allocate auxiliary arrays");

   gauss_dble((double*)(vmat),2*n*n);

   for (i=0;i<n;i++)
   {
      for (k=0;k<2;k++)
      {
         for (j=0;j<i;j++)
            vproject_dble(n,0,vmat+i*n,vmat+j*n);
      }

      (void)(vnormalize_dble(n,0,vmat+i*n));
   }

   for (i=0;i<n;i++)
   {
      for (j=0;j<n;j++)
      {
         wmat[i*n+j].re=vmat[j*n+i].re;
         wmat[i*n+j].im=-vmat[j*n+i].im;
      }
   }
}


static double wt_rotate_dble(void)
{
   int my_rank,nmax,n,i,ib;
   double wt1,wt2,wdt,wtav;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   for (i=0;i<5;i++)
   {
      random_sd(VOLUME,psd[i],1.0);
      ppk[i]=psd[i];
   }

   gram_schmidt(5,ppk);
   for (i=0;i<5;i++)
      assign_sd2sd(VOLUME,psd[i],psd[i+5]);
   random_unitary(5);
   nmax=1;

   for (ib=0;ib<1;nmax*=2)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      for (n=0;n<nmax;n++)
      {
         if ((n&0xf)==0x0)
         {
            for (i=0;i<5;i++)
               assign_sd2sd(VOLUME,psd[i+5],psd[i]);
         }

         rotate_dble(VOLUME,5,ppk,vmat);
         rotate_dble(VOLUME,5,ppk,wmat);
      }

      wt2=MPI_Wtime();
      wdt=wt2-wt1;

      MPI_Reduce(&wdt,&wtav,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      if (my_rank==0)
      {
         wtav/=(double)(NPROC);

         if (wtav>2.0)
            ib=1;

         wtav/=(double)(2*nmax);
      }

      MPI_Bcast(&ib,1,MPI_INT,0,MPI_COMM_WORLD);
   }

   MPI_Bcast(&wtav,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   return wtav;
}


int main(int argc,char *argv[])
{
   int my_rank,icom,nflds;
   double wdt;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("time2.log","w",stdout);

      printf("\n");
      printf("Timing of the salg_dble routines\n");
      printf("--------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      if (NPROC>1)
         printf("There are %d MPI processes\n",NPROC);
      else
         printf("There is 1 MPI process\n");

      if ((VOLUME*sizeof(double))<(64*1024))
         printf("The local size of a quark field is %d KB\n",
                (int)((24*VOLUME*sizeof(double))/(1024)));
      else
         printf("The local size of a quark field is %d MB\n",
                (int)((24*VOLUME*sizeof(double))/(1024*1024)));

#if (defined x64)
#if (defined AVX)
      printf("Using AVX instructions\n");
#else
      printf("Using SSE3 instructions and 16 xmm registers\n");
#endif
#if (defined P3)
      printf("Assuming SSE prefetch instructions fetch 32 bytes\n");
#elif (defined PM)
      printf("Assuming SSE prefetch instructions fetch 64 bytes\n");
#elif (defined P4)
      printf("Assuming SSE prefetch instructions fetch 128 bytes\n");
#else
      printf("SSE prefetch instructions are not used\n");
#endif
#endif
      printf("\n");
   }

   icom=1;
   start_ranlux(0,12345);
   geometry();

   nflds=(int)((4*1024*1024)/(VOLUME*sizeof(double)))+1;
   if ((nflds%2)==1)
      nflds+=1;
   if (nflds<10)
      nflds=10;
   alloc_wsd(nflds);
   psd=reserve_wsd(nflds);

   wdt=1.0e6*wt_spinor_prod_dble(nflds,icom)/(double)(VOLUME);

   if (my_rank==0)
   {
      printf("Function spinor_prod_dble:\n");
      printf("Time per lattice point: %4.3f micro sec\n",wdt);
      printf("%d Mflops [%d bit arithmetic]\n\n",
             (int)(96.0/wdt),(int)(sizeof(spinor_dble))/3);
   }

   wdt=1.0e6*wt_norm_square_dble(nflds,icom)/(double)(VOLUME);

   if (my_rank==0)
   {
      printf("Function norm_square_dble:\n");
      printf("Time per lattice point: %4.3f micro sec\n",wdt);
      printf("%d Mflops [%d bit arithmetic]\n\n",
             (int)(48.0/wdt),(int)(sizeof(spinor_dble))/3);
   }

   wdt=1.0e6*wt_normalize_dble(nflds,icom)/(double)(VOLUME);

   if (my_rank==0)
   {
      printf("Function normalize_dble:\n");
      printf("Time per lattice point: %4.3f micro sec\n",wdt);
      printf("%d Mflops [%d bit arithmetic]\n\n",
             (int)(72.0/wdt),(int)(sizeof(spinor_dble))/3);
   }

   wdt=1.0e6*wt_mulc_spinor_add_dble(nflds)/(double)(VOLUME);

   if (my_rank==0)
   {
      printf("Function mulc_spinor_add_dble:\n");
      printf("Time per lattice point: %4.3f micro sec\n",wdt);
      printf("%d Mflops [%d bit arithmetic]\n\n",
             (int)(96.0/wdt),(int)(sizeof(spinor_dble))/3);
   }

   wdt=1.0e6*wt_project_dble(nflds,icom)/(double)(VOLUME);

   if (my_rank==0)
   {
      printf("Function project_dble:\n");
      printf("Time per lattice point: %4.3f micro sec\n",wdt);
      printf("%d Mflops [%d bit arithmetic]\n\n",
             (int)(192.0/wdt),(int)(sizeof(spinor_dble))/3);
   }

   wdt=1.0e6*wt_rotate_dble()/(double)(25*VOLUME);

   if (my_rank==0)
   {
      printf("Function rotate_dble (n=5 fields):\n");
      printf("Time per lattice point: %4.3f*n^2 micro sec\n",wdt);
      printf("%d Mflops [%d bit arithmetic]\n\n",
             (int)(91.2/wdt),(int)(sizeof(spinor_dble))/3);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
