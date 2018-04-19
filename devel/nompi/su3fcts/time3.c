
/*******************************************************************************
*
* File time3.c
*
* Copyright (C) 2005, 2009, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of su3xsu3, su3dagxsu3, ...
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "random.h"
#include "utils.h"
#include "su3fcts.h"


int main(void)
{
   int n,count;
   double t1,t2,dt;
   su3_dble *u,*v,*w;
   u3_alg_dble *X;

   printf("\n");
   printf("Timing of su3xsu3, su3dagxsu3, ...\n");
   printf("----------------------------------\n\n");

#if (defined AVX)
#if (defined FMA3)
   printf("Using AVX and FMA3 instructions\n");
#else
   printf("Using AVX instructions\n");
#endif
#elif (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n");
#endif

   printf("Measurement made with all data in cache\n\n");

   u=amalloc(3*sizeof(su3_dble),4);
   X=amalloc(sizeof(u3_alg_dble),3);
   error((u==NULL)||(X==NULL),1,"main [time3.c]",
         "Unable to allocate auxiliary array");
   v=u+1;
   w=u+2;

   rlxd_init(1,23456);
   random_su3_dble(u);
   random_su3_dble(v);
   ranlxd((double*)(&(*X).c1),9);

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         su3xsu3(u,v,w);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6/(double)(n);

   printf("The times per application are:\n");
   printf("su3xsu3:       %4.3f nsec [%d Mflops]\n",1.0e3*dt,(int)(198.0/dt));

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         su3dagxsu3(u,v,w);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6/(double)(n);

   printf("su3dagxsu3:    %4.3f nsec [%d Mflops]\n",1.0e3*dt,(int)(198.0/dt));

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         su3xsu3dag(u,v,w);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6/(double)(n);

   printf("su3xsu3dag:    %4.3f nsec [%d Mflops]\n",1.0e3*dt,(int)(198.0/dt));

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         su3dagxsu3dag(u,v,w);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6/(double)(n);

   printf("su3dagxsu3dag: %4.3f nsec [%d Mflops]\n",1.0e3*dt,(int)(198.0/dt));

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         su3xu3alg(u,X,v);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6/(double)(n);

   printf("su3xu3alg:     %4.3f nsec [%d Mflops]\n",1.0e3*dt,(int)(198.0/dt));

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         su3dagxu3alg(u,X,v);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6/(double)(n);

   printf("su3dagxu3alg:  %4.3f nsec [%d Mflops]\n",1.0e3*dt,(int)(198.0/dt));

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         u3algxsu3(X,u,v);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6/(double)(n);

   printf("u3algxsu3:     %4.3f nsec [%d Mflops]\n",1.0e3*dt,(int)(198.0/dt));

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         u3algxsu3dag(X,u,v);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6/(double)(n);

   printf("u3algxsu3dag:  %4.3f nsec [%d Mflops]\n\n",1.0e3*dt,(int)(198.0/dt));

   exit(0);
}
