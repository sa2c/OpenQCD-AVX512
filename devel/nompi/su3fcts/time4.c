
/*******************************************************************************
*
* File time4.c
*
* Copyright (C) 2005, 2009, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of prod2su3alg, prod2u3alg and rotate_su3alg
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
   su3_dble *u,*v;
   su3_alg_dble *X;
   u3_alg_dble *Y;

   printf("\n");
   printf("Timing of prod2su3alg, prod2u3alg and rotate_su3alg\n");
   printf("---------------------------------------------------\n\n");

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

   u=amalloc(2*sizeof(*u),4);
   X=amalloc(sizeof(*X),4);
   Y=amalloc(sizeof(*Y),4);
   error((u==NULL)||(X==NULL)||(Y==NULL),1,
         "main [time4.c]","Unable to allocate auxiliary variables");
   v=u+1;

   rlxd_init(1,23456);
   random_su3_dble(u);
   random_su3_dble(v);
   ranlxd((double*)(X),8);
   ranlxd((double*)(Y),9);

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         prod2su3alg(u,v,X);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6/(double)(n);

   printf("The times per application are:\n");
   printf("prod2su3alg:   %4.3f nsec [%d Mflops]\n",1.0e3*dt,(int)(212.0/dt));

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         prod2u3alg(u,v,Y);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6/(double)(n);

   printf("prod2u3alg:    %4.3f nsec [%d Mflops]\n",1.0e3*dt,(int)(207.0/dt));

   n=(int)(1.0e6);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         rotate_su3alg(u,X);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6/(double)(n);

   printf("rotate_su3alg: %4.3f nsec [%d Mflops]\n\n",1.0e3*dt,(int)(274.0/dt));
   exit(0);
}
