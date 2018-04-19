
/*******************************************************************************
*
* File time2.c
*
* Copyright (C) 2007, 2009, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of cmat_vec_dble, cmat_mul_dble and cmat_inv_dble
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "random.h"
#include "su3.h"
#include "utils.h"
#include "linalg.h"


int main(void)
{
   int ie,ir,nm,n,count;
   double k,t1,t2,dt;
   complex_dble *a,*b,*c,*v,*w;

   printf("\n");
   printf("Timing of cmat_vec_dble, cmat_mul_dble and cmat_inv_dble\n");
   printf("--------------------------------------------------------\n\n");

#if (defined AVX)
#if (defined FMA3)
   printf("Using AVX and FMA3 instructions\n\n");
#else
   printf("Using AVX instructions\n\n");
#endif
#elif (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n\n");
#endif

   printf("Measurement made with all data in cache\n\n");

   printf("Matrix size: ");
   ir=scanf(" %d",&nm);

   error((ir!=1)||(nm<1),1,"main [time2.c]",
         "Read error or improper matrix size");

   a=amalloc((3*nm*nm+2*nm)*sizeof(*a),6);
   error(a==NULL,1,"main [time2.c]","Unable to allocate auxiliary arrays");

   rlxd_init(1,23456);
   ranlxd((double*)(a),6*nm*nm+4*nm);

   b=a+nm*nm;
   c=b+nm*nm;
   v=c+nm*nm;
   w=v+nm;

   n=(int)(1.0e7)/(nm*nm);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         cmat_vec_dble(nm,a,v,w);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6f/(double)(n);

   printf("\n");
   printf("Time per call of cmat_vec_dble():\n");
   printf("%.2e micro sec (%d Mflops)\n\n",
          dt,(int)((double)(nm*(6+(nm-1)*8))/dt));

   n=(int)(1.0e7)/(nm*nm*nm);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         cmat_mul_dble(nm,a,b,c);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6f/(double)(n);

   printf("Time per call of cmat_mul_dble():\n");
   printf("%.2e micro sec (%d Mflops)\n\n",
          dt,(int)((double)(nm*nm*(6+(nm-1)*8))/dt));

   for (n=0;n<nm;n++)
      a[n*nm+n].re=(double)(10*nm*nm);

   ie=cmat_inv_dble(nm,a,b,&k);
   error(ie!=0,1,"main [time2.c]","Matrix is not safely invertible");

   n=(int)(1.0e7)/(nm*nm*nm);
   dt=0.0;

   while (dt<2.0)
   {
      t1=(double)clock();
      for (count=0;count<n;count++)
         ie=cmat_inv_dble(nm,a,b,&k);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }

   dt*=2.0e6f/(double)(n);

   printf("Time per call of cmat_inv_dble():\n");
   printf("%.2e micro sec (~%d Mflops)\n\n",
          dt,(int)((double)(12*nm*nm*nm)/dt));

   exit(0);
}
