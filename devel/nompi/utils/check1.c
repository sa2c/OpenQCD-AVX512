
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2007, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Test of the endianness and byte swapping programs
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "random.h"
#include "utils.h"

#define N 18000

static int istd[N],istds[N];
static double dstd[N],dstds[N];
static su3_dble ufld[N/18];


static void set_u2v(su3_dble *u,double *v)
{
   v[ 0]=(*u).c11.re;
   v[ 1]=(*u).c11.im;
   v[ 2]=(*u).c12.re;
   v[ 3]=(*u).c12.im;
   v[ 4]=(*u).c13.re;
   v[ 5]=(*u).c13.im;

   v[ 6]=(*u).c21.re;
   v[ 7]=(*u).c21.im;
   v[ 8]=(*u).c22.re;
   v[ 9]=(*u).c22.im;
   v[10]=(*u).c23.re;
   v[11]=(*u).c23.im;

   v[12]=(*u).c31.re;
   v[13]=(*u).c31.im;
   v[14]=(*u).c32.re;
   v[15]=(*u).c32.im;
   v[16]=(*u).c33.re;
   v[17]=(*u).c33.im;
}


static void set_v2u(double *v,su3_dble *u)
{
   (*u).c11.re=v[ 0];
   (*u).c11.im=v[ 1];
   (*u).c12.re=v[ 2];
   (*u).c12.im=v[ 3];
   (*u).c13.re=v[ 4];
   (*u).c13.im=v[ 5];

   (*u).c21.re=v[ 6];
   (*u).c21.im=v[ 7];
   (*u).c22.re=v[ 8];
   (*u).c22.im=v[ 9];
   (*u).c23.re=v[10];
   (*u).c23.im=v[11];

   (*u).c31.re=v[12];
   (*u).c31.im=v[13];
   (*u).c32.re=v[14];
   (*u).c32.im=v[15];
   (*u).c33.re=v[16];
   (*u).c33.im=v[17];   
}


int main(void)
{
   int ie,k,it;
   stdint_t i[2];
   double d[2];
   char *ci[2],*cd[2];
   
   printf("\n");
   printf("Test of the endianness and byte swapping programs\n");
   printf("-------------------------------------------------\n\n");

   printf("sizeof(stdint_t) = %d\n",(int)(sizeof(stdint_t)));
   printf("sizeof(double) = %d\n",(int)(sizeof(double)));

   ie=endianness();
   if (ie==LITTLE_ENDIAN)
      printf("The machine is little endian\n\n");
   else if (ie==BIG_ENDIAN)
      printf("The machine is big endian\n\n");
   else
      printf("The machine has unknown endianness\n\n");

   ci[0]=(char*)(i);
   ci[1]=(char*)(i+1);

   ci[0][0]='A';
   ci[0][1]='B';
   ci[0][2]='C';
   ci[0][3]='D';   
   
   ci[1][0]='1';
   ci[1][1]='2';
   ci[1][2]='3';
   ci[1][3]='4';   

   printf("Byte swapping integers:\n");
   printf("%.4s, %.4s  ->  ",ci[0],ci[1]);
   bswap_int(2,i);
   printf("%.4s, %.4s\n\n",ci[0],ci[1]);

   cd[0]=(char*)(d);
   cd[1]=(char*)(d+1);

   cd[0][0]='A';
   cd[0][1]='B';
   cd[0][2]='C';
   cd[0][3]='D';
   cd[0][4]='E';
   cd[0][5]='F';
   cd[0][6]='G';
   cd[0][7]='H';    

   cd[1][0]='1';
   cd[1][1]='2';
   cd[1][2]='3';
   cd[1][3]='4';
   cd[1][4]='5';
   cd[1][5]='6';
   cd[1][6]='7';
   cd[1][7]='8';    

   printf("Byte swapping double precision numbers:\n");
   printf("%.8s, %.8s  ->  ",cd[0],cd[1]);
   bswap_double(2,d);
   printf("%.8s, %.8s\n\n",cd[0],cd[1]);

   gauss_dble(dstd,N);

   for (k=0;k<N;k++)
   {
      dstds[k]=dstd[k];
      istd[k]=(int)(dstd[k]);
      istds[k]=istd[k];
   }

   bswap_int(N,istd);
   bswap_int(N,istd);

   it=0;

   for (k=0;k<N;k++)
      if ((istds[k]-istd[k])!=0)
         it+=1;

   printf("2x(byte swap) %d integer random numbers:\n",N);
   printf("Number of failures = %d\n\n",it);
   
   bswap_double(N,dstd);
   bswap_double(N,dstd);

   it=0;

   for (k=0;k<N;k++)
      if ((dstds[k]-dstd[k])!=0.0)
         it+=1;

   printf("2x(byte swap) %d double-precision random numbers:\n",N);
   printf("Number of failures = %d\n\n",it);

   gauss_dble(dstd,N);

   for (k=0;k<(N/18);k++)
      set_v2u(dstd+18*k,ufld+k);

   bswap_double(N,ufld);

   for (k=0;k<(N/18);k++)
      set_u2v(ufld+k,dstds+18*k);   

   bswap_double(N,dstds);   
   it=0;

   for (k=0;k<N;k++)
      if ((dstds[k]-dstd[k])!=0.0)
         it+=1;

   printf("Byte swap %d random su3_dble matrices:\n",N/18);
   printf("Number of failures = %d\n\n",it);
   
   exit(0);
}
