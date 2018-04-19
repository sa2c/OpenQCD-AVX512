
/*******************************************************************************
*
* File time3.c
*
* Copyright (C) 2005, 2009, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of inv_pauli_dble() and det_pauli_dble()
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "random.h"
#include "utils.h"
#include "random.h"
#include "linalg.h"
#include "sw_term.h"


int main(void)
{
   int n,count,itest;
   double t1,t2,dt,mu;
   pauli_dble *m;
   
   printf("\n");
   printf("Timing of inv_pauli_dble() and det_pauli_dble()\n");
   printf("-----------------------------------------------\n\n");

#if (defined x64)
   printf("Using SSE3 instructions and up to 16 xmm registers\n\n");
#endif
   
   printf("Measurement made with all data in cache\n\n");
   
   m=amalloc(2*sizeof(*m),4);
   error(m==NULL,1,"main [time3.c]",
         "Unable to allocate auxiliary arrays");
   
   rlxd_init(1,23456);
   ranlxd((*m).u,36);
   mu=0.1234;

   for (n=0;n<6;n++)
      (*m).u[n]=1.0;

   for (n=6;n<36;n++)
      (*m).u[n]=0.01*((*m).u[n]-0.5);

   n=(int)(1.0e5);
   dt=0.0;
   itest=0;

   while (dt<2.0)
   {   
      t1=(double)clock();
      for (count=0;count<n;count++)
         itest=inv_pauli_dble(mu,m,m+1);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }   

   dt*=2.0e6/(double)(n);
   error(itest!=0,1,"main [time3.c]","Inversion was not safe");

   printf("Time per call of inv_pauli_dble():\n");
   printf("%.3f micro sec\n\n",dt);   

   n=(int)(1.0e5);
   dt=0.0;
   itest=0;

   while (dt<2.0)
   {   
      t1=(double)clock();
      for (count=0;count<n;count++)
         det_pauli_dble(mu,m);
      t2=(double)clock();
      dt=(t2-t1)/(double)(CLOCKS_PER_SEC);
      n*=2;
   }   

   dt*=2.0e6/(double)(n);

   printf("Time per call of det_pauli_dble():\n");
   printf("%.3f micro sec\n\n",dt);   

   exit(0);
}
