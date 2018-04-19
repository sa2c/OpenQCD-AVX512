
/*******************************************************************************
*
* File time2.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of ranlxd and gauss_dble
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "random.h"

#define NRLX 100
#define NGSS 24
#define NLOOPS 100000


int main(void)
{
   int k,level;
   float t1,t2,dt;
   double r[NRLX];

   printf("\n");
   printf("Timing of ranlxd ");
   printf("(average time per random number in microsec)\n\n");

   for (level=1;level<=2;level++)
   {
      rlxd_init(level,1);

      t1=(float)clock();
      for (k=1;k<=NLOOPS;k++) 
         ranlxd(r,NRLX);
      t2=(float)clock();
      
      dt=(t2-t1)/(float)(CLOCKS_PER_SEC);
      dt*=1.0e6f/(float)(NRLX*NLOOPS);      
      
      printf("%4.3f (level %1d)  ",dt,level);
   }

   printf("\n\n");
   printf("Timing of gauss_dble ");
   printf("(average time per random number in microsec)\n\n");

   for (level=1;level<=2;level++)
   {
      rlxd_init(level,1);

      t1=(float)clock();
      for (k=1;k<=NLOOPS;k++) 
         gauss_dble(r,NGSS);
      t2=(float)clock();
      
      dt=(t2-t1)/(float)(CLOCKS_PER_SEC);
      dt*=1.0e6f/(float)(NGSS*NLOOPS);      
      
      printf("%4.3f (level %1d)  ",dt,level);
   }
   
   printf("\n\n");
   exit(0);
}


