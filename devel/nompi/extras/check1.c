
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2010 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the modified Bessel function I0(x) [program i0m()]
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "extras.h"


int main(void)
{
   double x,y;
   
   printf("\n");
   printf("Modified Bessel function I0(x) [program i0m()]\n");
   printf("----------------------------------------------\n\n");

   printf("Print selected values:\n\n");
   
   for (;;)
   {
      printf("Specify x: ");

      if (scanf("%lf",&x)==1)
      {
         y=i0m(x);
         printf("x = %.4e, exp(-x)*I0(x) = %.15e\n\n",x,y);
      }
      else
      {
         printf("No value specified, program stopped\n\n");
         break;
      }
   }
   
   exit(0);
}
