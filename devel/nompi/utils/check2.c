
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Test of the program fdigits()
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"


int main(void)
{
   int n;
   double x;
   
   printf("\n");
   printf("Test of the program fdigits()\n");
   printf("-----------------------------\n\n");

   while (1)
   {
      printf("x = ");
      scanf("%lf",&x);
      n=fdigits(x);
      printf("    %.*f\n\n",n,x);
   }
   
   exit(0);
}
