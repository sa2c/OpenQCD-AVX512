
/*******************************************************************************
*
* File table1.c
*
* Copyright (C) 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Table of the relative error of the Zolotarev rational approximation to the
* function f(x)=1/|x| (suitable for plotting, for example)
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "ratfcts.h"

static int n;
static double eps,delta,As,*ars;


static void alloc_ars(void)
{
   ars=amalloc(2*n*sizeof(double),3);

   error(ars==NULL,1,"alloc_ars [table1.c]",
         "Unable to allocate coefficient array");
}


static double zolotarev_sign(double x)
{
   int r;
   double y,p;

   y=x*x;
   p=1.0;

   for (r=0;r<n;r++)
      p*=((y+ars[2*r])/(y+ars[2*r+1]));

   return x*As*p;
}


int main(void)
{
   int k,kmax,ir;
   double a,b,dx,x;
   FILE *dat=NULL;

   printf("\n");
   printf("Relative error of the Zolotorave approximation to 1/|x|\n");
   printf("-------------------------------------------------------\n\n");

   printf("Range [a,b] of |x| = ");

   ir=scanf("[");
   ir+=scanf("%lf",&a);
   ir+=scanf(",");
   ir+=scanf("%lf",&b);
   ir+=scanf("]");

   printf("Degree n = ");
   ir+=scanf("%d",&n);

   error(ir!=3,1,"main [table1.c]","Read error or invalid input data");   
   error((a<=0.0)||(b<=a)||(n<1),1,"main [table1.c]",
         "Improper range or degree");

   eps=a/b;
   alloc_ars();
   zolotarev(n,eps*eps,&As,ars,&delta);
   printf("Approximation error = %.1e\n\n",delta);
 
   dat=fopen("table1.dat","w");
   error(dat==NULL,1,"main [table1.c]","Unable to open output file");

   fprintf(dat,"#\n");
   fprintf(dat,"# Relative error of the Zolotarev rational approximation to\n");
   fprintf(dat,"# the function f(x)=1/|x|\n");
   fprintf(dat,"#\n");
   fprintf(dat,"# Approximation range %.2e<=|x|<=%.2e, degree = %d\n",a,b,n);   
   fprintf(dat,"# Approximation error = %.1e\n",delta);      
   fprintf(dat,"#\n");
   fprintf(dat,"#        x          error\n");

   kmax=100;
   dx=0.5*eps/(double)(kmax);
   
   for (k=0;k<kmax;k++)
   {
      x=0.5*eps+(double)(k)*dx;
      fprintf(dat,"  %.6e  % .6e\n",x*b,zolotarev_sign(x)-1.0);
   }

   kmax=1000;
   dx=(10.0*eps)/(double)(kmax);

   for (k=0;k<kmax;k++)
   {
      x=eps+(double)(k)*dx;
      fprintf(dat,"  %.6e  % .6e\n",x*b,zolotarev_sign(x)-1.0);
   }

   kmax=n*100;
   dx=(1.05-11.0*eps)/(double)(kmax);

   for (k=0;k<kmax;k++)
   {
      x=11.0*eps+(double)(k)*dx;
      fprintf(dat,"  %.6e  % .6e\n",x*b,zolotarev_sign(x)-1.0);
   }
   
   fclose(dat);

   printf("Printed data to file table1.dat\n\n");
   exit(0);
}
