
/*******************************************************************************
*
* File pchi_square.c
*
* Copyright (C) 2005, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Chi-square probability distribution
*
* The externally accessible function is
*
*   double pchi_square(double chi_square,int nu)
*     For chi_square>=0 and nu=1,2,...,1000 the program returns an 
*     approximation for P(chi_square|nu) which deviates from the exact
*     distribution by less than 10^(-8) [10^(-9) if nu=1]
*
* Notes:
*
* See the notes
*
*   M. Luescher: Statistical tests
*   
* for a detailed description.
*
*******************************************************************************/

#define PCHI_SQUARE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "extras.h"

static int init=0;
static double c0,c1,c2,c3,c4,c5,lng[40];
static double xd0,xd1,xd2,xd3,xd4,pi;


static void define_constants(void)
{
   int n;
   double x;

   xd0=0.0;
   xd1=1.0;
   xd2=2.0;
   xd3=3.0;
   xd4=4.0;

   pi=xd4*atan(xd1);

   c0=-xd1/xd2;
   c1=log(xd2*pi)/xd2;
   c2= xd1/12.0;
   c3=-xd1/360.0;
   c4= xd1/1260.0;
   c5=-xd1/1680.0;
           
   lng[1]=log(pi)/xd2;
   lng[2]=xd0;

   for (n=3;n<40;++n)
   {           
      x=(double)(n-2);
      lng[n]=lng[n-2]+log(x/xd2);
   }

   init=1;
}
  

static double ln_gamma(int k)
{
   double y,z,zm1,zm2;

   if (k<40)
      return lng[k];

   z=(double)k;
   z=z/xd2;
   zm1=xd1/z;
   zm2=zm1*zm1;

   y=c5;
   y=y*zm2+c4;
   y=y*zm2+c3;
   y=y*zm2+c2;
   y=y*zm1+c1;

   return (z+c0)*log(z)-z+y;
}


static double pchi1(double chi_square)
{
   double x,y,z,a,p;

   x=chi_square;

   if (x<=1.0e-18)
      return xd0;
   if (x>=40.0)
      return xd1;

   z=x/xd2;
   a=xd2*sqrt(z/pi)*exp(-z);
        
   y=xd0;
   z=xd3;
   p=x/z;

   for (;a>(xd1-p)*1.0e-9;)
   {
      y+=a;
      a*=p;
      z+=xd2;
      p=x/z;
   }
   return y;
}


static double pchi2(double chi_square,int nu)
{
   double x,y,z,xnu,lna,a,p;

   x=chi_square;
   xnu=(double)nu;

   if (x<=1.0e-18)
      return xd0;

   if (x<=xnu)
   {
      z=x/xd2;
      lna=(xnu/xd2)*log(z)-z-ln_gamma(nu+2);

      z=xnu+xd2;
      p=x/z;
      y=xd0;

      if ((lna-log(xd1-p))<-18.5)
         return y;
        
      a=exp(lna);

      for (;a>(xd1-p)*1.0e-8;)
      {
         y+=a;
         a*=p;
         z+=xd2;
         p=x/z;
      }
      return y;
   }
   else
   {
      z=x/xd2;
      lna=((xnu/xd2-xd1)*log(z)-z)-ln_gamma(nu);

      z=xnu-xd2;
      p=z/x;
      if (nu%2==1)
         y=pchi1(x);
      else
         y=xd1;

      if ((lna-log(xd1-p))<-18.5)
         return y;
        
      a=exp(lna);

      for (;(z>=xd0)&&(a>(xd1-p)*9.0e-9);)
      {
         y-=a;
         a*=p;
         z-=xd2;
         p=z/x;
      }
      return y;
   }
}


double pchi_square(double chi_square,int nu)
{
   if (init==0)
      define_constants();

   if ((nu<1)||(nu>1000)||(chi_square<xd0))
   {
      printf("Error in pchi_square: argument out of range\n");
      printf("Program aborted\n\n");
      exit(0);
   }

   if (nu==1)
      return pchi1(chi_square);
   else
      return pchi2(chi_square,nu);
}

