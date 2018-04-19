
/*******************************************************************************
*
* File fsolve.c
*
* Copyright (C) 2008, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* General purpose equation solver and function minimizers
*
* The externally accessible functions are
*
*   double inverse_fct(double y,double x1,double x2,double (*f)(double x),
*                      double omega1,double omega2)   
*     Finds a solution x of the equation f(x)=y in the interval [x1,x2]
*     to an absolute precision omega1 or a relative precision omega2
*     (whichever is reached first). The points x1,x2 must be such that
*     f(x1) and f(x2) have different sign 
*
*   double minimize_fct(double x0,double x1,double x2,double (*f)(double x),
*                       double omega1,double omega2)
*     Finds a local minimum x of f(x) in the interval [x0,x2] to an
*     absolute precision omega1 or a relative precision omega2 (whichever
*     is reached first). The point x1 is taken as an initial guess of the
*     position of the minimum (x0<x1<x2)
*
*   void powell(int n,double *x0,double *x1,double *x2,
*               double (*f)(int n,double *x),int imx,
*               double omega1,double omega2,double *xmin,int *status)
*     Finds a local minimum xmin of a given function f that depends on
*     on a vector x[0],..,x[n-1] of n variables. The minimum is searched
*     for in the hypercube x0[j]<x[j]<x2[j], j=0,..,n-1, starting from
*     x=x1 (which must be in the hypercube). At most imx iterations of
*     Powell's direction set method are applied to find the minimum.
*      The program terminates if the coordinates of the position xmin
*     changed by less than omega1 or less than omega2*xmin in the last
*     iteration. On output status reports the total number of iterations
*     that were required or a negative number if the program failed (-1
*     if the algorithm did not converge, -2 if the minimum could not be
*     bracketed)
*
* Notes:
*
* The program inverse_fct() uses a slightly modified secant method, while
* minimize_fct() proceeds according to a golden-ratio bisection method.
*
* The meaningful levels of precision depend on the function f() and the
* machine precision. In particular, rounding becomes important in the
* minimization routine at a precision roughly equal to the square root
* of the machine precision.
*
* Powell's method is described in Chapter 10 of
* 
*   W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery,
*   Numerical Recipes in FORTRAN, 2nd Edition 
*   (Cambridge University Press, Cambridge, 1992)
*
* The variant recommended in this book is implemented here with some small
* modifications
*
*******************************************************************************/

#define FSOLVE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "utils.h"
#include "extras.h"

static int nsv=0,isv;
static double *osv,*psv,**vsv,*xsv;
static double (*fsv)(int n,double *x);


static int relative_sign(double f1,double f2)
{
   if (((f1>=0.0)&&(f2<=0.0))||((f1<=0.0)&&(f2>=0.0)))
      return 1;
   else
      return 0;
}


double inverse_fct(double x1,double x2,double (*f)(double x),double y,
		   double omega1,double omega2)
{
   double x3,f1,f2,f3,dx;
   double lambda,eps;

   f1=f(x1)-y;
   f2=f(x2)-y;
      
   error((x1>x2)||(relative_sign(f1,f2)==0),1,"inverse_fct [fsolve.c]",
         "Improper bracket [x1,x2]");

   eps=0.1;
   omega2*=0.5;
   dx=x2-x1;   

   while ((dx>omega1)&&(dx>(omega2*(x1+x2))))
   {
      if (fabs(f1)<fabs(f2))
      {
         lambda=f1/(f1-f2);
         if (lambda<eps)
            lambda=eps;
         x3=x1+dx*lambda;
      }
      else
      {
         lambda=f2/(f2-f1);
         if (lambda<eps)
            lambda=eps;
         x3=x2-dx*lambda;
      }

      f3=f(x3)-y;

      if (relative_sign(f1,f3)==1)
      {
         x2=x3;
         f2=f3;
      }
      else
      {
         x1=x3;
         f1=f3;
      }

      dx=x2-x1;
   }

   if (fabs(f1)<fabs(f2))
      return x1;
   else
      return x2;
}


static int find_bracket(double *x0,double *x1,double *x2,double (*f)(double x))
{
   int ic;
   double y0,y1,y2,d0,d2;
   double f0,f1,f2;

   y0=(*x0);
   y1=(*x1);
   y2=(*x2);
   d0=0.05*(y0-y1);
   d2=0.05*(y2-y1);

   (*x0)=y1+d0;
   (*x2)=y1+d2;   

   f0=f(*x0);
   f1=f(*x1);
   f2=f(*x2);   

   if ((f0>f1)&&(f2>f1))
      return 0;

   if ((f(y0)>f1)&&(f(y2)>f1))
   {
      (*x0)=y0;
      (*x2)=y2;
      return 0;
   }
   
   for (ic=1;ic<20;ic++)
   {
      if (f1>f2)
      {
         (*x2)=(*x1);
         f2=f1;
      }

      (*x1)=(*x0);
      f1=f0;

      (*x0)+=d0;

      if ((*x0)<y0)
         (*x0)=y0;
      
      f0=f(*x0);

      if ((f0>f1)&&(f2>f1))
         return 0;         
   }
   
   (*x0)=y1+d0;
   (*x1)=y1;
   (*x2)=y1+d2;  

   f0=f(*x0);
   f1=f(*x1);
   f2=f(*x2);  
   
   for (ic=1;ic<20;ic++)
   {
      if (f0<f1)
      {
         (*x0)=(*x1);
         f0=f1;
      }
      
      (*x1)=(*x2);
      f1=f2;
      
      (*x2)+=d2;

      if ((*x2)>y2)
         (*x2)=y2;
      
      f2=f(*x2);

      if ((f0>f1)&&(f2>f1))
         return 0;
   }   

   return 1;
}


static double mini_fct(double x0,double x1,double x2,double (*f)(double x),
                       double omega1,double omega2)
{
   double s,x3,f1,f2,dx;

   omega2*=0.5;
   s=0.5*(3.0-sqrt(5.0));
   x3=x2;
   dx=x3-x0;
   f1=f(x1);
   
   if ((x1-x0)<(x3-x1))
   {
      x2=x1+s*(x3-x1);
      f2=f(x2);
   }
   else
   {
      x2=x1;
      f2=f1;
      x1=x2-s*(x2-x0);
      f1=f(x1);
   }

   while ((dx>omega1)&&(dx>(omega2*(fabs(x0)+fabs(x3)))))
   {
      if (f1<f2)
      {
         x3=x2;
         x2=x1;
         f2=f1;
         x1=x2-s*(x2-x0);
         f1=f(x1);
      }
      else
      {
         x0=x1;
         x1=x2;
         f1=f2;
         x2=x1+s*(x3-x1);
         f2=f(x2);
      }

      dx=x3-x0;
   }

   if (f1<f2)
      return x1;
   else
      return x2;
}


double minimize_fct(double x0,double x1,double x2,double (*f)(double x),
		    double omega1,double omega2)
{
   error((x1<=x0)||(x1>=x2),1,"minimize_fct [fsolve.c]",
         "Improper input values x0,x1,x2");

   error(find_bracket(&x0,&x1,&x2,f),1,"minimize_fct [fsolve.c]",
         "Unable to bracket minimum");   

   return mini_fct(x0,x1,x2,f,omega1,omega2);
}


static void alloc_arrays(int n)
{
   int i,j;
   
   if (nsv!=n)
   {
      if (nsv!=0)
      {
         afree(osv);
         afree(vsv);
      }

      if (n>0)
      {
         osv=amalloc(n*(n+3)*sizeof(*psv),3);
         vsv=amalloc(n*sizeof(*vsv),3);

         error((osv==NULL)||(vsv==NULL),1,"alloc_arrays [fsolve.c]",
               "Unable to allocate auxiliary arrays");

         psv=osv+n;
         vsv[0]=psv+n;

         for (i=1;i<n;i++)
            vsv[i]=vsv[i-1]+n;

         xsv=vsv[n-1]+n;
         nsv=n;
      }
      else
      {
         osv=NULL;
         vsv=NULL;
         xsv=NULL;
         nsv=0;
      }
   }
   
   for (i=0;i<n;i++)
   {
      for (j=0;j<n;j++)
         vsv[i][j]=0.0;

      vsv[i][i]=1.0;
   }
}


static void find_bnds(double *x0,double *x2,double *r0,double *r2,
                      double *rom1,double *rom2)
{
   int j,k;
   double *v,vmax,va,pa;

   k=0;
   v=vsv[isv];
   vmax=fabs(v[0]);

   for (j=1;j<nsv;j++)
   {
      va=fabs(v[j]);
      
      if (vmax<va)
      {
         k=j;
         vmax=va;
      }
   }

   (*rom1)=1.0/vmax;
   (*rom2)=fabs(psv[k])*(*rom1);
   
   if (v[k]>0.0)
   {
      (*r0)=(x0[k]-psv[k])/v[k];
      (*r2)=(x2[k]-psv[k])/v[k];   
   }
   else
   {
      (*r0)=(x2[k]-psv[k])/v[k]; 
      (*r2)=(x0[k]-psv[k])/v[k];
   }

   for (j=0;j<nsv;j++)
   {
      pa=fabs(psv[j]);
      va=fabs(v[j]);

      if (((*rom2)*va)>pa)
         (*rom2)=pa/va;

      if (v[j]>0.0)
      {
         if ((psv[j]+(*r0)*v[j])<x0[j])
            (*r0)=(x0[j]-psv[j])/v[j];
         if ((psv[j]+(*r2)*v[j])>x2[j])
            (*r2)=(x2[j]-psv[j])/v[j];         
      }
      else
      {
         if ((psv[j]+(*r0)*v[j])>x2[j])
            (*r0)=(x2[j]-psv[j])/v[j];
         if ((psv[j]+(*r2)*v[j])<x0[j])
            (*r2)=(x0[j]-psv[j])/v[j];         
      }         
   }
}


static double fline(double r)
{
   int j;

   for (j=0;j<nsv;j++)
      xsv[j]=psv[j]+r*vsv[isv][j];

   return fsv(nsv,xsv);
}


void powell(int n,double *x0,double *x1,double *x2,
            double (*f)(int n,double *x),int imx,double omega1,double omega2,
            double *xmin,int *status)
{
   int i,j,k,ifn,io1,io2;
   double r0,r1,r2,del;
   double rom1,rom2;
   double fo,fp,fe;

   ifn=0;

   for (j=0;j<n;j++)
   {
      if ((x0[j]>=x1[j])||(x2[j]<=x1[j]))
         ifn=1;
   }

   error(ifn,1,"powell [fsolve.c]","Improper parameter arrays x0,x1,x2");
   error((imx<4)||((omega1<=0.0)&&(omega2<=0.0)),1,"powell [fsolve.c]",
         "Improper parameters imx,omega1 or omega2");
   
   fsv=f;   
   alloc_arrays(n);

   for (j=0;j<n;j++)
   {
      osv[j]=x1[j];
      psv[j]=x1[j];
   }

   fo=f(n,osv);
   fp=fo;

   for (i=0;i<imx;i++)
   {
      del=0.0;
      k=0;

      for (isv=0;isv<n;isv++)
      {
         r1=0.0;         
         find_bnds(x0,x2,&r0,&r2,&rom1,&rom2);
         ifn=find_bracket(&r0,&r1,&r2,fline);

         if (ifn==0)
         {
            r1=mini_fct(r0,r1,r2,fline,rom1*omega1,rom2*omega2);

            for (j=0;j<n;j++)
               psv[j]+=r1*vsv[isv][j];

            fe=f(n,psv);
            r0=fabs(fe-fp);

            if (r0>del)
            {
               del=r0;
               k=isv;
            }

            fp=fe;
         }
         else if (i>=2)
         {
            for (j=0;j<n;j++)
               xmin[j]=psv[j];

            (*status)=-2;
            return;
         }
      }

      ifn=0;
      
      for (j=0;j<n;j++)
      {
         xsv[j]=psv[j]+(psv[j]-osv[j]);

         if ((xsv[j]<=x0[j])||(xsv[j]>=x2[j]))
            ifn=1;
      }

      if (ifn==0)
      {
         fe=f(n,xsv);
         r0=fe-fo;
         r1=fo-fp-del;
         r2=2.0*(fo-2.0*fp+fe)*r1*r1-del*r0*r0;

         if ((r0<(-4.0*DBL_EPSILON*fabs(fo)))&&(r2<0.0))
         {
            for (j=0;j<n;j++)
            {
               if (k>0)
                  vsv[k][j]=vsv[0][j];
               vsv[0][j]=psv[j]-osv[j];
            }
         }
      }

      io1=1;
      io2=1;
      
      for (j=0;j<n;j++)
      {
         r0=fabs(psv[j]-osv[j]);

         if (r0>omega1)
            io1=0;
         if (r0>(omega2*psv[j]))
            io2=0;

         osv[j]=psv[j];
      }

      fo=fp;

      if ((i>=3)&&((io1==1)||(io2==1)))
         break;
   }

   for (j=0;j<n;j++)         
      xmin[j]=psv[j];

   if (i<imx)
      (*status)=i+1;
   else
      (*status)=-1;
}
