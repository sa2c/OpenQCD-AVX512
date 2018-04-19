
/*******************************************************************************
*
* File chebyshev.c
*
* Copyright (C) 2005, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Chebyshev approximation and integration
*
* The externally accessible functions are
*
*   int cheby_fit(double a,double b,double (*f)(double x),
*                                   int nmax,double eps,double c[])
*     Computes the coefficients c[0],...,c[n], with n<=nmax being the 
*     value returned by the program and eps the desired absolute precision 
*     of the approximation
*
*   double cheby_val(double a,double b,int n,double c[],double x)
*     Computes the value of the Chebyshev approximation at x, assuming
*     the coefficients c_k are stored in the array c[0],...,c[n]
*
*   double cheby_int(double a,double b,double (*f)(double x),
*                                      int nmax,double eps)
*     Computes the definite integral of f(x) in the range a<=x<=b to an
*     absolute precision eps, using Chebyshev polynomials of degree n<=nmax
*
* Notes:
*
* For the numerical approximation and integration of a given function f(x),
* using the Chebyshev polynomials
*
*   T_k(z)=cos(k*theta),  z=cos(theta),  -1<=z<=1,
*
* the function is assumed to be defined in the range a<=x<=b and to be
* available as a function program. The approximation is then of the form
*
*   f(x)=sum{c_k*T_k(z),k=0..n}+r(x),   z=(a+b-2*x)/(a-b)
*
*   |r(x)|<eps for all x
*
* A detailed description is given in the notes
*
*   M. Luescher: Chebyshev approximation and integration
*
*******************************************************************************/

#define CHEBYSHEV_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "extras.h"

static int max_degree;
static double *alist,*clist,*flist;


static void allocate_arrays(int nmax)
{
   for (max_degree=16;max_degree<=nmax;)
      max_degree*=2;

   alist=malloc((max_degree+1)*sizeof(double));
   clist=malloc((max_degree*2)*sizeof(double));
   flist=malloc((max_degree+1)*sizeof(double));
}


static void free_arrays(void)
{
   free(alist);
   free(clist);
   free(flist);
}


static void update_clist(int n)
{
   int k,kmin,kmax,dk;
   double pi,x,dx;

   pi=4.0*atan(1.0);
   dx=pi/(double)(max_degree);

   kmin=0;
   kmax=2*max_degree;
   dk=max_degree/n;
   
   if (n>32)
   {
      kmin=dk;
      dk*=2;
   }
   
   for (k=kmin;k<kmax;k+=dk)
   {
      x=(double)(k)*dx;
      clist[k]=cos(x);
   }
}  


static void update_flist(int n,double a,double b,double (*f)(double x))
{
   int k,kmin,kmax,dk;
   double x;

   kmin=0;
   kmax=max_degree;
   dk=max_degree/n;

   if (n>32)
   {
      kmin=dk;
      kmax-=dk;
      dk*=2;
   }

   for (k=kmin;k<=kmax;k+=dk)
   {
      x=0.5*(a+b-(a-b)*clist[k]);
      flist[k]=(*f)(x);
   }
}


static void compute_alist(int n)
{
   int i,k,dk;
   double sum,r;
   
   dk=max_degree/n;
   r=2.0/(double)n;

   for (i=0;i<=n;++i)
   {
      sum=0.5*(flist[0]+flist[max_degree]);
      if (i%2==1)
         sum-=flist[max_degree];

      for (k=dk;k<max_degree;k+=dk)
         sum+=flist[k]*clist[(i*k)%(2*max_degree)];

      alist[i]=r*sum;
   }
}


static void compute_blist(int n,double a,double b)
{
   int i,k,dk;
   double sum,r;
   
   dk=max_degree/n;
   r=2.0/(double)n;

   for (i=0;i<=n;++i)
   {
      if (i%2==0)
      {
         sum=0.5*(flist[0]+flist[max_degree]);

         for (k=dk;k<max_degree;k+=dk)
            sum+=flist[k]*clist[(i*k)%(2*max_degree)];

         alist[i]=(a-b)*r*sum/(double)(i*i-1);
      }
      else
      {
         alist[i]=0.0;
      }
   }

   alist[0]*=0.5;
}


static int test_convergence(int n)
{
   int i,k,kmax;
   double m[4],a;

   kmax=n/4;

   for (i=0;i<4;++i)
   {
      m[i]=0.0;
      
      for (k=0;k<kmax;++k)
      {
         a=fabs(alist[i*kmax+k]);
         if (a>m[i])
            m[i]=a;
      }
   }

   if ((m[0]>=1.0e2*m[1])&&(m[0]>=1.0e4*m[2])&&(m[0]>=1.0e6*m[3]))   
      return(0);

   return(1);
}


static double abs_error(int n)
{
   int k,kmin;
   double err;

   kmin=n/2+1;
   err=0.0;

   for (k=0;k<kmin;++k)
      err+=fabs(alist[k]);
   err*=(DBL_EPSILON);
   
   for (k=kmin;k<n;++k)
      err+=fabs(alist[k]);

   return(2.0*err);
}


static int economize(int n,double eps,double err,double c[])
{
   int k;
   double r;

   r=err;

   for (k=n;k>=1;--k)
   {
      r+=fabs(c[k]);
      if (r>=eps)
         break;
   }

   return(k);
}


int cheby_fit(double a,double b,double (*f)(double x),
              int nmax,double eps,double c[])
{
   int n,k,itest;
   double err;
   
   if ((a>=b)||(nmax<16)||(eps<=0.0))
   {
      printf("Error in cheby_fit\n");
      printf("Arguments out of range\n");
      printf("Program aborted\n\n");
      exit(0);
   }

   itest=1;
   err=eps;
   allocate_arrays(nmax);

   for (n=32;n<=max_degree;n*=2)
   {
      update_clist(n);
      update_flist(n,a,b,f);
      compute_alist(n);

      itest=test_convergence(n);
      err=abs_error(n);

      if ((itest==0)&&(err<eps))
      {
         n/=2;
         c[0]=0.5*alist[0];
         for (k=1;k<=n;++k)
            c[k]=alist[k];
         break;
      }
   }   

   free_arrays();

   if ((itest!=0)||(err>=eps))
   {
      printf("Error in cheby_fit\n");
      printf("Specified accuracy has not been reached\n");
      printf("Program aborted\n\n");
      exit(0);
   }

   n=economize(n,eps,err,c); 
   return(n);
}


double cheby_val(double a,double b,int n,double c[],double x)
{
   int k;
   double u,v,w,z;

   if ((n<0)||(a>=b)||(x>b)||(x<a))
   {
      printf("Error in cheby_val\n");
      printf("Arguments out of range\n");
      printf("Program aborted\n\n");
      exit(0);
   }

   if (n==0)
      return(c[0]);

   z=2.0*(a+b-2.0*x)/(a-b);
   u=c[n];
   v=c[n-1];

   for (k=n-2;k>=0;--k)
   {
      w=z*u+v;
      v=c[k]-u;
      u=w;
   }

   return(0.5*z*u+v);
}


double cheby_int(double a,double b,double (*f)(double x),
                 int nmax,double eps)
{
   int n,k,itest;
   double err,sum;
   
   if ((a>=b)||(nmax<16)||(eps<=0.0))
   {
      printf("Error in cheby_int\n");
      printf("Arguments out of range\n");
      printf("Program aborted\n\n");
      exit(0);
   }

   itest=1;
   err=eps;
   sum=0.0;
   allocate_arrays(nmax);

   for (n=32;n<=max_degree;n*=2)
   {
      update_clist(n);
      update_flist(n,a,b,f);
      compute_blist(n,a,b);

      itest=test_convergence(n);
      err=abs_error(n);

      if ((itest==0)&&(err<eps))
      {
         n/=2;
         for (k=n;k>=0;k-=2)
            sum+=alist[k];
         break;
      }
   }   

   free_arrays();

   if ((itest!=0)||(err>=eps))
   {
      printf("Error in cheby_int\n");
      printf("Specified accuracy has not been reached\n");
      printf("Program aborted\n\n");
      exit(0);
   }

   return(sum);
}
