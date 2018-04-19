
/*******************************************************************************
*
* File stat.c
*
* Copyright (C) 2005, 2011 Martin Luescher, Leonardo Giusti
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Collection of simple statistical analysis programs
*
* The externally accessible functions are 
*
*   double average(int n,double *a)
*     Returns the average of the array elements a[0],..,a[n-1]
*
*   double sigma0(int n,double *a)
*     Returns the naive statistical error of the average of the array
*     elements a[0],..,a[n-1]
*
*   double auto_corr(int n,double *a,int tmax,double *g)
*     Computes the normalized autocorrelation function g[t] at time
*     separations t=0,..,tmax-1 of the sequence a[0],..,a[n-1] and 
*     returns the value of the (unnormalized) autocorrelation function
*     at t=0. The inequality tmax<=n must be respected
*
*   void sigma_auto_corr(int n,double *a,int tmax,int lambda,double *eg)
*     Computes the statistical error eg[t] at time t=0,..,tmax-1 of the 
*     normalized autocorrelation function of the sequence a[0],..,a[n-1].
*     The choice of the summation cutoff lambda is not critical, but it
*     should be set to a value not smaller than a few times the integrated
*     autocorrelation time of the sequence (see the notes below). The
*     inequality 2*tmax+lambda-1<=n must be respected
*
*   double tauint(int n,double *a,int tmax,int lambda,int *w,double *sigma)
*     Returns an estimate of the integrated autocorrelation time of the
*     sequence a[0],..,a[n-1]. On exit the summation window determined by
*     the program is assigned to *w and an estimate of the statistical 
*     error on the calculated autocorrelation time is assigned to *sigma.
*     The parameter tmax sets an upper limit on the summation window and
*     the summation cutoff lambda should be set to a value not smaller than
*     a few times the integrated autocorrelation time (see the notes below).
*     The inequality 2*tmax+lambda-1<=n must be respected
*
*   double print_auto(int n,double *a)
*     Prints a table of the approximate integrated auto-correlation time
*     tau(w)=1/2+sum_{t=1}^w g[t] and the associated statistical error
*     sigma(w)=sigma0*sqrt{2*tau(w)}, where g[t] denotes the normalized
*     autocorrelation function of the sequence a[0],..,a[n-1]. On exit 
*     the program returns the average of the array elements
*
*   double jack_err(int nx,int n,double **a,double (*f)(int nx,double *x),
*                   int bmax,double *sig)
*     Computes the standard estimate of an arbitrary function f() of
*     nx primary stochastic variables x[k], k=0,..,nx-1, for a given 
*     sequence a[k][0],..,a[k][n-1] of values of these. The associated
*     jackknife errors sig[bs-1] for bin size bs=1,..,bmax are also
*     computed. On exit the program returns the standard estimate of
*     the function f()
*
*   double print_jack(int nx,int n,double **a,double (*f)(int nx,double *x))
*     Prints a table of the jackknife errors calculated by the program
*     jack_err(), together with the estimated integrated autocorrelation
*     times, as a function of the bin size bs. On exit the program returns
*     the standard estimate of the function f()
*
* Notes:
*
* For a recent discussion of statistical error estimation see
*
*  Ulli Wolff, Monte Carlo errors with less errors,
*  Comput. Phys. Commun. 156 (2004) 143 [hep-lat/0306017]
*
* The jackknife procedure is explained in appendix B of this paper and
* an improved Madras-Sokal formula is derived [eq.(42)] for the statistical
* error of the integrated autocorrelation time. The program tauint() makes
* use of this formula, while print_auto() uses the unimproved formula
*
* The standard estimate of a function f() of a vector x[0],..,x[nx-1] of
* primary stochastic variables is obtained by setting the arguments x[k]
* to the ensemble averages <x[k]>
*
* The computation of the autocorrelation function and the integrated
* autocorrelation time follows the lines of appendix A of
*
*  M. L"uscher, Schwarz-preconditioned HMC algorithm for two-flavor
*  lattice QCD, Comput. Phys. Commun. 165 (2005) 199 [hep-lat/0409106] 
*
* In particular, the summation cutoff lambda is introduced there and
* the selection of the summation window *w is explained
*
* The programs in this module may be used in MPI programs, but should then
* only be called from the root process
*
*******************************************************************************/

#define STAT_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "utils.h"
#include "extras.h"


double average(int n,double *a)
{
   int i;
   double abar;

   error_root(n<1,1,"average [stat.c]",
              "Argument n is out of range (should be at least 1)");   
   
   abar=0.0;
   
   for (i=0;i<n;i++)
      abar+=a[i];

   abar/=(double)(n);
   
   return abar;   
}


double sigma0(int n,double *a)
{
   int i;
   double abar,var,xn;

   error_root(n<2,1,"sigma0 [stat.c]",
              "Argument n is out of range (should be at least 2)");
 
   abar=average(n,a);
   var=0.0;
   
   for (i=0;i<n;i++)
      var+=((a[i]-abar)*(a[i]-abar));

   xn=(double)(n);
   var/=(xn*(xn-1.0));
   
   return sqrt(var);
}


double auto_corr(int n,double *a,int tmax,double *g)
{
   int t,k;
   double abar,g0,var;

   error_root((n<2)||(tmax<1)||(tmax>n),1,"auto_corr [stat.c]",
              "Argument n or tmax is out of range"); 
   
   abar=average(n,a);
   g0=sigma0(n,a);  

   if (g0<=(10.0*DBL_EPSILON*fabs(abar)))
   {
      g0=0.0;

      for (t=0;t<tmax;t++)
         g[t]=0.0;  
   }
   else
   {
      g0=g0*g0*(double)(n-1);  
      g[0]=1.0;
    
      for (t=1;t<tmax;t++)
      {
         var=0.0;
      
         for (k=0;k<(n-t);k++)
            var+=((a[k]-abar)*(a[k+t]-abar));

         g[t]=var/(g0*(double)(n-t));
      }
   }

   return g0;
}


static void sigma_corr(int n,int tmax,int lambda,double *g,double *eg)
{
   int k,t;
   double sm,r;

   if (g[0]==0.0)
   {
      for (t=0;t<tmax;t++)
         eg[t]=0.0;
   }
   else
   {
      eg[0]=0.0;   

      for (t=1;t<tmax;t++)
      {
         sm=0.0;
    
         for (k=1;k<=(t+lambda);k++)
         {
            r=g[k+t]+g[abs(k-t)]-2.0*g[t]*g[k];
            sm+=(r*r);
         }

         eg[t]=sqrt(sm/(double)(n));
      }
   }
}


void sigma_auto_corr(int n,double *a,int tmax,int lambda,double *eg)
{
   int tmaxx;
   double *g;

   tmaxx=2*tmax+lambda-1;
   
   error_root((n<2)||(tmax<1)||(lambda<1)||(tmaxx>n),1,
              "sigma_auto_corr [stat.c]",
              "Argument n, tmax or lambda is out of range");
  
   g=amalloc(tmaxx*sizeof(*g),3);
   error_root(g==NULL,1,"sigma_auto_corr [stat.c]",
              "Unable to allocate auxiliary array");
              
   auto_corr(n,a,tmaxx,g);
   sigma_corr(n,tmax,lambda,g,eg);

   afree(g);
}


double tauint(int n,double *a,int tmax,int lambda,int *w,double *sigma)
{
   int t,tmaxx;
   double tau,g0;
   double *g,*eg;

   tmaxx=2*tmax+lambda-1;
   
   error_root((n<2)||(tmax<1)||(lambda<1)||(tmaxx>n),1,"tauint [stat.c]",
              "Argument n, tmax or lambda is out of range");
   
   g=amalloc(tmaxx*sizeof(*g),3);
   eg=amalloc(tmax*sizeof(*eg),3);

   error_root((g==NULL)||(eg==NULL),1,"tauint [stat.c]",
              "Unable to allocate auxiliary arrays");  

   g0=auto_corr(n,a,tmaxx,g);
   sigma_corr(n,tmax,lambda,g,eg);

   tau=0.5;
   (*w)=1;
   (*sigma)=0.0;

   if (g0!=0.0)
   {
      for (t=1;t<tmax;t++)
      {
         tau+=g[t];
      
         if (g[t]<=eg[t])
            break;
      }

      error_root(t==tmax,1,"tauint [stat.c]",
                 "The summation window was not properly determined");
      (*w)=t;
      (*sigma)=tau*sqrt((2.0/(double)(n))*
                        ((double)(2*(*w)+1)-3.0*tau+1.0/(4.0*tau)));
   }

   afree(g);
   afree(eg);

   return tau;
}


double print_auto(int n,double *a)
{
   int w,dw,iw,tmax;
   double *ga,*ta,sig0,abar,sig,err,g0;

   if (n<30)
      tmax=n/10+1;
   else
      tmax=n/30+3;
   
   ga=amalloc(tmax*sizeof(double),3);
   ta=amalloc(tmax*sizeof(double),3);   

   error_root((ga==NULL)||(ta==NULL),1,"print_auto [stat.c]",
              "Unable to allocate auxiliary arrays");

   sig0=sqrt(1.0-1.0/(double)(n))*sigma0(n,a);
   abar=average(n,a);
   g0=auto_corr(n,a,tmax,ga);

   if (g0==0.0)
      printf("Statistical error is negligible, no data to print\n");
   else
   {
      ta[0]=0.5;

      for (w=1;w<tmax;w++)
      {
	 ta[w]=ta[w-1]+ga[w];

         if (ta[w]<=0.0)
            tmax=w;
      }
      
      printf(" window w       tau_int             sigma\n\n");
      dw=1;
      iw=0;
      
      for (w=0;w<tmax;)
      {
         if (w==0)
            err=0.0;
         else
            err=sqrt((4.0*(double)(w)+2.0)/(double)(n));
         
         sig=sig0*sqrt(2.0*ta[w]);
         
         printf("    %d    \t%.1f (%.1f) \t%.1e (%.1e)\n",
                w,ta[w],ta[w]*err,sig,0.5*sig*err);

         if ((w>=8)&&(iw>=4))
         {
            iw=0;
            dw*=2;
         }

         w+=dw;
         iw+=1;
      }
   }
   
   printf("\n");
   afree(ga);
   afree(ta);

   return abar;
}


static double javg(int nx,int n,double **a,double (*f)(int nx,double *x))
{
   int i;
   double *x,fbar;

   x=amalloc(nx*sizeof(*x),3);
   error_root(x==NULL,1,"javg [stat.c]","Unable to allocate auxiliary array");

   for (i=0;i<nx;i++)
      x[i]=average(n,a[i]);

   fbar=f(nx,x);
   afree(x);

   return fbar;
}


static void jbin(int n,int bs,double *a,double *b)
{
   int nbin,i,j;
   double abar,sm,r;

   nbin=n/bs;
   r=1.0/(double)(bs);

   for (i=0;i<nbin;i++)
   {
      sm=0.0;

      for (j=0;j<bs;j++)
         sm+=a[i*bs+j];

      b[i]=r*sm;
   }
   
   abar=average(nbin,b);   
   r=(double)(nbin);
   r=1.0/sqrt(r*(r-1.0));

   for (i=0;i<nbin;i++)
       b[i]=abar+r*(abar-b[i]);
}


static double jerr(int nx,int nbin,double **b,double (*f)(int nx,double *x),
                   double fbar)
{
   int i,j;
   double *x,var,r;

   x=amalloc(nx*sizeof(*x),3);
   error_root(x==NULL,1,"jerr [stat.c]","Unable to allocate auxiliary array");

   var=0.0;

   for (i=0;i<nbin;i++)
   {
      for (j=0;j<nx;j++)
         x[j]=b[j][i];

      r=f(nx,x)-fbar;
      var+=(r*r);
   }

   afree(x);

   return sqrt(var);
}


double jack_err(int nx,int n,double **a,double (*f)(int nx,double *x),
                int bmax,double *sig)
{
   int i,bs,nbin;
   double **b,*p,fbar;

   error_root((nx<1)||(bmax<1)||((2*bmax)>n),1,"jack_err [stat.c]",
              "Argument nx,n or bmax is out of range");
   
   b=amalloc(nx*sizeof(*b),3);
   p=amalloc(nx*n*sizeof(*p),3);

   error_root((b==NULL)||(p==NULL),1,"jack_err [stat.c]",
              "Unable to allocate auxiliary arrays");
   
   for (i=0;i<nx;i++)
   {
      b[i]=p;
      p+=n;
   }

   fbar=javg(nx,n,a,f);

   for (bs=1;bs<=bmax;bs++)
   {
      for (i=0;i<nx;i++)
         jbin(n,bs,a[i],b[i]);

      nbin=n/bs;
      sig[bs-1]=jerr(nx,nbin,b,f,fbar);
   }
   
   afree(p);
   afree(b);

   return fbar;
}

   
double print_jack(int nx,int n,double **a,double (*f)(int nx,double *x))
{
   int bs,dbs,ibs,bmax;
   double fbar,tau,*sig;

   if (n<30)
      bmax=n/10+1;
   else
      bmax=n/30+3;

   sig=amalloc(bmax*sizeof(*sig),3);
   error_root(sig==NULL,1,"print_jack [stat.c]",
              "Unable to allocate auxiliary array");

   fbar=jack_err(nx,n,a,f,bmax,sig);
   
   printf(" bin size     tau_int            sigma\n\n");
   dbs=1;
   ibs=0;
      
   for (bs=1;bs<=bmax;)
   {
      tau=sig[bs-1]/sig[0];
      tau=0.5*(tau*tau);
      
      printf("    %d    \t%.1f \t\t%.1e\n",bs,tau,sig[bs-1]);

      if ((bs>=8)&&(ibs>=4))
      {
         ibs=0;
         dbs*=2;
      }

      bs+=dbs;
      ibs+=1;
   }

   printf("\n");
   afree(sig);
   
   return fbar;
}
