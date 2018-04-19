
/*******************************************************************************
*
* File ratfcts.c
*
* Copyright (C) 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Rational function coefficients data base
*
* The externally accessible functions are
*
*   ratfct_t ratfct(int *irat)
*     Returns a structure containing the coefficients of the rational
*     function specified by the integers irat[3] (see the notes).
*
* Notes:
*
* Currently only fractions of Zolotorev rational functions are supported.
* The parameters of the available Zolotarev rational functions can be
* retrieved from the parameter data base (see flags/rat_parms.c). 
*
* Zolotarev rational functions R(y) approximate 1/sqrt(y) in the range
* eps<=y<=1 for some value eps>0. The functions provided by this module
* instead approximate the function 1/|x| in a range ra<=|x|<=rb specified 
* in the parameter data base. The relation between x and y is
*
*   y=x^2/rb^2
*
* and thus eps=(ra/rb)^2. 
*
* The coefficients a[r], r=0,..,2*n-1, returned by the program zolotarev()
* are ordered such that
*
*   a[0]>a[1]>..>a[2*n-1]>0.
*
* For any given integers k,l satisfying k>=0 and k<=l<n, the subset
*
*   a[2*k],..,a[2*l+1]    
*
* of coefficients defines a fraction
*
*   r(y)=p(y)/q(y),
*
*   p(y)=(y+a[2*k])*(y+a[2*k+2])*..*(y+a[2*l]),
*
*   q(y)=(y+a[2*k+1])*(y+a[2*k+3])*..*(y+a[2*l+1]),
*
* of the Zolotarev rational function, which is completely specified by
*
*   irat[0]       Index of the Zolotarev rational function in the
*                 parameter data base.
*
*   irat[1]       Lower end k of the selected coefficient range.       
*
*   irat[2]       Upper end l of the selected coefficient range.
*
* This module administrates the coefficients of up to 32 such rational
* functions.
*
* The members of the ratfct_t structure are
*
*   np        Number of poles of r(y) (that is, np=l-k+1).
*
*   A         Normalization factor.
*
*   delta     Relative approximation error.
*
*   mu        Array of the coefficients rb*sqrt{a[2*k+2*j+1]}, j=0,..,np-1.
*
*   rmu       Array of the associated residues.
*
*   nu        Array of the coefficients rb*sqrt{a[2*k+2*j]}, j=0,..,np-1.
*
*   rnu       Array of the associated residues.
*
* The normalization factor A and the approximation error delta refer to the
* full Zolotarev rational function (corresponding to k=0,l=n-1), where
*
*   |1-A*|x|*r(y)|<=delta
*
* when x is in the approximation range.
*
* In all cases, the coefficient arrays are such that
*
*   r(y)=Prod{(x^2+nu[j]^2)/(x^2+mu[j]^2),j=0,..,np-1}
*
*       =1+Sum{rmu[j]/(x^2+mu[j]^2),j=0,..,np-1}
*
* Moreover, the partial fraction decomposition 
*
*   Prod{(x+i*mu[j])/(x+i*nu[j]),j=0,..,np-1)}
*
*       =1+i*Sum{rnu[j]/(x+i*nu[j]),j=0,..,np-1)}
*
* holds.
*   
* The program in this module may perform global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define RATFCTS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "flags.h"
#include "utils.h"
#include "ratfcts.h"

#define IRMAX 32

static int init=0,ns,irs,irats[IRMAX][3];
static double *ars;
static ratfct_t rats[IRMAX]={{0,0.0,1.0,NULL,NULL,NULL,NULL}};


static void init_rat(void)
{
   int ir;
   
   for (ir=0;ir<IRMAX;ir++)
   {
      if (ir>0)
         rats[ir]=rats[0];

      irats[ir][0]=0;
      irats[ir][1]=0;
      irats[ir][2]=0;
   }

   ns=0;
   irs=0;
   ars=NULL;
   init=1;
}


static int fnd_rat(int *irat)
{
   int ir;

   for (ir=0;ir<irs;ir++)
   {
      if ((irat[0]==irats[ir][0])&&(irat[1]==irats[ir][1])&&
          (irat[2]==irats[ir][2]))
         return ir;
   }

   return irs;
}


static void alloc_rat(int *irat)
{
   int n,np,k,l;
   double *mu;
   rat_parms_t rp;
   
   error(irs==IRMAX,1,"alloc_rat [ratfcts.c]",
         "Attempt to define more than %d rational functions",IRMAX);

   rp=rat_parms(irat[0]);
   n=rp.degree;
   k=irat[1];
   l=irat[2];
   np=l-k+1;
   
   error((k<0)||(l<k)||(l>=n),1,"alloc_rat [ratfcts.c]",
         "Improper coefficient range or undefined rational function");

   if (n>ns)
   {
      if (ns>0)
         free(ars);
      ars=malloc(2*n*sizeof(*ars));
      ns=n;
   }

   mu=malloc(4*np*sizeof(*mu));   

   error((ars==NULL)||(mu==NULL),1,"alloc_rat [ratfcts.c]",
         "Unable to allocate coefficient arrays");
   
   rats[irs].np=np;
   rats[irs].mu=mu;
   rats[irs].rmu=mu+np;
   rats[irs].nu=mu+2*np;
   rats[irs].rnu=mu+3*np;

   irats[irs][0]=irat[0];
   irats[irs][1]=irat[1];
   irats[irs][2]=irat[2];               
}


static void set_rat(int *irat)
{
   int n,np,k,l,i,j;
   double ra,rb,pmu,pnu;
   double eps,A,delta,*ar;
   double *mu,*nu,*rmu,*rnu;
   rat_parms_t rp;

   rp=rat_parms(irat[0]);
   n=rp.degree;
   k=irat[1];
   l=irat[2];
   np=l-k+1;

   ra=rp.range[0];
   rb=rp.range[1];
   eps=ra/rb;
   eps=eps*eps;

   zolotarev(n,eps,&A,ars,&delta);
   rats[irs].A=A/rb;
   rats[irs].delta=delta;

   ar=ars+2*k;
   mu=rats[irs].mu;
   nu=rats[irs].nu;
   rmu=rats[irs].rmu;
   rnu=rats[irs].rnu;

   for (i=0;i<np;i++)
   {
      mu[i]=rb*sqrt(ar[2*i+1]);
      nu[i]=rb*sqrt(ar[2*i]);
   }

    for (i=0;i<np;i++)
   {  
      pmu=1.0;
      pnu=1.0;

      for (j=0;j<np;j++)
      {
         if (j!=i)
         {
            pmu*=((ar[2*j]-ar[2*i+1])/(ar[2*j+1]-ar[2*i+1]));
            pnu*=((mu[j]-nu[i])/(nu[j]-nu[i]));
         }
      }

      rmu[i]=rb*rb*(ar[2*i]-ar[2*i+1])*pmu;
      rnu[i]=(mu[i]-nu[i])*pnu;
   }

   irs+=1;
}


ratfct_t ratfct(int *irat)
{
   int ir;
   
   if (init==0)
      init_rat();

   ir=fnd_rat(irat);

   if (ir==irs)
   {
      alloc_rat(irat);
      set_rat(irat);
   }

   return rats[ir];
}
