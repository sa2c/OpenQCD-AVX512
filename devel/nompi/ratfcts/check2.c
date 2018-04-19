
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2008, 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the Jacobi elliptic functions sn,cn,dn
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "random.h"
#include "utils.h"
#include "ratfcts.h"


static void sncndn_smallk(double u,double rk,double *sn,double *cn,double *dn)
{
   int n,nmax;
   double K,Kp,pi;
   double t,v,tn,vn,k,r;
   
   K=ellipticK(rk);
   Kp=ellipticK(1.0/rk);
   pi=4.0*atan(1.0);
   
   t=pi*Kp/K;
   v=(pi*u)/(2.0*K);

   nmax=(int)(-1.5*log(DBL_EPSILON)/t);
   if (nmax<1)
      nmax=1;

   (*sn)=0.0;
   (*cn)=0.0;
   (*dn)=0.0;

   for (n=nmax;n>=0;n--)
   {
      tn=(double)(n)*t;
      vn=(double)(2*n+1)*v;

      (*sn)+=(exp(-tn)*sin(vn)/(1.0-exp(-2.0*tn-t)));
      (*cn)+=(exp(-tn)*cos(vn)/(1.0+exp(-2.0*tn-t)));

      if (n>0)
      {
         vn=(double)(2*n)*v;
         (*dn)+=(exp(-tn)*cos(vn)/(1.0+exp(-2.0*tn)));
      }
   }

   k=rk/sqrt(1.0+rk*rk);
   r=(2.0*pi*exp(-0.5*t))/(k*K);

   (*sn)*=r;
   (*cn)*=r;
   (*dn)=(pi/(2.0*K))*(1.0+4.0*(*dn));
}


static void sncndn_landen(double u,double rk,double *sn,double *cn,double *dn)
{
   double k,kp,kt,ktp,r;

   kp=1.0/sqrt(1.0+rk*rk);
   k=rk*kp;
   kt=k/(1.0+kp);
   kt=(kt*kt);
   ktp=2.0*sqrt(kp)/(1.0+kp);

   sncndn(u/(1.0+kt),kt/ktp,sn,cn,dn);

   r=1.0/(1.0+kt*(*sn)*(*sn));

   (*sn)=(1.0+kt)*(*sn)*r;
   (*cn)=(*cn)*(*dn)*r;
   (*dn)=sqrt(kp*kp+k*k*(*cn)*(*cn));
}


int main(void)
{
   int n;
   double u,rk,K;
   double sn,cn,dn,snsk,cnsk,dnsk;
   double dmax_sn,dmax_cn,dmax_dn,dev;
   
   printf("\n");
   printf("Computation of the Jacobi elliptic functions sn,cn,dn\n");
   printf("-----------------------------------------------------\n\n");

   rlxd_init(1,1234);   

   dmax_sn=0.0;
   dmax_cn=0.0;
   dmax_dn=0.0;
   
   for (n=0;n<10000;n++)
   {
      ranlxd(&rk,1);
      rk*=0.1;
      K=ellipticK(rk);
      ranlxd(&u,1);
      u=K*(0.5-u);

      sncndn(u,rk,&sn,&cn,&dn);      
      sncndn_smallk(u,rk,&snsk,&cnsk,&dnsk);

      if (sn!=0.0)
      {
         dev=fabs(1.0-snsk/sn);

         if (dev>dmax_sn)
            dmax_sn=dev;
      }

      dev=fabs(1.0-cnsk/cn);

      if (dev>dmax_cn)
         dmax_cn=dev;

      dev=fabs(1.0-dnsk/dn);

      if (dev>dmax_dn)
         dmax_dn=dev;       
   }

   printf("-K/2<=u<=K/2, rk<=0.1:\n");
   printf("maximal relative error (sn,cn,dn) = (%.1e,%.1e,%.1e)\n\n",
          dmax_sn,dmax_cn,dmax_dn);

   dmax_sn=0.0;
   dmax_cn=0.0;
   dmax_dn=0.0;
   
   for (n=0;n<10000;n++)
   {
      ranlxd(&rk,1);
      rk*=0.1;
      K=ellipticK(rk);
      ranlxd(&u,1);
      u=16.0*K*(0.5-u);

      sncndn(u,rk,&sn,&cn,&dn);      
      sncndn_smallk(u,rk,&snsk,&cnsk,&dnsk);

      dev=fabs(snsk-sn);

      if (dev>dmax_sn)
         dmax_sn=dev;

      dev=fabs(cnsk-cn);

      if (dev>dmax_cn)
         dmax_cn=dev;

      dev=fabs(dnsk-dn);

      if (dev>dmax_dn)
         dmax_dn=dev;       
   }

   printf("-8*K<=u<=8K, rk<=0.1:\n");
   printf("maximal absolute error (sn,cn,dn) = (%.1e,%.1e,%.1e)\n\n",
          dmax_sn,dmax_cn,dmax_dn);

   dmax_sn=0.0;
   dmax_cn=0.0;
   dmax_dn=0.0;
   
   for (n=0;n<10000;n++)
   {
      ranlxd(&rk,1);
      rk=rk/(1.0-rk);
      K=ellipticK(rk);
      ranlxd(&u,1);
      u=K*(0.5-u);

      sncndn(u,rk,&sn,&cn,&dn);      
      sncndn_landen(u,rk,&snsk,&cnsk,&dnsk);

      if (sn!=0.0)
      {
         dev=fabs(1.0-snsk/sn);

         if (dev>dmax_sn)
            dmax_sn=dev;
      }

      dev=fabs(1.0-cnsk/cn);

      if (dev>dmax_cn)
         dmax_cn=dev;

      dev=fabs(1.0-dnsk/dn);

      if (dev>dmax_dn)
         dmax_dn=dev;       
   }

   printf("-K/2<=u<=K/2, Landen recursion:\n");
   printf("maximal relative error (sn,cn,dn) = (%.1e,%.1e,%.1e)\n\n",
          dmax_sn,dmax_cn,dmax_dn);

   dmax_sn=0.0;
   dmax_cn=0.0;
   dmax_dn=0.0;
   
   for (n=0;n<10000;n++)
   {
      ranlxd(&rk,1);
      rk=rk/(1.0-rk);
      K=ellipticK(rk);
      ranlxd(&u,1);
      u=16.0*K*(0.5-u);
      
      sncndn(u,rk,&sn,&cn,&dn);      
      sncndn_landen(u,rk,&snsk,&cnsk,&dnsk);

      dev=fabs(snsk-sn);

      if (dev>dmax_sn)
         dmax_sn=dev;

      dev=fabs(cnsk-cn);

      if (dev>dmax_cn)
         dmax_cn=dev;

      dev=fabs(dnsk-dn);

      if (dev>dmax_dn)
         dmax_dn=dev;       
   }

   printf("-8*K<=u<=8K, Landen recursion:\n");
   printf("maximal absolute error (sn,cn,dn) = (%.1e,%.1e,%.1e)\n\n",
          dmax_sn,dmax_cn,dmax_dn);
   printf("Print values at specfied u and k/k'\n\n");
   
   for (;;)
   {
      printf("u, k/k' =  ");
      
      if (scanf("%lf %lf",&u,&rk)==2)
      {
         sncndn(u,rk,&sn,&cn,&dn);
         printf("sn = %.16e,  cn = %.16e, dn = %.16e\n",sn,cn,dn);
      }
      else
      {
         printf("Invalid input values, program stopped\n\n");
         break;
      }
   }
   
   exit(0);
}
