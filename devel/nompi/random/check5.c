
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Statistical test of random_su3
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "utils.h"
#include "su3fcts.h"
#include "extras.h"


static void dev(su3 *u,double *d1,double *d2)
{
   int i;
   float *r;
   double d;
   complex det1,det2,det3,det;
   su3 v,w;

   _su3_dagger(v,(*u));
   _su3_times_su3(w,v,(*u));

   w.c11.re-=1.0f;
   w.c22.re-=1.0f;
   w.c33.re-=1.0f;

   *d1=0.0;
   r=(float*)(&w);

   for (i=0;i<18;i++)
   {
      d=fabs((double)(r[i]));
      if (d>(*d1))
         *d1=d;
   }

   det1.re=
      ((*u).c22.re*(*u).c33.re-(*u).c22.im*(*u).c33.im)-
      ((*u).c23.re*(*u).c32.re-(*u).c23.im*(*u).c32.im);
   det1.im=
      ((*u).c22.re*(*u).c33.im+(*u).c22.im*(*u).c33.re)-
      ((*u).c23.re*(*u).c32.im+(*u).c23.im*(*u).c32.re);
   det2.re=
      ((*u).c21.re*(*u).c33.re-(*u).c21.im*(*u).c33.im)-
      ((*u).c23.re*(*u).c31.re-(*u).c23.im*(*u).c31.im);
   det2.im=
      ((*u).c21.re*(*u).c33.im+(*u).c21.im*(*u).c33.re)-
      ((*u).c23.re*(*u).c31.im+(*u).c23.im*(*u).c31.re);    
   det3.re=
      ((*u).c21.re*(*u).c32.re-(*u).c21.im*(*u).c32.im)-
      ((*u).c22.re*(*u).c31.re-(*u).c22.im*(*u).c31.im);
   det3.im=
      ((*u).c21.re*(*u).c32.im+(*u).c21.im*(*u).c32.re)-
      ((*u).c22.re*(*u).c31.im+(*u).c22.im*(*u).c31.re);

   det.re=
      ((*u).c11.re*det1.re-(*u).c11.im*det1.im)-
      ((*u).c12.re*det2.re-(*u).c12.im*det2.im)+
      ((*u).c13.re*det3.re-(*u).c13.im*det3.im);
   det.im=
      ((*u).c11.re*det1.im+(*u).c11.im*det1.re)-
      ((*u).c12.re*det2.im+(*u).c12.im*det2.re)+
      ((*u).c13.re*det3.im+(*u).c13.im*det3.re);

   *d2=0.0;
   d=fabs((double)(det.re)-1.0);
   if (d>(*d2))
      *d2=d;
   d=fabs((double)(det.im));
   if (d>(*d2))
      *d2=d;
}


int main(void)
{
   int i,n;
   float *rw,*rz,wsq,zsq;
   double *a,abar,sig,d1,d2,dmax1,dmax2;
   complex wuz;
   su3_vector w,z,uz;
   su3 u;

   printf("\n");
   printf("Statistical test of random_su3\n");
   printf("------------------------------\n\n");

   dmax1=0.0;
   dmax2=0.0;
   
   for (i=0;i<10000;i++)
   {
      random_su3(&u);
      dev(&u,&d1,&d2);

      if (d1>dmax1)
         dmax1=d1;
      if (d2>dmax2)
         dmax2=d2;
   }

   printf("In 10000 trials:\n");
   printf("max |1-U^dag*U| = %.1e\n",dmax1);
   printf("max |1-det U|   = %.1e\n",dmax2);
   
   for (;;)
   {
      rw=(float*)(&w);
      rz=(float*)(&z);
      gauss(rw,6);
      gauss(rz,6);
      wsq=_vector_prod_re(w,w);
      zsq=_vector_prod_re(z,z);      

      printf("\n");      
      printf("Specify number of trials (0 exits): ");

      if (scanf("%d",&n)==1)
      {
         printf("\n");

         if (n<=0)
            exit(0);

         a=amalloc(n*sizeof(double),3);

         for (i=0;i<n;i++)
         {
            random_su3(&u);
            _su3_multiply(uz,u,z);
            wuz.re=_vector_prod_re(w,uz);
            wuz.im=_vector_prod_im(w,uz);         

            a[i]=(double)(wuz.re*wuz.re+wuz.im*wuz.im);
         }
      
         abar=average(n,a);
         sig=sigma0(n,a);
   
         printf("<|wbar*U*z|^2> = %1.4f [error %.1e]\n",abar,sig);
         printf("Exact          = %1.4f\n",wsq*zsq/3.0f);
      
         afree(a);
      }
      else
      {
         printf("Invalid input, program stopped\n\n");
         break;
      }         
   }
   
   exit(0);
}
