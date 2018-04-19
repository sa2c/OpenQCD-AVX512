
/*******************************************************************************
*
* File valg.c
*
* Copyright (C) 2007, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic linear algebra routines for single-precision complex fields
*
* The externally accessible functions are
*
*   complex vprod(int n,int icom,complex *v,complex *w)
*     Computes the scalar product of the n-vectors v and w.
*
*   float vnorm_square(int n,int icom,complex *v)
*     Computes the square of the norm of the n-vector v.
*
*   void mulc_vadd(int n,complex *v,complex *w,complex z)
*     Replaces the n-vector v by v+z*w.
*
*   void vproject(int n,int icom,complex *v,complex *w)
*     Replaces the n-vector v by v-(w,v)*w.
*
*   void vscale(int n,float r,complex_dble *v)
*     Replaces the n-vector v by r*v.
*
*   float vnormalize(int n,int icom,complex *v)
*     Normalizes the n-vector v to unity and returns the norm of the
*     input vector.
*
*   void vrotate(int n,int nv,complex **pv,complex *a)
*     Replaces the n-vectors vk=pv[k], k=0,..,nv-1, by the linear
*     combinations sum_{j=0}^{nv-1} vj*a[n*j+k].
*
* Notes:
*
* All these programs operate on complex n-vectors whose base addresses are
* passed through the arguments. The length n of the arrays is specified by
* the parameter n. Scalar products are globally summed if the parameter
* icom is equal to 1. In this case the calculated values are guaranteed to
* be exactly the same on all processes.
*
* The programs perform no communications except in the case of the scalar
* products if these are globally summed.
*
*******************************************************************************/

#define VALG_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "linalg.h"
#include "vflds.h"
#include "global.h"

static int nrot=0;
static complex *psi;


static void alloc_wrotate(int n)
{
   if (nrot>0)
      afree(psi);

   psi=amalloc(n*sizeof(*psi),ALIGN);
   error_loc(psi==NULL,1,"alloc_wrotate [valg.c]",
             "Unable to allocate workspace");
   nrot=n;
   set_v2zero(n,psi);
}


complex vprod(int n,int icom,complex *v,complex *w)
{
   complex z,*vm;
   complex_dble vd,wd;

   vd.re=0.0;
   vd.im=0.0;
   vm=v+n;

   for (;v<vm;v++)
   {
         vd.re+=(double)((*v).re*(*w).re+(*v).im*(*w).im);
         vd.im+=(double)((*v).re*(*w).im-(*v).im*(*w).re);
         w+=1;
   }

   if ((icom!=1)||(NPROC==1))
   {
      z.re=(float)(vd.re);
      z.im=(float)(vd.im);
   }
   else
   {
      MPI_Reduce(&vd.re,&wd.re,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&wd.re,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

      z.re=(float)(wd.re);
      z.im=(float)(wd.im);
   }

   return z;
}


float vnorm_square(int n,int icom,complex *v)
{
   complex *vm;
   double x,y;

   x=0.0;
   vm=v+n;

   for (;v<vm;v++)
      x+=(double)((*v).re*(*v).re+(*v).im*(*v).im);

   if ((icom!=1)||(NPROC==1))
      return (float)(x);
   else
   {
      MPI_Reduce(&x,&y,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&y,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      return (float)(y);
   }
}


void mulc_vadd(int n,complex *v,complex *w,complex z)
{
   complex *vm;

   vm=v+n;

   for (;v<vm;v++)
   {
      (*v).re+=(z.re*(*w).re-z.im*(*w).im);
      (*v).im+=(z.re*(*w).im+z.im*(*w).re);
      w+=1;
   }
}


void vproject(int n,int icom,complex *v,complex *w)
{
   complex z;

   z=vprod(n,icom,w,v);
   z.re=-z.re;
   z.im=-z.im;
   mulc_vadd(n,v,w,z);
}


void vscale(int n,float r,complex *v)
{
   complex *vm;

   vm=v+n;

   for (;v<vm;v++)
   {
      (*v).re*=r;
      (*v).im*=r;
   }
}


float vnormalize(int n,int icom,complex *v)
{
   float r;

   r=vnorm_square(n,icom,v);
   r=(float)(sqrt((double)(r)));

   if (r!=0.0f)
      vscale(n,1.0f/r,v);
   else
      error_loc(1,1,"vnormalize [valg.c]",
                "Vector field has vanishing norm");

   return r;
}


void vrotate(int n,int nv,complex **pv,complex *a)
{
   int i,k,j;
   complex s,*z,*vj;

   if (nv>nrot)
      alloc_wrotate(nv);

   for (i=0;i<n;i++)
   {
      for (k=0;k<nv;k++)
      {
         s.re=0.0f;
         s.im=0.0f;
         z=a+k;

         for (j=0;j<nv;j++)
         {
            vj=pv[j]+i;
            s.re+=((*z).re*(*vj).re-(*z).im*(*vj).im);
            s.im+=((*z).re*(*vj).im+(*z).im*(*vj).re);
            z+=nv;
         }

         psi[k].re=s.re;
         psi[k].im=s.im;
      }

      for (k=0;k<nv;k++)
      {
         pv[k][i].re=psi[k].re;
         pv[k][i].im=psi[k].im;
      }
   }
}
