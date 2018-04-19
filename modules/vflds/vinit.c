
/*******************************************************************************
*
* File vinit.c
*
* Copyright (C) 2007, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Generic initialization and assignment programs for complex single- and
* double-precision fields
*
* The externally accessible functions are
*
*   void set_v2zero(int n,complex *v)
*     Sets the single-precision field v to zero.
*
*   void set_vd2zero(int n,complex_dble *vd)
*     Sets the double-precision field vd to zero.
*
*   void random_v(int n,complex *v,float sigma)
*     Initializes the components of the single-precision field v to
*     (complex) random values z with distribution proportional to
*     exp{-|z|^2/sigma^2}.
*
*   void random_vd(int n,complex_dble *vd,double sigma)
*     Initializes the components of the double-precision field vd to
*     (complex) random values z with distribution proportional to
*     exp{-|z|^2/sigma^2}.
*
*   void assign_v2v(int n,complex *v,complex *w)
*     Assigns the single-precision field v to the single-precision
*     field w.
*
*   void assign_v2vd(int n,complex *v,complex_dble *wd)
*     Assigns the single-precision field v to the double-precision
*     field wd.
*
*   void assign_vd2v(int n,complex_dble *vd,complex *w)
*     Assigns the double-precision field vd to the single-precision
*     field w.
*
*   void assign_vd2vd(int n,complex_dble *vd,complex_dble *wd)
*     Assigns the double-precision field vd to the double-precision
*     field wd.
*
*   void add_v2vd(int n,complex *v,complex_dble *wd)
*     Adds the single-precision field v to the double-precision field
*     wd.
*
*   void diff_vd2v(int n,complex_dble *vd,complex_dble *wd,complex *w)
*     Assigns the difference vd-wd of the double-precision fields vd
*     and wd to the single-precision field w.
*
* Notes:
*
* All these programs operate on arrays of complex fields, whose base
* addresses are passed through the arguments. The length of the arrays
* is specified by the parameter n.
*
* Since no communications are performed, all programs in this file can be
* called locally.
*
*******************************************************************************/

#define VINIT_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "random.h"
#include "vflds.h"

static const complex v0={0.0f,0.0f};
static const complex_dble vd0={0.0,0.0};


void set_v2zero(int n,complex *v)
{
   complex *vm;

   vm=v+n;

   for (;v<vm;v++)
      (*v)=v0;
}


void set_vd2zero(int n,complex_dble *vd)
{
   complex_dble *vm;

   vm=vd+n;

   for (;vd<vm;vd++)
      (*vd)=vd0;
}


void random_v(int n,complex *v,float sigma)
{
   int m;
   float r[8];
   complex *vm;

   m=n&0x3;
   vm=v+(n-m);

   while (v<vm)
   {
      gauss(r,8);
      v[0].re=sigma*r[0];
      v[0].im=sigma*r[1];
      v[1].re=sigma*r[2];
      v[1].im=sigma*r[3];
      v[2].re=sigma*r[4];
      v[2].im=sigma*r[5];
      v[3].re=sigma*r[6];
      v[3].im=sigma*r[7];
      v+=4;
   }

   vm=v+m;

   while (v<vm)
   {
      gauss(r,2);
      v[0].re=sigma*r[0];
      v[0].im=sigma*r[1];
      v+=1;
   }
}


void random_vd(int n,complex_dble *vd,double sigma)
{
   int m;
   double rd[8];
   complex_dble *vm;

   m=n&0x3;
   vm=vd+(n-m);

   while (vd<vm)
   {
      gauss_dble(rd,8);
      vd[0].re=sigma*rd[0];
      vd[0].im=sigma*rd[1];
      vd[1].re=sigma*rd[2];
      vd[1].im=sigma*rd[3];
      vd[2].re=sigma*rd[4];
      vd[2].im=sigma*rd[5];
      vd[3].re=sigma*rd[6];
      vd[3].im=sigma*rd[7];
      vd+=4;
   }

   vm=vd+m;

   while (vd<vm)
   {
      gauss_dble(rd,2);
      vd[0].re=sigma*rd[0];
      vd[0].im=sigma*rd[1];
      vd+=1;
   }
}


void assign_v2v(int n,complex *v,complex *w)
{
   complex *vm;

   vm=v+n;

   for (;v<vm;v++)
   {
      (*w).re=(*v).re;
      (*w).im=(*v).im;
      w+=1;
   }
}


void assign_v2vd(int n,complex *v,complex_dble *wd)
{
   complex *vm;

   vm=v+n;

   for (;v<vm;v++)
   {
      (*wd).re=(double)((*v).re);
      (*wd).im=(double)((*v).im);
      wd+=1;
   }
}


void assign_vd2v(int n,complex_dble *vd,complex *w)
{
   complex_dble *vm;

   vm=vd+n;

   for (;vd<vm;vd++)
   {
      (*w).re=(float)((*vd).re);
      (*w).im=(float)((*vd).im);
      w+=1;
   }
}


void assign_vd2vd(int n,complex_dble *vd,complex_dble *wd)
{
   complex_dble *vm;

   vm=vd+n;

   for (;vd<vm;vd++)
   {
      (*wd).re=(*vd).re;
      (*wd).im=(*vd).im;
      wd+=1;
   }
}


void add_v2vd(int n,complex *v,complex_dble *wd)
{
   complex *vm;

   vm=v+n;

   for (;v<vm;v++)
   {
      (*wd).re+=(double)((*v).re);
      (*wd).im+=(double)((*v).im);
      wd+=1;
   }
}


void diff_vd2v(int n,complex_dble *vd,complex_dble *wd,complex *w)
{
   complex_dble *vm;

   vm=vd+n;

   for (;vd<vm;vd++)
   {
      (*w).re=(float)((*vd).re-(*wd).re);
      (*w).im=(float)((*vd).im-(*wd).im);
      w+=1;
      wd+=1;
   }
}
