
/*******************************************************************************
*
* File check8.c
*
* Copyright (C) 2009, 2010, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the programs cm3x3_zero(),...,cm3x3_lc2()
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "random.h"
#include "su3fcts.h"

#define NTEST 1000
#define SEED 376

static double rs,*rus,*rvs,*rws,*rzs;
static complex_dble *cs;
static su3_dble *us,*vs,*ws,*zs;
static su3_dble ud0={{0.0,0.0},{0.0,0.0},{0.0,0.0},
                     {0.0,0.0},{0.0,0.0},{0.0,0.0},
                     {0.0,0.0},{0.0,0.0},{0.0,0.0}};
static complex_dble trs ALIGNED16;


static void alloc_matrices(void)
{
   rus=amalloc(90*sizeof(*rus),4);
   cs=amalloc(3*sizeof(*cs),4);
   us=amalloc(5*sizeof(*us),4);

   error((rus==NULL)||(cs==NULL)||(us==NULL),1,
         "alloc_matrices [check8.c]","Unable to allocate matrices");

   rvs=rus+36;
   rws=rvs+18;
   rzs=rws+18;

   vs=us+2;
   ws=vs+1;
   zs=ws+1;
}


static void mat2vec(su3_dble *u,double *ru)
{
   ru[ 0]=(*u).c11.re;
   ru[ 1]=(*u).c11.im;
   ru[ 2]=(*u).c12.re;
   ru[ 3]=(*u).c12.im;
   ru[ 4]=(*u).c13.re;
   ru[ 5]=(*u).c13.im;

   ru[ 6]=(*u).c21.re;
   ru[ 7]=(*u).c21.im;
   ru[ 8]=(*u).c22.re;
   ru[ 9]=(*u).c22.im;
   ru[10]=(*u).c23.re;
   ru[11]=(*u).c23.im;

   ru[12]=(*u).c31.re;
   ru[13]=(*u).c31.im;
   ru[14]=(*u).c32.re;
   ru[15]=(*u).c32.im;
   ru[16]=(*u).c33.re;
   ru[17]=(*u).c33.im;
}


static void vec2mat(double *ru,su3_dble *u)
{
   (*u).c11.re=ru[ 0];
   (*u).c11.im=ru[ 1];
   (*u).c12.re=ru[ 2];
   (*u).c12.im=ru[ 3];
   (*u).c13.re=ru[ 4];
   (*u).c13.im=ru[ 5];

   (*u).c21.re=ru[ 6];
   (*u).c21.im=ru[ 7];
   (*u).c22.re=ru[ 8];
   (*u).c22.im=ru[ 9];
   (*u).c23.re=ru[10];
   (*u).c23.im=ru[11];

   (*u).c31.re=ru[12];
   (*u).c31.im=ru[13];
   (*u).c32.re=ru[14];
   (*u).c32.im=ru[15];
   (*u).c33.re=ru[16];
   (*u).c33.im=ru[17];
}


static void add_vec(double *ru,double *rv,double *rw)
{
   int i;

   for (i=0;i<18;i++)
      rw[i]=ru[i]+rv[i];
}


static void mulr_vec(double r,double *ru,double *rv)
{
   int i;

   for (i=0;i<18;i++)
      rv[i]=r*ru[i];
}


static void mulc_vec(complex_dble c,double *ru,double *rv)
{
   int i;

   for (i=0;i<18;i+=2)
   {
      rv[i  ]=c.re*ru[i  ]-c.im*ru[i+1];
      rv[i+1]=c.re*ru[i+1]+c.im*ru[i  ];
   }
}


static void dag_vec(double *ru,double *rv)
{
   int i,j;

   for (i=0;i<3;i++)
   {
      for (j=0;j<3;j++)
      {
         rv[6*i+2*j  ]= ru[6*j+2*i  ];
         rv[6*i+2*j+1]=-ru[6*j+2*i+1];
      }
   }
}


static void random_matrix(su3_dble *u)
{
   int i;
   double r[18];

   ranlxd(r,18);

   for (i=0;i<18;i++)
      r[i]=2.0*r[i]-1.0;

   vec2mat(r,u);
}


static void start_test(void)
{
   int i;
   double r[6];

   ranlxd(&rs,1);
   rs=2.0*rs-1.0;

   ranlxd(r,6);

   for (i=0;i<6;i++)
      r[i]=2.0*r[i]-1.0;

   cs[0].re=r[0];
   cs[0].im=r[1];
   cs[1].re=r[2];
   cs[1].im=r[3];
   cs[2].re=r[4];
   cs[2].im=r[5];

   random_matrix(us);
   random_matrix(us+1);
   random_matrix(vs);
   random_matrix(ws);
   random_matrix(zs);
}


static double dev_uv(su3_dble *u,su3_dble *v)
{
   int i;
   double r[18],dev,dmax;

   r[ 0]=(*u).c11.re-(*v).c11.re;
   r[ 1]=(*u).c11.im-(*v).c11.im;
   r[ 2]=(*u).c12.re-(*v).c12.re;
   r[ 3]=(*u).c12.im-(*v).c12.im;
   r[ 4]=(*u).c13.re-(*v).c13.re;
   r[ 5]=(*u).c13.im-(*v).c13.im;

   r[ 6]=(*u).c21.re-(*v).c21.re;
   r[ 7]=(*u).c21.im-(*v).c21.im;
   r[ 8]=(*u).c22.re-(*v).c22.re;
   r[ 9]=(*u).c22.im-(*v).c22.im;
   r[10]=(*u).c23.re-(*v).c23.re;
   r[11]=(*u).c23.im-(*v).c23.im;

   r[12]=(*u).c31.re-(*v).c31.re;
   r[13]=(*u).c31.im-(*v).c31.im;
   r[14]=(*u).c32.re-(*v).c32.re;
   r[15]=(*u).c32.im-(*v).c32.im;
   r[16]=(*u).c33.re-(*v).c33.re;
   r[17]=(*u).c33.im-(*v).c33.im;

   dmax=0.0;

   for (i=0;i<18;i++)
   {
      dev=fabs(r[i]);
      if (dev>dmax)
         dmax=dev;
   }

   return dmax;
}


int main(void)
{
   int i;
   double dev,dmax[15];

   printf("\n");
   printf("Check of the programs cm3x3_zero(),...,cm3x3_lc2()\n");
   printf("--------------------------------------------------\n\n");

   printf("Test performed on %d random matrices\n\n",NTEST);

   rlxd_init(1,SEED);
   alloc_matrices();

   for (i=0;i<15;i++)
      dmax[i]=0.0;

   for (i=0;i<NTEST;i++)
   {
      start_test();

      mat2vec(vs,rvs);
      cm3x3_zero(1,vs);
      dev=dev_uv(&ud0,vs);
      vec2mat(rvs,vs);
      if (dev>dmax[0])
         dmax[0]=dev;

      cm3x3_unity(1,vs);
      (*vs).c11.re-=1.0;
      (*vs).c22.re-=1.0;
      (*vs).c33.re-=1.0;
      dev=dev_uv(&ud0,vs);
      vec2mat(rvs,vs);
      if (dev>dmax[1])
         dmax[1]=dev;

      mat2vec(us,rus);
      mat2vec(ws,rws);
      vec2mat(rus,ws);
      cm3x3_assign(1,us,vs);
      dev=dev_uv(ws,vs);
      if (dev>dmax[2])
         dmax[2]=dev;
      dev=dev_uv(ws,us);
      if (dev>dmax[2])
         dmax[2]=dev;
      vec2mat(rvs,vs);
      vec2mat(rws,ws);

      mat2vec(us,rus);
      mat2vec(vs,rvs);
      mat2vec(ws,rws);
      cm3x3_swap(1,us,vs);
      vec2mat(rus,ws);
      dev=dev_uv(ws,vs);
      if (dev>dmax[14])
         dmax[14]=dev;
      vec2mat(rvs,ws);
      dev=dev_uv(ws,us);
      if (dev>dmax[14])
         dmax[14]=dev;
      vec2mat(rus,us);
      vec2mat(rvs,vs);
      vec2mat(rws,ws);

      mat2vec(us,rus);
      dag_vec(rus,rzs);
      vec2mat(rzs,ws);
      cm3x3_dagger(us,vs);
      dev=dev_uv(vs,ws);
      if (dev>dmax[3])
         dmax[3]=dev;
      vec2mat(rus,ws);
      dev=dev_uv(us,ws);
      if (dev>dmax[3])
         dmax[3]=dev;

      vec2mat(rzs,ws);
      cm3x3_dagger(us,us);
      dev=dev_uv(us,ws);
      if (dev>dmax[3])
      dmax[3]=dev;
      vec2mat(rus,us);
      vec2mat(rvs,vs);
      vec2mat(rws,ws);

      cm3x3_tr(us,us+1,&trs);
      su3xsu3(us,us+1,vs);
      trs.re-=((*vs).c11.re+(*vs).c22.re+(*vs).c33.re);
      trs.im-=((*vs).c11.im+(*vs).c22.im+(*vs).c33.im);
      dev=fabs(trs.re)+fabs(trs.im);
      vec2mat(rvs,vs);
      if (dev>dmax[4])
         dmax[4]=dev;

      cm3x3_tr(us,us,&trs);
      su3xsu3(us,us,vs);
      trs.re-=((*vs).c11.re+(*vs).c22.re+(*vs).c33.re);
      trs.im-=((*vs).c11.im+(*vs).c22.im+(*vs).c33.im);
      dev=fabs(trs.re)+fabs(trs.im);
      vec2mat(rvs,vs);
      if (dev>dmax[4])
         dmax[4]=dev;

      cm3x3_retr(us,us+1,&trs.re);
      su3xsu3(us,us+1,vs);
      trs.re-=((*vs).c11.re+(*vs).c22.re+(*vs).c33.re);
      dev=fabs(trs.re);
      vec2mat(rvs,vs);
      if (dev>dmax[5])
         dmax[5]=dev;

      cm3x3_retr(us,us,&trs.re);
      su3xsu3(us,us,vs);
      trs.re-=((*vs).c11.re+(*vs).c22.re+(*vs).c33.re);
      dev=fabs(trs.re);
      vec2mat(rvs,vs);
      if (dev>dmax[5])
         dmax[5]=dev;

      mat2vec(us,rus);
      mat2vec(vs,rvs);
      mat2vec(zs,rzs);
      add_vec(rus,rvs,rws);
      vec2mat(rws,zs);
      cm3x3_add(us,vs);
      dev=dev_uv(vs,zs);
      vec2mat(rvs,vs);
      vec2mat(rzs,zs);
      if (dev>dmax[6])
         dmax[6]=dev;

      mat2vec(us,rvs);
      add_vec(rus,rvs,rws);
      vec2mat(rws,zs);
      cm3x3_add(us,us);
      dev=dev_uv(us,zs);
      vec2mat(rus,us);
      vec2mat(rzs,zs);
      if (dev>dmax[6])
         dmax[6]=dev;

      mat2vec(ws,rws);
      mat2vec(zs,rzs);
      su3xsu3(us,vs,zs);
      mat2vec(zs,rvs);
      add_vec(rvs,rws,rus);
      vec2mat(rus,zs);
      cm3x3_mul_add(us,vs,ws);
      vec2mat(rws,ws);
      cm3x3_mul_add(us,vs,ws);
      dev=dev_uv(ws,zs);
      vec2mat(rws,ws);
      vec2mat(rzs,zs);
      if (dev>dmax[7])
         dmax[7]=dev;

      mat2vec(vs,rvs);
      mat2vec(zs,rzs);
      su3xsu3(us,vs,zs);
      mat2vec(zs,rws);
      add_vec(rws,rvs,rus);
      vec2mat(rus,zs);
      cm3x3_mul_add(us,vs,vs);
      vec2mat(rvs,vs);
      cm3x3_mul_add(us,vs,vs);
      dev=dev_uv(vs,zs);
      vec2mat(rvs,vs);
      vec2mat(rzs,zs);
      if (dev>dmax[7])
         dmax[7]=dev;

      mat2vec(ws,rws);
      mat2vec(zs,rzs);
      mat2vec(us,rus);
      mulr_vec(rs,rus,rvs);
      vec2mat(rvs,zs);
      cm3x3_mulr(&rs,us,ws);
      dev=dev_uv(ws,zs);
      if (dev>dmax[8])
         dmax[8]=dev;

      cm3x3_mulr(&rs,us,us);
      dev=dev_uv(us,zs);
      vec2mat(rus,us);
      vec2mat(rws,ws);
      vec2mat(rzs,zs);
      if (dev>dmax[8])
         dmax[8]=dev;

      mat2vec(us,rus);
      mat2vec(vs,rvs);
      mulr_vec(rs,rus,rws);
      add_vec(rws,rvs,rzs);
      mat2vec(ws,rws);
      vec2mat(rzs,ws);
      cm3x3_mulr_add(&rs,us,vs);
      dev=dev_uv(vs,ws);
      vec2mat(rvs,vs);
      vec2mat(rws,ws);
      if (dev>dmax[9])
         dmax[9]=dev;

      mulr_vec(rs,rus,rws);
      add_vec(rus,rws,rvs);
      mat2vec(ws,rws);
      vec2mat(rvs,ws);
      cm3x3_mulr_add(&rs,us,us);
      dev=dev_uv(us,ws);
      vec2mat(rus,us);
      vec2mat(rws,ws);
      if (dev>dmax[9])
         dmax[9]=dev;

      mat2vec(ws,rws);
      mat2vec(zs,rzs);
      mat2vec(us,rus);
      mulc_vec(cs[0],rus,rvs);
      vec2mat(rvs,zs);
      cm3x3_mulc(cs,us,ws);
      dev=dev_uv(ws,zs);
      if (dev>dmax[10])
         dmax[10]=dev;

      cm3x3_mulc(cs,us,us);
      dev=dev_uv(us,zs);
      vec2mat(rus,us);
      vec2mat(rws,ws);
      vec2mat(rzs,zs);
      if (dev>dmax[10])
         dmax[10]=dev;

      mat2vec(us,rus);
      mat2vec(vs,rvs);
      mulc_vec(cs[0],rus,rws);
      add_vec(rws,rvs,rzs);
      mat2vec(ws,rws);
      vec2mat(rzs,ws);
      cm3x3_mulc_add(cs,us,vs);
      dev=dev_uv(vs,ws);
      vec2mat(rvs,vs);
      vec2mat(rws,ws);
      if (dev>dmax[11])
         dmax[11]=dev;

      mulc_vec(cs[0],rus,rws);
      add_vec(rus,rws,rvs);
      mat2vec(ws,rws);
      vec2mat(rvs,ws);
      cm3x3_mulc_add(cs,us,us);
      dev=dev_uv(us,ws);
      vec2mat(rus,us);
      vec2mat(rws,ws);
      if (dev>dmax[11])
         dmax[11]=dev;

      mat2vec(us,rus);
      mat2vec(vs,rvs);
      mat2vec(ws,rws);
      mulc_vec(cs[1],rus,rzs);
      vec2mat(rzs,ws);
      (*ws).c11.re+=cs[0].re;
      (*ws).c11.im+=cs[0].im;
      (*ws).c22.re+=cs[0].re;
      (*ws).c22.im+=cs[0].im;
      (*ws).c33.re+=cs[0].re;
      (*ws).c33.im+=cs[0].im;
      cm3x3_lc1(cs,us,vs);
      dev=dev_uv(vs,ws);
      if (dev>dmax[12])
         dmax[12]=dev;

      cm3x3_lc1(cs,us,us);
      dev=dev_uv(us,ws);
      vec2mat(rus,us);
      vec2mat(rvs,vs);
      vec2mat(rws,ws);
      if (dev>dmax[12])
         dmax[12]=dev;

      mat2vec(us,rus);
      mat2vec(us+1,rus+18);
      mulc_vec(cs[1],rus,rvs);
      mulc_vec(cs[2],rus+18,rws);
      add_vec(rvs,rws,rzs);
      vec2mat(rzs,ws);
      (*ws).c11.re+=cs[0].re;
      (*ws).c11.im+=cs[0].im;
      (*ws).c22.re+=cs[0].re;
      (*ws).c22.im+=cs[0].im;
      (*ws).c33.re+=cs[0].re;
      (*ws).c33.im+=cs[0].im;
      cm3x3_lc2(cs,us,vs);
      dev=dev_uv(vs,ws);
      if (dev>dmax[13])
         dmax[13]=dev;

      cm3x3_lc2(cs,us,us);
      vec2mat(rus,us);
      cm3x3_lc2(cs,us,us);
      dev=dev_uv(us,ws);
      if (dev>dmax[13])
         dmax[13]=dev;
   }

   printf("Maximal deviation of cm3x3_zero()     = %.1e\n",dmax[0]);
   printf("Maximal deviation of cm3x3_unity()    = %.1e\n",dmax[1]);
   printf("Maximal deviation of cm3x3_assign()   = %.1e\n",dmax[2]);
   printf("Maximal deviation of cm3x3_swap()     = %.1e\n",dmax[14]);
   printf("Maximal deviation of cm3x3_dagger()   = %.1e\n",dmax[3]);
   printf("Maximal deviation of cm3x3_tr()       = %.1e\n",dmax[4]);
   printf("Maximal deviation of cm3x3_retr()     = %.1e\n",dmax[5]);
   printf("Maximal deviation of cm3x3_add()      = %.1e\n",dmax[6]);
   printf("Maximal deviation of cm3x3_mul_add()  = %.1e\n",dmax[7]);
   printf("Maximal deviation of cm3x3_mulr()     = %.1e\n",dmax[8]);
   printf("Maximal deviation of cm3x3_mulr_add() = %.1e\n",dmax[9]);
   printf("Maximal deviation of cm3x3_mulc()     = %.1e\n",dmax[10]);
   printf("Maximal deviation of cm3x3_mulc_add() = %.1e\n",dmax[11]);
   printf("Maximal deviation of cm3x3_lc1()      = %.1e\n",dmax[12]);
   printf("Maximal deviation of cm3x3_lc2()      = %.1e\n\n",dmax[13]);

   for (i=1;i<15;i++)
   {
      if (dmax[i]>dmax[0])
         dmax[0]=dmax[i];
   }

   printf("Maximal deviation (all tests)         = %.1e\n\n",dmax[0]);
   exit(0);
}
