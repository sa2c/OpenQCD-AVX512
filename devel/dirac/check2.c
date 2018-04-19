
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2005, 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Action of Dw() on plane waves.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "global.h"

static spinor rs ALIGNED16;
static const spinor sd0={{{0.0}}};


static su3_vector mul_cplx(complex z,su3_vector s)
{
   su3_vector r;

   r.c1.re=z.re*s.c1.re-z.im*s.c1.im;
   r.c1.im=z.im*s.c1.re+z.re*s.c1.im;
   r.c2.re=z.re*s.c2.re-z.im*s.c2.im;
   r.c2.im=z.im*s.c2.re+z.re*s.c2.im;
   r.c3.re=z.re*s.c3.re-z.im*s.c3.im;
   r.c3.im=z.im*s.c3.re+z.re*s.c3.im;

   return r;
}


static spinor mul_gamma(int mu,spinor s)
{
   spinor r;
   complex i,m_i,m_1;

   i.re=0.0f;
   i.im=1.0f;

   m_i.re=0.0f;
   m_i.im=-1.0f;

   m_1.re=-1.0f;
   m_1.im=0.0f;

   if (mu==0)
   {
      r.c1=mul_cplx(m_1,s.c3);
      r.c2=mul_cplx(m_1,s.c4);
      r.c3=mul_cplx(m_1,s.c1);
      r.c4=mul_cplx(m_1,s.c2);
   }
   else if (mu==1)
   {
      r.c1=mul_cplx(m_i,s.c4);
      r.c2=mul_cplx(m_i,s.c3);
      r.c3=mul_cplx(i,s.c2);
      r.c4=mul_cplx(i,s.c1);
   }
   else if (mu==2)
   {
      r.c1=mul_cplx(m_1,s.c4);
      r.c2=s.c3;
      r.c3=s.c2;
      r.c4=mul_cplx(m_1,s.c1);
   }
   else if (mu==3)
   {
      r.c1=mul_cplx(m_i,s.c3);
      r.c2=mul_cplx(i,s.c4);
      r.c3=mul_cplx(i,s.c1);
      r.c4=mul_cplx(m_i,s.c2);
   }
   else
   {
      r.c1=s.c1;
      r.c2=s.c2;
      r.c3=mul_cplx(m_1,s.c3);
      r.c4=mul_cplx(m_1,s.c4);
   }

   return r;
}


int main(int argc,char *argv[])
{
   int my_rank,bc;
   int n,i,ix,nu,x0,x1,x2,x3;
   int np[4],bo[4];
   float ran[4];
   float mu,pi,d,dmax;
   float mp,pt,pv,p[4],sp[4];
   double phi[2],phi_prime[2],theta[3];
   complex z;
   spinor **ps,s0,s1,s2,s3,s4;
   sw_parms_t swp;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      printf("\n");
      printf("Action of Dw() on plane waves\n");
      printf("-----------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      printf("For this test to pass, the calculated differences delta\n");
      printf("should be at most 1*10^(-5) or so\n\n");

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>]");
   }

   set_lat_parms(5.5,1.0,0,NULL,1.978);
   print_lat_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;
   set_bc_parms(bc,0.55,0.78,0.9012,1.2034,phi,phi_prime,theta);
   print_bc_parms(2);

   start_ranlux(0,12345);
   geometry();
   alloc_ws(3);
   ps=reserve_ws(3);

   swp=set_sw_parms(-0.0123);
   mu=0.0876f;

   if (my_rank==0)
      printf("m0 = %.4e, csw = %.4e, cF = %.4e, cF' = %.4e\n\n",
             swp.m0,swp.csw,swp.cF[0],swp.cF[1]);

   (void)udfld();
   set_ud_phase();
   sw_term(NO_PTS);
   assign_ud2u();
   assign_swd2sw();
   pi=(float)(4.0*atan(1.0));
   n=10;
   bo[0]=cpr[0]*L0;
   bo[1]=cpr[1]*L1;
   bo[2]=cpr[2]*L2;
   bo[3]=cpr[3]*L3;
   dmax=0.0f;

   for (i=0;i<n;i++)
   {
      ranlxs(ran,4);

      if (bc==0)
         np[0]=(int)(ran[0]*(float)(NPROC0*L0-1));
      else
         np[0]=(int)(ran[0]*(float)(NPROC0*L0));
      np[1]=(int)(ran[1]*(float)(NPROC1*L1));
      np[2]=(int)(ran[2]*(float)(NPROC2*L2));
      np[3]=(int)(ran[3]*(float)(NPROC3*L3));

      if (np[0]==0)
         np[0]=1;

      if (bc==0)
         p[0]=(float)(np[0])*pi/(float)(NPROC0*L0-1);
      else if (bc==3)
         p[0]=((float)(np[0])*2.0f*pi+pi)/(float)(NPROC0*L0);
      else
         p[0]=(float)(np[0])*pi/(float)(NPROC0*L0);
      p[1]=(float)(np[1])*2.0f*pi/(float)(NPROC1*L1);
      p[2]=(float)(np[2])*2.0f*pi/(float)(NPROC2*L2);
      p[3]=(float)(np[3])*2.0f*pi/(float)(NPROC3*L3);

      random_s(1,&rs,1.0f);

      MPI_Bcast(p,4,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(&rs,24,MPI_FLOAT,0,MPI_COMM_WORLD);

      sp[0]=(float)(sin((double)(p[0])));
      sp[1]=(float)(sin((double)(p[1])));
      sp[2]=(float)(sin((double)(p[2])));
      sp[3]=(float)(sin((double)(p[3])));

      mp=swp.m0;
      mp+=(float)(1.0-cos((double)(p[0])));
      mp+=(float)(1.0-cos((double)(p[1])));
      mp+=(float)(1.0-cos((double)(p[2])));
      mp+=(float)(1.0-cos((double)(p[3])));

      for (x0=0;x0<L0;x0++)
      {
         for (x1=0;x1<L1;x1++)
         {
            for (x2=0;x2<L2;x2++)
            {
               for (x3=0;x3<L3;x3++)
               {
                  ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

                  pt=p[0]*(float)(x0+bo[0]);
                  pv=p[1]*(float)(x1+bo[1])+p[2]*(float)(x2+bo[2])+
                     p[3]*(float)(x3+bo[3]);

                  if (bc!=3)
                  {
                     z.re=(float)(sin((double)(pt))*cos((double)(pv)));
                     z.im=(float)(sin((double)(pt))*sin((double)(pv)));
                  }
                  else
                  {
                     z.re=(float)(cos((double)(pt+pv)));
                     z.im=(float)(sin((double)(pt+pv)));
                  }

                  s0.c1=mul_cplx(z,rs.c1);
                  s0.c2=mul_cplx(z,rs.c2);
                  s0.c3=mul_cplx(z,rs.c3);
                  s0.c4=mul_cplx(z,rs.c4);

                  ps[0][ix]=s0;

                  z.re=cos(pt)*cos(pv);
                  z.im=cos(pt)*sin(pv);

                  s1.c1=mul_cplx(z,rs.c1);
                  s1.c2=mul_cplx(z,rs.c2);
                  s1.c3=mul_cplx(z,rs.c3);
                  s1.c4=mul_cplx(z,rs.c4);

                  z.re=mp;
                  z.im=0.0f;

                  if ((cpr[0]==0)&&(x0==1)&&(bc!=3))
                     z.re+=(float)(swp.cF[0]-1.0);

                  if ((cpr[0]==(NPROC0-1))&&(x0==(L0-2))&&(bc==0))
                     z.re+=(float)(swp.cF[1]-1.0);

                  if ((cpr[0]==(NPROC0-1))&&(x0==(L0-1))&&((bc==1)||(bc==2)))
                     z.re+=(float)(swp.cF[1]-1.0);

                  s2.c1=mul_cplx(z,s0.c1);
                  s2.c2=mul_cplx(z,s0.c2);
                  s2.c3=mul_cplx(z,s0.c3);
                  s2.c4=mul_cplx(z,s0.c4);

                  for (nu=0;nu<5;nu++)
                  {
                     if ((nu==0)&&(bc!=3))
                     {
                        s3=mul_gamma(0,s1);
                        z.re=sp[0];
                        z.im=0.0f;
                     }
                     else if (nu==4)
                     {
                        s3=mul_gamma(4,s0);
                        z.re=0.0f;
                        z.im=mu;
                     }
                     else
                     {
                        s3=mul_gamma(nu,s0);
                        z.re=0.0f;
                        z.im=sp[nu];
                     }

                     s4.c1=mul_cplx(z,s3.c1);
                     s4.c2=mul_cplx(z,s3.c2);
                     s4.c3=mul_cplx(z,s3.c3);
                     s4.c4=mul_cplx(z,s3.c4);

                     _vector_add_assign(s2.c1,s4.c1);
                     _vector_add_assign(s2.c2,s4.c2);
                     _vector_add_assign(s2.c3,s4.c3);
                     _vector_add_assign(s2.c4,s4.c4);
                  }

                  if (((cpr[0]==0)&&(x0==0)&&(bc!=3))||
                      ((cpr[0]==(NPROC0-1))&&(x0==(L0-1))&&(bc==0)))
                     ps[1][ix]=sd0;
                  else
                     ps[1][ix]=s2;
               }
            }
         }
      }

      Dw(mu,ps[0],ps[2]);

      mulr_spinor_add(VOLUME,ps[2],ps[1],-1.0f);
      d=norm_square(VOLUME,1,ps[2])/norm_square(VOLUME,1,ps[0]);
      d=(float)(sqrt((double)(d)));
      if (d>dmax)
         dmax=d;

      if (my_rank==0)
         printf("Normalized deviation = %.1e at p=(%d,%d,%d,%d)\n",
                d,np[0],np[1],np[2],np[3]);
   }

   if (my_rank==0)
   {
      printf("\n");
      printf("Maximal normalized deviation = %.1e\n\n",dmax);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
