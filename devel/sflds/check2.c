
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the communication programs in scom.c and sdcom.c.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "random.h"
#include "su3fcts.h"
#include "utils.h"
#include "sflds.h"
#include "linalg.h"
#include "lattice.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)
#define NFLDS 4

typedef union
{
   spinor_dble s;
   double r[24];
} spin_dble_t;

static double p[4];
static spinor_dble rs ALIGNED16;
static const spinor_dble sd0={{{0.0}}};


static int is_zero_dble(spinor_dble *s)
{
   int i,ie;
   spin_dble_t *sp;

   sp=(spin_dble_t*)(s);
   ie=1;

   for (i=0;i<24;i++)
      ie&=((*sp).r[i]==0.0);

   return ie;
}


static int check_int_bnd_dble(spinor_dble *s)
{
   int bc,ix,iy,t;
   int ie;

   bc=bc_type();
   ie=1;

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((ix<(VOLUME/2))&&
          (((t==0)&&(bc!=3))||((t==(N0-1))&&(bc==0))))
         ie&=is_zero_dble(s);
      else if ((ix>=(VOLUME/2))&&(t==(N0-1))&&((bc==1)||(bc==2)))
      {
         iy=iup[ix][0];
         ie&=is_zero_dble(s+iy-ix);
      }
      else
         ie&=(is_zero_dble(s)^0x1);

      s+=1;
   }

   return ie;
}


static int check_ext_bnd_dble(spinor_dble *s)
{
   int bc,ix,t;
   int ie;

   bc=bc_type();
   ie=1;

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((ix<(VOLUME/2))&&
          (((t==0)&&(bc!=3))||((t==(N0-1))&&(bc==0))))
         ie&=is_zero_dble(s);
      else
         ie&=(is_zero_dble(s)^0x1);

      s+=1;
   }

   return ie;
}


static su3_vector_dble mul_cplx(complex_dble z,su3_vector_dble *s)
{
   su3_vector_dble r;

   r.c1.re=z.re*(*s).c1.re-z.im*(*s).c1.im;
   r.c1.im=z.im*(*s).c1.re+z.re*(*s).c1.im;
   r.c2.re=z.re*(*s).c2.re-z.im*(*s).c2.im;
   r.c2.im=z.im*(*s).c2.re+z.re*(*s).c2.im;
   r.c3.re=z.re*(*s).c3.re-z.im*(*s).c3.im;
   r.c3.im=z.im*(*s).c3.re+z.re*(*s).c3.im;

   return r;
}


static spinor_dble mul_gamma(int mu,spinor_dble *s)
{
   spinor_dble r;
   complex_dble i,m_i,m_1;

   i.re=0.0;
   i.im=1.0;

   m_i.re=0.0;
   m_i.im=-1.0;

   m_1.re=-1.0;
   m_1.im=0.0;

   if (mu==0)
   {
      r.c1=mul_cplx(m_1,&((*s).c3));
      r.c2=mul_cplx(m_1,&((*s).c4));
      r.c3=mul_cplx(m_1,&((*s).c1));
      r.c4=mul_cplx(m_1,&((*s).c2));
   }
   else if (mu==1)
   {
      r.c1=mul_cplx(m_i,&((*s).c4));
      r.c2=mul_cplx(m_i,&((*s).c3));
      r.c3=mul_cplx(i,&((*s).c2));
      r.c4=mul_cplx(i,&((*s).c1));
   }
   else if (mu==2)
   {
      r.c1=mul_cplx(m_1,&((*s).c4));
      r.c2=(*s).c3;
      r.c3=(*s).c2;
      r.c4=mul_cplx(m_1,&((*s).c1));
   }
   else if (mu==3)
   {
      r.c1=mul_cplx(m_i,&((*s).c3));
      r.c2=mul_cplx(i,&((*s).c4));
      r.c3=mul_cplx(i,&((*s).c1));
      r.c4=mul_cplx(m_i,&((*s).c2));
   }
   else
   {
      r.c1=(*s).c1;
      r.c2=(*s).c2;
      r.c3=mul_cplx(m_1,&((*s).c3));
      r.c4=mul_cplx(m_1,&((*s).c4));
   }

   return r;
}


static spinor_dble theta(int ifc,spinor_dble *s)
{
   int i;
   spin_dble_t r,*sp;

   sp=(spin_dble_t*)(s);
   r.s=mul_gamma(ifc/2,s);

   if (ifc&0x1)
   {
      for (i=0;i<18;i++)
         r.r[i]=0.5*((*sp).r[i]-r.r[i]);
   }
   else
   {
      for (i=0;i<18;i++)
         r.r[i]=0.5*((*sp).r[i]+r.r[i]);
   }

   return r.s;
}


static void set_sd(spinor_dble *s)
{
   int bc,bo[4],np[4],k;
   int ix,x0,x1,x2,x3;
   float ran[4];
   double pi,pt,pv;
   complex_dble z;

   bc=bc_type();
   bo[0]=cpr[0]*L0;
   bo[1]=cpr[1]*L1;
   bo[2]=cpr[2]*L2;
   bo[3]=cpr[3]*L3;

   random_sd(1,&rs,1.0);
   MPI_Bcast(&rs,24,MPI_DOUBLE,0,MPI_COMM_WORLD);
   ranlxs(ran,4);
   MPI_Bcast(ran,4,MPI_FLOAT,0,MPI_COMM_WORLD);

   if (bc==0)
      np[0]=(int)(ran[0]*(float)(N0-1));
   else
      np[0]=(int)(ran[0]*(float)(N0));
   np[1]=(int)(ran[1]*(float)(N1));
   np[2]=(int)(ran[2]*(float)(N2));
   np[3]=(int)(ran[3]*(float)(N3));

   for (k=0;k<4;k++)
      if (np[k]==0)
         np[k]=1;

   pi=4.0*atan(1.0);

   if (bc==0)
      p[0]=pi*(double)(np[0])/(double)(N0-1);
   else if (bc==3)
      p[0]=2*pi*(double)(np[0])/(double)(N0);
   else
      p[0]=pi*(double)(np[0])/(double)(N0);

   p[1]=2*pi*(double)(np[1])/(double)(N1);
   p[2]=2*pi*(double)(np[2])/(double)(N2);
   p[3]=2*pi*(double)(np[3])/(double)(N3);

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

               pt=p[0]*(double)(x0+bo[0]);
               pv=p[1]*(double)(x1+bo[1])+
                  p[2]*(double)(x2+bo[2])+
                  p[3]*(double)(x3+bo[3]);

               if (bc==3)
               {
                  z.re=cos(pt+pv);
                  z.im=sin(pt+pv);
               }
               else
               {
                  z.re=sin(pt)*cos(pv);
                  z.im=sin(pt)*sin(pv);
               }

               s[ix].c1=mul_cplx(z,&(rs.c1));
               s[ix].c2=mul_cplx(z,&(rs.c2));
               s[ix].c3=mul_cplx(z,&(rs.c3));
               s[ix].c4=mul_cplx(z,&(rs.c4));
            }
         }
      }
   }

   bnd_sd2zero(ALL_PTS,s);
}


static void set_sd_bnd(int is,spinor_dble *s)
{
   int bc,bo[4],np[4],k;
   int x0,x1,x2,x3,x[4];
   int ix,iy,ifc,mu;
   float ran[4];
   double pi,pt,pv;
   complex_dble z;

   bc=bc_type();
   bo[0]=cpr[0]*L0;
   bo[1]=cpr[1]*L1;
   bo[2]=cpr[2]*L2;
   bo[3]=cpr[3]*L3;

   random_sd(1,&rs,1.0);
   MPI_Bcast(&rs,24,MPI_DOUBLE,0,MPI_COMM_WORLD);
   ranlxs(ran,4);
   MPI_Bcast(ran,4,MPI_FLOAT,0,MPI_COMM_WORLD);

   if (bc==0)
      np[0]=(int)(ran[0]*(float)(N0-1));
   else
      np[0]=(int)(ran[0]*(float)(N0));
   np[1]=(int)(ran[1]*(float)(N1));
   np[2]=(int)(ran[2]*(float)(N2));
   np[3]=(int)(ran[3]*(float)(N3));

   for (k=0;k<4;k++)
      if (np[k]==0)
         np[k]=1;

   pi=4.0*atan(1.0);

   if (bc==0)
      p[0]=pi*(double)(np[0])/(double)(N0-1);
   else if (bc==3)
      p[0]=2*pi*(double)(np[0])/(double)(N0);
   else
      p[0]=pi*(double)(np[0])/(double)(N0);

   p[1]=2*pi*(double)(np[1])/(double)(N1);
   p[2]=2*pi*(double)(np[2])/(double)(N2);
   p[3]=2*pi*(double)(np[3])/(double)(N3);

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               if ((x0+x1+x2+x3)&0x1)
               {
                  ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

                  for (ifc=0;ifc<8;ifc++)
                  {
                     mu=ifc/2;
                     x[0]=x0;
                     x[1]=x1;
                     x[2]=x2;
                     x[3]=x3;

                     if (ifc&0x1)
                     {
                        iy=iup[ix][mu];
                        x[mu]+=1;
                     }
                     else
                     {
                        iy=idn[ix][mu];
                        x[mu]-=1;
                     }

                     if ((iy>=VOLUME)&&
                         ((ifc>1)||
                          ((ifc==0)&&((cpr[0]>0)||(bc==3)))||
                          ((ifc==1)&&((cpr[0]<(NPROC0-1))||(bc==3)))))
                     {
                        pt=p[0]*(double)(x[0]+bo[0]);
                        pv=p[1]*(double)(x[1]+bo[1])+
                           p[2]*(double)(x[2]+bo[2])+
                           p[3]*(double)(x[3]+bo[3]);

                        if (bc==3)
                        {
                           z.re=cos(pt+pv);
                           z.im=sin(pt+pv);
                        }
                        else
                        {
                           z.re=sin(pt)*cos(pv);
                           z.im=sin(pt)*sin(pv);
                        }

                        s[iy].c1=mul_cplx(z,&(rs.c1));
                        s[iy].c2=mul_cplx(z,&(rs.c2));
                        s[iy].c3=mul_cplx(z,&(rs.c3));
                        s[iy].c4=mul_cplx(z,&(rs.c4));
                        s[iy]=theta(ifc^is,s+iy);
                     }
                  }
               }
            }
         }
      }
   }

   bnd_sd2zero(EVEN_PTS,s);
}


static double check_cpsd_int(int is,spinor_dble *s)
{
   int bc,bo[4];
   int x0,x1,x2,x3,x[4];
   int ix,iy,ifc,mu,i;
   double pt,pv,d,dmax;
   complex_dble z;
   spin_dble_t r,*sp;

   bc=bc_type();
   bo[0]=cpr[0]*L0;
   bo[1]=cpr[1]*L1;
   bo[2]=cpr[2]*L2;
   bo[3]=cpr[3]*L3;
   dmax=0.0;

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               if ((x0+x1+x2+x3)&0x1)
               {
                  ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

                  for (ifc=0;ifc<8;ifc++)
                  {
                     mu=ifc/2;
                     x[0]=x0;
                     x[1]=x1;
                     x[2]=x2;
                     x[3]=x3;

                     if (ifc&0x1)
                     {
                        iy=iup[ix][mu];
                        x[mu]+=1;
                     }
                     else
                     {
                        iy=idn[ix][mu];
                        x[mu]-=1;
                     }

                     if ((iy>=VOLUME)&&
                         ((ifc>1)||
                          ((ifc==0)&&((cpr[0]>0)||(bc==3)))||
                          ((ifc==1)&&((cpr[0]<(NPROC0-1))||(bc==3)))))
                     {
                        pt=p[0]*(double)(x[0]+bo[0]);
                        pv=p[1]*(double)(x[1]+bo[1])+
                           p[2]*(double)(x[2]+bo[2])+
                           p[3]*(double)(x[3]+bo[3]);

                        if (bc==3)
                        {
                           z.re=cos(pt+pv);
                           z.im=sin(pt+pv);
                        }
                        else
                        {
                           z.re=sin(pt)*cos(pv);
                           z.im=sin(pt)*sin(pv);
                        }

                        r.s.c1=mul_cplx(z,&(rs.c1));
                        r.s.c2=mul_cplx(z,&(rs.c2));
                        r.s.c3=mul_cplx(z,&(rs.c3));
                        r.s.c4=mul_cplx(z,&(rs.c4));
                        sp=(spin_dble_t*)(s+iy);

                        for (i=0;i<18;i++)
                           r.r[i]-=(*sp).r[i];

                        r.s=theta((ifc^0x1)^is,&(r.s));

                        for (i=0;i<18;i++)
                        {
                           d=fabs(r.r[i]);
                           if (d>dmax)
                              dmax=d;
                        }
                     }
                  }
               }
            }
         }
      }
   }

   if (NPROC>1)
   {
      d=dmax;
      MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return dmax;
}


static double check_cpsd_ext(int is,spinor_dble *s)
{
   int bc,bo[4];
   int x0,x1,x2,x3;
   int ix,iy,ifc,mu,i;
   double pt,pv,d,dmax;
   complex_dble z;
   spin_dble_t r,*sp;

   bc=bc_type();
   bo[0]=cpr[0]*L0;
   bo[1]=cpr[1]*L1;
   bo[2]=cpr[2]*L2;
   bo[3]=cpr[3]*L3;
   dmax=0.0;

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
               sp=(spin_dble_t*)(s+ix);

               if (((x0+x1+x2+x3)&0x1)==0)
               {
                  for (ifc=0;ifc<8;ifc++)
                  {
                     mu=ifc/2;

                     if (ifc&0x1)
                        iy=iup[ix][mu];
                     else
                        iy=idn[ix][mu];

                     if ((iy>=VOLUME)&&
                         ((ifc>1)||
                          ((ifc==0)&&((cpr[0]>0)||(bc==3)))||
                          ((ifc==1)&&((cpr[0]<(NPROC0-1))||(bc==3)))))
                     {
                        pt=p[0]*(double)(x0+bo[0]);
                        pv=p[1]*(double)(x1+bo[1])+
                           p[2]*(double)(x2+bo[2])+
                           p[3]*(double)(x3+bo[3]);

                        if (bc==3)
                        {
                           z.re=cos(pt+pv);
                           z.im=sin(pt+pv);
                        }
                        else
                        {
                           z.re=sin(pt)*cos(pv);
                           z.im=sin(pt)*sin(pv);
                        }

                        r.s.c1=mul_cplx(z,&(rs.c1));
                        r.s.c2=mul_cplx(z,&(rs.c2));
                        r.s.c3=mul_cplx(z,&(rs.c3));
                        r.s.c4=mul_cplx(z,&(rs.c4));
                        r.s=theta((ifc^0x1)^is,&(r.s));

                        for (i=0;i<18;i++)
                           (*sp).r[i]-=r.r[i];
                     }
                  }

                  for (i=0;i<18;i++)
                  {
                     d=fabs((*sp).r[i]);
                     if (d>dmax)
                        dmax=d;
                  }
               }
               else
               {
                  for (i=0;i<18;i++)
                  {
                     d=fabs((*sp).r[i]);
                     if (d>dmax)
                        dmax=d;
                  }
               }
            }
         }
      }
   }

   if (NPROC>1)
   {
      d=dmax;
      MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return dmax;
}


int main(int argc,char *argv[])
{
   int my_rank,bc,ie,is,k;
   double phi[2],phi_prime[2],theta[3];
   double d,dmax;
   spinor **ps;
   spinor_dble **psd;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      printf("\n");
      printf(" Check of the communication programs in scom.c and sdcom.c\n");
      printf("----------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>]");
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime,theta);
   print_bc_parms(0);

   start_ranlux(0,12345);
   geometry();
   alloc_ws(NFLDS);
   alloc_wsd(NFLDS);

   ps=reserve_ws(NFLDS);
   psd=reserve_wsd(NFLDS);
   dmax=0.0;

   for (is=0;is<2;is++)
   {
      for (k=0;k<NFLDS;k+=2)
      {
         random_sd(NSPIN,psd[k],1.0);
         assign_sd2s(NSPIN,psd[k],ps[k]);
         d=(double)(norm_square(NSPIN,1,ps[k]));
         cps_int_bnd(is,ps[k]);
         cpsd_int_bnd(is,psd[k]);
         assign_sd2s(NSPIN,psd[k],ps[k+1]);
         mulr_spinor_add(NSPIN,ps[k],ps[k+1],-1.0f);
         d=(double)(norm_square(NSPIN,1,ps[k]))/d;
         d=sqrt(d);
         if (d>dmax)
            dmax=d;

         random_sd(NSPIN,psd[k],1.0);
         assign_sd2s(NSPIN,psd[k],ps[k]);
         d=(double)(norm_square(NSPIN,1,ps[k]));
         cps_ext_bnd(is,ps[k]);
         cpsd_ext_bnd(is,psd[k]);
         assign_sd2s(NSPIN,psd[k],ps[k+1]);
         mulr_spinor_add(NSPIN,ps[k],ps[k+1],-1.0f);
         d=(double)(norm_square(NSPIN,1,ps[k]))/d;
         d=sqrt(d);
         if (d>dmax)
            dmax=d;
      }
   }

   if (my_rank==0)
   {
      printf("Maximal relative deviation single-/double-precision programs"
             " = %.1e\n",dmax);
      printf("Now checking double-precision programs:\n");
   }

   ie=1;

   for (is=0;is<2;is++)
   {
      for (k=0;k<NFLDS;k+=2)
      {
         random_sd(NSPIN,psd[k],1.0);
         cpsd_int_bnd(is,psd[k]);
         ie&=check_int_bnd_dble(psd[k]);

         random_sd(NSPIN,psd[k],1.0);
         cpsd_ext_bnd(is,psd[k]);
         ie&=check_ext_bnd_dble(psd[k]);
      }
   }

   error(ie!=1,1,"main [check2.c]",
         "Spinor fields vanish on an incorrect set of points");

   dmax=0.0;

   for (is=0;is<2;is++)
   {
      set_sd(psd[0]);
      random_sd(NSPIN-VOLUME,psd[0]+VOLUME,1.0);
      assign_sd2sd(VOLUME,psd[0],psd[1]);
      cpsd_int_bnd(is,psd[0]);
      mulr_spinor_add_dble(VOLUME,psd[1],psd[0],-1.0);
      d=norm_square_dble(VOLUME,1,psd[1]);
      if (d>dmax)
         dmax=d;
      d=check_cpsd_int(is,psd[0]);
      if (d>dmax)
         dmax=d;
   }

   if (my_rank==0)
      printf("Maximal deviation (cpsd_int_bnd) = %.1e\n",dmax);

   dmax=0.0;

   for (is=0;is<2;is++)
   {
      random_sd(NSPIN,psd[0],1.0);
      set_sd_bnd(is,psd[0]);
      assign_sd2sd(NSPIN,psd[0],psd[1]);
      cpsd_ext_bnd(is,psd[0]);
      mulr_spinor_add_dble(NSPIN-VOLUME,psd[1]+VOLUME,psd[0]+VOLUME,-1.0);
      d=norm_square_dble(NSPIN-VOLUME,1,psd[1]+VOLUME);
      if (d>dmax)
         dmax=d;
      mulr_spinor_add_dble(VOLUME,psd[0],psd[1],-1.0);
      d=check_cpsd_ext(is,psd[0]);
      if (d>dmax)
         dmax=d;
   }

   if (my_rank==0)
   {
      printf("Maximal deviation (cpsd_ext_bnd) = %.1e\n\n",dmax);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
