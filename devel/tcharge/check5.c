
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2010, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Yang-Mills action of constant abelian background fields.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "su3fcts.h"
#include "random.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "tcharge.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int bc,np[4],bo[4];
static double mt[4][4],inp[4],twopi;
static su3_dble ud0={{0.0}};


static double afld(int *x,int mu)
{
   int nu;
   double xt[4],phi;

   xt[0]=(double)(safe_mod(x[0],N0));
   xt[1]=(double)(safe_mod(x[1],N1));
   xt[2]=(double)(safe_mod(x[2],N2));
   xt[3]=(double)(safe_mod(x[3],N3));

   phi=0.0;

   for (nu=0;nu<mu;nu++)
      phi-=inp[nu]*mt[mu][nu]*xt[nu];

   phi*=inp[mu];

   if (safe_mod(x[mu],np[mu])==(np[mu]-1))
   {
      for (nu=(mu+1);nu<4;nu++)
         phi-=inp[nu]*mt[mu][nu]*xt[nu];
   }

   return twopi*phi;
}


static void ftplaq(int *x,int mu,int nu,double *ftp)
{
   double sm,om[3],*phi;
   bc_parms_t bcp;

   bcp=bc_parms();

   if ((x[0]==0)&&(mu==0)&&(bc==1))
   {
      sm=afld(x,mu);
      x[mu]+=1;
      sm+=afld(x,nu);
      x[mu]-=1;
      x[nu]+=1;
      sm-=afld(x,mu);
      x[nu]-=1;

      phi=bcp.phi[0];
      om[0]=sm-phi[0]*inp[nu];
      om[1]=sm-phi[1]*inp[nu];
      om[2]=-2.0*sm-phi[2]*inp[nu];
   }
   else if ((x[0]==(N0-1))&&(mu==0)&&((bc==1)||(bc==2)))
   {
      sm=afld(x,mu)-afld(x,nu);
      x[nu]+=1;
      sm-=afld(x,mu);
      x[nu]-=1;

      phi=bcp.phi[1];
      om[0]=sm+phi[0]*inp[nu];
      om[1]=sm+phi[1]*inp[nu];
      om[2]=-2.0*sm+phi[2]*inp[nu];
   }
   else
   {
      sm=afld(x,mu)-afld(x,nu);
      x[mu]+=1;
      sm+=afld(x,nu);
      x[mu]-=1;
      x[nu]+=1;
      sm-=afld(x,mu);
      x[nu]-=1;

      om[0]=sm;
      om[1]=sm;
      om[2]=-2.0*sm;
   }

   ftp[0]=sin(om[0]);
   ftp[1]=sin(om[1]);
   ftp[2]=sin(om[2]);
}


static double Abnd(void)
{
   int ib,x1,x2,x3,x[4];
   int mu,nu;
   double ft1,ft2,ft3,tr;
   double r0[3],r1[3],r2[3],r3[3];
   double aloc,aall;

   aloc=0.0;

   for (ib=0;ib<2;ib++)
   {
      if (ib==0)
         x[0]=1;
      else
         x[0]=N0-1;

      if (((ib==0)&&(bc==1)&&(cpr[0]==0))||
          ((ib==1)&&((bc==1)||(bc==2))&&(cpr[0]==(NPROC0-1))))
      {
         for (x1=0;x1<L1;x1++)
         {
            for (x2=0;x2<L2;x2++)
            {
               for (x3=0;x3<L3;x3++)
               {
                  x[1]=bo[1]+x1;
                  x[2]=bo[2]+x2;
                  x[3]=bo[3]+x3;

                  for (mu=0;mu<3;mu++)
                  {
                     for (nu=(mu+1);nu<4;nu++)
                     {
                        ftplaq(x,mu,nu,r0);
                        x[mu]-=1;
                        ftplaq(x,mu,nu,r1);
                        x[nu]-=1;
                        ftplaq(x,mu,nu,r2);
                        x[mu]+=1;
                        ftplaq(x,mu,nu,r3);
                        x[nu]+=1;

                        ft1=0.25*(r0[0]+r1[0]+r2[0]+r3[0]);
                        ft2=0.25*(r0[1]+r1[1]+r2[1]+r3[1]);
                        ft3=0.25*(r0[2]+r1[2]+r2[2]+r3[2]);

                        tr=(ft1+ft2+ft3)/3.0;
                        ft1-=tr;
                        ft2-=tr;
                        ft3-=tr;

                        aloc+=ft1*ft1+ft2*ft2+ft3*ft3;
                     }
                  }
               }
            }
         }
      }
   }

   if (NPROC>1)
   {
      MPI_Reduce(&aloc,&aall,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&aall,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      return aall;
   }
   else
      return aloc;
}


static double Amt(void)
{
   int mu,nu;
   double sm,pi;
   double xl[4],phi,ft1,ft2,ft3,tr;

   xl[0]=(double)(NPROC0*L0);
   xl[1]=(double)(NPROC1*L1);
   xl[2]=(double)(NPROC2*L2);
   xl[3]=(double)(NPROC3*L3);

   pi=4.0*atan(1.0);
   sm=0.0;

   for (mu=1;mu<4;mu++)
   {
      for (nu=0;nu<mu;nu++)
      {
         phi=2.0*pi*mt[mu][nu]/(xl[mu]*xl[nu]);

         ft1=sin(phi);
         ft2=ft1;
         ft3=-sin(2.0*phi);

         tr=(ft1+ft2+ft3)/3.0;

         ft1-=tr;
         ft2-=tr;
         ft3-=tr;

         sm+=ft1*ft1+ft2*ft2+ft3*ft3;
      }
   }

   if (bc==0)
      sm*=(double)((N0-2)*N1)*(double)(N2*N3);
   else if (bc==1)
   {
      sm*=(double)((N0-3)*N1)*(double)(N2*N3);
      sm+=Abnd();
   }
   else if (bc==2)
   {
      sm*=(double)((N0-2)*N1)*(double)(N2*N3);
      sm+=Abnd();
   }
   else
      sm*=(double)(N0*N1)*(double)(N2*N3);

   return sm;
}


static void choose_mt(void)
{
   int mu,nu;
   double r[6];

   ranlxd(r,6);
   MPI_Bcast(r,6,MPI_DOUBLE,0,MPI_COMM_WORLD);

   mt[0][1]=(double)((int)(3.0*r[0])-1);
   mt[0][2]=(double)((int)(3.0*r[1])-1);
   mt[0][3]=(double)((int)(3.0*r[2])-1);
   mt[1][2]=(double)((int)(3.0*r[3])-1);
   mt[1][3]=(double)((int)(3.0*r[4])-1);
   mt[2][3]=(double)((int)(3.0*r[5])-1);

   for (mu=0;mu<4;mu++)
   {
      mt[mu][mu]=0.0;

      for (nu=0;nu<mu;nu++)
         mt[mu][nu]=-mt[nu][mu];
   }
}


static void set_ud(void)
{
   int x[4];
   int x0,x1,x2,x3;
   int ix,ifc;
   double phi;
   su3_dble *udb,*u;

   udb=udfld();

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

               if (ix>=(VOLUME/2))
               {
                  x[0]=bo[0]+x0;
                  x[1]=bo[1]+x1;
                  x[2]=bo[2]+x2;
                  x[3]=bo[3]+x3;

                  u=udb+8*(ix-(VOLUME/2));

                  for (ifc=0;ifc<8;ifc++)
                  {
                     if (ifc&0x1)
                        x[ifc/2]-=1;

                     phi=afld(x,ifc/2);

                     if (ifc&0x1)
                        x[ifc/2]+=1;

                     (*u)=ud0;
                     (*u).c11.re=cos(phi);
                     (*u).c11.im=sin(phi);
                     (*u).c22.re=(*u).c11.re;
                     (*u).c22.im=(*u).c11.im;
                     (*u).c33.re=cos(-2.0*phi);
                     (*u).c33.im=sin(-2.0*phi);
                     u+=1;
                  }
               }
            }
         }
      }
   }

   set_bc();
   set_flags(UPDATED_UD);
}


int main(int argc,char *argv[])
{
   int my_rank,i;
   double phi[2],phi_prime[2],theta[3];
   double A1,A2,d,dmax;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check5.log","w",stdout);
      printf("\n");
      printf("Yang-Mills action of constant abelian background fields\n");
      printf("-------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check5.c]",
                    "Syntax: check5 [-bc <type>]");
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime,theta);
   print_bc_parms(0);

   start_ranlux(0,123);
   geometry();

   twopi=8.0*atan(1.0);

   np[0]=N0;
   np[1]=N1;
   np[2]=N2;
   np[3]=N3;

   bo[0]=cpr[0]*L0;
   bo[1]=cpr[1]*L1;
   bo[2]=cpr[2]*L2;
   bo[3]=cpr[3]*L3;

   inp[0]=1.0/(double)(np[0]);
   inp[1]=1.0/(double)(np[1]);
   inp[2]=1.0/(double)(np[2]);
   inp[3]=1.0/(double)(np[3]);

   dmax=0.0;

   for (i=0;i<10;i++)
   {
      choose_mt();
      set_ud();

      A1=Amt();
      A2=ym_action();

      if (my_rank==0)
         printf("Field no = %2d, A1 = %12.6e, A2 = %12.6e\n",i+1,A1,A2);

      d=fabs(A1-A2)/A1;
      if (d>dmax)
         dmax=d;
   }

   if (my_rank==0)
   {
      printf("\n");
      printf("Maximal relative deviation = %.1e\n\n",dmax);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
