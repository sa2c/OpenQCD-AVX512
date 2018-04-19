
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2012-2014, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Gauge action of constant Abelian background fields.
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
#include "forces.h"
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


static void Abnd(double *s0,double *s1)
{
   int k,x1,x2,x3,x[4];
   double r0,r1,r2,*cG,*phi;
   double rs0,rs1;
   bc_parms_t bcp;

   rs0=0.0;
   rs1=0.0;

   bcp=bc_parms();
   cG=bcp.cG;

   if ((bc==1)&&(cpr[0]==0))
   {
      x[0]=0;
      phi=bcp.phi[0];

      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               x[1]=bo[1]+x1;
               x[2]=bo[2]+x2;
               x[3]=bo[3]+x3;

               for (k=1;k<4;k++)
               {
                  r0=-afld(x,0);
                  x[0]+=1;
                  r2=r0-afld(x,0);
                  x[0]+=1;
                  r2-=afld(x,k);
                  x[0]-=1;
                  r0-=afld(x,k);
                  x[k]+=1;
                  r2+=afld(x,0);
                  r1=r0-afld(x,k);
                  x[0]-=1;
                  r0+=afld(x,0);
                  r2+=afld(x,0);
                  x[k]+=1;
                  r1+=afld(x,0);
                  x[k]-=2;

                  rs0-=cG[0]*(cos(phi[0]*inp[k]+r0)+
                              cos(phi[1]*inp[k]+r0)+
                              cos(phi[2]*inp[k]-2.0*r0)-3.0);

                  rs1-=0.5*(cos(2.0*(phi[0]*inp[k]+r0))+
                            cos(2.0*(phi[1]*inp[k]+r0))+
                            cos(2.0*(phi[2]*inp[k]-2.0*r0))-3.0);

                  rs1-=(cos(2.0*phi[0]*inp[k]+r1)+
                        cos(2.0*phi[1]*inp[k]+r1)+
                        cos(2.0*phi[2]*inp[k]-2.0*r1)-3.0);

                  rs1-=(cos(phi[0]*inp[k]+r2)+
                        cos(phi[1]*inp[k]+r2)+
                        cos(phi[2]*inp[k]-2.0*r2)-3.0);
               }
            }
         }
      }
   }

   if (((bc==1)||(bc==2))&&(cpr[0]==(NPROC0-1)))
   {
      x[0]=N0;
      phi=bcp.phi[1];

      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               x[1]=bo[1]+x1;
               x[2]=bo[2]+x2;
               x[3]=bo[3]+x3;

               for (k=1;k<4;k++)
               {
                  x[0]-=1;
                  r0=afld(x,0);
                  r2=r0;
                  r0-=afld(x,k);
                  x[0]-=1;
                  r2+=afld(x,0);
                  r2-=afld(x,k);
                  x[k]+=1;
                  r2-=afld(x,0);
                  x[0]+=1;
                  r1=r0-afld(x,k);
                  r0-=afld(x,0);
                  r2-=afld(x,0);
                  x[k]+=1;
                  r1-=afld(x,0);
                  x[0]+=1;
                  x[k]-=2;

                  rs0-=cG[1]*(cos(phi[0]*inp[k]+r0)+
                              cos(phi[1]*inp[k]+r0)+
                              cos(phi[2]*inp[k]-2.0*r0)-3.0);

                  rs1-=0.5*(cos(2.0*(phi[0]*inp[k]+r0))+
                            cos(2.0*(phi[1]*inp[k]+r0))+
                            cos(2.0*(phi[2]*inp[k]-2.0*r0))-3.0);

                  rs1-=(cos(2.0*phi[0]*inp[k]+r1)+
                        cos(2.0*phi[1]*inp[k]+r1)+
                        cos(2.0*phi[2]*inp[k]-2.0*r1)-3.0);

                  rs1-=(cos(phi[0]*inp[k]+r2)+
                        cos(phi[1]*inp[k]+r2)+
                        cos(phi[2]*inp[k]-2.0*r2)-3.0);
               }
            }
         }
      }
   }

   if (NPROC>1)
   {
      MPI_Reduce(&rs0,s0,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&rs1,s1,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(s0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(s1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
   {
      (*s0)=rs0;
      (*s1)=rs1;
   }
}


static double Amt(void)
{
   int mu,nu;
   double c0,c1,*cG;
   double smt0,smt1,sms0,sms1,pi;
   double xl[4],phi,n0,s0,s1,bs0,bs1;
   lat_parms_t lat;
   bc_parms_t bcp;

   lat=lat_parms();
   c0=lat.c0;
   c1=lat.c1;
   bcp=bc_parms();
   cG=bcp.cG;

   xl[0]=(double)(N0);
   xl[1]=(double)(N1);
   xl[2]=(double)(N2);
   xl[3]=(double)(N3);

   pi=4.0*atan(1.0);
   smt0=0.0;
   smt1=0.0;
   sms0=0.0;
   sms1=0.0;

   for (mu=1;mu<4;mu++)
   {
      for (nu=0;nu<mu;nu++)
      {
         phi=2.0*pi*mt[mu][nu]/(xl[mu]*xl[nu]);

         s0=3.0-2.0*cos(phi)-cos(2.0*phi);
         s1=3.0-2.0*cos(2.0*phi)-cos(4.0*phi);

         if (nu==0)
         {
            smt0+=s0;
            smt1+=s1;
         }
         else
         {
            sms0+=s0;
            sms1+=s1;
         }
      }
   }

   n0=(double)(N0);

   if (bc==0)
   {
      s0=(n0-1.0)*smt0+(n0-2.0+0.5*cG[0]+0.5*cG[1])*sms0;
      s1=(n0-1.0)*smt1+(n0-2.0)*smt1+(2.0*(n0-2.0)+cG[0]+cG[1])*sms1;
   }
   else if (bc==1)
   {
      s0=(n0-2.0)*smt0+(n0-1.0)*sms0;
      s1=(n0-2.0)*smt1+(n0-3.0)*smt1+2.0*(n0-1.0)*sms1;
   }
   else if (bc==2)
   {
      s0=(n0-1.0)*smt0+(n0-1.0+0.5*cG[0])*sms0;
      s1=(n0-1.0)*smt1+(n0-2.0)*smt1+(2.0*(n0-1.0)+cG[0])*sms1;
   }
   else
   {
      s0=n0*smt0+n0*sms0;
      s1=2.0*n0*smt1+2.0*n0*sms1;
   }

   s0*=(double)(N1*N2*N3);
   s1*=(double)(N1*N2*N3);

   Abnd(&bs0,&bs1);

   return (lat.beta/3.0)*(c0*(s0+bs0)+c1*(s1+bs1));
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

   set_flags(UPDATED_UD);
   set_flags(UNSET_UD_PHASE);
   set_bc();
}


int main(int argc,char *argv[])
{
   int my_rank,i;
   double A1,A2,d,dmax;
   double phi[2],phi_prime[2],theta[3];
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      printf("\n");
      printf("Gauge action of constant Abelian background fields\n");
      printf("--------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>]");
   }

   set_lat_parms(3.5,0.33,0,NULL,1.0);
   print_lat_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;
   set_bc_parms(bc,0.9012,1.2034,1.0,1.0,phi,phi_prime,theta);
   print_bc_parms(1);

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
      A2=action0(1);

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
