
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2009-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Gauge covariance of the Wilson flow.
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
#include "update.h"
#include "wflow.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int bc,nfc[8],ofs[8];
static const su3_dble ud0={{0.0}};
static su3_dble *g,*gbuf;
static su3_dble wd ALIGNED16;


static void pack_gbuf(void)
{
   int ifc,ib,ix;

   nfc[0]=FACE0/2;
   nfc[1]=FACE0/2;
   nfc[2]=FACE1/2;
   nfc[3]=FACE1/2;
   nfc[4]=FACE2/2;
   nfc[5]=FACE2/2;
   nfc[6]=FACE3/2;
   nfc[7]=FACE3/2;

   ofs[0]=0;
   ofs[1]=ofs[0]+nfc[0];
   ofs[2]=ofs[1]+nfc[1];
   ofs[3]=ofs[2]+nfc[2];
   ofs[4]=ofs[3]+nfc[3];
   ofs[5]=ofs[4]+nfc[4];
   ofs[6]=ofs[5]+nfc[5];
   ofs[7]=ofs[6]+nfc[6];

   for (ifc=0;ifc<8;ifc++)
   {
      for (ib=0;ib<nfc[ifc];ib++)
      {
         ix=map[ofs[ifc]+ib];
         gbuf[ofs[ifc]+ib]=g[ix];
      }
   }
}


static void send_gbuf(void)
{
   int ifc,np,saddr,raddr;
   int nbf,tag;
   su3_dble *sbuf,*rbuf;
   MPI_Status stat;

   np=cpr[0]+cpr[1]+cpr[2]+cpr[3];

   for (ifc=0;ifc<8;ifc++)
   {
      nbf=18*nfc[ifc];

      if (nbf>0)
      {
         tag=mpi_tag();
         saddr=npr[ifc^0x1];
         raddr=npr[ifc];
         sbuf=gbuf+ofs[ifc];
         rbuf=g+VOLUME+ofs[ifc];

         if (np&0x1)
         {
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         }
      }
   }
}


static void random_g(void)
{
   int ix,t;
   su3_dble unity,*gx;

   unity=ud0;
   unity.c11.re=1.0;
   unity.c22.re=1.0;
   unity.c33.re=1.0;
   gx=g;

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((t>0)||(bc!=1))
         random_su3_dble(gx);
      else
         (*gx)=unity;

      gx+=1;
   }

   if (BNDRY>0)
   {
      pack_gbuf();
      send_gbuf();
   }
}


static void transform_ud(void)
{
   int ix,iy,t,ifc;
   su3_dble *u;

   u=udfld();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         iy=iup[ix][0];
         su3xsu3dag(u,g+iy,&wd);
         su3xsu3(g+ix,&wd,u);
         u+=1;

         if (bc==3)
         {
            iy=idn[ix][0];
            su3xsu3dag(u,g+ix,&wd);
            su3xsu3(g+iy,&wd,u);
         }
         else if (bc!=0)
         {
            iy=idn[ix][0];
            su3xsu3(g+iy,u,&wd);
            (*u)=wd;
         }

         u+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
            {
               if (ifc&0x1)
               {
                  iy=idn[ix][ifc/2];
                  su3xsu3dag(u,g+ix,&wd);
                  su3xsu3(g+iy,&wd,u);
               }
               else
               {
                  iy=iup[ix][ifc/2];
                  su3xsu3dag(u,g+iy,&wd);
                  su3xsu3(g+ix,&wd,u);
               }
            }

            u+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc==3)
         {
            iy=iup[ix][0];
            su3xsu3dag(u,g+iy,&wd);
            su3xsu3(g+ix,&wd,u);
         }
         else if (bc!=0)
         {
            su3xsu3(g+ix,u,&wd);
            (*u)=wd;
         }

         u+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            if (ifc&0x1)
            {
               iy=idn[ix][ifc/2];
               su3xsu3dag(u,g+ix,&wd);
               su3xsu3(g+iy,&wd,u);
            }
            else
            {
               iy=iup[ix][ifc/2];
               su3xsu3dag(u,g+iy,&wd);
               su3xsu3(g+ix,&wd,u);
            }

            u+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            if (ifc&0x1)
            {
               iy=idn[ix][ifc/2];
               su3xsu3dag(u,g+ix,&wd);
               su3xsu3(g+iy,&wd,u);
            }
            else
            {
               iy=iup[ix][ifc/2];
               su3xsu3dag(u,g+iy,&wd);
               su3xsu3(g+ix,&wd,u);
            }

            u+=1;
         }
      }
   }

   set_flags(UPDATED_UD);
}


static double cmp_ud(su3_dble *u,su3_dble *v)
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

   for (i=0;i<18;i+=2)
   {
      dev=r[i]*r[i]+r[i+1]*r[i+1];
      if (dev>dmax)
         dmax=dev;
   }

   return dmax;
}


static double max_dev_ud(su3_dble *v)
{
   double d,dmax;
   su3_dble *u,*um;

   u=udfld();
   um=u+4*VOLUME;
   dmax=0.0;

   for (;u<um;u++)
   {
      d=cmp_ud(u,v);

      if (d>dmax)
         dmax=d;

      v+=1;
   }

   if (NPROC>1)
   {
      d=dmax;
      MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return sqrt(dmax);
}


int main(int argc,char *argv[])
{
   int my_rank,n;
   double phi[2],phi_prime[2],theta[3];
   double eps,dev;
   su3_dble *udb,**usv;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      fin=freopen("check2.in","r",stdin);

      printf("\n");
      printf("Gauge covariance of the Wilson flow\n");
      printf("-----------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("n","%d\n",&n);
      read_line("eps","%lf",&eps);
      fclose(fin);

      printf("n = %d\n",n);
      printf("eps = %.3e\n\n",eps);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>]");
   }

   MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&eps,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);

   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   theta[0]=0.0;
   theta[1]=0.0;
   theta[2]=0.0;
   set_bc_parms(bc,0.973,1.127,1.0,1.0,phi,phi_prime,theta);
   print_bc_parms(0);

   start_ranlux(0,1234);
   geometry();
   alloc_wud(2);
   alloc_wfd(1);
   usv=reserve_wud(2);
   udb=udfld();

   g=amalloc(NSPIN*sizeof(*g),4);

   if (BNDRY>0)
      gbuf=amalloc((BNDRY/2)*sizeof(*gbuf),4);

   error((g==NULL)||((BNDRY>0)&&(gbuf==NULL)),1,"main [check2.c]",
         "Unable to allocate auxiliary arrays");

   random_ud();
   random_g();
   cm3x3_assign(4*VOLUME,udb,usv[0]);
   fwd_euler(n,eps);
   transform_ud();
   cm3x3_assign(4*VOLUME,udb,usv[1]);
   cm3x3_assign(4*VOLUME,usv[0],udb);
   set_flags(UPDATED_UD);
   transform_ud();
   fwd_euler(n,eps);

   dev=max_dev_ud(usv[1]);

   if (my_rank==0)
   {
      printf("Maximal absolute deviation of U(x,mu) = %.1e\n\n",dev);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
