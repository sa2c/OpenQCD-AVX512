
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2005, 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Gauge covariance of Dw_dble().
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

#define N0 (NPROC0*L0)

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


static void transform_sd(spinor_dble *pk,spinor_dble *pl)
{
   int ix;
   su3_dble gx;
   spinor_dble r,s;

   for (ix=0;ix<VOLUME;ix++)
   {
      s=pk[ix];
      gx=g[ix];

      _su3_multiply(r.c1,gx,s.c1);
      _su3_multiply(r.c2,gx,s.c2);
      _su3_multiply(r.c3,gx,s.c3);
      _su3_multiply(r.c4,gx,s.c4);

      pl[ix]=r;
   }
}


int main(int argc,char *argv[])
{
   int my_rank,i;
   double phi[2],phi_prime[2],theta[3];
   double mu,d;
   spinor_dble **psd;
   sw_parms_t swp;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);
      printf("\n");
      printf("Gauge covariance of Dw_dble() (random fields)\n");
      printf("---------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check4.c]",
                    "Syntax: check1 [-bc <type>]");
   }

   set_lat_parms(5.5,1.0,0,NULL,1.978);
   print_lat_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   theta[0]=0.35;
   theta[1]=-1.25;
   theta[2]=0.78;
   set_bc_parms(bc,0.55,0.78,0.9012,1.2034,phi,phi_prime,theta);
   print_bc_parms(2);

   start_ranlux(0,12345);
   geometry();
   alloc_wsd(5);
   psd=reserve_wsd(5);

   g=amalloc(NSPIN*sizeof(*g),4);
   if (BNDRY!=0)
      gbuf=amalloc((BNDRY/2)*sizeof(*gbuf),4);

   error((g==NULL)||((BNDRY!=0)&&(gbuf==NULL)),1,"main [check4.c]",
         "Unable to allocate auxiliary arrays");

   swp=set_sw_parms(-0.0123);
   mu=0.0376;

   if (my_rank==0)
      printf("m0 = %.4e, csw = %.4e, cF = %.4e, cF' = %.4e\n\n",
             swp.m0,swp.csw,swp.cF[0],swp.cF[1]);

   random_g();
   random_ud();
   set_ud_phase();
   sw_term(NO_PTS);

   for (i=0;i<4;i++)
      random_sd(NSPIN,psd[i],1.0);

   assign_sd2sd(VOLUME,psd[0],psd[4]);
   bnd_sd2zero(ALL_PTS,psd[4]);
   Dw_dble(mu,psd[0],psd[1]);
   mulr_spinor_add_dble(VOLUME,psd[4],psd[0],-1.0);
   d=norm_square_dble(VOLUME,1,psd[4]);
   error(d!=0.0,1,"main [check4.c]","Dw_dble() changes the input field");

   Dw_dble(mu,psd[0],psd[4]);
   mulr_spinor_add_dble(VOLUME,psd[4],psd[1],-1.0);
   d=norm_square_dble(VOLUME,1,psd[4]);
   error(d!=0.0,1,"main [check4.c]","Action of Dw_dble() depends "
         "on the boundary values of the input field");

   assign_sd2sd(VOLUME,psd[1],psd[4]);
   bnd_sd2zero(ALL_PTS,psd[4]);
   mulr_spinor_add_dble(VOLUME,psd[4],psd[1],-1.0);
   d=norm_square_dble(VOLUME,1,psd[4]);
   error(d!=0.0,1,"main [check4.c]",
         "Dw_dble() does not vanish at global time 0 and NPROC0*L0-1 ");

   transform_sd(psd[0],psd[2]);
   transform_ud();
   sw_term(NO_PTS);
   Dw_dble(mu,psd[2],psd[3]);
   transform_sd(psd[1],psd[2]);

   mulr_spinor_add_dble(VOLUME,psd[3],psd[2],-1.0);
   d=norm_square_dble(VOLUME,1,psd[3])/norm_square_dble(VOLUME,1,psd[0]);

   if (my_rank==0)
   {
      printf("Normalized difference = %.2e\n",sqrt(d));
      printf("(should be around 1*10^(-15) or so)\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
