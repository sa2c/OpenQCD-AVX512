
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2012, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Gauge and translation invariance of the gauge action.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
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

static int bc,nfc[8],ofs[8];
static const su3_dble ud0={{0.0}};
static su3_dble *g,*gbuf;
static su3_dble wd ALIGNED16;


static double bnd_action(void)
{
   int i,j;
   double c0,c1,*cG,*phi;
   double s[3],d0[2],d1[2],act;
   lat_parms_t lat;
   bc_parms_t bcp;

   if ((bc==1)||(bc==2))
   {
      lat=lat_parms();
      bcp=bc_parms();

      s[0]=(double)(N1);
      s[1]=(double)(N2);
      s[2]=(double)(N3);

      for (i=0;i<2;i++)
      {
         d0[i]=0.0;
         d1[i]=0.0;
         phi=bcp.phi[i];

         for (j=0;j<3;j++)
         {
            d0[i]-=(cos(phi[0]/s[j])+cos(phi[1]/s[j])+
                    cos(phi[2]/s[j])-3.0);
            d1[i]-=(cos(2.0*phi[0]/s[j])+cos(2.0*phi[1]/s[j])+
                    cos(2.0*phi[2]/s[j])-3.0);
         }
      }

      c0=lat.c0;
      c1=lat.c1;
      cG=bcp.cG;

      act=c0*cG[1]*d0[1]+c1*d0[1]+c1*1.5*d1[1];

      if (bc==1)
         act+=(c0*cG[0]*d0[0]+c1*d0[0]+c1*1.5*d1[0]);

      return (lat.beta/3.0)*(double)(N1*N2*N3)*act;
   }
   else
      return 0.0;
}



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


static void random_vec(int *svec)
{
   int mu,bs[4];
   double r[4];

   bs[0]=NPROC0*L0;
   bs[1]=NPROC1*L1;
   bs[2]=NPROC2*L2;
   bs[3]=NPROC3*L3;

   ranlxd(r,4);

   for (mu=0;mu<4;mu++)
   {
      svec[mu]=(int)((double)(bs[mu])*r[mu]);
      if (svec[mu]>(bs[mu]/2))
         svec[mu]-=bs[mu];
   }

   MPI_Bcast(svec,4,MPI_INT,0,MPI_COMM_WORLD);
}


int main(int argc,char *argv[])
{
   int my_rank,n,s[4];
   double phi[2],phi_prime[2],theta[3],p1,p2;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);

      printf("\n");
      printf("Gauge and translation invariance of the gauge action\n");
      printf("----------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>]");
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

   start_ranlux(0,12345);
   geometry();

   g=amalloc(NSPIN*sizeof(*g),ALIGN);
   if (BNDRY!=0)
      gbuf=amalloc((BNDRY/2)*sizeof(*gbuf),ALIGN);

   error((g==NULL)||((BNDRY!=0)&&(gbuf==NULL)),1,"main [check1.c]",
         "Unable to allocate auxiliary arrays");

   p1=action0(1);
   p2=bnd_action();

   if (my_rank==0)
   {
      printf("Action after initialization = %.15e\n",p1);
      printf("Expected value              = %.15e\n\n",p2);
   }

   random_ud();
   p1=action0(1);
   random_g();
   transform_ud();
   p2=action0(1);

   if (my_rank==0)
   {
      printf("Random gauge field:\n");
      printf("Action = %.12e\n",p1);
      printf("Gauge invariance: relative difference = %.1e\n\n",
             fabs(1.0-p2/p1));
   }

   if (my_rank==0)
      printf("Translation invariance:\n");

   p1=action0(1);

   for (n=0;n<8;n++)
   {
      random_vec(s);
      if (bc!=3)
         s[0]=0;
      shift_ud(s);
      p2=action0(1);

      if (my_rank==0)
      {
         printf("s=(% d, % d,% d,% d), ",s[0],s[1],s[2],s[3]);
         printf("relative deviation = %.1e\n",fabs(1.0-p2/p1));
      }
   }

   if (my_rank==0)
   {
      printf("\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
