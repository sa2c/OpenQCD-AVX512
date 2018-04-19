
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2005, 2007-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the programs for the plaquette sums of the double-precision
* gauge field.
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
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int bc,nfc[8],ofs[8];
static double asl1[N0],asl2[N0];
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
   int my_rank,n,t,s[4];
   double phi[2],phi_prime[2],theta[3];
   double act1,nplaq1,nplaq2,p1,p2;
   double d1,d2,d3;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);

      printf("\n");
      printf("Plaquette sums of the double-precision gauge field\n");
      printf("--------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check4.c]",
                    "Syntax: check4 [-bc <type>]");
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

   start_ranlux(0,12345);
   geometry();

   g=amalloc(NSPIN*sizeof(*g),4);

   if (BNDRY>0)
      gbuf=amalloc((BNDRY/2)*sizeof(*gbuf),4);

   error((g==NULL)||((BNDRY>0)&&(gbuf==NULL)),1,"main [check4.c]",
         "Unable to allocate auxiliary arrays");

   p1=plaq_sum_dble(1);
   p2=plaq_wsum_dble(1);

   if (bc==0)
   {
      nplaq1=(double)((6*N0-3)*N1)*(double)(N2*N3);
      nplaq2=(double)((6*N0-6)*N1)*(double)(N2*N3);
   }
   else if (bc==3)
   {
      nplaq1=(double)(6*N0*N1)*(double)(N2*N3);
      nplaq2=nplaq1;
   }
   else
   {
      nplaq1=(double)((6*N0+3)*N1)*(double)(N2*N3);
      nplaq2=(double)(6*N0*N1)*(double)(N2*N3);
   }

   d1=0.0;
   d2=0.0;

   if (bc==1)
   {
      d1=cos(phi[0]/(double)(N1))+
         cos(phi[1]/(double)(N1))+
         cos((phi[0]+phi[1])/(double)(N1))+
         cos(phi[0]/(double)(N2))+
         cos(phi[1]/(double)(N2))+
         cos((phi[0]+phi[1])/(double)(N2))+
         cos(phi[0]/(double)(N3))+
         cos(phi[1]/(double)(N3))+
         cos((phi[0]+phi[1])/(double)(N3));

      d1=(d1-9.0)*(double)(N1*N2*N3);
   }

   if ((bc==1)||(bc==2))
   {
      d2=cos(phi_prime[0]/(double)(N1))+
         cos(phi_prime[1]/(double)(N1))+
         cos((phi_prime[0]+phi_prime[1])/(double)(N1))+
         cos(phi_prime[0]/(double)(N2))+
         cos(phi_prime[1]/(double)(N2))+
         cos((phi_prime[0]+phi_prime[1])/(double)(N2))+
         cos(phi_prime[0]/(double)(N3))+
         cos(phi_prime[1]/(double)(N3))+
         cos((phi_prime[0]+phi_prime[1])/(double)(N3));

      d2=(d2-9.0)*(double)(N1*N2*N3);
   }

   if (my_rank==0)
   {
      printf("After field initialization:\n");
      printf("Deviation from expected value (plaq_sum)  = %.1e\n",
             fabs(1.0-p1/(3.0*nplaq1+d1+d2)));
      printf("Deviation from expected value (plaq_wsum) = %.1e\n\n",
             fabs(1.0-p2/(3.0*nplaq2+d1+d2)));
   }

   print_flags();
   random_ud();

   p1=plaq_sum_dble(1);
   p2=plaq_wsum_dble(1);
   act1=plaq_action_slices(asl1);
   d1=act1;

   if ((bc==0)||(bc==3))
   {
      for (t=0;t<N0;t++)
         d1-=asl1[t];
   }

   if (my_rank==0)
   {
      printf("Comparison of plaq_wsum_dble() with plaq_action_slices():\n");
      printf("Absolute difference of total action = %.1e\n",
             fabs(3.0*nplaq2-0.5*act1-p2));
      if ((bc==0)||(bc==3))
         printf("Deviation from sum of action slices = %.1e\n\n",
                fabs(d1));
      else
         printf("\n");
   }

   random_g();
   transform_ud();
   d1=fabs(p1-plaq_sum_dble(1));
   d2=fabs(p2-plaq_wsum_dble(1));
   plaq_action_slices(asl2);
   d3=0.0;

   for (t=0;t<N0;t++)
      d3+=fabs(asl1[t]-asl2[t]);

   if (my_rank==0)
   {
      printf("Gauge invariance:\n");
      printf("Relative difference (plaq_sum_dble)  = %.1e\n",d1/fabs(p1));
      printf("Relative difference (plaq_wsum_dble) = %.1e\n",d2/fabs(p2));
      printf("Relative difference (action slices)  = %.1e\n\n",
             d3/((double)(N0)*asl2[1]));
   }

   if (my_rank==0)
      printf("Translation invariance:\n");

   random_ud();
   p1=plaq_sum_dble(1);
   p2=plaq_wsum_dble(1);
   plaq_action_slices(asl1);

   for (n=0;n<8;n++)
   {
      random_vec(s);
      if (bc!=3)
         s[0]=0;
      shift_ud(s);
      d1=fabs(p1-plaq_sum_dble(1));
      d2=fabs(p2-plaq_wsum_dble(1));
      plaq_action_slices(asl2);
      d3=0.0;

      for (t=0;t<N0;t++)
         d3+=fabs(asl1[safe_mod(t-s[0],N0)]-asl2[t]);

      for (t=0;t<N0;t++)
         asl1[t]=asl2[t];

      if (my_rank==0)
      {
         printf("s=(% 3d,% 3d,% 3d,% 3d):\n",s[0],s[1],s[2],s[3]);
         printf("Absolute deviation (plaq_sum_dble)  = %.1e\n",d1);
         printf("Absolute deviation (plaq_wsum_dble) = %.1e\n",d2);
         printf("Absolute difference (action slices) = %.1e\n\n",
                d2/(double)(N0));
      }
   }

   if (bc==1)
   {
      random_ud();
      p1=plaq_sum_dble(1);
      p2=plaq_wsum_dble(1);

      if (my_rank==0)
      {
         printf("\n");
         printf("Comparison of plaq_sum_dble() and plaq_wsum_dble():\n");
         printf("Absolute deviation = %.1e\n\n",
                fabs(p1-p2-9.0*(double)(N1*N2*N3)));
      }
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
