
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2012, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the program set_bstap().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "global.h"

static const int plns[6][2]={{0,1},{0,2},{0,3},{2,3},{3,1},{1,2}};
static int bc,nfc[8],ofs[8],hofs[8];
static double psum0[8],psum1[8];
static su3_dble *udb,*hdb;
static su3_dble wd1 ALIGNED16;
static su3_dble wd2 ALIGNED16;


static void set_ofs(void)
{
   nfc[0]=FACE0/2;
   nfc[1]=FACE0/2;
   nfc[2]=FACE1/2;
   nfc[3]=FACE1/2;
   nfc[4]=FACE2/2;
   nfc[5]=FACE2/2;
   nfc[6]=FACE3/2;
   nfc[7]=FACE3/2;

   ofs[0]=0;
   ofs[1]=ofs[0]+(FACE0/2);
   ofs[2]=ofs[1]+(FACE0/2);
   ofs[3]=ofs[2]+(FACE1/2);
   ofs[4]=ofs[3]+(FACE1/2);
   ofs[5]=ofs[4]+(FACE2/2);
   ofs[6]=ofs[5]+(FACE2/2);
   ofs[7]=ofs[6]+(FACE3/2);

   hofs[0]=0;
   hofs[1]=hofs[0]+3*FACE0;
   hofs[2]=hofs[1]+3*FACE0;
   hofs[3]=hofs[2]+3*FACE1;
   hofs[4]=hofs[3]+3*FACE1;
   hofs[5]=hofs[4]+3*FACE2;
   hofs[6]=hofs[5]+3*FACE2;
   hofs[7]=hofs[6]+3*FACE3;
}


static double plaq0(int n,int ix)
{
   int ip[4];
   double sm;

   plaq_uidx(n,ix,ip);

   su3xsu3(udb+ip[0],udb+ip[1],&wd1);
   su3dagxsu3dag(udb+ip[3],udb+ip[2],&wd2);
   cm3x3_retr(&wd1,&wd2,&sm);

   return sm;
}


static double plaq1(int iu,int ih)
{
   su3xsu3dag(udb+iu,hdb+ih,&wd1);

   return wd1.c11.re+wd1.c22.re+wd1.c33.re;
}


static void set_psum0(void)
{
   int ifc,n,ix,mu,nu;

   for (ifc=0;ifc<8;ifc++)
      psum0[ifc]=0.0;

   for (ix=0;ix<VOLUME;ix++)
   {
      for (n=0;n<6;n++)
      {
         mu=plns[n][0];
         nu=plns[n][1];

         if (iup[ix][mu]>=VOLUME)
            psum0[2*mu+1]+=plaq0(n,ix);

         if (iup[ix][nu]>=VOLUME)
            psum0[2*nu+1]+=plaq0(n,ix);

         if (idn[ix][mu]>=VOLUME)
            psum0[2*mu]+=plaq0(n,ix);

         if (idn[ix][nu]>=VOLUME)
            psum0[2*nu]+=plaq0(n,ix);
      }
   }
}


static void set_psum1(void)
{
   int ifc,n,ix,mu,nu,ip[4];
   int iy,ib,iu,ih;

   for (ifc=0;ifc<8;ifc++)
      psum1[ifc]=0.0;

   for (ix=0;ix<VOLUME;ix++)
   {
      for (n=0;n<6;n++)
      {
         mu=plns[n][0];
         nu=plns[n][1];

         if (iup[ix][mu]>=VOLUME)
         {
            plaq_uidx(n,ix,ip);
            iu=ip[1];

            ifc=2*mu+1;
            iy=iup[ix][mu]-VOLUME;

            if (iy<(BNDRY/2))
               ib=iy-ofs[ifc];
            else
               ib=iy-ofs[ifc]-(BNDRY/2)+nfc[ifc];

            ih=hofs[ifc]+3*ib+nu-(nu>mu);

            psum1[ifc]+=plaq1(iu,ih);
         }

         if (iup[ix][nu]>=VOLUME)
         {
            plaq_uidx(n,ix,ip);
            iu=ip[3];

            ifc=2*nu+1;
            iy=iup[ix][nu]-VOLUME;

            if (iy<(BNDRY/2))
               ib=iy-ofs[ifc];
            else
               ib=iy-ofs[ifc]-(BNDRY/2)+nfc[ifc];

            ih=hofs[ifc]+3*ib+mu-(mu>nu);

            psum1[ifc]+=plaq1(iu,ih);
         }

         if (idn[ix][mu]>=VOLUME)
         {
            plaq_uidx(n,ix,ip);
            iu=ip[2];

            ifc=2*mu;
            iy=idn[ix][mu]-VOLUME;

            if (iy<(BNDRY/2))
               ib=iy-ofs[ifc];
            else
               ib=iy-ofs[ifc]-(BNDRY/2)+nfc[ifc];

            ih=hofs[ifc]+3*ib+nu-(nu>mu);

            psum1[ifc]+=plaq1(iu,ih);
         }

         if (idn[ix][nu]>=VOLUME)
         {
            plaq_uidx(n,ix,ip);
            iu=ip[0];

            ifc=2*nu;
            iy=idn[ix][nu]-VOLUME;

            if (iy<(BNDRY/2))
               ib=iy-ofs[ifc];
            else
               ib=iy-ofs[ifc]-(BNDRY/2)+nfc[ifc];

            ih=hofs[ifc]+3*ib+mu-(mu>nu);

            psum1[ifc]+=plaq1(iu,ih);
         }
      }
   }
}


static void check_psums(void)
{
   int ifc,np;
   int saddr,raddr,nbf,tag;
   double sbuf,rbuf,dmy[8];
   MPI_Status stat;

   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;

   for (ifc=0;ifc<8;ifc++)
   {
      if (nfc[ifc]>0)
      {
         saddr=npr[ifc];
         raddr=npr[ifc^0x1];
         sbuf=psum0[ifc];
         nbf=1;
         tag=mpi_tag();

         if (np==0)
         {
            MPI_Send(&sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
            MPI_Recv(&rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(&rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
            MPI_Send(&sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         }

         if ((bc!=3)&&
             (((cpr[0]==0)&&(ifc==1))||((cpr[0]==(NPROC0-1))&&(ifc==0))))
            psum1[ifc^0x1]=0.0;
         else
            psum1[ifc^0x1]-=rbuf;
      }
   }

   for (ifc=0;ifc<8;ifc++)
      dmy[ifc]=fabs(psum0[ifc]);

   MPI_Reduce(dmy,psum0,8,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Bcast(psum0,8,MPI_DOUBLE,0,MPI_COMM_WORLD);

   for (ifc=0;ifc<8;ifc++)
      dmy[ifc]=fabs(psum1[ifc]);

   MPI_Reduce(dmy,psum1,8,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Bcast(psum1,8,MPI_DOUBLE,0,MPI_COMM_WORLD);
}


int main(int argc,char *argv[])
{
   int my_rank,ifc,ie;
   double phi[2],phi_prime[2],theta[3];
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check5.log","w",stdout);

      printf("\n");
      printf("Check of the program set_bstap()\n");
      printf("--------------------------------\n\n");

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

   start_ranlux(0,89103);
   geometry();

   print_flags();

   random_ud();
   set_bstap();

   print_flags();

   udb=udfld();
   hdb=bstap();
   set_ofs();
   set_psum0();
   set_psum1();
   check_psums();

   ie=check_bc(0.0);
   error_root(ie==0,1,"main [check5.c]","Boundary conditions changed");

   if (my_rank==0)
   {
      for (ifc=0;ifc<8;ifc++)
      {
         if (nfc[ifc]>0)
         {
            printf("ifc = %d, max|sum| = %.4e, maximal deviation = %.1e\n",
                   ifc,psum0[ifc],psum1[ifc]);
         }
      }

      printf("\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
