
/*******************************************************************************
*
* File bstap.c
*
* Copyright (C) 2012, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Allocation and computation of the boundary staple field.
*
* The externally accessible functions are
*
*   su3_dble *bstap(void)
*     Returns the base address of the boundary staple field. If it is
*     not already allocated, the field is allocated and initialized to
*     unity.
*
*   void set_bstap(void)
*     Computes the boundary staples and copies them to the neighbouring
*     MPI processes (see doc/gauge_actions.pdf).
*
* Notes:
*
* The boundary staple field has size 3*BNDRY and is logically divided into
* face segments. For the face with index ifc, the associated segment is
* at offset ofs[ifc] from the base address, where
*
*   ofs[0]=0
*   ofs[1]=ofs[0]+3*FACE0
*   ofs[2]=ofs[1]+3*FACE0
*   ofs[3]=ofs[2]+3*FACE1
*   ofs[4]=ofs[3]+3*FACE1
*   ofs[5]=ofs[4]+3*FACE2
*   ofs[6]=ofs[5]+3*FACE2
*   ofs[7]=ofs[6]+3*FACE3
*
* The ordering of the staples along the faces coincides with the ordering
* of the lattice points at the boundary (see main/README.global and also
* lattice/README.uidx for some further details). With open or SF boundary
* conditions, the program set_bstap() sets the time-like boundary staples
* at the boundaries of the lattice to zero.
*
* All these programs act globally and must be called on all MPI processes
* simultaneously.
*
*******************************************************************************/

#define BSTAP_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3fcts.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "global.h"

static const int plns[6][2]={{0,1},{0,2},{0,3},{2,3},{3,1},{1,2}};
static int bc,np,nfc[8],ofs[8],hofs[8],tags[8],nmu[8];
static su3_dble wd ALIGNED16;
static su3_dble *hdb=NULL;


static void set_ofs(void)
{
   int ifc;

   bc=bc_type();
   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;

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

   for (ifc=0;ifc<8;ifc++)
   {
      nmu[ifc]=cpr[ifc/2]&0x1;
      tags[ifc]=mpi_permanent_tag();
   }
}


static void alloc_hdb(void)
{
   int ifc,n;

   error(iup[0][0]==0,1,"alloc_hdb [bstap.c]",
         "Geometry arrays are not set");

   set_ofs();
   n=0;

   for (ifc=0;ifc<8;ifc+=2)
   {
      if (n<nfc[ifc])
         n=nfc[ifc];
   }

   n=3*(BNDRY+2*n);
   hdb=amalloc(n*sizeof(*hdb),ALIGN);
   error(hdb==NULL,1,"alloc_hdb [bstap.c]",
         "Unable to allocate the boundary staple field");
   cm3x3_unity(n,hdb);
}


su3_dble *bstap(void)
{
   if ((NPROC>1)&&(hdb==NULL))
      alloc_hdb();

   return hdb;
}


static void get_ofs(int mu,int nu,int ix,int *ip)
{
   int n,is;

   for (n=0;n<6;n++)
   {
      if (((plns[n][0]==mu)&&(plns[n][1]==nu))||
          ((plns[n][0]==nu)&&(plns[n][1]==mu)))
      {
         plaq_uidx(n,ix,ip);

         if (mu==plns[n][0])
         {
            is=ip[0];
            ip[0]=ip[2];
            ip[2]=is;

            is=ip[1];
            ip[1]=ip[3];
            ip[3]=is;
         }

         return;
      }
   }
}


static void get_staples(int ifc)
{
   int ib,ix,mu,nu,k,ip[4];
   su3_dble *udb,*sbuf;

   udb=udfld();
   sbuf=hdb+3*BNDRY;
   mu=ifc/2;

   for (ib=0;ib<(2*nfc[ifc]);ib++)
   {
      if (ib<(nfc[ifc]))
         ix=map[ofs[ifc]+ib];
      else
         ix=map[(BNDRY/2)+ofs[ifc]+ib-nfc[ifc]];

      for (k=0;k<3;k++)
      {
         nu=k+(k>=mu);
         get_ofs(mu,nu,ix,ip);

         if (ifc&0x1)
         {
            if ((mu>0)||(cpr[0]>0)||(bc==3))
            {
               su3xsu3dag(udb+ip[3],udb+ip[1],&wd);
               su3xsu3(udb+ip[2],&wd,sbuf);
            }
         }
         else
         {
            if ((mu>0)||(cpr[0]<(NPROC0-1))||(bc==3))
            {
               su3xsu3(udb+ip[0],udb+ip[1],&wd);
               su3dagxsu3(udb+ip[2],&wd,sbuf);
            }
         }

         sbuf+=1;
      }
   }
}


static void send_staples(int ifc,int tag)
{
   int saddr,raddr,nbf;
   su3_dble *sbuf,*rbuf;
   MPI_Status stat;

   saddr=npr[ifc^0x1];
   raddr=saddr;
   sbuf=hdb+3*BNDRY;
   rbuf=hdb+hofs[ifc^0x1];
   nbf=108*nfc[ifc];

   if ((ifc>1)||(bc==3)||
       ((ifc==1)&&(cpr[0]>0))||((ifc==0)&&(cpr[0]<(NPROC0-1))))
   {
      if (np==0)
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
   else
      cm3x3_zero(3*FACE0,rbuf);
}


void set_bstap(void)
{
   int ifc,sfc;

   if (query_flags(UDBUF_UP2DATE)!=1)
      copy_bnd_ud();

   if (NPROC>1)
   {
      if (hdb==NULL)
         alloc_hdb();

      for (ifc=0;ifc<8;ifc++)
      {
         sfc=ifc^nmu[ifc];

         if (nfc[sfc]>0)
         {
            get_staples(sfc);
            send_staples(sfc,tags[ifc]);
         }
      }
   }

   set_flags(SET_BSTAP);
}
