
/*******************************************************************************
*
* File ftcom.c
*
* Copyright (C) 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Communication of the field tensor components residing at the boundaries
* of the local lattices.
*
* The externally accessible functions are
*
*   void copy_bnd_ft(int n,u3_alg_dble *ft)
*     Fetches the boundary values the field ft from the neighbouring MPI
*     processes (see the notes). The boundary values at time NPROC0*L0
*     are fetched from the field at time 0 only in the case of periodic
*     boundary conditions.
*
*   void add_bnd_ft(int n,u3_alg_dble *ft)
*     Adds the boundary values of the field ft to the field on the
*     neighbouring MPI processes. The boundary values at time NPROC0*L0
*     are added to the field at time 0 only in the case of periodic
*     boundary conditions.
*
* Notes:
*
* Both communication programs assume that the field ft has the same size as
* the n-th component of the symmetric field tensor F_{mu nu}, where n=0,..,5
* labels the (mu,nu)-planes (0,1),(0,2),(0,3),(2,3),(3,1),(1,2). For further
* explanations, see the files lattice/README.ftidx and tcharge/ftensor.c.
*
* The programs in this module perform global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define FTCOM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "tcharge.h"
#include "global.h"

static const int plns[6][2]={{0,1},{0,2},{0,3},{2,3},{3,1},{1,2}};
static u3_alg_dble *ftbuf;
static ftidx_t *idx=NULL;


static void alloc_ftbuf(void)
{
   int n,nft,nbf;

   idx=ftidx();
   nbf=0;

   for (n=0;n<6;n++)
   {
      nft=idx[n].nft[0];
      if (nft>nbf)
         nbf=nft;

      nft=idx[n].nft[1];
      if (nft>nbf)
         nbf=nft;
   }

   ftbuf=amalloc(nbf*sizeof(*ftbuf),ALIGN);
   error(ftbuf==NULL,1,"alloc_ftbuf [ftcom.c]",
         "Unable to allocate communication buffers");
}


static void pack_buf(int n,int dir,u3_alg_dble *ft)
{
   int bc,mu,nft;
   int *ift,*ifm;
   u3_alg_dble *fb;

   nft=idx[n].nft[dir];

   if (nft>0)
   {
      bc=bc_type();
      mu=plns[n][dir];

      if ((mu>0)||(cpr[0]>0)||(bc==3))
      {
         ift=idx[n].ift[dir];
         ifm=ift+nft;
         fb=ftbuf;

         for (;ift<ifm;ift++)
         {
            fb[0]=ft[*ift];
            fb+=1;
         }
      }
   }
}


static void fwd_send(int n,int dir,u3_alg_dble *ft)
{
   int bc,mu,nft,nbf;
   int tag,saddr,raddr,np;
   u3_alg_dble *sbuf,*rbuf;
   MPI_Status stat;

   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   nft=idx[n].nft[dir];

   if (nft>0)
   {
      bc=bc_type();
      mu=plns[n][dir];
      tag=mpi_tag();
      saddr=npr[2*mu];
      raddr=npr[2*mu+1];
      sbuf=ftbuf;
      rbuf=ft+VOLUME;
      if (dir==1)
         rbuf+=idx[n].nft[0];
      nbf=9*nft;

      if (np==0)
      {
         if ((mu>0)||(cpr[0]>0)||(bc==3))
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         if ((mu>0)||(cpr[0]<(NPROC0-1))||(bc==3))
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
      }
      else
      {
         if ((mu>0)||(cpr[0]<(NPROC0-1))||(bc==3))
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         if ((mu>0)||(cpr[0]>0)||(bc==3))
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
      }
   }
}


void copy_bnd_ft(int n,u3_alg_dble *ft)
{
   if (NPROC>1)
   {
      if (idx==NULL)
         alloc_ftbuf();

      pack_buf(n,1,ft);
      fwd_send(n,1,ft);
      pack_buf(n,0,ft);
      fwd_send(n,0,ft);
   }
}


static void bck_send(int n,int dir,u3_alg_dble *ft)
{
   int bc,mu,nft,nbf;
   int tag,saddr,raddr,np;
   u3_alg_dble *sbuf,*rbuf;
   MPI_Status stat;

   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   nft=idx[n].nft[dir];

   if (nft>0)
   {
      bc=bc_type();
      mu=plns[n][dir];
      tag=mpi_tag();
      saddr=npr[2*mu+1];
      raddr=npr[2*mu];
      sbuf=ft+VOLUME;
      if (dir==1)
         sbuf+=idx[n].nft[0];
      rbuf=ftbuf;
      nbf=9*nft;

      if (np==0)
      {
         if ((mu>0)||(cpr[0]<(NPROC0-1))||(bc==3))
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         if ((mu>0)||(cpr[0]>0)||(bc==3))
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
      }
      else
      {
         if ((mu>0)||(cpr[0]>0)||(bc==3))
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         if ((mu>0)||(cpr[0]<(NPROC0-1))||(bc==3))
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
      }
   }
}


static void unpack_buf(int n,int dir,u3_alg_dble *ft)
{
   int bc,mu,nft;
   int *ift,*ifm;
   u3_alg_dble *f,*fb;

   nft=idx[n].nft[dir];

   if (nft>0)
   {
      bc=bc_type();
      mu=plns[n][dir];

      if ((mu>0)||(cpr[0]>0)||(bc==3))
      {
         ift=idx[n].ift[dir];
         ifm=ift+nft;
         fb=ftbuf;

         for (;ift<ifm;ift++)
         {
            f=ft+(*ift);

            (*f).c1+=(*fb).c1;
            (*f).c2+=(*fb).c2;
            (*f).c3+=(*fb).c3;
            (*f).c4+=(*fb).c4;
            (*f).c5+=(*fb).c5;
            (*f).c6+=(*fb).c6;
            (*f).c7+=(*fb).c7;
            (*f).c8+=(*fb).c8;
            (*f).c9+=(*fb).c9;

            fb+=1;
         }
      }
   }
}


void add_bnd_ft(int n,u3_alg_dble *ft)
{
   if (NPROC>1)
   {
      if (idx==NULL)
         alloc_ftbuf();

      bck_send(n,0,ft);
      unpack_buf(n,0,ft);
      bck_send(n,1,ft);
      unpack_buf(n,1,ft);
   }
}
