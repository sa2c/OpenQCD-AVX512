
/*******************************************************************************
*
* File sap_com.c
*
* Copyright (C) 2005, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* SAP communication program.
*
* The externally accessible functions are
*
*   void alloc_sap_bufs(void)
*     Allocates and initializes the buffers and index arrays needed for
*     the program sap_com().
*
*   void sap_com(int ic,spinor *r)
*     Subtracts the Weyl field b.bb.w[0] on the boundaries of all black
*     (if ic=0) or all white (if ic=1) blocks b of the SAP_BLOCKS grid
*     from the global spinor field r. Before subtraction, the Weyl fields
*     on the block faces in direction ifc are expanded to Dirac spinor
*     fields s satisfying theta[ifc]*s=0.
*
* Notes:
*
* The program alloc_sap_bufs() adds a single-precision Weyl field to the
* boundaries of the blocks in the SAP_BLOCKS grid. In memory these fields
* are arranged in particular way, but they are field arrays exactly like
* the ones on any block created by the allocation programs in the module
* block/block.c.
*
* alloc_sap_bufs() is called when the SAP_BLOCKS block grid is allocated.
* This program is not intended to be called from anywhere else and does
* nothing if called a second time.
*
* The operations performed by the program sap_com() are explained in some
* detail in README.sap_com. In the case of boundary conditions of type 0,
* 1 or 2, the Weyl fields residing at the exterior boundaries of the blocks
* at global time -1 and NPROC0*L0 are not subtracted from the field r.
*
* All programs in this module may involve communications and must be called
* simultaneously on all processes.
*
*******************************************************************************/

#define SAP_COM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "sflds.h"
#include "block.h"
#include "sap.h"
#include "global.h"

static int nb,nbh,isw,init=0;
static int bc,np,nmu[8],sflg[8];
static int nsbf[2][8],nlbf[2][8],*imb[2][8];
static weyl *snd_buf[2][8],*loc_buf[2][8],*rcv_buf[2][8];
static const weyl w0={{{0.0f}}};
static block_t *b0;
static MPI_Request snd_req[2][8],rcv_req[2][8];


static void set_nbf(void)
{
   int ifc,ibu,ibd;
   int *bo,*bs;
   block_t *b,*bm;
   bndry_t *bb;

   bc=bc_type();
   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;

   bs=(*b0).bs;
   ibu=((cpr[0]==(NPROC0-1))&&(bc!=3));
   ibd=((cpr[0]==0)&&(bc!=3));

   for (ifc=0;ifc<8;ifc++)
   {
      nmu[ifc]=cpr[ifc/2]&0x1;
      sflg[ifc]=((ifc>1)||
                 ((ifc==0)&&(cpr[0]!=0))||
                 ((ifc==1)&&(cpr[0]!=(NPROC0-1)))||
                 (bc==3));

      nlbf[0][ifc]=0;
      nsbf[0][ifc]=0;
      nlbf[1][ifc]=0;
      nsbf[1][ifc]=0;

      b=b0;
      bm=b+nbh;

      for (;b<bm;b++)
      {
         bo=(*b).bo;
         bb=(*b).bb;

         if ((bb[ifc].ibn)||
             ((ifc==0)&&(ibd)&&(bo[0]==0))||
             ((ifc==1)&&(ibu)&&((bo[0]+bs[0])==L0)))
            nsbf[isw][ifc]+=bb[ifc].vol;
         else
            nlbf[isw][ifc]+=bb[ifc].vol;
      }

      bm+=nbh;

      for (;b<bm;b++)
      {
         bo=(*b).bo;
         bb=(*b).bb;

         if ((bb[ifc].ibn)||
             ((ifc==0)&&(ibd)&&(bo[0]==0))||
             ((ifc==1)&&(ibu)&&((bo[0]+bs[0])==L0)))
            nsbf[isw^0x1][ifc]+=bb[ifc].vol;
         else
            nlbf[isw^0x1][ifc]+=bb[ifc].vol;
      }
   }
}


static void alloc_weyl(void)
{
   int n,ic,ifc;
   weyl *w,*wm;

   n=0;

   for (ifc=0;ifc<8;ifc++)
   {
      n+=(nlbf[0][ifc]+2*nsbf[0][ifc]);
      n+=(nlbf[1][ifc]+2*nsbf[1][ifc]);
   }

   w=amalloc(n*sizeof(*w),ALIGN);
   error(w==NULL,1,"alloc_weyl [sap_com.c]","Unable to allocate buffers");

   for (ic=0;ic<2;ic++)
   {
      for (ifc=0;ifc<8;ifc++)
      {
         snd_buf[ic][ifc]=w;
         w+=nsbf[ic][ifc];
         loc_buf[ic][ifc]=w;
         w+=nlbf[ic][ifc];
         rcv_buf[ic][ifc]=w;
         w+=nsbf[ic][ifc];
      }
   }

   w=snd_buf[0][0];
   wm=w+n;

   for (;w<wm;w++)
      (*w)=w0;
}


static void add_weyl(void)
{
   int ifc,ibd,ibu;
   int *bo,*bs;
   weyl **pw,*ws,*wl;
   block_t *b,*bm;
   bndry_t *bb;

   pw=malloc(8*nb*sizeof(*pw));
   error(pw==NULL,1,"add_weyl [sap_com.c]",
         "Unable to add the Weyl fields to the block boundaries");

   bs=(*b0).bs;
   ibd=((cpr[0]==0)&&(bc!=3));
   ibu=((cpr[0]==(NPROC0-1))&&(bc!=3));

   for (ifc=0;ifc<8;ifc++)
   {
      b=b0;
      bm=b+nbh;
      ws=snd_buf[isw][ifc];
      wl=loc_buf[isw][ifc];

      for (;b<bm;b++)
      {
         bo=(*b).bo;
         bb=(*b).bb;
         bb[ifc].nw=1;
         bb[ifc].w=pw;

         if ((bb[ifc].ibn)||
             ((ifc==0)&&(ibd)&&(bo[0]==0))||
             ((ifc==1)&&(ibu)&&((bo[0]+bs[0])==L0)))
         {
            (*pw)=ws;
            ws+=bb[ifc].vol;
         }
         else
         {
            (*pw)=wl;
            wl+=bb[ifc].vol;
         }

         pw+=1;
      }

      bm+=nbh;
      ws=snd_buf[isw^0x1][ifc];
      wl=loc_buf[isw^0x1][ifc];

      for (;b<bm;b++)
      {
         bo=(*b).bo;
         bb=(*b).bb;
         bb[ifc].nw=1;
         bb[ifc].w=pw;

         if ((bb[ifc].ibn)||
             ((ifc==0)&&(ibd)&&(bo[0]==0))||
             ((ifc==1)&&(ibu)&&((bo[0]+bs[0])==L0)))
         {
            (*pw)=ws;
            ws+=bb[ifc].vol;
         }
         else
         {
            (*pw)=wl;
            wl+=bb[ifc].vol;
         }

         pw+=1;
      }
   }
}


static void set_mpi_req(void)
{
   int ic,ifc;
   int saddr,raddr,tag,nbf;

   for (ic=0;ic<2;ic++)
   {
      for (ifc=0;ifc<8;ifc++)
      {
         saddr=npr[ifc];
         raddr=npr[ifc^0x1];
         tag=mpi_permanent_tag();
         nbf=12*nsbf[ic][ifc];

         MPI_Send_init(snd_buf[ic][ifc],nbf,MPI_FLOAT,
                       saddr,tag,MPI_COMM_WORLD,&snd_req[ic][ifc]);
         MPI_Recv_init(rcv_buf[ic][ifc],nbf,MPI_FLOAT,
                       raddr,tag,MPI_COMM_WORLD,&rcv_req[ic][ifc]);
      }
   }
}


static void alloc_imb(void)
{
   int n,ic,ifc,ibd,ibu;
   int *bo,*bs,*im;
   block_t *b,*bm;
   bndry_t *bb;

   n=0;

   for (ifc=0;ifc<8;ifc++)
   {
      n+=(nlbf[0][ifc]+nsbf[0][ifc]);
      n+=(nlbf[1][ifc]+nsbf[1][ifc]);
   }

   im=malloc(n*sizeof(*im));
   error(im==NULL,1,"alloc_imb [sap_com.c]",
         "Unable to allocate index arrays");

   bs=(*b0).bs;
   ibd=((cpr[0]==0)&&(bc!=3));
   ibu=((cpr[0]==(NPROC0-1))&&(bc!=3));

   for (ic=0;ic<2;ic++)
   {
      for (ifc=0;ifc<8;ifc++)
      {
         imb[ic][ifc]=im;

         if (ic^isw)
            b=b0+nbh;
         else
            b=b0;
         bm=b+nbh;

         for (;b<bm;b++)
         {
            bo=(*b).bo;
            bb=(*b).bb;

            if (!((bb[ifc].ibn)||
                  ((ifc==0)&&(ibd)&&(bo[0]==0))||
                  ((ifc==1)&&(ibu)&&((bo[0]+bs[0])==L0))))
            {
               for (n=0;n<bb[ifc].vol;n++)
                  im[n]=bb[ifc].imb[n];

               im+=bb[ifc].vol;
            }
         }

         if (ic^isw)
            b=b0;
         else
            b=b0+nbh;
         bm=b+nbh;

         for (;b<bm;b++)
         {
            bo=(*b).bo;
            bb=(*b).bb;

            if ((bb[ifc^0x1].ibn)||
                ((ifc==1)&&(ibd)&&(bo[0]==0))||
                ((ifc==0)&&(ibu)&&((bo[0]+bs[0])==L0)))
            {
               for (n=0;n<bb[ifc].vol;n++)
                  im[n]=(*b).imb[bb[ifc].map[n]];

               im+=bb[ifc].vol;
            }
         }
      }
   }
}


void alloc_sap_bufs(void)
{
   bndry_t *bb;

   if (init==1)
      return;

   b0=blk_list(SAP_BLOCKS,&nb,&isw);
   error(b0==NULL,1,"alloc_sap_bufs [sap_com.c]",
         "Block grid is not allocated");

   bb=(*b0).bb;
   error((bb==NULL)||((*bb).nw!=0),1,"alloc_sap_bufs [sap_com.c]",
         "Block boundary is not in the proper condition");

   nbh=nb/2;
   set_nbf();
   alloc_weyl();
   add_weyl();
   set_mpi_req();
   alloc_imb();

   init=1;
}


static void send_buf(int ic,int ifc,int eo)
{
   int io;

   io=ifc^nmu[ifc];

   if (sflg[io])
   {
      if (np==eo)
      {
         if (nsbf[ic][io])
            MPI_Start(&snd_req[ic][io]);
      }
      else
      {
         if (nsbf[ic][io^0x1])
            MPI_Start(&rcv_req[ic][io^0x1]);
      }
   }
}


static void wait_buf(int ic,int ifc,int eo)
{
   int io;
   MPI_Status stat;

   io=ifc^nmu[ifc];

   if (sflg[io])
   {
      if (np==eo)
      {
         if (nsbf[ic][io])
            MPI_Wait(&snd_req[ic][io],&stat);
      }
      else
      {
         if (nsbf[ic][io^0x1])
            MPI_Wait(&rcv_req[ic][io^0x1],&stat);
      }
   }
}


void sap_com(int ic,spinor *r)
{
   int ifc,io,nbf;

   if (init==0)
   {
      error_root(1,1,"sap_com [sap_com.c]",
                 "Communication buffers are not allocated");
      return;
   }

   send_buf(ic,0,0);
   send_buf(ic,0,1);

   for (ifc=0;ifc<8;ifc++)
   {
      wait_buf(ic,ifc,0);
      wait_buf(ic,ifc,1);

      if (ifc<7)
      {
         send_buf(ic,ifc+1,0);
         send_buf(ic,ifc+1,1);
      }

      io=(ifc^nmu[ifc])^0x1;

      if (sflg[io^0x1])
         nbf=nlbf[ic][io]+nsbf[ic][io];
      else
         nbf=nlbf[ic][io];

      sub_assign_w2s[io^0x1](imb[ic][io],nbf,loc_buf[ic][io],r);
   }
}
