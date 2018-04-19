
/*******************************************************************************
*
* File sdcom.c
*
* Copyright (C) 2005, 2008, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Communication functions for double-precision spinor fields.
*
*   void cpsd_int_bnd(int is,spinor_dble *sd)
*     Copies the spinors sd at the even interior boundary points of the
*     local lattice to the corresponding points on the neighbouring MPI
*     processes. Only half of the spinor components are copied, namely
*     theta[ifc^(is&0x1)]*sd, where ifc labels the faces of the local
*     lattice on the sending process.
*
*   void cpsd_ext_bnd(int is,spinor_dble *sd)
*     Copies the spinors sd at the even exterior boundary points of the
*     local lattice to the neighbouring MPI processes and *adds* them to
*     the field on the matching points of the target lattices. Only half
*     of the spinor components are copied, assuming the spinors sd satisfy
*     sd=theta[ifc^(is&0x1)]*sd, where ifc labels the faces of the local
*     lattice on the sending process.
*
* Notes:
*
* The spinor fields passed to cpsd_int_bnd() and cpsd_ext_bnd() must have at
* least NSPIN elements. They are interpreted as quark fields on the local
* lattice as described in main/README.global and doc/dirac.pdf. The projector
* theta[ifc] is defined at the top of the module sflds/Pbnd.c.
*
* If open, SF or open-SF boundary conditions are chosen, the programs do
* not copy any spinors in the time direction across the boundaries of the
* global lattice. The spinors on the even points at the boundaries are instead
* set zero as required by the boundary conditions (see doc/dirac.pdf). More
* precisely, cpsd_int_bnd() sets them to zero *before* copying any spinors in
* the space directions, while cpsd_ext_bnd() does the opposite.
*
* All these programs involve global communications and must be called on
* all processes simultaneously.
*
*******************************************************************************/

#define SDCOM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "sflds.h"
#include "global.h"

static int bc,np,nmu[8],nbf[8],ofs[8];
static int ns,sfc[8],rfc[8],sflg[8];
static int itags=0,tags[8];
static weyl_dble *wb=NULL,*snd_buf[8],*rcv_buf[8];
static const weyl_dble w0={{{0.0}}};
static MPI_Request snd_req[8],rcv_req[8];


static void get_tags(void)
{
   int i;

   if (itags==0)
   {
      for (i=0;i<8;i++)
         tags[i]=mpi_permanent_tag();

      itags=1;
   }
}


static void alloc_sdbufs(void)
{
   int n,ifc,tag,saddr,raddr;
   weyl_dble *w,*wm;

   error(iup[0][0]==0,1,"alloc_sdbufs [sdcom.c]",
         "Geometry arrays are not initialized");

   wb=amalloc(BNDRY*sizeof(*wb),ALIGN);
   error(wb==NULL,1,"alloc_sdbufs [sdcom.c]",
         "Unable to allocate communication buffers");

   w=wb;
   wm=wb+BNDRY;

   for (;w<wm;w++)
      (*w)=w0;

   bc=bc_type();
   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;

   nbf[0]=FACE0/2;
   nbf[1]=FACE0/2;
   nbf[2]=FACE1/2;
   nbf[3]=FACE1/2;
   nbf[4]=FACE2/2;
   nbf[5]=FACE2/2;
   nbf[6]=FACE3/2;
   nbf[7]=FACE3/2;

   get_tags();
   ofs[0]=0;
   ns=0;
   w=wb;

   for (ifc=0;ifc<8;ifc++)
   {
      nmu[ifc]=cpr[ifc/2]&0x1;

      if (ifc>0)
         ofs[ifc]=ofs[ifc-1]+nbf[ifc-1];

      if (nbf[ifc]>0)
      {
         sfc[ns]=ifc;
         ns+=1;

         snd_buf[ifc]=w;
         w+=nbf[ifc];
         rcv_buf[ifc]=w;
         w+=nbf[ifc];

         tag=tags[ifc];
         saddr=npr[ifc];
         raddr=npr[ifc^0x1];

         MPI_Send_init((double*)(snd_buf[ifc]),12*nbf[ifc],MPI_DOUBLE,saddr,
                       tag,MPI_COMM_WORLD,&snd_req[ifc]);
         MPI_Recv_init((double*)(rcv_buf[ifc]),12*nbf[ifc],MPI_DOUBLE,raddr,
                       tag,MPI_COMM_WORLD,&rcv_req[ifc]);
      }

      sflg[ifc]=((ifc>1)||
                 ((ifc==0)&&(cpr[0]!=0))||
                 ((ifc==1)&&(cpr[0]!=(NPROC0-1)))||
                 (bc==3));
   }

   for (n=0;n<ns;n+=2)
   {
      rfc[n]=sfc[ns-n-2];
      rfc[n+1]=sfc[ns-n-1];
   }
}

#if (defined x64)
#include "sse2.h"

static void zip_weyl(int vol,spinor_dble *pk,weyl_dble *pl)
{
   weyl_dble *pm;

   pm=pl+vol;

   for (;pl<pm;pl++)
   {
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5"
                            :
                            :
                            "m" ((*pk).c1.c1),
                            "m" ((*pk).c1.c2),
                            "m" ((*pk).c1.c3),
                            "m" ((*pk).c2.c1),
                            "m" ((*pk).c2.c2),
                            "m" ((*pk).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"
                            "movapd %%xmm5, %5"
                            :
                            "=m" ((*pl).c1.c1),
                            "=m" ((*pl).c1.c2),
                            "=m" ((*pl).c1.c3),
                            "=m" ((*pl).c2.c1),
                            "=m" ((*pl).c2.c2),
                            "=m" ((*pl).c2.c3));

      pk+=1;
   }
}


static void unzip_weyl(int vol,weyl_dble *pk,spinor_dble *pl)
{
   spinor_dble *pm;

   __asm__ __volatile__ ("xorpd %%xmm6, %%xmm6 \n\t"
                         "xorpd %%xmm7, %%xmm7"
                         :
                         :
                         :
                         "xmm6", "xmm7");

   pm=pl+vol;

   for (;pl<pm;pl++)
   {
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5"
                            :
                            :
                            "m" ((*pk).c1.c1),
                            "m" ((*pk).c1.c2),
                            "m" ((*pk).c1.c3),
                            "m" ((*pk).c2.c1),
                            "m" ((*pk).c2.c2),
                            "m" ((*pk).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      __asm__ __volatile__ ("addpd %%xmm0, %%xmm0 \n\t"
                            "addpd %%xmm1, %%xmm1 \n\t"
                            "addpd %%xmm2, %%xmm2 \n\t"
                            "addpd %%xmm3, %%xmm3 \n\t"
                            "addpd %%xmm4, %%xmm4 \n\t"
                            "addpd %%xmm5, %%xmm5"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"
                            "movapd %%xmm5, %5"
                            :
                            "=m" ((*pl).c1.c1),
                            "=m" ((*pl).c1.c2),
                            "=m" ((*pl).c1.c3),
                            "=m" ((*pl).c2.c1),
                            "=m" ((*pl).c2.c2),
                            "=m" ((*pl).c2.c3));

      __asm__ __volatile__ ("movapd %%xmm6, %0 \n\t"
                            "movapd %%xmm7, %1 \n\t"
                            "movapd %%xmm6, %2 \n\t"
                            "movapd %%xmm7, %3 \n\t"
                            "movapd %%xmm6, %4 \n\t"
                            "movapd %%xmm7, %5"
                            :
                            "=m" ((*pl).c3.c1),
                            "=m" ((*pl).c3.c2),
                            "=m" ((*pl).c3.c3),
                            "=m" ((*pl).c4.c1),
                            "=m" ((*pl).c4.c2),
                            "=m" ((*pl).c4.c3));

      pk+=1;
   }
}

#else

static void zip_weyl(int vol,spinor_dble *pk,weyl_dble *pl)
{
   weyl_dble *pm;

   pm=pl+vol;

   for (;pl<pm;pl++)
   {
      (*pl).c1=(*pk).c1;
      (*pl).c2=(*pk).c2;

      pk+=1;
   }
}


static void unzip_weyl(int vol,weyl_dble *pk,spinor_dble *pl)
{
   weyl_dble *pm;

   pm=pk+vol;

   for (;pk<pm;pk++)
   {
      _vector_add((*pl).c1,(*pk).c1,(*pk).c1);
      _vector_add((*pl).c2,(*pk).c2,(*pk).c2);
      (*pl).c3=w0.c1;
      (*pl).c4=w0.c1;

      pl+=1;
   }
}

#endif

static void send_bufs(int ifc,int eo)
{
   int io;

   io=(ifc^nmu[ifc]);

   if (sflg[io])
   {
      if (np==eo)
         MPI_Start(&snd_req[io]);
      else
         MPI_Start(&rcv_req[io^0x1]);
   }
}


static void wait_bufs(int ifc,int eo)
{
   int io;
   MPI_Status stat_snd,stat_rcv;

   io=(ifc^nmu[ifc]);

   if (sflg[io])
   {
      if (np==eo)
         MPI_Wait(&snd_req[io],&stat_snd);
      else
         MPI_Wait(&rcv_req[io^0x1],&stat_rcv);
   }
}


void cpsd_int_bnd(int is,spinor_dble *sd)
{
   int ifc,io;
   int n,m,eo;
   spinor_dble *sdb;

   if (NPROC0==1)
      bnd_sd2zero(EVEN_PTS,sd);

   if (NPROC==1)
      return;
   else if (wb==NULL)
   {
      alloc_sdbufs();
      bnd_sd2zero(NO_PTS,sd);
   }

   is&=0x1;
   m=0;
   eo=0;
   sdb=sd+VOLUME;

   for (n=0;n<ns;n++)
   {
      if (n>0)
         send_bufs(sfc[m],eo);

      ifc=sfc[n];
      io=ifc^nmu[ifc];

      if (sflg[io])
         assign_sd2wd[io^is](map+ofs[io^0x1],nbf[io],sd,snd_buf[io]);
      else
         bnd_sd2zero(EVEN_PTS,sd);

      if (n>0)
      {
         wait_bufs(sfc[m],eo);
         m+=eo;
         eo^=0x1;
      }
   }

   for (n=0;n<2;n++)
   {
      send_bufs(sfc[m],eo);
      wait_bufs(sfc[m],eo);
      m+=eo;
      eo^=0x1;
   }

   for (n=0;n<ns;n++)
   {
      if (m<ns)
         send_bufs(sfc[m],eo);

      ifc=sfc[n];
      io=(ifc^nmu[ifc])^0x1;

      if (sflg[io^0x1])
         unzip_weyl(nbf[io],rcv_buf[io],sdb+ofs[io^0x1]);
      else if ((io==0)&&((bc==1)||(bc==2)))
         set_sd2zero(nbf[io],sdb+ofs[io^0x1]);

      if (m<ns)
      {
         wait_bufs(sfc[m],eo);
         m+=eo;
         eo^=0x1;
      }
   }
}


void cpsd_ext_bnd(int is,spinor_dble *sd)
{
   int ifc,io;
   int n,m,eo;
   spinor_dble *sdb;

   if (NPROC==1)
   {
      bnd_sd2zero(EVEN_PTS,sd);
      return;
   }
   else if (wb==NULL)
   {
      alloc_sdbufs();
      bnd_sd2zero(NO_PTS,sd);
   }

   is&=0x1;
   m=0;
   eo=0;
   sdb=sd+VOLUME;

   for (n=0;n<ns;n++)
   {
      if (n>0)
         send_bufs(rfc[m],eo);

      ifc=rfc[n];
      io=ifc^nmu[ifc];

      if (sflg[io])
         zip_weyl(nbf[io],sdb+ofs[io],snd_buf[io]);

      if (n>0)
      {
         wait_bufs(rfc[m],eo);
         m+=eo;
         eo^=0x1;
      }
   }

   for (n=0;n<2;n++)
   {
      send_bufs(rfc[m],eo);
      wait_bufs(rfc[m],eo);
      m+=eo;
      eo^=0x1;
   }

   for (n=0;n<ns;n++)
   {
      if (m<ns)
         send_bufs(rfc[m],eo);

      ifc=rfc[n];
      io=(ifc^nmu[ifc])^0x1;

      if (sflg[io^0x1])
         add_assign_wd2sd[io^is](map+ofs[io],nbf[io],rcv_buf[io],sd);
      else
         bnd_sd2zero(EVEN_PTS,sd);

      if (m<ns)
      {
         wait_bufs(rfc[m],eo);
         m+=eo;
         eo^=0x1;
      }
   }

   if (NPROC0==1)
      bnd_sd2zero(EVEN_PTS,sd);
}
