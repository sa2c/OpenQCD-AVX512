
/*******************************************************************************
*
* File Aw_com.c
*
* Copyright (C) 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Communication functions needed for the computation of the little Dirac
* operator.
*
*   b2b_flds_t *b2b_flds(int n,int mu)
*     Extracts the spinor fields on the interior boundaries of the n'th
*     block of the DFL_BLOCKS grid and its neighbouring block in direction
*     mu. The spinors on the odd sites are multiplied by the link variables
*     in direction mu and -mu respectively. If the two blocks touch the
*     boundary of the local lattice, the fields extracted from the even
*     sites are copied to the neighbouring process. The program returns a
*     structure containing the extracted field arrays (see README.Aw_com
*     for detailed explanations).
*
*   void cpAoe_ext_bnd(void)
*     Copies the hopping terms Aoe and Aeo of the double-precision little
*     Dirac operator on the odd exterior boundary points of the local block
*     lattice to the neighbouring MPI processes and *adds* them to the hop-
*     ping terms on the matching blocks on the target lattices.
*
*   void cpAee_int_bnd(void)
*     Copies the even-even terms Aee of the double-precision little Dirac
*     operator on the (even) interior boundary points of the local block
*     lattice to the neighbouring MPI processes.
*
* Notes:
*
* The program b2b_flds() writes the extracted spinor fields to internally
* allocated field arrays. These are reused when the program is called
* the next time. The data in the field arrays returned by b2b_flds() are
* therefore preserved only up to the next call of the program.
*
*******************************************************************************/

#define AW_COM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "uflds.h"
#include "sflds.h"
#include "vflds.h"
#include "dfl.h"
#include "little.h"
#include "global.h"

typedef struct
{
   int *iud[2];
   int *ise[2];
   int *iso[2];
   su3_dble *ud[2];
   spinor_dble **sd[2];
   spinor_dble **snd_buf[2];
   b2b_flds_t b2b;
} bsd_t;

static int init_bsd=0,init_Aoe=0,init_Aee=0;
static int np,nmu[8];
static int Ns=0,nb,nbb,nbh,nbbh;
static int nbbe[8],nbbo[8],obbe[8],obbo[8];
static int (*inn)[8],*idx,*ipp,*mp;
static int nsnd,sfc[8];

static complex_dble **snd_buf_Aee[8];
static complex_dble **rcv_buf_Aoe[8],**rcv_buf_Aeo[8];
static bsd_t (*bsd)[4];

static MPI_Request snd_req_bsd[8],rcv_req_bsd[8];
static MPI_Request snd_req_Aee[8],rcv_req_Aee[8];
static MPI_Request snd_req_Aoe[8],rcv_req_Aoe[8];
static MPI_Request snd_req_Aeo[8],rcv_req_Aeo[8];


static void set_constants(void)
{
   int ifc;
   dfl_parms_t dfl;
   dfl_grid_t grd;

   dfl=dfl_parms();
   grd=dfl_geometry();

   Ns=dfl.Ns;
   nb=grd.nb;
   nbb=grd.nbb;
   nbh=nb/2;
   nbbh=nbb/2;

   for (ifc=0;ifc<8;ifc++)
   {
      nbbe[ifc]=grd.nbbe[ifc];
      nbbo[ifc]=grd.nbbo[ifc];
      obbe[ifc]=grd.obbe[ifc];
      obbo[ifc]=grd.obbo[ifc];
   }

   inn=grd.inn;
   idx=grd.idx;
   ipp=grd.ipp;
   mp=grd.map;

   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   nsnd=0;

   for (ifc=0;ifc<8;ifc++)
   {
      nmu[ifc]=cpr[ifc/2]&0x1;

      if (nbbe[ifc]+nbbo[ifc])
      {
         sfc[nsnd]=ifc;
         nsnd+=1;
      }
   }
}


static int fnd_nn(int n,int ifc)
{
   n=idx[n];
   n=inn[n][ifc];

   if (n>=nb)
      n=mp[n-nb];

   return idx[n];
}


static void set_snd_req_bsd(void)
{
   int ifc,vol,nbf;
   int tag,saddr,raddr;
   bsd_t *brd;
   b2b_flds_t *b2b;

   for (ifc=0;ifc<8;ifc++)
   {
      brd=bsd[0]+(ifc/2);
      b2b=&(*brd).b2b;
      vol=(*b2b).vol;

      nbf=24*Ns*vol;
      saddr=npr[ifc];
      raddr=npr[ifc^0x1];
      tag=mpi_permanent_tag();

      MPI_Send_init((*brd).snd_buf[(ifc&0x1)^0x1][0],nbf,
                    MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD,&snd_req_bsd[ifc]);
      MPI_Recv_init((*b2b).sde[(ifc&0x1)^0x1][0],nbf,
                    MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&rcv_req_bsd[ifc]);
   }
}


static void alloc_bsd(void)
{
   int nbs,isw,vbb,vbm;
   int n,m,mu,ifc,vol,*iud;
   int ix,iy,k;
   su3_dble *ud;
   spinor_dble **psd,*sd;
   block_t *b;
   bndry_t *bb;
   bsd_t *brd;

   if (Ns==0)
      set_constants();

   b=blk_list(DFL_BLOCKS,&nbs,&isw);
   error(nbs==0,1,"alloc_bsd [Aw_com.c]",
         "DFL_BLOCKS grid is not allocated");

   bb=(*b).bb;
   vbb=0;
   vbm=0;

   for (mu=0;mu<4;mu++)
   {
      vol=bb[2*mu].vol;
      vbb+=vol;

      if (vol>vbm)
         vbm=vol;
   }

   bsd=malloc(nb*sizeof(*bsd));
   iud=malloc(nb*vbb*sizeof(*iud));
   ud=amalloc(vbm*sizeof(*ud),ALIGN);
   psd=malloc(24*Ns*sizeof(*psd));
   sd=amalloc(3*Ns*vbm*sizeof(*sd),ALIGN);

   error((bsd==NULL)||(iud==NULL)||(ud==NULL)||(psd==NULL)||(sd==NULL),1,
         "alloc_bsd [Aw_com.c]","Unable to allocate buffers");

   set_sd2zero(3*Ns*vbm,sd);

   for (n=0;n<nb;n++)
   {
      brd=bsd[n];

      for (mu=0;mu<4;mu++)
      {
         ifc=2*mu+1;
         m=fnd_nn(n,ifc);
         bb=b[n].bb;
         vol=bb[ifc].vol/2;

         (*brd).iud[0]=iud;

         for (ix=0;ix<vol;ix++)
         {
            iy=bb[ifc].ipp[ix];
            (*iud)=8*(b[n].imb[iy]-VOLUME/2)+(ifc^0x1);
            iud+=1;
         }

         (*brd).iud[1]=iud;

         for (ix=0;ix<vol;ix++)
         {
            iy=bb[ifc].map[ix+vol];
            (*iud)=8*(b[m].imb[iy]-VOLUME/2)+ifc;
            iud+=1;
         }

         (*brd).sd[0]=b[n].sd+1;
         (*brd).sd[1]=b[m].sd+1;
         (*brd).ise[0]=bb[ifc].ipp+vol;
         (*brd).iso[0]=bb[ifc].ipp;
         (*brd).ise[1]=bb[ifc].map;
         (*brd).iso[1]=bb[ifc].map+vol;

         (*brd).ud[0]=ud;
         (*brd).ud[1]=ud+vol;

         if (n==0)
         {
            (*brd).b2b.vol=vol;
            (*brd).b2b.sde[0]=psd;
            (*brd).b2b.sdo[1]=psd+Ns;
            (*brd).b2b.sdo[0]=psd+2*Ns;
            (*brd).b2b.sde[1]=psd+3*Ns;
            (*brd).snd_buf[0]=psd+4*Ns;
            (*brd).snd_buf[1]=psd+5*Ns;

            for (k=0;k<(6*Ns);k++)
            {
               (*psd)=sd+k*vol;
               psd+=1;
            }
         }
         else
         {
            (*brd).b2b=bsd[0][mu].b2b;
            (*brd).snd_buf[0]=bsd[0][mu].snd_buf[0];
            (*brd).snd_buf[1]=bsd[0][mu].snd_buf[1];
         }

         (*brd).b2b.n[0]=n;
         (*brd).b2b.n[1]=m;
         (*brd).b2b.ibn=bb[ifc].ibn;

         brd+=1;
      }
   }

   set_snd_req_bsd();
   init_bsd=1;
}


static void send_bufs_bsd(int ifc,int eo)
{
   int io;

   io=(ifc^nmu[ifc]);

   if (np==eo)
      MPI_Start(&snd_req_bsd[io]);
   else
      MPI_Start(&rcv_req_bsd[io^0x1]);
}


static void wait_bufs_bsd(int ifc,int eo)
{
   int io;
   MPI_Status stat_snd,stat_rcv;

   io=(ifc^nmu[ifc]);

   if (np==eo)
      MPI_Wait(&snd_req_bsd[io],&stat_snd);
   else
      MPI_Wait(&rcv_req_bsd[io^0x1],&stat_rcv);
}


b2b_flds_t *b2b_flds(int n,int mu)
{
   int ifc,vol,ibn;
   int k,*imb;
   su3_dble *ud;
   spinor_dble **sd,**rd;
   bsd_t *brd;
   b2b_flds_t *b2b;

   if (init_bsd==0)
      alloc_bsd();

   ifc=2*mu+1;
   brd=bsd[n]+mu;
   b2b=&(*brd).b2b;
   vol=(*b2b).vol;
   ibn=(*b2b).ibn;

   imb=(*brd).ise[0];
   sd=(*brd).sd[0];

   if (ibn)
      rd=(*brd).snd_buf[0];
   else
      rd=(*b2b).sde[0];

   for (k=0;k<Ns;k++)
      gather_sd(vol,imb,sd[k],rd[k]);

   imb=(*brd).ise[1];
   sd=(*brd).sd[1];

   if (ibn)
      rd=(*brd).snd_buf[1];
   else
      rd=(*b2b).sde[1];

   for (k=0;k<Ns;k++)
      gather_sd(vol,imb,sd[k],rd[k]);

   if (ibn)
      send_bufs_bsd(ifc,0);

   gather_ud(vol,(*brd).iud[0],udfld(),(*brd).ud[0]);

   if (ibn)
   {
      wait_bufs_bsd(ifc,0);
      send_bufs_bsd(ifc,1);
   }

   imb=(*brd).iso[0];
   ud=(*brd).ud[0];
   sd=(*brd).sd[0];
   rd=(*b2b).sdo[0];

   for (k=0;k<Ns;k++)
      apply_udag2sd(vol,imb,ud,sd[k],rd[k]);

   if (ibn)
   {
      wait_bufs_bsd(ifc,1);
      send_bufs_bsd(ifc^0x1,0);
   }

   gather_ud(vol,(*brd).iud[1],udfld(),(*brd).ud[1]);

   if (ibn)
   {
      wait_bufs_bsd(ifc^0x1,0);
      send_bufs_bsd(ifc^0x1,1);
   }

   imb=(*brd).iso[1];
   ud=(*brd).ud[1];
   sd=(*brd).sd[1];
   rd=(*b2b).sdo[1];

   for (k=0;k<Ns;k++)
      apply_u2sd(vol,imb,ud,sd[k],rd[k]);

   if (ibn)
      wait_bufs_bsd(ifc^0x1,1);

   return b2b;
}


static void set_snd_req_Aoe(void)
{
   int ifc,nbf;
   int tag,saddr,raddr;
   Aw_dble_t Aw;

   Aw=Awop_dble();

   for (ifc=0;ifc<8;ifc++)
   {
      nbf=2*Ns*Ns*nbbo[ifc];
      saddr=npr[ifc];
      raddr=npr[ifc^0x1];
      tag=mpi_permanent_tag();

      MPI_Send_init(Aw.Aoe[8*nbh+obbo[ifc]-nbbh],nbf,
                    MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD,&snd_req_Aoe[ifc]);
      MPI_Recv_init(rcv_buf_Aoe[ifc][0],nbf,
                    MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&rcv_req_Aoe[ifc]);

      tag=mpi_permanent_tag();

      MPI_Send_init(Aw.Aeo[8*nbh+obbo[ifc]-nbbh],nbf,
                    MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD,&snd_req_Aeo[ifc]);
      MPI_Recv_init(rcv_buf_Aeo[ifc][0],nbf,
                    MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&rcv_req_Aeo[ifc]);
   }
}


static void alloc_Aoe(void)
{
   int n,k,ifc,nmat;
   complex_dble **pzd,*zd;

   if (Ns==0)
      set_constants();

   n=nbb;
   nmat=Ns*Ns;
   pzd=malloc(n*sizeof(*pzd));
   zd=amalloc(n*nmat*sizeof(*zd),ALIGN);

   error((pzd==NULL)||(zd==NULL),1,
         "alloc_Aoe [Aw_com.c]","Unable to allocate buffers");

   set_vd2zero(n*nmat,zd);

   for (k=0;k<n;k++)
      pzd[k]=zd+k*nmat;

   for (ifc=0;ifc<8;ifc++)
   {
      rcv_buf_Aoe[ifc]=pzd;
      pzd+=nbbe[ifc^0x1];
      rcv_buf_Aeo[ifc]=pzd;
      pzd+=nbbe[ifc^0x1];
   }

   set_snd_req_Aoe();
   init_Aoe=1;
}


static void add_mat(int ifc,int vol,int *imb,complex_dble **v,complex_dble **w)
{
   int nmat;
   complex_dble **vm,*rv,*wi,*wm;

   nmat=Ns*Ns;
   vm=v+vol;

   for (;v<vm;v++)
   {
      wi=w[8*((*imb)-nbh)+ifc];
      wm=wi+nmat;
      rv=(*v);

      for (;wi<wm;wi+=4)
      {
         wi[0].re+=rv[0].re;
         wi[0].im+=rv[0].im;
         wi[1].re+=rv[1].re;
         wi[1].im+=rv[1].im;
         wi[2].re+=rv[2].re;
         wi[2].im+=rv[2].im;
         wi[3].re+=rv[3].re;
         wi[3].im+=rv[3].im;

         rv+=4;
      }

      imb+=1;
   }
}


static void send_bufs_Aoe(int ifc,int eo)
{
   int io;

   io=(ifc^nmu[ifc]);

   if (np==eo)
   {
      if (nbbo[io])
         MPI_Start(&snd_req_Aoe[io]);
   }
   else
   {
      if (nbbe[io])
         MPI_Start(&rcv_req_Aoe[io^0x1]);
   }
}


static void wait_bufs_Aoe(int ifc,int eo)
{
   int io;
   MPI_Status stat_snd,stat_rcv;

   io=(ifc^nmu[ifc]);

   if (np==eo)
   {
      if (nbbo[io])
         MPI_Wait(&snd_req_Aoe[io],&stat_snd);
   }
   else
   {
      if (nbbe[io])
         MPI_Wait(&rcv_req_Aoe[io^0x1],&stat_rcv);
   }
}


static void send_bufs_Aeo(int ifc,int eo)
{
   int io;

   io=(ifc^nmu[ifc]);

   if (np==eo)
   {
      if (nbbo[io])
         MPI_Start(&snd_req_Aeo[io]);
   }
   else
   {
      if (nbbe[io])
         MPI_Start(&rcv_req_Aeo[io^0x1]);
   }
}


static void wait_bufs_Aeo(int ifc,int eo)
{
   int io;
   MPI_Status stat_snd,stat_rcv;

   io=(ifc^nmu[ifc]);

   if (np==eo)
   {
      if (nbbo[io])
         MPI_Wait(&snd_req_Aeo[io],&stat_snd);
   }
   else
   {
      if (nbbe[io])
         MPI_Wait(&rcv_req_Aeo[io^0x1],&stat_rcv);
   }
}


void cpAoe_ext_bnd(void)
{
   int ifc,io;
   int n,m,eo;
   Aw_dble_t Aw;

   if (NPROC==1)
      return;

   if (init_Aoe==0)
      alloc_Aoe();

   Aw=Awop_dble();
   m=0;
   eo=0;

   for (n=0;n<(nsnd+2);n++)
   {
      send_bufs_Aoe(sfc[m],eo);
      wait_bufs_Aoe(sfc[m],eo);
      send_bufs_Aeo(sfc[m],eo);
      wait_bufs_Aeo(sfc[m],eo);
      m+=eo;
      eo^=0x1;
   }

   for (n=0;n<nsnd;n++)
   {
      if (m<nsnd)
         send_bufs_Aoe(sfc[m],eo);

      ifc=sfc[n];
      io=(ifc^nmu[ifc])^0x1;

      add_mat(io^0x1,nbbe[io^0x1],ipp+obbe[io^0x1],rcv_buf_Aoe[io],Aw.Aoe);

      if (m<nsnd)
      {
         wait_bufs_Aoe(sfc[m],eo);
         send_bufs_Aeo(sfc[m],eo);
      }

      add_mat(io^0x1,nbbe[io^0x1],ipp+obbe[io^0x1],rcv_buf_Aeo[io],Aw.Aeo);

      if (m<nsnd)
      {
         wait_bufs_Aeo(sfc[m],eo);
         m+=eo;
         eo^=0x1;
      }
   }
}


static void set_snd_req_Aee(void)
{
   int ifc,nbf;
   int tag,saddr,raddr;
   Aw_dble_t Aw;

   Aw=Awophat_dble();

   for (ifc=0;ifc<8;ifc++)
   {
      nbf=2*Ns*Ns*nbbo[ifc];
      saddr=npr[ifc];
      raddr=npr[ifc^0x1];
      tag=mpi_permanent_tag();

      MPI_Send_init(snd_buf_Aee[ifc][0],nbf,
                    MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD,&snd_req_Aee[ifc]);
      MPI_Recv_init(Aw.Aee[nbh+obbe[ifc^0x1]],nbf,
                    MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&rcv_req_Aee[ifc]);
   }
}


static void alloc_Aee(void)
{
   int n,k,ifc,nmat;
   complex_dble **pzd,*zd;

   if (Ns==0)
      set_constants();

   n=nbbh;
   nmat=Ns*Ns;
   pzd=malloc(n*sizeof(*pzd));
   zd=amalloc(n*nmat*sizeof(*zd),ALIGN);

   error((pzd==NULL)||(zd==NULL),1,
         "alloc_Aee [Aw_com.c]","Unable to allocate buffers");

   set_vd2zero(n*nmat,zd);

   for (k=0;k<n;k++)
      pzd[k]=zd+k*nmat;

   for (ifc=0;ifc<8;ifc++)
   {
      snd_buf_Aee[ifc]=pzd;
      pzd+=nbbo[ifc];
   }

   set_snd_req_Aee();
   init_Aee=1;
}


static void get_mat(int vol,int *imb,complex_dble **v,complex_dble **w)
{
   int nmat;
   complex_dble *vi,*vm,**wm,*rw;

   nmat=Ns*Ns;
   wm=w+vol;

   for (;w<wm;w++)
   {
      vi=v[*imb];
      vm=vi+nmat;
      rw=(*w);

      for (;vi<vm;vi+=4)
      {
         rw[0].re=vi[0].re;
         rw[0].im=vi[0].im;
         rw[1].re=vi[1].re;
         rw[1].im=vi[1].im;
         rw[2].re=vi[2].re;
         rw[2].im=vi[2].im;
         rw[3].re=vi[3].re;
         rw[3].im=vi[3].im;

         rw+=4;
      }

      imb+=1;
   }
}


static void send_bufs_Aee(int ifc,int eo)
{
   int io;

   io=(ifc^nmu[ifc]);

   if (np==eo)
   {
      if (nbbo[io])
         MPI_Start(&snd_req_Aee[io]);
   }
   else
   {
      if (nbbe[io])
         MPI_Start(&rcv_req_Aee[io^0x1]);
   }
}


static void wait_bufs_Aee(int ifc,int eo)
{
   int io;
   MPI_Status stat_snd,stat_rcv;

   io=(ifc^nmu[ifc]);

   if (np==eo)
   {
      if (nbbo[io])
         MPI_Wait(&snd_req_Aee[io],&stat_snd);
   }
   else
   {
      if (nbbe[io])
         MPI_Wait(&rcv_req_Aee[io^0x1],&stat_rcv);
   }
}


void cpAee_int_bnd(void)
{
   int ifc,io;
   int n,m,eo;
   Aw_dble_t Aw;

   if (NPROC==1)
      return;

   if (init_Aee==0)
      alloc_Aee();

   Aw=Awophat_dble();
   m=0;
   eo=0;

   for (n=0;n<nsnd;n++)
   {
      if (n>0)
         send_bufs_Aee(sfc[m],eo);

      ifc=sfc[n];
      io=ifc^nmu[ifc];

      get_mat(nbbo[io],ipp+obbo[io],Aw.Aee,snd_buf_Aee[io]);

      if (n>0)
      {
         wait_bufs_Aee(sfc[m],eo);
         m+=eo;
         eo^=0x1;
      }
   }

   while (m<nsnd)
   {
      send_bufs_Aee(sfc[m],eo);
      wait_bufs_Aee(sfc[m],eo);
      m+=eo;
      eo^=0x1;
   }
}
