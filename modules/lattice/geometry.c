
/*******************************************************************************
*
* File geometry.c
*
* Copyright (C) 2005, 2008, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Programs related to the lattice and block geometry.
*
* The externally accessible functions are
*
*   int ipr_global(int *n)
*     This program returns the rank of the MPI process with Cartesian
*     coordinates n[0],..,n[3] in the process grid.
*
*   void ipt_global(int *x,int *ip,int *ix)
*     Given the Cartesian coordinates x[0],..,x[3] of a point on the full
*     lattice, this program finds the local lattice containing x. On exit
*     the rank of the associated MPI process is assigned to ip and the
*     local index of the point to ix.
*
*   int global_time(int ix)
*     Returns the (global) time coordinate of the lattice point with local
*     index ix.
*
*   void geometry(void)
*     Computes the global arrays cpr,npr describing the MPI process grid
*     and the index arrays ipt,iup,idn and map that characterize the lattice
*     geometry (see main/global.h).
*
*   void blk_geometry(block_t *b)
*     Computes the index arrays b.ipt,b.iup and b.idn that describe the
*     geometry of the block b.
*
*   void blk_imbed(block_t *b)
*     Computes the index arrays b.imb and b.ibp that describe the
*     embedding of the block b in the full lattice.
*
*   void bnd_geometry(block_t *b)
*     Computes the index arrays bb.ipp and bb.map that describe the
*     geometry of the exterior boundaries bb of the block b.
*
*   void bnd_imbed(block_t *b)
*     Computes the index arrays bb.imb that describe the embedding
*     of the exterior boundaries bb of the block b in the full lattice.
*
* Notes:
*
* See main/README.global for a description of the lattice geometry and
* block/README.block for explanations of the block structure.
*
* The programs geometry() and blk_geometry() may involve communications and
* must be called simultaneously on all processes. All other programs can be
* called locally.
*
*******************************************************************************/

#define GEOMETRY_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "global.h"

#define NPROC_BLK (NPROC0_BLK*NPROC1_BLK*NPROC2_BLK*NPROC3_BLK)
#define NBLK0 (NPROC0/NPROC0_BLK)
#define NBLK1 (NPROC1/NPROC1_BLK)
#define NBLK2 (NPROC2/NPROC2_BLK)
#define NBLK3 (NPROC3/NPROC3_BLK)

static int cbs[4],cbn[4],*cbix=NULL;
static int *tms=NULL;


int ipr_global(int *n)
{
   int ib,ip;
   int n0,n1,n2,n3;
   int nb0,nb1,nb2,nb3;
   int np0,np1,np2,np3;

   n0=safe_mod(n[0],NPROC0);
   n1=safe_mod(n[1],NPROC1);
   n2=safe_mod(n[2],NPROC2);
   n3=safe_mod(n[3],NPROC3);

   nb0=n0/NPROC0_BLK;
   nb1=n1/NPROC1_BLK;
   nb2=n2/NPROC2_BLK;
   nb3=n3/NPROC3_BLK;

   np0=n0%NPROC0_BLK;
   np1=n1%NPROC1_BLK;
   np2=n2%NPROC2_BLK;
   np3=n3%NPROC3_BLK;

   ib=nb0;
   ib=ib*NBLK1+nb1;
   ib=ib*NBLK2+nb2;
   ib=ib*NBLK3+nb3;

   ip=np0;
   ip=ip*NPROC1_BLK+np1;
   ip=ip*NPROC2_BLK+np2;
   ip=ip*NPROC3_BLK+np3;

   return ip+ib*NPROC_BLK;
}


void ipt_global(int *x,int *ip,int *ix)
{
   int x0,x1,x2,x3;
   int n[4];

   x0=safe_mod(x[0],NPROC0*L0);
   x1=safe_mod(x[1],NPROC1*L1);
   x2=safe_mod(x[2],NPROC2*L2);
   x3=safe_mod(x[3],NPROC3*L3);

   n[0]=x0/L0;
   n[1]=x1/L1;
   n[2]=x2/L2;
   n[3]=x3/L3;

   (*ip)=ipr_global(n);

   x0=x0%L0;
   x1=x1%L1;
   x2=x2%L2;
   x3=x3%L3;

   (*ix)=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
}


int global_time(int ix)
{
   if ((tms!=NULL)&&(ix>=0)&&(ix<VOLUME))
      return tms[ix];
   else
   {
      error_loc(1,1,"global_time [geometry.c]",
                "Time array is not set or the argument is out of range");
      return 0;
   }
}


static void set_cpr(void)
{
   int ib,ip;
   int np,nr;

   MPI_Comm_size(MPI_COMM_WORLD,&np);

   error(np!=NPROC,1,"set_cpr [geometry.c]",
         "Actual number of processes does not match NPROC");

   MPI_Comm_rank(MPI_COMM_WORLD,&nr);

   error((nr<0)||(nr>=NPROC),1,"set_cpr [geometry.c]",
         "Rank of process is out of range");

   ib=nr/NPROC_BLK;
   ip=nr%NPROC_BLK;

   cpr[3]=(ib%NBLK3)*NPROC3_BLK+(ip%NPROC3_BLK);
   ib/=NBLK3;
   ip/=NPROC3_BLK;

   cpr[2]=(ib%NBLK2)*NPROC2_BLK+(ip%NPROC2_BLK);
   ib/=NBLK2;
   ip/=NPROC2_BLK;

   cpr[1]=(ib%NBLK1)*NPROC1_BLK+(ip%NPROC1_BLK);
   ib/=NBLK1;
   ip/=NPROC1_BLK;

   cpr[0]=ib*NPROC0_BLK+ip;
}


static void set_npr(void)
{
   int mu,n[4];

   for (mu=0;mu<4;mu++)
      n[mu]=cpr[mu];

   for (mu=0;mu<4;mu++)
   {
      n[mu]-=1;
      npr[2*mu]=ipr_global(n);
      n[mu]+=2;
      npr[2*mu+1]=ipr_global(n);
      n[mu]-=1;
   }
}


static void cache_block(int *bs)
{
   int mu;

   cbs[0]=bs[0];
   cbn[0]=1;

   for (mu=1;mu<4;mu++)
   {
      if ((bs[mu]%4)==0)
         cbs[mu]=4;
      else if ((bs[mu]%3)==0)
         cbs[mu]=3;
      else if ((bs[mu]%2)==0)
         cbs[mu]=2;
      else
         cbs[mu]=1;

      cbn[mu]=bs[mu]/cbs[mu];
   }

   if (cbix!=NULL)
      free(cbix);

   cbix=malloc(cbs[0]*cbs[1]*cbs[2]*cbs[3]*sizeof(*cbix));
   error(cbix==NULL,1,"cache_block [geometry.c]",
         "Unable to allocate auxiliary array");
}


static void set_cbix(void)
{
   int x0,x1,x2,x3;
   int ig,iu,ib,is;

   ig=0;
   iu=0;

   for (x0=0;x0<cbs[0];x0++)
   {
      for (x1=0;x1<cbs[1];x1++)
      {
         for (x2=0;x2<cbs[2];x2++)
         {
            for (x3=0;x3<cbs[3];x3++)
            {
               ib=x3+cbs[3]*x2+cbs[2]*cbs[3]*x1+cbs[1]*cbs[2]*cbs[3]*x0;
               is=x0+x1+x2+x3;

               if ((is%2)==0)
               {
                  cbix[ib]=ig;
                  ig+=1;
               }
               else
               {
                  cbix[ib]=iu;
                  iu+=1;
               }
            }
         }
      }
   }
}


static int index(int *bo,int *bs,int x0,int x1,int x2,int x3)
{
   int y0,y1,y2,y3;
   int xb1,xb2,xb3;
   int xn1,xn2,xn3;
   int ib,in,is;

   y0=safe_mod(x0,bs[0]);
   y1=safe_mod(x1,bs[1]);
   y2=safe_mod(x2,bs[2]);
   y3=safe_mod(x3,bs[3]);

   xb1=y1%cbs[1];
   xb2=y2%cbs[2];
   xb3=y3%cbs[3];

   xn1=y1/cbs[1];
   xn2=y2/cbs[2];
   xn3=y3/cbs[3];

   ib=cbix[xb3+cbs[3]*xb2+cbs[2]*cbs[3]*xb1+cbs[1]*cbs[2]*cbs[3]*y0];
   in=xn3+cbn[3]*xn2+cbn[3]*cbn[2]*xn1;
   is=y0+y1+y2+y3;
   is+=(bo[0]+bo[1]+bo[2]+bo[3]);

   if ((is%2)!=0)
      ib+=((bs[0]*bs[1]*bs[2]*bs[3])/2);

   return ib+(cbs[0]*cbs[1]*cbs[2]*cbs[3]*in)/2;
}


static void set_tms(void)
{
   int ix,iy,x0;

   if (tms!=NULL)
      free(tms);

   tms=malloc(VOLUME*sizeof(*tms));
   error(tms==NULL,1,"set_tms [geometry.c]",
         "Unable to allocate time array");

   for (iy=0;iy<VOLUME;iy++)
   {
      x0=iy/(L1*L2*L3);
      ix=ipt[iy];

      tms[ix]=x0+cpr[0]*L0;
   }
}


void geometry(void)
{
   int x0,x1,x2,x3;
   int k,mu,ix,iy,iz,iw;
   int bo[4],bs[4],ifc[8];

   set_cpr();
   set_npr();

   bo[0]=0;
   bo[1]=0;
   bo[2]=0;
   bo[3]=0;

   bs[0]=L0;
   bs[1]=L1;
   bs[2]=L2;
   bs[3]=L3;

   cache_block(bs);
   set_cbix();

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               ix=index(bo,bs,x0,x1,x2,x3);
               iy=x3+L3*x2+L2*L3*x1+L1*L2*L3*x0;
               ipt[iy]=ix;

               iup[ix][0]=index(bo,bs,x0+1,x1,x2,x3);
               idn[ix][0]=index(bo,bs,x0-1,x1,x2,x3);

               iup[ix][1]=index(bo,bs,x0,x1+1,x2,x3);
               idn[ix][1]=index(bo,bs,x0,x1-1,x2,x3);

               iup[ix][2]=index(bo,bs,x0,x1,x2+1,x3);
               idn[ix][2]=index(bo,bs,x0,x1,x2-1,x3);

               iup[ix][3]=index(bo,bs,x0,x1,x2,x3+1);
               idn[ix][3]=index(bo,bs,x0,x1,x2,x3-1);

               if ((x0==(L0-1))&&(NPROC0>1))
                  iup[ix][0]=VOLUME;
               if ((x0==0)&&(NPROC0>1))
                  idn[ix][0]=VOLUME;

               if ((x1==(L1-1))&&(NPROC1>1))
                  iup[ix][1]=VOLUME;
               if ((x1==0)&&(NPROC1>1))
                  idn[ix][1]=VOLUME;

               if ((x2==(L2-1))&&(NPROC2>1))
                  iup[ix][2]=VOLUME;
               if ((x2==0)&&(NPROC2>1))
                  idn[ix][2]=VOLUME;

               if ((x3==(L3-1))&&(NPROC3>1))
                  iup[ix][3]=VOLUME;
               if ((x3==0)&&(NPROC3>1))
                  idn[ix][3]=VOLUME;
            }
         }
      }
   }

   ifc[0]=0;
   ifc[1]=ifc[0]+(FACE0/2);
   ifc[2]=ifc[1]+(FACE0/2);
   ifc[3]=ifc[2]+(FACE1/2);
   ifc[4]=ifc[3]+(FACE1/2);
   ifc[5]=ifc[4]+(FACE2/2);
   ifc[6]=ifc[5]+(FACE2/2);
   ifc[7]=ifc[6]+(FACE3/2);

   for (ix=0;ix<VOLUME;ix++)
   {
      if (ix==(VOLUME/2))
      {
         ifc[0]=(BNDRY/2);
         ifc[1]=ifc[0]+(FACE0/2);
         ifc[2]=ifc[1]+(FACE0/2);
         ifc[3]=ifc[2]+(FACE1/2);
         ifc[4]=ifc[3]+(FACE1/2);
         ifc[5]=ifc[4]+(FACE2/2);
         ifc[6]=ifc[5]+(FACE2/2);
         ifc[7]=ifc[6]+(FACE3/2);
      }

      iy=(ix+(VOLUME/2))%VOLUME;

      for (mu=0;mu<4;mu++)
      {
         if (idn[iy][mu]==VOLUME)
         {
            iz=ifc[2*mu];
            ifc[2*mu]+=1;

            idn[iy][mu]=VOLUME+iz;
            iw=iy;

            for (k=1;k<bs[mu];k++)
               iw=iup[iw][mu];

            map[iz]=iw;
         }

         if (iup[iy][mu]==VOLUME)
         {
            iz=ifc[2*mu+1];
            ifc[2*mu+1]+=1;

            iup[iy][mu]=VOLUME+iz;
            iw=iy;

            for (k=1;k<bs[mu];k++)
               iw=idn[iw][mu];

            map[iz]=iw;
         }
      }
   }

   set_tms();
   free(cbix);
   cbix=NULL;
}


void blk_geometry(block_t *b)
{
   int *bo,*bs;
   int x0,x1,x2,x3;
   int ix,iy;

   error(iup[0][0]==0,1,"blk_geometry [geometry.c]",
         "The global geometry arrays are not initialized");

   bo=(*b).bo;
   bs=(*b).bs;

   cache_block(bs);
   set_cbix();

   for (x0=0;x0<bs[0];x0++)
   {
      for (x1=0;x1<bs[1];x1++)
      {
         for (x2=0;x2<bs[2];x2++)
         {
            for (x3=0;x3<bs[3];x3++)
            {
               ix=index(bo,bs,x0,x1,x2,x3);
               iy=x3+bs[3]*x2+bs[2]*bs[3]*x1+bs[1]*bs[2]*bs[3]*x0;
               (*b).ipt[iy]=ix;

               if ((x0+1)<bs[0])
                  (*b).iup[ix][0]=index(bo,bs,x0+1,x1,x2,x3);
               else
                  (*b).iup[ix][0]=(*b).vol;

               if (x0>0)
                  (*b).idn[ix][0]=index(bo,bs,x0-1,x1,x2,x3);
               else
                  (*b).idn[ix][0]=(*b).vol;

               if ((x1+1)<bs[1])
                  (*b).iup[ix][1]=index(bo,bs,x0,x1+1,x2,x3);
               else
                  (*b).iup[ix][1]=(*b).vol;

               if (x1>0)
                  (*b).idn[ix][1]=index(bo,bs,x0,x1-1,x2,x3);
               else
                  (*b).idn[ix][1]=(*b).vol;

               if ((x2+1)<bs[2])
                  (*b).iup[ix][2]=index(bo,bs,x0,x1,x2+1,x3);
               else
                  (*b).iup[ix][2]=(*b).vol;

               if (x2>0)
                  (*b).idn[ix][2]=index(bo,bs,x0,x1,x2-1,x3);
               else
                  (*b).idn[ix][2]=(*b).vol;

               if ((x3+1)<bs[3])
                  (*b).iup[ix][3]=index(bo,bs,x0,x1,x2,x3+1);
               else
                  (*b).iup[ix][3]=(*b).vol;

               if (x3>0)
                  (*b).idn[ix][3]=index(bo,bs,x0,x1,x2,x3-1);
               else
                  (*b).idn[ix][3]=(*b).vol;
            }
         }
      }
   }

   (*b).ipt[(*b).vol]=(*b).ipt[0];

   free(cbix);
   cbix=NULL;
}


void blk_imbed(block_t *b)
{
   int *bo,*bs;
   int x0,x1,x2,x3;
   int ix,iy,ibd,ibu,*ibp;

   bo=(*b).bo;
   bs=(*b).bs;

   for (x0=0;x0<bs[0];x0++)
   {
      for (x1=0;x1<bs[1];x1++)
      {
         for (x2=0;x2<bs[2];x2++)
         {
            for (x3=0;x3<bs[3];x3++)
            {
               iy=x3+bs[3]*x2+bs[2]*bs[3]*x1+bs[1]*bs[2]*bs[3]*x0;
               ix=(*b).ipt[iy];

               iy=(bo[3]+x3)+L3*(bo[2]+x2)+L2*L3*(bo[1]+x1)+L1*L2*L3*(bo[0]+x0);
               (*b).imb[ix]=ipt[iy];
            }
         }
      }
   }

   (*b).imb[(*b).vol]=(*b).imb[0];

   ibd=((cpr[0]==0)&&((*b).bo[0]==0)&&(bc_type()!=3));
   ibu=((cpr[0]==(NPROC0-1))&&(((*b).bo[0]+(*b).bs[0])==L0)&&(bc_type()==0));
   ibp=(*b).ibp;

   for (ix=0;ix<(*b).vol;ix++)
   {
      if (((ibd)&&((*b).idn[ix][0]==(*b).vol))||
          ((ibu)&&((*b).iup[ix][0]==(*b).vol)))
      {
         (*ibp)=ix;
         ibp+=1;
      }
   }
}


void bnd_geometry(block_t *b)
{
   int ifc,mu,ix,iy,iw,iz;
   int vol,volh,*ipp[8],*map[8];
   bndry_t *bb;

   vol=(*b).vol;
   volh=vol/2;
   bb=(*b).bb;

   for (ifc=0;ifc<8;ifc++)
   {
      ipp[ifc]=bb[ifc].ipp;
      map[ifc]=bb[ifc].map;
   }

   for (ix=0;ix<vol;ix++)
   {
      if (ix<volh)
         iy=ix+volh;
      else
         iy=ix-volh;

      for (mu=0;mu<4;mu++)
      {
         if ((*b).iup[iy][mu]==vol)
         {
            ifc=2*mu+1;
            ipp[ifc][0]=iy;
            ipp[ifc]+=1;

            iw=iy;
            iz=iy;

            while (iw<vol)
            {
               iz=iw;
               iw=(*b).idn[iw][mu];
            }

            map[ifc][0]=iz;
            map[ifc]+=1;
         }

         if ((*b).idn[iy][mu]==vol)
         {
            ifc=2*mu;
            ipp[ifc][0]=iy;
            ipp[ifc]+=1;

            iw=iy;
            iz=iy;

            while (iw<vol)
            {
               iz=iw;
               iw=(*b).iup[iw][mu];
            }

            map[ifc][0]=iz;
            map[ifc]+=1;
         }
      }
   }

   for (ifc=0;ifc<8;ifc++)
   {
      vol=(*bb).vol;
      (*bb).ipp[vol]=(*bb).ipp[0];
      (*bb).map[vol]=(*bb).map[0];
      bb+=1;
   }
}


void bnd_imbed(block_t *b)
{
   int ifc,ix,iy;
   int vol,*ipp,*imb;
   bndry_t *bb;

   bb=(*b).bb;
   imb=(*b).imb;

   for (ifc=0;ifc<8;ifc++)
   {
      vol=(*bb).vol;
      ipp=(*bb).ipp;

      for (ix=0;ix<vol;ix++)
      {
         iy=imb[ipp[ix]];

         if (ifc&0x1)
            (*bb).imb[ix]=iup[iy][ifc/2];
         else
            (*bb).imb[ix]=idn[iy][ifc/2];
      }

      (*bb).imb[vol]=(*bb).imb[0];

      if ((*bb).imb[0]>=VOLUME)
         (*bb).ibn=1;
      else
         (*bb).ibn=0;

      bb+=1;
   }
}
