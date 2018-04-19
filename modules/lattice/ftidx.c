
/*******************************************************************************
*
* File ftidx.c
*
* Copyright (C) 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Labeling of the field tensor components on the faces of the local lattice.
*
* The externally accessible functions are
*
*   ftidx_t *ftidx(void)
*     Returns an array idx[6] of ftidx_t structures containing the offsets
*     of the field tensor components on the boundaries of the local lattice
*     (see the file README.ftidx).
*
*   void plaq_ftidx(int n,int ix,int *ip)
*     Calculates the offsets ip[4] of the field tensor components at the
*     corners of the (mu,nu)-plaquette at the point in the local lattice
*     with label ix. The indices (mu,nu) are determined by the parameter
*     n=0,..,5 (see the notes).
*
* Notes:
*
* For a detailed description of the layout of the field tensor array see the
* file README.ftidx in this directory.
*
* There are six planes
*
*  (mu,nu)={(0,1),(0,2),(0,3),(2,3),(3,1),(1,2)}
*
* which may labeled by an integer n running from 0 to 5. The corners in
* the (mu,nu)-plaquette at the point x are ordered such that
*
*   ip[0] -> F_{mu nu}(x)
*   ip[1] -> F_{mu nu}(x+mu)
*   ip[2] -> F_{mu nu}(x+nu)
*   ip[3] -> F_{mu nu}(x+mu+nu)
*
* In the program plaq_ftidx() it is taken for granted that 0<=ix<VOLUME.
*
*******************************************************************************/

#define FTIDX_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "utils.h"
#include "lattice.h"
#include "global.h"

static const int plns[6][2]={{0,1},{0,2},{0,3},{2,3},{3,1},{1,2}};
static int nfc[4],ofs[4],*cn[6][2],init=0;
static ftidx_t idx[6];


static void set_nft(void)
{
   int bs[4];
   int n,mu,nu;

   bs[0]=L0;
   bs[1]=L1;
   bs[2]=L2;
   bs[3]=L3;

   nfc[0]=FACE0;
   nfc[1]=FACE1;
   nfc[2]=FACE2;
   nfc[3]=FACE3;

   ofs[0]=VOLUME;
   ofs[1]=ofs[0]+FACE0;
   ofs[2]=ofs[1]+FACE1;
   ofs[3]=ofs[2]+FACE2;

   for (n=0;n<6;n++)
   {
      mu=plns[n][0];
      nu=plns[n][1];

      idx[n].nft[0]=nfc[mu];
      idx[n].nft[1]=nfc[nu];

      if (nfc[nu]>0)
         idx[n].nft[0]+=(nfc[mu]/bs[nu]);
   }
}


static void alloc_idx(void)
{
   int n,mu,nu;
   int np,*iw;

   set_nft();
   np=0;

   for (n=0;n<6;n++)
      np+=(idx[n].nft[0]+idx[n].nft[1]);

   if (BNDRY>0)
   {
      iw=malloc((np+9*(BNDRY/2))*sizeof(*iw));
      error(iw==NULL,1,"alloc_idx [ftidx.c]",
            "Unable to allocate index arrays");
   }
   else
      iw=NULL;

   for (n=0;n<6;n++)
   {
      idx[n].ift[0]=iw;
      iw+=idx[n].nft[0];

      idx[n].ift[1]=iw;
      iw+=idx[n].nft[1];
   }

   for (n=0;n<6;n++)
   {
      mu=plns[n][0];
      nu=plns[n][1];

      cn[n][0]=iw;
      iw+=3*nfc[mu];

      cn[n][1]=iw;
      iw+=3*nfc[nu];
   }
}


static int ibnd(int mu,int iy)
{
   if (iy>(VOLUME+(BNDRY/2)))
      return iy-ofs[mu]-BNDRY/2;
   else
      return iy-ofs[mu]-nfc[mu]/2;
}


static void set_idx(void)
{
   int n,mu,nu;
   int ix,iy,iw,iz;
   int iby,ibw,ibz;
   int *ift[2],*cnn[2],nft0,nfc0,icn;

   alloc_idx();

   for (n=0;n<6;n++)
   {
      mu=plns[n][0];
      nu=plns[n][1];

      ift[0]=idx[n].ift[0];
      ift[1]=idx[n].ift[1];
      cnn[0]=cn[n][0];
      cnn[1]=cn[n][1];

      nft0=idx[n].nft[0];
      nfc0=nfc[mu];
      icn=0;

      for (ix=0;ix<VOLUME;ix++)
      {
         iy=iup[ix][mu];
         iw=iup[ix][nu];

         if (iy>=VOLUME)
         {
            iby=ibnd(mu,iy);
            ift[0][iby]=map[iy-VOLUME];

            if (iw>=VOLUME)
            {
               ibw=ibnd(nu,iw);
               ift[1][ibw]=map[iw-VOLUME];

               iz=map[iy-VOLUME];
               iz=iup[iz][nu];
               ibz=ibnd(nu,iz);
               ift[0][nfc0+icn]=VOLUME+nft0+ibz;

               cnn[0][3*iby  ]=VOLUME+iby;
               cnn[0][3*iby+1]=VOLUME+nft0+ibw;
               cnn[0][3*iby+2]=VOLUME+nfc0+icn;

               cnn[1][3*ibw  ]=cnn[0][3*iby  ];
               cnn[1][3*ibw+1]=cnn[0][3*iby+1];
               cnn[1][3*ibw+2]=cnn[0][3*iby+2];

               icn+=1;
            }
            else
            {
               iz=iup[iw][mu];
               ibz=ibnd(mu,iz);

               cnn[0][3*iby  ]=VOLUME+iby;
               cnn[0][3*iby+1]=iw;
               cnn[0][3*iby+2]=VOLUME+ibz;
            }
         }
         else if (iw>=VOLUME)
         {
            ibw=ibnd(nu,iw);
            ift[1][ibw]=map[iw-VOLUME];

            iz=iup[iy][nu];
            ibz=ibnd(nu,iz);

            cnn[1][3*ibw  ]=iy;
            cnn[1][3*ibw+1]=VOLUME+nft0+ibw;
            cnn[1][3*ibw+2]=VOLUME+nft0+ibz;
         }
      }
   }

   init=1;
}


ftidx_t *ftidx(void)
{
   if (init==0)
      set_idx();

   return idx;
}


void plaq_ftidx(int n,int ix,int *ip)
{
   int mu,nu;
   int iy,iw,k;

   if (init==0)
      set_idx();

   mu=plns[n][0];
   nu=plns[n][1];

   iy=iup[ix][mu];
   iw=iup[ix][nu];
   ip[0]=ix;

   if (iy>=VOLUME)
   {
      k=3*ibnd(mu,iy);
      ip[1]=cn[n][0][k];
      ip[2]=cn[n][0][k+1];
      ip[3]=cn[n][0][k+2];
   }
   else if (iw>=VOLUME)
   {
      k=3*ibnd(nu,iw);
      ip[1]=cn[n][1][k];
      ip[2]=cn[n][1][k+1];
      ip[3]=cn[n][1][k+2];
   }
   else
   {
      ip[1]=iy;
      ip[2]=iw;
      ip[3]=iup[iy][nu];
   }
}
