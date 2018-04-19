
/*******************************************************************************
*
* File blk_grid.c
*
* Copyright (C) 2005, 2007, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Block grid allocation.
*
* The externally accessible functions are
*
*   void alloc_bgr(blk_grid_t grid)
*     Allocates the specified block grid. The block array and the block
*     fields are put in the static memory of this module and are properly
*     initialized.
*
*   block_t *blk_list(blk_grid_t grid,int *nb,int *isw)
*     Returns the pointer to the block array of the specified grid. The
*     number of blocks on the local lattice is assigned to nb and isw is
*     set to 0 or 1 depending on whether the first block is black or white
*     (by definition it is black on the first process). If the block grid
*     is not allocated, the program returns NULL and sets nb and isw to 0.
*
* Notes:
*
* The block sizes bs[4] and other parameters of the specified block grid
* are obtained from the parameter data base. These and the lattice sizes
* must be such that the lattice can be covered by non-overlapping blocks.
* Moreover, the number of blocks in each direction must be even and the
* local lattices must contain an even number of blocks. This ensures that
* the block grid can be chessboard-coloured and that the number of blocks
* in the local lattice is the same for both colours.
*
* On all processes, the blocks at a given position in the array of blocks
* returned by blk_list() have the same position in the local lattice. The
* blocks are ordered such that the first half of them have the same colour.
* For a given colour, the blocks are ordered according to their index
*
*   n[3]+nbl[3]*n[2]+nbl[2]*nbl[3]*n[1]+nbl[1]*nbl[2]*nbl[3]*n[0],
*
* where n[mu]=bo[mu]/bs[mu] are the Cartesian coordinates of the block in
* the block grid and nbl[mu] denotes the numbers of blocks in direction mu.
* All blocks have allocated boundaries and the protection flag set.
*
* The program alloc_bgr() involves communications and must be called on all
* processes simultaneously with the same parameters. A given block grid can
* be allocated only once.
*
*******************************************************************************/

#define BLK_GRID_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "sap.h"
#include "block.h"
#include "global.h"

typedef struct
{
   int nb,isw;
   block_t *b;
} bgrid_t;

static bgrid_t bgr[(int)(BLK_GRIDS)+1]={{0,0,NULL},{0,0,NULL},{0,0,NULL}};


static block_t *blks(int *bs,int iu,int iud,int ns,int nsd,
                     int iub,int iudb,int nw,int nwd,
                     int shf,int *nb,int *isw)
{
   int bo[4];
   int n0,n1,n2,n3,m0,m1,m2,m3;
   block_t *b,*rbe,*rbo;

   n0=L0/bs[0];
   n1=L1/bs[1];
   n2=L2/bs[2];
   n3=L3/bs[3];

   (*nb)=n0*n1*n2*n3;
   (*isw)=(cpr[0]*n0+cpr[1]*n1+cpr[2]*n2+cpr[3]*n3)&0x1;

   b=malloc((*nb)*sizeof(*b));
   error(b==NULL,1,"blks [blk_grid.c]","Unable to allocate block grid");

   rbe=b;
   rbo=b+(*nb)/2;

   for (m0=0;m0<n0;m0++)
   {
      for (m1=0;m1<n1;m1++)
      {
         for (m2=0;m2<n2;m2++)
         {
            for (m3=0;m3<n3;m3++)
            {
               bo[0]=m0*bs[0];
               bo[1]=m1*bs[1];
               bo[2]=m2*bs[2];
               bo[3]=m3*bs[3];

               if (((m0+m1+m2+m3)&0x1)==0)
               {
                  if (rbe==b)
                  {
                     alloc_blk(b,bo,bs,iu,iud,ns,nsd);
                     alloc_bnd(b,iub,iudb,nw,nwd);
                     (*b).shf=shf;
                  }
                  else
                     clone_blk(b,shf,bo,rbe);

                  rbe+=1;
               }
               else
               {
                  clone_blk(b,shf,bo,rbo);
                  rbo+=1;
               }
            }
         }
      }
   }

   return b;
}


void alloc_bgr(blk_grid_t grid)
{
   int iprms[1],igr,*bs;
   int iu,iud,ns,nsd,iub,iudb,nw,nwd,shf;
   sap_parms_t sap;
   dfl_parms_t dfl;

   igr=(int)(grid);

   if (NPROC>1)
   {
      iprms[0]=igr;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);

      error(iprms[0]!=igr,1,"alloc_bgr [blk_grid.c]",
            "Parameter is not global");
   }

   error(bgr[igr].b!=NULL,1,"alloc_bgr [blk_grid.c]",
            "Block grid is already allocated");

   bs=NULL;
   iu=0;
   iud=0;
   ns=0;
   nsd=0;
   iub=0;
   iudb=0;
   nw=0;
   nwd=0;
   shf=0x0;

   if (grid==SAP_BLOCKS)
   {
      sap=sap_parms();
      error_root(sap.ncy==0,1,"alloc_bgr [blk_grid.c]",
                 "SAP parameters are not set");

      bs=sap.bs;
      iu=1;
      ns=3;
      iub=1;
      shf=0x13;
   }
   else if (grid==DFL_BLOCKS)
   {
      dfl=dfl_parms();
      error_root(dfl.Ns==0,1,"alloc_bgr [blk_grid.c]",
                 "Deflation subspace parameters are not set");

      bs=dfl.bs;
      iud=1;
      ns=dfl.Ns+1;
      nsd=dfl.Ns+1;
      shf=0xb;
   }
   else
      error_root(1,1,"alloc_bgr [blk_grid.c]","Unknown block grid");

   bgr[igr].b=blks(bs,iu,iud,ns,nsd,iub,iudb,nw,nwd,shf,
                   &(bgr[igr].nb),&(bgr[grid].isw));

   if (grid==SAP_BLOCKS)
      alloc_sap_bufs();
}


block_t *blk_list(blk_grid_t grid,int *nb,int *isw)
{
   int igr;

   igr=(int)(grid);
   (*nb)=bgr[igr].nb;
   (*isw)=bgr[igr].isw;

   return bgr[igr].b;
}
