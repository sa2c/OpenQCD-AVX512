
/*******************************************************************************
*
* File map_s2blk.c
*
* Copyright (C) 2005, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Copying of the spinor fields to and from the blocks in a block grid
*
* The externally accessible functions are
*
*   void assign_s2sblk(blk_grid_t grid,int n,ptset_t set,spinor *s,int k)
*     Assigns the relevant part of the global single-precision spinor field s
*     to the single-precision field b.s[k] on the n'th block of the specified
*     block grid. Depending on the specified point set, the field on the even, 
*     odd or all points is copied.
*
*   void assign_sblk2s(blk_grid_t grid,int n,ptset_t set,int k,spinor *s)
*     Assigns the single-precision spinor field b.s[k] on the n'th block of
*     the specified block grid to the relevant part of the global single-
*     precision field s. Depending on the specified point set, the field on 
*     the even, odd or all points is copied.
*
*   void assign_s2sdblk(blk_grid_t grid,int n,ptset_t set,spinor *s,int k)
*     Assigns the relevant part of the global single-precision spinor field s
*     to the double-precision field b.sd[k] on the n'th block of the specified
*     block grid. Depending on the specified point set, the field on the even, 
*     odd or all points is copied.
*
*   void assign_sd2sdblk(blk_grid_t grid,int n,ptset_t set,
*                        spinor_dble *sd,int k)
*     Assigns the relevant part of the global double-precision spinor field sd
*     to the double-precision field b.sd[k] on the n'th block of the specified
*     block grid. Depending on the specified point set, the field on the even, 
*     odd or all points is copied.
*
*   void assign_sdblk2sd(blk_grid_t grid,int n,ptset_t set,
*                        int k,spinor_dble *sd)
*     Assigns the single-precision spinor field b.sd[k] on the n'th block of
*     the specified block grid to the relevant part of the global single-
*     precision field sd. Depending on the specified point set, the field on
*     the even, odd or all points is copied.
*
* Notes:
*
* Only the spinors residing on the blocks (but not those on the boundaries
* of the blocks) are copied. All these programs can be called locally.
*
*******************************************************************************/

#define MAP_S2BLK_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "block.h"
#include "global.h"

#if (defined x64)
#include "sse2.h"

void assign_s2sblk(blk_grid_t grid,int n,ptset_t set,spinor *s,int k)
{
   int nb,isw,vol,*imb;
   spinor *sb,*sm,*rs1,*rs2;
   block_t *b;

   b=blk_list(grid,&nb,&isw)+n;

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"assign_s2sblk [map_s2blk.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(k>=(*b).ns))
   {
      error_loc(1,1,"assign_s2sblk [map_s2blk.c]",
                "Block field number is out of range");
      return;
   }

   vol=(*b).vol;
   imb=(*b).imb;
   sb=(*b).s[k];
   sm=sb;

   if (set==ALL_PTS)
      sm+=vol;
   else if (set==EVEN_PTS)
      sm+=vol/2;
   else if (set==ODD_PTS)
   {
      imb+=vol/2;
      sb+=vol/2;
      sm+=vol;
   }

   rs2=s+(*imb);
   
   for (;sb<sm;sb++)
   {
      rs1=rs2;
      imb+=1;
      rs2=s+(*imb);
      _prefetch_spinor(rs2);      

      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %2, %%xmm1 \n\t"
                            "movaps %4, %%xmm2"
                            :
                            :
                            "m" ((*rs1).c1.c1),
                            "m" ((*rs1).c1.c2),
                            "m" ((*rs1).c1.c3),
                            "m" ((*rs1).c2.c1),
                            "m" ((*rs1).c2.c2),
                            "m" ((*rs1).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2");

      __asm__ __volatile__ ("movaps %0, %%xmm3 \n\t"
                            "movaps %2, %%xmm4 \n\t"
                            "movaps %4, %%xmm5"
                            :
                            :
                            "m" ((*rs1).c3.c1),
                            "m" ((*rs1).c3.c2),
                            "m" ((*rs1).c3.c3),
                            "m" ((*rs1).c4.c1),
                            "m" ((*rs1).c4.c2),
                            "m" ((*rs1).c4.c3)
                            :
                            "xmm3", "xmm4", "xmm5");      

      __asm__ __volatile__ ("movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm1, %2 \n\t"
                            "movaps %%xmm2, %4"
                            :
                            "=m" ((*sb).c1.c1),
                            "=m" ((*sb).c1.c2),
                            "=m" ((*sb).c1.c3),
                            "=m" ((*sb).c2.c1),
                            "=m" ((*sb).c2.c2),
                            "=m" ((*sb).c2.c3));

      __asm__ __volatile__ ("movaps %%xmm3, %0 \n\t"
                            "movaps %%xmm4, %2 \n\t"
                            "movaps %%xmm5, %4"
                            :
                            "=m" ((*sb).c3.c1),
                            "=m" ((*sb).c3.c2),
                            "=m" ((*sb).c3.c3),
                            "=m" ((*sb).c4.c1),
                            "=m" ((*sb).c4.c2),
                            "=m" ((*sb).c4.c3));      
   }
}


void assign_sblk2s(blk_grid_t grid,int n,ptset_t set,int k,spinor *s)
{
   int nb,isw,vol,*imb;
   spinor *sb,*sm,*rs1;
   block_t *b;

   b=blk_list(grid,&nb,&isw)+n;

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"assign_sblk2s [map_s2blk.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(k>=(*b).ns))
   {
      error_loc(1,1,"assign_sblk2s [map_s2blk.c]",
                "Block field number is out of range");
      return;
   }

   vol=(*b).vol;
   imb=(*b).imb;
   sb=(*b).s[k];
   sm=sb;

   if (set==ALL_PTS)
      sm+=vol;
   else if (set==EVEN_PTS)
      sm+=vol/2;
   else if (set==ODD_PTS)
   {
      imb+=vol/2;
      sb+=vol/2;
      sm+=vol;
   }

   for (;sb<sm;sb++)
   {
      rs1=s+(*imb);
      imb+=1;
                            
      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %2, %%xmm1 \n\t"
                            "movaps %4, %%xmm2"
                            :
                            :
                            "m" ((*sb).c1.c1),
                            "m" ((*sb).c1.c2),
                            "m" ((*sb).c1.c3),
                            "m" ((*sb).c2.c1),
                            "m" ((*sb).c2.c2),
                            "m" ((*sb).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2");

      __asm__ __volatile__ ("movaps %0, %%xmm3 \n\t"
                            "movaps %2, %%xmm4 \n\t"
                            "movaps %4, %%xmm5"
                            :
                            :
                            "m" ((*sb).c3.c1),
                            "m" ((*sb).c3.c2),
                            "m" ((*sb).c3.c3),
                            "m" ((*sb).c4.c1),
                            "m" ((*sb).c4.c2),
                            "m" ((*sb).c4.c3)
                            :
                            "xmm3", "xmm4", "xmm5");

      __asm__ __volatile__ ("movaps %%xmm0, %0 \n\t"
                            "movaps %%xmm1, %2 \n\t"
                            "movaps %%xmm2, %4"
                            :
                            "=m" ((*rs1).c1.c1),
                            "=m" ((*rs1).c1.c2),
                            "=m" ((*rs1).c1.c3),
                            "=m" ((*rs1).c2.c1),
                            "=m" ((*rs1).c2.c2),
                            "=m" ((*rs1).c2.c3));

      __asm__ __volatile__ ("movaps %%xmm3, %0 \n\t"
                            "movaps %%xmm4, %2 \n\t"
                            "movaps %%xmm5, %4"
                            :
                            "=m" ((*rs1).c3.c1),
                            "=m" ((*rs1).c3.c2),
                            "=m" ((*rs1).c3.c3),
                            "=m" ((*rs1).c4.c1),
                            "=m" ((*rs1).c4.c2),
                            "=m" ((*rs1).c4.c3));      
   }
}


void assign_s2sdblk(blk_grid_t grid,int n,ptset_t set,spinor *s,int k)
{
   int nb,isw,vol,*imb;
   spinor *rs1,*rs2;  
   spinor_dble *sb,*sm;
   block_t *b;

   b=blk_list(grid,&nb,&isw)+n;

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"assign_s2sdblk [map_s2blk.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(k>=(*b).nsd))
   {
      error_loc(1,1,"assign_s2sdblk [map_s2blk.c]",
                "Block field number is out of range");
      return;
   }

   vol=(*b).vol;
   imb=(*b).imb;
   sb=(*b).sd[k];
   sm=sb;

   if (set==ALL_PTS)
      sm+=vol;
   else if (set==EVEN_PTS)
      sm+=vol/2;
   else if (set==ODD_PTS)
   {
      imb+=vol/2;
      sb+=vol/2;
      sm+=vol;
   }

   rs2=s+(*imb);
   
   for (;sb<sm;sb++)
   {
      rs1=rs2;
      imb+=1;
      rs2=s+(*imb);
      _prefetch_spinor(rs2);

      __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                            "movaps %2, %%xmm2 \n\t"
                            "movaps %4, %%xmm4"
                            :
                            :
                            "m" ((*rs1).c1.c1),
                            "m" ((*rs1).c1.c2),
                            "m" ((*rs1).c1.c3),
                            "m" ((*rs1).c2.c1),
                            "m" ((*rs1).c2.c2),
                            "m" ((*rs1).c2.c3)
                            :
                            "xmm0", "xmm2", "xmm4");

      __asm__ __volatile__ ("movaps %0, %%xmm6 \n\t"
                            "movaps %2, %%xmm8 \n\t"
                            "movaps %4, %%xmm10"
                            :
                            :
                            "m" ((*rs1).c3.c1),
                            "m" ((*rs1).c3.c2),
                            "m" ((*rs1).c3.c3),
                            "m" ((*rs1).c4.c1),
                            "m" ((*rs1).c4.c2),
                            "m" ((*rs1).c4.c3)
                            :
                            "xmm6", "xmm8", "xmm10");      

      __asm__ __volatile__ ("movhlps %%xmm0, %%xmm1 \n\t"
                            "movhlps %%xmm2, %%xmm3 \n\t"
                            "movhlps %%xmm4, %%xmm5 \n\t"
                            "movhlps %%xmm6, %%xmm7 \n\t"
                            "movhlps %%xmm8, %%xmm9 \n\t"
                            "movhlps %%xmm10, %%xmm11"
                            :
                            :
                            :
                            "xmm1", "xmm3", "xmm5", "xmm7",
                            "xmm9", "xmm10");

      __asm__ __volatile__ ("cvtps2pd %%xmm0, %%xmm0 \n\t"
                            "cvtps2pd %%xmm1, %%xmm1 \n\t"
                            "cvtps2pd %%xmm2, %%xmm2 \n\t"
                            "cvtps2pd %%xmm3, %%xmm3 \n\t"
                            "cvtps2pd %%xmm4, %%xmm4 \n\t"
                            "cvtps2pd %%xmm5, %%xmm5 \n\t"
                            "cvtps2pd %%xmm6, %%xmm6 \n\t"
                            "cvtps2pd %%xmm7, %%xmm7 \n\t"
                            "cvtps2pd %%xmm8, %%xmm8 \n\t"
                            "cvtps2pd %%xmm9, %%xmm9 \n\t"
                            "cvtps2pd %%xmm10, %%xmm10 \n\t"
                            "cvtps2pd %%xmm11, %%xmm11"
                            :
                            :
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5", "xmm6", "xmm7",
                            "xmm8", "xmm9", "xmm10", "xmm11");
      
      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"
                            "movapd %%xmm5, %5"                            
                            :
                            "=m" ((*sb).c1.c1),
                            "=m" ((*sb).c1.c2),
                            "=m" ((*sb).c1.c3),
                            "=m" ((*sb).c2.c1),
                            "=m" ((*sb).c2.c2),
                            "=m" ((*sb).c2.c3));

      __asm__ __volatile__ ("movapd %%xmm6, %0 \n\t"
                            "movapd %%xmm7, %1 \n\t"
                            "movapd %%xmm8, %2 \n\t"
                            "movapd %%xmm9, %3 \n\t"
                            "movapd %%xmm10, %4 \n\t"
                            "movapd %%xmm11, %5"                            
                            :
                            "=m" ((*sb).c3.c1),
                            "=m" ((*sb).c3.c2),
                            "=m" ((*sb).c3.c3),
                            "=m" ((*sb).c4.c1),
                            "=m" ((*sb).c4.c2),
                            "=m" ((*sb).c4.c3));
   }
}


void assign_sd2sdblk(blk_grid_t grid,int n,ptset_t set,
                     spinor_dble *sd,int k)
{
   int nb,isw,vol,*imb;
   spinor_dble *sb,*sm,*rs1,*rs2;
   block_t *b;

   b=blk_list(grid,&nb,&isw)+n;

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"assign_sd2sdblk [map_s2blk.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(k>=(*b).nsd))
   {
      error_loc(1,1,"assign_sd2sdblk [map_s2blk.c]",
                "Block field number is out of range");
      return;
   }

   vol=(*b).vol;
   imb=(*b).imb;
   sb=(*b).sd[k];
   sm=sb;

   if (set==ALL_PTS)
      sm+=vol;
   else if (set==EVEN_PTS)
      sm+=vol/2;
   else if (set==ODD_PTS)
   {
      imb+=vol/2;
      sb+=vol/2;
      sm+=vol;
   }

   rs2=sd+(*imb);
   
   for (;sb<sm;sb++)
   {
      rs1=rs2;
      imb+=1;
      rs2=sd+(*imb);
      _prefetch_spinor_dble(rs2);      

      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5"
                            :
                            :
                            "m" ((*rs1).c1.c1),
                            "m" ((*rs1).c1.c2),
                            "m" ((*rs1).c1.c3),
                            "m" ((*rs1).c2.c1),
                            "m" ((*rs1).c2.c2),
                            "m" ((*rs1).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                            "movapd %1, %%xmm7 \n\t"
                            "movapd %2, %%xmm8 \n\t"
                            "movapd %3, %%xmm9 \n\t"
                            "movapd %4, %%xmm10 \n\t"
                            "movapd %5, %%xmm11"                            
                            :
                            :
                            "m" ((*rs1).c3.c1),
                            "m" ((*rs1).c3.c2),
                            "m" ((*rs1).c3.c3),
                            "m" ((*rs1).c4.c1),
                            "m" ((*rs1).c4.c2),
                            "m" ((*rs1).c4.c3)
                            :
                            "xmm6", "xmm7", "xmm8", "xmm9",
                            "xmm10", "xmm11");
      
      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"
                            "movapd %%xmm5, %5"
                            :
                            "=m" ((*sb).c1.c1),
                            "=m" ((*sb).c1.c2),
                            "=m" ((*sb).c1.c3),
                            "=m" ((*sb).c2.c1),
                            "=m" ((*sb).c2.c2),
                            "=m" ((*sb).c2.c3));

      __asm__ __volatile__ ("movapd %%xmm6, %0 \n\t"
                            "movapd %%xmm7, %1 \n\t"
                            "movapd %%xmm8, %2 \n\t"
                            "movapd %%xmm9, %3 \n\t"
                            "movapd %%xmm10, %4 \n\t"
                            "movapd %%xmm11, %5"
                            :
                            "=m" ((*sb).c3.c1),
                            "=m" ((*sb).c3.c2),
                            "=m" ((*sb).c3.c3),
                            "=m" ((*sb).c4.c1),
                            "=m" ((*sb).c4.c2),
                            "=m" ((*sb).c4.c3));      
   }
}


void assign_sdblk2sd(blk_grid_t grid,int n,ptset_t set,
                     int k,spinor_dble *sd)
{
   int nb,isw,vol,*imb;
   spinor_dble *sb,*sm,*rs1;
   block_t *b;

   b=blk_list(grid,&nb,&isw)+n;

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"assign_sdblk2sd [map_s2blk.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(k>=(*b).nsd))
   {
      error_loc(1,1,"assign_sdblk2sd [map_s2blk.c]",
                "Block field number is out of range");
      return;
   }

   vol=(*b).vol;
   imb=(*b).imb;
   sb=(*b).sd[k];
   sm=sb;

   if (set==ALL_PTS)
      sm+=vol;
   else if (set==EVEN_PTS)
      sm+=vol/2;
   else if (set==ODD_PTS)
   {
      imb+=vol/2;
      sb+=vol/2;
      sm+=vol;
   }

   for (;sb<sm;sb++)
   {
      rs1=sd+(*imb);
      imb+=1;
                            
      __asm__ __volatile__ ("movapd %0, %%xmm0 \n\t"
                            "movapd %1, %%xmm1 \n\t"
                            "movapd %2, %%xmm2 \n\t"
                            "movapd %3, %%xmm3 \n\t"
                            "movapd %4, %%xmm4 \n\t"
                            "movapd %5, %%xmm5"
                            :
                            :
                            "m" ((*sb).c1.c1),
                            "m" ((*sb).c1.c2),
                            "m" ((*sb).c1.c3),
                            "m" ((*sb).c2.c1),
                            "m" ((*sb).c2.c2),
                            "m" ((*sb).c2.c3)
                            :
                            "xmm0", "xmm1", "xmm2", "xmm3",
                            "xmm4", "xmm5");

      __asm__ __volatile__ ("movapd %0, %%xmm6 \n\t"
                            "movapd %1, %%xmm7 \n\t"
                            "movapd %2, %%xmm8 \n\t"
                            "movapd %3, %%xmm9 \n\t"
                            "movapd %4, %%xmm10 \n\t"
                            "movapd %5, %%xmm11"
                            :
                            :
                            "m" ((*sb).c3.c1),
                            "m" ((*sb).c3.c2),
                            "m" ((*sb).c3.c3),
                            "m" ((*sb).c4.c1),
                            "m" ((*sb).c4.c2),
                            "m" ((*sb).c4.c3)
                            :
                            "xmm6", "xmm7", "xmm8", "xmm9",
                            "xmm10", "xmm11");
      
      __asm__ __volatile__ ("movapd %%xmm0, %0 \n\t"
                            "movapd %%xmm1, %1 \n\t"
                            "movapd %%xmm2, %2 \n\t"
                            "movapd %%xmm3, %3 \n\t"
                            "movapd %%xmm4, %4 \n\t"
                            "movapd %%xmm5, %5"
                            :
                            "=m" ((*rs1).c1.c1),
                            "=m" ((*rs1).c1.c2),
                            "=m" ((*rs1).c1.c3),
                            "=m" ((*rs1).c2.c1),
                            "=m" ((*rs1).c2.c2),
                            "=m" ((*rs1).c2.c3));

      __asm__ __volatile__ ("movapd %%xmm6, %0 \n\t"
                            "movapd %%xmm7, %1 \n\t"
                            "movapd %%xmm8, %2 \n\t"
                            "movapd %%xmm9, %3 \n\t"
                            "movapd %%xmm10, %4 \n\t"
                            "movapd %%xmm11, %5"
                            :
                            "=m" ((*rs1).c3.c1),
                            "=m" ((*rs1).c3.c2),
                            "=m" ((*rs1).c3.c3),
                            "=m" ((*rs1).c4.c1),
                            "=m" ((*rs1).c4.c2),
                            "=m" ((*rs1).c4.c3));      
   }
}

#else

void assign_s2sblk(blk_grid_t grid,int n,ptset_t set,spinor *s,int k)
{
   int nb,isw,vol,*imb;
   spinor *sb,*sm;
   block_t *b;

   b=blk_list(grid,&nb,&isw)+n;

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"assign_s2sblk [map_s2blk.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(k>=(*b).ns))
   {
      error_loc(1,1,"assign_s2sblk [map_s2blk.c]",
                "Block field number is out of range");
      return;
   }

   vol=(*b).vol;
   imb=(*b).imb;
   sb=(*b).s[k];
   sm=sb;

   if (set==ALL_PTS)
      sm+=vol;
   else if (set==EVEN_PTS)
      sm+=vol/2;
   else if (set==ODD_PTS)
   {
      imb+=vol/2;
      sb+=vol/2;
      sm+=vol;
   }

   for (;sb<sm;sb++)
   {
      (*sb)=s[*imb];
      imb+=1;
   }
}


void assign_sblk2s(blk_grid_t grid,int n,ptset_t set,int k,spinor *s)
{
   int nb,isw,vol,*imb;
   spinor *sb,*sm;
   block_t *b;

   b=blk_list(grid,&nb,&isw)+n;

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"assign_sblk2s [map_s2blk.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(k>=(*b).ns))
   {
      error_loc(1,1,"assign_sblk2s [map_s2blk.c]",
                "Block field number is out of range");
      return;
   }

   vol=(*b).vol;
   imb=(*b).imb;
   sb=(*b).s[k];
   sm=sb;

   if (set==ALL_PTS)
      sm+=vol;
   else if (set==EVEN_PTS)
      sm+=vol/2;
   else if (set==ODD_PTS)
   {
      imb+=vol/2;
      sb+=vol/2;
      sm+=vol;
   }

   for (;sb<sm;sb++)
   {
      s[*imb]=(*sb);
      imb+=1;
   }
}


void assign_s2sdblk(blk_grid_t grid,int n,ptset_t set,spinor *s,int k)
{
   int nb,isw,vol,*imb;
   spinor *r;
   spinor_dble *sb,*sm;
   block_t *b;

   b=blk_list(grid,&nb,&isw)+n;

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"assign_s2sdblk [map_s2blk.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(k>=(*b).nsd))
   {
      error_loc(1,1,"assign_s2sdblk [map_s2blk.c]",
                "Block field number is out of range");
      return;
   }

   vol=(*b).vol;
   imb=(*b).imb;
   sb=(*b).sd[k];
   sm=sb;

   if (set==ALL_PTS)
      sm+=vol;
   else if (set==EVEN_PTS)
      sm+=vol/2;
   else if (set==ODD_PTS)
   {
      imb+=vol/2;
      sb+=vol/2;
      sm+=vol;
   }

   for (;sb<sm;sb++)
   {
      r=s+(*imb);

      (*sb).c1.c1.re=(double)((*r).c1.c1.re);
      (*sb).c1.c1.im=(double)((*r).c1.c1.im);
      (*sb).c1.c2.re=(double)((*r).c1.c2.re);
      (*sb).c1.c2.im=(double)((*r).c1.c2.im);
      (*sb).c1.c3.re=(double)((*r).c1.c3.re);
      (*sb).c1.c3.im=(double)((*r).c1.c3.im);

      (*sb).c2.c1.re=(double)((*r).c2.c1.re);
      (*sb).c2.c1.im=(double)((*r).c2.c1.im);
      (*sb).c2.c2.re=(double)((*r).c2.c2.re);
      (*sb).c2.c2.im=(double)((*r).c2.c2.im);
      (*sb).c2.c3.re=(double)((*r).c2.c3.re);
      (*sb).c2.c3.im=(double)((*r).c2.c3.im);

      (*sb).c3.c1.re=(double)((*r).c3.c1.re);
      (*sb).c3.c1.im=(double)((*r).c3.c1.im);
      (*sb).c3.c2.re=(double)((*r).c3.c2.re);
      (*sb).c3.c2.im=(double)((*r).c3.c2.im);
      (*sb).c3.c3.re=(double)((*r).c3.c3.re);
      (*sb).c3.c3.im=(double)((*r).c3.c3.im);

      (*sb).c4.c1.re=(double)((*r).c4.c1.re);
      (*sb).c4.c1.im=(double)((*r).c4.c1.im);
      (*sb).c4.c2.re=(double)((*r).c4.c2.re);
      (*sb).c4.c2.im=(double)((*r).c4.c2.im);
      (*sb).c4.c3.re=(double)((*r).c4.c3.re);
      (*sb).c4.c3.im=(double)((*r).c4.c3.im);      

      imb+=1;
   }
}


void assign_sd2sdblk(blk_grid_t grid,int n,ptset_t set,
                     spinor_dble *sd,int k)
{
   int nb,isw,vol,*imb;
   spinor_dble *sb,*sm;
   block_t *b;

   b=blk_list(grid,&nb,&isw)+n;

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"assign_sd2sdblk [map_s2blk.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(k>=(*b).nsd))
   {
      error_loc(1,1,"assign_sd2sdblk [map_s2blk.c]",
                "Block field number is out of range");
      return;
   }

   vol=(*b).vol;
   imb=(*b).imb;
   sb=(*b).sd[k];
   sm=sb;

   if (set==ALL_PTS)
      sm+=vol;
   else if (set==EVEN_PTS)
      sm+=vol/2;
   else if (set==ODD_PTS)
   {
      imb+=vol/2;
      sb+=vol/2;
      sm+=vol;
   }

   for (;sb<sm;sb++)
   {
      (*sb)=sd[*imb];
      imb+=1;
   }
}


void assign_sdblk2sd(blk_grid_t grid,int n,ptset_t set,
                     int k,spinor_dble *sd)
{
   int nb,isw,vol,*imb;
   spinor_dble *sb,*sm;
   block_t *b;

   b=blk_list(grid,&nb,&isw)+n;

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"assign_sdblk2sd [map_s2blk.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(k>=(*b).nsd))
   {
      error_loc(1,1,"assign_sdblk2sd [map_s2blk.c]",
                "Block field number is out of range");
      return;
   }

   vol=(*b).vol;
   imb=(*b).imb;
   sb=(*b).sd[k];
   sm=sb;

   if (set==ALL_PTS)
      sm+=vol;
   else if (set==EVEN_PTS)
      sm+=vol/2;
   else if (set==ODD_PTS)
   {
      imb+=vol/2;
      sb+=vol/2;
      sm+=vol;
   }

   for (;sb<sm;sb++)
   {
      sd[*imb]=(*sb);
      imb+=1;
   }
}

#endif
