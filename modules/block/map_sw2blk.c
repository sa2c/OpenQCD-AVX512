
/*******************************************************************************
*
* File map_sw2blk.c
*
* Copyright (C) 2005, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Copying of the SW fields to the blocks in a block grid
*
* The externally accessible functions are
*
*   int assign_swd2swbgr(blk_grid_t grid,ptset_t set)
*     Assigns the global double-precision SW field to the corresponding
*     single-precision fields in the specified grid. On the given point
*     set, the copied Pauli matrices are inverted before assignment and
*     the program returns 0 or 1 depending on whether the inversions were
*     safe or not.
*
*   int assign_swd2swdblk(blk_grid_t grid,int n,ptset_t set)
*     Assigns the global double-precision SW field to the corresponding
*     double-precision field on the n'th block of the specified grid. On
*     the given point set, the copied Pauli matrices are inverted before
*     assignment and the program returns 0 or 1 depending on whether the
*     inversions were safe or not.
*
* Notes:
*
* The possible point sets are defined in utils.h. Independently of the
* specified set, the source field is left unchanged, i.e. the inversions
* are performed "on the fly" (no inversions are performed if set=NO_PTS).
* Pauli matrix inversions are considered unsafe if the Frobenius condition
* number of the matrix exceeds 100 (see pauli_dble.c).
*
* In the case of assign_swd2swbgr(), the single-precision SW field on the
* block grid must not be shared. The program assign_swd2swdblk() instead
* requires that the double-precision SW field on the blocks is shared (an
* error occurs if these conditions are violated).
*
* Both programs in this module may involve communications and must be called
* on all processes simultaneously. The return value is globally the same in
* the case of assign_swd2swbgr(), but may depend on the rank of the process
* in the case of assign_swd2swdblk().
*
*******************************************************************************/

#define MAP_SW2BLK_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "sw_term.h"
#include "block.h"
#include "global.h"

static pauli_dble ms[2] ALIGNED16;


static int cp_swd2sw(block_t *b,ptset_t set)
{
   int *imb,ifail;
   pauli *pb,*pm;
   pauli_dble *swd,*p;

   swd=swdfld();
   pb=(*b).sw;
   pm=pb+(*b).vol;
   imb=(*b).imb;
   ifail=0;

   for (;pb<pm;pb+=2)
   {
      p=swd+2*(*imb);
      imb+=1;

      if ((set==ALL_PTS)||(set==EVEN_PTS))
      {
         ifail|=inv_pauli_dble(0.0,p,ms);
         ifail|=inv_pauli_dble(0.0,p+1,ms+1);
         assign_pauli(2,ms,pb);
      }
      else
         assign_pauli(2,p,pb);
   }

   pm+=(*b).vol;

   for (;pb<pm;pb+=2)
   {
      p=swd+2*(*imb);
      imb+=1;

      if ((set==ALL_PTS)||(set==ODD_PTS))
      {
         ifail|=inv_pauli_dble(0.0,p,ms);
         ifail|=inv_pauli_dble(0.0,p+1,ms+1);
         assign_pauli(2,ms,pb);
      }
      else
         assign_pauli(2,p,pb);
   }

   return ifail;
}


int assign_swd2swbgr(blk_grid_t grid,ptset_t set)
{
   int iprms[2],ie,io;
   int nb,isw,ifail,n;
   block_t *b,*bm;

   if (NPROC>1)
   {
      iprms[0]=(int)(grid);
      iprms[1]=(int)(set);

      MPI_Bcast(iprms,2,MPI_INT,0,MPI_COMM_WORLD);

      error((iprms[0]!=(int)(grid))||(iprms[1]!=(int)(set)),1,
            "assign_swd2swbgr [map_sw2blk.c]","Parameters are not global");
   }

   b=blk_list(grid,&nb,&isw);

   if (nb==0)
   {
      error_root(1,1,"assign_swd2swbgr [map_sw2blk.c]",
                 "Block grid is not allocated");
      return 0;
   }

   if (((*b).sw==NULL)||((*b).shf&0x4))
   {
      error_root(1,1,"assign_swd2swbgr [map_sw2blk.c]",
                 "SW field on the grid is not allocated or shared");
      return 0;
   }

   ie=query_flags(SWD_E_INVERTED);
   io=query_flags(SWD_O_INVERTED);

   error_root(((ie)&&((set==ALL_PTS)||(set==EVEN_PTS)))||
              ((io)&&((set==ALL_PTS)||(set==ODD_PTS))),1,
              "assign_swd2swbgr [map_sw2blk.c]",
              "Attempt to invert the SW field a second time");

   bm=b+nb;
   ifail=0;

   for (;b<bm;b++)
      ifail|=cp_swd2sw(b,set);

   set_grid_flags(grid,ASSIGNED_SWD2SWBGR);

   if ((set==ALL_PTS)||(set==EVEN_PTS))
      set_grid_flags(grid,INVERTED_SW_E);
   if ((set==ALL_PTS)||(set==ODD_PTS))
      set_grid_flags(grid,INVERTED_SW_O);

   if (set!=NO_PTS)
   {
      n=ifail;
      MPI_Reduce(&n,&ifail,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&ifail,1,MPI_INT,0,MPI_COMM_WORLD);
   }

   return ifail;
}


static int cp_swd2swd(block_t *b,ptset_t set)
{
   int *imb,ifail;
   pauli_dble *swd,*pb,*p,*pm;

   swd=swdfld();
   pb=(*b).swd;
   pm=pb+(*b).vol;
   imb=(*b).imb;
   ifail=0;

   for (;pb<pm;pb+=2)
   {
      p=swd+2*(*imb);
      imb+=1;

      if ((set==ALL_PTS)||(set==EVEN_PTS))
      {
         ifail|=inv_pauli_dble(0.0,p,pb);
         ifail|=inv_pauli_dble(0.0,p+1,pb+1);
      }
      else
      {
         pb[0]=p[0];
         pb[1]=p[1];
      }
   }

   pm+=(*b).vol;

   for (;pb<pm;pb+=2)
   {
      p=swd+2*(*imb);
      imb+=1;

      if ((set==ALL_PTS)||(set==ODD_PTS))
      {
         ifail|=inv_pauli_dble(0.0,p,pb);
         ifail|=inv_pauli_dble(0.0,p+1,pb+1);
      }
      else
      {
         pb[0]=p[0];
         pb[1]=p[1];
      }
   }

   return ifail;
}


int assign_swd2swdblk(blk_grid_t grid,int n,ptset_t set)
{
   int nb,isw,ie,io;
   block_t *b;

   b=blk_list(grid,&nb,&isw);

   if ((nb==0)||(n<0)||(n>=nb))
   {
      error_loc(1,1,"assign_swd2swdblk [map_sw2blk.c]",
                "Block grid is not allocated or block number out of range");
      return 0;
   }

   if (((*b).swd==NULL)||(!((*b).shf&0x8)))
   {
      error_loc(1,1,"assign_swd2swdblk [map_sw2blk.c]",
                "Block field is not allocated or not shared");
      return 0;
   }

   ie=query_flags(SWD_E_INVERTED);
   io=query_flags(SWD_O_INVERTED);

   if (((ie)&&((set==ALL_PTS)||(set==EVEN_PTS)))||
       ((io)&&((set==ALL_PTS)||(set==ODD_PTS))))
   {
      error_loc(1,1,"assign_swd2swdblk [map_sw2blk.c]",
                "Attempt to invert the SW field a second time");
      return 0;
   }

   return cp_swd2swd(b+n,set);
}
