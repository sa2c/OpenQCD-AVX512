
/*******************************************************************************
*
* File block.h
*
* Copyright (C) 2005, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef BLOCK_H
#define BLOCK_H

#ifndef SU3_H
#include "su3.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

typedef struct
{
   int ifc,ibn,vol,nw,nwd;
   int *ipp,*map,*imb;
   su3 *u;
   su3_dble *ud;
   weyl **w;
   weyl_dble **wd;
} bndry_t;

typedef struct
{
   int *bo,*bs,vol,vbb,nbp,ns,nsd,shf;
   int *ipt,*imb,*ibp;
   int (*iup)[4],(*idn)[4];
   su3 *u;
   su3_dble *ud;
   pauli *sw;
   pauli_dble *swd;
   spinor **s;
   spinor_dble **sd;
   bndry_t *bb;
} block_t;

typedef enum
{
   SAP_BLOCKS,DFL_BLOCKS,
   BLK_GRIDS
} blk_grid_t;

/* BLOCK_C */
extern void alloc_blk(block_t *b,int *bo,int *bs,
                      int iu,int iud,int ns,int nsd);
extern void alloc_bnd(block_t *b,int iu,int iud,int nw,int nwd);
extern void clone_blk(block_t *b,int shf,int *bo,block_t *c);
extern void free_blk(block_t *b);
extern int ipt_blk(block_t *b,int *x);

/* BLK_GRID_C */
extern void alloc_bgr(blk_grid_t grid);
extern block_t *blk_list(blk_grid_t grid,int *nb,int *isw);

/* MAP_U2BLK_C */
extern void assign_ud2ubgr(blk_grid_t grid);
extern void assign_ud2udblk(blk_grid_t grid,int n);

/* MAP_SW2BLK_C */
extern int assign_swd2swbgr(blk_grid_t grid,ptset_t set);
extern int assign_swd2swdblk(blk_grid_t grid,int n,ptset_t set);

/* MAP_S2BLK_C */
extern void assign_s2sblk(blk_grid_t grid,int n,ptset_t set,spinor *s,int k);
extern void assign_sblk2s(blk_grid_t grid,int n,ptset_t set,int k,spinor *s);
extern void assign_s2sdblk(blk_grid_t grid,int n,ptset_t set,spinor *s,int k);
extern void assign_sd2sdblk(blk_grid_t grid,int n,ptset_t set,
                            spinor_dble *sd,int k);
extern void assign_sdblk2sd(blk_grid_t grid,int n,ptset_t set,
                            int k,spinor_dble *sd);

#endif
