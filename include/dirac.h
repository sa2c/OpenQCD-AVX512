
/*******************************************************************************
*
* File dirac.h
*
* Copyright (C) 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef DIRAC_H
#define DIRAC_H

#ifndef SU3_H
#include "su3.h"
#endif

#ifndef BLOCK_H
#include "block.h"
#endif

/* DW_BND_C */
extern void Dw_bnd(blk_grid_t grid,int n,int k,int l);

/* DW_C */
extern void Dw(float mu,spinor *s,spinor *r);
extern void Dwee(float mu,spinor *s,spinor *r);
extern void Dwoo(float mu,spinor *s,spinor *r);
extern void Dweo(spinor *s,spinor *r);
extern void Dwoe(spinor *s,spinor *r);
extern void Dwhat(float mu,spinor *s,spinor *r);
extern void Dw_blk(blk_grid_t grid,int n,float mu,int k,int l);
extern void Dwee_blk(blk_grid_t grid,int n,float mu,int k,int l);
extern void Dwoo_blk(blk_grid_t grid,int n,float mu,int k,int l);
extern void Dwoe_blk(blk_grid_t grid,int n,int k,int l);
extern void Dweo_blk(blk_grid_t grid,int n,int k,int l);
extern void Dwhat_blk(blk_grid_t grid,int n,float mu,int k,int l);

/* DW_DBLE_C */
extern void Dw_dble(double mu,spinor_dble *s,spinor_dble *r);
extern void Dwee_dble(double mu,spinor_dble *s,spinor_dble *r);
extern void Dwoo_dble(double mu,spinor_dble *s,spinor_dble *r);
extern void Dweo_dble(spinor_dble *s,spinor_dble *r);
extern void Dwoe_dble(spinor_dble *s,spinor_dble *r);
extern void Dwhat_dble(double mu,spinor_dble *s,spinor_dble *r);
extern void Dw_blk_dble(blk_grid_t grid,int n,double mu,int k,int l);
extern void Dwee_blk_dble(blk_grid_t grid,int n,double mu,int k,int l);
extern void Dwoo_blk_dble(blk_grid_t grid,int n,double mu,int k,int l);
extern void Dwoe_blk_dble(blk_grid_t grid,int n,int k,int l);
extern void Dweo_blk_dble(blk_grid_t grid,int n,int k,int l);
extern void Dwhat_blk_dble(blk_grid_t grid,int n,double mu,int k,int l);

#endif
