
/*******************************************************************************
*
* File Dw_dble.c
*
* Copyright (C) 2005, 2011-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Application of the O(a)-improved Wilson-Dirac operator D (double-
* precision programs).
*
* The externally accessible functions are
*
*   void Dw_dble(double mu,spinor_dble *s,spinor_dble *r)
*     Depending on whether the twisted-mass flag is set or not, this
*     program applies D+i*mu*gamma_5*1e or D+i*mu*gamma_5 to the field
*     s and assigns the result to the field r.
*
*   void Dwee_dble(double mu,spinor_dble *s,spinor_dble *r)
*     Applies D_ee+i*mu*gamma_5 to the field s on the even points of the
*     lattice and assigns the result to the field r.
*
*   void Dwoo_dble(double mu,spinor_dble *s,spinor_dble *r)
*     Depending on whether the twisted-mass flag is set or not, this
*     program applies D_oo or D_oo+i*mu*gamma_5 to the field s on the
*     odd points of the lattice and assigns the result to the field r.
*
*   void Dwoe_dble(spinor_dble *s,spinor_dble *r)
*     Applies D_oe to the field s and assigns the result to the field r.
*
*   void Dweo_dble(spinor_dble *s,spinor_dble *r)
*     Applies D_eo to the field s and *subtracts* the result from the
*     field r.
*
*   void Dwhat_dble(double mu,spinor_dble *s,spinor_dble *r)
*     Applies Dhat+i*mu*gamma_5 to the field s and assigns the result to
*     the field r.
*
* The following programs operate on the fields in the n'th block b of the
* specified block grid:
*
*   void Dw_blk_dble(blk_grid_t grid,int n,double mu,int k,int l)
*     Depending on whether the twisted-mass flag is set or not, this
*     program applies D+i*mu*gamma_5*1e or D+i*mu*gamma_5 to the field
*     b.sd[k] and assigns the result to the field b.sd[l].
*
*   void Dwee_blk_dble(blk_grid_t grid,int n,double mu,int k,int l)
*     Applies D_ee+i*mu*gamma_5 to the field b.sd[k] on the even points and
*     assigns the result to the field b.sd[l].
*
*   void Dwoo_blk_dble(blk_grid_t grid,int n,double mu,int k,int l)
*     Depending on whether the twisted-mass flag is set or not, this
*     program applies D_oo or D_oo+i*mu*gamma_5 to the field b.sd[k] on
*     the odd points and assigns the result to the field b.sd[l].
*
*   void Dwoe_blk_dble(blk_grid_t grid,int n,int k,int l)
*     Applies D_oe to the field b.sd[k] and assigns the result to the field
*     b.sd[l].
*
*   void Dweo_blk_dble(blk_grid_t grid,int n,int k,int l)
*     Applies D_eo to the field b.sd[k] and *subtracts* the result from the
*     field b.sd[l].
*
*   void Dwhat_blk_dble(blk_grid_t grid,int n,double mu,int k,int l)
*     Applies Dhat+i*mu*gamma_5 to the field b.sd[k] and assigns the result
*     to the field b.sd[l].
*
* Notes:
*
* The notation and normalization conventions are specified in the notes
* "Implementation of the lattice Dirac operator" (file doc/dirac.pdf).
*
* In all these programs, it is assumed that the SW term is in the proper
* condition and that the global spinor fields have NSPIN elements. The
* programs check whether the twisted-mass flag (see flags/lat_parms.c) is
* set and turn off the twisted-mass term on the odd lattice sites if it is.
* The input and output fields may not coincide in the case of the programs
* Dw_dble(), Dwhat_dble(), Dw_blk_dble() and Dwhat_blk_dble().
*
* When the input and output fields are different, the input field is not
* changed except possibly at the points at global time 0 and NPROC0*L0-1,
* where both fields are set to zero if so required by the chosen boundary
* conditions. Depending on the operator considered, the fields are zeroed
* only on the even or odd points at these times.
*
* The programs Dw_dble(),..,Dwhat_dble() perform global operations and must
* be called simultaneously on all processes.
*
*******************************************************************************/

#define DW_DBLE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "flags.h"
#include "lattice.h"
#include "uflds.h"
#include "sflds.h"
#include "sw_term.h"
#include "dirac.h"
#include "global.h"

#define N0 (NPROC0*L0)

typedef union
{
   spinor_dble s;
   weyl_dble w[2];
} spin_t;

static double coe,ceo;
static const spinor_dble sd0={{{0.0,0.0},{0.0,0.0},{0.0,0.0}},
                              {{0.0,0.0},{0.0,0.0},{0.0,0.0}},
                              {{0.0,0.0},{0.0,0.0},{0.0,0.0}},
                              {{0.0,0.0},{0.0,0.0},{0.0,0.0}}};
static spin_t rs ALIGNED32;

#if (defined AVX)
#include "avx.h"

#define _load_cst(c) \
__asm__ __volatile__ ("vbroadcastsd %0, %%ymm15 \n\t" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm15")

#define _mul_cst() \
__asm__ __volatile__ ("vmulpd %%ymm15, %%ymm0, %%ymm0 \n\t" \
                      "vmulpd %%ymm15, %%ymm1, %%ymm1 \n\t" \
                      "vmulpd %%ymm15, %%ymm2, %%ymm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

#define _mul_cst_up() \
__asm__ __volatile__ ("vmulpd %%ymm15, %%ymm3, %%ymm3 \n\t" \
                      "vmulpd %%ymm15, %%ymm4, %%ymm4 \n\t" \
                      "vmulpd %%ymm15, %%ymm5, %%ymm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")


static void doe(int *piup,int *pidn,su3_dble *u,spinor_dble *pk)
{
   spinor_dble *sp,*sm;

/******************************* direction +0 *********************************/

   sp=pk+(*(piup++));

   _avx_pair_load_dble((*sp).c1,(*sp).c2);
   _avx_pair_load_up_dble((*sp).c3,(*sp).c4);

   sm=pk+(*(pidn++));
   _prefetch_spinor_dble(sm);

   _avx_vector_add_dble();
   sp=pk+(*(piup++));
   _prefetch_spinor_dble(sp);
   _avx_su3_multiply_pair_dble(*u);

   _avx_weyl_store_up_dble(rs.w[0]);
   _avx_weyl_store_up_dble(rs.w[1]);

/******************************* direction -0 *********************************/

   _avx_pair_load_dble((*sm).c1,(*sm).c2);
   _avx_pair_load_up_dble((*sm).c3,(*sm).c4);

   _avx_vector_sub_dble();
   sm=pk+(*(pidn++));
   _prefetch_spinor_dble(sm);
   u+=1;
   _avx_su3_inverse_multiply_pair_dble(*u);

   _avx_weyl_load_dble(rs.w[0]);
   _avx_vector_add_dble();
   _avx_weyl_store_dble(rs.w[0]);

   _avx_weyl_load_dble(rs.w[1]);
   _avx_vector_sub_dble();
   _avx_weyl_store_dble(rs.w[1]);

/******************************* direction +1 *********************************/

   _avx_pair_load_dble((*sp).c1,(*sp).c2);
   _avx_pair_load_up_dble((*sp).c4,(*sp).c3);

   _avx_vector_i_add_dble();
   sp=pk+(*(piup++));
   _prefetch_spinor_dble(sp);
   u+=1;
   _avx_su3_multiply_pair_dble(*u);

   _avx_weyl_load_dble(rs.w[0]);
   _avx_vector_add_dble();
   _avx_weyl_store_dble(rs.w[0]);

   _avx_weyl_load_dble(rs.w[1]);
   _avx_vector_xch_i_sub_dble();
   _avx_weyl_store_dble(rs.w[1]);

/******************************* direction -1 *********************************/

   _avx_pair_load_dble((*sm).c1,(*sm).c2);
   _avx_pair_load_up_dble((*sm).c4,(*sm).c3);

   _avx_vector_i_sub_dble();
   sm=pk+(*(pidn++));
   _prefetch_spinor_dble(sm);
   u+=1;
   _avx_su3_inverse_multiply_pair_dble(*u);

   _avx_weyl_load_dble(rs.w[0]);
   _avx_vector_add_dble();
   _avx_weyl_store_dble(rs.w[0]);

   _avx_weyl_load_dble(rs.w[1]);
   _avx_vector_xch_i_add_dble();
   _avx_weyl_store_dble(rs.w[1]);

/******************************* direction +2 *********************************/

   _avx_pair_load_dble((*sp).c1,(*sp).c2);
   _avx_pair_load_up_dble((*sp).c4,(*sp).c3);

   _avx_vector_addsub_dble();
   u+=1;
   _avx_su3_multiply_pair_dble(*u);
   sp=pk+(*(piup));
   _prefetch_spinor_dble(sp);
   _avx_weyl_load_dble(rs.w[0]);
   _avx_vector_add_dble();
   _avx_weyl_store_dble(rs.w[0]);

   _avx_weyl_load_dble(rs.w[1]);
   _avx_vector_xch_dble();
   _avx_vector_subadd_dble();
   _avx_weyl_store_dble(rs.w[1]);

/******************************* direction -2 *********************************/

   _avx_pair_load_dble((*sm).c1,(*sm).c2);
   _avx_pair_load_up_dble((*sm).c4,(*sm).c3);

   _avx_vector_subadd_dble();
   sm=pk+(*(pidn));
   _prefetch_spinor_dble(sm);
   u+=1;
   _avx_su3_inverse_multiply_pair_dble(*u);

   _avx_weyl_load_dble(rs.w[0]);
   _avx_vector_add_dble();
   _avx_weyl_store_dble(rs.w[0]);

   _avx_weyl_load_dble(rs.w[1]);
   _avx_vector_xch_dble();
   _avx_vector_addsub_dble();
   _avx_weyl_store_dble(rs.w[1]);

/******************************* direction +3 *********************************/

   _avx_pair_load_dble((*sp).c1,(*sp).c2);
   _avx_pair_load_up_dble((*sp).c3,(*sp).c4);

   _avx_vector_i_addsub_dble();
   u+=1;
   _avx_su3_multiply_pair_dble(*u);

   _avx_weyl_load_dble(rs.w[0]);
   _avx_vector_add_dble();
   _avx_weyl_store_dble(rs.w[0]);

   _avx_weyl_load_dble(rs.w[1]);
   _avx_vector_i_subadd_dble();
   _avx_weyl_store_dble(rs.w[1]);

/******************************* direction -3 *********************************/

   _avx_pair_load_dble((*sm).c1,(*sm).c2);
   _avx_pair_load_up_dble((*sm).c3,(*sm).c4);

   _avx_vector_i_subadd_dble();
   u+=1;
   _avx_su3_inverse_multiply_pair_dble(*u);

   _load_cst(coe);
   _avx_weyl_load_dble(rs.w[0]);
   _avx_vector_add_dble();
   _mul_cst();
   _avx_pair_store_dble(rs.s.c1,rs.s.c2);

   _avx_weyl_load_dble(rs.w[1]);
   _avx_vector_i_addsub_dble();
   _mul_cst();
   _avx_pair_store_dble(rs.s.c3,rs.s.c4);

   _avx_zeroupper();
}


static void deo(int *piup,int *pidn,su3_dble *u,spinor_dble *pl)
{
   spinor_dble *sp,*sm;

/******************************* direction +0 *********************************/

   sp=pl+(*(piup++));
   _prefetch_spinor_dble(sp);

   _load_cst(ceo);
   _avx_pair_load_dble(rs.s.c1,rs.s.c2);
   _avx_pair_load_up_dble(rs.s.c3,rs.s.c4);
   _mul_cst();
   _mul_cst_up();
   _avx_weyl_store_dble(rs.w[0]);
   _avx_weyl_store_up_dble(rs.w[1]);

   sm=pl+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _avx_vector_sub_dble();
   _avx_su3_inverse_multiply_pair_dble(*u);

   _avx_pair_load_dble((*sp).c1,(*sp).c2);
   _avx_vector_add_dble();
   _avx_pair_store_dble((*sp).c1,(*sp).c2);

   _avx_pair_load_dble((*sp).c3,(*sp).c4);
   _avx_vector_sub_dble();
   _avx_pair_store_dble((*sp).c3,(*sp).c4);

/******************************* direction -0 *********************************/

   _avx_weyl_load_dble(rs.w[0]);
   _avx_weyl_load_up_dble(rs.w[1]);

   sp=pl+(*(piup++));
   _prefetch_spinor_dble(sp);
   _avx_vector_add_dble();
   u+=1;
   _avx_su3_multiply_pair_dble(*u);

   _avx_pair_load_dble((*sm).c1,(*sm).c2);
   _avx_vector_add_dble();
   _avx_pair_store_dble((*sm).c1,(*sm).c2);

   _avx_pair_load_dble((*sm).c3,(*sm).c4);
   _avx_vector_add_dble();
   _avx_pair_store_dble((*sm).c3,(*sm).c4);

/******************************* direction +1 *********************************/

   _avx_weyl_load_dble(rs.w[0]);
   _avx_weyl_load_up_dble(rs.w[1]);

   sm=pl+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _avx_vector_xch_i_sub_dble();
   u+=1;
   _avx_su3_inverse_multiply_pair_dble(*u);

   _avx_pair_load_dble((*sp).c1,(*sp).c2);
   _avx_vector_add_dble();
   _avx_pair_store_dble((*sp).c1,(*sp).c2);

   _avx_pair_load_dble((*sp).c3,(*sp).c4);
   _avx_vector_xch_i_add_dble();
   _avx_pair_store_dble((*sp).c3,(*sp).c4);

/******************************* direction -1 *********************************/

   _avx_weyl_load_dble(rs.w[0]);
   _avx_weyl_load_up_dble(rs.w[1]);

   sp=pl+(*(piup++));
   _prefetch_spinor_dble(sp);
   _avx_vector_xch_i_add_dble();
   u+=1;
   _avx_su3_multiply_pair_dble(*u);

   _avx_pair_load_dble((*sm).c1,(*sm).c2);
   _avx_vector_add_dble();
   _avx_pair_store_dble((*sm).c1,(*sm).c2);

   _avx_pair_load_dble((*sm).c3,(*sm).c4);
   _avx_vector_xch_i_sub_dble();
   _avx_pair_store_dble((*sm).c3,(*sm).c4);

/******************************* direction +2 *********************************/

   _avx_weyl_load_dble(rs.w[0]);
   _avx_weyl_load_up_dble(rs.w[1]);

   sm=pl+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _avx_vector_xch_dble();
   _avx_vector_subadd_dble();
   u+=1;
   _avx_su3_inverse_multiply_pair_dble(*u);

   _avx_pair_load_dble((*sp).c1,(*sp).c2);
   _avx_vector_add_dble();
   _avx_pair_store_dble((*sp).c1,(*sp).c2);

   _avx_pair_load_dble((*sp).c3,(*sp).c4);
   _avx_vector_xch_dble();
   _avx_vector_addsub_dble();
   _avx_pair_store_dble((*sp).c3,(*sp).c4);

/******************************* direction -2 *********************************/

   _avx_weyl_load_dble(rs.w[0]);
   _avx_weyl_load_up_dble(rs.w[1]);

   sp=pl+(*(piup));
   _prefetch_spinor_dble(sp);
   _avx_vector_xch_dble();
   _avx_vector_addsub_dble();
   u+=1;
   _avx_su3_multiply_pair_dble(*u);

   _avx_pair_load_dble((*sm).c1,(*sm).c2);
   _avx_vector_add_dble();
   _avx_pair_store_dble((*sm).c1,(*sm).c2);

   _avx_pair_load_dble((*sm).c3,(*sm).c4);
   _avx_vector_xch_dble();
   _avx_vector_subadd_dble();
   _avx_pair_store_dble((*sm).c3,(*sm).c4);

/******************************* direction +3 *********************************/

   _avx_weyl_load_dble(rs.w[0]);
   _avx_weyl_load_up_dble(rs.w[1]);

   sm=pl+(*(pidn));
   _prefetch_spinor_dble(sm);
   _avx_vector_i_subadd_dble();
   u+=1;
   _avx_su3_inverse_multiply_pair_dble(*u);

   _avx_pair_load_dble((*sp).c1,(*sp).c2);
   _avx_vector_add_dble();
   _avx_pair_store_dble((*sp).c1,(*sp).c2);

   _avx_pair_load_dble((*sp).c3,(*sp).c4);
   _avx_vector_i_addsub_dble();
   _avx_pair_store_dble((*sp).c3,(*sp).c4);

/******************************* direction -3 *********************************/

   _avx_weyl_load_dble(rs.w[0]);
   _avx_weyl_load_up_dble(rs.w[1]);

   _avx_vector_i_addsub_dble();
   u+=1;
   _avx_su3_multiply_pair_dble(*u);

   _avx_pair_load_dble((*sm).c1,(*sm).c2);
   _avx_vector_add_dble();
   _avx_pair_store_dble((*sm).c1,(*sm).c2);

   _avx_pair_load_dble((*sm).c3,(*sm).c4);
   _avx_vector_i_subadd_dble();
   _avx_pair_store_dble((*sm).c3,(*sm).c4);

   _avx_zeroupper();
}

#elif (defined x64)
#include "sse2.h"

#define _load_cst(c) \
__asm__ __volatile__ ("movddup %0, %%xmm15" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm15")

#define _mul_cst() \
__asm__ __volatile__ ("mulpd %%xmm15, %%xmm0 \n\t" \
                      "mulpd %%xmm15, %%xmm1 \n\t" \
                      "mulpd %%xmm15, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

#define _mul_cst_up() \
__asm__ __volatile__ ("mulpd %%xmm15, %%xmm3 \n\t" \
                      "mulpd %%xmm15, %%xmm4 \n\t" \
                      "mulpd %%xmm15, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")


static void doe(int *piup,int *pidn,su3_dble *u,spinor_dble *pk)
{
   spinor_dble *sp,*sm;

/******************************* direction +0 *********************************/

   sp=pk+(*(piup++));

   _sse_load_dble((*sp).c1);
   _sse_load_up_dble((*sp).c3);

   sm=pk+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _sse_vector_add_dble();
   _sse_su3_multiply_dble(*u);
   _sse_store_up_dble(rs.s.c1);
   _sse_store_up_dble(rs.s.c3);

   _sse_load_dble((*sp).c2);
   _sse_load_up_dble((*sp).c4);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;

   _sse_vector_add_dble();
   _sse_su3_multiply_dble(*u);

   _sse_store_up_dble(rs.s.c2);
   _sse_store_up_dble(rs.s.c4);

/******************************* direction -0 *********************************/

   _sse_load_dble((*sm).c1);
   _sse_load_up_dble((*sm).c3);

   sp=pk+(*(piup++));
   _prefetch_spinor_dble(sp);
   _sse_vector_sub_dble();
   u+=1;
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble(rs.s.c1);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c1);

   _sse_load_dble(rs.s.c3);
   _sse_vector_sub_dble();
   _sse_store_dble(rs.s.c3);

   _sse_load_dble((*sm).c2);
   _sse_load_up_dble((*sm).c4);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;

   _sse_vector_sub_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble(rs.s.c2);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c2);

   _sse_load_dble(rs.s.c4);
   _sse_vector_sub_dble();
   _sse_store_dble(rs.s.c4);

/******************************* direction +1 *********************************/

   _sse_load_dble((*sp).c1);
   _sse_load_up_dble((*sp).c4);

   sm=pk+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _sse_vector_i_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble(rs.s.c1);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c1);

   _sse_load_dble(rs.s.c4);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble(rs.s.c4);

   _sse_load_dble((*sp).c2);
   _sse_load_up_dble((*sp).c3);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;

   _sse_vector_i_add_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble(rs.s.c2);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c2);

   _sse_load_dble(rs.s.c3);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble(rs.s.c3);

/******************************* direction -1 *********************************/

   _sse_load_dble((*sm).c1);
   _sse_load_up_dble((*sm).c4);

   sp=pk+(*(piup++));
   _prefetch_spinor_dble(sp);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   u+=1;
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble(rs.s.c1);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c1);

   _sse_load_dble(rs.s.c4);
   _sse_vector_i_add_dble();
   _sse_store_dble(rs.s.c4);

   _sse_load_dble((*sm).c2);
   _sse_load_up_dble((*sm).c3);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;

   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble(rs.s.c2);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c2);

   _sse_load_dble(rs.s.c3);
   _sse_vector_i_add_dble();
   _sse_store_dble(rs.s.c3);

/******************************* direction +2 *********************************/

   _sse_load_dble((*sp).c1);
   _sse_load_up_dble((*sp).c4);

   sm=pk+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _sse_vector_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble(rs.s.c1);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c1);

   _sse_load_dble(rs.s.c4);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c4);

   _sse_load_dble((*sp).c2);
   _sse_load_up_dble((*sp).c3);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;

   _sse_vector_sub_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble(rs.s.c2);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c2);

   _sse_load_dble(rs.s.c3);
   _sse_vector_sub_dble();
   _sse_store_dble(rs.s.c3);

/******************************* direction -2 *********************************/

   _sse_load_dble((*sm).c1);
   _sse_load_up_dble((*sm).c4);

   sp=pk+(*(piup));
   _prefetch_spinor_dble(sp);
   _sse_vector_sub_dble();
   u+=1;
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble(rs.s.c1);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c1);

   _sse_load_dble(rs.s.c4);
   _sse_vector_sub_dble();
   _sse_store_dble(rs.s.c4);

   _sse_load_dble((*sm).c2);
   _sse_load_up_dble((*sm).c3);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;

   _sse_vector_add_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble(rs.s.c2);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c2);

   _sse_load_dble(rs.s.c3);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c3);

/******************************* direction +3 *********************************/

   _sse_load_dble((*sp).c1);
   _sse_load_up_dble((*sp).c3);

   sm=pk+(*(pidn));
   _prefetch_spinor_dble(sm);
   _sse_vector_i_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble(rs.s.c1);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c1);

   _sse_load_dble(rs.s.c3);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble(rs.s.c3);

   _sse_load_dble((*sp).c2);
   _sse_load_up_dble((*sp).c4);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;

   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble(rs.s.c2);
   _sse_vector_add_dble();
   _sse_store_dble(rs.s.c2);

   _sse_load_dble(rs.s.c4);
   _sse_vector_i_add_dble();
   _sse_store_dble(rs.s.c4);

/******************************* direction -3 *********************************/

   _sse_load_dble((*sm).c1);
   _sse_load_up_dble((*sm).c3);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   u+=1;
   _sse_su3_inverse_multiply_dble(*u);

   _load_cst(coe);
   _sse_load_dble(rs.s.c1);
   _sse_vector_add_dble();
   _mul_cst();
   _sse_store_dble(rs.s.c1);

   _sse_load_dble(rs.s.c3);
   _sse_vector_i_add_dble();
   _mul_cst();
   _sse_store_dble(rs.s.c3);

   _sse_load_dble((*sm).c2);
   _sse_load_up_dble((*sm).c4);

   u+=1;
   _prefetch_su3_dble(u);
   u-=1;

   _sse_vector_i_add_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _load_cst(coe);
   _sse_load_dble(rs.s.c2);
   _sse_vector_add_dble();
   _mul_cst();
   _sse_store_dble(rs.s.c2);

   _sse_load_dble(rs.s.c4);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _mul_cst();
   _sse_store_dble(rs.s.c4);
}


static void deo(int *piup,int *pidn,su3_dble *u,spinor_dble *pl)
{
   spinor_dble *sp,*sm;

/******************************* direction +0 *********************************/

   sp=pl+(*(piup++));
   _prefetch_spinor_dble(sp);

   _load_cst(ceo);
   _sse_load_dble(rs.s.c1);
   _sse_load_up_dble(rs.s.c3);
   _mul_cst();
   _mul_cst_up();
   _sse_store_dble(rs.s.c1);
   _sse_store_up_dble(rs.s.c3);

   sm=pl+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _sse_vector_sub_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c1);
   _sse_vector_add_dble();
   _sse_store_dble((*sp).c1);

   _sse_load_dble((*sp).c3);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c3);

   _load_cst(ceo);
   _sse_load_dble(rs.s.c2);
   _sse_load_up_dble(rs.s.c4);
   _mul_cst();
   _mul_cst_up();
   _sse_store_dble(rs.s.c2);
   _sse_store_up_dble(rs.s.c4);

   _sse_vector_sub_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c2);
   _sse_vector_add_dble();
   _sse_store_dble((*sp).c2);

   _sse_load_dble((*sp).c4);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c4);

/******************************* direction -0 *********************************/

   _sse_load_dble(rs.s.c1);
   _sse_load_up_dble(rs.s.c3);

   sp=pl+(*(piup++));
   _prefetch_spinor_dble(sp);
   _sse_vector_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c1);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c1);

   _sse_load_dble((*sm).c3);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c3);

   _sse_load_dble(rs.s.c2);
   _sse_load_up_dble(rs.s.c4);

   _sse_vector_add_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c2);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c2);

   _sse_load_dble((*sm).c4);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c4);

/******************************* direction +1 *********************************/

   _sse_load_dble(rs.s.c1);
   _sse_load_up_dble(rs.s.c4);

   sm=pl+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   u+=1;
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c1);
   _sse_vector_add_dble();
   _sse_store_dble((*sp).c1);

   _sse_load_dble((*sp).c4);
   _sse_vector_i_add_dble();
   _sse_store_dble((*sp).c4);

   _sse_load_dble(rs.s.c2);
   _sse_load_up_dble(rs.s.c3);

   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c2);
   _sse_vector_add_dble();
   _sse_store_dble((*sp).c2);

   _sse_load_dble((*sp).c3);
   _sse_vector_i_add_dble();
   _sse_store_dble((*sp).c3);

/******************************* direction -1 *********************************/

   _sse_load_dble(rs.s.c1);
   _sse_load_up_dble(rs.s.c4);

   sp=pl+(*(piup++));
   _prefetch_spinor_dble(sp);
   _sse_vector_i_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c1);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c1);

   _sse_load_dble((*sm).c4);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c4);

   _sse_load_dble(rs.s.c2);
   _sse_load_up_dble(rs.s.c3);

   _sse_vector_i_add_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c2);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c2);

   _sse_load_dble((*sm).c3);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c3);

/******************************* direction +2 *********************************/

   _sse_load_dble(rs.s.c1);
   _sse_load_up_dble(rs.s.c4);

   sm=pl+(*(pidn++));
   _prefetch_spinor_dble(sm);
   _sse_vector_sub_dble();
   u+=1;
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c1);
   _sse_vector_add_dble();
   _sse_store_dble((*sp).c1);

   _sse_load_dble((*sp).c4);
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c4);

   _sse_load_dble(rs.s.c2);
   _sse_load_up_dble(rs.s.c3);

   _sse_vector_add_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c2);
   _sse_vector_add_dble();
   _sse_store_dble((*sp).c2);

   _sse_load_dble((*sp).c3);
   _sse_vector_add_dble();
   _sse_store_dble((*sp).c3);

/******************************* direction -2 *********************************/

   _sse_load_dble(rs.s.c1);
   _sse_load_up_dble(rs.s.c4);

   sp=pl+(*(piup));
   _prefetch_spinor_dble(sp);
   _sse_vector_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c1);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c1);

   _sse_load_dble((*sm).c4);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c4);

   _sse_load_dble(rs.s.c2);
   _sse_load_up_dble(rs.s.c3);

   _sse_vector_sub_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c2);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c2);

   _sse_load_dble((*sm).c3);
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c3);

/******************************* direction +3 *********************************/

   _sse_load_dble(rs.s.c1);
   _sse_load_up_dble(rs.s.c3);

   sm=pl+(*(pidn));
   _prefetch_spinor_dble(sm);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   u+=1;
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c1);
   _sse_vector_add_dble();
   _sse_store_dble((*sp).c1);

   _sse_load_dble((*sp).c3);
   _sse_vector_i_add_dble();
   _sse_store_dble((*sp).c3);

   _sse_load_dble(rs.s.c2);
   _sse_load_up_dble(rs.s.c4);

   _sse_vector_i_add_dble();
   _sse_su3_inverse_multiply_dble(*u);

   _sse_load_dble((*sp).c2);
   _sse_vector_add_dble();
   _sse_store_dble((*sp).c2);

   _sse_load_dble((*sp).c4);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble((*sp).c4);

/******************************* direction -3 *********************************/

   _sse_load_dble(rs.s.c1);
   _sse_load_up_dble(rs.s.c3);

   _sse_vector_i_add_dble();
   u+=1;
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c1);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c1);

   _sse_load_dble((*sm).c3);
   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_store_dble((*sm).c3);

   _sse_load_dble(rs.s.c2);
   _sse_load_up_dble(rs.s.c4);

   _sse_vector_i_mul_dble();
   _sse_vector_sub_dble();
   _sse_su3_multiply_dble(*u);

   _sse_load_dble((*sm).c2);
   _sse_vector_add_dble();
   _sse_store_dble((*sm).c2);

   _sse_load_dble((*sm).c4);
   _sse_vector_i_add_dble();
   _sse_store_dble((*sm).c4);
}

#else

#define _vector_mul_assign(r,c) \
   (r).c1.re*=(c); \
   (r).c1.im*=(c); \
   (r).c2.re*=(c); \
   (r).c2.im*=(c); \
   (r).c3.re*=(c); \
   (r).c3.im*=(c)


static void doe(int *piup,int *pidn,su3_dble *u,spinor_dble *pk)
{
   spinor_dble *sp,*sm;
   su3_vector_dble psi,chi;

/******************************* direction +0 *********************************/

   sp=pk+(*(piup++));

   _vector_add(psi,(*sp).c1,(*sp).c3);
   _su3_multiply(rs.s.c1,*u,psi);
   rs.s.c3=rs.s.c1;

   _vector_add(psi,(*sp).c2,(*sp).c4);
   _su3_multiply(rs.s.c2,*u,psi);
   rs.s.c4=rs.s.c2;

/******************************* direction -0 *********************************/

   sm=pk+(*(pidn++));
   u+=1;

   _vector_sub(psi,(*sm).c1,(*sm).c3);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c1,chi);
   _vector_sub_assign(rs.s.c3,chi);

   _vector_sub(psi,(*sm).c2,(*sm).c4);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c2,chi);
   _vector_sub_assign(rs.s.c4,chi);

/******************************* direction +1 *********************************/

   sp=pk+(*(piup++));
   u+=1;

   _vector_i_add(psi,(*sp).c1,(*sp).c4);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c1,chi);
   _vector_i_sub_assign(rs.s.c4,chi);

   _vector_i_add(psi,(*sp).c2,(*sp).c3);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c2,chi);
   _vector_i_sub_assign(rs.s.c3,chi);

/******************************* direction -1 *********************************/

   sm=pk+(*(pidn++));
   u+=1;

   _vector_i_sub(psi,(*sm).c1,(*sm).c4);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c1,chi);
   _vector_i_add_assign(rs.s.c4,chi);

   _vector_i_sub(psi,(*sm).c2,(*sm).c3);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c2,chi);
   _vector_i_add_assign(rs.s.c3,chi);

/******************************* direction +2 *********************************/

   sp=pk+(*(piup++));
   u+=1;

   _vector_add(psi,(*sp).c1,(*sp).c4);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c1,chi);
   _vector_add_assign(rs.s.c4,chi);

   _vector_sub(psi,(*sp).c2,(*sp).c3);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c2,chi);
   _vector_sub_assign(rs.s.c3,chi);

/******************************* direction -2 *********************************/

   sm=pk+(*(pidn++));
   u+=1;

   _vector_sub(psi,(*sm).c1,(*sm).c4);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c1,chi);
   _vector_sub_assign(rs.s.c4,chi);

   _vector_add(psi,(*sm).c2,(*sm).c3);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c2,chi);
   _vector_add_assign(rs.s.c3,chi);

/******************************* direction +3 *********************************/

   sp=pk+(*(piup));
   u+=1;

   _vector_i_add(psi,(*sp).c1,(*sp).c3);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c1,chi);
   _vector_i_sub_assign(rs.s.c3,chi);

   _vector_i_sub(psi,(*sp).c2,(*sp).c4);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c2,chi);
   _vector_i_add_assign(rs.s.c4,chi);

/******************************* direction -3 *********************************/

   sm=pk+(*(pidn));
   u+=1;

   _vector_i_sub(psi,(*sm).c1,(*sm).c3);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c1,chi);
   _vector_i_add_assign(rs.s.c3,chi);

   _vector_i_add(psi,(*sm).c2,(*sm).c4);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign(rs.s.c2,chi);
   _vector_i_sub_assign(rs.s.c4,chi);

   _vector_mul_assign(rs.s.c1,coe);
   _vector_mul_assign(rs.s.c2,coe);
   _vector_mul_assign(rs.s.c3,coe);
   _vector_mul_assign(rs.s.c4,coe);
}


static void deo(int *piup,int *pidn,su3_dble *u,spinor_dble *pl)
{
   spinor_dble *sp,*sm;
   su3_vector_dble psi,chi;

   _vector_mul_assign(rs.s.c1,ceo);
   _vector_mul_assign(rs.s.c2,ceo);
   _vector_mul_assign(rs.s.c3,ceo);
   _vector_mul_assign(rs.s.c4,ceo);

/******************************* direction +0 *********************************/

   sp=pl+(*(piup++));

   _vector_sub(psi,rs.s.c1,rs.s.c3);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign((*sp).c1,chi);
   _vector_sub_assign((*sp).c3,chi);

   _vector_sub(psi,rs.s.c2,rs.s.c4);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign((*sp).c2,chi);
   _vector_sub_assign((*sp).c4,chi);

/******************************* direction -0 *********************************/

   sm=pl+(*(pidn++));
   u+=1;

   _vector_add(psi,rs.s.c1,rs.s.c3);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign((*sm).c1,chi);
   _vector_add_assign((*sm).c3,chi);

   _vector_add(psi,rs.s.c2,rs.s.c4);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign((*sm).c2,chi);
   _vector_add_assign((*sm).c4,chi);

/******************************* direction +1 *********************************/

   sp=pl+(*(piup++));
   u+=1;

   _vector_i_sub(psi,rs.s.c1,rs.s.c4);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign((*sp).c1,chi);
   _vector_i_add_assign((*sp).c4,chi);

   _vector_i_sub(psi,rs.s.c2,rs.s.c3);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign((*sp).c2,chi);
   _vector_i_add_assign((*sp).c3,chi);

/******************************* direction -1 *********************************/

   sm=pl+(*(pidn++));
   u+=1;

   _vector_i_add(psi,rs.s.c1,rs.s.c4);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign((*sm).c1,chi);
   _vector_i_sub_assign((*sm).c4,chi);

   _vector_i_add(psi,rs.s.c2,rs.s.c3);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign((*sm).c2,chi);
   _vector_i_sub_assign((*sm).c3,chi);

/******************************* direction +2 *********************************/

   sp=pl+(*(piup++));
   u+=1;

   _vector_sub(psi,rs.s.c1,rs.s.c4);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign((*sp).c1,chi);
   _vector_sub_assign((*sp).c4,chi);

   _vector_add(psi,rs.s.c2,rs.s.c3);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign((*sp).c2,chi);
   _vector_add_assign((*sp).c3,chi);

/******************************* direction -2 *********************************/

   sm=pl+(*(pidn++));
   u+=1;

   _vector_add(psi,rs.s.c1,rs.s.c4);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign((*sm).c1,chi);
   _vector_add_assign((*sm).c4,chi);

   _vector_sub(psi,rs.s.c2,rs.s.c3);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign((*sm).c2,chi);
   _vector_sub_assign((*sm).c3,chi);

/******************************* direction +3 *********************************/

   sp=pl+(*(piup));
   u+=1;

   _vector_i_sub(psi,rs.s.c1,rs.s.c3);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign((*sp).c1,chi);
   _vector_i_add_assign((*sp).c3,chi);

   _vector_i_add(psi,rs.s.c2,rs.s.c4);
   _su3_inverse_multiply(chi,*u,psi);
   _vector_add_assign((*sp).c2,chi);
   _vector_i_sub_assign((*sp).c4,chi);

/******************************* direction -3 *********************************/

   sm=pl+(*(pidn));
   u+=1;

   _vector_i_add(psi,rs.s.c1,rs.s.c3);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign((*sm).c1,chi);
   _vector_i_sub_assign((*sm).c3,chi);

   _vector_i_sub(psi,rs.s.c2,rs.s.c4);
   _su3_multiply(chi,*u,psi);
   _vector_add_assign((*sm).c2,chi);
   _vector_i_add_assign((*sm).c4,chi);
}

#endif

void Dw_dble(double mu,spinor_dble *s,spinor_dble *r)
{
   int bc,ix,t;
   int *piup,*pidn;
   su3_dble *u,*um;
   pauli_dble *m;
   spin_t *so,*ro;
   tm_parms_t tm;

   cpsd_int_bnd(0x1,s);
   m=swdfld();
   apply_sw_dble(VOLUME/2,mu,m,s,r);
   set_sd2zero(BNDRY/2,r+VOLUME);
   tm=tm_parms();
   if (tm.eoflg==1)
      mu=0.0;

   coe=-0.5;
   ceo=-0.5;
   bc=bc_type();
   piup=iup[VOLUME/2];
   pidn=idn[VOLUME/2];

   so=(spin_t*)(s+(VOLUME/2));
   ro=(spin_t*)(r+(VOLUME/2));
   m+=VOLUME;
   u=udfld();
   um=u+4*VOLUME;

   if (((cpr[0]==0)&&(bc!=3))||((cpr[0]==(NPROC0-1))&&(bc==0)))
   {
      ix=VOLUME/2;

      for (;u<um;u+=8)
      {
         t=global_time(ix);
         ix+=1;

         if ((t>0)&&((t<(N0-1))||(bc!=0)))
         {
            doe(piup,pidn,u,s);

            mul_pauli_dble(mu,m,(*so).w,(*ro).w);
            mul_pauli_dble(-mu,m+1,(*so).w+1,(*ro).w+1);

            _vector_add_assign((*ro).s.c1,rs.s.c1);
            _vector_add_assign((*ro).s.c2,rs.s.c2);
            _vector_add_assign((*ro).s.c3,rs.s.c3);
            _vector_add_assign((*ro).s.c4,rs.s.c4);
            rs=(*so);

            deo(piup,pidn,u,r);
         }
         else
         {
            (*so).s=sd0;
            (*ro).s=sd0;
         }

         piup+=4;
         pidn+=4;
         so+=1;
         ro+=1;
         m+=2;
      }
   }
   else
   {
      for (;u<um;u+=8)
      {
         doe(piup,pidn,u,s);

         mul_pauli_dble(mu,m,(*so).w,(*ro).w);
         mul_pauli_dble(-mu,m+1,(*so).w+1,(*ro).w+1);

         _vector_add_assign((*ro).s.c1,rs.s.c1);
         _vector_add_assign((*ro).s.c2,rs.s.c2);
         _vector_add_assign((*ro).s.c3,rs.s.c3);
         _vector_add_assign((*ro).s.c4,rs.s.c4);
         rs=(*so);

         deo(piup,pidn,u,r);

         piup+=4;
         pidn+=4;
         so+=1;
         ro+=1;
         m+=2;
      }
   }

   cpsd_ext_bnd(0x1,r);
}


void Dwee_dble(double mu,spinor_dble *s,spinor_dble *r)
{
   int bc,ix,t;
   pauli_dble *m,*mm;
   spin_t *se,*re;

   bc=bc_type();
   m=swdfld();
   mm=m+VOLUME;
   se=(spin_t*)(s);
   re=(spin_t*)(r);

   if (((cpr[0]==0)&&(bc!=3))||((cpr[0]==(NPROC0-1))&&(bc==0)))
   {
      ix=0;

      for (;m<mm;m+=2)
      {
         t=global_time(ix);
         ix+=1;

         if ((t>0)&&((t<(N0-1))||(bc!=0)))
         {
            mul_pauli_dble(mu,m,(*se).w,(*re).w);
            mul_pauli_dble(-mu,m+1,(*se).w+1,(*re).w+1);
         }
         else
         {
            (*se).s=sd0;
            (*re).s=sd0;
         }

         se+=1;
         re+=1;
      }
   }
   else
   {
      for (;m<mm;m+=2)
      {
         mul_pauli_dble(mu,m,(*se).w,(*re).w);
         mul_pauli_dble(-mu,m+1,(*se).w+1,(*re).w+1);

         se+=1;
         re+=1;
      }
   }
}


void Dwoo_dble(double mu,spinor_dble *s,spinor_dble *r)
{
   int bc,ix,t;
   pauli_dble *m,*mm;
   spin_t *so,*ro;
   tm_parms_t tm;

   bc=bc_type();
   m=swdfld()+VOLUME;
   mm=m+VOLUME;
   so=(spin_t*)(s+(VOLUME/2));
   ro=(spin_t*)(r+(VOLUME/2));
   tm=tm_parms();
   if (tm.eoflg==1)
      mu=0.0;

   if (((cpr[0]==0)&&(bc!=3))||((cpr[0]==(NPROC0-1))&&(bc==0)))
   {
      ix=VOLUME/2;

      for (;m<mm;m+=2)
      {
         t=global_time(ix);
         ix+=1;

         if ((t>0)&&((t<(N0-1))||(bc!=0)))
         {
            mul_pauli_dble(mu,m,(*so).w,(*ro).w);
            mul_pauli_dble(-mu,m+1,(*so).w+1,(*ro).w+1);
         }
         else
         {
            (*so).s=sd0;
            (*ro).s=sd0;
         }

         so+=1;
         ro+=1;
      }
   }
   else
   {
      for (;m<mm;m+=2)
      {
         mul_pauli_dble(mu,m,(*so).w,(*ro).w);
         mul_pauli_dble(-mu,m+1,(*so).w+1,(*ro).w+1);

         so+=1;
         ro+=1;
      }
   }
}


void Dwoe_dble(spinor_dble *s,spinor_dble *r)
{
   int bc,ix,t;
   int *piup,*pidn;
   su3_dble *u,*um;
   spin_t *ro;

   cpsd_int_bnd(0x1,s);

   coe=-0.5;
   bc=bc_type();
   piup=iup[VOLUME/2];
   pidn=idn[VOLUME/2];

   ro=(spin_t*)(r+(VOLUME/2));
   u=udfld();
   um=u+4*VOLUME;

   if (((cpr[0]==0)&&(bc!=3))||((cpr[0]==(NPROC0-1))&&(bc==0)))
   {
      ix=VOLUME/2;

      for (;u<um;u+=8)
      {
         t=global_time(ix);
         ix+=1;

         if ((t>0)&&((t<(N0-1))||(bc!=0)))
         {
            doe(piup,pidn,u,s);
            (*ro)=rs;
         }
         else
            (*ro).s=sd0;

         piup+=4;
         pidn+=4;
         ro+=1;
      }
   }
   else
   {
      for (;u<um;u+=8)
      {
         doe(piup,pidn,u,s);
         (*ro)=rs;

         piup+=4;
         pidn+=4;
         ro+=1;
      }
   }
}


void Dweo_dble(spinor_dble *s,spinor_dble *r)
{
   int bc,ix,t;
   int *piup,*pidn;
   su3_dble *u,*um;
   spin_t *so;

   set_sd2zero(BNDRY/2,r+VOLUME);

   ceo=0.5;
   bc=bc_type();
   piup=iup[VOLUME/2];
   pidn=idn[VOLUME/2];

   so=(spin_t*)(s+(VOLUME/2));
   u=udfld();
   um=u+4*VOLUME;

   if (((cpr[0]==0)&&(bc!=3))||((cpr[0]==(NPROC0-1))&&(bc==0)))
   {
      ix=VOLUME/2;

      for (;u<um;u+=8)
      {
         t=global_time(ix);
         ix+=1;

         if ((t>0)&&((t<(N0-1))||(bc!=0)))
         {
            rs=(*so);
            deo(piup,pidn,u,r);
         }
         else
            (*so).s=sd0;

         piup+=4;
         pidn+=4;
         so+=1;
      }
   }
   else
   {
      for (;u<um;u+=8)
      {
         rs=(*so);
         deo(piup,pidn,u,r);

         piup+=4;
         pidn+=4;
         so+=1;
      }
   }

   cpsd_ext_bnd(0x1,r);
}


void Dwhat_dble(double mu,spinor_dble *s,spinor_dble *r)
{
   int bc,ix,t;
   int *piup,*pidn;
   su3_dble *u,*um;
   pauli_dble *m;

   cpsd_int_bnd(0x1,s);
   m=swdfld();
   apply_sw_dble(VOLUME/2,mu,m,s,r);
   set_sd2zero(BNDRY/2,r+VOLUME);

   coe=-0.5;
   ceo=0.5;
   bc=bc_type();
   piup=iup[VOLUME/2];
   pidn=idn[VOLUME/2];

   m+=VOLUME;
   u=udfld();
   um=u+4*VOLUME;

   if (((cpr[0]==0)&&(bc!=3))||((cpr[0]==(NPROC0-1))&&(bc==0)))
   {
      ix=VOLUME/2;

      for (;u<um;u+=8)
      {
         t=global_time(ix);
         ix+=1;

         if ((t>0)&&((t<(N0-1))||(bc!=0)))
         {
            doe(piup,pidn,u,s);

            mul_pauli_dble(0.0,m,rs.w,rs.w);
            mul_pauli_dble(0.0,m+1,rs.w+1,rs.w+1);

            deo(piup,pidn,u,r);
         }

         piup+=4;
         pidn+=4;
         m+=2;
      }
   }
   else
   {
      for (;u<um;u+=8)
      {
         doe(piup,pidn,u,s);

         mul_pauli_dble(0.0,m,rs.w,rs.w);
         mul_pauli_dble(0.0,m+1,rs.w+1,rs.w+1);

         deo(piup,pidn,u,r);

         piup+=4;
         pidn+=4;
         m+=2;
      }
   }

   cpsd_ext_bnd(0x1,r);
}


void Dw_blk_dble(blk_grid_t grid,int n,double mu,int k,int l)
{
   int nb,iswd,vol,volh,ibu,ibd;
   int *piup,*pidn,*ibp,*ibm;
   su3_dble *u,*um;
   pauli_dble *m;
   spinor_dble *s,*r;
   spin_t *so,*ro;
   block_t *b;
   tm_parms_t tm;

   b=blk_list(grid,&nb,&iswd);

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"Dw_blk_dble [Dw_dble.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(l<0)||(k==l)||(k>=(*b).nsd)||(l>=(*b).nsd)||((*b).ud==NULL))
   {
      error_loc(1,1,"Dw_blk_dbl [Dw_dble.c]",
                "Attempt to access unallocated memory space");
      return;
   }

   b+=n;
   vol=(*b).vol;
   volh=vol/2;
   s=(*b).sd[k];
   r=(*b).sd[l];
   so=(spin_t*)(s+volh);
   ro=(spin_t*)(r+volh);

   s[vol]=sd0;
   r[vol]=sd0;
   m=(*b).swd;
   apply_sw_dble(volh,mu,m,s,r);
   tm=tm_parms();
   if (tm.eoflg==1)
      mu=0.0;

   coe=-0.5;
   ceo=-0.5;
   piup=(*b).iup[volh];
   pidn=(*b).idn[volh];
   m+=vol;
   u=(*b).ud;
   um=u+4*vol;

   if ((*b).nbp)
   {
      ibp=(*b).ibp;
      ibm=ibp+(*b).nbp/2;

      for (;ibp<ibm;ibp++)
         s[*ibp]=sd0;

      ibu=((cpr[0]==(NPROC0-1))&&(((*b).bo[0]+(*b).bs[0])==L0)&&(bc_type()==0));
      ibd=((cpr[0]==0)&&((*b).bo[0]==0));

      for (;u<um;u+=8)
      {
         if (((pidn[0]<vol)||(!ibd))&&((piup[0]<vol)||(!ibu)))
         {
            doe(piup,pidn,u,s);

            mul_pauli_dble(mu,m,(*so).w,(*ro).w);
            mul_pauli_dble(-mu,m+1,(*so).w+1,(*ro).w+1);

            _vector_add_assign((*ro).s.c1,rs.s.c1);
            _vector_add_assign((*ro).s.c2,rs.s.c2);
            _vector_add_assign((*ro).s.c3,rs.s.c3);
            _vector_add_assign((*ro).s.c4,rs.s.c4);
            rs=(*so);

            deo(piup,pidn,u,r);
         }
         else
         {
            (*so).s=sd0;
            (*ro).s=sd0;
         }

         piup+=4;
         pidn+=4;
         ro+=1;
         so+=1;
         m+=2;
      }

      ibp=(*b).ibp;
      ibm=ibp+(*b).nbp/2;

      for (;ibp<ibm;ibp++)
         r[*ibp]=sd0;
   }
   else
   {
      for (;u<um;u+=8)
      {
         doe(piup,pidn,u,s);

         mul_pauli_dble(mu,m,(*so).w,(*ro).w);
         mul_pauli_dble(-mu,m+1,(*so).w+1,(*ro).w+1);

         _vector_add_assign((*ro).s.c1,rs.s.c1);
         _vector_add_assign((*ro).s.c2,rs.s.c2);
         _vector_add_assign((*ro).s.c3,rs.s.c3);
         _vector_add_assign((*ro).s.c4,rs.s.c4);
         rs=(*so);

         deo(piup,pidn,u,r);

         piup+=4;
         pidn+=4;
         ro+=1;
         so+=1;
         m+=2;
      }
   }
}


void Dwee_blk_dble(blk_grid_t grid,int n,double mu,int k,int l)
{
   int nb,iswd,vol,ibu,ibd;
   int *piup,*pidn;
   pauli_dble *m,*mm;
   spin_t *se,*re;
   block_t *b;

   b=blk_list(grid,&nb,&iswd);

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"Dwee_blk_dble [Dw_dble.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(l<0)||(k>=(*b).nsd)||(l>=(*b).nsd)||((*b).ud==NULL))
   {
      error_loc(1,1,"Dwee_blk_dbl [Dw_dble.c]",
                "Attempt to access unallocated memory space");
      return;
   }

   b+=n;
   vol=(*b).vol;
   se=(spin_t*)((*b).sd[k]);
   re=(spin_t*)((*b).sd[l]);
   m=(*b).swd;
   mm=m+vol;

   if ((*b).nbp)
   {
      piup=(*b).iup[0];
      pidn=(*b).idn[0];

      ibu=((cpr[0]==(NPROC0-1))&&(((*b).bo[0]+(*b).bs[0])==L0)&&(bc_type()==0));
      ibd=((cpr[0]==0)&&((*b).bo[0]==0));

      for (;m<mm;m+=2)
      {
         if (((pidn[0]<vol)||(!ibd))&&((piup[0]<vol)||(!ibu)))
         {
            mul_pauli_dble(mu,m,(*se).w,(*re).w);
            mul_pauli_dble(-mu,m+1,(*se).w+1,(*re).w+1);
         }
         else
         {
            (*se).s=sd0;
            (*re).s=sd0;
         }

         piup+=4;
         pidn+=4;
         se+=1;
         re+=1;
      }
   }
   else
   {
      for (;m<mm;m+=2)
      {
         mul_pauli_dble(mu,m,(*se).w,(*re).w);
         mul_pauli_dble(-mu,m+1,(*se).w+1,(*re).w+1);

         se+=1;
         re+=1;
      }
   }
}


void Dwoo_blk_dble(blk_grid_t grid,int n,double mu,int k,int l)
{
   int nb,iswd,vol,volh,ibu,ibd;
   int *piup,*pidn;
   pauli_dble *m,*mm;
   spin_t *so,*ro;
   block_t *b;
   tm_parms_t tm;

   b=blk_list(grid,&nb,&iswd);

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"Dwoo_blk_dble [Dw_dble.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(l<0)||(k>=(*b).nsd)||(l>=(*b).nsd)||((*b).ud==NULL))
   {
      error_loc(1,1,"Dwoo_blk_dbl [Dw_dble.c]",
                "Attempt to access unallocated memory space");
      return;
   }

   b+=n;
   vol=(*b).vol;
   volh=vol/2;
   so=(spin_t*)((*b).sd[k]+volh);
   ro=(spin_t*)((*b).sd[l]+volh);
   tm=tm_parms();
   if (tm.eoflg==1)
      mu=0.0;

   m=(*b).swd+vol;
   mm=m+vol;

   if ((*b).nbp)
   {
      piup=(*b).iup[volh];
      pidn=(*b).idn[volh];

      ibu=((cpr[0]==(NPROC0-1))&&(((*b).bo[0]+(*b).bs[0])==L0)&&(bc_type()==0));
      ibd=((cpr[0]==0)&&((*b).bo[0]==0));

      for (;m<mm;m+=2)
      {
         if (((pidn[0]<vol)||(!ibd))&&((piup[0]<vol)||(!ibu)))
         {
            mul_pauli_dble(mu,m,(*so).w,(*ro).w);
            mul_pauli_dble(-mu,m+1,(*so).w+1,(*ro).w+1);
         }
         else
         {
            (*so).s=sd0;
            (*ro).s=sd0;
         }

         piup+=4;
         pidn+=4;
         so+=1;
         ro+=1;
      }
   }
   else
   {
      for (;m<mm;m+=2)
      {
         mul_pauli_dble(mu,m,(*so).w,(*ro).w);
         mul_pauli_dble(-mu,m+1,(*so).w+1,(*ro).w+1);

         so+=1;
         ro+=1;
      }
   }
}


void Dwoe_blk_dble(blk_grid_t grid,int n,int k,int l)
{
   int nb,iswd,vol,volh,ibu,ibd;
   int *piup,*pidn,*ibp,*ibm;
   su3_dble *u,*um;
   spinor_dble *s;
   spin_t *ro;
   block_t *b;

   b=blk_list(grid,&nb,&iswd);

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"Dwoe_blk_dble [Dw_dble.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(l<0)||(k>=(*b).nsd)||(l>=(*b).nsd)||((*b).ud==NULL))
   {
      error_loc(1,1,"Dwoe_blk_dbl [Dw_dble.c]",
                "Attempt to access unallocated memory space");
      return;
   }

   b+=n;
   vol=(*b).vol;
   volh=vol/2;
   s=(*b).sd[k];
   ro=(spin_t*)((*b).sd[l]+volh);
   s[vol]=sd0;

   coe=-0.5;
   piup=(*b).iup[volh];
   pidn=(*b).idn[volh];
   u=(*b).ud;
   um=u+4*vol;

   if ((*b).nbp)
   {
      ibp=(*b).ibp;
      ibm=ibp+(*b).nbp/2;

      for (;ibp<ibm;ibp++)
         s[*ibp]=sd0;

      ibu=((cpr[0]==(NPROC0-1))&&(((*b).bo[0]+(*b).bs[0])==L0)&&(bc_type()==0));
      ibd=((cpr[0]==0)&&((*b).bo[0]==0));

      for (;u<um;u+=8)
      {
         if (((pidn[0]<vol)||(!ibd))&&((piup[0]<vol)||(!ibu)))
         {
            doe(piup,pidn,u,s);
            (*ro)=rs;
         }
         else
            (*ro).s=sd0;

         piup+=4;
         pidn+=4;
         ro+=1;
      }
   }
   else
   {
      for (;u<um;u+=8)
      {
         doe(piup,pidn,u,s);
         (*ro)=rs;

         piup+=4;
         pidn+=4;
         ro+=1;
      }
   }
}


void Dweo_blk_dble(blk_grid_t grid,int n,int k,int l)
{
   int nb,iswd,vol,volh,ibu,ibd;
   int *piup,*pidn,*ibp,*ibm;
   su3_dble *u,*um;
   spinor_dble *r;
   spin_t *so;
   block_t *b;

   b=blk_list(grid,&nb,&iswd);

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"Dweo_blk_dble [Dw_dble.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(l<0)||(k>=(*b).nsd)||(l>=(*b).nsd)||((*b).ud==NULL))
   {
      error_loc(1,1,"Dweo_blk_dbl [Dw_dble.c]",
                "Attempt to access unallocated memory space");
      return;
   }

   b+=n;
   vol=(*b).vol;
   volh=vol/2;
   so=(spin_t*)((*b).sd[k]+volh);
   r=(*b).sd[l];
   r[vol]=sd0;

   ceo=0.5;
   piup=(*b).iup[volh];
   pidn=(*b).idn[volh];
   u=(*b).ud;
   um=u+4*vol;

   if ((*b).nbp)
   {
      ibu=((cpr[0]==(NPROC0-1))&&(((*b).bo[0]+(*b).bs[0])==L0)&&(bc_type()==0));
      ibd=((cpr[0]==0)&&((*b).bo[0]==0));

      for (;u<um;u+=8)
      {
         if (((pidn[0]<vol)||(!ibd))&&((piup[0]<vol)||(!ibu)))
         {
            rs=(*so);
            deo(piup,pidn,u,r);
         }
         else
            (*so).s=sd0;

         piup+=4;
         pidn+=4;
         so+=1;
      }

      ibp=(*b).ibp;
      ibm=ibp+(*b).nbp/2;

      for (;ibp<ibm;ibp++)
         r[*ibp]=sd0;
   }
   else
   {
      for (;u<um;u+=8)
      {
         rs=(*so);
         deo(piup,pidn,u,r);

         piup+=4;
         pidn+=4;
         so+=1;
      }
   }
}


void Dwhat_blk_dble(blk_grid_t grid,int n,double mu,int k,int l)
{
   int nb,iswd,vol,volh,ibu,ibd;
   int *piup,*pidn,*ibp,*ibm;
   su3_dble *u,*um;
   pauli_dble *m;
   spinor_dble *s,*r;
   block_t *b;

   b=blk_list(grid,&nb,&iswd);

   if ((n<0)||(n>=nb))
   {
      error_loc(1,1,"Dwhat_blk_dble [Dw_dble.c]",
                "Block grid is not allocated or block number out of range");
      return;
   }

   if ((k<0)||(l<0)||(k==l)||(k>=(*b).nsd)||(l>=(*b).nsd)||((*b).ud==NULL))
   {
      error_loc(1,1,"Dwhat_blk_dbl [Dw_dble.c]",
                "Attempt to access unallocated memory space");
      return;
   }

   b+=n;
   vol=(*b).vol;
   volh=vol/2;
   s=(*b).sd[k];
   r=(*b).sd[l];

   s[vol]=sd0;
   r[vol]=sd0;
   m=(*b).swd;
   apply_sw_dble(volh,mu,m,s,r);

   coe=-0.5;
   ceo=0.5;
   piup=(*b).iup[volh];
   pidn=(*b).idn[volh];
   m+=vol;
   u=(*b).ud;
   um=u+4*vol;

   if ((*b).nbp)
   {
      ibp=(*b).ibp;
      ibm=ibp+(*b).nbp/2;

      for (;ibp<ibm;ibp++)
         s[*ibp]=sd0;

      ibu=((cpr[0]==(NPROC0-1))&&(((*b).bo[0]+(*b).bs[0])==L0)&&(bc_type()==0));
      ibd=((cpr[0]==0)&&((*b).bo[0]==0));

      for (;u<um;u+=8)
      {
         if (((pidn[0]<vol)||(!ibd))&&((piup[0]<vol)||(!ibu)))
         {
            doe(piup,pidn,u,s);

            mul_pauli_dble(0.0f,m,rs.w,rs.w);
            mul_pauli_dble(0.0f,m+1,rs.w+1,rs.w+1);

            deo(piup,pidn,u,r);
         }

         piup+=4;
         pidn+=4;
         m+=2;
      }

      ibp=(*b).ibp;
      ibm=ibp+(*b).nbp/2;

      for (;ibp<ibm;ibp++)
         r[*ibp]=sd0;
   }
   else
   {
      for (;u<um;u+=8)
      {
         doe(piup,pidn,u,s);

         mul_pauli_dble(0.0f,m,rs.w,rs.w);
         mul_pauli_dble(0.0f,m+1,rs.w+1,rs.w+1);

         deo(piup,pidn,u,r);

         piup+=4;
         pidn+=4;
         m+=2;
      }
   }
}
