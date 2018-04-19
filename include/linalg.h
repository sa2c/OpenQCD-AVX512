
/*******************************************************************************
*
* File linalg.h
*
* Copyright (C) 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef LINALG_H
#define LINALG_H

#ifndef SU3_H
#include "su3.h"
#endif

/* CMATRIX_C */
extern void cmat_vec(int n,complex *a,complex *v,complex *w);
extern void cmat_vec_assign(int n,complex *a,complex *v,complex *w);
extern void cmat_add(int n,complex *a,complex *b,complex *c);
extern void cmat_sub(int n,complex *a,complex *b,complex *c);
extern void cmat_mul(int n,complex *a,complex *b,complex *c);
extern void cmat_dag(int n,complex *a,complex *b);

/* CMATRIX_DBLE_C */
extern void cmat_vec_dble(int n,complex_dble *a,complex_dble *v,
                          complex_dble *w);
extern void cmat_vec_assign_dble(int n,complex_dble *a,complex_dble *v,
                                 complex_dble *w);
extern void cmat_add_dble(int n,complex_dble *a,complex_dble *b,
                          complex_dble *c);
extern void cmat_sub_dble(int n,complex_dble *a,complex_dble *b,
                          complex_dble *c);
extern void cmat_mul_dble(int n,complex_dble *a,complex_dble *b,
                          complex_dble *c);
extern void cmat_dag_dble(int n,complex_dble *a,complex_dble *b);
extern int cmat_inv_dble(int n,complex_dble *a,complex_dble *b,double *k);

/* LIEALG_C */
extern void random_alg(int vol,su3_alg_dble *X);
extern double norm_square_alg(int vol,int icom,su3_alg_dble *X);
extern double scalar_prod_alg(int vol,int icom,su3_alg_dble *X,su3_alg_dble *Y);
extern void set_alg2zero(int vol,su3_alg_dble *X);
extern void set_ualg2zero(int vol,u3_alg_dble *X);
extern void assign_alg2alg(int vol,su3_alg_dble *X,su3_alg_dble *Y);
extern void swap_alg(int vol,su3_alg_dble *X,su3_alg_dble *Y);
extern void muladd_assign_alg(int vol,double r,su3_alg_dble *X,su3_alg_dble *Y);

/* SALG_C */
extern complex spinor_prod(int vol,int icom,spinor *s,spinor *r);
extern float spinor_prod_re(int vol,int icom,spinor *s,spinor *r);
extern float norm_square(int vol,int icom,spinor *s);
extern void mulc_spinor_add(int vol,spinor *s,spinor *r,complex z);
extern void mulr_spinor_add(int vol,spinor *s,spinor *r,float c);
extern void project(int vol,int icom,spinor *s,spinor *r);
extern void scale(int vol,float c,spinor *s);
extern float normalize(int vol,int icom,spinor *s);
extern void rotate(int vol,int n,spinor **ppk,complex *v);
extern void mulg5(int vol,spinor *s);
extern void mulmg5(int vol,spinor *s);

/* SALG_DBLE_C */
extern complex_dble spinor_prod_dble(int vol,int icom,spinor_dble *s,
                                     spinor_dble *r);
extern double spinor_prod_re_dble(int vol,int icom,spinor_dble *s,
                                  spinor_dble *r);
extern complex_dble spinor_prod5_dble(int vol,int icom,spinor_dble *s,
                                      spinor_dble *r);
extern double norm_square_dble(int vol,int icom,spinor_dble *s);
extern void mulc_spinor_add_dble(int vol,spinor_dble *s,spinor_dble *r,
                                 complex_dble z);
extern void mulr_spinor_add_dble(int vol,spinor_dble *s,spinor_dble *r,
                                 double c);
extern void combine_spinor_dble(int vol,spinor_dble *s,spinor_dble *r,
                                double cs,double cr);
extern void project_dble(int vol,int icom,spinor_dble *s,spinor_dble *r);
extern void scale_dble(int vol,double c,spinor_dble *s);
extern double normalize_dble(int vol,int icom,spinor_dble *s);
extern void rotate_dble(int vol,int n,spinor_dble **ppk,complex_dble *v);
extern void mulg5_dble(int vol,spinor_dble *s);
extern void mulmg5_dble(int vol,spinor_dble *s);

/* VALG_C */
extern complex vprod(int n,int icom,complex *v,complex *w);
extern float vnorm_square(int n,int icom,complex *v);
extern void mulc_vadd(int n,complex *v,complex *w,complex z);
extern void vproject(int n,int icom,complex *v,complex *w);
extern void vscale(int n,float r,complex *v);
extern float vnormalize(int n,int icom,complex *v);
extern void vrotate(int n,int nv,complex **pv,complex *a);

/* VALG_DBLE_C */
extern complex_dble vprod_dble(int n,int icom,complex_dble *v,complex_dble *w);
extern double vnorm_square_dble(int n,int icom,complex_dble *v);
extern void mulc_vadd_dble(int n,complex_dble *v,complex_dble *w,
                           complex_dble z);
extern void vproject_dble(int n,int icom,complex_dble *v,complex_dble *w);
extern void vscale_dble(int n,double r,complex_dble *v);
extern double vnormalize_dble(int n,int icom,complex_dble *v);
extern void vrotate_dble(int n,int nv,complex_dble **pv,complex_dble *a);

#endif
