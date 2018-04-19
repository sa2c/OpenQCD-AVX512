
/*******************************************************************************
*
* File linsolv.h
*
* Copyright (C) 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef LINSOLV_H
#define LINSOLV_H

#ifndef SU3_H
#include "su3.h"
#endif

/* CGNE_C */
extern double cgne(int vol,int icom,void (*Dop)(spinor *s,spinor *r),
                   void (*Dop_dble)(spinor_dble *s,spinor_dble *r),
                   spinor **ws,spinor_dble **wsd,int nmx,double res,
                   spinor_dble *eta,spinor_dble *psi,int *status);

/* FGCR4VD_C */
extern double fgcr4vd(int vol,int icom,
                      void (*Dop)(complex_dble *v,complex_dble *w),
                      void (*Mop)(int k,complex *eta,complex *psi,complex *chi),
                      complex **wv,complex_dble **wvd,int nkv,int nmx,double res,
                      complex_dble *eta,complex_dble *psi,int *status);

/* FGCR_C */
extern double fgcr(int vol,int icom,
                   void (*Dop)(spinor_dble *s,spinor_dble *r),
                   void (*Mop)(int k,spinor *rho,spinor *phi,spinor *chi),
                   spinor **ws,spinor_dble **wsd,int nkv,int nmx,double res,
                   spinor_dble *eta,spinor_dble *psi,int *status);

/* MSCG_C */
extern void mscg(int vol,int icom,int nmu,double *mu,
                 void (*Dop_dble)(double mu,spinor_dble *s,spinor_dble *r),
                 spinor_dble **wsd,int nmx,double *res,
                 spinor_dble *eta,spinor_dble **psi,int *status);

#endif
