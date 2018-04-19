
/*******************************************************************************
*
* File dfl.h
*
* Copyright (C) 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef DFL_H
#define DFL_H

#ifndef SU3_H
#include "su3.h"
#endif

typedef struct
{
   int nb,nbb;
   int nbbe[8],nbbo[8];
   int obbe[8],obbo[8];
   int (*inn)[8];
   int *idx,*ipp,*map;
} dfl_grid_t;

/* DFL_GEOMETRY_C */
extern dfl_grid_t dfl_geometry(void);

/* DFL_MODES_C */
extern void dfl_modes(int *status);
extern void dfl_update(int nsm,int *status);
extern void dfl_modes2(int *status);
extern void dfl_update2(int nsm,int *status);

/* DFL_SAP_GCR_C */
extern double dfl_sap_gcr(int nkv,int nmx,double res,double mu,
                          spinor_dble *eta,spinor_dble *psi,int *status);
extern double dfl_sap_gcr2(int nkv,int nmx,double res,double mu,
                           spinor_dble *eta,spinor_dble *psi,int *status);

/* DFL_SUBSPACE_C */
extern void dfl_sd2vd(spinor_dble *sd,complex_dble *vd);
extern void dfl_vd2sd(complex_dble *vd,spinor_dble *sd);
extern void dfl_sub_vd2sd(complex_dble *vd,spinor_dble *sd);
extern void dfl_s2v(spinor *s,complex *v);
extern void dfl_v2s(complex *v,spinor *s);
extern void dfl_sub_v2s(complex *v,spinor *s);
extern void dfl_subspace(spinor **mds);

/* LTL_GCR */
extern double ltl_gcr(int nkv,int nmx,double res,double mu,
                      complex_dble *eta,complex_dble *psi,int *status);

#endif
