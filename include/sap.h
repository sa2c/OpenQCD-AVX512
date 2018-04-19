
/*******************************************************************************
*
* File sap.h
*
* Copyright (C) 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef SAP_H
#define SAP_H

#ifndef SU3_H
#include "su3.h"
#endif

/* BLK_SOLV_C */
extern void blk_mres(int n,float mu,int nmr);
extern void blk_eo_mres(int n,float mu,int nmr);

/* SAP_COM_C */
#if ((defined SAP_COM_C)||(defined BLK_GRID_C ))
extern void alloc_sap_bufs(void);
#endif
extern void sap_com(int ic,spinor *r);

/* SAP */
extern void sap(float mu,int isolv,int nmr,spinor *psi,spinor *eta);

/* SAP_GCR */
extern double sap_gcr(int nkv,int nmx,double res,double mu,
                      spinor_dble *eta,spinor_dble *psi,int *status);

#endif
