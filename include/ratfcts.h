
/*******************************************************************************
*
* File ratfcts.h
*
* Copyright (C) 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef RATFCTS_H
#define RATFCTS_H

typedef struct
{
   int np;
   double A,delta;
   double *mu,*rmu;
   double *nu,*rnu;
} ratfct_t;

/* ELLIPTIC_C */
extern double ellipticK(double rk);
extern void sncndn(double u,double rk,double *sn,double *cn,double *dn);

/* RATFCTS_C */
extern ratfct_t ratfct(int *irat);

/* ZOLOTAREV_C */
extern void zolotarev(int n,double eps,double *A,double *ar,double *delta);

#endif
