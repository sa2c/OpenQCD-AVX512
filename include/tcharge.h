
/*******************************************************************************
*
* File tcharge.h
*
* Copyright (C) 2010, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef TCHARGE_H
#define TCHARGE_H

#ifndef SU3_H
#include "su3.h"
#endif

/* FTCOM_C */
extern void copy_bnd_ft(int n,u3_alg_dble *ft);
extern void add_bnd_ft(int n,u3_alg_dble *ft);

/* FTENSOR_C */
extern u3_alg_dble **ftensor(void);

/* TCHARGE_C */
extern double tcharge(void);
extern double tcharge_slices(double *qsl);

/* YM_ACTION_C */
extern double ym_action(void);
extern double ym_action_slices(double *asl);

#endif
