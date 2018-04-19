
/*******************************************************************************
*
* File random.h
*
* Copyright (C) 2005, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef RANDOM_H
#define RANDOM_H

/* GAUSS_C */
extern void gauss(float r[],int n);
extern void gauss_dble(double r[],int n);

/* RANLUX_C */
extern void start_ranlux(int level,int seed);
extern void export_ranlux(int tag,char *out);
extern int import_ranlux(char *in);

/* RANLXS_C */
extern void ranlxs(float r[],int n);
extern void rlxs_init(int level,int seed);
extern int rlxs_size(void);
extern void rlxs_get(int state[]);
extern void rlxs_reset(int state[]);

/* RANLXD_C */
extern void ranlxd(double r[],int n);
extern void rlxd_init(int level,int seed);
extern int rlxd_size(void);
extern void rlxd_get(int state[]);
extern void rlxd_reset(int state[]);

#endif
