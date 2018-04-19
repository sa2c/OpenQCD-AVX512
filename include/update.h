
/*******************************************************************************
*
* File update.h
*
* Copyright (C) 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef UPDATE_H
#define UPDATE_H

#ifndef SU3_H
#include "su3.h"
#endif

typedef struct
{
   int iop;
   double eps;
} mdstep_t;

/* CHRONO */
extern void setup_chrono(void);
extern double mdtime(void);
extern void step_mdtime(double dt);
extern void add_chrono(int icr,spinor_dble *psi);
extern int get_chrono(int icr,spinor_dble *psi);
extern void reset_chrono(void);

/* COUNTERS */
extern void setup_counters(void);
extern void clear_counters(void);
extern void add2counter(char *type,int idx,int *status);
extern int get_count(char *type,int idx,int *status);
extern void print_avgstat(char *type,int idx);
extern void print_all_avgstat(void);

/* MDSTEPS_C */
extern void set_mdsteps(void);
extern mdstep_t *mdsteps(int *nop,int *itu);
extern void print_mdsteps(int ipr);

/* MDINT_C */
extern void run_mdint(void);

/* HMC_C */
extern void hmc_sanity_check(void);
extern void hmc_wsize(int *nwud,int *nws,int *nwsd,int *nwv,int *nwvd);
extern int run_hmc(double *act0,double *act1);

/* RWRAT_C */
extern double rwrat(int irp,int n,int *np,int *isp,double *sqn,int **status);

/* RWTM_C */
extern double rwtm1(double mu1,double mu2,int isp,double *sqn,int *status);
extern double rwtm2(double mu1,double mu2,int isp,double *sqn,int *status);

/* RWTMEO_C */
extern double rwtm1eo(double mu1,double mu2,int isp,double *sqn,int *status);
extern double rwtm2eo(double mu1,double mu2,int isp,double *sqn,int *status);

#endif
