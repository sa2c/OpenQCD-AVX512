
/*******************************************************************************
*
* File nompi/extras.h
*
* Copyright (C) 2009, 2010, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef EXTRAS_H 
#define EXTRAS_H

/* CHEBYSHEV_C */
extern int cheby_fit(double a,double b,double (*f)(double x),
                     int dmax,double eps,double c[]);
extern double cheby_int(double a,double b,double (*f)(double x),
                        int dmax,double eps);
extern double cheby_val(double a,double b,int n,double c[],double x);

/* FSOLVE_C */
extern double inverse_fct(double x1,double x2,double (*f)(double x),double y,
			  double omega1,double omega2);
extern double minimize_fct(double x0,double x1,double x2,double (*f)(double x),
			   double omega1,double omega2);
extern void powell(int n,double *x0,double *x1,double *x2,
                   double (*f)(int n,double *x),int imx,double omega1,
                   double omega2,double *xmin,int *status);

/* I0M_C */
extern double i0m(double x); 

/* KS_TEST_C */
extern void ks_test(int n,double f[],double *pkp,double *pkm);
extern void ks_prob(int n,double kp,double km,double *pp,double *pm);

/* PCHI_SQUARE_C */
extern double pchi_square(double chi_square,int nu);

/* STAT_C */
extern double average(int n,double *a);
extern double sigma0(int n,double *a);
extern double auto_corr(int n,double *a,int tmax,double *g);
extern void sigma_auto_corr(int n,double *a,int tmax,int lambda,double *eg);
extern double tauint(int n,double *a,int tmax,int lambda,int *w,double *sigma);
extern double print_auto(int n,double *a);
extern double jack_err(int nx,int n,double **a,double (*f)(int nx,double *x),
                       int bmax,double *sig);
extern double print_jack(int nx,int n,double **a,double (*f)(int nx,double *x));

#endif
