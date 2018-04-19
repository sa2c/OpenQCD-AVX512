
/*******************************************************************************
*
* File sflds.h
*
* Copyright (C) 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef SFLDS_H
#define SFLDS_H

#ifndef SU3_H
#include "su3.h"
#endif

/* PBND_C */
extern void (*assign_s2w[8])(int *imb,int vol,spinor *s,weyl *r);
extern void (*add_assign_w2s[8])(int *imb,int vol,weyl *s,spinor *r);
extern void (*sub_assign_w2s[8])(int *imb,int vol,weyl *s,spinor *r);
extern void (*mulg5_sub_assign_w2s[8])(int *imb,int vol,weyl *s,spinor *r);

/* PBND_DBLE_C */
extern void (*assign_sd2wd[8])(int *imb,int vol,spinor_dble *sd,
                               weyl_dble *rd);
extern void (*add_assign_wd2sd[8])(int *imb,int vol,weyl_dble *sd,
                                   spinor_dble *rd);
extern void (*sub_assign_wd2sd[8])(int *imb,int vol,weyl_dble *sd,
                                   spinor_dble *rd);
extern void (*mulg5_sub_assign_wd2sd[8])(int *imb,int vol,weyl_dble *sd,
                                         spinor_dble *rd);

/* SFLDS_C */
extern void set_s2zero(int vol,spinor *s);
extern void set_sd2zero(int vol,spinor_dble *sd);
extern void random_s(int vol,spinor *s,float sigma);
extern void random_sd(int vol,spinor_dble *sd,double sigma);
extern void assign_s2s(int vol,spinor *s,spinor *r);
extern void assign_s2sd(int vol,spinor *s,spinor_dble *rd);
extern void assign_sd2s(int vol,spinor_dble *sd,spinor *r);
extern void assign_sd2sd(int vol,spinor_dble *sd,spinor_dble *rd);
extern void diff_s2s(int vol,spinor *s,spinor *r);
extern void add_s2sd(int vol,spinor *s,spinor_dble *rd);
extern void diff_sd2s(int vol,spinor_dble *sd,spinor_dble *rd,spinor *r);

/* SCOM_C */
extern void cps_int_bnd(int is,spinor *s);
extern void cps_ext_bnd(int is,spinor *s);
   
/* SDCOM_C */
extern void cpsd_int_bnd(int is,spinor_dble *sd);
extern void cpsd_ext_bnd(int is,spinor_dble *sd);

#endif


