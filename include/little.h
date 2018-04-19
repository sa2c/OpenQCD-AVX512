
/*******************************************************************************
*
* File little.h
*
* Copyright (C) 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef LITTLE_H
#define LITTLE_H

#ifndef SU3_H
#include "su3.h"
#endif

typedef struct
{
   int Ns,nb;
   complex **Aee,**Aoo,**Aoe,**Aeo;
} Aw_t;

typedef struct
{
   int Ns,nb;
   complex_dble **Aee,**Aoo,**Aoe,**Aeo;
} Aw_dble_t;

typedef struct
{
   int n[2];
   int vol,ibn;
   spinor_dble **sde[2];
   spinor_dble **sdo[2];
} b2b_flds_t;

/* AW_COM_C */
extern b2b_flds_t *b2b_flds(int n,int mu);
extern void cpAoe_ext_bnd(void);
extern void cpAee_int_bnd(void);

/* AW_C */
extern void Aw(complex *v,complex *w);
extern void Aweeinv(complex *v,complex *w);
extern void Awooinv(complex *v,complex *w);
extern void Awoe(complex *v,complex *w);
extern void Aweo(complex *v,complex *w);
extern void Awhat(complex *v,complex *w);

/* AW_DBLE_C */
extern void Aw_dble(complex_dble *v,complex_dble *w);
extern void Aweeinv_dble(complex_dble *v,complex_dble *w);
extern void Awooinv_dble(complex_dble *v,complex_dble *w);
extern void Awoe_dble(complex_dble *v,complex_dble *w);
extern void Aweo_dble(complex_dble *v,complex_dble *w);
extern void Awhat_dble(complex_dble *v,complex_dble *w);

/* AW_GEN_C */
extern void gather_ud(int vol,int *imb,su3_dble *ud,su3_dble *vd);
extern void gather_sd(int vol,int *imb,spinor_dble *sd,spinor_dble *rd);
extern void apply_u2sd(int vol,int *imb,su3_dble *ud,spinor_dble *sd,
                       spinor_dble *rd);
extern void apply_udag2sd(int vol,int *imb,su3_dble *ud,spinor_dble *sd,
                          spinor_dble *rd);
extern void (*spinor_prod_gamma[])(int vol,spinor_dble *sd,spinor_dble *rd,
                                   complex_dble *sp);

/* AW_OPS_C */
extern Aw_t Awop(void);
extern Aw_t Awophat(void);
extern Aw_dble_t Awop_dble(void);
extern Aw_dble_t Awophat_dble(void);
extern void set_Aw(double mu);
extern int set_Awhat(double mu);

/* LTL_MODES_C */
extern int set_ltl_modes(void);
extern complex_dble *ltl_matrix(void);

#endif
