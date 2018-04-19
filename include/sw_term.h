
/*******************************************************************************
*
* File sw_term.h
*
* Copyright (C) 2005, 2009, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef SW_TERM_H
#define SW_TERM_H

#ifndef SU3_H
#include "su3.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

/* PAULI_C */
extern void mul_pauli(float mu,pauli *m,weyl *s,weyl *r);
extern void mul_pauli2(float mu,pauli *m,spinor *s,spinor *r);
extern void assign_pauli(int vol,pauli_dble *md,pauli *m);
extern void apply_sw(int vol,float mu,pauli *m,spinor *s,spinor *r);

/* PAULI_DBLE_C */
extern void mul_pauli_dble(double mu,pauli_dble *m,weyl_dble *s,weyl_dble *r);
extern int inv_pauli_dble(double mu,pauli_dble *m,pauli_dble *im);
extern complex_dble det_pauli_dble(double mu,pauli_dble *m);
extern void apply_sw_dble(int vol,double mu,pauli_dble *m,spinor_dble *s,
                          spinor_dble *r);
extern int apply_swinv_dble(int vol,double mu,pauli_dble *m,spinor_dble *s,
                            spinor_dble *r);

/* SWFLDS_C */
extern pauli *swfld(void);
extern pauli_dble *swdfld(void);
extern void assign_swd2sw(void);

/* SW_TERM_C */
extern int sw_term(ptset_t set);

#endif
