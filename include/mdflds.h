
/*******************************************************************************
*
* File mdflds.h
*
* Copyright (C) 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef MDFLDS_H
#define MDFLDS_H

#ifndef SU3_H
#include "su3.h"
#endif

typedef struct
{
   int npf;
   su3_alg_dble *mom,*frc;
   spinor_dble **pf;
} mdflds_t;

/* FCOM_C */
extern void copy_bnd_frc(void);
extern void add_bnd_frc(void);

/* MDFLDS_C */
extern mdflds_t *mdflds(void);
extern void set_frc2zero(void);
extern void bnd_mom2zero(void);
extern void random_mom(void);
extern double momentum_action(int icom);

#endif
