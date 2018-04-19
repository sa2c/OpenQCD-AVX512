
/*******************************************************************************
*
* File archive.h
*
* Copyright (C) 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef ARCHIVE_H
#define ARCHIVE_H

#ifndef SU3_H
#include "su3.h"
#endif

/* ARCHIVE_C */
extern void write_cnfg(char *out);
extern void read_cnfg(char *in);
extern void export_cnfg(char *out);
extern void import_cnfg(char *in);

/* MARCHIVE_C */
extern void write_mfld(char *out);
extern void read_mfld(char *in);
extern void export_mfld(char *out);
extern void import_mfld(char *in);

/* SARCHIVE_C */
extern void write_sfld(char *out,spinor_dble *sd);
extern void read_sfld(char *in,spinor_dble *sd);
extern void export_sfld(char *out,spinor_dble *sd);
extern void import_sfld(char *in,spinor_dble *sd);

#endif
