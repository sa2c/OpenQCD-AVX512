/*******************************************************************************
*
* File wflow.h
*
* Copyright (C) 2009, 2010, 2011, 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef WFLOW_H
#define WFLOW_H

/* WFLOW_C */
extern void fwd_euler(int n,double eps);
extern void fwd_rk2(int n,double eps);
extern void fwd_rk3(int n,double eps);

#endif
