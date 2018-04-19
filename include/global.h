
/*******************************************************************************
*
* File global.h
*
* Copyright (C) 2009, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Global parameters and arrays
*
*******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

#define NPROC0 2
#define NPROC1 1
#define NPROC2 1
#define NPROC3 1

#define L0 8
#define L1 8
#define L2 8
#define L3 8

#define NPROC0_BLK 1
#define NPROC1_BLK 1
#define NPROC2_BLK 1
#define NPROC3_BLK 1

#define NAME_SIZE 128

/****************************** do not change *********************************/

#if ((NPROC0<1)||(NPROC1<1)||(NPROC2<1)||(NPROC3<1)|| \
    ((NPROC0>1)&&((NPROC0%2)!=0))||((NPROC1>1)&&((NPROC1%2)!=0))|| \
    ((NPROC2>1)&&((NPROC2%2)!=0))||((NPROC3>1)&&((NPROC3%2)!=0)))
#error : The number of processes in each direction must be 1 or a multiple of 2
#endif

#if ((L0<4)||(L1<4)||(L2<4)||(L3<4)|| \
    ((L0%2)!=0)||((L1%2)!=0)||((L2%2)!=0)||((L3%2)!=0))
#error : The local lattice sizes must be even and not smaller than 4
#endif

#if ((NPROC0_BLK<1)||(NBROC0_BLK>NPROC0)||((NPROC0%NPROC0_BLK)!=0)|| \
     (NPROC1_BLK<1)||(NBROC1_BLK>NPROC1)||((NPROC1%NPROC1_BLK)!=0)|| \
     (NPROC2_BLK<1)||(NBROC2_BLK>NPROC2)||((NPROC2%NPROC2_BLK)!=0)|| \
     (NPROC3_BLK<1)||(NBROC3_BLK>NPROC3)||((NPROC3%NPROC3_BLK)!=0))   
#error : Improper processor block sizes NPROC0_BLK,..,NPROC3_BLK
#endif

#if (NAME_SIZE<128)
#error : NAME_SIZE must be greater or equal to 128
#endif

#define NPROC (NPROC0*NPROC1*NPROC2*NPROC3)
#define VOLUME (L0*L1*L2*L3)
#define FACE0 ((1-(NPROC0%2))*L1*L2*L3)
#define FACE1 ((1-(NPROC1%2))*L2*L3*L0)
#define FACE2 ((1-(NPROC2%2))*L3*L0*L1)
#define FACE3 ((1-(NPROC3%2))*L0*L1*L2)
#define BNDRY (2*(FACE0+FACE1+FACE2+FACE3))
#define NSPIN (VOLUME+(BNDRY/2))
#define ALIGN 6

#ifndef SU3_H
#include "su3.h"
#endif

#if defined MAIN_PROGRAM
  #define EXTERN
#else
  #define EXTERN extern
#endif

EXTERN int cpr[4];
EXTERN int npr[8];

EXTERN int ipt[VOLUME];
EXTERN int iup[VOLUME][4];
EXTERN int idn[VOLUME][4];
EXTERN int map[BNDRY+NPROC%2];

#undef EXTERN

#endif
