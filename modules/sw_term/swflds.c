
/*******************************************************************************
*
* File swflds.c
*
* Copyright (C) 2006, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Allocation and initialization of the global SW fields
*
* The externally accessible functions are
*
*   pauli *swfld(void)
*     Returns the base address of the single-precision SW field. If it
*     is not already allocated, the field is allocated and initialized
*     to unity.
*
*   pauli_dble *swdfld(void)
*     Returns the base address of the double-precision SW field. If it
*     is not already allocated, the field is allocated and initialized
*     to unity.
*
*   void assign_swd2sw(void)
*     Assigns the double-precision to the single-precision SW field.
*
* Notes:
*
* All these programs act globally and must be called simultaneously on all
* processes.
*
*******************************************************************************/

#define SWFLDS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "sw_term.h"
#include "global.h"

static pauli *swb=NULL;
static pauli_dble *swdb=NULL;


static void alloc_sw(void)
{
   int i;
   pauli *sw,*sm,unity;

   error_root(sizeof(pauli)!=(36*sizeof(float)),1,"alloc_sw [swflds.c]",
              "The pauli structures are not properly packed");

   swb=amalloc(2*VOLUME*sizeof(*swb),ALIGN);
   error(swb==NULL,1,"alloc_sw [swflds.c]",
         "Unable to allocate the global single-precision SW field");

   unity.u[0]=1.0f;
   unity.u[1]=1.0f;
   unity.u[2]=1.0f;
   unity.u[3]=1.0f;
   unity.u[4]=1.0f;
   unity.u[5]=1.0f;

   for (i=6;i<36;i++)
      unity.u[i]=0.0f;

   sw=swb;
   sm=sw+2*VOLUME;

   for (;sw<sm;sw++)
      (*sw)=unity;
}


pauli *swfld(void)
{
   if (swb==NULL)
      alloc_sw();

   return swb;
}


static void alloc_swd(void)
{
   int i;
   pauli_dble *sw,*sm,unity;

   error_root(sizeof(pauli_dble)!=(36*sizeof(double)),1,"alloc_swd [swflds.c]",
              "The pauli_dble structures are not properly packed");

   swdb=amalloc(2*VOLUME*sizeof(*swdb),ALIGN);
   error(swdb==NULL,1,"alloc_swd [swflds.c]",
         "Unable to allocate the global double-precision SW field");

   unity.u[0]=1.0;
   unity.u[1]=1.0;
   unity.u[2]=1.0;
   unity.u[3]=1.0;
   unity.u[4]=1.0;
   unity.u[5]=1.0;

   for (i=6;i<36;i++)
      unity.u[i]=0.0;

   sw=swdb;
   sm=sw+2*VOLUME;

   for (;sw<sm;sw++)
      (*sw)=unity;
}


pauli_dble *swdfld(void)
{
   if (swdb==NULL)
      alloc_swd();

   return swdb;
}


void assign_swd2sw(void)
{
   error(swdb==NULL,1,"assign_swd2sw [swflds.c]",
          "Attempt to access unallocated memory space");

   if (swb==NULL)
      alloc_sw();

   assign_pauli(2*VOLUME,swdb,swb);
   set_flags(ASSIGNED_SWD2SW);
}
