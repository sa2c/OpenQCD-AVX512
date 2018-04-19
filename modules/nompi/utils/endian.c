
/*******************************************************************************
*
* File endian.c
*
* Copyright (C) 2007, 2009, 2010 Bjoern Leder, Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Byte swapping programs
*
* The externally accessible functions are
*
*   int endianness(void)
*     Returns LITTLE_ENDIAN if the machine is little endian and BIG_ENDIAN
*     if it is big endian. Otherwise the return value is UNKNOWN_ENDIAN
*
*   void bswap_int(int n,void *a)
*     Inverts the byte order of the array elements a[0],..,a[n-1]
*     assuming these are 4 byte long 
*
*   void bswap_double(int n,void *a)
*     Inverts the byte order of the array elements a[0],..,a[n-1]
*     assuming these are 8 byte long 
*
* Notes:
*
* The integer types that are guaranteed to be 4 byte long are stdint_t
* and stduint_t. These are defined in the header file misc.h.
*
* On machines complying with the IEEE-754 standard, double precision
* floating-point numbers are 8 byte long. When the header file misc.h
* is read, the compiler checks whether the machine complies with the
* standard.

* The programs in this module do not involve any communications and can
* be called locally.
*
*******************************************************************************/

#define ENDIAN_C

#include <stdlib.h>
#include <stdio.h>
#include "utils.h"


int endianness(void)
{
   stduint_t i;
   unsigned char *b;

   i=0x04030201;
   b=(unsigned char*)(&i);

   if ((b[0]==1u)&&(b[1]==2u)&&(b[2]==3u)&&(b[3]==4u))
      return LITTLE_ENDIAN;
   else if ((b[0]==4u)&&(b[1]==3u)&&(b[2]==2u)&&(b[3]==1u))
      return BIG_ENDIAN;
   else return UNKNOWN_ENDIAN;
}


void bswap_int(int n,void *a)
{
   unsigned char *ba,*bam,bas;

   ba=(unsigned char*)(a);
   bam=ba+4*n;

   for (;ba<bam;ba+=4)
   {
      bas=ba[3];
      ba[3]=ba[0];
      ba[0]=bas;

      bas=ba[2];
      ba[2]=ba[1];
      ba[1]=bas;      
   }
}


void bswap_double(int n,void *a)
{
   unsigned char *ba,*bam,bas;

   ba=(unsigned char*)(a);
   bam=ba+8*n;

   for (;ba<bam;ba+=8)
   {
      bas=ba[7];
      ba[7]=ba[0];
      ba[0]=bas;

      bas=ba[6];
      ba[6]=ba[1];
      ba[1]=bas;      

      bas=ba[5];
      ba[5]=ba[2];
      ba[2]=bas;      

      bas=ba[4];
      ba[4]=ba[3];
      ba[3]=bas;
   }
}
