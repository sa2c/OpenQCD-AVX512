
/*******************************************************************************
*
* File utils.c
*
* Copyright (C) 2005, 2008, 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Collection of basic utility programs
*
* The externally accessible functions are
*
*   int safe_mod(int x,int y)
*     Returns x mod y, where y is assumed positive and x can have any
*     sign. The return value is in the interval [0,y)
*
*   void *amalloc(size_t size,int p)
*     Allocates an aligned memory area of "size" bytes, with a starting
*     address (the return value) that is an integer multiple of 2^p. A
*     NULL pointer is returned if the allocation was not successful
*
*   void afree(void *addr)
*     Frees the aligned memory area at address "addr" that was previously
*     allocated using amalloc. If the memory space at this address was
*     already freed using afree, or if the address does not match an
*     address previously returned by amalloc, the program does not do
*     anything
*
*   int mpi_permanent_tag(void)
*     Returns a new send tag that is guaranteed to be unique and which
*     is therefore suitable for use in permanent communication requests.
*     The available number of tags of this kind is 16384
*
*   int mpi_tag(void)
*     Returns a new send tag for use in non-permanent communications.
*     Note that the counter for these tags wraps around after 16384
*     tags have been delivered
*
*   void message(char *format,...)
*     Prints a message from process 0 to stdout. The usage and argument
*     list is the same as in the case of the printf function
*
*******************************************************************************/

#define UTILS_C

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "utils.h"
#include "global.h"

#define MAX_TAG 32767
#define MAX_PERMANENT_TAG MAX_TAG/2

static int pcmn_cnt=-1,cmn_cnt=MAX_TAG;

struct addr_t
{
   char *addr;
   char *true_addr;
   struct addr_t *next;
};

static struct addr_t *first=NULL;


int safe_mod(int x,int y)
{
   if (x>=0)
      return x%y;
   else
      return (y-(abs(x)%y))%y;
}


void *amalloc(size_t size,int p)
{
   int shift;
   char *true_addr,*addr;
   unsigned long mask;
   struct addr_t *new;

   if ((size<=0)||(p<0))
      return(NULL);

   shift=1<<p;
   mask=(unsigned long)(shift-1);

   true_addr=malloc(size+shift);
   new=malloc(sizeof(*first));

   if ((true_addr==NULL)||(new==NULL))
   {
      free(true_addr);
      free(new);
      return(NULL);
   }

   addr=(char*)(((unsigned long)(true_addr+shift))&(~mask));
   (*new).addr=addr;
   (*new).true_addr=true_addr;
   (*new).next=first;
   first=new;

   return (void*)(addr);
}


void afree(void *addr)
{
   struct addr_t *p,*q;

   q=NULL;

   for (p=first;p!=NULL;p=(*p).next)
   {
      if ((*p).addr==addr)
      {
         if (q!=NULL)
            (*q).next=(*p).next;
         else
            first=(*p).next;

         free((*p).true_addr);
         free(p);
         return;
      }

      q=p;
   }
}


int mpi_permanent_tag(void)
{
   if (pcmn_cnt<MAX_PERMANENT_TAG)
      pcmn_cnt+=1;
   else
      error_loc(1,1,"mpi_permanent_tag [utils.c]",
                "Requested more than 16384 tags");

   return pcmn_cnt;
}


int mpi_tag(void)
{
   if (cmn_cnt==MAX_TAG)
      cmn_cnt=MAX_PERMANENT_TAG;

   cmn_cnt+=1;

   return cmn_cnt;
}

void message(char *format,...)
{
   int my_rank;
   va_list args;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      va_start(args,format);
      vprintf(format,args);
      va_end(args);
   }
}
