
/*******************************************************************************
*
* File utils.c
*
* Copyright (C) 2005, 2008, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)

* Collection of basic utility programs
*
* The externally accessible functions are
*
*   int safe_mod(int x,int y)
*     Returns x mod y, where y is assumed positive and x can have any
*     sign. The return value is in the interval [0,y)
*
*   void *amalloc(size_t size,int p)
*     Allocates an aligned memory area of "size" bytes, with starting
*     address (the return value) that is an integer multiple of 2^p
*
*   void afree(void *addr)
*     Frees the aligned memory area at address "addr" that was
*     previously allocated using amalloc
*
*   void error(int test,int no,char *name,char *format,...)
*     Checks whether "test"=0 and if not aborts the program gracefully
*     with error number "no" after printing the "name" of the calling
*     program and an error message to stdout. The message is formed using
*     the "format" string and any additional arguments, exactly as in a
*     printf statement
*
*   void error_root(int test,int no,char *name,char *format,...)
*     Same as error(), provided for compatibility
*
*   void error_loc(int test,int no,char *name,char *format,...)
*     Same as error(), provided for compatibility
*
*   void message(char *format,...)
*     Same as printf(), provided for compatibility
*
*******************************************************************************/

#define UTILS_C

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include "utils.h"

struct addr_t
{
   char *addr;
   char *true_addr;
   struct addr_t *last,*next;
};

static struct addr_t *rpos=NULL;


int safe_mod(int x,int y)
{
   if (x>=0)
      return(x%y);
   else
      return((y-(abs(x)%y))%y);
}


void *amalloc(size_t size,int p)
{
   int shift;
   char *true_addr,*addr;
   unsigned long mask;
   struct addr_t *new,*rnxt;

   if ((size<=0)||(p<0))
      return(NULL);

   shift=1<<p;
   mask=(unsigned long)(shift-1);

   true_addr=malloc(size+shift);
   new=malloc(sizeof(*new));

   if ((true_addr==NULL)||(new==NULL))
   {
      free(true_addr);
      free(new);
      return NULL;
   }

   addr=(char*)(((unsigned long)(true_addr+shift))&(~mask));
   (*new).addr=addr;
   (*new).true_addr=true_addr;

   if (rpos!=NULL)
   {
      rnxt=(*rpos).next;

      (*new).next=rnxt;
      (*rpos).next=new;
      (*rnxt).last=new;
      (*new).last=rpos;
   }
   else
   {
      (*new).next=new;
      (*new).last=new;
   }

   rpos=new;

   return (void*)(addr);
}


void afree(void *addr)
{
   struct addr_t *p,*pn,*pl;

   if (rpos!=NULL)
   {
      p=rpos;

      for (;;)
      {
         if ((*p).addr==addr)
         {
            pn=(*p).next;
            pl=(*p).last;

            if (pn!=p)
            {
               (*pl).next=pn;
               (*pn).last=pl;
               rpos=pl;
            }
            else
               rpos=NULL;

            free((*p).true_addr);
            free(p);
            return;
         }

         p=(*p).next;
         if (p==rpos)
            return;
      }
   }
}


void error(int test,int no,char *name,char *format,...)
{
   va_list args;

   if (test!=0)
   {
      printf("\nError in %s:\n",name);
      va_start(args,format);
      vprintf(format,args);
      va_end(args);
      printf("\nProgram aborted\n\n");
      exit(no);
   }
}


void error_root(int test,int no,char *name,char *format,...)
{
   va_list args;

   if (test!=0)
   {
      printf("\nError in %s:\n",name);
      va_start(args,format);
      vprintf(format,args);
      va_end(args);
      printf("\nProgram aborted\n\n");
      exit(no);
   }
}


void error_loc(int test,int no,char *name,char *format,...)
{
   va_list args;

   if (test!=0)
   {
      printf("\nError in %s:\n",name);
      va_start(args,format);
      vprintf(format,args);
      va_end(args);
      printf("\nProgram aborted\n\n");
      exit(no);
   }
}


void message(char *format,...)
{
   va_list args;

   va_start(args,format);
   vprintf(format,args);
   va_end(args);
}
