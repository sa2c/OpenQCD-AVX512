
/*******************************************************************************
*
* File sse.h
*
* Copyright (C) 2005, 2008, 2009, 2011,  Martin Luescher, Filippo Palombi
*               2016
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Macros for Dirac spinors, SU(3) vectors and SU(3) matrices using inline
* assembly SSE3 instructions. The machine is assumed to comply with the
* x86-64 instruction set.
*
*******************************************************************************/

#ifndef SSE_H
#define SSE_H

typedef struct __attribute__ ((aligned (16)))
{
   int c1,c2,c3,c4;
} sse_int;

typedef struct __attribute__ ((aligned (16)))
{
   float c1,c2,c3,c4;
} sse_float;

typedef struct __attribute__ ((aligned (16)))
{
   sse_float c1,c2,c3;
} sse_vector;

static sse_float _sse_sgn12 __attribute__ ((unused)) ={-1.0f,-1.0f,1.0f,1.0f};
static sse_float _sse_sgn13 __attribute__ ((unused)) ={-1.0f,1.0f,-1.0f,1.0f};
static sse_float _sse_sgn14 __attribute__ ((unused)) ={-1.0f,1.0f,1.0f,-1.0f};
static sse_float _sse_sgn23 __attribute__ ((unused)) ={1.0f,-1.0f,-1.0f,1.0f};
static sse_float _sse_sgn24 __attribute__ ((unused)) ={1.0f,-1.0f,1.0f,-1.0f};
static sse_float _sse_sgn34 __attribute__ ((unused)) ={1.0f,1.0f,-1.0f,-1.0f};
static sse_float _sse_sgn   __attribute__ ((unused)) ={-1.0f,-1.0f,-1.0f,-1.0f};

/*******************************************************************************
*
* Prefetch macros
*
*******************************************************************************/

#if (defined P4)

#define _pfbase(addr) ((unsigned long)(addr)&(~0x7fL))

#define _prefetch_128b(addr) \
__asm__ __volatile__ ("prefetcht0 %0" \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))))

#define _prefetch_256b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1" \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L))))

#define _prefetch_384b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2" \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L))), \
                      "m" (*((char*)(_pfbase(addr)+0x100L))))

#define _prefetch_su3_alg_dble(addr) \
_prefetch_128b((addr))

#define _prefetch_weyl(addr) \
_prefetch_256b((addr))

#define _prefetch_spinor(addr) \
_prefetch_256b((addr))

#define _prefetch_su3(addr) \
_prefetch_256b((addr))

#define _prefetch_pauli(addr) \
_prefetch_256b((addr))

#define _prefetch_weyl_dble(addr) \
_prefetch_256b((addr))

#define _prefetch_spinor_dble(addr) \
_prefetch_256b((addr))

#define _prefetch_su3_dble(addr) \
_prefetch_256b((addr))

#define _prefetch_pauli_dble(addr) \
_prefetch_384b((addr))

#elif (defined PM)

#define _pfbase(addr) ((unsigned long)(addr)&(~0x3fL))

#define _prefetch_64b(addr) \
__asm__ __volatile__ ("prefetcht0 %0" \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))))

#define _prefetch_128b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1" \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))))

#define _prefetch_192b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2" \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L))))

#define _prefetch_320b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2 \n\t" \
                      "prefetcht0 %3 \n\t" \
                      "prefetcht0 %4" \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L))), \
                      "m" (*((char*)(_pfbase(addr)+0xc0L))), \
                      "m" (*((char*)(_pfbase(addr)+0x100L))))

#define _prefetch_su3_alg_dble(addr) \
_prefetch_64b((addr))

#define _prefetch_weyl(addr) \
_prefetch_64b((addr))

#define _prefetch_spinor(addr) \
_prefetch_128b((addr))

#define _prefetch_su3(addr) \
_prefetch_128b((addr))

#define _prefetch_pauli(addr) \
_prefetch_192b((addr))

#define _prefetch_weyl_dble(addr) \
_prefetch_128b((addr))

#define _prefetch_spinor_dble(addr) \
_prefetch_192b((addr))

#define _prefetch_su3_dble(addr) \
_prefetch_192b((addr))

#define _prefetch_pauli_dble(addr) \
_prefetch_320b((addr))

#elif (defined P3)

#define _pfbase(addr) ((unsigned long)(addr)&(~0x1fL))

#define _prefetch_64b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1" \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x20L))))

#define _prefetch_96b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2" \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x20L))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))))

#define _prefetch_160b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2 \n\t" \
                      "prefetcht0 %3 \n\t" \
                      "prefetcht0 %4" \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x20L))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))), \
                      "m" (*((char*)(_pfbase(addr)+0x60L))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L))))

#define _prefetch_192b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2 \n\t" \
                      "prefetcht0 %3 \n\t" \
                      "prefetcht0 %4 \n\t" \
                      "prefetcht0 %5" \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x20L))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))), \
                      "m" (*((char*)(_pfbase(addr)+0x60L))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L))), \
                      "m" (*((char*)(_pfbase(addr)+0xa0L))))

#define _prefetch_288b(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2 \n\t" \
                      "prefetcht0 %3 \n\t" \
                      "prefetcht0 %4" \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)))), \
                      "m" (*((char*)(_pfbase(addr)+0x20L))), \
                      "m" (*((char*)(_pfbase(addr)+0x40L))), \
                      "m" (*((char*)(_pfbase(addr)+0x60L))), \
                      "m" (*((char*)(_pfbase(addr)+0x80L)))); \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2 \n\t" \
                      "prefetcht0 %3" \
                      : \
                      : \
                      "m" (*((char*)(_pfbase(addr)+0xa0L))), \
                      "m" (*((char*)(_pfbase(addr)+0xc0L))), \
                      "m" (*((char*)(_pfbase(addr)+0xe0L))), \
                      "m" (*((char*)(_pfbase(addr)+0x100L))))

#define _prefetch_su3_alg_dble(addr) \
_prefetch_64b((addr))

#define _prefetch_weyl(addr) \
_prefetch_64b((addr))

#define _prefetch_spinor(addr) \
_prefetch_96b((addr))

#define _prefetch_su3(addr) \
_prefetch_96b((addr))

#define _prefetch_pauli(addr) \
_prefetch_160b((addr))

#define _prefetch_weyl_dble(addr) \
_prefetch_96b((addr))

#define _prefetch_spinor_dble(addr) \
_prefetch_192b((addr))

#define _prefetch_su3_dble(addr) \
_prefetch_160b((addr))

#define _prefetch_pauli_dble(addr) \
_prefetch_288b((addr))

#else

#define _prefetch_su3_alg_dble(addr)

#define _prefetch_weyl(addr)

#define _prefetch_spinor(addr)

#define _prefetch_su3(addr)

#define _prefetch_pauli(addr)

#define _prefetch_weyl_dble(addr)

#define _prefetch_spinor_dble(addr)

#define _prefetch_su3_dble(addr)

#define _prefetch_pauli_dble(addr)

#endif

/*******************************************************************************
*
* Macros for su3_vector data
*
* Most of these macros operate on pairs of su3 vectors that are stored
* in the low and high words of xmm0,xmm1,xmm2 or xmm3,xmm4,xmm5. For example,
*
* xmm0 -> sl.c1.re,sl.c1.im,sh.c1.re,sh.c1.im
* xmm1 -> sl.c2.re,sl.c2.im,sh.c2.re,sh.c2.im
* xmm2 -> sl.c3.re,sl.c3.im,sh.c3.re,sh.c3.im
*
* (where sl and sh are of type su3_vector). This can also be interpreted as
* an sse_vector s that is stored in these registers according to
*
* xmm0 -> s.c1.c1,s.c1.c2,s.c1.c3,s.c1.c4
* xmm1 -> s.c2.c1,s.c2.c2,s.c2.c3,s.c2.c4
* xmm2 -> s.c3.c1,s.c3.c2,s.c3.c3,s.c3.c4
*
* The load and store macros can be used to move data in either format
* from and to the xmm registers
*
*******************************************************************************/

/*
* Loads two su3 vectors sl and sh to the low and high words of xmm0,xmm1,xmm2
*/

#define _sse_pair_load(sl,sh) \
__asm__ __volatile__ ("movsd %0, %%xmm0 \n\t" \
                      "movsd %1, %%xmm1 \n\t" \
                      "movsd %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((sl).c1), \
                      "m" ((sl).c2), \
                      "m" ((sl).c3), \
                      "m" ((sh).c1), \
                      "m" ((sh).c2), \
                      "m" ((sh).c3) \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Loads two su3 vectors sl and sh to the low and high words of xmm3,xmm4,xmm5
*/

#define _sse_pair_load_up(sl,sh) \
__asm__ __volatile__ ("movsd %0, %%xmm3 \n\t" \
                      "movsd %1, %%xmm4 \n\t" \
                      "movsd %2, %%xmm5 \n\t" \
                      "movhps %3, %%xmm3 \n\t" \
                      "movhps %4, %%xmm4 \n\t" \
                      "movhps %5, %%xmm5" \
                      : \
                      : \
                      "m" ((sl).c1), \
                      "m" ((sl).c2), \
                      "m" ((sl).c3), \
                      "m" ((sh).c1), \
                      "m" ((sh).c2), \
                      "m" ((sh).c3) \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Stores the low and high words of xmm0,xmm1,xmm2 to the su3 vectors rl and rh
*/

#define _sse_pair_store(rl,rh) \
__asm__ __volatile__ ("movlps %%xmm0, %0 \n\t" \
                      "movlps %%xmm1, %1 \n\t" \
                      "movlps %%xmm2, %2 \n\t" \
                      "movhps %%xmm0, %3 \n\t" \
                      "movhps %%xmm1, %4 \n\t" \
                      "movhps %%xmm2, %5" \
                      : \
                      "=m" ((rl).c1), \
                      "=m" ((rl).c2), \
                      "=m" ((rl).c3), \
                      "=m" ((rh).c1), \
                      "=m" ((rh).c2), \
                      "=m" ((rh).c3))

/*
* Stores the low and high words of xmm3,xmm4,xmm5 to the su3 vectors rl and rh
*/

#define _sse_pair_store_up(rl,rh) \
__asm__ __volatile__ ("movlps %%xmm3, %0 \n\t" \
                      "movlps %%xmm4, %1 \n\t" \
                      "movlps %%xmm5, %2 \n\t" \
                      "movhps %%xmm3, %3 \n\t" \
                      "movhps %%xmm4, %4 \n\t" \
                      "movhps %%xmm5, %5" \
                      : \
                      "=m" ((rl).c1), \
                      "=m" ((rl).c2), \
                      "=m" ((rl).c3), \
                      "=m" ((rh).c1), \
                      "=m" ((rh).c2), \
                      "=m" ((rh).c3))

/*
* Loads the components of a Weyl spinor s to xmm0,xmm1,xmm2
*/

#define _sse_weyl_load(s) \
__asm__ __volatile__ ("movaps %0, %%xmm0 \n\t" \
                      "movaps %2, %%xmm1 \n\t" \
                      "movaps %4, %%xmm2" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Loads the components of a Weyl spinor s to xmm3,xmm4,xmm5
*/

#define _sse_weyl_load_up(s) \
__asm__ __volatile__ ("movaps %0, %%xmm3 \n\t" \
                      "movaps %2, %%xmm4 \n\t" \
                      "movaps %4, %%xmm5" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
* Stores xmm0,xmm1,xmm2 to the components of a Weyl spinor s
*/

#define _sse_weyl_store(s) \
__asm__ __volatile__ ("movaps %%xmm0, %0 \n\t" \
                      "movaps %%xmm1, %2 \n\t" \
                      "movaps %%xmm2, %4" \
                      : \
                      "=m" ((s).c1.c1), \
                      "=m" ((s).c1.c2), \
                      "=m" ((s).c1.c3), \
                      "=m" ((s).c2.c1), \
                      "=m" ((s).c2.c2), \
                      "=m" ((s).c2.c3))

/*
* Stores xmm3,xmm4,xmm5 to the components of a Weyl spinor s
*/

#define _sse_weyl_store_up(s) \
__asm__ __volatile__ ("movaps %%xmm3, %0 \n\t" \
                      "movaps %%xmm4, %2 \n\t" \
                      "movaps %%xmm5, %4" \
                      : \
                      "=m" ((s).c1.c1), \
                      "=m" ((s).c1.c2), \
                      "=m" ((s).c1.c3), \
                      "=m" ((s).c2.c1), \
                      "=m" ((s).c2.c2), \
                      "=m" ((s).c2.c3))

/*
* Adds xmm3,xmm4,xmm5 to xmm0,xmm1,xmm2
*/

#define _sse_vector_add() \
__asm__ __volatile__ ("addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Subtracts xmm3,xmm4,xmm5 from xmm0,xmm1,xmm2
*/

#define _sse_vector_sub() \
__asm__ __volatile__ ("subps %%xmm3, %%xmm0 \n\t" \
                      "subps %%xmm4, %%xmm1 \n\t" \
                      "subps %%xmm5, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2")

/*
* Multiplies the high words xmm3,xmm4,xmm5 with -1 and adds these registers
* to xmm0,xmm1,xmm2
*/

#define _sse_vector_addsub() \
__asm__ __volatile__ ("mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn34) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies the low words xmm3,xmm4,xmm5 with -1 and adds these registers
* to xmm0,xmm1,xmm2
*/

#define _sse_vector_subadd() \
__asm__ __volatile__ ("mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn12) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies xmm3,xmm4,xmm5 with i and adds them to xmm0,xmm1,xmm2
*/

#define _sse_vector_i_add() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "addsubps %%xmm3, %%xmm0 \n\t" \
                      "addsubps %%xmm4, %%xmm1 \n\t" \
                      "addsubps %%xmm5, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies xmm3,xmm4,xmm5 with i and subtracts them from xmm0,xmm1,xmm2
*/

#define _sse_vector_i_sub() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn24) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Exchanges the high and low words of xmm3,xmm4,xmm5, multiplies them with i
* and adds the result to xmm0,xmm1,xmm2
*/

#define _sse_vector_xch_i_add() \
__asm__ __volatile__ ("shufps $0x1b, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x1b, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x1b, %%xmm5, %%xmm5 \n\t" \
                      "addsubps %%xmm3, %%xmm0 \n\t" \
                      "addsubps %%xmm4, %%xmm1 \n\t" \
                      "addsubps %%xmm5, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Exchanges the high and low words of xmm3,xmm4,xmm5, multiplies them with i
* and subtracts the result from xmm0,xmm1,xmm2
*/

#define _sse_vector_xch_i_sub() \
__asm__ __volatile__ ("shufps $0x1b, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x1b, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x1b, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn24) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies the low and high words of xmm3,xmm4,xmm5 with i and -i
* respectively and adds these registers to xmm0,xmm1,xmm2
*/

#define _sse_vector_i_addsub() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn14) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Multiplies the low and high words of xmm3,xmm4,xmm5 with -i and i
* respectively and adds these registers to xmm0,xmm1,xmm2
*/

#define _sse_vector_i_subadd() \
__asm__ __volatile__ ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
                      "mulps %0, %%xmm3 \n\t" \
                      "mulps %0, %%xmm4 \n\t" \
                      "mulps %0, %%xmm5 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      "addps %%xmm4, %%xmm1 \n\t" \
                      "addps %%xmm5, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn23) \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm3", "xmm4", "xmm5")

/*
* Exchanges the high and low words in xmm3,xmm4,xmm5
*/

#define _sse_vector_xch() \
__asm__ __volatile__ ("shufps $0x4e, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x4e, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x4e, %%xmm5, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

/******************************************************************************
*
*  Action of su3 matrices on su3 vectors
*
******************************************************************************/

/*
* Multiplies a pair sl,sh of su3 vectors with an su3 matrix u,
* assuming sl and sh are in the low and high words of xmm0,xmm1,xmm2
*
* On output the result is in xmm3,xmm4,xmm5 and the registers
* xmm0,xmm1,xmm2 are changed
*/

#define _sse_su3_multiply(u) \
__asm__ __volatile__ ("movss %0, %%xmm3 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm4 \n\t" \
                      "movss %3, %%xmm7 \n\t" \
                      "movss %4, %%xmm5 \n\t" \
                      "movss %5, %%xmm8 \n\t" \
                      "shufps $0x0, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm5, %%xmm5 \n\t" \
                      "shufps $0x0, %%xmm8, %%xmm8" \
                      : \
                      : \
                      "m" ((u).c11.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c21.re), \
                      "m" ((u).c22.re), \
                      "m" ((u).c31.re), \
                      "m" ((u).c32.re) \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm4 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm5 \n\t" \
                      "mulps %%xmm1, %%xmm8 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "addps %%xmm8, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("movss %0, %%xmm9 \n\t" \
                      "movss %1, %%xmm10 \n\t" \
                      "movss %2, %%xmm11 \n\t" \
                      "movss %3, %%xmm6 \n\t" \
                      "movss %4, %%xmm7 \n\t" \
                      "movss %5, %%xmm8 \n\t" \
                      "shufps $0xb1, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0x0, %%xmm9, %%xmm9 \n\t" \
                      "shufps $0x0, %%xmm10, %%xmm10 \n\t" \
                      "shufps $0x0, %%xmm11, %%xmm11 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm8, %%xmm8" \
                      : \
                      : \
                      "m" ((u).c13.re), \
                      "m" ((u).c21.im), \
                      "m" ((u).c33.re), \
                      "m" ((u).c11.im), \
                      "m" ((u).c23.re), \
                      "m" ((u).c31.im) \
                      : \
                      "xmm0", "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("mulps %%xmm2, %%xmm9 \n\t" \
                      "mulps %%xmm0, %%xmm10 \n\t" \
                      "mulps %%xmm2, %%xmm11 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm8 \n\t" \
                      "addps %%xmm9, %%xmm3 \n\t" \
                      "addsubps %%xmm10, %%xmm4 \n\t" \
                      "addps %%xmm11, %%xmm5 \n\t" \
                      "addsubps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "addsubps %%xmm8, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8",  \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("movss %0, %%xmm9 \n\t" \
                      "movss %1, %%xmm10 \n\t" \
                      "movss %2, %%xmm11 \n\t" \
                      "movss %3, %%xmm6 \n\t" \
                      "movss %4, %%xmm7 \n\t" \
                      "movss %5, %%xmm8 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "shufps $0x0, %%xmm9, %%xmm9 \n\t" \
                      "shufps $0x0, %%xmm10, %%xmm10 \n\t" \
                      "shufps $0x0, %%xmm11, %%xmm11 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm8, %%xmm8" \
                      : \
                      : \
                      "m" ((u).c12.im), \
                      "m" ((u).c23.im), \
                      "m" ((u).c32.im), \
                      "m" ((u).c13.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c33.im) \
                      : \
                      "xmm1", "xmm2", "xmm6", "xmm7", \
                      "xmm8", "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("mulps %%xmm1, %%xmm9 \n\t" \
                      "mulps %%xmm2, %%xmm10 \n\t" \
                      "mulps %%xmm1, %%xmm11 \n\t" \
                      "mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm8 \n\t" \
                      "addsubps %%xmm9, %%xmm3 \n\t" \
                      "addsubps %%xmm10, %%xmm4 \n\t" \
                      "addsubps %%xmm11, %%xmm5 \n\t" \
                      "addsubps %%xmm6, %%xmm3 \n\t" \
                      "addsubps %%xmm7, %%xmm4 \n\t" \
                      "addsubps %%xmm8, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8",  \
                      "xmm9", "xmm10", "xmm11")

/*
* Multiplies a pair sl,sh of su3 vectors with an su3 matrix u^dagger,
* assuming sl and sh are in the low and high words of xmm0,xmm1,xmm2
*
* On output the result is in xmm3,xmm4,xmm5 and the registers
* xmm0,xmm1,xmm2 are changed
*/

#define _sse_su3_inverse_multiply(u) \
__asm__ __volatile__ ("movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm9 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "movss %3, %%xmm10 \n\t" \
                      "movss %4, %%xmm8 \n\t" \
                      "movss %5, %%xmm11 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm9, %%xmm9 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm10, %%xmm10 \n\t" \
                      "shufps $0x0, %%xmm8, %%xmm8 \n\t" \
                      "shufps $0x0, %%xmm11, %%xmm11" \
                      : \
                      : \
                      "m" ((u).c11.im), \
                      "m" ((u).c21.im), \
                      "m" ((u).c12.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c13.im), \
                      "m" ((u).c23.im) \
                      : \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm9 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm10 \n\t" \
                      "mulps %%xmm0, %%xmm8 \n\t" \
                      "mulps %%xmm1, %%xmm11 \n\t" \
                      "addps %%xmm6, %%xmm9 \n\t" \
                      "addps %%xmm7, %%xmm10 \n\t" \
                      "addps %%xmm8, %%xmm11" \
                      : \
                      : \
                      : \
                      "xmm6", "xmm7", "xmm8", \
                      "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("movss %0, %%xmm3 \n\t" \
                      "movss %1, %%xmm4 \n\t" \
                      "movss %2, %%xmm5 \n\t" \
                      "movss %3, %%xmm6 \n\t" \
                      "movss %4, %%xmm7 \n\t" \
                      "movss %5, %%xmm8 \n\t" \
                      "shufps $0xb1, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0x0, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x0, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x0, %%xmm5, %%xmm5 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm8, %%xmm8" \
                      : \
                      : \
                      "m" ((u).c11.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c13.re), \
                      "m" ((u).c31.im), \
                      "m" ((u).c32.im), \
                      "m" ((u).c33.im) \
                      : \
                      "xmm0", "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm0, %%xmm4 \n\t" \
                      "mulps %%xmm0, %%xmm5 \n\t" \
                      "mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm8 \n\t" \
                      "addsubps %%xmm9, %%xmm3 \n\t" \
                      "addsubps %%xmm10, %%xmm4 \n\t" \
                      "addsubps %%xmm11, %%xmm5 \n\t" \
                      "addsubps %%xmm6, %%xmm3 \n\t" \
                      "addsubps %%xmm7, %%xmm4 \n\t" \
                      "addsubps %%xmm8, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("movss %0, %%xmm9 \n\t" \
                      "movss %1, %%xmm10 \n\t" \
                      "movss %2, %%xmm11 \n\t" \
                      "movss %3, %%xmm6 \n\t" \
                      "movss %4, %%xmm7 \n\t" \
                      "movss %5, %%xmm8 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "shufps $0x0, %%xmm9, %%xmm9 \n\t" \
                      "shufps $0x0, %%xmm10, %%xmm10 \n\t" \
                      "shufps $0x0, %%xmm11, %%xmm11 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm8, %%xmm8" \
                      : \
                      : \
                      "m" ((u).c21.re), \
                      "m" ((u).c32.re), \
                      "m" ((u).c23.re), \
                      "m" ((u).c31.re), \
                      "m" ((u).c22.re), \
                      "m" ((u).c33.re) \
                      : \
                      "xmm1", "xmm2", "xmm6", "xmm7", \
                      "xmm8", "xmm9", "xmm10", "xmm11"); \
__asm__ __volatile__ ("mulps %%xmm1, %%xmm9 \n\t" \
                      "mulps %%xmm2, %%xmm10 \n\t" \
                      "mulps %%xmm1, %%xmm11 \n\t" \
                      "mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm8 \n\t" \
                      "addps %%xmm9, %%xmm3 \n\t" \
                      "addps %%xmm10, %%xmm4 \n\t" \
                      "addps %%xmm11, %%xmm5 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "addps %%xmm8, %%xmm5 \n\t" \
                      "shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0xb1, %%xmm5, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm6", "xmm7", "xmm8",  \
                      "xmm9", "xmm10", "xmm11")

/******************************************************************************
*
*  Macros for Dirac spinors
*
******************************************************************************/

/*
*  Loads the spinor s to the registers xmm0,..,xmm5 in linear order
*/

#define _sse_spinor_load(s) \
__asm__ __volatile__ ("movaps %0, %%xmm0 \n\t" \
                      "movaps %2, %%xmm1 \n\t" \
                      "movaps %4, %%xmm2" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm0", "xmm1", "xmm2"); \
__asm__ __volatile__ ("movaps %0, %%xmm3 \n\t" \
                      "movaps %2, %%xmm4 \n\t" \
                      "movaps %4, %%xmm5" \
                      : \
                      : \
                      "m" ((s).c3.c1), \
                      "m" ((s).c3.c2), \
                      "m" ((s).c3.c3), \
                      "m" ((s).c4.c1), \
                      "m" ((s).c4.c2), \
                      "m" ((s).c4.c3) \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
*  Loads the spinor s to the registers xmm6,..,xmm11 in linear order
*/

#define _sse_spinor_load_up(s) \
__asm__ __volatile__ ("movaps %0, %%xmm6 \n\t" \
                      "movaps %2, %%xmm7 \n\t" \
                      "movaps %4, %%xmm8" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm6", "xmm7", "xmm8"); \
__asm__ __volatile__ ("movaps %0, %%xmm9 \n\t" \
                      "movaps %2, %%xmm10 \n\t" \
                      "movaps %4, %%xmm11" \
                      : \
                      : \
                      "m" ((s).c3.c1), \
                      "m" ((s).c3.c2), \
                      "m" ((s).c3.c3), \
                      "m" ((s).c4.c1), \
                      "m" ((s).c4.c2), \
                      "m" ((s).c4.c3) \
                      : \
                      "xmm9", "xmm10", "xmm11")

/*
*  Stores the registers xmm0,..,xmm5 to the spinor s in linear order
*/

#define _sse_spinor_store(s) \
__asm__ __volatile__ ("movaps %%xmm0, %0 \n\t" \
                      "movaps %%xmm1, %2 \n\t" \
                      "movaps %%xmm2, %4" \
                      : \
                      "=m" ((s).c1.c1), \
                      "=m" ((s).c1.c2), \
                      "=m" ((s).c1.c3), \
                      "=m" ((s).c2.c1), \
                      "=m" ((s).c2.c2), \
                      "=m" ((s).c2.c3)); \
__asm__ __volatile__ ("movaps %%xmm3, %0 \n\t" \
                      "movaps %%xmm4, %2 \n\t" \
                      "movaps %%xmm5, %4" \
                      : \
                      "=m" ((s).c3.c1), \
                      "=m" ((s).c3.c2), \
                      "=m" ((s).c3.c3), \
                      "=m" ((s).c4.c1), \
                      "=m" ((s).c4.c2), \
                      "=m" ((s).c4.c3))

/*
*  Stores the registers xmm6,..,xmm11 to the spinor s in linear order
*/

#define _sse_spinor_store_up(s) \
__asm__ __volatile__ ("movaps %%xmm6, %0 \n\t" \
                      "movaps %%xmm7, %2 \n\t" \
                      "movaps %%xmm8, %4" \
                      : \
                      "=m" ((s).c1.c1), \
                      "=m" ((s).c1.c2), \
                      "=m" ((s).c1.c3), \
                      "=m" ((s).c2.c1), \
                      "=m" ((s).c2.c2), \
                      "=m" ((s).c2.c3)); \
__asm__ __volatile__ ("movaps %%xmm9, %0 \n\t" \
                      "movaps %%xmm10, %2 \n\t" \
                      "movaps %%xmm11, %4" \
                      : \
                      "=m" ((s).c3.c1), \
                      "=m" ((s).c3.c2), \
                      "=m" ((s).c3.c3), \
                      "=m" ((s).c4.c1), \
                      "=m" ((s).c4.c2), \
                      "=m" ((s).c4.c3))

/*
*  Loads (z.re,z.re,z.re,z.re) to xmm6 and (-z.im,z.im,-z.im,z.im) to xmm7
*/

#define _sse_load_cmplx(z) \
__asm__ __volatile__ ("movss %0, %%xmm6 \n\t" \
                      "movss %1, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
                      "mulps %2, %%xmm7" \
                      : \
                      : \
                      "m" ((z).re), \
                      "m" ((z).im), \
                      "m" (_sse_sgn13) \
                      : \
                      "xmm6", "xmm7")

/*
*  Multiplies the spinor s by the complex number z and assigns the result to
*  xmm0,..,xmm5, assuming z was loaded to xmm6,xmm7 using _sse_load_cmplx(z)
*/

#define _sse_mulc_spinor(s) \
__asm__ __volatile__ ("movaps %0, %%xmm0 \n\t" \
                      "movaps %2, %%xmm1 \n\t" \
                      "movaps %4, %%xmm2" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm0", "xmm1", "xmm2");  \
__asm__ __volatile__ ("movaps %%xmm0, %%xmm8 \n\t" \
                      "movaps %%xmm1, %%xmm9 \n\t" \
                      "movaps %%xmm2, %%xmm10 \n\t" \
                      "mulps %%xmm6, %%xmm0 \n\t" \
                      "mulps %%xmm6, %%xmm1 \n\t" \
                      "mulps %%xmm6, %%xmm2 \n\t" \
                      "shufps $0xb1, %%xmm8, %%xmm8 \n\t" \
                      "shufps $0xb1, %%xmm9, %%xmm9 \n\t" \
                      "shufps $0xb1, %%xmm10, %%xmm10 \n\t" \
                      "mulps %%xmm7, %%xmm8 \n\t" \
                      "mulps %%xmm7, %%xmm9 \n\t" \
                      "mulps %%xmm7, %%xmm10 \n\t" \
                      "addps %%xmm8, %%xmm0 \n\t" \
                      "addps %%xmm9, %%xmm1 \n\t" \
                      "addps %%xmm10, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm8", "xmm9", "xmm10"); \
__asm__ __volatile__ ("movaps %0, %%xmm3 \n\t" \
                      "movaps %2, %%xmm4 \n\t" \
                      "movaps %4, %%xmm5" \
                      : \
                      : \
                      "m" ((s).c3.c1), \
                      "m" ((s).c3.c2), \
                      "m" ((s).c3.c3), \
                      "m" ((s).c4.c1), \
                      "m" ((s).c4.c2), \
                      "m" ((s).c4.c3) \
                      : \
                      "xmm3", "xmm4", "xmm5");  \
__asm__ __volatile__ ("movaps %%xmm3, %%xmm11 \n\t" \
                      "movaps %%xmm4, %%xmm12 \n\t" \
                      "movaps %%xmm5, %%xmm13 \n\t" \
                      "mulps %%xmm6, %%xmm3 \n\t" \
                      "mulps %%xmm6, %%xmm4 \n\t" \
                      "mulps %%xmm6, %%xmm5 \n\t" \
                      "shufps $0xb1, %%xmm11, %%xmm11 \n\t" \
                      "shufps $0xb1, %%xmm12, %%xmm12 \n\t" \
                      "shufps $0xb1, %%xmm13, %%xmm13 \n\t" \
                      "mulps %%xmm7, %%xmm11 \n\t" \
                      "mulps %%xmm7, %%xmm12 \n\t" \
                      "mulps %%xmm7, %%xmm13 \n\t" \
                      "addps %%xmm11, %%xmm3 \n\t" \
                      "addps %%xmm12, %%xmm4 \n\t" \
                      "addps %%xmm13, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm11", "xmm12", "xmm13")


/*
*  Multiplies the spinor s by the complex number z and adds the result to
*  xmm0,..,xmm5, assuming z was loaded to xmm6,xmm7 using _sse_load_cmplx(z)
*/

#define _sse_mulc_spinor_add(s) \
__asm__ __volatile__ ("movaps %0, %%xmm8 \n\t" \
                      "movaps %2, %%xmm9 \n\t" \
                      "movaps %4, %%xmm10" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm8", "xmm9", "xmm10");  \
__asm__ __volatile__ ("movaps %%xmm8, %%xmm11 \n\t" \
                      "movaps %%xmm9, %%xmm12 \n\t" \
                      "movaps %%xmm10, %%xmm13 \n\t" \
                      "mulps %%xmm6, %%xmm8 \n\t" \
                      "mulps %%xmm6, %%xmm9 \n\t" \
                      "mulps %%xmm6, %%xmm10 \n\t" \
                      "shufps $0xb1, %%xmm11, %%xmm11 \n\t" \
                      "shufps $0xb1, %%xmm12, %%xmm12 \n\t" \
                      "shufps $0xb1, %%xmm13, %%xmm13 \n\t" \
                      "addps %%xmm8, %%xmm0 \n\t" \
                      "addps %%xmm9, %%xmm1 \n\t" \
                      "addps %%xmm10, %%xmm2 \n\t" \
                      "mulps %%xmm7, %%xmm11 \n\t" \
                      "mulps %%xmm7, %%xmm12 \n\t" \
                      "mulps %%xmm7, %%xmm13 \n\t" \
                      "addps %%xmm11, %%xmm0 \n\t" \
                      "addps %%xmm12, %%xmm1 \n\t" \
                      "addps %%xmm13, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm8", "xmm9", "xmm10", \
                      "xmm11", "xmm12", "xmm13"); \
__asm__ __volatile__ ("movaps %0, %%xmm8 \n\t" \
                      "movaps %2, %%xmm9 \n\t" \
                      "movaps %4, %%xmm10" \
                      : \
                      : \
                      "m" ((s).c3.c1), \
                      "m" ((s).c3.c2), \
                      "m" ((s).c3.c3), \
                      "m" ((s).c4.c1), \
                      "m" ((s).c4.c2), \
                      "m" ((s).c4.c3) \
                      : \
                      "xmm8", "xmm9", "xmm10");  \
__asm__ __volatile__ ("movaps %%xmm8, %%xmm11 \n\t" \
                      "movaps %%xmm9, %%xmm12 \n\t" \
                      "movaps %%xmm10, %%xmm13 \n\t" \
                      "mulps %%xmm6, %%xmm8 \n\t" \
                      "mulps %%xmm6, %%xmm9 \n\t" \
                      "mulps %%xmm6, %%xmm10 \n\t" \
                      "shufps $0xb1, %%xmm11, %%xmm11 \n\t" \
                      "shufps $0xb1, %%xmm12, %%xmm12 \n\t" \
                      "shufps $0xb1, %%xmm13, %%xmm13 \n\t" \
                      "addps %%xmm8, %%xmm3 \n\t" \
                      "addps %%xmm9, %%xmm4 \n\t" \
                      "addps %%xmm10, %%xmm5 \n\t" \
                      "mulps %%xmm7, %%xmm11 \n\t" \
                      "mulps %%xmm7, %%xmm12 \n\t" \
                      "mulps %%xmm7, %%xmm13 \n\t" \
                      "addps %%xmm11, %%xmm3 \n\t" \
                      "addps %%xmm12, %%xmm4 \n\t" \
                      "addps %%xmm13, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm8", "xmm9", "xmm10", \
                      "xmm11", "xmm12", "xmm13")

/*
*  Loads (c,c,c,c) to xmm6 and xmm7
*/

#define _sse_load_real(c) \
__asm__ __volatile__ ("movss %0, %%xmm6 \n\t" \
                      "movss %0, %%xmm7 \n\t" \
                      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x0, %%xmm7, %%xmm7" \
                      : \
                      : \
                      "m" (c) \
                      : \
                      "xmm6", "xmm7")

/*
*  Multiplies the spinor s by the real number c and assigns the result to
*  xmm0,..,xmm5, assuming c was loaded to xmm6,xmm7 using _sse_load_real(c)
*/

#define _sse_mulr_spinor(s) \
__asm__ __volatile__ ("movaps %0, %%xmm0 \n\t" \
                      "movaps %2, %%xmm1 \n\t" \
                      "movaps %4, %%xmm2" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm0", "xmm1", "xmm2");  \
__asm__ __volatile__ ("mulps %%xmm6, %%xmm0 \n\t" \
                      "mulps %%xmm7, %%xmm1 \n\t" \
                      "mulps %%xmm6, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2"); \
__asm__ __volatile__ ("movaps %0, %%xmm3 \n\t" \
                      "movaps %2, %%xmm4 \n\t" \
                      "movaps %4, %%xmm5" \
                      : \
                      : \
                      "m" ((s).c3.c1), \
                      "m" ((s).c3.c2), \
                      "m" ((s).c3.c3), \
                      "m" ((s).c4.c1), \
                      "m" ((s).c4.c2), \
                      "m" ((s).c4.c3) \
                      : \
                      "xmm3", "xmm4", "xmm5");  \
__asm__ __volatile__ ("mulps %%xmm7, %%xmm3 \n\t" \
                      "mulps %%xmm6, %%xmm4 \n\t" \
                      "mulps %%xmm7, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5")

/*
*  Multiplies the spinor s by the real number c and adds the result to
*  xmm0,..,xmm5, assuming c was loaded to xmm6,xmm7 using _sse_load_real(c)
*/

#define _sse_mulr_spinor_add(s) \
__asm__ __volatile__ ("movaps %0, %%xmm8 \n\t" \
                      "movaps %2, %%xmm9 \n\t" \
                      "movaps %4, %%xmm10" \
                      : \
                      : \
                      "m" ((s).c1.c1), \
                      "m" ((s).c1.c2), \
                      "m" ((s).c1.c3), \
                      "m" ((s).c2.c1), \
                      "m" ((s).c2.c2), \
                      "m" ((s).c2.c3) \
                      : \
                      "xmm8", "xmm9", "xmm10");  \
__asm__ __volatile__ ("mulps %%xmm6, %%xmm8 \n\t" \
                      "mulps %%xmm7, %%xmm9 \n\t" \
                      "mulps %%xmm6, %%xmm10 \n\t" \
                      "addps %%xmm8, %%xmm0 \n\t" \
                      "addps %%xmm9, %%xmm1 \n\t" \
                      "addps %%xmm10, %%xmm2" \
                      : \
                      : \
                      : \
                      "xmm0", "xmm1", "xmm2", \
                      "xmm8", "xmm9", "xmm10"); \
__asm__ __volatile__ ("movaps %0, %%xmm11 \n\t" \
                      "movaps %2, %%xmm12 \n\t" \
                      "movaps %4, %%xmm13" \
                      : \
                      : \
                      "m" ((s).c3.c1), \
                      "m" ((s).c3.c2), \
                      "m" ((s).c3.c3), \
                      "m" ((s).c4.c1), \
                      "m" ((s).c4.c2), \
                      "m" ((s).c4.c3) \
                      : \
                      "xmm11", "xmm12", "xmm13");  \
__asm__ __volatile__ ("mulps %%xmm7, %%xmm11 \n\t" \
                      "mulps %%xmm6, %%xmm12 \n\t" \
                      "mulps %%xmm7, %%xmm13 \n\t" \
                      "addps %%xmm11, %%xmm3 \n\t" \
                      "addps %%xmm12, %%xmm4 \n\t" \
                      "addps %%xmm13, %%xmm5" \
                      : \
                      : \
                      : \
                      "xmm3", "xmm4", "xmm5", \
                      "xmm11", "xmm12", "xmm13")

#endif
