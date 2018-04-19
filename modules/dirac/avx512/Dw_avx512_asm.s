# mark_description "Intel(R) C Intel(R) 64 Compiler for applications running on Intel(R) 64, Version 17.0.4.196 Build 20170411";
# mark_description "-I../../../include -I.. -I/cineca/prod/opt/compilers/intel/pe-xe-2017/binary/impi/2017.3.196/intel64/include";
# mark_description " -isystem /cineca/prod/opt/compilers/intel/pe-xe-2018/binary/impi/2018.1.163/include64/ -std=c89 -xCORE-AVX5";
# mark_description "12 -mtune=skylake -DAVX512 -O3 -Ddirac_counters -pedantic -fstrict-aliasing -Wno-long-long -Wstrict-prototyp";
# mark_description "es -S";
	.file "Dw_avx512.c"
	.text
..TXTST0:
# -- Begin  doe_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl doe_avx512
# --- doe_avx512(int *, int *, su3 *, spinor *, float, spin_t *)
doe_avx512:
# parameter 1: %rdi
# parameter 2: %rsi
# parameter 3: %rdx
# parameter 4: %rcx
# parameter 5: %xmm0
# parameter 6: %r8
..B1.1:                         # Preds ..B1.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_doe_avx512.1:
..L2:
                                                          #26.1
        pushq     %rbx                                          #26.1
	.cfi_def_cfa_offset 16
        movq      %rsp, %rbx                                    #26.1
	.cfi_def_cfa 3, 16
	.cfi_offset 3, -16
        andq      $-64, %rsp                                    #26.1
        pushq     %rbp                                          #26.1
        pushq     %rbp                                          #26.1
        movq      8(%rbx), %rbp                                 #26.1
        movq      %rbp, 8(%rsp)                                 #26.1
        movq      %rsp, %rbp                                    #26.1
	.cfi_escape 0x10, 0x06, 0x02, 0x76, 0x00
        movslq    (%rdi), %rax                                  #39.16
        movslq    (%rsi), %r11                                  #40.16
        movslq    4(%rdi), %r10                                 #41.17
        vmovss    %xmm0, -16(%rbp)                              #26.1
        lea       (%rax,%rax,2), %r9                            #39.8
        shlq      $5, %r9                                       #39.8
        lea       (%r11,%r11,2), %r11                           #40.8
        shlq      $5, %r11                                      #40.8
        lea       (%r10,%r10,2), %r10                           #41.9
        movslq    4(%rsi), %rax                                 #42.17
        shlq      $5, %r10                                      #41.9
        vmovups   (%rcx,%r9), %xmm29                            #44.3
        vmovups   16(%rcx,%r9), %xmm10                          #44.3
        vmovups   32(%rcx,%r9), %xmm8                           #44.3
        vmovups   48(%rcx,%r9), %xmm14                          #47.3
        vmovups   64(%rcx,%r9), %xmm12                          #47.3
        vmovups   80(%rcx,%r9), %xmm26                          #47.3
        vmovups   .L_2il0floatpacket.10(%rip), %zmm11           #47.3
        vmovups   .L_2il0floatpacket.11(%rip), %zmm6            #47.3
        movslq    8(%rdi), %r9                                  #51.16
        lea       (%rax,%rax,2), %rax                           #42.9
        shlq      $5, %rax                                      #42.9
        vmovaps   %zmm11, %zmm25                                #47.3
        lea       (%r9,%r9,2), %r9                              #51.8
        shlq      $5, %r9                                       #51.8
        vinsertf32x4 $1, (%rcx,%r11), %zmm29, %zmm21            #44.3
        vinsertf32x4 $2, (%rcx,%r10), %zmm21, %zmm22            #44.3
        vinsertf32x4 $3, (%rcx,%rax), %zmm22, %zmm19            #44.3
        vinsertf32x4 $1, 16(%rcx,%r11), %zmm10, %zmm16          #44.3
        vinsertf32x4 $1, 32(%rcx,%r11), %zmm8, %zmm9            #44.3
        vinsertf32x4 $1, 48(%rcx,%r11), %zmm14, %zmm0           #47.3
        vinsertf32x4 $1, 64(%rcx,%r11), %zmm12, %zmm7           #47.3
        vinsertf32x4 $1, 80(%rcx,%r11), %zmm26, %zmm29          #47.3
        vinsertf32x4 $2, 16(%rcx,%r10), %zmm16, %zmm17          #44.3
        vinsertf32x4 $2, 32(%rcx,%r10), %zmm9, %zmm18           #44.3
        vinsertf32x4 $2, 48(%rcx,%r10), %zmm0, %zmm15           #47.3
        vinsertf32x4 $2, 64(%rcx,%r10), %zmm7, %zmm24           #47.3
        vinsertf32x4 $2, 80(%rcx,%r10), %zmm29, %zmm5           #47.3
        vinsertf32x4 $3, 16(%rcx,%rax), %zmm17, %zmm20          #44.3
        vinsertf32x4 $3, 32(%rcx,%rax), %zmm18, %zmm13          #44.3
        vinsertf32x4 $3, 48(%rcx,%rax), %zmm15, %zmm22          #47.3
        vinsertf32x4 $3, 64(%rcx,%rax), %zmm24, %zmm16          #47.3
        vinsertf32x4 $3, 80(%rcx,%rax), %zmm5, %zmm28           #47.3
        vshufps   $228, %zmm20, %zmm19, %zmm27                  #44.3
        vshufps   $78, %zmm13, %zmm19, %zmm3                    #44.3
        vshufps   $228, %zmm13, %zmm20, %zmm4                   #44.3
        vpermi2ps %zmm16, %zmm22, %zmm25                        #47.3
        vpermt2ps %zmm28, %zmm6, %zmm22                         #47.3
        vpermt2ps %zmm28, %zmm11, %zmm16                        #47.3
        prefetcht0 (%rcx,%r9)                                   #52.3
        movslq    8(%rsi), %r10                                 #53.16
        lea       (%r10,%r10,2), %rax                           #53.8
        movl      $23055, %r10d                                 #64.3
        shlq      $5, %rax                                      #53.8
        kmovw     %r10d, %k1                                    #64.3
        movl      $42480, %r10d                                 #64.3
        kmovw     %r10d, %k2                                    #64.3
        movl      $38595, %r10d                                 #83.3
        kmovw     %r10d, %k3                                    #83.3
        movl      $26940, %r10d                                 #83.3
        kmovw     %r10d, %k4                                    #83.3
        prefetcht0 (%rcx,%rax)                                  #54.3
        movslq    12(%rdi), %rdi                                #55.16
        lea       (%rdi,%rdi,2), %r10                           #55.9
        shlq      $5, %r10                                      #55.9
        prefetcht0 (%rcx,%r10)                                  #56.3
        movslq    12(%rsi), %rsi                                #57.16
        lea       (%rsi,%rsi,2), %rdi                           #57.9
        shlq      $5, %rdi                                      #57.9
        prefetcht0 (%rcx,%rdi)                                  #58.3
        vmovups   .L_2il0floatpacket.12(%rip), %zmm11           #64.3
        vmovups   (%rdx), %zmm31                                #68.3
        vmovups   144(%rdx), %zmm0                              #68.3
        vmovups   .L_2il0floatpacket.15(%rip), %zmm12           #68.3
        vmovups   .L_2il0floatpacket.14(%rip), %zmm7            #68.3
        vmovups   .L_2il0floatpacket.18(%rip), %zmm15           #68.3
        vmovups   .L_2il0floatpacket.17(%rip), %zmm14           #68.3
        vmovups   .L_2il0floatpacket.16(%rip), %zmm13           #68.3
        vmovups   64(%rdx), %zmm5                               #68.3
        vmovups   208(%rdx), %zmm6                              #68.3
        vpermps   %zmm22, %zmm11, %zmm10                        #65.3
        vpermps   %zmm25, %zmm11, %zmm21                        #64.3
        vpermps   %zmm16, %zmm11, %zmm17                        #66.3
        vaddps    %zmm10, %zmm3, %zmm3{%k1}                     #65.3
        vaddps    %zmm21, %zmm27, %zmm27{%k1}                   #64.3
        vaddps    %zmm17, %zmm4, %zmm4{%k1}                     #66.3
        vsubps    %zmm10, %zmm3, %zmm3{%k2}                     #65.3
        vsubps    %zmm21, %zmm27, %zmm27{%k2}                   #64.3
        vsubps    %zmm17, %zmm4, %zmm4{%k2}                     #66.3
        vmovups   .L_2il0floatpacket.13(%rip), %zmm10           #68.3
        vmovups   .L_2il0floatpacket.20(%rip), %zmm17           #68.3
        vmovups   .L_2il0floatpacket.19(%rip), %zmm16           #68.3
        vmovups   .L_2il0floatpacket.27(%rip), %zmm22           #68.3
        vmovaps   %zmm31, %zmm26                                #68.3
        vpermt2ps 72(%rdx), %zmm10, %zmm26                      #68.3
        vmovaps   %zmm0, %zmm24                                 #68.3
        vpermt2ps 216(%rdx), %zmm10, %zmm24                     #68.3
        vmovaps   %zmm26, %zmm9                                 #68.3
        vpermt2ps %zmm24, %zmm12, %zmm9                         #68.3
        vmovaps   %zmm26, %zmm8                                 #68.3
        vpermt2ps %zmm24, %zmm7, %zmm8                          #68.3
        vmulps    %zmm9, %zmm3, %zmm28                          #68.3
        vmulps    %zmm8, %zmm27, %zmm25                         #68.3
        vmovups   .L_2il0floatpacket.23(%rip), %zmm9            #68.3
        vmovups   .L_2il0floatpacket.21(%rip), %zmm8            #68.3
        vpermt2ps 72(%rdx), %zmm9, %zmm31                       #68.3
        vpermt2ps 216(%rdx), %zmm9, %zmm0                       #68.3
        vpermilps $177, %zmm3, %zmm2                            #68.3
        vmulps    %zmm2, %zmm15, %zmm1                          #68.3
        vpermilps $177, %zmm27, %zmm20                          #68.3
        vmovaps   %zmm26, %zmm19                                #68.3
        vmulps    %zmm15, %zmm20, %zmm30                        #68.3
        vmovups   .L_2il0floatpacket.25(%rip), %zmm20           #68.3
        vpermt2ps %zmm24, %zmm14, %zmm19                        #68.3
        vmovaps   %zmm26, %zmm18                                #68.3
        vpermt2ps %zmm24, %zmm13, %zmm18                        #68.3
        vfmadd231ps %zmm27, %zmm19, %zmm28                      #68.3
        vmovups   .L_2il0floatpacket.24(%rip), %zmm19           #68.3
        vfmadd231ps %zmm3, %zmm18, %zmm25                       #68.3
        vmovups   .L_2il0floatpacket.22(%rip), %zmm18           #68.3
        vmovaps   %zmm26, %zmm29                                #68.3
        vpermt2ps %zmm24, %zmm17, %zmm29                        #68.3
        vmovaps   %zmm26, %zmm23                                #68.3
        vpermt2ps %zmm24, %zmm16, %zmm23                        #68.3
        vfmadd231ps %zmm1, %zmm29, %zmm28                       #68.3
        vfmadd231ps %zmm30, %zmm23, %zmm25                      #68.3
        vmovaps   %zmm26, %zmm21                                #68.3
        vpermt2ps %zmm24, %zmm18, %zmm26                        #68.3
        vpermt2ps %zmm24, %zmm8, %zmm21                         #68.3
        vfmadd231ps %zmm30, %zmm26, %zmm28                      #68.3
        vfmadd231ps %zmm1, %zmm21, %zmm25                       #68.3
        vmovups   .L_2il0floatpacket.26(%rip), %zmm21           #68.3
        vmovaps   %zmm31, %zmm26                                #68.3
        vpermt2ps %zmm0, %zmm20, %zmm26                         #68.3
        vmulps    %zmm26, %zmm27, %zmm26                        #68.3
        vmovaps   %zmm31, %zmm27                                #68.3
        vmovaps   %zmm31, %zmm24                                #68.3
        vpermt2ps %zmm0, %zmm21, %zmm27                         #68.3
        vpermt2ps %zmm0, %zmm19, %zmm24                         #68.3
        vfmadd231ps %zmm4, %zmm27, %zmm28                       #68.3
        vfmadd231ps %zmm4, %zmm24, %zmm25                       #68.3
        vpermilps $177, %zmm4, %zmm27                           #68.3
        vmovaps   %zmm31, %zmm24                                #68.3
        vmulps    %zmm27, %zmm15, %zmm23                        #68.3
        vmovups   .L_2il0floatpacket.30(%rip), %zmm27           #68.3
        vpermt2ps %zmm0, %zmm22, %zmm24                         #68.3
        vfmadd213ps %zmm26, %zmm24, %zmm3                       #68.3
        vmovups   .L_2il0floatpacket.28(%rip), %zmm24           #68.3
        vmovups   .L_2il0floatpacket.29(%rip), %zmm26           #68.3
        vmovaps   %zmm31, %zmm29                                #68.3
        vmovaps   %zmm31, %zmm2                                 #68.3
        vpermt2ps %zmm0, %zmm24, %zmm29                         #68.3
        vpermt2ps %zmm0, %zmm26, %zmm2                          #68.3
        vfmadd231ps %zmm23, %zmm29, %zmm25                      #68.3
        vmovups   .L_2il0floatpacket.31(%rip), %zmm29           #68.3
        vfmadd213ps %zmm3, %zmm2, %zmm30                        #68.3
        vmovups   .L_2il0floatpacket.32(%rip), %zmm2            #68.3
        vmovaps   %zmm31, %zmm3                                 #68.3
        vpermt2ps %zmm0, %zmm27, %zmm3                          #68.3
        vpermt2ps %zmm0, %zmm29, %zmm31                         #68.3
        vpermt2ps 136(%rdx), %zmm2, %zmm5                       #68.3
        vpermt2ps 280(%rdx), %zmm2, %zmm6                       #68.3
        vfmadd231ps %zmm23, %zmm3, %zmm28                       #68.3
        vmovups   .L_2il0floatpacket.33(%rip), %zmm3            #68.3
        vfmadd213ps %zmm30, %zmm31, %zmm1                       #68.3
        vmovups   .L_2il0floatpacket.35(%rip), %ymm0            #70.3
        vmovaps   %zmm5, %zmm31                                 #68.3
        vpermt2ps %zmm6, %zmm3, %zmm31                          #68.3
        vfmadd213ps %zmm1, %zmm31, %zmm4                        #68.3
        vmovups   .L_2il0floatpacket.34(%rip), %zmm1            #68.3
        vpermt2ps %zmm6, %zmm1, %zmm5                           #68.3
        vfmadd213ps %zmm4, %zmm5, %zmm23                        #68.3
        vmovups   .L_2il0floatpacket.36(%rip), %ymm4            #70.3
        vextractf64x4 $1, %zmm25, %ymm5                         #70.3
        vmovaps   %zmm25, %zmm6                                 #70.3
        vpermps   %ymm5, %ymm0, %ymm25                          #70.3
        vpermps   %ymm6, %ymm0, %ymm30                          #70.3
        vfmadd213ps %ymm5, %ymm4, %ymm25                        #70.3
        vmovups   .L_2il0floatpacket.37(%rip), %ymm5            #70.3
        vfmadd213ps %ymm30, %ymm4, %ymm6                        #70.3
        vmovups   .L_2il0floatpacket.38(%rip), %ymm30           #70.3
        vpermilps %ymm5, %ymm25, %ymm31                         #70.3
        vfmadd213ps %ymm6, %ymm30, %ymm31                       #70.3
        vmovups   %ymm31, -112(%rbp)                            #70.3[spill]
        vmovaps   %zmm28, %zmm31                                #71.3
        vextractf64x4 $1, %zmm28, %ymm28                        #71.3
        vpermps   %ymm31, %ymm0, %ymm6                          #71.3
        vfmadd213ps %ymm6, %ymm4, %ymm31                        #71.3
        vpermps   %ymm28, %ymm0, %ymm6                          #71.3
        vfmadd213ps %ymm28, %ymm4, %ymm6                        #71.3
        vpermilps %ymm5, %ymm6, %ymm25                          #71.3
        vfmadd213ps %ymm31, %ymm30, %ymm25                      #71.3
        vmovups   %ymm25, -80(%rbp)                             #71.3[spill]
        vmovaps   %zmm23, %zmm25                                #72.3
        vextractf64x4 $1, %zmm23, %ymm23                        #72.3
        vpermps   %ymm23, %ymm0, %ymm31                         #72.3
        vpermps   %ymm25, %ymm0, %ymm28                         #72.3
        vfmadd213ps %ymm23, %ymm4, %ymm31                       #72.3
        vfmadd213ps %ymm28, %ymm4, %ymm25                       #72.3
        vpermilps %ymm5, %ymm31, %ymm4                          #72.3
        vfmadd213ps %ymm25, %ymm30, %ymm4                       #72.3
        vmovups   16(%rcx,%r9), %xmm25                          #76.3
        vmovups   (%rcx,%r9), %xmm5                             #76.3
        vmovups   32(%rcx,%r9), %xmm30                          #76.3
        vmovups   %ymm4, -48(%rbp)                              #72.3[spill]
        vinsertf32x4 $1, 16(%rcx,%rax), %zmm25, %zmm23          #76.3
        vinsertf32x4 $2, 16(%rcx,%r10), %zmm23, %zmm31          #76.3
        vinsertf32x4 $1, (%rcx,%rax), %zmm5, %zmm6              #76.3
        vinsertf32x4 $2, (%rcx,%r10), %zmm6, %zmm28             #76.3
        vinsertf32x4 $3, 16(%rcx,%rdi), %zmm31, %zmm6           #76.3
        vmovups   48(%rcx,%r9), %xmm31                          #79.3
        vinsertf32x4 $3, (%rcx,%rdi), %zmm28, %zmm5             #76.3
        vshufps   $228, %zmm6, %zmm5, %zmm23                    #76.3
        vinsertf32x4 $1, 32(%rcx,%rax), %zmm30, %zmm4           #76.3
        vinsertf32x4 $2, 32(%rcx,%r10), %zmm4, %zmm28           #76.3
        vinsertf32x4 $3, 32(%rcx,%rdi), %zmm28, %zmm25          #76.3
        vshufps   $78, %zmm25, %zmm5, %zmm5                     #76.3
        vshufps   $228, %zmm25, %zmm6, %zmm6                    #76.3
        vmovups   64(%rcx,%r9), %xmm25                          #79.3
        vinsertf32x4 $1, 48(%rcx,%rax), %zmm31, %zmm30          #79.3
        vinsertf32x4 $2, 48(%rcx,%r10), %zmm30, %zmm4           #79.3
        vinsertf32x4 $3, 48(%rcx,%rdi), %zmm4, %zmm28           #79.3
        vmovups   80(%rcx,%r9), %xmm4                           #79.3
        vinsertf32x4 $1, 64(%rcx,%rax), %zmm25, %zmm31          #79.3
        vinsertf32x4 $2, 64(%rcx,%r10), %zmm31, %zmm30          #79.3
        vinsertf32x4 $3, 64(%rcx,%rdi), %zmm30, %zmm30          #79.3
        vinsertf32x4 $1, 80(%rcx,%rax), %zmm4, %zmm25           #79.3
        vinsertf32x4 $2, 80(%rcx,%r10), %zmm25, %zmm31          #79.3
        vmovups   .L_2il0floatpacket.39(%rip), %zmm25           #79.3
        vinsertf32x4 $3, 80(%rcx,%rdi), %zmm31, %zmm31          #79.3
        vmovaps   %zmm28, %zmm4                                 #79.3
        vpermt2ps %zmm30, %zmm25, %zmm4                         #79.3
        vpermt2ps %zmm31, %zmm25, %zmm30                        #79.3
        vpermps   %zmm4, %zmm11, %zmm25                         #83.3
        vpermps   %zmm30, %zmm11, %zmm30                        #85.3
        vmovups   496(%rdx), %zmm4                              #91.3
        vaddps    %zmm25, %zmm23, %zmm23{%k3}                   #83.3
        vaddps    %zmm30, %zmm6, %zmm6{%k3}                     #85.3
        vpermt2ps 568(%rdx), %zmm2, %zmm4                       #91.3
        vsubps    %zmm25, %zmm23, %zmm23{%k4}                   #83.3
        vsubps    %zmm30, %zmm6, %zmm6{%k4}                     #85.3
        vmovups   .L_2il0floatpacket.40(%rip), %zmm25           #79.3
        vpermt2ps %zmm31, %zmm25, %zmm28                        #79.3
        vmovups   288(%rdx), %zmm25                             #91.3
        vpermps   %zmm28, %zmm11, %zmm11                        #84.3
        vmovups   352(%rdx), %zmm28                             #91.3
        vaddps    %zmm11, %zmm5, %zmm5{%k3}                     #84.3
        vpermt2ps 424(%rdx), %zmm2, %zmm28                      #91.3
        vsubps    %zmm11, %zmm5, %zmm5{%k4}                     #84.3
        vmovups   432(%rdx), %zmm11                             #91.3
        vpermi2ps %zmm4, %zmm28, %zmm3                          #91.3
        vpermt2ps %zmm4, %zmm1, %zmm28                          #91.3
        vmovaps   %zmm25, %zmm2                                 #91.3
        vpermt2ps 360(%rdx), %zmm10, %zmm2                      #91.3
        vpermi2ps 504(%rdx), %zmm11, %zmm10                     #91.3
        vpermt2ps 504(%rdx), %zmm9, %zmm11                      #91.3
        vpermt2ps 360(%rdx), %zmm9, %zmm25                      #91.3
        vpermi2ps %zmm10, %zmm2, %zmm7                          #91.3
        vpermi2ps %zmm10, %zmm2, %zmm12                         #91.3
        vpermi2ps %zmm10, %zmm2, %zmm13                         #91.3
        vpermi2ps %zmm10, %zmm2, %zmm14                         #91.3
        vpermi2ps %zmm10, %zmm2, %zmm17                         #91.3
        vpermi2ps %zmm10, %zmm2, %zmm16                         #91.3
        vpermi2ps %zmm10, %zmm2, %zmm8                          #91.3
        vpermt2ps %zmm10, %zmm18, %zmm2                         #91.3
        vpermi2ps %zmm11, %zmm25, %zmm20                        #91.3
        vpermi2ps %zmm11, %zmm25, %zmm22                        #91.3
        vpermi2ps %zmm11, %zmm25, %zmm26                        #91.3
        vpermi2ps %zmm11, %zmm25, %zmm19                        #91.3
        vpermi2ps %zmm11, %zmm25, %zmm21                        #91.3
        vpermi2ps %zmm11, %zmm25, %zmm24                        #91.3
        vpermi2ps %zmm11, %zmm25, %zmm27                        #91.3
        vpermt2ps %zmm11, %zmm29, %zmm25                        #91.3
        vmulps    %zmm7, %zmm23, %zmm7                          #91.3
        vmulps    %zmm12, %zmm5, %zmm12                         #91.3
        vmovups   .L_2il0floatpacket.42(%rip), %ymm18           #96.3
        vfmadd231ps %zmm5, %zmm13, %zmm7                        #91.3
        vfmadd231ps %zmm23, %zmm14, %zmm12                      #91.3
        vpermilps $177, %zmm5, %zmm13                           #91.3
        vmulps    %zmm13, %zmm15, %zmm4                         #91.3
        vmovups   .L_2il0floatpacket.41(%rip), %ymm13           #96.3
        vfmadd231ps %zmm4, %zmm17, %zmm12                       #91.3
        vpermilps $177, %zmm23, %zmm1                           #91.3
        vmulps    %zmm15, %zmm1, %zmm1                          #91.3
        vfmadd231ps %zmm1, %zmm2, %zmm12                        #91.3
        vfmadd231ps %zmm1, %zmm16, %zmm7                        #91.3
        vmulps    %zmm20, %zmm23, %zmm2                         #91.3
        vfmadd231ps %zmm4, %zmm8, %zmm7                         #91.3
        vfmadd231ps %zmm6, %zmm21, %zmm12                       #91.3
        vmovups   .L_2il0floatpacket.43(%rip), %ymm21           #96.3
        vfmadd213ps %zmm2, %zmm22, %zmm5                        #91.3
        vfmadd231ps %zmm6, %zmm19, %zmm7                        #91.3
        vmovups   .L_2il0floatpacket.36(%rip), %ymm19           #96.3
        vmovups   .L_2il0floatpacket.44(%rip), %ymm22           #96.3
        vfmadd213ps %zmm5, %zmm26, %zmm1                        #91.3
        vbroadcastss -16(%rbp), %ymm26                          #94.10
        vfmadd213ps %zmm1, %zmm25, %zmm4                        #91.3
        vpermilps $177, %zmm6, %zmm8                            #91.3
        vmulps    %zmm8, %zmm15, %zmm15                         #91.3
        vfmadd213ps %zmm4, %zmm3, %zmm6                         #91.3
        vfmadd231ps %zmm15, %zmm24, %zmm7                       #91.3
        vfmadd231ps %zmm15, %zmm27, %zmm12                      #91.3
        vfmadd213ps %zmm6, %zmm28, %zmm15                       #91.3
        vpermps   %ymm7, %ymm0, %ymm3                           #96.3
        vpermps   %ymm12, %ymm0, %ymm10                         #97.3
        vpermps   %ymm15, %ymm0, %ymm17                         #98.3
        vfmadd213ps %ymm7, %ymm19, %ymm3                        #96.3
        vfmadd213ps %ymm12, %ymm19, %ymm10                      #97.3
        vfmadd213ps %ymm15, %ymm19, %ymm17                      #98.3
        vextractf64x4 $1, %zmm7, %ymm5                          #96.3
        vextractf64x4 $1, %zmm12, %ymm11                        #97.3
        vextractf64x4 $1, %zmm15, %ymm20                        #98.3
        vpermps   %ymm5, %ymm0, %ymm6                           #96.3
        vpermps   %ymm11, %ymm0, %ymm14                         #97.3
        vpermps   %ymm20, %ymm0, %ymm0                          #98.3
        vpermilps %ymm13, %ymm3, %ymm9                          #96.3
        vfmadd213ps %ymm5, %ymm19, %ymm6                        #96.3
        vfmadd213ps %ymm11, %ymm19, %ymm14                      #97.3
        vfmadd213ps %ymm20, %ymm19, %ymm0                       #98.3
        vfmadd213ps -112(%rbp), %ymm18, %ymm9                   #96.3[spill]
        vpermilps %ymm13, %ymm10, %ymm16                        #97.3
        vpermilps %ymm13, %ymm17, %ymm23                        #98.3
        vfmadd213ps -80(%rbp), %ymm18, %ymm16                   #97.3[spill]
        vfmadd213ps -48(%rbp), %ymm18, %ymm23                   #98.3[spill]
        vpermilps %ymm21, %ymm6, %ymm24                         #96.3
        vpermilps %ymm21, %ymm14, %ymm25                        #97.3
        vpermilps %ymm21, %ymm0, %ymm27                         #98.3
        vfmadd213ps %ymm9, %ymm22, %ymm24                       #96.3
        vfmadd213ps %ymm16, %ymm22, %ymm25                      #97.3
        vfmadd213ps %ymm23, %ymm22, %ymm27                      #98.3
        vmulps    %ymm26, %ymm24, %ymm29                        #100.8
        vmulps    %ymm25, %ymm26, %ymm31                        #101.8
        vmulps    %ymm27, %ymm26, %ymm0                         #102.8
        vshufps   $68, %ymm31, %ymm29, %ymm28                   #104.3
        vshufps   $228, %ymm29, %ymm0, %ymm30                   #104.3
        vshufps   $238, %ymm0, %ymm31, %ymm1                    #104.3
        vmovups   %xmm28, (%r8)                                 #104.3
        vmovups   %xmm30, 16(%r8)                               #104.3
        vmovups   %xmm1, 32(%r8)                                #104.3
        vextractf32x4 $1, %ymm28, 48(%r8)                       #104.3
        vextractf32x4 $1, %ymm30, 64(%r8)                       #104.3
        vextractf128 $1, %ymm1, 80(%r8)                         #104.3
        vzeroupper                                              #105.1
        movq      %rbp, %rsp                                    #105.1
        popq      %rbp                                          #105.1
	.cfi_restore 6
        movq      %rbx, %rsp                                    #105.1
        popq      %rbx                                          #105.1
	.cfi_def_cfa 7, 8
	.cfi_restore 3
        ret                                                     #105.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	doe_avx512,@function
	.size	doe_avx512,.-doe_avx512
	.data
# -- End  doe_avx512
	.text
# -- Begin  deo_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl deo_avx512
# --- deo_avx512(int *, int *, su3 *, spinor *, float, spin_t *)
deo_avx512:
# parameter 1: %rdi
# parameter 2: %rsi
# parameter 3: %rdx
# parameter 4: %rcx
# parameter 5: %xmm0
# parameter 6: %r8
..B2.1:                         # Preds ..B2.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_deo_avx512.11:
..L12:
                                                         #108.1
        pushq     %rbp                                          #108.1
	.cfi_def_cfa_offset 16
        movq      %rsp, %rbp                                    #108.1
	.cfi_def_cfa 6, 16
	.cfi_offset 6, -16
        movslq    (%rdi), %rax                                  #122.16
        vmovss    %xmm0, -16(%rbp)                              #108.1
        lea       (%rax,%rax,2), %r9                            #122.8
        shlq      $5, %r9                                       #122.8
        prefetcht0 (%rcx,%r9)                                   #123.3
        movl      $42255, %r10d                                 #165.3
        movslq    (%rsi), %r11                                  #124.16
        kmovw     %r10d, %k1                                    #165.3
        movl      $23280, %r10d                                 #165.3
        kmovw     %r10d, %k2                                    #165.3
        movl      $38595, %r10d                                 #183.3
        lea       (%r11,%r11,2), %rax                           #124.8
        shlq      $5, %rax                                      #124.8
        kmovw     %r10d, %k3                                    #183.3
        movl      $26940, %r10d                                 #183.3
        kmovw     %r10d, %k4                                    #183.3
        prefetcht0 (%rcx,%rax)                                  #125.3
        movslq    4(%rdi), %r10                                 #126.17
        lea       (%r10,%r10,2), %r11                           #126.9
        shlq      $5, %r11                                      #126.9
        prefetcht0 (%rcx,%r11)                                  #127.3
        movslq    4(%rsi), %r10                                 #128.17
        lea       (%r10,%r10,2), %r10                           #128.9
        shlq      $5, %r10                                      #128.9
        prefetcht0 (%rcx,%r10)                                  #129.3
        vmovups   (%r8), %xmm2                                  #131.3
        vmovups   16(%r8), %xmm6                                #131.3
        vmovups   32(%r8), %xmm4                                #131.3
        vbroadcastss -16(%rbp), %ymm27                          #133.10
        vmovups   .L_2il0floatpacket.35(%rip), %ymm23           #138.3
        vmovups   .L_2il0floatpacket.36(%rip), %ymm28           #138.3
        vmovups   .L_2il0floatpacket.47(%rip), %ymm18           #142.3
        vmovups   216(%rdx), %zmm9                              #150.3
        vmovups   .L_2il0floatpacket.15(%rip), %zmm26           #150.3
        vmovups   .L_2il0floatpacket.16(%rip), %zmm8            #150.3
        vinsertf128 $1, 48(%r8), %ymm2, %ymm15                  #131.3
        vinsertf128 $1, 64(%r8), %ymm6, %ymm13                  #131.3
        vshufps   $228, %ymm13, %ymm15, %ymm21                  #131.3
        vmulps    %ymm27, %ymm21, %ymm20                        #134.8
        vpermps   %ymm20, %ymm23, %ymm24                        #138.3
        vfmadd231ps %ymm20, %ymm28, %ymm24                      #138.3
        vinsertf128 $1, 80(%r8), %ymm4, %ymm1                   #131.3
        vshufps   $78, %ymm1, %ymm15, %ymm30                    #131.3
        vshufps   $228, %ymm1, %ymm13, %ymm31                   #131.3
        vmovups   .L_2il0floatpacket.45(%rip), %ymm13           #142.3
        vmovups   .L_2il0floatpacket.46(%rip), %ymm1            #142.3
        vmulps    %ymm30, %ymm27, %ymm25                        #135.8
        vmulps    %ymm31, %ymm27, %ymm14                        #136.8
        vpermps   %ymm20, %ymm13, %ymm15                        #142.3
        vpermps   %ymm20, %ymm1, %ymm21                         #142.3
        vfmadd213ps %ymm21, %ymm18, %ymm15                      #142.3
        vmovups   .L_2il0floatpacket.14(%rip), %zmm21           #150.3
        vpermps   %ymm25, %ymm1, %ymm16                         #143.3
        vpermps   %ymm14, %ymm1, %ymm12                         #144.3
        vmovaps   %zmm9, %zmm1                                  #150.3
        vpermps   %ymm14, %ymm23, %ymm29                        #140.3
        vpermps   %ymm14, %ymm13, %ymm30                        #144.3
        vfmadd231ps %ymm14, %ymm28, %ymm29                      #140.3
        vfmadd213ps %ymm12, %ymm18, %ymm30                      #144.3
        vpermps   %ymm25, %ymm23, %ymm22                        #139.3
        vpermps   %ymm25, %ymm13, %ymm17                        #143.3
        vfmadd231ps %ymm25, %ymm28, %ymm22                      #139.3
        vfmadd213ps %ymm16, %ymm18, %ymm17                      #143.3
        vmovups   .L_2il0floatpacket.25(%rip), %zmm28           #150.3
        vmovups   280(%rdx), %zmm13                             #150.3
        movslq    8(%rdi), %r8                                  #174.16
        vinsertf64x4 $1, %ymm15, %zmm24, %zmm10                 #142.3
        lea       (%r8,%r8,2), %r8                              #174.8
        vmovups   72(%rdx), %zmm15                              #150.3
        vmovups   .L_2il0floatpacket.13(%rip), %zmm24           #150.3
        vmovaps   %zmm15, %zmm7                                 #150.3
        vpermt2ps (%rdx), %zmm24, %zmm7                         #150.3
        vpermt2ps 144(%rdx), %zmm24, %zmm1                      #150.3
        vmovaps   %zmm7, %zmm5                                  #150.3
        vpermt2ps %zmm1, %zmm21, %zmm5                          #150.3
        vmulps    %zmm5, %zmm10, %zmm18                         #150.3
        vpermilps $177, %zmm10, %zmm31                          #150.3
        vmovaps   %zmm7, %zmm0                                  #150.3
        vmovaps   %zmm7, %zmm11                                 #150.3
        vpermt2ps %zmm1, %zmm26, %zmm0                          #150.3
        vpermt2ps %zmm1, %zmm8, %zmm11                          #150.3
        vmovaps   %zmm7, %zmm6                                  #150.3
        vmovaps   %zmm7, %zmm23                                 #150.3
        vmovaps   %zmm7, %zmm27                                 #150.3
        vmovaps   %zmm7, %zmm4                                  #150.3
        shlq      $5, %r8                                       #174.8
        vinsertf64x4 $1, %ymm30, %zmm29, %zmm12                 #144.3
        vmovups   .L_2il0floatpacket.18(%rip), %zmm30           #150.3
        vmovups   .L_2il0floatpacket.17(%rip), %zmm29           #150.3
        vmulps    %zmm30, %zmm31, %zmm5                         #150.3
        vmovups   .L_2il0floatpacket.27(%rip), %zmm31           #150.3
        vpermt2ps %zmm1, %zmm29, %zmm27                         #150.3
        vinsertf64x4 $1, %ymm17, %zmm22, %zmm3                  #143.3
        vmovups   .L_2il0floatpacket.21(%rip), %zmm22           #150.3
        vmovups   136(%rdx), %zmm17                             #150.3
        vfmadd231ps %zmm3, %zmm11, %zmm18                       #150.3
        vmulps    %zmm0, %zmm3, %zmm19                          #150.3
        vpermt2ps %zmm1, %zmm22, %zmm23                         #150.3
        vfmadd231ps %zmm10, %zmm27, %zmm19                      #150.3
        vmovups   .L_2il0floatpacket.26(%rip), %zmm27           #150.3
        vpermilps $177, %zmm3, %zmm2                            #150.3
        vmulps    %zmm2, %zmm30, %zmm16                         #150.3
        vmovups   .L_2il0floatpacket.19(%rip), %zmm2            #150.3
        vpermt2ps %zmm1, %zmm2, %zmm6                           #150.3
        vfmadd231ps %zmm5, %zmm6, %zmm18                        #150.3
        vmovups   .L_2il0floatpacket.20(%rip), %zmm6            #150.3
        vfmadd231ps %zmm16, %zmm23, %zmm18                      #150.3
        vmovups   .L_2il0floatpacket.23(%rip), %zmm23           #150.3
        vpermt2ps %zmm1, %zmm6, %zmm4                           #150.3
        vpermt2ps (%rdx), %zmm23, %zmm15                        #150.3
        vpermt2ps 144(%rdx), %zmm23, %zmm9                      #150.3
        vfmadd231ps %zmm16, %zmm4, %zmm19                       #150.3
        vmovups   .L_2il0floatpacket.22(%rip), %zmm4            #150.3
        vmovaps   %zmm15, %zmm0                                 #150.3
        vpermt2ps %zmm9, %zmm28, %zmm0                          #150.3
        vpermt2ps %zmm1, %zmm4, %zmm7                           #150.3
        vmovups   .L_2il0floatpacket.24(%rip), %zmm1            #150.3
        vmulps    %zmm0, %zmm10, %zmm0                          #150.3
        vfmadd231ps %zmm5, %zmm7, %zmm19                        #150.3
        vmovaps   %zmm15, %zmm10                                #150.3
        vmovaps   %zmm15, %zmm7                                 #150.3
        vpermt2ps %zmm9, %zmm27, %zmm10                         #150.3
        vpermt2ps %zmm9, %zmm1, %zmm7                           #150.3
        vfmadd231ps %zmm12, %zmm10, %zmm19                      #150.3
        vfmadd231ps %zmm12, %zmm7, %zmm18                       #150.3
        vpermilps $177, %zmm12, %zmm10                          #150.3
        vmovaps   %zmm15, %zmm7                                 #150.3
        vmulps    %zmm10, %zmm30, %zmm11                        #150.3
        vmovups   .L_2il0floatpacket.29(%rip), %zmm10           #150.3
        vpermt2ps %zmm9, %zmm31, %zmm7                          #150.3
        vfmadd213ps %zmm0, %zmm7, %zmm3                         #150.3
        vmovups   .L_2il0floatpacket.28(%rip), %zmm7            #150.3
        vmovaps   %zmm15, %zmm0                                 #150.3
        vpermt2ps %zmm9, %zmm7, %zmm0                           #150.3
        vfmadd231ps %zmm11, %zmm0, %zmm18                       #150.3
        vmovaps   %zmm15, %zmm0                                 #150.3
        vpermt2ps %zmm9, %zmm10, %zmm0                          #150.3
        vfmadd213ps %zmm3, %zmm0, %zmm5                         #150.3
        vmovups   .L_2il0floatpacket.30(%rip), %zmm3            #150.3
        vmovaps   %zmm15, %zmm0                                 #150.3
        vpermt2ps %zmm9, %zmm3, %zmm0                           #150.3
        vfmadd231ps %zmm11, %zmm0, %zmm19                       #150.3
        vmovups   .L_2il0floatpacket.31(%rip), %zmm0            #150.3
        vpermt2ps %zmm9, %zmm0, %zmm15                          #150.3
        vmovups   .L_2il0floatpacket.32(%rip), %zmm9            #150.3
        vfmadd213ps %zmm5, %zmm15, %zmm16                       #150.3
        vmovups   .L_2il0floatpacket.33(%rip), %zmm5            #150.3
        vpermt2ps 64(%rdx), %zmm9, %zmm17                       #150.3
        vpermt2ps 208(%rdx), %zmm9, %zmm13                      #150.3
        vmovaps   %zmm17, %zmm15                                #150.3
        vpermt2ps %zmm13, %zmm5, %zmm15                         #150.3
        vfmadd213ps %zmm16, %zmm15, %zmm12                      #150.3
        vmovups   .L_2il0floatpacket.34(%rip), %zmm16           #150.3
        vmovups   16(%rcx,%rax), %xmm15                         #152.3
        vpermt2ps %zmm13, %zmm16, %zmm17                        #150.3
        vfmadd213ps %zmm12, %zmm17, %zmm11                      #150.3
        vmovups   (%rcx,%rax), %xmm17                           #152.3
        vinsertf32x4 $1, (%rcx,%r9), %zmm17, %zmm13             #152.3
        vinsertf32x4 $2, (%rcx,%r10), %zmm13, %zmm12            #152.3
        vinsertf32x4 $3, (%rcx,%r11), %zmm12, %zmm13            #152.3
        vmovups   32(%rcx,%rax), %xmm12                         #152.3
        vinsertf32x4 $1, 16(%rcx,%r9), %zmm15, %zmm16           #152.3
        vinsertf32x4 $2, 16(%rcx,%r10), %zmm16, %zmm17          #152.3
        vinsertf32x4 $3, 16(%rcx,%r11), %zmm17, %zmm17          #152.3
        vinsertf32x4 $1, 32(%rcx,%r9), %zmm12, %zmm15           #152.3
        vinsertf32x4 $2, 32(%rcx,%r10), %zmm15, %zmm16          #152.3
        vshufps   $228, %zmm17, %zmm13, %zmm15                  #152.3
        vinsertf32x4 $3, 32(%rcx,%r11), %zmm16, %zmm12          #152.3
        vshufps   $78, %zmm12, %zmm13, %zmm16                   #152.3
        vshufps   $228, %zmm12, %zmm17, %zmm13                  #152.3
        vaddps    %zmm18, %zmm15, %zmm12                        #155.8
        vaddps    %zmm19, %zmm16, %zmm17                        #156.8
        vaddps    %zmm11, %zmm13, %zmm13                        #157.8
        vshufps   $68, %zmm17, %zmm12, %zmm15                   #158.3
        vshufps   $228, %zmm12, %zmm13, %zmm16                  #158.3
        vshufps   $238, %zmm13, %zmm17, %zmm17                  #158.3
        vmovups   %xmm15, (%rcx,%rax)                           #158.3
        vextractf32x4 $1, %zmm15, (%rcx,%r9)                    #158.3
        vextractf32x4 $2, %zmm15, (%rcx,%r10)                   #158.3
        vextractf32x4 $3, %zmm15, (%rcx,%r11)                   #158.3
        vmovups   %xmm16, 16(%rcx,%rax)                         #158.3
        vextractf32x4 $1, %zmm16, 16(%rcx,%r9)                  #158.3
        vextractf32x4 $2, %zmm16, 16(%rcx,%r10)                 #158.3
        vextractf32x4 $3, %zmm16, 16(%rcx,%r11)                 #158.3
        vmovups   %xmm17, 32(%rcx,%rax)                         #158.3
        vextractf32x4 $1, %zmm17, 32(%rcx,%r9)                  #158.3
        vextractf32x4 $2, %zmm17, 32(%rcx,%r10)                 #158.3
        vextractf32x4 $3, %zmm17, 32(%rcx,%r11)                 #158.3
        vmovups   48(%rcx,%rax), %xmm12                         #162.3
        vmovups   64(%rcx,%rax), %xmm13                         #162.3
        vmovups   80(%rcx,%rax), %xmm17                         #162.3
        vinsertf32x4 $1, 48(%rcx,%r9), %zmm12, %zmm15           #162.3
        vinsertf32x4 $2, 48(%rcx,%r10), %zmm15, %zmm16          #162.3
        vinsertf32x4 $3, 48(%rcx,%r11), %zmm16, %zmm12          #162.3
        vinsertf32x4 $1, 64(%rcx,%r9), %zmm13, %zmm15           #162.3
        vinsertf32x4 $2, 64(%rcx,%r10), %zmm15, %zmm16          #162.3
        vinsertf32x4 $3, 64(%rcx,%r11), %zmm16, %zmm13          #162.3
        vinsertf32x4 $1, 80(%rcx,%r9), %zmm17, %zmm15           #162.3
        vinsertf32x4 $2, 80(%rcx,%r10), %zmm15, %zmm16          #162.3
        vinsertf32x4 $3, 80(%rcx,%r11), %zmm16, %zmm17          #162.3
        vmovups   .L_2il0floatpacket.10(%rip), %zmm16           #162.3
        prefetcht0 (%rcx,%r8)                                   #175.3
        vmovaps   %zmm12, %zmm15                                #162.3
        vpermt2ps %zmm13, %zmm16, %zmm15                        #162.3
        vpermt2ps %zmm17, %zmm16, %zmm13                        #162.3
        vmovups   .L_2il0floatpacket.11(%rip), %zmm16           #162.3
        vpermt2ps %zmm17, %zmm16, %zmm12                        #162.3
        vmovups   .L_2il0floatpacket.12(%rip), %zmm16           #165.3
        vpermps   %zmm18, %zmm16, %zmm18                        #165.3
        vpermps   %zmm11, %zmm16, %zmm17                        #167.3
        vpermps   %zmm19, %zmm16, %zmm19                        #166.3
        vaddps    %zmm18, %zmm15, %zmm15{%k1}                   #165.3
        vaddps    %zmm17, %zmm13, %zmm13{%k1}                   #167.3
        vaddps    %zmm19, %zmm12, %zmm12{%k1}                   #166.3
        vsubps    %zmm18, %zmm15, %zmm15{%k2}                   #165.3
        vsubps    %zmm17, %zmm13, %zmm13{%k2}                   #167.3
        vsubps    %zmm19, %zmm12, %zmm12{%k2}                   #166.3
        vmovups   .L_2il0floatpacket.49(%rip), %zmm11           #168.3
        vmovups   .L_2il0floatpacket.48(%rip), %zmm17           #168.3
        vpermi2ps %zmm15, %zmm13, %zmm11                        #168.3
        vmovaps   %zmm15, %zmm18                                #168.3
        vmovups   .L_2il0floatpacket.50(%rip), %zmm15           #168.3
        vpermt2ps %zmm12, %zmm17, %zmm18                        #168.3
        vpermt2ps %zmm13, %zmm15, %zmm12                        #168.3
        vmovups   %xmm18, 48(%rcx,%rax)                         #168.3
        vextractf32x4 $1, %zmm18, 48(%rcx,%r9)                  #168.3
        vextractf32x4 $2, %zmm18, 48(%rcx,%r10)                 #168.3
        vextractf32x4 $3, %zmm18, 48(%rcx,%r11)                 #168.3
        vmovups   %xmm11, 64(%rcx,%rax)                         #168.3
        vextractf32x4 $1, %zmm11, 64(%rcx,%r9)                  #168.3
        vextractf32x4 $2, %zmm11, 64(%rcx,%r10)                 #168.3
        vextractf32x4 $3, %zmm11, 64(%rcx,%r11)                 #168.3
        vmovups   %xmm12, 80(%rcx,%rax)                         #168.3
        vextractf32x4 $1, %zmm12, 80(%rcx,%r9)                  #168.3
        vextractf32x4 $2, %zmm12, 80(%rcx,%r10)                 #168.3
        vextractf32x4 $3, %zmm12, 80(%rcx,%r11)                 #168.3
        movslq    8(%rsi), %r9                                  #176.16
        lea       (%r9,%r9,2), %rax                             #176.8
        shlq      $5, %rax                                      #176.8
        prefetcht0 (%rcx,%rax)                                  #177.3
        movslq    12(%rdi), %rdi                                #178.17
        lea       (%rdi,%rdi,2), %rdi                           #178.9
        shlq      $5, %rdi                                      #178.9
        prefetcht0 (%rcx,%rdi)                                  #179.3
        movslq    12(%rsi), %rsi                                #180.17
        lea       (%rsi,%rsi,2), %rsi                           #180.9
        shlq      $5, %rsi                                      #180.9
        prefetcht0 (%rcx,%rsi)                                  #181.3
        vmovups   .L_2il0floatpacket.51(%rip), %zmm13           #183.3
        vmovups   .L_2il0floatpacket.52(%rip), %zmm18           #183.3
        vmovups   424(%rdx), %zmm12                             #191.3
        vpermps   %zmm20, %zmm13, %zmm11                        #183.3
        vpermps   %zmm20, %zmm18, %zmm20                        #183.3
        vpermt2ps 352(%rdx), %zmm9, %zmm12                      #191.3
        vaddps    %zmm11, %zmm20, %zmm19{%k3}{z}                #183.3
        vsubps    %zmm11, %zmm20, %zmm19{%k4}                   #183.3
        vpermps   %zmm25, %zmm13, %zmm20                        #184.3
        vpermps   %zmm25, %zmm18, %zmm25                        #184.3
        vpermps   %zmm14, %zmm13, %zmm13                        #185.3
        vpermps   %zmm14, %zmm18, %zmm14                        #185.3
        vaddps    %zmm20, %zmm25, %zmm11{%k3}{z}                #184.3
        vmovups   568(%rdx), %zmm18                             #191.3
        vsubps    %zmm20, %zmm25, %zmm11{%k4}                   #184.3
        vaddps    %zmm13, %zmm14, %zmm20{%k3}{z}                #185.3
        vpermt2ps 496(%rdx), %zmm9, %zmm18                      #191.3
        vmovups   360(%rdx), %zmm25                             #191.3
        vsubps    %zmm13, %zmm14, %zmm20{%k4}                   #185.3
        vpermi2ps %zmm18, %zmm12, %zmm5                         #191.3
        vmovups   504(%rdx), %zmm14                             #191.3
        vmovaps   %zmm25, %zmm9                                 #191.3
        vpermt2ps 288(%rdx), %zmm24, %zmm9                      #191.3
        vpermi2ps 432(%rdx), %zmm14, %zmm24                     #191.3
        vpermt2ps 432(%rdx), %zmm23, %zmm14                     #191.3
        vpermt2ps 288(%rdx), %zmm23, %zmm25                     #191.3
        vpermi2ps %zmm24, %zmm9, %zmm21                         #191.3
        vpermi2ps %zmm24, %zmm9, %zmm26                         #191.3
        vpermi2ps %zmm24, %zmm9, %zmm8                          #191.3
        vpermi2ps %zmm24, %zmm9, %zmm29                         #191.3
        vpermi2ps %zmm24, %zmm9, %zmm6                          #191.3
        vpermi2ps %zmm24, %zmm9, %zmm2                          #191.3
        vpermi2ps %zmm24, %zmm9, %zmm22                         #191.3
        vpermt2ps %zmm24, %zmm4, %zmm9                          #191.3
        vpermi2ps %zmm14, %zmm25, %zmm28                        #191.3
        vpermi2ps %zmm14, %zmm25, %zmm31                        #191.3
        vpermi2ps %zmm14, %zmm25, %zmm1                         #191.3
        vpermi2ps %zmm14, %zmm25, %zmm10                        #191.3
        vpermi2ps %zmm14, %zmm25, %zmm27                        #191.3
        vpermi2ps %zmm14, %zmm25, %zmm7                         #191.3
        vpermi2ps %zmm14, %zmm25, %zmm3                         #191.3
        vpermt2ps %zmm14, %zmm0, %zmm25                         #191.3
        vmovups   (%rcx,%rax), %xmm0                            #193.3
        vmulps    %zmm21, %zmm19, %zmm21                        #191.3
        vmulps    %zmm26, %zmm11, %zmm13                        #191.3
        vmovups   16(%rcx,%rax), %xmm4                          #193.3
        vfmadd231ps %zmm11, %zmm8, %zmm21                       #191.3
        vfmadd231ps %zmm19, %zmm29, %zmm13                      #191.3
        vpermilps $177, %zmm11, %zmm26                          #191.3
        vmulps    %zmm26, %zmm30, %zmm26                        #191.3
        vpermilps $177, %zmm19, %zmm8                           #191.3
        vmulps    %zmm30, %zmm8, %zmm8                          #191.3
        vfmadd231ps %zmm26, %zmm6, %zmm13                       #191.3
        vfmadd231ps %zmm8, %zmm2, %zmm21                        #191.3
        vfmadd231ps %zmm8, %zmm9, %zmm13                        #191.3
        vmulps    %zmm28, %zmm19, %zmm9                         #191.3
        vfmadd231ps %zmm26, %zmm22, %zmm21                      #191.3
        vfmadd231ps %zmm20, %zmm27, %zmm13                      #191.3
        vfmadd213ps %zmm9, %zmm31, %zmm11                       #191.3
        vfmadd231ps %zmm20, %zmm1, %zmm21                       #191.3
        vfmadd213ps %zmm11, %zmm10, %zmm8                       #191.3
        vpermilps $177, %zmm20, %zmm28                          #191.3
        vmulps    %zmm28, %zmm30, %zmm1                         #191.3
        vfmadd213ps %zmm8, %zmm25, %zmm26                       #191.3
        vfmadd231ps %zmm1, %zmm3, %zmm13                        #191.3
        vmovups   .L_2il0floatpacket.34(%rip), %zmm3            #191.3
        vfmadd231ps %zmm1, %zmm7, %zmm21                        #191.3
        vmovups   32(%rcx,%rax), %xmm7                          #193.3
        vfmadd213ps %zmm26, %zmm5, %zmm20                       #191.3
        vpermt2ps %zmm18, %zmm3, %zmm12                         #191.3
        vfmadd213ps %zmm20, %zmm12, %zmm1                       #191.3
        vinsertf32x4 $1, (%rcx,%r8), %zmm0, %zmm2               #193.3
        vinsertf32x4 $2, (%rcx,%rsi), %zmm2, %zmm3              #193.3
        vinsertf32x4 $3, (%rcx,%rdi), %zmm3, %zmm2              #193.3
        vinsertf32x4 $1, 16(%rcx,%r8), %zmm4, %zmm5             #193.3
        vinsertf32x4 $2, 16(%rcx,%rsi), %zmm5, %zmm6            #193.3
        vinsertf32x4 $3, 16(%rcx,%rdi), %zmm6, %zmm3            #193.3
        vinsertf32x4 $1, 32(%rcx,%r8), %zmm7, %zmm8             #193.3
        vinsertf32x4 $2, 32(%rcx,%rsi), %zmm8, %zmm0            #193.3
                                # LOE rax rcx rbx rsi rdi r8 r12 r13 r14 r15 zmm0 zmm1 zmm2 zmm3 zmm13 zmm15 zmm16 zmm17 zmm21
..B2.4:                         # Preds ..B2.1
                                # Execution count [1.00e+00]
        vshufps   $228, %zmm3, %zmm2, %zmm4                     #193.3
        movl      $27075, %edx                                  #200.3
        vmovups   .L_2il0floatpacket.39(%rip), %zmm26           #199.3
        vmovups   .L_2il0floatpacket.40(%rip), %zmm25           #199.3
        vinsertf32x4 $3, 32(%rcx,%rdi), %zmm0, %zmm0            #193.3
        vaddps    %zmm21, %zmm4, %zmm5                          #194.8
        vpermps   %zmm21, %zmm16, %zmm21                        #200.3
        vshufps   $78, %zmm0, %zmm2, %zmm2                      #193.3
        vshufps   $228, %zmm0, %zmm3, %zmm3                     #193.3
        kmovw     %edx, %k1                                     #200.3
        vaddps    %zmm13, %zmm2, %zmm6                          #195.8
        vaddps    %zmm1, %zmm3, %zmm7                           #196.8
        vpermps   %zmm13, %zmm16, %zmm13                        #201.3
        vpermps   %zmm1, %zmm16, %zmm1                          #202.3
        vshufps   $68, %zmm6, %zmm5, %zmm8                      #197.3
        vshufps   $228, %zmm5, %zmm7, %zmm9                     #197.3
        vshufps   $238, %zmm7, %zmm6, %zmm10                    #197.3
        vmovups   .L_2il0floatpacket.53(%rip), %zmm16           #203.3
        vmovaps   %zmm26, %zmm28                                #199.3
        movl      $38460, %edx                                  #200.3
        kmovw     %edx, %k2                                     #200.3
        vmovups   %xmm8, (%rcx,%rax)                            #197.3
        vextractf32x4 $1, %zmm8, (%rcx,%r8)                     #197.3
        vextractf32x4 $2, %zmm8, (%rcx,%rsi)                    #197.3
        vextractf32x4 $3, %zmm8, (%rcx,%rdi)                    #197.3
        vmovups   %xmm9, 16(%rcx,%rax)                          #197.3
        vextractf32x4 $1, %zmm9, 16(%rcx,%r8)                   #197.3
        vextractf32x4 $2, %zmm9, 16(%rcx,%rsi)                  #197.3
        vextractf32x4 $3, %zmm9, 16(%rcx,%rdi)                  #197.3
        vmovups   %xmm10, 32(%rcx,%rax)                         #197.3
        vextractf32x4 $1, %zmm10, 32(%rcx,%r8)                  #197.3
        vextractf32x4 $2, %zmm10, 32(%rcx,%rsi)                 #197.3
        vextractf32x4 $3, %zmm10, 32(%rcx,%rdi)                 #197.3
        vmovups   48(%rcx,%rax), %xmm11                         #199.3
        vmovups   64(%rcx,%rax), %xmm18                         #199.3
        vmovups   80(%rcx,%rax), %xmm22                         #199.3
        vinsertf32x4 $1, 48(%rcx,%r8), %zmm11, %zmm12           #199.3
        vinsertf32x4 $1, 64(%rcx,%r8), %zmm18, %zmm19           #199.3
        vinsertf32x4 $1, 80(%rcx,%r8), %zmm22, %zmm23           #199.3
        vinsertf32x4 $2, 48(%rcx,%rsi), %zmm12, %zmm14          #199.3
        vinsertf32x4 $2, 64(%rcx,%rsi), %zmm19, %zmm20          #199.3
        vinsertf32x4 $2, 80(%rcx,%rsi), %zmm23, %zmm24          #199.3
        vinsertf32x4 $3, 48(%rcx,%rdi), %zmm14, %zmm30          #199.3
        vinsertf32x4 $3, 64(%rcx,%rdi), %zmm20, %zmm29          #199.3
        vinsertf32x4 $3, 80(%rcx,%rdi), %zmm24, %zmm27          #199.3
        vpermi2ps %zmm29, %zmm30, %zmm28                        #199.3
        vpermt2ps %zmm27, %zmm25, %zmm30                        #199.3
        vpermt2ps %zmm27, %zmm26, %zmm29                        #199.3
        vaddps    %zmm21, %zmm28, %zmm28{%k1}                   #200.3
        vaddps    %zmm13, %zmm30, %zmm30{%k1}                   #201.3
        vaddps    %zmm1, %zmm29, %zmm29{%k1}                    #202.3
        vsubps    %zmm21, %zmm28, %zmm28{%k2}                   #200.3
        vsubps    %zmm13, %zmm30, %zmm30{%k2}                   #201.3
        vsubps    %zmm1, %zmm29, %zmm29{%k2}                    #202.3
        vpermi2ps %zmm30, %zmm28, %zmm15                        #203.3
        vpermi2ps %zmm28, %zmm29, %zmm16                        #203.3
        vpermt2ps %zmm29, %zmm17, %zmm30                        #203.3
        vmovups   %xmm15, 48(%rcx,%rax)                         #203.3
        vextractf32x4 $1, %zmm15, 48(%rcx,%r8)                  #203.3
        vextractf32x4 $2, %zmm15, 48(%rcx,%rsi)                 #203.3
        vextractf32x4 $3, %zmm15, 48(%rcx,%rdi)                 #203.3
        vmovups   %xmm16, 64(%rcx,%rax)                         #203.3
        vextractf32x4 $1, %zmm16, 64(%rcx,%r8)                  #203.3
        vextractf32x4 $2, %zmm16, 64(%rcx,%rsi)                 #203.3
        vextractf32x4 $3, %zmm16, 64(%rcx,%rdi)                 #203.3
        vmovups   %xmm30, 80(%rcx,%rax)                         #203.3
        vextractf32x4 $1, %zmm30, 80(%rcx,%r8)                  #203.3
        vextractf32x4 $2, %zmm30, 80(%rcx,%rsi)                 #203.3
        vextractf32x4 $3, %zmm30, 80(%rcx,%rdi)                 #203.3
        vzeroupper                                              #204.1
        movq      %rbp, %rsp                                    #204.1
        popq      %rbp                                          #204.1
	.cfi_restore 6
        ret                                                     #204.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	deo_avx512,@function
	.size	deo_avx512,.-deo_avx512
	.data
# -- End  deo_avx512
	.section .rodata, "a"
	.align 64
	.align 64
.L_2il0floatpacket.10:
	.long	0x00000000,0x00000001,0x00000012,0x00000013,0x00000004,0x00000005,0x00000016,0x00000017,0x0000001a,0x0000001b,0x00000008,0x00000009,0x0000001e,0x0000001f,0x0000000c,0x0000000d
	.type	.L_2il0floatpacket.10,@object
	.size	.L_2il0floatpacket.10,64
	.align 64
.L_2il0floatpacket.11:
	.long	0x00000002,0x00000003,0x00000010,0x00000011,0x00000006,0x00000007,0x00000014,0x00000015,0x00000018,0x00000019,0x0000000a,0x0000000b,0x0000001c,0x0000001d,0x0000000e,0x0000000f
	.type	.L_2il0floatpacket.11,@object
	.size	.L_2il0floatpacket.11,64
	.align 64
.L_2il0floatpacket.12:
	.long	0x00000000,0x00000001,0x00000002,0x00000003,0x00000004,0x00000005,0x00000006,0x00000007,0x00000009,0x00000008,0x0000000b,0x0000000a,0x0000000d,0x0000000c,0x0000000f,0x0000000e
	.type	.L_2il0floatpacket.12,@object
	.size	.L_2il0floatpacket.12,64
	.align 64
.L_2il0floatpacket.13:
	.long	0x00000000,0x00000001,0x00000002,0x00000003,0x00000008,0x00000009,0x00000006,0x00000007,0x00000010,0x00000011,0x00000012,0x00000013,0x00000018,0x00000019,0x00000016,0x00000017
	.type	.L_2il0floatpacket.13,@object
	.size	.L_2il0floatpacket.13,64
	.align 64
.L_2il0floatpacket.14:
	.long	0x00000000,0x00000000,0x00000000,0x00000000,0x00000008,0x00000008,0x00000008,0x00000008,0x00000010,0x00000010,0x00000010,0x00000010,0x00000018,0x00000018,0x00000018,0x00000018
	.type	.L_2il0floatpacket.14,@object
	.size	.L_2il0floatpacket.14,64
	.align 64
.L_2il0floatpacket.15:
	.long	0x00000004,0x00000004,0x00000004,0x00000004,0x0000000c,0x0000000c,0x0000000c,0x0000000c,0x00000014,0x00000014,0x00000014,0x00000014,0x0000001c,0x0000001c,0x0000001c,0x0000001c
	.type	.L_2il0floatpacket.15,@object
	.size	.L_2il0floatpacket.15,64
	.align 64
.L_2il0floatpacket.16:
	.long	0x00000002,0x00000002,0x00000002,0x00000002,0x0000000e,0x0000000e,0x0000000e,0x0000000e,0x00000012,0x00000012,0x00000012,0x00000012,0x0000001e,0x0000001e,0x0000001e,0x0000001e
	.type	.L_2il0floatpacket.16,@object
	.size	.L_2il0floatpacket.16,64
	.align 64
.L_2il0floatpacket.17:
	.long	0x00000006,0x00000006,0x00000006,0x00000006,0x0000000a,0x0000000a,0x0000000a,0x0000000a,0x00000016,0x00000016,0x00000016,0x00000016,0x0000001a,0x0000001a,0x0000001a,0x0000001a
	.type	.L_2il0floatpacket.17,@object
	.size	.L_2il0floatpacket.17,64
	.align 64
.L_2il0floatpacket.18:
	.long	0xbf800000,0x3f800000,0xbf800000,0x3f800000,0x3f800000,0xbf800000,0x3f800000,0xbf800000,0xbf800000,0x3f800000,0xbf800000,0x3f800000,0x3f800000,0xbf800000,0x3f800000,0xbf800000
	.type	.L_2il0floatpacket.18,@object
	.size	.L_2il0floatpacket.18,64
	.align 64
.L_2il0floatpacket.19:
	.long	0x00000001,0x00000001,0x00000001,0x00000001,0x00000009,0x00000009,0x00000009,0x00000009,0x00000011,0x00000011,0x00000011,0x00000011,0x00000019,0x00000019,0x00000019,0x00000019
	.type	.L_2il0floatpacket.19,@object
	.size	.L_2il0floatpacket.19,64
	.align 64
.L_2il0floatpacket.20:
	.long	0x00000005,0x00000005,0x00000005,0x00000005,0x0000000d,0x0000000d,0x0000000d,0x0000000d,0x00000015,0x00000015,0x00000015,0x00000015,0x0000001d,0x0000001d,0x0000001d,0x0000001d
	.type	.L_2il0floatpacket.20,@object
	.size	.L_2il0floatpacket.20,64
	.align 64
.L_2il0floatpacket.21:
	.long	0x00000003,0x00000003,0x00000003,0x00000003,0x0000000f,0x0000000f,0x0000000f,0x0000000f,0x00000013,0x00000013,0x00000013,0x00000013,0x0000001f,0x0000001f,0x0000001f,0x0000001f
	.type	.L_2il0floatpacket.21,@object
	.size	.L_2il0floatpacket.21,64
	.align 64
.L_2il0floatpacket.22:
	.long	0x00000007,0x00000007,0x00000007,0x00000007,0x0000000b,0x0000000b,0x0000000b,0x0000000b,0x00000017,0x00000017,0x00000017,0x00000017,0x0000001b,0x0000001b,0x0000001b,0x0000001b
	.type	.L_2il0floatpacket.22,@object
	.size	.L_2il0floatpacket.22,64
	.align 64
.L_2il0floatpacket.23:
	.long	0x00000004,0x00000005,0x0000000c,0x0000000d,0x0000000a,0x0000000b,0x0000000e,0x0000000f,0x00000014,0x00000015,0x0000001c,0x0000001d,0x0000001a,0x0000001b,0x0000001e,0x0000001f
	.type	.L_2il0floatpacket.23,@object
	.size	.L_2il0floatpacket.23,64
	.align 64
.L_2il0floatpacket.24:
	.long	0x00000000,0x00000000,0x00000000,0x00000000,0x0000000a,0x0000000a,0x0000000a,0x0000000a,0x00000010,0x00000010,0x00000010,0x00000010,0x0000001a,0x0000001a,0x0000001a,0x0000001a
	.type	.L_2il0floatpacket.24,@object
	.size	.L_2il0floatpacket.24,64
	.align 64
.L_2il0floatpacket.25:
	.long	0x00000002,0x00000002,0x00000002,0x00000002,0x00000008,0x00000008,0x00000008,0x00000008,0x00000012,0x00000012,0x00000012,0x00000012,0x00000018,0x00000018,0x00000018,0x00000018
	.type	.L_2il0floatpacket.25,@object
	.size	.L_2il0floatpacket.25,64
	.align 64
.L_2il0floatpacket.26:
	.long	0x00000004,0x00000004,0x00000004,0x00000004,0x0000000e,0x0000000e,0x0000000e,0x0000000e,0x00000014,0x00000014,0x00000014,0x00000014,0x0000001e,0x0000001e,0x0000001e,0x0000001e
	.type	.L_2il0floatpacket.26,@object
	.size	.L_2il0floatpacket.26,64
	.align 64
.L_2il0floatpacket.27:
	.long	0x00000006,0x00000006,0x00000006,0x00000006,0x0000000c,0x0000000c,0x0000000c,0x0000000c,0x00000016,0x00000016,0x00000016,0x00000016,0x0000001c,0x0000001c,0x0000001c,0x0000001c
	.type	.L_2il0floatpacket.27,@object
	.size	.L_2il0floatpacket.27,64
	.align 64
.L_2il0floatpacket.28:
	.long	0x00000001,0x00000001,0x00000001,0x00000001,0x0000000b,0x0000000b,0x0000000b,0x0000000b,0x00000011,0x00000011,0x00000011,0x00000011,0x0000001b,0x0000001b,0x0000001b,0x0000001b
	.type	.L_2il0floatpacket.28,@object
	.size	.L_2il0floatpacket.28,64
	.align 64
.L_2il0floatpacket.29:
	.long	0x00000003,0x00000003,0x00000003,0x00000003,0x00000009,0x00000009,0x00000009,0x00000009,0x00000013,0x00000013,0x00000013,0x00000013,0x00000019,0x00000019,0x00000019,0x00000019
	.type	.L_2il0floatpacket.29,@object
	.size	.L_2il0floatpacket.29,64
	.align 64
.L_2il0floatpacket.30:
	.long	0x00000005,0x00000005,0x00000005,0x00000005,0x0000000f,0x0000000f,0x0000000f,0x0000000f,0x00000015,0x00000015,0x00000015,0x00000015,0x0000001f,0x0000001f,0x0000001f,0x0000001f
	.type	.L_2il0floatpacket.30,@object
	.size	.L_2il0floatpacket.30,64
	.align 64
.L_2il0floatpacket.31:
	.long	0x00000007,0x00000007,0x00000007,0x00000007,0x0000000d,0x0000000d,0x0000000d,0x0000000d,0x00000017,0x00000017,0x00000017,0x00000017,0x0000001d,0x0000001d,0x0000001d,0x0000001d
	.type	.L_2il0floatpacket.31,@object
	.size	.L_2il0floatpacket.31,64
	.align 64
.L_2il0floatpacket.32:
	.long	0x00000000,0x00000001,0x00000010,0x00000011,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000
	.type	.L_2il0floatpacket.32,@object
	.size	.L_2il0floatpacket.32,64
	.align 64
.L_2il0floatpacket.33:
	.long	0x00000000,0x00000000,0x00000000,0x00000000,0x00000002,0x00000002,0x00000002,0x00000002,0x00000010,0x00000010,0x00000010,0x00000010,0x00000012,0x00000012,0x00000012,0x00000012
	.type	.L_2il0floatpacket.33,@object
	.size	.L_2il0floatpacket.33,64
	.align 64
.L_2il0floatpacket.34:
	.long	0x00000001,0x00000001,0x00000001,0x00000001,0x00000003,0x00000003,0x00000003,0x00000003,0x00000011,0x00000011,0x00000011,0x00000011,0x00000013,0x00000013,0x00000013,0x00000013
	.type	.L_2il0floatpacket.34,@object
	.size	.L_2il0floatpacket.34,64
	.align 64
.L_2il0floatpacket.39:
	.long	0x00000012,0x00000013,0x00000000,0x00000001,0x00000016,0x00000017,0x00000004,0x00000005,0x00000008,0x00000009,0x0000001a,0x0000001b,0x0000000c,0x0000000d,0x0000001e,0x0000001f
	.type	.L_2il0floatpacket.39,@object
	.size	.L_2il0floatpacket.39,64
	.align 64
.L_2il0floatpacket.40:
	.long	0x00000010,0x00000011,0x00000002,0x00000003,0x00000014,0x00000015,0x00000006,0x00000007,0x0000000a,0x0000000b,0x00000018,0x00000019,0x0000000e,0x0000000f,0x0000001c,0x0000001d
	.type	.L_2il0floatpacket.40,@object
	.size	.L_2il0floatpacket.40,64
	.align 64
.L_2il0floatpacket.48:
	.long	0x00000000,0x00000001,0x00000010,0x00000011,0x00000004,0x00000005,0x00000014,0x00000015,0x0000000a,0x0000000b,0x0000001a,0x0000001b,0x0000000e,0x0000000f,0x0000001e,0x0000001f
	.type	.L_2il0floatpacket.48,@object
	.size	.L_2il0floatpacket.48,64
	.align 64
.L_2il0floatpacket.49:
	.long	0x00000000,0x00000001,0x00000012,0x00000013,0x00000004,0x00000005,0x00000016,0x00000017,0x0000000a,0x0000000b,0x00000018,0x00000019,0x0000000e,0x0000000f,0x0000001c,0x0000001d
	.type	.L_2il0floatpacket.49,@object
	.size	.L_2il0floatpacket.49,64
	.align 64
.L_2il0floatpacket.50:
	.long	0x00000002,0x00000003,0x00000012,0x00000013,0x00000006,0x00000007,0x00000016,0x00000017,0x00000008,0x00000009,0x00000018,0x00000019,0x0000000c,0x0000000d,0x0000001c,0x0000001d
	.type	.L_2il0floatpacket.50,@object
	.size	.L_2il0floatpacket.50,64
	.align 64
.L_2il0floatpacket.51:
	.long	0x00000006,0x00000007,0x00000004,0x00000005,0x00000006,0x00000007,0x00000004,0x00000005,0x00000005,0x00000004,0x00000007,0x00000006,0x00000005,0x00000004,0x00000007,0x00000006
	.type	.L_2il0floatpacket.51,@object
	.size	.L_2il0floatpacket.51,64
	.align 64
.L_2il0floatpacket.52:
	.long	0x00000000,0x00000001,0x00000002,0x00000003,0x00000000,0x00000001,0x00000002,0x00000003,0x00000000,0x00000001,0x00000002,0x00000003,0x00000000,0x00000001,0x00000002,0x00000003
	.type	.L_2il0floatpacket.52,@object
	.size	.L_2il0floatpacket.52,64
	.align 64
.L_2il0floatpacket.53:
	.long	0x00000002,0x00000003,0x00000010,0x00000011,0x00000006,0x00000007,0x00000014,0x00000015,0x00000008,0x00000009,0x0000001a,0x0000001b,0x0000000c,0x0000000d,0x0000001e,0x0000001f
	.type	.L_2il0floatpacket.53,@object
	.size	.L_2il0floatpacket.53,64
	.align 32
.L_2il0floatpacket.35:
	.long	0x00000004,0x00000005,0x00000006,0x00000007,0x00000000,0x00000001,0x00000002,0x00000003
	.type	.L_2il0floatpacket.35,@object
	.size	.L_2il0floatpacket.35,32
	.align 32
.L_2il0floatpacket.36:
	.long	0x3f800000,0x3f800000,0x3f800000,0x3f800000,0xbf800000,0xbf800000,0xbf800000,0xbf800000
	.type	.L_2il0floatpacket.36,@object
	.size	.L_2il0floatpacket.36,32
	.align 32
.L_2il0floatpacket.37:
	.long	0x00000000,0x00000001,0x00000002,0x00000003,0x00000003,0x00000002,0x00000001,0x00000000
	.type	.L_2il0floatpacket.37,@object
	.size	.L_2il0floatpacket.37,32
	.align 32
.L_2il0floatpacket.38:
	.long	0x3f800000,0x3f800000,0x3f800000,0x3f800000,0xbf800000,0x3f800000,0xbf800000,0x3f800000
	.type	.L_2il0floatpacket.38,@object
	.size	.L_2il0floatpacket.38,32
	.align 32
.L_2il0floatpacket.41:
	.long	0x00000000,0x00000001,0x00000002,0x00000003,0x00000002,0x00000003,0x00000000,0x00000001
	.type	.L_2il0floatpacket.41,@object
	.size	.L_2il0floatpacket.41,32
	.align 32
.L_2il0floatpacket.42:
	.long	0x3f800000,0x3f800000,0x3f800000,0x3f800000,0x3f800000,0x3f800000,0xbf800000,0xbf800000
	.type	.L_2il0floatpacket.42,@object
	.size	.L_2il0floatpacket.42,32
	.align 32
.L_2il0floatpacket.43:
	.long	0x00000000,0x00000001,0x00000002,0x00000003,0x00000001,0x00000000,0x00000003,0x00000002
	.type	.L_2il0floatpacket.43,@object
	.size	.L_2il0floatpacket.43,32
	.align 32
.L_2il0floatpacket.44:
	.long	0x3f800000,0x3f800000,0x3f800000,0x3f800000,0xbf800000,0x3f800000,0x3f800000,0xbf800000
	.type	.L_2il0floatpacket.44,@object
	.size	.L_2il0floatpacket.44,32
	.align 32
.L_2il0floatpacket.45:
	.long	0x00000007,0x00000006,0x00000005,0x00000004,0x00000007,0x00000006,0x00000005,0x00000004
	.type	.L_2il0floatpacket.45,@object
	.size	.L_2il0floatpacket.45,32
	.align 32
.L_2il0floatpacket.46:
	.long	0x00000000,0x00000001,0x00000002,0x00000003,0x00000000,0x00000001,0x00000002,0x00000003
	.type	.L_2il0floatpacket.46,@object
	.size	.L_2il0floatpacket.46,32
	.align 32
.L_2il0floatpacket.47:
	.long	0xbf800000,0x3f800000,0xbf800000,0x3f800000,0x3f800000,0xbf800000,0x3f800000,0xbf800000
	.type	.L_2il0floatpacket.47,@object
	.size	.L_2il0floatpacket.47,32
	.data
	.section .note.GNU-stack, ""
// -- Begin DWARF2 SEGMENT .eh_frame
	.section .eh_frame,"a",@progbits
.eh_frame_seg:
	.align 8
# End
