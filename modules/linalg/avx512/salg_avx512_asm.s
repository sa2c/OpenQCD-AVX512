# mark_description "Intel(R) C Intel(R) 64 Compiler for applications running on Intel(R) 64, Version 17.0.4.196 Build 20170411";
# mark_description "-I../../../include -I.. -I/cineca/prod/opt/compilers/intel/pe-xe-2017/binary/impi/2017.3.196/intel64/include";
# mark_description " -isystem /cineca/prod/opt/compilers/intel/pe-xe-2018/binary/impi/2018.1.163/include64/ -std=c89 -xCORE-AVX5";
# mark_description "12 -mtune=skylake -DAVX512 -O3 -Ddirac_counters -pedantic -fstrict-aliasing -Wno-long-long -Wstrict-prototyp";
# mark_description "es -S";
	.file "salg_avx512.c"
	.text
..TXTST0:
# -- Begin  mulc_spinor_add_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl mulc_spinor_add_avx512
# --- mulc_spinor_add_avx512(int, spinor *, const spinor *, complex)
mulc_spinor_add_avx512:
# parameter 1: %edi
# parameter 2: %rsi
# parameter 3: %rdx
# parameter 4: %xmm0
..B1.1:                         # Preds ..B1.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_mulc_spinor_add_avx512.1:
..L2:
                                                          #9.1
        vshufps   $1, %xmm0, %xmm0, %xmm1                       #9.1
        vbroadcastss %xmm1, %xmm3                               #19.8
        vbroadcastss %xmm0, %xmm2                               #18.8
        movslq    %edi, %rdi                                    #9.1
        vshuff32x4 $0, %zmm3, %zmm3, %zmm4                      #21.8
        vmulps    .L_2il0floatpacket.3(%rip), %zmm4, %zmm0      #23.8
        lea       (%rdi,%rdi,2), %rax                           #16.8
        shlq      $5, %rax                                      #16.8
        addq      %rsi, %rax                                    #16.8
        vshuff32x4 $0, %zmm2, %zmm2, %zmm1                      #20.8
        cmpq      %rax, %rsi                                    #25.14
        jae       ..B1.5        # Prob 10%                      #25.14
                                # LOE rax rdx rbx rbp rsi r12 r13 r14 r15 zmm0 zmm1
..B1.3:                         # Preds ..B1.1 ..B1.3
                                # Execution count [5.00e+00]
        vmovups   (%rdx), %zmm3                                 #26.30
        vmulps    %zmm3, %zmm0, %zmm2                           #27.10
        vpermilps $177, %zmm2, %zmm4                            #28.10
        vfmadd231ps %zmm1, %zmm3, %zmm4                         #29.10
        vaddps    (%rsi), %zmm4, %zmm5                          #31.10
        vmovups   %zmm5, (%rsi)                                 #32.26
        vmovups   64(%rdx), %zmm7                               #34.30
        vmulps    %zmm7, %zmm0, %zmm6                           #35.10
        vpermilps $177, %zmm6, %zmm8                            #36.10
        vfmadd231ps %zmm1, %zmm7, %zmm8                         #37.10
        vaddps    64(%rsi), %zmm8, %zmm9                        #39.10
        vmovups   %zmm9, 64(%rsi)                               #40.26
        vmovups   128(%rdx), %zmm11                             #42.30
        addq      $192, %rdx                                    #50.5
        vmulps    %zmm11, %zmm0, %zmm10                         #43.10
        vpermilps $177, %zmm10, %zmm12                          #44.10
        vfmadd231ps %zmm1, %zmm11, %zmm12                       #45.10
        vaddps    128(%rsi), %zmm12, %zmm13                     #47.10
        vmovups   %zmm13, 128(%rsi)                             #48.26
        addq      $192, %rsi                                    #25.18
        cmpq      %rax, %rsi                                    #25.14
        jb        ..B1.3        # Prob 82%                      #25.14
                                # LOE rax rdx rbx rbp rsi r12 r13 r14 r15 zmm0 zmm1
..B1.5:                         # Preds ..B1.3 ..B1.1
                                # Execution count [1.00e+00]
        vzeroupper                                              #52.1
        ret                                                     #52.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	mulc_spinor_add_avx512,@function
	.size	mulc_spinor_add_avx512,.-mulc_spinor_add_avx512
	.data
# -- End  mulc_spinor_add_avx512
	.text
# -- Begin  spinor_prod_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl spinor_prod_avx512
# --- spinor_prod_avx512(int, const spinor *, const spinor *)
spinor_prod_avx512:
# parameter 1: %edi
# parameter 2: %rsi
# parameter 3: %rdx
..B2.1:                         # Preds ..B2.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_spinor_prod_avx512.4:
..L5:
                                                          #65.1
        subq      $136, %rsp                                    #65.1
	.cfi_def_cfa_offset 144
        movslq    %edi, %rdi                                    #65.1
        vxorpd    %xmm0, %xmm0, %xmm0                           #72.3
        vmovapd   %xmm0, %xmm1                                  #73.3
        lea       (%rdi,%rdi,2), %rax                           #74.8
        shlq      $5, %rax                                      #74.8
        addq      %rsi, %rax                                    #74.8
        cmpq      %rax, %rsi                                    #77.14
        jae       ..B2.9        # Prob 10%                      #77.14
                                # LOE rax rdx rbx rbp rsi r12 r13 r14 r15 xmm0 xmm1
..B2.2:                         # Preds ..B2.1
                                # Execution count [9.00e-01]
        vmovups   .L_2il0floatpacket.3(%rip), %zmm5             #102.29
        vpxord    %zmm4, %zmm4, %zmm4                           #83.10
        vmovaps   %zmm4, %zmm2                                  #83.10
                                # LOE rax rdx rbx rbp rsi r12 r13 r14 r15 xmm0 xmm1 zmm2 zmm4 zmm5
..B2.3:                         # Preds ..B2.7 ..B2.2
                                # Execution count [5.00e+00]
        vmovaps   %zmm4, %zmm3                                  #84.10
        lea       768(%rsi), %rcx                               #78.11
        cmpq      %rax, %rcx                                    #79.5
        vmovaps   %zmm2, %zmm6                                  #83.10
        cmovae    %rax, %rcx                                    #79.5
        vmovaps   %zmm3, %zmm2                                  #84.10
        cmpq      %rcx, %rsi                                    #86.16
        jae       ..B2.7        # Prob 10%                      #86.16
                                # LOE rax rdx rcx rbx rbp rsi r12 r13 r14 r15 xmm0 xmm1 zmm2 zmm3 zmm4 zmm5 zmm6
..B2.5:                         # Preds ..B2.3 ..B2.5
                                # Execution count [2.50e+01]
        vmovups   (%rdx), %zmm7                                 #90.32
        vmovups   64(%rdx), %zmm8                               #91.32
        vmovups   128(%rdx), %zmm9                              #92.32
        vmovups   (%rsi), %zmm14                                #87.32
        vmovups   64(%rsi), %zmm16                              #88.32
        vmovups   128(%rsi), %zmm18                             #89.32
        vfmadd231ps %zmm14, %zmm7, %zmm6                        #94.12
        vpermilps $177, %zmm7, %zmm10                           #98.12
        addq      $192, %rsi                                    #86.21
        vmulps    %zmm5, %zmm10, %zmm13                         #103.12
        addq      $192, %rdx                                    #111.7
        vfmadd231ps %zmm16, %zmm8, %zmm6                        #95.12
        vfmadd231ps %zmm14, %zmm13, %zmm3                       #107.12
        vfmadd231ps %zmm18, %zmm9, %zmm6                        #96.12
        vpermilps $177, %zmm8, %zmm11                           #99.12
        vmulps    %zmm11, %zmm5, %zmm15                         #104.12
        vpermilps $177, %zmm9, %zmm12                           #100.12
        vmulps    %zmm12, %zmm5, %zmm17                         #105.12
        vfmadd231ps %zmm16, %zmm15, %zmm3                       #108.12
        vfmadd231ps %zmm18, %zmm17, %zmm3                       #109.12
        cmpq      %rcx, %rsi                                    #86.16
        jb        ..B2.5        # Prob 82%                      #86.16
                                # LOE rax rdx rcx rbx rbp rsi r12 r13 r14 r15 xmm0 xmm1 zmm2 zmm3 zmm4 zmm5 zmm6
..B2.7:                         # Preds ..B2.5 ..B2.3
                                # Execution count [5.00e+00]
        vmovups   %zmm6, (%rsp)                                 #114.19
        vmovups   %zmm3, 64(%rsp)                               #115.19
        vmovss    (%rsp), %xmm6                                 #114.19
        vmovss    8(%rsp), %xmm7                                #114.19
        vmovss    16(%rsp), %xmm10                              #114.19
        vmovss    24(%rsp), %xmm11                              #114.19
        vmovss    32(%rsp), %xmm16                              #114.19
        vmovss    40(%rsp), %xmm17                              #114.19
        vmovss    64(%rsp), %xmm3                               #115.19
        vmovss    48(%rsp), %xmm20                              #114.19
        vmovss    56(%rsp), %xmm21                              #114.19
        vmovss    72(%rsp), %xmm29                              #115.19
        vaddss    4(%rsp), %xmm6, %xmm8                         #114.5
        vaddss    12(%rsp), %xmm7, %xmm9                        #114.5
        vaddss    20(%rsp), %xmm10, %xmm12                      #114.5
        vaddss    28(%rsp), %xmm11, %xmm13                      #114.5
        vaddss    36(%rsp), %xmm16, %xmm18                      #114.5
        vaddss    44(%rsp), %xmm17, %xmm19                      #114.5
        vaddss    %xmm9, %xmm8, %xmm14                          #114.5
        vaddss    68(%rsp), %xmm3, %xmm30                       #115.5
        vaddss    %xmm13, %xmm12, %xmm15                        #114.5
        vaddss    52(%rsp), %xmm20, %xmm22                      #114.5
        vaddss    60(%rsp), %xmm21, %xmm23                      #114.5
        vaddss    76(%rsp), %xmm29, %xmm31                      #115.5
        vaddss    %xmm19, %xmm18, %xmm24                        #114.5
        vaddss    %xmm15, %xmm14, %xmm26                        #114.5
        vaddss    %xmm23, %xmm22, %xmm25                        #114.5
        vaddss    %xmm31, %xmm30, %xmm9                         #115.5
        vaddss    %xmm25, %xmm24, %xmm27                        #114.5
        vmovss    80(%rsp), %xmm3                               #115.19
        vaddss    %xmm27, %xmm26, %xmm28                        #114.5
        vaddss    84(%rsp), %xmm3, %xmm7                        #115.5
        vcvtss2sd %xmm28, %xmm28, %xmm28                        #114.19
        vmovss    88(%rsp), %xmm6                               #115.19
        vaddsd    %xmm28, %xmm0, %xmm0                          #114.5
        vaddss    92(%rsp), %xmm6, %xmm8                        #115.5
        vmovss    96(%rsp), %xmm11                              #115.19
        vaddss    %xmm8, %xmm7, %xmm10                          #115.5
        vaddss    100(%rsp), %xmm11, %xmm13                     #115.5
        vaddss    %xmm10, %xmm9, %xmm21                         #115.5
        vmovss    104(%rsp), %xmm12                             #115.19
        vmovss    112(%rsp), %xmm15                             #115.19
        vmovss    120(%rsp), %xmm16                             #115.19
        vaddss    108(%rsp), %xmm12, %xmm14                     #115.5
        vaddss    116(%rsp), %xmm15, %xmm17                     #115.5
        vaddss    124(%rsp), %xmm16, %xmm18                     #115.5
        vaddss    %xmm14, %xmm13, %xmm19                        #115.5
        vaddss    %xmm18, %xmm17, %xmm20                        #115.5
        vaddss    %xmm20, %xmm19, %xmm22                        #115.5
        vaddss    %xmm22, %xmm21, %xmm23                        #115.5
        vcvtss2sd %xmm23, %xmm23, %xmm23                        #115.19
        vaddsd    %xmm23, %xmm1, %xmm1                          #115.5
        cmpq      %rax, %rsi                                    #77.14
        jb        ..B2.3        # Prob 82%                      #77.14
                                # LOE rax rdx rbx rbp rsi r12 r13 r14 r15 xmm0 xmm1 zmm2 zmm4 zmm5
..B2.9:                         # Preds ..B2.7 ..B2.1
                                # Execution count [1.00e+00]
        vzeroupper                                              #122.10
        addq      $136, %rsp                                    #122.10
	.cfi_def_cfa_offset 8
        ret                                                     #122.10
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	spinor_prod_avx512,@function
	.size	spinor_prod_avx512,.-spinor_prod_avx512
	.data
# -- End  spinor_prod_avx512
	.section .rodata, "a"
	.align 64
	.align 64
.L_2il0floatpacket.3:
	.long	0x3f800000,0xbf800000,0x3f800000,0xbf800000,0x3f800000,0xbf800000,0x3f800000,0xbf800000,0x3f800000,0xbf800000,0x3f800000,0xbf800000,0x3f800000,0xbf800000,0x3f800000,0xbf800000
	.type	.L_2il0floatpacket.3,@object
	.size	.L_2il0floatpacket.3,64
	.data
	.section .note.GNU-stack, ""
// -- Begin DWARF2 SEGMENT .eh_frame
	.section .eh_frame,"a",@progbits
.eh_frame_seg:
	.align 8
# End
