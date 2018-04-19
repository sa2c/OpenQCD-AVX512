# mark_description "Intel(R) C Intel(R) 64 Compiler for applications running on Intel(R) 64, Version 17.0.4.196 Build 20170411";
# mark_description "-I../../../include -I.. -I/cineca/prod/opt/compilers/intel/pe-xe-2017/binary/impi/2017.3.196/intel64/include";
# mark_description " -isystem /cineca/prod/opt/compilers/intel/pe-xe-2018/binary/impi/2018.1.163/include64/ -std=c89 -xCORE-AVX5";
# mark_description "12 -mtune=skylake -DAVX512 -O3 -Ddirac_counters -pedantic -fstrict-aliasing -Wno-long-long -Wstrict-prototyp";
# mark_description "es -S";
	.file "salg_dble_avx512.c"
	.text
..TXTST0:
# -- Begin  spinor_prod_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl spinor_prod_dble_avx512
# --- spinor_prod_dble_avx512(const spinor_dble *, const spinor_dble *, const spinor_dble *)
spinor_prod_dble_avx512:
# parameter 1: %rdi
# parameter 2: %rsi
# parameter 3: %rdx
..B1.1:                         # Preds ..B1.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_spinor_prod_dble_avx512.1:
..L2:
                                                          #16.1
        subq      $136, %rsp                                    #16.1
	.cfi_def_cfa_offset 144
        vpxord    %zmm2, %zmm2, %zmm2                           #21.8
        vmovaps   %zmm2, %zmm1                                  #22.8
        cmpq      %rsi, %rdi                                    #23.14
        jae       ..B1.5        # Prob 10%                      #23.14
                                # LOE rdx rbx rbp rsi rdi r12 r13 r14 r15 zmm1 zmm2
..B1.2:                         # Preds ..B1.1
                                # Execution count [9.00e-01]
        vmovups   .L_2il0floatpacket.3(%rip), %zmm0             #36.27
                                # LOE rdx rbx rbp rsi rdi r12 r13 r14 r15 zmm0 zmm1 zmm2
..B1.3:                         # Preds ..B1.3 ..B1.2
                                # Execution count [5.00e+00]
        vmovups   (%rdx), %zmm3                                 #27.30
        vmovups   64(%rdx), %zmm4                               #28.30
        vmovups   128(%rdx), %zmm5                              #29.30
        vmovups   (%rdi), %zmm10                                #24.30
        vmovups   64(%rdi), %zmm12                              #25.30
        vmovups   128(%rdi), %zmm14                             #26.30
        vfmadd231pd %zmm10, %zmm3, %zmm2                        #30.10
        vpermilpd $85, %zmm3, %zmm6                             #33.10
        addq      $192, %rdi                                    #23.19
        vmulpd    %zmm0, %zmm6, %zmm9                           #37.10
        addq      $192, %rdx                                    #43.5
        vfmadd231pd %zmm12, %zmm4, %zmm2                        #31.10
        vfmadd231pd %zmm10, %zmm9, %zmm1                        #40.10
        vfmadd231pd %zmm14, %zmm5, %zmm2                        #32.10
        vpermilpd $85, %zmm4, %zmm7                             #34.10
        vmulpd    %zmm7, %zmm0, %zmm11                          #38.10
        vpermilpd $85, %zmm5, %zmm8                             #35.10
        vmulpd    %zmm8, %zmm0, %zmm13                          #39.10
        vfmadd231pd %zmm12, %zmm11, %zmm1                       #41.10
        vfmadd231pd %zmm14, %zmm13, %zmm1                       #42.10
        cmpq      %rsi, %rdi                                    #23.14
        jb        ..B1.3        # Prob 82%                      #23.14
                                # LOE rdx rbx rbp rsi rdi r12 r13 r14 r15 zmm0 zmm1 zmm2
..B1.5:                         # Preds ..B1.3 ..B1.1
                                # Execution count [1.00e+00]
        vmovups   %zmm2, (%rsp)                                 #45.10
        vmovups   %zmm1, 64(%rsp)                               #46.10
        vmovsd    (%rsp), %xmm2                                 #45.10
        vmovsd    16(%rsp), %xmm3                               #45.10
        vmovsd    32(%rsp), %xmm6                               #45.10
        vmovsd    48(%rsp), %xmm7                               #45.10
        vmovsd    64(%rsp), %xmm1                               #46.10
        vmovsd    80(%rsp), %xmm12                              #46.10
        vmovsd    96(%rsp), %xmm15                              #46.10
        vmovsd    112(%rsp), %xmm16                             #46.10
        vaddsd    8(%rsp), %xmm2, %xmm4                         #45.3
        vaddsd    24(%rsp), %xmm3, %xmm5                        #45.3
        vaddsd    40(%rsp), %xmm6, %xmm8                        #45.3
        vaddsd    56(%rsp), %xmm7, %xmm9                        #45.3
        vaddsd    72(%rsp), %xmm1, %xmm13                       #46.3
        vaddsd    88(%rsp), %xmm12, %xmm14                      #46.3
        vaddsd    104(%rsp), %xmm15, %xmm17                     #46.3
        vaddsd    120(%rsp), %xmm16, %xmm18                     #46.3
        vaddsd    %xmm5, %xmm4, %xmm10                          #45.3
        vaddsd    %xmm9, %xmm8, %xmm11                          #45.3
        vaddsd    %xmm14, %xmm13, %xmm19                        #46.3
        vaddsd    %xmm18, %xmm17, %xmm20                        #46.3
        vaddsd    %xmm11, %xmm10, %xmm0                         #45.3
        vaddsd    %xmm20, %xmm19, %xmm1                         #46.3
        vzeroupper                                              #47.10
        addq      $136, %rsp                                    #47.10
	.cfi_def_cfa_offset 8
        ret                                                     #47.10
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	spinor_prod_dble_avx512,@function
	.size	spinor_prod_dble_avx512,.-spinor_prod_dble_avx512
	.data
# -- End  spinor_prod_dble_avx512
	.text
# -- Begin  spinor_prod_re_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl spinor_prod_re_dble_avx512
# --- spinor_prod_re_dble_avx512(const spinor_dble *, const spinor_dble *, const spinor_dble *)
spinor_prod_re_dble_avx512:
# parameter 1: %rdi
# parameter 2: %rsi
# parameter 3: %rdx
..B2.1:                         # Preds ..B2.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_spinor_prod_re_dble_avx512.6:
..L7:
                                                          #51.1
        vpxord    %zmm0, %zmm0, %zmm0                           #54.8
        cmpq      %rsi, %rdi                                    #56.14
        jae       ..B2.5        # Prob 10%                      #56.14
                                # LOE rdx rbx rbp rsi rdi r12 r13 r14 r15 zmm0
..B2.3:                         # Preds ..B2.1 ..B2.3
                                # Execution count [5.00e+00]
        vmovups   (%rdx), %zmm1                                 #60.30
        vmovups   64(%rdx), %zmm2                               #61.30
        vmovups   128(%rdx), %zmm3                              #62.30
        vfmadd231pd (%rdi), %zmm1, %zmm0                        #63.10
        vfmadd231pd 64(%rdi), %zmm2, %zmm0                      #64.10
        addq      $192, %rdx                                    #69.5
        vfmadd231pd 128(%rdi), %zmm3, %zmm0                     #65.10
        addq      $192, %rdi                                    #56.19
        cmpq      %rsi, %rdi                                    #56.14
        jb        ..B2.3        # Prob 82%                      #56.14
                                # LOE rdx rbx rbp rsi rdi r12 r13 r14 r15 zmm0
..B2.5:                         # Preds ..B2.3 ..B2.1
                                # Execution count [1.00e+00]
        vmovups   %zmm0, -72(%rsp)                              #71.7
        vmovsd    -72(%rsp), %xmm0                              #71.7
        vmovsd    -56(%rsp), %xmm1                              #71.7
        vmovsd    -40(%rsp), %xmm4                              #71.7
        vmovsd    -24(%rsp), %xmm5                              #71.7
        vaddsd    -64(%rsp), %xmm0, %xmm2                       #72.10
        vaddsd    -48(%rsp), %xmm1, %xmm3                       #72.10
        vaddsd    -32(%rsp), %xmm4, %xmm6                       #72.10
        vaddsd    -16(%rsp), %xmm5, %xmm7                       #72.10
        vaddsd    %xmm3, %xmm2, %xmm8                           #72.10
        vaddsd    %xmm7, %xmm6, %xmm9                           #72.10
        vaddsd    %xmm9, %xmm8, %xmm0                           #72.10
        vcvtsd2ss %xmm0, %xmm0, %xmm0                           #71.7
        vcvtss2sd %xmm0, %xmm0, %xmm0                           #72.10
        vzeroupper                                              #72.10
        ret                                                     #72.10
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	spinor_prod_re_dble_avx512,@function
	.size	spinor_prod_re_dble_avx512,.-spinor_prod_re_dble_avx512
	.data
# -- End  spinor_prod_re_dble_avx512
	.text
# -- Begin  spinor_prod5_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl spinor_prod5_dble_avx512
# --- spinor_prod5_dble_avx512(const spinor_dble *, const spinor_dble *, const spinor_dble *)
spinor_prod5_dble_avx512:
# parameter 1: %rdi
# parameter 2: %rsi
# parameter 3: %rdx
..B3.1:                         # Preds ..B3.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_spinor_prod5_dble_avx512.9:
..L10:
                                                         #76.1
        subq      $136, %rsp                                    #76.1
	.cfi_def_cfa_offset 144
        vpxord    %zmm3, %zmm3, %zmm3                           #80.8
        vmovaps   %zmm3, %zmm2                                  #81.8
        cmpq      %rsi, %rdi                                    #82.14
        jae       ..B3.5        # Prob 10%                      #82.14
                                # LOE rdx rbx rbp rsi rdi r12 r13 r14 r15 zmm2 zmm3
..B3.2:                         # Preds ..B3.1
                                # Execution count [9.00e-01]
        vmovups   .L_2il0floatpacket.4(%rip), %zmm1             #89.27
        vmovups   .L_2il0floatpacket.3(%rip), %zmm0             #97.27
                                # LOE rdx rbx rbp rsi rdi r12 r13 r14 r15 zmm0 zmm1 zmm2 zmm3
..B3.3:                         # Preds ..B3.3 ..B3.2
                                # Execution count [5.00e+00]
        vmovups   (%rdx), %zmm4                                 #86.30
        vmovups   64(%rdx), %zmm5                               #87.30
        vmovups   128(%rdx), %zmm6                              #88.30
        vmovups   (%rdi), %zmm11                                #83.30
        vmovups   128(%rdi), %zmm15                             #85.30
        vmulpd    64(%rdi), %zmm1, %zmm13                       #90.10
        vfmadd231pd %zmm11, %zmm4, %zmm3                        #91.10
        vpermilpd $85, %zmm4, %zmm7                             #94.10
        addq      $192, %rdi                                    #82.19
        vmulpd    %zmm0, %zmm7, %zmm10                          #98.10
        addq      $192, %rdx                                    #104.5
        vfmadd231pd %zmm13, %zmm5, %zmm3                        #92.10
        vfmadd231pd %zmm11, %zmm10, %zmm2                       #101.10
        vfnmadd231pd %zmm15, %zmm6, %zmm3                       #93.10
        vpermilpd $85, %zmm5, %zmm8                             #95.10
        vmulpd    %zmm8, %zmm0, %zmm12                          #99.10
        vpermilpd $85, %zmm6, %zmm9                             #96.10
        vmulpd    %zmm9, %zmm0, %zmm14                          #100.10
        vfmadd213pd %zmm2, %zmm12, %zmm13                       #102.10
        vmovaps   %zmm13, %zmm2                                 #103.10
        vfnmadd231pd %zmm15, %zmm14, %zmm2                      #103.10
        cmpq      %rsi, %rdi                                    #82.14
        jb        ..B3.3        # Prob 82%                      #82.14
                                # LOE rdx rbx rbp rsi rdi r12 r13 r14 r15 zmm0 zmm1 zmm2 zmm3
..B3.5:                         # Preds ..B3.3 ..B3.1
                                # Execution count [1.00e+00]
        vmovups   %zmm3, (%rsp)                                 #106.10
        vmovups   %zmm2, 64(%rsp)                               #107.10
        vmovsd    (%rsp), %xmm1                                 #106.10
        vmovsd    16(%rsp), %xmm3                               #106.10
        vmovsd    32(%rsp), %xmm6                               #106.10
        vmovsd    48(%rsp), %xmm7                               #106.10
        vmovsd    64(%rsp), %xmm2                               #107.10
        vmovsd    80(%rsp), %xmm12                              #107.10
        vmovsd    96(%rsp), %xmm15                              #107.10
        vmovsd    112(%rsp), %xmm16                             #107.10
        vaddsd    8(%rsp), %xmm1, %xmm4                         #106.3
        vaddsd    24(%rsp), %xmm3, %xmm5                        #106.3
        vaddsd    40(%rsp), %xmm6, %xmm8                        #106.3
        vaddsd    56(%rsp), %xmm7, %xmm9                        #106.3
        vaddsd    72(%rsp), %xmm2, %xmm13                       #107.3
        vaddsd    88(%rsp), %xmm12, %xmm14                      #107.3
        vaddsd    104(%rsp), %xmm15, %xmm17                     #107.3
        vaddsd    120(%rsp), %xmm16, %xmm18                     #107.3
        vaddsd    %xmm5, %xmm4, %xmm10                          #106.3
        vaddsd    %xmm9, %xmm8, %xmm11                          #106.3
        vaddsd    %xmm14, %xmm13, %xmm19                        #107.3
        vaddsd    %xmm18, %xmm17, %xmm20                        #107.3
        vaddsd    %xmm11, %xmm10, %xmm0                         #106.3
        vaddsd    %xmm20, %xmm19, %xmm1                         #107.3
        vzeroupper                                              #108.10
        addq      $136, %rsp                                    #108.10
	.cfi_def_cfa_offset 8
        ret                                                     #108.10
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	spinor_prod5_dble_avx512,@function
	.size	spinor_prod5_dble_avx512,.-spinor_prod5_dble_avx512
	.data
# -- End  spinor_prod5_dble_avx512
	.text
# -- Begin  norm_square_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl norm_square_dble_avx512
# --- norm_square_dble_avx512(const spinor_dble *, const spinor_dble *)
norm_square_dble_avx512:
# parameter 1: %rdi
# parameter 2: %rsi
..B4.1:                         # Preds ..B4.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_norm_square_dble_avx512.14:
..L15:
                                                         #112.1
        vpxord    %zmm0, %zmm0, %zmm0                           #114.9
        cmpq      %rsi, %rdi                                    #115.14
        jae       ..B4.5        # Prob 10%                      #115.14
                                # LOE rbx rbp rsi rdi r12 r13 r14 r15 zmm0
..B4.3:                         # Preds ..B4.1 ..B4.3
                                # Execution count [5.00e+00]
        vmovups   (%rdi), %zmm1                                 #116.30
        vmovups   64(%rdi), %zmm2                               #117.30
        vmovups   128(%rdi), %zmm3                              #118.30
        vfmadd231pd %zmm1, %zmm1, %zmm0                         #119.11
        vfmadd231pd %zmm2, %zmm2, %zmm0                         #120.11
        addq      $192, %rdi                                    #115.19
        vfmadd231pd %zmm3, %zmm3, %zmm0                         #121.11
        cmpq      %rsi, %rdi                                    #115.14
        jb        ..B4.3        # Prob 82%                      #115.14
                                # LOE rbx rbp rsi rdi r12 r13 r14 r15 zmm0
..B4.5:                         # Preds ..B4.3 ..B4.1
                                # Execution count [1.00e+00]
        vmovups   %zmm0, -72(%rsp)                              #123.10
        vmovsd    -72(%rsp), %xmm0                              #123.10
        vmovsd    -56(%rsp), %xmm1                              #123.10
        vmovsd    -40(%rsp), %xmm4                              #123.10
        vmovsd    -24(%rsp), %xmm5                              #123.10
        vaddsd    -64(%rsp), %xmm0, %xmm2                       #123.10
        vaddsd    -48(%rsp), %xmm1, %xmm3                       #123.10
        vaddsd    -32(%rsp), %xmm4, %xmm6                       #123.10
        vaddsd    -16(%rsp), %xmm5, %xmm7                       #123.10
        vaddsd    %xmm3, %xmm2, %xmm8                           #123.10
        vaddsd    %xmm7, %xmm6, %xmm9                           #123.10
        vaddsd    %xmm9, %xmm8, %xmm0                           #123.10
        vzeroupper                                              #123.10
        ret                                                     #123.10
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	norm_square_dble_avx512,@function
	.size	norm_square_dble_avx512,.-norm_square_dble_avx512
	.data
# -- End  norm_square_dble_avx512
	.text
# -- Begin  mulc_spinor_add_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl mulc_spinor_add_dble_avx512
# --- mulc_spinor_add_dble_avx512(int, spinor_dble *, const spinor_dble *, complex_dble)
mulc_spinor_add_dble_avx512:
# parameter 1: %edi
# parameter 2: %rsi
# parameter 3: %rdx
# parameter 4: %xmm0 %xmm1
..B5.1:                         # Preds ..B5.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_mulc_spinor_add_dble_avx512.17:
..L18:
                                                         #128.1
        vmovapd   %xmm1, %xmm2                                  #134.8
        movslq    %edi, %rdi                                    #128.1
        vbroadcastsd %xmm0, %zmm1                               #135.8
        vbroadcastsd %xmm2, %zmm0                               #136.8
        lea       (%rdi,%rdi,2), %rax                           #138.8
        shlq      $6, %rax                                      #138.8
        addq      %rsi, %rax                                    #138.8
        cmpq      %rax, %rsi                                    #140.14
        jae       ..B5.5        # Prob 10%                      #140.14
                                # LOE rax rdx rbx rbp rsi r12 r13 r14 r15 zmm0 zmm1
..B5.3:                         # Preds ..B5.1 ..B5.3
                                # Execution count [5.00e+00]
        vmovups   (%rdx), %zmm3                                 #141.30
        vmulpd    %zmm3, %zmm0, %zmm2                           #142.10
        vpermilpd $85, %zmm2, %zmm4                             #143.10
        vfmaddsub231pd %zmm1, %zmm3, %zmm4                      #144.10
        vaddpd    (%rsi), %zmm4, %zmm5                          #146.10
        vmovups   %zmm5, (%rsi)                                 #147.26
        vmovups   64(%rdx), %zmm7                               #149.30
        vmulpd    %zmm7, %zmm0, %zmm6                           #150.10
        vpermilpd $85, %zmm6, %zmm8                             #151.10
        vfmaddsub231pd %zmm1, %zmm7, %zmm8                      #152.10
        vaddpd    64(%rsi), %zmm8, %zmm9                        #154.10
        vmovups   %zmm9, 64(%rsi)                               #155.26
        vmovups   128(%rdx), %zmm11                             #157.30
        addq      $192, %rdx                                    #165.5
        vmulpd    %zmm11, %zmm0, %zmm10                         #158.10
        vpermilpd $85, %zmm10, %zmm12                           #159.10
        vfmaddsub231pd %zmm1, %zmm11, %zmm12                    #160.10
        vaddpd    128(%rsi), %zmm12, %zmm13                     #162.10
        vmovups   %zmm13, 128(%rsi)                             #163.26
        addq      $192, %rsi                                    #140.18
        cmpq      %rax, %rsi                                    #140.14
        jb        ..B5.3        # Prob 82%                      #140.14
                                # LOE rax rdx rbx rbp rsi r12 r13 r14 r15 zmm0 zmm1
..B5.5:                         # Preds ..B5.3 ..B5.1
                                # Execution count [1.00e+00]
        vzeroupper                                              #167.1
        ret                                                     #167.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	mulc_spinor_add_dble_avx512,@function
	.size	mulc_spinor_add_dble_avx512,.-mulc_spinor_add_dble_avx512
	.data
# -- End  mulc_spinor_add_dble_avx512
	.text
# -- Begin  mulr_spinor_add_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl mulr_spinor_add_dble_avx512
# --- mulr_spinor_add_dble_avx512(int, spinor_dble *, const spinor_dble *, double)
mulr_spinor_add_dble_avx512:
# parameter 1: %edi
# parameter 2: %rsi
# parameter 3: %rdx
# parameter 4: %xmm0
..B6.1:                         # Preds ..B6.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_mulr_spinor_add_dble_avx512.20:
..L21:
                                                         #171.1
        movslq    %edi, %rdi                                    #171.1
        vbroadcastsd %xmm0, %zmm0                               #177.8
        lea       (%rdi,%rdi,2), %rax                           #179.8
        shlq      $6, %rax                                      #179.8
        addq      %rsi, %rax                                    #179.8
        cmpq      %rax, %rsi                                    #181.14
        jae       ..B6.5        # Prob 10%                      #181.14
        .align    16,0x90
                                # LOE rax rdx rbx rbp rsi r12 r13 r14 r15 zmm0
..B6.3:                         # Preds ..B6.1 ..B6.3
                                # Execution count [5.00e+00]
        vmovups   (%rdx), %zmm1                                 #182.30
        vfmadd213pd (%rsi), %zmm0, %zmm1                        #185.10
        vmovups   %zmm1, (%rsi)                                 #186.26
        vmovups   64(%rdx), %zmm2                               #188.30
        vfmadd213pd 64(%rsi), %zmm0, %zmm2                      #191.10
        vmovups   %zmm2, 64(%rsi)                               #192.26
        vmovups   128(%rdx), %zmm3                              #194.30
        addq      $192, %rdx                                    #200.5
        vfmadd213pd 128(%rsi), %zmm0, %zmm3                     #197.10
        vmovups   %zmm3, 128(%rsi)                              #198.26
        addq      $192, %rsi                                    #181.18
        cmpq      %rax, %rsi                                    #181.14
        jb        ..B6.3        # Prob 82%                      #181.14
                                # LOE rax rdx rbx rbp rsi r12 r13 r14 r15 zmm0
..B6.5:                         # Preds ..B6.3 ..B6.1
                                # Execution count [1.00e+00]
        vzeroupper                                              #202.1
        ret                                                     #202.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	mulr_spinor_add_dble_avx512,@function
	.size	mulr_spinor_add_dble_avx512,.-mulr_spinor_add_dble_avx512
	.data
# -- End  mulr_spinor_add_dble_avx512
	.text
# -- Begin  combine_spinor_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl combine_spinor_dble_avx512
# --- combine_spinor_dble_avx512(int, spinor_dble *, const spinor_dble *, double, double)
combine_spinor_dble_avx512:
# parameter 1: %edi
# parameter 2: %rsi
# parameter 3: %rdx
# parameter 4: %xmm0
# parameter 5: %xmm1
..B7.1:                         # Preds ..B7.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_combine_spinor_dble_avx512.23:
..L24:
                                                         #206.1
        vmovapd   %xmm1, %xmm2                                  #212.11
        movslq    %edi, %rdi                                    #206.1
        vbroadcastsd %xmm0, %zmm1                               #213.9
        vbroadcastsd %xmm2, %zmm0                               #214.9
        lea       (%rdi,%rdi,2), %rax                           #216.8
        shlq      $6, %rax                                      #216.8
        addq      %rsi, %rax                                    #216.8
        cmpq      %rax, %rsi                                    #218.14
        jae       ..B7.5        # Prob 10%                      #218.14
        .align    16,0x90
                                # LOE rax rdx rbx rbp rsi r12 r13 r14 r15 zmm0 zmm1
..B7.3:                         # Preds ..B7.1 ..B7.3
                                # Execution count [5.00e+00]
        vmulpd    (%rdx), %zmm0, %zmm2                          #220.10
        vfmadd231pd (%rsi), %zmm1, %zmm2                        #222.10
        vmovups   %zmm2, (%rsi)                                 #223.26
        vmulpd    64(%rdx), %zmm0, %zmm3                        #226.10
        vfmadd231pd 64(%rsi), %zmm1, %zmm3                      #228.10
        vmovups   %zmm3, 64(%rsi)                               #229.26
        vmulpd    128(%rdx), %zmm0, %zmm4                       #232.10
        addq      $192, %rdx                                    #237.5
        vfmadd231pd 128(%rsi), %zmm1, %zmm4                     #234.10
        vmovups   %zmm4, 128(%rsi)                              #235.26
        addq      $192, %rsi                                    #218.18
        cmpq      %rax, %rsi                                    #218.14
        jb        ..B7.3        # Prob 82%                      #218.14
                                # LOE rax rdx rbx rbp rsi r12 r13 r14 r15 zmm0 zmm1
..B7.5:                         # Preds ..B7.3 ..B7.1
                                # Execution count [1.00e+00]
        vzeroupper                                              #239.1
        ret                                                     #239.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	combine_spinor_dble_avx512,@function
	.size	combine_spinor_dble_avx512,.-combine_spinor_dble_avx512
	.data
# -- End  combine_spinor_dble_avx512
	.text
# -- Begin  scale_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl scale_dble_avx512
# --- scale_dble_avx512(int, double, spinor_dble *)
scale_dble_avx512:
# parameter 1: %edi
# parameter 2: %xmm0
# parameter 3: %rsi
..B8.1:                         # Preds ..B8.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_scale_dble_avx512.26:
..L27:
                                                         #242.1
        movslq    %edi, %rdi                                    #242.1
        vbroadcastsd %xmm0, %zmm0                               #248.8
        lea       (%rdi,%rdi,2), %rax                           #250.8
        shlq      $6, %rax                                      #250.8
        addq      %rsi, %rax                                    #250.8
        cmpq      %rax, %rsi                                    #252.14
        jae       ..B8.5        # Prob 10%                      #252.14
                                # LOE rax rbx rbp rsi r12 r13 r14 r15 zmm0
..B8.3:                         # Preds ..B8.1 ..B8.3
                                # Execution count [5.00e+00]
        vmulpd    (%rsi), %zmm0, %zmm1                          #254.10
        vmulpd    64(%rsi), %zmm0, %zmm2                        #258.10
        vmulpd    128(%rsi), %zmm0, %zmm3                       #262.10
        vmovups   %zmm1, (%rsi)                                 #255.26
        vmovups   %zmm2, 64(%rsi)                               #259.26
        vmovups   %zmm3, 128(%rsi)                              #263.26
        addq      $192, %rsi                                    #252.18
        cmpq      %rax, %rsi                                    #252.14
        jb        ..B8.3        # Prob 82%                      #252.14
                                # LOE rax rbx rbp rsi r12 r13 r14 r15 zmm0
..B8.5:                         # Preds ..B8.3 ..B8.1
                                # Execution count [1.00e+00]
        vzeroupper                                              #265.1
        ret                                                     #265.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	scale_dble_avx512,@function
	.size	scale_dble_avx512,.-scale_dble_avx512
	.data
# -- End  scale_dble_avx512
	.text
# -- Begin  rotate_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl rotate_dble_avx512
# --- rotate_dble_avx512(int, int, spinor_dble **, spinor_dble *, const complex_dble *)
rotate_dble_avx512:
# parameter 1: %edi
# parameter 2: %esi
# parameter 3: %rdx
# parameter 4: %rcx
# parameter 5: %r8
..B9.1:                         # Preds ..B9.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_rotate_dble_avx512.29:
..L30:
                                                         #268.1
        xorl      %r10d, %r10d                                  #273.8
        movslq    %edi, %r9                                     #268.1
        xorl      %edi, %edi                                    #273.8
        testq     %r9, %r9                                      #273.19
        jle       ..B9.9        # Prob 10%                      #273.19
                                # LOE rdx rcx rbx rbp rdi r8 r9 r10 r12 r13 r14 r15 esi
..B9.2:                         # Preds ..B9.1
                                # Execution count [9.00e-01]
        movslq    %esi, %rsi                                    #268.1
        movq      %r9, %rax                                     #306.7
        vmovups   .L_2il0floatpacket.3(%rip), %zmm0             #286.12
        shlq      $4, %rax                                      #306.7
        lea       (%rsi,%rsi,2), %rsi                           #278.10
        shlq      $6, %rsi                                      #278.10
        movq      %r15, -24(%rsp)                               #306.7[spill]
        movq      %rbx, -16(%rsp)                               #306.7[spill]
	.cfi_offset 3, -24
	.cfi_offset 15, -32
                                # LOE rax rdx rcx rbp rsi rdi r8 r9 r10 r12 r13 r14 zmm0
..B9.3:                         # Preds ..B9.7 ..B9.2
                                # Execution count [5.00e+00]
        vmulpd    8(%r8){1to8}, %zmm0, %zmm6                    #287.10
        movq      %r8, %r11                                     #279.5
        vbroadcastsd (%r8), %zmm1                               #283.10
        movq      (%rdx), %rbx                                  #278.10
        vmovups   (%rsi,%rbx), %zmm2                            #289.30
        vmovups   64(%rsi,%rbx), %zmm5                          #294.30
        vmovups   128(%rsi,%rbx), %zmm8                         #299.30
        vmulpd    %zmm2, %zmm6, %zmm3                           #290.10
        vmulpd    %zmm5, %zmm6, %zmm4                           #295.10
        vmulpd    %zmm8, %zmm6, %zmm7                           #300.10
        vpermilpd $85, %zmm3, %zmm3                             #291.10
        movl      $1, %ebx                                      #304.10
        vfmadd231pd %zmm1, %zmm2, %zmm3                         #292.10
        vpermilpd $85, %zmm4, %zmm2                             #296.10
        vpermilpd $85, %zmm7, %zmm9                             #301.10
        vfmadd231pd %zmm1, %zmm5, %zmm2                         #297.10
        vfmadd213pd %zmm9, %zmm8, %zmm1                         #302.10
        cmpq      $1, %r9                                       #304.21
        jle       ..B9.7        # Prob 10%                      #304.21
                                # LOE rax rdx rcx rbx rbp rsi rdi r8 r9 r10 r11 r12 r13 r14 zmm0 zmm1 zmm2 zmm3
..B9.5:                         # Preds ..B9.3 ..B9.5
                                # Execution count [2.50e+01]
        addq      %rax, %r11                                    #306.7
        movq      (%rdx,%rbx,8), %r15                           #305.12
        incq      %rbx                                          #304.24
        vmulpd    8(%r11){1to8}, %zmm0, %zmm10                  #312.12
        vmovups   (%r15,%rsi), %zmm5                            #314.32
        vmovups   64(%r15,%rsi), %zmm8                          #320.32
        vmovups   128(%r15,%rsi), %zmm12                        #326.32
        vbroadcastsd (%r11), %zmm14                             #310.12
        vmulpd    %zmm5, %zmm10, %zmm4                          #315.12
        vmulpd    %zmm8, %zmm10, %zmm7                          #321.12
        vmulpd    %zmm12, %zmm10, %zmm11                        #327.12
        vpermilpd $85, %zmm4, %zmm6                             #316.12
        vpermilpd $85, %zmm7, %zmm9                             #322.12
        vpermilpd $85, %zmm11, %zmm13                           #328.12
        vfmadd231pd %zmm14, %zmm5, %zmm6                        #317.12
        vfmadd231pd %zmm14, %zmm8, %zmm9                        #323.12
        vfmadd213pd %zmm13, %zmm12, %zmm14                      #329.12
        vaddpd    %zmm3, %zmm6, %zmm3                           #318.12
        vaddpd    %zmm2, %zmm9, %zmm2                           #324.12
        vaddpd    %zmm1, %zmm14, %zmm1                          #330.12
        cmpq      %r9, %rbx                                     #304.21
        jl        ..B9.5        # Prob 82%                      #304.21
                                # LOE rax rdx rcx rbx rbp rsi rdi r8 r9 r10 r11 r12 r13 r14 zmm0 zmm1 zmm2 zmm3
..B9.7:                         # Preds ..B9.5 ..B9.3
                                # Execution count [5.00e+00]
        incq      %r10                                          #273.22
        addq      $16, %r8                                      #273.22
        vmovups   %zmm3, (%rdi,%rcx)                            #333.26
        vmovups   %zmm2, 64(%rdi,%rcx)                          #334.26
        vmovups   %zmm1, 128(%rdi,%rcx)                         #335.26
        addq      $192, %rdi                                    #273.22
        cmpq      %r9, %r10                                     #273.19
        jl        ..B9.3        # Prob 82%                      #273.19
                                # LOE rax rdx rcx rbp rsi rdi r8 r9 r10 r12 r13 r14 zmm0
..B9.8:                         # Preds ..B9.7
                                # Execution count [9.00e-01]
        movq      -24(%rsp), %r15                               #[spill]
	.cfi_restore 15
        movq      -16(%rsp), %rbx                               #[spill]
	.cfi_restore 3
                                # LOE rbx rbp r12 r13 r14 r15
..B9.9:                         # Preds ..B9.8 ..B9.1
                                # Execution count [1.00e+00]
        vzeroupper                                              #337.1
        ret                                                     #337.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	rotate_dble_avx512,@function
	.size	rotate_dble_avx512,.-rotate_dble_avx512
	.data
# -- End  rotate_dble_avx512
	.text
# -- Begin  mulg5_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl mulg5_dble_avx512
# --- mulg5_dble_avx512(int, spinor_dble *)
mulg5_dble_avx512:
# parameter 1: %edi
# parameter 2: %rsi
..B10.1:                        # Preds ..B10.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_mulg5_dble_avx512.36:
..L37:
                                                         #340.1
        movslq    %edi, %rdi                                    #340.1
        lea       (%rdi,%rdi,2), %rax                           #343.8
        shlq      $6, %rax                                      #343.8
        addq      %rsi, %rax                                    #343.8
        cmpq      %rax, %rsi                                    #345.14
        jae       ..B10.5       # Prob 10%                      #345.14
                                # LOE rax rbx rbp rsi r12 r13 r14 r15
..B10.2:                        # Preds ..B10.1
                                # Execution count [9.00e-01]
        vpxord    %zmm1, %zmm1, %zmm1                           #350.25
        vxorpd    %ymm0, %ymm0, %ymm0                           #354.25
                                # LOE rax rbx rbp rsi r12 r13 r14 r15 ymm0 zmm1
..B10.3:                        # Preds ..B10.3 ..B10.2
                                # Execution count [5.00e+00]
        vsubpd    96(%rsi), %zmm1, %zmm2                        #350.10
        vsubpd    160(%rsi), %ymm0, %ymm3                       #354.10
        vmovups   %zmm2, 96(%rsi)                               #351.26
        vmovupd   %ymm3, 160(%rsi)                              #355.26
        addq      $192, %rsi                                    #345.18
        cmpq      %rax, %rsi                                    #345.14
        jb        ..B10.3       # Prob 82%                      #345.14
                                # LOE rax rbx rbp rsi r12 r13 r14 r15 ymm0 zmm1
..B10.5:                        # Preds ..B10.3 ..B10.1
                                # Execution count [1.00e+00]
        vzeroupper                                              #357.1
        ret                                                     #357.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	mulg5_dble_avx512,@function
	.size	mulg5_dble_avx512,.-mulg5_dble_avx512
	.data
# -- End  mulg5_dble_avx512
	.text
# -- Begin  mulmg5_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl mulmg5_dble_avx512
# --- mulmg5_dble_avx512(int, spinor_dble *)
mulmg5_dble_avx512:
# parameter 1: %edi
# parameter 2: %rsi
..B11.1:                        # Preds ..B11.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_mulmg5_dble_avx512.39:
..L40:
                                                         #360.1
        movslq    %edi, %rdi                                    #360.1
        lea       (%rdi,%rdi,2), %rax                           #363.8
        shlq      $6, %rax                                      #363.8
        addq      %rsi, %rax                                    #363.8
        cmpq      %rax, %rsi                                    #365.14
        jae       ..B11.5       # Prob 10%                      #365.14
                                # LOE rax rbx rbp rsi r12 r13 r14 r15
..B11.2:                        # Preds ..B11.1
                                # Execution count [9.00e-01]
        vpxord    %zmm1, %zmm1, %zmm1                           #370.25
        vxorpd    %ymm0, %ymm0, %ymm0                           #374.25
                                # LOE rax rbx rbp rsi r12 r13 r14 r15 ymm0 zmm1
..B11.3:                        # Preds ..B11.3 ..B11.2
                                # Execution count [5.00e+00]
        vsubpd    (%rsi), %zmm1, %zmm2                          #370.10
        vsubpd    64(%rsi), %ymm0, %ymm3                        #374.10
        vmovups   %zmm2, (%rsi)                                 #371.26
        vmovupd   %ymm3, 64(%rsi)                               #375.26
        addq      $192, %rsi                                    #365.18
        cmpq      %rax, %rsi                                    #365.14
        jb        ..B11.3       # Prob 82%                      #365.14
                                # LOE rax rbx rbp rsi r12 r13 r14 r15 ymm0 zmm1
..B11.5:                        # Preds ..B11.3 ..B11.1
                                # Execution count [1.00e+00]
        vzeroupper                                              #377.1
        ret                                                     #377.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	mulmg5_dble_avx512,@function
	.size	mulmg5_dble_avx512,.-mulmg5_dble_avx512
	.data
# -- End  mulmg5_dble_avx512
	.section .rodata, "a"
	.align 64
	.align 64
.L_2il0floatpacket.3:
	.long	0x00000000,0x3ff00000,0x00000000,0xbff00000,0x00000000,0x3ff00000,0x00000000,0xbff00000,0x00000000,0x3ff00000,0x00000000,0xbff00000,0x00000000,0x3ff00000,0x00000000,0xbff00000
	.type	.L_2il0floatpacket.3,@object
	.size	.L_2il0floatpacket.3,64
	.align 64
.L_2il0floatpacket.4:
	.long	0x00000000,0x3ff00000,0x00000000,0x3ff00000,0x00000000,0x3ff00000,0x00000000,0x3ff00000,0x00000000,0xbff00000,0x00000000,0xbff00000,0x00000000,0xbff00000,0x00000000,0xbff00000
	.type	.L_2il0floatpacket.4,@object
	.size	.L_2il0floatpacket.4,64
	.data
	.section .note.GNU-stack, ""
// -- Begin DWARF2 SEGMENT .eh_frame
	.section .eh_frame,"a",@progbits
.eh_frame_seg:
	.align 8
# End
