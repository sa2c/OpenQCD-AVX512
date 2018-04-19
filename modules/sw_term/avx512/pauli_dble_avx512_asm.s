# mark_description "Intel(R) C Intel(R) 64 Compiler for applications running on Intel(R) 64, Version 17.0.4.196 Build 20170411";
# mark_description "-I../../../include -I.. -I/cineca/prod/opt/compilers/intel/pe-xe-2017/binary/impi/2017.3.196/intel64/include";
# mark_description " -isystem /cineca/prod/opt/compilers/intel/pe-xe-2018/binary/impi/2018.1.163/include64/ -std=c89 -xCORE-AVX5";
# mark_description "12 -mtune=skylake -DAVX512 -O3 -Ddirac_counters -pedantic -fstrict-aliasing -Wno-long-long -Wstrict-prototyp";
# mark_description "es -S";
	.file "pauli_dble_avx512.c"
	.text
..TXTST0:
# -- Begin  mul_pauli2_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl mul_pauli2_dble_avx512
# --- mul_pauli2_dble_avx512(double, pauli_dble *, weyl_dble *, weyl_dble *)
mul_pauli2_dble_avx512:
# parameter 1: %xmm0
# parameter 2: %rdi
# parameter 3: %rsi
# parameter 4: %rdx
..B1.1:                         # Preds ..B1.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_mul_pauli2_dble_avx512.1:
..L2:
                                                          #20.1
        pushq     %rbp                                          #20.1
	.cfi_def_cfa_offset 16
        movq      %rsp, %rbp                                    #20.1
	.cfi_def_cfa 6, 16
	.cfi_offset 6, -16
        movl      $86, %eax                                     #81.8
        vmovups   .L_2il0floatpacket.9(%rip), %zmm17            #49.9
        vmovups   .L_2il0floatpacket.10(%rip), %zmm9            #51.9
        vmovups   (%rsi), %zmm22                                #31.28
        vmovups   64(%rsi), %zmm16                              #32.28
        vmovups   128(%rsi), %zmm14                             #33.28
        vmovups   .L_2il0floatpacket.13(%rip), %zmm5            #70.9
        vmovups   (%rdi), %zmm19                                #35.27
        vmovups   64(%rdi), %zmm29                              #36.27
        vmovups   288(%rdi), %zmm21                             #40.27
        vmovups   352(%rdi), %zmm6                              #41.27
        vmovups   .L_2il0floatpacket.14(%rip), %zmm23           #73.9
        vmovsd    %xmm0, -16(%rbp)                              #20.1
        vmovups   .L_2il0floatpacket.16(%rip), %zmm8            #78.9
        vmovups   -16(%rbp), %zmm13                             #46.27
        vmovups   .L_2il0floatpacket.11(%rip), %zmm11           #54.8
        vmovups   .L_2il0floatpacket.12(%rip), %zmm20           #57.9
        vmovups   128(%rdi), %zmm30                             #37.27
        vmovups   416(%rdi), %zmm10                             #42.27
        vmovups   480(%rdi), %zmm7                              #43.27
        vmovups   192(%rdi), %zmm27                             #38.27
        vmovups   256(%rdi), %zmm28                             #39.27
        vpermi2pd %zmm6, %zmm29, %zmm23                         #73.9
        vpermi2pd %zmm14, %zmm22, %zmm11                        #54.8
        kmovw     %eax, %k1                                     #81.8
        vmovaps   %zmm17, %zmm2                                 #50.8
        movl      $169, %eax                                    #82.8
        vmovaps   %zmm9, %zmm1                                  #52.8
        vpermi2pd %zmm16, %zmm22, %zmm2                         #50.8
        vpermi2pd %zmm16, %zmm22, %zmm1                         #52.8
        vpermi2pd %zmm14, %zmm16, %zmm17                        #56.8
        vpermt2pd %zmm14, %zmm9, %zmm16                         #60.8
        vpermt2pd %zmm14, %zmm20, %zmm22                        #58.8
        vmovups   .L_2il0floatpacket.15(%rip), %zmm9            #75.9
        kmovw     %eax, %k2                                     #82.8
        vmovaps   %zmm5, %zmm25                                 #71.9
        movl      $106, %eax                                    #93.8
        vpermi2pd %zmm21, %zmm19, %zmm25                        #71.9
        kmovw     %eax, %k3                                     #93.8
        vmovaps   %zmm9, %zmm31                                 #76.10
        movl      $149, %eax                                    #94.8
        vmovaps   %zmm8, %zmm26                                 #79.10
        vpermi2pd %zmm23, %zmm25, %zmm31                        #76.10
        vpermi2pd %zmm23, %zmm13, %zmm26                        #79.10
        kmovw     %eax, %k4                                     #94.8
        vmulpd    %zmm2, %zmm31, %zmm24                         #77.8
        vmovups   .L_2il0floatpacket.18(%rip), %zmm31           #88.10
        vpermilpd $85, %zmm2, %zmm4                             #62.9
        movl      $85, %eax                                     #147.8
        vmulpd    %zmm4, %zmm26, %zmm12                         #80.10
        vmovups   .L_2il0floatpacket.17(%rip), %zmm26           #84.9
        vaddpd    %zmm12, %zmm24, %zmm24{%k1}                   #81.8
        vpermt2pd %zmm21, %zmm26, %zmm19                        #85.9
        vpermi2pd %zmm10, %zmm30, %zmm26                        #98.9
        vsubpd    %zmm12, %zmm24, %zmm24{%k2}                   #82.8
        vpermi2pd %zmm19, %zmm23, %zmm31                        #88.10
        vpermi2pd %zmm26, %zmm13, %zmm8                         #104.10
        vmovups   .L_2il0floatpacket.19(%rip), %zmm12           #90.9
        vfmadd213pd %zmm24, %zmm22, %zmm31                      #89.8
        vmovaps   %zmm12, %zmm21                                #91.10
        vpermi2pd %zmm13, %zmm23, %zmm21                        #91.10
        vpermi2pd %zmm13, %zmm26, %zmm12                        #113.10
        vpermilpd $85, %zmm22, %zmm15                           #66.9
        vmulpd    %zmm15, %zmm21, %zmm0                         #92.10
        vmovups   .L_2il0floatpacket.20(%rip), %zmm21           #101.10
        vaddpd    %zmm0, %zmm31, %zmm31{%k3}                    #93.8
        vpermi2pd %zmm26, %zmm25, %zmm21                        #101.10
        vsubpd    %zmm0, %zmm31, %zmm31{%k4}                    #94.8
        vmulpd    %zmm1, %zmm21, %zmm21                         #102.8
        vpermilpd $85, %zmm1, %zmm3                             #63.9
        vmulpd    %zmm3, %zmm8, %zmm24                          #105.10
        vmovups   .L_2il0floatpacket.21(%rip), %zmm8            #110.10
        vaddpd    %zmm24, %zmm21, %zmm21{%k1}                   #106.8
        vpermi2pd %zmm19, %zmm26, %zmm8                         #110.10
        vsubpd    %zmm24, %zmm21, %zmm21{%k2}                   #107.8
        vmovups   .L_2il0floatpacket.23(%rip), %zmm24           #123.10
        vfmadd213pd %zmm21, %zmm17, %zmm8                       #111.8
        vpermilpd $85, %zmm17, %zmm18                           #65.9
        vmulpd    %zmm18, %zmm12, %zmm21                        #114.10
        vmovups   .L_2il0floatpacket.24(%rip), %zmm12           #126.10
        vaddpd    %zmm21, %zmm8, %zmm8{%k3}                     #115.8
        vsubpd    %zmm21, %zmm8, %zmm8{%k4}                     #116.8
        vmovups   .L_2il0floatpacket.22(%rip), %zmm21           #119.9
        vmovaps   %zmm21, %zmm0                                 #120.9
        vpermi2pd %zmm7, %zmm27, %zmm0                          #120.9
        vpermi2pd %zmm0, %zmm19, %zmm24                         #123.10
        vpermi2pd %zmm0, %zmm13, %zmm12                         #126.10
        vmulpd    %zmm11, %zmm24, %zmm24                        #124.8
        vpermilpd $85, %zmm11, %zmm14                           #64.9
        vmulpd    %zmm14, %zmm12, %zmm12                        #127.10
        vaddpd    %zmm12, %zmm24, %zmm24{%k1}                   #128.8
        kmovw     %eax, %k1                                     #147.8
        vsubpd    %zmm12, %zmm24, %zmm24{%k2}                   #129.8
        vmovups   .L_2il0floatpacket.25(%rip), %zmm12           #132.10
        vpermi2pd %zmm19, %zmm0, %zmm12                         #132.10
        movl      $170, %eax                                    #148.8
        vmovups   .L_2il0floatpacket.26(%rip), %zmm19           #135.10
        vfmadd213pd %zmm24, %zmm16, %zmm12                      #133.8
        vpermi2pd %zmm13, %zmm0, %zmm19                         #135.10
        kmovw     %eax, %k7                                     #148.8
        vpermilpd $85, %zmm16, %zmm20                           #67.9
        movl      $90, %eax                                     #166.8
        vmulpd    %zmm20, %zmm19, %zmm13                        #136.10
        vmovups   .L_2il0floatpacket.27(%rip), %zmm19           #141.9
        kmovw     %eax, %k5                                     #166.8
        vaddpd    %zmm13, %zmm12, %zmm12{%k3}                   #137.8
        vsubpd    %zmm13, %zmm12, %zmm12{%k4}                   #138.8
        movl      $165, %eax                                    #167.8
        kmovw     %eax, %k6                                     #167.8
        vmovaps   %zmm19, %zmm13                                #142.10
        movl      $240, %eax                                    #272.8
        vpermi2pd %zmm23, %zmm25, %zmm13                        #142.10
        kmovw     %eax, %k2                                     #272.8
        vfmadd213pd %zmm8, %zmm2, %zmm13                        #143.8
        vmovups   .L_2il0floatpacket.28(%rip), %zmm8            #144.9
        vmovaps   %zmm8, %zmm24                                 #145.10
        vpermi2pd %zmm23, %zmm25, %zmm24                        #145.10
        vmulpd    %zmm24, %zmm4, %zmm24                         #146.10
        vaddpd    %zmm24, %zmm13, %zmm13{%k1}                   #147.8
        vsubpd    %zmm24, %zmm13, %zmm13{%k7}                   #148.8
        vmovaps   %zmm9, %zmm24                                 #151.10
        vpermi2pd %zmm0, %zmm23, %zmm24                         #151.10
        vfmadd213pd %zmm31, %zmm17, %zmm24                      #152.8
        vmovups   .L_2il0floatpacket.29(%rip), %zmm31           #153.9
        vpermt2pd %zmm0, %zmm31, %zmm23                         #154.10
        vmulpd    %zmm23, %zmm18, %zmm23                        #155.10
        vaddpd    %zmm23, %zmm24, %zmm24{%k7}                   #156.8
        vsubpd    %zmm23, %zmm24, %zmm24{%k1}                   #157.8
        vmovaps   %zmm19, %zmm23                                #161.10
        vpermi2pd %zmm26, %zmm25, %zmm23                        #161.10
        vpermt2pd %zmm26, %zmm8, %zmm25                         #164.10
        vfmadd213pd %zmm24, %zmm1, %zmm23                       #162.8
        vmulpd    %zmm25, %zmm3, %zmm25                         #165.10
        vaddpd    %zmm25, %zmm23, %zmm23{%k5}                   #166.8
        vsubpd    %zmm25, %zmm23, %zmm23{%k6}                   #167.8
        vmovaps   %zmm9, %zmm25                                 #170.10
        vpermi2pd %zmm0, %zmm26, %zmm25                         #170.10
        vpermt2pd %zmm0, %zmm31, %zmm26                         #173.10
        vfmadd213pd %zmm13, %zmm22, %zmm25                      #171.8
        vmulpd    %zmm26, %zmm15, %zmm26                        #174.10
        vaddpd    %zmm26, %zmm25, %zmm25{%k5}                   #175.8
        vsubpd    %zmm26, %zmm25, %zmm25{%k6}                   #176.8
        vmovups   .L_2il0floatpacket.30(%rip), %zmm26           #178.9
        vmovaps   %zmm26, %zmm0                                 #179.10
        vpermi2pd %zmm6, %zmm29, %zmm0                          #179.10
        vpermi2pd %zmm10, %zmm30, %zmm26                        #189.10
        vfmadd213pd %zmm12, %zmm2, %zmm0                        #180.8
        vmovups   .L_2il0floatpacket.31(%rip), %zmm12           #181.9
        vmovaps   %zmm12, %zmm2                                 #182.10
        vpermi2pd %zmm6, %zmm29, %zmm2                          #182.10
        vpermi2pd %zmm10, %zmm30, %zmm12                        #192.10
        vpermt2pd %zmm6, %zmm5, %zmm29                          #200.9
        vmulpd    %zmm2, %zmm4, %zmm4                           #183.10
        vaddpd    %zmm4, %zmm0, %zmm0{%k1}                      #184.8
        vsubpd    %zmm4, %zmm0, %zmm0{%k7}                      #185.8
        vfmadd213pd %zmm0, %zmm1, %zmm26                        #190.8
        vmulpd    %zmm12, %zmm3, %zmm1                          #193.10
        vmovups   .L_2il0floatpacket.32(%rip), %zmm3            #201.9
        vaddpd    %zmm1, %zmm26, %zmm26{%k1}                    #194.8
        vpermt2pd %zmm7, %zmm3, %zmm27                          #202.9
        vmovups   .L_2il0floatpacket.33(%rip), %zmm7            #214.9
        vsubpd    %zmm1, %zmm26, %zmm26{%k7}                    #195.8
        vpermt2pd 544(%rdi), %zmm7, %zmm28                      #215.9
        vmovaps   %zmm31, %zmm5                                 #208.10
        vpermi2pd %zmm27, %zmm29, %zmm5                         #208.10
        vmulpd    %zmm5, %zmm14, %zmm6                          #209.10
        vmovaps   %zmm9, %zmm0                                  #205.10
        vpermi2pd %zmm27, %zmm29, %zmm0                         #205.10
        vfmadd213pd %zmm23, %zmm11, %zmm0                       #206.8
        vaddpd    %zmm6, %zmm0, %zmm0{%k5}                      #210.8
        vmovaps   %zmm19, %zmm1                                 #236.10
        vsubpd    %zmm6, %zmm0, %zmm0{%k6}                      #211.8
        vpermi2pd %zmm28, %zmm29, %zmm1                         #236.10
        vpermt2pd %zmm28, %zmm8, %zmm29                         #239.10
        vfmadd213pd %zmm0, %zmm16, %zmm1                        #237.8
        vmovups   .L_2il0floatpacket.34(%rip), %zmm0            #245.9
        vmulpd    %zmm29, %zmm20, %zmm29                        #240.10
        vpermt2pd %zmm10, %zmm0, %zmm30                         #246.9
        vaddpd    %zmm29, %zmm1, %zmm1{%k7}                     #241.8
        vmovaps   %zmm19, %zmm5                                 #218.10
        vpermi2pd %zmm28, %zmm27, %zmm5                         #218.10
        vpermi2pd %zmm27, %zmm30, %zmm19                        #249.10
        vsubpd    %zmm29, %zmm1, %zmm1{%k1}                     #242.8
        vfmadd213pd %zmm26, %zmm22, %zmm5                       #219.8
        vfmadd213pd %zmm25, %zmm11, %zmm19                      #250.8
        vmovaps   %zmm8, %zmm22                                 #221.10
        vpermi2pd %zmm28, %zmm27, %zmm22                        #221.10
        vpermi2pd %zmm27, %zmm30, %zmm8                         #252.10
        vmulpd    %zmm22, %zmm15, %zmm15                        #222.10
        vmulpd    %zmm8, %zmm14, %zmm11                         #253.10
        vaddpd    %zmm15, %zmm5, %zmm5{%k5}                     #223.8
        vaddpd    %zmm11, %zmm19, %zmm19{%k5}                   #254.8
        vsubpd    %zmm15, %zmm5, %zmm5{%k6}                     #224.8
        vsubpd    %zmm11, %zmm19, %zmm19{%k6}                   #255.8
        vmovaps   %zmm31, %zmm6                                 #230.10
        vmovaps   %zmm9, %zmm15                                 #227.10
        vpermi2pd %zmm28, %zmm27, %zmm6                         #230.10
        vpermi2pd %zmm28, %zmm30, %zmm9                         #258.10
        vpermt2pd %zmm28, %zmm31, %zmm30                        #261.10
        vpermi2pd %zmm28, %zmm27, %zmm15                        #227.10
        vmulpd    %zmm6, %zmm18, %zmm18                         #231.10
        vfmadd213pd %zmm19, %zmm16, %zmm9                       #259.8
        vfmadd213pd %zmm5, %zmm17, %zmm15                       #228.8
        vmovups   .L_2il0floatpacket.35(%rip), %zmm28           #268.9
        vmulpd    %zmm30, %zmm20, %zmm16                        #262.10
        vaddpd    %zmm18, %zmm15, %zmm15{%k5}                   #232.8
        vaddpd    %zmm16, %zmm9, %zmm9{%k7}                     #263.8
        vsubpd    %zmm18, %zmm15, %zmm15{%k6}                   #233.8
        vsubpd    %zmm16, %zmm9, %zmm9{%k1}                     #264.8
        vpermi2pd %zmm1, %zmm15, %zmm28                         #269.8
        vpermi2pd %zmm9, %zmm1, %zmm7                           #267.8
        vpermt2pd %zmm15, %zmm21, %zmm9                         #271.8
        vblendmpd %zmm28, %zmm7, %zmm10{%k2}                    #272.8
        vblendmpd %zmm7, %zmm9, %zmm27{%k2}                     #273.8
        vblendmpd %zmm9, %zmm28, %zmm30{%k2}                    #274.8
        vmovups   %zmm10, (%rdx)                                #276.22
        vmovups   %zmm27, 64(%rdx)                              #277.22
        vmovups   %zmm30, 128(%rdx)                             #278.22
        vzeroupper                                              #279.1
        movq      %rbp, %rsp                                    #279.1
        popq      %rbp                                          #279.1
	.cfi_restore 6
        ret                                                     #279.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	mul_pauli2_dble_avx512,@function
	.size	mul_pauli2_dble_avx512,.-mul_pauli2_dble_avx512
	.data
# -- End  mul_pauli2_dble_avx512
	.text
# -- Begin  fwd_house_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl fwd_house_avx512
# --- fwd_house_avx512(double, complex_dble *, complex_dble *, double *)
fwd_house_avx512:
# parameter 1: %xmm0
# parameter 2: %rdi
# parameter 3: %rsi
# parameter 4: %rdx
..B2.1:                         # Preds ..B2.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_fwd_house_avx512.8:
..L9:
                                                          #283.1
        pushq     %r12                                          #283.1
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
        pushq     %r13                                          #283.1
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
        pushq     %r14                                          #283.1
	.cfi_def_cfa_offset 32
	.cfi_offset 14, -32
        pushq     %r15                                          #283.1
	.cfi_def_cfa_offset 40
	.cfi_offset 15, -40
        pushq     %rbx                                          #283.1
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
        pushq     %rbp                                          #283.1
	.cfi_def_cfa_offset 56
	.cfi_offset 6, -56
        xorl      %eax, %eax                                    #288.3
        xorl      %r8d, %r8d                                    #290.3
        movq      %rdi, %r9                                     #283.1
        xorl      %r11d, %r11d                                  #290.3
        vmovapd   %xmm0, %xmm14                                 #283.1
        xorl      %r10d, %r10d                                  #290.3
        vxorpd    %xmm1, %xmm1, %xmm1                           #326.12
        vmovsd    .L_2il0floatpacket.38(%rip), %xmm11           #307.12
        xorl      %edi, %edi                                    #290.3
        vmovsd    .L_2il0floatpacket.36(%rip), %xmm0            #306.16
                                # LOE rdx rsi r8 r9 eax edi r10d r11d xmm0 xmm1 xmm11 xmm14
..B2.2:                         # Preds ..B2.35 ..B2.1
                                # Execution count [5.00e+00]
        movslq    %r10d, %r12                                   #292.29
        lea       1(%r8), %ecx                                  #295.10
        shlq      $4, %r12                                      #291.10
        vmovsd    8(%r9,%r12), %xmm3                            #292.29
        vmulsd    %xmm3, %xmm3, %xmm12                          #292.29
        vmovsd    (%r9,%r12), %xmm2                             #291.29
        vfmadd231sd %xmm2, %xmm2, %xmm12                        #291.5
        vsqrtsd   %xmm12, %xmm12, %xmm13                        #293.10
                                # LOE rdx rsi r8 r9 r12 eax ecx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..B2.3:                         # Preds ..B2.2
                                # Execution count [5.00e+00]
        xorl      %r13d, %r13d                                  #295.5
        lea       5(%r11), %r14d                                #295.5
        movl      %r14d, %ebp                                   #295.5
        movl      $1, %ebx                                      #295.5
        sarl      $2, %ebp                                      #295.5
        shrl      $29, %ebp                                     #295.5
        lea       5(%rbp,%r11), %r15d                           #295.5
        xorl      %ebp, %ebp                                    #296.7
        sarl      $3, %r15d                                     #295.5
        testl     %r15d, %r15d                                  #295.5
        jbe       ..B2.7        # Prob 10%                      #295.5
                                # LOE rdx rsi r8 r9 r12 eax ecx ebx ebp edi r10d r11d r13d r14d r15d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..B2.4:                         # Preds ..B2.3
                                # Execution count [1.56e-02]
        vxorpd    %xmm10, %xmm10, %xmm10                        #295.5
        vxorpd    %xmm9, %xmm9, %xmm9                           #295.5
        vxorpd    %xmm8, %xmm8, %xmm8                           #295.5
        vxorpd    %xmm4, %xmm4, %xmm4                           #295.5
        vxorpd    %xmm7, %xmm7, %xmm7                           #295.5
        vxorpd    %xmm6, %xmm6, %xmm6                           #295.5
        vxorpd    %xmm5, %xmm5, %xmm5                           #295.5
                                # LOE rdx rsi r8 r9 r12 eax ecx ebp edi r10d r11d r13d r14d r15d xmm0 xmm1 xmm2 xmm3 xmm4 xmm5 xmm6 xmm7 xmm8 xmm9 xmm10 xmm11 xmm12 xmm13 xmm14
..B2.5:                         # Preds ..B2.5 ..B2.4
                                # Execution count [3.12e+00]
        incl      %r13d                                         #295.5
        lea       (%r10,%rbp), %ebx                             #297.33
        movslq    %ebx, %rbx                                    #296.14
        addl      $48, %ebp                                     #295.5
        shlq      $4, %rbx                                      #297.33
        vmovsd    104(%r9,%rbx), %xmm15                         #297.14
        vmovsd    200(%r9,%rbx), %xmm18                         #297.14
        vmulsd    %xmm15, %xmm15, %xmm17                        #296.7
        vmulsd    %xmm18, %xmm18, %xmm20                        #296.7
        vmovsd    192(%r9,%rbx), %xmm19                         #296.14
        vmovsd    96(%r9,%rbx), %xmm16                          #296.14
        vfmadd231sd %xmm16, %xmm16, %xmm17                      #296.7
        vmovsd    296(%r9,%rbx), %xmm21                         #297.14
        vmovsd    392(%r9,%rbx), %xmm24                         #297.14
        vmovsd    488(%r9,%rbx), %xmm27                         #297.14
        vmovsd    584(%r9,%rbx), %xmm30                         #297.14
        vfmadd231sd %xmm19, %xmm19, %xmm20                      #296.7
        vaddsd    %xmm12, %xmm17, %xmm12                        #296.7
        vmulsd    %xmm21, %xmm21, %xmm23                        #296.7
        vmulsd    %xmm24, %xmm24, %xmm26                        #296.7
        vmulsd    %xmm27, %xmm27, %xmm29                        #296.7
        vaddsd    %xmm10, %xmm20, %xmm10                        #296.7
        vmulsd    %xmm30, %xmm30, %xmm15                        #296.7
        vmovsd    680(%r9,%rbx), %xmm16                         #297.14
        vmovsd    776(%r9,%rbx), %xmm19                         #297.14
        vmulsd    %xmm16, %xmm16, %xmm18                        #296.7
        vmulsd    %xmm19, %xmm19, %xmm21                        #296.7
        vmovsd    768(%r9,%rbx), %xmm20                         #296.14
        vmovsd    288(%r9,%rbx), %xmm22                         #296.14
        vmovsd    384(%r9,%rbx), %xmm25                         #296.14
        vmovsd    480(%r9,%rbx), %xmm28                         #296.14
        vmovsd    576(%r9,%rbx), %xmm31                         #296.14
        vmovsd    672(%r9,%rbx), %xmm17                         #296.14
        vfmadd231sd %xmm22, %xmm22, %xmm23                      #296.7
        vfmadd231sd %xmm25, %xmm25, %xmm26                      #296.7
        vfmadd231sd %xmm28, %xmm28, %xmm29                      #296.7
        vfmadd231sd %xmm31, %xmm31, %xmm15                      #296.7
        vfmadd231sd %xmm17, %xmm17, %xmm18                      #296.7
        vfmadd231sd %xmm20, %xmm20, %xmm21                      #296.7
        vaddsd    %xmm9, %xmm23, %xmm9                          #296.7
        vaddsd    %xmm8, %xmm26, %xmm8                          #296.7
        vaddsd    %xmm4, %xmm29, %xmm4                          #296.7
        vaddsd    %xmm7, %xmm15, %xmm7                          #296.7
        vaddsd    %xmm6, %xmm18, %xmm6                          #296.7
        vaddsd    %xmm5, %xmm21, %xmm5                          #296.7
        cmpl      %r15d, %r13d                                  #295.5
        jb        ..B2.5        # Prob 99%                      #295.5
                                # LOE rdx rsi r8 r9 r12 eax ecx ebp edi r10d r11d r13d r14d r15d xmm0 xmm1 xmm2 xmm3 xmm4 xmm5 xmm6 xmm7 xmm8 xmm9 xmm10 xmm11 xmm12 xmm13 xmm14
..B2.6:                         # Preds ..B2.5
                                # Execution count [4.50e+00]
        vaddsd    %xmm10, %xmm12, %xmm10                        #295.5
        vaddsd    %xmm8, %xmm9, %xmm8                           #295.5
        vaddsd    %xmm7, %xmm4, %xmm4                           #295.5
        vaddsd    %xmm5, %xmm6, %xmm5                           #295.5
        vaddsd    %xmm8, %xmm10, %xmm9                          #295.5
        vaddsd    %xmm5, %xmm4, %xmm6                           #295.5
        vaddsd    %xmm6, %xmm9, %xmm12                          #295.5
        lea       1(,%r13,8), %ebx                              #296.7
                                # LOE rdx rsi r8 r9 r12 eax ecx ebx edi r10d r11d r14d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..B2.7:                         # Preds ..B2.6 ..B2.3
                                # Execution count [5.00e+00]
        cmpl      %r14d, %ebx                                   #295.5
        ja        ..B2.23       # Prob 50%                      #295.5
                                # LOE rdx rsi r8 r9 r12 eax ecx ebx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..B2.8:                         # Preds ..B2.7
                                # Execution count [0.00e+00]
        lea       (%r8,%rbx), %ebp                              #295.5
        negl      %ebp                                          #295.5
        addl      $5, %ebp                                      #295.5
        jmp       *.2.10_2.switchtab.4(,%rbp,8)                 #295.5
                                # LOE rdx rsi r8 r9 r12 eax ecx ebx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..1.10_0.TAG.6:
..B2.10:                        # Preds ..B2.8
                                # Execution count [0.00e+00]
        lea       (%rbx,%rbx,2), %ebp                           #296.14
        lea       (%r10,%rbp,2), %r13d                          #297.33
        movslq    %r13d, %r13                                   #296.14
        shlq      $4, %r13                                      #297.33
        lea       584(%r9,%r13), %r14                           #297.14
        vmovsd    (%r14), %xmm4                                 #297.14
        vmulsd    %xmm4, %xmm4, %xmm6                           #297.33
        vmovsd    -8(%r14), %xmm5                               #296.14
        vfmadd231sd %xmm5, %xmm5, %xmm6                         #296.7
        vaddsd    %xmm6, %xmm12, %xmm12                         #296.7
                                # LOE rdx rsi r8 r9 r12 eax ecx ebx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..1.10_0.TAG.5:
..B2.12:                        # Preds ..B2.8 ..B2.10
                                # Execution count [0.00e+00]
        lea       (%rbx,%rbx,2), %ebp                           #296.14
        lea       (%r10,%rbp,2), %r13d                          #297.33
        movslq    %r13d, %r13                                   #296.14
        shlq      $4, %r13                                      #297.33
        lea       488(%r9,%r13), %r14                           #297.14
        vmovsd    (%r14), %xmm4                                 #297.14
        vmulsd    %xmm4, %xmm4, %xmm6                           #297.33
        vmovsd    -8(%r14), %xmm5                               #296.14
        vfmadd231sd %xmm5, %xmm5, %xmm6                         #296.7
        vaddsd    %xmm6, %xmm12, %xmm12                         #296.7
                                # LOE rdx rsi r8 r9 r12 eax ecx ebx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..1.10_0.TAG.4:
..B2.14:                        # Preds ..B2.8 ..B2.12
                                # Execution count [0.00e+00]
        lea       (%rbx,%rbx,2), %ebp                           #296.14
        lea       (%r10,%rbp,2), %r13d                          #297.33
        movslq    %r13d, %r13                                   #296.14
        shlq      $4, %r13                                      #297.33
        lea       392(%r9,%r13), %r14                           #297.14
        vmovsd    (%r14), %xmm4                                 #297.14
        vmulsd    %xmm4, %xmm4, %xmm6                           #297.33
        vmovsd    -8(%r14), %xmm5                               #296.14
        vfmadd231sd %xmm5, %xmm5, %xmm6                         #296.7
        vaddsd    %xmm6, %xmm12, %xmm12                         #296.7
                                # LOE rdx rsi r8 r9 r12 eax ecx ebx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..1.10_0.TAG.3:
..B2.16:                        # Preds ..B2.8 ..B2.14
                                # Execution count [0.00e+00]
        lea       (%rbx,%rbx,2), %ebp                           #296.14
        lea       (%r10,%rbp,2), %r13d                          #297.33
        movslq    %r13d, %r13                                   #296.14
        shlq      $4, %r13                                      #297.33
        lea       296(%r9,%r13), %r14                           #297.14
        vmovsd    (%r14), %xmm4                                 #297.14
        vmulsd    %xmm4, %xmm4, %xmm6                           #297.33
        vmovsd    -8(%r14), %xmm5                               #296.14
        vfmadd231sd %xmm5, %xmm5, %xmm6                         #296.7
        vaddsd    %xmm6, %xmm12, %xmm12                         #296.7
                                # LOE rdx rsi r8 r9 r12 eax ecx ebx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..1.10_0.TAG.2:
..B2.18:                        # Preds ..B2.8 ..B2.16
                                # Execution count [0.00e+00]
        lea       (%rbx,%rbx,2), %ebp                           #296.14
        lea       (%r10,%rbp,2), %r13d                          #297.33
        movslq    %r13d, %r13                                   #296.14
        shlq      $4, %r13                                      #297.33
        lea       200(%r9,%r13), %r14                           #297.14
        vmovsd    (%r14), %xmm4                                 #297.14
        vmulsd    %xmm4, %xmm4, %xmm6                           #297.33
        vmovsd    -8(%r14), %xmm5                               #296.14
        vfmadd231sd %xmm5, %xmm5, %xmm6                         #296.7
        vaddsd    %xmm6, %xmm12, %xmm12                         #296.7
                                # LOE rdx rsi r8 r9 r12 eax ecx ebx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..1.10_0.TAG.1:
..B2.20:                        # Preds ..B2.8 ..B2.18
                                # Execution count [0.00e+00]
        lea       (%rbx,%rbx,2), %ebp                           #296.14
        lea       (%r10,%rbp,2), %r13d                          #297.33
        movslq    %r13d, %r13                                   #296.14
        shlq      $4, %r13                                      #297.33
        vmovsd    104(%r9,%r13), %xmm4                          #297.14
        vmulsd    %xmm4, %xmm4, %xmm6                           #297.33
        vmovsd    96(%r9,%r13), %xmm5                           #296.14
        vfmadd231sd %xmm5, %xmm5, %xmm6                         #296.7
        vaddsd    %xmm6, %xmm12, %xmm12                         #296.7
                                # LOE rdx rsi r8 r9 r12 eax ecx ebx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..1.10_0.TAG.0:
..B2.22:                        # Preds ..B2.8 ..B2.20
                                # Execution count [4.50e+00]
        lea       (%rbx,%rbx,2), %ebx                           #296.14
        lea       (%r10,%rbx,2), %ebp                           #297.33
        movslq    %ebp, %rbp                                    #296.14
        shlq      $4, %rbp                                      #297.33
        vmovsd    8(%r9,%rbp), %xmm4                            #297.14
        vmulsd    %xmm4, %xmm4, %xmm6                           #297.33
        vmovsd    (%r9,%rbp), %xmm5                             #296.14
        vfmadd231sd %xmm5, %xmm5, %xmm6                         #296.7
        vaddsd    %xmm6, %xmm12, %xmm12                         #296.7
                                # LOE rdx rsi r8 r9 r12 eax ecx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..B2.23:                        # Preds ..B2.22 ..B2.7
                                # Execution count [5.00e+00]
        vcomisd   %xmm14, %xmm12                                #299.15
        jb        ..B2.25       # Prob 50%                      #299.15
                                # LOE rdx rsi r8 r9 r12 eax ecx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..B2.24:                        # Preds ..B2.23
                                # Execution count [2.50e+00]
        vsqrtsd   %xmm12, %xmm12, %xmm12                        #300.12
        jmp       ..B2.26       # Prob 100%                     #300.12
                                # LOE rdx rsi r8 r9 r12 eax ecx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..B2.25:                        # Preds ..B2.23
                                # Execution count [2.50e+00]
        vmovapd   %xmm11, %xmm12                                #303.7
        movl      $1, %eax                                      #302.7
                                # LOE rdx rsi r8 r9 r12 eax ecx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..B2.26:                        # Preds ..B2.24 ..B2.25
                                # Execution count [5.00e+00]
        vmulsd    %xmm0, %xmm12, %xmm4                          #306.30
        vcomisd   %xmm4, %xmm13                                 #306.30
        jb        ..B2.28       # Prob 50%                      #306.30
                                # LOE rdx rsi r8 r9 r12 eax ecx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm11 xmm12 xmm13 xmm14
..B2.27:                        # Preds ..B2.26
                                # Execution count [2.50e+00]
        vdivsd    %xmm13, %xmm11, %xmm4                         #307.18
        vmulsd    %xmm4, %xmm2, %xmm5                           #308.19
        vmulsd    %xmm3, %xmm4, %xmm4                           #309.19
        jmp       ..B2.29       # Prob 100%                     #309.19
                                # LOE rdx rsi r8 r9 r12 eax ecx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm4 xmm5 xmm11 xmm12 xmm13 xmm14
..B2.28:                        # Preds ..B2.26
                                # Execution count [2.50e+00]
        vmovapd   %xmm11, %xmm5                                 #311.7
        vxorpd    %xmm4, %xmm4, %xmm4                           #312.7
                                # LOE rdx rsi r8 r9 r12 eax ecx edi r10d r11d xmm0 xmm1 xmm2 xmm3 xmm4 xmm5 xmm11 xmm12 xmm13 xmm14
..B2.29:                        # Preds ..B2.27 ..B2.28
                                # Execution count [6.63e-01]
        vfmadd231sd %xmm12, %xmm4, %xmm3                        #316.5
        xorl      %ebp, %ebp                                    #323.5
        vfmadd231sd %xmm12, %xmm5, %xmm2                        #315.5
        vmovsd    %xmm3, 8(%r9,%r12)                            #316.5
        vaddsd    %xmm13, %xmm12, %xmm3                         #318.28
        vmulsd    %xmm3, %xmm12, %xmm12                         #318.28
        lea       6(%r11), %r13d                                #331.23
        vmulsd    %xmm3, %xmm5, %xmm5                           #320.5
        vmulsd    %xmm3, %xmm4, %xmm4                           #321.28
        vdivsd    %xmm12, %xmm11, %xmm6                         #318.28
        vmovsd    %xmm2, (%r9,%r12)                             #315.5
        movq      %r8, %rbx                                     #320.5
        vmulsd    %xmm6, %xmm5, %xmm2                           #320.5
        vmulsd    %xmm6, %xmm4, %xmm8                           #321.33
        movslq    %edi, %r12                                    #328.12
        shlq      $4, %rbx                                      #320.5
        addq      %r8, %r12                                     #328.12
        shlq      $4, %r12                                      #328.12
        vxorpd    .L_2il0floatpacket.37(%rip), %xmm2, %xmm7     #320.5
        vmovsd    %xmm6, (%rdx,%r8,8)                           #319.5
        negq      %r8                                           #323.27
        vmovsd    %xmm6, -24(%rsp)                              #318.5
        addq      $5, %r8                                       #323.27
        vmovsd    %xmm7, (%rsi,%rbx)                            #320.5
        vmovsd    %xmm8, 8(%rsi,%rbx)                           #321.5
        lea       (%r9,%r12), %rbx                              #328.12
        vmovddup  -24(%rsp), %xmm2                              #343.28
        lea       16(%r12,%r9), %r12                            #329.12
        movq      %rdx, -16(%rsp)                               #331.23[spill]
                                # LOE rbx rbp rsi r8 r9 r12 eax ecx edi r10d r11d r13d xmm0 xmm1 xmm2 xmm11 xmm14
..B2.30:                        # Preds ..B2.34 ..B2.29
                                # Execution count [2.12e+01]
        vmovapd   %xmm1, %xmm3                                  #326.12
        movq      %rbx, %r15                                    #328.7
        movq      %r12, %r14                                    #329.7
        xorl      %edx, %edx                                    #331.7
                                # LOE rbx rbp rsi r8 r9 r12 r14 r15 eax edx ecx edi r10d r11d r13d xmm0 xmm1 xmm2 xmm3 xmm11 xmm14
..B2.31:                        # Preds ..B2.31 ..B2.30
                                # Execution count [1.18e+02]
        vmovupd   (%r14), %xmm5                                 #334.27
        incl      %edx                                          #331.7
        vmulpd    8(%r15){1to2}, %xmm5, %xmm4                   #335.14
        vpermilpd $1, %xmm4, %xmm6                              #336.14
        addq      $96, %r14                                     #340.9
        vfmsubadd231pd (%r15){1to2}, %xmm5, %xmm6               #337.14
        addq      $96, %r15                                     #339.9
        vaddpd    %xmm3, %xmm6, %xmm3                           #338.14
        cmpl      %r13d, %edx                                   #331.7
        jb        ..B2.31       # Prob 82%                      #331.7
                                # LOE rbx rbp rsi r8 r9 r12 r14 r15 eax edx ecx edi r10d r11d r13d xmm0 xmm1 xmm2 xmm3 xmm11 xmm14
..B2.32:                        # Preds ..B2.31
                                # Execution count [2.25e+01]
        vmulpd    %xmm2, %xmm3, %xmm3                           #344.12
        movq      %rbx, %r15                                    #347.7
        movq      %r12, %r14                                    #348.7
        xorl      %edx, %edx                                    #349.7
                                # LOE rbx rbp rsi r8 r9 r12 r14 r15 eax edx ecx edi r10d r11d r13d xmm0 xmm1 xmm2 xmm3 xmm11 xmm14
..B2.33:                        # Preds ..B2.33 ..B2.32
                                # Execution count [1.25e+02]
        vmulpd    8(%r15){1to2}, %xmm3, %xmm4                   #353.14
        vpermilpd $1, %xmm4, %xmm6                              #354.14
        incl      %edx                                          #349.7
        vfmaddsub231pd (%r15){1to2}, %xmm3, %xmm6               #355.14
        addq      $96, %r15                                     #358.9
        vmovupd   (%r14), %xmm5                                 #352.27
        vsubpd    %xmm6, %xmm5, %xmm7                           #356.14
        vmovupd   %xmm7, (%r14)                                 #357.25
        addq      $96, %r14                                     #359.9
        cmpl      %r13d, %edx                                   #349.7
        jb        ..B2.33       # Prob 82%                      #349.7
                                # LOE rbx rbp rsi r8 r9 r12 r14 r15 eax edx ecx edi r10d r11d r13d xmm0 xmm1 xmm2 xmm3 xmm11 xmm14
..B2.34:                        # Preds ..B2.33
                                # Execution count [2.50e+01]
        incq      %rbp                                          #323.5
        addq      $16, %r12                                     #323.5
        cmpq      %r8, %rbp                                     #323.5
        jb        ..B2.30       # Prob 81%                      #323.5
                                # LOE rbx rbp rsi r8 r9 r12 eax ecx edi r10d r11d r13d xmm0 xmm1 xmm2 xmm11 xmm14
..B2.35:                        # Preds ..B2.34
                                # Execution count [5.00e+00]
        decl      %r11d                                         #295.10
        addl      $7, %r10d                                     #295.10
        addl      $6, %edi                                      #295.10
        movl      %ecx, %r8d                                    #290.3
        movq      -16(%rsp), %rdx                               #[spill]
        cmpl      $5, %ecx                                      #290.3
        jb        ..B2.2        # Prob 79%                      #290.3
                                # LOE rdx rsi r8 r9 eax edi r10d r11d xmm0 xmm1 xmm11 xmm14
..B2.36:                        # Preds ..B2.35
                                # Execution count [1.00e+00]
        vmovsd    568(%r9), %xmm2                               #364.44
        vmulsd    %xmm2, %xmm2, %xmm0                           #364.44
        vmovsd    560(%r9), %xmm1                               #364.8
        vfmadd231sd %xmm1, %xmm1, %xmm0                         #364.3
        vcomisd   %xmm14, %xmm0                                 #366.13
        jb        ..B2.38       # Prob 50%                      #366.13
                                # LOE rbx rbp rsi r12 r13 r14 r15 eax xmm0 xmm1 xmm2 xmm11
..B2.37:                        # Preds ..B2.36
                                # Execution count [5.00e-01]
        vdivsd    %xmm0, %xmm11, %xmm11                         #367.16
        jmp       ..B2.39       # Prob 100%                     #367.16
                                # LOE rbx rbp rsi r12 r13 r14 r15 eax xmm1 xmm2 xmm11
..B2.38:                        # Preds ..B2.36
                                # Execution count [5.00e-01]
        movl      $1, %eax                                      #369.5
                                # LOE rbx rbp rsi r12 r13 r14 r15 eax xmm1 xmm2 xmm11
..B2.39:                        # Preds ..B2.37 ..B2.38
                                # Execution count [1.00e+00]
        vmulsd    %xmm11, %xmm1, %xmm0                          #373.19
        vmulsd    %xmm2, %xmm11, %xmm1                          #374.3
        vxorpd    .L_2il0floatpacket.37(%rip), %xmm1, %xmm2     #374.3
        vmovsd    %xmm0, 80(%rsi)                               #373.3
        vmovsd    %xmm2, 88(%rsi)                               #374.3
	.cfi_restore 6
        popq      %rbp                                          #376.10
	.cfi_def_cfa_offset 48
	.cfi_restore 3
        popq      %rbx                                          #376.10
	.cfi_def_cfa_offset 40
	.cfi_restore 15
        popq      %r15                                          #376.10
	.cfi_def_cfa_offset 32
	.cfi_restore 14
        popq      %r14                                          #376.10
	.cfi_def_cfa_offset 24
	.cfi_restore 13
        popq      %r13                                          #376.10
	.cfi_def_cfa_offset 16
	.cfi_restore 12
        popq      %r12                                          #376.10
	.cfi_def_cfa_offset 8
        ret                                                     #376.10
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	fwd_house_avx512,@function
	.size	fwd_house_avx512,.-fwd_house_avx512
	.section .rodata, "a"
	.align 64
	.align 8
.2.10_2.switchtab.4:
	.quad	..1.10_0.TAG.0
	.quad	..1.10_0.TAG.1
	.quad	..1.10_0.TAG.2
	.quad	..1.10_0.TAG.3
	.quad	..1.10_0.TAG.4
	.quad	..1.10_0.TAG.5
	.quad	..1.10_0.TAG.6
	.data
# -- End  fwd_house_avx512
	.text
# -- Begin  solv_sys_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl solv_sys_avx512
# --- solv_sys_avx512(complex_dble *, complex_dble *)
solv_sys_avx512:
# parameter 1: %rdi
# parameter 2: %rsi
..B3.1:                         # Preds ..B3.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_solv_sys_avx512.35:
..L36:
                                                         #381.1
        pushq     %r12                                          #381.1
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
        pushq     %r13                                          #381.1
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
        pushq     %r14                                          #381.1
	.cfi_def_cfa_offset 32
	.cfi_offset 14, -32
        pushq     %r15                                          #381.1
	.cfi_def_cfa_offset 40
	.cfi_offset 15, -40
        pushq     %rbx                                          #381.1
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
        pushq     %rbp                                          #381.1
	.cfi_def_cfa_offset 56
	.cfi_offset 6, -56
        movl      $5, %edx                                      #386.8
        vxorpd    %xmm0, %xmm0, %xmm0                           #410.24
        movl      $80, %eax                                     #386.8
                                # LOE rax rdx rsi rdi xmm0
..B3.2:                         # Preds ..B3.10 ..B3.1
                                # Execution count [5.00e+00]
        lea       -1(%rdx), %r13d                               #387.19
        movslq    %r13d, %r14                                   #387.10
        lea       -3(%rdx,%rdx,2), %ebp                         #387.10
        movq      %r14, %r12                                    #405.28
        addl      %ebp, %ebp                                    #387.10
        shlq      $4, %r12                                      #405.28
        movslq    %ebp, %rbp                                    #387.10
        addq      %rsi, %r12                                    #381.1
        shlq      $4, %rbp                                      #388.28
        testl     %r13d, %r13d                                  #387.28
        js        ..B3.10       # Prob 2%                       #387.28
                                # LOE rax rdx rbp rsi rdi r12 r14 r13d xmm0
..B3.3:                         # Preds ..B3.2
                                # Execution count [4.90e+00]
        lea       -1(%rdx), %r11                                #395.21
        movq      %r11, %rbx                                    #395.12
        lea       (%rdi,%rax), %r8                              #388.28
        shlq      $4, %rbx                                      #395.12
        lea       (%rbp,%r8), %r9                               #388.28
                                # LOE rax rdx rbx rbp rsi rdi r8 r9 r11 r12 r14 r13d xmm0
..B3.4:                         # Preds ..B3.8 ..B3.3
                                # Execution count [2.72e+01]
        vmovupd   (%rax,%rsi), %xmm2                            #390.25
        movq      %r11, %rcx                                    #395.12
        vmulpd    8(%r9){1to2}, %xmm2, %xmm1                    #391.12
        vpermilpd $1, %xmm1, %xmm1                              #392.12
        vfmaddsub231pd (%r9){1to2}, %xmm2, %xmm1                #393.12
        cmpq      %r14, %r11                                    #395.29
        jle       ..B3.8        # Prob 10%                      #395.29
                                # LOE rax rdx rcx rbx rbp rsi rdi r8 r9 r11 r12 r14 r13d xmm0 xmm1
..B3.5:                         # Preds ..B3.4
                                # Execution count [2.45e+01]
        lea       (%rdi,%rbp), %r10                             #396.30
        addq      %rbx, %r10                                    #396.30
                                # LOE rax rdx rcx rbx rbp rsi rdi r8 r9 r10 r11 r12 r14 r13d xmm0 xmm1
..B3.6:                         # Preds ..B3.6 ..B3.5
                                # Execution count [1.36e+02]
        lea       (%rcx,%rcx,2), %r15d                          #398.34
        addl      %r15d, %r15d                                  #398.34
        decq      %rcx                                          #395.32
        movslq    %r15d, %r15                                   #398.27
        shlq      $4, %r15                                      #398.27
        vmovupd   (%r8,%r15), %xmm3                             #398.27
        vmulpd    8(%r10){1to2}, %xmm3, %xmm2                   #399.14
        vpermilpd $1, %xmm2, %xmm4                              #400.14
        vfmaddsub231pd (%r10){1to2}, %xmm3, %xmm4               #401.14
        addq      $-16, %r10                                    #395.32
        vaddpd    %xmm4, %xmm1, %xmm1                           #402.14
        cmpq      %r14, %rcx                                    #395.29
        jg        ..B3.6        # Prob 82%                      #395.29
                                # LOE rax rdx rcx rbx rbp rsi rdi r8 r9 r10 r11 r12 r14 r13d xmm0 xmm1
..B3.8:                         # Preds ..B3.6 ..B3.4
                                # Execution count [2.72e+01]
        vmulpd    8(%r12){1to2}, %xmm1, %xmm2                   #407.12
        vpermilpd $1, %xmm2, %xmm3                              #408.12
        addq      $-96, %rbp                                    #387.31
        vfmaddsub231pd (%r12){1to2}, %xmm1, %xmm3               #409.12
        decq      %r14                                          #387.31
        vsubpd    %xmm3, %xmm0, %xmm1                           #410.12
        vmovupd   %xmm1, (%r9)                                  #411.23
        addq      $-96, %r9                                     #387.31
        addq      $-16, %r12                                    #387.31
        decl      %r13d                                         #387.31
        jns       ..B3.4        # Prob 82%                      #387.28
                                # LOE rax rdx rbx rbp rsi rdi r8 r9 r11 r12 r14 r13d xmm0
..B3.10:                        # Preds ..B3.8 ..B3.2
                                # Execution count [5.00e+00]
        .byte     15                                            #386.22
        .byte     31                                            #386.22
        .byte     128                                           #386.22
        .byte     0                                             #386.22
        .byte     0                                             #386.22
        .byte     0                                             #386.22
        .byte     0                                             #386.22
        addq      $-16, %rax                                    #386.22
        decq      %rdx                                          #386.22
        jg        ..B3.2        # Prob 80%                      #386.19
                                # LOE rax rdx rsi rdi xmm0
..B3.11:                        # Preds ..B3.10
                                # Execution count [1.00e+00]
	.cfi_restore 6
        popq      %rbp                                          #414.1
	.cfi_def_cfa_offset 48
	.cfi_restore 3
        popq      %rbx                                          #414.1
	.cfi_def_cfa_offset 40
	.cfi_restore 15
        popq      %r15                                          #414.1
	.cfi_def_cfa_offset 32
	.cfi_restore 14
        popq      %r14                                          #414.1
	.cfi_def_cfa_offset 24
	.cfi_restore 13
        popq      %r13                                          #414.1
	.cfi_def_cfa_offset 16
	.cfi_restore 12
        popq      %r12                                          #414.1
	.cfi_def_cfa_offset 8
        ret                                                     #414.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	solv_sys_avx512,@function
	.size	solv_sys_avx512,.-solv_sys_avx512
	.data
# -- End  solv_sys_avx512
	.text
# -- Begin  bck_house_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl bck_house_avx512
# --- bck_house_avx512(complex_dble *, complex_dble *, double *)
bck_house_avx512:
# parameter 1: %rdi
# parameter 2: %rsi
# parameter 3: %rdx
..B4.1:                         # Preds ..B4.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_bck_house_avx512.62:
..L63:
                                                         #417.1
        pushq     %r12                                          #417.1
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
        pushq     %r13                                          #417.1
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
        pushq     %r14                                          #417.1
	.cfi_def_cfa_offset 32
	.cfi_offset 14, -32
        pushq     %r15                                          #417.1
	.cfi_def_cfa_offset 40
	.cfi_offset 15, -40
        pushq     %rbx                                          #417.1
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
        pushq     %rbp                                          #417.1
	.cfi_def_cfa_offset 56
	.cfi_offset 6, -56
        movq      %rsi, %r8                                     #417.1
        movq      %rdx, %r9                                     #417.1
        xorl      %edx, %edx                                    #424.3
        xorl      %esi, %esi                                    #424.3
        vxorpd    %xmm0, %xmm0, %xmm0                           #441.12
        movq      80(%r8), %rax                                 #421.15
        movq      88(%r8), %rcx                                 #422.15
        movq      %rax, 560(%rdi)                               #421.3
        xorl      %eax, %eax                                    #436.26
        movq      %rcx, 568(%rdi)                               #422.3
        xorl      %ecx, %ecx                                    #424.3
                                # LOE rax rdx rdi r8 r9 ecx esi xmm0
..B4.2:                         # Preds ..B4.15 ..B4.1
                                # Execution count [5.00e+00]
        movl      %edx, %r12d                                   #425.12
        movq      %r8, %r11                                     #425.12
        movq      %r12, %rbp                                    #425.12
        movslq    %esi, %r15                                    #427.16
        shlq      $4, %rbp                                      #425.12
        shlq      $4, %r15                                      #427.16
        subq      %rbp, %r11                                    #425.12
        movq      448(%rdi,%r15), %r13                          #427.16
        movq      64(%r11), %rbx                                #425.12
        movq      %r13, 64(%r11)                                #427.5
        lea       1(%rdx), %r13d                                #432.5
        movq      %rbx, 448(%rdi,%r15)                          #429.5
        lea       5(%rcx), %ebx                                 #432.10
        movq      456(%rdi,%r15), %r10                          #428.16
        movq      72(%r11), %r14                                #426.12
        movq      %r10, 72(%r11)                                #428.5
        movq      %r14, 456(%rdi,%r15)                          #430.5
        cmpl      $6, %ebx                                      #432.27
        jge       ..B4.9        # Prob 50%                      #432.27
                                # LOE rax rdx rbp rdi r8 r9 r11 r12 ecx esi r13d xmm0
..B4.3:                         # Preds ..B4.2
                                # Execution count [5.00e+00]
        xorl      %r14d, %r14d                                  #432.5
        lea       1(%rdx), %r15d                                #432.5
        shrl      $1, %r15d                                     #432.5
        movl      $1, %ebx                                      #432.5
        xorl      %r10d, %r10d                                  #433.7
        testl     %r15d, %r15d                                  #432.5
        jbe       ..B4.7        # Prob 9%                       #432.5
                                # LOE rax rdx rbp rdi r8 r9 r11 r12 ecx ebx esi r10d r13d r14d r15d xmm0
..B4.4:                         # Preds ..B4.3
                                # Execution count [4.50e+00]
        movq      %r8, -24(%rsp)                                #[spill]
        movq      %r9, -16(%rsp)                                #[spill]
        .align    16,0x90
                                # LOE rax rdx rbp rdi r11 r12 ecx esi r10d r13d r14d r15d xmm0
..B4.5:                         # Preds ..B4.5 ..B4.4
                                # Execution count [1.25e+01]
        lea       (%rsi,%r10), %ebx                             #434.18
        addl      $12, %r10d                                    #432.5
        movslq    %ebx, %rbx                                    #434.18
        lea       (%r14,%r14), %r8d                             #434.7
        movslq    %r8d, %r8                                     #434.7
        incl      %r14d                                         #432.5
        shlq      $4, %rbx                                      #434.18
        shlq      $4, %r8                                       #434.7
        movq      552(%rdi,%rbx), %r9                           #434.18
        movq      %r9, 88(%r11,%r8)                             #434.7
        movq      544(%rdi,%rbx), %r9                           #433.18
        movq      %r9, 80(%r11,%r8)                             #433.7
        movq      648(%rdi,%rbx), %r9                           #434.18
        movq      %rax, 552(%rdi,%rbx)                          #436.7
        movq      %r9, 104(%r11,%r8)                            #434.7
        movq      640(%rdi,%rbx), %r9                           #433.18
        movq      %rax, 544(%rdi,%rbx)                          #435.7
        movq      %r9, 96(%r11,%r8)                             #433.7
        movq      %rax, 648(%rdi,%rbx)                          #436.7
        movq      %rax, 640(%rdi,%rbx)                          #435.7
        cmpl      %r15d, %r14d                                  #432.5
        jb        ..B4.5        # Prob 63%                      #432.5
                                # LOE rax rdx rbp rdi r11 r12 ecx esi r10d r13d r14d r15d xmm0
..B4.6:                         # Preds ..B4.5
                                # Execution count [4.50e+00]
        movq      -24(%rsp), %r8                                #[spill]
        lea       1(%r14,%r14), %ebx                            #433.7
        movq      -16(%rsp), %r9                                #[spill]
                                # LOE rax rdx rbp rdi r8 r9 r11 r12 ecx ebx esi r13d xmm0
..B4.7:                         # Preds ..B4.6 ..B4.3
                                # Execution count [5.00e+00]
        lea       -1(%rbx), %r10d                               #432.5
        cmpl      %r13d, %r10d                                  #432.5
        jae       ..B4.9        # Prob 9%                       #432.5
                                # LOE rax rdx rbp rdi r8 r9 r11 r12 ecx ebx esi r13d xmm0
..B4.8:                         # Preds ..B4.7
                                # Execution count [4.50e+00]
        movslq    %ebx, %r10                                    #434.7
        lea       (%rbx,%rbx,2), %ebx                           #434.18
        subq      %r12, %r10                                    #434.7
        lea       (%rsi,%rbx,2), %r14d                          #434.18
        movslq    %r14d, %r14                                   #434.18
        shlq      $4, %r14                                      #434.18
        shlq      $4, %r10                                      #434.7
        movq      456(%rdi,%r14), %r15                          #434.18
        movq      %r15, 72(%r8,%r10)                            #434.7
        movq      448(%rdi,%r14), %r15                          #433.18
        movq      %r15, 64(%r8,%r10)                            #433.7
        movq      %rax, 456(%rdi,%r14)                          #436.7
        movq      %rax, 448(%rdi,%r14)                          #435.7
                                # LOE rax rdx rbp rdi r8 r9 r11 r12 ecx esi r13d xmm0
..B4.9:                         # Preds ..B4.2 ..B4.8 ..B4.7
                                # Execution count [3.96e-01]
        shlq      $3, %r12                                      #453.28
        negq      %rbp                                          #444.30
        negq      %r12                                          #453.28
        addq      %r9, %r12                                     #453.28
        addq      %rdi, %rbp                                    #444.30
        addq      $2, %rdx                                      #443.23
        xorb      %bl, %bl                                      #439.5
                                # LOE rax rdx rbp rdi r8 r9 r11 r12 ecx esi r13d bl xmm0
..B4.10:                        # Preds ..B4.14 ..B4.9
                                # Execution count [2.54e+01]
        movq      %rax, %r14                                    #443.7
        vmovapd   %xmm0, %xmm1                                  #441.12
        movq      %r14, %r10                                    #443.7
                                # LOE rax rdx rbp rdi r8 r9 r10 r11 r12 r14 ecx esi r13d bl xmm0 xmm1
..B4.11:                        # Preds ..B4.11 ..B4.10
                                # Execution count [1.41e+02]
        vmovupd   64(%r10,%r11), %xmm3                          #446.27
        incq      %r14                                          #443.7
        vmulpd    72(%r10,%rbp){1to2}, %xmm3, %xmm2             #447.14
        vpermilpd $1, %xmm2, %xmm4                              #448.14
        vfmaddsub231pd 64(%r10,%rbp){1to2}, %xmm3, %xmm4        #449.14
        addq      $16, %r10                                     #443.7
        vaddpd    %xmm1, %xmm4, %xmm1                           #450.14
        cmpq      %rdx, %r14                                    #443.7
        jb        ..B4.11       # Prob 82%                      #443.7
                                # LOE rax rdx rbp rdi r8 r9 r10 r11 r12 r14 ecx esi r13d bl xmm0 xmm1
..B4.12:                        # Preds ..B4.11
                                # Execution count [2.70e+01]
        movq      %rax, %r15                                    #456.7
        lea       64(%rbp), %r10                                #456.7
        vmulpd    32(%r12){1to2}, %xmm1, %xmm1                  #454.12
        movq      %r15, %r14                                    #456.7
                                # LOE rax rdx rbp rdi r8 r9 r10 r11 r12 r14 r15 ecx esi r13d bl xmm0 xmm1
..B4.13:                        # Preds ..B4.13 ..B4.12
                                # Execution count [1.50e+02]
        vmulpd    72(%r14,%r11){1to2}, %xmm1, %xmm2             #459.14
        vpermilpd $1, %xmm2, %xmm4                              #460.14
        incq      %r15                                          #456.7
        vfmsubadd231pd 64(%r14,%r11){1to2}, %xmm1, %xmm4        #461.14
        addq      $16, %r14                                     #456.7
        vmovupd   (%r10), %xmm3                                 #463.28
        vsubpd    %xmm4, %xmm3, %xmm5                           #464.14
        vmovupd   %xmm5, (%r10)                                 #465.25
        addq      $16, %r10                                     #456.7
        cmpq      %rdx, %r15                                    #456.7
        jb        ..B4.13       # Prob 82%                      #456.7
                                # LOE rax rdx rbp rdi r8 r9 r10 r11 r12 r14 r15 ecx esi r13d bl xmm0 xmm1
..B4.14:                        # Preds ..B4.13
                                # Execution count [3.00e+01]
        incb      %bl                                           #439.5
        addq      $96, %rbp                                     #439.5
        cmpb      $6, %bl                                       #439.5
        jb        ..B4.10       # Prob 83%                      #439.5
                                # LOE rax rdx rbp rdi r8 r9 r11 r12 ecx esi r13d bl xmm0
..B4.15:                        # Preds ..B4.14
                                # Execution count [5.00e+00]
        addl      $-7, %esi                                     #432.5
        decl      %ecx                                          #432.5
        movl      %r13d, %edx                                   #424.3
        cmpl      $5, %r13d                                     #424.3
        jb        ..B4.2        # Prob 79%                      #424.3
                                # LOE rax rdx rdi r8 r9 ecx esi xmm0
..B4.16:                        # Preds ..B4.15
                                # Execution count [1.00e+00]
	.cfi_restore 6
        popq      %rbp                                          #469.1
	.cfi_def_cfa_offset 48
	.cfi_restore 3
        popq      %rbx                                          #469.1
	.cfi_def_cfa_offset 40
	.cfi_restore 15
        popq      %r15                                          #469.1
	.cfi_def_cfa_offset 32
	.cfi_restore 14
        popq      %r14                                          #469.1
	.cfi_def_cfa_offset 24
	.cfi_restore 13
        popq      %r13                                          #469.1
	.cfi_def_cfa_offset 16
	.cfi_restore 12
        popq      %r12                                          #469.1
	.cfi_def_cfa_offset 8
        ret                                                     #469.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	bck_house_avx512,@function
	.size	bck_house_avx512,.-bck_house_avx512
	.data
# -- End  bck_house_avx512
	.section .rodata, "a"
	.space 8, 0x00 	# pad
	.align 64
.L_2il0floatpacket.9:
	.long	0x00000000,0x00000000,0x00000001,0x00000000,0x0000000c,0x00000000,0x0000000d,0x00000000,0x00000000,0x00000000,0x00000001,0x00000000,0x0000000c,0x00000000,0x0000000d,0x00000000
	.type	.L_2il0floatpacket.9,@object
	.size	.L_2il0floatpacket.9,64
	.align 64
.L_2il0floatpacket.10:
	.long	0x00000002,0x00000000,0x00000003,0x00000000,0x0000000e,0x00000000,0x0000000f,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000,0x0000000e,0x00000000,0x0000000f,0x00000000
	.type	.L_2il0floatpacket.10,@object
	.size	.L_2il0floatpacket.10,64
	.align 64
.L_2il0floatpacket.11:
	.long	0x00000004,0x00000000,0x00000005,0x00000000,0x00000008,0x00000000,0x00000009,0x00000000,0x00000004,0x00000000,0x00000005,0x00000000,0x00000008,0x00000000,0x00000009,0x00000000
	.type	.L_2il0floatpacket.11,@object
	.size	.L_2il0floatpacket.11,64
	.align 64
.L_2il0floatpacket.12:
	.long	0x00000006,0x00000000,0x00000007,0x00000000,0x0000000a,0x00000000,0x0000000b,0x00000000,0x00000006,0x00000000,0x00000007,0x00000000,0x0000000a,0x00000000,0x0000000b,0x00000000
	.type	.L_2il0floatpacket.12,@object
	.size	.L_2il0floatpacket.12,64
	.align 64
.L_2il0floatpacket.13:
	.long	0x00000000,0x00000000,0x00000001,0x00000000,0x00000008,0x00000000,0x00000009,0x00000000,0x00000006,0x00000000,0x00000007,0x00000000,0x0000000e,0x00000000,0x0000000f,0x00000000
	.type	.L_2il0floatpacket.13,@object
	.size	.L_2il0floatpacket.13,64
	.align 64
.L_2il0floatpacket.14:
	.long	0x00000004,0x00000000,0x00000005,0x00000000,0x0000000c,0x00000000,0x0000000d,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000,0x0000000a,0x00000000,0x0000000b,0x00000000
	.type	.L_2il0floatpacket.14,@object
	.size	.L_2il0floatpacket.14,64
	.align 64
.L_2il0floatpacket.15:
	.long	0x00000000,0x00000000,0x00000000,0x00000000,0x00000002,0x00000000,0x00000002,0x00000000,0x0000000c,0x00000000,0x0000000c,0x00000000,0x0000000e,0x00000000,0x0000000e,0x00000000
	.type	.L_2il0floatpacket.15,@object
	.size	.L_2il0floatpacket.15,64
	.align 64
.L_2il0floatpacket.16:
	.long	0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x0000000d,0x00000000,0x0000000d,0x00000000,0x0000000f,0x00000000,0x0000000f,0x00000000
	.type	.L_2il0floatpacket.16,@object
	.size	.L_2il0floatpacket.16,64
	.align 64
.L_2il0floatpacket.17:
	.long	0x00000002,0x00000000,0x00000003,0x00000000,0x0000000a,0x00000000,0x0000000b,0x00000000,0x00000004,0x00000000,0x00000005,0x00000000,0x0000000c,0x00000000,0x0000000d,0x00000000
	.type	.L_2il0floatpacket.17,@object
	.size	.L_2il0floatpacket.17,64
	.align 64
.L_2il0floatpacket.18:
	.long	0x00000004,0x00000000,0x00000004,0x00000000,0x00000006,0x00000000,0x00000006,0x00000000,0x00000009,0x00000000,0x00000009,0x00000000,0x0000000b,0x00000000,0x0000000b,0x00000000
	.type	.L_2il0floatpacket.18,@object
	.size	.L_2il0floatpacket.18,64
	.align 64
.L_2il0floatpacket.19:
	.long	0x00000005,0x00000000,0x00000005,0x00000000,0x00000007,0x00000000,0x00000007,0x00000000,0x00000008,0x00000000,0x00000008,0x00000000,0x00000008,0x00000000,0x00000008,0x00000000
	.type	.L_2il0floatpacket.19,@object
	.size	.L_2il0floatpacket.19,64
	.align 64
.L_2il0floatpacket.20:
	.long	0x00000001,0x00000000,0x00000001,0x00000000,0x00000003,0x00000000,0x00000003,0x00000000,0x0000000c,0x00000000,0x0000000c,0x00000000,0x0000000e,0x00000000,0x0000000e,0x00000000
	.type	.L_2il0floatpacket.20,@object
	.size	.L_2il0floatpacket.20,64
	.align 64
.L_2il0floatpacket.21:
	.long	0x00000004,0x00000000,0x00000004,0x00000000,0x00000006,0x00000000,0x00000006,0x00000000,0x0000000c,0x00000000,0x0000000c,0x00000000,0x0000000e,0x00000000,0x0000000e,0x00000000
	.type	.L_2il0floatpacket.21,@object
	.size	.L_2il0floatpacket.21,64
	.align 64
.L_2il0floatpacket.22:
	.long	0x00000004,0x00000000,0x00000005,0x00000000,0x0000000c,0x00000000,0x0000000d,0x00000000,0x00000006,0x00000000,0x00000007,0x00000000,0x0000000e,0x00000000,0x0000000f,0x00000000
	.type	.L_2il0floatpacket.22,@object
	.size	.L_2il0floatpacket.22,64
	.align 64
.L_2il0floatpacket.23:
	.long	0x00000000,0x00000000,0x00000000,0x00000000,0x00000002,0x00000000,0x00000002,0x00000000,0x00000008,0x00000000,0x00000008,0x00000000,0x0000000a,0x00000000,0x0000000a,0x00000000
	.type	.L_2il0floatpacket.23,@object
	.size	.L_2il0floatpacket.23,64
	.align 64
.L_2il0floatpacket.24:
	.long	0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000009,0x00000000,0x00000009,0x00000000,0x0000000b,0x00000000,0x0000000b,0x00000000
	.type	.L_2il0floatpacket.24,@object
	.size	.L_2il0floatpacket.24,64
	.align 64
.L_2il0floatpacket.25:
	.long	0x00000000,0x00000000,0x00000000,0x00000000,0x00000002,0x00000000,0x00000002,0x00000000,0x0000000d,0x00000000,0x0000000d,0x00000000,0x0000000f,0x00000000,0x0000000f,0x00000000
	.type	.L_2il0floatpacket.25,@object
	.size	.L_2il0floatpacket.25,64
	.align 64
.L_2il0floatpacket.26:
	.long	0x00000001,0x00000000,0x00000001,0x00000000,0x00000003,0x00000000,0x00000003,0x00000000,0x00000008,0x00000000,0x00000008,0x00000000,0x00000008,0x00000000,0x00000008,0x00000000
	.type	.L_2il0floatpacket.26,@object
	.size	.L_2il0floatpacket.26,64
	.align 64
.L_2il0floatpacket.27:
	.long	0x00000004,0x00000000,0x00000004,0x00000000,0x00000006,0x00000000,0x00000006,0x00000000,0x00000008,0x00000000,0x00000008,0x00000000,0x0000000a,0x00000000,0x0000000a,0x00000000
	.type	.L_2il0floatpacket.27,@object
	.size	.L_2il0floatpacket.27,64
	.align 64
.L_2il0floatpacket.28:
	.long	0x00000005,0x00000000,0x00000005,0x00000000,0x00000007,0x00000000,0x00000007,0x00000000,0x00000009,0x00000000,0x00000009,0x00000000,0x0000000b,0x00000000,0x0000000b,0x00000000
	.type	.L_2il0floatpacket.28,@object
	.size	.L_2il0floatpacket.28,64
	.align 64
.L_2il0floatpacket.29:
	.long	0x00000001,0x00000000,0x00000001,0x00000000,0x00000003,0x00000000,0x00000003,0x00000000,0x0000000d,0x00000000,0x0000000d,0x00000000,0x0000000f,0x00000000,0x0000000f,0x00000000
	.type	.L_2il0floatpacket.29,@object
	.size	.L_2il0floatpacket.29,64
	.align 64
.L_2il0floatpacket.30:
	.long	0x00000000,0x00000000,0x00000000,0x00000000,0x00000008,0x00000000,0x00000008,0x00000000,0x00000006,0x00000000,0x00000006,0x00000000,0x0000000e,0x00000000,0x0000000e,0x00000000
	.type	.L_2il0floatpacket.30,@object
	.size	.L_2il0floatpacket.30,64
	.align 64
.L_2il0floatpacket.31:
	.long	0x00000001,0x00000000,0x00000001,0x00000000,0x00000009,0x00000000,0x00000009,0x00000000,0x00000007,0x00000000,0x00000007,0x00000000,0x0000000f,0x00000000,0x0000000f,0x00000000
	.type	.L_2il0floatpacket.31,@object
	.size	.L_2il0floatpacket.31,64
	.align 64
.L_2il0floatpacket.32:
	.long	0x00000002,0x00000000,0x00000003,0x00000000,0x0000000a,0x00000000,0x0000000b,0x00000000,0x00000000,0x00000000,0x00000001,0x00000000,0x00000008,0x00000000,0x00000009,0x00000000
	.type	.L_2il0floatpacket.32,@object
	.size	.L_2il0floatpacket.32,64
	.align 64
.L_2il0floatpacket.33:
	.long	0x00000000,0x00000000,0x00000001,0x00000000,0x00000008,0x00000000,0x00000009,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000,0x0000000a,0x00000000,0x0000000b,0x00000000
	.type	.L_2il0floatpacket.33,@object
	.size	.L_2il0floatpacket.33,64
	.align 64
.L_2il0floatpacket.34:
	.long	0x00000006,0x00000000,0x00000007,0x00000000,0x0000000e,0x00000000,0x0000000f,0x00000000,0x00000000,0x00000000,0x00000001,0x00000000,0x00000008,0x00000000,0x00000009,0x00000000
	.type	.L_2il0floatpacket.34,@object
	.size	.L_2il0floatpacket.34,64
	.align 64
.L_2il0floatpacket.35:
	.long	0x00000002,0x00000000,0x00000003,0x00000000,0x0000000e,0x00000000,0x0000000f,0x00000000,0x00000000,0x00000000,0x00000001,0x00000000,0x0000000c,0x00000000,0x0000000d,0x00000000
	.type	.L_2il0floatpacket.35,@object
	.size	.L_2il0floatpacket.35,64
	.align 16
.L_2il0floatpacket.37:
	.long	0x00000000,0x80000000,0x00000000,0x00000000
	.type	.L_2il0floatpacket.37,@object
	.size	.L_2il0floatpacket.37,16
	.align 8
.L_2il0floatpacket.36:
	.long	0x00000000,0x3cb00000
	.type	.L_2il0floatpacket.36,@object
	.size	.L_2il0floatpacket.36,8
	.align 8
.L_2il0floatpacket.38:
	.long	0x00000000,0x3ff00000
	.type	.L_2il0floatpacket.38,@object
	.size	.L_2il0floatpacket.38,8
	.data
	.section .note.GNU-stack, ""
// -- Begin DWARF2 SEGMENT .eh_frame
	.section .eh_frame,"a",@progbits
.eh_frame_seg:
	.align 8
# End
