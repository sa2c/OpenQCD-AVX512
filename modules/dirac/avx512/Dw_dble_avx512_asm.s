# mark_description "Intel(R) C Intel(R) 64 Compiler for applications running on Intel(R) 64, Version 17.0.4.196 Build 20170411";
# mark_description "-I../../../include -I.. -I/cineca/prod/opt/compilers/intel/pe-xe-2017/binary/impi/2017.3.196/intel64/include";
# mark_description " -isystem /cineca/prod/opt/compilers/intel/pe-xe-2018/binary/impi/2018.1.163/include64/ -std=c89 -xCORE-AVX5";
# mark_description "12 -mtune=skylake -DAVX512 -O3 -Ddirac_counters -pedantic -fstrict-aliasing -Wno-long-long -Wstrict-prototyp";
# mark_description "es -S";
	.file "Dw_dble_avx512.c"
	.text
..TXTST0:
# -- Begin  doe_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl doe_dble_avx512
# --- doe_dble_avx512(const int *, const int *, const su3_dble *, const spinor_dble *, double, spin_t *)
doe_dble_avx512:
# parameter 1: %rdi
# parameter 2: %rsi
# parameter 3: %rdx
# parameter 4: %rcx
# parameter 5: %xmm0
# parameter 6: %r8
..B1.1:                         # Preds ..B1.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_doe_dble_avx512.1:
..L2:
                                                          #27.1
        pushq     %rbp                                          #27.1
	.cfi_def_cfa_offset 16
        movq      %rsp, %rbp                                    #27.1
	.cfi_def_cfa 6, 16
	.cfi_offset 6, -16
        movslq    (%rdi), %rax                                  #42.16
        movslq    (%rsi), %r9                                   #43.16
        vmovups   .L_2il0floatpacket.14(%rip), %zmm13           #45.3
        vmovups   .L_2il0floatpacket.15(%rip), %zmm15           #45.3
        vmovups   .L_2il0floatpacket.16(%rip), %zmm16           #45.3
        vmovups   .L_2il0floatpacket.17(%rip), %zmm19           #45.3
        vmovups   .L_2il0floatpacket.18(%rip), %zmm18           #45.3
        vmovsd    %xmm0, -16(%rbp)                              #27.1
        vmovaps   %zmm13, %zmm28                                #45.3
        lea       (%rax,%rax,2), %r11                           #42.8
        shlq      $6, %r11                                      #42.8
        lea       (%r9,%r9,2), %r10                             #43.8
        shlq      $6, %r10                                      #43.8
        movl      $15, %eax                                     #48.8
        vmovaps   %zmm13, %zmm27                                #46.3
        kmovw     %eax, %k4                                     #48.8
        movl      $240, %eax                                    #49.8
        kmovw     %eax, %k3                                     #49.8
        vmovups   (%rcx,%r10), %zmm30                           #45.3
        vmovups   (%rcx,%r11), %zmm24                           #45.3
        vmovups   96(%rcx,%r10), %zmm29                         #46.3
        vmovups   96(%rcx,%r11), %zmm23                         #46.3
        vpermi2pd %zmm24, %zmm30, %zmm28                        #45.3
        vpermt2pd 64(%rcx,%r11), %zmm16, %zmm24                 #45.3
        vpermt2pd 64(%rcx,%r10), %zmm15, %zmm30                 #45.3
        vpermi2pd %zmm23, %zmm29, %zmm27                        #46.3
        vpermt2pd 160(%rcx,%r10), %zmm15, %zmm29                #46.3
        vpermt2pd 160(%rcx,%r11), %zmm16, %zmm23                #46.3
        vaddpd    %zmm27, %zmm28, %zmm7{%k4}{z}                 #48.8
        movslq    4(%rdi), %rax                                 #55.16
        vmovaps   %zmm19, %zmm26                                #45.3
        vmovaps   %zmm19, %zmm25                                #46.3
        vpermi2pd %zmm30, %zmm24, %zmm26                        #45.3
        lea       (%rax,%rax,2), %rax                           #55.8
        vpermt2pd %zmm30, %zmm18, %zmm24                        #45.3
        vpermi2pd %zmm29, %zmm23, %zmm25                        #46.3
        vpermt2pd %zmm29, %zmm18, %zmm23                        #46.3
        vsubpd    %zmm27, %zmm28, %zmm7{%k3}                    #49.8
        vaddpd    %zmm25, %zmm26, %zmm4{%k4}{z}                 #50.8
        vaddpd    %zmm23, %zmm24, %zmm14{%k4}{z}                #52.8
        vsubpd    %zmm25, %zmm26, %zmm4{%k3}                    #51.8
        vsubpd    %zmm23, %zmm24, %zmm14{%k3}                   #53.8
        shlq      $6, %rax                                      #55.8
        prefetcht0 (%rcx,%rax)                                  #56.3
        movslq    4(%rsi), %r9                                  #57.16
        vpermilpd $85, %zmm7, %zmm22                            #62.3
        vpermilpd $85, %zmm4, %zmm10                            #62.3
        vpermilpd $85, %zmm14, %zmm29                           #62.3
        lea       (%r9,%r9,2), %r10                             #57.8
        movl      $90, %r9d                                     #74.8
        kmovw     %r9d, %k1                                     #74.8
        movl      $165, %r9d                                    #75.8
        kmovw     %r9d, %k2                                     #75.8
        shlq      $6, %r10                                      #57.8
        movl      $175, %r9d                                    #92.3
        kmovw     %r9d, %k5                                     #92.3
        movl      $80, %r9d                                     #92.3
        kmovw     %r9d, %k6                                     #92.3
        movl      $60, %r9d                                     #102.8
        kmovw     %r9d, %k7                                     #102.8
        prefetcht0 (%rcx,%r10)                                  #58.3
        vmovups   .L_2il0floatpacket.19(%rip), %zmm30           #62.3
        vmovups   (%rdx), %zmm11                                #62.3
        vmovups   .L_2il0floatpacket.25(%rip), %zmm24           #62.3
        vmovups   64(%rdx), %zmm17                              #62.3
        vmovups   128(%rdx), %zmm8                              #62.3
        vmulpd    %zmm10, %zmm30, %zmm21                        #62.3
        vmulpd    %zmm29, %zmm30, %zmm10                        #62.3
        vmulpd    %zmm30, %zmm22, %zmm6                         #62.3
        vmovups   .L_2il0floatpacket.20(%rip), %zmm29           #62.3
        vmovups   .L_2il0floatpacket.27(%rip), %zmm22           #62.3
        vmovaps   %zmm11, %zmm28                                #62.3
        vpermt2pd 144(%rdx), %zmm29, %zmm28                     #62.3
        vmulpd    %zmm7, %zmm28, %zmm27                         #62.3
        vmovups   .L_2il0floatpacket.21(%rip), %zmm28           #62.3
        vmovaps   %zmm11, %zmm26                                #62.3
        vpermt2pd 144(%rdx), %zmm28, %zmm26                     #62.3
        vfmadd213pd %zmm27, %zmm6, %zmm26                       #62.3
        vmovups   .L_2il0floatpacket.22(%rip), %zmm27           #62.3
        vmovaps   %zmm11, %zmm25                                #62.3
        vpermt2pd 144(%rdx), %zmm27, %zmm25                     #62.3
        vfmadd213pd %zmm26, %zmm4, %zmm25                       #62.3
        vmovups   .L_2il0floatpacket.23(%rip), %zmm26           #62.3
        vmovaps   %zmm11, %zmm9                                 #62.3
        vpermt2pd 144(%rdx), %zmm26, %zmm9                      #62.3
        vfmadd213pd %zmm25, %zmm21, %zmm9                       #62.3
        vmovups   .L_2il0floatpacket.24(%rip), %zmm25           #62.3
        vmovaps   %zmm11, %zmm23                                #62.3
        vpermt2pd 208(%rdx), %zmm25, %zmm23                     #62.3
        vfmadd213pd %zmm9, %zmm14, %zmm23                       #62.3
        vmovaps   %zmm11, %zmm9                                 #62.3
        vpermt2pd 208(%rdx), %zmm24, %zmm9                      #62.3
        vfmadd213pd %zmm23, %zmm10, %zmm9                       #62.3
        vmovups   .L_2il0floatpacket.26(%rip), %zmm23           #62.3
        vmovaps   %zmm17, %zmm1                                 #62.3
        vmovaps   %zmm11, %zmm3                                 #62.3
        vpermt2pd 144(%rdx), %zmm25, %zmm1                      #62.3
        vpermt2pd 144(%rdx), %zmm23, %zmm3                      #62.3
        vpermt2pd 144(%rdx), %zmm22, %zmm11                     #62.3
        vmulpd    %zmm1, %zmm7, %zmm31                          #62.3
        vmulpd    %zmm3, %zmm7, %zmm12                          #62.3
        vmovups   96(%rcx,%r10), %zmm3                          #71.3
        vfmadd213pd %zmm12, %zmm6, %zmm11                       #62.3
        vmovaps   %zmm17, %zmm7                                 #62.3
        vpermt2pd 144(%rdx), %zmm24, %zmm7                      #62.3
        vfmadd213pd %zmm31, %zmm6, %zmm7                        #62.3
        vmovaps   %zmm17, %zmm6                                 #62.3
        vmovaps   %zmm17, %zmm5                                 #62.3
        vpermt2pd 208(%rdx), %zmm23, %zmm6                      #62.3
        vpermt2pd 208(%rdx), %zmm29, %zmm5                      #62.3
        vfmadd213pd %zmm7, %zmm4, %zmm6                         #62.3
        vfmadd213pd %zmm11, %zmm4, %zmm5                        #62.3
        vmovups   96(%rcx,%rax), %zmm7                          #71.3
        vmovaps   %zmm17, %zmm2                                 #62.3
        vmovaps   %zmm17, %zmm20                                #62.3
        vmovaps   %zmm17, %zmm11                                #62.3
        vpermt2pd 208(%rdx), %zmm22, %zmm17                     #62.3
        vpermt2pd 208(%rdx), %zmm28, %zmm2                      #62.3
        vpermt2pd 208(%rdx), %zmm27, %zmm20                     #62.3
        vpermt2pd 208(%rdx), %zmm26, %zmm11                     #62.3
        vfmadd213pd %zmm6, %zmm21, %zmm17                       #62.3
        vfmadd213pd %zmm5, %zmm21, %zmm2                        #62.3
        vmovaps   %zmm8, %zmm21                                 #62.3
        vpermt2pd 272(%rdx), %zmm29, %zmm21                     #62.3
        vpermt2pd 272(%rdx), %zmm28, %zmm8                      #62.3
        vfmadd213pd %zmm2, %zmm14, %zmm20                       #62.3
        vfmadd213pd %zmm17, %zmm14, %zmm21                      #62.3
        vfmadd213pd %zmm20, %zmm10, %zmm11                      #62.3
        vfmadd213pd %zmm21, %zmm10, %zmm8                       #62.3
        vmovups   .L_2il0floatpacket.28(%rip), %zmm21           #64.3
        vmovups   .L_2il0floatpacket.31(%rip), %zmm20           #71.3
        vpermpd   %zmm9, %zmm21, %zmm17                         #64.3
        vpermpd   %zmm11, %zmm21, %zmm14                        #65.3
        vpermpd   %zmm8, %zmm21, %zmm4                          #66.3
        vaddpd    %zmm9, %zmm17, %zmm10{%k4}{z}                 #64.3
        vsubpd    %zmm9, %zmm17, %zmm10{%k3}                    #64.3
        vaddpd    %zmm11, %zmm14, %zmm9{%k4}{z}                 #65.3
        vmovups   .L_2il0floatpacket.29(%rip), %zmm17           #71.3
        vsubpd    %zmm11, %zmm14, %zmm9{%k3}                    #65.3
        vaddpd    %zmm8, %zmm4, %zmm11{%k4}{z}                  #66.3
        vmovups   .L_2il0floatpacket.30(%rip), %zmm14           #71.3
        vsubpd    %zmm8, %zmm4, %zmm11{%k3}                     #66.3
        vmovups   (%rcx,%r10), %zmm8                            #70.3
        vmovups   (%rcx,%rax), %zmm4                            #70.3
        vmovaps   %zmm3, %zmm12                                 #71.3
        vmovaps   %zmm13, %zmm5                                 #70.3
        vpermt2pd %zmm7, %zmm17, %zmm12                         #71.3
        vpermt2pd 160(%rcx,%r10), %zmm14, %zmm3                 #71.3
        vpermt2pd 160(%rcx,%rax), %zmm20, %zmm7                 #71.3
        vpermi2pd %zmm4, %zmm8, %zmm5                           #70.3
        vpermt2pd 64(%rcx,%r10), %zmm15, %zmm8                  #70.3
        vpermt2pd 64(%rcx,%rax), %zmm16, %zmm4                  #70.3
        vmovaps   %zmm19, %zmm0                                 #71.3
        vmovaps   %zmm19, %zmm6                                 #70.3
        vpermi2pd %zmm3, %zmm7, %zmm0                           #71.3
        vpermi2pd %zmm8, %zmm4, %zmm6                           #70.3
        vpermt2pd %zmm8, %zmm18, %zmm4                          #70.3
        vpermt2pd %zmm3, %zmm18, %zmm7                          #71.3
        vpermilpd $85, %zmm12, %zmm1                            #73.8
        vaddpd    %zmm1, %zmm5, %zmm2{%k1}{z}                   #74.8
        vpermilpd $85, %zmm0, %zmm8                             #76.8
        movslq    8(%rdi), %r11                                 #83.16
        vsubpd    %zmm1, %zmm5, %zmm2{%k2}                      #75.8
        vaddpd    %zmm8, %zmm6, %zmm5{%k1}{z}                   #77.8
        vsubpd    %zmm8, %zmm6, %zmm5{%k2}                      #78.8
        lea       (%r11,%r11,2), %r10                           #83.8
        vpermilpd $85, %zmm7, %zmm6                             #79.8
        shlq      $6, %r10                                      #83.8
        vaddpd    %zmm6, %zmm4, %zmm8{%k1}{z}                   #80.8
        vsubpd    %zmm6, %zmm4, %zmm8{%k2}                      #81.8
        prefetcht0 (%rcx,%r10)                                  #84.3
        movslq    8(%rsi), %rax                                 #85.16
        vpermilpd $85, %zmm2, %zmm4                             #90.3
        vpermilpd $85, %zmm5, %zmm12                            #90.3
        vpermilpd $85, %zmm8, %zmm1                             #90.3
        lea       (%rax,%rax,2), %r9                            #85.8
        shlq      $6, %r9                                       #85.8
        movl      $63, %eax                                     #117.3
        kmovw     %eax, %k1                                     #117.3
        movl      $192, %eax                                    #117.3
        kmovw     %eax, %k2                                     #117.3
        vmulpd    %zmm30, %zmm4, %zmm3                          #90.3
        vmulpd    %zmm12, %zmm30, %zmm4                         #90.3
        vmulpd    %zmm1, %zmm30, %zmm7                          #90.3
        prefetcht0 (%rcx,%r9)                                   #86.3
        movl      $195, %eax                                    #101.8
        vmovups   288(%rdx), %zmm1                              #90.3
        vmovups   352(%rdx), %zmm12                             #90.3
        vmovups   416(%rdx), %zmm6                              #90.3
        vmovaps   %zmm1, %zmm31                                 #90.3
        vpermt2pd 432(%rdx), %zmm29, %zmm31                     #90.3
        vmulpd    %zmm2, %zmm31, %zmm0                          #90.3
        vmovaps   %zmm1, %zmm31                                 #90.3
        vpermt2pd 432(%rdx), %zmm28, %zmm31                     #90.3
        vfmadd213pd %zmm0, %zmm3, %zmm31                        #90.3
        vmovaps   %zmm1, %zmm0                                  #90.3
        vpermt2pd 432(%rdx), %zmm27, %zmm0                      #90.3
        vfmadd213pd %zmm31, %zmm5, %zmm0                        #90.3
        vmovaps   %zmm1, %zmm31                                 #90.3
        vpermt2pd 432(%rdx), %zmm26, %zmm31                     #90.3
        vfmadd213pd %zmm0, %zmm4, %zmm31                        #90.3
        vmovaps   %zmm1, %zmm0                                  #90.3
        vpermt2pd 496(%rdx), %zmm25, %zmm0                      #90.3
        vfmadd213pd %zmm31, %zmm8, %zmm0                        #90.3
        vmovaps   %zmm1, %zmm31                                 #90.3
        vpermt2pd 496(%rdx), %zmm24, %zmm31                     #90.3
        vfmadd213pd %zmm0, %zmm7, %zmm31                        #90.3
        vmovaps   %zmm1, %zmm0                                  #90.3
        vpermt2pd 432(%rdx), %zmm23, %zmm0                      #90.3
        vpermt2pd 432(%rdx), %zmm22, %zmm1                      #90.3
        vmulpd    %zmm0, %zmm2, %zmm0                           #90.3
        vfmadd213pd %zmm0, %zmm3, %zmm1                         #90.3
        vmovaps   %zmm12, %zmm0                                 #90.3
        vpermt2pd 496(%rdx), %zmm29, %zmm0                      #90.3
        vfmadd213pd %zmm1, %zmm5, %zmm0                         #90.3
        vmovaps   %zmm12, %zmm1                                 #90.3
        vpermt2pd 496(%rdx), %zmm28, %zmm1                      #90.3
        vfmadd213pd %zmm0, %zmm4, %zmm1                         #90.3
        vmovaps   %zmm12, %zmm0                                 #90.3
        vpermt2pd 496(%rdx), %zmm27, %zmm0                      #90.3
        vfmadd213pd %zmm1, %zmm8, %zmm0                         #90.3
        vmovaps   %zmm12, %zmm1                                 #90.3
        vpermt2pd 496(%rdx), %zmm26, %zmm1                      #90.3
        vfmadd213pd %zmm0, %zmm7, %zmm1                         #90.3
        vmovaps   %zmm12, %zmm0                                 #90.3
        vpermt2pd 432(%rdx), %zmm25, %zmm0                      #90.3
        vmulpd    %zmm0, %zmm2, %zmm2                           #90.3
        vmovaps   %zmm12, %zmm0                                 #90.3
        vpermt2pd 432(%rdx), %zmm24, %zmm0                      #90.3
        vfmadd213pd %zmm2, %zmm3, %zmm0                         #90.3
        vmovups   .L_2il0floatpacket.32(%rip), %zmm2            #92.3
        vmovaps   %zmm12, %zmm3                                 #90.3
        vpermt2pd 496(%rdx), %zmm23, %zmm3                      #90.3
        vpermt2pd 496(%rdx), %zmm22, %zmm12                     #90.3
        vfmadd213pd %zmm0, %zmm5, %zmm3                         #90.3
        vfmadd213pd %zmm3, %zmm4, %zmm12                        #90.3
        vpermpd   %zmm31, %zmm21, %zmm4                         #92.3
        vmovaps   %zmm6, %zmm5                                  #90.3
        vpermt2pd 560(%rdx), %zmm29, %zmm5                      #90.3
        vpermt2pd 560(%rdx), %zmm28, %zmm6                      #90.3
        vaddpd    %zmm31, %zmm4, %zmm4{%k4}                     #92.3
        vfmadd213pd %zmm12, %zmm8, %zmm5                        #90.3
        vsubpd    %zmm4, %zmm31, %zmm4{%k3}                     #92.3
        vpermpd   %zmm1, %zmm21, %zmm8                          #93.3
        vfmadd213pd %zmm5, %zmm7, %zmm6                         #90.3
        vmovups   96(%rcx,%r9), %zmm31                          #99.3
        vpermpd   %zmm4, %zmm2, %zmm7                           #92.3
        vpermpd   %zmm6, %zmm21, %zmm0                          #94.3
        vaddpd    %zmm1, %zmm8, %zmm8{%k4}                      #93.3
        vaddpd    %zmm7, %zmm10, %zmm10{%k5}                    #92.3
        vaddpd    %zmm6, %zmm0, %zmm0{%k4}                      #94.3
        vsubpd    %zmm8, %zmm1, %zmm8{%k3}                      #93.3
        vsubpd    %zmm7, %zmm10, %zmm10{%k6}                    #92.3
        vsubpd    %zmm0, %zmm6, %zmm0{%k3}                      #94.3
        vpermpd   %zmm8, %zmm2, %zmm12                          #93.3
        vmovups   (%rcx,%r9), %zmm1                             #98.3
        vmovups   (%rcx,%r10), %zmm7                            #98.3
        vmovups   96(%rcx,%r10), %zmm8                          #99.3
        vpermpd   %zmm0, %zmm2, %zmm6                           #94.3
        vaddpd    %zmm12, %zmm9, %zmm9{%k5}                     #93.3
        vpermi2pd %zmm8, %zmm31, %zmm17                         #99.3
        vaddpd    %zmm6, %zmm11, %zmm11{%k5}                    #94.3
        vpermt2pd 160(%rcx,%r10), %zmm20, %zmm8                 #99.3
        vpermt2pd 160(%rcx,%r9), %zmm14, %zmm31                 #99.3
        vsubpd    %zmm6, %zmm11, %zmm11{%k6}                    #94.3
        vsubpd    %zmm12, %zmm9, %zmm9{%k6}                     #93.3
        kmovw     %eax, %k5                                     #101.8
        vmovaps   %zmm13, %zmm3                                 #98.3
        vpermi2pd %zmm7, %zmm1, %zmm3                           #98.3
        vpermt2pd 64(%rcx,%r9), %zmm15, %zmm1                   #98.3
        vpermt2pd 64(%rcx,%r10), %zmm16, %zmm7                  #98.3
        vaddpd    %zmm17, %zmm3, %zmm2{%k5}{z}                  #101.8
        movslq    12(%rdi), %rdi                                #108.15
        vmovaps   %zmm19, %zmm6                                 #98.3
        vmovaps   %zmm19, %zmm20                                #99.3
        vpermi2pd %zmm1, %zmm7, %zmm6                           #98.3
        lea       (%rdi,%rdi,2), %r9                            #108.8
        vpermt2pd %zmm1, %zmm18, %zmm7                          #98.3
        vpermi2pd %zmm31, %zmm8, %zmm20                         #99.3
        vpermt2pd %zmm31, %zmm18, %zmm8                         #99.3
        vsubpd    %zmm17, %zmm3, %zmm2{%k7}                     #102.8
        vaddpd    %zmm20, %zmm6, %zmm3{%k5}{z}                  #103.8
        vaddpd    %zmm8, %zmm7, %zmm4{%k5}{z}                   #105.8
        vsubpd    %zmm20, %zmm6, %zmm3{%k7}                     #104.8
        vsubpd    %zmm8, %zmm7, %zmm4{%k7}                      #106.8
        shlq      $6, %r9                                       #108.8
        prefetcht0 (%rcx,%r9)                                   #109.3
        movslq    12(%rsi), %rsi                                #110.15
        vpermilpd $85, %zmm2, %zmm6                             #115.3
        vpermilpd $85, %zmm3, %zmm17                            #115.3
        vpermilpd $85, %zmm4, %zmm14                            #115.3
        lea       (%rsi,%rsi,2), %rax                           #110.8
        shlq      $6, %rax                                      #110.8
        movl      $150, %esi                                    #127.8
        kmovw     %esi, %k5                                     #127.8
        movl      $105, %esi                                    #128.8
        kmovw     %esi, %k6                                     #128.8
        vmulpd    %zmm30, %zmm6, %zmm5                          #115.3
        vmulpd    %zmm17, %zmm30, %zmm12                        #115.3
        vmulpd    %zmm14, %zmm30, %zmm8                         #115.3
        prefetcht0 (%rcx,%rax)                                  #111.3
        vmovups   576(%rdx), %zmm20                             #115.3
        vmovups   640(%rdx), %zmm7                              #115.3
        vmovups   704(%rdx), %zmm6                              #115.3
        vmovaps   %zmm20, %zmm17                                #115.3
        vpermt2pd 720(%rdx), %zmm29, %zmm17                     #115.3
        vmulpd    %zmm2, %zmm17, %zmm14                         #115.3
        vmovaps   %zmm20, %zmm0                                 #115.3
        vpermt2pd 720(%rdx), %zmm28, %zmm0                      #115.3
        vfmadd213pd %zmm14, %zmm5, %zmm0                        #115.3
        vmovaps   %zmm20, %zmm1                                 #115.3
        vpermt2pd 720(%rdx), %zmm27, %zmm1                      #115.3
        vfmadd213pd %zmm0, %zmm3, %zmm1                         #115.3
        vmovaps   %zmm20, %zmm31                                #115.3
        vmovaps   %zmm20, %zmm0                                 #115.3
        vpermt2pd 720(%rdx), %zmm26, %zmm31                     #115.3
        vpermt2pd 720(%rdx), %zmm23, %zmm0                      #115.3
        vfmadd213pd %zmm1, %zmm12, %zmm31                       #115.3
        vmulpd    %zmm0, %zmm2, %zmm1                           #115.3
        vmovaps   %zmm20, %zmm17                                #115.3
        vpermt2pd 784(%rdx), %zmm25, %zmm17                     #115.3
        vfmadd213pd %zmm31, %zmm4, %zmm17                       #115.3
        vmovaps   %zmm20, %zmm14                                #115.3
        vpermt2pd 720(%rdx), %zmm22, %zmm20                     #115.3
        vpermt2pd 784(%rdx), %zmm24, %zmm14                     #115.3
        vfmadd213pd %zmm1, %zmm5, %zmm20                        #115.3
        vfmadd213pd %zmm17, %zmm8, %zmm14                       #115.3
        vmovaps   %zmm7, %zmm17                                 #115.3
        vpermt2pd 784(%rdx), %zmm29, %zmm17                     #115.3
        vfmadd213pd %zmm20, %zmm3, %zmm17                       #115.3
        vmovaps   %zmm7, %zmm20                                 #115.3
        vpermt2pd 784(%rdx), %zmm28, %zmm20                     #115.3
        vmovaps   %zmm7, %zmm1                                  #115.3
        vpermt2pd 720(%rdx), %zmm25, %zmm1                      #115.3
        vfmadd213pd %zmm17, %zmm12, %zmm20                      #115.3
        vmulpd    %zmm1, %zmm2, %zmm2                           #115.3
        vmovaps   %zmm7, %zmm0                                  #115.3
        vpermt2pd 784(%rdx), %zmm27, %zmm0                      #115.3
        vfmadd213pd %zmm20, %zmm4, %zmm0                        #115.3
        vmovaps   %zmm7, %zmm17                                 #115.3
        vpermt2pd 784(%rdx), %zmm26, %zmm17                     #115.3
        vfmadd213pd %zmm0, %zmm8, %zmm17                        #115.3
        vmovaps   %zmm7, %zmm0                                  #115.3
        vpermt2pd 720(%rdx), %zmm24, %zmm0                      #115.3
        vfmadd213pd %zmm2, %zmm5, %zmm0                         #115.3
        vmovaps   %zmm7, %zmm5                                  #115.3
        vpermt2pd 784(%rdx), %zmm23, %zmm5                      #115.3
        vpermt2pd 784(%rdx), %zmm22, %zmm7                      #115.3
        vfmadd213pd %zmm0, %zmm3, %zmm5                         #115.3
        vmovups   .L_2il0floatpacket.33(%rip), %zmm0            #117.3
        vfmadd213pd %zmm5, %zmm12, %zmm7                        #115.3
        vmovups   96(%rcx,%rax), %zmm5                          #124.3
        vmovaps   %zmm6, %zmm3                                  #115.3
        vpermt2pd 848(%rdx), %zmm29, %zmm3                      #115.3
        vpermt2pd 848(%rdx), %zmm28, %zmm6                      #115.3
        vfmadd213pd %zmm7, %zmm4, %zmm3                         #115.3
        vpermpd   %zmm14, %zmm21, %zmm4                         #117.3
        vfmadd213pd %zmm3, %zmm8, %zmm6                         #115.3
        vmovups   (%rcx,%rax), %zmm3                            #123.3
        vaddpd    %zmm14, %zmm4, %zmm4{%k4}                     #117.3
        vpermpd   %zmm6, %zmm21, %zmm1                          #119.3
        vsubpd    %zmm4, %zmm14, %zmm4{%k3}                     #117.3
        vaddpd    %zmm6, %zmm1, %zmm1{%k4}                      #119.3
        vpermpd   %zmm4, %zmm0, %zmm12                          #117.3
        vpermpd   %zmm17, %zmm21, %zmm4                         #118.3
        vsubpd    %zmm1, %zmm6, %zmm1{%k3}                      #119.3
        vaddpd    %zmm12, %zmm10, %zmm10{%k1}                   #117.3
        vaddpd    %zmm17, %zmm4, %zmm4{%k4}                     #118.3
        vpermpd   %zmm1, %zmm0, %zmm2                           #119.3
        vsubpd    %zmm12, %zmm10, %zmm10{%k2}                   #117.3
        vsubpd    %zmm4, %zmm17, %zmm4{%k3}                     #118.3
        vaddpd    %zmm2, %zmm11, %zmm11{%k1}                    #119.3
        vpermpd   %zmm4, %zmm0, %zmm12                          #118.3
        vsubpd    %zmm2, %zmm11, %zmm11{%k2}                    #119.3
        vmovups   (%rcx,%r9), %zmm4                             #123.3
        vmovups   96(%rcx,%r9), %zmm2                           #124.3
        vaddpd    %zmm12, %zmm9, %zmm9{%k1}                     #118.3
        vmovaps   %zmm13, %zmm6                                 #123.3
        vpermi2pd %zmm4, %zmm3, %zmm6                           #123.3
        vpermt2pd 64(%rcx,%r9), %zmm16, %zmm4                   #123.3
        vpermt2pd 64(%rcx,%rax), %zmm15, %zmm3                  #123.3
        vpermi2pd %zmm2, %zmm5, %zmm13                          #124.3
        vpermt2pd 160(%rcx,%rax), %zmm15, %zmm5                 #124.3
        vpermt2pd 160(%rcx,%r9), %zmm16, %zmm2                  #124.3
        vsubpd    %zmm12, %zmm9, %zmm9{%k2}                     #118.3
        vmovaps   %zmm19, %zmm0                                 #123.3
        vpermi2pd %zmm3, %zmm4, %zmm0                           #123.3
        vpermt2pd %zmm3, %zmm18, %zmm4                          #123.3
        vpermi2pd %zmm5, %zmm2, %zmm19                          #124.3
        vpermt2pd %zmm5, %zmm18, %zmm2                          #124.3
        vpermilpd $85, %zmm13, %zmm18                           #126.8
        vaddpd    %zmm18, %zmm6, %zmm13{%k5}{z}                 #127.8
        vpermilpd $85, %zmm19, %zmm1                            #129.8
        vsubpd    %zmm18, %zmm6, %zmm13{%k6}                    #128.8
        vaddpd    %zmm1, %zmm0, %zmm12{%k5}{z}                  #130.8
                                # LOE rdx rbx r8 r12 r13 r14 r15 zmm0 zmm1 zmm2 zmm4 zmm9 zmm10 zmm11 zmm12 zmm13 zmm21 zmm22 zmm23 zmm24 zmm25 zmm26 zmm27 zmm28 zmm29 zmm30 k3 k4 k5 k6
..B1.4:                         # Preds ..B1.1
                                # Execution count [1.00e+00]
        vpermilpd $85, %zmm2, %zmm14                            #132.8
        movl      $111, %eax                                    #140.3
        vmovups   864(%rdx), %zmm3                              #138.3
        vmovups   1008(%rdx), %zmm2                             #138.3
        vmovups   928(%rdx), %zmm5                              #138.3
        vmovups   992(%rdx), %zmm7                              #138.3
        vmovups   1136(%rdx), %zmm6                             #138.3
        vsubpd    %zmm1, %zmm0, %zmm12{%k6}                     #131.8
        vaddpd    %zmm14, %zmm4, %zmm8{%k5}{z}                  #133.8
        kmovw     %eax, %k1                                     #140.3
        vsubpd    %zmm14, %zmm4, %zmm8{%k6}                     #134.8
        vmovups   1072(%rdx), %zmm4                             #138.3
        vmovaps   %zmm29, %zmm18                                #138.3
        movl      $144, %eax                                    #140.3
        vpermi2pd %zmm2, %zmm3, %zmm18                          #138.3
        kmovw     %eax, %k2                                     #140.3
        vmulpd    %zmm13, %zmm18, %zmm19                        #138.3
        vpermilpd $85, %zmm13, %zmm15                           #138.3
        vpermilpd $85, %zmm12, %zmm16                           #138.3
        vmulpd    %zmm30, %zmm15, %zmm1                         #138.3
        vmulpd    %zmm16, %zmm30, %zmm0                         #138.3
        vmovaps   %zmm23, %zmm16                                #138.3
        vmovaps   %zmm25, %zmm15                                #138.3
        vpermi2pd %zmm2, %zmm3, %zmm16                          #138.3
        vpermi2pd %zmm2, %zmm5, %zmm25                          #138.3
        vpermi2pd %zmm4, %zmm3, %zmm15                          #138.3
        vpermi2pd %zmm4, %zmm5, %zmm23                          #138.3
        vpermilpd $85, %zmm8, %zmm17                            #138.3
        vmulpd    %zmm17, %zmm30, %zmm30                        #138.3
        vmulpd    %zmm16, %zmm13, %zmm17                        #138.3
        vmulpd    %zmm25, %zmm13, %zmm13                        #138.3
        vmovaps   %zmm28, %zmm20                                #138.3
        vpermi2pd %zmm2, %zmm3, %zmm20                          #138.3
        vfmadd213pd %zmm19, %zmm1, %zmm20                       #138.3
        vmovaps   %zmm27, %zmm31                                #138.3
        vmovaps   %zmm26, %zmm14                                #138.3
        vmovaps   %zmm24, %zmm19                                #138.3
        vpermi2pd %zmm2, %zmm3, %zmm31                          #138.3
        vpermi2pd %zmm2, %zmm3, %zmm14                          #138.3
        vpermi2pd %zmm4, %zmm3, %zmm19                          #138.3
        vpermt2pd %zmm2, %zmm22, %zmm3                          #138.3
        vpermi2pd %zmm2, %zmm5, %zmm24                          #138.3
        vpermi2pd %zmm4, %zmm5, %zmm27                          #138.3
        vpermi2pd %zmm4, %zmm5, %zmm26                          #138.3
        vfmadd213pd %zmm17, %zmm1, %zmm3                        #138.3
        vfmadd213pd %zmm13, %zmm1, %zmm24                       #138.3
        vfmadd213pd %zmm20, %zmm12, %zmm31                      #138.3
        vfmadd213pd %zmm24, %zmm12, %zmm23                      #138.3
        vfmadd213pd %zmm31, %zmm0, %zmm14                       #138.3
        vmovups   .L_2il0floatpacket.34(%rip), %zmm24           #140.3
        vfmadd213pd %zmm14, %zmm8, %zmm15                       #138.3
        vmovaps   %zmm29, %zmm18                                #138.3
        vpermi2pd %zmm4, %zmm5, %zmm18                          #138.3
        vpermi2pd %zmm6, %zmm7, %zmm29                          #138.3
        vpermt2pd %zmm6, %zmm28, %zmm7                          #138.3
        vfmadd213pd %zmm3, %zmm12, %zmm18                       #138.3
        vfmadd213pd %zmm15, %zmm30, %zmm19                      #138.3
        vmovaps   %zmm28, %zmm3                                 #138.3
        vpermi2pd %zmm4, %zmm5, %zmm3                           #138.3
        vpermt2pd %zmm4, %zmm22, %zmm5                          #138.3
        vpermpd   %zmm19, %zmm21, %zmm12                        #140.3
        vfmadd213pd %zmm18, %zmm0, %zmm3                        #138.3
        vfmadd213pd %zmm23, %zmm0, %zmm5                        #138.3
        vmovups   .L_2il0floatpacket.35(%rip), %zmm28           #150.3
        vaddpd    %zmm19, %zmm12, %zmm12{%k4}                   #140.3
        vfmadd213pd %zmm3, %zmm8, %zmm27                        #138.3
        vfmadd213pd %zmm5, %zmm8, %zmm29                        #138.3
        vsubpd    %zmm12, %zmm19, %zmm12{%k3}                   #140.3
        vfmadd213pd %zmm27, %zmm30, %zmm26                      #138.3
        vfmadd213pd %zmm29, %zmm30, %zmm7                       #138.3
        vmovups   .L_2il0floatpacket.36(%rip), %zmm29           #150.3
        vpermpd   %zmm26, %zmm21, %zmm23                        #141.3
        vpermpd   %zmm7, %zmm21, %zmm21                         #142.3
        vpermpd   %zmm12, %zmm24, %zmm22                        #140.3
        vaddpd    %zmm26, %zmm23, %zmm23{%k4}                   #141.3
        vaddpd    %zmm7, %zmm21, %zmm21{%k4}                    #142.3
        vaddpd    %zmm22, %zmm10, %zmm10{%k1}                   #140.3
        vsubpd    %zmm23, %zmm26, %zmm23{%k3}                   #141.3
        vsubpd    %zmm21, %zmm7, %zmm21{%k3}                    #142.3
        vsubpd    %zmm22, %zmm10, %zmm10{%k2}                   #140.3
        vpermpd   %zmm23, %zmm24, %zmm26                        #141.3
        vpermpd   %zmm21, %zmm24, %zmm25                        #142.3
        vbroadcastsd -16(%rbp), %zmm27                          #145.10
        vaddpd    %zmm26, %zmm9, %zmm9{%k1}                     #141.3
        vaddpd    %zmm25, %zmm11, %zmm11{%k1}                   #142.3
        vmulpd    %zmm10, %zmm27, %zmm10                        #146.8
        vsubpd    %zmm26, %zmm9, %zmm9{%k2}                     #141.3
        vsubpd    %zmm25, %zmm11, %zmm11{%k2}                   #142.3
        vmulpd    %zmm9, %zmm27, %zmm0                          #147.8
        vmulpd    %zmm11, %zmm27, %zmm9                         #148.8
        vmovups   .L_2il0floatpacket.37(%rip), %zmm11           #150.3
        vpermi2pd %zmm0, %zmm10, %zmm28                         #150.3
        vpermi2pd %zmm10, %zmm9, %zmm29                         #150.3
        vpermt2pd %zmm9, %zmm11, %zmm0                          #150.3
        vmovupd   %ymm28, (%r8)                                 #150.3
        vmovupd   %ymm29, 32(%r8)                               #150.3
        vmovupd   %ymm0, 64(%r8)                                #150.3
        vextractf64x4 $1, %zmm28, 96(%r8)                       #150.3
        vextractf64x4 $1, %zmm29, 128(%r8)                      #150.3
        vextractf64x4 $1, %zmm0, 160(%r8)                       #150.3
        vzeroupper                                              #151.1
        movq      %rbp, %rsp                                    #151.1
        popq      %rbp                                          #151.1
	.cfi_restore 6
        ret                                                     #151.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	doe_dble_avx512,@function
	.size	doe_dble_avx512,.-doe_dble_avx512
	.data
# -- End  doe_dble_avx512
	.text
# -- Begin  deo_dble_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl deo_dble_avx512
# --- deo_dble_avx512(const int *, const int *, const su3_dble *, spinor_dble *, double, spin_t *)
deo_dble_avx512:
# parameter 1: %rdi
# parameter 2: %rsi
# parameter 3: %rdx
# parameter 4: %rcx
# parameter 5: %xmm0
# parameter 6: %r8
..B2.1:                         # Preds ..B2.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_deo_dble_avx512.8:
..L9:
                                                          #154.1
        pushq     %rbp                                          #154.1
	.cfi_def_cfa_offset 16
        movq      %rsp, %rbp                                    #154.1
	.cfi_def_cfa 6, 16
	.cfi_offset 6, -16
        andq      $-64, %rsp                                    #154.1
        movslq    (%rdi), %rax                                  #168.16
        lea       (%rax,%rax,2), %r11                           #168.8
        shlq      $6, %r11                                      #168.8
        prefetcht0 (%rcx,%r11)                                  #169.3
        movl      $15, %eax                                     #181.3
        movslq    (%rsi), %r9                                   #170.16
        kmovw     %eax, %k5                                     #181.3
        movl      $240, %eax                                    #181.3
        kmovw     %eax, %k6                                     #181.3
        vbroadcastsd %xmm0, %zmm24                              #176.10
        movl      $90, %eax                                     #198.3
        lea       (%r9,%r9,2), %r10                             #170.8
        shlq      $6, %r10                                      #170.8
        kmovw     %eax, %k1                                     #198.3
        movl      $165, %eax                                    #198.3
        kmovw     %eax, %k2                                     #198.3
        movl      $195, %eax                                    #215.3
        kmovw     %eax, %k4                                     #215.3
        movl      $60, %eax                                     #215.3
        kmovw     %eax, %k3                                     #215.3
        prefetcht0 (%rcx,%r10)                                  #171.3
        vmovups   96(%r8), %zmm27                               #173.3
        vmovups   (%r8), %zmm23                                 #173.3
        vmovups   .L_2il0floatpacket.14(%rip), %zmm26           #173.3
        vmovups   .L_2il0floatpacket.15(%rip), %zmm30           #173.3
        vmovups   .L_2il0floatpacket.16(%rip), %zmm29           #173.3
        vmovups   .L_2il0floatpacket.17(%rip), %zmm25           #173.3
        vmovups   .L_2il0floatpacket.18(%rip), %zmm28           #173.3
        vmovups   .L_2il0floatpacket.28(%rip), %zmm20           #181.3
        vmovups   144(%rdx), %zmm18                             #187.3
        vmovups   208(%rdx), %zmm13                             #187.3
        vmovups   272(%rdx), %zmm1                              #187.3
        vmovups   (%rcx,%r10), %zmm3                            #189.3
        vpermi2pd %zmm23, %zmm27, %zmm26                        #173.3
        vpermt2pd 160(%r8), %zmm30, %zmm27                      #173.3
        vpermt2pd 64(%r8), %zmm29, %zmm23                       #173.3
        vmulpd    %zmm26, %zmm24, %zmm30                        #177.8
        vpermi2pd %zmm27, %zmm23, %zmm25                        #173.3
        vpermt2pd %zmm27, %zmm28, %zmm23                        #173.3
        vpermpd   %zmm30, %zmm20, %zmm22                        #181.3
        vmulpd    %zmm25, %zmm24, %zmm29                        #178.8
        vmulpd    %zmm23, %zmm24, %zmm28                        #179.8
        vmovups   .L_2il0floatpacket.19(%rip), %zmm27           #187.3
        vpermpd   %zmm29, %zmm20, %zmm21                        #182.3
        vaddpd    %zmm30, %zmm22, %zmm6{%k5}{z}                 #181.3
        vpermpd   %zmm28, %zmm20, %zmm19                        #183.3
        vaddpd    %zmm29, %zmm21, %zmm10{%k5}{z}                #182.3
        vsubpd    %zmm30, %zmm22, %zmm6{%k6}                    #181.3
        vaddpd    %zmm28, %zmm19, %zmm12{%k5}{z}                #183.3
        vsubpd    %zmm29, %zmm21, %zmm10{%k6}                   #182.3
        vsubpd    %zmm28, %zmm19, %zmm12{%k6}                   #183.3
        vmovups   .L_2il0floatpacket.27(%rip), %zmm19           #187.3
        vpermilpd $85, %zmm10, %zmm26                           #187.3
        vmulpd    %zmm26, %zmm27, %zmm16                        #187.3
        vmovups   .L_2il0floatpacket.20(%rip), %zmm26           #187.3
        vmovaps   %zmm18, %zmm24                                #187.3
        vpermt2pd (%rdx), %zmm26, %zmm24                        #187.3
        vpermilpd $85, %zmm6, %zmm15                            #187.3
        vmulpd    %zmm27, %zmm15, %zmm17                        #187.3
        vmulpd    %zmm6, %zmm24, %zmm23                         #187.3
        vmovups   .L_2il0floatpacket.22(%rip), %zmm24           #187.3
        vpermilpd $85, %zmm12, %zmm25                           #187.3
        vmulpd    %zmm25, %zmm27, %zmm15                        #187.3
        vmovups   .L_2il0floatpacket.21(%rip), %zmm25           #187.3
        vmovaps   %zmm18, %zmm22                                #187.3
        vpermt2pd (%rdx), %zmm25, %zmm22                        #187.3
        vfmadd213pd %zmm23, %zmm17, %zmm22                      #187.3
        vmovups   .L_2il0floatpacket.23(%rip), %zmm23           #187.3
        vmovaps   %zmm18, %zmm14                                #187.3
        vpermt2pd (%rdx), %zmm24, %zmm14                        #187.3
        vfmadd213pd %zmm22, %zmm10, %zmm14                      #187.3
        vmovups   .L_2il0floatpacket.24(%rip), %zmm22           #187.3
        vmovaps   %zmm18, %zmm21                                #187.3
        vpermt2pd (%rdx), %zmm23, %zmm21                        #187.3
        vfmadd213pd %zmm14, %zmm16, %zmm21                      #187.3
        vmovaps   %zmm18, %zmm20                                #187.3
        vpermt2pd 64(%rdx), %zmm22, %zmm20                      #187.3
        vfmadd213pd %zmm21, %zmm12, %zmm20                      #187.3
        vmovups   .L_2il0floatpacket.25(%rip), %zmm21           #187.3
        vmovaps   %zmm18, %zmm14                                #187.3
        vpermt2pd 64(%rdx), %zmm21, %zmm14                      #187.3
        vfmadd213pd %zmm20, %zmm15, %zmm14                      #187.3
        vmovups   .L_2il0floatpacket.26(%rip), %zmm20           #187.3
        vmovaps   %zmm18, %zmm11                                #187.3
        vpermt2pd (%rdx), %zmm20, %zmm11                        #187.3
        vpermt2pd (%rdx), %zmm19, %zmm18                        #187.3
        vmulpd    %zmm11, %zmm6, %zmm31                         #187.3
        vfmadd213pd %zmm31, %zmm17, %zmm18                      #187.3
        vmovaps   %zmm13, %zmm8                                 #187.3
        vmovaps   %zmm13, %zmm2                                 #187.3
        vpermt2pd 64(%rdx), %zmm26, %zmm8                       #187.3
        vpermt2pd (%rdx), %zmm22, %zmm2                         #187.3
        vfmadd213pd %zmm18, %zmm10, %zmm8                       #187.3
        vmulpd    %zmm2, %zmm6, %zmm18                          #187.3
        vmovaps   %zmm13, %zmm5                                 #187.3
        vpermt2pd (%rdx), %zmm21, %zmm5                         #187.3
        vfmadd213pd %zmm18, %zmm17, %zmm5                       #187.3
        vmovups   .L_2il0floatpacket.38(%rip), %zmm18           #189.3
        vmovaps   %zmm13, %zmm17                                #187.3
        vpermt2pd 64(%rdx), %zmm20, %zmm17                      #187.3
        vfmadd213pd %zmm5, %zmm10, %zmm17                       #187.3
        vmovups   64(%rcx,%r10), %zmm10                         #189.3
        vmovaps   %zmm13, %zmm9                                 #187.3
        vmovaps   %zmm13, %zmm7                                 #187.3
        vmovaps   %zmm13, %zmm11                                #187.3
        vpermt2pd 64(%rdx), %zmm19, %zmm13                      #187.3
        vpermt2pd 64(%rdx), %zmm25, %zmm9                       #187.3
        vpermt2pd 64(%rdx), %zmm24, %zmm7                       #187.3
        vpermt2pd 64(%rdx), %zmm23, %zmm11                      #187.3
        vpermt2pd 64(%rcx,%r11), %zmm18, %zmm10                 #189.3
        vfmadd213pd %zmm17, %zmm16, %zmm13                      #187.3
        vfmadd213pd %zmm8, %zmm16, %zmm9                        #187.3
        vmovups   96(%rcx,%r10), %zmm8                          #190.3
        vmovups   .L_2il0floatpacket.39(%rip), %zmm17           #189.3
        vfmadd213pd %zmm9, %zmm12, %zmm7                        #187.3
        vmovups   160(%rcx,%r10), %zmm9                         #190.3
        vfmadd213pd %zmm7, %zmm15, %zmm11                       #187.3
        vpermt2pd 160(%rcx,%r11), %zmm18, %zmm9                 #190.3
        vmovaps   %zmm1, %zmm16                                 #187.3
        vpermt2pd 128(%rdx), %zmm26, %zmm16                     #187.3
        vpermt2pd 128(%rdx), %zmm25, %zmm1                      #187.3
        vfmadd213pd %zmm13, %zmm12, %zmm16                      #187.3
        vfmadd213pd %zmm16, %zmm15, %zmm1                       #187.3
        vmovups   .L_2il0floatpacket.36(%rip), %zmm15           #189.3
        vmovups   .L_2il0floatpacket.35(%rip), %zmm16           #189.3
        vmovaps   %zmm1, %zmm12                                 #189.3
        vmovaps   %zmm14, %zmm13                                #189.3
        vpermt2pd %zmm14, %zmm15, %zmm12                        #189.3
        vpermt2pd %zmm11, %zmm16, %zmm13                        #189.3
        vmovups   .L_2il0floatpacket.37(%rip), %zmm14           #189.3
        vmovaps   %zmm8, %zmm0                                  #190.3
        vpermt2pd %zmm1, %zmm14, %zmm11                         #189.3
        vpermt2pd 96(%rcx,%r11), %zmm18, %zmm0                  #190.3
        vpermt2pd 96(%rcx,%r11), %zmm17, %zmm8                  #190.3
        vaddpd    %zmm11, %zmm9, %zmm9{%k5}                     #190.3
        vaddpd    %zmm13, %zmm0, %zmm0{%k5}                     #190.3
        vaddpd    %zmm12, %zmm8, %zmm8{%k5}                     #190.3
        vaddpd    %zmm11, %zmm10, %zmm31                        #189.3
        vsubpd    %zmm13, %zmm0, %zmm0{%k6}                     #190.3
        vsubpd    %zmm12, %zmm8, %zmm8{%k6}                     #190.3
        vsubpd    %zmm11, %zmm9, %zmm9{%k6}                     #190.3
        vmovaps   %zmm3, %zmm4                                  #189.3
        vpermt2pd (%rcx,%r11), %zmm18, %zmm4                    #189.3
        vpermt2pd (%rcx,%r11), %zmm17, %zmm3                    #189.3
        vaddpd    %zmm13, %zmm4, %zmm6                          #189.3
        vaddpd    %zmm12, %zmm3, %zmm1                          #189.3
        movslq    4(%rdi), %rax                                 #193.16
        vmovupd   %ymm6, (%rcx,%r10)                            #189.3
        vmovupd   %ymm1, 32(%rcx,%r10)                          #189.3
        vmovupd   %ymm31, 64(%rcx,%r10)                         #189.3
        vextractf64x4 $1, %zmm6, (%rcx,%r11)                    #189.3
        vextractf64x4 $1, %zmm1, 32(%rcx,%r11)                  #189.3
        vextractf64x4 $1, %zmm31, 64(%rcx,%r11)                 #189.3
        vmovupd   %ymm0, 96(%rcx,%r10)                          #190.3
        vmovupd   %ymm8, 128(%rcx,%r10)                         #190.3
        vmovupd   %ymm9, 160(%rcx,%r10)                         #190.3
        vextractf64x4 $1, %zmm0, 96(%rcx,%r11)                  #190.3
        vextractf64x4 $1, %zmm8, 128(%rcx,%r11)                 #190.3
        vextractf64x4 $1, %zmm9, 160(%rcx,%r11)                 #190.3
        lea       (%rax,%rax,2), %r10                           #193.8
        shlq      $6, %r10                                      #193.8
        prefetcht0 (%rcx,%r10)                                  #194.3
        movslq    4(%rsi), %r8                                  #195.16
        lea       (%r8,%r8,2), %r9                              #195.8
        shlq      $6, %r9                                       #195.8
        prefetcht0 (%rcx,%r9)                                   #196.3
        vmovups   .L_2il0floatpacket.40(%rip), %zmm7            #198.3
        vmovups   .L_2il0floatpacket.41(%rip), %zmm4            #198.3
        vmovups   432(%rdx), %zmm6                              #204.3
        vmovups   496(%rdx), %zmm8                              #204.3
        vmovups   560(%rdx), %zmm9                              #204.3
        vpermpd   %zmm30, %zmm7, %zmm12                         #198.3
        vpermpd   %zmm30, %zmm4, %zmm13                         #198.3
        vpermpd   %zmm29, %zmm7, %zmm11                         #199.3
        vpermpd   %zmm28, %zmm7, %zmm3                          #200.3
        vaddpd    %zmm12, %zmm13, %zmm5{%k1}{z}                 #198.3
        vsubpd    %zmm12, %zmm13, %zmm5{%k2}                    #198.3
        vpermpd   %zmm29, %zmm4, %zmm12                         #199.3
        vaddpd    %zmm11, %zmm12, %zmm2{%k1}{z}                 #199.3
        vsubpd    %zmm11, %zmm12, %zmm2{%k2}                    #199.3
        vpermpd   %zmm28, %zmm4, %zmm11                         #200.3
        vaddpd    %zmm3, %zmm11, %zmm7{%k1}{z}                  #200.3
        vsubpd    %zmm3, %zmm11, %zmm7{%k2}                     #200.3
        vmovaps   %zmm6, %zmm3                                  #204.3
        vpermt2pd 288(%rdx), %zmm26, %zmm3                      #204.3
        vpermilpd $85, %zmm5, %zmm10                            #204.3
        vmulpd    %zmm27, %zmm10, %zmm0                         #204.3
        vmulpd    %zmm5, %zmm3, %zmm10                          #204.3
        vpermilpd $85, %zmm2, %zmm1                             #204.3
        vpermilpd $85, %zmm7, %zmm4                             #204.3
        vmulpd    %zmm1, %zmm27, %zmm31                         #204.3
        vmulpd    %zmm4, %zmm27, %zmm1                          #204.3
        vmovaps   %zmm6, %zmm4                                  #204.3
        vpermt2pd 288(%rdx), %zmm25, %zmm4                      #204.3
        vfmadd213pd %zmm10, %zmm0, %zmm4                        #204.3
        vmovaps   %zmm6, %zmm10                                 #204.3
        vpermt2pd 288(%rdx), %zmm24, %zmm10                     #204.3
        vfmadd213pd %zmm4, %zmm2, %zmm10                        #204.3
        vmovaps   %zmm6, %zmm3                                  #204.3
        vpermt2pd 288(%rdx), %zmm23, %zmm3                      #204.3
        vfmadd213pd %zmm10, %zmm31, %zmm3                       #204.3
        vmovaps   %zmm6, %zmm4                                  #204.3
        vpermt2pd 352(%rdx), %zmm22, %zmm4                      #204.3
        vfmadd213pd %zmm3, %zmm7, %zmm4                         #204.3
        vmovaps   %zmm6, %zmm10                                 #204.3
        vmovaps   %zmm6, %zmm3                                  #204.3
        vpermt2pd 352(%rdx), %zmm21, %zmm10                     #204.3
        vpermt2pd 288(%rdx), %zmm20, %zmm3                      #204.3
        vpermt2pd 288(%rdx), %zmm19, %zmm6                      #204.3
        vfmadd213pd %zmm4, %zmm1, %zmm10                        #204.3
        vmulpd    %zmm3, %zmm5, %zmm4                           #204.3
        vfmadd213pd %zmm4, %zmm0, %zmm6                         #204.3
        vmovaps   %zmm8, %zmm3                                  #204.3
        vpermt2pd 352(%rdx), %zmm26, %zmm3                      #204.3
        vfmadd213pd %zmm6, %zmm2, %zmm3                         #204.3
        vmovaps   %zmm8, %zmm6                                  #204.3
        vpermt2pd 352(%rdx), %zmm25, %zmm6                      #204.3
        vfmadd213pd %zmm3, %zmm31, %zmm6                        #204.3
        vmovaps   %zmm8, %zmm3                                  #204.3
        vpermt2pd 288(%rdx), %zmm22, %zmm3                      #204.3
        vmovaps   %zmm8, %zmm4                                  #204.3
        vpermt2pd 352(%rdx), %zmm24, %zmm4                      #204.3
        vmulpd    %zmm3, %zmm5, %zmm5                           #204.3
        vfmadd213pd %zmm6, %zmm7, %zmm4                         #204.3
        vmovaps   %zmm8, %zmm6                                  #204.3
        vpermt2pd 352(%rdx), %zmm23, %zmm6                      #204.3
        vfmadd213pd %zmm4, %zmm1, %zmm6                         #204.3
        vmovaps   %zmm8, %zmm4                                  #204.3
        vpermt2pd 288(%rdx), %zmm21, %zmm4                      #204.3
        vfmadd213pd %zmm5, %zmm0, %zmm4                         #204.3
        vmovaps   %zmm8, %zmm0                                  #204.3
        vpermt2pd 352(%rdx), %zmm20, %zmm0                      #204.3
        vpermt2pd 352(%rdx), %zmm19, %zmm8                      #204.3
        vfmadd213pd %zmm4, %zmm2, %zmm0                         #204.3
        vfmadd213pd %zmm0, %zmm31, %zmm8                        #204.3
        vmovups   (%rcx,%r9), %zmm0                             #206.3
        vmovaps   %zmm9, %zmm2                                  #204.3
        vpermt2pd 416(%rdx), %zmm26, %zmm2                      #204.3
        vpermt2pd 416(%rdx), %zmm25, %zmm9                      #204.3
        vfmadd213pd %zmm8, %zmm7, %zmm2                         #204.3
        vmovups   64(%rcx,%r9), %zmm7                           #206.3
        vfmadd213pd %zmm2, %zmm1, %zmm9                         #204.3
        vpermt2pd 64(%rcx,%r10), %zmm18, %zmm7                  #206.3
        vmovaps   %zmm9, %zmm8                                  #206.3
        vmovaps   %zmm0, %zmm1                                  #206.3
        vpermt2pd (%rcx,%r10), %zmm17, %zmm0                    #206.3
        vpermt2pd %zmm10, %zmm15, %zmm8                         #206.3
        vpermt2pd (%rcx,%r10), %zmm18, %zmm1                    #206.3
        vaddpd    %zmm8, %zmm0, %zmm4                           #206.3
        vmovups   .L_2il0floatpacket.43(%rip), %zmm0            #207.3
        vmovaps   %zmm10, %zmm31                                #206.3
        vmovaps   %zmm6, %zmm5                                  #206.3
        vpermt2pd %zmm6, %zmm16, %zmm31                         #206.3
        vpermt2pd %zmm9, %zmm14, %zmm5                          #206.3
        vpermi2pd %zmm10, %zmm9, %zmm0                          #207.3
        vaddpd    %zmm31, %zmm1, %zmm2                          #206.3
        vaddpd    %zmm5, %zmm7, %zmm3                           #206.3
        vmovups   96(%rcx,%r9), %zmm7                           #207.3
        vmovups   .L_2il0floatpacket.42(%rip), %zmm1            #207.3
        vmovaps   %zmm10, %zmm31                                #207.3
        vmovups   .L_2il0floatpacket.44(%rip), %zmm10           #207.3
        vmovups   %zmm1, -64(%rsp)                              #207.3[spill]
        vpermt2pd %zmm6, %zmm1, %zmm31                          #207.3
        vpermt2pd %zmm9, %zmm10, %zmm6                          #207.3
        vmovaps   %zmm7, %zmm8                                  #207.3
        vpermt2pd 96(%rcx,%r10), %zmm18, %zmm8                  #207.3
        vpermt2pd 96(%rcx,%r10), %zmm17, %zmm7                  #207.3
        vaddpd    %zmm31, %zmm8, %zmm8{%k2}                     #207.3
        vaddpd    %zmm0, %zmm7, %zmm7{%k2}                      #207.3
        vsubpd    %zmm31, %zmm8, %zmm8{%k1}                     #207.3
        vsubpd    %zmm0, %zmm7, %zmm7{%k1}                      #207.3
        movslq    8(%rdi), %r11                                 #210.16
        lea       (%r11,%r11,2), %r8                            #210.8
        vmovupd   %ymm2, (%rcx,%r9)                             #206.3
        vmovupd   %ymm4, 32(%rcx,%r9)                           #206.3
        shlq      $6, %r8                                       #210.8
        vmovupd   %ymm3, 64(%rcx,%r9)                           #206.3
        vextractf64x4 $1, %zmm2, (%rcx,%r10)                    #206.3
        vmovups   160(%rcx,%r9), %zmm2                          #207.3
        vextractf64x4 $1, %zmm4, 32(%rcx,%r10)                  #206.3
        vextractf64x4 $1, %zmm3, 64(%rcx,%r10)                  #206.3
        vpermt2pd 160(%rcx,%r10), %zmm18, %zmm2                 #207.3
        vaddpd    %zmm6, %zmm2, %zmm2{%k2}                      #207.3
        vsubpd    %zmm6, %zmm2, %zmm2{%k1}                      #207.3
        vmovupd   %ymm8, 96(%rcx,%r9)                           #207.3
        vmovupd   %ymm7, 128(%rcx,%r9)                          #207.3
        vmovupd   %ymm2, 160(%rcx,%r9)                          #207.3
        vextractf64x4 $1, %zmm8, 96(%rcx,%r10)                  #207.3
        vextractf64x4 $1, %zmm7, 128(%rcx,%r10)                 #207.3
        vextractf64x4 $1, %zmm2, 160(%rcx,%r10)                 #207.3
        prefetcht0 (%rcx,%r8)                                   #211.3
        movslq    8(%rsi), %rax                                 #212.16
        lea       (%rax,%rax,2), %rax                           #212.8
        shlq      $6, %rax                                      #212.8
        prefetcht0 (%rcx,%rax)                                  #213.3
        vmovups   .L_2il0floatpacket.45(%rip), %zmm1            #215.3
        vpermpd   %zmm30, %zmm1, %zmm6                          #215.3
        vpermpd   %zmm29, %zmm1, %zmm9                          #216.3
        vaddpd    %zmm6, %zmm13, %zmm2{%k4}{z}                  #215.3
        vaddpd    %zmm9, %zmm12, %zmm5{%k4}{z}                  #216.3
        vsubpd    %zmm6, %zmm13, %zmm2{%k3}                     #215.3
        vpermpd   %zmm28, %zmm1, %zmm6                          #217.3
        vsubpd    %zmm9, %zmm12, %zmm5{%k3}                     #216.3
        vmovups   720(%rdx), %zmm1                              #221.3
        vmovups   784(%rdx), %zmm9                              #221.3
        vaddpd    %zmm6, %zmm11, %zmm8{%k4}{z}                  #217.3
        vpermilpd $85, %zmm2, %zmm31                            #221.3
        vmulpd    %zmm27, %zmm31, %zmm3                         #221.3
        vsubpd    %zmm6, %zmm11, %zmm8{%k3}                     #217.3
        vmovups   848(%rdx), %zmm6                              #221.3
        vmovaps   %zmm1, %zmm31                                 #221.3
        vpermt2pd 576(%rdx), %zmm26, %zmm31                     #221.3
        vpermilpd $85, %zmm5, %zmm0                             #221.3
        vmulpd    %zmm0, %zmm27, %zmm4                          #221.3
        vmulpd    %zmm2, %zmm31, %zmm0                          #221.3
        vmovaps   %zmm1, %zmm31                                 #221.3
        vpermt2pd 576(%rdx), %zmm25, %zmm31                     #221.3
        vfmadd213pd %zmm0, %zmm3, %zmm31                        #221.3
        vmovaps   %zmm1, %zmm0                                  #221.3
        vpermt2pd 576(%rdx), %zmm24, %zmm0                      #221.3
        vfmadd213pd %zmm31, %zmm5, %zmm0                        #221.3
        vmovaps   %zmm1, %zmm31                                 #221.3
        vpermt2pd 576(%rdx), %zmm23, %zmm31                     #221.3
        vpermilpd $85, %zmm8, %zmm7                             #221.3
        vmulpd    %zmm7, %zmm27, %zmm7                          #221.3
        vfmadd213pd %zmm0, %zmm4, %zmm31                        #221.3
        vmovaps   %zmm1, %zmm0                                  #221.3
        vpermt2pd 640(%rdx), %zmm22, %zmm0                      #221.3
        vfmadd213pd %zmm31, %zmm8, %zmm0                        #221.3
        vmovaps   %zmm1, %zmm31                                 #221.3
        vpermt2pd 640(%rdx), %zmm21, %zmm31                     #221.3
        vfmadd213pd %zmm0, %zmm7, %zmm31                        #221.3
        vmovaps   %zmm1, %zmm0                                  #221.3
        vpermt2pd 576(%rdx), %zmm20, %zmm0                      #221.3
        vpermt2pd 576(%rdx), %zmm19, %zmm1                      #221.3
        vmulpd    %zmm0, %zmm2, %zmm0                           #221.3
        vfmadd213pd %zmm0, %zmm3, %zmm1                         #221.3
        vmovaps   %zmm9, %zmm0                                  #221.3
        vpermt2pd 640(%rdx), %zmm26, %zmm0                      #221.3
        vfmadd213pd %zmm1, %zmm5, %zmm0                         #221.3
        vmovaps   %zmm9, %zmm1                                  #221.3
        vpermt2pd 640(%rdx), %zmm25, %zmm1                      #221.3
        vfmadd213pd %zmm0, %zmm4, %zmm1                         #221.3
        vmovaps   %zmm9, %zmm0                                  #221.3
        vpermt2pd 640(%rdx), %zmm24, %zmm0                      #221.3
        vfmadd213pd %zmm1, %zmm8, %zmm0                         #221.3
        vmovaps   %zmm9, %zmm1                                  #221.3
        vpermt2pd 640(%rdx), %zmm23, %zmm1                      #221.3
        vfmadd213pd %zmm0, %zmm7, %zmm1                         #221.3
        vmovaps   %zmm9, %zmm0                                  #221.3
        vpermt2pd 576(%rdx), %zmm22, %zmm0                      #221.3
        vmulpd    %zmm0, %zmm2, %zmm0                           #221.3
        vmovaps   %zmm9, %zmm2                                  #221.3
        vpermt2pd 576(%rdx), %zmm21, %zmm2                      #221.3
        vfmadd213pd %zmm0, %zmm3, %zmm2                         #221.3
        vmovaps   %zmm9, %zmm3                                  #221.3
        vpermt2pd 640(%rdx), %zmm20, %zmm3                      #221.3
        vpermt2pd 640(%rdx), %zmm19, %zmm9                      #221.3
        vfmadd213pd %zmm2, %zmm5, %zmm3                         #221.3
        vfmadd213pd %zmm3, %zmm4, %zmm9                         #221.3
        vmovaps   %zmm6, %zmm5                                  #221.3
        vpermt2pd 704(%rdx), %zmm26, %zmm5                      #221.3
        vpermt2pd 704(%rdx), %zmm25, %zmm6                      #221.3
        vfmadd213pd %zmm9, %zmm8, %zmm5                         #221.3
        vmovups   (%rcx,%rax), %zmm8                            #223.3
        vfmadd213pd %zmm5, %zmm7, %zmm6                         #221.3
        vmovups   64(%rcx,%rax), %zmm5                          #223.3
        vmovaps   %zmm8, %zmm2                                  #223.3
        vmovaps   %zmm31, %zmm4                                 #223.3
        vmovaps   %zmm6, %zmm3                                  #223.3
        vpermt2pd (%rcx,%r8), %zmm18, %zmm2                     #223.3
        vpermt2pd (%rcx,%r8), %zmm17, %zmm8                     #223.3
        vpermt2pd %zmm1, %zmm16, %zmm4                          #223.3
        vpermt2pd %zmm31, %zmm15, %zmm3                         #223.3
        vpermt2pd 64(%rcx,%r8), %zmm18, %zmm5                   #223.3
        vaddpd    %zmm4, %zmm2, %zmm0                           #223.3
        vaddpd    %zmm3, %zmm8, %zmm9                           #223.3
        vmovaps   %zmm1, %zmm7                                  #223.3
        vpermt2pd %zmm6, %zmm14, %zmm7                          #223.3
        vaddpd    %zmm7, %zmm5, %zmm5                           #223.3
        vmovupd   %ymm0, (%rcx,%rax)                            #223.3
        vmovupd   %ymm9, 32(%rcx,%rax)                          #223.3
                                # LOE rax rdx rcx rbx rsi rdi r8 r12 r13 r14 r15 zmm0 zmm1 zmm5 zmm6 zmm9 zmm10 zmm11 zmm12 zmm13 zmm14 zmm15 zmm16 zmm17 zmm18 zmm19 zmm20 zmm21 zmm22 zmm23 zmm24 zmm25 zmm26 zmm27 zmm28 zmm29 zmm30 zmm31 k1 k2 k3 k4 k5 k6
..B2.4:                         # Preds ..B2.1
                                # Execution count [1.00e+00]
        vmovups   .L_2il0floatpacket.46(%rip), %zmm7            #224.3
        vmovups   96(%rcx,%r8), %zmm8                           #224.3
        vpermpd   %zmm31, %zmm7, %zmm4                          #224.3
        vpermpd   %zmm1, %zmm7, %zmm2                           #224.3
        vmovups   160(%rcx,%rax), %zmm31                        #224.3
        vmovaps   %zmm18, %zmm1                                 #224.3
        vmovaps   %zmm15, %zmm3                                 #224.3
        vpermt2pd 160(%rcx,%r8), %zmm18, %zmm31                 #224.3
        vmovupd   %ymm5, 64(%rcx,%rax)                          #223.3
        vextractf64x4 $1, %zmm9, 32(%rcx,%r8)                   #223.3
        vextractf64x4 $1, %zmm0, (%rcx,%r8)                     #223.3
        vextractf64x4 $1, %zmm5, 64(%rcx,%r8)                   #223.3
        vpermpd   %zmm6, %zmm7, %zmm9                           #224.3
        vmovups   96(%rcx,%rax), %zmm6                          #224.3
        vpermi2pd %zmm4, %zmm9, %zmm3                           #224.3
        vpermi2pd %zmm8, %zmm6, %zmm1                           #224.3
        vpermt2pd %zmm8, %zmm17, %zmm6                          #224.3
        vmovaps   %zmm16, %zmm0                                 #224.3
        vpermi2pd %zmm2, %zmm4, %zmm0                           #224.3
        vpermt2pd %zmm9, %zmm14, %zmm2                          #224.3
        vaddpd    %zmm3, %zmm6, %zmm6{%k3}                      #224.3
        vaddpd    %zmm0, %zmm1, %zmm1{%k6}                      #224.3
        vaddpd    %zmm2, %zmm31, %zmm31{%k5}                    #224.3
        vsubpd    %zmm3, %zmm6, %zmm6{%k4}                      #224.3
        vsubpd    %zmm0, %zmm1, %zmm1{%k5}                      #224.3
        vsubpd    %zmm2, %zmm31, %zmm31{%k6}                    #224.3
        vmovupd   %ymm1, 96(%rcx,%rax)                          #224.3
        vmovupd   %ymm6, 128(%rcx,%rax)                         #224.3
        vmovupd   %ymm31, 160(%rcx,%rax)                        #224.3
        vextractf64x4 $1, %zmm1, 96(%rcx,%r8)                   #224.3
        vextractf64x4 $1, %zmm6, 128(%rcx,%r8)                  #224.3
        vextractf64x4 $1, %zmm31, 160(%rcx,%r8)                 #224.3
        movslq    12(%rdi), %rax                                #228.16
        lea       (%rax,%rax,2), %r8                            #228.8
        shlq      $6, %r8                                       #228.8
        prefetcht0 (%rcx,%r8)                                   #229.3
        movl      $150, %eax                                    #233.3
        movslq    12(%rsi), %rsi                                #230.16
        kmovw     %eax, %k4                                     #233.3
        movl      $105, %eax                                    #233.3
        kmovw     %eax, %k3                                     #233.3
        lea       (%rsi,%rsi,2), %rdi                           #230.8
        shlq      $6, %rdi                                      #230.8
        prefetcht0 (%rcx,%rdi)                                  #231.3
        vmovups   .L_2il0floatpacket.47(%rip), %zmm5            #233.3
        vmovups   1072(%rdx), %zmm1                             #239.3
        vmovups   928(%rdx), %zmm6                              #239.3
        vmovups   992(%rdx), %zmm31                             #239.3
        vmovups   1136(%rdx), %zmm0                             #239.3
        vpermpd   %zmm30, %zmm5, %zmm4                          #233.3
        vpermpd   %zmm28, %zmm5, %zmm28                         #235.3
        vaddpd    %zmm4, %zmm13, %zmm30{%k4}{z}                 #233.3
        vaddpd    %zmm28, %zmm11, %zmm2{%k4}{z}                 #235.3
        vsubpd    %zmm4, %zmm13, %zmm30{%k3}                    #233.3
        vpermpd   %zmm29, %zmm5, %zmm13                         #234.3
        vsubpd    %zmm28, %zmm11, %zmm2{%k3}                    #235.3
        vmovups   1008(%rdx), %zmm29                            #239.3
        vaddpd    %zmm13, %zmm12, %zmm3{%k4}{z}                 #234.3
        vsubpd    %zmm13, %zmm12, %zmm3{%k3}                    #234.3
        vmovups   864(%rdx), %zmm13                             #239.3
        vpermilpd $85, %zmm30, %zmm12                           #239.3
        vpermilpd $85, %zmm3, %zmm11                            #239.3
        vpermilpd $85, %zmm2, %zmm7                             #239.3
        vmulpd    %zmm27, %zmm12, %zmm28                        #239.3
        vmulpd    %zmm11, %zmm27, %zmm12                        #239.3
        vmulpd    %zmm7, %zmm27, %zmm11                         #239.3
        vmovaps   %zmm26, %zmm27                                #239.3
        vpermi2pd %zmm13, %zmm29, %zmm27                        #239.3
        vmulpd    %zmm30, %zmm27, %zmm8                         #239.3
        vmovaps   %zmm25, %zmm27                                #239.3
        vpermi2pd %zmm13, %zmm29, %zmm27                        #239.3
        vfmadd213pd %zmm8, %zmm28, %zmm27                       #239.3
        vmovaps   %zmm24, %zmm4                                 #239.3
        vpermi2pd %zmm13, %zmm29, %zmm4                         #239.3
        vpermi2pd %zmm6, %zmm1, %zmm24                          #239.3
        vfmadd213pd %zmm27, %zmm3, %zmm4                        #239.3
        vmovaps   %zmm23, %zmm5                                 #239.3
        vmovaps   %zmm20, %zmm9                                 #239.3
        vpermi2pd %zmm13, %zmm29, %zmm5                         #239.3
        vpermi2pd %zmm13, %zmm29, %zmm9                         #239.3
        vpermi2pd %zmm6, %zmm1, %zmm20                          #239.3
        vpermi2pd %zmm6, %zmm1, %zmm23                          #239.3
        vfmadd213pd %zmm4, %zmm12, %zmm5                        #239.3
        vmulpd    %zmm9, %zmm30, %zmm4                          #239.3
        vmovaps   %zmm22, %zmm7                                 #239.3
        vpermi2pd %zmm13, %zmm1, %zmm22                         #239.3
        vpermi2pd %zmm6, %zmm29, %zmm7                          #239.3
        vmulpd    %zmm22, %zmm30, %zmm22                        #239.3
        vfmadd213pd %zmm5, %zmm2, %zmm7                         #239.3
        vmovaps   %zmm21, %zmm27                                #239.3
        vpermi2pd %zmm6, %zmm29, %zmm27                         #239.3
        vpermt2pd %zmm13, %zmm19, %zmm29                        #239.3
        vpermi2pd %zmm13, %zmm1, %zmm21                         #239.3
        vfmadd213pd %zmm7, %zmm11, %zmm27                       #239.3
        vfmadd213pd %zmm4, %zmm28, %zmm29                       #239.3
        vfmadd213pd %zmm22, %zmm28, %zmm21                      #239.3
        vmovaps   %zmm26, %zmm5                                 #239.3
        vpermi2pd %zmm6, %zmm1, %zmm5                           #239.3
        vpermi2pd %zmm31, %zmm0, %zmm26                         #239.3
        vpermt2pd %zmm31, %zmm25, %zmm0                         #239.3
        vfmadd213pd %zmm29, %zmm3, %zmm5                        #239.3
        vfmadd213pd %zmm21, %zmm3, %zmm20                       #239.3
        vmovups   (%rcx,%rdi), %zmm21                           #241.3
        vmovaps   %zmm25, %zmm29                                #239.3
        vpermi2pd %zmm6, %zmm1, %zmm29                          #239.3
        vpermt2pd %zmm6, %zmm19, %zmm1                          #239.3
        vmovups   (%rcx,%r8), %zmm19                            #241.3
        vfmadd213pd %zmm5, %zmm12, %zmm29                       #239.3
        vfmadd213pd %zmm20, %zmm12, %zmm1                       #239.3
        vmovups   96(%rcx,%r8), %zmm25                          #242.3
        vfmadd213pd %zmm29, %zmm2, %zmm24                       #239.3
        vfmadd213pd %zmm1, %zmm2, %zmm26                        #239.3
        vmovups   160(%rcx,%rdi), %zmm2                         #242.3
        vmovups   96(%rcx,%rdi), %zmm1                          #242.3
        vfmadd213pd %zmm24, %zmm11, %zmm23                      #239.3
        vmovups   64(%rcx,%rdi), %zmm24                         #241.3
        vfmadd213pd %zmm26, %zmm11, %zmm0                       #239.3
        vpermt2pd 160(%rcx,%r8), %zmm18, %zmm2                  #242.3
        vpermt2pd 64(%rcx,%r8), %zmm18, %zmm24                  #241.3
        vpermi2pd %zmm23, %zmm27, %zmm16                        #241.3
        vpermi2pd %zmm0, %zmm23, %zmm14                         #241.3
        vpermi2pd %zmm23, %zmm27, %zmm10                        #242.3
        vpermi2pd %zmm27, %zmm0, %zmm15                         #241.3
        vaddpd    %zmm14, %zmm24, %zmm14                        #241.3
        vmovaps   %zmm18, %zmm20                                #241.3
        vmovaps   %zmm18, %zmm26                                #242.3
        vpermi2pd %zmm19, %zmm21, %zmm20                        #241.3
        vpermt2pd %zmm19, %zmm17, %zmm21                        #241.3
        vpermi2pd %zmm25, %zmm1, %zmm26                         #242.3
        vpermt2pd %zmm25, %zmm17, %zmm1                         #242.3
        vaddpd    %zmm16, %zmm20, %zmm16                        #241.3
        vaddpd    %zmm10, %zmm26, %zmm26{%k2}                   #242.3
        vaddpd    %zmm15, %zmm21, %zmm15                        #241.3
        vmovups   .L_2il0floatpacket.48(%rip), %zmm18           #242.3
        vmovups   -64(%rsp), %zmm17                             #242.3[spill]
        vsubpd    %zmm10, %zmm26, %zmm26{%k1}                   #242.3
        vpermi2pd %zmm27, %zmm0, %zmm18                         #242.3
        vpermt2pd %zmm0, %zmm17, %zmm23                         #242.3
        vaddpd    %zmm18, %zmm1, %zmm1{%k3}                     #242.3
        vaddpd    %zmm23, %zmm2, %zmm2{%k1}                     #242.3
        vsubpd    %zmm18, %zmm1, %zmm1{%k4}                     #242.3
        vsubpd    %zmm23, %zmm2, %zmm2{%k2}                     #242.3
        vmovupd   %ymm16, (%rcx,%rdi)                           #241.3
        vmovupd   %ymm15, 32(%rcx,%rdi)                         #241.3
        vmovupd   %ymm14, 64(%rcx,%rdi)                         #241.3
        vextractf64x4 $1, %zmm16, (%rcx,%r8)                    #241.3
        vextractf64x4 $1, %zmm15, 32(%rcx,%r8)                  #241.3
        vextractf64x4 $1, %zmm14, 64(%rcx,%r8)                  #241.3
        vmovupd   %ymm26, 96(%rcx,%rdi)                         #242.3
        vmovupd   %ymm1, 128(%rcx,%rdi)                         #242.3
        vmovupd   %ymm2, 160(%rcx,%rdi)                         #242.3
        vextractf64x4 $1, %zmm26, 96(%rcx,%r8)                  #242.3
        vextractf64x4 $1, %zmm1, 128(%rcx,%r8)                  #242.3
        vextractf64x4 $1, %zmm2, 160(%rcx,%r8)                  #242.3
        vzeroupper                                              #243.1
        movq      %rbp, %rsp                                    #243.1
        popq      %rbp                                          #243.1
	.cfi_def_cfa 7, 8
	.cfi_restore 6
        ret                                                     #243.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	deo_dble_avx512,@function
	.size	deo_dble_avx512,.-deo_dble_avx512
	.data
# -- End  deo_dble_avx512
	.section .rodata, "a"
	.align 64
	.align 64
.L_2il0floatpacket.14:
	.long	0x00000008,0x00000000,0x00000009,0x00000000,0x0000000e,0x00000000,0x0000000f,0x00000000,0x00000000,0x00000000,0x00000001,0x00000000,0x00000006,0x00000000,0x00000007,0x00000000
	.type	.L_2il0floatpacket.14,@object
	.size	.L_2il0floatpacket.14,64
	.align 64
.L_2il0floatpacket.15:
	.long	0x00000004,0x00000000,0x00000005,0x00000000,0x0000000a,0x00000000,0x0000000b,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000,0x00000008,0x00000000,0x00000009,0x00000000
	.type	.L_2il0floatpacket.15,@object
	.size	.L_2il0floatpacket.15,64
	.align 64
.L_2il0floatpacket.16:
	.long	0x00000002,0x00000000,0x00000003,0x00000000,0x00000008,0x00000000,0x00000009,0x00000000,0x00000004,0x00000000,0x00000005,0x00000000,0x0000000a,0x00000000,0x0000000b,0x00000000
	.type	.L_2il0floatpacket.16,@object
	.size	.L_2il0floatpacket.16,64
	.align 64
.L_2il0floatpacket.17:
	.long	0x00000000,0x00000000,0x00000001,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000,0x0000000c,0x00000000,0x0000000d,0x00000000,0x0000000e,0x00000000,0x0000000f,0x00000000
	.type	.L_2il0floatpacket.17,@object
	.size	.L_2il0floatpacket.17,64
	.align 64
.L_2il0floatpacket.18:
	.long	0x00000004,0x00000000,0x00000005,0x00000000,0x00000006,0x00000000,0x00000007,0x00000000,0x00000008,0x00000000,0x00000009,0x00000000,0x0000000a,0x00000000,0x0000000b,0x00000000
	.type	.L_2il0floatpacket.18,@object
	.size	.L_2il0floatpacket.18,64
	.align 64
.L_2il0floatpacket.19:
	.long	0x00000000,0xbff00000,0x00000000,0x3ff00000,0x00000000,0xbff00000,0x00000000,0x3ff00000,0x00000000,0x3ff00000,0x00000000,0xbff00000,0x00000000,0x3ff00000,0x00000000,0xbff00000
	.type	.L_2il0floatpacket.19,@object
	.size	.L_2il0floatpacket.19,64
	.align 64
.L_2il0floatpacket.20:
	.long	0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000008,0x00000000,0x00000008,0x00000000,0x00000008,0x00000000,0x00000008,0x00000000
	.type	.L_2il0floatpacket.20,@object
	.size	.L_2il0floatpacket.20,64
	.align 64
.L_2il0floatpacket.21:
	.long	0x00000001,0x00000000,0x00000001,0x00000000,0x00000001,0x00000000,0x00000001,0x00000000,0x00000009,0x00000000,0x00000009,0x00000000,0x00000009,0x00000000,0x00000009,0x00000000
	.type	.L_2il0floatpacket.21,@object
	.size	.L_2il0floatpacket.21,64
	.align 64
.L_2il0floatpacket.22:
	.long	0x00000002,0x00000000,0x00000002,0x00000000,0x00000002,0x00000000,0x00000002,0x00000000,0x0000000e,0x00000000,0x0000000e,0x00000000,0x0000000e,0x00000000,0x0000000e,0x00000000
	.type	.L_2il0floatpacket.22,@object
	.size	.L_2il0floatpacket.22,64
	.align 64
.L_2il0floatpacket.23:
	.long	0x00000003,0x00000000,0x00000003,0x00000000,0x00000003,0x00000000,0x00000003,0x00000000,0x0000000f,0x00000000,0x0000000f,0x00000000,0x0000000f,0x00000000,0x0000000f,0x00000000
	.type	.L_2il0floatpacket.23,@object
	.size	.L_2il0floatpacket.23,64
	.align 64
.L_2il0floatpacket.24:
	.long	0x00000004,0x00000000,0x00000004,0x00000000,0x00000004,0x00000000,0x00000004,0x00000000,0x0000000c,0x00000000,0x0000000c,0x00000000,0x0000000c,0x00000000,0x0000000c,0x00000000
	.type	.L_2il0floatpacket.24,@object
	.size	.L_2il0floatpacket.24,64
	.align 64
.L_2il0floatpacket.25:
	.long	0x00000005,0x00000000,0x00000005,0x00000000,0x00000005,0x00000000,0x00000005,0x00000000,0x0000000d,0x00000000,0x0000000d,0x00000000,0x0000000d,0x00000000,0x0000000d,0x00000000
	.type	.L_2il0floatpacket.25,@object
	.size	.L_2il0floatpacket.25,64
	.align 64
.L_2il0floatpacket.26:
	.long	0x00000006,0x00000000,0x00000006,0x00000000,0x00000006,0x00000000,0x00000006,0x00000000,0x0000000a,0x00000000,0x0000000a,0x00000000,0x0000000a,0x00000000,0x0000000a,0x00000000
	.type	.L_2il0floatpacket.26,@object
	.size	.L_2il0floatpacket.26,64
	.align 64
.L_2il0floatpacket.27:
	.long	0x00000007,0x00000000,0x00000007,0x00000000,0x00000007,0x00000000,0x00000007,0x00000000,0x0000000b,0x00000000,0x0000000b,0x00000000,0x0000000b,0x00000000,0x0000000b,0x00000000
	.type	.L_2il0floatpacket.27,@object
	.size	.L_2il0floatpacket.27,64
	.align 64
.L_2il0floatpacket.28:
	.long	0x00000004,0x00000000,0x00000005,0x00000000,0x00000006,0x00000000,0x00000007,0x00000000,0x00000000,0x00000000,0x00000001,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000
	.type	.L_2il0floatpacket.28,@object
	.size	.L_2il0floatpacket.28,64
	.align 64
.L_2il0floatpacket.29:
	.long	0x0000000e,0x00000000,0x0000000f,0x00000000,0x00000008,0x00000000,0x00000009,0x00000000,0x00000006,0x00000000,0x00000007,0x00000000,0x00000000,0x00000000,0x00000001,0x00000000
	.type	.L_2il0floatpacket.29,@object
	.size	.L_2il0floatpacket.29,64
	.align 64
.L_2il0floatpacket.30:
	.long	0x0000000a,0x00000000,0x0000000b,0x00000000,0x00000004,0x00000000,0x00000005,0x00000000,0x00000008,0x00000000,0x00000009,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000
	.type	.L_2il0floatpacket.30,@object
	.size	.L_2il0floatpacket.30,64
	.align 64
.L_2il0floatpacket.31:
	.long	0x00000008,0x00000000,0x00000009,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000,0x0000000a,0x00000000,0x0000000b,0x00000000,0x00000004,0x00000000,0x00000005,0x00000000
	.type	.L_2il0floatpacket.31,@object
	.size	.L_2il0floatpacket.31,64
	.align 64
.L_2il0floatpacket.32:
	.long	0x00000000,0x00000000,0x00000001,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000,0x00000007,0x00000000,0x00000006,0x00000000,0x00000005,0x00000000,0x00000004,0x00000000
	.type	.L_2il0floatpacket.32,@object
	.size	.L_2il0floatpacket.32,64
	.align 64
.L_2il0floatpacket.33:
	.long	0x00000000,0x00000000,0x00000001,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000,0x00000006,0x00000000,0x00000007,0x00000000,0x00000004,0x00000000,0x00000005,0x00000000
	.type	.L_2il0floatpacket.33,@object
	.size	.L_2il0floatpacket.33,64
	.align 64
.L_2il0floatpacket.34:
	.long	0x00000000,0x00000000,0x00000001,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000,0x00000005,0x00000000,0x00000004,0x00000000,0x00000007,0x00000000,0x00000006,0x00000000
	.type	.L_2il0floatpacket.34,@object
	.size	.L_2il0floatpacket.34,64
	.align 64
.L_2il0floatpacket.35:
	.long	0x00000000,0x00000000,0x00000001,0x00000000,0x00000008,0x00000000,0x00000009,0x00000000,0x00000004,0x00000000,0x00000005,0x00000000,0x0000000c,0x00000000,0x0000000d,0x00000000
	.type	.L_2il0floatpacket.35,@object
	.size	.L_2il0floatpacket.35,64
	.align 64
.L_2il0floatpacket.36:
	.long	0x00000000,0x00000000,0x00000001,0x00000000,0x0000000a,0x00000000,0x0000000b,0x00000000,0x00000004,0x00000000,0x00000005,0x00000000,0x0000000e,0x00000000,0x0000000f,0x00000000
	.type	.L_2il0floatpacket.36,@object
	.size	.L_2il0floatpacket.36,64
	.align 64
.L_2il0floatpacket.37:
	.long	0x00000002,0x00000000,0x00000003,0x00000000,0x0000000a,0x00000000,0x0000000b,0x00000000,0x00000006,0x00000000,0x00000007,0x00000000,0x0000000e,0x00000000,0x0000000f,0x00000000
	.type	.L_2il0floatpacket.37,@object
	.size	.L_2il0floatpacket.37,64
	.align 64
.L_2il0floatpacket.38:
	.long	0x00000000,0x00000000,0x00000001,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000,0x00000008,0x00000000,0x00000009,0x00000000,0x0000000a,0x00000000,0x0000000b,0x00000000
	.type	.L_2il0floatpacket.38,@object
	.size	.L_2il0floatpacket.38,64
	.align 64
.L_2il0floatpacket.39:
	.long	0x00000004,0x00000000,0x00000005,0x00000000,0x00000006,0x00000000,0x00000007,0x00000000,0x0000000c,0x00000000,0x0000000d,0x00000000,0x0000000e,0x00000000,0x0000000f,0x00000000
	.type	.L_2il0floatpacket.39,@object
	.size	.L_2il0floatpacket.39,64
	.align 64
.L_2il0floatpacket.40:
	.long	0x00000007,0x00000000,0x00000006,0x00000000,0x00000005,0x00000000,0x00000004,0x00000000,0x00000007,0x00000000,0x00000006,0x00000000,0x00000005,0x00000000,0x00000004,0x00000000
	.type	.L_2il0floatpacket.40,@object
	.size	.L_2il0floatpacket.40,64
	.align 64
.L_2il0floatpacket.41:
	.long	0x00000000,0x00000000,0x00000001,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000,0x00000000,0x00000000,0x00000001,0x00000000,0x00000002,0x00000000,0x00000003,0x00000000
	.type	.L_2il0floatpacket.41,@object
	.size	.L_2il0floatpacket.41,64
	.align 64
.L_2il0floatpacket.42:
	.long	0x00000003,0x00000000,0x00000002,0x00000000,0x0000000b,0x00000000,0x0000000a,0x00000000,0x00000007,0x00000000,0x00000006,0x00000000,0x0000000f,0x00000000,0x0000000e,0x00000000
	.type	.L_2il0floatpacket.42,@object
	.size	.L_2il0floatpacket.42,64
	.align 64
.L_2il0floatpacket.43:
	.long	0x00000003,0x00000000,0x00000002,0x00000000,0x00000009,0x00000000,0x00000008,0x00000000,0x00000007,0x00000000,0x00000006,0x00000000,0x0000000d,0x00000000,0x0000000c,0x00000000
	.type	.L_2il0floatpacket.43,@object
	.size	.L_2il0floatpacket.43,64
	.align 64
.L_2il0floatpacket.44:
	.long	0x00000001,0x00000000,0x00000000,0x00000000,0x00000009,0x00000000,0x00000008,0x00000000,0x00000005,0x00000000,0x00000004,0x00000000,0x0000000d,0x00000000,0x0000000c,0x00000000
	.type	.L_2il0floatpacket.44,@object
	.size	.L_2il0floatpacket.44,64
	.align 64
.L_2il0floatpacket.45:
	.long	0x00000006,0x00000000,0x00000007,0x00000000,0x00000004,0x00000000,0x00000005,0x00000000,0x00000006,0x00000000,0x00000007,0x00000000,0x00000004,0x00000000,0x00000005,0x00000000
	.type	.L_2il0floatpacket.45,@object
	.size	.L_2il0floatpacket.45,64
	.align 64
.L_2il0floatpacket.46:
	.long	0x00000002,0x00000000,0x00000003,0x00000000,0x00000000,0x00000000,0x00000001,0x00000000,0x00000006,0x00000000,0x00000007,0x00000000,0x00000004,0x00000000,0x00000005,0x00000000
	.type	.L_2il0floatpacket.46,@object
	.size	.L_2il0floatpacket.46,64
	.align 64
.L_2il0floatpacket.47:
	.long	0x00000005,0x00000000,0x00000004,0x00000000,0x00000007,0x00000000,0x00000006,0x00000000,0x00000005,0x00000000,0x00000004,0x00000000,0x00000007,0x00000000,0x00000006,0x00000000
	.type	.L_2il0floatpacket.47,@object
	.size	.L_2il0floatpacket.47,64
	.align 64
.L_2il0floatpacket.48:
	.long	0x00000001,0x00000000,0x00000000,0x00000000,0x0000000b,0x00000000,0x0000000a,0x00000000,0x00000005,0x00000000,0x00000004,0x00000000,0x0000000f,0x00000000,0x0000000e,0x00000000
	.type	.L_2il0floatpacket.48,@object
	.size	.L_2il0floatpacket.48,64
	.data
	.section .note.GNU-stack, ""
// -- Begin DWARF2 SEGMENT .eh_frame
	.section .eh_frame,"a",@progbits
.eh_frame_seg:
	.align 8
# End
