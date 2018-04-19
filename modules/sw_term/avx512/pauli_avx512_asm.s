# mark_description "Intel(R) C Intel(R) 64 Compiler for applications running on Intel(R) 64, Version 17.0.4.196 Build 20170411";
# mark_description "-I../../../include -I.. -I/cineca/prod/opt/compilers/intel/pe-xe-2017/binary/impi/2017.3.196/intel64/include";
# mark_description " -isystem /cineca/prod/opt/compilers/intel/pe-xe-2018/binary/impi/2018.1.163/include64/ -std=c89 -xCORE-AVX5";
# mark_description "12 -mtune=skylake -DAVX512 -O3 -Ddirac_counters -pedantic -fstrict-aliasing -Wno-long-long -Wstrict-prototyp";
# mark_description "es -S";
	.file "pauli_avx512.c"
	.text
..TXTST0:
# -- Begin  mul_pauli2_avx512
	.text
# mark_begin;
       .align    16,0x90
	.globl mul_pauli2_avx512
# --- mul_pauli2_avx512(float, pauli *, spinor *, spinor *)
mul_pauli2_avx512:
# parameter 1: %xmm0
# parameter 2: %rdi
# parameter 3: %rsi
# parameter 4: %rdx
..B1.1:                         # Preds ..B1.0
                                # Execution count [1.00e+00]
	.cfi_startproc
..___tag_value_mul_pauli2_avx512.1:
..L2:
                                                          #15.1
        movl      $42410, %eax                                  #82.9
        vmovups   (%rdi), %zmm8                                 #51.27
        vmovups   64(%rdi), %zmm7                               #52.27
        vmovups   144(%rdi), %zmm29                             #54.27
        vmovups   208(%rdi), %zmm3                              #55.27
        vmovups   .L_2il0floatpacket.11(%rip), %zmm25           #69.9
        vmovups   .L_2il0floatpacket.12(%rip), %zmm31           #71.9
        vmovups   .L_2il0floatpacket.13(%rip), %zmm16           #74.9
        vbroadcastss %xmm0, %xmm14                              #64.9
        vmovups   64(%rsi), %ymm12                              #43.9
        vmovups   (%rsi), %zmm11                                #42.29
        vmovups   .L_2il0floatpacket.8(%rip), %zmm9             #45.9
        vmovups   .L_2il0floatpacket.14(%rip), %zmm17           #80.10
        vmovups   .L_2il0floatpacket.15(%rip), %zmm19           #88.10
        vmovups   .L_2il0floatpacket.17(%rip), %zmm22           #101.9
        vmovups   .L_2il0floatpacket.9(%rip), %zmm10            #47.9
        vmovups   .L_2il0floatpacket.10(%rip), %zmm13           #49.9
        vmovups   .L_2il0floatpacket.16(%rip), %zmm20           #93.10
        vmovups   .L_2il0floatpacket.18(%rip), %zmm24           #107.10
        vmovups   .L_2il0floatpacket.19(%rip), %zmm28           #116.9
        vmovups   128(%rdi), %zmm2                              #53.27
        vpermi2ps %zmm29, %zmm8, %zmm25                         #69.9
        vpermi2ps %zmm3, %zmm7, %zmm31                          #71.9
        vbroadcastss %xmm14, %zmm26                             #65.9
        vpermt2ps %zmm29, %zmm28, %zmm8                         #116.9
        vpermi2ps %zmm31, %zmm25, %zmm16                        #74.9
        vpermi2ps %zmm25, %zmm31, %zmm19                        #88.10
        vpermt2ps %zmm31, %zmm22, %zmm25                        #101.9
        vpermi2ps %zmm16, %zmm26, %zmm17                        #80.10
        vpermi2ps %zmm26, %zmm31, %zmm20                        #93.10
        vpermt2ps %zmm25, %zmm24, %zmm26                        #107.10
        vmovups   .L_2il0floatpacket.22(%rip), %zmm28           #137.9
        vshufps   $170, %zmm31, %zmm8, %zmm30                   #118.10
        vshufps   $255, %zmm31, %zmm8, %zmm24                   #121.10
        kmovw     %eax, %k2                                     #82.9
        vpermilps $160, %zmm16, %zmm15                          #75.10
        movl      $23125, %eax                                  #83.9
        kmovw     %eax, %k3                                     #83.9
        movl      $21925, %eax                                  #95.9
        kmovw     %eax, %k4                                     #95.9
        vpermilps $160, %zmm25, %zmm23                          #102.10
        movl      $43610, %eax                                  #96.9
        kmovw     %eax, %k5                                     #96.9
        movl      $26022, %eax                                  #109.9
        kmovw     %eax, %k6                                     #109.9
        movl      $39513, %eax                                  #110.9
        kmovw     %eax, %k7                                     #110.9
        movl      $43690, %eax                                  #123.9
        vpermi2ps %zmm12, %zmm11, %zmm9                         #45.9
        vpermi2ps %zmm12, %zmm11, %zmm10                        #47.9
        vpermt2ps %zmm12, %zmm13, %zmm11                        #49.9
        vmulps    %zmm15, %zmm9, %zmm1                          #76.9
        vmulps    %zmm19, %zmm10, %zmm0                         #89.9
        vmulps    %zmm23, %zmm11, %zmm12                        #103.9
        vmovups   .L_2il0floatpacket.21(%rip), %zmm15           #133.9
        vpermilps $177, %zmm9, %zmm5                            #59.10
        vmulps    %zmm5, %zmm17, %zmm18                         #81.10
        vmovups   .L_2il0floatpacket.23(%rip), %zmm17           #158.9
        vpermi2ps %zmm3, %zmm7, %zmm15                          #133.9
        vaddps    %zmm18, %zmm1, %zmm1{%k2}                     #82.9
        vpermi2ps %zmm15, %zmm8, %zmm28                         #137.9
        vsubps    %zmm18, %zmm1, %zmm1{%k3}                     #83.9
        vpermi2ps %zmm15, %zmm8, %zmm17                         #158.9
        kmovw     %eax, %k3                                     #123.9
        vfmadd231ps %zmm10, %zmm30, %zmm1                       #119.9
        vpermilps $177, %zmm11, %zmm4                           #61.10
        movl      $21845, %eax                                  #124.9
        vmulps    %zmm4, %zmm26, %zmm27                         #108.10
        vmovups   .L_2il0floatpacket.20(%rip), %zmm26           #131.10
        kmovw     %eax, %k2                                     #124.9
        vpermt2ps 272(%rdi), %zmm26, %zmm2                      #131.10
        vaddps    %zmm27, %zmm12, %zmm12{%k6}                   #109.9
        vpermilps $177, %zmm10, %zmm6                           #60.10
        movl      $61680, %eax                                  #139.10
        vmulps    %zmm6, %zmm20, %zmm21                         #94.10
        vmulps    %zmm24, %zmm6, %zmm25                         #122.10
        vmovups   .L_2il0floatpacket.24(%rip), %zmm20           #169.9
        vsubps    %zmm27, %zmm12, %zmm12{%k7}                   #110.9
        vaddps    %zmm21, %zmm0, %zmm0{%k4}                     #95.9
        vpermt2ps %zmm3, %zmm20, %zmm7                          #169.9
        vaddps    %zmm25, %zmm1, %zmm1{%k3}                     #123.9
        vsubps    %zmm21, %zmm0, %zmm0{%k5}                     #96.9
        vmovups   .L_2il0floatpacket.25(%rip), %zmm3            #172.9
        vmovups   .L_2il0floatpacket.26(%rip), %zmm21           #195.9
        vsubps    %zmm25, %zmm1, %zmm1{%k2}                     #124.9
        vpermi2ps %zmm8, %zmm7, %zmm3                           #172.9
        kmovw     %eax, %k1                                     #139.10
        vpermilps $245, %zmm28, %zmm29                          #142.9
        movl      $3855, %eax                                   #148.9
        vshufps   $244, %zmm2, %zmm29, %zmm29{%k1}              #143.10
        vshufps   $68, %zmm2, %zmm7, %zmm7{%k1}                 #183.9
        kmovw     %eax, %k4                                     #148.9
        vmulps    %zmm29, %zmm5, %zmm30                         #144.10
        vpermilps $160, %zmm28, %zmm27                          #138.9
        movl      $42405, %eax                                  #154.9
        vshufps   $164, %zmm2, %zmm27, %zmm27{%k1}              #139.10
        vmovaps   %zmm15, %zmm13                                #148.9
        vshufps   $228, %zmm15, %zmm8, %zmm13{%k4}              #148.9
        vfmadd231ps %zmm9, %zmm27, %zmm0                        #140.9
        vpermilps $245, %zmm13, %zmm14                          #152.10
        vmulps    %zmm14, %zmm5, %zmm5                          #153.10
        vaddps    %zmm30, %zmm0, %zmm0{%k2}                     #145.9
        kmovw     %eax, %k2                                     #154.9
        vsubps    %zmm30, %zmm0, %zmm0{%k3}                     #146.9
        vpermilps $160, %zmm13, %zmm31                          #149.10
        movl      $23130, %eax                                  #155.9
        vpermilps $160, %zmm3, %zmm8                            #173.10
        vfmadd213ps %zmm12, %zmm31, %zmm9                       #150.9
        kmovw     %eax, %k3                                     #155.9
        vaddps    %zmm5, %zmm9, %zmm9{%k2}                      #154.9
        vpermilps $160, %zmm17, %zmm16                          #159.10
        movl      $42662, %eax                                  #164.9
        vpermilps $10, %zmm2, %zmm8{%k1}                        #174.10
        vfmadd231ps %zmm11, %zmm16, %zmm1                       #160.9
        vfmadd213ps %zmm0, %zmm8, %zmm11                        #175.9
        vsubps    %zmm5, %zmm9, %zmm9{%k3}                      #155.9
        kmovw     %eax, %k5                                     #164.9
        vpermilps $245, %zmm3, %zmm0                            #177.10
        movl      $22873, %eax                                  #165.9
        vpermilps $245, %zmm17, %zmm18                          #162.10
        vpermilps $95, %zmm2, %zmm0{%k1}                        #178.10
        vpermilps $160, %zmm7, %zmm2                            #184.10
        vpermilps $245, %zmm7, %zmm7                            #187.10
        vmulps    %zmm18, %zmm4, %zmm19                         #163.10
        vmulps    %zmm7, %zmm6, %zmm6                           #188.10
        vmulps    %zmm0, %zmm4, %zmm4                           #179.10
        vfmadd213ps %zmm9, %zmm2, %zmm10                        #185.9
        vmovups   .L_2il0floatpacket.27(%rip), %zmm9            #197.9
        vaddps    %zmm19, %zmm1, %zmm1{%k5}                     #164.9
        vaddps    %zmm6, %zmm10, %zmm10{%k2}                    #189.9
        kmovw     %eax, %k6                                     #165.9
        vsubps    %zmm6, %zmm10, %zmm10{%k3}                    #190.9
        vsubps    %zmm19, %zmm1, %zmm1{%k6}                     #165.9
        movl      $25957, %eax                                  #180.9
        kmovw     %eax, %k7                                     #180.9
        movl      $39578, %eax                                  #181.9
        kmovw     %eax, %k5                                     #181.9
        vaddps    %zmm4, %zmm11, %zmm11{%k7}                    #180.9
        vpermi2ps %zmm10, %zmm1, %zmm21                         #195.9
        vpermt2ps %zmm10, %zmm9, %zmm1                          #197.9
        vsubps    %zmm4, %zmm11, %zmm11{%k5}                    #181.9
        vaddps    %zmm1, %zmm21, %zmm1                          #198.9
        vextractf32x4 $1, %zmm11, %xmm10                        #204.11
        vextractf32x4 $2, %zmm11, %xmm0                         #211.11
        vextractf32x4 $3, %zmm11, %xmm23                        #212.11
        vaddps    %xmm11, %xmm10, %xmm22                        #205.11
        vmovups   %xmm22, 32(%rdx)                              #206.21
        vmovups   %ymm1, (%rdx)                                 #201.24
        vextractf64x4 $1, %zmm1, 48(%rdx)                       #209.24
        vaddps    %xmm0, %xmm23, %xmm1                          #213.11
        vmovups   %xmm1, 80(%rdx)                               #214.21
        vzeroupper                                              #215.1
        ret                                                     #215.1
        .align    16,0x90
                                # LOE
	.cfi_endproc
# mark_end;
	.type	mul_pauli2_avx512,@function
	.size	mul_pauli2_avx512,.-mul_pauli2_avx512
	.data
# -- End  mul_pauli2_avx512
	.section .rodata, "a"
	.align 64
	.align 64
.L_2il0floatpacket.8:
	.long	0x00000000,0x00000001,0x00000002,0x00000003,0x00000006,0x00000007,0x00000008,0x00000009,0x0000000c,0x0000000d,0x0000000e,0x0000000f,0x00000012,0x00000013,0x00000014,0x00000015
	.type	.L_2il0floatpacket.8,@object
	.size	.L_2il0floatpacket.8,64
	.align 64
.L_2il0floatpacket.9:
	.long	0x00000002,0x00000003,0x00000004,0x00000005,0x00000008,0x00000009,0x0000000a,0x0000000b,0x0000000e,0x0000000f,0x00000010,0x00000011,0x00000014,0x00000015,0x00000016,0x00000017
	.type	.L_2il0floatpacket.9,@object
	.size	.L_2il0floatpacket.9,64
	.align 64
.L_2il0floatpacket.10:
	.long	0x00000004,0x00000005,0x00000000,0x00000001,0x0000000a,0x0000000b,0x00000006,0x00000007,0x00000010,0x00000011,0x0000000c,0x0000000d,0x00000016,0x00000017,0x00000012,0x00000013
	.type	.L_2il0floatpacket.10,@object
	.size	.L_2il0floatpacket.10,64
	.align 64
.L_2il0floatpacket.11:
	.long	0x00000000,0x00000001,0x0000000a,0x0000000b,0x00000004,0x00000005,0x00000002,0x00000003,0x00000010,0x00000011,0x0000001a,0x0000001b,0x00000014,0x00000015,0x00000012,0x00000013
	.type	.L_2il0floatpacket.11,@object
	.size	.L_2il0floatpacket.11,64
	.align 64
.L_2il0floatpacket.12:
	.long	0x00000004,0x00000005,0x00000000,0x00000001,0x0000000c,0x0000000d,0x00000006,0x00000007,0x00000014,0x00000015,0x00000010,0x00000011,0x0000001c,0x0000001d,0x00000016,0x00000017
	.type	.L_2il0floatpacket.12,@object
	.size	.L_2il0floatpacket.12,64
	.align 64
.L_2il0floatpacket.13:
	.long	0x00000000,0x00000000,0x00000001,0x00000001,0x00000002,0x00000003,0x00000010,0x00000011,0x00000008,0x00000008,0x00000009,0x00000009,0x0000000a,0x0000000b,0x00000018,0x00000019
	.type	.L_2il0floatpacket.13,@object
	.size	.L_2il0floatpacket.13,64
	.align 64
.L_2il0floatpacket.14:
	.long	0x00000000,0x00000001,0x00000002,0x00000003,0x00000015,0x00000015,0x00000017,0x00000017,0x00000008,0x00000009,0x0000000a,0x0000000b,0x0000001d,0x0000001d,0x0000001f,0x0000001f
	.type	.L_2il0floatpacket.14,@object
	.size	.L_2il0floatpacket.14,64
	.align 64
.L_2il0floatpacket.15:
	.long	0x00000000,0x00000000,0x00000004,0x00000004,0x00000014,0x00000014,0x00000015,0x00000015,0x00000008,0x00000008,0x0000000c,0x0000000c,0x0000001c,0x0000001c,0x0000001d,0x0000001d
	.type	.L_2il0floatpacket.15,@object
	.size	.L_2il0floatpacket.15,64
	.align 64
.L_2il0floatpacket.16:
	.long	0x00000001,0x00000001,0x00000005,0x00000005,0x00000014,0x00000015,0x00000016,0x00000017,0x00000009,0x00000009,0x0000000d,0x0000000d,0x0000001c,0x0000001d,0x0000001e,0x0000001f
	.type	.L_2il0floatpacket.16,@object
	.size	.L_2il0floatpacket.16,64
	.align 64
.L_2il0floatpacket.17:
	.long	0x00000006,0x00000006,0x00000002,0x00000003,0x00000014,0x00000015,0x00000007,0x00000007,0x0000000e,0x0000000e,0x0000000a,0x0000000b,0x0000001c,0x0000001d,0x0000000f,0x0000000f
	.type	.L_2il0floatpacket.17,@object
	.size	.L_2il0floatpacket.17,64
	.align 64
.L_2il0floatpacket.18:
	.long	0x00000000,0x00000001,0x00000013,0x00000013,0x00000015,0x00000015,0x00000006,0x00000007,0x00000008,0x00000009,0x0000001b,0x0000001b,0x0000001d,0x0000001d,0x0000000e,0x0000000f
	.type	.L_2il0floatpacket.18,@object
	.size	.L_2il0floatpacket.18,64
	.align 64
.L_2il0floatpacket.19:
	.long	0x00000008,0x00000009,0x00000006,0x00000007,0x0000000e,0x0000000f,0x0000000c,0x0000000d,0x00000018,0x00000019,0x00000016,0x00000017,0x0000001e,0x0000001f,0x0000001c,0x0000001d
	.type	.L_2il0floatpacket.19,@object
	.size	.L_2il0floatpacket.19,64
	.align 64
.L_2il0floatpacket.20:
	.long	0x00000000,0x00000001,0x00000002,0x00000003,0x00000000,0x00000001,0x00000002,0x00000003,0x00000010,0x00000011,0x00000012,0x00000013,0x00000010,0x00000011,0x00000012,0x00000013
	.type	.L_2il0floatpacket.20,@object
	.size	.L_2il0floatpacket.20,64
	.align 64
.L_2il0floatpacket.21:
	.long	0x00000006,0x00000007,0x00000002,0x00000003,0x00000008,0x00000009,0x0000000e,0x0000000f,0x00000016,0x00000017,0x00000012,0x00000013,0x00000018,0x00000019,0x0000001e,0x0000001f
	.type	.L_2il0floatpacket.21,@object
	.size	.L_2il0floatpacket.21,64
	.align 64
.L_2il0floatpacket.22:
	.long	0x00000006,0x00000007,0x00000010,0x00000011,0x00000016,0x00000017,0x00000000,0x00000000,0x0000000e,0x0000000f,0x00000018,0x00000019,0x0000001e,0x0000001f,0x00000000,0x00000000
	.type	.L_2il0floatpacket.22,@object
	.size	.L_2il0floatpacket.22,64
	.align 64
.L_2il0floatpacket.23:
	.long	0x00000000,0x00000001,0x00000002,0x00000003,0x00000004,0x00000005,0x00000012,0x00000013,0x00000008,0x00000009,0x0000000a,0x0000000b,0x0000000c,0x0000000d,0x0000001a,0x0000001b
	.type	.L_2il0floatpacket.23,@object
	.size	.L_2il0floatpacket.23,64
	.align 64
.L_2il0floatpacket.24:
	.long	0x00000000,0x00000001,0x00000008,0x00000009,0x0000000a,0x0000000b,0xffffffff,0xffffffff,0x00000010,0x00000011,0x00000018,0x00000019,0x0000001a,0x0000001b,0xffffffff,0xffffffff
	.type	.L_2il0floatpacket.24,@object
	.size	.L_2il0floatpacket.24,64
	.align 64
.L_2il0floatpacket.25:
	.long	0x00000004,0x00000005,0x00000014,0x00000015,0x00000000,0x00000000,0x00000000,0x00000000,0x0000000c,0x0000000d,0x0000001c,0x0000001d,0x00000000,0x00000000,0x00000000,0x00000000
	.type	.L_2il0floatpacket.25,@object
	.size	.L_2il0floatpacket.25,64
	.align 64
.L_2il0floatpacket.26:
	.long	0x00000000,0x00000001,0x00000002,0x00000003,0x00000010,0x00000011,0x00000012,0x00000013,0x00000008,0x00000009,0x0000000a,0x0000000b,0x00000018,0x00000019,0x0000001a,0x0000001b
	.type	.L_2il0floatpacket.26,@object
	.size	.L_2il0floatpacket.26,64
	.align 64
.L_2il0floatpacket.27:
	.long	0x00000004,0x00000005,0x00000006,0x00000007,0x00000014,0x00000015,0x00000016,0x00000017,0x0000000c,0x0000000d,0x0000000e,0x0000000f,0x0000001c,0x0000001d,0x0000001e,0x0000001f
	.type	.L_2il0floatpacket.27,@object
	.size	.L_2il0floatpacket.27,64
	.data
	.section .note.GNU-stack, ""
// -- Begin DWARF2 SEGMENT .eh_frame
	.section .eh_frame,"a",@progbits
.eh_frame_seg:
	.align 8
# End
