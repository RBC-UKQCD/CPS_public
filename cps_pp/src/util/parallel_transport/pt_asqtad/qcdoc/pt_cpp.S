# Branch hints suppressed.
.set r0,0; .set SP,1; .set RTOC,2; .set r3,3; .set r4,4
.set r5,5; .set r6,6; .set r7,7; .set r8,8; .set r9,9
.set r10,10; .set r11,11; .set r12,12; .set r13,13; .set r14,14
.set r15,15; .set r16,16; .set r17,17; .set r18,18; .set r19,19
.set r20,20; .set r21,21; .set r22,22; .set r23,23; .set r24,24
.set r25,25; .set r26,26; .set r27,27; .set r28,28; .set r29,29
.set r30,30; .set r31,31
.set fp0,0; .set fp1,1; .set fp2,2; .set fp3,3; .set fp4,4
.set fp5,5; .set fp6,6; .set fp7,7; .set fp8,8; .set fp9,9
.set fp10,10; .set fp11,11; .set fp12,12; .set fp13,13; .set fp14,14
.set fp15,15; .set fp16,16; .set fp17,17; .set fp18,18; .set fp19,19
.set fp20,20; .set fp21,21; .set fp22,22; .set fp23,23; .set fp24,24
.set fp25,25; .set fp26,26; .set fp27,27; .set fp28,28; .set fp29,29
.set fp30,30; .set fp31,31
.set MQ,0; .set XER,1; .set FROM_RTCU,4; .set FROM_RTCL,5; .set FROM_DEC,6
.set LR,8; .set CTR,9; .set TID,17; .set DSISR,18; .set DAR,19; .set TO_RTCU,20
.set TO_RTCL,21; .set TO_DEC,22; .set SDR_0,24; .set SDR_1,25; .set SRR_0,26
.set SRR_1,27
.set BO_dCTR_NZERO_AND_NOT,0; .set BO_dCTR_NZERO_AND_NOT_1,1
.set BO_dCTR_ZERO_AND_NOT,2; .set BO_dCTR_ZERO_AND_NOT_1,3
.set BO_IF_NOT,4; .set BO_IF_NOT_1,5; .set BO_IF_NOT_2,6
.set BO_IF_NOT_3,7; .set BO_dCTR_NZERO_AND,8; .set BO_dCTR_NZERO_AND_1,9
.set BO_dCTR_ZERO_AND,10; .set BO_dCTR_ZERO_AND_1,11; .set BO_IF,12
.set BO_IF_1,13; .set BO_IF_2,14; .set BO_IF_3,15; .set BO_dCTR_NZERO,16
.set BO_dCTR_NZERO_1,17; .set BO_dCTR_ZERO,18; .set BO_dCTR_ZERO_1,19
.set BO_ALWAYS,20; .set BO_ALWAYS_1,21; .set BO_ALWAYS_2,22
.set BO_ALWAYS_3,23; .set BO_dCTR_NZERO_8,24; .set BO_dCTR_NZERO_9,25
.set BO_dCTR_ZERO_8,26; .set BO_dCTR_ZERO_9,27; .set BO_ALWAYS_8,28
.set BO_ALWAYS_9,29; .set BO_ALWAYS_10,30; .set BO_ALWAYS_11,31
.set CR0_LT,0; .set CR0_GT,1; .set CR0_EQ,2; .set CR0_SO,3
.set CR1_FX,4; .set CR1_FEX,5; .set CR1_VX,6; .set CR1_OX,7
.set CR2_LT,8; .set CR2_GT,9; .set CR2_EQ,10; .set CR2_SO,11
.set CR3_LT,12; .set CR3_GT,13; .set CR3_EQ,14; .set CR3_SO,15
.set CR4_LT,16; .set CR4_GT,17; .set CR4_EQ,18; .set CR4_SO,19
.set CR5_LT,20; .set CR5_GT,21; .set CR5_EQ,22; .set CR5_SO,23
.set CR6_LT,24; .set CR6_GT,25; .set CR6_EQ,26; .set CR6_SO,27
.set CR7_LT,28; .set CR7_GT,29; .set CR7_EQ,30; .set CR7_SO,31
.set TO_LT,16; .set TO_GT,8; .set TO_EQ,4; .set TO_LLT,2; .set TO_LGT,1

	.globl	cmm_agg_cpp            


# .text section
	.file	"pt_cpp.c"              


cmm_agg_cpp:                           # 0x00000000 (H.11.NO_SYMBOL)
	cmpi	0,0,r3,0
	stmw	r17,-60(SP)
	bc	BO_IF_NOT,CR0_GT,__L270
	addi	r4,r5,-144
	mulli	r8,r3,152
	addi	r9,r5,-160
	addi	r0,r0,3
	addi	r3,r5,-8
	add	r17,r8,r9
__L24:                                  # 0x00000024 (H.11.NO_SYMBOL+0x24)
	ori	r5,r3,0x0000
	lwz	r8,148(r4)
	mtspr	CTR,r0
	lwz	r9,144(r4)
	addi	r4,r4,152
	mulli	r8,r8,3
	mulli	r9,r9,3
	add	r8,r7,r8
	addi	r8,r8,-8
	add	r9,r6,r9
	addi	r10,r9,136
	addi	r11,r9,128
	addi	r12,r9,88
	addi	r31,r9,80
	addi	r30,r9,40
	addi	r29,r9,32
	addi	r28,r9,120
	addi	r27,r9,112
	addi	r26,r9,72
	addi	r25,r9,64
	addi	r24,r9,24
	addi	r23,r9,16
	addi	r22,r9,104
	addi	r21,r9,96
	addi	r20,r9,56
	addi	r19,r9,48
	addi	r18,r9,8
__L90:                                  # 0x00000090 (H.11.NO_SYMBOL+0x90)
	lfd	fp0,0(r18)
	lfd	fp1,24(r5)
	lfd	fp4,0(r9)
	fmul	fp2,fp0,fp1
	lfd	fp3,16(r5)
	lfd	fp5,0(r21)
	lfd	fp6,48(r5)
	fmsub	fp2,fp4,fp3,fp2
	lfd	fp0,0(r19)
	lfd	fp1,32(r5)
	fmadd	fp0,fp0,fp1,fp2
	lfd	fp3,0(r20)
	lfd	fp4,40(r5)
	fnmsub	fp0,fp3,fp4,fp0
	lfd	fp1,0(r22)
	lfd	fp2,56(r5)
	fmadd	fp0,fp5,fp6,fp0
	fnmsub	fp0,fp1,fp2,fp0
	stfd	fp0,8(r8)
	lfd	fp1,16(r5)
	lfd	fp4,0(r9)
	lfd	fp3,24(r5)
	lfd	fp0,0(r18)
	lfd	fp5,0(r21)
	fmul	fp2,fp1,fp0
	lfd	fp6,56(r5)
	fmadd	fp2,fp4,fp3,fp2
	lfd	fp0,0(r19)
	lfd	fp1,40(r5)
	fmadd	fp0,fp0,fp1,fp2
	lfd	fp3,0(r20)
	lfd	fp4,32(r5)
	fmadd	fp0,fp4,fp3,fp0
	lfd	fp1,0(r22)
	lfd	fp2,48(r5)
	fmadd	fp0,fp5,fp6,fp0
	fmadd	fp0,fp2,fp1,fp0
	stfd	fp0,16(r8)
	lfd	fp1,24(r5)
	lfd	fp4,0(r23)
	lfd	fp3,16(r5)
	lfd	fp0,0(r24)
	lfd	fp5,0(r27)
	fmul	fp2,fp1,fp0
	lfd	fp6,48(r5)
	fmsub	fp2,fp3,fp4,fp2
	lfd	fp0,0(r25)
	lfd	fp1,32(r5)
	fmadd	fp0,fp1,fp0,fp2
	lfd	fp3,0(r26)
	lfd	fp4,40(r5)
	fnmsub	fp0,fp4,fp3,fp0
	lfd	fp1,0(r28)
	lfd	fp2,56(r5)
	fmadd	fp0,fp6,fp5,fp0
	fnmsub	fp0,fp2,fp1,fp0
	stfd	fp0,24(r8)
	lfd	fp1,16(r5)
	lfd	fp4,0(r23)
	lfd	fp3,24(r5)
	lfd	fp0,0(r24)
	lfd	fp5,0(r27)
	fmul	fp2,fp1,fp0
	lfd	fp6,56(r5)
	fmadd	fp2,fp3,fp4,fp2
	lfd	fp0,0(r25)
	lfd	fp1,40(r5)
	fmadd	fp0,fp1,fp0,fp2
	lfd	fp3,0(r26)
	lfd	fp4,32(r5)
	fmadd	fp0,fp4,fp3,fp0
	lfd	fp1,0(r28)
	lfd	fp2,48(r5)
	fmadd	fp0,fp6,fp5,fp0
	fmadd	fp0,fp2,fp1,fp0
	stfd	fp0,32(r8)
	lfd	fp1,24(r5)
	lfd	fp4,0(r29)
	lfd	fp3,16(r5)
	lfd	fp0,0(r30)
	lfd	fp5,0(r11)
	fmul	fp2,fp1,fp0
	lfd	fp6,48(r5)
	fmsub	fp2,fp3,fp4,fp2
	lfd	fp0,0(r31)
	lfd	fp1,32(r5)
	fmadd	fp0,fp1,fp0,fp2
	lfd	fp3,0(r12)
	lfd	fp4,40(r5)
	fnmsub	fp0,fp4,fp3,fp0
	lfd	fp1,0(r10)
	lfd	fp2,56(r5)
	fmadd	fp0,fp6,fp5,fp0
	fnmsub	fp0,fp2,fp1,fp0
	stfd	fp0,40(r8)
	lfd	fp1,16(r5)
	lfd	fp4,0(r29)
	lfd	fp3,24(r5)
	lfd	fp0,0(r30)
	lfd	fp5,0(r12)
	fmul	fp2,fp1,fp0
	lfd	fp6,0(r11)
	fmadd	fp2,fp3,fp4,fp2
	lfd	fp0,0(r31)
	lfd	fp1,40(r5)
	fmadd	fp0,fp1,fp0,fp2
	lfd	fp3,32(r5)
	lfd	fp4,56(r5)
	fmadd	fp0,fp3,fp5,fp0
	lfd	fp1,0(r10)
	lfdu	fp2,48(r5)
	fmadd	fp0,fp4,fp6,fp0
	fmadd	fp0,fp2,fp1,fp0
	stfdu	fp0,48(r8)
	bc	BO_dCTR_NZERO,CR0_LT,__L90
	cmp	0,0,r3,r17
	addi	r3,r3,152
	bc	BO_IF,CR0_LT,__L24
	lmw	r17,-60(SP)
	bclr	BO_ALWAYS,CR0_LT
__L270:                                 # 0x00000270 (H.11.NO_SYMBOL+0x270)
	bclr	BO_ALWAYS,CR0_LT
