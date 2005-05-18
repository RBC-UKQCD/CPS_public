/* 
 * Modified by Chulwoo Jung from double precision
 * BAGEL machine generated output,
 * Written by Peter Boyle
 * This software is provided for NON-COMMERCIAL use only,
 * and may not be redistributed without permission.
 * It is provided as is, and is not guaranteed fit for any purpose.
 */  
.text
.align 2
.globl asq_vaxmy_vxdot_s
asq_vaxmy_vxdot_s:

	or    %r11 , %r4 , %r4
	la   %r1,	-736(%r1)
	or    %r12 , %r5 , %r5
	stw   %r14,	368(%r1)
	stw   %r15,	376(%r1)
	or    %r14 , %r6 , %r6
	stw   %r16,	384(%r1)
	or    %r15 , %r3 , %r3
	stw   %r17,	392(%r1)
	stw   %r18,	400(%r1)
	or    %r17 , %r7 , %r7
	stw   %r19,	408(%r1)
	stw   %r20,	416(%r1)
	stw   %r21,	424(%r1)
	stw   %r22,	432(%r1)
	stw   %r23,	440(%r1)
	stw   %r24,	448(%r1)
	stw   %r25,	456(%r1)
	stw   %r26,	464(%r1)
	stw   %r27,	472(%r1)
	stw   %r28,	480(%r1)
	stw   %r29,	488(%r1)
	stw   %r30,	496(%r1)
	stw   %r31,	504(%r1)
	stfd  14,	112(%r1)
	stfd  15,	120(%r1)
	stfd  16,	128(%r1)
	stfd  17,	136(%r1)
	stfd  18,	144(%r1)
	stfd  19,	152(%r1)
	stfd  20,	160(%r1)
	stfd  21,	168(%r1)
	stfd  22,	176(%r1)
	stfd  23,	184(%r1)
	stfd  24,	192(%r1)
	stfd  25,	200(%r1)
	stfd  26,	208(%r1)
	stfd  27,	216(%r1)
	stfd  28,	224(%r1)
	stfd  29,	232(%r1)
	stfd  30,	240(%r1)
	stfd  31,	248(%r1)
	lfs   0,	0(%r15)
	li   %r18,0
	stw   %r18,	0(%r17)
	stw   %r18,	4(%r17)
asq_vaxmy_vxdot_s_lab0:
	or.  %r14 , %r14 , %r14
	lfs   19,	0(%r17)
	lfs   20,	0(%r17)
	lfs   21,	0(%r17)
	lfs   22,	0(%r17)
	lfs   23,	0(%r17)
	lfs   24,	0(%r17)
	li   %r19,32
	li   %r20,64
	bf gt,  asq_vaxmy_vxdot_s_lab1
	lfs   1,	0(%r11)
	lfs   7,	0(%r12)
	lfs   2,	4(%r11)
	lfs   8,	4(%r12)
	lfs   3,	8(%r11)
	lfs   9,	8(%r12)
	lfs   4,	12(%r11)
	lfs   10,	12(%r12)
	mtctr 	%r14
	b       asq_vaxmy_vxdot_s_lab2
.align 4
asq_vaxmy_vxdot_s_lab2:
	addi  %r16 , %r12 , 0
	lfs   5,	16(%r11)
	fmsubs  13 , 0 , 1 , 7
	lfs   11,	16(%r12)
	fmsubs  14 , 0 , 2 , 8
	lfs   6,	20(%r11)
	fmsubs  15 , 0 , 3 , 9
	lfs   12,	20(%r12)
	fmsubs  16 , 0 , 4 , 10
	la   %r11,	24(%r11)
	la   %r12,	24(%r12)
	fmadds  19 , 13 , 1 , 19
	stfs  13,	0(%r16)
	fmadds  20 , 14 , 2 , 20
	stfs  14,	4(%r16)
	fmadds  21 , 15 , 3 , 21
	stfs  15,	8(%r16)
	fmadds  22 , 16 , 4 , 22
	stfs  16,	12(%r16)
	fmsubs  17 , 0 , 5 , 11
	dcbt 	%r19,%r11
	fmsubs  18 , 0 , 6 , 12
	lfs   1,	0(%r11)
	lfs   7,	0(%r12)
	lfs   2,	4(%r11)
	lfs   8,	4(%r12)
	fmadds  23 , 17 , 5 , 23
	stfs  17,	16(%r16)
	fmadds  24 , 18 , 6 , 24
	stfs  18,	20(%r16)
	dcbt 	%r20,%r11
	lfs   3,	8(%r11)
	lfs   9,	8(%r12)
	lfs   4,	12(%r11)
	lfs   10,	12(%r12)
	la   %r14,	-1(%r14)
	dcbt 	%r19,%r12
	dcbt 	%r20,%r12
	bdnz    asq_vaxmy_vxdot_s_lab2
	fadds   19 , 19 , 20
	fadds   21 , 21 , 22
	fadds   23 , 23 , 24
	fadds   19 , 19 , 21
	fadds   19 , 19 , 23
	stfs  19,	0(%r17)
asq_vaxmy_vxdot_s_lab1:
	lwz   %r14,	368(%r1)
	lwz   %r15,	376(%r1)
	lwz   %r16,	384(%r1)
	lwz   %r17,	392(%r1)
	lwz   %r18,	400(%r1)
	lwz   %r19,	408(%r1)
	lwz   %r20,	416(%r1)
	lwz   %r21,	424(%r1)
	lwz   %r22,	432(%r1)
	lwz   %r23,	440(%r1)
	lwz   %r24,	448(%r1)
	lwz   %r25,	456(%r1)
	lwz   %r26,	464(%r1)
	lwz   %r27,	472(%r1)
	lwz   %r28,	480(%r1)
	lwz   %r29,	488(%r1)
	lwz   %r30,	496(%r1)
	lwz   %r31,	504(%r1)
	lfd   14,	112(%r1)
	lfd   15,	120(%r1)
	lfd   16,	128(%r1)
	lfd   17,	136(%r1)
	lfd   18,	144(%r1)
	lfd   19,	152(%r1)
	lfd   20,	160(%r1)
	lfd   21,	168(%r1)
	lfd   22,	176(%r1)
	lfd   23,	184(%r1)
	lfd   24,	192(%r1)
	lfd   25,	200(%r1)
	lfd   26,	208(%r1)
	lfd   27,	216(%r1)
	lfd   28,	224(%r1)
	lfd   29,	232(%r1)
	lfd   30,	240(%r1)
	lfd   31,	248(%r1)
	la   %r1,	736(%r1)
	blr
