	.section	".text"
	.align 2
	.globl p4_dsum_s
	.type	p4_dsum_s,@function
p4_dsum_s:
	la   %r1,	-864(%r1)
	stw	%r14,	368(%r1)
	stw	%r15,	376(%r1)
	stw	%r16,	384(%r1)
	stw	%r17,	392(%r1)
	stw	%r18,	400(%r1)
	stw	%r19,	408(%r1)
	stw	%r20,	416(%r1)
	stw	%r21,	424(%r1)
	stw	%r22,	432(%r1)
	stw	%r23,	440(%r1)
	stw	%r24,	448(%r1)
	stw	%r25,	456(%r1)
	stw	%r26,	464(%r1)
	stw	%r27,	472(%r1)
	stw	%r28,	480(%r1)
	stw	%r29,	488(%r1)
	stw	%r30,	496(%r1)
	stw	%r31,	504(%r1)
	li   16,	0
	li   17,	32
	li   18,	64
	li   19,	96
	li   20,	128
	li   21,	160
	li   22,	192
	li   23,	224
	li   24,	256
	li   25,	288
	li   26,	320
	li   27,	352
	li   28,	384
	li   29,	416
	li   30,	448
	li   31,	480
	stfd	14,	112(%r1)
	stfd	15,	120(%r1)
	stfd	16,	128(%r1)
	stfd	17,	136(%r1)
	stfd	18,	144(%r1)
	stfd	19,	152(%r1)
	stfd	20,	160(%r1)
	stfd	21,	168(%r1)
	stfd	22,	176(%r1)
	stfd	23,	184(%r1)
	stfd	24,	192(%r1)
	stfd	25,	200(%r1)
	stfd	26,	208(%r1)
	stfd	27,	216(%r1)
	stfd	28,	224(%r1)
	stfd	29,	232(%r1)
	stfd	30,	240(%r1)
	stfd	31,	248(%r1)
	lfs	   30,	0(%r6)
	lfs	   31,	0(%r7)
	addi %r11, %r5, -24
	addi %r10, %r4, -256
	addi %r8, %r4, 256
	addi %r9, %r4, 0
	dcbt	 %r16,	 %r4
	dcbt	 %r17,	 %r4
	dcbt	 %r18,	 %r4
	dcbt	 %r19,	 %r4
	dcbt	 %r16,	 %r10
	dcbt	 %r17,	 %r10
	dcbt	 %r18,	 %r10
	dcbt	 %r19,	 %r10








	or.	 3,	3,	3
	bf gt,	p4_dsum_s_lab0
	lfs	   12,	0(%r4)
	lfs	   13,	4(%r4)
	lfs	   14,	8(%r4)
	lfs	   15,	12(%r4)
	lfs	   16,	16(%r4)
	lfs	   17,	20(%r4)
	lfs	   0,	0(%r5)
	lfs	   1,	4(%r5)
	lfs	   2,	8(%r5)
	lfs	   3,	12(%r5)
	lfs	   4,	16(%r5)
	lfs	   5,	20(%r5)
	subic.	 3,	3,	1
	bf gt,	p4_dsum_s_lab2
	b  p4_dsum_s_lab1
p4_dsum_s_lab1:
	lfs	   18,	32(%r4)
	dcbt	 %r8,	 %r16
	fmadds	0,	30,	12,	0
	lfs	   19,	36(%r4)

	fmadds	1,	30,	13,	1
	lfs	   20,	40(%r4)

	fmadds	2,	30,	14,	2
	lfs	   21,	44(%r4)

	fmadds	3,	30,	15,	3
	lfs	   22,	48(%r4)
	addi %r5, %r5, 24 
	fmadds	4,	30,	16,	4
	lfs	   23,	52(%r4)
	addi %r11, %r11, 24 
	fmadds	5,	30,	17,	5
	lfs	   24,	64(%r4)
	dcbt	 %r8,	 %r17
	fmadds	0,	31,	18,	0
	lfs	   25,	68(%r4)

	fmadds	1,	31,	19,	1
	lfs	   26,	72(%r4)

	fmadds	2,	31,	20,	2
	lfs	   27,	76(%r4)

	fmadds	3,	31,	21,	3
	lfs	   28,	80(%r4)
	addi %r9, %r9, 256 
	fmadds	4,	31,	22,	4
	lfs	   29,	84(%r4)
	addi %r10, %r10, 256 
	fmadds	5,	31,	23,	5
	lfs	   12,	96(%r4)
	dcbt	 %r8,	 %r18
	fmadds	0,	30,	24,	0
	lfs	   13,	100(%r4)

	fmadds	1,	30,	25,	1
	lfs	   14,	104(%r4)

	fmadds	2,	30,	26,	2
	lfs	   15,	108(%r4)
	dcbt	 %r5,	 %r16
	fmadds	3,	30,	27,	3
	lfs	   16,	112(%r4)
	dcbt	 %r5,	 %r17
	fmadds	4,	30,	28,	4
	lfs	   17,	116(%r4)

	fmadds	5,	30,	29,	5
	lfs	   18,	128(%r4)
	dcbt	 %r8,	 %r19
	fmadds	0,	31,	12,	0
	lfs	   19,	132(%r4)

	fmadds	1,	31,	13,	1
	lfs	   20,	136(%r4)

	fmadds	2,	31,	14,	2
	lfs	   21,	140(%r4)
	dcbt	 %r8,	 %r20
	fmadds	3,	31,	15,	3
	lfs	   22,	144(%r4)

	fmadds	4,	31,	16,	4
	lfs	   23,	148(%r4)

	fmadds	5,	31,	17,	5
	lfs	   24,	160(%r4)
	dcbt	 %r9,	 %r16
	subic.	 3,	3,	1
	lfs	   25,	164(%r4)

	fmadds	6,	30,	18,	0
	lfs	   26,	168(%r4)

	fmadds	7,	30,	19,	1
	lfs	   27,	172(%r4)

	fmadds	8,	30,	20,	2
	lfs	   28,	176(%r4)

	fmadds	9,	30,	21,	3
	lfs	   29,	180(%r4)
	dcbt	 %r9,	 %r17
	fmadds	10,	30,	22,	4
	lfs	   12,	192(%r4)

	fmadds	11,	30,	23,	5
	lfs	   13,	196(%r4)

	fmadds	6,	31,	24,	6
	lfs	   14,	200(%r4)

	fmadds	7,	31,	25,	7
	lfs	   15,	204(%r4)

	fmadds	8,	31,	26,	8
	lfs	   16,	208(%r4)
	dcbt	 %r9,	 %r18
	fmadds	9,	31,	27,	9
	lfs	   17,	212(%r4)

	fmadds	10,	31,	28,	10
	lfs	   18,	224(%r10)

	fmadds	11,	31,	29,	11
	lfs	   19,	228(%r10)

	fmadds	6,	30,	12,	6
	lfs	   20,	232(%r10)

	fmadds	7,	30,	13,	7
	lfs	   21,	236(%r10)
	dcbt	 %r9,	 %r19
	fmadds	8,	30,	14,	8
	lfs	   22,	240(%r10)
	addi %r4, %r4, 256 
	fmadds	9,	30,	15,	9
	lfs	   23,	244(%r10)
	addi %r8, %r8, 256 
	fmadds	10,	30,	16,	10
	lfs	   12,	256(%r10)
	dcbt	 %r9,	 %r20
	fmadds	11,	30,	17,	11
	lfs	   13,	260(%r10)

	fmadds	6,	31,	18,	6
	lfs	   14,	264(%r10)

	fmadds	7,	31,	19,	7
	lfs	   15,	268(%r10)

	fmadds	8,	31,	20,	8
	lfs	   16,	272(%r10)

	fmadds	9,	31,	21,	9
	lfs	   17,	276(%r10)

	fmadds	10,	31,	22,	10
	stfs	   6,	0(%r11)
	stfs	   7,	4(%r11)
	fmadds	11,	31,	23,	11
	stfs	   8,	8(%r11)
	stfs	   9,	12(%r11)
	stfs	   10,	16(%r11)
	stfs	   11,	20(%r11)
	lfs	   0,	0(%r5)
	lfs	   1,	4(%r5)
	lfs	   2,	8(%r5)
	lfs	   3,	12(%r5)
	lfs	   4,	16(%r5)
	lfs	   5,	20(%r5)
	bf gt,	p4_dsum_s_lab2
	b  p4_dsum_s_lab1
p4_dsum_s_lab2:
	fmadds	0,	30,	12,	0
	lfs	   18,	32(%r4)
	fmadds	1,	30,	13,	1
	lfs	   19,	36(%r4)
	fmadds	2,	30,	14,	2
	lfs	   20,	40(%r4)
	fmadds	3,	30,	15,	3
	lfs	   21,	44(%r4)
	fmadds	4,	30,	16,	4
	lfs	   22,	48(%r4)
	fmadds	5,	30,	17,	5
	lfs	   23,	52(%r4)
	fmadds	0,	31,	18,	0
	addi %r11, %r11, 24 
	fmadds	1,	31,	19,	1
	addi %r10, %r10, 256 
	fmadds	2,	31,	20,	2
	lfs	   24,	64(%r4)
	fmadds	3,	31,	21,	3
	lfs	   25,	68(%r4)
	fmadds	4,	31,	22,	4
	lfs	   26,	72(%r4)
	fmadds	5,	31,	23,	5
	lfs	   27,	76(%r4)
	fmadds	0,	30,	24,	0
	lfs	   28,	80(%r4)
	fmadds	1,	30,	25,	1
	lfs	   29,	84(%r4)
	fmadds	2,	30,	26,	2
	lfs	   12,	96(%r4)
	fmadds	3,	30,	27,	3
	lfs	   13,	100(%r4)
	fmadds	4,	30,	28,	4
	lfs	   14,	104(%r4)
	fmadds	5,	30,	29,	5
	lfs	   15,	108(%r4)
	fmadds	0,	31,	12,	0
	lfs	   16,	112(%r4)
	fmadds	1,	31,	13,	1
	lfs	   17,	116(%r4)
	fmadds	2,	31,	14,	2
	lfs	   18,	128(%r4)
	fmadds	3,	31,	15,	3
	lfs	   19,	132(%r4)
	fmadds	4,	31,	16,	4
	lfs	   20,	136(%r4)
	fmadds	5,	31,	17,	5
	lfs	   21,	140(%r4)
	subic.	 3,	3,	1
	lfs	   22,	144(%r4)
	fmadds	6,	30,	18,	0
	lfs	   23,	148(%r4)
	fmadds	7,	30,	19,	1
	lfs	   24,	160(%r4)
	fmadds	8,	30,	20,	2
	lfs	   25,	164(%r4)
	fmadds	9,	30,	21,	3
	lfs	   26,	168(%r4)
	fmadds	10,	30,	22,	4
	lfs	   27,	172(%r4)
	fmadds	11,	30,	23,	5
	lfs	   28,	176(%r4)
	fmadds	6,	31,	24,	6
	lfs	   29,	180(%r4)
	fmadds	7,	31,	25,	7
	lfs	   12,	192(%r4)
	fmadds	8,	31,	26,	8
	lfs	   13,	196(%r4)
	fmadds	9,	31,	27,	9
	lfs	   14,	200(%r4)
	fmadds	10,	31,	28,	10
	lfs	   15,	204(%r4)
	fmadds	11,	31,	29,	11
	lfs	   16,	208(%r4)
	fmadds	6,	30,	12,	6
	lfs	   17,	212(%r4)
	fmadds	7,	30,	13,	7
	lfs	   18,	224(%r10)
	fmadds	8,	30,	14,	8
	lfs	   19,	228(%r10)
	fmadds	9,	30,	15,	9
	lfs	   20,	232(%r10)
	fmadds	10,	30,	16,	10
	lfs	   21,	236(%r10)
	fmadds	11,	30,	17,	11
	lfs	   22,	240(%r10)
	fmadds	6,	31,	18,	6
	lfs	   23,	244(%r10)
	fmadds	7,	31,	19,	7

	fmadds	8,	31,	20,	8

	fmadds	9,	31,	21,	9

	fmadds	10,	31,	22,	10

	fmadds	11,	31,	23,	11


	stfs	   6,	0(%r11)
	stfs	   7,	4(%r11)
	stfs	   8,	8(%r11)
	stfs	   9,	12(%r11)
	stfs	   10,	16(%r11)
	stfs	   11,	20(%r11)
p4_dsum_s_lab0:
	lwz	%r14,	368(%r1)
	lwz	%r15,	376(%r1)
	lwz	%r16,	384(%r1)
	lwz	%r17,	392(%r1)
	lwz	%r18,	400(%r1)
	lwz	%r19,	408(%r1)
	lwz	%r20,	416(%r1)
	lwz	%r21,	424(%r1)
	lwz	%r22,	432(%r1)
	lwz	%r23,	440(%r1)
	lwz	%r24,	448(%r1)
	lwz	%r25,	456(%r1)
	lwz	%r26,	464(%r1)
	lwz	%r27,	472(%r1)
	lwz	%r28,	480(%r1)
	lwz	%r29,	488(%r1)
	lwz	%r30,	496(%r1)
	lwz	%r31,	504(%r1)
	lfd	14,	112(%r1)
	lfd	15,	120(%r1)
	lfd	16,	128(%r1)
	lfd	17,	136(%r1)
	lfd	18,	144(%r1)
	lfd	19,	152(%r1)
	lfd	20,	160(%r1)
	lfd	21,	168(%r1)
	lfd	22,	176(%r1)
	lfd	23,	184(%r1)
	lfd	24,	192(%r1)
	lfd	25,	200(%r1)
	lfd	26,	208(%r1)
	lfd	27,	216(%r1)
	lfd	28,	224(%r1)
	lfd	29,	232(%r1)
	lfd	30,	240(%r1)
	lfd	31,	248(%r1)
	la   %r1,	864(%r1)
	blr
