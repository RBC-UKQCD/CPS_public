	.section	".text"
	.align 2
	.globl asq_dsum_s
	.type	asq_dsum_s,@function
asq_dsum_s:
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
	addi %r6, %r4, -512
	addi %r7, %r4, 256
	addi %r8, %r4, 0
	dcbt	 %r16,	 %r4
	dcbt	 %r17,	 %r4
	dcbt	 %r18,	 %r4
	dcbt	 %r19,	 %r4
	or.	 3,	3,	3
	addi %r5, %r5, -24 
	bf gt,	lab0
	lfs	   6,	0(%r4)
	lfs	   7,	4(%r4)
	lfs	   8,	8(%r4)
	lfs	   9,	12(%r4)
	lfs	   10,	16(%r4)
	lfs	   11,	20(%r4)
	lfs	   12,	32(%r4)
	lfs	   13,	36(%r4)
	lfs	   14,	40(%r4)
	lfs	   15,	44(%r4)
	lfs	   16,	48(%r4)
	lfs	   17,	52(%r4)
	subic.	 3,	3,	1
	bf gt,	lab2
	b  lab1
lab1:
	lfs	   18,	64(%r4)
	dcbt	 %r7,	 %r16
	fadds	   0,	6,	12
	lfs	   19,	68(%r4)

	fadds	   1,	7,	13
	lfs	   20,	72(%r4)

	fadds	   2,	8,	14
	lfs	   21,	76(%r4)
	dcbt	 %r7,	 %r17
	fadds	   3,	9,	15
	lfs	   22,	80(%r4)

	fadds	   4,	10,	16
	lfs	   23,	84(%r4)

	fadds	   5,	11,	17
	lfs	   24,	96(%r4)
	dcbt	 %r7,	 %r18
	fadds	   0,	0,	18
	lfs	   25,	100(%r4)

	fadds	   1,	1,	19
	lfs	   26,	104(%r4)

	fadds	   2,	2,	20
	lfs	   27,	108(%r4)
	dcbt	 %r7,	 %r19
	fadds	   3,	3,	21
	lfs	   28,	112(%r4)

	fadds	   4,	4,	22
	lfs	   29,	116(%r4)

	fadds	   5,	5,	23
	lfs	   6,	128(%r4)

	fadds	   0,	0,	24
	lfs	   7,	132(%r4)

	fadds	   1,	1,	25
	lfs	   8,	136(%r4)

	fadds	   2,	2,	26
	lfs	   9,	140(%r4)

	fadds	   3,	3,	27
	lfs	   10,	144(%r4)

	fadds	   4,	4,	28
	lfs	   11,	148(%r4)

	fadds	   5,	5,	29
	lfs	   12,	160(%r4)

	fadds	   0,	0,	6
	lfs	   13,	164(%r4)

	fadds	   1,	1,	7
	lfs	   14,	168(%r4)

	fadds	   2,	2,	8
	lfs	   15,	172(%r4)

	fadds	   3,	3,	9
	lfs	   16,	176(%r4)
	addi %r8, %r8, 512 
	fadds	   4,	4,	10
	lfs	   17,	180(%r4)
	addi %r6, %r6, 512 
	fadds	   5,	5,	11
	lfs	   18,	192(%r4)
	dcbt	 %r7,	 %r20
	fadds	   0,	0,	12
	lfs	   19,	196(%r4)

	fadds	   1,	1,	13
	lfs	   20,	200(%r4)

	fadds	   2,	2,	14
	lfs	   21,	204(%r4)

	fadds	   3,	3,	15
	lfs	   22,	208(%r4)

	fadds	   4,	4,	16
	lfs	   23,	212(%r4)

	fadds	   5,	5,	17
	lfs	   24,	224(%r4)
	dcbt	 %r7,	 %r21
	fadds	   0,	0,	18
	lfs	   25,	228(%r4)

	fadds	   1,	1,	19
	lfs	   26,	232(%r4)

	fadds	   2,	2,	20
	lfs	   27,	236(%r4)

	fadds	   3,	3,	21
	lfs	   28,	240(%r4)

	fadds	   4,	4,	22
	lfs	   29,	244(%r4)

	fadds	   5,	5,	23
	lfs	   6,	256(%r4)
	dcbt	 %r7,	 %r22
	fadds	   0,	0,	24
	lfs	   7,	260(%r4)

	fadds	   1,	1,	25
	lfs	   8,	264(%r4)

	fadds	   2,	2,	26
	lfs	   9,	268(%r4)

	fadds	   3,	3,	27
	lfs	   10,	272(%r4)

	fadds	   4,	4,	28
	lfs	   11,	276(%r4)

	fadds	   5,	5,	29
	lfs	   12,	288(%r4)
	dcbt	 %r7,	 %r23
	subic.	 3,	3,	1
	lfs	   13,	292(%r4)

	fadds	   0,	0,	6
	lfs	   14,	296(%r4)

	fadds	   1,	1,	7
	lfs	   15,	300(%r4)

	fadds	   2,	2,	8
	lfs	   16,	304(%r4)

	fadds	   3,	3,	9
	lfs	   17,	308(%r4)

	fadds	   4,	4,	10
	addi %r5, %r5, 24 
	dcbt	 %r8,	 %r16
	fadds	   5,	5,	11
	lfs	   18,	320(%r6)

	fadds	   0,	0,	12
	lfs	   19,	324(%r6)

	fadds	   1,	1,	13
	lfs	   20,	328(%r6)
	dcbt	 %r8,	 %r17
	fadds	   2,	2,	14
	lfs	   21,	332(%r6)

	fadds	   3,	3,	15
	lfs	   22,	336(%r6)

	fadds	   4,	4,	16
	lfs	   23,	340(%r6)
	dcbt	 %r8,	 %r18
	fadds	   5,	5,	17
	lfs	   24,	352(%r6)

	fadds	   0,	0,	18
	lfs	   25,	356(%r6)

	fadds	   1,	1,	19
	lfs	   26,	360(%r6)
	dcbt	 %r8,	 %r19
	fadds	   2,	2,	20
	lfs	   27,	364(%r6)

	fadds	   3,	3,	21
	lfs	   28,	368(%r6)

	fadds	   4,	4,	22
	lfs	   29,	372(%r6)

	fadds	   5,	5,	23
	lfs	   6,	384(%r6)

	fadds	   0,	0,	24
	lfs	   7,	388(%r6)

	fadds	   1,	1,	25
	lfs	   8,	392(%r6)

	fadds	   2,	2,	26
	lfs	   9,	396(%r6)

	fadds	   3,	3,	27
	lfs	   10,	400(%r6)

	fadds	   4,	4,	28
	lfs	   11,	404(%r6)

	fadds	   5,	5,	29
	lfs	   12,	416(%r6)

	fadds	   0,	0,	6
	lfs	   13,	420(%r6)

	fadds	   1,	1,	7
	lfs	   14,	424(%r6)

	fadds	   2,	2,	8
	lfs	   15,	428(%r6)
	addi %r4, %r4, 512 
	fadds	   3,	3,	9
	lfs	   16,	432(%r6)
	addi %r7, %r7, 512 
	fadds	   4,	4,	10
	lfs	   17,	436(%r6)
	dcbt	 %r8,	 %r20
	fadds	   5,	5,	11
	lfs	   18,	448(%r6)

	fadds	   0,	0,	12
	lfs	   19,	452(%r6)

	fadds	   1,	1,	13
	lfs	   20,	456(%r6)

	fadds	   2,	2,	14
	lfs	   21,	460(%r6)

	fadds	   3,	3,	15
	lfs	   22,	464(%r6)

	fadds	   4,	4,	16
	lfs	   23,	468(%r6)
	dcbt	 %r8,	 %r21
	fadds	   5,	5,	17
	lfs	   24,	480(%r6)

	fadds	   0,	0,	18
	lfs	   25,	484(%r6)

	fadds	   1,	1,	19
	lfs	   26,	488(%r6)

	fadds	   2,	2,	20
	lfs	   27,	492(%r6)

	fadds	   3,	3,	21
	lfs	   28,	496(%r6)

	fadds	   4,	4,	22
	lfs	   29,	500(%r6)
	dcbt	 %r8,	 %r22
	fadds	   5,	5,	23
	lfs	   6,	512(%r6)

	fadds	   0,	0,	24
	lfs	   7,	516(%r6)

	fadds	   1,	1,	25
	lfs	   8,	520(%r6)

	fadds	   2,	2,	26
	lfs	   9,	524(%r6)

	fadds	   3,	3,	27
	lfs	   10,	528(%r6)

	fadds	   4,	4,	28
	lfs	   11,	532(%r6)
	dcbt	 %r8,	 %r23
	fadds	   5,	5,	29
#<<<<<< Break 1

	lfs	   12,	544(%r6)

	lfs	   13,	548(%r6)

	lfs	   14,	552(%r6)

	lfs	   15,	556(%r6)

	lfs	   16,	560(%r6)
	lfs	   17,	564(%r6)
	stfs	   0,	0(%r5)
	stfs	   1,	4(%r5)
	stfs	   2,	8(%r5)
	stfs	   3,	12(%r5)
	stfs	   4,	16(%r5)
	stfs	   5,	20(%r5)
	bf gt,	lab2
	b  lab1
lab2:
	fadds	   0,	6,	12
	lfs	   18,	64(%r4)
	fadds	   1,	7,	13
	lfs	   19,	68(%r4)
	fadds	   2,	8,	14
	lfs	   20,	72(%r4)
	fadds	   3,	9,	15
	lfs	   21,	76(%r4)
	fadds	   4,	10,	16
	lfs	   22,	80(%r4)
	fadds	   5,	11,	17
	lfs	   23,	84(%r4)
	fadds	   0,	0,	18
	lfs	   24,	96(%r4)
	fadds	   1,	1,	19
	lfs	   25,	100(%r4)
	fadds	   2,	2,	20
	lfs	   26,	104(%r4)
	fadds	   3,	3,	21
	lfs	   27,	108(%r4)
	fadds	   4,	4,	22
	lfs	   28,	112(%r4)
	fadds	   5,	5,	23
	lfs	   29,	116(%r4)
	fadds	   0,	0,	24
	lfs	   6,	128(%r4)
	fadds	   1,	1,	25
	lfs	   7,	132(%r4)
	fadds	   2,	2,	26
	lfs	   8,	136(%r4)
	fadds	   3,	3,	27
	lfs	   9,	140(%r4)
	fadds	   4,	4,	28
	lfs	   10,	144(%r4)
	fadds	   5,	5,	29
	lfs	   11,	148(%r4)
	fadds	   0,	0,	6
	lfs	   12,	160(%r4)
	fadds	   1,	1,	7
	lfs	   13,	164(%r4)
	fadds	   2,	2,	8
	lfs	   14,	168(%r4)
	fadds	   3,	3,	9
	lfs	   15,	172(%r4)
	fadds	   4,	4,	10
	lfs	   16,	176(%r4)
	fadds	   5,	5,	11
	lfs	   17,	180(%r4)
	fadds	   0,	0,	12
	addi %r6, %r6, 512 
	fadds	   1,	1,	13
	lfs	   18,	192(%r4)
	fadds	   2,	2,	14
	lfs	   19,	196(%r4)
	fadds	   3,	3,	15
	lfs	   20,	200(%r4)
	fadds	   4,	4,	16
	lfs	   21,	204(%r4)
	fadds	   5,	5,	17
	lfs	   22,	208(%r4)
	fadds	   0,	0,	18
	lfs	   23,	212(%r4)
	fadds	   1,	1,	19
	lfs	   24,	224(%r4)
	fadds	   2,	2,	20
	lfs	   25,	228(%r4)
	fadds	   3,	3,	21
	lfs	   26,	232(%r4)
	fadds	   4,	4,	22
	lfs	   27,	236(%r4)
	fadds	   5,	5,	23
	lfs	   28,	240(%r4)
	fadds	   0,	0,	24
	lfs	   29,	244(%r4)
	fadds	   1,	1,	25
	lfs	   6,	256(%r4)
	fadds	   2,	2,	26
	lfs	   7,	260(%r4)
	fadds	   3,	3,	27
	lfs	   8,	264(%r4)
	fadds	   4,	4,	28
	lfs	   9,	268(%r4)
	fadds	   5,	5,	29
	lfs	   10,	272(%r4)
	fadds	   0,	0,	6
	lfs	   11,	276(%r4)
	fadds	   1,	1,	7
	lfs	   12,	288(%r4)
	fadds	   2,	2,	8
	lfs	   13,	292(%r4)
	fadds	   3,	3,	9
	lfs	   14,	296(%r4)
	fadds	   4,	4,	10
	lfs	   15,	300(%r4)
	fadds	   5,	5,	11
	lfs	   16,	304(%r4)
	fadds	   0,	0,	12
	lfs	   17,	308(%r4)
	fadds	   1,	1,	13
	addi %r5, %r5, 24 
	fadds	   2,	2,	14
	lfs	   18,	320(%r6)
	fadds	   3,	3,	15
	lfs	   19,	324(%r6)
	fadds	   4,	4,	16
	lfs	   20,	328(%r6)
	fadds	   5,	5,	17
	lfs	   21,	332(%r6)
	fadds	   0,	0,	18
	lfs	   22,	336(%r6)
	fadds	   1,	1,	19
	lfs	   23,	340(%r6)
	fadds	   2,	2,	20
	lfs	   24,	352(%r6)
	fadds	   3,	3,	21
	lfs	   25,	356(%r6)
	fadds	   4,	4,	22
	lfs	   26,	360(%r6)
	fadds	   5,	5,	23
	lfs	   27,	364(%r6)
	fadds	   0,	0,	24
	lfs	   28,	368(%r6)
	fadds	   1,	1,	25
	lfs	   29,	372(%r6)
	fadds	   2,	2,	26
	lfs	   6,	384(%r6)
	fadds	   3,	3,	27
	lfs	   7,	388(%r6)
	fadds	   4,	4,	28
	lfs	   8,	392(%r6)
	fadds	   5,	5,	29
	lfs	   9,	396(%r6)
	fadds	   0,	0,	6
	lfs	   10,	400(%r6)
	fadds	   1,	1,	7
	lfs	   11,	404(%r6)
	fadds	   2,	2,	8
	lfs	   12,	416(%r6)
	fadds	   3,	3,	9
	lfs	   13,	420(%r6)
	fadds	   4,	4,	10
	lfs	   14,	424(%r6)
	fadds	   5,	5,	11
	lfs	   15,	428(%r6)
	fadds	   0,	0,	12
	lfs	   16,	432(%r6)
	fadds	   1,	1,	13
	lfs	   17,	436(%r6)
	fadds	   2,	2,	14
	lfs	   18,	448(%r6)
	fadds	   3,	3,	15
	lfs	   19,	452(%r6)
	fadds	   4,	4,	16
	lfs	   20,	456(%r6)
	fadds	   5,	5,	17
	lfs	   21,	460(%r6)
	fadds	   0,	0,	18
	lfs	   22,	464(%r6)
	fadds	   1,	1,	19
	lfs	   23,	468(%r6)
	fadds	   2,	2,	20
	lfs	   24,	480(%r6)
	fadds	   3,	3,	21
	lfs	   25,	484(%r6)
	fadds	   4,	4,	22
	lfs	   26,	488(%r6)
	fadds	   5,	5,	23
	lfs	   27,	492(%r6)
	fadds	   0,	0,	24
	lfs	   28,	496(%r6)
	fadds	   1,	1,	25
	lfs	   29,	500(%r6)
	fadds	   2,	2,	26

	fadds	   3,	3,	27

	fadds	   4,	4,	28

	fadds	   5,	5,	29









	stfs	   0,	0(%r5)
	stfs	   1,	4(%r5)
	stfs	   2,	8(%r5)
	stfs	   3,	12(%r5)
	stfs	   4,	16(%r5)
	stfs	   5,	20(%r5)
lab0:
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
