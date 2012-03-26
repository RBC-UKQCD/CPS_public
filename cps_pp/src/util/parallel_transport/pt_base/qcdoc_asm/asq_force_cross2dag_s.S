	.section	".text"
	.align 2
	.globl asq_force_cross2dag_s
	.type	asq_force_cross2dag_s,@function
asq_force_cross2dag_s:
	la   %r1,	-864(%r1)
	std	%r14,	368(%r1)
	std	%r15,	376(%r1)
	std	%r16,	384(%r1)
	std	%r17,	392(%r1)
	std	%r18,	400(%r1)
	std	%r19,	408(%r1)
	std	%r20,	416(%r1)
	std	%r21,	424(%r1)
	std	%r22,	432(%r1)
	std	%r23,	440(%r1)
	std	%r24,	448(%r1)
	std	%r25,	456(%r1)
	std	%r26,	464(%r1)
	std	%r27,	472(%r1)
	std	%r28,	480(%r1)
	std	%r29,	488(%r1)
	std	%r30,	496(%r1)
	std	%r31,	504(%r1)
	li   10,	48
	li   11,	80
	li   12,	112
	li   13,	128
	li   14,	160
	li   15,	192
	li   16,	224
	li   17,	256
	li   18,	288
	li   19,	320
	li   20,	352
	li   21,	384
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
	addi %r5, %r5, -144 
	lfs	   0,	0(%r7)
	lfs	   1,	0(%r3)
	lfs	   2,	4(%r3)
	lfs	   3,	8(%r3)
	lfs	   4,	12(%r3)
	lfs	   5,	16(%r3)
	lfs	   6,	20(%r3)
	or.	 6,	6,	6
	subic.	 6,	6,	1
	lfs	   7,	0(%r4)
	lfs	   8,	4(%r4)
	lfs	   9,	8(%r4)
	lfs	   10,	12(%r4)
	lfs	   11,	16(%r4)
	lfs	   12,	20(%r4)
	bf gt,	lab0
	b  lab1
lab1:
	lfs	   11,	16(%r4)
	addi %r5, %r5, 144 
	lfs	   12,	20(%r4)
	fmuls	1,	0,	1
	lfs	   13,	0(%r5)
	fmuls	2,	0,	2
	lfs	   14,	4(%r5)
	fmuls	3,	0,	3
	lfs	   15,	8(%r5)
	fmuls	4,	0,	4
	lfs	   16,	12(%r5)
	fmuls	5,	0,	5
	lfs	   17,	16(%r5)
	fmuls	6,	0,	6
	lfs	   18,	20(%r5)
#<<<<<< Break 0
	lfs	   19,	24(%r5)
	fmadd	13,	1,	7,	13
	lfs	   20,	28(%r5)
	fmadd	14,	2,	7,	14
	lfs	   21,	32(%r5)
	fmadd	15,	1,	9,	15
	dcbt	 %r10,	 %r3
	fmadd	16,	2,	9,	16
	lfs	   22,	36(%r5)
	fmadd	17,	1,	11,	17
	lfs	   23,	40(%r5)
	fmadd	18,	2,	11,	18
	lfs	   24,	44(%r5)
	fmadd	13,	2,	8,	13
	dcbt	 %r11,	 %r3
	fnmsub	14,	1,	8,	14
	fmadd	15,	2,	10,	15
	fnmsub	16,	1,	10,	16
	fmadd	17,	2,	12,	17
	fnmsub	18,	1,	12,	18
#<<<<<< Break 1
	fmadd	19,	3,	7,	19
	lfs	   25,	48(%r5)
	fmadd	20,	4,	7,	20
	stfs	   13,	0(%r5)
	fmadd	21,	3,	9,	21
	lfs	   26,	52(%r5)
	fmadd	22,	4,	9,	22
	stfs	   14,	4(%r5)
	fmadd	23,	3,	11,	23
	lfs	   27,	56(%r5)
	fmadd	24,	4,	11,	24
	stfs	   15,	8(%r5)
	fmadd	19,	4,	8,	19
	dcbt	 %r12,	 %r3
	fnmsub	20,	3,	8,	20
	stfs	   16,	12(%r5)
	fmadd	21,	4,	10,	21
	lfs	   28,	60(%r5)
	fnmsub	22,	3,	10,	22
	stfs	   17,	16(%r5)
	fmadd	23,	4,	12,	23
	lfs	   29,	64(%r5)
	fnmsub	24,	3,	12,	24
	stfs	   18,	20(%r5)
	lfs	   30,	68(%r5)
	dcbt	 %r10,	 %r4
	lfs	   1,	24(%r3)
	lfs	   2,	28(%r3)
#<<<<<< Break 2
	lfs	   3,	32(%r3)
	fmadd	25,	5,	7,	25
	stfs	   19,	24(%r5)
	fmadd	26,	6,	7,	26
	dcbt	 %r11,	 %r4
	fmadd	27,	5,	9,	27
	stfs	   20,	28(%r5)
	fmadd	28,	6,	9,	28
	lfs	   4,	36(%r3)
	fmadd	29,	5,	11,	29
	stfs	   21,	32(%r5)
	fmadd	30,	6,	11,	30
	lfs	   13,	72(%r5)
	fmadd	25,	6,	8,	25
	stfs	   22,	36(%r5)
	fnmsub	26,	5,	8,	26
	lfs	   14,	76(%r5)
	fmadd	27,	6,	10,	27
	stfs	   23,	40(%r5)
	fnmsub	28,	5,	10,	28
	dcbt	 %r12,	 %r4
	fmadd	29,	6,	12,	29
	stfs	   24,	44(%r5)
	fnmsub	30,	5,	12,	30
	lfs	   15,	80(%r5)
	lfs	   16,	84(%r5)
	lfs	   17,	88(%r5)
	dcbt	 %r13,	 %r5
	lfs	   18,	92(%r5)
#<<<<<< Break 3
	lfs	   5,	40(%r3)
	fmuls	1,	0,	1
	stfs	   25,	48(%r5)
	fmuls	2,	0,	2
	lfs	   6,	44(%r3)
	fmuls	3,	0,	3
	stfs	   26,	52(%r5)
	fmuls	4,	0,	4
	dcbt	 %r14,	 %r5
	fmuls	5,	0,	5
	stfs	   27,	56(%r5)
	fmuls	6,	0,	6
	lfs	   7,	24(%r4)
	stfs	   28,	60(%r5)
	lfs	   8,	28(%r4)
	stfs	   29,	64(%r5)
	lfs	   9,	32(%r4)
	stfs	   30,	68(%r5)
	dcbt	 %r15,	 %r5
	lfs	   10,	36(%r4)
	lfs	   11,	40(%r4)
	lfs	   12,	44(%r4)
	dcbt	 %r16,	 %r5
#<<<<<< Break 7
	lfs	   19,	96(%r5)
	fmadd	13,	1,	7,	13
	lfs	   20,	100(%r5)
	fmadd	14,	2,	7,	14
	lfs	   21,	104(%r5)
	fmadd	15,	1,	9,	15
	dcbt	 %r17,	 %r5
	fmadd	16,	2,	9,	16
	lfs	   22,	108(%r5)
	fmadd	17,	1,	11,	17
	lfs	   23,	112(%r5)
	fmadd	18,	2,	11,	18
	lfs	   24,	116(%r5)
	fmadd	13,	2,	8,	13
	dcbt	 %r18,	 %r5
	fnmsub	14,	1,	8,	14
	fmadd	15,	2,	10,	15
	fnmsub	16,	1,	10,	16
	fmadd	17,	2,	12,	17
	fnmsub	18,	1,	12,	18
#<<<<<< Break 4
	fmadd	19,	3,	7,	19
	lfs	   25,	120(%r5)
	fmadd	20,	4,	7,	20
	stfs	   13,	72(%r5)
	fmadd	21,	3,	9,	21
	lfs	   26,	124(%r5)
	fmadd	22,	4,	9,	22
	stfs	   14,	76(%r5)
	fmadd	23,	3,	11,	23
	lfs	   27,	128(%r5)
	fmadd	24,	4,	11,	24
	stfs	   15,	80(%r5)
	fmadd	19,	4,	8,	19
	lfs	   28,	132(%r5)
	fnmsub	20,	3,	8,	20
	stfs	   16,	84(%r5)
	fmadd	21,	4,	10,	21
	lfs	   29,	136(%r5)
	fnmsub	22,	3,	10,	22
	stfs	   17,	88(%r5)
	fmadd	23,	4,	12,	23
	lfs	   30,	140(%r5)
	fnmsub	24,	3,	12,	24
	stfs	   18,	92(%r5)
#<<<<<< Break 5
	stfs	   19,	96(%r5)
	fmadd	25,	5,	7,	25
	lfs	   1,	48(%r3)
	fmadd	26,	6,	7,	26
	stfs	   20,	100(%r5)
	fmadd	27,	5,	9,	27
	lfs	   2,	52(%r3)
	fmadd	28,	6,	9,	28
	stfs	   21,	104(%r5)
	fmadd	29,	5,	11,	29
	lfs	   3,	56(%r3)
	fmadd	30,	6,	11,	30
	stfs	   22,	108(%r5)
	fmadd	25,	6,	8,	25
	lfs	   4,	60(%r3)
	fnmsub	26,	5,	8,	26
	stfs	   23,	112(%r5)
	fmadd	27,	6,	10,	27
	stfs	   24,	116(%r5)
	fnmsub	28,	5,	10,	28
	fmadd	29,	6,	12,	29
	fnmsub	30,	5,	12,	30
#<<<<<< Break 6
	addi %r4, %r4, 48 
	stfs	   25,	120(%r5)
	lfs	   5,	64(%r3)
	stfs	   26,	124(%r5)
	lfs	   6,	68(%r3)
	stfs	   27,	128(%r5)
#<<<<<< Break 8
	lfs	   7,	0(%r4)
	addi %r3, %r3, 48 
	stfs	   28,	132(%r5)
	subic.	 6,	6,	1
	lfs	   8,	4(%r4)
	stfs	   29,	136(%r5)
	lfs	   9,	8(%r4)
	stfs	   30,	140(%r5)
	lfs	   10,	12(%r4)
	bf gt,	lab0
	b  lab1
lab0:
	lfs	   11,	16(%r4)
	addi %r5, %r5, 144 
	lfs	   12,	20(%r4)
	fmuls	1,	0,	1
	lfs	   13,	0(%r5)
	fmuls	2,	0,	2
	lfs	   14,	4(%r5)
	fmuls	3,	0,	3
	lfs	   15,	8(%r5)
	fmuls	4,	0,	4
	lfs	   16,	12(%r5)
	fmuls	5,	0,	5
	lfs	   17,	16(%r5)
	fmuls	6,	0,	6
	lfs	   18,	20(%r5)
#<<<<<< Break 0
	lfs	   19,	24(%r5)
	fmadd	13,	1,	7,	13
	lfs	   20,	28(%r5)
	fmadd	14,	2,	7,	14
	lfs	   21,	32(%r5)
	fmadd	15,	1,	9,	15
	lfs	   22,	36(%r5)
	fmadd	16,	2,	9,	16
	lfs	   23,	40(%r5)
	fmadd	17,	1,	11,	17
	lfs	   24,	44(%r5)
	fmadd	18,	2,	11,	18
	fmadd	13,	2,	8,	13
	fnmsub	14,	1,	8,	14
	fmadd	15,	2,	10,	15
	fnmsub	16,	1,	10,	16
	fmadd	17,	2,	12,	17
	fnmsub	18,	1,	12,	18
#<<<<<< Break 1
	fmadd	19,	3,	7,	19
	lfs	   25,	48(%r5)
	fmadd	20,	4,	7,	20
	stfs	   13,	0(%r5)
	fmadd	21,	3,	9,	21
	lfs	   26,	52(%r5)
	fmadd	22,	4,	9,	22
	stfs	   14,	4(%r5)
	fmadd	23,	3,	11,	23
	lfs	   27,	56(%r5)
	fmadd	24,	4,	11,	24
	stfs	   15,	8(%r5)
	fmadd	19,	4,	8,	19
	lfs	   28,	60(%r5)
	fnmsub	20,	3,	8,	20
	stfs	   16,	12(%r5)
	fmadd	21,	4,	10,	21
	lfs	   29,	64(%r5)
	fnmsub	22,	3,	10,	22
	stfs	   17,	16(%r5)
	fmadd	23,	4,	12,	23
	lfs	   30,	68(%r5)
	fnmsub	24,	3,	12,	24
	stfs	   18,	20(%r5)
	lfs	   1,	24(%r3)
	lfs	   2,	28(%r3)
#<<<<<< Break 2
	lfs	   3,	32(%r3)
	fmadd	25,	5,	7,	25
	stfs	   19,	24(%r5)
	fmadd	26,	6,	7,	26
	lfs	   4,	36(%r3)
	fmadd	27,	5,	9,	27
	stfs	   20,	28(%r5)
	fmadd	28,	6,	9,	28
	lfs	   13,	72(%r5)
	fmadd	29,	5,	11,	29
	stfs	   21,	32(%r5)
	fmadd	30,	6,	11,	30
	lfs	   14,	76(%r5)
	fmadd	25,	6,	8,	25
	stfs	   22,	36(%r5)
	fnmsub	26,	5,	8,	26
	lfs	   15,	80(%r5)
	fmadd	27,	6,	10,	27
	stfs	   23,	40(%r5)
	fnmsub	28,	5,	10,	28
	lfs	   16,	84(%r5)
	fmadd	29,	6,	12,	29
	stfs	   24,	44(%r5)
	fnmsub	30,	5,	12,	30
	lfs	   17,	88(%r5)
	lfs	   18,	92(%r5)
#<<<<<< Break 3
	lfs	   5,	40(%r3)
	fmuls	1,	0,	1
	stfs	   25,	48(%r5)
	fmuls	2,	0,	2
	lfs	   6,	44(%r3)
	fmuls	3,	0,	3
	stfs	   26,	52(%r5)
	fmuls	4,	0,	4
	lfs	   7,	24(%r4)
	fmuls	5,	0,	5
	stfs	   27,	56(%r5)
	fmuls	6,	0,	6
	lfs	   8,	28(%r4)
	stfs	   28,	60(%r5)
	lfs	   9,	32(%r4)
	stfs	   29,	64(%r5)
	lfs	   10,	36(%r4)
	stfs	   30,	68(%r5)
	lfs	   11,	40(%r4)
	lfs	   12,	44(%r4)
#<<<<<< Break 7
	lfs	   19,	96(%r5)
	fmadd	13,	1,	7,	13
	lfs	   20,	100(%r5)
	fmadd	14,	2,	7,	14
	lfs	   21,	104(%r5)
	fmadd	15,	1,	9,	15
	lfs	   22,	108(%r5)
	fmadd	16,	2,	9,	16
	lfs	   23,	112(%r5)
	fmadd	17,	1,	11,	17
	lfs	   24,	116(%r5)
	fmadd	18,	2,	11,	18
	fmadd	13,	2,	8,	13
	fnmsub	14,	1,	8,	14
	fmadd	15,	2,	10,	15
	fnmsub	16,	1,	10,	16
	fmadd	17,	2,	12,	17
	fnmsub	18,	1,	12,	18
#<<<<<< Break 4
	fmadd	19,	3,	7,	19
	lfs	   25,	120(%r5)
	fmadd	20,	4,	7,	20
	stfs	   13,	72(%r5)
	fmadd	21,	3,	9,	21
	lfs	   26,	124(%r5)
	fmadd	22,	4,	9,	22
	stfs	   14,	76(%r5)
	fmadd	23,	3,	11,	23
	lfs	   27,	128(%r5)
	fmadd	24,	4,	11,	24
	stfs	   15,	80(%r5)
	fmadd	19,	4,	8,	19
	lfs	   28,	132(%r5)
	fnmsub	20,	3,	8,	20
	stfs	   16,	84(%r5)
	fmadd	21,	4,	10,	21
	lfs	   29,	136(%r5)
	fnmsub	22,	3,	10,	22
	stfs	   17,	88(%r5)
	fmadd	23,	4,	12,	23
	lfs	   30,	140(%r5)
	fnmsub	24,	3,	12,	24
	stfs	   18,	92(%r5)
#<<<<<< Break 5
	stfs	   19,	96(%r5)
	fmadd	25,	5,	7,	25
	stfs	   20,	100(%r5)
	fmadd	26,	6,	7,	26
	stfs	   21,	104(%r5)
	fmadd	27,	5,	9,	27
	stfs	   22,	108(%r5)
	fmadd	28,	6,	9,	28
	stfs	   23,	112(%r5)
	fmadd	29,	5,	11,	29
	stfs	   24,	116(%r5)
	fmadd	30,	6,	11,	30
	fmadd	25,	6,	8,	25
	fnmsub	26,	5,	8,	26
	fmadd	27,	6,	10,	27
	fnmsub	28,	5,	10,	28
	fmadd	29,	6,	12,	29
	fnmsub	30,	5,	12,	30
#<<<<<< Break 6
	addi %r4, %r4, 48 
	stfs	   25,	120(%r5)
#<<<<<< Break 8
	stfs	   26,	124(%r5)
	addi %r3, %r3, 48 
	stfs	   27,	128(%r5)
	subic.	 6,	6,	1
	stfs	   28,	132(%r5)
	stfs	   29,	136(%r5)
	stfs	   30,	140(%r5)
	ld	%r14,	368(%r1)
	ld	%r15,	376(%r1)
	ld	%r16,	384(%r1)
	ld	%r17,	392(%r1)
	ld	%r18,	400(%r1)
	ld	%r19,	408(%r1)
	ld	%r20,	416(%r1)
	ld	%r21,	424(%r1)
	ld	%r22,	432(%r1)
	ld	%r23,	440(%r1)
	ld	%r24,	448(%r1)
	ld	%r25,	456(%r1)
	ld	%r26,	464(%r1)
	ld	%r27,	472(%r1)
	ld	%r28,	480(%r1)
	ld	%r29,	488(%r1)
	ld	%r30,	496(%r1)
	ld	%r31,	504(%r1)
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
