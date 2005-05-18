	.section	".text"
	.align 2
	.globl asq_cmv_4_s
	.type	asq_cmv_4_s,@function
asq_cmv_4_s:
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
	li   27,	8
	li   28,	40
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
	lwz	%r12, 0(%r5)
	lwz	%r14, 4(%r5)
	lwz	%r11, 80(%r5)
	lwz	%r13, 84(%r5)
	addi	%r3,	%r3,	1
	srawi.	3,	3,	1
	bf gt,	lab0
	lfs	   12,	8(%r5)
	lfs	   13,	12(%r5)
	lfs	   18,	32(%r5)
	lfs	   19,	36(%r5)
	lfs	   24,	56(%r5)
	lfs	   25,	60(%r5)
	lwz	%r30, 160(%r5)
	lwz	%r9, 164(%r5)
	lwz	%r31, 240(%r5)
	lwz	%r10, 244(%r5)
	lfs	   14,	16(%r5)
	lfs	   15,	20(%r5)
	lfs	   20,	40(%r5)
	lfs	   21,	44(%r5)
	lfs	   26,	64(%r5)
	lfs	   27,	68(%r5)
	add	%r12, %r12, %r6
	add	%r11, %r11, %r6
	add	%r14, %r14, %r7
	add	%r13, %r13, %r7
	add	%r29, %r5, %r16
	lfs	   0,	0 ( % r12)
	lfs	   1,	4 ( % r12)
	lfs	   2,	8 ( % r12)
	lfs	   3,	12 ( % r12)
	mr	 %r10,	 %r13
	dcbt	 %r16,	 %r5
	dcbt	 %r17,	 %r5
	dcbt	 %r18,	 %r5
	dcbt	 %r19,	 %r5
	dcbt	 %r20,	 %r5
	dcbt	 %r21,	 %r5





	addi %r29, %r29, 176 
	subic.	 3,	3,	1
	bf gt,	lab2
	b  lab1
lab1:
	fmuls	 6, 12,  0
	dcbt	 %r16,	 %r29
	fmuls	 7, 12,  1
	lwz	%r9, 164(%r5)
	fmuls	 8, 18,  0

	fmuls	 9, 18,  1
	lfs	   20,	40(%r5)
	fmuls	10, 24,  0

	fmuls	11, 24,  1
	lfs	   21,	44(%r5)
	fnmsubs	 6, 13,  1,  6

	fmadds	 7, 13,  0,  7
	lfs	   26,	64(%r5)
	fnmsubs	 8, 19,  1,  8

	add	%r9, %r9, %r7
	lfs	   27,	68(%r5)
	fmadds	 9, 19,  0,  9

	fnmsubs	10, 25,  1, 10
	lfs	   4,	16 ( % r12)
	fmadds	11, 25,  0, 11

	mr	 %r13,	 %r10
	lfs	   5,	20 ( % r12)
	fmadds	 6, 14,  2,  6
	dcbt	 %r17,	 %r29
	fmadds	 7, 14,  3,  7
	lfs	   16,	24(%r5)
	fmadds	 8, 20,  2,  8

	fmadds	 9, 20,  3,  9
	lfs	   17,	28(%r5)
	fmadds	10, 26,  2, 10

	fmadds	11, 26,  3, 11
	lfs	   22,	48(%r5)
	fnmsubs	 6, 15,  3,  6

	fmadds	 7, 15,  2,  7
	lfs	   23,	52(%r5)
	fnmsubs	 8, 21,  3,  8

	fmadds	 9, 21,  2,  9
	lfs	   28,	72(%r5)
	fnmsubs	10, 27,  3, 10

	fmadds	11, 27,  2, 11
	lfs	   29,	76(%r5)
	fmadds	 6, 16,  4,  6

	fmadds	 7, 16,  5,  7
	lfs	   12,	88(%r5)
	fmadds	 8, 22,  4,  8
	dcbt	 %r18,	 %r29
	fmadds	 9, 22,  5,  9
	lfs	   13,	92(%r5)
	fmadds	10, 28,  4, 10

	fmadds	11, 28,  5, 11
	lfs	   18,	112(%r5)
	fnmsubs	 6, 17,  5,  6

	fmadds	 7, 17,  4,  7
	lfs	   19,	116(%r5)
	fnmsubs	 8, 23,  5,  8

	fmadds	 9, 23,  4,  9
	lfs	   24,	136(%r5)
	fnmsubs	10, 29,  5, 10

	fmadds	11, 29,  4, 11
	stfs	   6,	0(%r14)
	subic.	 3,	3,	1

	lfs	   25,	140(%r5)
	stfs	   7,	4(%r14)
	lfs	   14,	96(%r5)
	stfs	   8,	8(%r14)
	lfs	   15,	100(%r5)
	stfs	   9,	12(%r14)
	stfs	   10,	16(%r14)
	stfs	   11,	20(%r14)
	fmuls	 6, 12,  0
	dcbt	 %r19,	 %r29
	fmuls	 7, 12,  1

	fmuls	 8, 18,  0

	fmuls	 9, 18,  1

	fmuls	10, 24,  0

	fmuls	11, 24,  1
	lwz	%r10, 244(%r5)
	fnmsubs	 6, 13,  1,  6

	fmadds	 7, 13,  0,  7
	lfs	   20,	120(%r5)
	fnmsubs	 8, 19,  1,  8

	add	%r10, %r10, %r7
	lfs	   21,	124(%r5)
	fmadds	 9, 19,  0,  9

	fnmsubs	10, 25,  1, 10
	lfs	   26,	144(%r5)
	fmadds	11, 25,  0, 11

	mr	 %r14,	 %r9
	lfs	   27,	148(%r5)
	fmadds	 6, 14,  2,  6

	fmadds	 7, 14,  3,  7
	lfs	   16,	104(%r5)
	fmadds	 8, 20,  2,  8

	fmadds	 9, 20,  3,  9
	lfs	   17,	108(%r5)
	fmadds	10, 26,  2, 10

	fmadds	11, 26,  3, 11
	lfs	   22,	128(%r5)
	fnmsubs	 6, 15,  3,  6
	dcbt	 %r20,	 %r29
	fmadds	 7, 15,  2,  7
	lfs	   23,	132(%r5)
	fnmsubs	 8, 21,  3,  8

	fmadds	 9, 21,  2,  9
	lfs	   28,	152(%r5)
	fnmsubs	10, 27,  3, 10

	fmadds	11, 27,  2, 11
	lfs	   29,	156(%r5)
	fmadds	 6, 16,  4,  6

	fmadds	 7, 16,  5,  7
	lfs	   12,	168(%r5)
	fmadds	 8, 22,  4,  8

	fmadds	 9, 22,  5,  9
	lfs	   13,	172(%r5)
	fmadds	10, 28,  4, 10

	fmadds	11, 28,  5, 11
	lfs	   18,	192(%r5)
	fnmsubs	 6, 17,  5,  6

	fmadds	 7, 17,  4,  7
	lfs	   19,	196(%r5)
	fnmsubs	 8, 23,  5,  8

	fmadds	 9, 23,  4,  9
	lfs	   24,	216(%r5)
	fnmsubs	10, 29,  5, 10

	fmadds	11, 29,  4, 11
	stfs	   6,	0(%r13)



	lfs	   25,	220(%r5)

	stfs	   7,	4(%r13)

	lfs	   14,	176(%r5)
	addi %r29, %r29, 160 
	stfs	   8,	8(%r13)

	lfs	   15,	180(%r5)

	stfs	   9,	12(%r13)

	stfs	   10,	16(%r13)

	stfs	   11,	20(%r13)

	addi %r5, %r5, 160 
	bf gt,	lab2
	fmuls	 6, 12,  0
	dcbt	 %r16,	 %r29
	fmuls	 7, 12,  1

	fmuls	 8, 18,  0

	fmuls	 9, 18,  1

	fmuls	10, 24,  0

	fmuls	11, 24,  1
	lwz	%r30, 160(%r5)
	fnmsubs	 6, 13,  1,  6

	fmadds	 7, 13,  0,  7
	lwz	%r9, 164(%r5)
	fnmsubs	 8, 19,  1,  8

	add	%r30, %r30, %r6
	lfs	   20,	40(%r5)
	add	%r9, %r9, %r7
	dcbt	 %r17,	 %r29
	fmadds	 9, 19,  0,  9
	lfs	   21,	44(%r5)
	fnmsubs	10, 25,  1, 10

	fmadds	11, 25,  0, 11
	lfs	   26,	64(%r5)
	mr	 %r13,	 %r10

	fmadds	 6, 14,  2,  6
	lfs	   27,	68(%r5)
	fmadds	 7, 14,  3,  7

	fmadds	 8, 20,  2,  8
	lfs	   16,	24(%r5)
	fmadds	 9, 20,  3,  9

	fmadds	10, 26,  2, 10
	lfs	   17,	28(%r5)
	fmadds	11, 26,  3, 11
	dcbt	 %r18,	 %r29
	fnmsubs	 6, 15,  3,  6
	lfs	   22,	48(%r5)
	fmadds	 7, 15,  2,  7

	fnmsubs	 8, 21,  3,  8
	lfs	   23,	52(%r5)
	fmadds	 9, 21,  2,  9

	fnmsubs	10, 27,  3, 10
	lfs	   28,	72(%r5)
	fmadds	11, 27,  2, 11

	fmadds	 6, 16,  4,  6
	lfs	   29,	76(%r5)
	fmadds	 7, 16,  5,  7

	fmadds	 8, 22,  4,  8
	lfs	   12,	88(%r5)
	fmadds	 9, 22,  5,  9

	fmadds	10, 28,  4, 10
	lfs	   13,	92(%r5)
	fmadds	11, 28,  5, 11
	dcbt	 %r16,	 %r30
	fnmsubs	 6, 17,  5,  6
	lfs	   18,	112(%r5)
	fmadds	 7, 17,  4,  7
	dcbt	 %r17,	 %r30
	fnmsubs	 8, 23,  5,  8
	lfs	   19,	116(%r5)
	fmadds	 9, 23,  4,  9
	dcbt	 %r28,	 %r30
	fnmsubs	10, 29,  5, 10
	lfs	   24,	136(%r5)
	fmadds	11, 29,  4, 11
	stfs	   6,	0(%r14)
	mr	 %r12,	 %r30
	lfs	   25,	140(%r5)
	subic.	 3,	3,	1
	stfs	   7,	4(%r14)
	lfs	   14,	96(%r5)
	stfs	   8,	8(%r14)
	lfs	   15,	100(%r5)
	stfs	   9,	12(%r14)
	stfs	   10,	16(%r14)
	stfs	   11,	20(%r14)
	fmuls	 6, 12,  0
	dcbt	 %r19,	 %r29
	fmuls	 7, 12,  1

	fmuls	 8, 18,  0

	fmuls	 9, 18,  1

	fmuls	10, 24,  0

	fmuls	11, 24,  1
	lwz	%r10, 244(%r5)
	fnmsubs	 6, 13,  1,  6

	fmadds	 7, 13,  0,  7
	lfs	   20,	120(%r5)
	fnmsubs	 8, 19,  1,  8

	add	%r10, %r10, %r7
	lfs	   21,	124(%r5)
	fmadds	 9, 19,  0,  9

	fnmsubs	10, 25,  1, 10
	lfs	   26,	144(%r5)
	fmadds	11, 25,  0, 11

	mr	 %r14,	 %r9
	lfs	   27,	148(%r5)
	fmadds	 6, 14,  2,  6
	dcbt	 %r20,	 %r29
	fmadds	 7, 14,  3,  7
	lfs	   16,	104(%r5)
	fmadds	 8, 20,  2,  8

	fmadds	 9, 20,  3,  9
	lfs	   17,	108(%r5)
	fmadds	10, 26,  2, 10

	fmadds	11, 26,  3, 11
	lfs	   22,	128(%r5)
	fnmsubs	 6, 15,  3,  6

	fmadds	 7, 15,  2,  7
	lfs	   23,	132(%r5)
	fnmsubs	 8, 21,  3,  8

	fmadds	 9, 21,  2,  9
	lfs	   28,	152(%r5)
	fnmsubs	10, 27,  3, 10

	fmadds	11, 27,  2, 11
	lfs	   29,	156(%r5)
	fmadds	 6, 16,  4,  6

	fmadds	 7, 16,  5,  7
	lfs	   0,	0 ( % r12)
	fmadds	 8, 22,  4,  8

	fmadds	 9, 22,  5,  9
	lfs	   1,	4 ( % r12)
	fmadds	10, 28,  4, 10
	lfs	   12,	168(%r5)
	fmadds	11, 28,  5, 11
	lfs	   13,	172(%r5)
	fnmsubs	 6, 17,  5,  6
	lfs	   18,	192(%r5)
	fmadds	 7, 17,  4,  7
	lfs	   19,	196(%r5)
	fnmsubs	 8, 23,  5,  8
	lfs	   24,	216(%r5)
	fmadds	 9, 23,  4,  9
	lfs	   25,	220(%r5)
	fnmsubs	10, 29,  5, 10
	lfs	   2,	8 ( % r12)
	fmadds	11, 29,  4, 11
	lfs	   3,	12 ( % r12)

	lfs	   14,	176(%r5)

	lfs	   15,	180(%r5)

	stfs	   6,	0(%r13)

	stfs	   7,	4(%r13)
	addi %r29, %r29, 160 
	stfs	   8,	8(%r13)

	stfs	   9,	12(%r13)

	stfs	   10,	16(%r13)

	stfs	   11,	20(%r13)

	stfs	   11,	24(%r13)

	stfs	   11,	28(%r13)
	addi %r5, %r5, 160 
	bf gt,	lab2
	b  lab1
lab2:
	fmuls	 6, 12,  0

	fmuls	 7, 12,  1

	fmuls	 8, 18,  0

	fmuls	 9, 18,  1

	fmuls	10, 24,  0

	fmuls	11, 24,  1
	lwz	%r30, 160(%r5)
	fnmsubs	 6, 13,  1,  6

	fmadds	 7, 13,  0,  7
	lwz	%r9, 164(%r5)
	fnmsubs	 8, 19,  1,  8

	add	%r30, %r30, %r6
	lfs	   20,	40(%r5)
	add	%r9, %r9, %r7

	fmadds	 9, 19,  0,  9
	lfs	   21,	44(%r5)
	fnmsubs	10, 25,  1, 10

	fmadds	11, 25,  0, 11
	lfs	   26,	64(%r5)
	mr	 %r13,	 %r10

	fmadds	 6, 14,  2,  6
	lfs	   27,	68(%r5)
	fmadds	 7, 14,  3,  7

	fmadds	 8, 20,  2,  8
	lfs	   16,	24(%r5)
	fmadds	 9, 20,  3,  9

	fmadds	10, 26,  2, 10
	lfs	   17,	28(%r5)
	fmadds	11, 26,  3, 11

	fnmsubs	 6, 15,  3,  6
	lfs	   22,	48(%r5)
	fmadds	 7, 15,  2,  7

	fnmsubs	 8, 21,  3,  8
	lfs	   23,	52(%r5)
	fmadds	 9, 21,  2,  9

	fnmsubs	10, 27,  3, 10
	lfs	   28,	72(%r5)
	fmadds	11, 27,  2, 11

	fmadds	 6, 16,  4,  6
	lfs	   29,	76(%r5)
	fmadds	 7, 16,  5,  7

	fmadds	 8, 22,  4,  8
	lfs	   12,	88(%r5)
	fmadds	 9, 22,  5,  9

	fmadds	10, 28,  4, 10
	lfs	   13,	92(%r5)
	fmadds	11, 28,  5, 11

	fnmsubs	 6, 17,  5,  6
	lfs	   18,	112(%r5)
	fmadds	 7, 17,  4,  7

	fnmsubs	 8, 23,  5,  8
	lfs	   19,	116(%r5)
	fmadds	 9, 23,  4,  9

	fnmsubs	10, 29,  5, 10
	lfs	   24,	136(%r5)
	fmadds	11, 29,  4, 11

	mr	 %r12,	 %r30
	stfs	   6,	0(%r14)
	subic.	 3,	3,	1

	lfs	   25,	140(%r5)

	stfs	   7,	4(%r14)

	lfs	   14,	96(%r5)

	stfs	   8,	8(%r14)

	lfs	   15,	100(%r5)

	stfs	   9,	12(%r14)

	stfs	   10,	16(%r14)

	stfs	   11,	20(%r14)


	fmuls	 6, 12,  0

	fmuls	 7, 12,  1

	fmuls	 8, 18,  0

	fmuls	 9, 18,  1

	fmuls	10, 24,  0

	fmuls	11, 24,  1
	lwz	%r10, 244(%r5)
	fnmsubs	 6, 13,  1,  6

	fmadds	 7, 13,  0,  7
	lfs	   20,	120(%r5)
	fnmsubs	 8, 19,  1,  8

	add	%r10, %r10, %r7
	lfs	   21,	124(%r5)
	fmadds	 9, 19,  0,  9

	fnmsubs	10, 25,  1, 10
	lfs	   26,	144(%r5)
	fmadds	11, 25,  0, 11

	mr	 %r14,	 %r9
	lfs	   27,	148(%r5)
	fmadds	 6, 14,  2,  6

	fmadds	 7, 14,  3,  7
	lfs	   16,	104(%r5)
	fmadds	 8, 20,  2,  8

	fmadds	 9, 20,  3,  9
	lfs	   17,	108(%r5)
	fmadds	10, 26,  2, 10

	fmadds	11, 26,  3, 11
	lfs	   22,	128(%r5)
	fnmsubs	 6, 15,  3,  6

	fmadds	 7, 15,  2,  7
	lfs	   23,	132(%r5)
	fnmsubs	 8, 21,  3,  8

	fmadds	 9, 21,  2,  9
	lfs	   28,	152(%r5)
	fnmsubs	10, 27,  3, 10

	fmadds	11, 27,  2, 11
	lfs	   29,	156(%r5)
	fmadds	 6, 16,  4,  6

	fmadds	 7, 16,  5,  7

	fmadds	 8, 22,  4,  8

	fmadds	 9, 22,  5,  9

	fmadds	10, 28,  4, 10

	fmadds	11, 28,  5, 11

	fnmsubs	 6, 17,  5,  6

	fmadds	 7, 17,  4,  7

	fnmsubs	 8, 23,  5,  8

	fmadds	 9, 23,  4,  9

	fnmsubs	10, 29,  5, 10

	fmadds	11, 29,  4, 11









	addi %r29, %r29, 160 











	addi %r5, %r5, 160 



	stfs	   6,	0(%r13)

	stfs	   7,	4(%r13)

	stfs	   8,	8(%r13)

	stfs	   9,	12(%r13)

	stfs	   10,	16(%r13)
	stfs	   11,	20(%r13)
	stfs	   11,	24(%r13)
	stfs	   11,	28(%r13)
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
