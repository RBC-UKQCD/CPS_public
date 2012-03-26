	.section	".text"
	.align 2
	.globl pt_cmv
	.type	pt_cmv,@function
pt_cmv:
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
	li   9,	0
	li   10,	32
	li   11,	64
	li   12,	96
	li   13,	128
	li   14,	160
	li   15,	192
	li   16,	224
	li   17,	256
	dcbt	 %r9,	 %r4
	dcbt	 %r10,	 %r4
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
	lwz	%r19, 0(%r4)
	lwz	%r23, 4(%r4)
	lwz	%r27, 8(%r4)
	dcbt	 %r11,	 %r4
	lwz	%r20, 12(%r4)
	lwz	%r24, 16(%r4)
	lwz	%r28, 20(%r4)
	add	%r27, %r27, %r5
	add	%r19, %r19, %r6
	add	%r23, %r23, %r7
	add	%r20, %r20, %r6
	add	%r24, %r24, %r7
	add	%r28, %r28, %r5
	dcbt	 %r12,	 %r4
	or.	 3,	3,	3
	mr	 %r26,	 %r24
	bf gt,	lab0
	lfd	   18,	0(%r27)
	lfd	   19,	8(%r27)
	lfd	   22,	48(%r27)
	lfd	   23,	56(%r27)
	lfd	   26,	96(%r27)
	lfd	   27,	104(%r27)
	add	%r18, %r5, %r9
	lfd	   0,	0 ( %r19)
	lfd	   1,	8 ( %r19)
	mr	 %r26,	 %r24
	subic.	 3,	3,	1
	bf gt,	lab2
	b  lab1
lab1:
#<<<<<< Break 0
	lwz	%r21, 24(%r4)
	fmul	 6, 18,  0
	lwz	%r25, 28(%r4)
	dcbt	 %r13,	 %r30
	lwz	%r29, 32(%r4)
	fmul	 7, 18,  1
	add	%r21, %r21, %r6
	lfd	   20,	16(%r27)
	add	%r25, %r25, %r7
	fmul	 8, 22,  0
	add	%r29, %r29, %r5

	fmul	 9, 22,  1
	stfd	   12,	0(%r24)
	fmul	10, 26,  0

	fmul	11, 26,  1
	lfd	   21,	24(%r27)
	fnmsub	 6, 19,  1,  6
	stfd	   13,	8(%r24)
	fmadd	 7, 19,  0,  7
	lfd	   24,	64(%r27)
	fnmsub	 8, 23,  1,  8
	stfd	   14,	16(%r24)
	fmadd	 9, 23,  0,  9
	lfd	   25,	72(%r27)
	fnmsub	10, 27,  1, 10
	stfd	   15,	24(%r24)
	fmadd	11, 27,  0, 11
	lfd	   28,	112(%r27)
	stfd	   16,	32(%r24)
	lfd	   29,	120(%r27)
	stfd	   17,	40(%r24)
	lfd	   2,	16 ( %r19)
	lfd	   3,	24 ( %r19)
#<<<<<< Break 1
	lfd	   18,	32(%r27)
	fmadd	 6, 20,  2,  6
	dcbt	 %r9,	 %r29
	fmadd	 7, 20,  3,  7
	lfd	   19,	40(%r27)
	fmadd	 8, 24,  2,  8

	fmadd	 9, 24,  3,  9
	lfd	   22,	80(%r27)
	fmadd	10, 28,  2, 10

	fmadd	11, 28,  3, 11
	lfd	   23,	88(%r27)
	fnmsub	 6, 21,  3,  6
	dcbt	 %r10,	 %r29
	fmadd	 7, 21,  2,  7
	lfd	   26,	128(%r27)
	fnmsub	 8, 25,  3,  8

	fmadd	 9, 25,  2,  9
	lfd	   27,	136(%r27)
	fnmsub	10, 29,  3, 10

	fmadd	11, 29,  2, 11
	lfd	   4,	32 ( %r19)
	dcbt	 %r11,	 %r29
	lfd	   5,	40 ( %r19)

#<<<<<< Break 2
	lfd	   20,	0(%r28)
	mr	 %r24,	 %r26
	fmadd	 6, 18,  4,  6

	fmadd	 7, 18,  5,  7
	lfd	   21,	8(%r28)
	fmadd	 8, 22,  4,  8
	dcbt	 %r9,	 %r21
	fmadd	 9, 22,  5,  9
	lfd	   24,	48(%r28)
	fmadd	10, 26,  4, 10
	dcbt	 %r10,	 %r21
	fmadd	11, 26,  5, 11
	lfd	   25,	56(%r28)
	fnmsub	 6, 19,  5,  6

	fmadd	 7, 19,  4,  7
	lfd	   28,	96(%r28)
	fnmsub	 8, 23,  5,  8
	dcbt	 %r12,	 %r29
	fmadd	 9, 23,  4,  9
	lfd	   29,	104(%r28)
	fnmsub	10, 27,  5, 10

	fmadd	11, 27,  4, 11
	lfd	   0,	0 ( %r20)

	lfd	   1,	8 ( %r20)
	dcbt	 %r13,	 %r29
#<<<<<< Break 3
	mr	 %r19,	 %r21
	mr	 %r27,	 %r29
	subic.	 3,	3,	1
#<<<<<< Break 4
	lwz	%r22, 36(%r4)
	fmul	12, 20,  0
	lwz	%r26, 40(%r4)
	lfd	   18,	16(%r28)
	lwz	%r30, 44(%r4)
	fmul	13, 20,  1
	add	%r22, %r22, %r6

	add	%r26, %r26, %r7
	fmul	14, 24,  0
	add	%r30, %r30, %r5
	stfd	   6,	0(%r23)
	fmul	15, 24,  1

	fmul	16, 28,  0
	lfd	   19,	24(%r28)
	fmul	17, 28,  1
	dcbt	 %r13,	 %r4
	fnmsub	12, 21,  1, 12
	stfd	   7,	8(%r23)
	fmadd	13, 21,  0, 13
	dcbt	 %r14,	 %r4
	fnmsub	14, 25,  1, 14
	lfd	   22,	64(%r28)
	fmadd	15, 25,  0, 15

	fnmsub	16, 29,  1, 16
	stfd	   8,	16(%r23)
	fmadd	17, 29,  0, 17
	lfd	   23,	72(%r28)
	stfd	   9,	24(%r23)
	lfd	   26,	112(%r28)
	stfd	   10,	32(%r23)
	lfd	   27,	120(%r28)
	stfd	   11,	40(%r23)
	lfd	   2,	16 ( %r20)
	lfd	   3,	24 ( %r20)
#<<<<<< Break 5
	lfd	   20,	32(%r28)
	fmadd	12, 18,  2, 12
	dcbt	 %r9,	 %r30
	fmadd	13, 18,  3, 13
	lfd	   21,	40(%r28)
	fmadd	14, 22,  2, 14

	fmadd	15, 22,  3, 15
	lfd	   24,	80(%r28)
	fmadd	16, 26,  2, 16

	fmadd	17, 26,  3, 17
	lfd	   25,	88(%r28)
	fnmsub	12, 19,  3, 12
	dcbt	 %r10,	 %r30
	fmadd	13, 19,  2, 13
	lfd	   28,	128(%r28)
	fnmsub	14, 23,  3, 14

	fmadd	15, 23,  2, 15
	lfd	   29,	136(%r28)
	fnmsub	16, 27,  3, 16

	fmadd	17, 27,  2, 17
	lfd	   4,	32 ( %r20)
	dcbt	 %r11,	 %r30
	lfd	   5,	40 ( %r20)

#<<<<<< Break 6
	lfd	   18,	0(%r29)
	mr	 %r23,	 %r25
	fmadd	12, 20,  4, 12
	addi %r4, %r4, 24 

	mr	 %r20,	 %r22
	fmadd	13, 20,  5, 13
	mr	 %r28,	 %r30
	lfd	   19,	8(%r29)

	fmadd	14, 24,  4, 14

	dcbt	 %r9,	 %r22

	fmadd	15, 24,  5, 15

	lfd	   22,	48(%r29)

	fmadd	16, 28,  4, 16

	dcbt	 %r10,	 %r22

	fmadd	17, 28,  5, 17
	lfd	   23,	56(%r29)
	fnmsub	12, 21,  5, 12

	fmadd	13, 21,  4, 13
	lfd	   26,	96(%r29)
	fnmsub	14, 25,  5, 14
	dcbt	 %r12,	 %r30
	fmadd	15, 25,  4, 15
	lfd	   27,	104(%r29)
	fnmsub	16, 29,  5, 16

	fmadd	17, 29,  4, 17
	lfd	   0,	0 ( %r21)
#<<<<<< Break 7

	lfd	   1,	8 ( %r21)
	dcbt	 %r13,	 %r30


	bf gt,	lab2
	b  lab1
lab2:
#<<<<<< Break 0
	fmul	 6, 18,  0
	lfd	   20,	16(%r27)
	fmul	 7, 18,  1
	stfd	   12,	0(%r24)
	fmul	 8, 22,  0
	lfd	   21,	24(%r27)
	fmul	 9, 22,  1
	stfd	   13,	8(%r24)
	fmul	10, 26,  0
	lfd	   24,	64(%r27)
	fmul	11, 26,  1
	stfd	   14,	16(%r24)
	fnmsub	 6, 19,  1,  6
	lfd	   25,	72(%r27)
	fmadd	 7, 19,  0,  7
	stfd	   15,	24(%r24)
	fnmsub	 8, 23,  1,  8
	lfd	   28,	112(%r27)
	fmadd	 9, 23,  0,  9
	stfd	   16,	32(%r24)
	fnmsub	10, 27,  1, 10
	lfd	   29,	120(%r27)
	fmadd	11, 27,  0, 11
	stfd	   17,	40(%r24)
	lfd	   2,	16 ( %r19)
	lfd	   3,	24 ( %r19)
#<<<<<< Break 1
	lfd	   18,	32(%r27)
	mr	 %r20,	 %r22
	fmadd	 6, 20,  2,  6
	mr	 %r28,	 %r30
	lfd	   19,	40(%r27)
	fmadd	 7, 20,  3,  7
	lfd	   22,	80(%r27)
	fmadd	 8, 24,  2,  8
	lfd	   23,	88(%r27)
	fmadd	 9, 24,  3,  9
	lfd	   26,	128(%r27)
	fmadd	10, 28,  2, 10
	lfd	   27,	136(%r27)
	fmadd	11, 28,  3, 11
	lfd	   4,	32 ( %r19)
	fnmsub	 6, 21,  3,  6
	lfd	   5,	40 ( %r19)
	fmadd	 7, 21,  2,  7
	fnmsub	 8, 25,  3,  8
	fmadd	 9, 25,  2,  9
	fnmsub	10, 29,  3, 10
	fmadd	11, 29,  2, 11
#<<<<<< Break 2
	fmadd	 6, 18,  4,  6
	mr	 %r24,	 %r26
	lfd	   20,	0(%r28)
	fmadd	 7, 18,  5,  7
	lfd	   21,	8(%r28)
	fmadd	 8, 22,  4,  8
	lfd	   24,	48(%r28)
	fmadd	 9, 22,  5,  9
	lfd	   25,	56(%r28)
	fmadd	10, 26,  4, 10
	lfd	   28,	96(%r28)
	fmadd	11, 26,  5, 11
	lfd	   29,	104(%r28)
	fnmsub	 6, 19,  5,  6
	lfd	   0,	0 ( %r20)
	fmadd	 7, 19,  4,  7
	lfd	   1,	8 ( %r20)
	fnmsub	 8, 23,  5,  8
	fmadd	 9, 23,  4,  9
	fnmsub	10, 27,  5, 10
	fmadd	11, 27,  4, 11
#<<<<<< Break 3
#<<<<<< Break 4
	lfd	   18,	16(%r28)
	fmul	12, 20,  0
	stfd	   6,	0(%r23)
	fmul	13, 20,  1
	lfd	   19,	24(%r28)
	fmul	14, 24,  0
	stfd	   7,	8(%r23)
	fmul	15, 24,  1
	lfd	   22,	64(%r28)
	fmul	16, 28,  0
	stfd	   8,	16(%r23)
	fmul	17, 28,  1
	lfd	   23,	72(%r28)
	fnmsub	12, 21,  1, 12
	stfd	   9,	24(%r23)
	fmadd	13, 21,  0, 13
	lfd	   26,	112(%r28)
	fnmsub	14, 25,  1, 14
	stfd	   10,	32(%r23)
	fmadd	15, 25,  0, 15
	lfd	   27,	120(%r28)
	fnmsub	16, 29,  1, 16
	stfd	   11,	40(%r23)
	fmadd	17, 29,  0, 17
	lfd	   2,	16 ( %r20)
	lfd	   3,	24 ( %r20)
#<<<<<< Break 5
	lfd	   20,	32(%r28)
	fmadd	12, 18,  2, 12
	lfd	   21,	40(%r28)
	fmadd	13, 18,  3, 13
	lfd	   24,	80(%r28)
	fmadd	14, 22,  2, 14
	lfd	   25,	88(%r28)
	fmadd	15, 22,  3, 15
	lfd	   28,	128(%r28)
	fmadd	16, 26,  2, 16
	lfd	   29,	136(%r28)
	fmadd	17, 26,  3, 17
	lfd	   4,	32 ( %r20)
	fnmsub	12, 19,  3, 12
	lfd	   5,	40 ( %r20)
	fmadd	13, 19,  2, 13
	fnmsub	14, 23,  3, 14
	fmadd	15, 23,  2, 15
	fnmsub	16, 27,  3, 16
	fmadd	17, 27,  2, 17
#<<<<<< Break 6
	fmadd	12, 20,  4, 12

	fmadd	13, 20,  5, 13

	fmadd	14, 24,  4, 14

	fmadd	15, 24,  5, 15

	fmadd	16, 28,  4, 16

	fmadd	17, 28,  5, 17

	fnmsub	12, 21,  5, 12

	fmadd	13, 21,  4, 13
	fnmsub	14, 25,  5, 14
	fmadd	15, 25,  4, 15
	fnmsub	16, 29,  5, 16
	fmadd	17, 29,  4, 17
#<<<<<< Break 7
	stfd	   12,	0(%r24)
	stfd	   13,	8(%r24)
	stfd	   14,	16(%r24)
	stfd	   15,	24(%r24)
	stfd	   16,	32(%r24)
	stfd	   17,	40(%r24)
lab0:
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
