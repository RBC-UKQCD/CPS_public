	.section	".text"
	.align 2
	.globl pt_cmv2_dag_pad
	.type	pt_cmv2_dag_pad,@function
pt_cmv2_dag_pad:
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
	lfd	   20,	16(%r27)
	lfd	   21,	24(%r27)
	lfd	   22,	32(%r27)
	lfd	   23,	40(%r27)
	add	%r18, %r5, %r9
	lfd	   0,	0 ( % r19)
	lfd	   1,	8 ( % r19)
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
	lfd	   24,	48(%r27)
	add	%r25, %r25, %r7
	fmul	 8, 20,  0
	add	%r29, %r29, %r5

	fmul	 9, 20,  1
	stfd	   12,	0(%r24)
	fmul	10, 22,  0

	fmul	11, 22,  1
	lfd	   25,	56(%r27)
	fmadd	 6, 19,  1,  6
	stfd	   13,	8(%r24)
	fnmsub	 7, 19,  0,  7
	lfd	   26,	64(%r27)
	fmadd	 8, 21,  1,  8
	stfd	   14,	16(%r24)
	fnmsub	 9, 21,  0,  9
	lfd	   27,	72(%r27)
	fmadd	10, 23,  1, 10
	stfd	   15,	24(%r24)
	fnmsub	11, 23,  0, 11
	lfd	   28,	80(%r27)
	stfd	   16,	32(%r24)
	lfd	   29,	88(%r27)
	stfd	   17,	40(%r24)
	lfd	   2,	16 ( % r19)
	stfd	   30,	48(%r24)
	lfd	   3,	24 ( % r19)
	stfd	   30,	56(%r24)
#<<<<<< Break 1
	lfd	   18,	96(%r27)
	fmadd	 6, 24,  2,  6
	dcbt	 %r9,	 %r29
	fmadd	 7, 24,  3,  7
	lfd	   19,	104(%r27)
	fmadd	 8, 26,  2,  8

	fmadd	 9, 26,  3,  9
	lfd	   20,	112(%r27)
	fmadd	10, 28,  2, 10

	fmadd	11, 28,  3, 11
	lfd	   21,	120(%r27)
	fmadd	 6, 25,  3,  6
	dcbt	 %r10,	 %r29
	fnmsub	 7, 25,  2,  7
	lfd	   22,	128(%r27)
	fmadd	 8, 27,  3,  8

	fnmsub	 9, 27,  2,  9
	lfd	   23,	136(%r27)
	fmadd	10, 29,  3, 10

	fnmsub	11, 29,  2, 11
	lfd	   4,	32 ( % r19)
	dcbt	 %r11,	 %r29
	lfd	   5,	40 ( % r19)

#<<<<<< Break 2
	lfd	   24,	0(%r28)
	mr	 %r24,	 %r26
	fmadd	 6, 18,  4,  6

	fmadd	 7, 18,  5,  7
	lfd	   25,	8(%r28)
	fmadd	 8, 20,  4,  8
	dcbt	 %r9,	 %r21
	fmadd	 9, 20,  5,  9
	lfd	   26,	16(%r28)
	fmadd	10, 22,  4, 10
	dcbt	 %r10,	 %r21
	fmadd	11, 22,  5, 11
	lfd	   27,	24(%r28)
	fmadd	 6, 19,  5,  6

	fnmsub	 7, 19,  4,  7
	lfd	   28,	32(%r28)
	fmadd	 8, 21,  5,  8
	dcbt	 %r12,	 %r29
	fnmsub	 9, 21,  4,  9
	lfd	   29,	40(%r28)
	fmadd	10, 23,  5, 10

	fnmsub	11, 23,  4, 11
	lfd	   0,	0 ( % r20)

	lfd	   1,	8 ( % r20)
	dcbt	 %r13,	 %r29
#<<<<<< Break 3
	mr	 %r19,	 %r21
	mr	 %r27,	 %r29
	subic.	 3,	3,	1
#<<<<<< Break 4
	lwz	%r22, 36(%r4)
	fmul	12, 24,  0
	lwz	%r26, 40(%r4)
	lfd	   18,	48(%r28)
	lwz	%r30, 44(%r4)
	fmul	13, 24,  1
	add	%r22, %r22, %r6

	add	%r26, %r26, %r7
	fmul	14, 26,  0
	add	%r30, %r30, %r5
	stfd	   6,	0(%r23)
	fmul	15, 26,  1

	fmul	16, 28,  0
	lfd	   19,	56(%r28)
	fmul	17, 28,  1
	stfd	   7,	8(%r23)
	fmadd	12, 25,  1, 12
	lfd	   20,	64(%r28)
	fnmsub	13, 25,  0, 13
	stfd	   8,	16(%r23)
	fmadd	14, 27,  1, 14
	lfd	   21,	72(%r28)
	fnmsub	15, 27,  0, 15
	stfd	   9,	24(%r23)
	fmadd	16, 29,  1, 16
	lfd	   22,	80(%r28)
	fnmsub	17, 29,  0, 17
	stfd	   10,	32(%r23)
	lfd	   23,	88(%r28)
	stfd	   11,	40(%r23)
	lfd	   2,	16 ( % r20)
	stfd	   30,	48(%r23)
	lfd	   3,	24 ( % r20)
	stfd	   30,	56(%r23)
#<<<<<< Break 5
	lfd	   24,	96(%r28)
	fmadd	12, 18,  2, 12
	dcbt	 %r9,	 %r30
	fmadd	13, 18,  3, 13
	lfd	   25,	104(%r28)
	fmadd	14, 20,  2, 14

	fmadd	15, 20,  3, 15
	lfd	   26,	112(%r28)
	fmadd	16, 22,  2, 16

	fmadd	17, 22,  3, 17
	lfd	   27,	120(%r28)
	fmadd	12, 19,  3, 12
	dcbt	 %r10,	 %r30
	fnmsub	13, 19,  2, 13
	lfd	   28,	128(%r28)
	fmadd	14, 21,  3, 14

	fnmsub	15, 21,  2, 15
	lfd	   29,	136(%r28)
	fmadd	16, 23,  3, 16

	fnmsub	17, 23,  2, 17
	lfd	   4,	32 ( % r20)
	dcbt	 %r11,	 %r30
	lfd	   5,	40 ( % r20)

#<<<<<< Break 6
	lfd	   18,	0(%r29)
	mr	 %r23,	 %r25
	fmadd	12, 24,  4, 12
	addi %r4, %r4, 24 

	mr	 %r20,	 %r22
	fmadd	13, 24,  5, 13
	mr	 %r28,	 %r30
	lfd	   19,	8(%r29)

	fmadd	14, 26,  4, 14

	dcbt	 %r9,	 %r22

	fmadd	15, 26,  5, 15

	lfd	   20,	16(%r29)

	fmadd	16, 28,  4, 16

	dcbt	 %r10,	 %r22

	fmadd	17, 28,  5, 17
	lfd	   21,	24(%r29)
	fmadd	12, 25,  5, 12

	fnmsub	13, 25,  4, 13
	lfd	   22,	32(%r29)
	fmadd	14, 27,  5, 14
	dcbt	 %r12,	 %r30
	fnmsub	15, 27,  4, 15
	lfd	   23,	40(%r29)
	fmadd	16, 29,  5, 16

	fnmsub	17, 29,  4, 17
	lfd	   0,	0 ( % r21)
#<<<<<< Break 7

	lfd	   1,	8 ( % r21)
	dcbt	 %r13,	 %r30


	bf gt,	lab2
	b  lab1
lab2:
#<<<<<< Break 0
	fmul	 6, 18,  0
	lfd	   24,	48(%r27)
	fmul	 7, 18,  1
	stfd	   12,	0(%r24)
	fmul	 8, 20,  0
	lfd	   25,	56(%r27)
	fmul	 9, 20,  1
	stfd	   13,	8(%r24)
	fmul	10, 22,  0
	lfd	   26,	64(%r27)
	fmul	11, 22,  1
	stfd	   14,	16(%r24)
	fmadd	 6, 19,  1,  6
	lfd	   27,	72(%r27)
	fnmsub	 7, 19,  0,  7
	stfd	   15,	24(%r24)
	fmadd	 8, 21,  1,  8
	lfd	   28,	80(%r27)
	fnmsub	 9, 21,  0,  9
	stfd	   16,	32(%r24)
	fmadd	10, 23,  1, 10
	lfd	   29,	88(%r27)
	fnmsub	11, 23,  0, 11
	stfd	   17,	40(%r24)
	lfd	   2,	16 ( % r19)
	stfd	   30,	48(%r24)
	lfd	   3,	24 ( % r19)
	stfd	   30,	56(%r24)
#<<<<<< Break 1
	lfd	   18,	96(%r27)
	mr	 %r20,	 %r22
	fmadd	 6, 24,  2,  6
	mr	 %r28,	 %r30
	lfd	   19,	104(%r27)
	fmadd	 7, 24,  3,  7
	lfd	   20,	112(%r27)
	fmadd	 8, 26,  2,  8
	lfd	   21,	120(%r27)
	fmadd	 9, 26,  3,  9
	lfd	   22,	128(%r27)
	fmadd	10, 28,  2, 10
	lfd	   23,	136(%r27)
	fmadd	11, 28,  3, 11
	lfd	   4,	32 ( % r19)
	fmadd	 6, 25,  3,  6
	lfd	   5,	40 ( % r19)
	fnmsub	 7, 25,  2,  7
	fmadd	 8, 27,  3,  8
	fnmsub	 9, 27,  2,  9
	fmadd	10, 29,  3, 10
	fnmsub	11, 29,  2, 11
#<<<<<< Break 2
	fmadd	 6, 18,  4,  6
	mr	 %r24,	 %r26
	lfd	   24,	0(%r28)
	fmadd	 7, 18,  5,  7
	lfd	   25,	8(%r28)
	fmadd	 8, 20,  4,  8
	lfd	   26,	16(%r28)
	fmadd	 9, 20,  5,  9
	lfd	   27,	24(%r28)
	fmadd	10, 22,  4, 10
	lfd	   28,	32(%r28)
	fmadd	11, 22,  5, 11
	lfd	   29,	40(%r28)
	fmadd	 6, 19,  5,  6
	lfd	   0,	0 ( % r20)
	fnmsub	 7, 19,  4,  7
	lfd	   1,	8 ( % r20)
	fmadd	 8, 21,  5,  8
	fnmsub	 9, 21,  4,  9
	fmadd	10, 23,  5, 10
	fnmsub	11, 23,  4, 11
#<<<<<< Break 3
#<<<<<< Break 4
	lfd	   18,	48(%r28)
	fmul	12, 24,  0
	stfd	   6,	0(%r23)
	fmul	13, 24,  1
	lfd	   19,	56(%r28)
	fmul	14, 26,  0
	stfd	   7,	8(%r23)
	fmul	15, 26,  1
	lfd	   20,	64(%r28)
	fmul	16, 28,  0
	stfd	   8,	16(%r23)
	fmul	17, 28,  1
	lfd	   21,	72(%r28)
	fmadd	12, 25,  1, 12
	stfd	   9,	24(%r23)
	fnmsub	13, 25,  0, 13
	lfd	   22,	80(%r28)
	fmadd	14, 27,  1, 14
	stfd	   10,	32(%r23)
	fnmsub	15, 27,  0, 15
	lfd	   23,	88(%r28)
	fmadd	16, 29,  1, 16
	stfd	   11,	40(%r23)
	fnmsub	17, 29,  0, 17
	lfd	   2,	16 ( % r20)
	stfd	   30,	48(%r23)
	lfd	   3,	24 ( % r20)
	stfd	   30,	56(%r23)
#<<<<<< Break 5
	lfd	   24,	96(%r28)
	fmadd	12, 18,  2, 12
	lfd	   25,	104(%r28)
	fmadd	13, 18,  3, 13
	lfd	   26,	112(%r28)
	fmadd	14, 20,  2, 14
	lfd	   27,	120(%r28)
	fmadd	15, 20,  3, 15
	lfd	   28,	128(%r28)
	fmadd	16, 22,  2, 16
	lfd	   29,	136(%r28)
	fmadd	17, 22,  3, 17
	lfd	   4,	32 ( % r20)
	fmadd	12, 19,  3, 12
	lfd	   5,	40 ( % r20)
	fnmsub	13, 19,  2, 13
	fmadd	14, 21,  3, 14
	fnmsub	15, 21,  2, 15
	fmadd	16, 23,  3, 16
	fnmsub	17, 23,  2, 17
#<<<<<< Break 6
	fmadd	12, 24,  4, 12

	fmadd	13, 24,  5, 13

	fmadd	14, 26,  4, 14

	fmadd	15, 26,  5, 15

	fmadd	16, 28,  4, 16

	fmadd	17, 28,  5, 17

	fmadd	12, 25,  5, 12

	fnmsub	13, 25,  4, 13
	fmadd	14, 27,  5, 14
	fnmsub	15, 27,  4, 15
	fmadd	16, 29,  5, 16
	fnmsub	17, 29,  4, 17
#<<<<<< Break 7
	stfd	   12,	0(%r24)
	stfd	   13,	8(%r24)
	stfd	   14,	16(%r24)
	stfd	   15,	24(%r24)
	stfd	   16,	32(%r24)
	stfd	   17,	40(%r24)
	stfd	   30,	48(%r24)
	stfd	   30,	56(%r24)
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
