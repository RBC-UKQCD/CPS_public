	.section	".text"
	.align 2
	.globl dirac_cmv_jcw_agg
	.type	dirac_cmv_jcw_agg,@function
dirac_cmv_jcw_agg:
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
	stfd	14,	112(%r1)
	li   17,	32
	stfd	15,	120(%r1)
	li   18,	64
	stfd	16,	128(%r1)
	li   19,	96
	stfd	17,	136(%r1)
	li   20,	128
	li   21,	160
	stfd	18,	144(%r1)
	li   22,	192
	li   23,	224
	stfd	19,	152(%r1)
	li   24,	256
	li   25,	288
	stfd	20,	160(%r1)
	li   26,	320
	li   27,	16
	stfd	21,	168(%r1)
	li   28,	48
	stfd	22,	176(%r1)
	stfd	23,	184(%r1)
	lwz	%r12, 0(%r5)
	stfd	24,	192(%r1)
	lwz	%r14, 4(%r5)
	stfd	25,	200(%r1)
	lwz	%r11, 152(%r5)
	stfd	26,	208(%r1)
	stfd	27,	216(%r1)
	lwz	%r13, 156(%r5)
	stfd	28,	224(%r1)
	stfd	29,	232(%r1)
	stfd	30,	240(%r1)
	stfd	31,	248(%r1)
	addi	%r3,	%r3,	1
	srawi.	3,	3,	1
	bf gt,	lab0
	lfd	   12,	8(%r5)
	lfd	   13,	16(%r5)
	lfd	   14,	24(%r5)
	lfd	   15,	32(%r5)
	lfd	   18,	56(%r5)
	lfd	   19,	64(%r5)
	lwz	%r30, 304(%r5)
	lwz	%r9, 308(%r5)
	lwz	%r31, 456(%r5)
	lwz	%r10, 460(%r5)
	lfd	   20,	72(%r5)
	lfd	   21,	80(%r5)
	lfd	   24,	104(%r5)
	lfd	   25,	112(%r5)
	lfd	   26,	120(%r5)
	lfd	   27,	128(%r5)
	add	%r12, %r12, %r6
	add	%r11, %r11, %r6
	add	%r14, %r14, %r7
	add	%r13, %r13, %r7
	add	%r29, %r5, %r16
	lfd	   0,	0 ( % r12)
	lfd	   1,	8 ( % r12)
	lfd	   2,	16 ( % r12)
	lfd	   3,	24 ( % r12)
	mr	 %r10,	 %r13
	addi %r29, %r29, 432
	subic.	 3,	3,	1
	bf gt,	lab2
	b  dirac_cmv_jcw_agg_lab1
dirac_cmv_jcw_agg_lab1:
	fmul	 6, 12,  0
	dcbt	 %r16,	 %r29
	fmul	 7, 12,  1
	lwz	%r30, 304(%r5)
	fmul	 8, 18,  0

	fmul	 9, 18,  1
	lwz	%r9, 308(%r5)
	fmul	10, 24,  0

	fmul	11, 24,  1
	lfd	   20,	72(%r5)
	fnmsub	 6, 13,  1,  6
	dcbt	 %r17,	 %r29
	fmadd	 7, 13,  0,  7
	lfd	   21,	80(%r5)
	fnmsub	 8, 19,  1,  8

	add	%r30, %r30, %r6
	lfd	   26,	120(%r5)
	add	%r9, %r9, %r7

	fmadd	 9, 19,  0,  9
	lfd	   27,	128(%r5)
	fnmsub	10, 25,  1, 10
	dcbt	 %r18,	 %r29
	fmadd	11, 25,  0, 11
	lfd	   4,	32 ( % r12)
	mr	 %r13,	 %r10

	fmadd	 6, 14,  2,  6
	lfd	   5,	40 ( % r12)
	fmadd	 7, 14,  3,  7

	fmadd	 8, 20,  2,  8
	lfd	   16,	40(%r5)
	fmadd	 9, 20,  3,  9
	dcbt	 %r16,	 %r30
	fmadd	10, 26,  2, 10
	lfd	   17,	48(%r5)
	fmadd	11, 26,  3, 11
	dcbt	 %r17,	 %r30
	fnmsub	 6, 15,  3,  6
	lfd	   22,	88(%r5)
	fmadd	 7, 15,  2,  7

	fnmsub	 8, 21,  3,  8
	lfd	   23,	96(%r5)
	fmadd	 9, 21,  2,  9
	dcbt	 %r19,	 %r29
	fnmsub	10, 27,  3, 10
	lfd	   28,	136(%r5)
	fmadd	11, 27,  2, 11

	fmadd	 6, 16,  4,  6
	lfd	   29,	144(%r5)
	fmadd	 7, 16,  5,  7

	fmadd	 8, 22,  4,  8
	lfd	   0,	0 ( % r11)
	fmadd	 9, 22,  5,  9
	dcbt	 %r20,	 %r29
	fmadd	10, 28,  4, 10
	lfd	   1,	8 ( % r11)
	fmadd	11, 28,  5, 11

	fnmsub	 6, 17,  5,  6
	lfd	   12,	160(%r5)
	fmadd	 7, 17,  4,  7

	fnmsub	 8, 23,  5,  8
	lfd	   13,	168(%r5)
	fmadd	 9, 23,  4,  9
	lfd	   18,	208(%r5)
	fnmsub	10, 29,  5, 10
	lfd	   19,	216(%r5)
	fmadd	11, 29,  4, 11
	lfd	   24,	256(%r5)
	mr	 %r12,	 %r30
	stfd	   6,	0(%r14)
	subic.	 3,	3,	1
	lfd	   25,	264(%r5)
	stfd	   7,	8(%r14)
	lfd	   2,	16 ( % r11)
	stfd	   8,	16(%r14)
	lfd	   3,	24 ( % r11)
	stfd	   9,	24(%r14)
	lfd	   14,	176(%r5)
	stfd	   10,	32(%r14)
	lfd	   15,	184(%r5)
	stfd	   11,	40(%r14)
	stfd	   30,	48(%r14)
	stfd	   30,	56(%r14)
	fmul	 6, 12,  0
	dcbt	 %r21,	 %r29
	fmul	 7, 12,  1
	lwz	%r31, 456(%r5)
	fmul	 8, 18,  0

	fmul	 9, 18,  1
	lwz	%r10, 460(%r5)
	fmul	10, 24,  0

	fmul	11, 24,  1
	lfd	   20,	224(%r5)
	fnmsub	 6, 13,  1,  6
	dcbt	 %r22,	 %r29
	fmadd	 7, 13,  0,  7
	lfd	   21,	232(%r5)
	fnmsub	 8, 19,  1,  8

	add	%r31, %r31, %r6
	lfd	   26,	272(%r5)
	add	%r10, %r10, %r7

	fmadd	 9, 19,  0,  9
	lfd	   27,	280(%r5)
	fnmsub	10, 25,  1, 10
	dcbt	 %r23,	 %r29
	fmadd	11, 25,  0, 11
	lfd	   4,	32 ( % r11)
	mr	 %r14,	 %r9

	fmadd	 6, 14,  2,  6
	lfd	   5,	40 ( % r11)
	fmadd	 7, 14,  3,  7

	fmadd	 8, 20,  2,  8
	lfd	   16,	192(%r5)
	fmadd	 9, 20,  3,  9
	dcbt	 %r27,	 %r31
	fmadd	10, 26,  2, 10
	lfd	   17,	200(%r5)
	fmadd	11, 26,  3, 11

	fnmsub	 6, 15,  3,  6
	lfd	   22,	240(%r5)
	fmadd	 7, 15,  2,  7

	fnmsub	 8, 21,  3,  8
	lfd	   23,	248(%r5)
	fmadd	 9, 21,  2,  9
	dcbt	 %r24,	 %r29
	fnmsub	10, 27,  3, 10
	lfd	   28,	288(%r5)
	fmadd	11, 27,  2, 11

	fmadd	 6, 16,  4,  6
	lfd	   29,	296(%r5)
	fmadd	 7, 16,  5,  7

	fmadd	 8, 22,  4,  8
	lfd	   0,	0 ( % r12)
	fmadd	 9, 22,  5,  9
	dcbt	 %r25,	 %r29
	fmadd	10, 28,  4, 10
	lfd	   1,	8 ( % r12)
	fmadd	11, 28,  5, 11

	fnmsub	 6, 17,  5,  6
	lfd	   12,	312(%r5)
	fmadd	 7, 17,  4,  7

	fnmsub	 8, 23,  5,  8
	lfd	   13,	320(%r5)
	fmadd	 9, 23,  4,  9
	lfd	   18,	360(%r5)
	fnmsub	10, 29,  5, 10
	lfd	   19,	368(%r5)
	fmadd	11, 29,  4, 11
	lfd	   24,	408(%r5)
	mr	 %r11,	 %r31
	stfd	   6,	0(%r13)

	lfd	   25,	416(%r5)

	stfd	   7,	8(%r13)

	lfd	   2,	16 ( % r12)

	stfd	   8,	16(%r13)
	addi %r29, %r29, 304
	lfd	   3,	24 ( % r12)

	stfd	   9,	24(%r13)

	lfd	   14,	328(%r5)

	stfd	   10,	32(%r13)

	lfd	   15,	336(%r5)

	stfd	   11,	40(%r13)
	addi %r5, %r5, 304
	stfd	   30,	48(%r13)
	stfd	   30,	56(%r13)
	bf gt,	lab2
	b  dirac_cmv_jcw_agg_lab1
lab2:
	fmul	 6, 12,  0

	fmul	 7, 12,  1
	lwz	%r30, 304(%r5)
	fmul	 8, 18,  0

	fmul	 9, 18,  1
	lwz	%r9, 308(%r5)
	fmul	10, 24,  0

	fmul	11, 24,  1
	lfd	   20,	72(%r5)
	fnmsub	 6, 13,  1,  6

	fmadd	 7, 13,  0,  7
	lfd	   21,	80(%r5)
	fnmsub	 8, 19,  1,  8

	add	%r30, %r30, %r6
	lfd	   26,	120(%r5)
	add	%r9, %r9, %r7

	fmadd	 9, 19,  0,  9
	lfd	   27,	128(%r5)
	fnmsub	10, 25,  1, 10

	fmadd	11, 25,  0, 11
	lfd	   4,	32 ( % r12)
	mr	 %r13,	 %r10

	fmadd	 6, 14,  2,  6
	lfd	   5,	40 ( % r12)
	fmadd	 7, 14,  3,  7

	fmadd	 8, 20,  2,  8
	lfd	   16,	40(%r5)
	fmadd	 9, 20,  3,  9

	fmadd	10, 26,  2, 10
	lfd	   17,	48(%r5)
	fmadd	11, 26,  3, 11

	fnmsub	 6, 15,  3,  6
	lfd	   22,	88(%r5)
	fmadd	 7, 15,  2,  7

	fnmsub	 8, 21,  3,  8
	lfd	   23,	96(%r5)
	fmadd	 9, 21,  2,  9

	fnmsub	10, 27,  3, 10
	lfd	   28,	136(%r5)
	fmadd	11, 27,  2, 11

	fmadd	 6, 16,  4,  6
	lfd	   29,	144(%r5)
	fmadd	 7, 16,  5,  7

	fmadd	 8, 22,  4,  8
	lfd	   0,	0 ( % r11)
	fmadd	 9, 22,  5,  9

	fmadd	10, 28,  4, 10
	lfd	   1,	8 ( % r11)
	fmadd	11, 28,  5, 11

	fnmsub	 6, 17,  5,  6
	lfd	   12,	160(%r5)
	fmadd	 7, 17,  4,  7

	fnmsub	 8, 23,  5,  8
	lfd	   13,	168(%r5)
	fmadd	 9, 23,  4,  9

	fnmsub	10, 29,  5, 10
	lfd	   18,	208(%r5)
	fmadd	11, 29,  4, 11

	mr	 %r12,	 %r30
	lfd	   19,	216(%r5)
	subic.	 3,	3,	1

	lfd	   24,	256(%r5)
	stfd	   6,	0(%r14)
	lfd	   25,	264(%r5)
	stfd	   7,	8(%r14)
	lfd	   2,	16 ( % r11)
	stfd	   8,	16(%r14)
	lfd	   3,	24 ( % r11)
	stfd	   9,	24(%r14)
	lfd	   14,	176(%r5)
	stfd	   10,	32(%r14)
	lfd	   15,	184(%r5)
	stfd	   11,	40(%r14)
	stfd	   30,	48(%r14)
	stfd	   30,	56(%r14)
	fmul	 6, 12,  0

	fmul	 7, 12,  1
	lwz	%r31, 456(%r5)
	fmul	 8, 18,  0

	fmul	 9, 18,  1
	lwz	%r10, 460(%r5)
	fmul	10, 24,  0

	fmul	11, 24,  1
	lfd	   20,	224(%r5)
	fnmsub	 6, 13,  1,  6

	fmadd	 7, 13,  0,  7
	lfd	   21,	232(%r5)
	fnmsub	 8, 19,  1,  8

	add	%r31, %r31, %r6
	lfd	   26,	272(%r5)
	add	%r10, %r10, %r7

	fmadd	 9, 19,  0,  9
	lfd	   27,	280(%r5)
	fnmsub	10, 25,  1, 10

	fmadd	11, 25,  0, 11
	lfd	   4,	32 ( % r11)
	mr	 %r14,	 %r9

	fmadd	 6, 14,  2,  6
	lfd	   5,	40 ( % r11)
	fmadd	 7, 14,  3,  7

	fmadd	 8, 20,  2,  8
	lfd	   16,	192(%r5)
	fmadd	 9, 20,  3,  9

	fmadd	10, 26,  2, 10
	lfd	   17,	200(%r5)
	fmadd	11, 26,  3, 11

	fnmsub	 6, 15,  3,  6
	lfd	   22,	240(%r5)
	fmadd	 7, 15,  2,  7

	fnmsub	 8, 21,  3,  8
	lfd	   23,	248(%r5)
	fmadd	 9, 21,  2,  9

	fnmsub	10, 27,  3, 10
	lfd	   28,	288(%r5)
	fmadd	11, 27,  2, 11

	fmadd	 6, 16,  4,  6
	lfd	   29,	296(%r5)
	fmadd	 7, 16,  5,  7

	fmadd	 8, 22,  4,  8

	fmadd	 9, 22,  5,  9

	fmadd	10, 28,  4, 10

	fmadd	11, 28,  5, 11

	fnmsub	 6, 17,  5,  6

	fmadd	 7, 17,  4,  7

	fnmsub	 8, 23,  5,  8

	fmadd	 9, 23,  4,  9

	fnmsub	10, 29,  5, 10

	fmadd	11, 29,  4, 11

	mr	 %r11,	 %r31









	addi %r29, %r29, 304


	stfd	   6,	0(%r13)

	stfd	   7,	8(%r13)

	stfd	   8,	16(%r13)

	stfd	   9,	24(%r13)

	stfd	   10,	32(%r13)
	addi %r5, %r5, 304
	stfd	   11,	40(%r13)
	stfd	   30,	48(%r13)
	stfd	   30,	56(%r13)
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
