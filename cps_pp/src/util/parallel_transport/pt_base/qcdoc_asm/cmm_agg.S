	.section	".text"
	.align 2
	.globl cmm_agg
	.type	cmm_agg,@function
cmm_agg:
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
	li   19,	152
	li   20,	184
	li   21,	216
	li   22,	248
	li   23,	280
	li   24,	312
	li   25,	0
	li   26,	32
	li   27,	64
	li   28,	96
	li   29,	128
	li   30,	160
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
	mr	 %r7,	 %r3
	mr	 %r18,	 %r5
	addi %r13, %r3, -152
	lwz	%r8, 0(%r7)
	lwz	%r10, 4(%r7)
	mulli   9,8,3
	mulli   11,10,3
	add	%r8, %r9, %r4
	add	%r12, %r11, %r5
	add	%r18, %r11, %r5
	or.	 6,	6,	6
	lfd	0,	8(%r7)
	lfd	1,	16(%r7)
	lfd	4,	0(%r8)
	lfd	5,	8(%r8)
	lfd	6,	16(%r8)
	lfd	7,	24(%r8)
	lfd	16,	32(%r8)
	lfd	17,	40(%r8)
	bf gt,	lab0
	b  lab1
lab1:
#<<<<<< Break 0
	lfd	2,	24(%r7)
	addi %r13, %r13, 304 
	lfd	3,	32(%r7)
	fmul	20,  0,  4
	lfd	8,	48(%r8)
	fmul	21,  0,  5
	lfd	9,	56(%r8)
	fmul	22,  0,  6
	lfd	10,	64(%r8)
	fmul	23,  0,  7
	lfd	11,	72(%r8)
	fmul	24,  0, 16
	stfd	26,	96(%r18)
	fmul	25,  0, 17
	lfd	18,	80(%r8)
	fnmsub	20,  1,  5, 20
	lfd	19,	88(%r8)
	fmadd	21,  1,  4, 21
	stfd	27,	104(%r18)
	fnmsub	22,  1,  7, 22
	fmadd	23,  1,  6, 23
	fnmsub	24,  1, 17, 24
	fmadd	25,  1, 16, 25
#<<<<<< Break 1
	fmadd	20,  2,  8, 20
	lwz	%r14, 0(%r13)
	fmadd	21,  2,  9, 21
	lwz	%r16, 4(%r13)
	fmadd	22,  2, 10, 22
	lfd	0,	40(%r7)
	fmadd	23,  2, 11, 23
	stfd	28,	112(%r18)
	fmadd	24,  2, 18, 24
	lfd	1,	48(%r7)
	fmadd	25,  2, 19, 25
	dcbt	 %r19,	 %r7
	fnmsub	20,  3,  9, 20
	lfd	12,	96(%r8)
	fmadd	21,  3,  8, 21
	lfd	13,	104(%r8)
	fnmsub	22,  3, 11, 22
	lfd	14,	112(%r8)
	fmadd	23,  3, 10, 23
	lfd	15,	120(%r8)
	fnmsub	24,  3, 19, 24
	stfd	29,	120(%r18)
	fmadd	25,  3, 18, 25
	dcbt	 %r20,	 %r7
	lfd	16,	128(%r8)
	lfd	17,	136(%r8)
#<<<<<< Break 2
	lfd	2,	56(%r7)
	mulli   15,14,3
	lfd	3,	64(%r7)
	mulli   17,16,3
	stfd	30,	128(%r18)
	fmadd	20,  0, 12, 20
	dcbt	 %r21,	 %r7
	fmadd	21,  0, 13, 21
	lfd	18,	32(%r8)
	fmadd	22,  0, 14, 22
	lfd	19,	40(%r8)
	fmadd	23,  0, 15, 23
	stfd	31,	136(%r18)
	fmadd	24,  0, 16, 24
	fmadd	25,  0, 17, 25
	fnmsub	20,  1, 13, 20
	fmadd	21,  1, 12, 21
	fnmsub	22,  1, 15, 22
	fmadd	23,  1, 14, 23
	fnmsub	24,  1, 17, 24
	fmadd	25,  1, 16, 25
#<<<<<< Break 3
	add	%r14, %r15, %r4
	lfd	0,	72(%r7)
	add	%r18, %r17, %r5
	lfd	1,	80(%r7)
	fmul	26,  2,  4
	stfd	20,	0(%r12)
	fmul	27,  2,  5
	dcbt	 %r22,	 %r7
	fmul	28,  2,  6
	lfd	16,	80(%r8)
	fmul	29,  2,  7
	lfd	17,	88(%r8)
	fmul	30,  2, 18
	dcbt	 %r23,	 %r7
	fmul	31,  2, 19
	dcbt	 %r24,	 %r7
	fnmsub	26,  3,  5, 26
	stfd	21,	8(%r12)
	fmadd	27,  3,  4, 27
	fnmsub	28,  3,  7, 28
	fmadd	29,  3,  6, 29
	fnmsub	30,  3, 19, 30
	fmadd	31,  3, 18, 31
#<<<<<< Break 4
	fmadd	26,  0,  8, 26
	stfd	22,	16(%r12)
	fmadd	27,  0,  9, 27
	lfd	2,	88(%r7)
	fmadd	28,  0, 10, 28
	dcbt	 %r25,	 %r14
	fmadd	29,  0, 11, 29
	lfd	3,	96(%r7)
	fmadd	30,  0, 16, 30
	lfd	18,	128(%r8)
	fmadd	31,  0, 17, 31
	lfd	19,	136(%r8)
	fnmsub	26,  1,  9, 26
	stfd	23,	24(%r12)
	fmadd	27,  1,  8, 27
	fnmsub	28,  1, 11, 28
	fmadd	29,  1, 10, 29
	fnmsub	30,  1, 17, 30
	fmadd	31,  1, 16, 31
#<<<<<< Break 5
	fmadd	26,  2, 12, 26
	stfd	24,	32(%r12)
	fmadd	27,  2, 13, 27
	lfd	0,	104(%r7)
	fmadd	28,  2, 14, 28
	dcbt	 %r26,	 %r14
	fmadd	29,  2, 15, 29
	lfd	1,	112(%r7)
	fmadd	30,  2, 18, 30
	lfd	16,	32(%r8)
	fmadd	31,  2, 19, 31
	lfd	17,	40(%r8)
	fnmsub	26,  3, 13, 26
	stfd	25,	40(%r12)
	fmadd	27,  3, 12, 27
	fnmsub	28,  3, 15, 28
	fmadd	29,  3, 14, 29
	fnmsub	30,  3, 19, 30
	fmadd	31,  3, 18, 31
#<<<<<< Break 6
	fmul	20,  0,  4
	stfd	26,	48(%r12)
	fmul	21,  0,  5
	lfd	2,	120(%r7)
	fmul	22,  0,  6
	dcbt	 %r27,	 %r14
	fmul	23,  0,  7
	lfd	3,	128(%r7)
	fmul	24,  0, 16
	lfd	18,	80(%r8)
	fmul	25,  0, 17
	lfd	19,	88(%r8)
	fnmsub	20,  1,  5, 20
	stfd	27,	56(%r12)
	fmadd	21,  1,  4, 21
	fnmsub	22,  1,  7, 22
	fmadd	23,  1,  6, 23
	fnmsub	24,  1, 17, 24
	fmadd	25,  1, 16, 25
#<<<<<< Break 7
	fmadd	20,  2,  8, 20
	stfd	28,	64(%r12)
	fmadd	21,  2,  9, 21
	lfd	0,	136(%r7)
	fmadd	22,  2, 10, 22
	dcbt	 %r28,	 %r14
	fmadd	23,  2, 11, 23
	lfd	1,	144(%r7)
	fmadd	24,  2, 18, 24
	lfd	16,	128(%r8)
	fmadd	25,  2, 19, 25
	lfd	17,	136(%r8)
	fnmsub	20,  3,  9, 20
	stfd	29,	72(%r12)
	fmadd	21,  3,  8, 21
	fnmsub	22,  3, 11, 22
	fmadd	23,  3, 10, 23
	fnmsub	24,  3, 19, 24
	fmadd	25,  3, 18, 25
#<<<<<< Break 8
	fmadd	20,  0, 12, 20
	stfd	30,	80(%r12)
	fmadd	21,  0, 13, 21
	lfd	2,	8(%r13)
	fmadd	22,  0, 14, 22
	dcbt	 %r29,	 %r14
	fmadd	23,  0, 15, 23
	lfd	3,	16(%r13)
	fmadd	24,  0, 16, 24
	lfd	4,	0(%r14)
	fmadd	25,  0, 17, 25
	lfd	5,	8(%r14)
	fnmsub	20,  1, 13, 20
	lfd	6,	16(%r14)
	fmadd	21,  1, 12, 21
	stfd	31,	88(%r12)
	fnmsub	22,  1, 15, 22
	dcbt	 %r30,	 %r14
	fmadd	23,  1, 14, 23
	lfd	7,	24(%r14)
	fnmsub	24,  1, 17, 24
	lfd	18,	32(%r14)
	fmadd	25,  1, 16, 25
	lfd	19,	40(%r14)
#<<<<<< Break 9
	lfd	0,	24(%r13)
	addi %r7, %r7, 304 
	stfd	20,	96(%r12)
	fmul	26,  2,  4
	lfd	1,	32(%r13)
	fmul	27,  2,  5
	lfd	8,	48(%r14)
	fmul	28,  2,  6
	lfd	9,	56(%r14)
	fmul	29,  2,  7
	lfd	10,	64(%r14)
	fmul	30,  2, 18
	lfd	11,	72(%r14)
	fmul	31,  2, 19
	lfd	16,	80(%r14)
	fnmsub	26,  3,  5, 26
	stfd	21,	104(%r12)
	fmadd	27,  3,  4, 27
	lfd	17,	88(%r14)
	fnmsub	28,  3,  7, 28
	fmadd	29,  3,  6, 29
	fnmsub	30,  3, 19, 30
	fmadd	31,  3, 18, 31
#<<<<<< Break 10
	fmadd	26,  0,  8, 26
	lwz	%r8, 0(%r7)
	fmadd	27,  0,  9, 27
	lwz	%r10, 4(%r7)
	fmadd	28,  0, 10, 28
	lfd	2,	40(%r13)
	fmadd	29,  0, 11, 29
	lfd	3,	48(%r13)
	fmadd	30,  0, 16, 30
	stfd	22,	112(%r12)
	fmadd	31,  0, 17, 31
	dcbt	 %r19,	 %r13
	fnmsub	26,  1,  9, 26
	lfd	12,	96(%r14)
	fmadd	27,  1,  8, 27
	lfd	13,	104(%r14)
	fnmsub	28,  1, 11, 28
	lfd	14,	112(%r14)
	fmadd	29,  1, 10, 29
	lfd	15,	120(%r14)
	fnmsub	30,  1, 17, 30
	dcbt	 %r20,	 %r13
	fmadd	31,  1, 16, 31
	stfd	23,	120(%r12)
	lfd	18,	128(%r14)
	lfd	19,	136(%r14)
#<<<<<< Break 11
	lfd	0,	56(%r13)
	mulli   9,8,3
	lfd	1,	64(%r13)
	mulli   11,10,3
	dcbt	 %r21,	 %r13
	fmadd	26,  2, 12, 26
	stfd	24,	128(%r12)
	fmadd	27,  2, 13, 27
	lfd	16,	32(%r14)
	fmadd	28,  2, 14, 28
	lfd	17,	40(%r14)
	fmadd	29,  2, 15, 29
	stfd	25,	136(%r12)
	fmadd	30,  2, 18, 30
	fmadd	31,  2, 19, 31
	fnmsub	26,  3, 13, 26
	fmadd	27,  3, 12, 27
	fnmsub	28,  3, 15, 28
	fmadd	29,  3, 14, 29
	fnmsub	30,  3, 19, 30
	fmadd	31,  3, 18, 31
#<<<<<< Break 12
	add	%r8, %r9, %r4
	lfd	2,	72(%r13)
	add	%r12, %r11, %r5
	lfd	3,	80(%r13)
	fmul	20,  0,  4
	dcbt	 %r22,	 %r13
	fmul	21,  0,  5
	stfd	26,	0(%r18)
	fmul	22,  0,  6
	lfd	18,	80(%r14)
	fmul	23,  0,  7
	lfd	19,	88(%r14)
	fmul	24,  0, 16
	dcbt	 %r23,	 %r13
	fmul	25,  0, 17
	dcbt	 %r24,	 %r13
	fnmsub	20,  1,  5, 20
	stfd	27,	8(%r18)
	fmadd	21,  1,  4, 21
	fnmsub	22,  1,  7, 22
	fmadd	23,  1,  6, 23
	fnmsub	24,  1, 17, 24
	fmadd	25,  1, 16, 25
#<<<<<< Break 13
	fmadd	20,  2,  8, 20
	lfd	0,	88(%r13)
	fmadd	21,  2,  9, 21
	stfd	28,	16(%r18)
	fmadd	22,  2, 10, 22
	dcbt	 %r25,	 %r8
	fmadd	23,  2, 11, 23
	lfd	1,	96(%r13)
	fmadd	24,  2, 18, 24
	lfd	16,	128(%r14)
	fmadd	25,  2, 19, 25
	lfd	17,	136(%r14)
	fnmsub	20,  3,  9, 20
	stfd	29,	24(%r18)
	fmadd	21,  3,  8, 21
	fnmsub	22,  3, 11, 22
	fmadd	23,  3, 10, 23
	fnmsub	24,  3, 19, 24
	fmadd	25,  3, 18, 25
#<<<<<< Break 14
	fmadd	20,  0, 12, 20
	lfd	2,	104(%r13)
	fmadd	21,  0, 13, 21
	stfd	30,	32(%r18)
	fmadd	22,  0, 14, 22
	dcbt	 %r26,	 %r8
	fmadd	23,  0, 15, 23
	lfd	3,	112(%r13)
	fmadd	24,  0, 16, 24
	lfd	18,	32(%r14)
	fmadd	25,  0, 17, 25
	lfd	19,	40(%r14)
	fnmsub	20,  1, 13, 20
	stfd	31,	40(%r18)
	fmadd	21,  1, 12, 21
	fnmsub	22,  1, 15, 22
	fmadd	23,  1, 14, 23
	fnmsub	24,  1, 17, 24
	fmadd	25,  1, 16, 25
#<<<<<< Break 15
	fmul	26,  2,  4
	lfd	0,	120(%r13)
	fmul	27,  2,  5
	stfd	20,	48(%r18)
	fmul	28,  2,  6
	dcbt	 %r27,	 %r8
	fmul	29,  2,  7
	lfd	1,	128(%r13)
	fmul	30,  2, 18
	lfd	16,	80(%r14)
	fmul	31,  2, 19
	lfd	17,	88(%r14)
	fnmsub	26,  3,  5, 26
	stfd	21,	56(%r18)
	fmadd	27,  3,  4, 27
	fnmsub	28,  3,  7, 28
	fmadd	29,  3,  6, 29
	fnmsub	30,  3, 19, 30
	fmadd	31,  3, 18, 31
#<<<<<< Break 16
	fmadd	26,  0,  8, 26
	lfd	2,	136(%r13)
	fmadd	27,  0,  9, 27
	stfd	22,	64(%r18)
	fmadd	28,  0, 10, 28
	dcbt	 %r28,	 %r8
	fmadd	29,  0, 11, 29
	lfd	3,	144(%r13)
	fmadd	30,  0, 16, 30
	lfd	18,	128(%r14)
	fmadd	31,  0, 17, 31
	lfd	19,	136(%r14)
	fnmsub	26,  1,  9, 26
	stfd	23,	72(%r18)
	fmadd	27,  1,  8, 27
	fnmsub	28,  1, 11, 28
	fmadd	29,  1, 10, 29
	fnmsub	30,  1, 17, 30
	fmadd	31,  1, 16, 31
#<<<<<< Break 17
	subic.	 6,	6,	1
	lfd	0,	8(%r7)
	fmadd	26,  2, 12, 26
	stfd	24,	80(%r18)
	fmadd	27,  2, 13, 27
	dcbt	 %r29,	 %r8
	fmadd	28,  2, 14, 28
	lfd	1,	16(%r7)
	fmadd	29,  2, 15, 29
	lfd	4,	0(%r8)
	fmadd	30,  2, 18, 30
	lfd	5,	8(%r8)
	fmadd	31,  2, 19, 31
	lfd	6,	16(%r8)
	fnmsub	26,  3, 13, 26
	dcbt	 %r30,	 %r8
	fmadd	27,  3, 12, 27
	stfd	25,	88(%r18)
	fnmsub	28,  3, 15, 28
	lfd	7,	24(%r8)
	fmadd	29,  3, 14, 29
	lfd	16,	32(%r8)
	fnmsub	30,  3, 19, 30
	lfd	17,	40(%r8)
	fmadd	31,  3, 18, 31
#<<<<<< Break 18
	bf gt,	lab0
	b  lab1
lab0:
	stfd	26,	96(%r18)
	stfd	27,	104(%r18)
	stfd	28,	112(%r18)
	stfd	29,	120(%r18)
	stfd	30,	128(%r18)
	stfd	31,	136(%r18)
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
