	.section	".text"
	.align 2
	.globl asq_cross2dag
	.type	asq_cross2dag,@function
asq_cross2dag:
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
	li   10,	96
	li   11,	128
	li   12,	160
	li   13,	192
	li   14,	224
	li   15,	256
	li   16,	288
	li   17,	320
	li   18,	352
	li   19,	384
	li   20,	416
	li   21,	448
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
	addi %r5, %r5, -288 
	lfd	   0,	0(%r7)
	lfd	   1,	0(%r3)
	lfd	   2,	8(%r3)
	lfd	   3,	16(%r3)
	lfd	   4,	24(%r3)
	lfd	   5,	32(%r3)
	lfd	   6,	40(%r3)
	or.	 6,	6,	6
	subic.	 6,	6,	1
	lfd	   7,	0(%r4)
	lfd	   8,	8(%r4)
	lfd	   9,	16(%r4)
	lfd	   10,	24(%r4)
	lfd	   11,	32(%r4)
	lfd	   12,	40(%r4)
	bf gt,	lab0
	b  lab1
lab1:
	lfd	   11,	32(%r4)
	addi %r5, %r5, 288 
	lfd	   12,	40(%r4)
	fmul	1,	0,	1
	lfd	   13,	0(%r5)
	fmul	2,	0,	2
	lfd	   14,	8(%r5)
	fmul	3,	0,	3
	lfd	   15,	16(%r5)
	fmul	4,	0,	4
	lfd	   16,	24(%r5)
	fmul	5,	0,	5
	lfd	   17,	32(%r5)
	fmul	6,	0,	6
	lfd	   18,	40(%r5)
#<<<<<< Break 0
	lfd	   19,	48(%r5)
	fmadd	13,	1,	7,	13
	lfd	   20,	56(%r5)
	fmadd	14,	2,	7,	14
	lfd	   21,	64(%r5)
	fmadd	15,	1,	9,	15
	dcbt	 %r10,	 %r3
	fmadd	16,	2,	9,	16
	lfd	   22,	72(%r5)
	fmadd	17,	1,	11,	17
	lfd	   23,	80(%r5)
	fmadd	18,	2,	11,	18
	lfd	   24,	88(%r5)
	fmadd	13,	2,	8,	13
	dcbt	 %r11,	 %r3
	fnmsub	14,	1,	8,	14
	fmadd	15,	2,	10,	15
	fnmsub	16,	1,	10,	16
	fmadd	17,	2,	12,	17
	fnmsub	18,	1,	12,	18
#<<<<<< Break 1
	fmadd	19,	3,	7,	19
	lfd	   25,	96(%r5)
	fmadd	20,	4,	7,	20
	stfd	   13,	0(%r5)
	fmadd	21,	3,	9,	21
	lfd	   26,	104(%r5)
	fmadd	22,	4,	9,	22
	stfd	   14,	8(%r5)
	fmadd	23,	3,	11,	23
	lfd	   27,	112(%r5)
	fmadd	24,	4,	11,	24
	stfd	   15,	16(%r5)
	fmadd	19,	4,	8,	19
	dcbt	 %r12,	 %r3
	fnmsub	20,	3,	8,	20
	stfd	   16,	24(%r5)
	fmadd	21,	4,	10,	21
	lfd	   28,	120(%r5)
	fnmsub	22,	3,	10,	22
	stfd	   17,	32(%r5)
	fmadd	23,	4,	12,	23
	lfd	   29,	128(%r5)
	fnmsub	24,	3,	12,	24
	stfd	   18,	40(%r5)
	lfd	   30,	136(%r5)
	dcbt	 %r10,	 %r4
	lfd	   1,	48(%r3)
	lfd	   2,	56(%r3)
#<<<<<< Break 2
	lfd	   3,	64(%r3)
	fmadd	25,	5,	7,	25
	stfd	   19,	48(%r5)
	fmadd	26,	6,	7,	26
	dcbt	 %r11,	 %r4
	fmadd	27,	5,	9,	27
	stfd	   20,	56(%r5)
	fmadd	28,	6,	9,	28
	lfd	   4,	72(%r3)
	fmadd	29,	5,	11,	29
	stfd	   21,	64(%r5)
	fmadd	30,	6,	11,	30
	lfd	   13,	144(%r5)
	fmadd	25,	6,	8,	25
	stfd	   22,	72(%r5)
	fnmsub	26,	5,	8,	26
	lfd	   14,	152(%r5)
	fmadd	27,	6,	10,	27
	stfd	   23,	80(%r5)
	fnmsub	28,	5,	10,	28
	dcbt	 %r12,	 %r4
	fmadd	29,	6,	12,	29
	stfd	   24,	88(%r5)
	fnmsub	30,	5,	12,	30
	lfd	   15,	160(%r5)
	lfd	   16,	168(%r5)
	lfd	   17,	176(%r5)
	dcbt	 %r13,	 %r5
	lfd	   18,	184(%r5)
#<<<<<< Break 3
	lfd	   5,	80(%r3)
	fmul	1,	0,	1
	stfd	   25,	96(%r5)
	fmul	2,	0,	2
	lfd	   6,	88(%r3)
	fmul	3,	0,	3
	stfd	   26,	104(%r5)
	fmul	4,	0,	4
	dcbt	 %r14,	 %r5
	fmul	5,	0,	5
	stfd	   27,	112(%r5)
	fmul	6,	0,	6
	lfd	   7,	48(%r4)
	stfd	   28,	120(%r5)
	lfd	   8,	56(%r4)
	stfd	   29,	128(%r5)
	lfd	   9,	64(%r4)
	stfd	   30,	136(%r5)
	dcbt	 %r15,	 %r5
	lfd	   10,	72(%r4)
	lfd	   11,	80(%r4)
	lfd	   12,	88(%r4)
	dcbt	 %r16,	 %r5
#<<<<<< Break 7
	lfd	   19,	192(%r5)
	fmadd	13,	1,	7,	13
	lfd	   20,	200(%r5)
	fmadd	14,	2,	7,	14
	lfd	   21,	208(%r5)
	fmadd	15,	1,	9,	15
	dcbt	 %r17,	 %r5
	fmadd	16,	2,	9,	16
	lfd	   22,	216(%r5)
	fmadd	17,	1,	11,	17
	lfd	   23,	224(%r5)
	fmadd	18,	2,	11,	18
	lfd	   24,	232(%r5)
	fmadd	13,	2,	8,	13
	dcbt	 %r18,	 %r5
	fnmsub	14,	1,	8,	14
	fmadd	15,	2,	10,	15
	fnmsub	16,	1,	10,	16
	fmadd	17,	2,	12,	17
	fnmsub	18,	1,	12,	18
#<<<<<< Break 4
	fmadd	19,	3,	7,	19
	lfd	   25,	240(%r5)
	fmadd	20,	4,	7,	20
	stfd	   13,	144(%r5)
	fmadd	21,	3,	9,	21
	lfd	   26,	248(%r5)
	fmadd	22,	4,	9,	22
	stfd	   14,	152(%r5)
	fmadd	23,	3,	11,	23
	lfd	   27,	256(%r5)
	fmadd	24,	4,	11,	24
	stfd	   15,	160(%r5)
	fmadd	19,	4,	8,	19
	dcbt	 %r19,	 %r5
	fnmsub	20,	3,	8,	20
	stfd	   16,	168(%r5)
	fmadd	21,	4,	10,	21
	lfd	   28,	264(%r5)
	fnmsub	22,	3,	10,	22
	stfd	   17,	176(%r5)
	fmadd	23,	4,	12,	23
	lfd	   29,	272(%r5)
	fnmsub	24,	3,	12,	24
	stfd	   18,	184(%r5)
	lfd	   30,	280(%r5)
	dcbt	 %r20,	 %r5
#<<<<<< Break 5
	lfd	   1,	96(%r3)
	fmadd	25,	5,	7,	25
	stfd	   19,	192(%r5)
	fmadd	26,	6,	7,	26
	lfd	   2,	104(%r3)
	fmadd	27,	5,	9,	27
	stfd	   20,	200(%r5)
	fmadd	28,	6,	9,	28
	lfd	   3,	112(%r3)
	fmadd	29,	5,	11,	29
	stfd	   21,	208(%r5)
	fmadd	30,	6,	11,	30
	dcbt	 %r21,	 %r5
	fmadd	25,	6,	8,	25
	stfd	   22,	216(%r5)
	fnmsub	26,	5,	8,	26
	lfd	   4,	120(%r3)
	fmadd	27,	6,	10,	27
	stfd	   23,	224(%r5)
	fnmsub	28,	5,	10,	28
	stfd	   24,	232(%r5)
	fmadd	29,	6,	12,	29
	fnmsub	30,	5,	12,	30
#<<<<<< Break 6
	addi %r4, %r4, 96 
	stfd	   25,	240(%r5)
	lfd	   5,	128(%r3)
	stfd	   26,	248(%r5)
	lfd	   6,	136(%r3)
	stfd	   27,	256(%r5)
#<<<<<< Break 8
	lfd	   7,	0(%r4)
	addi %r3, %r3, 96 
	stfd	   28,	264(%r5)
	subic.	 6,	6,	1
	lfd	   8,	8(%r4)
	stfd	   29,	272(%r5)
	lfd	   9,	16(%r4)
	stfd	   30,	280(%r5)
	lfd	   10,	24(%r4)
	bf gt,	lab0
	b  lab1
lab0:
	lfd	   11,	32(%r4)
	addi %r5, %r5, 288 
	lfd	   12,	40(%r4)
	fmul	1,	0,	1
	lfd	   13,	0(%r5)
	fmul	2,	0,	2
	lfd	   14,	8(%r5)
	fmul	3,	0,	3
	lfd	   15,	16(%r5)
	fmul	4,	0,	4
	lfd	   16,	24(%r5)
	fmul	5,	0,	5
	lfd	   17,	32(%r5)
	fmul	6,	0,	6
	lfd	   18,	40(%r5)
#<<<<<< Break 0
	lfd	   19,	48(%r5)
	fmadd	13,	1,	7,	13
	lfd	   20,	56(%r5)
	fmadd	14,	2,	7,	14
	lfd	   21,	64(%r5)
	fmadd	15,	1,	9,	15
	lfd	   22,	72(%r5)
	fmadd	16,	2,	9,	16
	lfd	   23,	80(%r5)
	fmadd	17,	1,	11,	17
	lfd	   24,	88(%r5)
	fmadd	18,	2,	11,	18
	fmadd	13,	2,	8,	13
	fnmsub	14,	1,	8,	14
	fmadd	15,	2,	10,	15
	fnmsub	16,	1,	10,	16
	fmadd	17,	2,	12,	17
	fnmsub	18,	1,	12,	18
#<<<<<< Break 1
	fmadd	19,	3,	7,	19
	lfd	   25,	96(%r5)
	fmadd	20,	4,	7,	20
	stfd	   13,	0(%r5)
	fmadd	21,	3,	9,	21
	lfd	   26,	104(%r5)
	fmadd	22,	4,	9,	22
	stfd	   14,	8(%r5)
	fmadd	23,	3,	11,	23
	lfd	   27,	112(%r5)
	fmadd	24,	4,	11,	24
	stfd	   15,	16(%r5)
	fmadd	19,	4,	8,	19
	lfd	   28,	120(%r5)
	fnmsub	20,	3,	8,	20
	stfd	   16,	24(%r5)
	fmadd	21,	4,	10,	21
	lfd	   29,	128(%r5)
	fnmsub	22,	3,	10,	22
	stfd	   17,	32(%r5)
	fmadd	23,	4,	12,	23
	lfd	   30,	136(%r5)
	fnmsub	24,	3,	12,	24
	stfd	   18,	40(%r5)
	lfd	   1,	48(%r3)
	lfd	   2,	56(%r3)
#<<<<<< Break 2
	lfd	   3,	64(%r3)
	fmadd	25,	5,	7,	25
	stfd	   19,	48(%r5)
	fmadd	26,	6,	7,	26
	lfd	   4,	72(%r3)
	fmadd	27,	5,	9,	27
	stfd	   20,	56(%r5)
	fmadd	28,	6,	9,	28
	lfd	   13,	144(%r5)
	fmadd	29,	5,	11,	29
	stfd	   21,	64(%r5)
	fmadd	30,	6,	11,	30
	lfd	   14,	152(%r5)
	fmadd	25,	6,	8,	25
	stfd	   22,	72(%r5)
	fnmsub	26,	5,	8,	26
	lfd	   15,	160(%r5)
	fmadd	27,	6,	10,	27
	stfd	   23,	80(%r5)
	fnmsub	28,	5,	10,	28
	lfd	   16,	168(%r5)
	fmadd	29,	6,	12,	29
	stfd	   24,	88(%r5)
	fnmsub	30,	5,	12,	30
	lfd	   17,	176(%r5)
	lfd	   18,	184(%r5)
#<<<<<< Break 3
	lfd	   5,	80(%r3)
	fmul	1,	0,	1
	stfd	   25,	96(%r5)
	fmul	2,	0,	2
	lfd	   6,	88(%r3)
	fmul	3,	0,	3
	stfd	   26,	104(%r5)
	fmul	4,	0,	4
	lfd	   7,	48(%r4)
	fmul	5,	0,	5
	stfd	   27,	112(%r5)
	fmul	6,	0,	6
	lfd	   8,	56(%r4)
	stfd	   28,	120(%r5)
	lfd	   9,	64(%r4)
	stfd	   29,	128(%r5)
	lfd	   10,	72(%r4)
	stfd	   30,	136(%r5)
	lfd	   11,	80(%r4)
	lfd	   12,	88(%r4)
#<<<<<< Break 7
	lfd	   19,	192(%r5)
	fmadd	13,	1,	7,	13
	lfd	   20,	200(%r5)
	fmadd	14,	2,	7,	14
	lfd	   21,	208(%r5)
	fmadd	15,	1,	9,	15
	lfd	   22,	216(%r5)
	fmadd	16,	2,	9,	16
	lfd	   23,	224(%r5)
	fmadd	17,	1,	11,	17
	lfd	   24,	232(%r5)
	fmadd	18,	2,	11,	18
	fmadd	13,	2,	8,	13
	fnmsub	14,	1,	8,	14
	fmadd	15,	2,	10,	15
	fnmsub	16,	1,	10,	16
	fmadd	17,	2,	12,	17
	fnmsub	18,	1,	12,	18
#<<<<<< Break 4
	fmadd	19,	3,	7,	19
	lfd	   25,	240(%r5)
	fmadd	20,	4,	7,	20
	stfd	   13,	144(%r5)
	fmadd	21,	3,	9,	21
	lfd	   26,	248(%r5)
	fmadd	22,	4,	9,	22
	stfd	   14,	152(%r5)
	fmadd	23,	3,	11,	23
	lfd	   27,	256(%r5)
	fmadd	24,	4,	11,	24
	stfd	   15,	160(%r5)
	fmadd	19,	4,	8,	19
	lfd	   28,	264(%r5)
	fnmsub	20,	3,	8,	20
	stfd	   16,	168(%r5)
	fmadd	21,	4,	10,	21
	lfd	   29,	272(%r5)
	fnmsub	22,	3,	10,	22
	stfd	   17,	176(%r5)
	fmadd	23,	4,	12,	23
	lfd	   30,	280(%r5)
	fnmsub	24,	3,	12,	24
	stfd	   18,	184(%r5)
#<<<<<< Break 5
	stfd	   19,	192(%r5)
	fmadd	25,	5,	7,	25
	stfd	   20,	200(%r5)
	fmadd	26,	6,	7,	26
	stfd	   21,	208(%r5)
	fmadd	27,	5,	9,	27
	stfd	   22,	216(%r5)
	fmadd	28,	6,	9,	28
	stfd	   23,	224(%r5)
	fmadd	29,	5,	11,	29
	stfd	   24,	232(%r5)
	fmadd	30,	6,	11,	30
	fmadd	25,	6,	8,	25
	fnmsub	26,	5,	8,	26
	fmadd	27,	6,	10,	27
	fnmsub	28,	5,	10,	28
	fmadd	29,	6,	12,	29
	fnmsub	30,	5,	12,	30
#<<<<<< Break 6
	addi %r4, %r4, 96 
	stfd	   25,	240(%r5)
#<<<<<< Break 8
	stfd	   26,	248(%r5)
	addi %r3, %r3, 96 
	stfd	   27,	256(%r5)
	subic.	 6,	6,	1
	stfd	   28,	264(%r5)
	stfd	   29,	272(%r5)
	stfd	   30,	280(%r5)
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
