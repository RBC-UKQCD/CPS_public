	.section	".text"
	.align 2
	.globl pt_copy_pad
	.type	pt_copy_pad,@function
pt_copy_pad:
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
	lwz	%r18, 0(%r4)
	lwz	%r22, 4(%r4)
	dcbt	 %r11,	 %r4
	lwz	%r19, 12(%r4)
	lwz	%r23, 16(%r4)
	add	%r18, %r18, %r5
	add	%r22, %r22, %r6
	add	%r19, %r19, %r5
	add	%r23, %r23, %r6
	dcbt	 %r12,	 %r4
	or.	 3,	3,	3
	mr	 %r25,	 %r23
	bf gt,	lab0
	lfd	   0,	0 ( % r18)
	lfd	   1,	8 ( % r18)
	mr	 %r25,	 %r23
	subic.	 3,	3,	1
	bf gt,	lab2
	b  lab1
lab1:
	lfd	   2,	16 ( % r18)
#<<<<<< Break 0
	lwz	%r20, 24(%r4)
	stfd	   6,	0(%r23)
	lwz	%r24, 28(%r4)
	lfd	   3,	24 ( % r18)
	add	%r20, %r20, %r5
	stfd	   7,	8(%r23)
	add	%r24, %r24, %r6
	stfd	   8,	16(%r23)
	stfd	   9,	24(%r23)
	stfd	   10,	32(%r23)
	stfd	   11,	40(%r23)
	stfd	   11,	48(%r23)
	stfd	   11,	56(%r23)
#<<<<<< Break 1
	lfd	   4,	32 ( % r18)
	dcbt	 %r9,	 %r20
	lfd	   5,	40 ( % r18)
	dcbt	 %r10,	 %r20
#<<<<<< Break 2
	lfd	   6,	0 ( % r19)
	mr	 %r23,	 %r25

	lfd	   7,	8 ( % r19)
	dcbt	 %r13,	 %r4
#<<<<<< Break 3
	mr	 %r18,	 %r20
	subic.	 3,	3,	1
#<<<<<< Break 4
	lwz	%r21, 36(%r4)
	lfd	   8,	16 ( % r19)
	lwz	%r25, 40(%r4)
	dcbt	 %r14,	 %r4
	add	%r21, %r21, %r5
	stfd	   0,	0(%r22)
	add	%r25, %r25, %r6

	lfd	   9,	24 ( % r19)
	stfd	   1,	8(%r22)
	stfd	   2,	16(%r22)
	stfd	   3,	24(%r22)
	stfd	   4,	32(%r22)
	stfd	   5,	40(%r22)
	stfd	   5,	48(%r22)
	stfd	   5,	56(%r22)
#<<<<<< Break 5
	lfd	   10,	32 ( % r19)
	dcbt	 %r9,	 %r21
	lfd	   11,	40 ( % r19)
	dcbt	 %r10,	 %r21
#<<<<<< Break 6
	lfd	   0,	0 ( % r20)
	mr	 %r22,	 %r24

	addi %r4, %r4, 24 
	lfd	   1,	8 ( % r20)
	mr	 %r19,	 %r21







	bf gt,	lab2
	b  lab1
lab2:
	lfd	   2,	16 ( % r18)
#<<<<<< Break 0
	stfd	   6,	0(%r23)
	lfd	   3,	24 ( % r18)
	stfd	   7,	8(%r23)
	stfd	   8,	16(%r23)
	stfd	   9,	24(%r23)
	stfd	   10,	32(%r23)
	stfd	   11,	40(%r23)
	stfd	   11,	48(%r23)
	stfd	   11,	56(%r23)
#<<<<<< Break 1
	lfd	   4,	32 ( % r18)
	mr	 %r19,	 %r21
	lfd	   5,	40 ( % r18)
#<<<<<< Break 2
	lfd	   6,	0 ( % r19)
	mr	 %r23,	 %r25
	lfd	   7,	8 ( % r19)
#<<<<<< Break 3
#<<<<<< Break 4
	lfd	   8,	16 ( % r19)
	stfd	   0,	0(%r22)
	lfd	   9,	24 ( % r19)
	stfd	   1,	8(%r22)
	stfd	   2,	16(%r22)
	stfd	   3,	24(%r22)
	stfd	   4,	32(%r22)
	stfd	   5,	40(%r22)
	stfd	   5,	48(%r22)
	stfd	   5,	56(%r22)
#<<<<<< Break 5
	lfd	   10,	32 ( % r19)
	lfd	   11,	40 ( % r19)
#<<<<<< Break 6







	stfd	   6,	0(%r23)
	stfd	   7,	8(%r23)
	stfd	   8,	16(%r23)
	stfd	   9,	24(%r23)
	stfd	   10,	32(%r23)
	stfd	   11,	40(%r23)
	stfd	   11,	48(%r23)
	stfd	   11,	56(%r23)
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
