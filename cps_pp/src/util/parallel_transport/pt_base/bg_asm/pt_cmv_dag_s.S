	.section	".text"
	.align 2
	.globl pt_cmv_dag_s
	.type	pt_cmv_dag_s,@function
pt_cmv_dag_s:
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
	bf gt,	pt_cmv_dag_s_lab0
	lfs	   18,	0(%r27)
	lfs	   19,	4(%r27)
	lfs	   20,	8(%r27)
	lfs	   21,	12(%r27)
	lfs	   22,	16(%r27)
	lfs	   23,	20(%r27)
	add	%r18, %r5, %r9
	lfs	   0,	0 ( %r19)
	lfs	   1,	4 ( %r19)
	mr	 %r26,	 %r24
	subic.	 3,	3,	1
	bf gt,	pt_cmv_dag_s_lab2
	b  pt_cmv_dag_s_lab1
pt_cmv_dag_s_lab1:
#<<<<<< Break 0
	lwz	%r21, 24(%r4)
	fmuls	 6, 18,  0
	lwz	%r25, 28(%r4)
	dcbt	 %r11,	 %r30
	lwz	%r29, 32(%r4)
	fmuls	 7, 18,  1
	add	%r21, %r21, %r6
	lfs	   24,	24(%r27)
	add	%r25, %r25, %r7
	fmuls	 8, 20,  0
	add	%r29, %r29, %r5

	fmuls	 9, 20,  1
	stfs	   12,	0(%r24)
	fmuls	10, 22,  0

	fmuls	11, 22,  1
	lfs	   25,	28(%r27)
	fmadds	 6, 19,  1,  6
	stfs	   13,	4(%r24)
	fnmsubs	 7, 19,  0,  7
	lfs	   26,	32(%r27)
	fmadds	 8, 21,  1,  8
	stfs	   14,	8(%r24)
	fnmsubs	 9, 21,  0,  9
	lfs	   27,	36(%r27)
	fmadds	10, 23,  1, 10
	stfs	   15,	12(%r24)
	fnmsubs	11, 23,  0, 11
	lfs	   28,	40(%r27)
	stfs	   16,	16(%r24)
	lfs	   29,	44(%r27)
	stfs	   17,	20(%r24)
	lfs	   2,	8 ( %r19)
	lfs	   3,	12 ( %r19)
#<<<<<< Break 1
	lfs	   18,	48(%r27)
	fmadds	 6, 24,  2,  6
	dcbt	 %r9,	 %r29
	fmadds	 7, 24,  3,  7
	lfs	   19,	52(%r27)
	fmadds	 8, 26,  2,  8

	fmadds	 9, 26,  3,  9
	lfs	   20,	56(%r27)
	fmadds	10, 28,  2, 10

	fmadds	11, 28,  3, 11
	lfs	   21,	60(%r27)
	fmadds	 6, 25,  3,  6

	fnmsubs	 7, 25,  2,  7
	lfs	   22,	64(%r27)
	fmadds	 8, 27,  3,  8
	dcbt	 %r10,	 %r29
	fnmsubs	 9, 27,  2,  9
	lfs	   23,	68(%r27)
	fmadds	10, 29,  3, 10

	fnmsubs	11, 29,  2, 11
	lfs	   4,	16 ( %r19)

	lfs	   5,	20 ( %r19)

#<<<<<< Break 2
	lfs	   24,	0(%r28)
	mr	 %r24,	 %r26
	fmadds	 6, 18,  4,  6
	dcbt	 %r9,	 %r21
	fmadds	 7, 18,  5,  7
	lfs	   25,	4(%r28)
	fmadds	 8, 20,  4,  8
	dcbt	 %r10,	 %r21
	fmadds	 9, 20,  5,  9
	lfs	   26,	8(%r28)
	fmadds	10, 22,  4, 10

	fmadds	11, 22,  5, 11
	lfs	   27,	12(%r28)
	fmadds	 6, 19,  5,  6
	dcbt	 %r11,	 %r29
	fnmsubs	 7, 19,  4,  7
	lfs	   28,	16(%r28)
	fmadds	 8, 21,  5,  8

	fnmsubs	 9, 21,  4,  9
	lfs	   29,	20(%r28)
	fmadds	10, 23,  5, 10

	fnmsubs	11, 23,  4, 11
	lfs	   0,	0 ( %r20)

	lfs	   1,	4 ( %r20)

#<<<<<< Break 3
	mr	 %r19,	 %r21
	mr	 %r27,	 %r29
	subic.	 3,	3,	1
#<<<<<< Break 4
	lwz	%r22, 36(%r4)
	fmuls	12, 24,  0
	lwz	%r26, 40(%r4)
	lfs	   18,	24(%r28)
	lwz	%r30, 44(%r4)
	fmuls	13, 24,  1
	add	%r22, %r22, %r6
	dcbt	 %r13,	 %r4
	add	%r26, %r26, %r7
	fmuls	14, 26,  0
	add	%r30, %r30, %r5
	stfs	   6,	0(%r23)
	fmuls	15, 26,  1
	dcbt	 %r14,	 %r4
	fmuls	16, 28,  0
	lfs	   19,	28(%r28)
	fmuls	17, 28,  1

	fmadds	12, 25,  1, 12
	stfs	   7,	4(%r23)
	fnmsubs	13, 25,  0, 13
	lfs	   20,	32(%r28)
	fmadds	14, 27,  1, 14
	stfs	   8,	8(%r23)
	fnmsubs	15, 27,  0, 15
	lfs	   21,	36(%r28)
	fmadds	16, 29,  1, 16
	stfs	   9,	12(%r23)
	fnmsubs	17, 29,  0, 17
	lfs	   22,	40(%r28)
	stfs	   10,	16(%r23)
	lfs	   23,	44(%r28)
	stfs	   11,	20(%r23)
	lfs	   2,	8 ( %r20)
	lfs	   3,	12 ( %r20)
#<<<<<< Break 5
	lfs	   24,	48(%r28)
	fmadds	12, 18,  2, 12
	dcbt	 %r9,	 %r30
	fmadds	13, 18,  3, 13
	lfs	   25,	52(%r28)
	fmadds	14, 20,  2, 14

	fmadds	15, 20,  3, 15
	lfs	   26,	56(%r28)
	fmadds	16, 22,  2, 16

	fmadds	17, 22,  3, 17
	lfs	   27,	60(%r28)
	fmadds	12, 19,  3, 12

	fnmsubs	13, 19,  2, 13
	lfs	   28,	64(%r28)
	fmadds	14, 21,  3, 14
	dcbt	 %r10,	 %r30
	fnmsubs	15, 21,  2, 15
	lfs	   29,	68(%r28)
	fmadds	16, 23,  3, 16

	fnmsubs	17, 23,  2, 17
	lfs	   4,	16 ( %r20)

	lfs	   5,	20 ( %r20)

#<<<<<< Break 6
	lfs	   18,	0(%r29)
	mr	 %r23,	 %r25
	fmadds	12, 24,  4, 12
	addi %r4, %r4, 24 
	dcbt	 %r9,	 %r22
	mr	 %r20,	 %r22
	fmadds	13, 24,  5, 13
	mr	 %r28,	 %r30
	lfs	   19,	4(%r29)

	fmadds	14, 26,  4, 14

	dcbt	 %r10,	 %r22

	fmadds	15, 26,  5, 15

	lfs	   20,	8(%r29)

	fmadds	16, 28,  4, 16



	fmadds	17, 28,  5, 17
	lfs	   21,	12(%r29)
	fmadds	12, 25,  5, 12
	dcbt	 %r11,	 %r30
	fnmsubs	13, 25,  4, 13
	lfs	   22,	16(%r29)
	fmadds	14, 27,  5, 14

	fnmsubs	15, 27,  4, 15
	lfs	   23,	20(%r29)
	fmadds	16, 29,  5, 16

	fnmsubs	17, 29,  4, 17
	lfs	   0,	0 ( %r21)
#<<<<<< Break 7

	lfs	   1,	4 ( %r21)

	bf gt,	pt_cmv_dag_s_lab2
	b  pt_cmv_dag_s_lab1
pt_cmv_dag_s_lab2:
#<<<<<< Break 0
	fmuls	 6, 18,  0
	lfs	   24,	24(%r27)
	fmuls	 7, 18,  1
	stfs	   12,	0(%r24)
	fmuls	 8, 20,  0
	lfs	   25,	28(%r27)
	fmuls	 9, 20,  1
	stfs	   13,	4(%r24)
	fmuls	10, 22,  0
	lfs	   26,	32(%r27)
	fmuls	11, 22,  1
	stfs	   14,	8(%r24)
	fmadds	 6, 19,  1,  6
	lfs	   27,	36(%r27)
	fnmsubs	 7, 19,  0,  7
	stfs	   15,	12(%r24)
	fmadds	 8, 21,  1,  8
	lfs	   28,	40(%r27)
	fnmsubs	 9, 21,  0,  9
	stfs	   16,	16(%r24)
	fmadds	10, 23,  1, 10
	lfs	   29,	44(%r27)
	fnmsubs	11, 23,  0, 11
	stfs	   17,	20(%r24)
	lfs	   2,	8 ( %r19)
	lfs	   3,	12 ( %r19)
#<<<<<< Break 1
	lfs	   18,	48(%r27)
	mr	 %r20,	 %r22
	fmadds	 6, 24,  2,  6
	mr	 %r28,	 %r30
	lfs	   19,	52(%r27)
	fmadds	 7, 24,  3,  7
	lfs	   20,	56(%r27)
	fmadds	 8, 26,  2,  8
	lfs	   21,	60(%r27)
	fmadds	 9, 26,  3,  9
	lfs	   22,	64(%r27)
	fmadds	10, 28,  2, 10
	lfs	   23,	68(%r27)
	fmadds	11, 28,  3, 11
	lfs	   4,	16 ( %r19)
	fmadds	 6, 25,  3,  6
	lfs	   5,	20 ( %r19)
	fnmsubs	 7, 25,  2,  7
	fmadds	 8, 27,  3,  8
	fnmsubs	 9, 27,  2,  9
	fmadds	10, 29,  3, 10
	fnmsubs	11, 29,  2, 11
#<<<<<< Break 2
	fmadds	 6, 18,  4,  6
	mr	 %r24,	 %r26
	lfs	   24,	0(%r28)
	fmadds	 7, 18,  5,  7
	lfs	   25,	4(%r28)
	fmadds	 8, 20,  4,  8
	lfs	   26,	8(%r28)
	fmadds	 9, 20,  5,  9
	lfs	   27,	12(%r28)
	fmadds	10, 22,  4, 10
	lfs	   28,	16(%r28)
	fmadds	11, 22,  5, 11
	lfs	   29,	20(%r28)
	fmadds	 6, 19,  5,  6
	lfs	   0,	0 ( %r20)
	fnmsubs	 7, 19,  4,  7
	lfs	   1,	4 ( %r20)
	fmadds	 8, 21,  5,  8
	fnmsubs	 9, 21,  4,  9
	fmadds	10, 23,  5, 10
	fnmsubs	11, 23,  4, 11
#<<<<<< Break 3
#<<<<<< Break 4
	lfs	   18,	24(%r28)
	fmuls	12, 24,  0
	stfs	   6,	0(%r23)
	fmuls	13, 24,  1
	lfs	   19,	28(%r28)
	fmuls	14, 26,  0
	stfs	   7,	4(%r23)
	fmuls	15, 26,  1
	lfs	   20,	32(%r28)
	fmuls	16, 28,  0
	stfs	   8,	8(%r23)
	fmuls	17, 28,  1
	lfs	   21,	36(%r28)
	fmadds	12, 25,  1, 12
	stfs	   9,	12(%r23)
	fnmsubs	13, 25,  0, 13
	lfs	   22,	40(%r28)
	fmadds	14, 27,  1, 14
	stfs	   10,	16(%r23)
	fnmsubs	15, 27,  0, 15
	lfs	   23,	44(%r28)
	fmadds	16, 29,  1, 16
	stfs	   11,	20(%r23)
	fnmsubs	17, 29,  0, 17
	lfs	   2,	8 ( %r20)
	lfs	   3,	12 ( %r20)
#<<<<<< Break 5
	lfs	   24,	48(%r28)
	fmadds	12, 18,  2, 12
	lfs	   25,	52(%r28)
	fmadds	13, 18,  3, 13
	lfs	   26,	56(%r28)
	fmadds	14, 20,  2, 14
	lfs	   27,	60(%r28)
	fmadds	15, 20,  3, 15
	lfs	   28,	64(%r28)
	fmadds	16, 22,  2, 16
	lfs	   29,	68(%r28)
	fmadds	17, 22,  3, 17
	lfs	   4,	16 ( %r20)
	fmadds	12, 19,  3, 12
	lfs	   5,	20 ( %r20)
	fnmsubs	13, 19,  2, 13
	fmadds	14, 21,  3, 14
	fnmsubs	15, 21,  2, 15
	fmadds	16, 23,  3, 16
	fnmsubs	17, 23,  2, 17
#<<<<<< Break 6
	fmadds	12, 24,  4, 12

	fmadds	13, 24,  5, 13

	fmadds	14, 26,  4, 14

	fmadds	15, 26,  5, 15

	fmadds	16, 28,  4, 16

	fmadds	17, 28,  5, 17

	fmadds	12, 25,  5, 12

	fnmsubs	13, 25,  4, 13
	fmadds	14, 27,  5, 14
	fnmsubs	15, 27,  4, 15
	fmadds	16, 29,  5, 16
	fnmsubs	17, 29,  4, 17
#<<<<<< Break 7
	stfs	   12,	0(%r24)
	stfs	   13,	4(%r24)
	stfs	   14,	8(%r24)
	stfs	   15,	12(%r24)
	stfs	   16,	16(%r24)
	stfs	   17,	20(%r24)
pt_cmv_dag_s_lab0:
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
