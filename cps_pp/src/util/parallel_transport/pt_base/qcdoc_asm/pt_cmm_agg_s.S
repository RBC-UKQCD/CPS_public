	.section	".text"
	.align 2
	.globl pt_cmm_agg_s
	.type	pt_cmm_agg_s,@function
pt_cmm_agg_s:
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
	li   19,	80
	li   20,	112
	li   21,	144
	li   22,	176
	li   23,	208
	li   24,	240
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
	addi %r13, %r3, -80
	lwz	%r8, 0(%r7)
	lwz	%r10, 4(%r7)
	mulli   9,8,3
	mulli   11,10,3
	add	%r8, %r9, %r4
	add	%r12, %r11, %r5
	add	%r18, %r11, %r5
	or.	 6,	6,	6
	subic.	 6,	6,	1
	lfs	   0,	8(%r7)
	lfs	   1,	12(%r7)
	lfs	   4,	0(%r8)
	lfs	   5,	4(%r8)
	lfs	   6,	8(%r8)
	lfs	   7,	12(%r8)
	lfs	   16,	16(%r8)
	lfs	   17,	20(%r8)
	bf gt,	lab0
	b  lab1
lab1:
#<<<<<< Break 0
	lfs	   2,	16(%r7)
	addi %r13, %r13, 160 
	lfs	   3,	20(%r7)
	fmuls	20,  0,  4
	lfs	   8,	24(%r8)
	fmuls	21,  0,  5
	lfs	   9,	28(%r8)
	fmuls	22,  0,  6
	lfs	   10,	32(%r8)
	fmuls	23,  0,  7
	lfs	   11,	36(%r8)
	fmuls	24,  0, 16
	stfs	   26,	48(%r18)
	fmuls	25,  0, 17
	lfs	   18,	40(%r8)
	fnmsubs	20,  1,  5, 20
	lfs	   19,	44(%r8)
	fmadds	21,  1,  4, 21
	stfs	   27,	52(%r18)
	fnmsubs	22,  1,  7, 22
	fmadds	23,  1,  6, 23
	fnmsubs	24,  1, 17, 24
	fmadds	25,  1, 16, 25
#<<<<<< Break 1
	fmadds	20,  2,  8, 20
	lwz	%r14, 0(%r13)
	fmadds	21,  2,  9, 21
	lwz	%r16, 4(%r13)
	fmadds	22,  2, 10, 22
	lfs	   0,	24(%r7)
	fmadds	23,  2, 11, 23
	stfs	   28,	56(%r18)
	fmadds	24,  2, 18, 24
	lfs	   1,	28(%r7)
	fmadds	25,  2, 19, 25
	dcbt	 %r19,	 %r7
	fnmsubs	20,  3,  9, 20
	lfs	   12,	48(%r8)
	fmadds	21,  3,  8, 21
	lfs	   13,	52(%r8)
	fnmsubs	22,  3, 11, 22
	lfs	   14,	56(%r8)
	fmadds	23,  3, 10, 23
	lfs	   15,	60(%r8)
	fnmsubs	24,  3, 19, 24
	stfs	   29,	60(%r18)
	fmadds	25,  3, 18, 25

	lfs	   16,	64(%r8)
	lfs	   17,	68(%r8)
#<<<<<< Break 2
	lfs	   2,	32(%r7)
	mulli   15,14,3
	lfs	   3,	36(%r7)
	mulli   17,16,3
	stfs	   30,	64(%r18)
	fmadds	20,  0, 12, 20

	fmadds	21,  0, 13, 21
	lfs	   18,	16(%r8)
	fmadds	22,  0, 14, 22
	lfs	   19,	20(%r8)
	fmadds	23,  0, 15, 23
	stfs	   31,	68(%r18)
	fmadds	24,  0, 16, 24
	fmadds	25,  0, 17, 25
	fnmsubs	20,  1, 13, 20
	fmadds	21,  1, 12, 21
	fnmsubs	22,  1, 15, 22
	fmadds	23,  1, 14, 23
	fnmsubs	24,  1, 17, 24
	fmadds	25,  1, 16, 25
#<<<<<< Break 3
	add	%r14, %r15, %r4
	lfs	   0,	40(%r7)
	add	%r18, %r17, %r5
	lfs	   1,	44(%r7)
	fmuls	26,  2,  4
	stfs	   20,	0(%r12)
	fmuls	27,  2,  5
	dcbt	 %r20,	 %r7
	fmuls	28,  2,  6
	lfs	   16,	40(%r8)
	fmuls	29,  2,  7
	lfs	   17,	44(%r8)
	fmuls	30,  2, 18

	fmuls	31,  2, 19

	fnmsubs	26,  3,  5, 26
	dcbt	 %r21,	 %r7
	fmadds	27,  3,  4, 27
	stfs	   21,	4(%r12)
	fnmsubs	28,  3,  7, 28

	fmadds	29,  3,  6, 29

	fnmsubs	30,  3, 19, 30
	dcbt	 %r22,	 %r7
	fmadds	31,  3, 18, 31


#<<<<<< Break 4
	stfs	   22,	8(%r12)
	fmadds	26,  0,  8, 26
	lfs	   2,	48(%r7)
	fmadds	27,  0,  9, 27
	dcbt	 %r25,	 %r14
	fmadds	28,  0, 10, 28
	lfs	   3,	52(%r7)
	fmadds	29,  0, 11, 29
	lfs	   18,	64(%r8)
	fmadds	30,  0, 16, 30
	lfs	   19,	68(%r8)
	fmadds	31,  0, 17, 31
	stfs	   23,	12(%r12)
	fnmsubs	26,  1,  9, 26
	fmadds	27,  1,  8, 27
	fnmsubs	28,  1, 11, 28
	fmadds	29,  1, 10, 29
	fnmsubs	30,  1, 17, 30
	fmadds	31,  1, 16, 31
#<<<<<< Break 5
	fmadds	26,  2, 12, 26
	stfs	   24,	16(%r12)
	fmadds	27,  2, 13, 27
	lfs	   0,	56(%r7)
	fmadds	28,  2, 14, 28

	fmadds	29,  2, 15, 29
	lfs	   1,	60(%r7)
	fmadds	30,  2, 18, 30
	lfs	   16,	16(%r8)
	fmadds	31,  2, 19, 31
	lfs	   17,	20(%r8)
	fnmsubs	26,  3, 13, 26
	stfs	   25,	20(%r12)
	fmadds	27,  3, 12, 27
	fnmsubs	28,  3, 15, 28
	fmadds	29,  3, 14, 29
	fnmsubs	30,  3, 19, 30
	fmadds	31,  3, 18, 31
#<<<<<< Break 6
	fmuls	20,  0,  4
	stfs	   26,	24(%r12)
	fmuls	21,  0,  5
	lfs	   2,	64(%r7)
	fmuls	22,  0,  6

	fmuls	23,  0,  7
	lfs	   3,	68(%r7)
	fmuls	24,  0, 16
	lfs	   18,	40(%r8)
	fmuls	25,  0, 17
	lfs	   19,	44(%r8)
	fnmsubs	20,  1,  5, 20
	stfs	   27,	28(%r12)
	fmadds	21,  1,  4, 21
	fnmsubs	22,  1,  7, 22
	fmadds	23,  1,  6, 23
	fnmsubs	24,  1, 17, 24
	fmadds	25,  1, 16, 25
#<<<<<< Break 7
	fmadds	20,  2,  8, 20
	stfs	   28,	32(%r12)
	fmadds	21,  2,  9, 21
	lfs	   0,	72(%r7)
	fmadds	22,  2, 10, 22
	dcbt	 %r26,	 %r14
	fmadds	23,  2, 11, 23
	lfs	   1,	76(%r7)
	fmadds	24,  2, 18, 24
	lfs	   16,	64(%r8)
	fmadds	25,  2, 19, 25
	lfs	   17,	68(%r8)
	fnmsubs	20,  3,  9, 20
	stfs	   29,	36(%r12)
	fmadds	21,  3,  8, 21
	fnmsubs	22,  3, 11, 22
	fmadds	23,  3, 10, 23
	fnmsubs	24,  3, 19, 24
	fmadds	25,  3, 18, 25
#<<<<<< Break 8
	fmadds	20,  0, 12, 20
	stfs	   30,	40(%r12)
	fmadds	21,  0, 13, 21
	lfs	   2,	8(%r13)
	fmadds	22,  0, 14, 22

	fmadds	23,  0, 15, 23
	lfs	   3,	12(%r13)
	fmadds	24,  0, 16, 24
	lfs	   4,	0(%r14)
	fmadds	25,  0, 17, 25
	lfs	   5,	4(%r14)
	fnmsubs	20,  1, 13, 20
	lfs	   6,	8(%r14)
	fmadds	21,  1, 12, 21
	stfs	   31,	44(%r12)
	fnmsubs	22,  1, 15, 22

	fmadds	23,  1, 14, 23
	lfs	   7,	12(%r14)
	fnmsubs	24,  1, 17, 24
	lfs	   18,	16(%r14)
	fmadds	25,  1, 16, 25
	lfs	   19,	20(%r14)
#<<<<<< Break 9
	lfs	   0,	16(%r13)
	addi %r7, %r7, 160 
	stfs	   20,	48(%r12)
	fmuls	26,  2,  4
	dcbt	 %r27,	 %r14
	fmuls	27,  2,  5
	lfs	   1,	20(%r13)
	fmuls	28,  2,  6
	lfs	   8,	24(%r14)
	fmuls	29,  2,  7
	lfs	   9,	28(%r14)
	fmuls	30,  2, 18
	lfs	   10,	32(%r14)
	fmuls	31,  2, 19

	fnmsubs	26,  3,  5, 26
	stfs	   21,	52(%r12)
	fmadds	27,  3,  4, 27
	lfs	   11,	36(%r14)
	fnmsubs	28,  3,  7, 28
	lfs	   16,	40(%r14)
	fmadds	29,  3,  6, 29
	lfs	   17,	44(%r14)
	fnmsubs	30,  3, 19, 30

	fmadds	31,  3, 18, 31
	dcbt	 %r28,	 %r14


#<<<<<< Break 10
	dcbt	 %r19,	 %r13
	fmadds	26,  0,  8, 26
	lwz	%r8, 0(%r7)
	fmadds	27,  0,  9, 27
	lwz	%r10, 4(%r7)
	fmadds	28,  0, 10, 28
	lfs	   2,	24(%r13)
	fmadds	29,  0, 11, 29
	lfs	   3,	28(%r13)
	fmadds	30,  0, 16, 30

	fmadds	31,  0, 17, 31
	stfs	   22,	56(%r12)
	fnmsubs	26,  1,  9, 26
	lfs	   12,	48(%r14)
	fmadds	27,  1,  8, 27
	lfs	   13,	52(%r14)
	fnmsubs	28,  1, 11, 28
	lfs	   14,	56(%r14)
	fmadds	29,  1, 10, 29
	lfs	   15,	60(%r14)
	fnmsubs	30,  1, 17, 30

	fmadds	31,  1, 16, 31
	lfs	   18,	64(%r14)
	stfs	   23,	60(%r12)
	lfs	   19,	68(%r14)
#<<<<<< Break 11
	lfs	   0,	32(%r13)
	mulli   9,8,3
	lfs	   1,	36(%r13)
	mulli   11,10,3
	dcbt	 %r20,	 %r13
	fmadds	26,  2, 12, 26
	lfs	   16,	16(%r14)
	fmadds	27,  2, 13, 27
	stfs	   24,	64(%r12)
	fmadds	28,  2, 14, 28
	lfs	   17,	20(%r14)
	fmadds	29,  2, 15, 29
	stfs	   25,	68(%r12)
	fmadds	30,  2, 18, 30
	fmadds	31,  2, 19, 31
	fnmsubs	26,  3, 13, 26
	fmadds	27,  3, 12, 27
	fnmsubs	28,  3, 15, 28
	fmadds	29,  3, 14, 29
	fnmsubs	30,  3, 19, 30
	fmadds	31,  3, 18, 31
#<<<<<< Break 12
	add	%r8, %r9, %r4
	lfs	   2,	40(%r13)
	add	%r12, %r11, %r5
	lfs	   3,	44(%r13)
	fmuls	20,  0,  4

	fmuls	21,  0,  5
	lfs	   18,	40(%r14)
	fmuls	22,  0,  6
	stfs	   26,	0(%r18)
	fmuls	23,  0,  7
	lfs	   19,	44(%r14)
	fmuls	24,  0, 16

	fmuls	25,  0, 17
	dcbt	 %r21,	 %r13
	fnmsubs	20,  1,  5, 20

	fmadds	21,  1,  4, 21

	fnmsubs	22,  1,  7, 22
	dcbt	 %r22,	 %r13
	fmadds	23,  1,  6, 23
	stfs	   27,	4(%r18)
	fnmsubs	24,  1, 17, 24

	fmadds	25,  1, 16, 25

#<<<<<< Break 13
	lfs	   0,	48(%r13)
	fmadds	20,  2,  8, 20
	dcbt	 %r25,	 %r8
	fmadds	21,  2,  9, 21
	lfs	   1,	52(%r13)
	fmadds	22,  2, 10, 22
	stfs	   28,	8(%r18)
	fmadds	23,  2, 11, 23
	lfs	   16,	64(%r14)
	fmadds	24,  2, 18, 24
	lfs	   17,	68(%r14)
	fmadds	25,  2, 19, 25
	stfs	   29,	12(%r18)
	fnmsubs	20,  3,  9, 20
	fmadds	21,  3,  8, 21
	fnmsubs	22,  3, 11, 22
	fmadds	23,  3, 10, 23
	fnmsubs	24,  3, 19, 24
	fmadds	25,  3, 18, 25
#<<<<<< Break 14
	fmadds	20,  0, 12, 20
	lfs	   2,	56(%r13)
	fmadds	21,  0, 13, 21

	fmadds	22,  0, 14, 22
	lfs	   3,	60(%r13)
	fmadds	23,  0, 15, 23
	stfs	   30,	16(%r18)
	fmadds	24,  0, 16, 24
	lfs	   18,	16(%r14)
	fmadds	25,  0, 17, 25
	lfs	   19,	20(%r14)
	fnmsubs	20,  1, 13, 20
	stfs	   31,	20(%r18)
	fmadds	21,  1, 12, 21
	fnmsubs	22,  1, 15, 22
	fmadds	23,  1, 14, 23
	fnmsubs	24,  1, 17, 24
	fmadds	25,  1, 16, 25
#<<<<<< Break 15
	fmuls	26,  2,  4
	lfs	   0,	64(%r13)
	fmuls	27,  2,  5

	fmuls	28,  2,  6
	lfs	   1,	68(%r13)
	fmuls	29,  2,  7
	stfs	   20,	24(%r18)
	fmuls	30,  2, 18
	lfs	   16,	40(%r14)
	fmuls	31,  2, 19
	lfs	   17,	44(%r14)
	fnmsubs	26,  3,  5, 26
	stfs	   21,	28(%r18)
	fmadds	27,  3,  4, 27
	fnmsubs	28,  3,  7, 28
	fmadds	29,  3,  6, 29
	fnmsubs	30,  3, 19, 30
	fmadds	31,  3, 18, 31
#<<<<<< Break 16
	fmadds	26,  0,  8, 26
	lfs	   2,	72(%r13)
	fmadds	27,  0,  9, 27
	dcbt	 %r26,	 %r8
	fmadds	28,  0, 10, 28
	lfs	   3,	76(%r13)
	fmadds	29,  0, 11, 29
	stfs	   22,	32(%r18)
	fmadds	30,  0, 16, 30
	lfs	   18,	64(%r14)
	fmadds	31,  0, 17, 31
	lfs	   19,	68(%r14)
	fnmsubs	26,  1,  9, 26
	stfs	   23,	36(%r18)
	fmadds	27,  1,  8, 27
	fnmsubs	28,  1, 11, 28
	fmadds	29,  1, 10, 29
	fnmsubs	30,  1, 17, 30
	fmadds	31,  1, 16, 31
#<<<<<< Break 17
	subic.	 6,	6,	1
	lfs	   0,	8(%r7)
	fmadds	26,  2, 12, 26

	fmadds	27,  2, 13, 27
	lfs	   1,	12(%r7)
	fmadds	28,  2, 14, 28
	stfs	   24,	40(%r18)
	fmadds	29,  2, 15, 29
	lfs	   4,	0(%r8)
	fmadds	30,  2, 18, 30
	lfs	   5,	4(%r8)
	fmadds	31,  2, 19, 31
	lfs	   6,	8(%r8)
	fnmsubs	26,  3, 13, 26

	fmadds	27,  3, 12, 27
	lfs	   7,	12(%r8)
	fnmsubs	28,  3, 15, 28
	lfs	   16,	16(%r8)
	fmadds	29,  3, 14, 29
	stfs	   25,	44(%r18)
	fnmsubs	30,  3, 19, 30
	lfs	   17,	20(%r8)
	fmadds	31,  3, 18, 31
#<<<<<< Break 18
	dcbt	 %r27,	 %r8


	dcbt	 %r28,	 %r8


	bf gt,	lab0
	b  lab1
lab0:
#<<<<<< Break 0
	lfs	   2,	16(%r7)
	addi %r13, %r13, 160 
	lfs	   3,	20(%r7)
	fmuls	20,  0,  4
	lfs	   8,	24(%r8)
	fmuls	21,  0,  5
	lfs	   9,	28(%r8)
	fmuls	22,  0,  6
	lfs	   10,	32(%r8)
	fmuls	23,  0,  7
	lfs	   11,	36(%r8)
	fmuls	24,  0, 16
	stfs	   26,	48(%r18)
	fmuls	25,  0, 17
	lfs	   18,	40(%r8)
	fnmsubs	20,  1,  5, 20
	lfs	   19,	44(%r8)
	fmadds	21,  1,  4, 21
	stfs	   27,	52(%r18)
	fnmsubs	22,  1,  7, 22
	fmadds	23,  1,  6, 23
	fnmsubs	24,  1, 17, 24
	fmadds	25,  1, 16, 25
#<<<<<< Break 1
	fmadds	20,  2,  8, 20
	lwz	%r14, 0(%r13)
	fmadds	21,  2,  9, 21
	lwz	%r16, 4(%r13)
	fmadds	22,  2, 10, 22
	lfs	   0,	24(%r7)
	fmadds	23,  2, 11, 23
	stfs	   28,	56(%r18)
	fmadds	24,  2, 18, 24
	lfs	   1,	28(%r7)
	fmadds	25,  2, 19, 25
	dcbt	 %r19,	 %r7
	fnmsubs	20,  3,  9, 20
	lfs	   12,	48(%r8)
	fmadds	21,  3,  8, 21
	lfs	   13,	52(%r8)
	fnmsubs	22,  3, 11, 22
	lfs	   14,	56(%r8)
	fmadds	23,  3, 10, 23
	lfs	   15,	60(%r8)
	fnmsubs	24,  3, 19, 24
	stfs	   29,	60(%r18)
	fmadds	25,  3, 18, 25

	lfs	   16,	64(%r8)
	lfs	   17,	68(%r8)
#<<<<<< Break 2
	lfs	   2,	32(%r7)
	mulli   15,14,3
	lfs	   3,	36(%r7)
	mulli   17,16,3
	stfs	   30,	64(%r18)
	fmadds	20,  0, 12, 20

	fmadds	21,  0, 13, 21
	lfs	   18,	16(%r8)
	fmadds	22,  0, 14, 22
	lfs	   19,	20(%r8)
	fmadds	23,  0, 15, 23
	stfs	   31,	68(%r18)
	fmadds	24,  0, 16, 24
	fmadds	25,  0, 17, 25
	fnmsubs	20,  1, 13, 20
	fmadds	21,  1, 12, 21
	fnmsubs	22,  1, 15, 22
	fmadds	23,  1, 14, 23
	fnmsubs	24,  1, 17, 24
	fmadds	25,  1, 16, 25
#<<<<<< Break 3
	add	%r14, %r15, %r4
	lfs	   0,	40(%r7)
	add	%r18, %r17, %r5
	lfs	   1,	44(%r7)
	fmuls	26,  2,  4
	stfs	   20,	0(%r12)
	fmuls	27,  2,  5
	dcbt	 %r20,	 %r7
	fmuls	28,  2,  6
	lfs	   16,	40(%r8)
	fmuls	29,  2,  7
	lfs	   17,	44(%r8)
	fmuls	30,  2, 18

	fmuls	31,  2, 19

	fnmsubs	26,  3,  5, 26
	dcbt	 %r21,	 %r7
	fmadds	27,  3,  4, 27
	stfs	   21,	4(%r12)
	fnmsubs	28,  3,  7, 28

	fmadds	29,  3,  6, 29

	fnmsubs	30,  3, 19, 30
	dcbt	 %r22,	 %r7
	fmadds	31,  3, 18, 31


#<<<<<< Break 4
	stfs	   22,	8(%r12)
	fmadds	26,  0,  8, 26
	lfs	   2,	48(%r7)
	fmadds	27,  0,  9, 27
	dcbt	 %r25,	 %r14
	fmadds	28,  0, 10, 28
	lfs	   3,	52(%r7)
	fmadds	29,  0, 11, 29
	lfs	   18,	64(%r8)
	fmadds	30,  0, 16, 30
	lfs	   19,	68(%r8)
	fmadds	31,  0, 17, 31
	stfs	   23,	12(%r12)
	fnmsubs	26,  1,  9, 26
	fmadds	27,  1,  8, 27
	fnmsubs	28,  1, 11, 28
	fmadds	29,  1, 10, 29
	fnmsubs	30,  1, 17, 30
	fmadds	31,  1, 16, 31
#<<<<<< Break 5
	fmadds	26,  2, 12, 26
	stfs	   24,	16(%r12)
	fmadds	27,  2, 13, 27
	lfs	   0,	56(%r7)
	fmadds	28,  2, 14, 28

	fmadds	29,  2, 15, 29
	lfs	   1,	60(%r7)
	fmadds	30,  2, 18, 30
	lfs	   16,	16(%r8)
	fmadds	31,  2, 19, 31
	lfs	   17,	20(%r8)
	fnmsubs	26,  3, 13, 26
	stfs	   25,	20(%r12)
	fmadds	27,  3, 12, 27
	fnmsubs	28,  3, 15, 28
	fmadds	29,  3, 14, 29
	fnmsubs	30,  3, 19, 30
	fmadds	31,  3, 18, 31
#<<<<<< Break 6
	fmuls	20,  0,  4
	stfs	   26,	24(%r12)
	fmuls	21,  0,  5
	lfs	   2,	64(%r7)
	fmuls	22,  0,  6

	fmuls	23,  0,  7
	lfs	   3,	68(%r7)
	fmuls	24,  0, 16
	lfs	   18,	40(%r8)
	fmuls	25,  0, 17
	lfs	   19,	44(%r8)
	fnmsubs	20,  1,  5, 20
	stfs	   27,	28(%r12)
	fmadds	21,  1,  4, 21
	fnmsubs	22,  1,  7, 22
	fmadds	23,  1,  6, 23
	fnmsubs	24,  1, 17, 24
	fmadds	25,  1, 16, 25
#<<<<<< Break 7
	fmadds	20,  2,  8, 20
	stfs	   28,	32(%r12)
	fmadds	21,  2,  9, 21
	lfs	   0,	72(%r7)
	fmadds	22,  2, 10, 22
	dcbt	 %r26,	 %r14
	fmadds	23,  2, 11, 23
	lfs	   1,	76(%r7)
	fmadds	24,  2, 18, 24
	lfs	   16,	64(%r8)
	fmadds	25,  2, 19, 25
	lfs	   17,	68(%r8)
	fnmsubs	20,  3,  9, 20
	stfs	   29,	36(%r12)
	fmadds	21,  3,  8, 21
	fnmsubs	22,  3, 11, 22
	fmadds	23,  3, 10, 23
	fnmsubs	24,  3, 19, 24
	fmadds	25,  3, 18, 25
#<<<<<< Break 8
	fmadds	20,  0, 12, 20
	stfs	   30,	40(%r12)
	fmadds	21,  0, 13, 21
	lfs	   2,	8(%r13)
	fmadds	22,  0, 14, 22

	fmadds	23,  0, 15, 23
	lfs	   3,	12(%r13)
	fmadds	24,  0, 16, 24
	lfs	   4,	0(%r14)
	fmadds	25,  0, 17, 25
	lfs	   5,	4(%r14)
	fnmsubs	20,  1, 13, 20
	lfs	   6,	8(%r14)
	fmadds	21,  1, 12, 21
	stfs	   31,	44(%r12)
	fnmsubs	22,  1, 15, 22

	fmadds	23,  1, 14, 23
	lfs	   7,	12(%r14)
	fnmsubs	24,  1, 17, 24
	lfs	   18,	16(%r14)
	fmadds	25,  1, 16, 25
	lfs	   19,	20(%r14)
#<<<<<< Break 9
	lfs	   0,	16(%r13)
	addi %r7, %r7, 160 
	stfs	   20,	48(%r12)
	fmuls	26,  2,  4
	dcbt	 %r27,	 %r14
	fmuls	27,  2,  5
	lfs	   1,	20(%r13)
	fmuls	28,  2,  6
	lfs	   8,	24(%r14)
	fmuls	29,  2,  7
	lfs	   9,	28(%r14)
	fmuls	30,  2, 18
	lfs	   10,	32(%r14)
	fmuls	31,  2, 19

	fnmsubs	26,  3,  5, 26
	stfs	   21,	52(%r12)
	fmadds	27,  3,  4, 27
	lfs	   11,	36(%r14)
	fnmsubs	28,  3,  7, 28
	lfs	   16,	40(%r14)
	fmadds	29,  3,  6, 29
	lfs	   17,	44(%r14)
	fnmsubs	30,  3, 19, 30

	fmadds	31,  3, 18, 31
	dcbt	 %r28,	 %r14


#<<<<<< Break 10
	dcbt	 %r19,	 %r13
	fmadds	26,  0,  8, 26
	lwz	%r8, 0(%r7)
	fmadds	27,  0,  9, 27
	lwz	%r10, 4(%r7)
	fmadds	28,  0, 10, 28
	lfs	   2,	24(%r13)
	fmadds	29,  0, 11, 29
	lfs	   3,	28(%r13)
	fmadds	30,  0, 16, 30

	fmadds	31,  0, 17, 31
	stfs	   22,	56(%r12)
	fnmsubs	26,  1,  9, 26
	lfs	   12,	48(%r14)
	fmadds	27,  1,  8, 27
	lfs	   13,	52(%r14)
	fnmsubs	28,  1, 11, 28
	lfs	   14,	56(%r14)
	fmadds	29,  1, 10, 29
	lfs	   15,	60(%r14)
	fnmsubs	30,  1, 17, 30

	fmadds	31,  1, 16, 31
	lfs	   18,	64(%r14)
	stfs	   23,	60(%r12)
	lfs	   19,	68(%r14)
#<<<<<< Break 11
	lfs	   0,	32(%r13)
	mulli   9,8,3
	lfs	   1,	36(%r13)
	mulli   11,10,3
	dcbt	 %r20,	 %r13
	fmadds	26,  2, 12, 26
	lfs	   16,	16(%r14)
	fmadds	27,  2, 13, 27
	stfs	   24,	64(%r12)
	fmadds	28,  2, 14, 28
	lfs	   17,	20(%r14)
	fmadds	29,  2, 15, 29
	stfs	   25,	68(%r12)
	fmadds	30,  2, 18, 30
	fmadds	31,  2, 19, 31
	fnmsubs	26,  3, 13, 26
	fmadds	27,  3, 12, 27
	fnmsubs	28,  3, 15, 28
	fmadds	29,  3, 14, 29
	fnmsubs	30,  3, 19, 30
	fmadds	31,  3, 18, 31
#<<<<<< Break 12
	add	%r8, %r9, %r4
	lfs	   2,	40(%r13)
	add	%r12, %r11, %r5
	lfs	   3,	44(%r13)
	fmuls	20,  0,  4

	fmuls	21,  0,  5
	lfs	   18,	40(%r14)
	fmuls	22,  0,  6
	stfs	   26,	0(%r18)
	fmuls	23,  0,  7
	lfs	   19,	44(%r14)
	fmuls	24,  0, 16

	fmuls	25,  0, 17
	dcbt	 %r21,	 %r13
	fnmsubs	20,  1,  5, 20

	fmadds	21,  1,  4, 21

	fnmsubs	22,  1,  7, 22
	dcbt	 %r22,	 %r13
	fmadds	23,  1,  6, 23
	stfs	   27,	4(%r18)
	fnmsubs	24,  1, 17, 24

	fmadds	25,  1, 16, 25

#<<<<<< Break 13
	lfs	   0,	48(%r13)
	fmadds	20,  2,  8, 20
	dcbt	 %r25,	 %r8
	fmadds	21,  2,  9, 21
	lfs	   1,	52(%r13)
	fmadds	22,  2, 10, 22
	stfs	   28,	8(%r18)
	fmadds	23,  2, 11, 23
	lfs	   16,	64(%r14)
	fmadds	24,  2, 18, 24
	lfs	   17,	68(%r14)
	fmadds	25,  2, 19, 25
	stfs	   29,	12(%r18)
	fnmsubs	20,  3,  9, 20
	fmadds	21,  3,  8, 21
	fnmsubs	22,  3, 11, 22
	fmadds	23,  3, 10, 23
	fnmsubs	24,  3, 19, 24
	fmadds	25,  3, 18, 25
#<<<<<< Break 14
	fmadds	20,  0, 12, 20
	lfs	   2,	56(%r13)
	fmadds	21,  0, 13, 21

	fmadds	22,  0, 14, 22
	lfs	   3,	60(%r13)
	fmadds	23,  0, 15, 23
	stfs	   30,	16(%r18)
	fmadds	24,  0, 16, 24
	lfs	   18,	16(%r14)
	fmadds	25,  0, 17, 25
	lfs	   19,	20(%r14)
	fnmsubs	20,  1, 13, 20
	stfs	   31,	20(%r18)
	fmadds	21,  1, 12, 21
	fnmsubs	22,  1, 15, 22
	fmadds	23,  1, 14, 23
	fnmsubs	24,  1, 17, 24
	fmadds	25,  1, 16, 25
#<<<<<< Break 15
	fmuls	26,  2,  4
	lfs	   0,	64(%r13)
	fmuls	27,  2,  5

	fmuls	28,  2,  6
	lfs	   1,	68(%r13)
	fmuls	29,  2,  7
	stfs	   20,	24(%r18)
	fmuls	30,  2, 18
	lfs	   16,	40(%r14)
	fmuls	31,  2, 19
	lfs	   17,	44(%r14)
	fnmsubs	26,  3,  5, 26
	stfs	   21,	28(%r18)
	fmadds	27,  3,  4, 27
	fnmsubs	28,  3,  7, 28
	fmadds	29,  3,  6, 29
	fnmsubs	30,  3, 19, 30
	fmadds	31,  3, 18, 31
#<<<<<< Break 16
	fmadds	26,  0,  8, 26
	lfs	   2,	72(%r13)
	fmadds	27,  0,  9, 27
	dcbt	 %r26,	 %r8
	fmadds	28,  0, 10, 28
	lfs	   3,	76(%r13)
	fmadds	29,  0, 11, 29
	stfs	   22,	32(%r18)
	fmadds	30,  0, 16, 30
	lfs	   18,	64(%r14)
	fmadds	31,  0, 17, 31
	lfs	   19,	68(%r14)
	fnmsubs	26,  1,  9, 26
	stfs	   23,	36(%r18)
	fmadds	27,  1,  8, 27
	fnmsubs	28,  1, 11, 28
	fmadds	29,  1, 10, 29
	fnmsubs	30,  1, 17, 30
	fmadds	31,  1, 16, 31
#<<<<<< Break 17
	subic.	 6,	6,	1
	stfs	   24,	40(%r18)
	fmadds	26,  2, 12, 26
	stfs	   25,	44(%r18)
	fmadds	27,  2, 13, 27
	fmadds	28,  2, 14, 28
	fmadds	29,  2, 15, 29
	fmadds	30,  2, 18, 30
	fmadds	31,  2, 19, 31
	fnmsubs	26,  3, 13, 26
	fmadds	27,  3, 12, 27
	fnmsubs	28,  3, 15, 28
	fmadds	29,  3, 14, 29
	fnmsubs	30,  3, 19, 30
	fmadds	31,  3, 18, 31
#<<<<<< Break 18


	dcbt	 %r27,	 %r8


	dcbt	 %r28,	 %r8


	stfs	   26,	48(%r18)
	stfs	   27,	52(%r18)
	stfs	   28,	56(%r18)
	stfs	   29,	60(%r18)
	stfs	   30,	64(%r18)
	stfs	   31,	68(%r18)
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
