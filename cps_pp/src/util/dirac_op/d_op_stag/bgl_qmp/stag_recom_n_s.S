	.section	".text"
	.align 2
	.globl stag_recom_n_s
	.type	stag_recom_n_s,@function
stag_recom_n_s:
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
	li   27,	352
	li   28,	384
	li   29,	416
	li   30,	448
	li   31,	480
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
	lwz	%r8, 0(%r5)
	lwz	%r9, 4(%r5)
	lwz	%r10, 8(%r5)
	lwz	%r11, 12(%r5)
	addi %r8, %r8, -24
	addi %r9, %r9, -24
	addi %r10, %r10, -24
	addi %r11, %r11, -24
	lfs	   0,	0(%r4)
	lfs	   1,	4(%r4)
	lfs	   2,	8(%r4)
	lfs	   3,	12(%r4)
	lfs	   4,	16(%r4)
	lfs	   5,	20(%r4)
	lfs	   6,	32(%r4)
	lfs	   7,	36(%r4)
	lfs	   8,	40(%r4)
	lfs	   9,	44(%r4)
	lfs	   10,	48(%r4)
	lfs	   11,	52(%r4)
	lfs	   12,	64(%r4)
	lfs	   13,	68(%r4)
	lfs	   14,	72(%r4)
	lfs	   15,	76(%r4)
	lfs	   16,	80(%r4)
	lfs	   17,	84(%r4)
	addi %r7, %r4, 256
	dcbt	 %r16,	 %r4
	dcbt	 %r17,	 %r4
	dcbt	 %r18,	 %r4
	dcbt	 %r19,	 %r4
	or.	 3,	3,	3
	dcbt	 %r20,	 %r4
	dcbt	 %r21,	 %r4
	dcbt	 %r22,	 %r4
	dcbt	 %r23,	 %r4
	bf gt,	stag_recom_n_s_lab0
	dcbt	 %r24,	 %r4
	dcbt	 %r25,	 %r4
	dcbt	 %r26,	 %r4
	dcbt	 %r27,	 %r4
	subic.	 3,	3,	1
	dcbt	 %r28,	 %r4
	dcbt	 %r29,	 %r4
	dcbt	 %r30,	 %r4
	dcbt	 %r31,	 %r4
	bf gt,	stag_recom_n_s_lab2
	b  stag_recom_n_s_lab1
stag_recom_n_s_lab1:
	lfs	   18,	96(%r4)
	addi %r7, %r7, 256 
	fsubs	   0,	0,	6
	lfs	   19,	100(%r4)
	lfs	   20,	104(%r4)
	fsubs	   1,	1,	7
	lfs	   21,	108(%r4)
	lfs	   22,	112(%r4)
	fsubs	   2,	2,	8
	lfs	   23,	116(%r4)
	fsubs	   3,	3,	9
	fsubs	   4,	4,	10
	fsubs	   5,	5,	11
#<<<<<< Break 1
	fsubs	   12,	12,	18
	lfs	   6,	128(%r4)
	lfs	   7,	132(%r4)
	fsubs	   13,	13,	19
	lfs	   8,	136(%r4)
	lfs	   9,	140(%r4)
	fsubs	   14,	14,	20
	lfs	   10,	144(%r4)
	lfs	   11,	148(%r4)
	fsubs	   15,	15,	21
	lfs	   24,	160(%r4)
	lfs	   25,	164(%r4)
	fsubs	   16,	16,	22
	lfs	   26,	168(%r4)
	lfs	   27,	172(%r4)
	fsubs	   17,	17,	23
	lfs	   28,	176(%r4)
	lfs	   29,	180(%r4)
#<<<<<< Break 2
	lfs	   18,	192(%r4)
	addi %r8, %r8, 24 
	fsubs	   6,	6,	24
	lfs	   19,	196(%r4)
	addi %r9, %r9, 24 
	fsubs	   7,	7,	25
	lfs	   20,	200(%r4)
	addi %r10, %r10, 24 
	fsubs	   8,	8,	26
	lfs	   21,	204(%r4)
	addi %r11, %r11, 24 
	fsubs	   9,	9,	27
	lfs	   22,	208(%r4)
	lfs	   23,	212(%r4)
	fsubs	   10,	10,	28
	lfs	   30,	224(%r4)
	lfs	   31,	228(%r4)
	fsubs	   11,	11,	29
#<<<<<< Break 3
	lfs	   24,	232(%r4)
	fsubs	   18,	18,	30
	lfs	   25,	236(%r4)
	lfs	   26,	240(%r4)
	fsubs	   19,	19,	31
	lfs	   27,	244(%r4)
#<<<<<< Break 4
	fsubs	   20,	20,	24
	subic.	 3,	3,	1
	dcbt	 %r7,	 %r16
	fsubs	   21,	21,	25


	fsubs	   22,	22,	26
	fsubs	   23,	23,	27
#<<<<<< Break 5
	dcbt	 %r7,	 %r17
	fadds	   24,	0,	12


	fadds	   25,	1,	13
	dcbt	 %r7,	 %r18

	fadds	   26,	2,	14

	dcbt	 %r7,	 %r19
	fadds	   27,	3,	15


	fadds	   28,	4,	16
	dcbt	 %r7,	 %r20

	fadds	   29,	5,	17

	dcbt	 %r7,	 %r21
	fadds	   24,	24,	6


	fadds	   25,	25,	7
	fadds	   26,	26,	8
	fadds	   27,	27,	9
	fadds	   28,	28,	10
	fadds	   29,	29,	11
	fadds	   24,	24,	18
	fadds	   25,	25,	19
	fadds	   26,	26,	20
	fadds	   27,	27,	21
	fadds	   28,	28,	22
	fadds	   29,	29,	23
#<<<<<< Break 6
	dcbt	 %r7,	 %r22
	fsubs	   0,	24,	0
	addi %r4, %r4, 256 

	fsubs	   1,	25,	1

	dcbt	 %r7,	 %r23
	fsubs	   2,	26,	2


	fsubs	   3,	27,	3
	fsubs	   4,	28,	4
	fsubs	   5,	29,	5
#<<<<<< Break 7
	stfs	   0,	0(%r8)
	fsubs	   6,	24,	6
	stfs	   1,	4(%r8)
	stfs	   2,	8(%r8)
	fsubs	   7,	25,	7
	stfs	   3,	12(%r8)
	stfs	   4,	16(%r8)
	fsubs	   8,	26,	8
	stfs	   5,	20(%r8)
	stfs	   5,	24(%r8)
	fsubs	   9,	27,	9
	stfs	   5,	28(%r8)
	fsubs	   10,	28,	10
	fsubs	   11,	29,	11
#<<<<<< Break 8
	fsubs	   12,	24,	12
	stfs	   6,	0(%r10)
	lfs	   0,	0(%r4)
	fsubs	   13,	25,	13
	stfs	   7,	4(%r10)
	lfs	   1,	4(%r4)
	fsubs	   14,	26,	14
	stfs	   8,	8(%r10)
	lfs	   2,	8(%r4)
	fsubs	   15,	27,	15
	stfs	   9,	12(%r10)
	lfs	   3,	12(%r4)
	fsubs	   16,	28,	16
	stfs	   10,	16(%r10)
	lfs	   4,	16(%r4)
	fsubs	   17,	29,	17
	stfs	   11,	20(%r10)
	lfs	   5,	20(%r4)
	stfs	   11,	24(%r10)
	stfs	   11,	28(%r10)
#<<<<<< Break 9
	stfs	   12,	0(%r9)
	lfs	   6,	32(%r4)
	fsubs	   24,	24,	18
	stfs	   13,	4(%r9)
	lfs	   7,	36(%r4)
	fsubs	   25,	25,	19
	stfs	   14,	8(%r9)
	lfs	   8,	40(%r4)
	fsubs	   26,	26,	20
	stfs	   15,	12(%r9)
	lfs	   9,	44(%r4)
	fsubs	   27,	27,	21
	stfs	   16,	16(%r9)
	lfs	   10,	48(%r4)
	fsubs	   28,	28,	22
	stfs	   17,	20(%r9)
	lfs	   11,	52(%r4)
	fsubs	   29,	29,	23
	stfs	   17,	24(%r9)
	stfs	   17,	28(%r9)
#<<<<<< Break 10
	stfs	   24,	0(%r11)
	lfs	   12,	64(%r4)
	stfs	   25,	4(%r11)
	lfs	   13,	68(%r4)
	stfs	   26,	8(%r11)
	lfs	   14,	72(%r4)
	stfs	   27,	12(%r11)
	lfs	   15,	76(%r4)
	stfs	   28,	16(%r11)
	lfs	   16,	80(%r4)
	stfs	   29,	20(%r11)
	lfs	   17,	84(%r4)
	stfs	   29,	24(%r11)
	stfs	   29,	28(%r11)
	bf gt,	stag_recom_n_s_lab2
	b  stag_recom_n_s_lab1
stag_recom_n_s_lab2:
	lfs	   18,	96(%r4)
	lfs	   19,	100(%r4)
	fsubs	   0,	0,	6
	lfs	   20,	104(%r4)
	lfs	   21,	108(%r4)
	fsubs	   1,	1,	7
	lfs	   22,	112(%r4)
	lfs	   23,	116(%r4)
	fsubs	   2,	2,	8
	fsubs	   3,	3,	9
	fsubs	   4,	4,	10
	fsubs	   5,	5,	11
#<<<<<< Break 1
	lfs	   6,	128(%r4)
	fsubs	   12,	12,	18
	lfs	   7,	132(%r4)
	lfs	   8,	136(%r4)
	fsubs	   13,	13,	19
	lfs	   9,	140(%r4)
	lfs	   10,	144(%r4)
	fsubs	   14,	14,	20
	lfs	   11,	148(%r4)
	lfs	   24,	160(%r4)
	fsubs	   15,	15,	21
	lfs	   25,	164(%r4)
	lfs	   26,	168(%r4)
	fsubs	   16,	16,	22
	lfs	   27,	172(%r4)
	lfs	   28,	176(%r4)
	fsubs	   17,	17,	23
	lfs	   29,	180(%r4)
#<<<<<< Break 2
	fsubs	   6,	6,	24
	lfs	   18,	192(%r4)
	addi %r8, %r8, 24 
	fsubs	   7,	7,	25
	lfs	   19,	196(%r4)
	addi %r9, %r9, 24 
	fsubs	   8,	8,	26
	lfs	   20,	200(%r4)
	addi %r10, %r10, 24 
	fsubs	   9,	9,	27
	lfs	   21,	204(%r4)
	addi %r11, %r11, 24 
	fsubs	   10,	10,	28
	lfs	   22,	208(%r4)
	lfs	   23,	212(%r4)
	fsubs	   11,	11,	29
	lfs	   30,	224(%r4)
	lfs	   31,	228(%r4)
#<<<<<< Break 3
	lfs	   24,	232(%r4)
	lfs	   25,	236(%r4)
	fsubs	   18,	18,	30
	lfs	   26,	240(%r4)
	lfs	   27,	244(%r4)
	fsubs	   19,	19,	31
#<<<<<< Break 4
	fsubs	   20,	20,	24
	fsubs	   21,	21,	25
	fsubs	   22,	22,	26
	fsubs	   23,	23,	27
#<<<<<< Break 5
	fadds	   24,	0,	12
	fadds	   25,	1,	13
	fadds	   26,	2,	14
	fadds	   27,	3,	15
	fadds	   28,	4,	16
	fadds	   29,	5,	17
	fadds	   24,	24,	6
	fadds	   25,	25,	7
	fadds	   26,	26,	8
	fadds	   27,	27,	9
	fadds	   28,	28,	10
	fadds	   29,	29,	11
	fadds	   24,	24,	18
	fadds	   25,	25,	19
	fadds	   26,	26,	20
	fadds	   27,	27,	21
	fadds	   28,	28,	22
	fadds	   29,	29,	23
#<<<<<< Break 6
	fsubs	   0,	24,	0
	fsubs	   1,	25,	1
	fsubs	   2,	26,	2
	fsubs	   3,	27,	3
	fsubs	   4,	28,	4
	fsubs	   5,	29,	5
#<<<<<< Break 7
	fsubs	   6,	24,	6
	stfs	   0,	0(%r8)
	stfs	   1,	4(%r8)
	fsubs	   7,	25,	7
	stfs	   2,	8(%r8)
	stfs	   3,	12(%r8)
	fsubs	   8,	26,	8
	stfs	   4,	16(%r8)
	stfs	   5,	20(%r8)
	fsubs	   9,	27,	9
	fsubs	   10,	28,	10
	fsubs	   11,	29,	11
#<<<<<< Break 8
	stfs	   6,	0(%r10)
	fsubs	   12,	24,	12
	stfs	   7,	4(%r10)
	stfs	   8,	8(%r10)
	fsubs	   13,	25,	13
	stfs	   9,	12(%r10)
	stfs	   10,	16(%r10)
	fsubs	   14,	26,	14
	stfs	   11,	20(%r10)
	fsubs	   15,	27,	15
	fsubs	   16,	28,	16
	fsubs	   17,	29,	17
#<<<<<< Break 9
	fsubs	   24,	24,	18
	stfs	   12,	0(%r9)
	stfs	   13,	4(%r9)
	fsubs	   25,	25,	19
	stfs	   14,	8(%r9)
	stfs	   15,	12(%r9)
	fsubs	   26,	26,	20
	stfs	   16,	16(%r9)
	stfs	   17,	20(%r9)
	fsubs	   27,	27,	21
	fsubs	   28,	28,	22
	fsubs	   29,	29,	23
#<<<<<< Break 10
	stfs	   24,	0(%r11)
	stfs	   25,	4(%r11)
	stfs	   26,	8(%r11)
	stfs	   27,	12(%r11)
	stfs	   28,	16(%r11)
	stfs	   29,	20(%r11)
stag_recom_n_s_lab0:
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
