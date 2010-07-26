	.section	".text"
	.align 2
	.globl stag_recom
	.type	stag_recom,@function
stag_recom:
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
	addi %r8, %r8, -48
	addi %r9, %r9, -48
	addi %r10, %r10, -48
	addi %r11, %r11, -48
	lfd	   0,	0(%r4)
	lfd	   1,	8(%r4)
	lfd	   2,	16(%r4)
	lfd	   3,	24(%r4)
	lfd	   4,	32(%r4)
	lfd	   5,	40(%r4)
	lfd	   6,	64(%r4)
	lfd	   7,	72(%r4)
	lfd	   8,	80(%r4)
	lfd	   9,	88(%r4)
	lfd	   10,	96(%r4)
	lfd	   11,	104(%r4)
	lfd	   12,	128(%r4)
	lfd	   13,	136(%r4)
	lfd	   14,	144(%r4)
	lfd	   15,	152(%r4)
	lfd	   16,	160(%r4)
	lfd	   17,	168(%r4)
	addi %r7, %r4, 64
	dcbt	 %r16,	 %r4
	dcbt	 %r17,	 %r4
	dcbt	 %r18,	 %r4
	dcbt	 %r19,	 %r4
	or.	 3,	3,	3
	dcbt	 %r20,	 %r4
	dcbt	 %r21,	 %r4
	dcbt	 %r22,	 %r4
	dcbt	 %r23,	 %r4
	bf gt,	stag_recom_lab0
	dcbt	 %r24,	 %r4
	dcbt	 %r25,	 %r4
	dcbt	 %r26,	 %r4
	dcbt	 %r27,	 %r4
	subic.	 3,	3,	1
	dcbt	 %r28,	 %r4
	dcbt	 %r29,	 %r4
	dcbt	 %r30,	 %r4
	dcbt	 %r31,	 %r4
	bf gt,	stag_recom_lab2
	b  stag_recom_lab1
stag_recom_lab1:
	lfd	   18,	192(%r4)
	addi %r7, %r7, 512 
	fadd	   0,	0,	6
	lfd	   19,	200(%r4)
	lfd	   20,	208(%r4)
	fadd	   1,	1,	7
	lfd	   21,	216(%r4)
	lfd	   22,	224(%r4)
	fadd	   2,	2,	8
	lfd	   23,	232(%r4)
	fadd	   3,	3,	9
	fadd	   4,	4,	10
	fadd	   5,	5,	11
#<<<<<< Break 1
	fadd	   12,	12,	18
	lfd	   6,	256(%r4)
	lfd	   7,	264(%r4)
	fadd	   13,	13,	19
	lfd	   8,	272(%r4)
	lfd	   9,	280(%r4)
	fadd	   14,	14,	20
	lfd	   10,	288(%r4)
	lfd	   11,	296(%r4)
	fadd	   15,	15,	21
	lfd	   24,	320(%r4)
	lfd	   25,	328(%r4)
	fadd	   16,	16,	22
	lfd	   26,	336(%r4)
	lfd	   27,	344(%r4)
	fadd	   17,	17,	23
	lfd	   28,	352(%r4)
	lfd	   29,	360(%r4)
#<<<<<< Break 2
	lfd	   18,	384(%r4)
	addi %r8, %r8, 48 
	fadd	   6,	6,	24
	lfd	   19,	392(%r4)
	addi %r9, %r9, 48 
	fadd	   7,	7,	25
	lfd	   20,	400(%r4)
	addi %r10, %r10, 48 
	fadd	   8,	8,	26
	lfd	   21,	408(%r4)
	addi %r11, %r11, 48 
	fadd	   9,	9,	27
	lfd	   22,	416(%r4)
	lfd	   23,	424(%r4)
	fadd	   10,	10,	28
	lfd	   30,	448(%r4)
	lfd	   31,	456(%r4)
	fadd	   11,	11,	29
#<<<<<< Break 3
	lfd	   24,	464(%r4)
	fadd	   18,	18,	30
	lfd	   25,	472(%r4)
	lfd	   26,	480(%r4)
	fadd	   19,	19,	31
	lfd	   27,	488(%r4)
#<<<<<< Break 4
	fadd	   20,	20,	24
	subic.	 3,	3,	1
	dcbt	 %r7,	 %r16
	fadd	   21,	21,	25


	fadd	   22,	22,	26
	dcbt	 %r7,	 %r17

	fadd	   23,	23,	27

#<<<<<< Break 5
	fadd	   24,	0,	12
	dcbt	 %r7,	 %r18

	fadd	   25,	1,	13

	dcbt	 %r7,	 %r19
	fadd	   26,	2,	14


	fadd	   27,	3,	15
	dcbt	 %r7,	 %r20

	fadd	   28,	4,	16

	dcbt	 %r7,	 %r21
	fadd	   29,	5,	17


	fadd	   24,	24,	6
	dcbt	 %r7,	 %r22

	fadd	   25,	25,	7

	dcbt	 %r7,	 %r23
	fadd	   26,	26,	8


	fadd	   27,	27,	9
	dcbt	 %r7,	 %r24

	fadd	   28,	28,	10

	dcbt	 %r7,	 %r25
	fadd	   29,	29,	11


	fadd	   24,	24,	18
	dcbt	 %r7,	 %r26

	fadd	   25,	25,	19

	dcbt	 %r7,	 %r27
	fadd	   26,	26,	20


	fadd	   27,	27,	21
	dcbt	 %r7,	 %r28

	fadd	   28,	28,	22

	fadd	   29,	29,	23
#<<<<<< Break 6
	fsub	   0,	24,	0
	dcbt	 %r7,	 %r29
	addi %r4, %r4, 512 
	fsub	   1,	25,	1


	fsub	   2,	26,	2
	dcbt	 %r7,	 %r30

	fsub	   3,	27,	3

	dcbt	 %r7,	 %r31
	fsub	   4,	28,	4


	fsub	   5,	29,	5
#<<<<<< Break 7
	stfd	   0,	0(%r8)
	fsub	   6,	24,	6
	stfd	   1,	8(%r8)
	stfd	   2,	16(%r8)
	fsub	   7,	25,	7
	stfd	   3,	24(%r8)
	stfd	   4,	32(%r8)
	fsub	   8,	26,	8
	stfd	   5,	40(%r8)
	stfd	   5,	48(%r8)
	fsub	   9,	27,	9
	stfd	   5,	56(%r8)
	fsub	   10,	28,	10
	fsub	   11,	29,	11
#<<<<<< Break 8
	fsub	   12,	24,	12
	stfd	   6,	0(%r10)
	lfd	   0,	0(%r4)
	fsub	   13,	25,	13
	stfd	   7,	8(%r10)
	lfd	   1,	8(%r4)
	fsub	   14,	26,	14
	stfd	   8,	16(%r10)
	lfd	   2,	16(%r4)
	fsub	   15,	27,	15
	stfd	   9,	24(%r10)
	lfd	   3,	24(%r4)
	fsub	   16,	28,	16
	stfd	   10,	32(%r10)
	lfd	   4,	32(%r4)
	fsub	   17,	29,	17
	stfd	   11,	40(%r10)
	lfd	   5,	40(%r4)
	stfd	   11,	48(%r10)
	stfd	   11,	56(%r10)
#<<<<<< Break 9
	stfd	   12,	0(%r9)
	lfd	   6,	64(%r4)
	fsub	   24,	24,	18
	stfd	   13,	8(%r9)
	lfd	   7,	72(%r4)
	fsub	   25,	25,	19
	stfd	   14,	16(%r9)
	lfd	   8,	80(%r4)
	fsub	   26,	26,	20
	stfd	   15,	24(%r9)
	lfd	   9,	88(%r4)
	fsub	   27,	27,	21
	stfd	   16,	32(%r9)
	lfd	   10,	96(%r4)
	fsub	   28,	28,	22
	stfd	   17,	40(%r9)
	lfd	   11,	104(%r4)
	fsub	   29,	29,	23
	stfd	   17,	48(%r9)
	stfd	   17,	56(%r9)
#<<<<<< Break 10
	stfd	   24,	0(%r11)
	lfd	   12,	128(%r4)
	stfd	   25,	8(%r11)
	lfd	   13,	136(%r4)
	stfd	   26,	16(%r11)
	lfd	   14,	144(%r4)
	stfd	   27,	24(%r11)
	lfd	   15,	152(%r4)
	stfd	   28,	32(%r11)
	lfd	   16,	160(%r4)
	stfd	   29,	40(%r11)
	lfd	   17,	168(%r4)
	stfd	   29,	48(%r11)
	stfd	   29,	56(%r11)
	bf gt,	stag_recom_lab2
	b  stag_recom_lab1
stag_recom_lab2:
	lfd	   18,	192(%r4)
	lfd	   19,	200(%r4)
	fadd	   0,	0,	6
	lfd	   20,	208(%r4)
	lfd	   21,	216(%r4)
	fadd	   1,	1,	7
	lfd	   22,	224(%r4)
	lfd	   23,	232(%r4)
	fadd	   2,	2,	8
	fadd	   3,	3,	9
	fadd	   4,	4,	10
	fadd	   5,	5,	11
#<<<<<< Break 1
	lfd	   6,	256(%r4)
	fadd	   12,	12,	18
	lfd	   7,	264(%r4)
	lfd	   8,	272(%r4)
	fadd	   13,	13,	19
	lfd	   9,	280(%r4)
	lfd	   10,	288(%r4)
	fadd	   14,	14,	20
	lfd	   11,	296(%r4)
	lfd	   24,	320(%r4)
	fadd	   15,	15,	21
	lfd	   25,	328(%r4)
	lfd	   26,	336(%r4)
	fadd	   16,	16,	22
	lfd	   27,	344(%r4)
	lfd	   28,	352(%r4)
	fadd	   17,	17,	23
	lfd	   29,	360(%r4)
#<<<<<< Break 2
	fadd	   6,	6,	24
	lfd	   18,	384(%r4)
	addi %r8, %r8, 48 
	fadd	   7,	7,	25
	lfd	   19,	392(%r4)
	addi %r9, %r9, 48 
	fadd	   8,	8,	26
	lfd	   20,	400(%r4)
	addi %r10, %r10, 48 
	fadd	   9,	9,	27
	lfd	   21,	408(%r4)
	addi %r11, %r11, 48 
	fadd	   10,	10,	28
	lfd	   22,	416(%r4)
	lfd	   23,	424(%r4)
	fadd	   11,	11,	29
	lfd	   30,	448(%r4)
	lfd	   31,	456(%r4)
#<<<<<< Break 3
	lfd	   24,	464(%r4)
	lfd	   25,	472(%r4)
	fadd	   18,	18,	30
	lfd	   26,	480(%r4)
	lfd	   27,	488(%r4)
	fadd	   19,	19,	31
#<<<<<< Break 4
	fadd	   20,	20,	24
	fadd	   21,	21,	25
	fadd	   22,	22,	26
	fadd	   23,	23,	27
#<<<<<< Break 5
	fadd	   24,	0,	12
	fadd	   25,	1,	13
	fadd	   26,	2,	14
	fadd	   27,	3,	15
	fadd	   28,	4,	16
	fadd	   29,	5,	17
	fadd	   24,	24,	6
	fadd	   25,	25,	7
	fadd	   26,	26,	8
	fadd	   27,	27,	9
	fadd	   28,	28,	10
	fadd	   29,	29,	11
	fadd	   24,	24,	18
	fadd	   25,	25,	19
	fadd	   26,	26,	20
	fadd	   27,	27,	21
	fadd	   28,	28,	22
	fadd	   29,	29,	23
#<<<<<< Break 6
	fsub	   0,	24,	0
	fsub	   1,	25,	1
	fsub	   2,	26,	2
	fsub	   3,	27,	3
	fsub	   4,	28,	4
	fsub	   5,	29,	5
#<<<<<< Break 7
	fsub	   6,	24,	6
	stfd	   0,	0(%r8)
	stfd	   1,	8(%r8)
	fsub	   7,	25,	7
	stfd	   2,	16(%r8)
	stfd	   3,	24(%r8)
	fsub	   8,	26,	8
	stfd	   4,	32(%r8)
	stfd	   5,	40(%r8)
	fsub	   9,	27,	9
	fsub	   10,	28,	10
	fsub	   11,	29,	11
#<<<<<< Break 8
	stfd	   6,	0(%r10)
	fsub	   12,	24,	12
	stfd	   7,	8(%r10)
	stfd	   8,	16(%r10)
	fsub	   13,	25,	13
	stfd	   9,	24(%r10)
	stfd	   10,	32(%r10)
	fsub	   14,	26,	14
	stfd	   11,	40(%r10)
	fsub	   15,	27,	15
	fsub	   16,	28,	16
	fsub	   17,	29,	17
#<<<<<< Break 9
	fsub	   24,	24,	18
	stfd	   12,	0(%r9)
	stfd	   13,	8(%r9)
	fsub	   25,	25,	19
	stfd	   14,	16(%r9)
	stfd	   15,	24(%r9)
	fsub	   26,	26,	20
	stfd	   16,	32(%r9)
	stfd	   17,	40(%r9)
	fsub	   27,	27,	21
	fsub	   28,	28,	22
	fsub	   29,	29,	23
#<<<<<< Break 10
	stfd	   24,	0(%r11)
	stfd	   25,	8(%r11)
	stfd	   26,	16(%r11)
	stfd	   27,	24(%r11)
	stfd	   28,	32(%r11)
	stfd	   29,	40(%r11)
stag_recom_lab0:
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
