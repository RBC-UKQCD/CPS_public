	.section	".text"
	.align 2
	.globl asq_dsum
	.type	asq_dsum,@function
asq_dsum:
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
	addi %r6, %r4, -1024
	addi %r7, %r4, 384
	addi %r8, %r4, -128
	dcbt	 %r16,	 %r4
	dcbt	 %r17,	 %r4
	dcbt	 %r18,	 %r4
	dcbt	 %r19,	 %r4
	or.	 3,	3,	3
	addi %r5, %r5, -48 
	bf gt,	lab0
	lfd	   6,	0(%r4)
	lfd	   7,	8(%r4)
	lfd	   8,	16(%r4)
	lfd	   9,	24(%r4)
	lfd	   10,	32(%r4)
	lfd	   11,	40(%r4)
	lfd	   12,	64(%r4)
	lfd	   13,	72(%r4)
	lfd	   14,	80(%r4)
	lfd	   15,	88(%r4)
	lfd	   16,	96(%r4)
	lfd	   17,	104(%r4)
	subic.	 3,	3,	1
	bf gt,	lab2
	b  lab1
lab1:
	lfd	   18,	128(%r4)
	dcbt	 %r7,	 %r16
	fadd	   0,	6,	12
	lfd	   19,	136(%r4)

	fadd	   1,	7,	13
	lfd	   20,	144(%r4)

	fadd	   2,	8,	14
	lfd	   21,	152(%r4)
	dcbt	 %r7,	 %r17
	fadd	   3,	9,	15
	lfd	   22,	160(%r4)

	fadd	   4,	10,	16
	lfd	   23,	168(%r4)

	fadd	   5,	11,	17
	lfd	   24,	192(%r4)
	dcbt	 %r7,	 %r18
	fadd	   0,	0,	18
	lfd	   25,	200(%r4)

	fadd	   1,	1,	19
	lfd	   26,	208(%r4)

	fadd	   2,	2,	20
	lfd	   27,	216(%r4)
	dcbt	 %r7,	 %r19
	fadd	   3,	3,	21
	lfd	   28,	224(%r4)

	fadd	   4,	4,	22
	lfd	   29,	232(%r4)

	fadd	   5,	5,	23
	lfd	   6,	256(%r4)
	dcbt	 %r7,	 %r20
	fadd	   0,	0,	24
	lfd	   7,	264(%r4)

	fadd	   1,	1,	25
	lfd	   8,	272(%r4)

	fadd	   2,	2,	26
	lfd	   9,	280(%r4)
	dcbt	 %r7,	 %r21
	fadd	   3,	3,	27
	lfd	   10,	288(%r4)

	fadd	   4,	4,	28
	lfd	   11,	296(%r4)

	fadd	   5,	5,	29
	lfd	   12,	320(%r4)
	dcbt	 %r7,	 %r22
	fadd	   0,	0,	6
	lfd	   13,	328(%r4)

	fadd	   1,	1,	7
	lfd	   14,	336(%r4)

	fadd	   2,	2,	8
	lfd	   15,	344(%r4)
	dcbt	 %r7,	 %r23
	fadd	   3,	3,	9
	lfd	   16,	352(%r4)
	addi %r8, %r8, 1024 
	fadd	   4,	4,	10
	lfd	   17,	360(%r4)
	addi %r6, %r6, 1024 
	fadd	   5,	5,	11
	lfd	   18,	384(%r4)
	dcbt	 %r7,	 %r24
	fadd	   0,	0,	12
	lfd	   19,	392(%r4)

	fadd	   1,	1,	13
	lfd	   20,	400(%r4)

	fadd	   2,	2,	14
	lfd	   21,	408(%r4)
	dcbt	 %r7,	 %r25
	fadd	   3,	3,	15
	lfd	   22,	416(%r4)

	fadd	   4,	4,	16
	lfd	   23,	424(%r4)

	fadd	   5,	5,	17
	lfd	   24,	448(%r4)
	dcbt	 %r7,	 %r26
	fadd	   0,	0,	18
	lfd	   25,	456(%r4)

	fadd	   1,	1,	19
	lfd	   26,	464(%r4)

	fadd	   2,	2,	20
	lfd	   27,	472(%r4)
	dcbt	 %r7,	 %r27
	fadd	   3,	3,	21
	lfd	   28,	480(%r4)

	fadd	   4,	4,	22
	lfd	   29,	488(%r4)

	fadd	   5,	5,	23
	lfd	   6,	512(%r4)
	dcbt	 %r7,	 %r28
	fadd	   0,	0,	24
	lfd	   7,	520(%r4)

	fadd	   1,	1,	25
	lfd	   8,	528(%r4)

	fadd	   2,	2,	26
	lfd	   9,	536(%r4)
	dcbt	 %r7,	 %r29
	fadd	   3,	3,	27
	lfd	   10,	544(%r4)

	fadd	   4,	4,	28
	lfd	   11,	552(%r4)

	fadd	   5,	5,	29
	lfd	   12,	576(%r4)
	dcbt	 %r7,	 %r30
	subic.	 3,	3,	1
	lfd	   13,	584(%r4)

	fadd	   0,	0,	6
	lfd	   14,	592(%r4)

	fadd	   1,	1,	7
	lfd	   15,	600(%r4)
	dcbt	 %r7,	 %r31
	fadd	   2,	2,	8
	lfd	   16,	608(%r4)

	fadd	   3,	3,	9
	lfd	   17,	616(%r4)

	fadd	   4,	4,	10
	addi %r5, %r5, 48 
	dcbt	 %r8,	 %r16
	fadd	   5,	5,	11
	lfd	   18,	640(%r6)

	fadd	   0,	0,	12
	lfd	   19,	648(%r6)

	fadd	   1,	1,	13
	lfd	   20,	656(%r6)
	dcbt	 %r8,	 %r17
	fadd	   2,	2,	14
	lfd	   21,	664(%r6)

	fadd	   3,	3,	15
	lfd	   22,	672(%r6)

	fadd	   4,	4,	16
	lfd	   23,	680(%r6)
	dcbt	 %r8,	 %r18
	fadd	   5,	5,	17
	lfd	   24,	704(%r6)

	fadd	   0,	0,	18
	lfd	   25,	712(%r6)

	fadd	   1,	1,	19
	lfd	   26,	720(%r6)
	dcbt	 %r8,	 %r19
	fadd	   2,	2,	20
	lfd	   27,	728(%r6)

	fadd	   3,	3,	21
	lfd	   28,	736(%r6)

	fadd	   4,	4,	22
	lfd	   29,	744(%r6)
	dcbt	 %r8,	 %r20
	fadd	   5,	5,	23
	lfd	   6,	768(%r6)

	fadd	   0,	0,	24
	lfd	   7,	776(%r6)

	fadd	   1,	1,	25
	lfd	   8,	784(%r6)
	dcbt	 %r8,	 %r21
	fadd	   2,	2,	26
	lfd	   9,	792(%r6)

	fadd	   3,	3,	27
	lfd	   10,	800(%r6)

	fadd	   4,	4,	28
	lfd	   11,	808(%r6)
	dcbt	 %r8,	 %r22
	fadd	   5,	5,	29
	lfd	   12,	832(%r6)

	fadd	   0,	0,	6
	lfd	   13,	840(%r6)

	fadd	   1,	1,	7
	lfd	   14,	848(%r6)
	dcbt	 %r8,	 %r23
	fadd	   2,	2,	8
	lfd	   15,	856(%r6)

	fadd	   3,	3,	9
	lfd	   16,	864(%r6)

	fadd	   4,	4,	10
	lfd	   17,	872(%r6)
	addi %r4, %r4, 1024 
	fadd	   5,	5,	11
	lfd	   18,	896(%r6)
	addi %r7, %r7, 1024 
	fadd	   0,	0,	12
	lfd	   19,	904(%r6)
	dcbt	 %r8,	 %r24
	fadd	   1,	1,	13
	lfd	   20,	912(%r6)

	fadd	   2,	2,	14
	lfd	   21,	920(%r6)

	fadd	   3,	3,	15
	lfd	   22,	928(%r6)
	dcbt	 %r8,	 %r25
	fadd	   4,	4,	16
	lfd	   23,	936(%r6)

	fadd	   5,	5,	17
	lfd	   24,	960(%r6)

	fadd	   0,	0,	18
	lfd	   25,	968(%r6)
	dcbt	 %r8,	 %r26
	fadd	   1,	1,	19
	lfd	   26,	976(%r6)

	fadd	   2,	2,	20
	lfd	   27,	984(%r6)

	fadd	   3,	3,	21
	lfd	   28,	992(%r6)
	dcbt	 %r8,	 %r27
	fadd	   4,	4,	22
	lfd	   29,	1000(%r6)

	fadd	   5,	5,	23
	lfd	   6,	1024(%r6)

	fadd	   0,	0,	24
	lfd	   7,	1032(%r6)
	dcbt	 %r8,	 %r28
	fadd	   1,	1,	25
	lfd	   8,	1040(%r6)

	fadd	   2,	2,	26
	lfd	   9,	1048(%r6)

	fadd	   3,	3,	27
	lfd	   10,	1056(%r6)
	dcbt	 %r8,	 %r29
	fadd	   4,	4,	28
	lfd	   11,	1064(%r6)

	fadd	   5,	5,	29
#<<<<<< Break 1

	lfd	   12,	1088(%r6)
	dcbt	 %r8,	 %r30
	lfd	   13,	1096(%r6)

	lfd	   14,	1104(%r6)

	lfd	   15,	1112(%r6)
	dcbt	 %r8,	 %r31
	lfd	   16,	1120(%r6)

	lfd	   17,	1128(%r6)

	stfd	   0,	0(%r5)
	stfd	   1,	8(%r5)
	stfd	   2,	16(%r5)
	stfd	   3,	24(%r5)
	stfd	   4,	32(%r5)
	stfd	   5,	40(%r5)
	bf gt,	lab2
	b  lab1
lab2:
	fadd	   0,	6,	12
	lfd	   18,	128(%r4)
	fadd	   1,	7,	13
	lfd	   19,	136(%r4)
	fadd	   2,	8,	14
	lfd	   20,	144(%r4)
	fadd	   3,	9,	15
	lfd	   21,	152(%r4)
	fadd	   4,	10,	16
	lfd	   22,	160(%r4)
	fadd	   5,	11,	17
	lfd	   23,	168(%r4)
	fadd	   0,	0,	18
	lfd	   24,	192(%r4)
	fadd	   1,	1,	19
	lfd	   25,	200(%r4)
	fadd	   2,	2,	20
	lfd	   26,	208(%r4)
	fadd	   3,	3,	21
	lfd	   27,	216(%r4)
	fadd	   4,	4,	22
	lfd	   28,	224(%r4)
	fadd	   5,	5,	23
	lfd	   29,	232(%r4)
	fadd	   0,	0,	24
	lfd	   6,	256(%r4)
	fadd	   1,	1,	25
	lfd	   7,	264(%r4)
	fadd	   2,	2,	26
	lfd	   8,	272(%r4)
	fadd	   3,	3,	27
	lfd	   9,	280(%r4)
	fadd	   4,	4,	28
	lfd	   10,	288(%r4)
	fadd	   5,	5,	29
	lfd	   11,	296(%r4)
	fadd	   0,	0,	6
	lfd	   12,	320(%r4)
	fadd	   1,	1,	7
	lfd	   13,	328(%r4)
	fadd	   2,	2,	8
	lfd	   14,	336(%r4)
	fadd	   3,	3,	9
	lfd	   15,	344(%r4)
	fadd	   4,	4,	10
	lfd	   16,	352(%r4)
	fadd	   5,	5,	11
	lfd	   17,	360(%r4)
	fadd	   0,	0,	12
	addi %r6, %r6, 1024 
	fadd	   1,	1,	13
	lfd	   18,	384(%r4)
	fadd	   2,	2,	14
	lfd	   19,	392(%r4)
	fadd	   3,	3,	15
	lfd	   20,	400(%r4)
	fadd	   4,	4,	16
	lfd	   21,	408(%r4)
	fadd	   5,	5,	17
	lfd	   22,	416(%r4)
	fadd	   0,	0,	18
	lfd	   23,	424(%r4)
	fadd	   1,	1,	19
	lfd	   24,	448(%r4)
	fadd	   2,	2,	20
	lfd	   25,	456(%r4)
	fadd	   3,	3,	21
	lfd	   26,	464(%r4)
	fadd	   4,	4,	22
	lfd	   27,	472(%r4)
	fadd	   5,	5,	23
	lfd	   28,	480(%r4)
	fadd	   0,	0,	24
	lfd	   29,	488(%r4)
	fadd	   1,	1,	25
	lfd	   6,	512(%r4)
	fadd	   2,	2,	26
	lfd	   7,	520(%r4)
	fadd	   3,	3,	27
	lfd	   8,	528(%r4)
	fadd	   4,	4,	28
	lfd	   9,	536(%r4)
	fadd	   5,	5,	29
	lfd	   10,	544(%r4)
	fadd	   0,	0,	6
	lfd	   11,	552(%r4)
	fadd	   1,	1,	7
	lfd	   12,	576(%r4)
	fadd	   2,	2,	8
	lfd	   13,	584(%r4)
	fadd	   3,	3,	9
	lfd	   14,	592(%r4)
	fadd	   4,	4,	10
	lfd	   15,	600(%r4)
	fadd	   5,	5,	11
	lfd	   16,	608(%r4)
	fadd	   0,	0,	12
	lfd	   17,	616(%r4)
	fadd	   1,	1,	13
	addi %r5, %r5, 48 
	fadd	   2,	2,	14
	lfd	   18,	640(%r6)
	fadd	   3,	3,	15
	lfd	   19,	648(%r6)
	fadd	   4,	4,	16
	lfd	   20,	656(%r6)
	fadd	   5,	5,	17
	lfd	   21,	664(%r6)
	fadd	   0,	0,	18
	lfd	   22,	672(%r6)
	fadd	   1,	1,	19
	lfd	   23,	680(%r6)
	fadd	   2,	2,	20
	lfd	   24,	704(%r6)
	fadd	   3,	3,	21
	lfd	   25,	712(%r6)
	fadd	   4,	4,	22
	lfd	   26,	720(%r6)
	fadd	   5,	5,	23
	lfd	   27,	728(%r6)
	fadd	   0,	0,	24
	lfd	   28,	736(%r6)
	fadd	   1,	1,	25
	lfd	   29,	744(%r6)
	fadd	   2,	2,	26
	lfd	   6,	768(%r6)
	fadd	   3,	3,	27
	lfd	   7,	776(%r6)
	fadd	   4,	4,	28
	lfd	   8,	784(%r6)
	fadd	   5,	5,	29
	lfd	   9,	792(%r6)
	fadd	   0,	0,	6
	lfd	   10,	800(%r6)
	fadd	   1,	1,	7
	lfd	   11,	808(%r6)
	fadd	   2,	2,	8
	lfd	   12,	832(%r6)
	fadd	   3,	3,	9
	lfd	   13,	840(%r6)
	fadd	   4,	4,	10
	lfd	   14,	848(%r6)
	fadd	   5,	5,	11
	lfd	   15,	856(%r6)
	fadd	   0,	0,	12
	lfd	   16,	864(%r6)
	fadd	   1,	1,	13
	lfd	   17,	872(%r6)
	fadd	   2,	2,	14
	lfd	   18,	896(%r6)
	fadd	   3,	3,	15
	lfd	   19,	904(%r6)
	fadd	   4,	4,	16
	lfd	   20,	912(%r6)
	fadd	   5,	5,	17
	lfd	   21,	920(%r6)
	fadd	   0,	0,	18
	lfd	   22,	928(%r6)
	fadd	   1,	1,	19
	lfd	   23,	936(%r6)
	fadd	   2,	2,	20
	lfd	   24,	960(%r6)
	fadd	   3,	3,	21
	lfd	   25,	968(%r6)
	fadd	   4,	4,	22
	lfd	   26,	976(%r6)
	fadd	   5,	5,	23
	lfd	   27,	984(%r6)
	fadd	   0,	0,	24
	lfd	   28,	992(%r6)
	fadd	   1,	1,	25
	lfd	   29,	1000(%r6)
	fadd	   2,	2,	26

	fadd	   3,	3,	27

	fadd	   4,	4,	28

	fadd	   5,	5,	29









	stfd	   0,	0(%r5)
	stfd	   1,	8(%r5)
	stfd	   2,	16(%r5)
	stfd	   3,	24(%r5)
	stfd	   4,	32(%r5)
	stfd	   5,	40(%r5)
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
