	.section	".text"
	.align 2
	.globl dirac_sum2_64
	.type	dirac_sum2_64,@function
dirac_sum2_64:
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
	li   27,	16
	li   28,	48
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
	lfd	   6,	0(%r5)
	lfd	   7,	8(%r5)
	lfd	   8,	16(%r5)
	lfd	   9,	24(%r5)
	lfd	   10,	32(%r5)
	lfd	   11,	40(%r5)
	or.	 3,	3,	3
	lfd	   12,	64(%r5)
	lfd	   13,	72(%r5)
	lfd	   14,	80(%r5)
	lfd	   15,	88(%r5)
	lfd	   16,	96(%r5)
	lfd	   17,	104(%r5)
	bf gt,	lab0
	subic.	 3,	3,	1
	bf gt,	lab2
	b  dirac_sum_64_lab1
dirac_sum_64_lab1:
	fadd	   0,	6,	12
	lfd	   18,	128(%r5)
	fadd	   1,	7,	13
	lfd	   19,	136(%r5)
	fadd	   2,	8,	14
	lfd	   20,	144(%r5)
	fadd	   3,	9,	15
	lfd	   21,	152(%r5)
	fadd	   4,	10,	16
	lfd	   22,	160(%r5)
	fadd	   5,	11,	17
	lfd	   23,	168(%r5)
	subic.	 3,	3,	1
	lfd	   24,	192(%r5)
	fadd	   0,	0,	18
	lfd	   25,	200(%r5)
	fadd	   1,	1,	19
	lfd	   26,	208(%r5)
	fadd	   2,	2,	20
	lfd	   27,	216(%r5)
	fadd	   3,	3,	21
	lfd	   28,	224(%r5)
	fadd	   4,	4,	22
	lfd	   29,	232(%r5)
	fadd	   5,	5,	23
	lfd	   6,	256(%r5)
	fadd	   0,	0,	24
	lfd	   7,	264(%r5)
	fadd	   1,	1,	25
	lfd	   8,	272(%r5)
	fadd	   2,	2,	26
	lfd	   9,	280(%r5)
	fadd	   3,	3,	27
	lfd	   10,	288(%r5)
	fadd	   4,	4,	28
	lfd	   11,	296(%r5)
	fadd	   5,	5,	29
	lfd	   12,	320(%r5)
	fadd	   0,	0,	6
	lfd	   13,	328(%r5)
	fadd	   1,	1,	7
	lfd	   14,	336(%r5)
	fadd	   2,	2,	8
	lfd	   15,	344(%r5)
	fadd	   3,	3,	9
	lfd	   16,	352(%r5)
	fadd	   4,	4,	10
	lfd	   17,	360(%r5)
	fadd	   5,	5,	11
	lfd	   18,	384(%r5)
	fadd	   0,	0,	12
	lfd	   19,	392(%r5)
	fadd	   1,	1,	13
	lfd	   20,	400(%r5)
	fadd	   2,	2,	14
	lfd	   21,	408(%r5)
	fadd	   3,	3,	15
	lfd	   22,	416(%r5)
	fadd	   4,	4,	16
	lfd	   23,	424(%r5)
	fadd	   5,	5,	17
	lfd	   24,	448(%r5)
	fadd	   0,	0,	18
	lfd	   25,	456(%r5)
	fadd	   1,	1,	19
	lfd	   26,	464(%r5)
	fadd	   2,	2,	20
	lfd	   27,	472(%r5)
	fadd	   3,	3,	21
	lfd	   28,	480(%r5)
	fadd	   4,	4,	22
	lfd	   29,	488(%r5)
	fadd	   5,	5,	23
	lfd	   6,	512(%r5)
	fadd	   0,	0,	24
	lfd	   7,	520(%r5)
	fadd	   1,	1,	25
	lfd	   8,	528(%r5)
	fadd	   2,	2,	26
	lfd	   9,	536(%r5)
	fadd	   3,	3,	27
	lfd	   10,	544(%r5)
	fadd	   4,	4,	28
	lfd	   11,	552(%r5)
	fadd	   5,	5,	29
	lfd	   12,	576(%r5)
	fadd	   0,	0,	6
	lfd	   13,	584(%r5)
	fadd	   1,	1,	7
	lfd	   14,	592(%r5)
	fadd	   2,	2,	8
	lfd	   15,	600(%r5)
	fadd	   3,	3,	9
	lfd	   16,	608(%r5)
	fadd	   4,	4,	10
	lfd	   17,	616(%r5)
	fadd	   5,	5,	11
	lfd	   18,	640(%r5)
	fadd	   0,	0,	12
	lfd	   19,	648(%r5)
	fadd	   1,	1,	13
	lfd	   20,	656(%r5)
	fadd	   2,	2,	14
	lfd	   21,	664(%r5)
	fadd	   3,	3,	15
	lfd	   22,	672(%r5)
	fadd	   4,	4,	16
	lfd	   23,	680(%r5)
	fadd	   5,	5,	17
	lfd	   24,	704(%r5)
	fadd	   0,	0,	18
	lfd	   25,	712(%r5)
	fadd	   1,	1,	19
	lfd	   26,	720(%r5)
	fadd	   2,	2,	20
	lfd	   27,	728(%r5)
	fadd	   3,	3,	21
	lfd	   28,	736(%r5)
	fadd	   4,	4,	22
	lfd	   29,	744(%r5)
	fadd	   5,	5,	23
	lfd	   6,	768(%r5)
	fadd	   0,	0,	24
	lfd	   7,	776(%r5)
	fadd	   1,	1,	25
	lfd	   8,	784(%r5)
	fadd	   2,	2,	26
	lfd	   9,	792(%r5)
	fadd	   3,	3,	27
	lfd	   10,	800(%r5)
	fadd	   4,	4,	28
	lfd	   11,	808(%r5)
	fadd	   5,	5,	29
	lfd	   12,	832(%r5)
	fadd	   0,	0,	6
	lfd	   13,	840(%r5)
	fadd	   1,	1,	7
	lfd	   14,	848(%r5)
	fadd	   2,	2,	8
	lfd	   15,	856(%r5)
	fadd	   3,	3,	9
	lfd	   16,	864(%r5)
	fadd	   4,	4,	10
	lfd	   17,	872(%r5)
	fadd	   5,	5,	11
	lfd	   18,	896(%r5)
	fadd	   0,	0,	12
	lfd	   19,	904(%r5)
	fadd	   1,	1,	13
	lfd	   20,	912(%r5)
	fadd	   2,	2,	14
	lfd	   21,	920(%r5)
	fadd	   3,	3,	15
	lfd	   22,	928(%r5)
	fadd	   4,	4,	16
	lfd	   23,	936(%r5)
	fadd	   5,	5,	17
	lfd	   24,	960(%r5)
	fadd	   0,	0,	18
	lfd	   25,	968(%r5)
	fadd	   1,	1,	19
	lfd	   26,	976(%r5)
	fadd	   2,	2,	20
	lfd	   27,	984(%r5)
	fadd	   3,	3,	21
	lfd	   28,	992(%r5)
	fadd	   4,	4,	22
	lfd	   29,	1000(%r5)
	fadd	   5,	5,	23
	lfd	   6,	1024(%r5)
	fadd	   0,	0,	24
	lfd	   7,	1032(%r5)
	fadd	   1,	1,	25
	lfd	   8,	1040(%r5)
	fadd	   2,	2,	26
	lfd	   9,	1048(%r5)
	fadd	   3,	3,	27
	lfd	   10,	1056(%r5)
	fadd	   4,	4,	28
	lfd	   11,	1064(%r5)
	fadd	   5,	5,	29
	lfd	   12,	1088(%r5)
	lfd	   13,	1096(%r5)
	lfd	   14,	1104(%r5)
	lfd	   15,	1112(%r5)
	lfd	   16,	1120(%r5)
	lfd	   17,	1128(%r5)
	stfd	   0,	0(%r6)
	stfd	   1,	8(%r6)
	addi %r5, %r5, 1024
	stfd	   2,	16(%r6)
	stfd	   3,	24(%r6)
	stfd	   4,	32(%r6)
	stfd	   5,	40(%r6)
	addi %r6, %r6, 48
	bf gt,	lab2
	b  dirac_sum_64_lab1
lab2:
	fadd	   0,	6,	12
	lfd	   18,	128(%r5)
	fadd	   1,	7,	13
	lfd	   19,	136(%r5)
	fadd	   2,	8,	14
	lfd	   20,	144(%r5)
	fadd	   3,	9,	15
	lfd	   21,	152(%r5)
	fadd	   4,	10,	16
	lfd	   22,	160(%r5)
	fadd	   5,	11,	17
	lfd	   23,	168(%r5)
	fadd	   0,	0,	18
	lfd	   24,	192(%r5)
	fadd	   1,	1,	19
	lfd	   25,	200(%r5)
	fadd	   2,	2,	20
	lfd	   26,	208(%r5)
	fadd	   3,	3,	21
	lfd	   27,	216(%r5)
	fadd	   4,	4,	22
	lfd	   28,	224(%r5)
	fadd	   5,	5,	23
	lfd	   29,	232(%r5)
	fadd	   0,	0,	24
	lfd	   6,	256(%r5)
	fadd	   1,	1,	25
	lfd	   7,	264(%r5)
	fadd	   2,	2,	26
	lfd	   8,	272(%r5)
	fadd	   3,	3,	27
	lfd	   9,	280(%r5)
	fadd	   4,	4,	28
	lfd	   10,	288(%r5)
	fadd	   5,	5,	29
	lfd	   11,	296(%r5)
	fadd	   0,	0,	6
	lfd	   12,	320(%r5)
	fadd	   1,	1,	7
	lfd	   13,	328(%r5)
	fadd	   2,	2,	8
	lfd	   14,	336(%r5)
	fadd	   3,	3,	9
	lfd	   15,	344(%r5)
	fadd	   4,	4,	10
	lfd	   16,	352(%r5)
	fadd	   5,	5,	11
	lfd	   17,	360(%r5)
	fadd	   0,	0,	12
	lfd	   18,	384(%r5)
	fadd	   1,	1,	13
	lfd	   19,	392(%r5)
	fadd	   2,	2,	14
	lfd	   20,	400(%r5)
	fadd	   3,	3,	15
	lfd	   21,	408(%r5)
	fadd	   4,	4,	16
	lfd	   22,	416(%r5)
	fadd	   5,	5,	17
	lfd	   23,	424(%r5)
	fadd	   0,	0,	18
	lfd	   24,	448(%r5)
	fadd	   1,	1,	19
	lfd	   25,	456(%r5)
	fadd	   2,	2,	20
	lfd	   26,	464(%r5)
	fadd	   3,	3,	21
	lfd	   27,	472(%r5)
	fadd	   4,	4,	22
	lfd	   28,	480(%r5)
	fadd	   5,	5,	23
	lfd	   29,	488(%r5)
	fadd	   0,	0,	24
	lfd	   6,	512(%r5)
	fadd	   1,	1,	25
	lfd	   7,	520(%r5)
	fadd	   2,	2,	26
	lfd	   8,	528(%r5)
	fadd	   3,	3,	27
	lfd	   9,	536(%r5)
	fadd	   4,	4,	28
	lfd	   10,	544(%r5)
	fadd	   5,	5,	29
	lfd	   11,	552(%r5)
	fadd	   0,	0,	6
	lfd	   12,	576(%r5)
	fadd	   1,	1,	7
	lfd	   13,	584(%r5)
	fadd	   2,	2,	8
	lfd	   14,	592(%r5)
	fadd	   3,	3,	9
	lfd	   15,	600(%r5)
	fadd	   4,	4,	10
	lfd	   16,	608(%r5)
	fadd	   5,	5,	11
	lfd	   17,	616(%r5)
	fadd	   0,	0,	12
	lfd	   18,	640(%r5)
	fadd	   1,	1,	13
	lfd	   19,	648(%r5)
	fadd	   2,	2,	14
	lfd	   20,	656(%r5)
	fadd	   3,	3,	15
	lfd	   21,	664(%r5)
	fadd	   4,	4,	16
	lfd	   22,	672(%r5)
	fadd	   5,	5,	17
	lfd	   23,	680(%r5)
	fadd	   0,	0,	18
	lfd	   24,	704(%r5)
	fadd	   1,	1,	19
	lfd	   25,	712(%r5)
	fadd	   2,	2,	20
	lfd	   26,	720(%r5)
	fadd	   3,	3,	21
	lfd	   27,	728(%r5)
	fadd	   4,	4,	22
	lfd	   28,	736(%r5)
	fadd	   5,	5,	23
	lfd	   29,	744(%r5)
	fadd	   0,	0,	24
	lfd	   6,	768(%r5)
	fadd	   1,	1,	25
	lfd	   7,	776(%r5)
	fadd	   2,	2,	26
	lfd	   8,	784(%r5)
	fadd	   3,	3,	27
	lfd	   9,	792(%r5)
	fadd	   4,	4,	28
	lfd	   10,	800(%r5)
	fadd	   5,	5,	29
	lfd	   11,	808(%r5)
	fadd	   0,	0,	6
	lfd	   12,	832(%r5)
	fadd	   1,	1,	7
	lfd	   13,	840(%r5)
	fadd	   2,	2,	8
	lfd	   14,	848(%r5)
	fadd	   3,	3,	9
	lfd	   15,	856(%r5)
	fadd	   4,	4,	10
	lfd	   16,	864(%r5)
	fadd	   5,	5,	11
	lfd	   17,	872(%r5)
	fadd	   0,	0,	12
	lfd	   18,	896(%r5)
	fadd	   1,	1,	13
	lfd	   19,	904(%r5)
	fadd	   2,	2,	14
	lfd	   20,	912(%r5)
	fadd	   3,	3,	15
	lfd	   21,	920(%r5)
	fadd	   4,	4,	16
	lfd	   22,	928(%r5)
	fadd	   5,	5,	17
	lfd	   23,	936(%r5)
	fadd	   0,	0,	18
	lfd	   24,	960(%r5)
	fadd	   1,	1,	19
	lfd	   25,	968(%r5)
	fadd	   2,	2,	20
	lfd	   26,	976(%r5)
	fadd	   3,	3,	21
	lfd	   27,	984(%r5)
	fadd	   4,	4,	22
	lfd	   28,	992(%r5)
	fadd	   5,	5,	23
	lfd	   29,	1000(%r5)
	fadd	   0,	0,	24

	fadd	   1,	1,	25

	fadd	   2,	2,	26

	fadd	   3,	3,	27

	fadd	   4,	4,	28

	fadd	   5,	5,	29







	stfd	   0,	0(%r6)
	stfd	   1,	8(%r6)
	stfd	   2,	16(%r6)
	stfd	   3,	24(%r6)
	stfd	   4,	32(%r6)
	stfd	   5,	40(%r6)
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
