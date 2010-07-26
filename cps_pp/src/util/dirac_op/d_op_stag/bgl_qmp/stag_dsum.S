	.section	".text"
	.align 2
	.globl stag_dsum
	.type	stag_dsum,@function
stag_dsum:
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
	lfd	   30,	0(%r6)
	lfd	   31,	0(%r7)
	addi %r11, %r5, -48
	addi %r10, %r4, -512
	addi %r8, %r4, 384
	addi %r9, %r4, 128
	dcbt	 %r16,	 %r4
	dcbt	 %r17,	 %r4
	dcbt	 %r18,	 %r4
	dcbt	 %r19,	 %r4
	dcbt	 %r20,	 %r4
	dcbt	 %r21,	 %r4
	dcbt	 %r22,	 %r4
	dcbt	 %r23,	 %r4
	dcbt	 %r16,	 %r10
	dcbt	 %r17,	 %r10
	dcbt	 %r18,	 %r10
	dcbt	 %r19,	 %r10
	dcbt	 %r20,	 %r10
	dcbt	 %r21,	 %r10
	dcbt	 %r22,	 %r10
	dcbt	 %r23,	 %r10
	or.	 3,	3,	3
	bf gt,	stag_dsum_lab0
	lfd	   12,	0(%r4)
	lfd	   13,	8(%r4)
	lfd	   14,	16(%r4)
	lfd	   15,	24(%r4)
	lfd	   16,	32(%r4)
	lfd	   17,	40(%r4)
	lfd	   0,	0(%r5)
	lfd	   1,	8(%r5)
	lfd	   2,	16(%r5)
	lfd	   3,	24(%r5)
	lfd	   4,	32(%r5)
	lfd	   5,	40(%r5)
	subic.	 3,	3,	1
	bf gt,	stag_dsum_lab2
	b  stag_dsum_lab1
stag_dsum_lab1:
	lfd	   18,	64(%r4)
	dcbt	 %r8,	 %r16
	fmadd	0,	30,	12,	0
	lfd	   19,	72(%r4)

	fmadd	1,	30,	13,	1
	lfd	   20,	80(%r4)

	fmadd	2,	30,	14,	2
	lfd	   21,	88(%r4)
	dcbt	 %r8,	 %r17
	fmadd	3,	30,	15,	3
	lfd	   22,	96(%r4)
	addi %r5, %r5, 48 
	fmadd	4,	30,	16,	4
	lfd	   23,	104(%r4)
	addi %r11, %r11, 48 
	fmadd	5,	30,	17,	5
	lfd	   24,	128(%r4)
	dcbt	 %r8,	 %r18
	fmadd	0,	31,	18,	0
	lfd	   25,	136(%r4)

	fmadd	1,	31,	19,	1
	lfd	   26,	144(%r4)

	fmadd	2,	31,	20,	2
	lfd	   27,	152(%r4)
	dcbt	 %r8,	 %r19
	fmadd	3,	31,	21,	3
	lfd	   28,	160(%r4)
	addi %r9, %r9, 512 
	fmadd	4,	31,	22,	4
	lfd	   29,	168(%r4)
	addi %r10, %r10, 512 
	fmadd	5,	31,	23,	5
	lfd	   12,	192(%r4)
	dcbt	 %r8,	 %r20
	fmadd	0,	30,	24,	0
	lfd	   13,	200(%r4)
	dcbt	 %r8,	 %r21
	fmadd	1,	30,	25,	1
	lfd	   14,	208(%r4)

	fmadd	2,	30,	26,	2
	lfd	   15,	216(%r4)
	dcbt	 %r5,	 %r16
	fmadd	3,	30,	27,	3
	lfd	   16,	224(%r4)
	dcbt	 %r5,	 %r17
	fmadd	4,	30,	28,	4
	lfd	   17,	232(%r4)

	fmadd	5,	30,	29,	5
	lfd	   18,	256(%r4)
	dcbt	 %r8,	 %r22
	fmadd	0,	31,	12,	0
	lfd	   19,	264(%r4)

	fmadd	1,	31,	13,	1
	lfd	   20,	272(%r4)

	fmadd	2,	31,	14,	2
	lfd	   21,	280(%r4)
	dcbt	 %r8,	 %r23
	fmadd	3,	31,	15,	3
	lfd	   22,	288(%r4)

	fmadd	4,	31,	16,	4
	lfd	   23,	296(%r4)

	fmadd	5,	31,	17,	5
	lfd	   24,	320(%r4)
	dcbt	 %r9,	 %r16
	subic.	 3,	3,	1
	lfd	   25,	328(%r4)

	fmadd	6,	30,	18,	0
	lfd	   26,	336(%r4)

	fmadd	7,	30,	19,	1
	lfd	   27,	344(%r4)
	dcbt	 %r9,	 %r17
	fmadd	8,	30,	20,	2
	lfd	   28,	352(%r4)

	fmadd	9,	30,	21,	3
	lfd	   29,	360(%r4)

	fmadd	10,	30,	22,	4
	lfd	   12,	384(%r4)
	dcbt	 %r9,	 %r18
	fmadd	11,	30,	23,	5
	lfd	   13,	392(%r4)

	fmadd	6,	31,	24,	6
	lfd	   14,	400(%r4)

	fmadd	7,	31,	25,	7
	lfd	   15,	408(%r4)
	dcbt	 %r9,	 %r19
	fmadd	8,	31,	26,	8
	lfd	   16,	416(%r4)

	fmadd	9,	31,	27,	9
	lfd	   17,	424(%r4)

	fmadd	10,	31,	28,	10
	lfd	   18,	448(%r10)
	dcbt	 %r9,	 %r20
	fmadd	11,	31,	29,	11
	lfd	   19,	456(%r10)

	fmadd	6,	30,	12,	6
	lfd	   20,	464(%r10)

	fmadd	7,	30,	13,	7
	lfd	   21,	472(%r10)
	dcbt	 %r9,	 %r21
	fmadd	8,	30,	14,	8
	lfd	   22,	480(%r10)
	addi %r4, %r4, 512 
	fmadd	9,	30,	15,	9
	lfd	   23,	488(%r10)
	addi %r8, %r8, 512 
	fmadd	10,	30,	16,	10
	lfd	   12,	512(%r10)
	dcbt	 %r9,	 %r22
	fmadd	11,	30,	17,	11
	lfd	   13,	520(%r10)

	fmadd	6,	31,	18,	6
	lfd	   14,	528(%r10)

	fmadd	7,	31,	19,	7
	lfd	   15,	536(%r10)
	dcbt	 %r9,	 %r23
	fmadd	8,	31,	20,	8
	lfd	   16,	544(%r10)

	fmadd	9,	31,	21,	9
	lfd	   17,	552(%r10)

	fmadd	10,	31,	22,	10
	stfd	   6,	0(%r11)
	stfd	   7,	8(%r11)
	fmadd	11,	31,	23,	11
	stfd	   8,	16(%r11)
	stfd	   9,	24(%r11)
	stfd	   10,	32(%r11)
	stfd	   11,	40(%r11)
	lfd	   0,	0(%r5)
	lfd	   1,	8(%r5)
	lfd	   2,	16(%r5)
	lfd	   3,	24(%r5)
	lfd	   4,	32(%r5)
	lfd	   5,	40(%r5)
	bf gt,	stag_dsum_lab2
	b  stag_dsum_lab1
stag_dsum_lab2:
	fmadd	0,	30,	12,	0
	lfd	   18,	64(%r4)
	fmadd	1,	30,	13,	1
	lfd	   19,	72(%r4)
	fmadd	2,	30,	14,	2
	lfd	   20,	80(%r4)
	fmadd	3,	30,	15,	3
	lfd	   21,	88(%r4)
	fmadd	4,	30,	16,	4
	lfd	   22,	96(%r4)
	fmadd	5,	30,	17,	5
	lfd	   23,	104(%r4)
	fmadd	0,	31,	18,	0
	addi %r11, %r11, 48 
	fmadd	1,	31,	19,	1
	addi %r10, %r10, 512 
	fmadd	2,	31,	20,	2
	lfd	   24,	128(%r4)
	fmadd	3,	31,	21,	3
	lfd	   25,	136(%r4)
	fmadd	4,	31,	22,	4
	lfd	   26,	144(%r4)
	fmadd	5,	31,	23,	5
	lfd	   27,	152(%r4)
	fmadd	0,	30,	24,	0
	lfd	   28,	160(%r4)
	fmadd	1,	30,	25,	1
	lfd	   29,	168(%r4)
	fmadd	2,	30,	26,	2
	lfd	   12,	192(%r4)
	fmadd	3,	30,	27,	3
	lfd	   13,	200(%r4)
	fmadd	4,	30,	28,	4
	lfd	   14,	208(%r4)
	fmadd	5,	30,	29,	5
	lfd	   15,	216(%r4)
	fmadd	0,	31,	12,	0
	lfd	   16,	224(%r4)
	fmadd	1,	31,	13,	1
	lfd	   17,	232(%r4)
	fmadd	2,	31,	14,	2
	lfd	   18,	256(%r4)
	fmadd	3,	31,	15,	3
	lfd	   19,	264(%r4)
	fmadd	4,	31,	16,	4
	lfd	   20,	272(%r4)
	fmadd	5,	31,	17,	5
	lfd	   21,	280(%r4)
	subic.	 3,	3,	1
	lfd	   22,	288(%r4)
	fmadd	6,	30,	18,	0
	lfd	   23,	296(%r4)
	fmadd	7,	30,	19,	1
	lfd	   24,	320(%r4)
	fmadd	8,	30,	20,	2
	lfd	   25,	328(%r4)
	fmadd	9,	30,	21,	3
	lfd	   26,	336(%r4)
	fmadd	10,	30,	22,	4
	lfd	   27,	344(%r4)
	fmadd	11,	30,	23,	5
	lfd	   28,	352(%r4)
	fmadd	6,	31,	24,	6
	lfd	   29,	360(%r4)
	fmadd	7,	31,	25,	7
	lfd	   12,	384(%r4)
	fmadd	8,	31,	26,	8
	lfd	   13,	392(%r4)
	fmadd	9,	31,	27,	9
	lfd	   14,	400(%r4)
	fmadd	10,	31,	28,	10
	lfd	   15,	408(%r4)
	fmadd	11,	31,	29,	11
	lfd	   16,	416(%r4)
	fmadd	6,	30,	12,	6
	lfd	   17,	424(%r4)
	fmadd	7,	30,	13,	7
	lfd	   18,	448(%r10)
	fmadd	8,	30,	14,	8
	lfd	   19,	456(%r10)
	fmadd	9,	30,	15,	9
	lfd	   20,	464(%r10)
	fmadd	10,	30,	16,	10
	lfd	   21,	472(%r10)
	fmadd	11,	30,	17,	11
	lfd	   22,	480(%r10)
	fmadd	6,	31,	18,	6
	lfd	   23,	488(%r10)
	fmadd	7,	31,	19,	7

	fmadd	8,	31,	20,	8

	fmadd	9,	31,	21,	9

	fmadd	10,	31,	22,	10

	fmadd	11,	31,	23,	11


	stfd	   6,	0(%r11)
	stfd	   7,	8(%r11)
	stfd	   8,	16(%r11)
	stfd	   9,	24(%r11)
	stfd	   10,	32(%r11)
	stfd	   11,	40(%r11)
stag_dsum_lab0:
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
