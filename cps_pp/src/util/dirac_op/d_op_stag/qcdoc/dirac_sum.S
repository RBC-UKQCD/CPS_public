	.section	".text"
	.align 2
	.globl dirac_sum
	.type	 dirac_sum,@function
dirac_sum:
	la   %r1,	-864(%r1)
	
	stw   %r14,	368(%r1)
	stw   %r15,	376(%r1)
	stw   %r16,	384(%r1)
	stw   %r17,	392(%r1)
	stw   %r18,	400(%r1)
	stw   %r19,	408(%r1)
	stw   %r20,	416(%r1)
	stw   %r21,	424(%r1)
	stw   %r22,	432(%r1)

	stw   %r23,	440(%r1)
	li 14, 0
	stw   %r24,	448(%r1)
	li 15, 32
	stw   %r25,	456(%r1)
	stw   %r26,	464(%r1)
	dcbt %r14, %r4	
	stw   %r27,	472(%r1)
	dcbt %r15, %r4
	stw   %r28,	480(%r1)
	stw   %r29,	488(%r1)
	stw   %r30,	496(%r1)
	stw   %r31,	504(%r1)
	
	stfd  14,	112(%r1)
	stfd  15,	120(%r1)
	stfd  16,	128(%r1)
	stfd  17,	136(%r1)
	stfd  18,	144(%r1)
	stfd  19,	152(%r1)
	stfd  20,	160(%r1)
	stfd  21,	168(%r1)
	stfd  22,	176(%r1)
	stfd  23,	184(%r1)
	stfd  24,	192(%r1)
	stfd  25,	200(%r1)
	stfd  26,	208(%r1)
	stfd  27,	216(%r1)
	stfd  28,	224(%r1)
	stfd  29,	232(%r1)
	stfd  30,	240(%r1)
	stfd  31,	248(%r1)
	
	mr %r0, %r3
	or. %r0, %r0, %r0
	mtctr %r0
	bf gt, lab0

	lwz %r10, 0(%r4)
	lwz %r11, 4(%r4)
	lwz %r12, 8(%r4)
	lwz %r13, 12(%r4)
	lwz %r14, 16(%r4)
	lwz %r15, 20(%r4)
	add %r20, %r6, %r10	
	lwz %r16, 24(%r4)
	add %r21, %r5, %r11	
	lwz %r17, 28(%r4)
	add %r22, %r5, %r12	
	lwz %r18, 32(%r4)
	add %r23, %r5, %r13
	lfd  6, 0(%r21)	
	add %r24, %r5, %r14
	lfd  10, 32(%r21)	
	add %r25, %r5, %r15
	lfd  12, 0(%r22)	
	add %r26, %r5, %r16
	lfd  16, 32(%r22)	
	add %r27, %r5, %r17
	lfd  18, 0(%r23)	
	add %r28, %r5, %r18
	lfd  22, 32(%r23)
	lfd  7, 8(%r21)
	lfd  8, 16(%r21)
	lfd  9, 24(%r21)
	lfd  24, 0(%r24)
	lfd  11, 40(%r21)
	lfd  28, 32(%r24)	
	lfd  13, 8(%r22)
	lfd  14, 16(%r22)
	lfd  15, 24(%r22)
	lfd  17, 40(%r22)
	lfd  19, 8(%r23)
	lfd  20, 16(%r23)
	lfd  21, 24(%r23)
	lfd  23, 40(%r23)
	lfd  25, 8(%r24)
	lfd  26, 16(%r24)
	lfd  27, 24(%r24)
	lfd  29, 40(%r24)

	subic. %r0, %r0, 1
	beq lab2
	b lab1
.align 4
lab1:
	subic. %r0, %r0, 1	
	lfd 30, 0(%r25)
	addi %r4, %r4, 36
	lwz %r10, 0(%r4)
	fadd 0, 6,  12	
	lwz %r11, 4(%r4)
	fadd 1, 7,  13	
	lwz %r12, 8(%r4)
	fadd 2, 8,  14	
	lwz %r13, 12(%r4)
	fadd 3, 9,  15	
	lwz %r14, 16(%r4)
	fadd 4, 10, 16	
	lwz %r15, 20(%r4)
	fadd 5, 11, 17
	lwz %r16, 24(%r4)	
	lfd 31, 32(%r25)			
	fadd 0, 0,  18	
	lwz %r17, 28(%r4)
	fadd 1, 1,  19	
	lwz %r18, 32(%r4)
	lfd  12, 0(%r26)		
	fadd 2, 2,  20
	lfd  7, 8(%r25)
	fadd 3, 3,  21	
	lfd  8, 16(%r25)	
	lfd  9, 24(%r25)	
	lfd  16, 32(%r26)	
	lfd  11, 40(%r25)	
	fadd 4, 4,  22
	lfd  13, 8(%r26)	
	fadd 5, 5,  23	
	lfd  18, 0(%r27)
	fadd 0, 0,  24
	lfd  14, 16(%r26)
	lfd  15, 24(%r26)
	fadd 1, 1,  25
	lfd  17, 40(%r26)
	add %r21, %r5, %r11
	lfd  19, 8(%r27)
	fadd 2, 2,  26
	lfd  22, 32(%r27)
	add %r22, %r5, %r12
	lfd  20, 16(%r27)	
	fadd 3, 3,  27
	lfd  24, 0(%r28)	
	add %r23, %r5, %r13
	fadd 4, 4,  28
	lfd  21, 24(%r27)
	add %r24, %r5, %r14
	fadd 5, 5,  29
	lfd  28, 32(%r28)	
	fadd 0, 0,  30
	lfd  23, 40(%r27)			
	fadd 1, 1,  7
	lfd  25, 8(%r28)		
	fadd 2, 2,  8
	lfd  6, 0(%r21)
	fadd 3, 3,  9
	lfd  26, 16(%r28)
	fadd 4, 4,  31
	lfd  27, 24(%r28)	
	fadd 5, 5,  11
	lfd  10, 32(%r21)	
	fadd 0, 0,  12
	lfd  29, 40(%r28)		
	fadd 1, 1,  13
	lfd  7, 8(%r21)
	fadd 2, 2,  14
	lfd  12, 0(%r22)
	fadd 3, 3,  15
	lfd  8, 16(%r21)
	fadd 4, 4,  16
	lfd  9, 24(%r21)
	fadd 5, 5,  17
	lfd  16, 32(%r22)
	fadd 0, 0,  18
	lfd  11, 40(%r21)
	fadd 1, 1,  19
	lfd  13, 8(%r22)
	fadd 2, 2,  20
	lfd  18, 0(%r23)
	fadd 3, 3,  21
	lfd  14, 16(%r22)
	fadd 4, 4,  22
	lfd  15, 24(%r22)
	fadd 5, 5,  23
	lfd  22, 32(%r23)
	fadd 0, 0,  24
	beq lab3				
	lfd  17, 40(%r22)
	fadd 1, 1,  25
	lfd  19, 8(%r23)
	fadd 2, 2,  26
	lfd  24, 0(%r24)
	fadd 3, 3,  27
	lfd  20, 16(%r23)
	fadd 4, 4,  28
	lfd  21, 24(%r23)
	fadd 5, 5,  29
	lfd  28, 32(%r24)
	add %r25, %r5, %r15	
	lfd  23, 40(%r23)
	add %r26, %r5, %r16
	lfd  25, 8(%r24)	
	add %r27, %r5, %r17
	stfd 0, 0(%r20)	
	add %r28, %r5, %r18
	stfd 1, 8(%r20)
	stfd 2, 16(%r20)
	lfd  26, 16(%r24)
	lfd  27, 24(%r24)
	lfd  29, 40(%r24)
	stfd 3, 24(%r20)
	stfd 4, 32(%r20)
	stfd 5, 40(%r20)
	add %r20, %r6, %r10

	bdnz lab1
	b lab0

lab2:
	lfd 30, 0(%r25)
	fadd 0, 6,  12	
	fadd 1, 7,  13	
	fadd 2, 8,  14	
	fadd 3, 9,  15	
	fadd 4, 10, 16	
	fadd 5, 11, 17
	lfd 31, 32(%r25)			
	fadd 0, 0,  18	
	fadd 1, 1,  19	
	lfd  12, 0(%r26)		
	fadd 2, 2,  20
	lfd  7, 8(%r25)
	fadd 3, 3,  21	
	lfd  8, 16(%r25)	
	lfd  9, 24(%r25)	
	lfd  16, 32(%r26)	
	lfd  11, 40(%r25)	
	fadd 4, 4,  22
	lfd  13, 8(%r26)	
	fadd 5, 5,  23	
	lfd  18, 0(%r27)
	fadd 0, 0,  24
	lfd  14, 16(%r26)
	lfd  15, 24(%r26)
	fadd 1, 1,  25
	lfd  17, 40(%r26)
	lfd  19, 8(%r27)
	fadd 2, 2,  26
	lfd  22, 32(%r27)
	lfd  20, 16(%r27)	
	fadd 3, 3,  27
	lfd  24, 0(%r28)	
	fadd 4, 4,  28
	lfd  21, 24(%r27)
	fadd 5, 5,  29
	lfd  28, 32(%r28)	
	fadd 0, 0,  30
	lfd  23, 40(%r27)			
	fadd 1, 1,  7
	lfd  25, 8(%r28)		
	fadd 2, 2,  8
	fadd 3, 3,  9
	lfd  26, 16(%r28)
	fadd 4, 4,  31
	lfd  27, 24(%r28)	
	fadd 5, 5,  11
	fadd 0, 0,  12
	lfd  29, 40(%r28)		
	fadd 1, 1,  13
	fadd 2, 2,  14
	fadd 3, 3,  15
	fadd 4, 4,  16
	fadd 5, 5,  17
	fadd 0, 0,  18
	fadd 1, 1,  19
	fadd 2, 2,  20
	fadd 3, 3,  21
	fadd 4, 4,  22
	fadd 5, 5,  23
	fadd 0, 0,  24
	fadd 1, 1,  25
	fadd 2, 2,  26
	fadd 3, 3,  27
	fadd 4, 4,  28
	fadd 5, 5,  29
	stfd 0, 0(%r20)	
	stfd 1, 8(%r20)
	stfd 2, 16(%r20)
	stfd 3, 24(%r20)
	stfd 4, 32(%r20)
	stfd 5, 40(%r20)

	b lab0	
	
lab3:	
	lfd  17, 40(%r22)
	fadd 1, 1,  25
	lfd  19, 8(%r23)
	fadd 2, 2,  26
	lfd  24, 0(%r24)
	fadd 3, 3,  27
	lfd  20, 16(%r23)
	fadd 4, 4,  28
	lfd  21, 24(%r23)
	fadd 5, 5,  29
	lfd  28, 32(%r24)
	add %r25, %r5, %r15	
	lfd  23, 40(%r23)
	add %r26, %r5, %r16
	lfd  25, 8(%r24)	
	add %r27, %r5, %r17
	stfd 0, 0(%r20)	
	add %r28, %r5, %r18
	stfd 1, 8(%r20)
	stfd 2, 16(%r20)
	lfd  26, 16(%r24)
	lfd  27, 24(%r24)
	lfd  29, 40(%r24)
	stfd 3, 24(%r20)
	stfd 4, 32(%r20)
	stfd 5, 40(%r20)
	add %r20, %r6, %r10

	b lab2

	
lab0:
	
	lwz   %r14,	368(%r1)
	lwz   %r15,	376(%r1)
	lwz   %r16,	384(%r1)
	lwz   %r17,	392(%r1)
	lwz   %r18,	400(%r1)
	lwz   %r19,	408(%r1)
	lwz   %r20,	416(%r1)
	lwz   %r21,	424(%r1)
	lwz   %r22,	432(%r1)
	
	lwz   %r23,	440(%r1)
	lwz   %r24,	448(%r1)
	lwz   %r25,	456(%r1)
	lwz   %r26,	464(%r1)
	lwz   %r27,	472(%r1)
	lwz   %r28,	480(%r1)
	lwz   %r29,	488(%r1)
	lwz   %r30,	496(%r1)
	lwz   %r31,	504(%r1)
	
	lfd   14,	112(%r1)
	lfd   15,	120(%r1)
	lfd   16,	128(%r1)
	lfd   17,	136(%r1)
	lfd   18,	144(%r1)
	lfd   19,	152(%r1)
	lfd   20,	160(%r1)
	lfd   21,	168(%r1)
	lfd   22,	176(%r1)
	lfd   23,	184(%r1)
	lfd   24,	192(%r1)
	lfd   25,	200(%r1)
	lfd   26,	208(%r1)
	lfd   27,	216(%r1)
	lfd   28,	224(%r1)
	lfd   29,	232(%r1)
	lfd   30,	240(%r1)
	lfd   31,	248(%r1)
	
	la   %r1,	864(%r1)
	blr
.Lfe1:
        .size    dirac_sum,.Lfe1-dirac_sum
	.align 2
	.globl dirac_sum_acc
	.type	 dirac_sum_acc,@function
dirac_sum_acc:
	la   %r1,	-864(%r1)
	stw   %r14,	368(%r1)
	stw   %r15,	376(%r1)
	stw   %r16,	384(%r1)
	stw   %r17,	392(%r1)
	stw   %r18,	400(%r1)
	stw   %r19,	408(%r1)
	stw   %r20,	416(%r1)
	stw   %r21,	424(%r1)
	stw   %r22,	432(%r1)
	stw   %r23,	440(%r1)
	stw   %r24,	448(%r1)
	stw   %r25,	456(%r1)
	stw   %r26,	464(%r1)
	stw   %r27,	472(%r1)
	stw   %r28,	480(%r1)
	stw   %r29,	488(%r1)
	stw   %r30,	496(%r1)
	stw   %r31,	504(%r1)
	stfd  14,	112(%r1)
	stfd  15,	120(%r1)
	stfd  16,	128(%r1)
	stfd  17,	136(%r1)
	stfd  18,	144(%r1)
	stfd  19,	152(%r1)
	stfd  20,	160(%r1)
	stfd  21,	168(%r1)
	stfd  22,	176(%r1)
	stfd  23,	184(%r1)
	stfd  24,	192(%r1)
	stfd  25,	200(%r1)
	stfd  26,	208(%r1)
	stfd  27,	216(%r1)
	stfd  28,	224(%r1)
	stfd  29,	232(%r1)
	stfd  30,	240(%r1)
	stfd  31,	248(%r1)
	mr %r0, %r3
	or. %r0, %r0, %r0
	mtctr %r0
	bf gt, .Lab0

	lwz %r10, 0(%r4)
	lwz %r11, 4(%r4)
	lwz %r12, 8(%r4)
	lwz %r13, 12(%r4)
	lwz %r14, 16(%r4)
	lwz %r15, 20(%r4)
	lwz %r16, 24(%r4)
	lwz %r17, 28(%r4)
	lwz %r18, 32(%r4)
	add %r20, %r6, %r10
	add %r21, %r5, %r11
	add %r22, %r5, %r12
	add %r23, %r5, %r13
	add %r24, %r5, %r14
	add %r25, %r5, %r15
	add %r26, %r5, %r16
	add %r27, %r5, %r17
	add %r28, %r5, %r18
	lfd  0, 0(%r20)
	lfd  1, 8(%r20)
	lfd  2, 16(%r20)
	lfd  3, 24(%r20)
	lfd  4, 32(%r20)
	lfd  5, 40(%r20)			
	lfd  6, 0(%r21)
	lfd  7, 8(%r21)
	lfd  8, 16(%r21)
	lfd  9, 24(%r21)
	lfd  10, 32(%r21)
	lfd  11, 40(%r21)
	
	lfd  12, 0(%r22)
	lfd  13, 8(%r22)
	lfd  14, 16(%r22)
	lfd  15, 24(%r22)
	lfd  16, 32(%r22)
	lfd  17, 40(%r22)
	
	lfd  18, 0(%r23)
	lfd  19, 8(%r23)
	lfd  20, 16(%r23)
	lfd  21, 24(%r23)
	lfd  22, 32(%r23)
	lfd  23, 40(%r23)
	
	lfd  24, 0(%r24)
	lfd  25, 8(%r24)
	lfd  26, 16(%r24)
	lfd  27, 24(%r24)
	lfd  28, 32(%r24)
	lfd  29, 40(%r24)

	subic. %r0, %r0, 1
	beq .Lab2
	b .Lab1
.align 4
.Lab1:
	fadd 0, 0,  6		
	subic. %r0, %r0, 1	
	lfd 30, 0(%r25)
	addi %r4, %r4, 36
	fadd 1, 1,  7		
	lwz %r10, 0(%r4)
	fadd 2, 2,  8	
	fadd 3, 3,  9	
	fadd 4, 4,  10	
	fadd 5, 5,  11							
	fadd 0, 0,  12	
	lwz %r11, 4(%r4)
	fadd 1, 1,  13	
	lwz %r12, 8(%r4)
	fadd 2, 2,  14	
	lwz %r13, 12(%r4)
	fadd 3, 3,  15	
	lwz %r14, 16(%r4)
	fadd 4, 4, 16	
	lwz %r15, 20(%r4)
	fadd 5, 5, 17
	lwz %r16, 24(%r4)	
	lfd 31, 32(%r25)			
	fadd 0, 0,  18	
	lwz %r17, 28(%r4)
	fadd 1, 1,  19	
	lwz %r18, 32(%r4)
	lfd  12, 0(%r26)		
	fadd 2, 2,  20
	lfd  7, 8(%r25)
	fadd 3, 3,  21	
	lfd  8, 16(%r25)	
	lfd  9, 24(%r25)	
	lfd  16, 32(%r26)	
	lfd  11, 40(%r25)	
	fadd 4, 4,  22
	lfd  13, 8(%r26)	
	fadd 5, 5,  23	
	lfd  18, 0(%r27)
	fadd 0, 0,  24
	lfd  14, 16(%r26)
	lfd  15, 24(%r26)
	fadd 1, 1,  25
	lfd  17, 40(%r26)
	add %r21, %r5, %r11
	lfd  19, 8(%r27)
	fadd 2, 2,  26
	lfd  22, 32(%r27)
	add %r22, %r5, %r12
	lfd  20, 16(%r27)	
	fadd 3, 3,  27
	lfd  24, 0(%r28)	
	add %r23, %r5, %r13
	fadd 4, 4,  28
	lfd  21, 24(%r27)
	add %r24, %r5, %r14
	fadd 5, 5,  29
	lfd  28, 32(%r28)	
	fadd 0, 0,  30
	lfd  23, 40(%r27)			
	fadd 1, 1,  7
	lfd  25, 8(%r28)		
	fadd 2, 2,  8
	lfd  6, 0(%r21)
	fadd 3, 3,  9
	lfd  26, 16(%r28)
	fadd 4, 4,  31
	lfd  27, 24(%r28)	
	fadd 5, 5,  11
	lfd  10, 32(%r21)	
	fadd 0, 0,  12
	lfd  29, 40(%r28)		
	fadd 1, 1,  13
	lfd  7, 8(%r21)
	fadd 2, 2,  14
	lfd  12, 0(%r22)
	fadd 3, 3,  15
	lfd  8, 16(%r21)
	fadd 4, 4,  16
	lfd  9, 24(%r21)
	fadd 5, 5,  17
	lfd  16, 32(%r22)
	fadd 0, 0,  18
	lfd  11, 40(%r21)
	fadd 1, 1,  19
	lfd  13, 8(%r22)
	fadd 2, 2,  20
	lfd  18, 0(%r23)
	fadd 3, 3,  21
	lfd  14, 16(%r22)
	fadd 4, 4,  22
	lfd  15, 24(%r22)
	fadd 5, 5,  23
	lfd  22, 32(%r23)
	fadd 0, 0,  24
	beq .Lab3				
	lfd  17, 40(%r22)
	fadd 1, 1,  25
	lfd  19, 8(%r23)
	fadd 2, 2,  26
	lfd  24, 0(%r24)
	fadd 3, 3,  27
	lfd  20, 16(%r23)
	fadd 4, 4,  28
	lfd  21, 24(%r23)
	fadd 5, 5,  29
	lfd  28, 32(%r24)
	add %r25, %r5, %r15	
	lfd  23, 40(%r23)
	add %r26, %r5, %r16
	lfd  25, 8(%r24)	
	add %r27, %r5, %r17
	stfd 0, 0(%r20)	
	add %r28, %r5, %r18
	stfd 1, 8(%r20)
	stfd 2, 16(%r20)
	lfd  26, 16(%r24)
	lfd  27, 24(%r24)
	lfd  29, 40(%r24)
	stfd 3, 24(%r20)
	stfd 4, 32(%r20)
	stfd 5, 40(%r20)
	add %r20, %r6, %r10
	lfd  0, 0(%r20) /* can be done better */
	lfd  1, 8(%r20)
	lfd  2, 16(%r20)
	lfd  3, 24(%r20)
	lfd  4, 32(%r20)
	lfd  5, 40(%r20)	
	bdnz .Lab1
	b .Lab0

.Lab2:
	lfd 30, 0(%r25)
	fadd 0, 0,  6	
	fadd 1, 1,  7	
	fadd 2, 2,  8	
	fadd 3, 3,  9	
	fadd 4, 4, 10	
	fadd 5, 5, 11	
	fadd 0, 0,  12	
	fadd 1, 1,  13	
	fadd 2, 2,  14	
	fadd 3, 3,  15	
	fadd 4, 4, 16	
	fadd 5, 5, 17
	lfd 31, 32(%r25)			
	fadd 0, 0,  18	
	fadd 1, 1,  19	
	lfd  12, 0(%r26)		
	fadd 2, 2,  20
	lfd  7, 8(%r25)
	fadd 3, 3,  21	
	lfd  8, 16(%r25)	
	lfd  9, 24(%r25)	
	lfd  16, 32(%r26)	
	lfd  11, 40(%r25)	
	fadd 4, 4,  22
	lfd  13, 8(%r26)	
	fadd 5, 5,  23	
	lfd  18, 0(%r27)
	fadd 0, 0,  24
	lfd  14, 16(%r26)
	lfd  15, 24(%r26)
	fadd 1, 1,  25
	lfd  17, 40(%r26)
	lfd  19, 8(%r27)
	fadd 2, 2,  26
	lfd  22, 32(%r27)
	lfd  20, 16(%r27)	
	fadd 3, 3,  27
	lfd  24, 0(%r28)	
	fadd 4, 4,  28
	lfd  21, 24(%r27)
	fadd 5, 5,  29
	lfd  28, 32(%r28)	
	fadd 0, 0,  30
	lfd  23, 40(%r27)			
	fadd 1, 1,  7
	lfd  25, 8(%r28)		
	fadd 2, 2,  8
	fadd 3, 3,  9
	lfd  26, 16(%r28)
	fadd 4, 4,  31
	lfd  27, 24(%r28)	
	fadd 5, 5,  11
	fadd 0, 0,  12
	lfd  29, 40(%r28)		
	fadd 1, 1,  13
	fadd 2, 2,  14
	fadd 3, 3,  15
	fadd 4, 4,  16
	fadd 5, 5,  17
	fadd 0, 0,  18
	fadd 1, 1,  19
	fadd 2, 2,  20
	fadd 3, 3,  21
	fadd 4, 4,  22
	fadd 5, 5,  23
	fadd 0, 0,  24
	fadd 1, 1,  25
	fadd 2, 2,  26
	fadd 3, 3,  27
	fadd 4, 4,  28
	fadd 5, 5,  29
	stfd 0, 0(%r20)	
	stfd 1, 8(%r20)
	stfd 2, 16(%r20)
	stfd 3, 24(%r20)
	stfd 4, 32(%r20)
	stfd 5, 40(%r20)

	b .Lab0	
	
.Lab3:	
	lfd  17, 40(%r22)
	fadd 1, 1,  25
	lfd  19, 8(%r23)
	fadd 2, 2,  26
	lfd  24, 0(%r24)
	fadd 3, 3,  27
	lfd  20, 16(%r23)
	fadd 4, 4,  28
	lfd  21, 24(%r23)
	fadd 5, 5,  29
	lfd  28, 32(%r24)
	add %r25, %r5, %r15	
	lfd  23, 40(%r23)
	add %r26, %r5, %r16
	lfd  25, 8(%r24)	
	add %r27, %r5, %r17
	stfd 0, 0(%r20)	
	add %r28, %r5, %r18
	stfd 1, 8(%r20)
	stfd 2, 16(%r20)
	lfd  26, 16(%r24)
	lfd  27, 24(%r24)
	lfd  29, 40(%r24)
	stfd 3, 24(%r20)
	stfd 4, 32(%r20)
	stfd 5, 40(%r20)
	add %r20, %r6, %r10

	b .Lab2

	
.Lab0:
	lwz   %r14,	368(%r1)
	lwz   %r15,	376(%r1)
	lwz   %r16,	384(%r1)
	lwz   %r17,	392(%r1)
	lwz   %r18,	400(%r1)
	lwz   %r19,	408(%r1)
	lwz   %r20,	416(%r1)
	lwz   %r21,	424(%r1)
	lwz   %r22,	432(%r1)
	lwz   %r23,	440(%r1)
	lwz   %r24,	448(%r1)
	lwz   %r25,	456(%r1)
	lwz   %r26,	464(%r1)
	lwz   %r27,	472(%r1)
	lwz   %r28,	480(%r1)
	lwz   %r29,	488(%r1)
	lwz   %r30,	496(%r1)
	lwz   %r31,	504(%r1)
	lfd   14,	112(%r1)
	lfd   15,	120(%r1)
	lfd   16,	128(%r1)
	lfd   17,	136(%r1)
	lfd   18,	144(%r1)
	lfd   19,	152(%r1)
	lfd   20,	160(%r1)
	lfd   21,	168(%r1)
	lfd   22,	176(%r1)
	lfd   23,	184(%r1)
	lfd   24,	192(%r1)
	lfd   25,	200(%r1)
	lfd   26,	208(%r1)
	lfd   27,	216(%r1)
	lfd   28,	224(%r1)
	lfd   29,	232(%r1)
	lfd   30,	240(%r1)
	lfd   31,	248(%r1)
	la   %r1,	864(%r1)
	blr
.Lfe2:
        .size    dirac_sum_acc,.Lfe2-dirac_sum_acc
