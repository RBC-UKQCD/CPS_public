        .section        ".text"
        .align 2
        .globl pt_copy_buffer
        .type    pt_copy_buffer,@function
pt_copy_buffer:
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
	li 10, 0
        srawi 0,3,31
        srwi 0,0,30
        add 3,3,0
	dcbt 10, 6 
        srawi. 0,3,2
        mtctr 0
        bclr 4,1
        mr 3,6

.L1:
	lwz 6, 0(3)
	lwz 7, 4(3)
	lwz 8, 8(3)
	slwi 6,6,3
	lwz 9, 12(3)	
	slwi 7,7,3	
	add 10, 4, 6
	slwi 8,8,3	
	add 11, 4, 7
	slwi 9,9,3	
	lfd 0, 0(10)
	add 12, 4, 8
	lfd 4, 32(10)	
	add 13, 4, 9
	lfd 6, 0(11)	
	lfd 1, 8(10)		
	lfd 2, 16(10)
	lfd 3, 24(10)
	lfd 5, 40(10)
	lfd 7, 8(11)
	lfd 8, 16(11)
	stfd 0, 0(5)
	stfd 1, 8(5)
	stfd 2, 16(5)	
	lfd 9, 24(11)
	lfd 10, 32(11)
	lfd 11, 40(11)
	stfd 3, 24(5)
	stfd 4, 32(5)
	stfd 5, 40(5)
	lfd 0, 0(12)
	lfd 1, 8(12)
	lfd 2, 16(12)	
	stfd 6, 48(5)
	stfd 7, 56(5)
	stfd 8, 64(5)
	lfd 3, 24(12)
	lfd 4, 32(12)
	lfd 5, 40(12)	
	stfd 9, 72(5)
	stfd 10, 80(5)
	stfd 11, 88(5)

	lfd 6, 0(13)
	lfd 7, 8(13)
	lfd 8, 16(13)
	stfd 0, 96(5)
	stfd 1, 104(5)
	stfd 2, 112(5)	
	lfd 9, 24(13)
	lfd 10, 32(13)
	lfd 11, 40(13)
	stfd 3, 120(5)
	stfd 4, 128(5)
	stfd 5, 136(5)
	
	stfd 6, 144(5)
	stfd 7, 152(5)
	stfd 8, 160(5)
	stfd 9, 168(5)
	stfd 10, 176(5)
	stfd 11, 184(5)

	addi 5, 5, 192
	addi 3, 3, 16
	
	bdnz .L1 
lab0:
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
