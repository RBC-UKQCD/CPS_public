        .section        ".text"
        .align 2
        .globl copy_buffer
        .type    copy_buffer,@function
copy_buffer:
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
	blr
