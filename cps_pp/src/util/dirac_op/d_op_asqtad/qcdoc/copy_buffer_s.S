        .section        ".text"
        .align 2
        .globl copy_buffer_s
        .type    copy_buffer,@function
copy_buffer_s:
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
	slwi 6,6,2
	lwz 9, 12(3)	
	slwi 7,7,2	
	add 10, 4, 6
	slwi 8,8,2	
	add 11, 4, 7
	slwi 9,9,2	
	lfs 0, 0(10)
	add 12, 4, 8
	lfs 4, 16(10)	
	add 13, 4, 9
	lfs 6, 0(11)	
	lfs 1, 4(10)		
	lfs 2, 8(10)
	lfs 3, 12(10)
	lfs 5, 20(10)
	lfs 7, 4(11)
	lfs 8, 8(11)
	stfs 0, 0(5)
	stfs 1, 4(5)
	stfs 2, 8(5)	
	lfs 9, 12(11)
	lfs 10, 16(11)
	lfs 11, 20(11)
	stfs 3, 12(5)
	stfs 4, 16(5)
	stfs 5, 20(5)
	lfs 0, 0(12)
	lfs 1, 4(12)
	lfs 2, 8(12)	
	stfs 6, 24(5)
	stfs 7, 28(5)
	stfs 8, 32(5)
	lfs 3, 12(12)
	lfs 4, 16(12)
	lfs 5, 20(12)	
	stfs 9, 36(5)
	stfs 10, 40(5)
	stfs 11, 44(5)

	lfs 6, 0(13)
	lfs 7, 4(13)
	lfs 8, 8(13)
	stfs 0, 48(5)
	stfs 1, 52(5)
	stfs 2, 56(5)	
	lfs 9, 12(13)
	lfs 10, 16(13)
	lfs 11, 20(13)
	stfs 3, 60(5)
	stfs 4, 64(5)
	stfs 5, 68(5)
	
	stfs 6, 72(5)
	stfs 7, 76(5)
	stfs 8, 80(5)
	stfs 9, 84(5)
	stfs 10, 88(5)
	stfs 11, 92(5)

	addi 5, 5, 96
	addi 3, 3, 16
	
	bdnz .L1 
	blr
