	.section	".text"
	.align 2
	.globl flush_cache_spinor
	.type	 flush_cache_spinor,@function
/* flush 6 * %r3 cache lines + 2 */
flush_cache_spinor:
	la   %r1, -64(%r1)
	mr %r0, %r3
	or. %r0, %r0, %r0
	mtctr %r0
	bf gt, lab0
	
        li 7, 0
	li 8, 64
	li 9, 96
	li 10, 128
	dcbf %r7, %r4
	li 11, 160
	li 12, 192
	li 7, 32
	
lab1:
	dcbf %r7, %r4
	dcbf %r8, %r4
	dcbf %r9, %r4
	addi %r7, %r7, 192	
	dcbf %r10, %r4
	addi %r8, %r8, 192	
	dcbf %r11, %r4
	addi %r9, %r9, 192	
	dcbf %r12, %r4
	addi %r10, %r10, 192
	addi %r11, %r11, 192
	addi %r12, %r12, 192
	bdnz lab1
	dcbf %r7, %r4
lab0:
	la   %r1,	64(%r1)
	blr
.Lfe1:
        .size    flush_cache_spinor,.Lfe1-flush_cache_spinor
