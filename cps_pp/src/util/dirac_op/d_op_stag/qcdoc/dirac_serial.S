        .section        ".text"
        .align 2
        .globl dirac_cmv
        .type    dirac_cmv,@function
dirac_cmv:
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
	li 15, 288	
	stw   %r24,	448(%r1)
	li 16, 16			
	stw   %r25,	456(%r1)
	li 18, 96
	lwz %r12, 0(%r4)
	li 19, 0		
	lwz %r14, 4(%r4)
	li 20, 32	
	lwz %r11, 8(%r4)
	li 21, 64	
	lwz %r13, 12(%r4)
	li 22, 96	
	stw   %r26,	464(%r1)
	li 23, 128
	stw   %r27,	472(%r1)
	li 24, 160
	stw   %r28,	480(%r1)
	li 25, 192	
	stw   %r29,	488(%r1)
	li 26, 224	
	stw   %r30,	496(%r1)
	li 27, 256	
	stw   %r31,	504(%r1)
	add %r12, %r12, %r6 /*fp0 pointer*/
	add %r11, %r11, %r6 /*fp0 pointer*/	
	
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

	
	
        srwi 0, 3, 31
	dcbt %r19, %r5	
        add 0, 0, 3
	dcbt %r20, %r5	
        srawi. 3, 0, 1
	dcbt %r21, %r5	
	bf gt,   lab0
        dcbt %r19, %r12
	add %r14, %r14, %r7 /*fp1 pointer*/
	dcbt %r20, %r12
	add %r13, %r13, %r7 /*fp1 pointer*/
	dcbt %r19, %r11
	add %r29, %r15, %r5
	dcbt %r20, %r11	
	dcbt %r22, %r5

	dcbt %r23, %r5	

	dcbt %r24, %r5	

	dcbt %r25, %r5

	dcbt %r26, %r5

	dcbt %r27, %r5

        lfd 12, 0(%r5)         /* uu */
	lfd 13, 8(%r5)         /* uu+1 */	
	lfd 18, 48(%r5)        /* uu+6 */
        lfd 19, 56(%r5)        /* uu+7 */			
	lfd 0,	0(%r12)		/* fp0 */
	lfd 1,	8(%r12)		/* fp0+1 */
	lfd 2,	16(%r12)	/* fp0+2 */
	lfd 3,	24(%r12)	/* fp0+3 */			
        lfd 24, 96(%r5)        /* uu+12 */
	lfd 25, 104(%r5)       /* uu+13 */
	subic. 3, 3, 1
	bf gt, lab3
	lwz %r30, 16(%r4)
	lwz %r9, 20(%r4)
	lwz %r31, 24(%r4)
	lwz %r10, 28(%r4)
	add %r30, %r30, %r6 /*fp0 pointer*/
	add %r9 , %r9 , %r7 /*fp1 pointer*/	
	add %r31, %r31, %r6 /*fp0 pointer*/
	add %r10, %r10, %r7 /*fp1 pointer*/	
	subic. 3, 3, 1
	bf gt, lab2
        mtctr 3
	b lab1
.align 4
lab1:
	dcbt %r19, %r30	
	fmul 6,  0, 12   /*  *(fp1) += *fp0 * *(uu) */
	fmul 7,  1, 12   /*  *(fp1+1) += *(fp0+1) * *(uu) */
	dcbt %r20, %r30				
	addi %r4, %r4, 16	/*increment offsets pointer*/		
	fmul 8,  0, 18   /*  *(fp1+2) += *fp0 * *(uu+6) */
        lfd 14, 16(%r5)        /* uu+2 */				
	fmul 9,  1, 18   /*  *(fp1+3) += *(fp0+1) * *(uu+6) */
        lfd 20, 64(%r5)        /* uu+8 */		
	fmul 10, 0, 24  /*  *(fp1+4) += *fp0 * *(uu+12) */
        lfd 26, 112(%r5)       /* uu+14 */		
	fmul 11, 1, 24  /*  *(fp1+5) += *(fp0+1) * *(uu+12) */
        lfd 15, 24(%r5)        /* uu+3 */	
	fnmsub 6,  1, 13, 6   /* *(fp1) -= *(fp0+1) * *(uu+1) */

	fmadd 7,  0, 13, 7   /*  *(fp1+1) += *fp0 * *(uu+1) */
	dcbt %r19, %r31  
	fnmsub 8,  1, 19, 8   /*  *(fp1+2) -= *(fp0+1) * *(uu+7) */
        lfd 21, 72(%r5)        /* uu+9 */			
	fmadd 9,  0, 19, 9   /*  *(fp1+3) += *fp0 * *(uu+7) */
        lfd 27, 120(%r5)       /* uu+15 */		
	fnmsub 10, 1, 25, 10  /*  *(fp1+4) -= *(fp0+1) * *(uu+13) */
	lfd 4,	32(%r12)	/* fp0+4 */				
	fmadd 11, 0, 25, 11  /*  *(fp1+5) += *fp0 * *(uu+13) */
	lfd 5,	40(%r12)	/* fp0+5 */
/* now switch fpr's 6-12, 7-18, 8-24, 9-13, 10-19, 11-25 */	
	fmadd 12,  2, 14, 6   /*  *(fp1) += *(fp0+2) * *(uu+2) */

	fmadd 18,  3, 14, 7   /*  *(fp1+1) += *(fp0+3) * *(uu+2) */
	dcbt %r20, %r31  
	mr %r12, %r30	
	fmadd 24,  2, 20, 8   /*  *(fp1+2) +=*(fp0+2) * *(uu+8) */
        lfd 16, 32(%r5)        /* uu+4 */		
	fmadd 13,  3, 20, 9   /*  *(fp1+3) += *(fp0+3) * *(uu+8) */
        lfd 22, 80(%r5)        /* uu+10 */		
	fmadd 19, 2, 26, 10  /*  *(fp1+4) +=*(fp0+2) * *(uu+14) */
        lfd 28, 128(%r5)       /* uu+16 */			
	fmadd 25, 3, 26, 11  /*  *(fp1+5) += *(fp0+3) * *(uu+14) */
        lfd 17, 40(%r5)        /* uu+5 */		
/* switch 12-30, 18-31 */	
	fnmsub 30,  3, 15, 12   /*  *(fp1) -= *(fp0+3)* *(uu+3) */

	fmadd 31,  2, 15, 18   /*  *(fp1+1) += *(fp0+2) * *(uu+3) */
        dcbt %r19, %r29				
	fnmsub 24,  3, 21, 24   /*  *(fp1+2) -= *(fp0+3)* *(uu+9) */
        lfd 23, 88(%r5)        /* uu+11 */		
	fmadd 13,  2, 21, 13   /*  *(fp1+3) += *(fp0+2) * *(uu+9) */
	lfd 29, 136(%r5)       /* uu+17 */
/* done loading for current iteration */
       	fnmsub 19, 3, 27, 19  /* *(fp1+4) -= *(fp0+3)* *(uu+15) */
/* start loading for next iteration */
	lfd 0,	0(%r11)		/* fp0 */ 
	fmadd 25, 2, 27, 25  /*  *(fp1+5) += *(fp0+2) * *(uu+15) */
	lfd 1,	8(%r11)	/* fp0+1 */ 
/* now switch fpr's 30-14, 31-20, 24-26, 13-15, 19-21, 25-27 */
	fmadd 14,  4, 16, 30   /*  *(fp1) += *(fp0+4) * *(uu+4) */

	fmadd 20,  5, 16, 31   /*  *(fp1+1) += *(fp0+5) * *(uu+4) */
	dcbt %r20, %r29			
	fmadd 26,  4, 22, 24   /* *(fp1+2) += *(fp0+4) * *(uu+10) */
	lfd 12, 144(%r5)    /* uu */	 		
	fmadd 15,  5, 22, 13   /* *(fp1+3) += *(fp0+5) * *(uu+10) */
	lfd 18, 192(%r5)        /* uu+6 */			
	fmadd 21, 4, 28, 19   /*  *(fp1+4) += *(fp0+4) * *(uu+16) */
	lfd 24, 240(%r5)        /* uu+12 */		
	fmadd 27, 5, 28, 25   /*  *(fp1+5) += *(fp0+5) * *(uu+16) */
	lfd 13, 152(%r5)         /* uu+1 */	
	fnmsub 14,  5, 17, 14   /*  *(fp1) -= *(fp0+5) * *(uu+5) */

	fmadd 20,  4, 17, 20   /*  *(fp1+1) += *(fp0+4) * *(uu+5) */
	dcbt %r21, %r29				
	fnmsub 26,  5, 23, 26  /* *(fp1+2) -= *(fp0+5) * *(uu+11) */
	lfd 19, 200(%r5)        /* uu+7 */	
	fmadd 15,  4, 23, 15   /* *(fp1+3) += *(fp0+4) * *(uu+11) */
        lfd 25, 248(%r5)       /* uu+13 */		
	fnmsub 21, 5, 29, 21  /* *(fp1+4) -= *(fp0+5) * *(uu+17) */
	lfd 2,	16(%r11)	/* fp0+2 */		
	fmadd 27, 4, 29, 27  /*  *(fp1+5) += *(fp0+4) * *(uu+17) */
	stfd  14,  0(%r14)   /* fp1 */	
	fmul 6,  0, 12   /*  *(fp1) += *fp0 * *(uu) */
	stfd  20,  8(%r14)   /* fp1+1 */	
	fmul 7,  1, 12   /*  *(fp1+1) += *(fp0+1) * *(uu) */
	stfd  26,  16(%r14)  /* fp1+2 */	
	fmul 8,  0, 18   /*  *(fp1+2) += *fp0 * *(uu+6) */
	dcbt %r22, %r29				
	fmul 9,  1, 18   /*  *(fp1+3) += *(fp0+1) * *(uu+6) */	

	fmul 10, 0, 24  /*  *(fp1+4) += *fp0 * *(uu+12) */
	stfd  15,  24(%r14)  /* fp1+3 */
	stfd  21, 32(%r14)  /* fp1+4 */
	stfd  27, 40(%r14)  /* fp1+5 */
	lwz %r30, 16(%r4)  /*cache miss half of the time*/
	lfd 3,	24(%r11)	/* fp0+3 */			
        lfd 14, 160(%r5)        /* uu+2 */
	mr %r14, %r9	
	fmul 11, 1, 24  /*  *(fp1+5) += *(fp0+1) * *(uu+12) */
	
	fnmsub 6,  1, 13, 6   /* *(fp1) -= *(fp0+1) * *(uu+1) */	
	dcbt %r23, %r29
	fmadd 7,  0, 13, 7   /*  *(fp1+1) += *fp0 * *(uu+1) */
        lfd 20, 208(%r5)        /* uu+8 */	
	fnmsub 8,  1, 19, 8   /*  *(fp1+2) -= *(fp0+1) * *(uu+7) */
        lfd 26, 256(%r5)       /* uu+14 */	
	fmadd 9,  0, 19, 9   /*  *(fp1+3) += *fp0 * *(uu+7) */
        lfd 15, 168(%r5)        /* uu+3 */
	fnmsub 10, 1, 25, 10  /*  *(fp1+4) -= *(fp0+1) * *(uu+13) */
        lfd 21, 216(%r5)        /* uu+9 */
	fmadd 11, 0, 25, 11  /*  *(fp1+5) += *fp0 * *(uu+13) */

/* now switch fpr's 6-12, 7-18, 8-24, 9-13, 10-19, 11-25 */	
	fmadd 12,  2, 14, 6   /*  *(fp1) += *(fp0+2) * *(uu+2) */
	dcbt %r24, %r29					
	fmadd 18,  3, 14, 7   /*  *(fp1+1) += *(fp0+3) * *(uu+2) */
        lfd 27, 264(%r5)       /* uu+15 */	
	fmadd 24,  2, 20, 8   /*  *(fp1+2) +=*(fp0+2) * *(uu+8) */
        lfd 16, 176(%r5)        /* uu+4 */	
	fmadd 13,  3, 20, 9   /*  *(fp1+3) += *(fp0+3) * *(uu+8) */
	lfd 5,	40(%r11)	/* fp0+5 */	
	fmadd 19, 2, 26, 10  /*  *(fp1+4) +=*(fp0+2) * *(uu+14) */
	lfd 4,	32(%r11)	/* fp0+4 */
	fmadd 25, 3, 26, 11  /*  *(fp1+5) += *(fp0+3) * *(uu+14) */
	mr %r8, %r10 /* save fp1 pointer*/
	mr %r11, %r31 	
	fnmsub 12,  3, 15, 12   /*  *(fp1) -= *(fp0+3)* *(uu+3) */
	dcbt %r25, %r29
	fmadd 18,  2, 15, 18   /*  *(fp1+1) += *(fp0+2) * *(uu+3) */
        lfd 22, 224(%r5)        /* uu+10 */	
	fnmsub 24,  3, 21, 24   /*  *(fp1+2) -= *(fp0+3)* *(uu+9) */
        lfd 28, 272(%r5)       /* uu+16 */
	fmadd 13,  2, 21, 13   /*  *(fp1+3) += *(fp0+2) * *(uu+9) */
        lfd 17, 184(%r5)        /* uu+5 */		
       	fnmsub 19, 3, 27, 19  /* *(fp1+4) -= *(fp0+3)* *(uu+15) */
        lfd 23, 232(%r5)        /* uu+11 */		
	fmadd 25, 2, 27, 25  /*  *(fp1+5) += *(fp0+2) * *(uu+15) */
        lwz %r9, 20(%r4)
	lwz %r31, 24(%r4)
	lwz %r10, 28(%r4)
/* now switch fpr's 12-14, 18-20, 24-26, 13-15, 19-21, 25-27 */
	fmadd 14,  4, 16, 12   /*  *(fp1) += *(fp0+4) * *(uu+4) */
	dcbt %r26, %r29					
	fmadd 20,  5, 16, 18   /*  *(fp1+1) += *(fp0+5) * *(uu+4) */
        lfd 29, 280(%r5)       /* uu+17 */
/* done loading for current iteration */
	fmadd 26,  4, 22, 24   /* *(fp1+2) += *(fp0+4) * *(uu+10) */
/* start loading for next iteration */
	lfd 0,	0(%r12)		/* fp0 */	
	fmadd 15,  5, 22, 13   /* *(fp1+3) += *(fp0+5) * *(uu+10) */
	lfd 1,	8(%r12)		/* fp0+1 */
	fmadd 21, 4, 28, 19   /*  *(fp1+4) += *(fp0+4) * *(uu+16) */
	lfd 2,	16(%r12)	/* fp0+2 */
	fmadd 27, 5, 28, 25   /*  *(fp1+5) += *(fp0+5) * *(uu+16) */
	fnmsub 14,  5, 17, 14   /*  *(fp1) -= *(fp0+5) * *(uu+5) */
	dcbt %r27, %r29
	fmadd 20,  4, 17, 20   /*  *(fp1+1) += *(fp0+4) * *(uu+5) */
	lfd 3,	24(%r12)	/* fp0+3 */	
	fnmsub 26,  5, 23, 26  /* *(fp1+2) -= *(fp0+5) * *(uu+11) */
	lfd 12, 288(%r5)      /* uu */	
	fmadd 15,  4, 23, 15   /* *(fp1+3) += *(fp0+4) * *(uu+11) */
	lfd 18, 336(%r5)        /* uu+6 */			
	fnmsub 21, 5, 29, 21  /* *(fp1+4) -= *(fp0+5) * *(uu+17) */
	lfd 24, 384(%r5)        /* uu+12 */	
	fmadd 27, 4, 29, 27  /*  *(fp1+5) += *(fp0+4) * *(uu+17) */
	stfd  14,  0(%r13)
	add %r9, %r9, %r7
	stfd  20,  8(%r13)	
	add %r10, %r10, %r7
	stfd  26,  16(%r13)	
	add %r30, %r30, %r6 /*fp0 pointer*/
	lfd 13, 296(%r5)         /* uu+1 */
	add %r31, %r31, %r6 /*fp0+6 pointer*/	
	lfd 19, 344(%r5)        /* uu+7 */
        lfd 25, 392(%r5)       /* uu+13 */	
	add  %r5, %r15, %r5
	stfd  15,  24(%r13)
	add  %r29, %r15, %r29
	stfd  21, 32(%r13)
	stfd  27, 40(%r13)
	mr %r13, %r8

	bdnz	lab1
        b lab2
.align 4
lab2:	
	fmul 6,  0, 12   /*  *(fp1) += *fp0 * *(uu) */
	fmul 7,  1, 12   /*  *(fp1+1) += *(fp0+1) * *(uu) */
	dcbt %r19, %r30
	fmul 8,  0, 18   /*  *(fp1+2) += *fp0 * *(uu+6) */
        lfd 14, 16(%r5)        /* uu+2 */				
	fmul 9,  1, 18   /*  *(fp1+3) += *(fp0+1) * *(uu+6) */
        lfd 20, 64(%r5)        /* uu+8 */		
	fmul 10, 0, 24  /*  *(fp1+4) += *fp0 * *(uu+12) */
        lfd 26, 112(%r5)       /* uu+14 */		
	fmul 11, 1, 24  /*  *(fp1+5) += *(fp0+1) * *(uu+12) */
        lfd 15, 24(%r5)        /* uu+3 */	
	fnmsub 6,  1, 13, 6   /* *(fp1) -= *(fp0+1) * *(uu+1) */

	fmadd 7,  0, 13, 7   /*  *(fp1+1) += *fp0 * *(uu+1) */
	dcbt %r20, %r30				
	fnmsub 8,  1, 19, 8   /*  *(fp1+2) -= *(fp0+1) * *(uu+7) */
        lfd 21, 72(%r5)        /* uu+9 */			
	fmadd 9,  0, 19, 9   /*  *(fp1+3) += *fp0 * *(uu+7) */
        lfd 27, 120(%r5)       /* uu+15 */		
	fnmsub 10, 1, 25, 10  /*  *(fp1+4) -= *(fp0+1) * *(uu+13) */
	lfd 4,	32(%r12)	/* fp0+4 */				
	fmadd 11, 0, 25, 11  /*  *(fp1+5) += *fp0 * *(uu+13) */
	lfd 5,	40(%r12)	/* fp0+5 */
/* now switch fpr's 6-12, 7-18, 8-24, 9-13, 10-19, 11-25 */	
	fmadd 12,  2, 14, 6   /*  *(fp1) += *(fp0+2) * *(uu+2) */

	fmadd 18,  3, 14, 7   /*  *(fp1+1) += *(fp0+3) * *(uu+2) */
	dcbt %r16, %r31  /*should work most of the time */
	mr %r12, %r30	
	fmadd 24,  2, 20, 8   /*  *(fp1+2) +=*(fp0+2) * *(uu+8) */
        lfd 16, 32(%r5)        /* uu+4 */		
	fmadd 13,  3, 20, 9   /*  *(fp1+3) += *(fp0+3) * *(uu+8) */
        lfd 22, 80(%r5)        /* uu+10 */		
	fmadd 19, 2, 26, 10  /*  *(fp1+4) +=*(fp0+2) * *(uu+14) */
        lfd 28, 128(%r5)       /* uu+16 */			
	fmadd 25, 3, 26, 11  /*  *(fp1+5) += *(fp0+3) * *(uu+14) */
        lfd 17, 40(%r5)        /* uu+5 */		
/* switch 12-30, 18-31 */	
	fnmsub 30,  3, 15, 12   /*  *(fp1) -= *(fp0+3)* *(uu+3) */

	fmadd 31,  2, 15, 18   /*  *(fp1+1) += *(fp0+2) * *(uu+3) */
        dcbt %r19, %r29				
	fnmsub 24,  3, 21, 24   /*  *(fp1+2) -= *(fp0+3)* *(uu+9) */
        lfd 23, 88(%r5)        /* uu+11 */		
	fmadd 13,  2, 21, 13   /*  *(fp1+3) += *(fp0+2) * *(uu+9) */
	lfd 29, 136(%r5)       /* uu+17 */
/* done loading for current iteration */
       	fnmsub 19, 3, 27, 19  /* *(fp1+4) -= *(fp0+3)* *(uu+15) */
/* start loading for next iteration */
	lfd 0,	0(%r11)		/* fp0 */ /*sometimes a cache miss*/
	fmadd 25, 2, 27, 25  /*  *(fp1+5) += *(fp0+2) * *(uu+15) */
	lfd 1,	8(%r11)	/* fp0+1 */ /*sometimes a cache miss*/
/* now switch fpr's 30-14, 31-20, 24-26, 13-15, 19-21, 25-27 */
	fmadd 14,  4, 16, 30   /*  *(fp1) += *(fp0+4) * *(uu+4) */

	fmadd 20,  5, 16, 31   /*  *(fp1+1) += *(fp0+5) * *(uu+4) */
	dcbt %r20, %r29			
	fmadd 26,  4, 22, 24   /* *(fp1+2) += *(fp0+4) * *(uu+10) */
	lfd 12, 144(%r5)    /* uu */	 		
	fmadd 15,  5, 22, 13   /* *(fp1+3) += *(fp0+5) * *(uu+10) */
	lfd 18, 192(%r5)        /* uu+6 */			
	fmadd 21, 4, 28, 19   /*  *(fp1+4) += *(fp0+4) * *(uu+16) */
	lfd 24, 240(%r5)        /* uu+12 */		
	fmadd 27, 5, 28, 25   /*  *(fp1+5) += *(fp0+5) * *(uu+16) */
	lfd 13, 152(%r5)         /* uu+1 */	
	fnmsub 14,  5, 17, 14   /*  *(fp1) -= *(fp0+5) * *(uu+5) */

	fmadd 20,  4, 17, 20   /*  *(fp1+1) += *(fp0+4) * *(uu+5) */
	dcbt %r21, %r29				
	fnmsub 26,  5, 23, 26  /* *(fp1+2) -= *(fp0+5) * *(uu+11) */
	lfd 19, 200(%r5)        /* uu+7 */	
	fmadd 15,  4, 23, 15   /* *(fp1+3) += *(fp0+4) * *(uu+11) */
        lfd 25, 248(%r5)       /* uu+13 */		
	fnmsub 21, 5, 29, 21  /* *(fp1+4) -= *(fp0+5) * *(uu+17) */
	lfd 2,	16(%r11)	/* fp0+2 */		
	fmadd 27, 4, 29, 27  /*  *(fp1+5) += *(fp0+4) * *(uu+17) */
	stfd  14,  0(%r14)   /* fp1 */	
	fmul 6,  0, 12   /*  *(fp1) += *fp0 * *(uu) */
	stfd  20,  8(%r14)   /* fp1+1 */	
	fmul 7,  1, 12   /*  *(fp1+1) += *(fp0+1) * *(uu) */
	stfd  26,  16(%r14)  /* fp1+2 */	
	fmul 8,  0, 18   /*  *(fp1+2) += *fp0 * *(uu+6) */
	dcbt %r22, %r29				
	fmul 9,  1, 18   /*  *(fp1+3) += *(fp0+1) * *(uu+6) */	

	fmul 10, 0, 24  /*  *(fp1+4) += *fp0 * *(uu+12) */
	stfd  15,  24(%r14)  /* fp1+3 */
	stfd  21, 32(%r14)  /* fp1+4 */
	stfd  27, 40(%r14)  /* fp1+5 */
	lfd 3,	24(%r11)	/* fp0+3 */			
        lfd 14, 160(%r5)        /* uu+2 */
	mr %r14, %r9	
	fmul 11, 1, 24  /*  *(fp1+5) += *(fp0+1) * *(uu+12) */
	
	fnmsub 6,  1, 13, 6   /* *(fp1) -= *(fp0+1) * *(uu+1) */	
	dcbt %r23, %r29
	fmadd 7,  0, 13, 7   /*  *(fp1+1) += *fp0 * *(uu+1) */
        lfd 20, 208(%r5)        /* uu+8 */	
	fnmsub 8,  1, 19, 8   /*  *(fp1+2) -= *(fp0+1) * *(uu+7) */
        lfd 26, 256(%r5)       /* uu+14 */	
	fmadd 9,  0, 19, 9   /*  *(fp1+3) += *fp0 * *(uu+7) */
        lfd 15, 168(%r5)        /* uu+3 */
	fnmsub 10, 1, 25, 10  /*  *(fp1+4) -= *(fp0+1) * *(uu+13) */
        lfd 21, 216(%r5)        /* uu+9 */
	fmadd 11, 0, 25, 11  /*  *(fp1+5) += *fp0 * *(uu+13) */

/* now switch fpr's 6-12, 7-18, 8-24, 9-13, 10-19, 11-25 */	
	fmadd 12,  2, 14, 6   /*  *(fp1) += *(fp0+2) * *(uu+2) */
	dcbt %r24, %r29					
	fmadd 18,  3, 14, 7   /*  *(fp1+1) += *(fp0+3) * *(uu+2) */
        lfd 27, 264(%r5)       /* uu+15 */	
	fmadd 24,  2, 20, 8   /*  *(fp1+2) +=*(fp0+2) * *(uu+8) */
        lfd 16, 176(%r5)        /* uu+4 */	
	fmadd 13,  3, 20, 9   /*  *(fp1+3) += *(fp0+3) * *(uu+8) */
	lfd 5,	40(%r11)	/* fp0+5 */	
	fmadd 19, 2, 26, 10  /*  *(fp1+4) +=*(fp0+2) * *(uu+14) */
	lfd 4,	32(%r11)	/* fp0+4 */
	fmadd 25, 3, 26, 11  /*  *(fp1+5) += *(fp0+3) * *(uu+14) */
	mr %r11, %r31 	
	fnmsub 12,  3, 15, 12   /*  *(fp1) -= *(fp0+3)* *(uu+3) */
	dcbt %r25, %r29
	fmadd 18,  2, 15, 18   /*  *(fp1+1) += *(fp0+2) * *(uu+3) */
        lfd 22, 224(%r5)        /* uu+10 */	
	fnmsub 24,  3, 21, 24   /*  *(fp1+2) -= *(fp0+3)* *(uu+9) */
        lfd 28, 272(%r5)       /* uu+16 */
	fmadd 13,  2, 21, 13   /*  *(fp1+3) += *(fp0+2) * *(uu+9) */
        lfd 17, 184(%r5)        /* uu+5 */		
       	fnmsub 19, 3, 27, 19  /* *(fp1+4) -= *(fp0+3)* *(uu+15) */
        lfd 23, 232(%r5)        /* uu+11 */		
	fmadd 25, 2, 27, 25  /*  *(fp1+5) += *(fp0+2) * *(uu+15) */
/* now switch fpr's 12-14, 18-20, 24-26, 13-15, 19-21, 25-27 */
	fmadd 14,  4, 16, 12   /*  *(fp1) += *(fp0+4) * *(uu+4) */
	dcbt %r26, %r29					
	fmadd 20,  5, 16, 18   /*  *(fp1+1) += *(fp0+5) * *(uu+4) */
        lfd 29, 280(%r5)       /* uu+17 */
/* done loading for current iteration */
	fmadd 26,  4, 22, 24   /* *(fp1+2) += *(fp0+4) * *(uu+10) */
/* start loading for next iteration */
	lfd 0,	0(%r12)		/* fp0 */	
	fmadd 15,  5, 22, 13   /* *(fp1+3) += *(fp0+5) * *(uu+10) */
	lfd 1,	8(%r12)		/* fp0+1 */
	fmadd 21, 4, 28, 19   /*  *(fp1+4) += *(fp0+4) * *(uu+16) */
	lfd 2,	16(%r12)	/* fp0+2 */
	fmadd 27, 5, 28, 25   /*  *(fp1+5) += *(fp0+5) * *(uu+16) */
	fnmsub 14,  5, 17, 14   /*  *(fp1) -= *(fp0+5) * *(uu+5) */
	dcbt %r27, %r29
	fmadd 20,  4, 17, 20   /*  *(fp1+1) += *(fp0+4) * *(uu+5) */
	lfd 3,	24(%r12)	/* fp0+3 */	
	fnmsub 26,  5, 23, 26  /* *(fp1+2) -= *(fp0+5) * *(uu+11) */
	lfd 12, 288(%r5)      /* uu */	
	fmadd 15,  4, 23, 15   /* *(fp1+3) += *(fp0+4) * *(uu+11) */
	lfd 18, 336(%r5)        /* uu+6 */			
	fnmsub 21, 5, 29, 21  /* *(fp1+4) -= *(fp0+5) * *(uu+17) */
	lfd 24, 384(%r5)        /* uu+12 */	
	fmadd 27, 4, 29, 27  /*  *(fp1+5) += *(fp0+4) * *(uu+17) */
	stfd  14,  0(%r13)
	stfd  20,  8(%r13)	
	stfd  26,  16(%r13)	
	lfd 13, 296(%r5)         /* uu+1 */
	lfd 19, 344(%r5)        /* uu+7 */
        lfd 25, 392(%r5)       /* uu+13 */	
	add  %r5, %r15, %r5
	stfd  15,  24(%r13)
	add  %r29, %r15, %r29
	stfd  21, 32(%r13)
	stfd  27, 40(%r13)
	or %r13, %r10, %r10
        b lab3
	
/* last 2 U*a multiplies */
.align 4
lab3:		
	fmul 6,  0, 12   /*  *(fp1) += *fp0 * *(uu) */
	fmul 7,  1, 12   /*  *(fp1+1) += *(fp0+1) * *(uu) */
	fmul 8,  0, 18   /*  *(fp1+2) += *fp0 * *(uu+6) */
        lfd 14, 16(%r5)        /* uu+2 */				
	fmul 9,  1, 18   /*  *(fp1+3) += *(fp0+1) * *(uu+6) */
        lfd 20, 64(%r5)        /* uu+8 */		
	fmul 10, 0, 24  /*  *(fp1+4) += *fp0 * *(uu+12) */
        lfd 26, 112(%r5)       /* uu+14 */		
	fmul 11, 1, 24  /*  *(fp1+5) += *(fp0+1) * *(uu+12) */
        lfd 15, 24(%r5)        /* uu+3 */	
	fnmsub 6,  1, 13, 6   /* *(fp1) -= *(fp0+1) * *(uu+1) */

	fmadd 7,  0, 13, 7   /*  *(fp1+1) += *fp0 * *(uu+1) */

	fnmsub 8,  1, 19, 8   /*  *(fp1+2) -= *(fp0+1) * *(uu+7) */
        lfd 21, 72(%r5)        /* uu+9 */			
	fmadd 9,  0, 19, 9   /*  *(fp1+3) += *fp0 * *(uu+7) */
        lfd 27, 120(%r5)       /* uu+15 */		
	fnmsub 10, 1, 25, 10  /*  *(fp1+4) -= *(fp0+1) * *(uu+13) */
	lfd 4,	32(%r12)	/* fp0+4 */				
	fmadd 11, 0, 25, 11  /*  *(fp1+5) += *fp0 * *(uu+13) */
	lfd 5,	40(%r12)	/* fp0+5 */
/* now switch fpr's 6-12, 7-18, 8-24, 9-13, 10-19, 11-25 */	
	fmadd 12,  2, 14, 6   /*  *(fp1) += *(fp0+2) * *(uu+2) */

	fmadd 18,  3, 14, 7   /*  *(fp1+1) += *(fp0+3) * *(uu+2) */

	fmadd 24,  2, 20, 8   /*  *(fp1+2) +=*(fp0+2) * *(uu+8) */
        lfd 16, 32(%r5)        /* uu+4 */		
	fmadd 13,  3, 20, 9   /*  *(fp1+3) += *(fp0+3) * *(uu+8) */
        lfd 22, 80(%r5)        /* uu+10 */		
	fmadd 19, 2, 26, 10  /*  *(fp1+4) +=*(fp0+2) * *(uu+14) */
        lfd 28, 128(%r5)       /* uu+16 */			
	fmadd 25, 3, 26, 11  /*  *(fp1+5) += *(fp0+3) * *(uu+14) */
        lfd 17, 40(%r5)        /* uu+5 */		
/* switch 12-30, 18-31 */	
	fnmsub 30,  3, 15, 12   /*  *(fp1) -= *(fp0+3)* *(uu+3) */

	fmadd 31,  2, 15, 18   /*  *(fp1+1) += *(fp0+2) * *(uu+3) */

	fnmsub 24,  3, 21, 24   /*  *(fp1+2) -= *(fp0+3)* *(uu+9) */
        lfd 23, 88(%r5)        /* uu+11 */		
	fmadd 13,  2, 21, 13   /*  *(fp1+3) += *(fp0+2) * *(uu+9) */
	lfd 29, 136(%r5)       /* uu+17 */
/* done loading for current iteration */
       	fnmsub 19, 3, 27, 19  /* *(fp1+4) -= *(fp0+3)* *(uu+15) */
/* start loading for next iteration */
	lfd 0,	0(%r11)		/* fp0 */ /*sometimes a cache miss*/
	fmadd 25, 2, 27, 25  /*  *(fp1+5) += *(fp0+2) * *(uu+15) */
	lfd 1,	8(%r11)	/* fp0+1 */ /*sometimes a cache miss*/
/* now switch fpr's 30-14, 31-20, 24-26, 13-15, 19-21, 25-27 */
	fmadd 14,  4, 16, 30   /*  *(fp1) += *(fp0+4) * *(uu+4) */

	fmadd 20,  5, 16, 31   /*  *(fp1+1) += *(fp0+5) * *(uu+4) */

	fmadd 26,  4, 22, 24   /* *(fp1+2) += *(fp0+4) * *(uu+10) */
	lfd 12, 144(%r5)    /* uu */	 		
	fmadd 15,  5, 22, 13   /* *(fp1+3) += *(fp0+5) * *(uu+10) */
	lfd 18, 192(%r5)        /* uu+6 */			
	fmadd 21, 4, 28, 19   /*  *(fp1+4) += *(fp0+4) * *(uu+16) */
	lfd 24, 240(%r5)        /* uu+12 */		
	fmadd 27, 5, 28, 25   /*  *(fp1+5) += *(fp0+5) * *(uu+16) */
	lfd 13, 152(%r5)         /* uu+1 */	
	fnmsub 14,  5, 17, 14   /*  *(fp1) -= *(fp0+5) * *(uu+5) */

	fmadd 20,  4, 17, 20   /*  *(fp1+1) += *(fp0+4) * *(uu+5) */

	fnmsub 26,  5, 23, 26  /* *(fp1+2) -= *(fp0+5) * *(uu+11) */
	lfd 19, 200(%r5)        /* uu+7 */	
	fmadd 15,  4, 23, 15   /* *(fp1+3) += *(fp0+4) * *(uu+11) */
        lfd 25, 248(%r5)       /* uu+13 */		
	fnmsub 21, 5, 29, 21  /* *(fp1+4) -= *(fp0+5) * *(uu+17) */
	lfd 2,	16(%r11)	/* fp0+2 */		
	fmadd 27, 4, 29, 27  /*  *(fp1+5) += *(fp0+4) * *(uu+17) */
	stfd  14,  0(%r14)   /* fp1 */	
	fmul 6,  0, 12   /*  *(fp1) += *fp0 * *(uu) */
	stfd  20,  8(%r14)   /* fp1+1 */	
	fmul 7,  1, 12   /*  *(fp1+1) += *(fp0+1) * *(uu) */
	stfd  26,  16(%r14)  /* fp1+2 */	
	fmul 8,  0, 18   /*  *(fp1+2) += *fp0 * *(uu+6) */
	stfd  15,  24(%r14)  /* fp1+3 */
	fmul 9,  1, 18   /*  *(fp1+3) += *(fp0+1) * *(uu+6) */	
	stfd  21, 32(%r14)  /* fp1+4 */
	fmul 10, 0, 24  /*  *(fp1+4) += *fp0 * *(uu+12) */
	stfd  27, 40(%r14)  /* fp1+5 */
	fmul 11, 1, 24  /*  *(fp1+5) += *(fp0+1) * *(uu+12) */
	lfd 3,	24(%r11)	/* fp0+3 */				
	fnmsub 6,  1, 13, 6   /* *(fp1) -= *(fp0+1) * *(uu+1) */	
        lfd 14, 160(%r5)        /* uu+2 */
	fmadd 7,  0, 13, 7   /*  *(fp1+1) += *fp0 * *(uu+1) */
        lfd 20, 208(%r5)        /* uu+8 */	
	fnmsub 8,  1, 19, 8   /*  *(fp1+2) -= *(fp0+1) * *(uu+7) */
        lfd 26, 256(%r5)       /* uu+14 */	
	fmadd 9,  0, 19, 9   /*  *(fp1+3) += *fp0 * *(uu+7) */
        lfd 15, 168(%r5)        /* uu+3 */
	fnmsub 10, 1, 25, 10  /*  *(fp1+4) -= *(fp0+1) * *(uu+13) */
        lfd 21, 216(%r5)        /* uu+9 */
	fmadd 11, 0, 25, 11  /*  *(fp1+5) += *fp0 * *(uu+13) */

/* now switch fpr's 6-12, 7-18, 8-24, 9-13, 10-19, 11-25 */	
	fmadd 12,  2, 14, 6   /*  *(fp1) += *(fp0+2) * *(uu+2) */

	fmadd 18,  3, 14, 7   /*  *(fp1+1) += *(fp0+3) * *(uu+2) */
        lfd 27, 264(%r5)       /* uu+15 */	
	fmadd 24,  2, 20, 8   /*  *(fp1+2) +=*(fp0+2) * *(uu+8) */
        lfd 16, 176(%r5)        /* uu+4 */	
	fmadd 13,  3, 20, 9   /*  *(fp1+3) += *(fp0+3) * *(uu+8) */
	lfd 5,	40(%r11)	/* fp0+5 */	
	fmadd 19, 2, 26, 10  /*  *(fp1+4) +=*(fp0+2) * *(uu+14) */
	lfd 4,	32(%r11)	/* fp0+4 */
	fmadd 25, 3, 26, 11  /*  *(fp1+5) += *(fp0+3) * *(uu+14) */

	fnmsub 12,  3, 15, 12   /*  *(fp1) -= *(fp0+3)* *(uu+3) */

	fmadd 18,  2, 15, 18   /*  *(fp1+1) += *(fp0+2) * *(uu+3) */

        lfd 22, 224(%r5)        /* uu+10 */	
	fnmsub 24,  3, 21, 24   /*  *(fp1+2) -= *(fp0+3)* *(uu+9) */
        lfd 28, 272(%r5)       /* uu+16 */
	fmadd 13,  2, 21, 13   /*  *(fp1+3) += *(fp0+2) * *(uu+9) */
        lfd 17, 184(%r5)        /* uu+5 */		
       	fnmsub 19, 3, 27, 19  /* *(fp1+4) -= *(fp0+3)* *(uu+15) */
        lfd 23, 232(%r5)        /* uu+11 */		
	fmadd 25, 2, 27, 25  /*  *(fp1+5) += *(fp0+2) * *(uu+15) */
/* now switch fpr's 12-14, 18-20, 24-26, 13-15, 19-21, 25-27 */
	fmadd 14,  4, 16, 12   /*  *(fp1) += *(fp0+4) * *(uu+4) */

	fmadd 20,  5, 16, 18   /*  *(fp1+1) += *(fp0+5) * *(uu+4) */
        lfd 29, 280(%r5)       /* uu+17 */
/* done loading for current iteration */
	fmadd 26,  4, 22, 24   /* *(fp1+2) += *(fp0+4) * *(uu+10) */
/* start loading for next iteration */

	fmadd 15,  5, 22, 13   /* *(fp1+3) += *(fp0+5) * *(uu+10) */

	fmadd 21, 4, 28, 19   /*  *(fp1+4) += *(fp0+4) * *(uu+16) */

	fmadd 27, 5, 28, 25   /*  *(fp1+5) += *(fp0+5) * *(uu+16) */

	fnmsub 14,  5, 17, 14   /*  *(fp1) -= *(fp0+5) * *(uu+5) */

	fmadd 20,  4, 17, 20   /*  *(fp1+1) += *(fp0+4) * *(uu+5) */

	fnmsub 26,  5, 23, 26  /* *(fp1+2) -= *(fp0+5) * *(uu+11) */

	fmadd 15,  4, 23, 15   /* *(fp1+3) += *(fp0+4) * *(uu+11) */

	fnmsub 21, 5, 29, 21  /* *(fp1+4) -= *(fp0+5) * *(uu+17) */

	fmadd 27, 4, 29, 27  /*  *(fp1+5) += *(fp0+4) * *(uu+17) */
	stfd  14,  0(%r13)
	stfd  20,  8(%r13)	
	stfd  26,  16(%r13)	
	stfd  15,  24(%r13)
	stfd  21, 32(%r13)
	stfd  27, 40(%r13)

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
        .size    dirac_cmv,.Lfe1-dirac_cmv
