        .section        ".text"
        .align 2
        .globl save_reg
        .type    save_reg,@function
save_reg:
        stw   %r14,     0(%r3)
        stw   %r15,     4(%r3)
        stw   %r16,     8(%r3)
        stw   %r17,     12(%r3)
        stw   %r18,     16(%r3)
        stw   %r19,     20(%r3)
        stw   %r20,     24(%r3)
        stw   %r21,     28(%r3)
        stw   %r22,     32(%r3)
        stfd  14,       0(%r4)
        stfd  15,       8(%r4)
        stfd  16,       16(%r4)
        stfd  17,       24(%r4)
        stfd  18,       32(%r4)
        stfd  19,       40(%r4)
        stfd  20,       48(%r4)
        stfd  21,       56(%r4)
        stfd  22,       64(%r4)
        stfd  23,       72(%r4)
        stfd  24,       80(%r4)
        stfd  25,       88(%r4)
        stfd  26,       96(%r4)
        stfd  27,       104(%r4)
        stfd  28,       112(%r4)
        stfd  29,       120(%r4)
        stfd  30,       128(%r4)
        stfd  31,       136(%r4)

	blr

        .section        ".text"
        .align 2
        .globl restore_reg
        .type    restore_reg,@function
restore_reg:
        lwz   %r14,     0(%r3)
        lwz   %r15,     4(%r3)
        lwz   %r16,     8(%r3)
        lwz   %r17,     12(%r3)
        lwz   %r18,     16(%r3)
        lwz   %r19,     20(%r3)
        lwz   %r20,     24(%r3)
        lwz   %r21,     28(%r3)
        lwz   %r22,     32(%r3)
        lfd  14,       0(%r4)
        lfd  18,       32(%r4)
        lfd  22,       64(%r4)	
        lfd  15,       8(%r4)
        lfd  16,       16(%r4)
        lfd  17,       24(%r4)
        lfd  19,       40(%r4)
        lfd  20,       48(%r4)
        lfd  21,       56(%r4)
        lfd  23,       72(%r4)
        lfd  24,       80(%r4)
        lfd  25,       88(%r4)
        lfd  26,       96(%r4)
        lfd  27,       104(%r4)
        lfd  28,       112(%r4)
        lfd  29,       120(%r4)
        lfd  30,       128(%r4)	
        lfd  31,       136(%r4)

	blr
