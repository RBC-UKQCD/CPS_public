/*
 * PAB, calling these routines allows us to reduce
 * the register save and restore overhead in the wilson
 * kernels. This is a non-trivial gain on 2*2*2*2
 * the C-routine calling this needs to be compiled
 * with -fcall-used-%r14...31 -fcall-used-%fr14...31
 * In particular doing the whole wilson_mdagm without save
 * and restore would save something like 1000 cycles, or 2 microseconds.
 */
.data
.align 3
wfm_abi_save_area:
.space 1024
	
.text
.align 2
.globl wfm_abi_save
wfm_abi_save:
	lis %r3,wfm_abi_save_area@h
	ori %r3,%r3,wfm_abi_save_area@l
	stw   %r14,	368(%r3)
	stw   %r15,	376(%r3)
	stw   %r16,	384(%r3)
	stw   %r17,	392(%r3)
	stw   %r18,	400(%r3)
	stw   %r19,	408(%r3)
	stw   %r20,	416(%r3)
	stw   %r21,	424(%r3)
	stw   %r22,	432(%r3)
	stw   %r23,	440(%r3)
	stw   %r24,	448(%r3)
	stw   %r25,	456(%r3)
	stw   %r26,	464(%r3)
	stw   %r27,	472(%r3)
	stw   %r28,	480(%r3)
	stw   %r29,	488(%r3)
	stw   %r30,	496(%r3)
	stw   %r31,	504(%r3)
	stfd  14,	112(%r3)
	stfd  15,	120(%r3)
	stfd  16,	128(%r3)
	stfd  17,	136(%r3)
	stfd  18,	144(%r3)
	stfd  19,	152(%r3)
	stfd  20,	160(%r3)
	stfd  21,	168(%r3)
	stfd  22,	176(%r3)
	stfd  23,	184(%r3)
	stfd  24,	192(%r3)
	stfd  25,	200(%r3)
	stfd  26,	208(%r3)
	stfd  27,	216(%r3)
	stfd  28,	224(%r3)
	stfd  29,	232(%r3)
	stfd  30,	240(%r3)
	stfd  31,	248(%r3)
	blr
	
.align 2
.globl wfm_abi_restore
wfm_abi_restore:	
	lis %r3,wfm_abi_save_area@h
	ori %r3,%r3,wfm_abi_save_area@l
	lwz   %r14,	368(%r3)
	lwz   %r15,	376(%r3)
	lwz   %r16,	384(%r3)
	lwz   %r17,	392(%r3)
	lwz   %r18,	400(%r3)
	lwz   %r19,	408(%r3)
	lwz   %r20,	416(%r3)
	lwz   %r21,	424(%r3)
	lwz   %r22,	432(%r3)
	lwz   %r23,	440(%r3)
	lwz   %r24,	448(%r3)
	lwz   %r25,	456(%r3)
	lwz   %r26,	464(%r3)
	lwz   %r27,	472(%r3)
	lwz   %r28,	480(%r3)
	lwz   %r29,	488(%r3)
	lwz   %r30,	496(%r3)
	lwz   %r31,	504(%r3)
	lfd   14,	112(%r3)
	lfd   15,	120(%r3)
	lfd   16,	128(%r3)
	lfd   17,	136(%r3)
	lfd   18,	144(%r3)
	lfd   19,	152(%r3)
	lfd   20,	160(%r3)
	lfd   21,	168(%r3)
	lfd   22,	176(%r3)
	lfd   23,	184(%r3)
	lfd   24,	192(%r3)
	lfd   25,	200(%r3)
	lfd   26,	208(%r3)
	lfd   27,	216(%r3)
	lfd   28,	224(%r3)
	lfd   29,	232(%r3)
	lfd   30,	240(%r3)
	lfd   31,	248(%r3)
	blr



