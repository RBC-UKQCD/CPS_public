/* 
 * BAGEL machine generated output.  
 * Copyright UKQCD Collaboration  
 * This software is provided for NON-COMMERCIAL use only,  
 * and may not be redistributed without permission.  
 * It is provided as is, and is not guaranteed fit for any purpose.
 * Written by Peter Boyle
 */  
#define CACHELINE 32 
#define Isize 4 
#define Bdy -3 
#define Nahead 1 
#define SHFT_PREF_IMM 4 
#define zeromulti 0 
#define onemulti 4 
#define twomulti 8 
#define threemulti 12 
#define fourmulti 16 
#define fivemulti 20 
#define minusone -1 
#define ZERO 0 
#define VEC_ATOM 48 
#define MAT_ATOM 144 
#define VEC_IMM00 0 
#define VEC_IMM01 8 
#define VEC_IMM10 16 
#define VEC_IMM11 24 
#define VEC_IMM20 32 
#define VEC_IMM21 40 
#define MAT_IMM000 0 
#define MAT_IMM001 8 
#define MAT_IMM010 16 
#define MAT_IMM011 24 
#define MAT_IMM020 32 
#define MAT_IMM021 40 
#define MAT_IMM100 48 
#define MAT_IMM101 56 
#define MAT_IMM110 64 
#define MAT_IMM111 72 
#define MAT_IMM120 80 
#define MAT_IMM121 88 
#define MAT_IMM200 96 
#define MAT_IMM201 104 
#define MAT_IMM210 112 
#define MAT_IMM211 120 
#define MAT_IMM220 128 
#define MAT_IMM221 136 
#define frame_down -736 
#define frame_up 736 
#define bias 736 
#define S_I_14 368 
#define S_I_15 376 
#define S_I_16 384 
#define S_I_17 392 
#define S_I_18 400 
#define S_I_19 408 
#define S_I_20 416 
#define S_I_21 424 
#define S_I_22 432 
#define S_I_23 440 
#define S_I_24 448 
#define S_I_25 456 
#define S_I_26 464 
#define S_I_27 472 
#define S_I_28 480 
#define S_I_29 488 
#define S_I_30 496 
#define S_I_31 504 
#define S_F_14 112 
#define S_F_15 120 
#define S_F_16 128 
#define S_F_17 136 
#define S_F_18 144 
#define S_F_19 152 
#define S_F_20 160 
#define S_F_21 168 
#define S_F_22 176 
#define S_F_23 184 
#define S_F_24 192 
#define S_F_25 200 
#define S_F_26 208 
#define S_F_27 216 
#define S_F_28 224 
#define S_F_29 232 
#define S_F_30 240 
#define S_F_31 248 
#define PRE0 0 
#define PRE1 32 
#define PRE2 64 
#define PRE3 96 
#define PRE4 128 
#define minus1 -1 

#define A 0
#define vec100 1
#define vec101 2
#define vec110 3
#define vec111 4
#define vec120 5
#define vec121 6
#define vec200 7
#define vec201 8
#define vec210 9
#define vec211 10
#define vec220 11
#define vec221 12
#define mat000 13
#define mat001 14
#define mat010 15
#define mat011 16
#define mat020 17
#define mat021 18
#define mat100 19
#define mat101 20
#define mat110 21
#define mat111 22
#define mat120 23
#define mat121 24
#define mat200 25
#define mat201 26
#define mat210 27
#define mat211 28
#define mat220 29
#define mat221 30
#define GPR0 %r0
#define SP %r1
#define GPR2 %r2
#define Arg0 %r3
#define Arg1 %r4
#define Arg2 %r5
#define Arg3 %r6
#define Arg4 %r7
#define Arg5 %r8
#define Arg6 %r9
#define predicate %r10
#define prefscratch %r11
#define vec1ptr %r12
#define vec2ptr %r13
#define matptr %r14
#define counter %r15
#define Aptr %r16
#define outptr %r17
#define src1_idx %r18
#define dest_idx %r19
#define src2_idx %r20
#define vec1tmp %r21
#define vec2tmp %r22
#define LOOKUP0 %r23
#define PRER0 %r24
#define PRER1 %r25
#define LOOKUP1 %r26
#define LOOKUP2 %r27
#define PRER2 %r28
#define PRER3 %r29
#define PRER4 %r30
.text
.align 2
.globl cross_look
cross_look:

	or    vec1ptr , Arg2 , Arg2
	   la   SP,	frame_down(SP)
	or    vec2ptr , Arg3 , Arg3
	   std   matptr,	S_I_14(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   std   counter,	S_I_15(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	or    matptr , Arg0 , Arg0
			/* Last Dependency:	   std   matptr,	S_I_14(SP)*/
	   std   Aptr,	S_I_16(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	or    counter , Arg4 , Arg4
			/* Last Dependency:	   std   counter,	S_I_15(SP)*/
	   std   outptr,	S_I_17(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	or    Aptr , Arg1 , Arg1
			/* Last Dependency:	   std   Aptr,	S_I_16(SP)*/
	   std   src1_idx,	S_I_18(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	or.  counter , counter , counter
			/* Last Dependency:	or    counter , Arg4 , Arg4*/
	   std   dest_idx,	S_I_19(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   std   src2_idx,	S_I_20(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	or    dest_idx , Arg6 , Arg6
			/* Last Dependency:	   std   dest_idx,	S_I_19(SP)*/
	   std   vec1tmp,	S_I_21(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	or    src2_idx , Arg5 , Arg5
			/* Last Dependency:	   std   src2_idx,	S_I_20(SP)*/
	   std   vec2tmp,	S_I_22(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	or    src1_idx , dest_idx , dest_idx
			/* Last Dependency:	or    dest_idx , Arg6 , Arg6*/
	   std   LOOKUP0,	S_I_23(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   std   PRER0,	S_I_24(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	or    LOOKUP0 , vec1ptr , vec1ptr
			/* Last Dependency:	   std   LOOKUP0,	S_I_23(SP)*/
	   std   PRER1,	S_I_25(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   std   LOOKUP1,	S_I_26(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   std   LOOKUP2,	S_I_27(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	or    LOOKUP1 , vec2ptr , vec2ptr
			/* Last Dependency:	   std   LOOKUP1,	S_I_26(SP)*/
	   std   PRER2,	S_I_28(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	or    LOOKUP2 , matptr , matptr
			/* Last Dependency:	   std   LOOKUP2,	S_I_27(SP)*/
	   std   PRER3,	S_I_29(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   std   PRER4,	S_I_30(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   std   %r31,	S_I_31(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat001,	S_F_14(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat010,	S_F_15(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat011,	S_F_16(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat020,	S_F_17(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat021,	S_F_18(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat100,	S_F_19(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat101,	S_F_20(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat110,	S_F_21(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat111,	S_F_22(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat120,	S_F_23(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat121,	S_F_24(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat200,	S_F_25(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat201,	S_F_26(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat210,	S_F_27(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat211,	S_F_28(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat220,	S_F_29(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  mat221,	S_F_30(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   stfd  31,	S_F_31(SP)
			/* Last Dependency:	   la   SP,	frame_down(SP)*/
	   lfd   A,	ZERO(Aptr)
			/* Last Dependency:	or    Aptr , Arg1 , Arg1*/
	li   PRER0,PRE0
			/* Last Dependency:	   std   PRER0,	S_I_24(SP)*/
	li   PRER1,PRE1
			/* Last Dependency:	   std   PRER1,	S_I_25(SP)*/
	   lwz   prefscratch,	PRE0(src1_idx)
			/* Last Dependency:	or    src1_idx , dest_idx , dest_idx*/
	li   PRER2,PRE2
			/* Last Dependency:	   std   PRER2,	S_I_28(SP)*/
	mulli  prefscratch , prefscratch , VEC_ATOM
			/* Last Dependency:	   lwz   prefscratch,	PRE0(src1_idx)*/
	li   PRER3,PRE3
			/* Last Dependency:	   std   PRER3,	S_I_29(SP)*/
	add   vec1ptr , LOOKUP0 , prefscratch
			/* Last Dependency:	mulli  prefscratch , prefscratch , VEC_ATOM*/
	li   PRER4,PRE4
			/* Last Dependency:	   std   PRER4,	S_I_30(SP)*/
	   lwz   prefscratch,	PRE0(src2_idx)
			/* Last Dependency:	add   vec1ptr , LOOKUP0 , prefscratch*/
	mulli  prefscratch , prefscratch , VEC_ATOM
			/* Last Dependency:	   lwz   prefscratch,	PRE0(src2_idx)*/
	add   vec2ptr , LOOKUP1 , prefscratch
			/* Last Dependency:	mulli  prefscratch , prefscratch , VEC_ATOM*/
	   lwz   prefscratch,	PRE0(dest_idx)
			/* Last Dependency:	add   vec2ptr , LOOKUP1 , prefscratch*/
	mulli  prefscratch , prefscratch , MAT_ATOM
			/* Last Dependency:	   lwz   prefscratch,	PRE0(dest_idx)*/
	add   matptr , LOOKUP2 , prefscratch
			/* Last Dependency:	mulli  prefscratch , prefscratch , MAT_ATOM*/
	bf gt,  cross_look_lab0
			/* Last Dependency:	add   matptr , LOOKUP2 , prefscratch*/
	   lfd   vec100,	VEC_IMM00(vec1ptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/
	   lfd   vec101,	VEC_IMM01(vec1ptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/
	   lfd   vec200,	VEC_IMM00(vec2ptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/
	   lfd   vec201,	VEC_IMM01(vec2ptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/
	   lfd   vec210,	VEC_IMM10(vec2ptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/
	   lfd   vec211,	VEC_IMM11(vec2ptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/
	   lfd   vec220,	VEC_IMM20(vec2ptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/

	fmul   vec100 , A , vec100
			/* Last Dependency:	   lfd   vec100,	VEC_IMM00(vec1ptr)*/
	   lfd   vec221,	VEC_IMM21(vec2ptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/

	fmul   vec101 , A , vec101
			/* Last Dependency:	   lfd   vec101,	VEC_IMM01(vec1ptr)*/
	   lfd   mat000,	MAT_IMM000(matptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/
	   lfd   mat001,	MAT_IMM001(matptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/
	   lfd   mat010,	MAT_IMM010(matptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/
	   lfd   mat011,	MAT_IMM011(matptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/
	   lfd   mat020,	MAT_IMM020(matptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/
	   lfd   mat021,	MAT_IMM021(matptr)
			/* Last Dependency:	bf gt,  cross_look_lab0*/
	mtctr 	counter
			/* Last Dependency:	bf gt,  cross_look_lab0*/
	b       cross_look_lab1
.align 4
			/* Last Dependency:	b       cross_look_lab1*/
cross_look_lab1:
	addi  outptr , matptr , ZERO
	addi  vec1tmp , vec1ptr , ZERO
	addi  vec2tmp , vec2ptr , ZERO
	   lfd   vec110,	VEC_IMM10(vec1tmp)
			/* Last Dependency:	addi  vec1tmp , vec1ptr , ZERO*/

	fmadd  mat000 , vec200 , vec100 , mat000
	   lfd   vec111,	VEC_IMM11(vec1tmp)
			/* Last Dependency:	addi  vec1tmp , vec1ptr , ZERO*/

	fmadd  mat001 , vec200 , vec101 , mat001
	   lfd   mat100,	MAT_IMM100(matptr)

	fmadd  mat010 , vec210 , vec100 , mat010
	   lfd   mat101,	MAT_IMM101(matptr)

	fmadd  mat011 , vec210 , vec101 , mat011
	   lfd   mat110,	MAT_IMM110(matptr)

	fmadd  mat020 , vec220 , vec100 , mat020
			/* Last Dependency:	   lfd   mat020,	MAT_IMM020(matptr)*/
	   lfd   mat111,	MAT_IMM111(matptr)

	fmadd  mat000 , vec201 , vec101 , mat000
			/* Last Dependency:	fmadd  mat000 , vec200 , vec100 , mat000*/
	   lfd   mat120,	MAT_IMM120(matptr)

	fnmsub mat001 , vec201 , vec100 , mat001
			/* Last Dependency:	fmadd  mat001 , vec200 , vec101 , mat001*/
	   lfd   mat121,	MAT_IMM121(matptr)

	fmadd  mat010 , vec211 , vec101 , mat010
			/* Last Dependency:	fmadd  mat010 , vec210 , vec100 , mat010*/
	   la   src1_idx,	Isize(src1_idx)

	fnmsub mat011 , vec211 , vec100 , mat011
			/* Last Dependency:	fmadd  mat011 , vec210 , vec101 , mat011*/
	or.  predicate , predicate , predicate

	fmadd  mat021 , vec220 , vec101 , mat021
			/* Last Dependency:	   lfd   mat021,	MAT_IMM021(matptr)*/
	   lwz   prefscratch,	PRE0(src1_idx)
			/* Last Dependency:	   la   src1_idx,	Isize(src1_idx)*/

	fmadd  mat020 , vec221 , vec101 , mat020
			/* Last Dependency:	fmadd  mat020 , vec220 , vec100 , mat020*/
	   stfd  mat000,	MAT_IMM000(outptr)
			/* Last Dependency:	fmadd  mat000 , vec201 , vec101 , mat000*/
	mulli  prefscratch , prefscratch , VEC_ATOM
			/* Last Dependency:	   lwz   prefscratch,	PRE0(src1_idx)*/
	   stfd  mat001,	MAT_IMM001(outptr)
			/* Last Dependency:	fnmsub mat001 , vec201 , vec100 , mat001*/
	add   vec1ptr , LOOKUP0 , prefscratch
			/* Last Dependency:	mulli  prefscratch , prefscratch , VEC_ATOM*/
	   stfd  mat010,	MAT_IMM010(outptr)
			/* Last Dependency:	fmadd  mat010 , vec211 , vec101 , mat010*/
	   stfd  mat011,	MAT_IMM011(outptr)
			/* Last Dependency:	fnmsub mat011 , vec211 , vec100 , mat011*/

	fnmsub mat021 , vec221 , vec100 , mat021
			/* Last Dependency:	fmadd  mat021 , vec220 , vec101 , mat021*/
	   lwz   prefscratch,	SHFT_PREF_IMM(src1_idx)
			/* Last Dependency:	add   vec1ptr , LOOKUP0 , prefscratch*/
	mulli  prefscratch , prefscratch , VEC_ATOM
			/* Last Dependency:	   lwz   prefscratch,	SHFT_PREF_IMM(src1_idx)*/
	add   prefscratch , LOOKUP0 , prefscratch
			/* Last Dependency:	mulli  prefscratch , prefscratch , VEC_ATOM*/
	bf lt,  cross_look_lab2
			/* Last Dependency:	add   prefscratch , LOOKUP0 , prefscratch*/
	or    prefscratch , SP , SP
			/* Last Dependency:	bf lt,  cross_look_lab2*/
cross_look_lab2:
			/* Last Dependency:	or    prefscratch , SP , SP*/

	fmul   vec110 , A , vec110
	dcbt 	PRER0,prefscratch

	fmul   vec111 , A , vec111
	   stfd  mat020,	MAT_IMM020(outptr)
	or.  predicate , predicate , predicate
	   stfd  mat021,	MAT_IMM021(outptr)
	   lfd   vec120,	VEC_IMM20(vec1tmp)
	   lfd   vec121,	VEC_IMM21(vec1tmp)

	fmadd  mat100 , vec200 , vec110 , mat100
			/* Last Dependency:	fmul   vec110 , A , vec110*/
	   lfd   mat200,	MAT_IMM200(outptr)

	fmadd  mat101 , vec200 , vec111 , mat101
			/* Last Dependency:	fmul   vec111 , A , vec111*/
	dcbt 	PRER1,prefscratch
			/* Last Dependency:	Stalled by Cache unit GPR0*/
	   lfd   mat201,	MAT_IMM201(outptr)
	   lfd   mat210,	MAT_IMM210(outptr)
	   lfd   mat211,	MAT_IMM211(outptr)

	fmul   vec120 , A , vec120
			/* Last Dependency:	   lfd   vec120,	VEC_IMM20(vec1tmp)*/
	   lfd   mat220,	MAT_IMM220(outptr)

	fmul   vec121 , A , vec121
			/* Last Dependency:	   lfd   vec121,	VEC_IMM21(vec1tmp)*/
	   lfd   mat221,	MAT_IMM221(outptr)

	fmadd  mat100 , vec201 , vec111 , mat100
			/* Last Dependency:	fmadd  mat100 , vec200 , vec110 , mat100*/
	dcbt 	PRER1,src1_idx
			/* Last Dependency:	Stalled by Cache unit GPR0*/

	fnmsub mat101 , vec201 , vec110 , mat101
			/* Last Dependency:	fmadd  mat101 , vec200 , vec111 , mat101*/
	   la   src2_idx,	Isize(src2_idx)
	   lwz   prefscratch,	PRE0(src2_idx)
			/* Last Dependency:	   la   src2_idx,	Isize(src2_idx)*/
	mulli  prefscratch , prefscratch , VEC_ATOM
			/* Last Dependency:	   lwz   prefscratch,	PRE0(src2_idx)*/
	add   vec2ptr , LOOKUP1 , prefscratch
			/* Last Dependency:	mulli  prefscratch , prefscratch , VEC_ATOM*/
	   lwz   prefscratch,	SHFT_PREF_IMM(src2_idx)
			/* Last Dependency:	add   vec2ptr , LOOKUP1 , prefscratch*/
	mulli  prefscratch , prefscratch , VEC_ATOM
			/* Last Dependency:	   lwz   prefscratch,	SHFT_PREF_IMM(src2_idx)*/
	add   prefscratch , LOOKUP1 , prefscratch
			/* Last Dependency:	mulli  prefscratch , prefscratch , VEC_ATOM*/
	bf lt,  cross_look_lab3
			/* Last Dependency:	add   prefscratch , LOOKUP1 , prefscratch*/
	or    prefscratch , SP , SP
			/* Last Dependency:	bf lt,  cross_look_lab3*/
cross_look_lab3:
			/* Last Dependency:	or    prefscratch , SP , SP*/
	or.  predicate , predicate , predicate
	dcbt 	PRER0,prefscratch
	   stfd  mat100,	MAT_IMM100(outptr)
	   stfd  mat101,	MAT_IMM101(outptr)
	   la   dest_idx,	Isize(dest_idx)
	dcbt 	PRER1,prefscratch
			/* Last Dependency:	Stalled by Cache unit GPR0*/
	   lwz   prefscratch,	PRE0(dest_idx)
	mulli  prefscratch , prefscratch , MAT_ATOM
			/* Last Dependency:	   lwz   prefscratch,	PRE0(dest_idx)*/
	add   matptr , LOOKUP2 , prefscratch
			/* Last Dependency:	mulli  prefscratch , prefscratch , MAT_ATOM*/
	dcbt 	PRER1,src2_idx
			/* Last Dependency:	Stalled by Cache unit GPR0*/
	   lwz   prefscratch,	SHFT_PREF_IMM(dest_idx)
			/* Last Dependency:	add   matptr , LOOKUP2 , prefscratch*/
	mulli  prefscratch , prefscratch , MAT_ATOM
			/* Last Dependency:	   lwz   prefscratch,	SHFT_PREF_IMM(dest_idx)*/
	add   prefscratch , LOOKUP2 , prefscratch
			/* Last Dependency:	mulli  prefscratch , prefscratch , MAT_ATOM*/
	bf lt,  cross_look_lab4
			/* Last Dependency:	add   prefscratch , LOOKUP2 , prefscratch*/
	or    prefscratch , SP , SP
			/* Last Dependency:	bf lt,  cross_look_lab4*/
cross_look_lab4:
			/* Last Dependency:	or    prefscratch , SP , SP*/

	fmadd  mat110 , vec210 , vec110 , mat110
	dcbt 	PRER0,prefscratch

	fmadd  mat111 , vec210 , vec111 , mat111
	   lfd   vec100,	VEC_IMM00(vec1ptr)

	fmadd  mat120 , vec220 , vec110 , mat120
	   lfd   vec101,	VEC_IMM01(vec1ptr)

	fmadd  mat121 , vec220 , vec111 , mat121
	   lfd   mat000,	MAT_IMM000(matptr)

	fmadd  mat200 , vec200 , vec120 , mat200
	   lfd   mat001,	MAT_IMM001(matptr)

	fmadd  mat110 , vec211 , vec111 , mat110
			/* Last Dependency:	fmadd  mat110 , vec210 , vec110 , mat110*/
	   lfd   mat010,	MAT_IMM010(matptr)

	fnmsub mat111 , vec211 , vec110 , mat111
			/* Last Dependency:	fmadd  mat111 , vec210 , vec111 , mat111*/
	dcbt 	PRER1,prefscratch
			/* Last Dependency:	Stalled by Cache unit GPR0*/

	fmadd  mat120 , vec221 , vec111 , mat120
			/* Last Dependency:	fmadd  mat120 , vec220 , vec110 , mat120*/
	   lfd   mat011,	MAT_IMM011(matptr)

	fnmsub mat121 , vec221 , vec110 , mat121
			/* Last Dependency:	fmadd  mat121 , vec220 , vec111 , mat121*/
	   lfd   mat020,	MAT_IMM020(matptr)

	fmadd  mat201 , vec200 , vec121 , mat201
	   lfd   mat021,	MAT_IMM021(matptr)

	fmadd  mat200 , vec201 , vec121 , mat200
			/* Last Dependency:	fmadd  mat200 , vec200 , vec120 , mat200*/
	   stfd  mat110,	MAT_IMM110(outptr)
			/* Last Dependency:	fmadd  mat110 , vec211 , vec111 , mat110*/

	fmadd  mat210 , vec210 , vec120 , mat210
	   stfd  mat111,	MAT_IMM111(outptr)
			/* Last Dependency:	fnmsub mat111 , vec211 , vec110 , mat111*/

	fmadd  mat211 , vec210 , vec121 , mat211
	dcbt 	PRER2,prefscratch
			/* Last Dependency:	Stalled by Cache unit GPR0*/

	fmadd  mat220 , vec220 , vec120 , mat220
	   stfd  mat120,	MAT_IMM120(outptr)
			/* Last Dependency:	fmadd  mat120 , vec221 , vec111 , mat120*/

	fnmsub mat201 , vec201 , vec120 , mat201
			/* Last Dependency:	fmadd  mat201 , vec200 , vec121 , mat201*/
	   stfd  mat121,	MAT_IMM121(outptr)
			/* Last Dependency:	fnmsub mat121 , vec221 , vec110 , mat121*/

	fmadd  mat221 , vec220 , vec121 , mat221
	   stfd  mat200,	MAT_IMM200(outptr)
			/* Last Dependency:	fmadd  mat200 , vec201 , vec121 , mat200*/

	fmadd  mat210 , vec211 , vec121 , mat210
			/* Last Dependency:	fmadd  mat210 , vec210 , vec120 , mat210*/
	   lfd   vec200,	VEC_IMM00(vec2ptr)
			/* Last Dependency:	fmadd  mat201 , vec200 , vec121 , mat201*/

	fnmsub mat211 , vec211 , vec120 , mat211
			/* Last Dependency:	fmadd  mat211 , vec210 , vec121 , mat211*/
	   lfd   vec201,	VEC_IMM01(vec2ptr)
			/* Last Dependency:	fnmsub mat201 , vec201 , vec120 , mat201*/

	fmadd  mat220 , vec221 , vec121 , mat220
			/* Last Dependency:	fmadd  mat220 , vec220 , vec120 , mat220*/
	dcbt 	PRER3,prefscratch
			/* Last Dependency:	Stalled by Cache unit GPR0*/

	fmul   vec100 , A , vec100
			/* Last Dependency:	   lfd   vec100,	VEC_IMM00(vec1ptr)*/
	   stfd  mat201,	MAT_IMM201(outptr)
			/* Last Dependency:	fnmsub mat201 , vec201 , vec120 , mat201*/

	fnmsub mat221 , vec221 , vec120 , mat221
			/* Last Dependency:	fmadd  mat221 , vec220 , vec121 , mat221*/
	   lfd   vec210,	VEC_IMM10(vec2ptr)
			/* Last Dependency:	fmadd  mat211 , vec210 , vec121 , mat211*/

	fmul   vec101 , A , vec101
			/* Last Dependency:	   lfd   vec101,	VEC_IMM01(vec1ptr)*/
	   stfd  mat210,	MAT_IMM210(outptr)
			/* Last Dependency:	fmadd  mat210 , vec211 , vec121 , mat210*/
	   stfd  mat211,	MAT_IMM211(outptr)
			/* Last Dependency:	fnmsub mat211 , vec211 , vec120 , mat211*/
	   stfd  mat220,	MAT_IMM220(outptr)
			/* Last Dependency:	fmadd  mat220 , vec221 , vec121 , mat220*/
	dcbt 	PRER4,prefscratch
			/* Last Dependency:	Stalled by Cache unit GPR0*/
	   stfd  mat221,	MAT_IMM221(outptr)
			/* Last Dependency:	fnmsub mat221 , vec221 , vec120 , mat221*/
	   lfd   vec211,	VEC_IMM11(vec2ptr)
			/* Last Dependency:	fnmsub mat211 , vec211 , vec120 , mat211*/
	   lfd   vec220,	VEC_IMM20(vec2ptr)
			/* Last Dependency:	fmadd  mat221 , vec220 , vec121 , mat221*/
	   lfd   vec221,	VEC_IMM21(vec2ptr)
			/* Last Dependency:	fnmsub mat221 , vec221 , vec120 , mat221*/
	   la   counter,	minus1(counter)
	dcbt 	PRER1,dest_idx
			/* Last Dependency:	Stalled by Cache unit GPR0*/
	bdnz    cross_look_lab1
cross_look_lab0:
			/* Last Dependency:	bdnz    cross_look_lab1*/
	   ld   matptr,	S_I_14(SP)
	   ld   counter,	S_I_15(SP)
	   ld   Aptr,	S_I_16(SP)
	   ld   outptr,	S_I_17(SP)
	   ld   src1_idx,	S_I_18(SP)
	   ld   dest_idx,	S_I_19(SP)
	   ld   src2_idx,	S_I_20(SP)
	   ld   vec1tmp,	S_I_21(SP)
	   ld   vec2tmp,	S_I_22(SP)
	   ld   LOOKUP0,	S_I_23(SP)
	   ld   PRER0,	S_I_24(SP)
	   ld   PRER1,	S_I_25(SP)
	   ld   LOOKUP1,	S_I_26(SP)
	   ld   LOOKUP2,	S_I_27(SP)
	   ld   PRER2,	S_I_28(SP)
	   ld   PRER3,	S_I_29(SP)
	   ld   PRER4,	S_I_30(SP)
	   ld   %r31,	S_I_31(SP)
	   lfd   mat001,	S_F_14(SP)
	   lfd   mat010,	S_F_15(SP)
	   lfd   mat011,	S_F_16(SP)
	   lfd   mat020,	S_F_17(SP)
	   lfd   mat021,	S_F_18(SP)
	   lfd   mat100,	S_F_19(SP)
	   lfd   mat101,	S_F_20(SP)
	   lfd   mat110,	S_F_21(SP)
	   lfd   mat111,	S_F_22(SP)
	   lfd   mat120,	S_F_23(SP)
	   lfd   mat121,	S_F_24(SP)
	   lfd   mat200,	S_F_25(SP)
	   lfd   mat201,	S_F_26(SP)
	   lfd   mat210,	S_F_27(SP)
	   lfd   mat211,	S_F_28(SP)
	   lfd   mat220,	S_F_29(SP)
	   lfd   mat221,	S_F_30(SP)
	   lfd   31,	S_F_31(SP)
	   la   SP,	frame_up(SP)
			/* Last Dependency:	   lfd   31,	S_F_31(SP)*/
	blr
			/* Last Dependency:	   la   SP,	frame_up(SP)*/
