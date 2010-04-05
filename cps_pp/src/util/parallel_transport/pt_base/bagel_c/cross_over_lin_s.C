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
#define VEC_ATOM 24 
#define MAT_ATOM 72 
#define VEC_IMM00 0 
#define VEC_IMM01 4 
#define VEC_IMM10 8 
#define VEC_IMM11 12 
#define VEC_IMM20 16 
#define VEC_IMM21 20 
#define MAT_IMM000 0 
#define MAT_IMM001 4 
#define MAT_IMM010 8 
#define MAT_IMM011 12 
#define MAT_IMM020 16 
#define MAT_IMM021 20 
#define MAT_IMM100 24 
#define MAT_IMM101 28 
#define MAT_IMM110 32 
#define MAT_IMM111 36 
#define MAT_IMM120 40 
#define MAT_IMM121 44 
#define MAT_IMM200 48 
#define MAT_IMM201 52 
#define MAT_IMM210 56 
#define MAT_IMM211 60 
#define MAT_IMM220 64 
#define MAT_IMM221 68 
#define bias 0 
#define PRE0 0 
#define PRE1 32 
#define PRE2 64 
#define minus1 -1 

typedef float Float;

#ifdef __cplusplus
extern "C" {
#endif
 void cross_over_lin_s (unsigned long Arg0,unsigned long Arg1,unsigned long Arg2,unsigned long Arg3,unsigned long Arg4,unsigned long Arg5,unsigned long Arg6)
{

  /*Ensure natural alignment*/  int iStackArea[0+1];

  char *StackArea=(char *)iStackArea;

  unsigned long predicate;
  unsigned long prefscratch;
  unsigned long vec1ptr;
  unsigned long vec2ptr;
  unsigned long matptr;
  unsigned long counter;
  unsigned long Aptr;
  unsigned long outptr;
  unsigned long src1_idx;
  unsigned long dest_idx;
  unsigned long src2_idx;
  unsigned long vec1tmp;
  unsigned long vec2tmp;
  unsigned long LOOKUP0;
  unsigned long PRER0;
  unsigned long PRER1;
  unsigned long LOOKUP1;
  unsigned long PRER2;
  unsigned long R25;
  unsigned long R26;
  unsigned long R27;
  unsigned long R28;
  unsigned long R29;
  unsigned long R31;
  unsigned long StackPointer = (unsigned long) StackArea ;
  Float A;
  Float vec100;
  Float vec101;
  Float vec110;
  Float vec111;
  Float vec120;
  Float vec121;
  Float vec200;
  Float vec201;
  Float vec210;
  Float vec211;
  Float vec220;
  Float vec221;
  Float mat000;
  Float mat001;
  Float mat010;
  Float mat011;
  Float mat020;
  Float mat021;
  Float mat100;
  Float mat101;
  Float mat110;
  Float mat111;
  Float mat120;
  Float mat121;
  Float mat200;
  Float mat201;
  Float mat210;
  Float mat211;
  Float mat220;
  Float mat221;
  Float F31;
  matptr = Arg0 | Arg0 ; 
  Aptr = Arg1 | Arg1 ; 
  vec1ptr = Arg2 | Arg2 ; 
  vec2ptr = Arg3 | Arg3 ; 
  counter = Arg4 | Arg4 ; 
  src2_idx = Arg5 | Arg5 ; 
  dest_idx = Arg6 | Arg6 ; 
  src1_idx = dest_idx | dest_idx ; 
  A = *( (Float *) ( (ZERO) + (Aptr) ) ) ; 
  LOOKUP0 = vec1ptr | vec1ptr ; 
PRER0=PRE0;
PRER1=PRE1;
  prefscratch = * ( (unsigned long *) ( (PRE0)+(src1_idx) ) )  ;
  prefscratch = prefscratch * VEC_ATOM ; 
  vec1ptr = LOOKUP0 + prefscratch ; 
  LOOKUP1 = matptr | matptr ; 
PRER2=PRE2;
  prefscratch = * ( (unsigned long *) ( (PRE0)+(dest_idx) ) )  ;
  prefscratch = prefscratch * MAT_ATOM ; 
  matptr = LOOKUP1 + prefscratch ; 
  if ( counter <= 0 ) goto cross_over_lin_s_lab0;
  vec100 = *( (Float *) ( (VEC_IMM00) + (vec1ptr) ) ) ; 
  vec101 = *( (Float *) ( (VEC_IMM01) + (vec1ptr) ) ) ; 
  vec100 = A * vec100  ;
  vec101 = A * vec101  ;
  vec200 = *( (Float *) ( (VEC_IMM00) + (vec2ptr) ) ) ; 
  vec201 = *( (Float *) ( (VEC_IMM01) + (vec2ptr) ) ) ; 
  vec210 = *( (Float *) ( (VEC_IMM10) + (vec2ptr) ) ) ; 
  vec211 = *( (Float *) ( (VEC_IMM11) + (vec2ptr) ) ) ; 
  vec220 = *( (Float *) ( (VEC_IMM20) + (vec2ptr) ) ) ; 
  vec221 = *( (Float *) ( (VEC_IMM21) + (vec2ptr) ) ) ; 
cross_over_lin_s_lab1:
/*pragma_dcbt_space 5*/
  outptr = matptr + ZERO ; 
  vec1tmp = vec1ptr + ZERO ; 
  vec2tmp = vec2ptr + ZERO ; 
  vec110 = *( (Float *) ( (VEC_IMM10) + (vec1tmp) ) ) ; 
  vec111 = *( (Float *) ( (VEC_IMM11) + (vec1tmp) ) ) ; 
  mat000 = vec200 * vec100  ;
  mat001 = vec200 * vec101  ;
  mat000 = vec201 * vec101 + mat000 ;
  mat001 = - ( vec201 * vec100 - mat001 );
  *( (Float *) ( (MAT_IMM000) + (outptr) ) ) = mat000; 
  *( (Float *) ( (MAT_IMM001) + (outptr) ) ) = mat001; 
  mat010 = vec210 * vec100  ;
  mat011 = vec210 * vec101  ;
  mat010 = vec211 * vec101 + mat010 ;
  mat011 = - ( vec211 * vec100 - mat011 );
  *( (Float *) ( (MAT_IMM010) + (outptr) ) ) = mat010; 
  *( (Float *) ( (MAT_IMM011) + (outptr) ) ) = mat011; 
  mat020 = vec220 * vec100  ;
  mat021 = vec220 * vec101  ;
  mat020 = vec221 * vec101 + mat020 ;
  mat021 = - ( vec221 * vec100 - mat021 );
  src1_idx =   ( Isize ) + ( src1_idx )  ;
  prefscratch = * ( (unsigned long *) ( (PRE0)+(src1_idx) ) )  ;
  prefscratch = prefscratch * VEC_ATOM ; 
  vec1ptr = LOOKUP0 + prefscratch ; 
  prefscratch = * ( (unsigned long *) ( (SHFT_PREF_IMM)+(src1_idx) ) )  ;
  prefscratch = prefscratch * VEC_ATOM ; 
  prefscratch = LOOKUP0 + prefscratch ; 
  if ( predicate >= 0 ) goto cross_over_lin_s_lab2;
  prefscratch = StackPointer | StackPointer ; 
cross_over_lin_s_lab2:
/*preload PRER0 prefscratch */
/*preload PRER1 src1_idx */
  *( (Float *) ( (MAT_IMM020) + (outptr) ) ) = mat020; 
  *( (Float *) ( (MAT_IMM021) + (outptr) ) ) = mat021; 
  vec110 = A * vec110  ;
  vec111 = A * vec111  ;
  vec120 = *( (Float *) ( (VEC_IMM20) + (vec1tmp) ) ) ; 
  vec121 = *( (Float *) ( (VEC_IMM21) + (vec1tmp) ) ) ; 
  vec120 = A * vec120  ;
  vec121 = A * vec121  ;
  mat100 = vec200 * vec110  ;
  mat101 = vec200 * vec111  ;
  mat100 = vec201 * vec111 + mat100 ;
  mat101 = - ( vec201 * vec110 - mat101 );
  vec2ptr =   ( VEC_ATOM ) + ( vec2ptr )  ;
/*preload PRER0 vec2ptr */
  *( (Float *) ( (MAT_IMM100) + (outptr) ) ) = mat100; 
  *( (Float *) ( (MAT_IMM101) + (outptr) ) ) = mat101; 
  dest_idx =   ( Isize ) + ( dest_idx )  ;
  prefscratch = * ( (unsigned long *) ( (PRE0)+(dest_idx) ) )  ;
  prefscratch = prefscratch * MAT_ATOM ; 
  matptr = LOOKUP1 + prefscratch ; 
  prefscratch = * ( (unsigned long *) ( (SHFT_PREF_IMM)+(dest_idx) ) )  ;
  prefscratch = prefscratch * MAT_ATOM ; 
  prefscratch = LOOKUP1 + prefscratch ; 
  if ( predicate >= 0 ) goto cross_over_lin_s_lab3;
  prefscratch = StackPointer | StackPointer ; 
cross_over_lin_s_lab3:
/*preload PRER0 prefscratch */
/*preload PRER1 prefscratch */
/*preload PRER2 prefscratch */
/*preload PRER1 dest_idx */
  mat110 = vec210 * vec110  ;
  mat111 = vec210 * vec111  ;
  mat110 = vec211 * vec111 + mat110 ;
  mat111 = - ( vec211 * vec110 - mat111 );
  *( (Float *) ( (MAT_IMM110) + (outptr) ) ) = mat110; 
  *( (Float *) ( (MAT_IMM111) + (outptr) ) ) = mat111; 
  mat120 = vec220 * vec110  ;
  mat121 = vec220 * vec111  ;
  mat120 = vec221 * vec111 + mat120 ;
  mat121 = - ( vec221 * vec110 - mat121 );
  *( (Float *) ( (MAT_IMM120) + (outptr) ) ) = mat120; 
  *( (Float *) ( (MAT_IMM121) + (outptr) ) ) = mat121; 
  vec100 = *( (Float *) ( (VEC_IMM00) + (vec1ptr) ) ) ; 
  vec101 = *( (Float *) ( (VEC_IMM01) + (vec1ptr) ) ) ; 
  mat200 = vec200 * vec120  ;
  mat201 = vec200 * vec121  ;
  mat200 = vec201 * vec121 + mat200 ;
  mat201 = - ( vec201 * vec120 - mat201 );
  *( (Float *) ( (MAT_IMM200) + (outptr) ) ) = mat200; 
  *( (Float *) ( (MAT_IMM201) + (outptr) ) ) = mat201; 
  mat210 = vec210 * vec120  ;
  mat211 = vec210 * vec121  ;
  mat210 = vec211 * vec121 + mat210 ;
  mat211 = - ( vec211 * vec120 - mat211 );
  *( (Float *) ( (MAT_IMM210) + (outptr) ) ) = mat210; 
  *( (Float *) ( (MAT_IMM211) + (outptr) ) ) = mat211; 
  mat220 = vec220 * vec120  ;
  mat221 = vec220 * vec121  ;
  mat220 = vec221 * vec121 + mat220 ;
  mat221 = - ( vec221 * vec120 - mat221 );
  *( (Float *) ( (MAT_IMM220) + (outptr) ) ) = mat220; 
  *( (Float *) ( (MAT_IMM221) + (outptr) ) ) = mat221; 
  vec100 = A * vec100  ;
  vec101 = A * vec101  ;
  vec200 = *( (Float *) ( (VEC_IMM00) + (vec2ptr) ) ) ; 
  vec201 = *( (Float *) ( (VEC_IMM01) + (vec2ptr) ) ) ; 
  vec210 = *( (Float *) ( (VEC_IMM10) + (vec2ptr) ) ) ; 
  vec211 = *( (Float *) ( (VEC_IMM11) + (vec2ptr) ) ) ; 
  vec220 = *( (Float *) ( (VEC_IMM20) + (vec2ptr) ) ) ; 
  vec221 = *( (Float *) ( (VEC_IMM21) + (vec2ptr) ) ) ; 
  counter =   ( minus1 ) + ( counter )  ;
  if ( counter >  0 ) goto cross_over_lin_s_lab1;
cross_over_lin_s_lab0:

  return;
}
#ifdef __cplusplus
}
#endif
