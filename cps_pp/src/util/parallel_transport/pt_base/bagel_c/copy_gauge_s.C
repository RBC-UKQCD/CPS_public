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
#define UIMM0 0 
#define UIMM1 4 
#define UIMM2 8 
#define UIMM3 12 
#define UIMM4 16 
#define UIMM5 20 
#define UIMM6 24 
#define UIMM7 28 
#define UIMM8 32 
#define UIMM9 36 
#define UIMM10 40 
#define UIMM11 44 
#define UIMM12 48 
#define UIMM13 52 
#define UIMM14 56 
#define UIMM15 60 
#define UIMM16 64 
#define UIMM17 68 
#define UIMM18 72 
#define UIMM19 76 
#define ZERO_IMM 0 
#define MAT_IMM 72 
#define GAUGE_AGG_IMM 80 
#define bias 0 
#define PRE0 0 
#define PRE1 32 
#define PRE2 64 
#define PRE3 96 
#define PRE4 128 
#define minus1 -1 

typedef float Float;

#ifdef __cplusplus
extern "C" {
#endif
 void copy_gauge_s (unsigned long Arg0,unsigned long Arg1,unsigned long Arg2,unsigned long Arg3)
{

  /*Ensure natural alignment*/  int iStackArea[0+1];

  char *StackArea=(char *)iStackArea;

  unsigned long predicate;
  unsigned long prefscratch;
  unsigned long Uptr;
  unsigned long Vptr;
  unsigned long counter;
  unsigned long lenptr;
  unsigned long V_idx;
  unsigned long LOOKUP0;
  unsigned long PRER0;
  unsigned long PRER1;
  unsigned long PRER2;
  unsigned long PRER3;
  unsigned long PRER4;
  unsigned long R17;
  unsigned long R18;
  unsigned long R19;
  unsigned long R20;
  unsigned long R21;
  unsigned long R22;
  unsigned long R23;
  unsigned long R24;
  unsigned long R25;
  unsigned long R26;
  unsigned long R27;
  unsigned long R28;
  unsigned long R29;
  unsigned long R31;
  unsigned long StackPointer = (unsigned long) StackArea ;
  Float U000;
  Float U001;
  Float U010;
  Float U011;
  Float U020;
  Float U021;
  Float U100;
  Float U101;
  Float U110;
  Float U111;
  Float U120;
  Float U121;
  Float U200;
  Float U201;
  Float U210;
  Float U211;
  Float U220;
  Float U221;
  Float F18;
  Float F19;
  Float F20;
  Float F21;
  Float F22;
  Float F23;
  Float F24;
  Float F25;
  Float F26;
  Float F27;
  Float F28;
  Float F29;
  Float F30;
  Float F31;
  Vptr = Arg0 | Arg0 ; 
  Uptr = Arg1 | Arg1 ; 
  lenptr = Arg2 | Arg2 ; 
  V_idx = Arg3 | Arg3 ; 
  counter = * ( (unsigned long *) ( (ZERO_IMM)+(lenptr) ) )  ;
  if ( counter <= 0 ) goto copy_gauge_s_lab0;
  LOOKUP0 = Vptr | Vptr ; 
PRER0=PRE0;
PRER1=PRE1;
PRER2=PRE2;
  prefscratch = * ( (unsigned long *) ( (PRE0)+(V_idx) ) )  ;
  prefscratch = prefscratch * MAT_IMM ; 
  Vptr = LOOKUP0 + prefscratch ; 
PRER3=PRE3;
PRER4=PRE4;
/*pragma_store_lim 1*/
/*pragma_dcbt_post 2*/
copy_gauge_s_lab1:
  U000 = *( (Float *) ( (UIMM2) + (Uptr) ) ) ; 
  U001 = *( (Float *) ( (UIMM3) + (Uptr) ) ) ; 
  U010 = *( (Float *) ( (UIMM4) + (Uptr) ) ) ; 
  U011 = *( (Float *) ( (UIMM5) + (Uptr) ) ) ; 
  U020 = *( (Float *) ( (UIMM6) + (Uptr) ) ) ; 
  U021 = *( (Float *) ( (UIMM7) + (Uptr) ) ) ; 
/*preload PRER2 Uptr */
/*preload PRER3 Uptr */
/*preload PRER4 Uptr */
  *( (Float *) ( (UIMM0) + (Vptr) ) ) = U000; 
  *( (Float *) ( (UIMM1) + (Vptr) ) ) = U001; 
  *( (Float *) ( (UIMM2) + (Vptr) ) ) = U010; 
  *( (Float *) ( (UIMM3) + (Vptr) ) ) = U011; 
  *( (Float *) ( (UIMM4) + (Vptr) ) ) = U020; 
  *( (Float *) ( (UIMM5) + (Vptr) ) ) = U021; 
  U100 = *( (Float *) ( (UIMM8) + (Uptr) ) ) ; 
  U101 = *( (Float *) ( (UIMM9) + (Uptr) ) ) ; 
  U110 = *( (Float *) ( (UIMM10) + (Uptr) ) ) ; 
  U111 = *( (Float *) ( (UIMM11) + (Uptr) ) ) ; 
  U120 = *( (Float *) ( (UIMM12) + (Uptr) ) ) ; 
  U121 = *( (Float *) ( (UIMM13) + (Uptr) ) ) ; 
  *( (Float *) ( (UIMM6) + (Vptr) ) ) = U100; 
  *( (Float *) ( (UIMM7) + (Vptr) ) ) = U101; 
  *( (Float *) ( (UIMM8) + (Vptr) ) ) = U110; 
  *( (Float *) ( (UIMM9) + (Vptr) ) ) = U111; 
  *( (Float *) ( (UIMM10) + (Vptr) ) ) = U120; 
  *( (Float *) ( (UIMM11) + (Vptr) ) ) = U121; 
  U200 = *( (Float *) ( (UIMM14) + (Uptr) ) ) ; 
  U201 = *( (Float *) ( (UIMM15) + (Uptr) ) ) ; 
  U210 = *( (Float *) ( (UIMM16) + (Uptr) ) ) ; 
  U211 = *( (Float *) ( (UIMM17) + (Uptr) ) ) ; 
  U220 = *( (Float *) ( (UIMM18) + (Uptr) ) ) ; 
  U221 = *( (Float *) ( (UIMM19) + (Uptr) ) ) ; 
  *( (Float *) ( (UIMM12) + (Vptr) ) ) = U200; 
  *( (Float *) ( (UIMM13) + (Vptr) ) ) = U201; 
  *( (Float *) ( (UIMM14) + (Vptr) ) ) = U210; 
  *( (Float *) ( (UIMM15) + (Vptr) ) ) = U211; 
  *( (Float *) ( (UIMM16) + (Vptr) ) ) = U220; 
  *( (Float *) ( (UIMM17) + (Vptr) ) ) = U221; 
  V_idx =   ( Isize ) + ( V_idx )  ;
  prefscratch = * ( (unsigned long *) ( (PRE0)+(V_idx) ) )  ;
  prefscratch = prefscratch * MAT_IMM ; 
  Vptr = LOOKUP0 + prefscratch ; 
  Uptr =   ( GAUGE_AGG_IMM ) + ( Uptr )  ;
  counter =   ( minus1 ) + ( counter )  ;
  if ( counter >  0 ) goto copy_gauge_s_lab1;
copy_gauge_s_lab0:

  return;
}
#ifdef __cplusplus
}
#endif
