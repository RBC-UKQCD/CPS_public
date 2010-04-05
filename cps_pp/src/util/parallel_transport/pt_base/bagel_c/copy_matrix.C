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
#define UIMM000 0 
#define UIMM001 8 
#define UIMM010 16 
#define UIMM011 24 
#define UIMM020 32 
#define UIMM021 40 
#define UIMM100 48 
#define UIMM101 56 
#define UIMM110 64 
#define UIMM111 72 
#define UIMM120 80 
#define UIMM121 88 
#define UIMM200 96 
#define UIMM201 104 
#define UIMM210 112 
#define UIMM211 120 
#define UIMM220 128 
#define UIMM221 136 
#define ZERO_IMM 0 
#define MAT_IMM 144 
#define bias 0 
#define PRE0 0 
#define PRE1 32 
#define PRE2 64 
#define PRE3 96 
#define PRE4 128 
#define minus1 -1 

typedef double Float;

#ifdef __cplusplus
extern "C" {
#endif
 void copy_matrix (unsigned long Arg0,unsigned long Arg1,unsigned long Arg2,unsigned long Arg3,unsigned long Arg4)
{

  /*Ensure natural alignment*/  int iStackArea[0+1];

  char *StackArea=(char *)iStackArea;

  unsigned long predicate;
  unsigned long prefscratch;
  unsigned long Uptr;
  unsigned long Vptr;
  unsigned long counter;
  unsigned long lenptr;
  unsigned long U_idx;
  unsigned long V_idx;
  unsigned long LOOKUP0;
  unsigned long PRER0;
  unsigned long PRER1;
  unsigned long PRER2;
  unsigned long PRER3;
  unsigned long PRER4;
  unsigned long LOOKUP1;
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
  U_idx = Arg4 | Arg4 ; 
  counter = * ( (unsigned long *) ( (ZERO_IMM)+(lenptr) ) )  ;
  if ( counter <= 0 ) goto copy_matrix_lab0;
  LOOKUP0 = Uptr | Uptr ; 
PRER0=PRE0;
PRER1=PRE1;
PRER2=PRE2;
PRER3=PRE3;
PRER4=PRE4;
  prefscratch = * ( (unsigned long *) ( (PRE0)+(U_idx) ) )  ;
  prefscratch = prefscratch * MAT_IMM ; 
  Uptr = LOOKUP0 + prefscratch ; 
  LOOKUP1 = Vptr | Vptr ; 
  prefscratch = * ( (unsigned long *) ( (PRE0)+(V_idx) ) )  ;
  prefscratch = prefscratch * MAT_IMM ; 
  Vptr = LOOKUP1 + prefscratch ; 
/*pragma_store_lim 1*/
/*pragma_dcbt_post 2*/
copy_matrix_lab1:
  U000 = *( (Float *) ( (UIMM000) + (Uptr) ) ) ; 
  U001 = *( (Float *) ( (UIMM001) + (Uptr) ) ) ; 
  U010 = *( (Float *) ( (UIMM010) + (Uptr) ) ) ; 
  U011 = *( (Float *) ( (UIMM011) + (Uptr) ) ) ; 
  U020 = *( (Float *) ( (UIMM020) + (Uptr) ) ) ; 
  U021 = *( (Float *) ( (UIMM021) + (Uptr) ) ) ; 
  prefscratch = * ( (unsigned long *) ( (SHFT_PREF_IMM)+(U_idx) ) )  ;
  prefscratch = prefscratch * MAT_IMM ; 
  prefscratch = LOOKUP0 + prefscratch ; 
  if ( predicate >= 0 ) goto copy_matrix_lab2;
  prefscratch = StackPointer | StackPointer ; 
copy_matrix_lab2:
/*preload PRER0 prefscratch */
/*preload PRER1 prefscratch */
/*preload PRER2 prefscratch */
/*preload PRER3 prefscratch */
/*preload PRER4 prefscratch */
/*preload PRER1 U_idx */
  *( (Float *) ( (UIMM000) + (Vptr) ) ) = U000; 
  *( (Float *) ( (UIMM001) + (Vptr) ) ) = U001; 
  *( (Float *) ( (UIMM010) + (Vptr) ) ) = U010; 
  *( (Float *) ( (UIMM011) + (Vptr) ) ) = U011; 
  *( (Float *) ( (UIMM020) + (Vptr) ) ) = U020; 
  *( (Float *) ( (UIMM021) + (Vptr) ) ) = U021; 
  U100 = *( (Float *) ( (UIMM100) + (Uptr) ) ) ; 
  U101 = *( (Float *) ( (UIMM101) + (Uptr) ) ) ; 
  U110 = *( (Float *) ( (UIMM110) + (Uptr) ) ) ; 
  U111 = *( (Float *) ( (UIMM111) + (Uptr) ) ) ; 
  U120 = *( (Float *) ( (UIMM120) + (Uptr) ) ) ; 
  U121 = *( (Float *) ( (UIMM121) + (Uptr) ) ) ; 
  *( (Float *) ( (UIMM100) + (Vptr) ) ) = U100; 
  *( (Float *) ( (UIMM101) + (Vptr) ) ) = U101; 
  *( (Float *) ( (UIMM110) + (Vptr) ) ) = U110; 
  *( (Float *) ( (UIMM111) + (Vptr) ) ) = U111; 
  *( (Float *) ( (UIMM120) + (Vptr) ) ) = U120; 
  *( (Float *) ( (UIMM121) + (Vptr) ) ) = U121; 
  U200 = *( (Float *) ( (UIMM200) + (Uptr) ) ) ; 
  U201 = *( (Float *) ( (UIMM201) + (Uptr) ) ) ; 
  U210 = *( (Float *) ( (UIMM210) + (Uptr) ) ) ; 
  U211 = *( (Float *) ( (UIMM211) + (Uptr) ) ) ; 
  U220 = *( (Float *) ( (UIMM220) + (Uptr) ) ) ; 
  U221 = *( (Float *) ( (UIMM221) + (Uptr) ) ) ; 
  *( (Float *) ( (UIMM200) + (Vptr) ) ) = U200; 
  *( (Float *) ( (UIMM201) + (Vptr) ) ) = U201; 
  *( (Float *) ( (UIMM210) + (Vptr) ) ) = U210; 
  *( (Float *) ( (UIMM211) + (Vptr) ) ) = U211; 
  *( (Float *) ( (UIMM220) + (Vptr) ) ) = U220; 
  *( (Float *) ( (UIMM221) + (Vptr) ) ) = U221; 
  V_idx =   ( Isize ) + ( V_idx )  ;
  prefscratch = * ( (unsigned long *) ( (PRE0)+(V_idx) ) )  ;
  prefscratch = prefscratch * MAT_IMM ; 
  Vptr = LOOKUP1 + prefscratch ; 
  U_idx =   ( Isize ) + ( U_idx )  ;
  prefscratch = * ( (unsigned long *) ( (PRE0)+(U_idx) ) )  ;
  prefscratch = prefscratch * MAT_IMM ; 
  Uptr = LOOKUP0 + prefscratch ; 
  counter =   ( minus1 ) + ( counter )  ;
  if ( counter >  0 ) goto copy_matrix_lab1;
copy_matrix_lab0:

  return;
}
#ifdef __cplusplus
}
#endif
