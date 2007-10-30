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
#define UIMM00 0 
#define UIMM01 8 
#define UIMM10 16 
#define UIMM11 24 
#define UIMM20 32 
#define UIMM21 40 
#define VIMM00 48 
#define VIMM01 56 
#define VIMM10 64 
#define VIMM11 72 
#define VIMM20 80 
#define VIMM21 88 
#define WIMM00 96 
#define WIMM01 104 
#define WIMM10 112 
#define WIMM11 120 
#define WIMM20 128 
#define WIMM21 136 
#define XIMM00 144 
#define XIMM01 152 
#define XIMM10 160 
#define XIMM11 168 
#define XIMM20 176 
#define XIMM21 184 
#define GIMM000 0 
#define GIMM001 8 
#define GIMM010 16 
#define GIMM011 24 
#define GIMM020 32 
#define GIMM021 40 
#define GIMM100 48 
#define GIMM101 56 
#define GIMM110 64 
#define GIMM111 72 
#define GIMM120 80 
#define GIMM121 88 
#define GIMM200 96 
#define GIMM201 104 
#define GIMM210 112 
#define GIMM211 120 
#define GIMM220 128 
#define GIMM221 136 
#define CHIIMM000 0 
#define CHIIMM001 4 
#define CHIIMM010 8 
#define CHIIMM011 12 
#define CHIIMM100 16 
#define CHIIMM101 20 
#define CHIIMM110 24 
#define CHIIMM111 28 
#define CHIIMM200 32 
#define CHIIMM201 36 
#define CHIIMM210 40 
#define CHIIMM211 44 
#define ZERO_IMM 0 
#define PSI_IMM 192 
#define CHI_IMM 48 
#define PAD_CHI_IMM 64 
#define MAT_IMM 144 
#define Ndim 4 
#define Ndimm1 3 
#define hIsize 4 
#define hbitbucket 64 
#define hstk0 112 
#define hstk1 160 
#define hstk2 208 
#define hstk3 256 
#define bias 0 
#define CONST0 0 
#define CONST1 32 
#define CONST2 64 
#define CONST3 96 
#define CONST4 128 
#define minus1 -1 
#define minus1 -1 

typedef double Float;

#ifdef __cplusplus
extern "C" {
#endif
 void s_dec_hsu3 (unsigned long Arg0,unsigned long Arg1,unsigned long Arg2,unsigned long Arg3)
{

  /*Ensure natural alignment*/  int iStackArea[304+1];

  char *StackArea=(char *)iStackArea;

  unsigned long predicate;
  unsigned long prefscratch;
  unsigned long psi;
  unsigned long Umu;
  unsigned long Chiin;
  unsigned long Chiout;
  unsigned long mem;
  unsigned long length;
  unsigned long Chiplus;
  unsigned long Chitmp1;
  unsigned long Chitmp2;
  unsigned long Chiminus;
  unsigned long pref;
  unsigned long mu;
  unsigned long tab;
  unsigned long CONSTR0;
  unsigned long CONSTR1;
  unsigned long CONSTR2;
  unsigned long CONSTR3;
  unsigned long CONSTR4;
  unsigned long R24;
  unsigned long R25;
  unsigned long R26;
  unsigned long R27;
  unsigned long R28;
  unsigned long R29;
  unsigned long R31;
  unsigned long StackPointer = (unsigned long) StackArea ;
  Float U00;
  Float U01;
  Float U10;
  Float U11;
  Float U20;
  Float U21;
  Float V00;
  Float V01;
  Float V10;
  Float V11;
  Float V20;
  Float V21;
  Float W00;
  Float W01;
  Float W10;
  Float W11;
  Float W20;
  Float W21;
  Float X00;
  Float X01;
  Float X10;
  Float X11;
  Float X20;
  Float X21;
  Float C0;
  Float C1;
  Float D0;
  Float D1;
  Float E0;
  Float E1;
  Float F0;
  Float F1;
  Float F32;
  Float F33;
  Float F34;
  Float F35;
  Float F36;
  Float F37;
  Float F38;
  Float F39;
  Float F40;
  Float F41;
  Float F42;
  Float F43;
  Float F44;
  Float F45;
  Float F46;
  Float F47;
  Float F48;
  Float F49;
  Float F50;
  Float F51;
  Float F52;
  Float F53;
  Float F54;
  Float F55;
  Float F56;
  Float F57;
  Float F58;
  Float F59;
  Float F60;
  Float F61;
  Float F62;
  Float F63;
  mem = StackPointer + bias ; 
CONSTR0=CONST0;
CONSTR1=CONST1;
CONSTR2=CONST2;
CONSTR3=CONST3;
CONSTR4=CONST4;
  psi = Arg0 | Arg0 ; 
  Umu = Arg1 | Arg1 ; 
  length = Arg2 | Arg2 ; 
  tab = Arg3 | Arg3 ; 
  length = * ( (unsigned long *) ( (ZERO_IMM)+(length) ) )  ;
  if ( length <= 0 ) goto s_dec_hsu3_lab0;
  Chiout = * ( (unsigned long *) ( (zeromulti)+(tab) ) )  ;
  mem = StackPointer + bias ; 
  U00 = *( (Float *) ( (UIMM00) + (psi) ) ) ; 
  U01 = *( (Float *) ( (UIMM01) + (psi) ) ) ; 
  V00 = *( (Float *) ( (VIMM00) + (psi) ) ) ; 
  V01 = *( (Float *) ( (VIMM01) + (psi) ) ) ; 
  W00 = *( (Float *) ( (WIMM00) + (psi) ) ) ; 
  W01 = *( (Float *) ( (WIMM01) + (psi) ) ) ; 
  X00 = *( (Float *) ( (XIMM00) + (psi) ) ) ; 
  X01 = *( (Float *) ( (XIMM01) + (psi) ) ) ; 
  U10 = *( (Float *) ( (UIMM10) + (psi) ) ) ; 
  U11 = *( (Float *) ( (UIMM11) + (psi) ) ) ; 
  V10 = *( (Float *) ( (VIMM10) + (psi) ) ) ; 
  V11 = *( (Float *) ( (VIMM11) + (psi) ) ) ; 
  W10 = *( (Float *) ( (WIMM10) + (psi) ) ) ; 
  W11 = *( (Float *) ( (WIMM11) + (psi) ) ) ; 
  X10 = *( (Float *) ( (XIMM10) + (psi) ) ) ; 
  X11 = *( (Float *) ( (XIMM11) + (psi) ) ) ; 
s_dec_hsu3_lab1:
/*pragma_load_lim 1*/
/*pragma_dcbt_space 8*/
  U20 = *( (Float *) ( (UIMM20) + (psi) ) ) ; 
  U21 = *( (Float *) ( (UIMM21) + (psi) ) ) ; 
  V20 = *( (Float *) ( (VIMM20) + (psi) ) ) ; 
  V21 = *( (Float *) ( (VIMM21) + (psi) ) ) ; 
  W20 = *( (Float *) ( (WIMM20) + (psi) ) ) ; 
  W21 = *( (Float *) ( (WIMM21) + (psi) ) ) ; 
  X20 = *( (Float *) ( (XIMM20) + (psi) ) ) ; 
  X21 = *( (Float *) ( (XIMM21) + (psi) ) ) ; 
  Chiplus = mem + hstk0 ; 
  Chiminus = Chiout + ZERO_IMM ; 
  tab =   ( Isize ) + ( tab )  ;
  Chiout = * ( (unsigned long *) ( (zeromulti)+(tab) ) )  ;
  Chitmp2 = mem + hstk1 ; 
  Chitmp1 = Chiout + ZERO_IMM ; 
  tab =   ( Isize ) + ( tab )  ;
  Chiout = * ( (unsigned long *) ( (zeromulti)+(tab) ) )  ;
  C0 = U00 + X01 ; 
  *( (float *) ( (CHIIMM000) + (Chiminus) ) ) = C0; 
  C1 = U01 - X00 ; 
  *( (float *) ( (CHIIMM001) + (Chiminus) ) ) = C1; 
  D0 = U00 - X01 ; 
  *( (float *) ( (CHIIMM000) + (Chiplus) ) ) = D0; 
  D1 = U01 + X00 ; 
  *( (float *) ( (CHIIMM001) + (Chiplus) ) ) = D1; 
  E0 = V00 + W01 ; 
  *( (float *) ( (CHIIMM010) + (Chiminus) ) ) = E0; 
  E1 = V01 - W00 ; 
  *( (float *) ( (CHIIMM011) + (Chiminus) ) ) = E1; 
  F0 = V00 - W01 ; 
  *( (float *) ( (CHIIMM010) + (Chiplus) ) ) = F0; 
  F1 = V01 + W00 ; 
  *( (float *) ( (CHIIMM011) + (Chiplus) ) ) = F1; 
  C0 = U10 + X11 ; 
  *( (float *) ( (CHIIMM100) + (Chiminus) ) ) = C0; 
  C1 = U11 - X10 ; 
  *( (float *) ( (CHIIMM101) + (Chiminus) ) ) = C1; 
  D0 = U10 - X11 ; 
  *( (float *) ( (CHIIMM100) + (Chiplus) ) ) = D0; 
  D1 = U11 + X10 ; 
  *( (float *) ( (CHIIMM101) + (Chiplus) ) ) = D1; 
  E0 = V10 + W11 ; 
  *( (float *) ( (CHIIMM110) + (Chiminus) ) ) = E0; 
  E1 = V11 - W10 ; 
  *( (float *) ( (CHIIMM111) + (Chiminus) ) ) = E1; 
  F0 = V10 - W11 ; 
  *( (float *) ( (CHIIMM110) + (Chiplus) ) ) = F0; 
  F1 = V11 + W10 ; 
  *( (float *) ( (CHIIMM111) + (Chiplus) ) ) = F1; 
  C0 = U20 + X21 ; 
  *( (float *) ( (CHIIMM200) + (Chiminus) ) ) = C0; 
  C1 = U21 - X20 ; 
  *( (float *) ( (CHIIMM201) + (Chiminus) ) ) = C1; 
  D0 = U20 - X21 ; 
  *( (float *) ( (CHIIMM200) + (Chiplus) ) ) = D0; 
  D1 = U21 + X20 ; 
  *( (float *) ( (CHIIMM201) + (Chiplus) ) ) = D1; 
  E0 = V20 + W21 ; 
  *( (float *) ( (CHIIMM210) + (Chiminus) ) ) = E0; 
  E1 = V21 - W20 ; 
  *( (float *) ( (CHIIMM211) + (Chiminus) ) ) = E1; 
  F0 = V20 - W21 ; 
  *( (float *) ( (CHIIMM210) + (Chiplus) ) ) = F0; 
  F1 = V21 + W20 ; 
  *( (float *) ( (CHIIMM211) + (Chiplus) ) ) = F1; 
  C0 = U00 + X00 ; 
  *( (float *) ( (CHIIMM000) + (Chitmp1) ) ) = C0; 
  C1 = U01 + X01 ; 
  *( (float *) ( (CHIIMM001) + (Chitmp1) ) ) = C1; 
  D0 = U00 - X00 ; 
  *( (float *) ( (CHIIMM000) + (Chitmp2) ) ) = D0; 
  D1 = U01 - X01 ; 
  *( (float *) ( (CHIIMM001) + (Chitmp2) ) ) = D1; 
  E0 = V00 - W00 ; 
  *( (float *) ( (CHIIMM010) + (Chitmp1) ) ) = E0; 
  E1 = V01 - W01 ; 
  *( (float *) ( (CHIIMM011) + (Chitmp1) ) ) = E1; 
  F0 = V00 + W00 ; 
  *( (float *) ( (CHIIMM010) + (Chitmp2) ) ) = F0; 
  F1 = V01 + W01 ; 
  *( (float *) ( (CHIIMM011) + (Chitmp2) ) ) = F1; 
  C0 = U10 + X10 ; 
  *( (float *) ( (CHIIMM100) + (Chitmp1) ) ) = C0; 
  C1 = U11 + X11 ; 
  *( (float *) ( (CHIIMM101) + (Chitmp1) ) ) = C1; 
  D0 = U10 - X10 ; 
  *( (float *) ( (CHIIMM100) + (Chitmp2) ) ) = D0; 
  D1 = U11 - X11 ; 
  *( (float *) ( (CHIIMM101) + (Chitmp2) ) ) = D1; 
  E0 = V10 - W10 ; 
  *( (float *) ( (CHIIMM110) + (Chitmp1) ) ) = E0; 
  E1 = V11 - W11 ; 
  *( (float *) ( (CHIIMM111) + (Chitmp1) ) ) = E1; 
  F0 = V10 + W10 ; 
  *( (float *) ( (CHIIMM110) + (Chitmp2) ) ) = F0; 
  F1 = V11 + W11 ; 
  *( (float *) ( (CHIIMM111) + (Chitmp2) ) ) = F1; 
  C0 = U20 + X20 ; 
  *( (float *) ( (CHIIMM200) + (Chitmp1) ) ) = C0; 
  C1 = U21 + X21 ; 
  *( (float *) ( (CHIIMM201) + (Chitmp1) ) ) = C1; 
  D0 = U20 - X20 ; 
  *( (float *) ( (CHIIMM200) + (Chitmp2) ) ) = D0; 
  D1 = U21 - X21 ; 
  *( (float *) ( (CHIIMM201) + (Chitmp2) ) ) = D1; 
  E0 = V20 - W20 ; 
  *( (float *) ( (CHIIMM210) + (Chitmp1) ) ) = E0; 
  E1 = V21 - W21 ; 
  *( (float *) ( (CHIIMM211) + (Chitmp1) ) ) = E1; 
  F0 = V20 + W20 ; 
  *( (float *) ( (CHIIMM210) + (Chitmp2) ) ) = F0; 
  F1 = V21 + W21 ; 
  *( (float *) ( (CHIIMM211) + (Chitmp2) ) ) = F1; 
/*pragma_dcbt_space 4*/
/*preload CONSTR0 Umu */
/*preload CONSTR1 Umu */
/*preload CONSTR2 Umu */
/*preload CONSTR3 Umu */
/*preload CONSTR4 Umu */
  Chiplus = mem + hstk2 ; 
  Chiminus = Chiout + ZERO_IMM ; 
  tab =   ( Isize ) + ( tab )  ;
  Chiout = * ( (unsigned long *) ( (zeromulti)+(tab) ) )  ;
  C0 = U00 + W01 ; 
  *( (float *) ( (CHIIMM000) + (Chiminus) ) ) = C0; 
  C1 = U01 - W00 ; 
  *( (float *) ( (CHIIMM001) + (Chiminus) ) ) = C1; 
  D0 = U00 - W01 ; 
  *( (float *) ( (CHIIMM000) + (Chiplus) ) ) = D0; 
  D1 = U01 + W00 ; 
  *( (float *) ( (CHIIMM001) + (Chiplus) ) ) = D1; 
  E0 = V00 - X01 ; 
  *( (float *) ( (CHIIMM010) + (Chiminus) ) ) = E0; 
  E1 = V01 + X00 ; 
  *( (float *) ( (CHIIMM011) + (Chiminus) ) ) = E1; 
  F0 = V00 + X01 ; 
  *( (float *) ( (CHIIMM010) + (Chiplus) ) ) = F0; 
  F1 = V01 - X00 ; 
  *( (float *) ( (CHIIMM011) + (Chiplus) ) ) = F1; 
  C0 = U10 + W11 ; 
  *( (float *) ( (CHIIMM100) + (Chiminus) ) ) = C0; 
  C1 = U11 - W10 ; 
  *( (float *) ( (CHIIMM101) + (Chiminus) ) ) = C1; 
  D0 = U10 - W11 ; 
  *( (float *) ( (CHIIMM100) + (Chiplus) ) ) = D0; 
  D1 = U11 + W10 ; 
  *( (float *) ( (CHIIMM101) + (Chiplus) ) ) = D1; 
  E0 = V10 - X11 ; 
  *( (float *) ( (CHIIMM110) + (Chiminus) ) ) = E0; 
  E1 = V11 + X10 ; 
  *( (float *) ( (CHIIMM111) + (Chiminus) ) ) = E1; 
  F0 = V10 + X11 ; 
  *( (float *) ( (CHIIMM110) + (Chiplus) ) ) = F0; 
  F1 = V11 - X10 ; 
  *( (float *) ( (CHIIMM111) + (Chiplus) ) ) = F1; 
  C0 = U20 + W21 ; 
  *( (float *) ( (CHIIMM200) + (Chiminus) ) ) = C0; 
  C1 = U21 - W20 ; 
  *( (float *) ( (CHIIMM201) + (Chiminus) ) ) = C1; 
  D0 = U20 - W21 ; 
  *( (float *) ( (CHIIMM200) + (Chiplus) ) ) = D0; 
  D1 = U21 + W20 ; 
  *( (float *) ( (CHIIMM201) + (Chiplus) ) ) = D1; 
  E0 = V20 - X21 ; 
  *( (float *) ( (CHIIMM210) + (Chiminus) ) ) = E0; 
  E1 = V21 + X20 ; 
  *( (float *) ( (CHIIMM211) + (Chiminus) ) ) = E1; 
  F0 = V20 + X21 ; 
  *( (float *) ( (CHIIMM210) + (Chiplus) ) ) = F0; 
  F1 = V21 - X20 ; 
  *( (float *) ( (CHIIMM211) + (Chiplus) ) ) = F1; 
  Chiplus = mem + hstk3 ; 
  Chiminus = Chiout + ZERO_IMM ; 
  tab =   ( Isize ) + ( tab )  ;
  Chiout = * ( (unsigned long *) ( (zeromulti)+(tab) ) )  ;
  C0 = U00 - W00 ; 
  *( (float *) ( (CHIIMM000) + (Chiminus) ) ) = C0; 
  C1 = U01 - W01 ; 
  *( (float *) ( (CHIIMM001) + (Chiminus) ) ) = C1; 
  D0 = U00 + W00 ; 
  *( (float *) ( (CHIIMM000) + (Chiplus) ) ) = D0; 
  D1 = U01 + W01 ; 
  *( (float *) ( (CHIIMM001) + (Chiplus) ) ) = D1; 
  E0 = V00 - X00 ; 
  *( (float *) ( (CHIIMM010) + (Chiminus) ) ) = E0; 
  E1 = V01 - X01 ; 
  *( (float *) ( (CHIIMM011) + (Chiminus) ) ) = E1; 
  F0 = V00 + X00 ; 
  *( (float *) ( (CHIIMM010) + (Chiplus) ) ) = F0; 
  F1 = V01 + X01 ; 
  *( (float *) ( (CHIIMM011) + (Chiplus) ) ) = F1; 
/*pragma_store_lim 1*/
  C0 = U10 - W10 ; 
  *( (float *) ( (CHIIMM100) + (Chiminus) ) ) = C0; 
  C1 = U11 - W11 ; 
  *( (float *) ( (CHIIMM101) + (Chiminus) ) ) = C1; 
  D0 = U10 + W10 ; 
  *( (float *) ( (CHIIMM100) + (Chiplus) ) ) = D0; 
  D1 = U11 + W11 ; 
  *( (float *) ( (CHIIMM101) + (Chiplus) ) ) = D1; 
  F0 = V10 + X10 ; 
  *( (float *) ( (CHIIMM110) + (Chiplus) ) ) = F0; 
  F1 = V11 + X11 ; 
  *( (float *) ( (CHIIMM111) + (Chiplus) ) ) = F1; 
  C0 = U20 + W20 ; 
  *( (float *) ( (CHIIMM200) + (Chiplus) ) ) = C0; 
  C1 = U21 + W21 ; 
  *( (float *) ( (CHIIMM201) + (Chiplus) ) ) = C1; 
  D0 = V20 + X20 ; 
  *( (float *) ( (CHIIMM210) + (Chiplus) ) ) = D0; 
  D1 = V21 + X21 ; 
  *( (float *) ( (CHIIMM211) + (Chiplus) ) ) = D1; 
  F0 = V10 - X10 ; 
  F1 = V11 - X11 ; 
  C0 = U20 - W20 ; 
  C1 = U21 - W21 ; 
  D0 = V20 - X20 ; 
  D1 = V21 - X21 ; 
mu=Ndimm1;
  Chiin = mem + hstk0 ; 
  Chitmp1 = Chiminus + ZERO_IMM ; 
/*pragma_load_lim 10*/
  U00 = *( (float *) ( (CHIIMM000) + (Chiin) ) ) ; 
  U10 = *( (float *) ( (CHIIMM010) + (Chiin) ) ) ; 
  U01 = *( (float *) ( (CHIIMM001) + (Chiin) ) ) ; 
  U11 = *( (float *) ( (CHIIMM011) + (Chiin) ) ) ; 
  W00 = *( (Float *) ( (GIMM000) + (Umu) ) ) ; 
  X00 = *( (Float *) ( (GIMM100) + (Umu) ) ) ; 
  W01 = *( (Float *) ( (GIMM001) + (Umu) ) ) ; 
  X01 = *( (Float *) ( (GIMM101) + (Umu) ) ) ; 
  U20 = *( (float *) ( (CHIIMM100) + (Chiin) ) ) ; 
  V00 = *( (float *) ( (CHIIMM110) + (Chiin) ) ) ; 
  U21 = *( (float *) ( (CHIIMM101) + (Chiin) ) ) ; 
  V01 = *( (float *) ( (CHIIMM111) + (Chiin) ) ) ; 
/*pragma_dcbt_space 4*/
/*pragma_dcbt_post 1*/
s_dec_hsu3_lab2:
  *( (float *) ( (CHIIMM110) + (Chitmp1) ) ) = F0; 
  *( (float *) ( (CHIIMM111) + (Chitmp1) ) ) = F1; 
  *( (float *) ( (CHIIMM200) + (Chitmp1) ) ) = C0; 
  *( (float *) ( (CHIIMM201) + (Chitmp1) ) ) = C1; 
  *( (float *) ( (CHIIMM210) + (Chitmp1) ) ) = D0; 
  *( (float *) ( (CHIIMM211) + (Chitmp1) ) ) = D1; 
  W10 = *( (Float *) ( (GIMM010) + (Umu) ) ) ; 
  X10 = *( (Float *) ( (GIMM110) + (Umu) ) ) ; 
  W11 = *( (Float *) ( (GIMM011) + (Umu) ) ) ; 
  X11 = *( (Float *) ( (GIMM111) + (Umu) ) ) ; 
/*pragma_load_lim 1*/
  V10 = *( (float *) ( (CHIIMM200) + (Chiin) ) ) ; 
  V20 = *( (float *) ( (CHIIMM210) + (Chiin) ) ) ; 
  V11 = *( (float *) ( (CHIIMM201) + (Chiin) ) ) ; 
  V21 = *( (float *) ( (CHIIMM211) + (Chiin) ) ) ; 
  W20 = *( (Float *) ( (GIMM020) + (Umu) ) ) ; 
  X20 = *( (Float *) ( (GIMM120) + (Umu) ) ) ; 
  W21 = *( (Float *) ( (GIMM021) + (Umu) ) ) ; 
  X21 = *( (Float *) ( (GIMM121) + (Umu) ) ) ; 
  C0 = W00 * U00  ;
  C1 = W00 * U01  ;
  D0 = W00 * U10  ;
  D1 = W00 * U11  ;
  E0 = X00 * U00  ;
  E1 = X00 * U01  ;
  C0 = W01 * U01 + C0 ;
  C1 = - ( W01 * U00 - C1 );
  D0 = W01 * U11 + D0 ;
  D1 = - ( W01 * U10 - D1 );
  E0 = X01 * U01 + E0 ;
  E1 = - ( X01 * U00 - E1 );
  C0 = W10 * U20 + C0 ;
  C1 = W10 * U21 + C1 ;
  D0 = W10 * V00 + D0 ;
  D1 = W10 * V01 + D1 ;
  E0 = X10 * U20 + E0 ;
  E1 = X10 * U21 + E1 ;
  C0 = W11 * U21 + C0 ;
  C1 = - ( W11 * U20 - C1 );
  D0 = W11 * V01 + D0 ;
  D1 = - ( W11 * V00 - D1 );
  E0 = X11 * U21 + E0 ;
  E1 = - ( X11 * U20 - E1 );
  C0 = W20 * V10 + C0 ;
  C1 = W20 * V11 + C1 ;
  D0 = W20 * V20 + D0 ;
  D1 = W20 * V21 + D1 ;
  E0 = X20 * V10 + E0 ;
  E1 = X20 * V11 + E1 ;
  C0 = W21 * V11 + C0 ;
  C1 = - ( W21 * V10 - C1 );
  D0 = W21 * V21 + D0 ;
  D1 = - ( W21 * V20 - D1 );
  E0 = X21 * V11 + E0 ;
  E1 = - ( X21 * V10 - E1 );
  *( (float *) ( (CHIIMM000) + (Chiout) ) ) = C0; 
  *( (float *) ( (CHIIMM001) + (Chiout) ) ) = C1; 
  *( (float *) ( (CHIIMM010) + (Chiout) ) ) = D0; 
  *( (float *) ( (CHIIMM011) + (Chiout) ) ) = D1; 
  *( (float *) ( (CHIIMM100) + (Chiout) ) ) = E0; 
  *( (float *) ( (CHIIMM101) + (Chiout) ) ) = E1; 
/*pragma_dcbt_space 4*/
  W00 = *( (Float *) ( (GIMM200) + (Umu) ) ) ; 
  W01 = *( (Float *) ( (GIMM201) + (Umu) ) ) ; 
  W10 = *( (Float *) ( (GIMM210) + (Umu) ) ) ; 
  W11 = *( (Float *) ( (GIMM211) + (Umu) ) ) ; 
  W20 = *( (Float *) ( (GIMM220) + (Umu) ) ) ; 
  W21 = *( (Float *) ( (GIMM221) + (Umu) ) ) ; 
  pref = Umu + MAT_IMM ; 
/*preload CONSTR0 pref */
/*preload CONSTR1 pref */
/*preload CONSTR2 pref */
  pref = pref + CHI_IMM ; 
/*preload CONSTR0 pref */
/*preload CONSTR1 pref */
  Umu = Umu + MAT_IMM ; 
  F0 = X00 * U10  ;
  F1 = X00 * U11  ;
  C0 = W00 * U00  ;
  C1 = W00 * U01  ;
  D0 = W00 * U10  ;
  D1 = W00 * U11  ;
  F0 = X01 * U11 + F0 ;
  F1 = - ( X01 * U10 - F1 );
  C0 = W01 * U01 + C0 ;
  C1 = - ( W01 * U00 - C1 );
  D0 = W01 * U11 + D0 ;
  D1 = - ( W01 * U10 - D1 );
  F0 = X10 * V00 + F0 ;
  F1 = X10 * V01 + F1 ;
  C0 = W10 * U20 + C0 ;
  C1 = W10 * U21 + C1 ;
  D0 = W10 * V00 + D0 ;
  D1 = W10 * V01 + D1 ;
  F0 = X11 * V01 + F0 ;
  F1 = - ( X11 * V00 - F1 );
  C0 = W11 * U21 + C0 ;
  C1 = - ( W11 * U20 - C1 );
  D0 = W11 * V01 + D0 ;
  D1 = - ( W11 * V00 - D1 );
  F0 = X20 * V20 + F0 ;
  F1 = X20 * V21 + F1 ;
  C0 = W20 * V10 + C0 ;
  C1 = W20 * V11 + C1 ;
  D0 = W20 * V20 + D0 ;
  D1 = W20 * V21 + D1 ;
  F0 = X21 * V21 + F0 ;
  F1 = - ( X21 * V20 - F1 );
  C0 = W21 * V11 + C0 ;
  C1 = - ( W21 * V10 - C1 );
  D0 = W21 * V21 + D0 ;
  D1 = - ( W21 * V20 - D1 );
  Chitmp1 = Chiout + ZERO_IMM ; 
  Chiin = Chiin + CHI_IMM ; 
  tab =   ( Isize ) + ( tab )  ;
  Chiout = * ( (unsigned long *) ( (zeromulti)+(tab) ) )  ;
  U00 = *( (float *) ( (CHIIMM000) + (Chiin) ) ) ; 
  U10 = *( (float *) ( (CHIIMM010) + (Chiin) ) ) ; 
  U01 = *( (float *) ( (CHIIMM001) + (Chiin) ) ) ; 
  U11 = *( (float *) ( (CHIIMM011) + (Chiin) ) ) ; 
  W00 = *( (Float *) ( (GIMM000) + (Umu) ) ) ; 
  X00 = *( (Float *) ( (GIMM100) + (Umu) ) ) ; 
  W01 = *( (Float *) ( (GIMM001) + (Umu) ) ) ; 
  X01 = *( (Float *) ( (GIMM101) + (Umu) ) ) ; 
  U20 = *( (float *) ( (CHIIMM100) + (Chiin) ) ) ; 
  V00 = *( (float *) ( (CHIIMM110) + (Chiin) ) ) ; 
  U21 = *( (float *) ( (CHIIMM101) + (Chiin) ) ) ; 
  V01 = *( (float *) ( (CHIIMM111) + (Chiin) ) ) ; 
  mu =   ( minus1 ) + ( mu )  ;
  if ( mu >  0 ) goto s_dec_hsu3_lab2;
s_dec_hsu3_lab3:
/*pragma_dcbt_space 4*/
/*pragma_dcbt_post 1*/
  *( (float *) ( (CHIIMM110) + (Chitmp1) ) ) = F0; 
  *( (float *) ( (CHIIMM111) + (Chitmp1) ) ) = F1; 
  *( (float *) ( (CHIIMM200) + (Chitmp1) ) ) = C0; 
  *( (float *) ( (CHIIMM201) + (Chitmp1) ) ) = C1; 
  *( (float *) ( (CHIIMM210) + (Chitmp1) ) ) = D0; 
  *( (float *) ( (CHIIMM211) + (Chitmp1) ) ) = D1; 
  W10 = *( (Float *) ( (GIMM010) + (Umu) ) ) ; 
  X10 = *( (Float *) ( (GIMM110) + (Umu) ) ) ; 
  W11 = *( (Float *) ( (GIMM011) + (Umu) ) ) ; 
  X11 = *( (Float *) ( (GIMM111) + (Umu) ) ) ; 
/*pragma_dcbt_space 4*/
/*pragma_dcbt_post 0*/
  psi = psi + PSI_IMM ; 
/*preload CONSTR0 psi */
/*preload CONSTR1 psi */
/*preload CONSTR2 psi */
  pref = psi + CHI_IMM ; 
/*preload CONSTR0 pref */
/*preload CONSTR1 pref */
/*preload CONSTR2 pref */
/*pragma_load_lim 1*/
  V10 = *( (float *) ( (CHIIMM200) + (Chiin) ) ) ; 
  V20 = *( (float *) ( (CHIIMM210) + (Chiin) ) ) ; 
  V11 = *( (float *) ( (CHIIMM201) + (Chiin) ) ) ; 
  V21 = *( (float *) ( (CHIIMM211) + (Chiin) ) ) ; 
  W20 = *( (Float *) ( (GIMM020) + (Umu) ) ) ; 
  X20 = *( (Float *) ( (GIMM120) + (Umu) ) ) ; 
  W21 = *( (Float *) ( (GIMM021) + (Umu) ) ) ; 
  X21 = *( (Float *) ( (GIMM121) + (Umu) ) ) ; 
  C0 = W00 * U00  ;
  C1 = W00 * U01  ;
  D0 = W00 * U10  ;
  D1 = W00 * U11  ;
  E0 = X00 * U00  ;
  E1 = X00 * U01  ;
  C0 = W01 * U01 + C0 ;
  C1 = - ( W01 * U00 - C1 );
  D0 = W01 * U11 + D0 ;
  D1 = - ( W01 * U10 - D1 );
  E0 = X01 * U01 + E0 ;
  E1 = - ( X01 * U00 - E1 );
  C0 = W10 * U20 + C0 ;
  C1 = W10 * U21 + C1 ;
  D0 = W10 * V00 + D0 ;
  D1 = W10 * V01 + D1 ;
  E0 = X10 * U20 + E0 ;
  E1 = X10 * U21 + E1 ;
  C0 = W11 * U21 + C0 ;
  C1 = - ( W11 * U20 - C1 );
  D0 = W11 * V01 + D0 ;
  D1 = - ( W11 * V00 - D1 );
  E0 = X11 * U21 + E0 ;
  E1 = - ( X11 * U20 - E1 );
  C0 = W20 * V10 + C0 ;
  C1 = W20 * V11 + C1 ;
  D0 = W20 * V20 + D0 ;
  D1 = W20 * V21 + D1 ;
  E0 = X20 * V10 + E0 ;
  E1 = X20 * V11 + E1 ;
  C0 = W21 * V11 + C0 ;
  C1 = - ( W21 * V10 - C1 );
  D0 = W21 * V21 + D0 ;
  D1 = - ( W21 * V20 - D1 );
  E0 = X21 * V11 + E0 ;
  E1 = - ( X21 * V10 - E1 );
  *( (float *) ( (CHIIMM000) + (Chiout) ) ) = C0; 
  *( (float *) ( (CHIIMM001) + (Chiout) ) ) = C1; 
  *( (float *) ( (CHIIMM010) + (Chiout) ) ) = D0; 
  *( (float *) ( (CHIIMM011) + (Chiout) ) ) = D1; 
  *( (float *) ( (CHIIMM100) + (Chiout) ) ) = E0; 
  *( (float *) ( (CHIIMM101) + (Chiout) ) ) = E1; 
/*pragma_dcbt_space 4*/
  W00 = *( (Float *) ( (GIMM200) + (Umu) ) ) ; 
  W01 = *( (Float *) ( (GIMM201) + (Umu) ) ) ; 
  W10 = *( (Float *) ( (GIMM210) + (Umu) ) ) ; 
  W11 = *( (Float *) ( (GIMM211) + (Umu) ) ) ; 
  W20 = *( (Float *) ( (GIMM220) + (Umu) ) ) ; 
  W21 = *( (Float *) ( (GIMM221) + (Umu) ) ) ; 
  Umu = Umu + MAT_IMM ; 
  F0 = X00 * U10  ;
  F1 = X00 * U11  ;
  C0 = W00 * U00  ;
  C1 = W00 * U01  ;
  D0 = W00 * U10  ;
  D1 = W00 * U11  ;
  F0 = X01 * U11 + F0 ;
  F1 = - ( X01 * U10 - F1 );
  C0 = W01 * U01 + C0 ;
  C1 = - ( W01 * U00 - C1 );
  D0 = W01 * U11 + D0 ;
  D1 = - ( W01 * U10 - D1 );
  F0 = X10 * V00 + F0 ;
  F1 = X10 * V01 + F1 ;
  C0 = W10 * U20 + C0 ;
  C1 = W10 * U21 + C1 ;
  D0 = W10 * V00 + D0 ;
  D1 = W10 * V01 + D1 ;
  F0 = X11 * V01 + F0 ;
  F1 = - ( X11 * V00 - F1 );
  C0 = W11 * U21 + C0 ;
  C1 = - ( W11 * U20 - C1 );
  D0 = W11 * V01 + D0 ;
  D1 = - ( W11 * V00 - D1 );
  F0 = X20 * V20 + F0 ;
  F1 = X20 * V21 + F1 ;
  C0 = W20 * V10 + C0 ;
  C1 = W20 * V11 + C1 ;
  D0 = W20 * V20 + D0 ;
  D1 = W20 * V21 + D1 ;
  F0 = X21 * V21 + F0 ;
  F1 = - ( X21 * V20 - F1 );
  C0 = W21 * V11 + C0 ;
  C1 = - ( W21 * V10 - C1 );
  D0 = W21 * V21 + D0 ;
  D1 = - ( W21 * V20 - D1 );
  Chitmp1 = Chiout + ZERO_IMM ; 
  Chiin = Chiin + CHI_IMM ; 
  tab =   ( Isize ) + ( tab )  ;
  Chiout = * ( (unsigned long *) ( (zeromulti)+(tab) ) )  ;
/*pragma_store_lim 1*/
s_dec_hsu3_lab4:
/*pragma_load_lim 10*/
  *( (float *) ( (CHIIMM110) + (Chitmp1) ) ) = F0; 
  *( (float *) ( (CHIIMM111) + (Chitmp1) ) ) = F1; 
  *( (float *) ( (CHIIMM200) + (Chitmp1) ) ) = C0; 
  *( (float *) ( (CHIIMM201) + (Chitmp1) ) ) = C1; 
  *( (float *) ( (CHIIMM210) + (Chitmp1) ) ) = D0; 
  *( (float *) ( (CHIIMM211) + (Chitmp1) ) ) = D1; 
  U00 = *( (Float *) ( (UIMM00) + (psi) ) ) ; 
  U01 = *( (Float *) ( (UIMM01) + (psi) ) ) ; 
  V00 = *( (Float *) ( (VIMM00) + (psi) ) ) ; 
  V01 = *( (Float *) ( (VIMM01) + (psi) ) ) ; 
  W00 = *( (Float *) ( (WIMM00) + (psi) ) ) ; 
  W01 = *( (Float *) ( (WIMM01) + (psi) ) ) ; 
  X00 = *( (Float *) ( (XIMM00) + (psi) ) ) ; 
  X01 = *( (Float *) ( (XIMM01) + (psi) ) ) ; 
  U10 = *( (Float *) ( (UIMM10) + (psi) ) ) ; 
  U11 = *( (Float *) ( (UIMM11) + (psi) ) ) ; 
  V10 = *( (Float *) ( (VIMM10) + (psi) ) ) ; 
  V11 = *( (Float *) ( (VIMM11) + (psi) ) ) ; 
  W10 = *( (Float *) ( (WIMM10) + (psi) ) ) ; 
  W11 = *( (Float *) ( (WIMM11) + (psi) ) ) ; 
  X10 = *( (Float *) ( (XIMM10) + (psi) ) ) ; 
  X11 = *( (Float *) ( (XIMM11) + (psi) ) ) ; 
  length =   ( minus1 ) + ( length )  ;
  if ( length >  0 ) goto s_dec_hsu3_lab1;
s_dec_hsu3_lab0:

  return;
}
#ifdef __cplusplus
}
#endif
