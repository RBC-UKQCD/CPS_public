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
#define PSI_IMM000 0 
#define PSI_IMM001 8 
#define PSI_IMM010 16 
#define PSI_IMM011 24 
#define PSI_IMM020 32 
#define PSI_IMM021 40 
#define PSI_IMM100 48 
#define PSI_IMM101 56 
#define PSI_IMM110 64 
#define PSI_IMM111 72 
#define PSI_IMM120 80 
#define PSI_IMM121 88 
#define PSI_IMM200 96 
#define PSI_IMM201 104 
#define PSI_IMM210 112 
#define PSI_IMM211 120 
#define PSI_IMM220 128 
#define PSI_IMM221 136 
#define PSI_IMM300 144 
#define PSI_IMM301 152 
#define PSI_IMM310 160 
#define PSI_IMM311 168 
#define PSI_IMM320 176 
#define PSI_IMM321 184 
#define CHI_IMM000 0 
#define CHI_IMM001 4 
#define CHI_IMM010 8 
#define CHI_IMM011 12 
#define CHI_IMM100 16 
#define CHI_IMM101 20 
#define CHI_IMM110 24 
#define CHI_IMM111 28 
#define CHI_IMM200 32 
#define CHI_IMM201 36 
#define CHI_IMM210 40 
#define CHI_IMM211 44 
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
#define ZERO_IMM 0 
#define PSI_ATOM 192 
#define CHI_ATOM 48 
#define PAD_CHI_ATOM 64 
#define MAT_IMM 144 
#define Ndim 4 
#define Ndimm1 3 
#define hbitbucket 64 
#define hstk0 112 
#define hstk1 160 
#define hstk2 208 
#define hstk3 256 
#define Isize 4 
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
 void s_rec_su3 (unsigned long Arg0,unsigned long Arg1,unsigned long Arg2,unsigned long Arg3)
{

  /*Ensure natural alignment*/  int iStackArea[544+1];

  char *StackArea=(char *)iStackArea;

  unsigned long predicate;
  unsigned long prefscratch;
  unsigned long psi;
  unsigned long Umu;
  unsigned long Ufetch;
  unsigned long Chiin;
  unsigned long Chiout;
  unsigned long Chifetch;
  unsigned long Chiplus0;
  unsigned long Chiplus1;
  unsigned long Chiplus2;
  unsigned long Chiplus3;
  unsigned long Chiminus0;
  unsigned long Chiminus1;
  unsigned long Chiminus2;
  unsigned long Chiminus3;
  unsigned long mu;
  unsigned long Chidrain;
  unsigned long pref;
  unsigned long mem;
  unsigned long length;
  unsigned long CONSTR0;
  unsigned long CONSTR1;
  unsigned long CONSTR2;
  unsigned long CONSTR3;
  unsigned long CONSTR4;
  unsigned long R31;
  unsigned long StackPointer = (unsigned long) StackArea ;
  Float PSI000;
  Float PSI001;
  Float PSI010;
  Float PSI011;
  Float PSI020;
  Float PSI021;
  Float PSI030;
  Float PSI031;
  Float PSI100;
  Float PSI101;
  Float PSI110;
  Float PSI111;
  Float PSI120;
  Float PSI121;
  Float PSI130;
  Float PSI131;
  Float PSI200;
  Float PSI201;
  Float PSI210;
  Float PSI211;
  Float PSI220;
  Float PSI221;
  Float PSI230;
  Float PSI231;
  Float Atmp000;
  Float Atmp001;
  Float Atmp010;
  Float Atmp011;
  Float Btmp000;
  Float Btmp001;
  Float Btmp010;
  Float Btmp011;
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
  psi = Arg0 | Arg0 ; 
  Umu = Arg1 | Arg1 ; 
  Chiin = Arg2 | Arg2 ; 
  length = Arg3 | Arg3 ; 
  length = * ( (unsigned long *) ( (ZERO_IMM)+(length) ) )  ;
  if ( length <= 0 ) goto s_rec_su3_lab0;
CONSTR0=CONST0;
CONSTR1=CONST1;
CONSTR2=CONST2;
CONSTR3=CONST3;
CONSTR4=CONST4;
/*pragma_dcbt_space 5*/
/*pragma_dcbt_post 1*/
  PSI020 = *( (Float *) ( (GIMM000) + (Umu) ) ) ; 
  PSI110 = *( (Float *) ( (GIMM010) + (Umu) ) ) ; 
  PSI021 = *( (Float *) ( (GIMM001) + (Umu) ) ) ; 
  PSI111 = *( (Float *) ( (GIMM011) + (Umu) ) ) ; 
  PSI210 = *( (Float *) ( (GIMM100) + (Umu) ) ) ; 
  PSI200 = *( (Float *) ( (GIMM110) + (Umu) ) ) ; 
  PSI211 = *( (Float *) ( (GIMM101) + (Umu) ) ) ; 
  PSI201 = *( (Float *) ( (GIMM111) + (Umu) ) ) ; 
/*preload CONSTR0 Chiin */
/*preload CONSTR1 Chiin */
  Atmp000 = *( (float *) ( (CHI_IMM000) + (Chiin) ) ) ; 
  Atmp001 = *( (float *) ( (CHI_IMM001) + (Chiin) ) ) ; 
  Atmp010 = *( (float *) ( (CHI_IMM010) + (Chiin) ) ) ; 
  Atmp011 = *( (float *) ( (CHI_IMM011) + (Chiin) ) ) ; 
  Btmp000 = *( (float *) ( (CHI_IMM100) + (Chiin) ) ) ; 
  Btmp001 = *( (float *) ( (CHI_IMM101) + (Chiin) ) ) ; 
  Btmp010 = *( (float *) ( (CHI_IMM110) + (Chiin) ) ) ; 
  Btmp011 = *( (float *) ( (CHI_IMM111) + (Chiin) ) ) ; 
  Chidrain = mem + hbitbucket ; 
s_rec_su3_lab1:
  Chiout = mem + hstk0 ; 
mu=Ndimm1;
s_rec_su3_lab2:
/*pragma_dcbt_space 5*/
/*pragma_store_lim 1*/
/*pragma_load_lim 2*/
  *( (float *) ( (CHI_IMM110) + (Chidrain) ) ) = PSI030; 
  *( (float *) ( (CHI_IMM111) + (Chidrain) ) ) = PSI031; 
  *( (float *) ( (CHI_IMM200) + (Chidrain) ) ) = PSI130; 
  *( (float *) ( (CHI_IMM201) + (Chidrain) ) ) = PSI131; 
  *( (float *) ( (CHI_IMM210) + (Chidrain) ) ) = PSI230; 
  *( (float *) ( (CHI_IMM211) + (Chidrain) ) ) = PSI231; 
  PSI000 = *( (float *) ( (CHI_IMM200) + (Chiin) ) ) ; 
  PSI001 = *( (float *) ( (CHI_IMM201) + (Chiin) ) ) ; 
  PSI010 = *( (float *) ( (CHI_IMM210) + (Chiin) ) ) ; 
  PSI011 = *( (float *) ( (CHI_IMM211) + (Chiin) ) ) ; 
  PSI100 = *( (Float *) ( (GIMM200) + (Umu) ) ) ; 
  PSI120 = *( (Float *) ( (GIMM210) + (Umu) ) ) ; 
  PSI101 = *( (Float *) ( (GIMM201) + (Umu) ) ) ; 
  PSI121 = *( (Float *) ( (GIMM211) + (Umu) ) ) ; 
  Chifetch = Chiin + PAD_CHI_ATOM ; 
/*preload CONSTR0 Chifetch */
/*preload CONSTR1 Chifetch */
  Ufetch = Umu + MAT_IMM ; 
/*preload CONSTR0 Ufetch */
/*preload CONSTR1 Ufetch */
/*preload CONSTR2 Ufetch */
/*preload CONSTR3 Ufetch */
/*preload CONSTR4 Ufetch */
  PSI130 = PSI020 * Atmp000  ;
  PSI131 = PSI020 * Atmp001  ;
  PSI130 = - ( PSI021 * Atmp001 - PSI130 );
  PSI131 = PSI021 * Atmp000 + PSI131 ;
  PSI230 = PSI020 * Atmp010  ;
  PSI231 = PSI020 * Atmp011  ;
  PSI230 = - ( PSI021 * Atmp011 - PSI230 );
  PSI231 = PSI021 * Atmp010 + PSI231 ;
  PSI220 = PSI110 * Atmp000  ;
  PSI221 = PSI110 * Atmp001  ;
  PSI220 = - ( PSI111 * Atmp001 - PSI220 );
  PSI221 = PSI111 * Atmp000 + PSI221 ;
  PSI130 = PSI210 * Btmp000 + PSI130 ;
  PSI131 = PSI210 * Btmp001 + PSI131 ;
  PSI230 = PSI210 * Btmp010 + PSI230 ;
  PSI231 = PSI210 * Btmp011 + PSI231 ;
  PSI220 = PSI200 * Btmp000 + PSI220 ;
  PSI221 = PSI200 * Btmp001 + PSI221 ;
  PSI130 = - ( PSI211 * Btmp001 - PSI130 );
  PSI131 = PSI211 * Btmp000 + PSI131 ;
  PSI230 = - ( PSI211 * Btmp011 - PSI230 );
  PSI231 = PSI211 * Btmp010 + PSI231 ;
  PSI220 = - ( PSI201 * Btmp001 - PSI220 );
  PSI221 = PSI201 * Btmp000 + PSI221 ;
  PSI130 = PSI100 * PSI000 + PSI130 ;
  PSI131 = PSI100 * PSI001 + PSI131 ;
  PSI230 = PSI100 * PSI010 + PSI230 ;
  PSI231 = PSI100 * PSI011 + PSI231 ;
  PSI220 = PSI120 * PSI000 + PSI220 ;
  PSI221 = PSI120 * PSI001 + PSI221 ;
  PSI130 = - ( PSI101 * PSI001 - PSI130 );
  PSI131 = PSI101 * PSI000 + PSI131 ;
  PSI230 = - ( PSI101 * PSI011 - PSI230 );
  PSI231 = PSI101 * PSI010 + PSI231 ;
  PSI220 = - ( PSI121 * PSI001 - PSI220 );
  PSI221 = PSI121 * PSI000 + PSI221 ;
  *( (float *) ( (CHI_IMM000) + (Chiout) ) ) = PSI130; 
  *( (float *) ( (CHI_IMM001) + (Chiout) ) ) = PSI131; 
  *( (float *) ( (CHI_IMM010) + (Chiout) ) ) = PSI230; 
  *( (float *) ( (CHI_IMM011) + (Chiout) ) ) = PSI231; 
  *( (float *) ( (CHI_IMM100) + (Chiout) ) ) = PSI220; 
  *( (float *) ( (CHI_IMM101) + (Chiout) ) ) = PSI221; 
  PSI020 = *( (Float *) ( (GIMM020) + (Umu) ) ) ; 
  PSI021 = *( (Float *) ( (GIMM021) + (Umu) ) ) ; 
  PSI210 = *( (Float *) ( (GIMM120) + (Umu) ) ) ; 
  PSI211 = *( (Float *) ( (GIMM121) + (Umu) ) ) ; 
  PSI100 = *( (Float *) ( (GIMM220) + (Umu) ) ) ; 
  PSI101 = *( (Float *) ( (GIMM221) + (Umu) ) ) ; 
  Umu = Umu + MAT_IMM ; 
  PSI030 = PSI110 * Atmp010  ;
  PSI031 = PSI110 * Atmp011  ;
  PSI030 = - ( PSI111 * Atmp011 - PSI030 );
  PSI031 = PSI111 * Atmp010 + PSI031 ;
  PSI130 = PSI020 * Atmp000  ;
  PSI131 = PSI020 * Atmp001  ;
  PSI130 = - ( PSI021 * Atmp001 - PSI130 );
  PSI131 = PSI021 * Atmp000 + PSI131 ;
  PSI230 = PSI020 * Atmp010  ;
  PSI231 = PSI020 * Atmp011  ;
  PSI230 = - ( PSI021 * Atmp011 - PSI230 );
  PSI231 = PSI021 * Atmp010 + PSI231 ;
  PSI030 = PSI200 * Btmp010 + PSI030 ;
  PSI031 = PSI200 * Btmp011 + PSI031 ;
  PSI130 = PSI210 * Btmp000 + PSI130 ;
  PSI131 = PSI210 * Btmp001 + PSI131 ;
  PSI230 = PSI210 * Btmp010 + PSI230 ;
  PSI231 = PSI210 * Btmp011 + PSI231 ;
  PSI030 = - ( PSI201 * Btmp011 - PSI030 );
  PSI031 = PSI201 * Btmp010 + PSI031 ;
  PSI130 = - ( PSI211 * Btmp001 - PSI130 );
  PSI131 = PSI211 * Btmp000 + PSI131 ;
  PSI230 = - ( PSI211 * Btmp011 - PSI230 );
  PSI231 = PSI211 * Btmp010 + PSI231 ;
  PSI030 = PSI120 * PSI010 + PSI030 ;
  PSI031 = PSI120 * PSI011 + PSI031 ;
  PSI130 = PSI100 * PSI000 + PSI130 ;
  PSI131 = PSI100 * PSI001 + PSI131 ;
  PSI230 = PSI100 * PSI010 + PSI230 ;
  PSI231 = PSI100 * PSI011 + PSI231 ;
  PSI030 = - ( PSI121 * PSI011 - PSI030 );
  PSI031 = PSI121 * PSI010 + PSI031 ;
  PSI130 = - ( PSI101 * PSI001 - PSI130 );
  PSI131 = PSI101 * PSI000 + PSI131 ;
  PSI230 = - ( PSI101 * PSI011 - PSI230 );
  PSI231 = PSI101 * PSI010 + PSI231 ;
  Chiin = Chiin + PAD_CHI_ATOM ; 
  Chidrain = Chiout + ZERO_IMM ; 
  Chiout = Chiout + CHI_ATOM ; 
  PSI020 = *( (Float *) ( (GIMM000) + (Umu) ) ) ; 
  PSI110 = *( (Float *) ( (GIMM010) + (Umu) ) ) ; 
  PSI021 = *( (Float *) ( (GIMM001) + (Umu) ) ) ; 
  PSI111 = *( (Float *) ( (GIMM011) + (Umu) ) ) ; 
  PSI210 = *( (Float *) ( (GIMM100) + (Umu) ) ) ; 
  PSI200 = *( (Float *) ( (GIMM110) + (Umu) ) ) ; 
  PSI211 = *( (Float *) ( (GIMM101) + (Umu) ) ) ; 
  PSI201 = *( (Float *) ( (GIMM111) + (Umu) ) ) ; 
  Atmp000 = *( (float *) ( (CHI_IMM000) + (Chiin) ) ) ; 
  Atmp001 = *( (float *) ( (CHI_IMM001) + (Chiin) ) ) ; 
  Atmp010 = *( (float *) ( (CHI_IMM010) + (Chiin) ) ) ; 
  Atmp011 = *( (float *) ( (CHI_IMM011) + (Chiin) ) ) ; 
  Btmp000 = *( (float *) ( (CHI_IMM100) + (Chiin) ) ) ; 
  Btmp001 = *( (float *) ( (CHI_IMM101) + (Chiin) ) ) ; 
  Btmp010 = *( (float *) ( (CHI_IMM110) + (Chiin) ) ) ; 
  Btmp011 = *( (float *) ( (CHI_IMM111) + (Chiin) ) ) ; 
  mu =   ( minus1 ) + ( mu )  ;
  if ( mu >  0 ) goto s_rec_su3_lab2;
s_rec_su3_lab3:
/*pragma_store_lim 2*/
/*pragma_dcbt_space 5*/
/*pragma_dcbt_post 1*/
/*pragma_dcbt_pre 0*/
/*pragma_load_lim 2*/
  *( (float *) ( (CHI_IMM110) + (Chidrain) ) ) = PSI030; 
  *( (float *) ( (CHI_IMM111) + (Chidrain) ) ) = PSI031; 
  *( (float *) ( (CHI_IMM200) + (Chidrain) ) ) = PSI130; 
  *( (float *) ( (CHI_IMM201) + (Chidrain) ) ) = PSI131; 
  *( (float *) ( (CHI_IMM210) + (Chidrain) ) ) = PSI230; 
  *( (float *) ( (CHI_IMM211) + (Chidrain) ) ) = PSI231; 
  PSI000 = *( (float *) ( (CHI_IMM200) + (Chiin) ) ) ; 
  PSI001 = *( (float *) ( (CHI_IMM201) + (Chiin) ) ) ; 
  PSI010 = *( (float *) ( (CHI_IMM210) + (Chiin) ) ) ; 
  PSI011 = *( (float *) ( (CHI_IMM211) + (Chiin) ) ) ; 
  PSI100 = *( (Float *) ( (GIMM200) + (Umu) ) ) ; 
  PSI120 = *( (Float *) ( (GIMM210) + (Umu) ) ) ; 
  PSI101 = *( (Float *) ( (GIMM201) + (Umu) ) ) ; 
  PSI121 = *( (Float *) ( (GIMM211) + (Umu) ) ) ; 
/*pragma_dcbt_space 3*/
  Chifetch = Chiin + PAD_CHI_ATOM ; 
/*preload CONSTR0 Chifetch */
/*preload CONSTR1 Chifetch */
  Chifetch = Chifetch + PAD_CHI_ATOM ; 
/*preload CONSTR0 Chifetch */
/*preload CONSTR1 Chifetch */
  Chifetch = Chifetch + PAD_CHI_ATOM ; 
/*preload CONSTR0 Chifetch */
/*preload CONSTR1 Chifetch */
  Chifetch = Chifetch + PAD_CHI_ATOM ; 
/*preload CONSTR0 Chifetch */
/*preload CONSTR1 Chifetch */
  PSI130 = PSI020 * Atmp000  ;
  PSI131 = PSI020 * Atmp001  ;
  PSI130 = - ( PSI021 * Atmp001 - PSI130 );
  PSI131 = PSI021 * Atmp000 + PSI131 ;
  PSI230 = PSI020 * Atmp010  ;
  PSI231 = PSI020 * Atmp011  ;
  PSI230 = - ( PSI021 * Atmp011 - PSI230 );
  PSI231 = PSI021 * Atmp010 + PSI231 ;
  PSI220 = PSI110 * Atmp000  ;
  PSI221 = PSI110 * Atmp001  ;
  PSI220 = - ( PSI111 * Atmp001 - PSI220 );
  PSI221 = PSI111 * Atmp000 + PSI221 ;
  PSI130 = PSI210 * Btmp000 + PSI130 ;
  PSI131 = PSI210 * Btmp001 + PSI131 ;
  PSI230 = PSI210 * Btmp010 + PSI230 ;
  PSI231 = PSI210 * Btmp011 + PSI231 ;
  PSI220 = PSI200 * Btmp000 + PSI220 ;
  PSI221 = PSI200 * Btmp001 + PSI221 ;
  PSI130 = - ( PSI211 * Btmp001 - PSI130 );
  PSI131 = PSI211 * Btmp000 + PSI131 ;
  PSI230 = - ( PSI211 * Btmp011 - PSI230 );
  PSI231 = PSI211 * Btmp010 + PSI231 ;
  PSI220 = - ( PSI201 * Btmp001 - PSI220 );
  PSI221 = PSI201 * Btmp000 + PSI221 ;
  PSI130 = PSI100 * PSI000 + PSI130 ;
  PSI131 = PSI100 * PSI001 + PSI131 ;
  PSI230 = PSI100 * PSI010 + PSI230 ;
  PSI231 = PSI100 * PSI011 + PSI231 ;
  PSI220 = PSI120 * PSI000 + PSI220 ;
  PSI221 = PSI120 * PSI001 + PSI221 ;
  PSI130 = - ( PSI101 * PSI001 - PSI130 );
  PSI131 = PSI101 * PSI000 + PSI131 ;
  PSI230 = - ( PSI101 * PSI011 - PSI230 );
  PSI231 = PSI101 * PSI010 + PSI231 ;
  PSI220 = - ( PSI121 * PSI001 - PSI220 );
  PSI221 = PSI121 * PSI000 + PSI221 ;
  *( (float *) ( (CHI_IMM000) + (Chiout) ) ) = PSI130; 
  *( (float *) ( (CHI_IMM001) + (Chiout) ) ) = PSI131; 
  *( (float *) ( (CHI_IMM010) + (Chiout) ) ) = PSI230; 
  *( (float *) ( (CHI_IMM011) + (Chiout) ) ) = PSI231; 
  *( (float *) ( (CHI_IMM100) + (Chiout) ) ) = PSI220; 
  *( (float *) ( (CHI_IMM101) + (Chiout) ) ) = PSI221; 
  PSI020 = *( (Float *) ( (GIMM020) + (Umu) ) ) ; 
  PSI021 = *( (Float *) ( (GIMM021) + (Umu) ) ) ; 
  PSI210 = *( (Float *) ( (GIMM120) + (Umu) ) ) ; 
  PSI211 = *( (Float *) ( (GIMM121) + (Umu) ) ) ; 
  PSI100 = *( (Float *) ( (GIMM220) + (Umu) ) ) ; 
  PSI101 = *( (Float *) ( (GIMM221) + (Umu) ) ) ; 
  Umu = Umu + MAT_IMM ; 
  PSI030 = PSI110 * Atmp010  ;
  PSI031 = PSI110 * Atmp011  ;
  PSI030 = - ( PSI111 * Atmp011 - PSI030 );
  PSI031 = PSI111 * Atmp010 + PSI031 ;
  PSI130 = PSI020 * Atmp000  ;
  PSI131 = PSI020 * Atmp001  ;
  PSI130 = - ( PSI021 * Atmp001 - PSI130 );
  PSI131 = PSI021 * Atmp000 + PSI131 ;
  PSI230 = PSI020 * Atmp010  ;
  PSI231 = PSI020 * Atmp011  ;
  PSI230 = - ( PSI021 * Atmp011 - PSI230 );
  PSI231 = PSI021 * Atmp010 + PSI231 ;
  PSI030 = PSI200 * Btmp010 + PSI030 ;
  PSI031 = PSI200 * Btmp011 + PSI031 ;
  PSI130 = PSI210 * Btmp000 + PSI130 ;
  PSI131 = PSI210 * Btmp001 + PSI131 ;
  PSI230 = PSI210 * Btmp010 + PSI230 ;
  PSI231 = PSI210 * Btmp011 + PSI231 ;
  PSI030 = - ( PSI201 * Btmp011 - PSI030 );
  PSI031 = PSI201 * Btmp010 + PSI031 ;
  PSI130 = - ( PSI211 * Btmp001 - PSI130 );
  PSI131 = PSI211 * Btmp000 + PSI131 ;
  PSI230 = - ( PSI211 * Btmp011 - PSI230 );
  PSI231 = PSI211 * Btmp010 + PSI231 ;
  PSI030 = PSI120 * PSI010 + PSI030 ;
  PSI031 = PSI120 * PSI011 + PSI031 ;
  PSI130 = PSI100 * PSI000 + PSI130 ;
  PSI131 = PSI100 * PSI001 + PSI131 ;
  PSI230 = PSI100 * PSI010 + PSI230 ;
  PSI231 = PSI100 * PSI011 + PSI231 ;
  PSI030 = - ( PSI121 * PSI011 - PSI030 );
  PSI031 = PSI121 * PSI010 + PSI031 ;
  PSI130 = - ( PSI101 * PSI001 - PSI130 );
  PSI131 = PSI101 * PSI000 + PSI131 ;
  PSI230 = - ( PSI101 * PSI011 - PSI230 );
  PSI231 = PSI101 * PSI010 + PSI231 ;
  Chiin = Chiin + PAD_CHI_ATOM ; 
  Chidrain = Chiout + ZERO_IMM ; 
  Chiout = Chiout + CHI_ATOM ; 
/*pragma_store_inorder 1*/
  Chiminus0 = mem + hstk0 ; 
/*pragma_load_lim 2*/
  Atmp000 = *( (float *) ( (CHI_IMM000) + (Chiminus0) ) ) ; 
  Btmp000 = *( (float *) ( (CHI_IMM000) + (Chiin) ) ) ; 
  Chiplus0 = Chiin + ZERO_IMM ; 
  Atmp001 = *( (float *) ( (CHI_IMM001) + (Chiminus0) ) ) ; 
  Btmp001 = *( (float *) ( (CHI_IMM001) + (Chiin) ) ) ; 
  Chiplus0 = Chiin + ZERO_IMM ; 
  Atmp010 = *( (float *) ( (CHI_IMM010) + (Chiminus0) ) ) ; 
  Btmp010 = *( (float *) ( (CHI_IMM010) + (Chiin) ) ) ; 
  Chiplus0 = Chiin + ZERO_IMM ; 
  Atmp011 = *( (float *) ( (CHI_IMM011) + (Chiminus0) ) ) ; 
  Btmp011 = *( (float *) ( (CHI_IMM011) + (Chiin) ) ) ; 
  Chiplus0 = Chiin + ZERO_IMM ; 
/*pragma_dcbt_post 1*/
  *( (float *) ( (CHI_IMM110) + (Chidrain) ) ) = PSI030; 
  *( (float *) ( (CHI_IMM111) + (Chidrain) ) ) = PSI031; 
  *( (float *) ( (CHI_IMM200) + (Chidrain) ) ) = PSI130; 
  *( (float *) ( (CHI_IMM201) + (Chidrain) ) ) = PSI131; 
  *( (float *) ( (CHI_IMM210) + (Chidrain) ) ) = PSI230; 
  *( (float *) ( (CHI_IMM211) + (Chidrain) ) ) = PSI231; 
/*pragma_load_lim 1*/
  Chiminus1 = mem + hstk1 ; 
  Chiplus1 = Chiin + PAD_CHI_ATOM ; 
  PSI000 = Btmp000 + Atmp000 ; 
  PSI001 = Btmp001 + Atmp001 ; 
  PSI010 = Btmp010 + Atmp010 ; 
  PSI011 = Btmp011 + Atmp011 ; 
  PSI020 = Btmp011 - Atmp011 ; 
  PSI021 = Atmp010 - Btmp010 ; 
  PSI030 = Btmp001 - Atmp001 ; 
  PSI031 = Atmp000 - Btmp000 ; 
  PSI200 = *( (float *) ( (CHI_IMM000) + (Chiminus1) ) ) ; 
  PSI220 = *( (float *) ( (CHI_IMM000) + (Chiplus1) ) ) ; 
  PSI201 = *( (float *) ( (CHI_IMM001) + (Chiminus1) ) ) ; 
  PSI221 = *( (float *) ( (CHI_IMM001) + (Chiplus1) ) ) ; 
  PSI210 = *( (float *) ( (CHI_IMM010) + (Chiminus1) ) ) ; 
  PSI230 = *( (float *) ( (CHI_IMM010) + (Chiplus1) ) ) ; 
  PSI211 = *( (float *) ( (CHI_IMM011) + (Chiminus1) ) ) ; 
  PSI231 = *( (float *) ( (CHI_IMM011) + (Chiplus1) ) ) ; 
  Chiminus2 = mem + hstk2 ; 
  Chiminus3 = mem + hstk3 ; 
  Chiplus2 = Chiplus1 + PAD_CHI_ATOM ; 
  Chiplus3 = Chiplus2 + PAD_CHI_ATOM ; 
  PSI000 = PSI000 + PSI220 ; 
  PSI001 = PSI001 + PSI221 ; 
  PSI010 = PSI010 + PSI230 ; 
  PSI011 = PSI011 + PSI231 ; 
  PSI000 = PSI000 + PSI200 ; 
  PSI001 = PSI001 + PSI201 ; 
  PSI010 = PSI010 + PSI210 ; 
  PSI011 = PSI011 + PSI211 ; 
  PSI020 = PSI020 + PSI230 ; 
  PSI021 = PSI021 + PSI231 ; 
  PSI020 = PSI020 - PSI210 ; 
  PSI021 = PSI021 - PSI211 ; 
  PSI030 = PSI030 - PSI220 ; 
  PSI031 = PSI031 - PSI221 ; 
  PSI030 = PSI030 + PSI200 ; 
  PSI031 = PSI031 + PSI201 ; 
  Atmp000 = *( (float *) ( (CHI_IMM000) + (Chiminus2) ) ) ; 
  Btmp000 = *( (float *) ( (CHI_IMM000) + (Chiplus2) ) ) ; 
  Atmp001 = *( (float *) ( (CHI_IMM001) + (Chiminus2) ) ) ; 
  Btmp001 = *( (float *) ( (CHI_IMM001) + (Chiplus2) ) ) ; 
  Atmp010 = *( (float *) ( (CHI_IMM010) + (Chiminus2) ) ) ; 
  Btmp010 = *( (float *) ( (CHI_IMM010) + (Chiplus2) ) ) ; 
  Atmp011 = *( (float *) ( (CHI_IMM011) + (Chiminus2) ) ) ; 
  Btmp011 = *( (float *) ( (CHI_IMM011) + (Chiplus2) ) ) ; 
  PSI000 = PSI000 + Btmp000 ; 
  PSI001 = PSI001 + Btmp001 ; 
  PSI010 = PSI010 + Btmp010 ; 
  PSI011 = PSI011 + Btmp011 ; 
  PSI000 = PSI000 + Atmp000 ; 
  PSI001 = PSI001 + Atmp001 ; 
  PSI010 = PSI010 + Atmp010 ; 
  PSI011 = PSI011 + Atmp011 ; 
  PSI020 = PSI020 + Btmp001 ; 
  PSI021 = PSI021 - Btmp000 ; 
  PSI020 = PSI020 - Atmp001 ; 
  PSI021 = PSI021 + Atmp000 ; 
  PSI030 = PSI030 - Btmp011 ; 
  PSI031 = PSI031 + Btmp010 ; 
  PSI030 = PSI030 + Atmp011 ; 
  PSI031 = PSI031 - Atmp010 ; 
/*pragma_load_lim 2*/
  PSI200 = *( (float *) ( (CHI_IMM000) + (Chiminus3) ) ) ; 
  PSI220 = *( (float *) ( (CHI_IMM000) + (Chiplus3) ) ) ; 
  PSI201 = *( (float *) ( (CHI_IMM001) + (Chiminus3) ) ) ; 
  PSI221 = *( (float *) ( (CHI_IMM001) + (Chiplus3) ) ) ; 
  PSI210 = *( (float *) ( (CHI_IMM010) + (Chiminus3) ) ) ; 
  PSI230 = *( (float *) ( (CHI_IMM010) + (Chiplus3) ) ) ; 
  PSI211 = *( (float *) ( (CHI_IMM011) + (Chiminus3) ) ) ; 
  PSI231 = *( (float *) ( (CHI_IMM011) + (Chiplus3) ) ) ; 
  PSI000 = PSI000 + PSI220 ; 
  PSI001 = PSI001 + PSI221 ; 
  PSI010 = PSI010 + PSI230 ; 
  PSI011 = PSI011 + PSI231 ; 
  PSI020 = PSI020 + PSI220 ; 
  PSI021 = PSI021 + PSI221 ; 
  PSI030 = PSI030 + PSI230 ; 
  PSI031 = PSI031 + PSI231 ; 
  PSI000 = PSI000 + PSI200 ; 
  PSI001 = PSI001 + PSI201 ; 
  PSI010 = PSI010 + PSI210 ; 
  PSI011 = PSI011 + PSI211 ; 
  PSI020 = PSI020 - PSI200 ; 
  PSI021 = PSI021 - PSI201 ; 
  PSI030 = PSI030 - PSI210 ; 
  PSI031 = PSI031 - PSI211 ; 
  Atmp000 = *( (float *) ( (CHI_IMM100) + (Chiminus0) ) ) ; 
  Btmp000 = *( (float *) ( (CHI_IMM100) + (Chiplus0) ) ) ; 
  Atmp001 = *( (float *) ( (CHI_IMM101) + (Chiminus0) ) ) ; 
  Btmp001 = *( (float *) ( (CHI_IMM101) + (Chiplus0) ) ) ; 
  Atmp010 = *( (float *) ( (CHI_IMM110) + (Chiminus0) ) ) ; 
  Btmp010 = *( (float *) ( (CHI_IMM110) + (Chiplus0) ) ) ; 
  Atmp011 = *( (float *) ( (CHI_IMM111) + (Chiminus0) ) ) ; 
  Btmp011 = *( (float *) ( (CHI_IMM111) + (Chiplus0) ) ) ; 
  *( (Float *) ( (PSI_IMM000) + (psi) ) ) = PSI000; 
  *( (Float *) ( (PSI_IMM001) + (psi) ) ) = PSI001; 
/*pragma_load_lim 1*/
  PSI100 = Btmp000 + Atmp000 ; 
  PSI101 = Btmp001 + Atmp001 ; 
  PSI110 = Btmp010 + Atmp010 ; 
  PSI111 = Btmp011 + Atmp011 ; 
  PSI120 = Btmp011 - Atmp011 ; 
  PSI121 = Atmp010 - Btmp010 ; 
  PSI130 = Btmp001 - Atmp001 ; 
  PSI131 = Atmp000 - Btmp000 ; 
  PSI200 = *( (float *) ( (CHI_IMM100) + (Chiminus1) ) ) ; 
  PSI220 = *( (float *) ( (CHI_IMM100) + (Chiplus1) ) ) ; 
  PSI201 = *( (float *) ( (CHI_IMM101) + (Chiminus1) ) ) ; 
  PSI221 = *( (float *) ( (CHI_IMM101) + (Chiplus1) ) ) ; 
  PSI210 = *( (float *) ( (CHI_IMM110) + (Chiminus1) ) ) ; 
  PSI230 = *( (float *) ( (CHI_IMM110) + (Chiplus1) ) ) ; 
  PSI211 = *( (float *) ( (CHI_IMM111) + (Chiminus1) ) ) ; 
  PSI231 = *( (float *) ( (CHI_IMM111) + (Chiplus1) ) ) ; 
  PSI100 = PSI100 + PSI220 ; 
  PSI101 = PSI101 + PSI221 ; 
  PSI110 = PSI110 + PSI230 ; 
  PSI111 = PSI111 + PSI231 ; 
  PSI100 = PSI100 + PSI200 ; 
  PSI101 = PSI101 + PSI201 ; 
  PSI110 = PSI110 + PSI210 ; 
  PSI111 = PSI111 + PSI211 ; 
  PSI120 = PSI120 + PSI230 ; 
  PSI121 = PSI121 + PSI231 ; 
  PSI120 = PSI120 - PSI210 ; 
  PSI121 = PSI121 - PSI211 ; 
  PSI130 = PSI130 - PSI220 ; 
  PSI131 = PSI131 - PSI221 ; 
  PSI130 = PSI130 + PSI200 ; 
  PSI131 = PSI131 + PSI201 ; 
  Atmp000 = *( (float *) ( (CHI_IMM100) + (Chiminus2) ) ) ; 
  Btmp000 = *( (float *) ( (CHI_IMM100) + (Chiplus2) ) ) ; 
  Atmp001 = *( (float *) ( (CHI_IMM101) + (Chiminus2) ) ) ; 
  Btmp001 = *( (float *) ( (CHI_IMM101) + (Chiplus2) ) ) ; 
  Atmp010 = *( (float *) ( (CHI_IMM110) + (Chiminus2) ) ) ; 
  Btmp010 = *( (float *) ( (CHI_IMM110) + (Chiplus2) ) ) ; 
  Atmp011 = *( (float *) ( (CHI_IMM111) + (Chiminus2) ) ) ; 
  Btmp011 = *( (float *) ( (CHI_IMM111) + (Chiplus2) ) ) ; 
  PSI100 = PSI100 + Btmp000 ; 
  PSI101 = PSI101 + Btmp001 ; 
  PSI110 = PSI110 + Btmp010 ; 
  PSI111 = PSI111 + Btmp011 ; 
  PSI100 = PSI100 + Atmp000 ; 
  PSI101 = PSI101 + Atmp001 ; 
  PSI110 = PSI110 + Atmp010 ; 
  PSI111 = PSI111 + Atmp011 ; 
  PSI120 = PSI120 + Btmp001 ; 
  PSI121 = PSI121 - Btmp000 ; 
  PSI120 = PSI120 - Atmp001 ; 
  PSI121 = PSI121 + Atmp000 ; 
  PSI130 = PSI130 - Btmp011 ; 
  PSI131 = PSI131 + Btmp010 ; 
  PSI130 = PSI130 + Atmp011 ; 
  PSI131 = PSI131 - Atmp010 ; 
/*pragma_load_lim 2*/
  PSI200 = *( (float *) ( (CHI_IMM100) + (Chiminus3) ) ) ; 
  PSI220 = *( (float *) ( (CHI_IMM100) + (Chiplus3) ) ) ; 
  PSI201 = *( (float *) ( (CHI_IMM101) + (Chiminus3) ) ) ; 
  PSI221 = *( (float *) ( (CHI_IMM101) + (Chiplus3) ) ) ; 
  PSI210 = *( (float *) ( (CHI_IMM110) + (Chiminus3) ) ) ; 
  PSI230 = *( (float *) ( (CHI_IMM110) + (Chiplus3) ) ) ; 
  PSI211 = *( (float *) ( (CHI_IMM111) + (Chiminus3) ) ) ; 
  PSI231 = *( (float *) ( (CHI_IMM111) + (Chiplus3) ) ) ; 
  PSI100 = PSI100 + PSI220 ; 
  PSI101 = PSI101 + PSI221 ; 
  PSI110 = PSI110 + PSI230 ; 
  PSI111 = PSI111 + PSI231 ; 
  PSI120 = PSI120 + PSI220 ; 
  PSI121 = PSI121 + PSI221 ; 
  PSI130 = PSI130 + PSI230 ; 
  PSI131 = PSI131 + PSI231 ; 
  PSI100 = PSI100 + PSI200 ; 
  PSI101 = PSI101 + PSI201 ; 
  PSI110 = PSI110 + PSI210 ; 
  PSI111 = PSI111 + PSI211 ; 
  PSI120 = PSI120 - PSI200 ; 
  PSI121 = PSI121 - PSI201 ; 
  PSI130 = PSI130 - PSI210 ; 
  PSI131 = PSI131 - PSI211 ; 
  Atmp000 = *( (float *) ( (CHI_IMM200) + (Chiminus0) ) ) ; 
  Btmp000 = *( (float *) ( (CHI_IMM200) + (Chiplus0) ) ) ; 
  Atmp001 = *( (float *) ( (CHI_IMM201) + (Chiminus0) ) ) ; 
  Btmp001 = *( (float *) ( (CHI_IMM201) + (Chiplus0) ) ) ; 
  Atmp010 = *( (float *) ( (CHI_IMM210) + (Chiminus0) ) ) ; 
  Btmp010 = *( (float *) ( (CHI_IMM210) + (Chiplus0) ) ) ; 
  Atmp011 = *( (float *) ( (CHI_IMM211) + (Chiminus0) ) ) ; 
  Btmp011 = *( (float *) ( (CHI_IMM211) + (Chiplus0) ) ) ; 
  *( (Float *) ( (PSI_IMM010) + (psi) ) ) = PSI100; 
  *( (Float *) ( (PSI_IMM011) + (psi) ) ) = PSI101; 
/*pragma_load_lim 1*/
/*pragma_dcbt_post 0*/
/*pragma_dcbt_space 1*/
  Ufetch = Umu + ZERO_IMM ; 
/*preload CONSTR1 Ufetch */
/*preload CONSTR2 Ufetch */
/*preload CONSTR3 Ufetch */
/*preload CONSTR4 Ufetch */
  PSI200 = Btmp000 + Atmp000 ; 
  PSI201 = Btmp001 + Atmp001 ; 
  PSI210 = Btmp010 + Atmp010 ; 
  PSI211 = Btmp011 + Atmp011 ; 
  PSI220 = Btmp011 - Atmp011 ; 
  PSI221 = Atmp010 - Btmp010 ; 
  PSI230 = Btmp001 - Atmp001 ; 
  PSI231 = Atmp000 - Btmp000 ; 
  PSI000 = *( (float *) ( (CHI_IMM200) + (Chiminus1) ) ) ; 
  Btmp000 = *( (float *) ( (CHI_IMM200) + (Chiplus1) ) ) ; 
  PSI001 = *( (float *) ( (CHI_IMM201) + (Chiminus1) ) ) ; 
  Btmp001 = *( (float *) ( (CHI_IMM201) + (Chiplus1) ) ) ; 
  PSI100 = *( (float *) ( (CHI_IMM210) + (Chiminus1) ) ) ; 
  Btmp010 = *( (float *) ( (CHI_IMM210) + (Chiplus1) ) ) ; 
  PSI101 = *( (float *) ( (CHI_IMM211) + (Chiminus1) ) ) ; 
  Btmp011 = *( (float *) ( (CHI_IMM211) + (Chiplus1) ) ) ; 
  PSI200 = PSI200 + Btmp000 ; 
  PSI201 = PSI201 + Btmp001 ; 
  PSI210 = PSI210 + Btmp010 ; 
  PSI211 = PSI211 + Btmp011 ; 
  PSI200 = PSI200 + PSI000 ; 
  PSI201 = PSI201 + PSI001 ; 
  PSI210 = PSI210 + PSI100 ; 
  PSI211 = PSI211 + PSI101 ; 
  PSI220 = PSI220 + Btmp010 ; 
  PSI221 = PSI221 + Btmp011 ; 
  PSI220 = PSI220 - PSI100 ; 
  PSI221 = PSI221 - PSI101 ; 
  PSI230 = PSI230 - Btmp000 ; 
  PSI231 = PSI231 - Btmp001 ; 
  PSI230 = PSI230 + PSI000 ; 
  PSI231 = PSI231 + PSI001 ; 
  Atmp000 = *( (float *) ( (CHI_IMM200) + (Chiminus2) ) ) ; 
  Btmp000 = *( (float *) ( (CHI_IMM200) + (Chiplus2) ) ) ; 
  Atmp001 = *( (float *) ( (CHI_IMM201) + (Chiminus2) ) ) ; 
  Btmp001 = *( (float *) ( (CHI_IMM201) + (Chiplus2) ) ) ; 
  Atmp010 = *( (float *) ( (CHI_IMM210) + (Chiminus2) ) ) ; 
  Btmp010 = *( (float *) ( (CHI_IMM210) + (Chiplus2) ) ) ; 
  Atmp011 = *( (float *) ( (CHI_IMM211) + (Chiminus2) ) ) ; 
  Btmp011 = *( (float *) ( (CHI_IMM211) + (Chiplus2) ) ) ; 
  PSI200 = PSI200 + Btmp000 ; 
  PSI201 = PSI201 + Btmp001 ; 
  PSI210 = PSI210 + Btmp010 ; 
  PSI211 = PSI211 + Btmp011 ; 
  PSI200 = PSI200 + Atmp000 ; 
  PSI201 = PSI201 + Atmp001 ; 
  PSI210 = PSI210 + Atmp010 ; 
  PSI211 = PSI211 + Atmp011 ; 
  PSI220 = PSI220 + Btmp001 ; 
  PSI221 = PSI221 - Btmp000 ; 
  PSI220 = PSI220 - Atmp001 ; 
  PSI221 = PSI221 + Atmp000 ; 
  PSI230 = PSI230 - Btmp011 ; 
  PSI231 = PSI231 + Btmp010 ; 
  PSI230 = PSI230 + Atmp011 ; 
  PSI231 = PSI231 - Atmp010 ; 
/*pragma_load_lim 2*/
  PSI000 = *( (float *) ( (CHI_IMM200) + (Chiminus3) ) ) ; 
  Btmp000 = *( (float *) ( (CHI_IMM200) + (Chiplus3) ) ) ; 
  PSI001 = *( (float *) ( (CHI_IMM201) + (Chiminus3) ) ) ; 
  Btmp001 = *( (float *) ( (CHI_IMM201) + (Chiplus3) ) ) ; 
  PSI100 = *( (float *) ( (CHI_IMM210) + (Chiminus3) ) ) ; 
  Btmp010 = *( (float *) ( (CHI_IMM210) + (Chiplus3) ) ) ; 
  PSI101 = *( (float *) ( (CHI_IMM211) + (Chiminus3) ) ) ; 
  Btmp011 = *( (float *) ( (CHI_IMM211) + (Chiplus3) ) ) ; 
  PSI200 = PSI200 + Btmp000 ; 
  PSI201 = PSI201 + Btmp001 ; 
  PSI210 = PSI210 + Btmp010 ; 
  PSI211 = PSI211 + Btmp011 ; 
  PSI220 = PSI220 + Btmp000 ; 
  PSI221 = PSI221 + Btmp001 ; 
  PSI230 = PSI230 + Btmp010 ; 
  PSI231 = PSI231 + Btmp011 ; 
  PSI200 = PSI200 + PSI000 ; 
  PSI201 = PSI201 + PSI001 ; 
  PSI210 = PSI210 + PSI100 ; 
  PSI211 = PSI211 + PSI101 ; 
  PSI220 = PSI220 - PSI000 ; 
  PSI221 = PSI221 - PSI001 ; 
  PSI230 = PSI230 - PSI100 ; 
  PSI231 = PSI231 - PSI101 ; 
  *( (Float *) ( (PSI_IMM020) + (psi) ) ) = PSI200; 
  *( (Float *) ( (PSI_IMM021) + (psi) ) ) = PSI201; 
/*pragma_store_lim 2*/
/*pragma_dcbt_space 8*/
  *( (Float *) ( (PSI_IMM100) + (psi) ) ) = PSI010; 
  *( (Float *) ( (PSI_IMM101) + (psi) ) ) = PSI011; 
  *( (Float *) ( (PSI_IMM110) + (psi) ) ) = PSI110; 
  *( (Float *) ( (PSI_IMM111) + (psi) ) ) = PSI111; 
  *( (Float *) ( (PSI_IMM120) + (psi) ) ) = PSI210; 
  *( (Float *) ( (PSI_IMM121) + (psi) ) ) = PSI211; 
  *( (Float *) ( (PSI_IMM200) + (psi) ) ) = PSI020; 
  *( (Float *) ( (PSI_IMM201) + (psi) ) ) = PSI021; 
  *( (Float *) ( (PSI_IMM210) + (psi) ) ) = PSI120; 
  *( (Float *) ( (PSI_IMM211) + (psi) ) ) = PSI121; 
  *( (Float *) ( (PSI_IMM220) + (psi) ) ) = PSI220; 
  *( (Float *) ( (PSI_IMM221) + (psi) ) ) = PSI221; 
  Chidrain = mem + hbitbucket ; 
  *( (Float *) ( (PSI_IMM300) + (psi) ) ) = PSI030; 
  *( (Float *) ( (PSI_IMM301) + (psi) ) ) = PSI031; 
  *( (Float *) ( (PSI_IMM310) + (psi) ) ) = PSI130; 
  *( (Float *) ( (PSI_IMM311) + (psi) ) ) = PSI131; 
  *( (Float *) ( (PSI_IMM320) + (psi) ) ) = PSI230; 
  *( (Float *) ( (PSI_IMM321) + (psi) ) ) = PSI231; 
  psi = psi + PSI_ATOM ; 
  Chiplus3 = Chiplus3 + ZERO_IMM ; 
  Chiin = Chiplus3 + PAD_CHI_ATOM ; 
/*pragma_dcbt_space 0*/
/*preload CONSTR0 Chiin */
/*preload CONSTR1 Chiin */
  PSI020 = *( (Float *) ( (GIMM000) + (Umu) ) ) ; 
  PSI110 = *( (Float *) ( (GIMM010) + (Umu) ) ) ; 
  PSI021 = *( (Float *) ( (GIMM001) + (Umu) ) ) ; 
  PSI111 = *( (Float *) ( (GIMM011) + (Umu) ) ) ; 
  PSI210 = *( (Float *) ( (GIMM100) + (Umu) ) ) ; 
  PSI200 = *( (Float *) ( (GIMM110) + (Umu) ) ) ; 
  PSI211 = *( (Float *) ( (GIMM101) + (Umu) ) ) ; 
  PSI201 = *( (Float *) ( (GIMM111) + (Umu) ) ) ; 
  Atmp000 = *( (float *) ( (CHI_IMM000) + (Chiin) ) ) ; 
  Atmp001 = *( (float *) ( (CHI_IMM001) + (Chiin) ) ) ; 
  Atmp010 = *( (float *) ( (CHI_IMM010) + (Chiin) ) ) ; 
  Atmp011 = *( (float *) ( (CHI_IMM011) + (Chiin) ) ) ; 
  Btmp000 = *( (float *) ( (CHI_IMM100) + (Chiin) ) ) ; 
  Btmp001 = *( (float *) ( (CHI_IMM101) + (Chiin) ) ) ; 
  Btmp010 = *( (float *) ( (CHI_IMM110) + (Chiin) ) ) ; 
  Btmp011 = *( (float *) ( (CHI_IMM111) + (Chiin) ) ) ; 
  length =   ( minus1 ) + ( length )  ;
  if ( length >  0 ) goto s_rec_su3_lab1;
  *( (float *) ( (CHI_IMM110) + (Chidrain) ) ) = PSI030; 
  *( (float *) ( (CHI_IMM111) + (Chidrain) ) ) = PSI031; 
  *( (float *) ( (CHI_IMM200) + (Chidrain) ) ) = PSI130; 
  *( (float *) ( (CHI_IMM201) + (Chidrain) ) ) = PSI131; 
  *( (float *) ( (CHI_IMM210) + (Chidrain) ) ) = PSI230; 
  *( (float *) ( (CHI_IMM211) + (Chidrain) ) ) = PSI231; 
s_rec_su3_lab0:

  return;
}
#ifdef __cplusplus
}
#endif
