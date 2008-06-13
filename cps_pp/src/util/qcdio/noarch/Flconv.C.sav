//
//   Flconv.C
//               convert   IEEE32BIG <->  TIDSP



#include <util/qcdio.h>
#include <stdlib.h> // exit()
#include <util/Flconv.h>
//
// switch between big-endian and little endian
//

void byterevn(type32 w[], int n)
{
  register type32 old,newv;
  int j;
  
  for(j=0; j<n; j++)
    {
      old = w[j];
      newv = old >> 24 & 0x000000ff;
      newv |= old >> 8 & 0x0000ff00;
      newv |= old << 8 & 0x00ff0000;
      newv |= old << 24 & 0xff000000;
      w[j] = newv;
    }
} 

// TMS320C30-IEEE Floating-Point Format Converter
// c.f Texas Instruments Application Report  SPR400
//    N.B. This report has several typos.

// x-th bit of tmp 
#define tmpB(x)  ( tmp>>(x)&1 )
//#define DEB(x)  printf("%s\n",(x))
#define DEB(x)  


void ti2ieee(type32 *ti, int Num){
  register type32 tmp, tisave;
  
  type32 sign(0);
  type32 expo(0);
  type32 sfct(0);

  type32 sign_ti;
  type32 expo_ti;
  type32 sfct_ti;
  type32 EXP80_81, EXP7F, MANT0;

  int i;

  for(i=0;i<Num;i++){
    tmp = *(ti+i);

    tisave = tmp;

    // divide TIFP into sign_ti/expo_ti/sfct_ti
    // TIFP format :
    //  31             24     23     22                     0
    // |    expo          | sign  |         sfct             |
    sign_ti= tmpB(23);
    expo_ti= tmp >> 24;
    sfct_ti= tmp & 0x007fffff;

    EXP80_81 = (expo_ti == 0x80) || (expo_ti == 0x81);
    EXP7F = expo_ti == 0x7f;
    MANT0 =  (sfct_ti == 0) ;

//    printf("expo_ti = %d EXP80_81 = %d",expo_ti,EXP80_81);


    // CASE6  Positive number >= 2^{-126}
    //  = !EXP80_81 & !ti(23)
    if( !EXP80_81 && !sign_ti ){
      DEB("C6");
      sign= sign_ti;
      expo= expo_ti + 0x7f;
      sfct= sfct_ti; 
    }
   // CASE7  Positive number [2^{-127},  (2-2^{-23}) 2^{-127}]
    else if( EXP80_81 && tmpB(24) && !sign_ti ){
      DEB("C7");
      sign= sign_ti;
      expo= 0;
      sfct= sfct_ti/2+0x400000; 
    }
   // CASE8  zero
    else if( EXP80_81 && !tmpB(24) ){
      DEB("C8");
      sign= 0;
      expo= 0;
      sfct= 0;
    }
   // CASE9  Negative number [(-2-2^{-23}) 2^{-127}, (-1-2^{-23}) 2^{-127}]
    else if( EXP80_81 && tmpB(23) && !MANT0 ){
      DEB("C9");
      sign= sign_ti;
      expo= 0;
      sfct= (0x800000-sfct_ti)/2+0x400000; 
    }
   // CASE10  Negative number [- 2^{127}, -2^{-126}]
    else if( !(EXP80_81 && !tmpB(24))  && !EXP7F && tmpB(23) && MANT0 ){
      DEB("C10");
      sign= sign_ti;
      expo= expo_ti+0x80;
      sfct= 0;
    }
   // CASE11  Negative number (- 2^{128}, -2^{-126})
    else if( !EXP80_81 && sign_ti && !MANT0 ){
      DEB("C11");
      sign= sign_ti;
      expo= expo_ti+0x7f;
      sfct= 0x800000-sfct_ti;
    }
   // CASE12  Negative number - 2^{128}
    else if( EXP7F && sign_ti && MANT0 ){
      DEB("C12");
      sign= sign_ti;
      expo= 0xff;
      sfct= 0;
    }
    else {
      DEB("UNSUPPORTED");
      printf("unsupported bit pattern %x\n",(int)tmp);
      exit(-13);
    }


    // IEEE FP format :
    //    31     30           23   22                 0
    // | sign  |     expo        |     sfct             |
    sign = sign << 31;
    expo = (expo << 23) & 0x7f800000;
    sfct = sfct & 0x7fffff;

    *(ti+i) = sign | expo | sfct;
  }

}

