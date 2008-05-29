#include <config.h>
#include <string.h>
#include <iostream>
#include <util/data_types.h>
#include <util/fpconv.h>
#include <util/qioarg.h>
#include <sys/types.h>

CPS_START_NAMESPACE
using namespace std;

const char * FP_FORMAT_NAME[] = { 
  "n/a",
  "AUTOMATIC",
  "TIDSP32",
  "IEEE32",
  "IEEE32BIG",
  "IEEE32LITTLE",
  "IEEE64",
  "IEEE64BIG",
  "IEEE64LITTLE"
};

const int FP_FORMAT_ENTRIES = sizeof(FP_FORMAT_NAME)/sizeof(FP_FORMAT_NAME[0]);


FPConv::FPConv() 
  : fileFormat(FP_UNKNOWN), cname("FPConv"), sim_qcdsp(0) { 
  testHostFormat();

}

FPConv::~FPConv() {
}


const char * FPConv::name(const enum FP_FORMAT format) {
  return FP_FORMAT_NAME[int(format)];
}

char * FPConv::file2host(char * hbuf, const char * fdat, const int fdat_len)const {
  const char * fname = "file2host()";
  // trivial case
  if(hostFormat == fileFormat) {
    memcpy(hbuf,fdat,fdat_len*size(hostFormat));
    return hbuf;
  }
  
  // need conversion
  // 1. adjust endian
  if(big_endian(hostFormat) != big_endian(fileFormat)) {
    if(size(fileFormat) == 8) 
      byterevn64((type64*)fdat,fdat_len);
    else
      byterevn((type32*)fdat,fdat_len);
  }

  // 2. distinguish by host format
  if(hostFormat == FP_TIDSP32) {          // host  ti
    if(size(fileFormat) == 4) {
      ieee2ti((type32*)fdat, fdat_len);
      copy32((type32*)hbuf, (type32*)fdat, fdat_len);
    }
    else {
      ERR.NotImplemented(cname,fname,"Read IEEE64 on QCDSP not implemented!\n");
      //      return 0;
    }
  }
  else if(size(hostFormat) == 4) {     // host  ieee32
    if(fileFormat==FP_TIDSP32) {
      ti2ieee((type32*)fdat, fdat_len);
      copy32((type32*)hbuf, (type32*)fdat, fdat_len);
    }
    else if(size(fileFormat) == 4)
      copy32((type32*)hbuf, (type32*)fdat, fdat_len);
    else
      conv64to32((type32*)hbuf, (type64*)fdat, fdat_len);
  }
  else  {                              // host  ieee64
    if(fileFormat == FP_TIDSP32) {
      ti2ieee((type32*)fdat, fdat_len);
      conv32to64((type64*)hbuf, (type32*)fdat, fdat_len);
    }
    else if(size(fileFormat) == 4)
      conv32to64((type64*)hbuf, (type32*)fdat, fdat_len);
    else
      copy64((type64*)hbuf, (type64*)fdat, fdat_len);
  }

  return hbuf;
}

char * FPConv::host2file(char *fbuf, const char * hdat, const int hdat_len) const{
  const char * fname = "host2file";
  // trivial case
  if(hostFormat == fileFormat) {
    memcpy(fbuf,hdat,hdat_len*size(fileFormat));
    return fbuf;
  }
  
  // need conversion
  // 1. distinguish by host format
  if(hostFormat == FP_TIDSP32) {          // host  ti
    if(size(fileFormat) == 4) {
      copy32((type32*)fbuf, (type32*)hdat, hdat_len);
      ti2ieee((type32*)fbuf, hdat_len);
    }
    else {
      ERR.NotImplemented(cname,fname, "Write IEEE64 on QCDSP : not implemented!\n");
      //      return 0;
    }
  }
  else if(size(hostFormat) == 4) {     // host  ieee32
    if(fileFormat==FP_TIDSP32) {
      copy32((type32*)fbuf, (type32*)hdat, hdat_len);
      ieee2ti((type32*)fbuf, hdat_len);
    }
    else if(size(fileFormat) == 4)
      copy32((type32*)fbuf, (type32*)hdat, hdat_len);
    else
      conv32to64((type64*)fbuf, (type32*)hdat, hdat_len);
  }
  else  {                              // host  ieee64
    if(fileFormat == FP_TIDSP32) {
      conv64to32((type32*)fbuf, (type64*)hdat, hdat_len);
      ieee2ti((type32*)fbuf, hdat_len);
    }
    else if(size(fileFormat) == 4)
      conv64to32((type32*)fbuf, (type64*)hdat, hdat_len);
    else
      copy64((type64*)fbuf, (type64*)hdat, hdat_len);
  }

  // 2. adjust endian
  if(big_endian(hostFormat) != big_endian(fileFormat)) {
    if(size(fileFormat) == 8) 
      byterevn64((type64*)fbuf,hdat_len);
    else
      byterevn((type32*)fbuf,hdat_len);
  }

  return fbuf;
}


void FPConv::byterevn(type32 w[], int n) const{
  /*  char * buf = (char*)w;
  cout << "First 16 bytes: ";
  for(int i=0;i<16;i++) cout << hex << (unsigned int)buf[i] << " ";
  cout << dec << endl;
  */
  //  cout << "Byte reverse 32 bits" << endl;

  register type32 oldv, newv;

  for(int i=0;i<n;i++) {
    oldv = w[i];
    newv = 0;
    for(int j=0;j<4;j++) {
      newv = (newv << 8) | (oldv & 0xff);
      oldv >>= 8;
    }
    w[i] = newv;
  }
  /*
  cout << "First 16 bytes: ";
  for(int i=0;i<16;i++) cout << hex << (unsigned int)buf[i] << " ";
  cout << dec << endl;
  */
}

void FPConv::byterevn64(type64 w[], int n) const{
  /*
  char * buf = (char*)w;
  cout << "First 16 bytes: ";
  for(int i=0;i<16;i++) cout << hex << (unsigned int)buf[i] << " ";
  cout << dec << endl;
  */
  //  cout << "Byte reverse 64 bits" << endl;

  register type64 oldv, newv;

  for(int i=0;i<n;i++) {
    oldv = w[i];
    newv = 0;
    for(int j=0;j<8;j++) {
      newv = (newv << 8) | (oldv & 0xff);
      oldv >>= 8;
    }
    w[i] = newv;
  }
  /*
  cout << "First 16 bytes: ";
  for(int i=0;i<16;i++) cout << hex << (unsigned int)buf[i] << " ";
  cout << dec << endl;
  */
}

void FPConv::conv64to32(type32 tgt[], type64 src[], int n) const{
  //  cout << "Conversion 64 ==> 32" << endl;
  double *s = (double*)src;
  float  *t = (float*)tgt;
  for(int i=0;i<n;i++)  *t++ = *s++;
}

void FPConv::conv32to64(type64 tgt[], type32 src[], int n) const{
  //  cout << "Convesion 32 ==> 64 " << endl;
  float  *s = (float*) src;
  double *t = (double*)tgt;
  for(int i=0;i<n;i++)  *t++ = *s++;
}

void FPConv::copy64(type64 tgt[], type64 src[], int n) const{
  double *s = (double*)src;
  double *t = (double*)tgt;
  for(int i=0;i<n;i++)  *t++ = *s++;
}

void FPConv::copy32(type32 tgt[], type32 src[], int n) const{
  float *s = (float*)src;
  float *t = (float*)tgt;
  for(int i=0;i<n;i++)  *t++ = *s++;
}

enum FP_FORMAT  FPConv::testHostFormat() { // test the type of CPS::Float
  const char * fname = "testHostFormat()";
  // 1. endian
  char end_check[4] = {1,0,0,0};
//  unsigned long *lp = (unsigned long *)end_check;
  uint32_t *lp = (uint32_t *)end_check;
  int host_big;

  VRB.Flow(cname,fname,"size(uint32_t)=%d lp=%x",sizeof(uint32_t),*lp);
  if ( *lp == 0x1 ) { 
    //    cout << "Host is little-endian\n";
    host_big = 0;
  } else {
    //    cout << "Host is big-endian\n";
    host_big = 1;
  }

  // 2. pi test
  if(sizeof(Float) == 8) {  // 64 bits
    if(host_big)  hostFormat = FP_IEEE64BIG;
    else          hostFormat = FP_IEEE64LITTLE;
  }
  else {  // 32 bits
    union { 
      float pinum;
      char pichar[4];
    }cpspi;

    cpspi.pinum = FPConv_PI;
    if(host_big) {
      hostFormat = FP_IEEE32BIG;
      for(int i=0;i<4;i++) {
	if(cpspi.pichar[i] != FPConv_ieee32pi_big[i]) {
	  hostFormat = FP_TIDSP32;
	  break;
	}
      }
    }
    else {
      hostFormat = FP_IEEE32LITTLE;
      for(int i=0;i<4;i++) {
	if(cpspi.pichar[i] != FPConv_ieee32pi_big[3-i]) {
	  hostFormat = FP_TIDSP32;
	  break;
	}
      }
    }

  } // end of 32 bits
    VRB.Flow(cname,fname, "Host FP Format : %s\n", name(hostFormat) );

  return hostFormat;
}

enum FP_FORMAT  FPConv::setFileFormat(const enum FP_FORMAT dataFormat) {
  const char * fname = "setFileFormat";
  fileFormat = dataFormat;
  if(dataFormat == FP_AUTOMATIC) 
    fileFormat = hostFormat;
  if(fileFormat == FP_UNKNOWN) {
    VRB.Flow(cname,fname,"Floating point format cannot be FP_UNKNOWN!\n");
  }
  return fileFormat;
}

enum FP_FORMAT  FPConv::setFileFormat(const char * desc) {
  const char * fname = "setFileFormat()";
  fileFormat = FP_UNKNOWN;
  for(int i=1;i<FP_FORMAT_ENTRIES;i++) {
    if(!strcmp(FP_FORMAT_NAME[i],desc)) {
      fileFormat = FP_FORMAT(i);
      break;
    }
  }
  if(fileFormat == FP_UNKNOWN) {
    VRB.Flow(cname,fname, "Floating point format \"%s\" not recognized!\n",desc);
  }
  else if(fileFormat == FP_AUTOMATIC)  {
    fileFormat = hostFormat;
  }

  return fileFormat;
}

unsigned int FPConv::checksum(char * data, const int data_len,
			       const enum FP_FORMAT dataFormat) const{
  const char * fname = "checksum()";
  // checksum always done on 32-bits

  enum FP_FORMAT chkFormat = dataFormat;
  if(dataFormat == FP_AUTOMATIC)  chkFormat = fileFormat;

  if(chkFormat == FP_UNKNOWN) {
    VRB.Flow(cname,fname, "checksum data format UNKNOWN!\n");
    return 0;
  }

  int csumcnt = data_len;
  if(size(chkFormat) == 8)  csumcnt *= 2; // always sum as 32 bits data

  int need_byterevn = 0;
  if(sim_qcdsp && !big_endian(chkFormat)) 
    need_byterevn = 1;
  if(!sim_qcdsp && big_endian(hostFormat) != big_endian(chkFormat))
    need_byterevn = 1;

  uint32_t *buf = (uint32_t*)data;
//  VRB.Result(cname,fname,"data = %x %x",*buf,*(buf+1));
  if(need_byterevn)
    byterevn((type32*)data, csumcnt);
//  VRB.Result(cname,fname,"data = %x %x",*buf,*(buf+1));

  uint32_t s = 0;
  for(int i=0;i<csumcnt;i++)  {
    s += *buf;
//	if(i%1000<2)
//       VRB.Result(cname,fname,"i=%d s=%x",i,s);
    buf++;
  }

  if(need_byterevn)
    byterevn((type32*)data, csumcnt);

//  VRB.Result(cname,fname," need_byterevn=%d s=%x",need_byterevn,s);
  return s;
}

int FPConv::size(const enum FP_FORMAT datatype) const {
  switch(datatype) {
  case FP_TIDSP32:
  case FP_IEEE32:
  case FP_IEEE32BIG:
  case FP_IEEE32LITTLE:
    return 4;
  case FP_IEEE64:
  case FP_IEEE64BIG:
  case FP_IEEE64LITTLE:
    return 8;
  default:
    return 0;
  }
}

bool FPConv::big_endian(const enum FP_FORMAT datatype) const {
  switch(datatype) {
  case FP_TIDSP32:
  case FP_IEEE32:
  case FP_IEEE32BIG:
  case FP_IEEE64:
  case FP_IEEE64BIG:
    return true;
  case FP_IEEE32LITTLE:
  case FP_IEEE64LITTLE:
    return false;
  default:
    return false;
  }
}


// as this is no longer integrated in QCDOC release of CPS
// i copy the code here to do  TIDSP32 ==> IEEE32

// TMS320C30-IEEE Floating-Point Format Converter
// c.f Texas Instruments Application Report  SPR400
//    N.B. This report has several typos.

// x-th bit of tmp 
#define tmpB(x)  ( tmp>>(x)&1 )
//#define DEB(x)  printf("%s\n",(x))
#define DEB(x)  


void FPConv::ti2ieee(type32 *ti, int Num) const{
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


void FPConv::ieee2ti(type32 *ti, int Num) const{
  const char * fname = "ieee2ti()";
  ERR.NotImplemented(cname,fname);
  exit(13);
}


CPS_END_NAMESPACE
