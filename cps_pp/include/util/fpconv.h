#ifndef __FP_CONV__
#define __FP_CONV__

#include <config.h>

CPS_START_NAMESPACE
using namespace std;

typedef unsigned long int type32 ;
typedef unsigned long long type64;


enum FP_FORMAT {
  FP_UNKNOWN = 0,
  FP_AUTOMATIC, // used as arguments to function, 
                //will be substituted by hostFormat or fileFormat depending on situation
  FP_TIDSP32,
  FP_IEEE32,
  FP_IEEE32BIG,
  FP_IEEE32LITTLE,
  FP_IEEE64,
  FP_IEEE64BIG,
  FP_IEEE64LITTLE
};

const double             FPConv_PI       = 3.14159265358979323846264338327950288319716939937510;
const unsigned long      FPConv_ieee32pi = 0x40490fdb;
const unsigned long long FPConv_ieee64pi = 0x400921fb54442d18;

class FPConv {
 protected:
  void ti2ieee (type32 *,int);
  void ieee2ti (type32 *,int);

  void byterevn(type32 w[], int n);
  void byterevn64(type64 w[], int n);

  void conv64to32(type32 tgt[], type64 src[], int n);
  void conv32to64(type64 tgt[], type32 src[], int n);
  void copy64(type64 tgt[], type64 src[], int n);
  void copy32(type32 tgt[], type32 src[], int n);
  enum FP_FORMAT  testHostFormat();

 public:
  enum FP_FORMAT  hostFormat, fileFormat; // open to public so I can broadcast

  FPConv();
  ~FPConv();
  char * file2host(char * hbuf, const char *fdat, const int fdat_len);
  char * host2file(char * fbuf, const char *hdat, const int hdat_len);

  enum FP_FORMAT  setFileFormat(const char * desc);
  enum FP_FORMAT  setFileFormat(const enum FP_FORMAT dataFormat);
  unsigned long checksum(char * data, const int data_len, 
			 const enum FP_FORMAT dataFormat = FP_AUTOMATIC); // default use fileFormat

  int size(const enum FP_FORMAT datatype);
  bool big_endian(const enum FP_FORMAT datatype);
  inline int fileFpSize() { return size(fileFormat); }
  inline int hostFpSize() { return size(hostFormat); }
  const char * name(const enum FP_FORMAT format);
};

CPS_END_NAMESPACE

#endif
