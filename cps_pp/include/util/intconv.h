#ifndef __INT_CONV__
#define __INT_CONV__

#include <config.h>

CPS_START_NAMESPACE
using namespace std;

typedef unsigned long int type32;
typedef unsigned long long type64; // maybe later we'll have 64-bit int...


enum INT_FORMAT {
  INT_UNKNOWN = 0,
  INT_AUTOMATIC, // used as arguments to function, 
                //will be substituted by hostFormat or fileFormat depending on situation
  INT_32BIG,
  INT_32LITTLE
};


class IntConv {
 protected:
  void byterevn(type32 w[], int n);

  void copy32(type32 tgt[], type32 src[], int n);
  enum INT_FORMAT  testHostFormat();

 public:
  enum INT_FORMAT  hostFormat, fileFormat; // open to public so I can broadcast them
  enum INT_FORMAT  setFileFormat(const char * desc);
  enum INT_FORMAT  setFileFormat(const enum INT_FORMAT dataFormat);

  IntConv();
  virtual ~IntConv();
  char * file2host(char * hbuf, const char *fdat, const int fdat_len);
  char * host2file(char * fbuf, const char *hdat, const int hdat_len);
  unsigned int checksum(char * data, const int data_len, 
			const enum INT_FORMAT dataFormat = INT_AUTOMATIC); // default use fileFormat
  unsigned int posDepCsum(char * data, const int data_len, 
			  const enum INT_FORMAT dataFormat = INT_AUTOMATIC);

  // functions to describe int format
  int size(const enum INT_FORMAT datatype);
  inline int fileIntSize() { return size(fileFormat); }
  inline int hostIntSize() { return size(hostFormat); }
  bool big_endian(const enum INT_FORMAT datatype);
  const char * name(const enum INT_FORMAT format);
};

CPS_END_NAMESPACE

#endif
