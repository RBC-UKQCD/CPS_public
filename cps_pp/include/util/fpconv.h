/*!
  \file
  \brief
  
  $Id: fpconv.h,v 1.10 2012-08-10 14:05:33 chulwoo Exp $
 */

#ifndef __FP_CONV__
#define __FP_CONV__

#include <config.h>
#include <stdint.h>


CPS_START_NAMESPACE

typedef uint32_t type32 ;
typedef uint64_t type64;


class QioArg;

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

const double             FPConv_PI       = 3.14159265358979323846264338327950288419716939937510;
//const unsigned long      FPConv_ieee32pi = 0x40490fdbL;
//const unsigned long long FPConv_ieee64pi = 0x400921fb54442d18LL;
const char FPConv_ieee32pi_big[4] = { '\x40', '\x49', '\x0f', '\xdb' };
const char FPConv_ieee64pi_big[8] = { '\x40', '\x09', '\x21', '\xfb',
				      '\x54', '\x44', '\x2d', '\x18' };


class QioArg;  // forward declaration

class DataConversion {
 public:
  virtual char * host2file(char * fbuf, const char * hbuf, const int dat_len) const = 0;
  virtual char * file2host(char * hbuf, const char * fbuf, const int dat_len) const = 0;
  virtual int fileDataSize() const = 0;
  virtual int hostDataSize() const = 0;
  virtual unsigned int checksum(char * data, const int data_len) const = 0;
  virtual unsigned int posDepCsum(char * data, const int data_len, const int dimension,
				  const QioArg & qio_arg, 
				  const int siteid, const int global_id) const = 0;
  DataConversion() { }
  virtual ~DataConversion() { }
};


class FPConv : public DataConversion
{
 private:
  char * cname;
 protected:
  void ti2ieee (type32 *,int) const;
  void ieee2ti (type32 *,int) const;

  void byterevn(type32 w[], int n) const;
  void byterevn64(type64 w[], int n) const;

  void conv64to32(type32 tgt[], type64 src[], int n) const;
  void conv32to64(type64 tgt[], type32 src[], int n) const;
  void copy64(type64 tgt[], type64 src[], int n) const;
  void copy32(type32 tgt[], type32 src[], int n) const;
  enum FP_FORMAT  testHostFormat();

 public:
  enum FP_FORMAT  hostFormat, fileFormat; // open to public so I can broadcast

  FPConv();
  virtual ~FPConv();
  char * file2host(char * hbuf, const char *fdat, const int fdat_len) const;
  char * host2file(char * fbuf, const char *hdat, const int hdat_len) const;

  enum FP_FORMAT  setFileFormat(const char * desc);
  enum FP_FORMAT  setFileFormat(const enum FP_FORMAT dataFormat);
  unsigned int checksum(char * data, const int data_len, 
			const enum FP_FORMAT dataFormat) const;

  int size(const enum FP_FORMAT datatype) const ;
  bool big_endian(const enum FP_FORMAT datatype) const;
  inline int fileFpSize() const { return size(fileFormat); }
  inline int hostFpSize() const { return size(hostFormat); }
  static const char * name(const enum FP_FORMAT format);

  int fileDataSize() const { return fileFpSize(); }
  int hostDataSize() const { return hostFpSize(); }
  unsigned int checksum(char * data, const int data_len) const {
    return checksum(data, data_len, FP_AUTOMATIC);
  }
  // not used in FPConv
  unsigned int posDepCsum(char * data, const int data_len, const int dimension,
			  const QioArg & qio_arg, 
			  const int site_id, const int global_id) const { return 0; }

 private:
  bool sim_qcdsp;
 public:
  void SimQCDSP(int sim) { sim_qcdsp = sim?true:false; }
};

CPS_END_NAMESPACE

#endif
