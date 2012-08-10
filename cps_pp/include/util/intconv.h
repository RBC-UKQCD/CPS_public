#ifndef __INT_CONV__
#define __INT_CONV__

#include <config.h>
#include <util/fpconv.h> // need base class DataConversion

CPS_START_NAMESPACE

typedef uint32_t type32;
typedef uint64_t type64; // maybe later we'll have 64-bit int...


enum INT_FORMAT {
  //INT_UNKNOWN = 0,
  INT_UNKNOWN,
  INT_AUTOMATIC, // used as arguments to function, 
                //will be substituted by hostFormat or fileFormat depending on situation
  INT_32BIG,
  INT_32LITTLE
};


class DataConversion;  // forward declaration 
class QioArg;


class IntConv : public DataConversion
{
 private:
  char * cname;
 protected:
  void byterevn(type32 w[], int n) const;

  void copy32(type32 tgt[], type32 src[], int n) const;
  enum INT_FORMAT  testHostFormat();

 public:
  enum INT_FORMAT  hostFormat, fileFormat; // open to public so I can broadcast them
  enum INT_FORMAT  setFileFormat(const char * desc);
  enum INT_FORMAT  setFileFormat(const enum INT_FORMAT dataFormat);

  IntConv();
  virtual ~IntConv();
  char * file2host(char * hbuf, const char *fdat, const int fdat_len) const;
  char * host2file(char * fbuf, const char *hdat, const int hdat_len) const;
  unsigned int checksum(char * data, const int data_len, 
			const enum INT_FORMAT dataFormat) const; // default use fileFormat
  unsigned int posDepCsum(char * data, const int data_len, const int dimension,
			  const QioArg & qio_arg, const int siteid, const int global_id,
			  const enum INT_FORMAT dataFormat) const;

  // functions to describe int format
  int size(const enum INT_FORMAT datatype) const;
  inline int fileIntSize() const { return size(fileFormat); }
  inline int hostIntSize() const { return size(hostFormat); }
  bool big_endian(const enum INT_FORMAT datatype) const;
  static const char * name(const enum INT_FORMAT format);

  int fileDataSize() const { return fileIntSize(); }
  int hostDataSize() const { return hostIntSize(); }
  unsigned int checksum(char * data, const int data_len) const {
    return checksum(data, data_len, INT_AUTOMATIC);
  }
  unsigned int posDepCsum(char * data, const int data_len, const int dimension,
			  const QioArg & qio_arg, 
			  const int siteid, const int global_id) const {
    return posDepCsum(data, data_len, dimension, qio_arg, siteid, global_id, INT_AUTOMATIC);
  }
};

CPS_END_NAMESPACE

#endif
