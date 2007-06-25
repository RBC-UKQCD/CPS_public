#ifndef __QIOREADPROP__
#define __QIOREADPROP__

#include <util/qio_general.h>

CPS_START_NAMESPACE
using namespace std;


class qio_readPropagator: private qio_init {

 private:

  char *cname;

 public:

  qio_readPropagator(int argc, char *argv[]): 
    qio_init(argc, argv), cname("qio_readPropagator"), readProps(0), readSources(0), readSourceType(QIO_UNKNOWN_SOURCE){}


  // general reader
  qio_readPropagator(char *infile, void *prop, void *source, 
		     const int max_Float_prop, const int max_Float_source, 
		     int argc, char *argv[], int volFormat=QIO_UNKNOWN): 
    qio_init(argc, argv), cname("qio_readPropagator"), readProps(0), readSources(0), readSourceType(QIO_UNKNOWN_SOURCE)
    {read(infile, prop, source, max_Float_prop, max_Float_source, volFormat);}


  // read scS_12sink
  qio_readPropagator(char *infile, void *prop, void *source, 
		     int argc, char *argv[], int volFormat=QIO_UNKNOWN): 
    qio_init(argc, argv), cname("qio_readPropagator"), readProps(0), readSources(0), readSourceType(QIO_UNKNOWN_SOURCE)
    {read_ScS_12sink(infile, prop, source, volFormat);}

  // read 12pairs
   qio_readPropagator(char *infile, const QIO_PROP_SOURCE_TYPES sType, void *prop, void *source, 
		     int argc, char *argv[], int volFormat=QIO_UNKNOWN): 
    qio_init(argc, argv), cname("qio_readPropagator"), readProps(0), readSources(0), readSourceType(QIO_UNKNOWN_SOURCE)
     { read_12pairs(infile, prop, source, sType, volFormat); }

  virtual ~qio_readPropagator(){
    #ifdef DEBUG_Init
    printf("finished qio_readPropagator\n");
    #endif //DEBUG_Init
  }

  void read(char *infile, void *prop, void *source, int max_Float_prop, int max_Float_source, int volFormat=QIO_UNKNOWN);

  void read_ScS_12sink(char *infile, void *prop, void *source, int volFormat=QIO_UNKNOWN);

  void read_12pairs(char *infile, void *prop, void *source, const QIO_PROP_SOURCE_TYPES sType, int volFormat=QIO_UNKNOWN);

  int readProps;
  int readSources;
  QIO_PROP_SOURCE_TYPES readSourceType;  // 0=UNKNOWN, 1=Scalar, 2=full

 private:

  QIO_Reader *qio_Input;

  void qio_openInput(char *filename, QIO_String *xml_file_in, int volFormat);

  void qio_closeInput()
    { QIO_close_read(qio_Input);}

  void qio_communicateTchunks(void *start, const int sizepersite);

  int qio_readNextPropagatorRecord(int datacount, Float *rprop);
  
  int qio_readNextPropagatorRecordSingle(int datacount, Float *rprop);

  
  int qio_readNextSourcePairRecord(int datacount, Float *rprop, int &readSpin, int &readColor);

  int qio_readNextPropPairRecord(int datacount, Float *rprop, const int readSpin, const int readColor);

};



CPS_END_NAMESPACE
#endif // __QIOREADRPROP__



