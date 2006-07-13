#ifndef __QIOWRITELAT__
#define __QIOWRITELAT__

#include <util/qio_general.h>
#include <util/fpconv.h>


CPS_START_NAMESPACE
using namespace std;


class qio_writeLattice: private qio_init {

 private:

  char *cname;

 public:

  qio_writeLattice( int argc, char *argv[]): qio_init(argc, argv), cname("qio_writeLattice"){}

  qio_writeLattice(char *outfile, Lattice &lat, int argc, char *argv[], FP_FORMAT floatFormat=FP_AUTOMATIC):
   qio_init(argc, argv), cname("qio_writeLattice")
    { write(outfile, lat, floatFormat);}

  virtual ~qio_writeLattice(){ 
    #ifdef DEBUG_Init
    printf("finished: qio_writeLattice\n");
    #endif // DEBUG_Init
  }

  void write(char *outfile, Lattice &lat, FP_FORMAT floatFormat=FP_AUTOMATIC);
 
 private:

  QIO_Writer *qio_Output;

  void qio_openOutput(char *filename, char *stringLFN, char *xml_write_file);

  void qio_closeOutput()
    { QIO_close_write(qio_Output);}

    

};



CPS_END_NAMESPACE
#endif // __QIOWRITELAT__
