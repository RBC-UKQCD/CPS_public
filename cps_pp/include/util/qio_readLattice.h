#ifndef __QIOREADLAT__
#define __QIOREADLAT__

#include <util/qio_general.h>

CPS_START_NAMESPACE
using namespace std;


class qio_readLattice: private qio_init {

 private:

  char *cname;

 public:

  qio_readLattice(int argc, char *argv[]): qio_init(argc, argv), cname("qio_readLattice"){}

  qio_readLattice(char *infile, Lattice &lat, int argc, char *argv[], int volFormat=QIO_UNKNOWN): qio_init(argc, argv), cname("qio_readLattice")
   {read(infile, lat, volFormat);}

  virtual ~qio_readLattice(){
    #ifdef DEBUG_Init
    printf("finished qio_readLattice\n");
    #endif //DEBUG_Init
  }

  void read(char *infile, Lattice &lat, int volFormat=QIO_UNKNOWN);

 private:

  QIO_Reader *qio_Input;

  void qio_openInput(char *filename, int volFormat);

  void qio_closeInput()
    { QIO_close_read(qio_Input);}

  void qio_communicateTchunks(void *start);

};



CPS_END_NAMESPACE
#endif // __QIOREADLAT__


