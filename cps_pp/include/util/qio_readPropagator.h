#ifndef __QIOREADPROP__
#define __QIOREADPROP__

#include <util/qio_general.h>

CPS_START_NAMESPACE
using namespace std;


class qio_readPropagator: private qio_init {

 private:

  char *cname;

 public:

  // version with argc, argv

  qio_readPropagator(int argc, char *argv[]): 
    qio_init(argc, argv), cname("qio_readPropagator"), readProps(0), readSources(0), readSourceType(QIO_UNKNOWN_SOURCE), lastColor(-1), lastSpin(0){}


  //! read propagator, auto-detect format
  /*!
    \param infile file to read from
    \param prop   allocated memory for propagator
    \param source allocated memory for source
    \param max_Float_prop maximal number of Floats allocated for propagator
    \param max_Float_source maximal number of Floats allocated for source
  */
  qio_readPropagator(char *infile, void *prop, void *source, 
		     const int max_Float_prop, const int max_Float_source, 
		     int argc, char *argv[], int volFormat=QIO_UNKNOWN): 
    qio_init(argc, argv), cname("qio_readPropagator"), readProps(0), readSources(0), readSourceType(QIO_UNKNOWN_SOURCE), lastColor(-1), lastSpin(0)
    {read(infile, prop, source, max_Float_prop, max_Float_source, volFormat);}


  //! read propagator, format scalar source + 12 sinks
  /*!
    \param infile file to read from
    \param prop   allocated memory for propagator
    \param source allocated memory for source
  */
  qio_readPropagator(char *infile, void *prop, void *source, 
		     int argc, char *argv[], int volFormat=QIO_UNKNOWN): 
    qio_init(argc, argv), cname("qio_readPropagator"), readProps(0), readSources(0), readSourceType(QIO_UNKNOWN_SOURCE), lastColor(-1), lastSpin(0)
    {read_ScS_12sink(infile, prop, source, volFormat);}

  //! read propagator 12 pairs of sink and source
  /*!
    \param infile file to read from
    \param sType  source type (full, scalar)
    \param prop   allocated memory for propagator
    \param source allocated memory for source
   */
  qio_readPropagator(char *infile, const QIO_PROP_SOURCE_TYPES sType, void *prop, void *source, 
		     int argc, char *argv[], int volFormat=QIO_UNKNOWN): 
    qio_init(argc, argv), cname("qio_readPropagator"), readProps(0), readSources(0), readSourceType(QIO_UNKNOWN_SOURCE), lastColor(-1), lastSpin(0)
    { read_12pairs(infile, prop, source, sType, volFormat); }
  
  // versions w/o argc/argv, using GJP

  qio_readPropagator(): 
    qio_init(GJP.argc(), GJP.argv()), cname("qio_readPropagator"), readProps(0), readSources(0), readSourceType(QIO_UNKNOWN_SOURCE), lastColor(-1), lastSpin(0){}

  //! read propagator, auto-detect format
  /*!
    \param infile file to read from
    \param prop   allocated memory for propagator
    \param source allocated memory for source
    \param max_Float_prop maximal number of Floats allocated for propagator
    \param max_Float_source maximal number of Floats allocated for source
  */
  qio_readPropagator(char *infile, void *prop, void *source, 
		     const int max_Float_prop, const int max_Float_source, 
		     int volFormat=QIO_UNKNOWN): 
    qio_init(GJP.argc(), GJP.argv()), cname("qio_readPropagator"), readProps(0), readSources(0), readSourceType(QIO_UNKNOWN_SOURCE), lastColor(-1), lastSpin(0)
    {read(infile, prop, source, max_Float_prop, max_Float_source, volFormat);}

  qio_readPropagator(char *infile,  
		     int spin, int color,
		     void *prop, void *source,
		     const int max_Float_prop, const int max_Float_source, 
		     int volFormat=QIO_UNKNOWN): 
    qio_init(GJP.argc(), GJP.argv()), cname("qio_readPropagator"), readProps(0), readSources(0), readSourceType(QIO_UNKNOWN_SOURCE), lastColor(-1), lastSpin(0)
      {read(infile, spin, color, prop, source, max_Float_prop, max_Float_source, volFormat);}

  //! read propagator, format scalar source + 12 sinks
  /*!
    \param infile file to read from
    \param prop   allocated memory for propagator
    \param source allocated memory for source
  */
  qio_readPropagator(char *infile, void *prop, void *source, 
		     int volFormat=QIO_UNKNOWN): 
    qio_init(GJP.argc(), GJP.argv()), cname("qio_readPropagator"), readProps(0), readSources(0), readSourceType(QIO_UNKNOWN_SOURCE), lastColor(-1), lastSpin(0)
    {read_ScS_12sink(infile, prop, source, volFormat);}

  //! read propagator 12 pairs of sink and source
  /*!
    \param infile file to read from
    \param sType  source type (full, scalar)
    \param prop   allocated memory for propagator
    \param source allocated memory for source
  */
  qio_readPropagator(char *infile, const QIO_PROP_SOURCE_TYPES sType, void *prop, void *source, 
		     int volFormat=QIO_UNKNOWN): 
    qio_init(GJP.argc(), GJP.argv()), cname("qio_readPropagator"), readProps(0), readSources(0), readSourceType(QIO_UNKNOWN_SOURCE), lastColor(-1), lastSpin(0)
    { read_12pairs(infile, prop, source, sType, volFormat); }
  


  virtual ~qio_readPropagator(){
    #ifdef DEBUG_Init
    printf("finished qio_readPropagator\n");
    #endif //DEBUG_Init
  }

  //! read propagator, auto-detect format
  /*!
    \param infile file to read from
    \param prop   allocated memory for propagator
    \param source allocated memory for source
    \param max_Float_prop maximal number of Floats allocated for propagator
    \param max_Float_source maximal number of Floats allocated for source
  */
  void read(char *infile, void *prop, void *source, int max_Float_prop, int max_Float_source, int volFormat=QIO_UNKNOWN);
  void read(char *infile, int spin, int color, void *prop, void *source, int max_Float_prop, int max_Float_source, int volFormat=QIO_UNKNOWN);

  //! read propagator, format scalar source + 12 sinks
  /*!
    \param infile file to read from
    \param prop   allocated memory for propagator
    \param source allocated memory for source
  */
  void read_ScS_12sink(char *infile, void *prop, void *source, int volFormat=QIO_UNKNOWN);

  //! read propagator 12 pairs of sink and source
  /*!
    \param infile file to read from
    \param sType  source type (full, scalar)
    \param prop   allocated memory for propagator
    \param source allocated memory for source
  */
  void read_12pairs(char *infile, void *prop, void *source, const QIO_PROP_SOURCE_TYPES sType, int volFormat=QIO_UNKNOWN);
  void read_pair(char *infile, int spin, int color, void *prop, void *source, const QIO_PROP_SOURCE_TYPES sType, int volFormat=QIO_UNKNOWN);

  //! number of sinks read
  int readProps;
  //! number of sources read
  int readSources;
  //! type of source read
  QIO_PROP_SOURCE_TYPES readSourceType;  // 0=UNKNOWN, 1=Scalar, 2=full

 private:

  QIO_Reader *qio_Input;

  // in case the data does not specify spin, color: we assume spin(<4), color(<3) ordering
  int lastColor, lastSpin;
  void qio_guessSpinColor(int &spin, int &color);
  
  

  void qio_openInput(char *filename, QIO_String *xml_file_in, int volFormat);

  void qio_closeInput()
    { QIO_close_read(qio_Input); }

  void qio_communicateTchunks(void *start, const int sizepersite);

  int qio_readNextPropagatorRecord(Float *rprop);
  
  int qio_readNextPropagatorRecordSingle(Float *rprop);

  int qio_readNextPropPairRecord(Float *rprop, int &readSpin, int &readColor);

  int qio_readTmpSourceRecord(const QIO_PROP_SOURCE_TYPES sType, Float *rsource);

  public:
  int hyper_n;
  int hyper_lower[4];
  int hyper_upper[4];

};



CPS_END_NAMESPACE
#endif // __QIOREADRPROP__



