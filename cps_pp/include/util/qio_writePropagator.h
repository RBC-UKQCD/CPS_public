#ifndef __QIOWRITEPROP__
#define __QIOWRITEPROP__

#include <util/qio_general.h>
#include <util/fpconv.h>


CPS_START_NAMESPACE
using namespace std;


class qio_writePropagator: private qio_init {

 private:

  char *cname;

 public:

  // version with argc/argv

  qio_writePropagator( int argc, char *argv[]): qio_init(argc, argv), cname("qio_writePropagator"), source_hypercube(0)
    {initHeader();}


  //! write propagator, format scalar source + 12 sinks
  /*!
    \param outfile file to write to
    \param prop    pointer to sinks
    \param scSource pointer to scalar source
   */
  qio_writePropagator(char *outfile, const void *prop, const void *scSource, 
		      int argc, char *argv[], int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(argc, argv), cname("qio_writePropagator"), source_hypercube(0)
    { initHeader(); write_ScS_12sink(outfile, prop, scSource, volFormat, floatFormat);}

   
  //! write propagator, format scalar source + 12 sinks, set add. info
  /*!
    \param outfile file to write to
    \param prop    pointer to sinks
    \param scSource pointer to scalar source
    \param ensemble_id ensemble id
    \param ensemble_label ensemble label
    \param traj trajectory number
   */
  qio_writePropagator(char *outfile, const void *prop, const void *scSource, 
		      const char * ensemble_id, const char * ensemble_label, const int traj,
		      int argc, char *argv[], int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(argc, argv), cname("qio_writePropagator"), source_hypercube(0)
    { setHeader(ensemble_id, ensemble_label, traj); write_ScS_12sink(outfile, prop, scSource, volFormat, floatFormat);}

  //! write propagator, format scalar source + 12 sinks, set add. info
  /*!
    \param outfile file to write to
    \param prop    pointer to sinks
    \param scSource pointer to scalar source
    \param ensemble_id ensemble id
    \param ensemble_label ensemble label
    \param traj trajectory number
    \param propagator_type descr. of prop. type
    \param source_type descr. of source type
   */
  qio_writePropagator(char *outfile, const void *prop, const void *scSource, 
		      const char * ensemble_id, const char * ensemble_label, const int traj,
		      const char * propagator_type, const char * source_type,
		      int argc, char *argv[], int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(argc, argv), cname("qio_writePropagator"), source_hypercube(0)
    { setHeader(ensemble_id, ensemble_label, traj, propagator_type, source_type); 
    write_ScS_12sink(outfile, prop, scSource, volFormat, floatFormat);}

  
  
  //! write propagator, format 12 pairs of source and sink
  /*!
    \param outfile file to write to
    \param sType source type (scalar, full)
    \param prop    pointer to sinks
    \param source pointer to scalar source
   */
  qio_writePropagator(char *outfile, const QIO_PROP_SOURCE_TYPES sType, const void *prop, const void *source, 
		      int argc, char *argv[], int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(argc, argv), cname("qio_writePropagator"), source_hypercube(0)
    { initHeader(); write_12pairs(outfile, sType, prop, source, volFormat, floatFormat);}

  //! write propagator, format 12 pairs of source and sink, set add. info
  /*!
    \param outfile file to write to
    \param sType source type (scalar, full)
    \param prop    pointer to sinks
    \param source pointer to scalar source
    \param ensemble_id ensemble id
    \param ensemble_label ensemble label
    \param traj trajectory number
   */
  qio_writePropagator(char *outfile, const QIO_PROP_SOURCE_TYPES sType, const void *prop, const void *source, 
		      const char * ensemble_id, const char * ensemble_label, const int traj, 
		      int argc, char *argv[], int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(argc, argv), cname("qio_writePropagator"), source_hypercube(0)
    { setHeader(ensemble_id, ensemble_label, traj); write_12pairs(outfile, sType, prop, source, volFormat, floatFormat);}
  

  //! write propagator, format 12 pairs of source and sink, set add. info
  /*!
    \param outfile file to write to
    \param sType source type (scalar, full)
    \param prop    pointer to sinks
    \param source pointer to scalar source
    \param ensemble_id ensemble id
    \param ensemble_label ensemble label
    \param traj trajectory number
    \param propagator_type descr. of prop. type
    \param source_type descr. of source type
   */
  qio_writePropagator(char *outfile, const QIO_PROP_SOURCE_TYPES sType, const void *prop, const void *source, 
		      const char * ensemble_id, const char * ensemble_label, const int traj,
		      const char * propagator_type, const char * source_type,
		      int argc, char *argv[], int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(argc, argv), cname("qio_writePropagator"), source_hypercube(0)
    { setHeader(ensemble_id, ensemble_label, traj, propagator_type, source_type); 
    write_12pairs(outfile, sType, prop, source, volFormat, floatFormat);}
  
  // version w/o argc/argv, using GJP

  qio_writePropagator(): qio_init(GJP.argc(), GJP.argv()), cname("qio_writePropagator"), source_hypercube(0)
    {initHeader();}
  
  
  //! write propagator, format scalar source + 12 sinks
  /*!
    \param outfile file to write to
    \param prop    pointer to sinks
    \param scSource pointer to scalar source
  */
  qio_writePropagator(char *outfile, const void *prop, const void *scSource, 
		      int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(GJP.argc(), GJP.argv()), cname("qio_writePropagator"), source_hypercube(0)
    { initHeader(); write_ScS_12sink(outfile, prop, scSource, volFormat, floatFormat);}

   
  //! write propagator, format scalar source + 12 sinks, set add. info
  /*!
    \param outfile file to write to
    \param prop    pointer to sinks
    \param scSource pointer to scalar source
    \param ensemble_id ensemble id
    \param ensemble_label ensemble label
    \param traj trajectory number
   */
  qio_writePropagator(char *outfile, const void *prop, const void *scSource, 
		      const char * ensemble_id, const char * ensemble_label, const int traj,
		      int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(GJP.argc(), GJP.argv()), cname("qio_writePropagator"), source_hypercube(0)
    { setHeader(ensemble_id, ensemble_label, traj); write_ScS_12sink(outfile, prop, scSource, volFormat, floatFormat);}

  //! write propagator, format scalar source + 12 sinks, set add. info
  /*!
    \param outfile file to write to
    \param prop    pointer to sinks
    \param scSource pointer to scalar source
    \param ensemble_id ensemble id
    \param ensemble_label ensemble label
    \param traj trajectory number
    \param propagator_type descr. of prop. type
    \param source_type descr. of source type
   */
  qio_writePropagator(char *outfile, const void *prop, const void *scSource, 
		      const char * ensemble_id, const char * ensemble_label, const int traj,
		      const char * propagator_type, const char * source_type,
		      int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(GJP.argc(), GJP.argv()), cname("qio_writePropagator"), source_hypercube(0)
    { setHeader(ensemble_id, ensemble_label, traj, propagator_type, source_type); 
    write_ScS_12sink(outfile, prop, scSource, volFormat, floatFormat);}

  
  
  //! write propagator, format 12 pairs of source and sink
  /*!
    \param outfile file to write to
    \param sType source type (scalar, full)
    \param prop    pointer to sinks
    \param source pointer to scalar source
   */
  qio_writePropagator(char *outfile, const QIO_PROP_SOURCE_TYPES sType, const void *prop, const void *source, 
		      int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(GJP.argc(), GJP.argv()), cname("qio_writePropagator"), source_hypercube(0)
    { initHeader(); write_12pairs(outfile, sType, prop, source, volFormat, floatFormat);}

  //! write propagator, format 12 pairs of source and sink, set add. info
  /*!
    \param outfile file to write to
    \param sType source type (scalar, full)
    \param prop    pointer to sinks
    \param source pointer to scalar source
    \param ensemble_id ensemble id
    \param ensemble_label ensemble label
    \param traj trajectory number
   */
  qio_writePropagator(char *outfile, const QIO_PROP_SOURCE_TYPES sType, const void *prop, const void *source, 
		      const char * ensemble_id, const char * ensemble_label, const int traj, 
		      int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(GJP.argc(), GJP.argv()), cname("qio_writePropagator"), source_hypercube(0)
    { setHeader(ensemble_id, ensemble_label, traj); write_12pairs(outfile, sType, prop, source, volFormat, floatFormat);}
  
  //! write propagator, format 12 pairs of source and sink, set add. info
  /*!
    \param outfile file to write to
    \param sType source type (scalar, full)
    \param prop    pointer to sinks
    \param source pointer to scalar source
    \param ensemble_id ensemble id
    \param ensemble_label ensemble label
    \param traj trajectory number
    \param propagator_type descr. of prop. type
    \param source_type descr. of source type
   */
  qio_writePropagator(char *outfile, const QIO_PROP_SOURCE_TYPES sType, const void *prop, const void *source, 
		      const char * ensemble_id, const char * ensemble_label, const int traj,
		      const char * propagator_type, const char * source_type,
		      int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(GJP.argc(), GJP.argv()), cname("qio_writePropagator"), source_hypercube(0)
    { setHeader(ensemble_id, ensemble_label, traj, propagator_type, source_type); 
    write_12pairs(outfile, sType, prop, source, volFormat, floatFormat);}
  

  

  virtual ~qio_writePropagator(){ 
    #ifdef DEBUG_Init
    printf("finished: qio_writePropagator\n");
    #endif // DEBUG_Init
  }


  //! write propagator, format scalar source + 12 sinks
  /*!
    \param outfile file to write to
    \param prop    pointer to sinks
    \param scSource pointer to scalar source
  */  
  void write_ScS_12sink(char *outfile, const void *prop, const void *scSource, 
			int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC);

  //! write propagator, format 12 pairs of source and sink
  /*!
    \param outfile file to write to
    \param sType source type (scalar, full)
    \param prop    pointer to sinks
    \param source pointer to scalar source
   */  
  void write_12pairs(char *outfile, const QIO_PROP_SOURCE_TYPES sType, const void *prop, const void *source, 
		     int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC);
  void write_pair(char *outfile, const QIO_PROP_SOURCE_TYPES sType, const void *prop, const void *source,
		  int spin,
		  int color,
		  int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC);

  //! set additional info
  /*!
    \param ensemble_id ensemble id
    \param ensemble_label ensemble label
    \param traj trajectory number
  */
  void setHeader(const char * ensemble_id, const char * ensemble_label, const int traj);

  //! set additional info
  /*!
    \param propagator_type descr. of prop. type
    \param source_type descr. of source type
   */
  void setHeader(const char * propagator_type, const char * source_type);

   //! set additional info
  /*!
    \param ensemble_id ensemble id
    \param ensemble_label ensemble label
    \param traj trajectory number
    \param propagator_type descr. of prop. type
    \param source_type descr. of source type
  */
  void setHeader(const char * ensemble_id, const char * ensemble_label, const int traj, 
		 const char * propagator_type, const char * source_type)
    {
      setHeader(ensemble_id, ensemble_label, traj);
      setHeader(propagator_type, source_type);
    }

  //! set tslice (only writing this one and only tslice for the source) 
  void setSourceTslice( const int tslice);

  //! set tslices (only writing from tslices from t_start to t_end for the source)
  void setSourceTslices( const int t_start, const int t_end);

  //! only write a hypercube with start-, end-coordinates for the source
  void setSourceHypercube( const int start[4], const int end[4]);

  //! write a complete source
  void setSourceComplete(void)
    { source_hypercube=0;}

  //! writing a complete (0) or hypercube/tslices source (1)
  int SourceHypercube(void)
    { return source_hypercube;}

 
 private:

  QIO_Writer *qio_Output;

  void qio_openOutput(char *filename, QIO_String *record_file, int volFormat);

  void qio_closeOutput()
    { QIO_close_write(qio_Output);}


  void initHeader()
    { 
      header_traj = -1; 
      strcpy(header_ensemble_id, "not specified" );
      strcpy(header_ensemble_label, "not specified");
      strcpy(header_propagator_type, "not specified");
      strcpy(header_source_type, "not specified");
    }

  int header_traj;
  char header_ensemble_id[MAX_HEADER_LINE];
  char header_ensemble_label[MAX_HEADER_LINE];
  char header_propagator_type[MAX_HEADER_LINE];
  char header_source_type[MAX_HEADER_LINE];

  int source_hypercube; 
  int source_start[QIO_RW_DIMENSION];
  int source_end[QIO_RW_DIMENSION];
  


};



CPS_END_NAMESPACE
#endif // __QIOWRITEPROP__
