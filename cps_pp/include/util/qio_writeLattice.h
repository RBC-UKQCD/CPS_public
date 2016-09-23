#ifndef __QIOWRITELAT__
#define __QIOWRITELAT__

#include <util/qio_general.h>
#include <util/fpconv.h>


CPS_START_NAMESPACE

class qio_writeLattice: private qio_init {

 private:

  char *cname;

 public:

  // versions with argc/argv

  qio_writeLattice( int argc, char *argv[]): qio_init(argc, argv), cname("qio_writeLattice")
    {initHeader();}

  //! write a lattice configuration using QIO
  /*!
    \param outfile file to write to
    \param lat pointer to lattice
   */
  qio_writeLattice(char *outfile, Lattice &lat, int argc, char *argv[], int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
   qio_init(argc, argv), cname("qio_writeLattice")
    { initHeader(); write(outfile, "", lat, volFormat, floatFormat);}

  //! write a lattice configuration using QIO, set add. info
  /*!
    \param outfile file to write to
    \param lat pointer to lattice
    \param ildgLFN ILDG logical file name
   */
  qio_writeLattice(char *outfile, Lattice &lat, char *ildgLFN, int argc, char *argv[], int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
   qio_init(argc, argv), cname("qio_writeLattice")
    { initHeader(); write(outfile, ildgLFN, lat, volFormat, floatFormat);}
 
  //! write a lattice configuration using QIO, set add. info
  /*!
    \param outfile file to write to
    \param lat pointer to lattice
    \param ildgLFN ILDG logical file name
    \param ensemble_id ensemble ID
    \param ensemble_label ensemble label
    \param traj trajectory number
   */
  qio_writeLattice(char *outfile, Lattice &lat, char *ildgLFN, const char * ensemble_id, const char * ensemble_label, const int traj,
		   int argc, char *argv[], int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(argc, argv), cname("qio_writeLattice")
    { setHeader(ensemble_id, ensemble_label, traj); write(outfile, ildgLFN, lat, volFormat, floatFormat);}
  
  // versions w/o argc/argv, using GJP
  
  //! write a lattice configuration using QIO, set add. info
  /*!
    \param outfile file to write to
    \param lat pointer to lattice
    \param ildgLFN ILDG logical file name
    \param ensemble_id ensemble ID
    \param ensemble_label ensemble label
    \param traj trajectory number
   */
  qio_writeLattice(): qio_init(GJP.argc(), GJP.argv()), cname("qio_writeLattice")
    {initHeader();}

  //! write a lattice configuration using QIO
  /*!
    \param outfile file to write to
    \param lat pointer to lattice
  */
  qio_writeLattice(char *outfile, Lattice &lat, int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
   qio_init(GJP.argc(), GJP.argv()), cname("qio_writeLattice")
    { initHeader(); write(outfile, "", lat, volFormat, floatFormat);}

  //! write a lattice configuration using QIO, set add. info
  /*!
    \param outfile file to write to
    \param lat pointer to lattice
    \param ildgLFN ILDG logical file name
  */
  qio_writeLattice(char *outfile, Lattice &lat, char *ildgLFN, int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
   qio_init(GJP.argc(), GJP.argv()), cname("qio_writeLattice")
    { initHeader(); write(outfile, ildgLFN, lat, volFormat, floatFormat);}
 
  //! write a lattice configuration using QIO, set add. info
  /*!
    \param outfile file to write to
    \param lat pointer to lattice
    \param ildgLFN ILDG logical file name
    \param ensemble_id ensemble ID
    \param ensemble_label ensemble label
    \param traj trajectory number
   */
  qio_writeLattice(char *outfile, Lattice &lat, char *ildgLFN, const char * ensemble_id, const char * ensemble_label, const int traj,
		   int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC):
   qio_init(GJP.argc(), GJP.argv()), cname("qio_writeLattice")
    { setHeader(ensemble_id, ensemble_label, traj); write(outfile, ildgLFN, lat, volFormat, floatFormat);}
  


  virtual ~qio_writeLattice(){ 
    #ifdef DEBUG_Init
      std::printf("finished: qio_writeLattice\n");
    #endif // DEBUG_Init
  }

  //! write a lattice configuration using QIO
  /*!
    \param outfile file to write to
    \param lat pointer to lattice
    \param ildgLFN ILDG logical file name
  */
  void write(char *outfile, char *ildgLFN, Lattice &lat, int volFormat=QIO_VOLFMT, FP_FORMAT floatFormat=FP_AUTOMATIC);

  //! set additional info
  /*!
    \param ensemble_id ensemble ID
    \param ensemble_label ensemble label
    \param traj trajectory number
   */
  void setHeader(const char * ensemble_id, const char * ensemble_label, const int traj);

 
 private:

  QIO_Writer *qio_Output;

  void qio_openOutput(char *filename, char *stringLFN, char *xml_write_file, int volFormat);

  void qio_closeOutput()
    { QIO_close_write(qio_Output);}


  void initHeader()
    { 
      header_traj = -1; 
      strcpy(header_ensemble_id, "not specified" );
      strcpy(header_ensemble_label, "not specified");
    }

  int header_traj;
  char header_ensemble_id[MAX_HEADER_LINE];
  char header_ensemble_label[MAX_HEADER_LINE];
    


};



CPS_END_NAMESPACE
#endif // __QIOWRITELAT__
