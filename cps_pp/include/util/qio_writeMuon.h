#ifndef __QIOWRITEMUON__
#define __QIOWRITEMUON__

#include <util/fpconv.h>


using namespace std;

#ifndef USE_QIO
CPS_START_NAMESPACE
class qio_writeMuon {

 private:

  char *cname;

 public:

 qio_writeMuon( int argc, char *argv[]): cname("qio_writeMuon") { }
  
  virtual ~qio_writeMuon() { }
  
  void write(char *outfile, Lattice &lat, Rcomplex *muon, int count,
	      char *ildgLFN="added ildgLFN", 
	      int volFormat=0, 
	      FP_FORMAT floatFormat=FP_AUTOMATIC){
      ERR.NotImplemented(cname,"wrte()");
	}
  
  void setHeader(const char * ensemble_id, const char * ensemble_label, const int traj);
  
  
 private:
  

  void qio_openOutput(char *filename, char *stringLFN, char *xml_write_file, int volFormat){
      ERR.NotImplemented(cname,"qio_openOutput()");
	}

  void qio_closeOutput(){
      ERR.NotImplemented(cname,"qio_closeOutput()");
	}

  void initHeader() { }



};
CPS_END_NAMESPACE
#else


#include <util/qio_general.h>
CPS_START_NAMESPACE
class qio_writeMuon: private qio_init {

 private:

  char *cname;

 public:

 qio_writeMuon( int argc, char *argv[]): qio_init(argc, argv), cname("qio_writeMuon")
    {
      initHeader();
    }
  
  virtual ~qio_writeMuon()
    { 
#ifdef DEBUG_Init
      printf("finished: qio_writeMuon\n");
#endif // DEBUG_Init
    }
  
  void write(char *outfile, Lattice &lat, Rcomplex *muon, int count,
	      char *ildgLFN="added ildgLFN", 
	      int volFormat=QIO_VOLFMT, 
	      FP_FORMAT floatFormat=FP_AUTOMATIC);
  
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
#endif



#endif // __QIOWRITEMUON__
