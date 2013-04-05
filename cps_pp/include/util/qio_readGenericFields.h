#ifndef __QIO_READENERICFIELDS_H__
#define __QIO_READGENERICFIELDS_H__

#include <util/qio_general.h>
#include <util/gjp.h>
#include <util/fpconv.h>

CPS_START_NAMESPACE
using namespace std;

// QIO for Generic Fields
class qio_readGenericFields: private qio_init {

 private:

  char *cname;

public:

  // version with argc/argv

  qio_readGenericFields( int argc, char *argv[]):
    qio_init(argc, argv),
    cname( (char*) "qio_readGenericFields" )
    { }

  qio_readGenericFields(): qio_init(GJP.argc(), GJP.argv()), cname((char*)"qio_readGenericFields")
  {}
  
  //! read field
  /*!
    \param outfile file to read to
    \param field    pointer to field
   */
  qio_readGenericFields(char *outfile,
			 const int n_fields, const int f_size_per_site,  void *field,
			 int argc, char *argv[], int volFormat=QIO_VOLFMT,
			 FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(argc, argv),
    cname((char*) "qio_readGenericFields")
    {
      read_genericfields(outfile, n_fields, f_size_per_site, field,  volFormat, floatFormat);
    }


  //! read generic field, plus  infos
  /*!
    \param outfile file to read to
    \param field    pointer to field
    \param ensemble_id ensemble id
    \param ensemble_label ensemble label
    \param traj trajectory number
   */
  qio_readGenericFields(char *outfile,
			 const int n_fields, const int f_size_per_site,  void *field,
			 const char * ensemble_id,
			 const char * ensemble_label,
			 const int traj,
         		 const char * field_type_label,
			 int argc, char *argv[],
			 int volFormat=QIO_VOLFMT,
			 FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(argc, argv), cname((char*) "qio_readGenericFields")
  {
    read_genericfields(outfile, n_fields, f_size_per_site, field,  volFormat, floatFormat);
  }

  virtual ~qio_readGenericFields(){ 
    #ifdef DEBUG_Init
    printf("finished: qio_readGenericFields\n");
    #endif // DEBUG_Init
  }


  //! read field
  /*!
    \param infile file to read to
     \param f_size_per_site  degree on a site  in unit of Float.
     \param n_fields  number of fields
     \param field   pointer to field

     The field should store data in memory in the following format :

        [ 1 st field ]  [ 2 nd field ] ...   [ (n_fields-1q)-th field ]

    where  [ n-th field ]  is
    [ f_size_per_site  Floats for (0,0,0,0) ]  [ f_size_per_site  Floats for (1,0,0,0) ]  [ f_size_per_site  Floats for (2,0,0,0) ] ....    [f_size_per_site Floats for (Nx-1, Ny-1, Nz-1, Nt-1) ]

    To save the number of io, we rearrange the file format as follows :
    
    [ n_fields* f_size_per_site  Floats for (0,0,0,0) ]  [ n_fields* f_size_per_site  Floats for (1,0,0,0) ]  [ n_fields* f_size_per_site  Floats for (2,0,0,0) ] ....    [n_fields* f_size_per_site Floats for (Nx-1, Ny-1, Nz-1, Nt-1) ]


    the most fastest changing index is the f_size_per_site degree in one field, then the index for the field, 0 ... n_field-1

    This rearrangement requires the non-local memory access, but I hope the benefit of n_field times smaller number of I/O will supersede the slow down.
    
    CONFESSION:  I didn't fully understand all of the QIO's sophisticated specifications in
    http://usqcd.jlab.org/usqcd-docs/qio/qio_2p3.pdf

  */  

  void read_genericfields(char *infile,
			    const int n_fields, const int f_size_per_site,  void *field,
			    int volFormat=QIO_VOLFMT,
			    FP_FORMAT floatFormat=FP_AUTOMATIC);

  //! read field lives on one of even odd sites, see qio_readgenericfields(...) for arguments and formats.
  /*!
  \param outfile   file to read to
  \param f_size_per_site  degree on a site  in unit of Float.  e.g. for even/odd wilson fermion  f_size_per_site is still 24 not 12
  */
  void read_genericfields_eo(
			      char *infile,
			      const int n_fields, int f_size_per_site,  void *field,
			      int volFormat=QIO_VOLFMT ,
			      FP_FORMAT floatFormat=FP_AUTOMATIC )
  {
    const char * fname = "read_genericfields_eo(...)";
    if(GJP.XnodeSites() %2 != 0 )
      ERR.NotImplemented(cname, fname, "X-direction needs to be even length");
    
    read_genericfields(infile,  n_fields, f_size_per_site/2, field,
			volFormat,  floatFormat);
  }
  
private:

  int header_traj;
  char header_ensemble_id[MAX_HEADER_LINE];
  char header_ensemble_label[MAX_HEADER_LINE];
  char header_field_type_label[MAX_HEADER_LINE];

  QIO_Reader *qio_Input;

  void qio_openInput(char *filename, QIO_String *record_file, int volFormat);

  void qio_closeInput()
  { QIO_close_read(qio_Input); }



};


CPS_END_NAMESPACE
#endif // __QIO_READGENERICFIELDS_H__
