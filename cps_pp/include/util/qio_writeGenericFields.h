#ifndef __QIO_WRITEGENERICFIELDS_H__
#define __QIO_WRITEGENERICFIELDS_H__

#include <util/qio_general.h>
#include <util/gjp.h>
#include <util/fpconv.h>

CPS_START_NAMESPACE
using namespace std;

// QIO for Generic Fields
class qio_writeGenericFields: private qio_init {

 private:

  char *cname;

public:

  // version with argc/argv

  qio_writeGenericFields( int argc, char *argv[]):
    qio_init(argc, argv),
    cname( (char*) "qio_writeGenericFields" )
    {
      // if( GJP. XnodeSites() %2 !=0 )
      // ERR.NotImplemented(cname, "qio_writeGenericFields(i,c)",
      // "# of sites in x-direction needs to be even");
      
      initHeader();
    }

  qio_writeGenericFields(): qio_init(GJP.argc(), GJP.argv()), cname("qio_writeGenericFields")
  {initHeader();}
  
  //! write field
  /*!
    \param outfile file to write to
    \param field    pointer to field
   */
  qio_writeGenericFields(char *outfile,
			 const int n_fields, const int f_size_per_site,  void *field,
			 int argc, char *argv[], int volFormat=QIO_VOLFMT,
			 FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(argc, argv),
    cname((char*) "qio_writeGenericFields")
    {
      //      if( GJP. XnodeSites() %2 !=0 )
      //ERR.NotImplemented(cname, "qio_writeGenericFields(i,c)",
      // "# of sites in x-direction needs to be even");
      initHeader();
      write_genericfields(outfile, n_fields, f_size_per_site, field,  volFormat, floatFormat);
    }


  //! write generic field, plus  infos
  /*!
    \param outfile file to write to
    \param field    pointer to field
    \param ensemble_id ensemble id
    \param ensemble_label ensemble label
    \param traj trajectory number
   */
  qio_writeGenericFields(char *outfile,
			 const int n_fields, const int f_size_per_site,  void *field,
			 const char * ensemble_id,
			 const char * ensemble_label,
			 const int traj,
         		 const char * field_type_label,
			 int argc, char *argv[],
			 int volFormat=QIO_VOLFMT,
			 FP_FORMAT floatFormat=FP_AUTOMATIC):
    qio_init(argc, argv), cname((char*) "qio_writeGenericFields")
  { setHeader(ensemble_id, ensemble_label, traj, field_type_label);
    write_genericfields(outfile, n_fields, f_size_per_site, field,  volFormat, floatFormat);
  }

  virtual ~qio_writeGenericFields(){ 
    #ifdef DEBUG_Init
    printf("finished: qio_writeGenericFields\n");
    #endif // DEBUG_Init
  }


  //! write field
  /*!
    \param outfile file to write to
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

  void write_genericfields(char *outfile,
			    const int n_fields, const int f_size_per_site,  void *field,
			    int volFormat=QIO_VOLFMT,
			    FP_FORMAT floatFormat=FP_AUTOMATIC);

  //! write field lives on one of even odd sites, see qio_writegenericfields(...) for arguments and formats.
  /*!
  \param outfile   file to write to
  \param f_size_per_site  degree on a site  in unit of Float.  e.g. for even/odd wilson fermion  f_size_per_site is still 24 not 12
  */
  void write_genericfields_eo(
			      char *outfile,
			      const int n_fields, int f_size_per_site,  void *field,
			      int volFormat=QIO_VOLFMT ,
			      FP_FORMAT floatFormat=FP_AUTOMATIC )
  {
    const char * fname = "write_genericfields_eo(...)";
    if(GJP.XnodeSites() %2 != 0 )
      ERR.NotImplemented(cname, fname, "X-direction needs to be even length");
    
    write_genericfields(outfile,  n_fields, f_size_per_site/2, field,
			volFormat,  floatFormat);
  }
  
  //! set additional info
  /*!
    \param ensemble_id ensemble id
    \param ensemble_label ensemble label
    \param traj trajectory number
    \param field_type descr. of field. type
  */
  void setHeader(const char * ensemble_id, const char * ensemble_label, const int traj, 
		 const char * field_type_label);

private:

  int header_traj;
  char header_ensemble_id[MAX_HEADER_LINE];
  char header_ensemble_label[MAX_HEADER_LINE];
  char header_field_type_label[MAX_HEADER_LINE];

  QIO_Writer *qio_Output;

  void qio_openOutput(char *filename, QIO_String *record_file, int volFormat);

  void qio_closeOutput()
    { QIO_close_write(qio_Output);}


  void initHeader()
    { 
      header_traj = -1; 
      strcpy(header_ensemble_id, "not specified" );
      strcpy(header_ensemble_label, "not specified");
      strcpy(header_field_type_label, "not specified");
    }



};


CPS_END_NAMESPACE
#endif // __QIO_WRITEGENERICFIELDS_H__
