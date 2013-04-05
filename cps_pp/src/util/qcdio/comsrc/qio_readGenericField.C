#ifdef USE_QIO
#include <config.h>
#include <util/qio_readGenericFields.h>
#include <util/time_cps.h>



CPS_START_NAMESPACE
using namespace std;

#define PROFILE

// qio-factory functions 

struct qio_genfield_glb_type {  int precision; int n_fields; int f_size_per_site; int n_sites; };
static qio_genfield_glb_type qio_genfield_glb;

void qio_genfield_put_glb(char *buf_, size_t site_index, int count, void *arg_)
{
  moveMem( arg_, buf_, sizeof(qio_genfield_glb_type));
}

// global variable for fanctor  (could make a separated record and  read from it, but then we would need if statement...)
static int qio_glb_genfield_n_fields, qio_glb_genfield_f_size_per_site, qio_glb_n_sites;
void qio_putGenField(char *buf_, size_t site_index, int count, void *arg)
{
 /*printf(" called  with count %i\n",count);*/\
  
  Float *field = (Float*) arg;

  const int n_field = qio_genfield_glb. n_fields;
  const int f_size = qio_genfield_glb. f_size_per_site;
  const int n_sites = qio_genfield_glb. n_sites;

  Float *buf = (Float*) buf_;
  
  /*
    The field should store data in memory in the following format :

    [ 1 st field ]  [ 2 nd field ] ...   [ (n_fields-1q)-th field ]

    where  [ n-th field ]  is
    [ f_size_per_site  Floats for (0,0,0,0) ]  [ f_size_per_site  Floats for (1,0,0,0) ]  [ f_size_per_site  Floats for (2,0,0,0) ] ....    [f_size_per_site Floats for (Nx-1, Ny-1, Nz-1, Nt-1) ]

    To save the number of io, we rearrange the file format as follows :
    
    [ n_fields* f_size_per_site  Floats for (0,0,0,0) ]  [ n_fields* f_size_per_site  Floats for (1,0,0,0) ]  [ n_fields* f_size_per_site  Floats for (2,0,0,0) ] ....    [n_fields* f_size_per_site Floats for (Nx-1, Ny-1, Nz-1, Nt-1) ]


    the most fastest changing index is the f_size_per_site degree in one field, then the index for the field, 0 ... n_field-1

    This rearrangement requires the non-local memory access, but I hope the benefit of n_field times smaller number of I/O will supersede the slow down.
    
  */

  for(int field_i=0; field_i<n_field;++field_i) {
 
    moveMem(  field+ f_size* (site_index + n_sites * field_i   ), 
	      buf  + f_size* ( field_i ),  
	      sizeof(Float)* f_size);
  }
  
  /*printf("  %i wrote %f %f\n",ii,*(array + 2*ii), *(array+2*ii+1));*/\  
}




// now start class-functions...

void qio_readGenericFields::qio_openInput(char *filename, QIO_String *record_file, int volFormat)
{

  char * fname = "qio_openInput(...)";
  VRB.Func(cname,fname);

  //#define DEBUG_openInput
#ifdef DEBUG_openInput
  printf("called qio_openInput with filename: %s volform=%d\n",filename,volFormat);
#endif // DEBUG_openInput

  
  // QIO_PARALLEL ?
  const int serpar(QIO_SERPAR);
  //const int serpar(QIO_PARALLEL);

  QIO_Reader *infile;
  QIO_Iflag iflag;
  iflag.serpar = serpar;
  iflag.volfmt = volFormat;
  
  VRB.Flow(cname,fname,"open input file %s\n",filename);
  infile = QIO_open_read(record_file, filename, &layout, pointer_fs, &iflag);


#ifdef DEBUG_openInput
  printf("done QIO_open_read\n");
  #endif // DEBUG_openInput

  #ifdef DEBUG_openInput
  printf("finishing qio_openInput, filename: %s\n",filename);
  #endif // DEBUG_openInput

  qio_Input = infile;

  VRB.FuncEnd(cname,fname);

}


void qio_readGenericFields::read_genericfields(
						 char *infile,
						 int n_fields, int f_size_per_site, void *field,
						 int volFormat, FP_FORMAT floatFormat)
 {

#ifdef PROFILE
   time_elapse();
#endif

   const char * fname = "read_genericfields(...)";
   VRB.Func(cname,fname);
   VRB.Result(cname,fname,"reading generic fields     from %s \n",infile);
   
   int return_val(0);

   qio_setLayout();
   qio_setFilesystem();
   
#if 0
   //detect output format

  int SingleDouble(1);
  // 1=double (standard), else single
  
  switch( floatFormat )
    {
    case FP_TIDSP32 :      
    case FP_IEEE32 :       
    case FP_IEEE32BIG :    
    case FP_IEEE32LITTLE : SingleDouble=0; break;
    default: SingleDouble=1;
    }

  if (SingleDouble) VRB.Flow(cname,fname," input-precision: DOUBLE\n");
  else VRB.Flow(cname,fname," output-precision: SINGLE\n");

#endif
  
  
  char xml_info_field[7*(MAX_HEADER_LINE+10)];

  QIO_String *qio_file_header = QIO_string_create();


  qio_openInput(infile, qio_file_header, volFormat);

  //  Fist , store the n_fields and f_size_per_site information as global data

  QIO_String *no_string = QIO_string_create();
  // QIO_string_set( no_string, "");

  QIO_RecordInfo* glb_rec_info
    = QIO_create_record_info(QIO_GLOBAL, NULL, NULL, 0, (char*)"", (char*)"I", 0, 0, 0, 0 );

  return_val += 
    QIO_read( qio_Input, glb_rec_info, no_string, qio_genfield_put_glb,
	      sizeof(qio_genfield_glb_type), sizeof(int),
	      (void*) &qio_genfield_glb);

  n_fields = qio_genfield_glb.n_fields;
  f_size_per_site = qio_genfield_glb. f_size_per_site;



  // In the file, the number of global sites are saved
  if ( qio_genfield_glb. n_sites != GJP. VolSites() )
    ERR.General(cname,fname,"read VolSites is wrong %d vs %d", qio_genfield_glb. n_sites, GJP.VolSites() );

  VRB.Flow(cname,fname, "header n_fields = %d, f_size_per_site = %d", n_fields, f_size_per_site );
  

  // For the actual I/O, we will use the number of local sites
  qio_genfield_glb. n_sites =  GJP.VolNodeSites();
  
  if( qio_genfield_glb. precision != 1 )
    ERR.NotImplemented(cname,fname, "Only double precision is supported for I/O");

  int SingleDouble(1);
    // 1=double (standard), else single

  // Now create the record info for fireld itself
  
  QIO_RecordInfo *record_field 
    =   QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, (char*)"", (char*)"", 0, 0, 0, 0);

  
  // Do the actual work
  if(SingleDouble)
    {
      //output in double-precision

      return_val += QIO_read( qio_Input, record_field, no_string,
			       qio_putGenField,
			       n_fields* f_size_per_site*sizeof(Float),
			       sizeof(Float), field);

    }
  else
    {
      ERR.NotImplemented(cname,fname,"Sorry for lack of support for single precision sotre at this moment");
    }


  if ( (return_val == 0) ) 
    VRB.Result(cname,fname,"QIO_read successfull...\n");
  else
    ERR.General(cname,fname,"ERROR QIO: QIO_read(s) returned %i\n",return_val);


  // clean-up
  QIO_destroy_record_info( glb_rec_info );
  QIO_destroy_record_info(record_field);
  QIO_string_destroy(qio_file_header);
  QIO_string_destroy(no_string);

  qio_closeInput();

#ifdef PROFILE
  if(!UniqueID()) printf("%s %e sec\n",cname, time_elapse());
#endif


  VRB.FuncEnd(cname,fname);

}


CPS_END_NAMESPACE
#endif
