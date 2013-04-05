#ifdef USE_QIO
#include <config.h>
#include <util/qio_writeGenericFields.h>
#include <util/time_cps.h>


#define PROFILE

CPS_START_NAMESPACE
using namespace std;


// qio-factory functions 

// global variable for fanctor  (could make a separated record and  read from it, but then we would need if statement...)
struct qio_genfield_glb_type {  int precision; int n_fields; int f_size_per_site; int n_sites; };
static qio_genfield_glb_type qio_genfield_glb;

void qio_genfield_get_glb(char *buf_, size_t site_index, int count, void *arg_)
{
  moveMem(buf_, arg_, sizeof(qio_genfield_glb_type)) ;

  qio_genfield_glb_type* a;

#if 0
  printf("qio_genfield_get_glb: writing %d %d %d %d\n",
	 a-> precision,
	 a-> n_fields,
	 a-> f_size_per_site,
	 a-> n_sites );
#endif
  
    //  int *arg=(int*)arg_;   int *buf=(int*)buf_;
    //  buf[0]=arg[0]; buf[1]=arg[1];
}

void qio_getGenField(char *buf_, size_t site_index, int count, void *arg)
{
 /*printf(" called  with count %i\n",count);*/\
  
  const Float *field = (Float*) arg;

  const int n_field = qio_genfield_glb. n_fields;
  const int f_size = qio_genfield_glb. f_size_per_site;
  const int n_sites = qio_genfield_glb. n_sites;

#if 0
  printf("qio_geGenField: %d %d %d %d %d %x\n", n_field, f_size, n_sites,count, site_index, buf_);
#endif
  
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

  for(int field_i=0; field_i<n_field;++field_i){
    moveMem(  buf  + f_size* ( field_i   ),  
	      field+ f_size* (site_index + n_sites * field_i   ), 
	      sizeof(Float)* f_size);
  }
  //printf("\n");
  
  /*printf("  %i wrote %f %f\n",ii,*(array + 2*ii), *(array+2*ii+1));*/\  
}


// now start class-functions...

void qio_writeGenericFields::qio_openOutput(char *filename, QIO_String *record_file, int volFormat)
{

  const char * fname = "qio_openOutput(...)";
  VRB.Func(cname,fname);

  const int volfmt(volFormat);
  const int serpar(QIO_SERPAR);
  const int ildgstyle = QIO_ILDGNO;
  
  VRB.Result(cname,fname,"volfmt=%d serpar=%d",volfmt,serpar);

  QIO_Writer *outfile;
  QIO_Oflag oflag;


  VRB.Flow(cname,fname,"open output file %s \n",filename);

  oflag.serpar=serpar;
  oflag.ildgstyle=ildgstyle;

  oflag.ildgLFN = NULL;

  oflag.mode = QIO_TRUNC;

  outfile = QIO_open_write(record_file, filename, volfmt, &layout, pointer_fs, &oflag);

  qio_Output = outfile;

  VRB.FuncEnd(cname,fname);

}


void qio_writeGenericFields::write_genericfields(
						 char *outfile,
						 int n_fields, int f_size_per_site, void *field,
						 int volFormat, FP_FORMAT floatFormat)
 {
#ifdef PROFILE
   time_elapse();
#endif

   const char * fname = "write_genericfields(...)";
   VRB.Func(cname,fname);
   VRB.Result(cname,fname,"writing generic fields   to %s \n",outfile);
   
   int return_val(0);

   qio_setLayout();
   qio_setFilesystem();

   
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

  qio_genfield_glb.precision=SingleDouble;
  if (SingleDouble)  VRB.Flow(cname,fname," output-precision: DOUBLE\n");
  else VRB.Flow(cname,fname," output-precision: SINGLE\n");

  char xml_info_field[7*(MAX_HEADER_LINE+10)];

  sprintf(xml_info_field, "DATATYPE = 4D_GENERICFIELD(n_fields = %i, f_size_per_site = %i\nFIELDTYPE = %s\nENSEMBLE_ID = %s\nENSEMBLE_LABEL = %s\nSEQUENCE_NUMBER = %i", n_fields, f_size_per_site,
 	  header_field_type_label, header_ensemble_id, header_ensemble_label, header_traj);

  //check for length...
  if( (strlen(xml_info_field)*sizeof(char)) > QIO_INFO_STRING_MAX){
    ERR.General(cname,fname," xml_info_file too large: %i\n %s\n", (strlen(xml_info_field)*sizeof(char)), xml_info_field);
    exit(-12);
  }
  QIO_String *qio_file_header = QIO_string_create();
  QIO_string_set( qio_file_header, xml_info_field );

  qio_openOutput(outfile, qio_file_header, volFormat);

  //  Fist , store the n_fields and f_size_per_site information as global data

  QIO_String *no_string = QIO_string_create();
  QIO_string_set( no_string, "");

  QIO_RecordInfo* glb_rec_info
    = QIO_create_record_info(QIO_GLOBAL, NULL, NULL, 0, (char*)"", (char*)"I", 0, 0, sizeof(qio_genfield_glb_type), 1);

  qio_genfield_glb.n_fields = n_fields;
  qio_genfield_glb. f_size_per_site = f_size_per_site;

  // To the file, we record the number of global sites
  qio_genfield_glb. n_sites = GJP. VolSites();
  
  return_val += 
    QIO_write( qio_Output, glb_rec_info, no_string, qio_genfield_get_glb, sizeof(qio_genfield_glb_type), sizeof(int), (void*) &qio_genfield_glb);
  
  // For the actual I/O, we use the number of local sites
  qio_genfield_glb. n_sites = GJP. VolNodeSites();

  // Now create the record info for fireld itself

  // printf("record :%d %d = %d\n", n_fields, f_size_per_site, n_fields * f_size_per_site);
  
  QIO_RecordInfo *record_field 
    =   QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, (char*)"GenericField", (char*)"D", 0, 0,
			       n_fields * f_size_per_site*sizeof(Float), 1);
  // color=0 (not used), spin=0 (not used),  sizepersite, there will be num_fields
  
  // Do the actual work
  if(SingleDouble)
    {
      //output in double-precision

      return_val += QIO_write( qio_Output, record_field, no_string,
			       qio_getGenField,
			       n_fields* f_size_per_site*sizeof(Float),
			       sizeof(Float), field);

    }
  else
    {
      ERR.NotImplemented(cname,fname,"Sorry for lack of support for single precision sotre at this moment");
    }

  if ( (return_val == 0) ) 
    VRB.Result(cname,fname,"QIO_write successfull...\n");
  else
    ERR.General(cname,fname,"ERROR QIO: QIO_write(s) returned %i\n",return_val);


  // clean-up
  QIO_destroy_record_info( glb_rec_info );
  QIO_destroy_record_info(record_field);
  QIO_string_destroy(qio_file_header);
  QIO_string_destroy(no_string);

  qio_closeOutput();

#ifdef PROFILE
  if(!UniqueID()) printf("%s %g sec\n",cname,time_elapse());
#endif

  VRB.FuncEnd(cname,fname);

}



void qio_writeGenericFields::setHeader(const char * ensemble_id, const char * ensemble_label, const int traj,
						     		 const char * field_type_label)
{

  const char * fname = "setHeader(...)";

  VRB.Func(cname,fname);

  if(strlen(ensemble_id) > MAX_HEADER_LINE)
    {
     ERR.General(cname,fname,"ERROR: length of ensemble_id exceeds maximum length!\n");
	  exit(-1); 
    }
  else
    strcpy(header_ensemble_id, ensemble_id);

  if(strlen(ensemble_label) > MAX_HEADER_LINE)
    {
       ERR.General(cname,fname,"ERROR: length of ensemble_label exceeds maximum length!\n");
	  exit(-1); 
    }
  else
    strcpy(header_ensemble_label, ensemble_label);

  header_traj = traj;
  
  if(strlen(field_type_label) > MAX_HEADER_LINE)
    {
       ERR.General(cname,fname,"ERROR: length of field_type_label exceeds maximum length!\n");
	  exit(-1); 
    }
  else
    strcpy(header_field_type_label, field_type_label);


  VRB.FuncEnd(cname,fname);

}





CPS_END_NAMESPACE
#endif
