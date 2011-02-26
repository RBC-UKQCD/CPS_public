#ifdef USE_QIO
#include <config.h>
#include <util/qio_writePropagator.h>


// adapted functions to Qpropw-class
// 
// assuming ordering of propagator object:
// 
// [volume][sink_spin][sink_color][source_spin][source_color][ReIm]
// (volume slowest, ReIm fastest)
//
// the same ordering is applied for the source field (if full source/pairs otherwise unneeded indices removed)
//
// ees, 02/05/08



// adjust datacount issue ( 1 instead of 12 or 2, not counting complex numbers...)
//
// ees, 02/29/08

// to do:  add functions for writing timeslice, hypercubic(?) sources


CPS_START_NAMESPACE
using namespace std;


// qio-factory functions outside class!!


// qio-factory functions for full propagator or full source

#define QIO_GETPROP_BASIC(SOURCE_SPIN, SOURCE_COLOR)\
\
{\
 /*printf(" called getprop with count %i\n",count);*/\
\
 Float *prop = (Float *)arg; \
\
 Float *array = (Float *)buf; \
\
 Float *src_array = prop + index*2*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX + 2*SOURCE_COLOR + 2*QIO_PROP_COLOR_MAX*SOURCE_SPIN ;\
\
 for( int ii(0); ii < QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX; ++ii){ \
   *(array + 2*ii) = *(src_array + 2*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*ii ); \
\
   *(array + 2*ii + 1 ) = *(src_array + 2*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*ii + 1 ); } \
\
}

#define QIO_GETPROP_SINGLE_BASIC(SOURCE_SPIN, SOURCE_COLOR)\
\
{\
 /*printf(" called getprop single with index %i, s %i c %i\n",index, SOURCE_SPIN, SOURCE_COLOR);*/ \
\
 float *array = (float *)buf; \
\
 Float *prop = (double *)arg; \
\
 Float *src_array = prop + index*2*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX + 2*SOURCE_COLOR + 2*QIO_PROP_COLOR_MAX*SOURCE_SPIN;\
\
 for( int ii(0); ii < QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX; ++ii){ \
\
   *(array + 2*ii) = *(src_array + 2*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*ii); \
\
   *(array + 2*ii + 1 ) = *(src_array + 2*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*ii + 1) ; \
\
  /*printf("  %i wrote %f %f\n",ii,*(array + 2*ii), *(array+2*ii+1));*/\
\
}\
\
}


//for use in function call...
#define qio_getProp_SpinColor(SPIN, COLOR) qio_getProp_Spin_ ## SPIN ## _Color_ ## COLOR

void qio_getProp_SpinColor(0,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_BASIC(0,0)

void qio_getProp_SpinColor(0,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_BASIC(0,1)

void qio_getProp_SpinColor(0,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_BASIC(0,2)

void qio_getProp_SpinColor(1,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_BASIC(1,0)

void qio_getProp_SpinColor(1,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_BASIC(1,1)

void qio_getProp_SpinColor(1,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_BASIC(1,2)

void qio_getProp_SpinColor(2,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_BASIC(2,0)

void qio_getProp_SpinColor(2,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_BASIC(2,1)

void qio_getProp_SpinColor(2,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_BASIC(2,2)

void qio_getProp_SpinColor(3,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_BASIC(3,0)

void qio_getProp_SpinColor(3,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_BASIC(3,1)

void qio_getProp_SpinColor(3,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_BASIC(3,2)

  // single prec.

#define qio_getPropSingle_SpinColor(SPIN, COLOR) qio_getPropSingle_Spin_ ## SPIN ## _Color_ ## COLOR

void qio_getPropSingle_SpinColor(0,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_SINGLE_BASIC(0,0)

void qio_getPropSingle_SpinColor(0,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_SINGLE_BASIC(0,1)

void qio_getPropSingle_SpinColor(0,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_SINGLE_BASIC(0,2)

void qio_getPropSingle_SpinColor(1,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_SINGLE_BASIC(1,0)

void qio_getPropSingle_SpinColor(1,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_SINGLE_BASIC(1,1)

void qio_getPropSingle_SpinColor(1,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_SINGLE_BASIC(1,2)

void qio_getPropSingle_SpinColor(2,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_SINGLE_BASIC(2,0)

void qio_getPropSingle_SpinColor(2,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_SINGLE_BASIC(2,1)

void qio_getPropSingle_SpinColor(2,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_SINGLE_BASIC(2,2)

void qio_getPropSingle_SpinColor(3,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_SINGLE_BASIC(3,0)

void qio_getPropSingle_SpinColor(3,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_SINGLE_BASIC(3,1)

void qio_getPropSingle_SpinColor(3,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETPROP_SINGLE_BASIC(3,2)


  // calls for ScalarSource

#define QIO_GETSCSRC_BASIC(SOURCE_SPIN, SOURCE_COLOR)\
\
{\
 /*printf(" called getScSrc with count %i\n",count);*/\
\
 Float *source = (Float *)arg; \
\
 Float *array = (Float *)buf; \
\
 Float *src_array = source + index*2*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX + 2*SOURCE_COLOR + 2*QIO_PROP_COLOR_MAX*SOURCE_SPIN ;\
\
 for( int ii(0); ii < 2; ++ii) \
   *(array + ii) = *(src_array + ii); \
\
}

#define QIO_GETSCSRC_SINGLE_BASIC(SOURCE_SPIN, SOURCE_COLOR)\
\
{\
 /*printf(" called getScSrc with count %i\n",count);*/\
\
 Float *source = (double *)arg; \
\
 float *array = (float *)buf; \
\
 Float *src_array = source + index*2*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX + 2*SOURCE_COLOR + 2*QIO_PROP_COLOR_MAX*SOURCE_SPIN ;\
\
 for( int ii(0); ii < 2; ++ii) \
   *(array + ii) = *(src_array + ii); \
\
}

//for use in function call...
#define qio_getScSrc_SpinColor(SPIN, COLOR) qio_getScSrc_Spin_ ## SPIN ## _Color_ ## COLOR

void qio_getScSrc_SpinColor(0,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_BASIC(0,0)

void qio_getScSrc_SpinColor(0,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_BASIC(0,1)

void qio_getScSrc_SpinColor(0,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_BASIC(0,2)

void qio_getScSrc_SpinColor(1,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_BASIC(1,0)

void qio_getScSrc_SpinColor(1,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_BASIC(1,1)

void qio_getScSrc_SpinColor(1,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_BASIC(1,2)

void qio_getScSrc_SpinColor(2,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_BASIC(2,0)

void qio_getScSrc_SpinColor(2,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_BASIC(2,1)

void qio_getScSrc_SpinColor(2,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_BASIC(2,2)

void qio_getScSrc_SpinColor(3,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_BASIC(3,0)

void qio_getScSrc_SpinColor(3,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_BASIC(3,1)

void qio_getScSrc_SpinColor(3,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_BASIC(3,2)

  // single prec.

#define qio_getScSrcSingle_SpinColor(SPIN, COLOR) qio_getScSrcSingle_Spin_ ## SPIN ## _Color_ ## COLOR

void qio_getScSrcSingle_SpinColor(0,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_SINGLE_BASIC(0,0)

void qio_getScSrcSingle_SpinColor(0,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_SINGLE_BASIC(0,1)

void qio_getScSrcSingle_SpinColor(0,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_SINGLE_BASIC(0,2)

void qio_getScSrcSingle_SpinColor(1,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_SINGLE_BASIC(1,0)

void qio_getScSrcSingle_SpinColor(1,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_SINGLE_BASIC(1,1)

void qio_getScSrcSingle_SpinColor(1,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_SINGLE_BASIC(1,2)

void qio_getScSrcSingle_SpinColor(2,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_SINGLE_BASIC(2,0)

void qio_getScSrcSingle_SpinColor(2,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_SINGLE_BASIC(2,1)

void qio_getScSrcSingle_SpinColor(2,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_SINGLE_BASIC(2,2)

void qio_getScSrcSingle_SpinColor(3,0) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_SINGLE_BASIC(3,0)

void qio_getScSrcSingle_SpinColor(3,1) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_SINGLE_BASIC(3,1)

void qio_getScSrcSingle_SpinColor(3,2) (char *buf, size_t index, int count, void *arg)
  QIO_GETSCSRC_SINGLE_BASIC(3,2)


  // calls for ScalarSource (cannot be replaced by qio_getScSource_SpinColor(0,0)


void qio_getScSource(char *buf, size_t index, int count, void *arg)
{

#ifdef DEBUG_GetProp
  printf("UID: %i, called qio_getScSource with index: %i, count: %i\n",UniqueID(), index, count);
#endif // DEBUG_GetProp

  Float *source = (Float *)arg;

  Float *array = (Float *) buf;

  Float *src_array = source + index*2;

  for (int ii(0); ii < 2; ++ii)
    *(array + ii)  = *(src_array + ii); 
  

#ifdef DEBUG_GetProp
  printf("UID: %i, finished qio_getScSource (called with index: %i, count: %i)\n",UniqueID(), index, count);
#endif // DEBUG_GetProp

}


void qio_getScSourceSingle(char *buf, size_t index, int count, void *arg)
{


#ifdef DEBUG_GetProp
  printf("UID: %i, called qio_ScSourceSingle with index: %i, count: %i\n",UniqueID(), index, count);
#endif // DEBUG_GetProp

  float *array = (float *) buf;

  Float *source = (double *)arg;

  Float *src_array = source + index * 2;

  for(int ii(0); ii < 2; ++ii)
    *(array + ii) = *(src_array + ii);

#ifdef DEBUG_GetProp
  printf("UID: %i, finished qio_getScSourceSingle (called with index: %i, count: %i)\n",UniqueID(), index, count);
#endif // DEBUG_GetProp

}




// now start class-functions...



void qio_writePropagator::qio_openOutput(char *filename, QIO_String *record_file, int volFormat)
{

  char * fname = "qio_openOutput(...)";
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

  outfile = QIO_open_write(record_file, filename, volfmt, &layout, NULL, &oflag);

  qio_Output = outfile;

  VRB.FuncEnd(cname,fname);

}




void qio_writePropagator::write_ScS_12sink(char *outfile, const void *prop, const void *scSource, int volFormat, FP_FORMAT floatFormat)
{

  const char * fname = "write_ScS_12sink(...)";

  VRB.Func(cname,fname);

  VRB.Result(cname,fname,"writing propagator format:\n    scalar source plus 12 sinks\n    to %s \n", outfile);

  int return_val(0);

  const int sizepersite(QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX);

  Float *wprop = (Float *)prop;

  Float *wsource = (Float *)scSource;


  QIO_RecordInfo *record_prop;
  record_prop = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0,0,0,0);

  QIO_RecordInfo *record_source;
  record_source = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0,0,0,0);

  QIO_String *record_xml_prop;
  record_xml_prop = QIO_string_create();

  QIO_String *record_xml_source;
  record_xml_source = QIO_string_create();

  QIO_String *record_file;
  record_file = QIO_string_create();

  qio_setLayout();

  //detect output format

  int SingleDouble(1);
  // 1=double (standard), else single

  switch( floatFormat )
    {
    case FP_TIDSP32 :      SingleDouble=0; break;
    case FP_IEEE32 :       SingleDouble=0; break;
    case FP_IEEE32BIG :    SingleDouble=0; break;
    case FP_IEEE32LITTLE : SingleDouble=0; break;

    default: SingleDouble=1;
    }

  if (SingleDouble) VRB.Flow(cname,fname," output-precision: DOUBLE\n");
  else VRB.Flow(cname,fname," output-precision: SINGLE\n");

  
  char xml_info_file[7*(MAX_HEADER_LINE+10)];
  char xml_info_prop[6*(MAX_HEADER_LINE+10)];
  char xml_info_source[6*(MAX_HEADER_LINE+10)];

  sprintf(xml_info_file, "DATATYPE = 4D_PROPAGATOR\nPROPAGATORTYPE = %s\nSOURCETYPE = %s\nENSEMBLE_ID = %s\nENSEMBLE_LABEL = %s\nSEQUENCE_NUMBER = %i",
  	  header_propagator_type, header_source_type, header_ensemble_id, header_ensemble_label, header_traj);

  sprintf(xml_info_prop, "DATATYPE = 4D_PROPAGATOR\nPROPAGATORTYPE = %s\nENSEMBLE_ID = %s\nENSEMBLE_LABEL = %s\nSEQUENCE_NUMBER = %i",
  	  header_propagator_type, header_ensemble_id, header_ensemble_label, header_traj);

  sprintf(xml_info_source, "DATATYPE = 4D_PROPAGATOR\nSOURCETYPE = %s\nENSEMBLE_ID = %s\nENSEMBLE_LABEL = %s\nSEQUENCE_NUMBER = %i",
  	  header_source_type, header_ensemble_id, header_ensemble_label, header_traj);

  
  //check for length...
  if( (strlen(xml_info_file)*sizeof(char)) > QIO_INFO_STRING_MAX){
    ERR.General(cname,fname," xml_info_file too large: %i\n %s\n", (strlen(xml_info_file)*sizeof(char)), xml_info_file);
    exit(-12);
  }
  if( (strlen(xml_info_prop)*sizeof(char)) > QIO_INFO_STRING_MAX){
    ERR.General(cname,fname," xml_info_prop too large: %i\n %s\n", (strlen(xml_info_prop)*sizeof(char)), xml_info_prop);
    exit(-12);
  }
  if( (strlen(xml_info_source)*sizeof(char)) > QIO_INFO_STRING_MAX){
    ERR.General(cname,fname," xml_info_source too large: %i\n %s\n", (strlen(xml_info_source)*sizeof(char)), xml_info_source);
    exit(-12);
  }

  QIO_USQCDPropFileInfo *filerecordinfo = QIO_create_usqcd_propfile_info(QIO_USQCDPROPFILETYPE_C1D12, xml_info_file);

  QIO_encode_usqcd_propfile_info(record_file, filerecordinfo);

  QIO_USQCDPropRecordInfo *userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0, 0, xml_info_prop);

  QIO_USQCDPropSourceInfo *userrecordinfo_source = QIO_create_usqcd_propsource_info(xml_info_source);


  qio_openOutput(outfile, record_file, volFormat);
  

  if(SingleDouble)
    {
      //output in double-precision

      // create the record info
      record_prop = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_D3_DiracFermion", "D", 3, 4, 2*sizepersite*sizeof(Float), 1);
      // color=3 (not used), spin=4 (not used), size 2 (complex), count per site

      if(source_hypercube)
	record_source = QIO_create_record_info(QIO_HYPER, source_start, source_end, QIO_RW_DIMENSION, "USQCD_D3_Complex", "D", 0, 0, 2*sizeof(Float), 1);
      else
	record_source = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_D3_Complex", "D", 0, 0, 2*sizeof(Float), 1);
      // color, spin not used, size 2(complex)

      return_val = 0 ;

      QIO_encode_usqcd_propsource_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSource, 2*sizeof(Float), sizeof(Float), wsource);
      

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 0);
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0, 0, xml_info_prop);
#endif

      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0, 1, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0, 2, xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 1);

      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(1, 0, xml_info_prop);      
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(1, 1, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(1, 2, xml_info_prop);      
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 2);

      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(2, 0, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(2, 1, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(2, 2, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 3);

      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(3, 0, xml_info_prop);      
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(3, 1, xml_info_prop);      
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(3, 2, xml_info_prop);      
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);


    }
  else
    {
      //output in single-precision

      // create the record info
      record_prop = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_F3_DiracFermion", "F", 3, 4, 2*sizepersite*sizeof(float), 1);
      // color=3 (not used), spin=4 (not used), size 2 (complex), count per side

      if(source_hypercube)
	record_source = QIO_create_record_info(QIO_HYPER, source_start, source_end, QIO_RW_DIMENSION, "USQCD_F3_Complex", "F", 0, 0, 2*sizeof(float), 1);
      else
	record_source = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_F3_Complex", "F", 0, 0, 2*sizeof(float), 1);
      // color, spin not used, size 2(complex)

      //now scalarSource + 12 sinks
      return_val = 0 ;

      QIO_encode_usqcd_propsource_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSourceSingle, 2*sizeof(float), sizeof(float), wsource);
      

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 0);

      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0, 0, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(0,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0, 1, xml_info_prop);
#endif
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(0,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0, 2, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(0,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 1);

      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(1, 0, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(1,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(1, 1, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(1,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(1, 2, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);      
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(1,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 2);
 
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(2, 0, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(2,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(2, 1, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);      
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(2,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(2, 2, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(2,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 3);

      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(3, 0, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(3,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(3, 1, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(3,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else      
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(3, 2, xml_info_prop);
#endif      
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);      
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(3,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
    }
  

  if ( (return_val == 0) ) 
    VRB.Result(cname,fname,"QIO_write successfull...\n");
  else
    ERR.General(cname,fname,"ERROR QIO: QIO_write(s) returned %i\n",return_val);


#ifdef PRINT_checksums
  uint32_t writeCheckA(QIO_get_writer_last_checksuma( output));
  uint32_t writeCheckB(QIO_get_writer_last_checksumb( output));
 
  printf("Checksums: a: %s -- b: %s \n",writeCheckA, writeCheckB);
#endif //PRINT_checksums

  // clean-up

  QIO_destroy_record_info(record_prop);
  QIO_destroy_record_info(record_source);

  
  QIO_string_destroy(record_xml_prop);
  QIO_string_destroy(record_xml_source);

  QIO_destroy_usqcd_propfile_info(filerecordinfo);
  QIO_destroy_usqcd_proprecord_info(userrecordinfo_prop);
  QIO_destroy_usqcd_propsource_info(userrecordinfo_source);

  QIO_string_destroy(record_file);

  qio_closeOutput();

  VRB.FuncEnd(cname,fname);

}


void qio_writePropagator::write_12pairs(char *outfile, const QIO_PROP_SOURCE_TYPES sType, const void *prop, const void *source, 
					int volFormat, FP_FORMAT floatFormat)
{

  const char * fname = "write_12pairs(...)";

  VRB.Func(cname,fname);

  VRB.Result(cname,fname,"writing propagator format:\n    12 pairs of %s source and sink\n    to %s \n",
	     ( sType==QIO_FULL_SOURCE ? "full" : ( sType==QIO_SCALAR_SOURCE ? "scalar" : "unknown") ),
	     outfile);



  int return_val(0);

  const int sizepersite(QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX);

  int source_size(0);

  Float *wprop = (Float *)prop;

  Float *wsource = (Float *)source;

  QIO_RecordInfo *record_prop;
  record_prop = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0,0,0,0);

  QIO_RecordInfo *record_source;
  record_source = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0,0,0,0);

  QIO_String *record_xml_prop;
  record_xml_prop = QIO_string_create();

  QIO_String *record_xml_source;
  record_xml_source = QIO_string_create();

  QIO_String *record_file;
  record_file = QIO_string_create();

  qio_setLayout();
  
  //detect output format

  int SingleDouble(1);
  // 1=double (standard), else single
  
  switch( floatFormat )
    {
    case FP_TIDSP32 :      SingleDouble=0; break;
    case FP_IEEE32 :       SingleDouble=0; break;
    case FP_IEEE32BIG :    SingleDouble=0; break;
    case FP_IEEE32LITTLE : SingleDouble=0; break;

    default: SingleDouble=1;
    }

  if (SingleDouble) VRB.Flow(cname,fname," output-precision: DOUBLE\n");
  else VRB.Flow(cname,fname," output-precision: SINGLE\n");

  char xml_info_file[7*(MAX_HEADER_LINE+10)];
  char xml_info_prop[6*(MAX_HEADER_LINE+10)];
  char xml_info_source[6*(MAX_HEADER_LINE+10)];



  sprintf(xml_info_file, "DATATYPE = 4D_PROPAGATOR\nPROPAGATORTYPE = %s\nSOURCETYPE = %s\nENSEMBLE_ID = %s\nENSEMBLE_LABEL = %s\nSEQUENCE_NUMBER = %i",
  	  header_propagator_type, header_source_type, header_ensemble_id, header_ensemble_label, header_traj);

  sprintf(xml_info_prop, "DATATYPE = 4D_PROPAGATOR\nPROPAGATORTYPE = %s\nENSEMBLE_ID = %s\nENSEMBLE_LABEL = %s\nSEQUENCE_NUMBER = %i",
  	  header_propagator_type, header_ensemble_id, header_ensemble_label, header_traj);

  sprintf(xml_info_source, "DATATYPE = 4D_PROPAGATOR\nSOURCETYPE = %s\nENSEMBLE_ID = %s\nENSEMBLE_LABEL = %s\nSEQUENCE_NUMBER = %i",
  	  header_source_type, header_ensemble_id, header_ensemble_label, header_traj);

  //check for length...
  if( (strlen(xml_info_file)*sizeof(char)) > QIO_INFO_STRING_MAX){
    ERR.General(cname,fname," xml_info_file too large: %i\n %s\n", (strlen(xml_info_file)*sizeof(char)), xml_info_file);
    exit(-12);
  }
  if( (strlen(xml_info_prop)*sizeof(char)) > QIO_INFO_STRING_MAX){
    ERR.General(cname,fname," xml_info_prop too large: %i\n %s\n", (strlen(xml_info_prop)*sizeof(char)), xml_info_prop);
    exit(-12);
  }
  if( (strlen(xml_info_source)*sizeof(char)) > QIO_INFO_STRING_MAX){
    ERR.General(cname,fname," xml_info_source too large: %i\n %s\n", (strlen(xml_info_source)*sizeof(char)), xml_info_source);
    exit(-12);
  }



  QIO_USQCDPropFileInfo *filerecordinfo;


  switch(sType)
    {
    case 1:
      filerecordinfo = QIO_create_usqcd_propfile_info(QIO_USQCDPROPFILETYPE_CD_PAIRS, xml_info_file);
      break;

    case 2:
      filerecordinfo = QIO_create_usqcd_propfile_info(QIO_USQCDPROPFILETYPE_DD_PAIRS, xml_info_file);
      break;

    default:
      ERR.General(cname,fname,"ERROR QIO: unknown source type %s\n",sType);
      filerecordinfo = QIO_create_usqcd_propfile_info(0,xml_info_file);


    }


  QIO_encode_usqcd_propfile_info(record_file, filerecordinfo);

  QIO_USQCDPropRecordInfo *userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0,0,xml_info_prop);

  QIO_USQCDPropSourceInfo *userrecordinfo_source = QIO_create_usqcd_propsource_info(xml_info_source);

  qio_openOutput(outfile, record_file, volFormat);


  if(SingleDouble)
    {
      //output in double-precision

      // create the record info
      record_prop = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_D3_DiracFermion", "D", 3, 4, 2*sizepersite*sizeof(Float), 1);
      // color=3 (not used), spin=4 (not used), size 2 (complex), count per site


      if( sType == QIO_SCALAR_SOURCE)  //1  scalar source
	{
	  if(source_hypercube)
	    record_source = QIO_create_record_info(QIO_HYPER, source_start, source_end, QIO_RW_DIMENSION, "USQCD_D3_Complex", "D", 0, 0, 2*sizeof(Float), 1);
	  else
	    record_source = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_D3_Complex", "D", 0, 0, 2*sizeof(Float), 1);
	  // color, spin not used, size 2(complex)
	  source_size = 1;
	}
      
      if( sType == QIO_FULL_SOURCE ) // 2)  full source
	{
	  if(source_hypercube)
	    record_source = QIO_create_record_info(QIO_HYPER, source_start, source_end, QIO_RW_DIMENSION, "USQCD_D3_DiracFermion", "D", 3, 4, 2*sizepersite*sizeof(Float), 1);
	  else
	    record_source = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_D3_DiracFermion", "D", 3, 4, 2*sizepersite*sizeof(Float), 1);
	  // color=3 (not used), spin=4 (not used), size 2(complex), count per site
	  source_size = sizepersite;
	}

      return_val = 0 ;



      // now 12 pairs of source, sink, ordered by spin, color.

      QIO_encode_usqcd_propsource_info(record_xml_source, userrecordinfo_source);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 0);
      
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0,0,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(0,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(0,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0,1,xml_info_prop);
#endif
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(0,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(0,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0,2,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(0,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      else
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(0,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 1);
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(1,0,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(1,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(1,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(1,1,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(1,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(1,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(1,2,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(1,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(1,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 2);
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(2,0,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(2,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(2,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(2,1,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(2,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(2,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(2,2,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(2,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(2,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 3);
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(3,0,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(3,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(3,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(3,1,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(3,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(3,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(3,2,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(3,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(3,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);


    }
  else
    {
      //output in single-precision

      // create the record info
      record_prop = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_F3_DiracFermion", "F", 3, 4, 2*sizepersite*sizeof(float), 1);
      // color=3 (not used), spin=4 (not used), size 2 (complex), count per site


      if( sType == QIO_SCALAR_SOURCE)  // scalar source
	{
	  if(source_hypercube)
	    record_source = QIO_create_record_info(QIO_HYPER, source_start, source_end, QIO_RW_DIMENSION, "USQCD_F3_Complex", "F", 0, 0, 2*sizeof(float), 1);
	  else
	    record_source = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_F3_Complex", "F", 0, 0, 2*sizeof(float), 1);
	  // color, spin not used, size 2(complex)
	  source_size = 1;
	}

      if( sType == QIO_FULL_SOURCE) // full source
	{
	  if(source_hypercube)
	    record_source = QIO_create_record_info(QIO_HYPER, source_start, source_end, QIO_RW_DIMENSION, "USQCD_F3_DiracFermion", "F", 3, 4, 2*sizepersite*sizeof(float), 1);
	  else
	    record_source = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_F3_DiracFermion", "F", 3, 4, 2*sizepersite*sizeof(float), 1);
	  // color=3 (not used), spin=4 (not used), size 2(complex), count per site
	  source_size = sizepersite;
	}

      return_val = 0 ;


      // now 12 pairs of source, sink, ordered by spin, color.


      QIO_encode_usqcd_propsource_info(record_xml_source, userrecordinfo_source);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 0);
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0,0,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(0,0), 2*source_size*sizeof(float), sizeof(float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(0,0), 2*source_size*sizeof(float), sizeof(float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(0,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0,1,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(0,1), 2*source_size*sizeof(float), sizeof(float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(0,1), 2*source_size*sizeof(float), sizeof(float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(0,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0,2,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(0,2), 2*source_size*sizeof(float), sizeof(float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(0,2), 2*source_size*sizeof(float), sizeof(float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(0,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 1);
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(1,0,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(1,0), 2*source_size*sizeof(float), sizeof(float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(1,0), 2*source_size*sizeof(float), sizeof(float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(1,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(1,1,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(1,1), 2*source_size*sizeof(float), sizeof(float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(1,1), 2*source_size*sizeof(float), sizeof(float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(1,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(1,2,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(1,2), 2*source_size*sizeof(float), sizeof(float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(1,2), 2*source_size*sizeof(float), sizeof(float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(1,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 2);
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(2,0,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(2,0), 2*source_size*sizeof(float), sizeof(float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(2,0), 2*source_size*sizeof(float), sizeof(float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(2,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(2,1,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(2,1), 2*source_size*sizeof(float), sizeof(float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(2,1), 2*source_size*sizeof(float), sizeof(float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(2,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(2,2,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(2,2), 2*source_size*sizeof(float), sizeof(float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(2,2), 2*source_size*sizeof(float), sizeof(float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(2,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, 3);
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 0);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(3,0,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(3,0), 2*source_size*sizeof(float), sizeof(float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(3,0), 2*source_size*sizeof(float), sizeof(float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(3,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 1);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(3,1,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(3,1), 2*source_size*sizeof(float), sizeof(float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(3,1), 2*source_size*sizeof(float), sizeof(float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(3,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, 2);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(3,2,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      if( sType == QIO_SCALAR_SOURCE)
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(3,2), 2*source_size*sizeof(float), sizeof(float), wsource);
      else
	return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(3,2), 2*source_size*sizeof(float), sizeof(float), wsource);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(3,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);

    }


  if ( (return_val == 0) ) 
    VRB.Result(cname,fname,"QIO_write successfull...\n");
  else
    ERR.General(cname,fname,"ERROR QIO: QIO_write(s) returned %i\n",return_val);


  // clean-up

  QIO_destroy_record_info(record_prop);
  QIO_destroy_record_info(record_source);

  
  QIO_string_destroy(record_xml_prop);
  QIO_string_destroy(record_xml_source);

  QIO_destroy_usqcd_propfile_info(filerecordinfo);
  QIO_destroy_usqcd_proprecord_info(userrecordinfo_prop);
  QIO_destroy_usqcd_propsource_info(userrecordinfo_source);

  QIO_string_destroy(record_file);

  qio_closeOutput();



  VRB.FuncEnd(cname,fname);

}





void qio_writePropagator::write_pair(char *outfile, const QIO_PROP_SOURCE_TYPES sType, const void *prop, const void *source,
				     int spin,
				     int color,
				     int volFormat, FP_FORMAT floatFormat)
{
  
  const char * fname = "write_pair(...)";

  VRB.Func(cname,fname);

  VRB.Result(cname,fname,"writing propagator format:\n  %s source %d %d and sink \n    to %s \n",
	     ( sType==QIO_FULL_SOURCE ? "full" : ( sType==QIO_SCALAR_SOURCE ? "scalar" : "unknown") ),
	     spin,
	     color,
	     outfile);

  int return_val(0);

  const int sizepersite(QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX);

  int source_size(0);

  Float *wprop = (Float *)prop;

  Float *wsource = (Float *)source;

  QIO_RecordInfo *record_prop;
  record_prop = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0,0,0,0);

  QIO_RecordInfo *record_source;
  record_source = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0,0,0,0);

  QIO_String *record_xml_prop;
  record_xml_prop = QIO_string_create();

  QIO_String *record_xml_source;
  record_xml_source = QIO_string_create();

  QIO_String *record_file;
  record_file = QIO_string_create();

  qio_setLayout();
  
  //detect output format

  int SingleDouble(1);
  // 1=double (standard), else single
  
  switch( floatFormat )
    {
    case FP_TIDSP32 :      SingleDouble=0; break;
    case FP_IEEE32 :       SingleDouble=0; break;
    case FP_IEEE32BIG :    SingleDouble=0; break;
    case FP_IEEE32LITTLE : SingleDouble=0; break;

    default: SingleDouble=1;
    }

  if (SingleDouble) VRB.Flow(cname,fname," output-precision: DOUBLE\n");
  else VRB.Flow(cname,fname," output-precision: SINGLE\n");

  char xml_info_file[7*(MAX_HEADER_LINE+10)];
  char xml_info_prop[6*(MAX_HEADER_LINE+10)];
  char xml_info_source[6*(MAX_HEADER_LINE+10)];



  sprintf(xml_info_file, "DATATYPE = 4D_PROPAGATOR\nPROPAGATORTYPE = %s\nSOURCETYPE = %s\nENSEMBLE_ID = %s\nENSEMBLE_LABEL = %s\nSEQUENCE_NUMBER = %i",
  	  header_propagator_type, header_source_type, header_ensemble_id, header_ensemble_label, header_traj);

  sprintf(xml_info_prop, "DATATYPE = 4D_PROPAGATOR\nPROPAGATORTYPE = %s\nENSEMBLE_ID = %s\nENSEMBLE_LABEL = %s\nSEQUENCE_NUMBER = %i",
  	  header_propagator_type, header_ensemble_id, header_ensemble_label, header_traj);

  sprintf(xml_info_source, "DATATYPE = 4D_PROPAGATOR\nSOURCETYPE = %s\nENSEMBLE_ID = %s\nENSEMBLE_LABEL = %s\nSEQUENCE_NUMBER = %i",
  	  header_source_type, header_ensemble_id, header_ensemble_label, header_traj);

  //check for length...
  if( (strlen(xml_info_file)*sizeof(char)) > QIO_INFO_STRING_MAX){
    ERR.General(cname,fname," xml_info_file too large: %i\n %s\n", (strlen(xml_info_file)*sizeof(char)), xml_info_file);
    exit(-12);
  }
  if( (strlen(xml_info_prop)*sizeof(char)) > QIO_INFO_STRING_MAX){
    ERR.General(cname,fname," xml_info_prop too large: %i\n %s\n", (strlen(xml_info_prop)*sizeof(char)), xml_info_prop);
    exit(-12);
  }
  if( (strlen(xml_info_source)*sizeof(char)) > QIO_INFO_STRING_MAX){
    ERR.General(cname,fname," xml_info_source too large: %i\n %s\n", (strlen(xml_info_source)*sizeof(char)), xml_info_source);
    exit(-12);
  }



  QIO_USQCDPropFileInfo *filerecordinfo;


  switch(sType)
    {
    case 1:
      filerecordinfo = QIO_create_usqcd_propfile_info(QIO_USQCDPROPFILETYPE_CD_PAIRS, xml_info_file);
      break;

    case 2:
      filerecordinfo = QIO_create_usqcd_propfile_info(QIO_USQCDPROPFILETYPE_DD_PAIRS, xml_info_file);
      break;

    default:
      ERR.General(cname,fname,"ERROR QIO: unknown source type %s\n",sType);
      filerecordinfo = QIO_create_usqcd_propfile_info(0,xml_info_file);


    }


  QIO_encode_usqcd_propfile_info(record_file, filerecordinfo);

  QIO_USQCDPropRecordInfo *userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(0,0,xml_info_prop);

  QIO_USQCDPropSourceInfo *userrecordinfo_source = QIO_create_usqcd_propsource_info(xml_info_source);

  qio_openOutput(outfile, record_file, volFormat);


  if(SingleDouble)
    {
      //output in double-precision

      // create the record info
      record_prop = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_D3_DiracFermion", "D", 3, 4, 2*sizepersite*sizeof(Float), 1);
      // color=3 (not used), spin=4 (not used), size 2 (complex), count per site


      if( sType == QIO_SCALAR_SOURCE)  //1  scalar source
	{
	  if(source_hypercube)
	    record_source = QIO_create_record_info(QIO_HYPER, source_start, source_end, QIO_RW_DIMENSION, "USQCD_D3_Complex", "D", 0, 0, 2*sizeof(Float), 1);
	  else
	    record_source = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_D3_Complex", "D", 0, 0, 2*sizeof(Float), 1);
	  // color, spin not used, size 2(complex)
	  source_size = 1;
	}
      
      if( sType == QIO_FULL_SOURCE ) // 2)  full source
	{
	  if(source_hypercube)
	    record_source = QIO_create_record_info(QIO_HYPER, source_start, source_end, QIO_RW_DIMENSION, "USQCD_D3_DiracFermion", "D", 3, 4, 2*sizepersite*sizeof(Float), 1);
	  else
	    record_source = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_D3_DiracFermion", "D", 3, 4, 2*sizepersite*sizeof(Float), 1);
	  // color=3 (not used), spin=4 (not used), size 2(complex), count per site
	  source_size = sizepersite;
	}

      return_val = 0 ;



      // now source, sink pair of spin, color.

      QIO_encode_usqcd_propsource_info(record_xml_source, userrecordinfo_source);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, spin);
      
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, color);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(spin,color,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);
      switch(spin){
      case 0:
	switch(color){
	case 0:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(0,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(0,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 1:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(0,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(0,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 2:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(0,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(0,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	default:
	  break;
	}
	break;
      case 1:
	switch(color){
	case 0:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(1,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(1,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 1:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(1,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(1,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 2:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(1,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(1,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	default:
	  break;
	}
	break;
      case 2:
	switch(color){
	case 0:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(2,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(2,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 1:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(2,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(2,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 2:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(2,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(2,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	default:
	  break;
	}
	break;
      case 3:
	switch(color){
	case 0:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(3,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(3,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 1:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(3,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(3,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 2:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrc_SpinColor(3,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(3,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	default:
	  break;
	}
	break;
      default:
	break;
      }
    }
  else
    {
      //output in single-precision

      // create the record info
      record_prop = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_F3_DiracFermion", "F", 3, 4, 2*sizepersite*sizeof(float), 1);
      // color=3 (not used), spin=4 (not used), size 2 (complex), count per site


      if( sType == QIO_SCALAR_SOURCE)  // scalar source
	{
	  if(source_hypercube)
	    record_source = QIO_create_record_info(QIO_HYPER, source_start, source_end, QIO_RW_DIMENSION, "USQCD_F3_Complex", "F", 0, 0, 2*sizeof(float), 1);
	  else
	    record_source = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_F3_Complex", "F", 0, 0, 2*sizeof(float), 1);
	  // color, spin not used, size 2(complex)
	  source_size = 1;
	}

      if( sType == QIO_FULL_SOURCE) // full source
	{
	  if(source_hypercube)
	    record_source = QIO_create_record_info(QIO_HYPER, source_start, source_end, QIO_RW_DIMENSION, "USQCD_F3_DiracFermion", "F", 3, 4, 2*sizepersite*sizeof(float), 1);
	  else
	    record_source = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0, "USQCD_F3_DiracFermion", "F", 3, 4, 2*sizepersite*sizeof(float), 1);
	  // color=3 (not used), spin=4 (not used), size 2(complex), count per site
	  source_size = sizepersite;
	}

      return_val = 0 ;


      // now source, sink, pair for spin, color.

      QIO_encode_usqcd_propsource_info(record_xml_source, userrecordinfo_source);

#ifdef SPINCOLORINSERT
      QIO_insert_usqcd_proprecord_spin( userrecordinfo_prop, spin);
      QIO_insert_usqcd_proprecord_color( userrecordinfo_prop, color);
#else
      userrecordinfo_prop = QIO_create_usqcd_proprecord_sc_info(spin,color,xml_info_prop);
#endif     
      QIO_encode_usqcd_proprecord_info(record_xml_prop, userrecordinfo_prop);

      switch(spin){
      case 0:
	switch(color){
	case 0:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(0,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(0,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(0,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 1:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(0,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(0,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(0,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 2:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(0,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(0,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(0,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	default:
	  break;
	}
	break;
      case 1:
	switch(color){
	case 0:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(1,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(1,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(1,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 1:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(1,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(1,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(1,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 2:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(1,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(1,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(1,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	default:
	  break;
	}
	break;
      case 2:
	switch(color){
	case 0:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(2,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(2,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(2,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 1:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(2,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(2,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(2,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 2:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(2,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(2,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(2,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	default:
	  break;
	}
	break;
      case 3:
	switch(color){
	case 0:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(3,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(3,0), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(3,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 1:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(3,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(3,1), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(3,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	case 2:
	  if( sType == QIO_SCALAR_SOURCE)
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSrcSingle_SpinColor(3,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  else
	    return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getPropSingle_SpinColor(3,2), 2*source_size*sizeof(Float), sizeof(Float), wsource);
	  return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(3,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
	  break;
	default:
	  break;
	}
	break;
      default:
	break;
      }
    }


  if ( (return_val == 0) ) 
    VRB.Result(cname,fname,"QIO_write successfull...\n");
  else
    ERR.General(cname,fname,"ERROR QIO: QIO_write(s) returned %i\n",return_val);


  // clean-up

  QIO_destroy_record_info(record_prop);
  QIO_destroy_record_info(record_source);

  
  QIO_string_destroy(record_xml_prop);
  QIO_string_destroy(record_xml_source);

  QIO_destroy_usqcd_propfile_info(filerecordinfo);
  QIO_destroy_usqcd_proprecord_info(userrecordinfo_prop);
  QIO_destroy_usqcd_propsource_info(userrecordinfo_source);

  QIO_string_destroy(record_file);

  qio_closeOutput();



  VRB.FuncEnd(cname,fname);

}










void qio_writePropagator::setHeader(const char * ensemble_id, const char * ensemble_label, const int traj)
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
  
  VRB.FuncEnd(cname,fname);

}


void qio_writePropagator::setHeader(const char * propagator_type, const char * source_type)
{

  const char * fname = "setHeader(...)";

  VRB.Func(cname,fname);

  if(strlen(propagator_type) > MAX_HEADER_LINE)
    {
     ERR.General(cname,fname,"ERROR: length of propagator_type exceeds maximum length!\n");
	  exit(-1); 
    }
  else
    strcpy(header_propagator_type, propagator_type);

  if(strlen(source_type) > MAX_HEADER_LINE)
    {
       ERR.General(cname,fname,"ERROR: length of source_type exceeds maximum length!\n");
	  exit(-1); 
    }
  else
    strcpy(header_source_type, source_type);


  
  VRB.FuncEnd(cname,fname);

}


void qio_writePropagator::setSourceTslice( const int tslice)
{

  const char * fname = "setSourceTslice(...)";

  VRB.Func(cname,fname);

  if( QIO_RW_DIMENSION > 4)
    ERR.General(cname, fname," not implemented for more than 4 dimensions");

  if( (tslice < 0)  || (tslice >= GJP.Sites(3)) )
    ERR.General(cname, fname," tslice %i out of range",tslice);

  source_hypercube=1;

  for(int ii(0); ii<3; ++ii){
    source_start[ii] = 0;
    source_end[ii] = GJP.Sites(ii)-1;
  }
  source_start[3] = tslice;
  source_end[3] = tslice;


  VRB.FuncEnd(cname,fname);
}

void qio_writePropagator::setSourceTslices( const int t_start, const int t_end)
{

  const char * fname = "setSourceTslices(...)";

  VRB.Func(cname,fname);

  if( QIO_RW_DIMENSION > 4)
    ERR.General(cname, fname," not implemented for more than 4 dimensions");

  if( (t_start < 0)  || (t_start > t_end) || t_end >= GJP.Sites(3) )
    ERR.General(cname, fname," t_start %i , t_end %i out of range",t_start, t_end);

  source_hypercube=1;

  for(int ii(0); ii<3; ++ii){
    source_start[ii] = 0;
    source_end[ii] = GJP.Sites(ii)-1;
  }
  source_start[3] = t_start;
  source_end[3] = t_end;


  VRB.FuncEnd(cname,fname);
}

void qio_writePropagator::setSourceHypercube( const int start[4], const int end[4])
{

  const char * fname = "setSourceHypercube(...)";

  VRB.Func(cname,fname);

  source_hypercube=1;

  for(int ii(0); ii<QIO_RW_DIMENSION; ++ii){

    if( ii >= 4)
      ERR.General(cname, fname," not implemented for more than 4 dimensions");

    if( (start[ii] < 0) || (start[ii] > end[ii]) || ( end[ii] >= GJP.Sites(ii) ) )
      ERR.General(cname, fname," range %i = %i, %i out of range",ii,start[ii],end[ii]);

    source_start[ii] = start[ii];
    source_end[ii] = end[ii];
  }


  VRB.FuncEnd(cname,fname);
}


CPS_END_NAMESPACE
#endif
