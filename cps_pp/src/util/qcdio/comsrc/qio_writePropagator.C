#ifdef USE_QIO
#include <config.h>
#include <util/qio_writePropagator.h>



CPS_START_NAMESPACE
using namespace std;


// qio-factory functions outside class!!



#define QIO_GETPROP_BASIC(SOURCE_SPIN, SOURCE_COLOR)\
\
{\
 Float *prop = (Float *)arg; \
\
 Float *array = (Float *)buf; \
\
 Float *src_array = prop + index*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX*count + SOURCE_SPIN*QIO_PROP_COLOR_MAX*count + SOURCE_COLOR*count; \
\
 for (int ii(0); ii < count; ++ii) \
    *(array + ii)  = *(src_array + ii); \
\
}

#define QIO_GETPROP_SINGLE_BASIC(SOURCE_SPIN, SOURCE_COLOR)\
\
{\
 float *array = (float *)buf; \
\
 double *prop = (double *)arg; \
\
 double *src_array = prop +  index*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX*count + SOURCE_SPIN*QIO_PROP_COLOR_MAX*count + SOURCE_COLOR*count; \
\
 for(int ii(0); ii < count; ++ii) \
    *(array + ii) = *(src_array + ii); \
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

  // these may be replaced by identical getProp calls! (just count differs...)

void qio_getScSource(char *buf, size_t index, int count, void *arg)
{

#ifdef DEBUG_GetProp
  printf("UID: %i, called qio_getScSource with index: %i, count: %i\n",UniqueID(), index, count);
#endif // DEBUG_GetProp

  Float *source = (Float *)arg;

  Float *array = (Float *) buf;

  Float *src_array = source + index*count;

  for (int ii(0); ii < count; ++ii)
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

  double *source = (double *)arg;

  double *src_array = source + index * count;

  for(int ii(0); ii < count; ++ii)
    *(array + ii) = *(src_array + ii);

#ifdef DEBUG_GetProp
  printf("UID: %i, finished qio_getScSourceSingle (called with index: %i, count: %i)\n",UniqueID(), index, count);
#endif // DEBUG_GetProp

}


// scalar Source for spin/color pairs





// end scalar Source for spin/color pairs




  // old calls for whole propagator...
  // may be used for ARBIRTARY paired format, if order is [record][volume][spin][color][C]  ( not what CPS is using!!!)
  // source can be treated the same way, just different size if scalar source!

void qio_getProp(char *buf, size_t index, int count, void *arg)
{

#ifdef DEBUG_GetProp
  printf("UID: %i, called qio_getProp with index: %i, count: %i\n",UniqueID(), index, count);
#endif // DEBUG_GetProp

  Float *prop = (Float *)arg;

  Float *array = (Float *) buf;

  Float *src_array = prop + index*count;

  for (int ii(0); ii < count; ++ii)
    *(array + ii)  = *(src_array + ii); 
  

#ifdef DEBUG_GetProp
  printf("UID: %i, finished qio_getProp (called with index: %i, count: %i)\n",UniqueID(), index, count);
#endif // DEBUG_GetProp

}


void qio_getPropSingle(char *buf, size_t index, int count, void *arg)
{


#ifdef DEBUG_GetProp
  printf("UID: %i, called qio_getPropSingle with index: %i, count: %i\n",UniqueID(), index, count);
#endif // DEBUG_GetProp

  float *array = (float *) buf;

  double *prop = (double *)arg;

  double *src_array = prop + index * count;

  for(int ii(0); ii < count; ++ii)
    *(array + ii) = *(src_array + ii);

#ifdef DEBUG_GetProp
  printf("UID: %i, finished qio_getPropSingle (called with index: %i, count: %i)\n",UniqueID(), index, count);
#endif // DEBUG_GetProp

}

// global functions can be taken from qio_writeLattice MOVE TO GENERAL ???



// now start class-functions...



void qio_writePropagator::qio_openOutput(char *filename, QIO_String *record_file, int volFormat)
{

  char * fname = "qio_openOutput(...)";
  VRB.Func(cname,fname);

  const int volfmt(volFormat);
  const int serpar(QIO_SERPAR);
  const int ildgstyle = QIO_ILDGNO;

  QIO_Writer *outfile;
  QIO_Oflag oflag;


  VRB.Flow(cname,fname,"open output file %s \n",filename);

  oflag.serpar=serpar;
  oflag.ildgstyle=ildgstyle;

  oflag.ildgLFN = NULL;

  oflag.mode = QIO_TRUNC;

  outfile = QIO_open_write(record_file, filename, volfmt, &layout, &oflag);

  qio_Output = outfile;

  VRB.FuncEnd(cname,fname);

}





// this was the old write, now using for scalarSource, 12 Sinks...
void qio_writePropagator::write_ScS_12sink(char *outfile, const void *prop, const void *scSource, int volFormat, FP_FORMAT floatFormat)
{

  const char * fname = "write_ScS_12sink(...)";

  VRB.Func(cname,fname);

  VRB.Result(cname,fname,"writing propagator format:\n    scalar source plus 12 sinks\n    to %s \n", outfile);

  int return_val(0);

  // sizepersite -> sizeperSiteSpinColor
  const int sizepersite(12);

  Float *wprop = (Float *)prop;

  Float *wsource = (Float *)scSource;


  QIO_RecordInfo *record_prop;
  record_prop = QIO_create_record_info(0, "", "", 0,0,0,0);

  QIO_RecordInfo *record_source;
  record_source = QIO_create_record_info(0, "", "", 0,0,0,0);

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

  

  CPS_QIO_PROP_FileRecordInfo *filerecordinfo = CPS_QIO_PROP_create_file_record_info("USQCD_DiracFermion_ScalarSource_TwelveSink", xml_info_file);

  CPS_QIO_PROP_encode_file_record_info(record_file, filerecordinfo);

  CPS_QIO_PROP_UserRecordInfo *userrecordinfo_prop = CPS_QIO_PROP_create_user_record_info(0, 0, xml_info_prop);

  CPS_QIO_SOURCE_UserRecordInfo *userrecordinfo_source = CPS_QIO_SOURCE_create_user_record_info(xml_info_source);


  qio_openOutput(outfile, record_file, volFormat);
  

  if(SingleDouble)
    {
      //output in double-precision

      // create the record info
      record_prop = QIO_create_record_info(QIO_FIELD, "USQCD_D3_DiracFermion", "D", 3, 4, sizeof(Float), 2*sizepersite);
      // color=3 (not used), spin=4 (not used), size 2 (complex), count per site

      record_source = QIO_create_record_info(QIO_FIELD, "USQCD_D3_Complex", "D", 0, 0, sizeof(Float), 2);
      // color, spin not used, size 2(complex)

      return_val = 0 ;

      CPS_QIO_SOURCE_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSource, 2*sizeof(Float), sizeof(Float), wsource);


      CPS_QIO_PROP_insert_userrecordinfo_spin( userrecordinfo_prop, 0);

      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 0);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 1);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 2);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      
      CPS_QIO_PROP_insert_userrecordinfo_spin( userrecordinfo_prop, 1);

      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 0);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 1);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 2);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

      CPS_QIO_PROP_insert_userrecordinfo_spin( userrecordinfo_prop, 2);

      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 0);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 1);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 2);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

      CPS_QIO_PROP_insert_userrecordinfo_spin( userrecordinfo_prop, 3);

      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 0);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 1);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 2);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);


    }
  else
    {
      //output in single-precision

      // create the record info
      record_prop = QIO_create_record_info(QIO_FIELD, "USQCD_F3_DiracFermion", "F", 3, 4, sizeof(float), 2*sizepersite);
      // color=3 (not used), spin=4 (not used), size 2 (complex), count per side

      record_source = QIO_create_record_info(QIO_FIELD, "USQCD_F3_Complex", "F", 0, 0, sizeof(float), 2);
      // color, spin not used, size 2(complex)

      //now scalarSource + 12 sinks
      return_val = 0 ;

      CPS_QIO_SOURCE_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getScSourceSingle, 2*sizeof(float), sizeof(float), wsource);

      CPS_QIO_PROP_insert_userrecordinfo_spin( userrecordinfo_prop, 0);

      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 0);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(0,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 1);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(0,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 2);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(0,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
      CPS_QIO_PROP_insert_userrecordinfo_spin( userrecordinfo_prop, 1);

      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 0);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(1,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 1);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(1,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 2);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);      
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(1,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
      CPS_QIO_PROP_insert_userrecordinfo_spin( userrecordinfo_prop, 3);
 
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 0);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(2,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 1);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);      
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(2,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 2);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(2,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      
      CPS_QIO_PROP_insert_userrecordinfo_spin( userrecordinfo_prop, 3);

      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 0);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(3,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 1);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getPropSingle_SpinColor(3,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_PROP_insert_userrecordinfo_color( userrecordinfo_prop, 2);
      CPS_QIO_PROP_encode_user_record_info(record_xml_prop, userrecordinfo_prop);      
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

  CPS_QIO_PROP_destroy_file_record_info(filerecordinfo);
  CPS_QIO_PROP_destroy_user_record_info(userrecordinfo_prop);
  CPS_QIO_SOURCE_destroy_user_record_info(userrecordinfo_source);

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

  const int sizepersite(12);

  int source_size(0);

  Float *wprop = (Float *)prop;

  Float *wsource = (Float *)source;

  QIO_RecordInfo *record_prop;
  record_prop = QIO_create_record_info(0, "", "", 0,0,0,0);

  QIO_RecordInfo *record_source;
  record_source = QIO_create_record_info(0, "", "", 0,0,0,0);

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


  CPS_QIO_PROP_FileRecordInfo *filerecordinfo;


  switch(sType)
    {
    case 1:
      filerecordinfo = CPS_QIO_PROP_create_file_record_info("USQCD_DiracFermion_ScalarSource_Sink_Pairs", xml_info_file);
      break;

    case 2:
      filerecordinfo = CPS_QIO_PROP_create_file_record_info("USQCD_DiracFermion_Source_Sink_Pairs", xml_info_file);
      break;

    default:
      ERR.General(cname,fname,"ERROR QIO: unknown source type %s\n",sType);

    }


  CPS_QIO_PROP_encode_file_record_info(record_file, filerecordinfo);

  CPS_QIO_PROP_PAIRS_UserRecordInfo *userrecordinfo_prop = CPS_QIO_PROP_PAIRS_create_user_record_info(xml_info_prop);

  CPS_QIO_SOURCE_PAIRS_UserRecordInfo *userrecordinfo_source = CPS_QIO_SOURCE_PAIRS_create_user_record_info(0,0,xml_info_source);

  qio_openOutput(outfile, record_file, volFormat);


  if(SingleDouble)
    {
      //output in double-precision

      // create the record info
      record_prop = QIO_create_record_info(QIO_FIELD, "USQCD_D3_DiracFermion", "D", 3, 4, sizeof(Float), 2*sizepersite);
      // color=3 (not used), spin=4 (not used), size 2 (complex), count per site


      if( sType == 1)  // scalar source
	{
	  record_source = QIO_create_record_info(QIO_FIELD, "USQCD_D3_Complex", "D", 0, 0, sizeof(Float), 2);
	  // color, spin not used, size 2(complex)
	  source_size = 2;
	}

      if( sType == 2) // full source
	{
	  record_source = QIO_create_record_info(QIO_FIELD, "USQCD_D3_DiracFermion", "D", 3, 4, sizeof(Float), 2*sizepersite);
	  // color=3 (not used), spin=4 (not used), size 2(complex), count per site
	  source_size = 2*sizepersite;
	}

      return_val = 0 ;



      // now 12 pairs of source, sink, ordered by spin, color.

	  
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_spin( userrecordinfo_source, 0);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 0);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(0,0), source_size*sizeof(Float), sizeof(Float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 1);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(0,1), source_size*sizeof(Float), sizeof(Float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 2);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(0,2), source_size*sizeof(Float), sizeof(Float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_spin( userrecordinfo_source, 1);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 0);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(1,0), source_size*sizeof(Float), sizeof(Float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 1);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(1,1), source_size*sizeof(Float), sizeof(Float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 2);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(1,2), source_size*sizeof(Float), sizeof(Float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_spin( userrecordinfo_source, 2);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 0);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(2,0), source_size*sizeof(Float), sizeof(Float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 1);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(2,1), source_size*sizeof(Float), sizeof(Float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 2);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(2,2), source_size*sizeof(Float), sizeof(Float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_spin( userrecordinfo_source, 3);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 0);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(3,0), source_size*sizeof(Float), sizeof(Float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,0), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 1);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(3,1), source_size*sizeof(Float), sizeof(Float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,1), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 2);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(3,2), source_size*sizeof(Float), sizeof(Float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,2), 2*sizepersite*sizeof(Float), sizeof(Float), wprop);


    }
  else
    {
      //output in single-precision

      // create the record info
      record_prop = QIO_create_record_info(QIO_FIELD, "USQCD_F3_DiracFermion", "F", 3, 4, sizeof(float), 2*sizepersite);
      // color=3 (not used), spin=4 (not used), size 2 (complex), count per site


      if( sType == 1)  // scalar source
	{
	  record_source = QIO_create_record_info(QIO_FIELD, "USQCD_F3_Complex", "F", 0, 0, sizeof(float), 2);
	  // color, spin not used, size 2(complex)
	  source_size = 2;
	}

      if( sType == 2) // full source
	{
	  record_source = QIO_create_record_info(QIO_FIELD, "USQCD_F3_DiracFermion", "F", 3, 4, sizeof(float), 2*sizepersite);
	  // color=3 (not used), spin=4 (not used), size 2(complex), count per site
	  source_size = 2*sizepersite;
	}

      return_val = 0 ;


      // now 12 pairs of source, sink, ordered by spin, color.

	  
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_spin( userrecordinfo_source, 0);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 0);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(0,0), source_size*sizeof(float), sizeof(float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 1);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(0,1), source_size*sizeof(float), sizeof(float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 2);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(0,2), source_size*sizeof(float), sizeof(float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(0,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_spin( userrecordinfo_source, 1);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 0);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(1,0), source_size*sizeof(float), sizeof(float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 1);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(1,1), source_size*sizeof(float), sizeof(float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 2);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(1,2), source_size*sizeof(float), sizeof(float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(1,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_spin( userrecordinfo_source, 2);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 0);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(2,0), source_size*sizeof(float), sizeof(float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 1);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(2,1), source_size*sizeof(float), sizeof(float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 2);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(2,2), source_size*sizeof(float), sizeof(float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(2,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_spin( userrecordinfo_source, 3);

      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 0);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(3,0), source_size*sizeof(float), sizeof(float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,0), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 1);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(3,1), source_size*sizeof(float), sizeof(float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,1), 2*sizepersite*sizeof(float), sizeof(float), wprop);
      CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( userrecordinfo_source, 2);
      CPS_QIO_SOURCE_PAIRS_encode_user_record_info(record_xml_source, userrecordinfo_source);
      return_val += QIO_write( qio_Output, record_source, record_xml_source, qio_getProp_SpinColor(3,2), source_size*sizeof(float), sizeof(float), wsource);
      CPS_QIO_PROP_PAIRS_encode_user_record_info(record_xml_prop, userrecordinfo_prop);
      return_val += QIO_write( qio_Output, record_prop, record_xml_prop, qio_getProp_SpinColor(3,2), 2*sizepersite*sizeof(float), sizeof(float), wprop);


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

  CPS_QIO_PROP_destroy_file_record_info(filerecordinfo);
  CPS_QIO_PROP_PAIRS_destroy_user_record_info(userrecordinfo_prop);
  CPS_QIO_SOURCE_PAIRS_destroy_user_record_info(userrecordinfo_source);

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


CPS_END_NAMESPACE
#endif
