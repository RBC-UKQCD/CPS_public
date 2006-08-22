#ifdef USE_QIO
#include <config.h>
#include <util/qio_writeLattice.h>

CPS_START_NAMESPACE
using namespace std;


// qio-factory functions outside class!!

void qio_getField(char *buf, size_t index, int count, void *arg)
{

#ifdef DEBUG_GetField
  printf("UID: %i, called qio_getField with index: %i, count: %i\n",UniqueID(), index, count);
#endif // DEBUG_GetField

  Matrix *lat = (Matrix *)arg;

  Matrix *mat0 = (Matrix *) buf;
  Matrix *mat1 = (Matrix *) buf + 1;
  Matrix *mat2 = (Matrix *) buf + 2;
  Matrix *mat3 = (Matrix *) buf + 3;

  Matrix *src_mat0 = lat + index * 4 ;
  Matrix *src_mat1 = lat + index * 4 + 1;
  Matrix *src_mat2 = lat + index * 4 + 2;
  Matrix *src_mat3 = lat + index * 4 + 3;

  *mat0 = *src_mat0;
  *mat1 = *src_mat1;
  *mat2 = *src_mat2;
  *mat3 = *src_mat3;

#ifdef DEBUG_GetField
  printf("UID: %i, finished qio_getField (called with index: %i, count: %i)\n",UniqueID(), index, count);
#endif // DEBUG_GetField

}


void qio_getFieldSingle(char *buf, size_t index, int count, void *arg)
{


#ifdef DEBUG_GetField
  printf("UID: %i, called qio_getFieldSingle with index: %i, count: %i\n",UniqueID(), index, count);
#endif // DEBUG_GetField

  float *mat = (float *) buf;

  double *lat = (double *)arg;

  double *src_mat = lat + index * 4 * 9 * 2 ;

  for(int ii(0); ii < 72; ++ii)
    *(mat + ii) = *(src_mat +ii);
  

#ifdef DEBUG_GetField
  printf("UID: %i, finished qio_getFieldSingle (called with index: %i, count: %i)\n",UniqueID(), index, count);
#endif // DEBUG_GetField

}

void qio_getGlobal(char *buf, size_t index, int count, void *arg)
{

#ifdef DEBUG_GlobalWrite
  printf("GetGlobal: UID: %i: called with index: %i, count: %i\n",UniqueID(),index, count);
#endif // DEBUG_GlobalWrite
    

  Float *data = (Float *) arg;
  Float *outp = (Float *) buf;

  for(int ii(0); ii < count; ++ ii)
    *(outp + ii) = *(data + ii);

#ifdef DEBUG_GlobalWrite
  printf("GetGlobal: UID: %i: wrote %f %f to GlobalData\n",UniqueID(),*(outp), *(outp+1));
#endif // DEBUG_GlobalWrite

}


void qio_getGlobalSingle(char *buf, size_t index, int count, void *arg)
{

#ifdef DEBUG_GlobalWrite
  printf("GetGlobalSingle: UID: %i: called with index: %i, count: %i\n",UniqueID(),index, count);
#endif // DEBUG_GlobalWrite
    

  Float *data = (Float *) arg;
  float *outp = (float *) buf;

  for(int ii(0); ii < count; ++ ii)
    *(outp + ii) = *(data + ii);

#ifdef DEBUG_GlobalWrite
  printf("GetGlobalSingle: UID: %i: wrote %f %f to GlobalData\n",UniqueID(),*(outp), *(outp+1));
#endif // DEBUG_GlobalWrite

}


// now start class-functions...


void qio_writeLattice::qio_openOutput(char *filename, char *stringLFN, char *xml_write_file)
{

  char * fname = "qio_openOutput(...)";
  VRB.Func(cname,fname);

  const int volfmt(QIO_VOLFMT);
  const int serpar(QIO_SERPAR);
  const int ildgstyle(QIO_ILDGSTYLE);

  QIO_String *xml_file_out;
  QIO_Writer *outfile;
  QIO_Oflag oflag;


  VRB.Flow(cname,fname,"open output file %s with xml-info: %s\n",filename, xml_write_file);

  xml_file_out = QIO_string_create();
  QIO_string_set(xml_file_out, xml_write_file);

  oflag.serpar=serpar;
  oflag.ildgstyle=ildgstyle;

  if(stringLFN != NULL)
    { 
      oflag.ildgLFN = QIO_string_create();
      QIO_string_set(oflag.ildgLFN, stringLFN);
    }
  else
    oflag.ildgLFN = NULL;

  oflag.mode = QIO_TRUNC;

  outfile = QIO_open_write(xml_file_out, filename, volfmt, &layout, &oflag);

  QIO_string_destroy(xml_file_out);
  
  qio_Output = outfile;

  VRB.FuncEnd(cname,fname);

}



void qio_writeLattice::write(char *outfile, Lattice &lat, FP_FORMAT floatFormat)
{

  const char * fname = "write(...)";

  VRB.Func(cname,fname);

  int return_val(0), return_val_globDat(0);


#ifdef DEBUG_GlobalWrite
  printf(" initializing global-data\n");
#endif

  const int globCount(2);  // number of global-data Floats

  Float globdat[globCount];

  globdat[0] = lat.SumReTrPlaq()/(18*GJP.VolSites()) ;

  Float ltrace(0.0);
  Matrix * lpoint = lat.GaugeField();
  if(GJP.SnodeCoor() == 0) 
    {
      for(int i=0;i<GJP.VolNodeSites()*4 ;i++)
	{
	  ltrace += (lpoint+i)->ReTr();
	}
      glb_sum_five(&ltrace) ;
    }
  else
    glb_sum_five(&ltrace);  // everyone has to participate in global ops

  globdat[1] = ltrace / (4*3*GJP.VolSites());

#ifdef DEBUG_GlobalWrite
  printf("UID: %i: loaded %f %f to GlobalData\n",UniqueID(),globdat[0], globdat[1]);
#endif // DEBUG_GlobalWrite



  QIO_RecordInfo *record;
  QIO_RecordInfo *record_globDat;
  QIO_String *record_xml;
  QIO_String *record_xml_globDat;

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

  
  char dummy[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Dummy QCDML</title>";

  record_xml = QIO_string_create();
  QIO_string_set(record_xml, dummy);

  char dummy2[] = "?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Plaquette, LinkTrace</title>";

  record_xml_globDat = QIO_string_create();
  QIO_string_set(record_xml_globDat, dummy2);

  char xml_outfile[] = "gauge_field";
 
  Matrix * wlat = lat.GaugeField();
  

  qio_openOutput(outfile, NULL, xml_outfile);
  // no idea what ildgLFN member of oflag should be, fill with NULL

  if(SingleDouble)
    {
      //output in double-precision

      // create the record info
      record = QIO_create_record_info(QIO_FIELD, "float_double", "D", 3, 4, 18*sizeof(Float), 4);
      // color=3 (not used), spin=4 (not used), 3*3*2*size (color*color*complex*size), 4 matrices per side

      record_globDat = QIO_create_record_info(QIO_GLOBAL, "float_double", "D",0,0, sizeof(Float),globCount);
      // store Plaquette and LinkTr

      // call write...  
      return_val_globDat = QIO_write( qio_Output, record_globDat, record_xml_globDat, qio_getGlobal, globCount*sizeof(Float), sizeof(Float), &globdat[0]);
      return_val = QIO_write( qio_Output, record, record_xml, qio_getField, 4*3*3*2*sizeof(Float), sizeof(Float), wlat);  

    }
  else
    {
      //output in single-precision

      // create the record info
      record = QIO_create_record_info(QIO_FIELD, "float_single", "F", 3, 4, 18*sizeof(float), 4);
      // color=3 (not used), spin=4 (not used), 3*3*2*size (color*color*complex*size), 4 matrices per side

      record_globDat = QIO_create_record_info(QIO_GLOBAL, "float_single", "F",0,0, sizeof(float),globCount);
      // store Plaquette and LinkTr

      return_val_globDat = QIO_write( qio_Output, record_globDat, record_xml_globDat, qio_getGlobalSingle, globCount*sizeof(float), sizeof(float), &globdat[0]);
      return_val = QIO_write( qio_Output, record, record_xml, qio_getFieldSingle, 4*3*3*2*sizeof(float), sizeof(float), wlat);  
    }

  if ( (return_val == 0) && (return_val_globDat == 0) )
    VRB.Result(cname,fname,"QIO_write successfull...\n");
  else
    ERR.General(cname,fname,"ERROR QIO: QIO_write returned %i %i (global, field)\n",return_val_globDat, return_val);


#ifdef PRINT_checksums
  uint32_t writeCheckA(QIO_get_writer_last_checksuma( output));
  uint32_t writeCheckB(QIO_get_writer_last_checksumb( output));
 
  printf("Checksums: a: %s -- b: %s \n",writeCheckA, writeCheckB);
#endif //PRINT_checksums

  // clean-up

  QIO_destroy_record_info(record);
  QIO_destroy_record_info(record_globDat);

  //different call
  //qio_closeOutput( output);
  qio_closeOutput();

  QIO_string_destroy(record_xml);
  QIO_string_destroy(record_xml_globDat);


  VRB.FuncEnd(cname,fname);

}
CPS_END_NAMESPACE
#endif
