#ifdef USE_QIO
#include <util/qio_readLattice.h>

CPS_START_NAMESPACE
using namespace std;


// qio-factory functions outside class!!


void qio_putField(char *buf, size_t index, int count, void *arg)
{
  #ifdef DEBUG_PutField
  printf("UID: %i, called qio_putField with index: %i, count: %i\n",UniqueID(), index, count);
  #endif // DEBUG_PutField

  Matrix *mat = (Matrix *)arg;
  
  Matrix *mat0 = mat + index * 4;
  Matrix *mat1 = mat + index * 4 + 1;
  Matrix *mat2 = mat + index * 4 + 2;
  Matrix *mat3 = mat + index * 4 + 3;
  
  Matrix *src_mat0 = (Matrix *) buf;
  Matrix *src_mat1 = (Matrix *) buf + 1;
  Matrix *src_mat2 = (Matrix *) buf + 2;
  Matrix *src_mat3 = (Matrix *) buf + 3;
  
  *mat0 = *src_mat0;
  *mat1 = *src_mat1;
  *mat2 = *src_mat2;
  *mat3 = *src_mat3;

  #ifdef DEBUG_PutField
  printf("UID: %i, finished qio_putField (called with index: %i, count: %i)\n",UniqueID(), index, count);
  #endif // DEBUG_PutField
   
}




void qio_putFieldSingle(char *buf, size_t index, int count, void *arg)
{
  #ifdef DEBUG_PutField
  printf("UID: %i, called qio_putFieldSingle with index: %i, count: %i\n",UniqueID(), index, count);
  #endif // DEBUG_PutField


  float *src_mat = (float *) buf;

  double *lat = (double *)arg;

  double *mat = lat + index * 4 * 9 * 2 ;    

  for(int ii(0); ii < 72; ++ii)
    *(mat + ii) = *(src_mat + ii);
  
  #ifdef DEBUG_PutField
  printf("UID: %i, finished qio_putFieldSingle (called with index: %i, count: %i)\n",UniqueID(), index, count);
  #endif // DEBUG_PutField
     
}


void qio_putGlobal(char *buf, size_t index, int count, void *arg)
{

#ifdef DEBUG_GlobalRead
  printf("PutGlobal: UID: %i: called with index: %i, count: %i\n",UniqueID(),index, count);
#endif // DEBUG_GlobalRead
    
  Float *data = (Float *) buf;
  Float *outp = (Float *) arg;

  for(int ii(0); ii < count; ++ ii)
    *(outp + ii) = *(data + ii);

#ifdef DEBUG_GlobalRead
  printf("PutGlobal: UID: %i: wrote %f %f to GlobalData\n",UniqueID(),*(outp), *(outp+1));
#endif // DEBUG_GlobalRead

}

void qio_putGlobalSingle(char *buf, size_t index, int count, void *arg)
{

#ifdef DEBUG_GlobalRead
  printf("PutGlobalSingle: UID: %i: called with index: %i, count: %i\n",UniqueID(),index, count);
#endif // DEBUG_GlobalRead
    
  float *data = (float *) buf;
  Float *outp = (Float *) arg;

  for(int ii(0); ii < count; ++ ii)
    *(outp + ii) = *(data + ii);

#ifdef DEBUG_GlobalRead
  printf("PutGlobal: UID: %i: wrote %f %f to GlobalData\n",UniqueID(),*(outp), *(outp+1));
#endif // DEBUG_GlobalRead

}


// now start class-functions...


void qio_readLattice::qio_openInput(char *filename, int volFormat)
{

  char * fname = "qio_openInput(...)";
  VRB.Func(cname,fname);

  #ifdef DEBUG_openInput
  printf("called qio_openInput with filename: %s\n",filename);
  #endif // DEBUG_openInput

  const int serpar(QIO_SERPAR);


  QIO_String *xml_file_in;
  QIO_Reader *infile;
  QIO_Iflag iflag;

  xml_file_in = QIO_string_create();

  iflag.serpar = serpar;
  iflag.volfmt = volFormat;
  
  VRB.Flow(cname,fname,"open input file %s\n",filename);

  infile = QIO_open_read(xml_file_in, filename, &layout, &iflag);

  #ifdef DEBUG_openInput
  printf("done QIO_open_read\n");
  #endif // DEBUG_openInput





  QIO_string_destroy(xml_file_in);

  #ifdef DEBUG_openInput
  printf("finishing qio_openInput, filename: %s\n",filename);
  #endif // DEBUG_openInput

  qio_Input = infile;

  VRB.FuncEnd(cname,fname);

} 




void qio_readLattice::read(char *infile, Lattice &lat, int volFormat)
{

  char * fname="read(...)";
  VRB.Func(cname,fname);


#ifdef DEBUG_ReadField
  printf("UID: %i, qio_readField called with filename: %s\n",UniqueID(),infile);
#endif // DEBUG_ReadField

  int return_val(0); 

  #ifdef DEBUG_GlobalRead
  printf(" initializing global-data\n");
  #endif



  
  QIO_RecordInfo *record;
  QIO_String *record_xml;
  CPS_QIO_UserRecordInfo *qio_UserRecordInfo;

  record_xml = QIO_string_create();


  qio_setLayout();

  // forces reading ???
  layout.latdim = 0;
  

  // create the record info
   record         = QIO_create_record_info(0, "", "", 0, 0, 0, 0);
   qio_UserRecordInfo = CPS_QIO_create_user_record_info("", "", "");

  
   //input = qio_openInput(infile, &layout);
   qio_openInput(infile, volFormat);


  #ifdef DEBUG_ReadField
   printf("UID: %i, qio_readField intermediate with filename: %s\n",UniqueID(),infile);
  #endif // DEBUG_ReadField

  #ifdef DO_readDimSize

  int readLatdim(QIO_get_reader_latdim( qio_Input));
  int *readLatsize =   QIO_get_reader_latsize(qio_Input);

  #ifdef DEBUG_ReadField
   printf("UID: %i, qio_readField intermediate with filenames: %s\n",UniqueID(),infile);
  #endif // DEBUG_ReadField

  VRB.Result(cname,fname,"read lattice-dimension: %i\n",readLatdim);
  for(int ii(0); ii < readLatdim; ++ii) VRB.Result(cname,fname,"  size[%i]: %i\n",ii,readLatsize[ii]);
  //VRB.Result(cname,fname," \n");

  #endif //DO_readDimSize

  #ifdef PRINT_checksums
  uint32_t readCheckA(QIO_get_reader_last_checksuma( input));
  uint32_t readCheckB(QIO_get_reader_last_checksumb( input));
  
  printf("Checksums: a: %s -- b: %s \n",readCheckA, readCheckB);
  #endif // PRINT_checksums

  #ifdef DEBUG_ReadField
   printf("UID: %i, qio_readField intermediate with filenames: %s\n",UniqueID(),infile);
  #endif // DEBUG_ReadField



  #ifdef DO_recordInfo
    return_val = QIO_read_record_info( qio_Input, record, record_xml);

    if( return_val == 0)
      VRB.Flow(cname,fname,"QIO_read_record_info successfull...\n");
    else
      ERR.General(cname,fname,"ERROR QIO: QIO_read_record_info returned %i\n",return_val);

    #ifdef DEBUG_ReadField
    printf("UID: %i, qio_readField intermediate with filenames: %s\n",UniqueID(),infile);
    #endif // DEBUG_ReadField


    int readGlobaldata(QIO_get_globaldata(record));
    char *readDatatype(QIO_get_datatype(record));
    char *readPrecision(QIO_get_precision(record));
    int readColors(QIO_get_colors(record));
    int readSpins(QIO_get_spins(record));
    int readTypesize(QIO_get_typesize(record));
    int readDatacount(QIO_get_datacount(record));
    char *readRecorddate(QIO_get_record_date(record));

    VRB.Result(cname,fname,"read values: \n  globaldata: %i\n  datatype: %s\n  precision: %s\n  colors: %i\n  spins: %i\n  typesize: %i\n  datacount: %i\n  date: %s\n",readGlobaldata, readDatatype, readPrecision, readColors, readSpins, readTypesize, readDatacount, readRecorddate);



   #else

    // we need the precision

   return_val = QIO_read_record_info( qio_Input, record, record_xml);

    if( return_val == 0)
      VRB.Flow(cname,fname,"QIO_read_record_info successfull...\n");
    else
      ERR.General(cname,fname,"ERROR QIO: QIO_read_record_info returned %i\n",return_val);


    char *readPrecision(QIO_get_precision(record));

   #endif //DO_recordInfo

   // cannot read directly into lat
   Matrix * rlat = lat.GaugeField();


   switch (*readPrecision)
     {
     case 'D':

     #ifdef DEBUG_ReadField
       printf("UID: %i, qio_readField calling QIO_read Double, intermediate with filenames: %s\n",UniqueID(),infile);
     #endif // DEBUG_ReadField

       return_val = QIO_read(qio_Input, record, record_xml, qio_putField, 4*3*3*2*sizeof(Float), sizeof(Float), rlat);
  

       if ( return_val == 0)
	 VRB.Result(cname,fname,"QIO_read (D) successfull...\n");
       else
	 ERR.General(cname,fname,"ERROR QIO: QIO_read (D) returned %i\n",return_val);

       break;

     case 'F':

       #ifdef DEBUG_ReadField
       printf("UID: %i, qio_readField calling QIO_read Single, intermediate with filenames: %s\n",UniqueID(),infile);
        #endif // DEBUG_ReadField

       return_val = QIO_read(qio_Input, record, record_xml, qio_putFieldSingle, 4*3*3*2*sizeof(float), sizeof(float), rlat);
  
       if ( return_val == 0)
	 VRB.Result(cname,fname,"QIO_read (F) successfull...\n");
       else
	 ERR.General(cname,fname,"ERROR QIO: QIO_read (F) returned %i\n",return_val);

       break;

     default:

       ERR.General(cname,fname,"ERROR: unrecognized precision: %s\n",*readPrecision);
       exit(-1);

     }


   #ifdef DEBUG_ReadField
   printf("UID: %i, qio_readField intermediate with filenames: %s\n",UniqueID(),infile);
   #endif // DEBUG_ReadField



   CPS_QIO_decode_user_record_info(qio_UserRecordInfo, record_xml);
  
   
  
  // now check the plaquette, linkTrace

  char *plaq_str, *linktr_str, *info_str;
  Float readPlaq, readLink;

   
  Float measPlaq(lat.SumReTrPlaq()/(18*GJP.VolSites())), measLink(0.0);
  Float diffPlaq, diffLink;

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

  measLink = ltrace / (4*3*GJP.VolSites());






  if( CPS_QIO_defined_plaq(qio_UserRecordInfo) )
    {
      plaq_str = CPS_QIO_get_plaq(qio_UserRecordInfo);
      readPlaq = atof(plaq_str);
    }
  else
    {
      plaq_str = "UNDEFINED";
      readPlaq = 0.0;
    }
    
  if( CPS_QIO_defined_linktr(qio_UserRecordInfo) )
    {
      linktr_str = CPS_QIO_get_linktr(qio_UserRecordInfo);
      readLink = atof(linktr_str);
    }
  else
    {
      linktr_str = "UNDEFINED";
      readLink = 0.0;
    }

  if( CPS_QIO_defined_info(qio_UserRecordInfo) )
    {
      info_str = CPS_QIO_get_info(qio_UserRecordInfo);
    }
  else
    {
      info_str = "UNDEFINED";
    }
  


  #ifdef DEBUG_ReadField
  VRB.Result(cname,fname," Plaquette read: %s (%12.12g), LinkTr read: %s (%12.12g)\n info-str:\n%s\n", plaq_str,readPlaq,linktr_str,readLink, info_str);
  #endif // DEBUG_ReadField



  diffPlaq = fabs(readPlaq - measPlaq);
  diffLink = fabs(readLink - measLink);

  VRB.Result(cname,fname," Plaquette: measured %12.12g <-> read %12.12g; delta: %12.12g\n LinkTrace: measured %12.12g <-> read %12.12g; delta %12.12g\n",
	 measPlaq, readPlaq, diffPlaq, 
	 measLink, readLink, diffLink);



  VRB.Result(cname, fname," additional information:\n%s", info_str);

  // tolerance depend on read-precision???

  if ( 
      ( ( fabs(diffLink/measLink) > Float(TOLERANCE) ) && ( strcmp(linktr_str, "UNDEFINED")!=0 ) ) 
      || 
      ( ( fabs(diffPlaq/measPlaq) > Float(TOLERANCE) ) && ( strcmp(plaq_str  , "UNDEFINED")!=0 ) ) 
      )
    {
      ERR.General(cname,fname,"ERROR: plaquette or linkTrace mismatch!!!\n");
      //      ERR.General(cname,fname,"%8.4g %8.4g %8.4g \n", (fabs(diffPlaq/measPlaq)), ( fabs(diffLink/measLink)), Float(TOLERANCE) );

      exit(-1);
    }
 
  


  // clean-up
  QIO_destroy_record_info(record);
  qio_closeInput();
  QIO_string_destroy(record_xml);
  CPS_QIO_destroy_user_record_info(qio_UserRecordInfo);




#ifdef DEBUG_ReadField
  printf("UID: %i, qio_readField finished with filenames: %s\n",UniqueID(),infile);
#endif // DEBUG_ReadField
  
  VRB.FuncEnd(cname,fname);


}

CPS_END_NAMESPACE
#endif
