#ifdef USE_QIO
#include <util/qio_readPropagator.h>

CPS_START_NAMESPACE
using namespace std;


// qio-factory functions outside class!!


#define QIO_PUTPROP_BASIC(SOURCE_SPIN, SOURCE_COLOR)\
\
{\
 Float *prop = (Float *)arg; \
\
 Float *array = prop + index*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX*count + SOURCE_SPIN*QIO_PROP_COLOR_MAX*count + SOURCE_COLOR*count; \
\
 Float *src_array = (Float *)buf;\
\
 for(int ii(0); ii < count; ++ii) \
  *(array + ii) = *(src_array + ii); \
\
}


#define QIO_PUTPROP_SINGLE_BASIC(SOURCE_SPIN, SOURCE_COLOR)\
\
{\
 float *src_array = (float *)buf; \
\
 double *prop = (double *)arg; \
\
 double *array = prop + index*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX*count + SOURCE_SPIN*QIO_PROP_COLOR_MAX*count + SOURCE_COLOR*count; \
\
 for(int ii(0); ii < count; ++ii) \
    *(array +ii) = *(src_array + ii); \
\
}

//for use in function call...
#define qio_putProp_SpinColor(SPIN, COLOR) qio_putProp_Spin_ ## SPIN ## _Color_ ## COLOR

void qio_putProp_SpinColor(0,0) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_BASIC(0,0)

void qio_putProp_SpinColor(0,1) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_BASIC(0,1)

void qio_putProp_SpinColor(0,2) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_BASIC(0,2)

void qio_putProp_SpinColor(1,0) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_BASIC(1,0)

void qio_putProp_SpinColor(1,1) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_BASIC(1,1)

void qio_putProp_SpinColor(1,2) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_BASIC(1,2)

void qio_putProp_SpinColor(2,0) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_BASIC(2,0)

void qio_putProp_SpinColor(2,1) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_BASIC(2,1)

void qio_putProp_SpinColor(2,2) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_BASIC(2,2)

void qio_putProp_SpinColor(3,0) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_BASIC(3,0)

void qio_putProp_SpinColor(3,1) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_BASIC(3,1)

void qio_putProp_SpinColor(3,2) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_BASIC(3,2)


  // single prec.

#define qio_putPropSingle_SpinColor(SPIN, COLOR) qio_putPropSingle_Spin_ ## SPIN ## _Color_ ## COLOR

void qio_putPropSingle_SpinColor(0,0) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_SINGLE_BASIC(0,0)

void qio_putPropSingle_SpinColor(0,1) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_SINGLE_BASIC(0,1)

void qio_putPropSingle_SpinColor(0,2) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_SINGLE_BASIC(0,2)

void qio_putPropSingle_SpinColor(1,0) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_SINGLE_BASIC(1,0)

void qio_putPropSingle_SpinColor(1,1) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_SINGLE_BASIC(1,1)

void qio_putPropSingle_SpinColor(1,2) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_SINGLE_BASIC(1,2)

void qio_putPropSingle_SpinColor(2,0) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_SINGLE_BASIC(2,0)

void qio_putPropSingle_SpinColor(2,1) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_SINGLE_BASIC(2,1)

void qio_putPropSingle_SpinColor(2,2) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_SINGLE_BASIC(2,2)

void qio_putPropSingle_SpinColor(3,0) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_SINGLE_BASIC(3,0)

void qio_putPropSingle_SpinColor(3,1) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_SINGLE_BASIC(3,1)

void qio_putPropSingle_SpinColor(3,2) (char *buf, size_t index, int count, void *arg)
  QIO_PUTPROP_SINGLE_BASIC(3,2)


 

// old functions for whole propagator in one record

void qio_putProp(char *buf, size_t index, int count, void *arg)
{


  #ifdef DEBUG_PutProp
  printf("UID: %i, called qio_putProp with index: %i, count: %i\n",UniqueID(), index, count);
  #endif // DEBUG_PutProp

  Float *prop = (Float *) arg;  

  Float *array = prop + index * count;

  Float *src_array = (Float *) buf;

  for(int ii(0); ii < count; ++ii)
    {
      #ifdef DEBUG_PutProp
      printf("UID: %i, qio_putProp value at %i before: %f\n",UniqueID(), ii, *(array +ii));
      #endif // DEBUG_PutProp
      *(array + ii) = *(src_array + ii);
      #ifdef DEBUG_PutProp
      printf("UID: %i, qio_putProp value at %i after : %f\n",UniqueID(), ii, *(array +ii));
      #endif // DEBUG_PutProp
    }

  #ifdef DEBUG_PutProp
  printf("UID: %i, finished qio_putProp (called with index: %i, count: %i)\n",UniqueID(), index, count);
  #endif // DEBUG_PutProp
   
  

}




void qio_putPropSingle(char *buf, size_t index, int count, void *arg)
{
  
  #ifdef DEBUG_PutProp
  printf("UID: %i, called qio_putPropSingle with index: %i, count: %i\n",UniqueID(), index, count);
  #endif // DEBUG_PutProp

  float *src_array = (float *)buf;

  double *prop = (double *)arg;

  double *array = prop + index*count;

  for(int ii(0); ii < count; ++ii)
    *(array +ii) = *(src_array + ii);

  #ifdef DEBUG_PutProp
  printf("UID: %i, finished qio_putPropSingle (called with index: %i, count: %i)\n",UniqueID(), index, count);
  #endif // DEBUG_PutProp
    
  
 
}


// calls for ScalarSource

void qio_putScSource(char *buf, size_t index, int count, void *arg)
{


  #ifdef DEBUG_PutProp
  printf("UID: %i, called qio_putScSource with index: %i, count: %i\n",UniqueID(), index, count);
  #endif // DEBUG_PutProp

  Float *source = (Float *) arg;  

  Float *array = source + index * count;

  Float *src_array = (Float *) buf;

  for(int ii(0); ii < count; ++ii)
    {
      #ifdef DEBUG_PutProp
      printf("UID: %i, qio_putScSource value at %i before: %f\n",UniqueID(), ii, *(array +ii));
      #endif // DEBUG_PutProp
      *(array + ii) = *(src_array + ii);
      #ifdef DEBUG_PutProp
      printf("UID: %i, qio_putScSource value at %i after : %f\n",UniqueID(), ii, *(array +ii));
      #endif // DEBUG_PutProp
    }

  #ifdef DEBUG_PutProp
  printf("UID: %i, finished qio_putScSource (called with index: %i, count: %i)\n",UniqueID(), index, count);
  #endif // DEBUG_PutProp
   
  

}




void qio_putScSourceSingle(char *buf, size_t index, int count, void *arg)
{
  
  #ifdef DEBUG_PutProp
  printf("UID: %i, called qio_putScSourceSingle with index: %i, count: %i\n",UniqueID(), index, count);
  #endif // DEBUG_PutProp

  float *src_array = (float *)buf;

  double *source = (double *)arg;

  double *array = source + index*count;

  for(int ii(0); ii < count; ++ii)
    *(array +ii) = *(src_array + ii);

  #ifdef DEBUG_PutProp
  printf("UID: %i, finished qio_putScSourceSingle (called with index: %i, count: %i)\n",UniqueID(), index, count);
  #endif // DEBUG_PutProp
    
  
 
}




// now start class-functions...


void qio_readPropagator::qio_openInput(char *filename, QIO_String *xml_file_in, int volFormat)
{

  char * fname = "qio_openInput(...)";
  VRB.Func(cname,fname);

  #ifdef DEBUG_openInput
  printf("called qio_openInput with filename: %s\n",filename);
  #endif // DEBUG_openInput

  const int serpar(QIO_SERPAR);

  QIO_Reader *infile;
  QIO_Iflag iflag;

  iflag.serpar = serpar;
  iflag.volfmt = volFormat;
  
  VRB.Flow(cname,fname,"open input file %s\n",filename);

  infile = QIO_open_read(xml_file_in, filename, &layout, &iflag);

  #ifdef DEBUG_openInput
  printf("done QIO_open_read\n");
  #endif // DEBUG_openInput


  #ifdef DEBUG_openInput
  printf("finishing qio_openInput, filename: %s\n",filename);
  #endif // DEBUG_openInput

  qio_Input = infile;

  VRB.FuncEnd(cname,fname);

} 


void qio_readPropagator::read(char *infile, void *prop, void *source, int max_Float_prop, int max_Float_source, int volFormat)
{

  char * fname="read(...)";
  VRB.Func(cname,fname);

  VRB.Result(cname,fname," ... not ready yet, doing nothing.... ;-)\n");

  // here we want to have a general read function, which deals with everything....

 
#ifdef DEBUG_ReadPropagator
  printf("UID: %i, qio_readPropagator called with filename: %s\n",UniqueID(),infile);
#endif // DEBUG_ReadPropagator

  VRB.FuncEnd(cname,fname);

}


void qio_readPropagator::read_12pairs(char *infile, void *prop, void *source, const QIO_PROP_SOURCE_TYPES sType, int volFormat)
{

  char * fname="read_12pairs(...)";
  VRB.Func(cname,fname);


  VRB.Result(cname,fname,"trying to read 12 pairs of %s source and sink\n    from %s \n",
	     ( sType==QIO_FULL_SOURCE ? "full" : ( sType==QIO_SCALAR_SOURCE ? "scalar" : "unknown") ),
	     infile );

  int return_val(0);

  const int sizepersite(12);

  int readColor, readSpin;

  QIO_RecordInfo *record_prop;
  QIO_String *record_xml_prop;

  QIO_RecordInfo *record_source;
  QIO_String *record_xml_source;

  QIO_String *record_file;

  record_xml_prop = QIO_string_create();
  record_xml_source = QIO_string_create();
  record_file = QIO_string_create();
  

  qio_setLayout();

  // create the record info
  record_prop         = QIO_create_record_info(0, "", "", 0, 0, 0, 0);
  record_source       = QIO_create_record_info(0, "", "", 0, 0, 0, 0);


  CPS_QIO_PROP_FileRecordInfo *qio_FileRecordInfo;

  qio_FileRecordInfo = CPS_QIO_PROP_create_file_record_info("", "");
  
  qio_openInput(infile, record_file, volFormat);

  CPS_QIO_PROP_decode_file_record_info( qio_FileRecordInfo, record_file);
  

  char *type = CPS_QIO_PROP_file_get_type(qio_FileRecordInfo);
  char *file_info = CPS_QIO_PROP_file_get_info(qio_FileRecordInfo); 


  char detectType('U');  // U: UNKNOWN, A: scalarSource+12sinks, B: SourceSinkPairs, C: ScalarSourceSinkPairs

  if(strcmp(type, "USQCD_DiracFermion_ScalarSource_TwelveSink") == 0)
    detectType = 'A';

  if(strcmp(type, "USQCD_DiracFermion_Source_Sink_Pairs") == 0)
    detectType = 'B';

  if(strcmp(type, "USQCD_DiracFermion_ScalarSource_Sink_Pairs") == 0)
    detectType = 'C';


  if( detectType == 'U' )
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: don't know how to handle type: %s\n", type);
      return;
    } 



  #ifdef DEBUG_ReadPropagator
   printf("UID: %i, qio_readField intermediate with filename: %s\n",UniqueID(),infile);
  #endif // DEBUG_ReadPropagator



   Float *rprop = (Float *)prop;

   Float *rsource = (Float *)source;


   char switchPrecision('U');

   char *readPrecision;

   int source_size(0);

   switch (detectType)
     {

     case 'A':
       // ScalarSource + 12 sinks....

       ERR.General(cname,fname,"ERROR: trying to read ScS_12sink format in pairs\n");
       exit(-1);
	
       break;


     case 'B':
       // pairs with full source

       if(sType == QIO_FULL_SOURCE)
	 source_size = 12;
       else
	 {
	   ERR.General(cname,fname,"ERROR: wrong sourceType (%s) for reading full source pairs\n", sType);
	   exit(-1);
	 }
	   


       break;

     case 'C':
       // pairs with scalar source

       if(sType == QIO_SCALAR_SOURCE)
	 source_size =1 ;
       else
	 {
	   ERR.General(cname,fname,"ERROR: wrong sourceType (%s) for reading scalar source pairs\n", sType);
	   exit(-1);
	 }

       break;


     default:

       ERR.General(cname,fname,"ERROR: unrecognized format-type...\n");
       exit(-1);

     

     }


      #ifdef DEBUG_ReadPropagator
   printf("UID: %i, qio_readField intermediate with filenames: %s\n",UniqueID(),infile);
   #endif // DEBUG_ReadPropagator


   // now do the reading

   for(int ii(0); ii < 12; ++ii)
     {

       // read next source record -> spin, color, precision

       return_val += qio_readNextSourcePairRecord(2*source_size, rsource, readSpin, readColor);

       
       // read next sink record (pass spin, color) -> precision
       
       return_val += qio_readNextPropPairRecord(2*sizepersite, rprop, readSpin, readColor);


     }



   if ( return_val == 0 )
     VRB.Result(cname,fname,"QIO_read  successfull...\n");
   else
     ERR.General(cname,fname,"ERROR QIO: QIO_read   returned %i \n",return_val );


   // communicate T-chunks in S-direction
   if( GJP.Snodes() > 1)
     { qio_communicateTchunks(rprop, 12*sizepersite); qio_communicateTchunks(rsource, 12*source_size);}


   readProps = 12;
   readSources = 12;
   readSourceType = sType;


  // clean-up
   QIO_destroy_record_info(record_prop);
   QIO_destroy_record_info(record_source);
   
   qio_closeInput();
   QIO_string_destroy(record_xml_prop);
   QIO_string_destroy(record_xml_source);


   QIO_string_destroy(record_file);


   CPS_QIO_PROP_destroy_file_record_info(qio_FileRecordInfo);


#ifdef DEBUG_ReadPropagator
  printf("UID: %i, qio_readField finished with filenames: %s\n",UniqueID(),infile);
#endif // DEBUG_ReadPropagator
  
  VRB.FuncEnd(cname,fname);




} 
  


void qio_readPropagator::read_ScS_12sink(char *infile, void *prop, void *source, int volFormat)
{

  char * fname="read_ScS_12sink(...)";
  VRB.Func(cname,fname);


  VRB.Result(cname,fname,"trying to read scalar source and 12 sinks\n    from %s \n", infile );
  

#ifdef DEBUG_ReadPropagator
  printf("UID: %i, qio_readPropagator called with filename: %s\n",UniqueID(),infile);
#endif // DEBUG_ReadPropagator


  const int sizepersite(12);

  int return_val(0); 

  QIO_RecordInfo *record_prop;
  QIO_String *record_xml_prop;

  QIO_RecordInfo *record_source;
  QIO_String *record_xml_source;

  QIO_String *record_file;

  record_xml_prop = QIO_string_create();
  record_xml_source = QIO_string_create();
  record_file = QIO_string_create();
  

  qio_setLayout();

  // create the record info
  record_prop         = QIO_create_record_info(0, "", "", 0, 0, 0, 0);
  record_source       = QIO_create_record_info(0, "", "", 0, 0, 0, 0);

  CPS_QIO_PROP_FileRecordInfo *qio_FileRecordInfo;

  qio_FileRecordInfo = CPS_QIO_PROP_create_file_record_info("", "");
  
  qio_openInput(infile, record_file, volFormat);

  CPS_QIO_PROP_decode_file_record_info( qio_FileRecordInfo, record_file);
  

  char *type = CPS_QIO_PROP_file_get_type(qio_FileRecordInfo);
  char *file_info = CPS_QIO_PROP_file_get_info(qio_FileRecordInfo); 


  char detectType('U');  // U: UNKNOWN, A: scalarSource+12sinks, B: SourceSinkPairs, C: ScalarSourceSinkPairs

  if(strcmp(type, "USQCD_DiracFermion_ScalarSource_TwelveSink") == 0)
    detectType = 'A';

  if(strcmp(type, "USQCD_DiracFermion_Source_Sink_Pairs") == 0)
    detectType = 'B';

  if(strcmp(type, "USQCD_DiracFermion_ScalarSource_Sink_Pairs") == 0)
    detectType = 'C';


  if( detectType == 'U' )
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: don't know how to handle type: %s\n", type);
      return;
    } 



  #ifdef DEBUG_ReadPropagator
   printf("UID: %i, qio_readField intermediate with filename: %s\n",UniqueID(),infile);
  #endif // DEBUG_ReadPropagator



   Float *rprop = (Float *)prop;

   Float *rsource = (Float *)source;


   char switchPrecision('U');

   char *readPrecision;

   switch (detectType)
     {
     case 'A':
       // ScalarSource + 12 sinks....

#ifdef DEBUG_ReadPropagator
       printf("UID: %i, qio_readField calling intermediate with filenames: %s\n",UniqueID(),infile);
#endif // DEBUG_ReadPropagator


	
       return_val = 0;

       // read source userrecordinfo to detect precision
       return_val += QIO_read_record_info(qio_Input, record_source, record_xml_source);

       readPrecision = QIO_get_precision(record_source);

       if(strcmp( readPrecision, "D" ) == 0)
	 switchPrecision='D';
       if(strcmp( readPrecision, "F" ) == 0)
	 switchPrecision='F';


       switch(switchPrecision)
	 {

	 case 'D':
	   
	   return_val += QIO_read(qio_Input, record_source, record_xml_source, qio_putScSource, 2*sizeof(Float), sizeof(Float), rsource);


	   for(int ii(0); ii < 12; ++ii)
	     return_val += qio_readNextPropagatorRecord(2*sizepersite, rprop);


	   if ( return_val == 0 )
	     VRB.Result(cname,fname,"QIO_read (D) successfull...\n");
	   else
	     ERR.General(cname,fname,"ERROR QIO: QIO_read (D)  returned %i \n",return_val );



	   break;

	 case 'F':

       	   return_val += QIO_read(qio_Input, record_source, record_xml_source, qio_putScSourceSingle, 2*sizeof(float), sizeof(float), rsource);

	   for(int ii(0); ii < 12; ++ii)
	     return_val += qio_readNextPropagatorRecordSingle(2*sizepersite, rprop);


	   if ( return_val == 0 )
	     VRB.Result(cname,fname,"QIO_read (F) successfull...\n");
	   else
	     ERR.General(cname,fname,"ERROR QIO: QIO_read (F)  returned %i \n",return_val );


	   
	   break;

	 default:
	   ERR.General(cname, fname,"ERROR QIO-Prop: don't know how to handle precision: %s\n", readPrecision);
	   return;
	   

	 }
	   


       readProps = 12;
       readSources = 1;
       readSourceType = QIO_SCALAR_SOURCE;

	
       break;


     case 'B':

       ERR.General(cname,fname,"ERROR: trying to read pair format in ScS_12sink\n");
       exit(-1);

       break;

     case 'C':

       ERR.General(cname,fname,"ERROR: trying to read pair format in ScS_12sink\n");
       exit(-1);

       break;

     default:

       ERR.General(cname,fname,"ERROR: unrecognized format-type...\n");
       exit(-1);

     }
   



   #ifdef DEBUG_ReadPropagator
   printf("UID: %i, qio_readField intermediate with filenames: %s\n",UniqueID(),infile);
   #endif // DEBUG_ReadPropagator




   // communicate T-chunks in S-direction
   if( GJP.Snodes() > 1)
     { qio_communicateTchunks(rprop, 12*sizepersite); qio_communicateTchunks(rsource, 1); }
   
  

  // clean-up
   QIO_destroy_record_info(record_prop);
   QIO_destroy_record_info(record_source);
   qio_closeInput();
   QIO_string_destroy(record_xml_prop);
   QIO_string_destroy(record_xml_source);


   QIO_string_destroy(record_file);


   CPS_QIO_PROP_destroy_file_record_info(qio_FileRecordInfo);


#ifdef DEBUG_ReadPropagator
  printf("UID: %i, qio_readField finished with filenames: %s\n",UniqueID(),infile);
#endif // DEBUG_ReadPropagator
  
  VRB.FuncEnd(cname,fname);


}


void qio_readPropagator::qio_communicateTchunks(void *start, const int sizepersite)
{


  char * fname="qio_communicateTchunks()";
  VRB.Func(cname,fname);


#ifdef DEBUG_commTchunk
  printf("UID: %i, qio_communicateTchunks called\n",UniqueID());
#endif // DEBUG_commTchunk


  Float *array = (Float *)start;


  //figure out chunk size
  
  int chunk_size = 2 * sizepersite * GJP.NodeSites(0) * GJP.NodeSites(1) * GJP.NodeSites(2) * GJP.NodeSites(3)/ GJP.Snodes();
  int data_size = chunk_size * sizeof(Float);

#ifdef DEBUG_commTchunk
  printf("array size: %i, Float size: %i, comm-size: %i\n",sizeof(Float), sizeof(Float), COMMS_DATASIZE);
  printf("UID: %i, qio_communicateTchunks chunk_size: %i, data_size: %i\n",UniqueID(), chunk_size, data_size);
#endif // DEBUG_commTchunk



  SCUDirArg* data_send    = new SCUDirArg();
  SCUDirArg* data_receive = new SCUDirArg();
  

  int send_chunk, receive_chunk;

  send_chunk = GJP.SnodeCoor();

  receive_chunk = send_chunk + GJP.Snodes() - 1;
  receive_chunk %= GJP.Snodes();

  //loop over number nodes in S-direction
  for( int nn(0); nn < (GJP.Snodes()-1); ++nn)
    {
      

#ifdef DEBUG_commTchunk
      printf("UID: %i, qio_communicateTchunks: S-coor:%i sends %i, receives %i\n",UniqueID(),GJP.SnodeCoor(),send_chunk, receive_chunk);
#endif // DEBUG_commTchunk



      // send address, receive address
      
      Float *send_array = array + send_chunk*chunk_size;
      Float *receive_array = array + receive_chunk*chunk_size;
      
      // sending up, receiving down 

#ifdef DEBUG_commTchunk
      printf("UID: %i, qio_communicateTchunks: S-coor:%i start send/receive at %i/%i...\n",UniqueID(),GJP.SnodeCoor(), send_array, receive_array);
#endif // DEBUG_commTchunk


      data_send->   Init( (void*)send_array,    SCU_SP, SCU_SEND, data_size );
      data_receive->Init( (void*)receive_array, SCU_SM, SCU_REC,  data_size );

#ifdef DEBUG_commTchunk
      printf("UID: %i, qio_communicateTchunks: S-coor:%i ...initialized...\n",UniqueID(),GJP.SnodeCoor());

      Float *test1_s= (Float *)send_array;
      Float *test1_r= (Float *)receive_array;
      Float *test2_s= (Float *)send_array + ( 2*sizepersite * GJP.NodeSites(0) * GJP.NodeSites(1) * GJP.NodeSites(2) * GJP.NodeSites(3)/ GJP.Snodes() ) -1;
      Float *test2_r= (Float *)receive_array + ( 2*sizepersite * GJP.NodeSites(0) * GJP.NodeSites(1) * GJP.NodeSites(2) * GJP.NodeSites(3)/ GJP.Snodes() ) -1;
      printf("UID: %i, qio_communicateTchunks: S-coor:%i  send: %f, receive: %f, end: %f %f\n",UniqueID(),GJP.SnodeCoor(), *test1_s, *test1_r, *test2_s, *test2_r);
#endif // DEBUG_commTchunk



      SCUTrans(data_send);
      SCUTrans(data_receive);
      SCUTransComplete();

#ifdef DEBUG_commTchunk
      printf("UID: %i, qio_communicateTchunks: S-coor:%i ...done!\n",UniqueID(),GJP.SnodeCoor());

      printf("UID: %i, qio_communicateTchunks: S-coor:%i  send: %f, receive: %f end: %f %f\n",UniqueID(),GJP.SnodeCoor(), *test1_s, *test1_r, *test2_s, *test2_r);
#endif // DEBUG_commTchunk


      // new send/receive

      send_chunk = receive_chunk;
      
      receive_chunk = send_chunk + GJP.Snodes() - 1;
      receive_chunk %= GJP.Snodes();



    }


#ifdef DEBUG_commTchunk
  printf("UID: %i, qio_communicateTchunks finished\n",UniqueID());
#endif // DEBUG_commTchunk
  
  VRB.FuncEnd(cname,fname);

}


int qio_readPropagator::qio_readNextPropagatorRecord(int datacount, Float *rprop)
{

  char * fname="qio_readNextPropagatorRecord()";
  VRB.Func(cname,fname);

  int return_val(0);

  QIO_RecordInfo *record;
  QIO_String *record_xml;
  
  record         = QIO_create_record_info(0, "", "", 0, 0, 0, 0);
  record_xml = QIO_string_create();

  CPS_QIO_PROP_UserRecordInfo *qio_UserRecordInfo;
  qio_UserRecordInfo = CPS_QIO_PROP_create_user_record_info(0, 0, "");

  return_val += QIO_read_record_info( qio_Input, record, record_xml);

  return_val += CPS_QIO_PROP_decode_user_record_info(qio_UserRecordInfo, record_xml);

  if(return_val != 0)
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: failed reading next propagator record info!");
 
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      CPS_QIO_PROP_destroy_user_record_info(qio_UserRecordInfo);
      return (1);
    }


  int readDatacount=QIO_get_datacount(record);
  int readSpin= CPS_QIO_PROP_user_get_spin(qio_UserRecordInfo);
  int readColor= CPS_QIO_PROP_user_get_color(qio_UserRecordInfo);

  //check datacount

  if(readDatacount != datacount)
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: datacount mismatch: read: %i, expected: %i\n",readDatacount, datacount);
      
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      CPS_QIO_PROP_destroy_user_record_info(qio_UserRecordInfo);
      return (1);
    }
      

  switch(readSpin)
    {
    case 0:
      {
	switch(readColor)
	  {
	  case 0:
	    return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,0), datacount*sizeof(Float), sizeof(Float), rprop);
	    break;
	  case 1:
	    return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,1), datacount*sizeof(Float), sizeof(Float), rprop);
	    break;
	  case 2:
	    return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,2), datacount*sizeof(Float), sizeof(Float), rprop);
	    break;
	  default:
	    ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	    ++return_val;
	  }
	
	break;
		
      case 1:
	{
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,0), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,1), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,2), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	}
	break;
	
      case 2:
	{
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,0), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,1), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,2), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	}
	break;
		
      case 3:
	{
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,0), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,1), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,2), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	}
	break;
	
      default:
	ERR.General(cname, fname, "ERROR QIO: spin out of range: %i\n", readSpin);
	++return_val;
      }
      
    }


  QIO_destroy_record_info(record);
  QIO_string_destroy(record_xml);
  CPS_QIO_PROP_destroy_user_record_info(qio_UserRecordInfo);
  
  return return_val;


}


int qio_readPropagator::qio_readNextPropagatorRecordSingle(int datacount, Float *rprop)
{

  char * fname="qio_readNextPropagatorRecordSingle()";
  VRB.Func(cname,fname);

  int return_val(0);

  QIO_RecordInfo *record;
  QIO_String *record_xml;
  
  record         = QIO_create_record_info(0, "", "", 0, 0, 0, 0);
  record_xml = QIO_string_create();

  CPS_QIO_PROP_UserRecordInfo *qio_UserRecordInfo;
  qio_UserRecordInfo = CPS_QIO_PROP_create_user_record_info(0, 0, "");

  return_val += QIO_read_record_info( qio_Input, record, record_xml);

  return_val += CPS_QIO_PROP_decode_user_record_info(qio_UserRecordInfo, record_xml);

  if(return_val != 0)
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: failed reading next propagator record info!");
 
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      CPS_QIO_PROP_destroy_user_record_info(qio_UserRecordInfo);
      return (1);
    }


  int readDatacount=QIO_get_datacount(record);
  int readSpin= CPS_QIO_PROP_user_get_spin(qio_UserRecordInfo);
  int readColor= CPS_QIO_PROP_user_get_color(qio_UserRecordInfo);

  //check datacount

  if(readDatacount != datacount)
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: datacount mismatch: read: %i, expected: %i\n",readDatacount, datacount);
      
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      CPS_QIO_PROP_destroy_user_record_info(qio_UserRecordInfo);
      return (1);
    }
      

  switch(readSpin)
    {
    case 0:
      {
	switch(readColor)
	  {
	  case 0:
	    return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,0), datacount*sizeof(float), sizeof(float), rprop);
	    break;
	  case 1:
	    return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,1), datacount*sizeof(float), sizeof(float), rprop);
	    break;
	  case 2:
	    return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,2), datacount*sizeof(float), sizeof(float), rprop);
	    break;
	  default:
	    ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	    ++return_val;
	  }
	
	break;
		
      case 1:
	{
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,0), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,1), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,2), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	}
	break;
	
      case 2:
	{
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,0), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,1), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,2), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	}
	break;
		
      case 3:
	{
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,0), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,1), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,2), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	}
	break;
	
      default:
	ERR.General(cname, fname, "ERROR QIO: spin out of range: %i\n", readSpin);
	++return_val;
      }
      
    }


  QIO_destroy_record_info(record);
  QIO_string_destroy(record_xml);
  CPS_QIO_PROP_destroy_user_record_info(qio_UserRecordInfo);
  
  return return_val;


}





int qio_readPropagator::qio_readNextSourcePairRecord(int datacount, Float *rprop, int &readSpin, int &readColor)
{

  char * fname="qio_readNextSourcePairRecord()";
  VRB.Func(cname,fname);

  int return_val(0);

  QIO_RecordInfo *record;
  QIO_String *record_xml;
  
  record         = QIO_create_record_info(0, "", "", 0, 0, 0, 0);
  record_xml = QIO_string_create();

  CPS_QIO_SOURCE_PAIRS_UserRecordInfo *qio_UserRecordInfo;
  qio_UserRecordInfo = CPS_QIO_SOURCE_PAIRS_create_user_record_info(0, 0, "");

  return_val += QIO_read_record_info( qio_Input, record, record_xml);

  return_val += CPS_QIO_SOURCE_PAIRS_decode_user_record_info(qio_UserRecordInfo, record_xml);

  if(return_val != 0)
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: failed reading next source record info!");
 
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      CPS_QIO_SOURCE_PAIRS_destroy_user_record_info(qio_UserRecordInfo);
      return (1);
    }


  int readDatacount=QIO_get_datacount(record);
  readSpin= CPS_QIO_SOURCE_PAIRS_user_get_spin(qio_UserRecordInfo);
  readColor= CPS_QIO_SOURCE_PAIRS_user_get_color(qio_UserRecordInfo);

  char *readPrecision = QIO_get_precision(record);
  //check datacount

  if(readDatacount != datacount)
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: datacount mismatch: read: %i, expected: %i\n",readDatacount, datacount);
      
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      CPS_QIO_SOURCE_PAIRS_destroy_user_record_info(qio_UserRecordInfo);
      return (1);
    }
      


  char switchPrecision('U');

  if(strcmp( readPrecision, "D" ) == 0)
    switchPrecision='D';
  if(strcmp( readPrecision, "F" ) == 0)
    switchPrecision='F';

  switch(switchPrecision)
    {


    case 'D':

      switch(readSpin)
	{
	case 0:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,0), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,1), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,2), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }

	  
	  break;
	  
	case 1:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,0), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,1), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,2), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	  
	  break;
	    
	case 2:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,0), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,1), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,2), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	  
	  break;
	    
	case 3:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,0), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,1), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,2), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	  
	  break;
	  
	default:
	  ERR.General(cname, fname, "ERROR QIO: spin out of range: %i\n", readSpin);
	  ++return_val;
	  
	}
      


      break;
      

    case 'F':

      switch(readSpin)
	{
	case 0:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,0), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,1), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,2), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	  
	  break;
	    
	case 1:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,0), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,1), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,2), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	  
	  break;
	  
	case 2:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,0), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,1), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,2), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	  
	  break;
	  
	case 3:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,0), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,1), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,2), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	  
	  break;
	  
	default:
	  ERR.General(cname, fname, "ERROR QIO: spin out of range: %i\n", readSpin);
	  ++return_val;
	}

    
      break;


    default:
      ERR.General(cname, fname, "ERROR QIO: precision out of range: %s\n", readPrecision);
      ++return_val;
      
      
    }




  QIO_destroy_record_info(record);
  QIO_string_destroy(record_xml);
  CPS_QIO_SOURCE_PAIRS_destroy_user_record_info(qio_UserRecordInfo);
 
  return return_val;
  
  
}





int qio_readPropagator::qio_readNextPropPairRecord(int datacount, Float *rprop, const int readSpin, const int readColor)
{

  char * fname="qio_readNextPropPairRecord()";
  VRB.Func(cname,fname);

  int return_val(0);

  QIO_RecordInfo *record;
  QIO_String *record_xml;
  
  record         = QIO_create_record_info(0, "", "", 0, 0, 0, 0);
  record_xml = QIO_string_create();

  CPS_QIO_PROP_PAIRS_UserRecordInfo *qio_UserRecordInfo;
  qio_UserRecordInfo = CPS_QIO_PROP_PAIRS_create_user_record_info("");

  return_val += QIO_read_record_info( qio_Input, record, record_xml);

  return_val += CPS_QIO_PROP_PAIRS_decode_user_record_info(qio_UserRecordInfo, record_xml);

  if(return_val != 0)
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: failed reading next source record info!");
 
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      CPS_QIO_PROP_PAIRS_destroy_user_record_info(qio_UserRecordInfo);
      return (1);
    }


  int readDatacount=QIO_get_datacount(record);
  
  char *readPrecision = QIO_get_precision(record);
  //check datacount

  if(readDatacount != datacount)
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: datacount mismatch: read: %i, expected: %i\n",readDatacount, datacount);
      
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      CPS_QIO_PROP_PAIRS_destroy_user_record_info(qio_UserRecordInfo);
      return (1);
    }
      


  char switchPrecision('U');

  if(strcmp( readPrecision, "D" ) == 0)
    switchPrecision='D';
  if(strcmp( readPrecision, "F" ) == 0)
    switchPrecision='F';

  switch(switchPrecision)
    {
      
    case 'D':

      switch(readSpin)
	{
	case 0:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,0), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,1), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,2), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	  
	
	  break;
	    
	case 1:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,0), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,1), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,2), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	  
	  break;
	    
	case 2:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,0), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,1), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,2), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	  
	  break;
	    
	case 3:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,0), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,1), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,2), datacount*sizeof(Float), sizeof(Float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	    
	  
	  break;
	
	default:
	  ERR.General(cname, fname, "ERROR QIO: spin out of range: %i\n", readSpin);
	  ++return_val;
	  
	}
	      
	  


      break;
      

    case 'F':

      switch(readSpin)
	{
	case 0:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,0), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,1), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,2), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	
	  
	  break;
	  
	case 1:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,0), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,1), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,2), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	  
	  break;
	    
	case 2:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,0), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,1), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,2), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	  
	  break;
	  
	case 3:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,0), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,1), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,2), datacount*sizeof(float), sizeof(float), rprop);
	      break;
	    default:
	      ERR.General(cname, fname, "ERROR QIO: color out of range: %i\n", readColor);
	      ++return_val;
	    }
	  
	  break;
	  
	default:
	  ERR.General(cname, fname, "ERROR QIO: spin out of range: %i\n", readSpin);
	  ++return_val;
	}

    




      break;
      

    default:
      ERR.General(cname, fname, "ERROR QIO: precision out of range: %s\n", readPrecision);
      ++return_val;
      
    }
      
  




  QIO_destroy_record_info(record);
  QIO_string_destroy(record_xml);
  CPS_QIO_PROP_PAIRS_destroy_user_record_info(qio_UserRecordInfo);
  
  return return_val;


}








CPS_END_NAMESPACE
#endif
