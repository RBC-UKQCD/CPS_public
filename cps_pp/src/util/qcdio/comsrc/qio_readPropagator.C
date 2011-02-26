#ifdef USE_QIO
#include <util/qio_readPropagator.h>



// adapted functions to Qpropw-class
// 
// assuming ordering of propagator object:
// 
// [volume][sink_spin][sink_color][source_spin][source_color][ReIm]
// (volume slowest, ReIm fastest)
//
// ees, 02/05/08


// adjust datacount issue ( 1 instead of 12 or 2, not counting complex numbers...)
//
// ees, 02/29/08

// to do: add general reader...

CPS_START_NAMESPACE
using namespace std;

#define IGNORE_QIO_BAD_XML


// qio-factory functions outside class!!

void qio_putTmpSource(char *buf, size_t index, int count, void *arg)
{
  Float *prop = (Float *)arg; 
  Float *array = prop + index*2*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX; 
  Float *src_array = (Float *)buf;
  for(int ii(0); ii < 2*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX; ++ii) 
    *(array + ii) = *(src_array + ii); 

}

void qio_putTmpSourceSingle(char *buf, size_t index, int count, void *arg)
{
  Float *prop = (double *)arg; 
  Float *array = prop + index*2*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX; 
  float *src_array = (float *)buf;
  for(int ii(0); ii < 2*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX; ++ii) 
    *(array + ii) = *(src_array + ii); 

}


#define QIO_PUTPROP_BASIC(SOURCE_SPIN, SOURCE_COLOR)\
\
{\
 /*printf(" called putprop with count %i\n",count);*/ \
\
 Float *prop = (Float *)arg; \
\
 Float *array = prop + index*2*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX + 2*SOURCE_COLOR + 2*QIO_PROP_COLOR_MAX*SOURCE_SPIN;\
\
 Float *src_array = (Float *)buf;\
\
 for(int ii(0); ii < QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX; ++ii){ \
  *(array + 2*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*ii )  = *(src_array + 2*ii );\
\
  *(array + 2*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*ii + 1 )  = *(src_array + 2*ii+ 1 ); }\
\
}


#define QIO_PUTPROP_SINGLE_BASIC(SOURCE_SPIN, SOURCE_COLOR)\
\
{\
 /*printf("  called putprop single with index %i, s %i c%i\n",index, SOURCE_SPIN, SOURCE_COLOR);*/\
\
 float *src_array = (float *)buf; \
\
 Float *prop = (double *)arg; \
\
 Float *array = prop + index * 2 * QIO_PROP_SPIN_MAX * QIO_PROP_COLOR_MAX * QIO_PROP_SPIN_MAX * QIO_PROP_COLOR_MAX + 2*SOURCE_COLOR + 2*QIO_PROP_COLOR_MAX*SOURCE_SPIN;\
\
 for(int ii(0); ii < QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX; ++ii){ \
  *(array + 2*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*ii )  = *(src_array + 2*ii );\
\
  *(array + 2*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*ii + 1 )  = *(src_array + 2*ii + 1 );\
\
  /*printf("  %i read %f %f\n",ii,*(src_array + 2*ii), *(src_array + 2*ii+1));*/\
 }\
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




// calls for ScalarSource

void qio_putScSource(char *buf, size_t index, int count, void *arg)
{


  #ifdef DEBUG_PutProp
  printf("UID: %i, called qio_putScSource with index: %i, count: %i\n",UniqueID(), index, count);
  #endif // DEBUG_PutProp

  Float *source = (Float *) arg;  

  Float *array = source + index * 2;

  Float *src_array = (Float *) buf;

  for(int ii(0); ii < 2; ++ii)
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

  Float *source = (double *)arg;

  Float *array = source + index*2;

  for(int ii(0); ii < 2; ++ii)
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

  infile = QIO_open_read(xml_file_in, filename, &layout, NULL, &iflag);

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
  
  VRB.Result(cname,fname,"trying to read propagator from %s\n",infile);
  
 
#ifdef DEBUG_ReadPropagator
  printf("UID: %i, qio_readPropagator called with filename: %s\n",UniqueID(),infile);
#endif // DEBUG_ReadPropagator
  
  int req_Float_prop(0), req_Float_source(0);
  
  QIO_String *record_file;
  
  record_file = QIO_string_create();
  
  qio_setLayout();

  QIO_USQCDPropFileInfo *qio_FileRecordInfo;

  qio_FileRecordInfo = QIO_create_usqcd_propfile_info(0, "");
  
  qio_openInput(infile, record_file, volFormat);

  QIO_decode_usqcd_propfile_info( qio_FileRecordInfo, record_file);

  int type = QIO_get_usqcd_propfile_type(qio_FileRecordInfo);
  char *file_info = QIO_get_usqcd_propfile_info(qio_FileRecordInfo); 


  VRB.Result(cname,fname," read file info:\n%s\n",file_info);


  char detectType('U');  // U: UNKNOWN, A: scalarSource+12sinks, B: SourceSinkPairs, C: ScalarSourceSinkPairs

  if( type == QIO_USQCDPROPFILETYPE_C1D12 )
    detectType = 'A';

  if( type == QIO_USQCDPROPFILETYPE_DD_PAIRS )
    detectType = 'B';

  if( type == QIO_USQCDPROPFILETYPE_CD_PAIRS )
    detectType = 'C';

  if( detectType == 'U' )
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: don't know how to handle type: %i\n", type);
      return;
    } 

  qio_closeInput();

  // check, if size fits and call read-function (if supported format)
  switch (detectType){

  case 'A':
    // ScalarSource + 12 sinks....
    VRB.Result(cname,fname," found prop. type, scalar source + 12 sinks\n");

    req_Float_prop = 288*GJP.VolNodeSites();
    req_Float_source = 2*GJP.VolNodeSites();

    if( (max_Float_prop >= req_Float_prop) && (max_Float_source >= req_Float_source) )
      read_ScS_12sink(infile, prop, source, volFormat);
    else 
      ERR.General(cname, fname,"not enough memory allocated to store source/sink\n need:  src: %i snk: %i Floats\n avail: src: %i snk: %i Floats\n NO PROPAGATOR WILL BE READ\n",
		  req_Float_prop, req_Float_source, max_Float_prop, max_Float_source);

    break;

  case 'B':
    // 12 pairs, full Source
    VRB.Result(cname,fname," found prop. type, 12 pairs of full source + sink\n");

    req_Float_prop = 288*GJP.VolNodeSites();
    req_Float_source = 288*GJP.VolNodeSites();

    if( (max_Float_prop >= req_Float_prop) && (max_Float_source >= req_Float_source) )
      read_12pairs(infile, prop, source,QIO_FULL_SOURCE, volFormat);
    else 
      ERR.General(cname, fname,"not enough memory allocated to store source/sink\n need:  src: %i snk: %i Floats\n avail: src: %i snk: %i Floats\n NO PROPAGATOR WILL BE READ\n",
		  req_Float_prop, req_Float_source, max_Float_prop, max_Float_source);



    break;

  case 'C':
    // 12 pairs, scalar Source
    VRB.Result(cname,fname," found prop. type, 12 pairs of scalar source + sink\n");

    req_Float_prop = 288*GJP.VolNodeSites();
    req_Float_source = 24*GJP.VolNodeSites();

    if( (max_Float_prop >= req_Float_prop) && (max_Float_source >= req_Float_source) )
      read_12pairs(infile, prop, source,QIO_SCALAR_SOURCE, volFormat);
    else 
      ERR.General(cname, fname,"not enough memory allocated to store source/sink\n need:  src: %i snk: %i Floats\n avail: src: %i snk: %i Floats\n NO PROPAGATOR WILL BE READ\n",
		  req_Float_prop, req_Float_source, max_Float_prop, max_Float_source);



    break;

  default:

    ERR.General(cname,fname,"ERROR: unrecognized or unsupported  propagator format-type\n");
    exit(-1);

  }

  // readProps, readSources, readSourceType are set in specific reader funcs

  // clean-up
  QIO_string_destroy(record_file);
  QIO_destroy_usqcd_propfile_info(qio_FileRecordInfo);

  VRB.FuncEnd(cname,fname);

}


void qio_readPropagator::read(char *infile, int spin, int color, void *prop, void *source, int max_Float_prop, int max_Float_source, int volFormat)
{

  char * fname="read(...)";
  VRB.Func(cname,fname);
  
  VRB.Result(cname,fname,"trying to read propagator from %s\n",infile);
  
 
#ifdef DEBUG_ReadPropagator
  printf("UID: %i, qio_readPropagator called with filename: %s\n",UniqueID(),infile);
#endif // DEBUG_ReadPropagator
  
  int req_Float_prop(0), req_Float_source(0);
  
  QIO_String *record_file;
  
  record_file = QIO_string_create();
  
  qio_setLayout();

  QIO_USQCDPropFileInfo *qio_FileRecordInfo;

  qio_FileRecordInfo = QIO_create_usqcd_propfile_info(0, "");
  
  qio_openInput(infile, record_file, volFormat);

  QIO_decode_usqcd_propfile_info( qio_FileRecordInfo, record_file);

  int type = QIO_get_usqcd_propfile_type(qio_FileRecordInfo);
  char *file_info = QIO_get_usqcd_propfile_info(qio_FileRecordInfo); 


  VRB.Result(cname,fname," read file info:\n%s\n",file_info);


  char detectType('U');  // U: UNKNOWN, A: scalarSource+12sinks, B: SourceSinkPairs, C: ScalarSourceSinkPairs

  if( type == QIO_USQCDPROPFILETYPE_C1D12 )
    detectType = 'A';

  if( type == QIO_USQCDPROPFILETYPE_DD_PAIRS )
    detectType = 'B';

  if( type == QIO_USQCDPROPFILETYPE_CD_PAIRS )
    detectType = 'C';

  if( detectType == 'U' )
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: don't know how to handle type: %i\n", type);
      return;
    } 

  qio_closeInput();

  // check, if size fits and call read-function (if supported format)
  switch (detectType){

  case 'B':
    // pairs, full Source
    VRB.Result(cname,fname," found prop. type, full source + sink pair\n");

    req_Float_prop = 24*GJP.VolNodeSites();
    req_Float_source = 24*GJP.VolNodeSites();

    if( (max_Float_prop >= req_Float_prop) && (max_Float_source >= req_Float_source) )
      read_pair(infile, spin, color, prop, source,QIO_FULL_SOURCE, volFormat);
    else 
      ERR.General(cname, fname,"not enough memory allocated to store source/sink\n need:  src: %i snk: %i Floats\n avail: src: %i snk: %i Floats\n NO PROPAGATOR WILL BE READ\n",
		  req_Float_prop, req_Float_source, max_Float_prop, max_Float_source);

    break;

  default:

    ERR.General(cname,fname,"ERROR: unrecognized or unsupported  propagator format-type\n");
    exit(-1);

  }

  // readProps, readSources, readSourceType are set in specific reader funcs

  // clean-up
  QIO_string_destroy(record_file);
  QIO_destroy_usqcd_propfile_info(qio_FileRecordInfo);

  VRB.FuncEnd(cname,fname);

}




void qio_readPropagator::read_pair(char *infile, int spin, int color, void *prop, void *source, const QIO_PROP_SOURCE_TYPES sType, int volFormat)
{

  char * fname="read_pairs(...)";
  VRB.Func(cname,fname);


  VRB.Result(cname,fname,"trying to read pair of %s source and sink\n    from %s \n",
	     ( sType==QIO_FULL_SOURCE ? "full" : ( sType==QIO_SCALAR_SOURCE ? "scalar" : "unknown") ),
	     infile );

  int return_val(0);

  const int sizepersite(QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX);

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
  record_prop         = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
  record_source       = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);


  QIO_USQCDPropFileInfo *qio_FileRecordInfo;

  qio_FileRecordInfo = QIO_create_usqcd_propfile_info(0, "");
  
  qio_openInput(infile, record_file, volFormat);

  QIO_decode_usqcd_propfile_info( qio_FileRecordInfo, record_file);
  

  int type = QIO_get_usqcd_propfile_type(qio_FileRecordInfo);
  //UNUSED  char *file_info = QIO_get_usqcd_propfile_info(qio_FileRecordInfo); 


  char detectType('U');  // U: UNKNOWN, A: scalarSource+12sinks, B: SourceSinkPairs, C: ScalarSourceSinkPairs

  //if(strcmp(type, "USQCD_DiracFermion_ScalarSource_TwelveSink") == 0)
  if( type == QIO_USQCDPROPFILETYPE_C1D12 )
    detectType = 'A';

  //if(strcmp(type, "USQCD_DiracFermion_Source_Sink_Pairs") == 0)
  if( type == QIO_USQCDPROPFILETYPE_DD_PAIRS )
    detectType = 'B';

  //if(strcmp(type, "USQCD_DiracFermion_ScalarSource_Sink_Pairs") == 0)
  if( type == QIO_USQCDPROPFILETYPE_CD_PAIRS )
    detectType = 'C';


  if( detectType == 'U' )
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: don't know how to handle type: %i\n", type);
      return;
    } 



  #ifdef DEBUG_ReadPropagator
   printf("UID: %i, qio_readField intermediate with filename: %s\n",UniqueID(),infile);
  #endif // DEBUG_ReadPropagator



   Float *rprop = (Float *)prop;

   Float *rsource = (Float *)source;


   int source_size(0);

   switch (detectType)
     {

     case 'B':
       // pairs with full source

       if(sType == QIO_FULL_SOURCE)
	 source_size = sizepersite;
       else
	 {
	   ERR.General(cname,fname,"ERROR: wrong sourceType (%s) for reading full source pairs\n", sType);
	   exit(-1);
	 }
	   


       break;

     default:

       ERR.General(cname,fname,"ERROR: unrecognized format-type...\n");
       exit(-1);

     }


   // now do the reading

   Float *tmpSrc = (Float*) smalloc(GJP.VolNodeSites()*2*source_size*sizeof(Float));

   // we have to introduce a tmp - source, since spin, color will only determined after reading the prop!
   
   // zero src 
   for(int jj(0); jj<GJP.VolNodeSites()*2*source_size;++jj)
     *(tmpSrc +jj) =0.0;
   
   //return_val += qio_readTmpSourceRecord(source_size, tmpSrc);
   return_val += qio_readTmpSourceRecord(sType, tmpSrc);
   
   // read next sink record (pass spin, color) -> precision
   return_val += qio_readNextPropPairRecord(rprop, spin, color);
   
   // now move tmpSrc to right position (spin, color)
   for(int site(0); site < GJP.VolNodeSites(); ++site){
     
     Float *putSource = rsource + site*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX*2*source_size + 2*spin*QIO_PROP_COLOR_MAX+2*color;
     for(int jj(0); jj < source_size; ++jj){
       *(putSource +2*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*jj) = *(tmpSrc + 2*source_size*site + 2*jj);
       *(putSource +2*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*jj + 1 ) = *(tmpSrc + 2*source_size*site + 2*jj + 1 );
     }
   }
   
   sfree(tmpSrc);

   if ( return_val == 0 )
     VRB.Result(cname,fname,"QIO_read  successfull...\n");
   else
     ERR.General(cname,fname,"ERROR QIO: QIO_read   returned %i \n",return_val );

   // communicate T-chunks in S-direction
   if( GJP.Snodes() > 1)
     { //qio_communicateTchunks(rprop, 12*sizepersite); qio_communicateTchunks(rsource, 12*source_size);
       ERR.General(cname,fname,"read prop does not work for spread Ls");
     }


   readProps = 1;
   readSources = 1;
   readSourceType = sType;
   

   // clean-up
   QIO_destroy_record_info(record_prop);
   QIO_destroy_record_info(record_source);
   
   qio_closeInput();
   QIO_string_destroy(record_xml_prop);
   QIO_string_destroy(record_xml_source);


   QIO_string_destroy(record_file);


   QIO_destroy_usqcd_propfile_info(qio_FileRecordInfo);
   
   
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

  const int sizepersite(QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX);

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
  record_prop         = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
  record_source       = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);


  QIO_USQCDPropFileInfo *qio_FileRecordInfo;

  qio_FileRecordInfo = QIO_create_usqcd_propfile_info(0, "");
  
  qio_openInput(infile, record_file, volFormat);

  QIO_decode_usqcd_propfile_info( qio_FileRecordInfo, record_file);
  

  int type = QIO_get_usqcd_propfile_type(qio_FileRecordInfo);
  //UNUSED  char *file_info = QIO_get_usqcd_propfile_info(qio_FileRecordInfo); 


  char detectType('U');  // U: UNKNOWN, A: scalarSource+12sinks, B: SourceSinkPairs, C: ScalarSourceSinkPairs

  //if(strcmp(type, "USQCD_DiracFermion_ScalarSource_TwelveSink") == 0)
  if( type == QIO_USQCDPROPFILETYPE_C1D12 )
    detectType = 'A';

  //if(strcmp(type, "USQCD_DiracFermion_Source_Sink_Pairs") == 0)
  if( type == QIO_USQCDPROPFILETYPE_DD_PAIRS )
    detectType = 'B';

  //if(strcmp(type, "USQCD_DiracFermion_ScalarSource_Sink_Pairs") == 0)
  if( type == QIO_USQCDPROPFILETYPE_CD_PAIRS )
    detectType = 'C';


  if( detectType == 'U' )
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: don't know how to handle type: %i\n", type);
      return;
    } 



  #ifdef DEBUG_ReadPropagator
   printf("UID: %i, qio_readField intermediate with filename: %s\n",UniqueID(),infile);
  #endif // DEBUG_ReadPropagator



   Float *rprop = (Float *)prop;

   Float *rsource = (Float *)source;


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
	 source_size = sizepersite;
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

   Float *tmpSrc = (Float*) smalloc(GJP.VolNodeSites()*2*source_size*sizeof(Float));

   for(int ii(0); ii < 12; ++ii){

     
     // we have to introduce a tmp - source, since spin, color will only determined after reading the prop!
   
     // zero src 
     for(int jj(0); jj<GJP.VolNodeSites()*2*source_size;++jj)
       *(tmpSrc +jj) =0.0;
  
     //return_val += qio_readTmpSourceRecord(source_size, tmpSrc);
     return_val += qio_readTmpSourceRecord(sType, tmpSrc);

     // read next sink record (pass spin, color) -> precision
     return_val += qio_readNextPropPairRecord(rprop, readSpin, readColor);

#ifdef DEBUG_ReadSpinColor
  printf("finished reading source spin %i color %i \n",readSpin, readColor);
#endif


     // now move tmpSrc to right position (spin, color)
     for(int site(0); site < GJP.VolNodeSites(); ++site){
	
       Float *putSource = rsource + site*QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX*2*source_size + 2*readSpin*QIO_PROP_COLOR_MAX+2*readColor;
       for(int jj(0); jj < source_size; ++jj){
	 *(putSource +2*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*jj) = *(tmpSrc + 2*source_size*site + 2*jj);
	 *(putSource +2*QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX*jj + 1 ) = *(tmpSrc + 2*source_size*site + 2*jj + 1 );
       }
       

     }
       

#ifdef DEBUG_PAIRRECORD
     // only works if S-direction is local

     printf(" checking for source-index %i %i...\n",readSpin, readColor);

     Float err_cnt(0.);

     for(int site(0); site < GJP.VolNodeSites(); ++site){
       for( int mm(0); mm < 4; ++mm)
	 for( int cc(0); cc < 3; ++cc){
	   
	   Float index = ( readColor + 3*readSpin + 12*cc + 36*mm)/1000.; 

	   Float *tmp_r = rprop +     2*readColor + 6*readSpin + 24*cc + 72*mm + 288*site;  
	   Float *tmp_i = rprop + 1 + 2*readColor + 6*readSpin + 24*cc + 72*mm + 288*site;  
	   

#ifdef DEBUG_PAIRRECORD_STD
	   if( (fabs(0.5*(*tmp_r - *tmp_i) - index) > 1e-5)){
#else
	   if( (fabs(*tmp_r - site - index) > 1e-5) || (fabs(*tmp_i - site +index) > 1e-5 ) ){
#endif
	     ++ err_cnt;
	     printf("mismatch propagator at %i %i %i %i %i, read: (%f, %f) expected: (%f, %f)\n",site, mm, cc, readSpin, readColor, *tmp_r, *tmp_i, (site+index), (site-index));
	   }
	 }

       for( int ii(0); ii < source_size; ++ii){

	 Float index = ( readColor + 3*readSpin + 12*ii)/1000.;

	 Float *tmp_r = rsource +     2*readColor + 6*readSpin + 24*ii + 24*source_size*site;  
	 Float *tmp_i = rsource + 1 + 2*readColor + 6*readSpin + 24*ii + 24*source_size*site;  
	   
#ifdef DEBUG_PAIRRECORD_STD
	 if( (fabs(0.5*(*tmp_r - *tmp_i) - index) > 1e-5)){
#else
	 if( (fabs(*tmp_r - site - index) > 1e-5) || (fabs(*tmp_i - site +index) > 1e-5 ) ){
#endif
	   ++ err_cnt;
	   printf("mismatch source at %i %i %i %i, read: (%f, %f) expected: (%f, %f)\n",site, ii, readSpin, readColor, *tmp_r, *tmp_i, (site+index), (site-index));
	 }
	 
       } 

     }

     glb_sum_five(&err_cnt);

      printf(" ...error-cnt: %f\n",err_cnt);

#endif      
 


   }

   sfree(tmpSrc);


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


   QIO_destroy_usqcd_propfile_info(qio_FileRecordInfo);


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


  const int sizepersite(QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX);

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
  record_prop         = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
  record_source       = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);

  QIO_USQCDPropFileInfo *qio_FileRecordInfo;

  qio_FileRecordInfo = QIO_create_usqcd_propfile_info(0, "");
  
  qio_openInput(infile, record_file, volFormat);

  QIO_decode_usqcd_propfile_info( qio_FileRecordInfo, record_file);
  
  int type = QIO_get_usqcd_propfile_type(qio_FileRecordInfo);
  //UNUSED char *file_info = QIO_get_usqcd_propfile_info(qio_FileRecordInfo); 


  char detectType('U');  // U: UNKNOWN, A: scalarSource+12sinks, B: SourceSinkPairs, C: ScalarSourceSinkPairs

  //if(strcmp(type, "USQCD_DiracFermion_ScalarSource_TwelveSink") == 0)
  if( type == QIO_USQCDPROPFILETYPE_C1D12 )
    detectType = 'A';

  //if(strcmp(type, "USQCD_DiracFermion_Source_Sink_Pairs") == 0)
  if( type == QIO_USQCDPROPFILETYPE_DD_PAIRS )
    detectType = 'B';

  //if(strcmp(type, "USQCD_DiracFermion_ScalarSource_Sink_Pairs") == 0)
  if( type == QIO_USQCDPROPFILETYPE_CD_PAIRS )
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


       // in case of a hyperfield source record, set scalar source to zero first

       for(int ii(0); ii < GJP.VolNodeSites()*2; ++ii)
	 *(rsource + ii ) = 0.0;
	
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

#ifdef DEBUG_PAIRRECORD
	   // only works if S-direction is local
	   {
	     printf(" checking for source ...\n");

	     Float err_cnt(0);

	     for(int site(0); site < GJP.VolNodeSites(); ++site){

	       Float *tmp_r = rsource + 2*site;  
	       Float *tmp_i = rsource + 1 + 2*site;  
	   

	       if( (fabs(*tmp_r - site ) > 1e-5) || (fabs(*tmp_i + site ) > 1e-5 ) ){
		 ++ err_cnt;
		 printf("mismatch source at %i, read: (%f, %f) expected: (%i, %i)\n",site, *tmp_r, *tmp_i, (site), (-site));
	       }

	     }
	     glb_sum_five(&err_cnt);
	     
	     printf(" ...error-cnt: %f\n",err_cnt);
	       
	   }
#endif //DEBUG_PAIRRECORD

	   for(int ii(0); ii < 12; ++ii){

	     return_val += qio_readNextPropagatorRecord(rprop);

#ifdef DEBUG_PAIRRECORD

	     printf(" checking for source-index, record  %i...\n",ii);

	     Float err_cnt(0.);
	     
	     for(int site(0); site < GJP.VolNodeSites(); ++site)
	       for( int mm(0); mm < 4; ++mm)
		 for( int cc(0); cc < 3; ++cc){
	   
		   Float index = ( ii + 12*cc + 36*mm)/1000.; 
		   
		   Float *tmp_r = rprop +     2*ii + 24*cc + 72*mm + 288*site;  
		   Float *tmp_i = rprop + 1 + 2*ii + 24*cc + 72*mm + 288*site;  
		   
		   
		   if( (fabs(*tmp_r - site - index) > 1e-5) || (fabs(*tmp_i - site +index) > 1e-5 ) ){
		     ++ err_cnt;
		     printf("mismatch propagator at %i %i %i %i, read: (%f, %f) expected: (%f, %f)\n",site, mm, cc, ii, *tmp_r, *tmp_i, (site+index), (site-index));
		   }
		 }
	     
	     glb_sum_five(&err_cnt);
	     
	     printf(" ...error-cnt: %f\n",err_cnt);
	     
#endif //DEBUG_PAIRRECORD
	   }

	   if ( return_val == 0 )
	     VRB.Result(cname,fname,"QIO_read (D) successfull...\n");
	   else
	     ERR.General(cname,fname,"ERROR QIO: QIO_read (D)  returned %i \n",return_val );



	   break;

	 case 'F':

       	   return_val += QIO_read(qio_Input, record_source, record_xml_source, qio_putScSourceSingle, 2*sizeof(float), sizeof(float), rsource); 

	   for(int ii(0); ii < 12; ++ii)
	     return_val += qio_readNextPropagatorRecordSingle(rprop);


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


   QIO_destroy_usqcd_propfile_info(qio_FileRecordInfo);


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
  
  int chunk_size = GJP.NodeSites(0) * GJP.NodeSites(1) * GJP.NodeSites(2) * GJP.NodeSites(3)/ GJP.Snodes();
  int data_size = 2*sizepersite;

#ifdef DEBUG_commTchunk
  printf("array size: %i, Float size: %i, comm-size: %i\n",sizeof(Float), sizeof(Float), COMMS_DATASIZE);
  printf("UID: %i, qio_communicateTchunks chunk_size: %i, data_size: %i\n",UniqueID(), chunk_size, data_size);
#endif // DEBUG_commTchunk

  int send_chunk, receive_chunk;

  send_chunk = GJP.SnodeCoor();

  receive_chunk = send_chunk + GJP.Snodes() - 1;
  receive_chunk %= GJP.Snodes();

  //loop over number nodes in S-direction
  for( int nn(0); nn < (GJP.Snodes()-1); ++nn){
    
    for( int mm(0); mm < chunk_size; ++mm){

#ifdef DEBUG_commTchunk
      printf("UID: %i, qio_communicateTchunks: S-coor:%i sends %i, receives %i\n",UniqueID(),GJP.SnodeCoor(),send_chunk, receive_chunk);
#endif // DEBUG_commTchunk



      // send address, receive address
      
      IFloat *send_array = (IFloat *) array + send_chunk*chunk_size*data_size + mm*data_size;
      IFloat *receive_array = (IFloat *) array + receive_chunk*chunk_size*data_size + mm*data_size;

      // sending up, receiving down 

#ifdef DEBUG_commTchunk
      printf("UID: %i, qio_communicateTchunks: S-coor:%i start send/receive at %i/%i...\n",UniqueID(),GJP.SnodeCoor(), send_array, receive_array);
#endif // DEBUG_commTchunk



#ifdef DEBUG_commTchunk
      printf("UID: %i, qio_communicateTchunks: S-coor:%i ...initialized...\n",UniqueID(),GJP.SnodeCoor());

      IFloat *test1_s= (IFloat *)send_array;
      IFloat *test1_r= (IFloat *)receive_array;
      IFloat *test2_s= (IFloat *)send_array + 2*sizepersite -1;
      IFloat *test2_r= (IFloat *)receive_array + 2*sizepersite  -1;
      printf("UID: %i, qio_communicateTchunks: S-coor:%i  send: %f, receive: %f, end: %f %f\n",UniqueID(),GJP.SnodeCoor(), *test1_s, *test1_r, *test2_s, *test2_r);
#endif // DEBUG_commTchunk
      

      getMinusData( receive_array, send_array, data_size, 4);

      
#ifdef DEBUG_commTchunk
      printf("UID: %i, qio_communicateTchunks: S-coor:%i ...done!\n",UniqueID(),GJP.SnodeCoor());

      printf("UID: %i, qio_communicateTchunks: S-coor:%i  send: %f, receive: %f end: %f %f\n",UniqueID(),GJP.SnodeCoor(), *test1_s, *test1_r, *test2_s, *test2_r);
#endif // DEBUG_commTchunk

    } // loop mm


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


int qio_readPropagator::qio_readNextPropagatorRecord(Float *rprop)
{

  char * fname="qio_readNextPropagatorRecord()";
  VRB.Func(cname,fname);

  int return_val(0), return_val_info(0), return_val_decode(0);

  int sizepersite(QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX);

  QIO_RecordInfo *record;
  QIO_String *record_xml;
  
  record         = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
  record_xml = QIO_string_create();

  QIO_USQCDPropRecordInfo *qio_UserRecordInfo;
  qio_UserRecordInfo = QIO_create_usqcd_proprecord_sc_info(0, 0, "");

  return_val_info += QIO_read_record_info( qio_Input, record, record_xml);

  return_val_decode += QIO_decode_usqcd_proprecord_info(qio_UserRecordInfo, record_xml);

#ifdef IGNORE_QIO_BAD_XML
  if(return_val_info == QIO_BAD_XML){
    
    VRB.Result(cname,fname,"WARNING QIO: reading info (prop record) returned %d (QIO_BAD_XML), IGNORING...\n",return_val_info);
    return_val_info=0;

  }

   if(return_val_decode == QIO_BAD_XML){
    
    VRB.Result(cname,fname,"WARNING QIO: decoding info (prop record) returned %d (QIO_BAD_XML), IGNORING...\n",return_val_decode);
    return_val_decode=0;

  }
#endif //IGNORE_QIO_BAD_XML



  if(return_val_info != 0 || return_val_decode != 0)
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: failed reading next propagator record info!  read record: %d, decode_usqcd_propsource_info: %d\n", return_val_decode, return_val_info);
 
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      QIO_destroy_usqcd_proprecord_info(qio_UserRecordInfo);
      return (1);
    }


  int readDatacount=QIO_get_datacount(record);
  int readDatasize=QIO_get_typesize(record);
  int readSpin, readColor;

#ifdef SPINCOLORCHECK
  if( (QIO_defined_usqcd_proprecord_spin(qio_UserRecordInfo)) && QIO_defined_usqcd_proprecord_color(qio_UserRecordInfo) ) {
#endif
    readSpin= QIO_get_usqcd_proprecord_spin(qio_UserRecordInfo);
    readColor= QIO_get_usqcd_proprecord_color(qio_UserRecordInfo);
    
    lastSpin = readSpin;
    lastColor = readColor;

#ifdef SPINCOLORCHECK
  }
  else{

    qio_guessSpinColor(readSpin, readColor);

    VRB.Result(cname,fname,"no spin/color defined, assuming normal ordering, next spin %i, color %i\n", readSpin, readColor );
  }
#endif

#ifdef DEBUG_ReadSpinColor
  printf("now reading source spin %i color %i \n",readSpin, readColor);
#endif


  //check datacount

  if( int(readDatacount*readDatasize) != int(sizepersite*2*sizeof(Float)) )
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: datacount/typesize mismatch:\n read: %i of size %i ( = %i),\n expected: 1 of size %i x %i = %i\n",
		  readDatacount, readDatasize, (readDatacount*readDatasize),
		  sizepersite, 2*sizeof(Float), (sizepersite*2*sizeof(Float)));
      
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      QIO_destroy_usqcd_proprecord_info(qio_UserRecordInfo);
      return (1);
    }
      

  switch(readSpin)
    {
    case 0:
      {
	switch(readColor)
	  {
	  case 0:
	    return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,0), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	    break;
	  case 1:
	    return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,1), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	    break;
	  case 2:
	    return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,2), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
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
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,0), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,1), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,2), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
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
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,0), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,1), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,2), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
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
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,0), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,1), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,2), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
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
  QIO_destroy_usqcd_proprecord_info(qio_UserRecordInfo);
  
  return return_val;


}



int qio_readPropagator::qio_readNextPropagatorRecordSingle(Float *rprop)
{

  char * fname="qio_readNextPropagatorRecordSingle()";
  VRB.Func(cname,fname);

  int return_val(0), return_val_info(0), return_val_decode(0);

  int sizepersite(QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX);

  QIO_RecordInfo *record;
  QIO_String *record_xml;
  
  record         = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
  record_xml = QIO_string_create();

  QIO_USQCDPropRecordInfo *qio_UserRecordInfo;
  qio_UserRecordInfo = QIO_create_usqcd_proprecord_sc_info(0, 0, "");

  return_val_info += QIO_read_record_info( qio_Input, record, record_xml);

  return_val_decode += QIO_decode_usqcd_proprecord_info(qio_UserRecordInfo, record_xml);

#ifdef IGNORE_QIO_BAD_XML
  if(return_val_info == QIO_BAD_XML){
    
    VRB.Result(cname,fname,"WARNING QIO: reading info (prop record) returned %d (QIO_BAD_XML), IGNORING...\n",return_val_info);
    return_val_info=0;

  }

   if(return_val_decode == QIO_BAD_XML){
    
    VRB.Result(cname,fname,"WARNING QIO: decoding info (prop record) returned %d (QIO_BAD_XML), IGNORING...\n",return_val_decode);
    return_val_decode=0;

  }
#endif //IGNORE_QIO_BAD_XML


  if(return_val_info != 0 || return_val_decode != 0)
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: failed reading next propagator record info! read record: %d, decode_usqcd_propsource_info: %d\n", return_val_decode, return_val_info);
 
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      QIO_destroy_usqcd_proprecord_info(qio_UserRecordInfo);
      return (1);
    }


  int readDatacount=QIO_get_datacount(record);
  int readDatasize=QIO_get_typesize(record);
  int readSpin, readColor;

#ifdef SPINCOLORCHECK
  if( (QIO_defined_usqcd_proprecord_spin(qio_UserRecordInfo)) && QIO_defined_usqcd_proprecord_color(qio_UserRecordInfo) ) {
#endif
    readSpin= QIO_get_usqcd_proprecord_spin(qio_UserRecordInfo);
    readColor= QIO_get_usqcd_proprecord_color(qio_UserRecordInfo);
    
    lastSpin = readSpin;
    lastColor = readColor;

#ifdef SPINCOLORCHECK
  }
  else{

    qio_guessSpinColor(readSpin, readColor);

    VRB.Result(cname,fname,"no spin/color defined, assuming normal ordering, next spin %i, color %i\n", readSpin, readColor );
  }
#endif

#ifdef DEBUG_ReadSpinColor
  printf("now reading source spin %i color %i \n",readSpin, readColor);
#endif


  //check datacount

  if( int(readDatacount*readDatasize) != int(sizepersite*2*sizeof(float)))
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: datacount/typesize mismatch:\n read: %i of size %i (= %i),\n expected: 1 of size %i x %i = %i\n",
		  readDatacount, readDatasize, (readDatacount*readDatasize),
		  sizepersite, 2*sizeof(float), (sizepersite*2*sizeof(float)) );
      
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      QIO_destroy_usqcd_proprecord_info(qio_UserRecordInfo);
      return (1);
    }
      

  switch(readSpin)
    {
    case 0:
      {
	switch(readColor)
	  {
	  case 0:
	    return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,0), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	    break;
	  case 1:
	    return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,1), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	    break;
	  case 2:
	    return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,2), 2*sizepersite*sizeof(float), sizeof(float), rprop);
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
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,0), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,1), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,2), 2*sizepersite*sizeof(float), sizeof(float), rprop);
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
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,0), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,1), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,2), 2*sizepersite*sizeof(float), sizeof(float), rprop);
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
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,0), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,1), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,2), 2*sizepersite*sizeof(float), sizeof(float), rprop);
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
  QIO_destroy_usqcd_proprecord_info(qio_UserRecordInfo);
  
  return return_val;


}




int qio_readPropagator::qio_readNextPropPairRecord(Float *rprop, int &readSpin, int &readColor)
{

  char * fname="qio_readNextPropPairRecord()";
  VRB.Func(cname,fname);

  int return_val(0), return_val_info(0), return_val_decode(0);

  int sizepersite(QIO_PROP_COLOR_MAX*QIO_PROP_SPIN_MAX);

  QIO_RecordInfo *record;
  QIO_String *record_xml;
  
  record         = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
  record_xml = QIO_string_create();

  QIO_USQCDPropRecordInfo *qio_UserRecordInfo;
  qio_UserRecordInfo = QIO_create_usqcd_proprecord_sc_info(0, 0, "");

  return_val_info += QIO_read_record_info( qio_Input, record, record_xml);

  return_val_decode += QIO_decode_usqcd_proprecord_info(qio_UserRecordInfo, record_xml);

#ifdef IGNORE_QIO_BAD_XML
  if(return_val_info == QIO_BAD_XML){
    
    VRB.Result(cname,fname,"WARNING QIO: reading info (prop record) returned %d (QIO_BAD_XML), IGNORING...\n",return_val_info);
    return_val_info=0;

  }

   if(return_val_decode == QIO_BAD_XML){
    
    VRB.Result(cname,fname,"WARNING QIO: decoding info (prop record) returned %d (QIO_BAD_XML), IGNORING...\n",return_val_decode);
    return_val_decode=0;

  }
#endif //IGNORE_QIO_BAD_XML


  if(return_val_decode != 0 || return_val_info != 0)
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: failed reading next prop record info! read record: %d, decode_usqcd_propsource_info: %d\n", return_val_decode, return_val_info);
 
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      QIO_destroy_usqcd_proprecord_info(qio_UserRecordInfo);
      return (1);
    }

  // defined ??, WORK!
  int readDatacount=QIO_get_datacount(record);
  int readDatasize=QIO_get_typesize(record);


#ifdef SPINCOLORCHECK
  if( (QIO_defined_usqcd_proprecord_spin(qio_UserRecordInfo)) && QIO_defined_usqcd_proprecord_color(qio_UserRecordInfo) ) {
#endif
    readSpin= QIO_get_usqcd_proprecord_spin(qio_UserRecordInfo);
    readColor= QIO_get_usqcd_proprecord_color(qio_UserRecordInfo);
    
    lastSpin = readSpin;
    lastColor = readColor;

#ifdef SPINCOLORCHECK
  }
  else{

    qio_guessSpinColor(readSpin, readColor);

    VRB.Result(cname,fname,"no spin/color defined, assuming normal ordering, next spin %i, color %i\n", readSpin, readColor );
  }
#endif


#ifdef DEBUG_ReadSpinColor
  printf("now reading source spin %i color %i \n",readSpin, readColor);
#endif


  
  char *readPrecision = QIO_get_precision(record);
      


  char switchPrecision('U');
  int basesize(0);

  if(strcmp( readPrecision, "D" ) == 0)
    {switchPrecision='D'; basesize=2*sizeof(Float);} 
  if(strcmp( readPrecision, "F" ) == 0)
    {switchPrecision='F'; basesize=2*sizeof(float);}

  //check datacount
  
  if( (readDatacount*readDatasize) != (sizepersite*basesize))
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: datacount/typesize mismatch:\n read: %i of size %i ( = %i),\n expected: 1 of size %i x %i = %i\n",
		  readDatacount, readDatasize,
		  sizepersite, basesize, (sizepersite*basesize));
      
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      QIO_destroy_usqcd_proprecord_info(qio_UserRecordInfo);
      return (1);
    }



  switch(switchPrecision)
    {
      
    case 'D':

      switch(readSpin)
	{
	case 0:
	  
	  switch(readColor)
	    {
	    case 0:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,0), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,1), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(0,2), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
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
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,0), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,1), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(1,2), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
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
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,0), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,1), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(2,2), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
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
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,0), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,1), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putProp_SpinColor(3,2), 2*sizepersite*sizeof(Float), sizeof(Float), rprop);
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
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,0), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,1), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(0,2), 2*sizepersite*sizeof(float), sizeof(float), rprop);
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
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,0), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,1), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(1,2), 2*sizepersite*sizeof(float), sizeof(float), rprop);
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
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,0), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,1), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(2,2), 2*sizepersite*sizeof(float), sizeof(float), rprop);
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
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,0), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 1:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,1), 2*sizepersite*sizeof(float), sizeof(float), rprop);
	      break;
	    case 2:
	      return_val += QIO_read(qio_Input, record, record_xml, qio_putPropSingle_SpinColor(3,2), 2*sizepersite*sizeof(float), sizeof(float), rprop);
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
  QIO_destroy_usqcd_proprecord_info(qio_UserRecordInfo);
  
  return return_val;


}

 

int qio_readPropagator::qio_readTmpSourceRecord(const QIO_PROP_SOURCE_TYPES sType, Float *rsource)
{

  char * fname="qio_readTmpSourceRecord()";
  VRB.Func(cname,fname);

  int return_val(0), return_val_info(0), return_val_decode(0);

  QIO_RecordInfo *record;
  QIO_String *record_xml;
  
  record         = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
  record_xml = QIO_string_create();

  QIO_USQCDPropSourceInfo *qio_UserRecordInfo;
  qio_UserRecordInfo = QIO_create_usqcd_propsource_info("");

  return_val_info += QIO_read_record_info( qio_Input, record, record_xml);

  return_val_decode += QIO_decode_usqcd_propsource_info(qio_UserRecordInfo, record_xml);


#ifdef IGNORE_QIO_BAD_XML
  if(return_val_info == QIO_BAD_XML){
    
    VRB.Result(cname,fname,"WARNING QIO: reading info (source record) returned %d (QIO_BAD_XML), IGNORING...\n",return_val_info);
    return_val_info=0;

  }

   if(return_val_decode == QIO_BAD_XML){
    
    VRB.Result(cname,fname,"WARNING QIO: decoding info (source record) returned %d (QIO_BAD_XML), IGNORING...\n",return_val_decode);
    return_val_decode=0;

  }
#endif //IGNORE_QIO_BAD_XML


  if(return_val_decode != 0 || return_val_info != 0 )
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: failed reading next source record info! read record: %d, decode_usqcd_propsource_info: %d\n", return_val_decode, return_val_info );
 
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      QIO_destroy_usqcd_propsource_info(qio_UserRecordInfo);
      return (1);
    }


  int sizepersite(0);

  if(sType == QIO_SCALAR_SOURCE)
	 sizepersite =1 ;
  if(sType == QIO_FULL_SOURCE)
	 sizepersite = QIO_PROP_SPIN_MAX*QIO_PROP_COLOR_MAX ;


  int readDatacount=QIO_get_datacount(record);
  int readDatasize=QIO_get_typesize(record);
  
  char *readPrecision = QIO_get_precision(record);
  int basesize(0);

  char switchPrecision('U');

  if(strcmp( readPrecision, "D" ) == 0)
    {switchPrecision='D'; basesize=2*sizeof(Float);}
  if(strcmp( readPrecision, "F" ) == 0)
    {switchPrecision='F'; basesize=2*sizeof(float);} 


  //check datacount

  if( int(readDatacount*readDatasize) != int(sizepersite*basesize))
    {
      ERR.General(cname, fname,"ERROR QIO-Prop: datacount/typesize mismatch:\n read: %i of size %i ( = %i),\n expected: 1 of size %i x %i = %i\n",
		  readDatacount, readDatasize, (readDatacount*readDatasize),
		  sizepersite, basesize, (sizepersite*basesize));
      
      
      QIO_destroy_record_info(record);
      QIO_string_destroy(record_xml);
      QIO_destroy_usqcd_propsource_info(qio_UserRecordInfo);
      return (1);
    }

  switch(switchPrecision)
    {


    case 'D':

      if( sType == QIO_SCALAR_SOURCE )
	return_val += QIO_read(qio_Input, record, record_xml, qio_putScSource, 2*sizepersite*sizeof(Float), sizeof(Float), rsource);
      else
	return_val += QIO_read(qio_Input, record, record_xml, qio_putTmpSource, 2*sizepersite*sizeof(Float), sizeof(Float), rsource);
      break;


    case 'F':

      
      if( sType == QIO_SCALAR_SOURCE )
	return_val += QIO_read(qio_Input, record, record_xml, qio_putScSourceSingle, 2*sizepersite*sizeof(float), sizeof(float), rsource);
      else
	return_val += QIO_read(qio_Input, record, record_xml, qio_putTmpSourceSingle, 2*sizepersite*sizeof(float), sizeof(float), rsource);
      break;

    default:
      ERR.General(cname, fname, "ERROR QIO: precision out of range: %s\n", readPrecision);
      ++return_val;
      
      
    }


  QIO_destroy_record_info(record);
  QIO_string_destroy(record_xml);
  QIO_destroy_usqcd_propsource_info(qio_UserRecordInfo);
 
  return return_val;
  
  
}



void qio_readPropagator::qio_guessSpinColor(int &spin, int &color)
{

#ifdef DEBUG_ReadSpinColor
  printf("qio_guessSpinColor has been called...\n");
#endif


  char * fname="qio_guessSpinColor()";
  VRB.Func(cname,fname);

  int tmpS, tmpC;

  tmpS = lastSpin;

  tmpC = lastColor + 1;
  if (tmpC == QIO_PROP_COLOR_MAX){
    tmpC = 0; 
    ++tmpS;
  }

  if ( tmpS > QIO_PROP_SPIN_MAX ){
    ERR.General(cname, fname, "ERROR QIO: SPIN OUT OF RANGE!!!\n");
    exit(-13);
  }
 
  spin = tmpS;
  color = tmpC;

  lastSpin = tmpS;
  lastColor = tmpC;

#ifdef DEBUG_ReadSpinColor
  printf("qio_guessSpinColor guessed spin %i color %i.\n",spin,color);
#endif

  
}






CPS_END_NAMESPACE
#endif
