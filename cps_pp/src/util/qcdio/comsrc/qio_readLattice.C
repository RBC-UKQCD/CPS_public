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

  infile = QIO_open_read(xml_file_in, filename, &layout, NULL, &iflag);

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
  QIO_USQCDLatticeInfo *qio_UserRecordInfo;

  record_xml = QIO_string_create();


  qio_setLayout();

  // forces reading ???  ERROR after v2.1.2 !!!
  //layout.latdim = 0;
  

  // create the record info
   record         = QIO_create_record_info(0,NULL,NULL,0,  "", "", 0, 0, 0, 0);
   qio_UserRecordInfo = QIO_create_usqcd_lattice_info("", "", "");

  
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


    int readGlobaldata(QIO_get_recordtype(record));
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



   QIO_decode_usqcd_lattice_info(qio_UserRecordInfo, record_xml);
  

   // communicate T-chunks in S-direction
   if( GJP.Snodes() > 1)
     qio_communicateTchunks(rlat);
   
  
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






  if( QIO_defined_plaq(qio_UserRecordInfo) )
    {
      plaq_str = QIO_get_plaq(qio_UserRecordInfo);
      readPlaq = atof(plaq_str);
    }
  else
    {
      plaq_str = "UNDEFINED";
      readPlaq = 0.0;
    }
    
  if( QIO_defined_linktr(qio_UserRecordInfo) )
    {
      linktr_str = QIO_get_linktr(qio_UserRecordInfo);
      readLink = atof(linktr_str);
    }
  else
    {
      linktr_str = "UNDEFINED";
      readLink = 0.0;
    }

  if( QIO_defined_info(qio_UserRecordInfo) )
    {
      info_str = QIO_get_info(qio_UserRecordInfo);
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
  QIO_destroy_usqcd_lattice_info(qio_UserRecordInfo);




#ifdef DEBUG_ReadField
  printf("UID: %i, qio_readField finished with filenames: %s\n",UniqueID(),infile);
#endif // DEBUG_ReadField
  
  VRB.FuncEnd(cname,fname);


}


void qio_readLattice::qio_communicateTchunks(void *start)
{


  char * fname="qio_communicateTchunks()";
  VRB.Func(cname,fname);




#ifdef DEBUG_commTchunk
  printf("UID: %i, qio_communicateTchunks called\n",UniqueID());
#endif // DEBUG_commTchunk


  Matrix *mat = (Matrix *)start;


  //figure out chunk size
  // using getPlusData or getMinusData we run into trouble with max. number of IFloats (currently 4096) when the local lattice size (XYZT/S) is > 56 ( already 3^4 is problematic!)


  int chunk_size = 4 * GJP.NodeSites(0) * GJP.NodeSites(1) * GJP.NodeSites(2) * GJP.NodeSites(3)/ GJP.Snodes();

  int data_size = 18; //  3x3 complex matrix = 18 floats ,add a loop over chunk_size

#ifdef DEBUG_commTchunk
  printf("matrix size: %i, Float size: %i, comm-size: %i\n",sizeof(Matrix), sizeof(Float), COMMS_DATASIZE);
  printf("UID: %i, qio_communicateTchunks chunk_size: %i, data_size: %i\n",UniqueID(), chunk_size, data_size);
#endif // DEBUG_commTchunk


  int send_chunk, receive_chunk;

  send_chunk = GJP.SnodeCoor();

  receive_chunk = send_chunk + GJP.Snodes() - 1;
  receive_chunk %= GJP.Snodes();

  //loop over number nodes in S-direction
  for( int nn(0); nn < (GJP.Snodes()-1); ++nn){

    for( int mm(0); mm < chunk_size; ++mm)
      {
      
	
#ifdef DEBUG_commTchunk
	printf("UID: %i, qio_communicateTchunks: S-coor:%i sends %i, receives %i\n",UniqueID(),GJP.SnodeCoor(),send_chunk, receive_chunk);
#endif // DEBUG_commTchunk



	// send address, receive address
            
	IFloat *send_mat = (IFloat *) mat + send_chunk*data_size*chunk_size + mm*data_size;
	IFloat *receive_mat = (IFloat *) mat + receive_chunk*data_size*chunk_size + mm*data_size;

	// sending up, receiving down 

#ifdef DEBUG_commTchunk
	printf("UID: %i, qio_communicateTchunks: S-coor:%i start send/receive at %i/%i...\n",UniqueID(),GJP.SnodeCoor(), send_mat, receive_mat);
#endif // DEBUG_commTchunk



#ifdef DEBUG_commTchunk
	printf("UID: %i, qio_communicateTchunks: S-coor:%i ...initialized...\n",UniqueID(),GJP.SnodeCoor());
	
	Float *test1_s= (Float *)send_mat;
	Float *test1_r= (Float *)receive_mat;
	Float *test2_s= (Float *)send_mat + 17;
	Float *test2_r= (Float *)receive_mat + 17;
	printf("UID: %i, qio_communicateTchunks: S-coor:%i  send: %f, receive: %f, end: %f %f\n",UniqueID(),GJP.SnodeCoor(), *test1_s, *test1_r, *test2_s, *test2_r);
#endif // DEBUG_commTchunk

	getMinusData( receive_mat, send_mat, data_size, 4 );
      
#ifdef DEBUG_commTchunk
	printf("UID: %i, qio_communicateTchunks: S-coor:%i ...done!\n",UniqueID(),GJP.SnodeCoor());
	
	printf("UID: %i, qio_communicateTchunks: S-coor:%i  send: %f, receive: %f end: %f %f\n",UniqueID(),GJP.SnodeCoor(), *test1_s, *test1_r, *test2_s, *test2_r);
#endif // DEBUG_commTchunk

      } // loop mm
    

    // new send/receive
    
    send_chunk = receive_chunk;
    
    receive_chunk = send_chunk + GJP.Snodes() - 1;
    receive_chunk %= GJP.Snodes();
    


  } // loop nn


#ifdef DEBUG_commTchunk
  printf("UID: %i, qio_communicateTchunks finished\n",UniqueID());
#endif // DEBUG_commTchunk
  


  
  VRB.FuncEnd(cname,fname);

}



CPS_END_NAMESPACE
#endif
