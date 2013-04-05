#ifdef USE_QIO
#include <config.h>
#include <util/qio_writeMuon.h>

CPS_START_NAMESPACE
using namespace std;


// qio-factory functions outside class!!

// get the vacuum polarization or muon line 
void qio_getField(char *buf, size_t index, int count, void *arg)
// arg is address of vacpol
// index is site on node?
{
  Rcomplex *start = (Rcomplex *)arg; // start of data
  Rcomplex *wrt_buf = (Rcomplex *)buf; // start of write buffer
  Rcomplex *src_buf = start + index*count; // start of actual site data
  for(int ii(0); ii< count; ++ii)
    *(wrt_buf + ii) = *(src_buf + ii); // copy data to buffer
}


void qio_getFieldSingle(char *buf, size_t index, int count, void *arg)
{


  //float *cc = (float *) buf;

  //double *muon = (double *)arg;
 
  // 2 is for re,im
  //double *src = muon + index * count*2 ;
  
  //for(int ii(0); ii < count*2; ++ii)
   // *(cc + ii) = *(src +ii);
  
#ifdef DEBUG_GetField
  printf("UID: %i, finished qio_getFieldSingle (called with index: %i, count: %i)\n",NodeID(), index, count);
#endif // DEBUG_GetField

}



// now start class-functions...

void qio_writeMuon::qio_openOutput(char *filename, char *stringLFN, 
				   char *xml_write_file, int volFormat)
{

  char * fname = "qio_openOutput(...)";
  VRB.Func(cname,fname);

  const int volfmt(volFormat);
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

  outfile = QIO_open_write(xml_file_out, filename, volfmt, &layout, NULL, &oflag);

  QIO_string_destroy(xml_file_out);
  
  qio_Output = outfile;

  VRB.FuncEnd(cname,fname);

}



void qio_writeMuon::write(char *outfile,
			  Lattice &lat,
			  Rcomplex *muon,
			  int count,
			  char *ildgLFN, 
			  int volFormat, 
			  FP_FORMAT floatFormat)
{

  const char * fname = "write(...)";

  VRB.Func(cname,fname);

  int return_val(0);

#ifdef DEBUG_GlobalWrite
  printf(" initializing global-data\n");
#endif


  Float plaq( lat.SumReTrPlaq()/(18*GJP.VolSites()) );

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

  //globdat[1] = ltrace / (4*3*GJP.VolSites());

  ltrace =  ltrace / (4*3*GJP.VolSites());

#ifdef DEBUG_GlobalWrite
  printf("UID: %i: loaded %f %f to GlobalData\n",NodeID(),plaq, ltrace);
#endif // DEBUG_GlobalWrite



  QIO_RecordInfo *record;
  //record = QIO_create_record_info(0,NULL,NULL,0, "", "", 0,0,0,0);


  QIO_String *record_xml;
  record_xml = QIO_string_create();

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

  
  //char xml_file[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Dummy QCDML</title>";
  char xml_file[] = QIO_XML_FILE_GAUGE ;

  char xml_plaq[20];
  char xml_linktr[20];


  if(SingleDouble)
    {
      sprintf( xml_plaq, "%.12g", plaq);
      sprintf( xml_linktr, "%.12g", ltrace); 
    }
  else
    {
      sprintf( xml_plaq, "%.8g", plaq);
      sprintf( xml_linktr, "%.8g", ltrace);
    }
  

  #ifdef DEBUG_GlobalWrite
  VRB.Result(cname,fname,"%.12g %s\n%.12g %s\n",plaq,xml_plaq,ltrace,xml_linktr);
  #endif // DEBUG_GlobalWrite


  // this can contain the add. header
  // char xml_info[] = "gauge_field";
 
  char xml_info[6*(MAX_HEADER_LINE+10)];



  sprintf(xml_info, "DATATYPE = 4D_SU3_GAUGE\nLINK_TRACE = %.10g\nPLAQUETTE = %.10g\nENSEMBLE_ID = %s\nENSEMBLE_LABEL = %s\nSEQUENCE_NUMBER = %i",
	  ltrace, plaq, header_ensemble_id, header_ensemble_label, header_traj);


  //check for length...
  if( (strlen(xml_info)*sizeof(char)) > QIO_INFO_STRING_MAX){
    ERR.General(cname,fname," xml_info too large: %i\n %s\n", (strlen(xml_info)*sizeof(char)), xml_info);
    exit(-12);
  }




  #ifdef DEBUG_GlobalWrite
  VRB.Result(cname, fname," adding plaq %s and linkTr %s to xml-record\n add. info:\n %s\n", xml_plaq, xml_linktr, xml_info);
  #endif //DEBUG_GlobalWrite

  QIO_USQCDLatticeInfo *userrecordinfo = QIO_create_usqcd_lattice_info(xml_plaq, xml_linktr, xml_info);

  QIO_encode_usqcd_lattice_info(record_xml,userrecordinfo); 


  Matrix * wlat = lat.GaugeField();
  
  //char ildg_lfn[] = "ildg-lfn n/a";



  if( strlen(ildgLFN) >  0 )

      qio_openOutput(outfile, ildgLFN, xml_file, volFormat);

  else
    {
      char ildg_lfn[] = "ildg-lfn n/a";
      qio_openOutput(outfile, ildg_lfn, xml_file, volFormat);
    }


  if(SingleDouble)
    {
      //output in double-precision

      // create the record info
      record = QIO_create_record_info(QIO_FIELD,NULL,NULL,0, "QDP_FD_ColorMatrix", "D", 3, 4, 2*sizeof(Float), count);
      // color=3 (not used), spin=4 (not used), size of Rcomplex, # Rcomplex
      // call write...  
      //return_val = QIO_write( qio_Output, record, record_xml, qio_getField, 
      //count*2*sizeof(Float), 2*sizeof(Float), muon);
      return_val = QIO_write( qio_Output, record, record_xml, qio_getField, 
			      count*2*sizeof(Float), sizeof(Float), muon);
    }
  else
    {
      //output in single-precision

      // create the record info
      record = QIO_create_record_info(QIO_FIELD,NULL,NULL,0, "QDP_F3_ColorMatrix", "F", 3, 4, 2*sizeof(Float), count);
      // color=3 (not used), spin=4 (not used), size of Rcomplex, # Rcomplex
      //return_val = QIO_write( qio_Output, record, record_xml, 
      //qio_getFieldSingle, count*2*sizeof(float), 
      //		      2*sizeof(float), muon);  
      return_val = QIO_write( qio_Output, record, record_xml, 
			      qio_getFieldSingle, count*2*sizeof(float), 
			      sizeof(float), muon);  
    }


  if ( return_val == 0 )
    VRB.Result(cname,fname,"QIO_write successfull...\n");
  else
    ERR.General(cname,fname,"ERROR QIO: QIO_write returned %d \n",return_val);


#ifdef PRINT_checksums
  uint32_t writeCheckA(QIO_get_writer_last_checksuma( output));
  uint32_t writeCheckB(QIO_get_writer_last_checksumb( output));
 
  printf("Checksums: a: %s -- b: %s \n",writeCheckA, writeCheckB);
#endif //PRINT_checksums

  // clean-up
  QIO_destroy_record_info(record);
  QIO_destroy_usqcd_lattice_info(userrecordinfo);  
  QIO_string_destroy(record_xml);

  qio_closeOutput();

  VRB.FuncEnd(cname,fname);

}


void qio_writeMuon::setHeader(const char * ensemble_id, const char * ensemble_label, const int traj)
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

CPS_END_NAMESPACE
#endif

