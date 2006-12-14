#include <config.h>
#include <util/latheader.h>
#include <util/qioarg.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <unistd.h>
using namespace std;

CPS_START_NAMESPACE

///////////////////////////////////////////////////////////////////
// GCFheaderPar class members
string elmSpacePar(string str)
{
  const int i0(str.find_first_not_of(" "));
  const int i1(str.find_last_not_of (" "));
  if(i1 - i0>0){ return(str.substr(i0,i1-i0+1)); }
  else         { return(str);  }
}

bool GCFheaderPar::add(string key_eq_value)
{
  const int eqp(key_eq_value.find("="));
  if( eqp  > 0  )
    {
      const string key( elmSpacePar( key_eq_value.substr(0,eqp) ) );
      const string val( elmSpacePar( key_eq_value.substr(eqp+1) ) );
      headerMap.insert(GCFHMapParT::value_type(key,val));
      return true;
    } 
  else 
    {
      return false;
    }
}


void GCFheaderPar::Show() const
{
  cout << endl;
  for (GCFHMapParT::const_iterator iter = headerMap.begin(); 
       iter != headerMap.end(); ++iter) 
    {
      cout << iter->first << ":" << iter->second << endl;
    }
  cout << endl;
};


string GCFheaderPar::asString( const string key ) 
{
  GCFHMapParT::const_iterator n(headerMap.find(key));

  if (n == headerMap.end()) {
    cout << "header::asString key " << key << " not found. use Default." << endl;
    prevFound = false;
    return string("");
  }
  else {
    prevFound = true;
    return ( n->second );
  }
}		  


int GCFheaderPar::asInt( const string key ) 
{
  GCFHMapParT::const_iterator n(headerMap.find(key));

  if (n == headerMap.end()) {
    cout << "header::asInt key "<<key<<" not found. use Default." << endl;
    prevFound = false;
    return int(0);
  }

  else {
    prevFound = true;
    int tmp;
    sscanf((n->second).c_str() , "%d ", &tmp);
    return ( tmp );
  }
}		  

unsigned int GCFheaderPar::asHex( const string key ) 
{
  GCFHMapParT::const_iterator n(headerMap.find(key));

  if (n == headerMap.end()) {
    cout << "header::asHex key "<<key<<" not found. use Default." << endl;
    prevFound = false;
    return 0;
  }

  else {
    prevFound = true;
    int tmp;
    sscanf((n->second).c_str() , "%x ", &tmp);
    return ( tmp );
  }
}		  

Float GCFheaderPar::asFloat( const string key ) 
{
  GCFHMapParT::const_iterator n(headerMap.find(key));

  if (n == headerMap.end()) {
    cout << "header::asFloat key " << key << " not found. use Default." << endl;
    prevFound = false;
    return Float(0.0);
  }
  else {
    prevFound = true;
    float tmp;
    sscanf((n->second).c_str() , "%f ", &tmp);
    return ( (Float) tmp );
  }
}



/////////////////////////////////////////////////////////////////////
// LatticeHeader members 
/////////////////////////////////////////////////////////////////////
void LatticeHeader::init(const QioArg & qio_arg, FP_FORMAT FileFormat, Float LinkTrace, Float Plaq) {
  hdr_version = "1.0";
  recon_row_3 = qio_arg.ReconRow3;
  storage_format = "1.0";

  for(int i=0;i<4;i++)
    dimension[i] = qio_arg.Nodes(i)*qio_arg.NodeSites(i);

  link_trace = LinkTrace;
  plaquette = Plaq;

  for(int i=0;i<4;i++) 
    boundary[i] = qio_arg.Bc(i);

  checksum = 0;

  //  ensemble_id = "unspecified";
  //  ensemble_label = "unspecified";
  //  sequence_number = 0;
  

#if TARGET == QCDOC 
    creator = "RBC"; // getlogin() not supported on QCDOC yet
    creator_hardware = "QCDOC";
#elif TARGET == BGL 
    creator = "RBC"; // getlogin() not supported on QCDOC yet
    creator_hardware = "BGL";
#elif TARGET == QCDSP
    creator = "RBC";
    creator_hardware = "CU-QCDSP ";
#else
//    creator = getlogin();
    creator = "RBC";
    creator_hardware = "CU-NOARCH ";
#endif

    char buf[256];
#if TARGET == QCDOC
    strcpy(buf,"one rack"); // gethostname() not supported on QCDOC yet
#else
    gethostname(buf,256);
#endif
    creator_hardware += buf;

  // default archive time is current time
    struct timeval tp;
    gettimeofday(&tp,NULL);
    time_t ptm = tp.tv_sec;

    cout << "time = " << ptm << endl;
    archive_date = asctime(localtime(&ptm));
    int i1( archive_date.find_last_not_of("\n"));
    archive_date = archive_date.substr(0,i1+1);

    creation_date = archive_date;

    floating_point = FileFormat;

}


void LatticeHeader::setHeader(const char * EnsembleId, const char * EnsembleLabel, const int SequenceNumber, const char *CreatorName, const char *CreatorHardware) {
  ensemble_id = EnsembleId;
  ensemble_label= EnsembleLabel;
  sequence_number = SequenceNumber;
  if (CreatorName) creator = CreatorName;
  if (CreatorHardware) creator_hardware = CreatorHardware;
}

void LatticeHeader::write(ostream & fout) {
  fout.seekp(0,ios::beg);
  fout << "BEGIN_HEADER" << endl;
  fout << "HDR_VERSION = " << hdr_version << endl;
  if(recon_row_3) 
    fout << "DATATYPE = 4D_SU3_GAUGE" << endl;
  else
    fout << "DATATYPE = 4D_SU3_GAUGE_3x3" << endl;
  fout << "STORAGE_FORMAT = " << storage_format << endl;

  for(int i=0;i<4;i++){
    fout << "DIMENSION_" << i+1 << " = " << dimension[i] << endl ;
  }
  // just to keep the space and write it later
  fout << "LINK_TRACE = " << setprecision(10) << link_trace  << endl;
  fout << "PLAQUETTE  = " << setprecision(10) << plaquette << endl;

  for(int i=0;i<4;i++){
    fout << "BOUNDARY_"<<i+1<<" = " <<
      (boundary[i]== BND_CND_APRD? "ANTIPERIODIC" :"PERIODIC") << endl;
  }

  fout << "CHECKSUM = ";
  // store checksum position
  csum_pos = fout.tellp();
  fout << hex << setw(8) << 0 << dec << endl;
  
  fout << "ENSEMBLE_ID = " << ensemble_id << endl;
  fout << "ENSEMBLE_LABEL = "<< ensemble_label << endl;
  fout << "SEQUENCE_NUMBER = " << sequence_number << endl;
  fout << "CREATOR = " << creator << endl;
  fout << "CREATOR_HARDWARE = " << creator_hardware << endl;
  fout << "CREATION_DATE = " << creation_date << endl;
  fout << "ARCHIVE_DATE = " << archive_date << endl;

  fout << "FLOATING_POINT = " << FPConv::name(floating_point) << endl;

  fout << "END_HEADER" << endl;
  
  data_start = fout.tellp();
}


void LatticeHeader::fillInChecksum(ostream & fout, unsigned int checksum) const{
  fout.seekp(csum_pos);
  fout << hex << setw(8) << checksum << dec;
}


void LatticeHeader::read(istream & fin) {
  
  string line;
  do {
    getline(fin,line); // read one line
    hd.add(line);
  } while( line.find("END_HEADER") == string::npos);
  
  data_start = fin.tellg();
  
  // interpret header
  hdr_version = hd.asString("HDR_VERSION");
  if(hd.asString("DATATYPE") == "4D_SU3_GAUGE")
    recon_row_3 = 1;
  else
    recon_row_3 = 0;
  storage_format = hd.asString("STORAGE_FORMAT");

  dimension[0] = hd.asInt("DIMENSION_1");
  dimension[1] = hd.asInt("DIMENSION_2");
  dimension[2] = hd.asInt("DIMENSION_3");
  dimension[3] = hd.asInt("DIMENSION_4");

  link_trace = hd.asFloat("LINK_TRACE");
  plaquette = hd.asFloat("PLAQUETTE");

  boundary[0] = (hd.asString("BOUNDARY_1")=="ANTIPERIODIC" ? BND_CND_APRD : BND_CND_PRD);
  boundary[1] = (hd.asString("BOUNDARY_2")=="ANTIPERIODIC" ? BND_CND_APRD : BND_CND_PRD);
  boundary[2] = (hd.asString("BOUNDARY_3")=="ANTIPERIODIC" ? BND_CND_APRD : BND_CND_PRD);
  boundary[3] = (hd.asString("BOUNDARY_4")=="ANTIPERIODIC" ? BND_CND_APRD : BND_CND_PRD);

  checksum = hd.asHex("CHECKSUM");
  //  sscanf(hd.asString("CHECKSUM").c_str(), "%x ", &checksum);

  ensemble_id = hd.asString("ENSEMBLE_ID");
  ensemble_label = hd.asString("ENSEMBLE_LABEL");
  sequence_number = hd.asInt("SEQUENCE_NUMBER");
  creator = hd.asString("CREATOR");
  creator_hardware = hd.asString("CREATOR_HARDWARE");
  creation_date = hd.asString("CREATION_DATE");
  archive_date = hd.asString("ARCHIVE_DATE");

  FPConv fp;
  floating_point = fp.setFileFormat(hd.asString("FLOATING_POINT").c_str());
}


/////////////////////////////////////////////////////////////////////
// LatRngHeader members 
/////////////////////////////////////////////////////////////////////
void LatRngHeader::init(const QioArg & qio_arg, INT_FORMAT FileFormat) {
  hdr_version = "1.0";
  datatype = "LATTICE_RNG_5D_4D";
  storage_format = "1.0";

  for(int i=0;i<5;i++)
    dimension[i] = qio_arg.Nodes(i)*qio_arg.NodeSites(i);

  checksum = 0;
  pos_dep_csum = 0;

  average = 0.0;
  variance = 1.0;
  

#if TARGET == QCDOC
  creator = "RBC"; // getlogin() not supported on QCDOC yet
  creator_hardware = "CU-QCDOC ";
#elif TARGET == BGL
  creator = "RBC"; // getlogin() not supported on QCDOC yet
  creator_hardware = "BGL";
#elif TARGET == QCDSP
  creator = "RBC";
  creator_hardware = "CU-QCDSP ";
#else
  creator = "RBC";
  creator_hardware = "CU-NOARCH ";
#endif
  
  char buf[256];
#if TARGET == QCDOC
  strcpy(buf,"one rack"); // gethostname() not supported on QCDOC yet
#else
  gethostname(buf,256);
#endif
  creator_hardware += buf;
  
  // default archive time is current time
  struct timeval tp;
  gettimeofday(&tp,NULL);
  time_t ptm = tp.tv_sec;
  
  cout << "time = " << ptm << endl;
  archive_date = asctime(localtime(&ptm));
  int i1( archive_date.find_last_not_of("\n"));
  archive_date = archive_date.substr(0,i1+1);
  
  creation_date = archive_date;
  
  int_format = FileFormat;
  
}


void LatRngHeader::write(ostream & fout) {
  fout.seekp(0,ios::beg);
  fout << "BEGIN_HEADER" << endl;
  fout << "HDR_VERSION = " << hdr_version << endl;
  fout << "DATATYPE = " << datatype << endl;
  fout << "STORAGE_FORMAT = " << storage_format << endl;

  for(int i=0;i<5;i++){
    fout << "DIMENSION_" << i+1 << " = " << dimension[i] << endl ;
  }
  // just to keep the space and write it later
  fout << "CHECKSUM = ";
  csum_pos = fout.tellp();
  fout << hex << setw(8) << 0 << dec << endl;

  // due to an error in the calculation of the old tag "POS_DEP_CSUM", 
  // now use a new tag "ORDER_CSUM" to deprecate old one
  fout << "ORDER_CSUM = ";
  pdcs_pos = fout.tellp();
  fout << hex << setw(8) << 0 << dec << endl;

  fout << "AVERAGE = ";
  avg_pos = fout.tellp();
  fout << setw(20) << left << ' ' << endl;  // 20 whitespaces, hopefully

  fout << "VARIANCE = ";
  var_pos = fout.tellp();
  fout<< setw(20) << left << ' ' << endl;

  fout << "CREATOR = " << creator << endl;
  fout << "CREATOR_HARDWARE = " << creator_hardware << endl;
  fout << "CREATION_DATE = " << creation_date << endl;
  fout << "ARCHIVE_DATE = " << archive_date << endl;

  fout << "INT_FORMAT = " << IntConv::name(int_format) << endl;
    
  fout << "END_HEADER" << endl;
  
  data_start = fout.tellp();
}


void LatRngHeader::fillInCheckInfo(ostream & fout,
				   unsigned int csum, unsigned int pdcs,
				   Float avg, Float var) const {
  fout.seekp(csum_pos,ios::beg);
  fout << hex << setw(8) << csum << dec;
  
  fout.seekp(pdcs_pos, ios::beg);
  fout <<hex << setw(8) << pdcs << dec;

  char numstr[100];

  fout.seekp(avg_pos, ios::beg);
  sprintf(numstr,"%-20.10lf",avg);
  fout << numstr;

  fout.seekp(var_pos, ios::beg);
  sprintf(numstr,"%-20.10lf",var);
  fout << numstr;

}


void LatRngHeader::read(istream & fin) {
  
  string line;
  do {
    getline(fin,line); // read one line
    hd.add(line);
  } while( line.find("END_HEADER") == string::npos);
  
  data_start = fin.tellg();
  
  // interpret header
  hdr_version = hd.asString("HDR_VERSION");
  datatype = hd.asString("DATATYPE");
  storage_format = hd.asString("STORAGE_FORMAT");

  dimension[0] = hd.asInt("DIMENSION_1");
  dimension[1] = hd.asInt("DIMENSION_2");
  dimension[2] = hd.asInt("DIMENSION_3");
  dimension[3] = hd.asInt("DIMENSION_4");
  dimension[4] = hd.asInt("DIMENSION_5");
  if(!hd.found()) dimension[4] = 1;  // s-dim default to 1

  checksum = hd.asHex("CHECKSUM");

  // due to an error in the calculation of the old tag "POS_DEP_CSUM", 
  // now use a new tag "ORDER_CSUM" and deprecate the old one
  pos_dep_csum = hd.asHex("ORDER_CSUM");
  //cout << "ORDER_CSUM = " << pos_dep_csum << endl;

  average = hd.asFloat("AVERAGE");
  variance = hd.asFloat("VARIANCE");

  creator = hd.asString("CREATOR");
  creator_hardware = hd.asString("CREATOR_HARDWARE");
  creation_date = hd.asString("CREATION_DATE");
  archive_date = hd.asString("ARCHIVE_DATE");

  IntConv intconv;
  int_format = intconv.setFileFormat(hd.asString("INT_FORMAT").c_str());

}




CPS_END_NAMESPACE
