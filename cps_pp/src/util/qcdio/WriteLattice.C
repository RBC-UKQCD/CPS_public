// Write the format of Gauge Connection Format
// of QCDSP {load,unload}_lattice format

#include <config.h>
#include <util/WriteLattice.h>
#include "Flconv.h"
CPS_START_NAMESPACE

// code from cdawson/kostas

WriteLattice::WriteLattice(char* file )
{
    filename=file;

    // default archive time is current time
    time_t ptm;
    time(&ptm);
    hd_archive_date = asctime(localtime(&ptm));
    int i1( hd_archive_date.find_last_not_of("\n"));
    hd_archive_date = hd_archive_date.substr(0,i1+1);

    hd_ensemble_id = "unspecified";
    hd_ensemble_label = "unspecified";
    hd_sequence_number = 0;
    hd_creator = "unspecified";
    hd_creator_hardware = "unspecified";
    hd_creation_date = "unspecified";

}

void WriteLattice::write()
{
  ofstream output(filename);
  cout << "FILENAME: " << filename << endl;
      
  GwilsonFnone lat;



  const int size_matrices(GJP.VolSites()*4);
  const int size_ints    ( size_matrices * ( 9 * 2 ) );
  const int size_chars   ( size_ints * sizeof(float) );

  Matrix *lpoint =  lat.GaugeField();


  //  calc Plaq and LinkTrace
  Float plaq = lat.SumReTrPlaq()/(18*GJP.VolSites()) ;
  Float ltrace(0.0);
  for(int i=0;i<size_matrices;i++){
     ltrace += (lpoint+i)->ReTr();
  }
  ltrace /= (4*3*GJP.VolSites());


  
  int xend=1;
  if(*(char *)&xend == 1){
    //Machine is Little Endian
    cout << "Machine is Little Endian\n";
    // reverse byte order
    byterevn((type32*)lpoint,size_ints);
  }

   // calc checksum
  unsigned int csum(0);
  unsigned int *ip = (unsigned int*)lpoint;
  for(int i=0; i < size_ints;i++)  csum += ip[i];

  // header 
  ostringstream header;

  header << "BEGIN_HEADER\n";
  header << "CHECKSUM = " << hex << csum << dec << "\n";
  header << "LINK_TRACE = "<< ltrace  << "\n";
  header << "PLAQUETTE  = "<< plaq << "\n";
  header << "DATATYPE = 4D_SU3_GAUGE_3x3\n";
  header << "HDR_VERSION = 1.0\n";
  header << "STORAGE_FORMAT = 1.0\n";
  for(int i=1;i<=4;i++){
    header << "DIMENSION_" << i << " = " << GJP.Nodes(i-1)*GJP.NodeSites(i-1) << "\n" ;
  }
  for(int i=1;i<=4;i++){
    header << "BOUNDARY_"<<i<<" = " <<
      ( (GJP.Bc(i-1) == BND_CND_APRD ) ? "ANTIPERIODIC\n" :"PERIODIC\n");
  }
  header << "ENSEMBLE_ID = " << hd_ensemble_id << endl;
  header << "ENSEMBLE_LABEL = "<< hd_ensemble_label << endl;
  header << "SEQUENCE_NUMBER = " << hd_sequence_number << endl;
  header << "CREATOR =" << hd_creator << endl;
  header << "CREATOR_HARDWARE = " << hd_creator_hardware << endl;
  header << "CREATION_DATE = " << hd_creation_date << endl;
  header << "ARCHIVE_DATE = " << hd_archive_date << endl;
  header << "FLOATING_POINT = IEEE32BIG\n";
  header << "END_HEADER\n";

  cout   << header.str();
  output << header.str();

  
  output.write((char*)lpoint,size_chars);
  output.close();


  if(*(char *)&xend == 1){
      //Machine is Little Endian
      cout << "Machie is Little Endian\n";
      // reverse byte order again
      byterevn((type32*)lpoint,size_ints);
  }

};

CPS_END_NAMESPACE
