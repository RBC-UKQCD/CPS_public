// A slight modification to WriteLattice class
// to enable parallel writing of Gauge Connection Format

// Write the format of Gauge Connection Format
// of QCDSP {load,unload}_lattice format

#include <config.h>
#include <util/WriteLatticePar.h>
#include <iomanip>
#include <fstream>
#include <iostream>

CPS_START_NAMESPACE
using namespace std;

// code from cdawson/kostas

void WriteLatticeParallel::writeHeader(ostream & fout, Float link_trace, Float plaq, const QioArg & wt_arg)
{
    // default archive time is current time
    time_t ptm;
    time(&ptm);
    cout << "time = " << ptm << endl;
    hd_archive_date = asctime(localtime(&ptm));
    int i1( hd_archive_date.find_last_not_of("\n"));
    hd_archive_date = hd_archive_date.substr(0,i1+1);

    hd_ensemble_id = "unspecified";
    hd_ensemble_label = "unspecified";
    hd_sequence_number = 0;
    hd_creator = "unspecified";
    hd_creator_hardware = "unspecified";
    hd_creation_date = "unspecified";

    cout << "header init'd" << endl;


  if(isRoot()) {
    fout.seekp(0,ios::beg);
    fout << "BEGIN_HEADER\n";
    fout << "HDR_VERSION = 1.0\n";
    if(recon_row_3) 
      fout << "DATATYPE = 4D_SU3_GAUGE\n";
    else
      fout << "DATATYPE = 4D_SU3_GAUGE_3x3\n";
    fout << "STORAGE_FORMAT = 1.0\n";
    for(int i=1;i<=4;i++){
      fout << "DIMENSION_" << i << " = " << wt_arg.Nodes(i-1)*wt_arg.NodeSites(i-1) << "\n" ;
    }
    // just to keep the space and write it later
    fout << "LINK_TRACE = " << setprecision(10) << link_trace  << "\n";
    fout << "PLAQUETTE  = " << setprecision(10) << plaq << "\n";
    for(int i=1;i<=4;i++){
      fout << "BOUNDARY_"<<i<<" = " <<
	( (wt_arg.Bc(i-1) == BND_CND_APRD ) ? "ANTIPERIODIC\n" :"PERIODIC\n");
    }
    fout << "CHECKSUM = ";

    // store checksum position
    csum_pos = fout.tellp();

    fout << hex << setw(8) << 0 << dec << endl;

    fout << "ENSEMBLE_ID = " << hd_ensemble_id << endl;
    fout << "ENSEMBLE_LABEL = "<< hd_ensemble_label << endl;
    fout << "SEQUENCE_NUMBER = " << hd_sequence_number << endl;
    fout << "CREATOR = " << hd_creator << endl;
    fout << "CREATOR_HARDWARE = " << hd_creator_hardware << endl;
    fout << "CREATION_DATE = " << hd_creation_date << endl;
    fout << "ARCHIVE_DATE = " << hd_archive_date << endl;
    fout << "FLOATING_POINT = " << fpconv.name(fpconv.fileFormat) << endl;
    fout << "END_HEADER" << endl;

    cout << "header written" << endl;

    data_start = fout.tellp();
  }    

  broadcastInt(&data_start); // from 0 to all
}

void WriteLatticeParallel::write(Lattice & lat, const QioArg & wt_arg)
{
  cout << endl << "Unloading lattice..." << endl << endl;

  // init
  int error = 0;

  const char * filename = wt_arg.FileName;
  recon_row_3 = wt_arg.ReconRow3;
  FP_FORMAT dataFormat = wt_arg.FileFpFormat;

  fpconv.setFileFormat(dataFormat);
  if(fpconv.fileFormat == FP_UNKNOWN) {
    cout << "Output Floating Point format UNKNOWN"<< endl;
    return;
  }

  //  calc Plaq and LinkTrace
  const int size_matrices(wt_arg.VolNodeSites()*4);
  Matrix *lpoint =  lat.GaugeField();
  cout << "Writing Gauge Field at Lattice::GaugeField() = " << hex << lpoint << dec << endl;

  Float plaq = lat.SumReTrPlaq()/(18*wt_arg.VolSites()) ;
  Float ltrace(0.0);
  if(wt_arg.Scoor() == 0) {
    for(int i=0;i<size_matrices;i++){
      ltrace += (lpoint+i)->ReTr();
    }
    ltrace = globalSumFloat(ltrace) / (4*3*wt_arg.VolSites());
  }
  else
    globalSumFloat(0.0);  // everyone has to participate in global ops


  // all open file, start writing
  ofstream output(filename);

  if(!output.good())
    {
      cout << "Could not open file:\n   "
	   << filename
           << "\nfor output.\nNot even a little bit.\n";
      error = 1;
    }

  if(synchronize(error) > 0)  return;

  writeHeader(output, ltrace, plaq, wt_arg);// write a header without CSUM


  const int Floats_per_site = recon_row_3 ? 4*12 : 4*18;
  const int chars_per_site = Floats_per_site * fpconv.fileFpSize();


  // TempBufAlloc is a mem allocator that prevents mem leak on "return"
  TempBufAlloc  fbuf(chars_per_site);

  int xbegin = wt_arg.XnodeSites() * wt_arg.Xcoor(), xend = wt_arg.XnodeSites() * (wt_arg.Xcoor()+1);
  int ybegin = wt_arg.YnodeSites() * wt_arg.Ycoor(), yend = wt_arg.YnodeSites() * (wt_arg.Ycoor()+1);
  int zbegin = wt_arg.ZnodeSites() * wt_arg.Zcoor(), zend = wt_arg.ZnodeSites() * (wt_arg.Zcoor()+1);
  int tbegin = wt_arg.TnodeSites() * wt_arg.Tcoor(), tend = wt_arg.TnodeSites() * (wt_arg.Tcoor()+1);

  int nx = wt_arg.XnodeSites() * wt_arg.Xnodes();
  int ny = wt_arg.YnodeSites() * wt_arg.Ynodes();
  int nz = wt_arg.ZnodeSites() * wt_arg.Znodes();
  //  int nt = wt_arg.TnodeSites() * wt_arg.Tnodes();  // not needed


  int tblk = nx*ny*nz*chars_per_site;
  int zblk = nx*ny*chars_per_site;
  int yblk = nx*chars_per_site;

  cout << endl;
  cout << "Trying to write " << wt_arg.VolNodeSites() * chars_per_site << " bytes"<<endl;

  // write in parallel manner, node 0 will assign & dispatch IO time slots
  int mat=0;
  unsigned int csum = 0;

  setConcurIONumber(wt_arg.ConcurIONumber);
  getIOTimeSlot();

  if(wt_arg.Scoor() == 0) {
    output.seekp(data_start,ios_base::beg);
    int jump = tbegin * tblk;
    for(int tr=tbegin;tr<tend;tr++) {
      jump += zbegin * zblk;
      for(int zr=zbegin;zr<zend;zr++) {
	jump += ybegin * yblk;
	for(int yr=ybegin;yr<yend;yr++) {
	  jump += xbegin * chars_per_site;
	  output.seekp(jump, ios_base::cur);

	  for(int xr=xbegin;xr<xend;xr++) {
	    for(int i=0;i<4;i++) {
	      fpconv.host2file((char*)fbuf + chars_per_site/4*i, (char*)&lpoint[mat++],
			       Floats_per_site/4);
	    }
	    csum += fpconv.checksum((char*)fbuf,Floats_per_site);
	    output.write(fbuf,chars_per_site);
	  }	  
	  
	  jump = (nx-xend) * chars_per_site;  // jump restarted from 0
	}
	jump += (ny-yend) * yblk;
      }
      jump += (nz-zend) * zblk;
    }
    if ( !output.good() ) { cout << "Output stream error!" << endl; error = 1; }
  }

  finishIOTimeSlot();
  //

  if(wt_arg.Scoor() == 0) {
    csum = globalSumUint(csum);
  }
  else
    globalSumUint(0);

  // fill in checksum
  if(isRoot()) {
    output.seekp(csum_pos,ios::beg);
    output << hex << setw(8) << csum << dec;
  }

  output.close();

  if(synchronize(error) != 0)  return;

  unload_good = true;
}


CPS_END_NAMESPACE
