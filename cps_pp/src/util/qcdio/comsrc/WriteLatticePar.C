// A slight modification to WriteLattice class
// to enable parallel writing of Gauge Connection Format

// Write the format of Gauge Connection Format
// of QCDSP {load,unload}_lattice format

#include <config.h>
#include <util/WriteLatticePar.h>

CPS_START_NAMESPACE

// code from cdawson/kostas

void WriteLatticeParallel::initHeader()
{
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

void WriteLatticeParallel::write(Lattice & lat, const QioArg & wt_arg)
{
  cout << endl << "Unloading lattice..." << endl << endl;

  int error = 0;

  const char * filename = wt_arg.FileName;
  int recon_row_3 = wt_arg.ReconRow3;
  FP_FORMAT dataFormat = wt_arg.FileFpFormat;

  fpconv.setFileFormat(dataFormat);
  if(fpconv.fileFormat == FP_UNKNOWN) {
    cout << "Output Floating Point format UNKNOWN"<< endl;
    return;
  }

  const int size_matrices(wt_arg.VolNodeSites()*4);
  const int size_ints  = recon_row_3 ? size_matrices * 12 : size_matrices * 18;
  const int size_chars   ( size_ints * fpconv.fileFpSize() );

  Matrix *lpoint =  lat.GaugeField();
  cout << "Writing Gauge Field at Lattice::GaugeField() = " << hex << lpoint << dec << endl;

  //  calc Plaq and LinkTrace
  Float plaq = lat.SumReTrPlaq()/(18*wt_arg.VolSites()) ;
  Float ltrace(0.0);
  for(int i=0;i<size_matrices;i++){
     ltrace += (lpoint+i)->ReTr();
  }
  ltrace = globalSumFloat(ltrace) / (4*3*wt_arg.VolSites());


  cout << "Transferring " << fpconv.name(fpconv.hostFormat) << " ==> " << fpconv.name(fpconv.fileFormat) << endl;
  int size_per_mat = recon_row_3 ? 12*fpconv.fileFpSize() : 18*fpconv.fileFpSize();
  char * fpoint = new char[size_chars];
  if(recon_row_3) {
    for(int mat=0;mat<size_matrices;mat++) {
      fpconv.host2file(fpoint + mat * size_per_mat,(char*)&lpoint[mat], 12);
    }
  }
  else {
    fpconv.host2file(fpoint,(char*)lpoint, size_ints);
  }

  // calc checksum
  unsigned long csum = globalSumUint(fpconv.checksum(fpoint,size_ints));

  // header 
  ostringstream header;    
  if(isRoot()) {
    header << "BEGIN_HEADER\n";
    header << "HDR_VERSION = 1.0\n";
    if(recon_row_3) 
      header << "DATATYPE = 4D_SU3_GAUGE\n";
    else
      header << "DATATYPE = 4D_SU3_GAUGE_3x3\n";
    header << "STORAGE_FORMAT = 1.0\n";
    for(int i=1;i<=4;i++){
      header << "DIMENSION_" << i << " = " << wt_arg.Nodes(i-1)*wt_arg.NodeSites(i-1) << "\n" ;
    }
    header << "LINK_TRACE = "<< ltrace  << "\n";
    header << "PLAQUETTE  = "<< plaq << "\n";
    for(int i=1;i<=4;i++){
      header << "BOUNDARY_"<<i<<" = " <<
	( (wt_arg.Bc(i-1) == BND_CND_APRD ) ? "ANTIPERIODIC\n" :"PERIODIC\n");
    }
    header << "CHECKSUM = " << hex << csum << dec << "\n";
    header << "ENSEMBLE_ID = " << hd_ensemble_id << endl;
    header << "ENSEMBLE_LABEL = "<< hd_ensemble_label << endl;
    header << "SEQUENCE_NUMBER = " << hd_sequence_number << endl;
    header << "CREATOR = " << hd_creator << endl;
    header << "CREATOR_HARDWARE = " << hd_creator_hardware << endl;
    header << "CREATION_DATE = " << hd_creation_date << endl;
    header << "ARCHIVE_DATE = " << hd_archive_date << endl;
    header << "FLOATING_POINT = " << fpconv.name(fpconv.fileFormat) << endl;
    header << "END_HEADER" << endl;
    
    cout   << header.str();
  }    


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


  int data_start;
  if(isRoot()) {
    output << header.str();
    data_start = output.tellp();
  }
  
  broadcastInt(&data_start); // from 0 to all

  int xbegin = wt_arg.XnodeSites * wt_arg.Xcoor, xend = wt_arg.XnodeSites * (wt_arg.Xcoor+1);
  int ybegin = wt_arg.YnodeSites * wt_arg.Ycoor, yend = wt_arg.YnodeSites * (wt_arg.Ycoor+1);
  int zbegin = wt_arg.ZnodeSites * wt_arg.Zcoor, zend = wt_arg.ZnodeSites * (wt_arg.Zcoor+1);
  int tbegin = wt_arg.TnodeSites * wt_arg.Tcoor, tend = wt_arg.TnodeSites * (wt_arg.Tcoor+1);

  int nx = wt_arg.XnodeSites * wt_arg.Xnodes;
  int ny = wt_arg.YnodeSites * wt_arg.Ynodes;
  int nz = wt_arg.ZnodeSites * wt_arg.Znodes;
  //  int nt = wt_arg.TnodeSites * wt_arg.Tnodes;  // not needed

  int size_per_site = 4*size_per_mat;

  int tblk = nx*ny*nz*size_per_site;
  int zblk = nx*ny*size_per_site;
  int yblk = nx*size_per_site;

  cout << endl;
  cout << "Trying to write " << wt_arg.VolNodeSites() * size_per_site << " bytes"<<endl;

  char * buf = fpoint;

  // write in parallel manner, node 0 will assign & dispatch IO time slots
  setConcurIONumber(wt_arg.ConcurIONumber);
  getIOTimeSlot();

  output.seekp(data_start,ios_base::beg);
  int jump = tbegin * tblk;
  for(int tr=tbegin;tr<tend;tr++) {
    jump += zbegin * zblk;
    for(int zr=zbegin;zr<zend;zr++) {
      jump += ybegin * yblk;
      for(int yr=ybegin;yr<yend;yr++) {
	jump += xbegin * size_per_site;
	output.seekp(jump, ios_base::cur);

	output.write(buf,(xend-xbegin) * size_per_site);
	buf += (xend-xbegin) * size_per_site;

	jump = (nx-xend) * size_per_site;  // jump restarted from 0
      }
      jump += (ny-yend) * yblk;
    }
    jump += (nz-zend) * zblk;
  }
  if ( !output.good() ) { cout << "Output stream error!" << endl; error = 1; }

  finishIOTimeSlot();
  //
  
  output.close();

  cout << "Actually wrote " << buf-fpoint << " bytes" << endl;

  if(synchronize(error) != 0)  return;

  delete[] fpoint;

  unload_good = true;
}


CPS_END_NAMESPACE
