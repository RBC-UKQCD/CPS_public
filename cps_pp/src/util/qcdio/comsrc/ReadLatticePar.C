#include <config.h>
#include <math.h>
#include <util/ReadLatticePar.h>
#include <util/time_cps.h>

#ifdef USE_QMP
#include <qmp.h>
#endif

CPS_START_NAMESPACE
using namespace std;

#define PROFILE
void ReadLatticeParallel::read(Lattice & lat, const QioArg & rd_arg)
{
  const char * fname = "read()";
  VRB.Func(cname,fname);
  
  char loginfo[100];
  sprintf(loginfo,"Load %s",rd_arg.FileName);
  startLogging(loginfo);

#ifdef PROFILE
  struct timeval start,end;
  gettimeofday(&start,NULL);
#endif

  io_good = false;
  int error = 0;

  if (isRoot()) { // commander, analyze file header
    
    ifstream input(rd_arg.FileName);
    if ( !input.good() )
      {
	//	VRB.Flow(cname,fname,"Could not open file [%s] for input.\n",rd_arg.FileName);
	//	VRB.Flow(cname,fname,"USER: maybe you should kill the process!!\n");
	error = 1;
      }

    if(!error) {
      hd.read(input);
      input.close();
    }
  }
  if(synchronize(error) != 0)
    ERR.FileR(cname, fname, rd_arg.FileName);
  log();

#ifndef USE_QMP
	int temp_start = hd.data_start;
	broadcastInt(&temp_start);
	hd.data_start = temp_start;
#else
  QMP_broadcast(&hd.data_start, sizeof(long));
  broadcastInt(&hd.recon_row_3);
  //  cout << "recon_row_3 = " << hd.recon_row_3 << endl;
#endif



  // check all conditions between FILE and GJP
  int nx = rd_arg.Xnodes() * rd_arg.XnodeSites();
  int ny = rd_arg.Ynodes() * rd_arg.YnodeSites();
  int nz = rd_arg.Znodes() * rd_arg.ZnodeSites();
  int nt = rd_arg.Tnodes() * rd_arg.TnodeSites();

  if(isRoot()) {
    if(hd.dimension[0] != nx || hd.dimension[1] != ny || hd.dimension[2] != nz || hd.dimension[3] != nt) {
      VRB.Flow(cname,fname,"Dimensions in file DISAGREE with GlobalJobParameter!\n");
      VRB.Flow(cname,fname,"In File: %d x %d x %d x %d\n",
	       hd.dimension[0],hd.dimension[1], hd.dimension[2], hd.dimension[3]);
      VRB.Flow(cname,fname,"In GJP:  %d x %d x %d x %d\n",nx, ny, nz, nt);
      error = 1;
    }

// Turned off the boundary check, as it is inconsistent  with NERSC convention.
// 04/03/05, CJ
#if 0
    if(hd.boundary[0] != rd_arg.Xbc() || hd.boundary[1] != rd_arg.Ybc() || hd.boundary[2] != rd_arg.Zbc() || hd.boundary[3] != rd_arg.Tbc()) {
      VRB.Flow(cname,fname,"Boundary conditions in file DISAGREE with GlobalJobParameter!\n");

      VRB.Flow(cname,fname,"In File: %s x %s x %s x %s\n",
	       hd.boundary[0]==BND_CND_PRD ? "P":"A",
	       hd.boundary[1]==BND_CND_PRD ? "P":"A",
	       hd.boundary[2]==BND_CND_PRD ? "P":"A",
	       hd.boundary[3]==BND_CND_PRD ? "P":"A");

      VRB.Flow(cname,fname,"In GJP:  %s x %s x %s x %s\n",
	       rd_arg.Xbc()==BND_CND_PRD ? "P":"A",
	       rd_arg.Ybc()==BND_CND_PRD ? "P":"A",
	       rd_arg.Zbc()==BND_CND_PRD ? "P":"A",
	       rd_arg.Tbc()==BND_CND_PRD ? "P":"A");

      error = 1;
    }
 #endif
  }

  if(synchronize(error) != 0)  
    ERR.General(cname, fname, "Wrong Parameters Specified\n");

  // see if file Floating Points is acceptable
  if(isRoot()) {
    fpconv.setFileFormat(hd.floating_point);
  }
  VRB.Flow(cname,fname,"FileFormat=%d",hd.floating_point);
  broadcastInt((int*)&fpconv.fileFormat);
  if(fpconv.fileFormat == FP_UNKNOWN) {
    ERR.General(cname,fname, "Data file Floating Point Format UNKNOWN\n");
  }
  
  VRB.Flow(cname,fname,"A copy of header info from file:\n");
  if(isRoot())  hd.show();

  int data_per_site = hd.recon_row_3 ? 4*12 : 4*18;

  // read lattice data, using parallel style or serial (all on node 0) style
  unsigned int csum;


#if TARGET != QCDOC
//  setSerial();
#endif

  log();

  VRB.Flow(cname,fname,"Reading configuation to address: %p\n", rd_arg.StartConfLoadAddr);
  if(parIO()) {
    ParallelIO pario(rd_arg);
    if(! pario.load((char*)rd_arg.StartConfLoadAddr, data_per_site, sizeof(Matrix)*4,
		    hd, fpconv, 4, &csum))  
      ERR.General(cname,fname,"Load Failed\n");  // failed to load
  }
#if 1
  else {
    SerialIO serio(rd_arg);
    if(! serio.load((char*)rd_arg.StartConfLoadAddr, data_per_site, sizeof(Matrix)*4,
		    hd, fpconv, 4, &csum))
      ERR.General(cname,fname,"Load Failed\n");  // failed to load
  }
#endif

  log();

//  printf("Node %d: lattice read csum=%x\n",UniqueID(),csum);
  //  cout << "loader finish, csum = " << hex << csum << dec << endl << endl;
  //  cout << "loader done" << endl << endl;


  // After reading...
  // STEP 1: checksum
  if(rd_arg.Scoor() == 0)
    csum = globalSumUint(csum);
  else
    globalSumUint(0);

  if(isRoot()) {
    if( hd.checksum != csum ) {
      VRB.Flow(cname,fname, "CheckSUM error !! Header: %x  Host calc: %x\n",hd.checksum,csum);
      
      printf("Node %d: CheckSUM error !! Header: %x  Host calc: %x\n",UniqueID(),hd.checksum,csum);
      error = 1;
    }
    else
      VRB.Flow(cname,fname,"CheckSUM is ok\n");
  }

  if(synchronize(error) != 0) 
    ERR.General(cname, fname, "Checksum error\n");


  // STEP 2: reconstruct row 3
  int size_matrices( rd_arg.XnodeSites() * rd_arg.YnodeSites() 
		     * rd_arg.ZnodeSites() * rd_arg.TnodeSites() * 4); 

  Matrix * lpoint = rd_arg.StartConfLoadAddr;

  if(hd.recon_row_3) {
    VRB.Flow(cname,fname,"Reconstructing row 3\n");
    for(int mat=0; mat<size_matrices; mat++) {
      Float * rec = (Float*)&lpoint[mat];
      // reconstruct the 3rd row
      rec[12] =  rec[2] * rec[10] - rec[3] * rec[11] - rec[4] * rec[8] + rec[5] * rec[9];
      rec[13] = -rec[2] * rec[11] - rec[3] * rec[10] + rec[4] * rec[9] + rec[5] * rec[8];
      rec[14] = -rec[0] * rec[10] + rec[1] * rec[11] + rec[4] * rec[6] - rec[5] * rec[7];
      rec[15] =  rec[0] * rec[11] + rec[1] * rec[10] - rec[4] * rec[7] - rec[5] * rec[6];
      rec[16] =  rec[0] * rec[ 8] - rec[1] * rec[ 9] - rec[2] * rec[6] + rec[3] * rec[7];
      rec[17] = -rec[0] * rec[ 9] - rec[1] * rec[ 8] + rec[2] * rec[7] + rec[3] * rec[6];
    }
  }

  // STEP 3: check plaq and linktrace
  if(lat.GaugeField() != lpoint) lat.GaugeField(lpoint);
  if(! CheckPlaqLinktrace(lat,rd_arg, hd.plaquette, hd.link_trace))  
    ERR.General(cname,fname,"Plaquette or Link trace check failed\n");

#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(cname,"read",0,&start,&end);
#endif

  io_good = true;

  log();
  finishLogging();

  lat.ClearSmeared();
  VRB.FuncEnd(cname,fname);
};


bool ReadLatticeParallel::CheckPlaqLinktrace(Lattice &lat, const QioArg & rd_arg,
					     const Float plaq_inheader, const Float linktrace_inheader) 
{
  const char * fname = "CheckPlaqLinktrace()";
  int error = 0;

  Float plaq = lat.SumReTrPlaq() / 18.0 / rd_arg.VolSites() ;
  Float devplaq(0.0);
  if(isRoot()) {
    devplaq = fabs(  (plaq - plaq_inheader) / plaq ) ;
    printf("%s::%s: plaquette::\n  calc: %0.8e  header: %0.8e   rel.dev.: %0.8e\n",
	     cname,fname,plaq, plaq_inheader, devplaq);
  }

  Float linktrace(0);
  int is;
  Matrix *m =  lat.GaugeField(); 
  for(is=0;is< rd_arg.VolNodeSites()*4; is++){
    linktrace += m->ReTr();
    m++;
  }
  if(rd_arg.Scoor() == 0) 
    linktrace = globalSumFloat(linktrace) / (rd_arg.VolSites()*12.0);
  else
    globalSumFloat(0.0);

  if(isRoot()) {
    Float devlinktrace =   
      fabs(  (linktrace - linktrace_inheader) / linktrace );

    printf("%s::%s: linktrace::\n  calc: %0.8e  header: %0.8e   rel.dev.: %0.8e\n",
	     cname,fname,linktrace, linktrace_inheader, devlinktrace);
  
    Float chkprec = rd_arg.CheckPrecision;

//  CJ:  turning off the link trace test, as some lattices have a very small link trace by accident.
//    if(devplaq > chkprec || devlinktrace > chkprec) {
//      VRB.Flow(cname,fname, "Plaquette and/or Link trace different from header\n");
    if(devplaq > chkprec) {
      VRB.Flow(cname,fname, "Plaquette different from header\n");
      error = 1;
    }
  }

  if(synchronize(error) != 0) return false;

  return true;
}

// added EES for getting header info consistent on every node
int ReadLatticeParallel::getSequenceNumber(){

  int tmp=hd.sequence_number;
  broadcastInt(&tmp);
  return tmp;

}




CPS_END_NAMESPACE
