
#include <config.h>

#include <stdlib.h>
#include <math.h>
#include <util/latrngio.h>
#include <util/qioarg.h>
#include <util/iostyle.h>
#include <util/intconv.h>
#include <util/random.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <sys/time.h>
#include <unistd.h>

CPS_START_NAMESPACE
using namespace std;


///////////////////////////////////////////////////////////////
// Lat Random Generator Read  functions 

LatRngRead::LatRngRead() 
  : LatRngIO()
{ }

LatRngRead::~LatRngRead() 
{ }

void LatRngRead::read(UGrandomGenerator * ugran, UGrandomGenerator * ugran_4d,
		      const QioArg & rd_arg) { 
  cout << endl << "Loading RNGs..." << endl << endl;

  io_good = false;
  int error = 0;

  if(isRoot()) {
    // all open file and check error
    ifstream input(rd_arg.FileName);
    if ( !input.good() ) {
      cout << "Could not open file:\n   "
	   << rd_arg.FileName
	   << "\nfor input.\n";
      error = 1;
    }

    if(!error) {
      hd.read(input);
      input.close();
    }
  }
  // first sync point
  // executed by all, sync and share error status information
  if(synchronize(error) != 0)   return;

  broadcastInt(&hd.data_start);


  if(isRoot()) {
    if(hd.datatype != "LATTICE_RNG_5D_4D"){  // need both 5d & 4d
      cout << "Invalid RNG type: " << hd.datatype << endl;
      error = 1;
    }
  }
  if(synchronize(error) != 0)  return;


  // check dimensions, etc
  int nx = rd_arg.Xnodes() * rd_arg.XnodeSites();
  int ny = rd_arg.Ynodes() * rd_arg.YnodeSites();
  int nz = rd_arg.Znodes() * rd_arg.ZnodeSites();
  int nt = rd_arg.Tnodes() * rd_arg.TnodeSites();
  int ns = rd_arg.Snodes() * rd_arg.SnodeSites();

  if(isRoot()) {
    cout << "Lattice Dimensions in File: " << hd.dimension[0] <<"x" << hd.dimension[1] <<"x" << hd.dimension[2] <<"x" << hd.dimension[3] <<"x" << hd.dimension[4]<< endl;
    if(hd.dimension[0] != nx || hd.dimension[1] != ny || hd.dimension[2] != nz || hd.dimension[3] != nt || hd.dimension[4]!=ns) {
      cout << "Dimensions in file DISAGREE with GlobalJobParameter!"<<endl;
      error = 1;
    }
  }

  if(synchronize(error) != 0)  return;

  // check floating point format
  if(isRoot()) {
    intconv.setFileFormat(hd.int_format);
  }
  broadcastInt((int*)&intconv.fileFormat);
  if(intconv.fileFormat == INT_UNKNOWN) {
    cout << "Data file Integer format UNKNOWN" << endl;
    return;
  }

  if(isRoot())  hd.show();

  // the UGrandomGenerator has multiple type members, need a universal
  // format to make it cross-platform
  int size_rng_ints = ugran[0].RNGints();
  int size_rng_chars = size_rng_ints * intconv.fileIntSize();

  QioArg rng_arg(rd_arg);
  rng_arg.cutHalf();

  unsigned int csum[2] = {0}, pos_dep_csum[2] = {0};
  Float RandSum[2]={0.0}, Rand2Sum[2]={0.0};

#if TARGET != QCDOC  // when on LINUX, only parallel (direct IO) mode is used
  setParallel();
#endif

  if(parIO()) {
    cout << endl << "================= Loading 5-D RNGs =================" << endl << endl;

    ParallelIO pario(rng_arg);
    if(! pario.load((char*)ugran, size_rng_ints, sizeof(UGrandomGenerator),
		    hd, intconv, 5,
		    &csum[0], &pos_dep_csum[0], &RandSum[0], &Rand2Sum[0])) return;

    hd.data_start += size_rng_chars * rng_arg.VolSites() * 
                     rng_arg.Snodes() * rng_arg.SnodeSites();
 
    cout << endl << "================= Loading 4-D RNGs =================" << endl << endl;
    if(! pario.load((char*)ugran_4d, size_rng_ints, sizeof(UGrandomGenerator),
		    hd, intconv, 4,
		    &csum[1], &pos_dep_csum[1], &RandSum[1], &Rand2Sum[1])) return;
  }
#if TARGET == QCDOC
  else {
    cout << endl << "================= Loading 5-D RNGs =================" << endl << endl;

    SerialIO serio(rng_arg);
    if(! serio.load((char*)ugran, size_rng_ints, sizeof(UGrandomGenerator),
		    hd, intconv, 5,
		    &csum[0], &pos_dep_csum[0], &RandSum[0], &Rand2Sum[0])) return;

    hd.data_start += size_rng_chars * rng_arg.VolSites() * 
                     rng_arg.Snodes() * rng_arg.SnodeSites();
    
    cout << endl << "================= Loading 4-D RNGs =================" << endl << endl;

    if(! serio.load((char*)ugran_4d, size_rng_ints, sizeof(UGrandomGenerator),
		    hd, intconv, 4,
		    &csum[1], &pos_dep_csum[1], &RandSum[1], &Rand2Sum[1])) return;
  }
#endif

  cout << endl << "================= Verification =================" << endl << endl;

  // Step 2.1: verify  checksum
  csum[0] += csum[1];
  csum[0] = globalSumUint(csum[0]);

  if(isRoot()) {
    if( hd.checksum != csum[0] ) {
      cout << "CheckSUM error !! Header:" << hex << hd.checksum << dec 
	   << " Host calc:" <<hex << csum[0] << dec << endl;
      error = 1;
    }
    else
      cout << "CheckSUM is ok\n";
  }

  if(synchronize(error) != 0) return;

  // Step 2.2: verify position-dep. Checksum
  pos_dep_csum[0] += pos_dep_csum[1];
  pos_dep_csum[0] = globalSumUint(pos_dep_csum[0]);
  if(isRoot()) {
    // pos_dep_csum could be absent
    if( hd.pos_dep_csum > 0 && hd.pos_dep_csum != pos_dep_csum[0] ) { 
      cout << "Position Dependent CheckSUM error!!" 
	   << " Header:" <<hex << hd.pos_dep_csum << dec
	   << " Host_calc:" <<hex << pos_dep_csum[0] << dec << endl;
      error = 1;
    }
    else
      cout << "Position Dependent CheckSUM is ok\n";
  }

  if(synchronize(error) != 0) return;


  // STEP 3: Verify Rand Average and Variance
  RandSum[0] += RandSum[1];
  Rand2Sum[0] += Rand2Sum[1];
  int total_rngs_4d = rng_arg.VolSites();
  int total_rngs_5d = total_rngs_4d * rng_arg.Snodes() * rng_arg.SnodeSites(); 
  Float RandAvg = globalSumFloat(RandSum[0]) / (total_rngs_5d + total_rngs_4d);
  Float RandVar = globalSumFloat(Rand2Sum[0]) / (total_rngs_5d + total_rngs_4d)
                  - RandAvg * RandAvg;
  if(isRoot()) {
  cout << "Average::  calc: " << RandAvg << "  header: " << hd.average 
       << "  rel.dev.: " << fabs((RandAvg-hd.average)/RandAvg) <<  endl;
  cout << "Variance:: calc: " << RandVar << "  header: " << hd.variance 
       << "  rel.dev.: " << fabs((RandVar-hd.variance)/RandVar) << endl;
  }

  cout << endl << "================= Loading Complete =================" << endl << endl;

  io_good = true;

}


///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// LatRngWrite functions

LatRngWrite::LatRngWrite() 
  : LatRngIO()
{ }

LatRngWrite::~LatRngWrite()
{ }



void LatRngWrite::write(UGrandomGenerator * ugran, UGrandomGenerator * ugran_4d,
			const QioArg & wt_arg) {
  cout << "LatRngWrite::write() start..." << endl;

  io_good = false;
  int error = 0;

  // some cross-platform issued to be discussed
  intconv.setFileFormat(wt_arg.FileIntFormat);
  cout << "Set File INT format : " << intconv.name(intconv.fileFormat) << endl;

  // size information
  int nx = wt_arg.Xnodes() * wt_arg.XnodeSites();
  int ny = wt_arg.Ynodes() * wt_arg.YnodeSites();
  int nz = wt_arg.Znodes() * wt_arg.ZnodeSites();
  int nt = wt_arg.Tnodes() * wt_arg.TnodeSites();
  int ns = wt_arg.Snodes() * wt_arg.SnodeSites();

  cout << "Lattice size = " << nx <<"x"<<ny<<"x"<<nz<<"x"<<nt<<"x"<<ns<<endl;

  // num of integers per RNG (NOTE: not all data in RNG are saved)
  int size_rng_ints = ugran[0].RNGints();
  cout << "size_rng_ints = " << size_rng_ints << endl;
  int size_rng_chars = size_rng_ints * intconv.fileIntSize();

#if TARGET != QCDOC
  setParallel();
#endif

  // start writing file
  ofstream output;

  if(parIO()) {
    output.open(wt_arg.FileName);
    if(!output.good())
      {
	cout << "Could not open file:\n   "
	     << wt_arg.FileName
	     << "\nfor output.\n";
	error = 1;
      }
  }
#if TARGET == QCDOC
  else {
    if(isRoot()) {
      output.open(wt_arg.FileName);
      if(!output.good())
	{
	  cout << "Could not open file:\n   "
	       << wt_arg.FileName
	       << "\nfor output.\n";
	  error = 1;
	}
    }    
  }  
#endif
  if(synchronize(error) > 0)  return;

  // write header
  if(isRoot()) {
    hd.init(wt_arg, intconv.fileFormat);
    hd.write(output);
  }
  broadcastInt(&hd.data_start); // from 0 to all

  unsigned int csum[2]={0}, pos_dep_csum[2] = {0};
  Float RandSum[2]={0.0}, Rand2Sum[2] = {0.0};

  QioArg rng_arg(wt_arg);
  rng_arg.cutHalf();

  if(parIO()) {
    cout << endl << "================= Unloading 5-D RNGs =================" << endl << endl;

    ParallelIO pario(rng_arg);
    if(! pario.store(output, (char*)ugran, size_rng_ints,
		     sizeof(UGrandomGenerator), hd, intconv, 5,
		     &csum[0], &pos_dep_csum[0], &RandSum[0], &Rand2Sum[0])) return;

    hd.data_start += size_rng_chars * rng_arg.VolSites() * 
                     rng_arg.Snodes() * rng_arg.SnodeSites();
 
    cout << endl << "================= Unloading 4-D RNGs =================" << endl << endl;

    if(! pario.store(output, (char*)ugran_4d, size_rng_ints, 
		     sizeof(UGrandomGenerator), hd, intconv, 4,
		     &csum[1], &pos_dep_csum[1], &RandSum[1], &Rand2Sum[1])) return;
  }
#if TARGET == QCDOC
  else {
    cout << endl << "================= Unloading 5-D RNGs =================" << endl << endl;

    SerialIO serio(rng_arg);
    if(! serio.store(output, (char*)ugran, size_rng_ints, 
		     sizeof(UGrandomGenerator), hd, intconv, 5,
		     &csum[0], &pos_dep_csum[0], &RandSum[0], &Rand2Sum[0])) return;

    hd.data_start += size_rng_chars * rng_arg.VolSites() * 
                     rng_arg.Snodes() * rng_arg.SnodeSites();

    cout << endl << "================= Unloading 4-D RNGs =================" << endl << endl;
    
    if(! serio.store(output, (char*)ugran_4d, size_rng_ints, 
		     sizeof(UGrandomGenerator), hd, intconv, 4,
		     &csum[1], &pos_dep_csum[1], &RandSum[1], &Rand2Sum[1])) return;
  }
#endif


  // fill in verification information
  csum[0] += csum[1];
  csum[0] = globalSumUint(csum[0]);
  pos_dep_csum[0] += pos_dep_csum[1];
  pos_dep_csum[0] = globalSumUint(pos_dep_csum[0]);

  RandSum[0] += RandSum[1];
  Rand2Sum[0] += Rand2Sum[1];
  int total_rngs_4d = rng_arg.VolSites();
  int total_rngs_5d = total_rngs_4d * rng_arg.Snodes() * rng_arg.SnodeSites(); 
  Float RandAvg = globalSumFloat(RandSum[0]) / (total_rngs_5d + total_rngs_4d);
  Float RandVar = globalSumFloat(Rand2Sum[0]) / (total_rngs_5d + total_rngs_4d)
                  - RandAvg * RandAvg;

  if(isRoot()) {
    hd.fillInCheckInfo(output, csum[0], pos_dep_csum[0], RandAvg, RandVar);
    if ( !output.good() ) { cout << "Output stream error!" << endl; error = 1; }
  }

  if(synchronize(error) != 0)  return;

  if(parIO()) 
    output.close();
  else
    if(isRoot()) output.close();
  
  cout << endl << "================= Unloading Complete =================" << endl << endl;

  io_good = true;
}



CPS_END_NAMESPACE
 
