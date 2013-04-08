#include <config.h>
#include <util/WriteU1LatticePar.h>
#include <util/iostyle.h>

using namespace std;
CPS_START_NAMESPACE

void WriteU1LatticeParallel::write(Lattice & lat, const QioArg & wt_arg)
{
  const char * fname = "write()";
  VRB.Func(cname,fname);

  char loginfo[100];
  sprintf(loginfo,"Unload %s",wt_arg.FileName);
  startLogging(loginfo);

#ifdef PROFILE
  struct timeval start,end;
  gettimeofday(&start,NULL);
#endif
   sync();

  // init
  int error = 0;
  io_good = false;

  // always zero for U1:
  recon_row_3 = 0;

  FP_FORMAT dataFormat = wt_arg.FileFpFormat;
  fpconv.setFileFormat(dataFormat);
  if(fpconv.fileFormat == FP_UNKNOWN) {
    ERR.General(cname,fname,"Output Floating Point format UNKNOWN\n");
  }

  // calc Plaq and LinkTrace
  const int size_matrices(wt_arg.VolNodeSites()*4);
  Float * lpoint = lat.U1GaugeField();
  VRB.Flow(cname,fname, "Writing U1 Gauge Field at Lattice::U1GaugeField() = %p\n", lpoint);

  Float plaq = lat.SumReU1Plaq()/ 6.0 / (wt_arg.VolSites()) ;
  //printf("U1 plaq= %e\n", plaq);fflush(stdout); exit(0);
  Float ltrace(0.0);
  if(wt_arg.Scoor() == 0) {
    for(int i=0;i<size_matrices;i++){
      ltrace += *(lpoint+i);
    }
    ltrace = globalSumFloat(ltrace) / (4*wt_arg.VolSites());
  }
  else
    globalSumFloat(0.0);  // everyone has to participate in global ops

  log();
  
  // write lattice data, in Parallel or Serial manner
  // determined by the template parameter "IoStyle" of this class
  int data_per_site = 4;
  const int chars_per_site = data_per_site * fpconv.fileFpSize();

  unsigned int csum = 0;

#if TARGET != QCDOC   // when not on QCDOC(like on LINUX), use serial IO mode
  setSerial();
#endif
  
  fstream output;

  if (isRoot()){
    output.open(wt_arg.FileName,ios::out);
	output.close();
  }
  Float temp=0.;
  glb_sum(&temp);
  sync();
  if(parIO()) {
    // all open file, start writing
    output.open(wt_arg.FileName);
    if(!output.good())    {
      //      VRB.Flow(cname,fname, "Could not open file: [%s] for output.\n",wt_arg.FileName);
      //      VRB.Flow(cname,fname,"USER: maybe you should kill the process\n");
      
      printf("Node %d:Could not open file: [%s] for output.\n",UniqueID(),wt_arg.FileName);
      error = 1;
    }
  }
  else {
    // only node 0 open file, start writing
    if(isRoot()) {
      FILE *fp = fopen(wt_arg.FileName,"w");
      fclose(fp);
      output.open(wt_arg.FileName);
      if(!output.good())    {
	//	VRB.Flow(cname,fname, "Could not open file: [%s] for output.\n",wt_arg.FileName);
	//	VRB.Flow(cname,fname,"USER: maybe you should kill the process\n");
	error = 1;
      }
    }
  }
  if (error)
    printf("Node %d: says opening %s failed\n",UniqueID(),wt_arg.FileName);
//  if(synchronize(error) > 0)  
//    ERR.FileW(cname,fname,wt_arg.FileName);
   error=0;

  // write header
  if(isRoot()){
    hd.init(wt_arg, fpconv.fileFormat, ltrace, plaq);
    hd.write(output);
  }
  if (error)
    printf("Node %d: says Writing header failed  %s failed\n",UniqueID(),wt_arg.FileName);
//  if(synchronize(error) > 0) 
//    ERR.General(cname,fname,"Writing header failed\n");
  log();

{
  int temp_start = hd.data_start;
  broadcastInt(&temp_start);
  hd.data_start = temp_start;
}

  if(parIO()) {
    ParallelIO pario(wt_arg);
    if(!pario.store(output, (char*)lpoint, data_per_site, sizeof(Float)*4,
		    hd, fpconv, 4, &csum)) 
      ERR.General(cname, fname, "Unload failed\n");

//    printf("Node %d: lattice write csum=%x\n",UniqueID(),csum);
    if(wt_arg.Scoor() == 0) 
      csum = globalSumUint(csum);
    else 
      globalSumUint(0);
  }
#if 1
  else {
    SerialIO serio(wt_arg);
    if(!serio.store(output, (char*)lpoint, data_per_site, sizeof(Float)*4,
		    hd, fpconv, 4, &csum)) 
      ERR.General(cname, fname, "Unload failed\n");
  }
#endif
//  output.flush();

  log();

  // after unloading, fill in checksum
  if(isRoot()) 
    hd.fillInChecksum(output,csum);

  if(parIO()) {
    output.close();
  }
  else {
    if(isRoot())
      output.close();
  }
  if(synchronize(error) != 0)  
    ERR.General(cname, fname, "Writing checksum failed\n");
  
#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(cname,"write",0,&start,&end);
#endif

  io_good = true;

  log();
  finishLogging();
  sync();
  VRB.FuncEnd(cname,fname);
}




CPS_END_NAMESPACE
