#include <config.h>
#include <util/latrngio.h>
#include <util/gjp.h>
#include <alg/do_arg.h>
#include <comms/sysfunc_cps.h>
#include <stdlib.h>


USING_NAMESPACE_CPS


int main(int argc, char ** argv) {
  if(argc<6) {
    cout << "Usage:" << endl<<"      qrun QCDOC.x  -[r|w]  <RNG.file>  <x/y/z sites>  <t sites>  <s sites>"<< endl;
    cout << "Eg,   qrun QCDOC.x -r  rng8x8x8x8x4.file   8  8  4"<< endl;
    cout << "      qrun QCDOC.x -w  rng4x4x4x32x32.file 4  32 32"<< endl;
    exit(1);
  }


  // init  GJP
  DoArg do_arg;

  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = SizeS();

  int nx = atoi(argv[3]);
  int ny = atoi(argv[3]);
  int nz = atoi(argv[3]);
  int nt = atoi(argv[4]);
  int ns = atoi(argv[5]);

  do_arg.x_node_sites = nx/do_arg.x_nodes;
  do_arg.y_node_sites = ny/do_arg.y_nodes;
  do_arg.z_node_sites = nz/do_arg.z_nodes;
  do_arg.t_node_sites = nt/do_arg.t_nodes;
  do_arg.s_node_sites = ns/do_arg.s_nodes;

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;

  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.beta = 5.3;
  do_arg.dwf_height = 0.9;

  do_arg.start_conf_kind = START_CONF_DISORD;
  GJP.Initialize(do_arg);

  cout << "GJP initialized ok" << endl;

  // both write or read, need to init LRG first
  LRG.Initialize(); 

  srand(time(NULL));
  int shuffletime = rand()%100+1;
  cout << "Shuffle for " << shuffletime << " times" << endl;
  
  for(int s=0;s<shuffletime;s++) {
    for(int x=0;x<GJP.XnodeSites();x++)
      for(int y=0;y<GJP.YnodeSites();y++)
	for(int z=0;z<GJP.ZnodeSites();z++)
	  for(int t=0;t<GJP.TnodeSites();t++)
	    for(int s=0;s<GJP.SnodeSites();s++) {
	      LRG.AssignGenerator(x,y,z,t,s);
	      LRG.Grand(FIVE_D);
	      LRG.Grand(FOUR_D);
	    }
    cout << "Shuffling #" << s+1 << " done!" << endl;
  }
  
  cout << "Shuffling done!"<< endl;

  LRG.setLogDir("logs");

  if(!strcmp(argv[1],"-r")) { // read
    cout << endl << "***********************  PARALLEL MODE  **************************" << endl;
    LRG.setParallel();
    if(LRG.Read(argv[2])) {
      cout << "Parallel Reading from file [" << argv[2] << "] success!" << endl;
    }
    else {
      cout << "Parallel Reading from file [" << argv[2] << "] failed!" << endl;
    }

    cout << endl << "***********************  SERIAL MODE  **************************" << endl;
    LRG.setSerial();
    if(LRG.Read(argv[2])) {
      cout << "Serial Reading from file [" << argv[2] << "] success!" << endl;
    }
    else {
      cout << "Serial Reading from file [" << argv[2] << "] failed!" << endl;
    }
  }
  else if(!strcmp(argv[1],"-w")) {  // write
    cout << endl << "***********************  PARALLEL MODE  **************************" << endl;
    LRG.setParallel();
    if(LRG.Write(argv[2])) {
      cout << "Parallel Writing to file [" << argv[2] << "] success!" << endl;
    }
    else {
      cout << "Parallel Writing to file [" << argv[2] << "] failed!" << endl;
    }

    char fname[256];
    strcpy(fname,argv[2]);
    strcat(fname,".serial");

    cout << endl << "***********************  SERIAL MODE  **************************" << endl;
    LRG.setSerial();
    if(LRG.Write(fname)) {
      cout << "Serial Writing to file [" << fname << "] success!" << endl;
    }
    else {
      cout << "Serial Writing to file [" << fname << "] failed!" << endl;
    }
  }
  else {
    cout << "Parameter format wrong!" << endl;
    cout << "For help, run this program w/o params"<< endl;
    exit(1);
  }


  exit(0);
}

  
