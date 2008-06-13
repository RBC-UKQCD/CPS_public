#include <config.h>
#include <stdlib.h>
#include <util/gjp.h>
#include <alg/common_arg.h>
#include <alg/alg_base.h>
#include <alg/alg_ghb.h>
#include <alg/do_arg.h>
#include <alg/ghb_arg.h>
//#include <sysfunc_cps.h>
#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>
#include <util/fpconv.h>
#include <comms/sysfunc_cps.h>


USING_NAMESPACE_CPS



int read_lattice(int argc, char ** argv) {
  cout << "ReadLatticeParallel class test!" << endl;

  GwilsonFnone lat;  // create a lattice at the GJP-specified address

  
  // load lattice, and use the above created lattice to check plaq and link trace
  ReadLatticeParallel  rd;
  rd.setLogDir("logs");
  //  rd.SimQCDSP(1);
  rd.read(lat,argv[2]);

  //// to have more control, use
  //
  // QioArg  rd_arg(char * filename, Float check_precision /* optional */ );
  // 
  //   (modify rd_arg members....)   // including setting ConcurIONumber // see  qioarg.h
  //
  // ReadLatticeParallel rd(lat, rd_arg);
  //

  
  if(!rd.good()) {
    cout << "Loading failed" << endl;
    exit(-13);
  }
  

  cout << "Load complete" << endl << endl;

  cout << "================== PARALLEL IO ======================" << endl;
  cout << "-----------------------------------------------------" << endl;
  cout << "================== SERIAL   IO ======================" << endl << endl;

  ReadLatticeSerial rd2;
  rd2.setLogDir("logs");
  //  rd2.SimQCDSP(1);
  rd2.read(lat,argv[2]);

  if(!rd2.good()) {
    cout << "Serial Loading failed" << endl;
    exit(-13);
  }

  cout << "Load Complete" << endl;

  return 0;
}


int write_lattice(int argc, char ** argv) {
  cout << "WriteLatticeParallel class test!" << endl;

  // generate lattice
  GwilsonFnone lat;
  int ITERATIONS = 10;
  
  CommonArg common_arg_ghb;
  common_arg_ghb.results = (void*)"ghb.dat";
  GhbArg ghb_arg;
  ghb_arg.num_iter=3;
  AlgGheatBath ghb(lat,&common_arg_ghb,&ghb_arg);

  for(int i=0;i<ITERATIONS;++i)  ghb.run();

  ////////////////////////////////////////////////////////////////////
  // start unloading lattice
  //WriteLatticeParallel wt(lat, argv[2]); 

  // more options, e.g.
  //   WriteLatticeParallel wt(lat, argv[2], FP_IEEE32BIG, 1);//  for NERSC format, see util/fpconv.h


  //// to have even more control, use
  //
  // QioArg  wt_arg(char * filename, FP_FORMAT fileFpFormat /* optional */, int recon_row_3 /* optional */ );
  // 
  //   (modify wt_arg members....)  // including setting ConcurIONumber // see  qioarg.h
  //
  // WriteLatticeParallel wt(lat, wt_arg);
  //


  // TWO-STEP writing, to set some labels in header
  WriteLatticeParallel wt;
  wt.setHeader("EnsTest","Testing Ensemble by Sam Oct,2004",1001);
  wt.setLogDir("logs");
  wt.write(lat,argv[2]);

  if(!wt.good()) {
    cout << "Unloading failed" << endl;
    exit(-13);
  }

  cout << "===================  PARALLEL UNLOADING ["<<argv[2]<<"] =========================" << endl;
  cout << "--------------------------------------------------------------------------------------" << endl;
  cout << "===================  SERIAL   UNLOADING ["<<argv[2]<<".serial] ==================" << endl << endl;

  WriteLatticeSerial wt2;
  wt2.setHeader("EnsSerial", "Serial unloading test Nov, 2004", 1002);
  wt2.setLogDir("logs");

  char filename2[256];
  strcpy(filename2,argv[2]);
  strcat(filename2,".serial");

  wt2.write(lat,filename2);
  if(!wt2.good()) {
    cout << "Unloading failed" << endl;
    exit(-13);
  }

  cout << "Unload complete" <<endl << endl;

  return 0;
}


int main(int argc, char ** argv) {
  if(argc<11) {
    cout << "Usage:" << endl<<"      qrun QCDOC.x  -[r|w]  <conf.dat>  <x sites> <y sites> <z sites> <t sites>  <Xbc> <Ybc> <Zbc> <Tbc>"<< endl;
    cout << "(use letter \'P\' or \'A\' for arguments of gauge BC's)" << endl;
    cout << "Eg,   qrun QCDOC.x -r  conf8x8x8x16.file   8 8 8 16  P P P P"<< endl;
    cout << "      qrun QCDOC.x -w  conf4x4x4x32.file   4 4 4 32  P P A A"<< endl;
    exit(1);
  }


  Start(&argc,&argv);
  // init  GJP
  DoArg do_arg;

  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = SizeS();

  int nx = atoi(argv[3]);
  int ny = atoi(argv[4]);
  int nz = atoi(argv[5]);
  int nt = atoi(argv[6]);

  do_arg.x_node_sites = nx/do_arg.x_nodes;
  do_arg.y_node_sites = ny/do_arg.y_nodes;
  do_arg.z_node_sites = nz/do_arg.z_nodes;
  do_arg.t_node_sites = nt/do_arg.t_nodes;
  do_arg.s_node_sites = 1;

  do_arg.x_bc = (argv[7][0]=='A' ? BND_CND_APRD : BND_CND_PRD);
  do_arg.y_bc = (argv[8][0]=='A' ? BND_CND_APRD : BND_CND_PRD);
  do_arg.z_bc = (argv[9][0]=='A' ? BND_CND_APRD : BND_CND_PRD);
  do_arg.t_bc = (argv[10][0]=='A' ? BND_CND_APRD : BND_CND_PRD);

  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.beta = 5.3;
  do_arg.dwf_height = 0.9;

  // start testing
  if(!strcmp(argv[1],"-w")) {
    do_arg.start_conf_kind = START_CONF_DISORD;
    GJP.Initialize(do_arg);

    cout << "Initialized ok" << endl;

    write_lattice(argc,argv);
  }
  else {
    do_arg.start_conf_load_addr =
      (unsigned long)smalloc(do_arg.x_node_sites * do_arg.y_node_sites *
		       do_arg.z_node_sites * do_arg.t_node_sites * 
		       4 * sizeof(Matrix));
#if 1
    do_arg.start_conf_kind = START_CONF_LOAD;
    GJP.Initialize(do_arg);

    cout << "Initialized ok" << endl;


  #if TARGET == QCDOC
    // *Shift() functions test
    QioArg qarg("NoFile");
    SerialIO  serio(qarg);
    if(serio.backForthTest()) 
      cout << "Back-Forth test done!" << endl;      
    else
      cout << "Back-Forth test failed!!" << endl;
    if(serio.rotateTest()) { 
      cout << "Rotation test done!" << endl;
    }
    else {
      cout << "Rotation test failed!!" << endl;
    }


    cout << "========================================================================"<<endl;

  #endif



    read_lattice(argc,argv);
#else
//   An equivalent way to load the lattice. ReadLatticePar is called
//   inside Lattice::Lattice()

    do_arg.start_conf_kind = START_CONF_FILE;
    do_arg.start_conf_filename  = argv[2];
    GJP.Initialize(do_arg);
    cout << "Initialized ok" << endl;
    GwilsonFnone lat;
    cout << "lattice loaded ok" << endl;
    const char *write_file_name = "test_out.lat";
//  Should dupilcate the NERSC format lattice except header
    WriteLatticeParallel(lat,write_file_name,FP_IEEE32BIG,1);
#endif

    sfree((Matrix *)do_arg.start_conf_load_addr);
  }


  exit(0);
}

  
