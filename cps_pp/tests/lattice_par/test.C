#include <util/gjp.h>
#include <alg/common_arg.h>
#include <alg/alg_base.h>
#include <alg/alg_ghb.h>
#include <alg/do_arg.h>
#include <alg/ghb_arg.h>
#include <sysfunc.h>
#include <qmp.h>
#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>
#include <util/fpconv.h>

USING_NAMESPACE_CPS


int read_lattice(int argc, char ** argv) {
  cout << "ReadLatticeParallel class test!" << endl;

  // load lattice
  ReadLatticeParallel  rd;
  rd.setConcurIONumber(4); // num of processors to do concurrent IO, any number, default is 8
  rd.read(argv[2]);

  cout << "Load complete" << endl << endl;

  // check its link trace and plaq trace

  // now re-set GJP from gauge file (to extract gauge field addr, etc)
  GJP.Initialize(rd.do_arg);

  GwilsonFnone latNodeLoad;
  latNodeLoad.GaugeField(rd.GaugeField());

  cout << endl << "Checking lattice data" << endl;

  rd.CheckPlaqLinktrace(latNodeLoad, 0.1);

  cout << "Verification done" << endl << endl;

  return 0;
}


int write_lattice(int argc, char ** argv) {
  cout << "WriteLatticeParallel class test!" << endl;

  // generate lattice
  GwilsonFnone lat;
  int ITERATIONS = 10;
  
  CommonArg common_arg_ghb;
  common_arg_ghb.results = (void*)"/host/lishu/qcdio/ghb.dat";
  GhbArg ghb_arg;
  ghb_arg.num_iter=3;
  AlgGheatBath ghb(lat,&common_arg_ghb,&ghb_arg);

  for(int i=0;i<ITERATIONS;++i)  ghb.run();

  // unload lattice
  WriteLatticeParallel  wt;

  wt.setConcurIONumber(4);
  if(argc < 6) {
    wt.write(lat,argv[2]); // default: output format = host format, output 3 rows
  }
  else { // specified output floating format
    FPConv fp;
    fp.setFileFormat(argv[5]);
    if(fp.fileFormat == FP_UNKNOWN)  fp.fileFormat = FP_AUTOMATIC;

    wt.write(lat,argv[2],fp.fileFormat,1);  // output format = user set, output 2 rows (recon_row_3 == 1)
  }

  cout << "Unload complete" <<endl << endl;

  return 0;
}


int main(int argc, char ** argv) {
  if(argc<5) {
    cout << "Usage:" << endl<<"      qrun QCDOC.x  -[r|w]  <conf.dat>  <x/y/z sites>  <t sites> "<< endl;
    exit(1);
  }

  QMP_init_msg_passing(&argc, &argv, QMP_SMP_ONE_ADDRESS);

  cout << "QMP init'd" << endl;

  // init  GJP
  DoArg do_arg;
  ParallelControl pc;

  do_arg.x_nodes = pc.Xnodes();
  do_arg.y_nodes = pc.Ynodes();
  do_arg.z_nodes = pc.Znodes();
  do_arg.t_nodes = pc.Tnodes();
  do_arg.s_nodes = 1;

  int nx = atoi(argv[3]);
  int ny = atoi(argv[3]);
  int nz = atoi(argv[3]);
  int nt = atoi(argv[4]);

  do_arg.x_node_sites = nx/do_arg.x_nodes;
  do_arg.y_node_sites = ny/do_arg.y_nodes;
  do_arg.z_node_sites = nz/do_arg.z_nodes;
  do_arg.t_node_sites = nt/do_arg.t_nodes;
  do_arg.s_node_sites = 1;

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  //  do_arg.t_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
//  do_arg.start_conf_kind = START_CONF_ORD;
//  do_arg.start_conf_kind = START_CONF_LOAD;
//  do_arg.start_conf_load_addr = (Matrix *)0x5f700;
  do_arg.start_seed_kind = START_SEED_FIXED;
  //  do_arg.colors = 3;
  do_arg.beta = 5.3;
  do_arg.dwf_height = 0.9;
  //  do_arg.verbose_level = 1;  //-100502;
//  do_arg.verbose_level = 100;
//  do_arg.verbose_level = 0;

  GJP.Initialize(do_arg);

  cout << "Initialized ok" << endl;


  // start testing
  if(!strcmp(argv[1],"-w")) {
    write_lattice(argc,argv);
  }
  else {
    read_lattice(argc,argv);
  }

  QMP_finalize_msg_passing();

  exit(0);
}

  
