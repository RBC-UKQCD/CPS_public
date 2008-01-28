#include <stdio.h>
#include <stdlib.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <alg/alg_fix_gauge.h>
#include <alg/alg_plaq.h>
#include <alg/alg_nuc3pt.h>
#include <util/command_line.h>
#include <util/ReadLatticePar.h>

USING_NAMESPACE_CPS

int main(int argc,char *argv[])
{
  char *cname = argv[0] ;
  char *fname = "main()" ;
  char *filename;
  filename = (char*) smalloc( 128*sizeof(char) );

  CommandLine::is(argc,argv);

  int serial_io = 0;
  int concur_io_number = 64;
  //------------------------------
  //set log directory
  //------------------------------
  char log_dir[255];
  sprintf(&log_dir[0],"IOLOG");

  // set lattice IO parameters
  //-------------------------------------------------
  //WriteLatticeParallel wl;
  //wl.setLogDir(&log_dir[0]);
  //wl.setHeader(ensemble_id,ensemble_label,sequence);
  

  // load defaults
  DoArg do_arg;
  // working directory
  const char* dir(CommandLine::arg());
  strcpy(filename,dir);
  strcat(filename,"/do_arg.default");
 
  do_arg.Decode(filename,"do_arg");
  // want to load in the file
  //do_arg.start_conf_kind       = START_CONF_FILE;
  //do_arg.start_seed_kind       = START_SEED_FIXED;
  // with this name
  // print out the inputs
  strcpy(filename, dir);
  strcat(filename,"/do_arg.dat");
  do_arg.Encode(filename,"do_arg");

  int conf_flag=0;
  if( do_arg.start_conf_kind == START_CONF_FILE ){
    conf_flag=1;
    do_arg.start_conf_kind = START_CONF_ORD;
  }

  GJP.Initialize(do_arg);


  // set the verbose level
  VRB.Level(2);

  GwilsonFdwf lattice;
  // enable the link buffer
  //lattice.EnableLinkBuffer(1000);
  
  const char* conf_file(CommandLine::arg());
  if( conf_flag ) do_arg.start_conf_kind = START_CONF_FILE;

  NoArg no_arg;
  CommonArg common_arg;

  int loop_s = CommandLine::arg_as_int();
  int loop_e = CommandLine::arg_as_int();
  int loop_skp = CommandLine::arg_as_int();


  // conf. loop start
  for(int loop=loop_s; loop<loop_e; loop+=loop_skp){

  sprintf(filename,"%s/nuc3pt.dat.%d",dir,loop);
  common_arg.set_filename(filename);
  FILE *fp;
  if( (fp = Fopen((char *)common_arg.results, "a")) == NULL ) {
    ERR.FileA(cname, fname, (char *)common_arg.results );
  }

#if TARGET == QCDOC
  Fprintf( fp, "Machine size in nodes:  (X,Y,Z,T) = (%d,%d,%d,%d)\n",
    SizeX(), SizeY(), SizeZ(), SizeT() );
#else
  Fprintf( fp, "Scalar Machine\n" ) ;
#endif
  Fprintf( fp, "\nCommand line arguments\n");
  Fprintf( fp, "\tprogram name:\t%s\n", argv[0] );

  if(do_arg.start_conf_kind == START_CONF_FILE){
     
    ReadLatticeParallel rd;
    rd.setLogDir(&log_dir[0]);
  
    if(serial_io) {
      //LRG.setSerial();
      rd.setSerial();
      //wl.setSerial();
    }

    char confname[100];
    sprintf(confname,"%s.%d",conf_file,loop);
    printf("%s %s\n",conf_file,confname);
    QioArg rd_arg(confname,1e-6);
    rd_arg.ConcurIONumber=concur_io_number;
    if(serial_io) rd.setSerial();

    cout<<"Start loading lattice"<<endl;
    rd.read(lattice,rd_arg);
    if(!rd.good()) ERR.General(cname,fname,"Failed to load lattice\n");
    cout<<"Loading latice finished"<<endl;
    Fprintf( fp, "\tlattice filename:\t%s\n",confname );
  }

  Fclose(fp);//output file

  AlgPlaq     plaq   (lattice,&common_arg,&no_arg);
  plaq.run();

  // args for the nucleon correlation functions:
  Nuc3ptArg nuc3pt_arg;
  // load defaults
  strcpy(filename, dir);
  strcat(filename,"/nuc3pt_arg.default");
  nuc3pt_arg.Decode(filename,"nuc3pt_arg");
  nuc3pt_arg.ensemble_id = loop;
  //we cannot calculate conserved current with 4D prop.
  if(nuc3pt_arg.calc_QProp==0) nuc3pt_arg.DoConserved=0;
  sprintf(filename,"%s/nuc3pt_arg.dat.%d",dir,loop);
  nuc3pt_arg.Encode(filename,"nuc3pt_arg");

  // compute gauge-fixing matrices for use with momentum source
  if(nuc3pt_arg.src_type == ( SourceType ) 3 ){ // BOX = 3
    FixGaugeArg fix_arg;
    fix_arg.fix_gauge_kind = ( FixGaugeType )  3; // Coulomb time slices
    fix_arg.hyperplane_start = nuc3pt_arg.t_source;
    fix_arg.hyperplane_step =  nuc3pt_arg.source_inc;
    fix_arg.hyperplane_num =  nuc3pt_arg.num_src;
    fix_arg.max_iter_num =  20000;
    fix_arg.stop_cond    =  1e-8;
    /* FIX the gauge */
    if(fix_arg.fix_gauge_kind != FIX_GAUGE_NONE){
       AlgFixGauge fix_gauge(lattice,&common_arg,&fix_arg);
       fix_gauge.run();
    }
    if( (fp = Fopen((char *)common_arg.results, "a")) == NULL ) {
      ERR.FileA(cname, fname, (char *)common_arg.results );
    }
    Fprintf( fp, "Coulomb gauge fixing matrices computed\n");
    Fprintf( fp, "start = %d  step = %d num = %d\n", 
             fix_arg.hyperplane_start,
             fix_arg.hyperplane_step,
             fix_arg.hyperplane_num);
    Fprintf( fp, "max it = %d  tol = %g\n", 
             fix_arg.max_iter_num,
             fix_arg.stop_cond);
    Fclose(fp);//output file
  } 

  // calculate the nucleon correlation functions
  AlgNuc3pt alg_nuc3pt(lattice,&common_arg,&nuc3pt_arg);
  alg_nuc3pt.run() ;

  } // loop end

  return 0;
} 








