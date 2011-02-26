
#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<sys/stat.h> //mkdir()
#include<sys/types.h>
#include<string.h>
#include<util/lattice.h>

#include<alg/common_arg.h>
#include<alg/do_arg.h>
#include<alg/no_arg.h>
#include<alg/array_arg.h> //FloatArray

#include<alg/meas_arg.h>
#include<alg/alg_meas.h>
#include<alg/alg_w_spect.h>
#include<alg/alg_plaq.h>
#include<alg/alg_pbp.h>
#include<alg/alg_fix_gauge.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>
#include<util/qioarg.h>

#undef QCDOC
#ifdef QCDOC
#include<qcdocos/gint.h> //SynchMach()
#endif

#define NUC3PT
#ifdef NUC3PT
#include<alg/alg_nuc3pt.h>
#endif

#undef USE_SCU_CHECKSUMS
#ifdef USE_SCU_CHECKSUMS
#include <qcdocos/scu_checksum.h>
#endif
#define MAX_FILENAME 256
//--------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

void CheckDirectory(char *);
void ReadGauge(MeasArg *meas);
void WriteGauge(MeasArg *meas);

void ReadRNG(MeasArg *meas);
void WriteRNG(MeasArg *meas);
void checkpoint(MeasArg *meas);
void run(Lattice &lat, MeasArg &meas_arg);

DoArg do_arg;
MeasArg meas_arg;
MeasArg meas_arg_save;
WspectArg w_spect_arg;
CgArg w_spect_cg_arg;

int main(int argc, char *argv[])
{ 
  char *cname=argv[0];
  char *fname="main()";

  Start(&argc, &argv);

  VRB.Result(cname,fname,"Starting the measurements.\n");

  if ( argc!=7) { 
    printf("Args: doarg-file measarg-file  initial-directory enable_reproduce repro_percent repro_offset\n");
    exit(-1);
  }
 
  chdir (argv[3]);

  if ( !do_arg.Decode(argv[1],"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}
  if ( !meas_arg.Decode(argv[2],"meas_arg")){printf("Bum meas_arg\n"); exit(-1);}
  int reproduce = atoi(argv[4]);
  double repro_percent = atof(argv[5]);
  int repro_offset = atoi(argv[6]);

  printf("Mapping the machine topology\n");

  /*Layout the lattice on the machine (without regard to even-odd)*/
  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = SizeS();
  VRB.Result(cname,fname,"Logical machine topology: %d x %d x %d x %d x %d\n",do_arg.x_nodes, do_arg.y_nodes, do_arg.z_nodes, do_arg.t_nodes, do_arg.s_nodes);

  do_arg.x_node_sites = do_arg.x_sites/do_arg.x_nodes; 
  do_arg.y_node_sites = do_arg.y_sites/do_arg.y_nodes;
  do_arg.z_node_sites = do_arg.z_sites/do_arg.z_nodes;
  do_arg.t_node_sites = do_arg.t_sites/do_arg.t_nodes;
  do_arg.s_node_sites = do_arg.s_sites/do_arg.s_nodes;

  if (do_arg.x_sites!=do_arg.x_node_sites*do_arg.x_nodes) {printf("Lattice does not fit\n");exit(-1);}
  if (do_arg.y_sites!=do_arg.y_node_sites*do_arg.y_nodes) {printf("Lattice does not fit\n");exit(-1);}
  if (do_arg.z_sites!=do_arg.z_node_sites*do_arg.z_nodes) {printf("Lattice does not fit\n");exit(-1);}
  if (do_arg.t_sites!=do_arg.t_node_sites*do_arg.t_nodes) {printf("Lattice does not fit\n");exit(-1);}
  if (do_arg.s_sites!=do_arg.s_node_sites*do_arg.s_nodes) {printf("Lattice does not fit\n");exit(-1);}

    //  chdir (meas_arg.WorkDirectory);
#ifdef USE_SCU_CHECKSUMS
  ScuChecksum::Initialise(meas_arg.HdwXCsum,meas_arg.HdwRCsum);
#endif

  //  do_arg.verbose_level=VERBOSE_RESULT_LEVEL;
  GJP.Initialize(do_arg);
  VRB.Level(do_arg.verbose_level);
  LRG.Initialize();
 
  { /*Force a lattice instantiation to read stuff*/

    Lattice &lattice = LatticeFactory::Create(F_CLASS_DWF,     G_CLASS_IMPR_RECT);
    //GnoneFnone lat;
    LatticeFactory::Destroy();
  }

  /************************************************
   * config loop
   ************************************************/
  for(int traj=meas_arg.TrajStart; 
      traj< meas_arg.TrajLessThanLimit; 
      traj +=meas_arg.TrajIncrement ) {

    meas_arg.TrajCur = traj;

    ReadGauge(&meas_arg);
    ReadRNG  (&meas_arg); 
   
    Lattice &lat = LatticeFactory::Create(meas_arg.Fermion, meas_arg.Gluon);


    // --------------------------------------------------------
    // Save the rng state to be used for reproducibility test
    // -------------------------------------------------------
    LRGState rng_state;
    rng_state.GetStates();    
    
    // ---
    // do the measurements
    // ---
    run(lat,meas_arg);
    
    //--------------------------------
    // Reproducibility test
    // -------------------------------
    if ( reproduce ) {
      int repro_period = (int) (100.0/repro_percent) * meas_arg.TrajIncrement;
      if ( traj % repro_period == repro_offset ) {
	
        printf("Reproducing Traj.%d\n",traj);
	
	// --------------------------------
	// Restore the random number seeds
	// --------------------------------
	rng_state.SetStates();
	
	for(int task = 0; task < meas_arg.TaskList.TaskList_len; task++){
	  MeasTask meas_task = meas_arg.TaskList.TaskList_val[task];
	  char tmp[128];
	  sprintf(tmp,"%s_repro",meas_task.OutputFilestem);	
	  strcpy(meas_task.OutputFilestem,tmp);
	}
	run(lat,meas_arg);
	//---------------------------------------
	// change meas_arg to its original form
	//---------------------------------------
	if ( !meas_arg.Decode(argv[2],"meas_arg")) {printf("Bum meas_arg\n"); exit(-1);}
      }

    }
  
    LatticeFactory::Destroy();
    
    WriteGauge(&meas_arg);
    
#ifdef USE_SCU_CHECKSUMS
    if ( ! ScuChecksum::CsumSwap() ) { 
      fprintf(stderr, "Checksum mismatch\n");
      exit(-1);
    }
#endif
    //checkpoint current measurement
    checkpoint(&meas_arg);
    meas_arg.TrajCur += meas_arg.TrajIncrement;

    //checkpoint for next measurement
    checkpoint(&meas_arg);
    WriteRNG  (&meas_arg);
    
  } /*End config loop*/
 
  End();
 
  return(0);
}

void CheckDirectory(char *dir){

  char test[50];
  sprintf(test,"%s/TEST",dir);
  ofstream ftest;
  if ( UniqueID() == 0 ) {
    ftest.open(test,ios::trunc);
    if(!ftest.good()){
      cout<<"Data directory "<<dir<<"does not exist. Trying to create one..."<<endl;
      int stat = mkdir(dir,0755);
      if(stat!=0) {
	cout<<"Could not create directory : "<<dir<<endl;
	exit(-1);
      }
    }
    
  }
  

}

void ReadGauge(MeasArg *meas)
{
  char *cname = "cps";
  char *fname = "ReadGauge";
  if ( meas->GaugeIO == MeasIOLoad ) { 
    GnoneFnone lat;
    char lat_file[256];
    sprintf(lat_file,"%s.%d",meas->GaugeStem,meas->TrajCur);
    QioArg rd_arg(lat_file,0.00001);
    rd_arg.ConcurIONumber=meas->IOconcurrency;

    ReadLatticeParallel rl;
    rl.read(lat,rd_arg);
    if(!rl.good()) 
      ERR.General(cname,fname,"Failed read lattice %s",lat_file);
    
  }
}
void WriteGauge(MeasArg *meas)
{ 
  char *cname = "cps";
  char *fname = "WriteGauge";
  if ( meas->GaugeIO == MeasIOSave ) { 
    GnoneFnone lat;
    char lat_file[256];
    sprintf(lat_file,"%s.%d",meas->GaugeStem,meas->TrajCur);
    QioArg wt_arg(lat_file,0.00001);
    wt_arg.ConcurIONumber=meas->IOconcurrency;

    WriteLatticeParallel wl;
    wl.setHeader("MeasurementCodeNotForGaugeProd","Dummy",meas->TrajCur);
    wl.write(lat,wt_arg);
    if(!wl.good()) 
      ERR.General(cname,fname,"Failed write lattice %s",lat_file);
  }
}
void ReadRNG(MeasArg *meas)
{
  char *cname = "cps";
  char *fname = "ReadRNG";
  if ( meas->RNGIO == MeasIOLoad ) { 
    char rng_file[256];
    sprintf(rng_file,"%s.%d",meas->RNGStem,meas->TrajCur);
    if ( !LRG.Read(rng_file,meas->IOconcurrency) ) 
      ERR.General(cname,fname,"Failed RNG file %s",rng_file);
  }
}
void WriteRNG(MeasArg *meas)
{
  char *cname = "cps";
  char *fname = "WriteRNG";
  if ( meas->RNGIO == MeasIOSave ) { 
    char rng_file[256];
    sprintf(rng_file,"%s.%d",meas->RNGStem,meas->TrajCur);
    if ( !LRG.Write(rng_file,meas->IOconcurrency) ) 
      ERR.General(cname,fname,"Failed RNG file %s",rng_file);
  }
}

void checkpoint(MeasArg *meas_arg)
{
  char *cname = "cps";
  char *fname = "checkpoint";
  char vml_file[256];
  int traj = meas_arg->TrajCur;
  char *work =meas_arg->WorkDirectory;
  CheckDirectory(work);

  meas_arg->TrajStart = traj;
  sprintf(vml_file,"%s/meas_arg.%d",work,traj); 
  if ( !meas_arg->Encode(vml_file,"meas_arg")){
    printf("bad meas_arg encode\n");
    exit(-1);
  }

  char start_conf[256];
  sprintf(start_conf,"%s.%d",meas_arg->GaugeStem,traj);
  strcpy(do_arg.start_conf_filename,start_conf);
  sprintf(start_conf,"%s.%d",meas_arg->RNGStem,traj);
  strcpy(do_arg.start_seed_filename,start_conf);
  sprintf(vml_file,"%s/do_arg.%d",work,traj); 
  if ( !do_arg.Encode(vml_file,"do_arg")){
    printf("bad do_arg encode\n");
    exit(-1);
  }
  sprintf(vml_file,"%s/w_spect_arg.%d",work,traj); 
  if ( !w_spect_arg.Encode(vml_file,"w_spect_arg")){
    printf("bad w_spect_arg encode\n");
    exit(-1);
  }
  sprintf(vml_file,"%s/w_spect_cg_arg.%d",work,traj); 
  if ( !w_spect_cg_arg.Encode(vml_file,"w_spect_cg_arg")){
    printf("bad w_spect_cg_arg encode\n");
    exit(-1);
  }
    
}


void run(Lattice &lat, MeasArg &meas_arg) {
  CommonArg common_arg;
  CommonArg common_arg_fg;
  FixGaugeArg w_spect_fg_arg;	  
 
  WspectOutput w_spect_output;
  w_spect_output.fold = BARYON_RAW;
  //flag to determine whether the lattice is gauge fixed      
  int gauge_fixed = 0; 
  
  for(int task = 0; task < meas_arg.TaskList.TaskList_len; task++){

      MeasTask meas_task = meas_arg.TaskList.TaskList_val[task];
      printf("task = %d\n",task); 

   
      FILE *truncate_it;
      /*************************************************
       * AlgPlaq
       *************************************************/
      if(meas_task.Measurement == MeasAlgPlaq){
        printf("Measuring plaq\n");
	NoArg noarg;
	CommonArg common_arg_plaq;
	CheckDirectory(meas_task.OutputFilestem);
	char plaq_file[256];
	sprintf(plaq_file,"%s/plaq.dat.%d",meas_task.OutputFilestem,meas_arg.TrajCur);
	truncate_it = Fopen(plaq_file,"w");
	Fclose(truncate_it);
	common_arg_plaq.set_filename(plaq_file);
	AlgPlaq plaq(lat,&common_arg_plaq,&noarg);
	plaq.run();
      } // end AlgPlaq

      /*************************************************
       * AlgPbp
       *************************************************/     
      if(meas_task.Measurement == MeasAlgPbp){
        printf("Measuring pbp\n");
	PbpArg pbp_arg;
	CommonArg common_arg_pbp;
	if ( !pbp_arg.Decode(meas_task.ArgFilename,"pbp_arg"))
	  {
	    printf("Bum pbp_arg\n"); exit(-1);
	  }

	CheckDirectory(meas_task.OutputFilestem);
	char pbp_file[256];
	sprintf(pbp_file,"%s/pbp.dat.%d",meas_task.OutputFilestem,meas_arg.TrajCur);
	truncate_it = Fopen(pbp_file,"w");
	Fclose(truncate_it);
	
	common_arg_pbp.set_filename(pbp_file);
	AlgPbp pbp(lat,&common_arg_pbp,&pbp_arg);
	pbp.run();
      } // end AlgPbp
      
      /*************************************************
       * AlgFixGauge
       *************************************************/
      if(meas_task.Measurement == MeasAlgFixGauge){

	printf("Gauge Fixing the Lattice\n");
	if ( !w_spect_fg_arg.Decode(meas_task.ArgFilename,"w_spect_fg_arg"))
	  {
	    printf("Bum w_spect_fg_arg\n"); exit(-1);
	  }
	AlgFixGauge fg(lat,&common_arg_fg,&w_spect_fg_arg);
	fg.run();	
	gauge_fixed = 1;
      }

      /*************************************************
       * AlgWspect
       *************************************************/
      
      if (meas_task.Measurement == MeasAlgWspect){
	char filenames[26][MAX_FILENAME];
	char suffix[MAX_FILENAME];
	char w_dir[MAX_FILENAME]; //directory to save AlgWspect data files

	strcpy(w_dir, meas_task.OutputFilestem);

	/********************************************/
	/*Test if the directory already exists*/
	/********************************************/
	CheckDirectory(w_dir);
	if ( !w_spect_arg.Decode(meas_task.ArgFilename,"w_spect_arg")){printf("Bum w_spect_arg\n"); exit(-1);}
	if ( !w_spect_cg_arg.Decode("w_spect_cg_arg.vml","w_spect_cg_arg")){printf("Bum w_spect_cg_arg\n"); exit(-1);}
	
	printf("Measuring hadron masses traj. %d\n",meas_arg.TrajCur);
	int traj = meas_arg.TrajCur;
	sprintf(filenames[0],"%s/scalar.dat.%d" ,w_dir,traj);
	sprintf(filenames[1],"%s/vector_x.dat.%d" ,w_dir,traj);
	sprintf(filenames[2],"%s/vector_y.dat.%d" ,w_dir,traj);
	sprintf(filenames[3],"%s/tensor_xy.dat.%d",w_dir,traj);
	sprintf(filenames[4], "%s/vector_z.dat.%d",w_dir,traj);
	sprintf(filenames[5], "%s/tensor_xz.dat.%d",w_dir,traj);
	sprintf(filenames[6], "%s/tensor_yz.dat.%d",w_dir,traj);
	sprintf(filenames[7], "%s/axial_vector_t.dat.%d",w_dir,traj);
	sprintf(filenames[8], "%s/vector_t.dat.%d",w_dir,traj);
	sprintf(filenames[9], "%s/tensor_xt.dat.%d",w_dir,traj);
	sprintf(filenames[10], "%s/tensor_yt.dat.%d",w_dir,traj);
	sprintf(filenames[11], "%s/axial_vector_z.dat.%d" ,w_dir,traj);
	sprintf(filenames[12], "%s/tensor_zt.dat.%d" ,w_dir,traj);
	sprintf(filenames[13], "%s/axial_vector_y.dat.%d" ,w_dir,traj);
	sprintf(filenames[14], "%s/axial_vector_x.dat.%d" ,w_dir,traj);
	sprintf(filenames[15], "%s/pseudoscalar.dat.%d",w_dir,traj);
	
	for(int ind=0;ind<16;ind++) filenames[ind][MAX_FILENAME-1]='\0';
	
	w_spect_output.meson_name00            = &(filenames[0][0]) ;
	w_spect_output.meson_name01            = &(filenames[1][0]) ;
	w_spect_output.meson_name02            = &(filenames[2][0]) ;
	w_spect_output.meson_name03            = &(filenames[3][0]) ;
	w_spect_output.meson_name04            = &(filenames[4][0]) ;
	w_spect_output.meson_name05            = &(filenames[5][0]) ;
	w_spect_output.meson_name06            = &(filenames[6][0]) ;
	w_spect_output.meson_name07            = &(filenames[7][0]) ;
	w_spect_output.meson_name08            = &(filenames[8][0]) ;
	w_spect_output.meson_name09            = &(filenames[9][0]) ;
	w_spect_output.meson_name10            = &(filenames[10][0]) ;
	w_spect_output.meson_name11            = &(filenames[11][0]) ;
	w_spect_output.meson_name12            = &(filenames[12][0]) ;
	w_spect_output.meson_name13            = &(filenames[13][0]) ;
	w_spect_output.meson_name14            = &(filenames[14][0]) ;
	w_spect_output.meson_name15            = &(filenames[15][0]) ;
	
	sprintf(filenames[16],"%s/cg.dat.%d" ,w_dir,traj);
	sprintf(filenames[17],"%s/mid_point.dat.%d",w_dir,traj);
	sprintf(filenames[18], "%s/a0_p.dat.%d",w_dir,traj);
	sprintf(filenames[19], "%s/nucleon.dat.%d",w_dir,traj);
	sprintf(filenames[20], "%s/nucleon_prime.dat.%d",w_dir,traj);
	sprintf(filenames[21], "%s/delta_x.dat.%d",w_dir,traj);
	sprintf(filenames[22], "%s/delta_y.dat.%d",w_dir,traj);
	sprintf(filenames[23], "%s/delta_z.dat.%d",w_dir,traj);
	sprintf(filenames[24], "%s/delta_t.dat.%d",w_dir,traj);
	sprintf(filenames[25], "%s/pbp.dat.%d",w_dir,traj);
	
	for(int ind=16;ind<26;ind++) filenames[ind][MAX_FILENAME-1]='\0';
	for(int ind=0;ind<26;ind++){
	  truncate_it = Fopen(filenames[ind],"w");
	  Fclose(truncate_it);
	}
	
	//-------------------------------------------------
	//cg,nucleon and other stuff
	//-------------------------------------------------
	w_spect_output.cg            = &(filenames[16][0]) ;
	w_spect_output.mid_point     = &(filenames[17][0]) ;
	w_spect_output.a0_p          = &(filenames[18][0]) ;
	w_spect_output.nucleon       = &(filenames[19][0]) ;
	w_spect_output.nucleon_prime = &(filenames[20][0]) ;
	w_spect_output.delta_x       = &(filenames[21][0]) ;
	w_spect_output.delta_y       = &(filenames[22][0]) ;
	w_spect_output.delta_z       = &(filenames[23][0]) ;
	w_spect_output.delta_t       = &(filenames[24][0]) ;
	w_spect_output.pbp           = &(filenames[25][0]) ;
	
	common_arg.results = (void *) &w_spect_output ;

	
	printf("start AlgWspect\n");
	AlgWspect ws(lat,&common_arg,&w_spect_arg,&w_spect_cg_arg);
	
	ws.SetCounter(traj, 0) ;
	
	ws.run();
	/*
	cout<<"====================================="<<endl;
	cout<<"Synchronizing machine..."<<endl;
	Gint::SynchMachine();
	cout<<"Machine Synched!"<<endl;
	cout<<"====================================="<<endl;
	*/
      }

      /*************************************************
       * AlgNuc3pt
       *************************************************/  
#ifdef NUC3PT     
      if(meas_task.Measurement == MeasAlgNuc3pt){
	printf("Measuring Nuc3pt\n");
	Nuc3ptArg nuc3pt_arg;
	CommonArg common_arg_nuc;
	if ( !nuc3pt_arg.Decode(meas_task.ArgFilename,"nuc3pt_arg"))
	  {
	    printf("Bum nuc3pt_arg\n"); exit(-1);
	  }
	
	CheckDirectory(meas_task.OutputFilestem);
	char nuc_file[256];
	sprintf(nuc_file,"%s/nuc3pt.dat.%d",meas_task.OutputFilestem,meas_arg.TrajCur);
	truncate_it = Fopen(nuc_file,"w");
	Fclose(truncate_it);
	
	common_arg_nuc.set_filename(nuc_file);
	AlgNuc3pt nuc3pt(lat,&common_arg_nuc,&nuc3pt_arg);
	nuc3pt.run();
      } // end AlgNuc3pt
#endif
    } //end task loop  
    
    if( gauge_fixed ){
      AlgFixGauge fg(lat,&common_arg_fg,&w_spect_fg_arg);
      fg.free();
    }



}
