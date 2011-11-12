#undef RELOAD
#define VOLFMT QIO_SINGLEFILE
#include<config.h>
#include<precision.h>

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
#include<alg/qpropw.h>
#include<alg/meson.h>
#include<alg/nuc2pt.h>
#include<alg/eigcg_arg.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>

#include <util/qio_readLattice.h>
#include <util/qio_writeLattice.h>

#include<util/qioarg.h>

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
void CalMesons(QPropW &q, char *out_dir, int traj, char *src_str="POINT");
void CalBaryons(QPropW &q, char *out_dir, int traj, char *src_str="POINT");

DoArg do_arg;
MeasArg meas_arg;
MeasArg meas_arg_save;
WspectArg w_spect_arg;
CgArg w_spect_cg_arg;

CPS_START_NAMESPACE
class EigCgDriver{
	static char *cname;
	int nev;
	int m;
	Float max_eig;
	int max_def_len;
	int restart_len;
	Float *restart;
	bool always_restart;

	int def_len;
	int vec_len;
	Vector **V;
	Float *M;
	float **U;
	Rcomplex *H;
	int do_rerun;
	double precision;
	public:
		EigCgDriver(){
			V=NULL;
			U=NULL;
			H=NULL;
		};
		~EigCgDriver(){};
		void init(Lattice &lat, EigCGArg &eigcg_arg){
			char *fname = "init()";
			nev=eigcg_arg.nev;
			m=eigcg_arg.m;
			max_eig=eigcg_arg.max_eig_cut;
			max_def_len=eigcg_arg.max_def_len;
			restart_len=eigcg_arg.restart_len;
			restart=eigcg_arg.restart;
			always_restart=eigcg_arg.always_restart;
			do_rerun=0;
			precision=1e-8;

			def_len=0;
			vec_len=GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
			V=new Vector*[m];
			for(int i=0;i<m;i++)
				V[i]=(Vector *)smalloc(cname,fname,"V",vec_len*sizeof(Float));
			M=new Float[m];
			U=new float*[max_def_len];
			for(int i=0;i<max_def_len;i++)
				U[i]=(float *)smalloc(cname,fname,"U",vec_len*sizeof(float));
			H=new Rcomplex[max_def_len*max_def_len];
		}
		int run_qpropw_acc( QPropW &qpropw ){
			qpropw.eig_Run(V,vec_len,M,max_eig,nev,m,U,H,max_def_len,def_len,restart,restart_len,always_restart,do_rerun,precision);
			return 1;
		}
		int run_qpropw( QPropW &qpropw ){
			qpropw.eig_Run(NULL,vec_len,M,max_eig,0,0,U,H,max_def_len,def_len,restart,restart_len,always_restart,do_rerun,precision);
			return 1;
		}
};

char *EigCgDriver::cname = "EigCgDriver";

CPS_END_NAMESPACE

int main(int argc, char *argv[])
{ 
//  char *cname=argv[0];
//  char *fname="main()";
  Start(&argc, &argv);

  if ( argc!=7) { 
    printf("Args: doarg-file measarg-file  initial-directory enable_reproduce repro_percent repro_offset\n");
    exit(-1);
  }
 
  chdir (argv[3]);

  if ( !do_arg.Decode(argv[1],"do_arg") ) 
    { 
      printf("Failedto decode do_arg\n"); exit(-1);
    }
  if ( !meas_arg.Decode(argv[2],"meas_arg"))
    {
      printf("Failed to decode meas_arg\n"); exit(-1);
    }
 
  
  int reproduce = atoi(argv[4]);
  double repro_percent = atof(argv[5]);
  int repro_offset = atoi(argv[6]);

  /*Layout the lattice on the machine (without regard to even-odd)*/
  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = SizeS();

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
  GJP.setArg(&argc, &argv);
  VRB.Level(do_arg.verbose_level);
  LRG.Initialize();
 
  { /*Force a lattice instantiation to read stuff*/
    GnoneFnone lat;
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
{
	//fix gauge
	VRB.Result("","main()","Fixing  Gauge\n");
  	CommonArg common_arg;
	FixGaugeArg fix_gauge_arg;
	fix_gauge_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
	fix_gauge_arg.hyperplane_start = 0;
	fix_gauge_arg.hyperplane_step = 1;
	fix_gauge_arg.hyperplane_num = GJP.TnodeSites()*GJP.Tnodes();
	fix_gauge_arg.stop_cond = 1e-8;
	fix_gauge_arg.max_iter_num = 10000;
	AlgFixGauge algfixgauge(lat,&common_arg,&fix_gauge_arg);
	algfixgauge.run();
	VRB.Result("","main()","Fixing  Gauge End\n");

	run(lat,meas_arg);
}
    
    //--------------------------------
    // Reproducibility test
    // -------------------------------
    if ( reproduce ) {
      int repro_period = (int) (100.0/repro_percent) * meas_arg.TrajIncrement;
      if ( traj % repro_period == repro_offset ) {
	
        VRB.Result("","main()","Reproducing Traj.%d\n",traj);
	
	// --------------------------------
	// Restore the random number seeds
	// --------------------------------
	rng_state.SetStates();
	
	for(unsigned int task = 0; task < meas_arg.TaskList.TaskList_len; task++){
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

  char test[256];
  sprintf(test,"%s/TEST",dir);
  ofstream ftest;
  if ( UniqueID() == 0 ) 
  {
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
//  char *cname = "cps";
//  char *fname = "checkpoint";
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
    
}


void run(Lattice &lat, MeasArg &meas_arg) {
  CommonArg common_arg;
  CommonArg common_arg_fg;
  FixGaugeArg w_spect_fg_arg;	  
 
  WspectOutput w_spect_output;
  w_spect_output.fold = BARYON_RAW;
  //flag to determine whether the lattice is gauge fixed      
  int gauge_fixed = 0; 
  
  for(unsigned int task = 0; task < meas_arg.TaskList.TaskList_len; task++){

      MeasTask meas_task = meas_arg.TaskList.TaskList_val[task];
      VRB.Result("","main()","task = %d\n",task); 

   
      FILE *truncate_it;
      /*************************************************
       * AlgPlaq
       *************************************************/
      if(meas_task.Measurement == MeasAlgPlaq){
        VRB.Result("","main()","Measuring plaq\n");
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
        VRB.Result("","main()","Measuring pbp\n");
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

	VRB.Result("","main()","Gauge Fixing the Lattice\n");
	if ( !w_spect_fg_arg.Decode(meas_task.ArgFilename,"w_spect_fg_arg"))
	  {
	    printf("Bum w_spect_fg_arg\n"); exit(-1);
	  }
	AlgFixGauge fg(lat,&common_arg_fg,&w_spect_fg_arg);
	fg.run();	
	gauge_fixed = 1;
      }


      /*************************************************
       * AlgNuc3pt
       *************************************************/  
#ifdef NUC3PT     
      if(meas_task.Measurement == MeasAlgNuc3pt){
	VRB.Result("","main()","Measuring Nuc3pt\n");
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
      //----------------------------------------------------//
      // Run QPropW to generate quark propagators           //
      //----------------------------------------------------//
      if(meas_task.Measurement == MeasAlgQPropW){
	VRB.Result("","main()","Measuring QPropW\n");
	QPropWArg qpropw_arg;
	QPropWBoxArg box_arg;
	EigCGArg eigcg_arg;
	CommonArg c_arg_qprop;
	
	if ( !qpropw_arg.Decode(meas_task.ArgFilename,"qpropw_arg"))
	  {
	    printf("Bum qpropw_arg\n"); exit(-1);
	  }
	if ( !eigcg_arg.Decode("eigcg_arg.vml","eigcg_arg"))
	  {
	    printf("Bum eigcg_arg\n"); exit(-1);
	  }
  
	CheckDirectory(meas_task.OutputFilestem);
	char qprop_file[256];
	char prop_dir[256];
	strcpy(prop_dir, qpropw_arg.file); 
	CheckDirectory(prop_dir);
	
	qpropw_arg.seqNum = meas_arg.TrajCur;
	int box_size = 16;
	int box_start = GJP.Sites(0)/2;
//	if (box_size > box_start) box_start = GJP.Sites(0)-box_size;
	int t_offset = GJP.Sites(3)/8;
 	char src_center[256];
	EigCgDriver eig_cg;
	eig_cg.init(lat,eigcg_arg);
	for (int i = 0; i<8; i++){
		if (i%2==0){ box_arg.box_start=0; }
		else { box_arg.box_start=box_start; }
		box_arg.box_end = box_arg.box_start + box_size;
		int b_i = box_arg.box_start;
		qpropw_arg.t = i*t_offset;
		sprintf(src_center,"x%dy%dz%dt%d",b_i, b_i, b_i, qpropw_arg.t);	
		sprintf(qprop_file,"%s/qpropw_%s.dat.%d", meas_task.OutputFilestem, src_center, meas_arg.TrajCur);
		
		c_arg_qprop.set_filename(qprop_file);
		FILE *truncate_it = Fopen(qprop_file,"w");
		Fclose(truncate_it);
		
		sprintf(qprop_file,"%sprop_%s.%d", prop_dir, src_center, meas_arg.TrajCur); 
		qpropw_arg.file = qprop_file;
		
		QPropWBoxSrc prop(lat, &qpropw_arg, &box_arg, &c_arg_qprop);
	  if ( i < 2 ) eig_cg.run_qpropw_acc(prop);
		else  eig_cg.run_qpropw(prop);
		
		CalMesons(prop, meas_task.OutputFilestem, meas_arg.TrajCur, &src_center[0]);
		CalBaryons(prop, meas_task.OutputFilestem, meas_arg.TrajCur, &src_center[0]);
	}

#ifdef RELOAD
        prop.ReLoad(qpropw_arg.file);
	CalMesons(prop, meas_task.OutputFilestem, meas_arg.TrajCur, &src_center[0]);
        CalBaryons(prop, meas_task.OutputFilestem, meas_arg.TrajCur, &src_center[0]);
#endif
	
      }
  } //end task loop  
    
  if( gauge_fixed ){
    AlgFixGauge fg(lat,&common_arg_fg,&w_spect_fg_arg);
    fg.free();
  }

   



}

void CalMesons(QPropW &q, char *out_dir, int traj, char *src_str){
 
  char file[256];
  sprintf(file, "%s/meson_%s.dat.%d", out_dir,src_str,traj);
  FILE *fp = Fopen(&file[0],"w");
 
  int mu, nu;
  //vector and scalar
  for ( mu = 0; mu < 4; mu++ ){
    Meson mes(mu, src_str);
    mes.Zero();
    mes.setMass(q.Mass(),q.Mass());
    mes.calcMeson(q,q);
    if(!UniqueID())
      mes.Print(fp);
  }

  //pseudoscalar
  {
    mu = -5;
    Meson mes(mu, src_str);
    mes.Zero();
    mes.setMass(q.Mass(),q.Mass());
    mes.calcMeson(q,q);
    if(!UniqueID())
    mes.Print(fp);
  }
  //tensor
  for ( mu = 0; mu < 4; mu++ ){
    for ( nu = mu+1; nu < 4; nu++ ){
      Meson mes1(mu, nu, src_str);
      mes1.Zero();
      mes1.setMass(q.Mass(),q.Mass());
      mes1.calcMeson(q,q);
      if(!UniqueID())
      mes1.Print(fp);

      Meson mes2(nu, mu, src_str);
      mes2.Zero();
      mes2.setMass(q.Mass(),q.Mass());
      mes2.calcMeson(q,q);
    if(!UniqueID())
      mes2.Print(fp);
    }
  }

  //axialvector
  nu = -5;
  for ( mu = 0; mu < 4; mu++ ){
    Meson mes(mu, nu, src_str);
    mes.Zero();
    mes.setMass(q.Mass(),q.Mass());
    mes.calcMeson(q,q);
    if(!UniqueID())
    mes.Print(fp);
  }
  Fclose(fp);

}

void CalBaryons(QPropW &q, char *out_dir, int traj, char *src_str){
  Nuc2pt nuc(NUC_G5C, POINT);
  nuc.Zero();
  nuc.calcNucleon(q);
  char file[256];
  sprintf(file, "%s/nucleon_%s.dat.%d", out_dir,src_str,traj);
  FILE *fp = Fopen(&file[0],"w");
  if(!UniqueID())
  nuc.Print(fp);

  Nuc2pt nuc2(NUC_C, POINT);
  nuc2.Zero();
  nuc2.calcNucleon(q);
    if(!UniqueID())
  nuc2.Print(fp);
  Fclose(fp);
}
