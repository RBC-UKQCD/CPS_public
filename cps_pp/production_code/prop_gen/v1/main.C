# define VOLFMT QIO_PARTFILE

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

#include <util/qcdio.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <comms/scu.h>
#include <comms/glb.h>

#include <util/lattice.h>
#include <util/time_cps.h>
#include <alg/do_arg.h>
#include <alg/array_arg.h>
#include <alg/no_arg.h>
#include <alg/common_arg.h>
#include <alg/hmd_arg.h>
#include <alg/alg_plaq.h>
#include <alg/alg_rnd_gauge.h>

#if(0==1)
 #include <ReadLattice.h>
 #include <WriteLattice.h>
#endif

#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

#include <util/command_line.h>
#include <sstream>

#include<unistd.h>
#include<config.h>

#include <alg/qpropw.h>
#include <alg/qpropw_arg.h>

#include <alg/alg_fix_gauge.h>
#include <alg/fix_gauge_arg.h>

#include <util/data_shift.h>

using namespace std;
USING_NAMESPACE_CPS

#if 0
const int num_masses=6;
const Float masses[num_masses]={0.002,0.004,0.006,0.008,0.025,0.03};
#else
//const int num_masses=2;
//const Float masses[num_masses]={0.025,0.03};
FloatArray mass_list;
static int num_masses;
static Float *masses;
#endif

void chkpt(const int num_nodes,int& chkpoint_no,Float dtime[],const int dtime_size);

inline Matrix operator * (const Matrix& m1, const Matrix& m2)
{ Matrix r; r.DotMEqual(m1,m2); return r; }

void rotate_gauge_explicit(Lattice &lat,int dir=3);

int IND(int x, int y, int z, int t, int l);

int main(int argc,char *argv[])
{

  char *cname="";
  char *fname="main()";
  Start(&argc,&argv);

  VRB.Result(cname,fname,"Starting next configuration.\n");

  //Initialize Timing
  const int dtime_size=200;
  Float dtime[dtime_size];
  dtime[0]=dclock();

  int chkpoint_no=0;

  CommandLine::is(argc,argv);

  DoArg do_arg;
  EvoArg evo_arg;
  QPropWArg qpropw_arg;


  /* get parameter from command line */

  /* 17 parameters: */
  /* x.x  
	  chkpoints (1 to turn checkpoints on, 0 to turn them off)
          DIRECTORY do_arg evo_arg mass_list
          lat_stand_in lat_stand_out
	  number (which random gauge to take)
	  t_src
	  label id seqNum
	  logdir
	  plaquette_gfix_info_dir
	  midprop_contractions_dir
	  configs_done_log
	  gauge_rotated_lats_dir
	  gauge_rotated_lats_saved_log
	  
 */

  int chkpoints(CommandLine::arg_as_int() );
 
  chdir(CommandLine::arg());
  
  
  if ( !do_arg.Decode(CommandLine::arg(),"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}
  if ( !evo_arg.Decode(CommandLine::arg(),"evo_arg") ) { printf("Bum evo_arg\n"); exit(-1);}
  if ( !mass_list.Decode(CommandLine::arg(),"mass_list") ) { printf("Bum mass_lilst\n"); exit(-1);}
  num_masses = mass_list.Floats.Floats_len;
  masses = mass_list.Floats.Floats_val;

  do_arg.gfix_chkb=1;
  do_arg.cg_reprod_freq=10;

  GJP.Initialize(do_arg);


//  const int num_nodes=do_arg.x_nodes*do_arg.y_nodes*do_arg.z_nodes*do_arg.t_nodes*do_arg.s_nodes;
  const int num_nodes=GJP.Xnodes()*GJP.Ynodes()*GJP.Znodes()*GJP.Tnodes()*GJP.Snodes();

  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime,dtime_size);

  char lat_stand_in[200], lat_stand_out[200];
  char label[200], id[200];
  char logdir[200], plaquette_gfix_info_dir[200], midprop_contractions_dir[200], configs_done_log[200];
  char gauge_rotated_lats_dir[200], gauge_rotated_lats_saved_log[200];

  sprintf(lat_stand_in, CommandLine::arg() );
  sprintf(lat_stand_out, CommandLine::arg() );
  
  int gauge_ran(CommandLine::arg_as_int() );
  int t_src(CommandLine::arg_as_int() );
  if ( t_src<0 || t_src>=GJP.Sites(3) ) {
    printf("t_src must be in the range 0 <= t_src < %d.\n",GJP.Sites(3));
    printf("Abort.\n");
    return 1;
  }
  if ( t_src % GJP.TnodeSites() ) {
    printf("t_src must be at the beginning of a node since we are shifting\n");
    printf("the lattice rather than putting a source at time t_src.\n");
    printf("Abort.\n");
    return 2;
  }

  sprintf(label, CommandLine::arg() );
  sprintf(id, CommandLine::arg() );
  int seqNum(CommandLine::arg_as_int() );

  sprintf(logdir, CommandLine::arg() );
  sprintf(plaquette_gfix_info_dir, CommandLine::arg() );
  sprintf(midprop_contractions_dir, CommandLine::arg() );
  sprintf(configs_done_log, CommandLine::arg() );
  sprintf(gauge_rotated_lats_dir, CommandLine::arg() );
  sprintf(gauge_rotated_lats_saved_log, CommandLine::arg() );
  
  VRB.Result(cname,fname,"This is configuration number %d.\n",seqNum);

  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime,dtime_size);

  const FP_FORMAT outformat ( 	FP_IEEE64BIG);

			     /* one of
				FP_UNKNOWN
				FP_AUTOMATIC
				FP_TIDSP32
				FP_IEEE32
				FP_IEEE32BIG
				FP_IEEE32LITTLE
				FP_IEEE64
				FP_IEEE64BIG
				FP_IEEE64LITTLE
			     */

  
  GwilsonFdwf lattice;


  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime,dtime_size);

  //If lat_stand_in is "none" then generate a random lattice, otherwise load
  //the file lat_stand_in.
  
  if (strcmp(lat_stand_in,"none")==0){
    VRB.Result(cname,fname,"ncreate random gauge - field (taking # %i) \n", gauge_ran);
    
    for(int ii(0); ii < gauge_ran; ++ii)
      lattice.SetGfieldDisOrd();
    
    VRB.Result(cname,fname,"Created random gauge field.\n");
  } 
  else {
    ReadLatticeParallel readLat;
    
    VRB.Result(cname,fname,"  reading: %s (NERSC-format)\n",lat_stand_in);

    readLat.read(lattice,lat_stand_in);

    VRB.Result(cname,fname,"Lattice read.\n");
    
  }

  //If lat_stand_out is specified then save the non-gauge rotated lattice
  //to this file.

  if (strcmp(lat_stand_out,"none")!=0) {
      WriteLatticeParallel writeLat;

      VRB.Result(cname,fname,"  writing: %s (NERSC-format)\n",lat_stand_out);

      writeLat.write(lattice, lat_stand_out, outformat);
  
      VRB.Result(cname,fname,"Lattice written.\n");
  }

  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime,dtime_size);

  //Shift lattice so that time slice t_src moves to t=0
  int t_shift=-t_src/GJP.TnodeSites();
  VRB.Result(cname,fname,"Shifting lattice by %d, command GDS.Set(0,0,0,%d).\n",-t_src,t_shift);
  GDS.Set(0,0,0,t_shift);
  lattice.Shift();
  VRB.Result(cname,fname,"Lattice shifted.\n");
  
  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime,dtime_size);

  FixGaugeArg fix_arg;
  fix_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
  fix_arg.hyperplane_start = 0;
  fix_arg.hyperplane_step = 1;
  fix_arg.hyperplane_num = GJP.Tnodes()*GJP.TnodeSites();
  fix_arg.stop_cond = 1e-8;
  fix_arg.max_iter_num = 10000;

  FILE *fp;

  CommonArg common_arg_gfix;
  char gfix_fname[200];
  sprintf(gfix_fname,"%s/gfix.%d.dat",plaquette_gfix_info_dir,seqNum);
  common_arg_gfix.results=gfix_fname;
  if( (fp = Fopen((char *)(gfix_fname),"w")) == NULL ){
    ERR.FileA("main","main", (char *)(gfix_fname) );
  }
  Fclose(fp);

  CommonArg common_arg_plaq;
  NoArg plaq_arg;
  char plaq_fname[200];
  sprintf(plaq_fname,"%s/plaq.%d.dat",plaquette_gfix_info_dir,seqNum);
  common_arg_plaq.results=plaq_fname;
  if( (fp = Fopen((char *)(plaq_fname),"w")) == NULL ){
    ERR.FileA("main","main", (char *)(plaq_fname) );
  }
  Fclose(fp);

  VRB.Result(cname,fname,"Calculating plaquette.\n");
  AlgPlaq plaq(lattice,&common_arg_plaq,&plaq_arg);
  plaq.run();
  VRB.Result(cname,fname,"Plaquette calculated.\n");

  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime,dtime_size);

  VRB.Result(cname,fname,"Calculating gauge fixing matrices.\n");
  AlgFixGauge fix_gauge(lattice,&common_arg_gfix,&fix_arg);
  fix_gauge.run();
  VRB.Result(cname,fname,"Gauge fixing matrices calculated.\n");

  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime,dtime_size);

  qpropw_arg.x=0;
  qpropw_arg.y=0;
  qpropw_arg.z=0;
  qpropw_arg.t=0; //Shift lattice instead of putting source at t!=0 (see below)
  qpropw_arg.gauge_fix_src=1;
  qpropw_arg.gauge_fix_snk=1;
  qpropw_arg.store_midprop=0;
  qpropw_arg.save_prop=1;
  qpropw_arg.do_half_fermion=0;	
  
  qpropw_arg.cg.max_num_iter=99999;
  qpropw_arg.cg.stop_rsd=1e-10;
  qpropw_arg.cg.true_rsd=1e-10;
  qpropw_arg.cg.RitzMatOper=MAT_HERM;
  qpropw_arg.cg.Inverter=CG;
  qpropw_arg.cg.bicgstab_n=0;

  //See which masses we need to calculate props for
  float mass_read_tmp;
  Float mass_read;
  int bc;
  int mass_done[num_masses][2];
  char logfile[200];
  sprintf(logfile,"%s/%d",logdir,seqNum);
  if( (fp = fopen((char *)(logfile),"r")) != NULL ){
//    ERR.FileA("main","main", (char *)(logfile) );
    for (int i=0; i<num_masses; i++)
      for (int j=0; j<2; j++)
        mass_done[i][j]=0;
    VRB.Result(cname,fname,"Finished masses and boundary conditions as read from log (config no. %d):\n",seqNum);
    while (fscanf(fp,"%f",&mass_read_tmp)>=0){
      mass_read=(Float) mass_read_tmp;
      fscanf(fp,"%d",&bc);
      VRB.Result(cname,fname,"f %d\n",mass_read,bc);
      for (int i=0; i<num_masses; i++) {
        if (fabs(masses[i]-mass_read)<=1e-6)
  	mass_done[i][bc]=1;
      }
    }
    fclose(fp);
  }
  VRB.Result(cname,fname,"Masses and boundary conditions to do (config. no. %d):\n",seqNum);
  for (int i=0; i<num_masses; i++)
    for (bc=0; bc<2; bc++)
      if (!mass_done[i][bc])
	VRB.Result(cname,fname,"%f %d",masses[i],bc);
  
  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime,dtime_size);

  //Loop over masses and calculate props
  CommonArg common_arg;
  char prop_name[200], midprop_contractions_fname[200], qio_out[200];
  VRB.Result(cname,fname,"Begin mass loop for configuration number %d.\n",seqNum);
  for (int i=0; i<num_masses; i++) {
    if (mass_done[i][0] && mass_done[i][1]){
      VRB.Result(cname,fname,"Nothing to do mass %f\n",masses[i]);
      continue;
    }
    qpropw_arg.cg.mass=masses[i];
    for (bc=0; bc<2; bc++) {
      if (!mass_done[i][bc]) {
	char bc_type[30], bc_label[30];
	if (bc) {
	  sprintf(bc_type,"antiperiodic");
	  sprintf(bc_label,"APRD");
	} else {
	  sprintf(bc_type,"periodic");
	  sprintf(bc_label,"PRD");
	}
	VRB.Result(cname,fname,"Inverting mass %f, %s boundary conditions, configuration number %d.\n",masses[i],bc_type,seqNum);
	//Checkpoint
	if (chkpoints)
	  chkpt(num_nodes,chkpoint_no,dtime,dtime_size);
	sprintf(prop_name,"32cube_b2.25_Iw_msea0.004_traj%d_m%0.3f_tshift%d_%s",seqNum,masses[i],-t_src,bc_label); //Note that t_src is w.r.t. the unshifted lattice in the file name, but the source time in the qio header will always say 0 because it is equal to qpropw_arg.t
	sprintf(midprop_contractions_fname,"%s/midprop_contractions.%s",midprop_contractions_dir,prop_name);
	if( (fp = Fopen((char *)(midprop_contractions_fname),"w")) == NULL ){
	  ERR.FileA("main","main", (char *)(midprop_contractions_fname) );
	}
	Fclose(fp);
#if TARGET==QCDOC
	sprintf(qio_out,"/pfs/qio_prop.%s",prop_name);
#else
	sprintf(qio_out,"%d",seqNum);
        mkdir(qio_out,0777);
	sprintf(qio_out,"%d/qio_prop.%s",seqNum,prop_name);
#endif
	common_arg.results=midprop_contractions_fname;
	qpropw_arg.file=qio_out;
	qpropw_arg.ensemble_id=id;
	qpropw_arg.ensemble_label=label;
	qpropw_arg.seqNum=seqNum;
	if (bc)
	  GJP.Tbc(BND_CND_APRD);
	else
	  GJP.Tbc(BND_CND_PRD);
	QPropWWallSrc propagator(lattice, &qpropw_arg, &common_arg);
	VRB.Result(cname,fname,"Finished inversion, mass %f, %s boundary conditions, configuration number %d.\n",masses[i],bc_type,seqNum);
	//Checkpoint
	if (chkpoints)
	  chkpt(num_nodes,chkpoint_no,dtime,dtime_size);
	if( (fp = Fopen((char *)(logfile),"a")) == NULL ){
	  ERR.FileA("main","main", (char *)(logfile) );
	}
	Fprintf(fp,"%0.3f %d\n",masses[i],bc);
	Fclose(fp);
      }
    }
  }
  VRB.Result(cname,fname,"Finished mass loop for configuration number %d.\n",seqNum);

  //Gauge rotate and save the lattice if this hasn't been done already
  //First check if the gauge rotated lattice has already been saved
  if( (fp = fopen((char *)(gauge_rotated_lats_saved_log),"r")) == NULL ){
    ERR.FileA("main","main", (char *)(gauge_rotated_lats_saved_log) );
  }
  int config_no_read;
  int config_saved=0;
  while (fscanf(fp,"%d",&config_no_read)>=0){
    if (config_no_read==seqNum){
      config_saved=1;
      break;
    }
  }
  fclose(fp);
  //If it hasn't then gauge rotate the lattice and save it
  if (!config_saved){
    VRB.Result(cname,fname,"Gauge rotating the lattice and saving it.\n");
    //Gauge rotate lattice
    VRB.Result(cname,fname,"Gauge rotating lattice.\n");
    rotate_gauge_explicit(lattice);
    VRB.Result(cname,fname,"Gauge rotated lattice.\n");
    //Shift the lattice back
    VRB.Result(cname,fname,"Shifting lattice back.  Shifting by %d, command GDS.Set(0,0,0,%d).\n",t_src,-t_shift);
    GDS.Set(0,0,0,-t_shift);
    lattice.Shift();
    VRB.Result(cname,fname,"Lattice shifted.\n");
    //Save it
    char rotated_lat_fname[200];
    sprintf(rotated_lat_fname,"%s/ckpoint_lat_gauge_rotated.%d",gauge_rotated_lats_dir,seqNum);
    VRB.Result(cname,fname,"Writing gauge rotated lattice.\n");
    WriteLatticeParallel writeLat;
    writeLat.write(lattice, rotated_lat_fname, outformat);
    VRB.Result(cname,fname,"Wrote gauge rotated lattice.\n");
    VRB.Result(cname,fname,"Finished gauge rotating and saving lattice.\n");
    //Checkpoint
    if (chkpoints)
      chkpt(num_nodes,chkpoint_no,dtime,dtime_size);
    //Record in log that this lattice has been saved
    if( (fp = Fopen((char *)(gauge_rotated_lats_saved_log),"a")) == NULL ){
    ERR.FileA("main","main", (char *)(gauge_rotated_lats_saved_log) );
    }
    Fprintf(fp,"%d\n",seqNum);
    Fclose(fp);
  }

  if( (fp = Fopen((char *)(configs_done_log),"a")) == NULL ){
    ERR.FileA("main","main", (char *)(configs_done_log ) );
  }
  Fprintf(fp,"%d\n",seqNum);
  Fclose(fp);


  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime,dtime_size);
  

  fix_gauge.free();

  VRB.Result(cname,fname,"Finished configuration number %d.\n",seqNum);
  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime,dtime_size);


  return 0;

}

//---------------------------------------------------------------------------
//Checksum, machine synchronization and timing checkpoint function
//---------------------------------------------------------------------------
void chkpt(const int num_nodes,int& chkpoint_no,Float dtime[],const int dtime_size)
{
  ERR.HdwCheck("main","chkpt");
  int dummy=0;
  Float test_val;
  Float* ptest_val=&test_val;
  char *cname="";
  char *fname="chkpt()";
  
  for (int ii=0; ii<num_nodes; ii++){
    test_val=0.0;
    if (UniqueID()==ii){
      test_val=1.0;
    }
    glb_sum_five(ptest_val);
    while(test_val==0.0)
      dummy=0; //Does absolutely nothing, here as a placeholder for while
  }
  VRB.Result(cname,fname,"Checkpoint no. %d reached.\n",chkpoint_no++);
  if ( chkpoint_no < dtime_size ) {
    dtime[chkpoint_no]=dclock();

    Float time_tmp=dtime[chkpoint_no]-dtime[chkpoint_no-1];
    int hr_tmp=time_tmp/3600.0;
    int min_tmp=(time_tmp-3600.0*hr_tmp)/60.0;
    Float sec_tmp=time_tmp-3600.0*hr_tmp-60.0*min_tmp;
    VRB.Result(cname,fname,"Time since last checkpoint: %d hours %d minutes %f seconds.\n",hr_tmp,min_tmp,sec_tmp);

    time_tmp=dtime[chkpoint_no]-dtime[0];
    hr_tmp=time_tmp/3600.0;
    min_tmp=(time_tmp-3600.0*hr_tmp)/60.0;
    sec_tmp=time_tmp-3600.0*hr_tmp-60.0*min_tmp;
    VRB.Result(cname,fname,"Time since beginning: %d hours %d minutes %f seconds.\n",hr_tmp, min_tmp,sec_tmp);
  }
}
//---------------------------------------------------------------------------


//-----------------------------------------------------------------------
// Rotate gauge of the lattice (from Min)
//-----------------------------------------------------------------------
void rotate_gauge_explicit(Lattice &lat,int dir)
{

  const char *cname="none";
  const char *fname="rotate_gauge_explicit";
  VRB.Func(cname,fname);

  if ((dir<0)||(dir>3)){
    VRB.Result(cname,fname,"Error:: direction should be 0,1,2,3\n");
    return;
  }

  int NX(GJP.XnodeSites());
  int NY(GJP.YnodeSites());
  int NZ(GJP.ZnodeSites());
  int NT(GJP.TnodeSites());


  //---------------------------------------------------------------
  // Set up the gauge fixing and other required conditions
  //---------------------------------------------------------------

  int num_nodes[4]
    = { GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes() } ;
  
  int node_sites[4]
    = { GJP.XnodeSites(), GJP.YnodeSites(), GJP.ZnodeSites(), GJP.TnodeSites() } ;
  
  int Size[4];
  for (int i=0; i<4; i++){
    Size[i] = num_nodes[i]*node_sites[i];
  }
  
  int *plns;
  plns = (int*) smalloc(Size[dir]*sizeof(int));
  
  for (int i=0; i<Size[dir]; i++){
    plns[i] = i;
  }
  
  int npln = Size[dir];
  
  FixGaugeType normdir;
  
  if (dir==3) normdir = FIX_GAUGE_COULOMB_T;
  else if (dir==0) normdir = FIX_GAUGE_COULOMB_X;
  else if (dir==1) normdir = FIX_GAUGE_COULOMB_Y;
  else normdir = FIX_GAUGE_COULOMB_Z;
      

  //----------------------------------------------------------------------
  //initialize the parameters need to gauge fixing ---------------------
  //----------------------------------------------------------------------
  Matrix *L;
  Matrix **Gp;
  int ii;
  
  int volume = NX*NY*NZ*NT;
  
  L = (Matrix*) smalloc(4*volume*sizeof(Matrix));
  Gp = (Matrix**) smalloc(npln*sizeof(Matrix*));
  
  for(ii=0; ii<npln; ii++)
    Gp[ii] = (Matrix*) smalloc(volume/node_sites[dir] * sizeof(Matrix));
  
  
  //-----------------------------------------------------------------------------
  //GAUGE FIXING MATRICES
  //-----------------------------------------------------------------------------

  for (int slice=0; slice<node_sites[dir]; slice++)
    for(int cnt=0; cnt<volume/node_sites[dir]; cnt++)
      Gp[slice][cnt]=lat.FixGaugePtr()[slice][cnt];
  
  //-------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------
  // TRY TO Transform the Lattice to the Coulomb Gauge and store it ---------------
  //-------------------------------------------------------------------------------

  int tt,xx,yy,zz,temp_ind;
  int slice_ind[3];
  
  //slice_ind[3] stores the 3 directions on the 'dir' slice with indices increasing 
  temp_ind = 0;
  for (int i=0; i<4; i++){
    if (i!=dir) {
      slice_ind[temp_ind] = i;
      temp_ind++;
    }
  }
  
  VRB.Result(cname,fname,"dir == %d \n",dir);
  VRB.Result(cname,fname,"slice index == %d %d %d\n",slice_ind[0],slice_ind[1],slice_ind[2]);
  
  // the dummy node_sites for each dummy dirction
  int NN[4];
  NN[0] = node_sites[slice_ind[0]];
  NN[1] = node_sites[slice_ind[1]];
  NN[2] = node_sites[slice_ind[2]];
  NN[3] = node_sites[dir];
  
  int s[4];
  for (int i=0; i<3; i++)
    s[i] = slice_ind[i];
  s[3] = dir;

  //---------------------------------------------------------------------------------
  //copy the old lattice config to matrix array L, transform L and then copy back
  //---------------------------------------------------------------------------------
  
  lat.CopyGaugeField(L);  
  int x[4];
  
  // xx yy zz tt are dummy position vector, as tt represents the gfixing direction  
  for (x[3]=0; x[3]<node_sites[3]; x[3]++)
    for (x[2]=0; x[2]<node_sites[2]; x[2]++)
      for (x[1]=0; x[1]<node_sites[1]; x[1]++)
	for (x[0]=0; x[0]<node_sites[0]; x[0]++)
	  {
	    xx = x[slice_ind[0]];
	    yy = x[slice_ind[1]];
	    zz = x[slice_ind[2]];
	    tt = x[dir];
	
	    Matrix g = Gp[tt][(zz*NN[1]+yy)*NN[0]+xx];
	    Matrix D; 
	    
	    Matrix gg ;
	    Matrix transmit;
	    
	    //----------------- T Direction ----------------------------
	    if (tt+1<NN[3]) gg = Gp[tt+1][(zz*NN[1]+yy)*NN[0]+xx];
	    else { 
	      transmit = Gp[0][(zz*NN[1]+yy)*NN[0]+xx];
	      getPlusData((IFloat *)&gg, (IFloat *)&transmit,
			  sizeof(Matrix)/sizeof(IFloat), dir) ;
	    }
	    D.Dagger(gg);
	    L[IND(x[0],x[1],x[2],x[3],dir)] = g*L[IND(x[0],x[1],x[2],x[3],dir)]*D;
	    
	    //----------------- Z Direction ----------------------------
	    if (zz+1<NN[2]) gg = Gp[tt][((zz+1)*NN[1]+yy)*NN[0]+xx]; 
	    else {
	      transmit = Gp[tt][((zz+1)%NN[2]*NN[1]+yy)*NN[0]+xx]; 
	      getPlusData((IFloat *)&gg, (IFloat *)&transmit,
			  sizeof(Matrix)/sizeof(IFloat), slice_ind[2]) ;
	    }
	    D.Dagger(gg);
	    L[IND(x[0],x[1],x[2],x[3],slice_ind[2])] = g*L[IND(x[0],x[1],x[2],x[3],slice_ind[2])]*D;
	    
	    //----------------- Y Direction ----------------------------
	    if (yy+1<NN[1]) gg = Gp[tt][(zz*NN[1]+(yy+1))*NN[0]+xx];
	    else {
	      transmit = Gp[tt][(zz*NN[1]+(yy+1)%NN[1])*NN[0]+xx];
	      getPlusData((IFloat *)&gg, (IFloat *)&transmit,
			  sizeof(Matrix)/sizeof(IFloat), slice_ind[1]) ;
	    }
	    D.Dagger(gg);
	    L[IND(x[0],x[1],x[2],x[3],slice_ind[1])] = g*L[IND(x[0],x[1],x[2],x[3],slice_ind[1])]*D;
	    
	    //----------------- X Direction ----------------------------
	    if (xx+1<NN[0]) gg = Gp[tt][(zz*NN[1]+yy)*NN[0]+(xx+1)];
	    else {
	      transmit = Gp[tt][(zz*NN[1]+yy)*NN[0]+(xx+1)%NN[0]];
	      getPlusData((IFloat *)&gg, (IFloat *)&transmit,
			  sizeof(Matrix)/sizeof(IFloat), slice_ind[0]) ;
	    }
	    D.Dagger(gg);
	    L[IND(x[0],x[1],x[2],x[3],slice_ind[0])] = g*L[IND(x[0],x[1],x[2],x[3],slice_ind[0])]*D;
	  }
  
  
  lat.GaugeField(L);
  
  //--------------------------------Free L and Gp -------------------------------------
  sfree(L);
  
  for (ii=0; ii<npln; ii++)
    sfree(Gp[ii]);
  
  sfree(Gp);
  sfree(plns);

}


int IND(int x, int y, int z, int t, int l)
{
  int NX(GJP.XnodeSites());
  int NY(GJP.YnodeSites());
  int NZ(GJP.ZnodeSites());
  int NT(GJP.TnodeSites()); 
  return ((((((t+NT)%NT)*NZ+((z+NZ)%NZ))*NY+((y+NY)%NY))*NX+((x+NX)%NX))*4+l);
}
