#include <config.h>
#include <stdio.h>
#include <stdlib.h>
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
#include <alg/no_arg.h>
#include <alg/common_arg.h>
#include <alg/hmd_arg.h>
#include <alg/alg_plaq.h>
#include <alg/alg_rnd_gauge.h>
#include <alg/threept_arg.h>
#include <alg/threept_prop_arg.h>
#include <alg/alg_threept.h>

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

void chkpt(const int num_nodes,int& chkpoint_no,const Float dtime_start,Float& dtime_last,Float& dtime_now);

inline Matrix operator * (const Matrix& m1, const Matrix& m2)
{ Matrix r; r.DotMEqual(m1,m2); return r; }

void rotate_gauge_explicit(Lattice &lat,int dir=3);

int IND(int x, int y, int z, int t, int l);

int main(int argc,char *argv[])
{

  Start(&argc,&argv);

  VRB.Result("","main()","Starting next configuration.\n");

  //Initialize Timing
  Float dtime_start=dclock();
  Float dtime_last=dtime_start;
  Float dtime_now=dtime_start;

  int chkpoint_no=0;

  CommandLine::is(argc,argv);

  const char *cname="none";
  const char *fname="main";

  DoArg do_arg;
  EvoArg evo_arg;
  ThreePtArg threept_arg;
  QPropWArg qpropw_arg;


  /* get parameter from command line */

  /* 23 parameters: */
  /* x.x  
	  chkpoints (1 to turn checkpoints on, 0 to turn them off)
          DIRECTORY do_arg evo_arg threept_arg
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
	  contractions_done_list
	  save_light_props
	  save_strange_props
	  save_gauge_rot_lats
	  overwrite_previous_config_props
	  
 */

  int chkpoints(CommandLine::arg_as_int() );
 
  chdir(CommandLine::arg());
  
  
  if ( !do_arg.Decode(CommandLine::arg(),"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}
  if ( !evo_arg.Decode(CommandLine::arg(),"evo_arg") ) { printf("Bum evo_arg\n"); exit(-1);}
  if ( !threept_arg.Decode(CommandLine::arg(),"threept_arg") ) { printf("Bum threept_arg\n"); exit(-1);}


  GJP.Initialize(do_arg);


  const int num_nodes=GJP.Xnodes()*GJP.Ynodes()*GJP.Znodes()*GJP.Tnodes()*GJP.Snodes();

  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

  char lat_stand_in[200], lat_stand_out[200];
  char label[200], id[200];
  char logdir[200], plaquette_gfix_info_dir[200], midprop_contractions_dir[200], configs_done_log[200];
  char gauge_rotated_lats_dir[200], gauge_rotated_lats_saved_log[200];
  char contractions_done_list[200];

  sprintf(lat_stand_in, CommandLine::arg() );
  sprintf(lat_stand_out, CommandLine::arg() );
  
  int gauge_ran(CommandLine::arg_as_int() );
  int t_src(CommandLine::arg_as_int() );
  if ( t_src<0 || t_src>=GJP.Tnodes()*GJP.TnodeSites() ) {
    printf("t_src must be in the range 0 <= t_src < %d.\n",GJP.Tnodes()*GJP.TnodeSites());
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
  
  sprintf(contractions_done_list, CommandLine::arg() );

  int save_light_props(CommandLine::arg_as_int() );
  int save_strange_props(CommandLine::arg_as_int() );
  int save_gauge_rot_lats(CommandLine::arg_as_int() );
  int overwrite_previous_config_props(CommandLine::arg_as_int() );
  
  VRB.Result("","main()","This is configuration number %d.\n",seqNum);

  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

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
    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

  //If lat_stand_in is "none" then generate a random lattice, otherwise load
  //the file lat_stand_in.
  
  if (strcmp(lat_stand_in,"none")==0){
    VRB.Result("","main()","create random gauge - field (taking # %i) \n", gauge_ran);
    
    for(int ii(0); ii < gauge_ran; ++ii)
      lattice.SetGfieldDisOrd();
    
    VRB.Result("","main()","Created random gauge field.\n");
  } 
  else {
    ReadLatticeParallel readLat;
    
    VRB.Result("","main()","  reading: %s (NERSC-format)\n",lat_stand_in);

    readLat.read(lattice,lat_stand_in);

    VRB.Result("","main()","Lattice read.\n");
    
  }

  //If lat_stand_out is specified then save the non-gauge rotated lattice
  //to this file.

  if (strcmp(lat_stand_out,"none")!=0) {
      WriteLatticeParallel writeLat;

      VRB.Result("","main()","  writing: %s (NERSC-format)\n",lat_stand_out);

      QioArg qio_arg_lat_stand;
      qio_arg_lat_stand.init(lat_stand_out,evo_arg.io_concurrency,0.01,outformat,INT_AUTOMATIC,1); //Have to do it this way to pass it the io_concurrency from evo_arg
      writeLat.write(lattice, qio_arg_lat_stand);
  
      VRB.Result("","main()","Lattice written.\n");
  }

  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

  //Shift lattice so that time slice t_src moves to t=0
  int t_shift=-t_src/GJP.TnodeSites();
  VRB.Result("","main()","Shifting lattice by %d, command GDS.Set(0,0,0,%d).\n",-t_src,t_shift);
  GDS.Set(0,0,0,t_shift);
  lattice.Shift();
  VRB.Result("","main()","Lattice shifted.\n");
  
  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

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
    ERR.FileW("main","main", (char *)(gfix_fname) );
  }
  Fclose(fp);

  CommonArg common_arg_plaq;
  NoArg plaq_arg;
  char plaq_fname[200];
  sprintf(plaq_fname,"%s/plaq.%d.dat",plaquette_gfix_info_dir,seqNum);
  common_arg_plaq.results=plaq_fname;
  if( (fp = Fopen((char *)(plaq_fname),"w")) == NULL ){
    ERR.FileW("main","main", (char *)(plaq_fname) );
  }
  Fclose(fp);

  VRB.Result("","main()","Calculating plaquette.\n");
  AlgPlaq plaq(lattice,&common_arg_plaq,&plaq_arg);
  plaq.run();
  VRB.Result("","main()","Plaquette calculated.\n");

  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

  VRB.Result("","main()","Calculating gauge fixing matrices.\n");
  AlgFixGauge fix_gauge(lattice,&common_arg_gfix,&fix_arg);
  fix_gauge.run();
  VRB.Result("","main()","Gauge fixing matrices calculated.\n");

  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

  qpropw_arg.x=0;
  qpropw_arg.y=0;
  qpropw_arg.z=0;
  qpropw_arg.t=0; //Will also need source times tK[i] for the kaon correlator
  qpropw_arg.gauge_fix_src=1;
  qpropw_arg.gauge_fix_snk=1;
  qpropw_arg.store_midprop=0;
  qpropw_arg.do_half_fermion=0;	
  
  qpropw_arg.cg.max_num_iter=99999;
  qpropw_arg.cg.stop_rsd=1e-8;
  qpropw_arg.cg.true_rsd=1e-8;
  qpropw_arg.cg.RitzMatOper=MAT_HERM;
  qpropw_arg.cg.Inverter=CG;
  qpropw_arg.cg.bicgstab_n=0;

  //See which masses we need to calculate props for
  int num_light=threept_arg.num_light;
  Float l_mass[num_light];
  int num_strange=threept_arg.num_strange;
  Float s_mass[num_strange];
  for (int i=0; i<num_light; i++)
    l_mass[i]=threept_arg.l_mass[i];
  for (int i=0; i<num_strange; i++)
    s_mass[i]=threept_arg.s_mass[i];

  int num_tK_tmp=threept_arg.num_tK;
  int num_tK=num_tK_tmp;
  for (int i=0; i<num_tK_tmp; i++)
    if (threept_arg.tK[i]==0)
      num_tK--;
  int tK[num_tK];
  int ii=0;
  for (int i=0; i<num_tK_tmp; i++) {
    if (threept_arg.tK[i]!=0) {
      tK[ii]=threept_arg.tK[i];
      ii++;
    } else {
      VRB.Result("","main()","threept_arg.tK[%d]=0.  Won't do repeat of kaon with source at 0.\n",i);
    }
  }
  //Now fix up threept_arg.tK so that there are no zeroes.
  for (int i=0; i<num_tK; i++)
    threept_arg.tK[i]=tK[i];
  threept_arg.num_tK=num_tK;

  char logfile[200];
  sprintf(logfile,"%s/%d",logdir,seqNum);
  if( (fp = fopen((char *)(logfile),"r")) == NULL ){
    ERR.FileR("main","main", (char *)(logfile) );
  }
  float mass_read_tmp;
  Float mass_read;
  int bc;
  int mom_num;
  int mom_dir;
  int tsrc_read;
  int l_mass_tpi_done[num_light][2][4][3]; //mass number, bc, momentum number (number of twists), momentum direction
  int s_mass_tpi_done[num_strange][2]; //mass number, bc
  int l_mass_tK_done[num_light][2][num_tK]; //mass number, bc, source time
  int s_mass_tK_done[num_strange][2][num_tK]; //mass number, bc, source time
  for (int i=0; i<num_light; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<4; k++)
	for (int l=0; l<3; l++)
	  l_mass_tpi_done[i][j][k][l]=0;
  for (int i=0; i<num_strange; i++)
    for (int j=0; j<2; j++)
      s_mass_tpi_done[i][j]=0;
  for (int i=0; i<num_light; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<num_tK; k++)
	l_mass_tK_done[i][j][k]=0;
  for (int i=0; i<num_strange; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<num_tK; k++)
	s_mass_tK_done[i][j][k]=0;
  VRB.Result("","main()","Finished masses, boundary conditions, momenta/twists, and source times as read from log (config no. %d):\n",seqNum);
  while (fscanf(fp,"%f",&mass_read_tmp)>=0){
    mass_read=(Float) mass_read_tmp;
    fscanf(fp,"%d",&bc);
    fscanf(fp,"%d",&mom_num);
    fscanf(fp,"%d",&mom_dir);
    fscanf(fp,"%d",&tsrc_read);
    VRB.Result("","main()","f %d %d %d %d ",mass_read,bc,mom_num,mom_dir,tsrc_read);
    if (bc<0 || bc>1) {
      VRB.Result("","main()","n");
      continue;
    }
    if (mom_num<0 || mom_num>3) {
      VRB.Result("","main()","n");
      continue;
    }
    if (mom_dir<0 || mom_dir>2) {
      VRB.Result("","main()","n");
      continue;
    }
    if ( tsrc_read<0 || tsrc_read>=GJP.Tnodes()*GJP.TnodeSites() ) {
      VRB.Result("","main()","n");
      continue;
    }
    for (int i=0; i<num_light; i++) {
      if (fabs(l_mass[i]-mass_read)<=1e-6) {
	if (tsrc_read==0) {
	  l_mass_tpi_done[i][bc][mom_num][mom_dir]=1;
	  if (mom_num==0 || mom_num==3) //no direction for 0 momentum or three twists,
	                                //set all directions to "done"
	    for (int tmp_dir=0; tmp_dir<3; tmp_dir++)
	      l_mass_tpi_done[i][bc][mom_num][tmp_dir]=1;
	} else { //tsrc_read if statement
	  for (int j=0; j<num_tK; j++)
	    if (tK[j]==tsrc_read)
	      l_mass_tK_done[i][bc][j]=1;
	} //tsrc_read if statement
	VRB.Result("","main()","light mass) ");
      } //mass_read if statement
    } //i for loop
    for (int i=0; i<num_strange; i++) {
      if (fabs(s_mass[i]-mass_read)<=1e-6) {
	if (mom_num==0 && mom_dir==0) {
	  if (tsrc_read==0)
	    s_mass_tpi_done[i][bc]=1;
	  else
	    for (int j=0; j<num_tK; j++)
	      if (tK[j]==tsrc_read)
		s_mass_tK_done[i][bc][j]=1;
	} //mom_num, mom_dir if statement
	VRB.Result("","main()","strange mass) ");
      } //mass_read if statement
    } //i for loop
    VRB.Result("","main()","n");
  } //fscanf while loop
  fclose(fp);
  VRB.Result("","main()","Masses, boundary conditions, momenta/twists and source times to do (config. no. %d):\n",seqNum);
  int all_masses_done=1;
  int do_first_mom=threept_arg.do_first_mom;
  int do_second_mom=threept_arg.do_second_mom;
  int do_third_mom=threept_arg.do_third_mom;
  int do_p_plus_a_kaon=threept_arg.do_p_plus_a_kaon;
  int do_kaon_at_walls=threept_arg.do_kaon_at_walls;
  int do_kaons_tK=threept_arg.do_kaons_tK;
  VRB.Result("","main()","Light masses\n");
  for (int i=0; i<num_light; i++) {
    for (bc=0; bc<2; bc++) {
      //Zero Momentum
      //always need these even for non-zero momentum calculations
      mom_num=0;
      mom_dir=0;
      if (!l_mass_tpi_done[i][bc][mom_num][mom_dir]) {
	VRB.Result("","main()","f %d %d %d 0\n",l_mass[i],bc,mom_num,mom_dir);
	all_masses_done=0;
      }
      //One Twist
      if (do_first_mom) {
	mom_num=1;
	for (mom_dir=0; mom_dir<3; mom_dir++)
	  if (!l_mass_tpi_done[i][bc][mom_num][mom_dir]) {
	    VRB.Result("","main()","f %d %d %d 0\n",l_mass[i],bc,mom_num,mom_dir);
	    all_masses_done=0;
	  }
      } //do_first_mom if statement
      //Two Twists
      if (do_second_mom) {
	mom_num=2;
	for (mom_dir=0; mom_dir<3; mom_dir++)
	  if (!l_mass_tpi_done[i][bc][mom_num][mom_dir]) {
	    VRB.Result("","main()","f %d %d %d 0\n",l_mass[i],bc,mom_num,mom_dir);
	    all_masses_done=0;
	  }
      } //do_second_mom if statement
      //Three Twists
      if (do_third_mom) {
	mom_num=3;
	mom_dir=0;
	if (!l_mass_tpi_done[i][bc][mom_num][mom_dir]) {
	  VRB.Result("","main()","f %d %d %d 0\n",l_mass[i],bc,mom_num,mom_dir);
	  all_masses_done=0;
	}
      } //do_third_mom if statement
      //Light quark propagators with source at tK[j]
      if (do_kaons_tK) {
	for (int j=0; j<num_tK; j++) {
	  if (!l_mass_tK_done[i][bc][j]) {
	    VRB.Result("","main()","f %d 0 0 %d\n",l_mass[i],bc,tK[j]);
	    all_masses_done=0;
	  } //l_mass_tK_done if statement
	} //j for loop
      } //do_kaons_tK if statement
    } //bc for loop
  } //i for loop
  VRB.Result("","main()","Strange masses\n");
  for (int i=0; i<num_strange; i++) {

    if (do_kaon_at_walls) {
      for (bc=0; bc<2; bc++) {
	if (!s_mass_tpi_done[i][bc]) {
	  VRB.Result("","main()","f %d 0 0 0\n",s_mass[i],bc);
	  all_masses_done=0;
	} //s_mass_tpi_done if statement
      } //bc for loop
    } //do_kaon_at_walls if statement

    int bc_num=1+do_p_plus_a_kaon;
    for (bc=0; bc<bc_num; bc++) {
      for (int j=0; j<num_tK; j++) {
	if (!s_mass_tK_done[i][bc][j]) {
	  VRB.Result("","main()","f %d 0 0 %d\n",s_mass[i],bc,tK[j]);
	  all_masses_done=0;
	} //s_mass_tK_done if statement
      } //j for loop
    } //bc for loop

  } //i for loop
  
  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

  //Check if contractions have been done for this configuration
  //and for each of the light quark masses
  if( (fp = fopen((char *)(contractions_done_list),"r")) == NULL ){
    ERR.FileR("main","main", (char *)(contractions_done_list) );
  }
  int config_no_read;
  int no_contractions_done=1;
  int l_mass_contracted[num_light];
  for (int i=0; i<num_light; i++)
    l_mass_contracted[i]=0;
  VRB.Result("","main()","Contractions done for light quark masses:\n");
  while (fscanf(fp,"%d",&config_no_read)>=0){
    fscanf(fp,"%f",&mass_read_tmp);
    mass_read=(Float) mass_read_tmp;
    if (config_no_read==seqNum){
      for (int i=0; i<num_light; i++) {
	if (fabs(l_mass[i]-mass_read)<=1e-6) {
	  l_mass_contracted[i]=1;
	  no_contractions_done=0;
	  VRB.Result("","main()","%f\n",l_mass[i]);
	}  //l_mass[i]-mass_read if statement
      } //i for loop
    } //config_no_read if statement
  } //fscanf while loop
  fclose(fp);
  //See if any light masses still need to be contracted.
  int all_contractions_done=1;
  VRB.Result("","main()","Light quark masses for which contractions must be done:\n");
  for (int i=0; i<num_light; i++) {
    if (!l_mass_contracted[i]) {
      VRB.Result("","main()","f\n",l_mass[i]);
      all_contractions_done=0;
    }
  }

  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);


  //Set the threept_arg correctly for it to be passed to alg_threept
  char wme_results[num_light][200], wme_results_mres_ZA[num_light][200], wme_results_pipi[num_light][200];
  for (int i_light=0; i_light<num_light; i_light++) {
    sprintf(wme_results[i_light],"%s.ml%0.4f.config%d",threept_arg.results,l_mass[i_light],seqNum);
    sprintf(wme_results_mres_ZA[i_light],"%s.ml%0.4f.config%d",threept_arg.results_mres_ZA,l_mass[i_light],seqNum);
    sprintf(wme_results_pipi[i_light],"%s.ml%0.4f.config%d",threept_arg.results_pipi,l_mass[i_light],seqNum);
  }
  threept_arg.t_src=0; //remember that the threept_arg.tK array contains
                       //the non-zero source times of some of the propagators
  threept_arg.t_snk=GJP.TnodeSites()*GJP.Tnodes();
  threept_arg.t_shift=0; //Shift of lattice already done, don't need to do it in alg_threept.
  threept_arg.ensemble_label=label;
  threept_arg.ensemble_id=id;
  threept_arg.seqNum=seqNum;

  //Delete the contents (if any) of the files that threept_arg will write to.
  for (int i_light=0; i_light<num_light; i_light++) {
    if (!l_mass_contracted[i_light]) {
      FILE *fp_tmp;
      if( (fp_tmp = Fopen((char *)(wme_results[i_light]),"w")) == NULL ){
	ERR.FileW("main","main", (char *)(wme_results[i_light]) );
      }
      Fclose(fp_tmp);
      if( (fp_tmp = Fopen((char *)(wme_results_mres_ZA[i_light]),"w")) == NULL ){
	ERR.FileW("main","main", (char *)(wme_results_mres_ZA[i_light]) );
      }
      Fclose(fp_tmp);
      if( (fp_tmp = Fopen((char *)(wme_results_pipi[i_light]),"w")) == NULL ){
	ERR.FileW("main","main", (char *)(wme_results_pipi[i_light]) );
      }
      Fclose(fp_tmp);
    }
  }
  


  if ( !all_masses_done || !all_contractions_done) {

    //Loop over masses and calculate props.  Reload props if already done.
    CommonArg common_arg;
    char prop_name[200], midprop_contractions_fname[200], qio_filename[200];
    char prop_name_w_config_no[200];
    VRB.Result("","main()","Begin mass loop for configuration number %d.\n",seqNum);
    QPropW* q_light_tpi[1][2][4][3];
    QPropW* q_light_tK[1][2][num_tK];
    QPropW* q_strange_tpi[num_strange][2];
    QPropW* q_strange_tK[num_strange][2][num_tK];

    //Set all pointers to null.
    for (int bc_tmp=0; bc_tmp<2; bc_tmp++) {
      for (int mom_num_tmp=0; mom_num_tmp<4; mom_num_tmp++)
	for (int mom_dir_tmp=0; mom_dir_tmp<3; mom_dir_tmp++)
	  q_light_tpi[0][bc_tmp][mom_num_tmp][mom_dir_tmp]=0;
      for (int j_tmp=0; j_tmp<num_tK; j_tmp++)
	q_light_tK[0][bc_tmp][j_tmp]=0;
    }
    for (int i_tmp=0; i_tmp<num_strange; i_tmp++)
      for (int bc_tmp=0; bc_tmp<2; bc_tmp++) {
	q_strange_tpi[i_tmp][bc_tmp]=0;
	for (int j_tmp=0; j_tmp<num_tK; j_tmp++)
	  q_strange_tK[i_tmp][bc_tmp][j_tmp]=0;
      }


    //Need to calculate strange quark propagators once and only once.
    //Use flags to ensure this despite the loop over light quark masses below.
    int s_mass_tpi_in_mem[num_strange][2]; //mass number, bc
    int s_mass_tK_in_mem[num_strange][2][num_tK]; //mass number, bc, source time
    for (int i_tmp=0; i_tmp<num_strange; i_tmp++)
      for (int bc_tmp=0; bc_tmp<2; bc_tmp++) {
	s_mass_tpi_in_mem[i_tmp][bc_tmp]=0;
	for (int j_tmp=0; j_tmp<num_tK; j_tmp++)
	  s_mass_tK_in_mem[i_tmp][bc_tmp][j_tmp]=0;
      }


    //Loop over light quark masses and for each light quark mass generate the
    //propagators and do the contractions.  Then delete light quark propagators
    //and move onto the next light quark mass.
    for (int i_light=0; i_light<num_light; i_light++) {
      
      for (int q_type=0; q_type<2; q_type++) {
	
	int num_q;
	if (q_type) {
	  num_q=num_strange;
	  qpropw_arg.save_prop=save_strange_props;
	} else {
	  num_q=1;
	  qpropw_arg.save_prop=save_light_props;
	}
	
	for (int src_type=0; src_type<2; src_type++) { //Whether doing q_(...)_tpi (src_type=0)
	                                               //or q_(...)_tK (src_type=1) propagators
	  
	  int n_mom;
	  if (q_type==0 && src_type==0)
	    n_mom=4;
	  else
	    n_mom=1;
	  
	  int bc_num=2-src_type*q_type*(1-do_p_plus_a_kaon);
	  int num_src_times;
	  if (src_type)
	    num_src_times=num_tK;
	  else
	    num_src_times=1;
	  
	  if (src_type==0 && q_type==1 && !do_kaon_at_walls) continue;
	  if (src_type==1 && q_type==0 && !do_kaons_tK) continue;
	  
	  for (int src_num=0; src_num<num_src_times; src_num++) {
	    
	    int src_time;
	    if (src_type)
	      src_time=tK[src_num];
	    else
	      src_time=0;
	    qpropw_arg.t=src_time;
	    
	    char source_info[30];
	    if (src_type)
	      sprintf(source_info,"tsrc%d_",src_time);
	    else
	      sprintf(source_info,"");
	    char source_string[100];
	    sprintf(source_string,"source at t=%d ",src_time);
	    
	    
	    for (int i=0; i<num_q; i++) {
	      
	      Float this_mass;
	      if (q_type)
		this_mass=s_mass[i];
	      else
		this_mass=l_mass[i_light];
	      qpropw_arg.cg.mass=this_mass;
	      
	      for (bc=0; bc<bc_num; bc++) {
		
		char bc_type[30], bc_label[30];
		if (bc) {
		  sprintf(bc_type,"antiperiodic");
		  sprintf(bc_label,"APRD");
		  GJP.Tbc(BND_CND_APRD);
		} else {
		  sprintf(bc_type,"periodic");
		  sprintf(bc_label,"PRD");
		  GJP.Tbc(BND_CND_PRD);
		}
		
		for (mom_num=0; mom_num<n_mom; mom_num++) {
		  
		  if ( !do_first_mom && mom_num==1 ) continue;
		  if ( !do_second_mom && mom_num==2 ) continue;
		  if ( !do_third_mom && mom_num==3 ) continue;
		  
		  int n_dir;
		  if (mom_num==0 || mom_num==3)
		    n_dir=1;
		  else
		    n_dir=3;
		  
		  for (mom_dir=0; mom_dir<n_dir; mom_dir++) {
		    
		    char mom_info[30];
		    char mom_string[30];
		    char err_string[100];
		    GJP.Xbc(BND_CND_PRD);
		    GJP.Ybc(BND_CND_PRD);
		    GJP.Zbc(BND_CND_PRD);
		    if (mom_num==0) {
		      sprintf(mom_info,"");
		      sprintf(mom_string,"");
		    } else if (mom_num==1) {
		      if (mom_dir==0) {
			sprintf(mom_info,"_Twist_px");
			sprintf(mom_string,"momentum Pi/L in x direction, ");
			GJP.Xbc(BND_CND_APRD);
		      } else if (mom_dir==1) {
			sprintf(mom_info,"_Twist_py");
			sprintf(mom_string,"momentum Pi/L in y direction, ");
			GJP.Ybc(BND_CND_APRD);
		      } else if (mom_dir==2) {
			sprintf(mom_info,"_Twist_pz");
			sprintf(mom_string,"momentum Pi/L in z direction, ");
			GJP.Zbc(BND_CND_APRD);
		      } else {
			sprintf(err_string,"mom_dir should be 0,1, or 2, but mom_dir=%d",mom_dir);
			ERR.General(cname,fname,err_string);
		      }
		    } else if (mom_num==2) {
		      if (mom_dir==0) {
			sprintf(mom_info,"_Twist_pypz");
			sprintf(mom_string,"momentum Pi/L in y and z direction, ");
			GJP.Ybc(BND_CND_APRD);
			GJP.Zbc(BND_CND_APRD);
		      } else if (mom_dir==1) {
			sprintf(mom_info,"_Twist_pxpz");
			sprintf(mom_string,"momentum Pi/L in x and z direction, ");
			GJP.Xbc(BND_CND_APRD);
			GJP.Zbc(BND_CND_APRD);
		      } else if (mom_dir==2) {
			sprintf(mom_info,"_Twist_pxpy");
			sprintf(mom_string,"momentum Pi/L in x and y direction, ");
			GJP.Xbc(BND_CND_APRD);
			GJP.Ybc(BND_CND_APRD);
		      } else {
			sprintf(err_string,"mom_dir should be 0,1, or 2, but mom_dir=%d",mom_dir);
			ERR.General(cname,fname,err_string);
		      }
		    } else if (mom_num==3) {
		      sprintf(mom_info,"_Twist_pxpypz");
		      sprintf(mom_string,"momentum Pi/L in x, y, and z direction, ");
		      GJP.Xbc(BND_CND_APRD);
		      GJP.Ybc(BND_CND_APRD);
		      GJP.Zbc(BND_CND_APRD);
		    } else {
		      sprintf(err_string,"mom_num should be 0,1,2, or 3, but mom_num=%d\n",mom_num);
		      ERR.General(cname,fname,err_string);
		    }
		    
		    sprintf(prop_name_w_config_no,"24cube_b0.87_DBW2_traj%d_m%0.4f_tshift%d_%s%s%s",seqNum,this_mass,-t_src,source_info,bc_label,mom_info); //Note that the variable t_src is w.r.t. the unshifted lattice in the file name, whereas the source time in source_info refers to the value of qpropw_arg.t which is what will appear in the qio header
		    if (overwrite_previous_config_props)
		      sprintf(prop_name,"24cube_b0.87_DBW2_m%0.4f_tshift%d_%s%s%s",this_mass,-t_src,source_info,bc_label,mom_info);
		    else
		      sprintf(prop_name,"%s",prop_name_w_config_no);
//		    sprintf(qio_filename,"/pfs/qio_prop.%s",prop_name);
		    sprintf(qio_filename,"%d/qio_prop.%s",seqNum,prop_name);
		    qpropw_arg.file=qio_filename;
		    qpropw_arg.ensemble_id=id;
		    qpropw_arg.ensemble_label=label;
		    qpropw_arg.seqNum=seqNum;

		    //Test whether a propagator is on disk to be reloaded. This
		    //does not test whether it is already in memory (which is 
		    //possible for strange quarks).
		    int testval;
		    if (q_type) {
		      if (src_type)
			testval=s_mass_tK_done[i][bc][src_num];
		      else
			testval=s_mass_tpi_done[i][bc];
		    } else {
		      if (src_type)
			testval=l_mass_tK_done[i_light][bc][src_num];
		      else
			testval=l_mass_tpi_done[i_light][bc][mom_num][mom_dir];
		    }
		    
		    if (testval) {
		      
		      if (!l_mass_contracted[i_light]) { //Only need to reload if actually doing contractions
			
			//Test to see if the propagator is already in memory.
			//(Can't be true for a light quark).
			int testval2;
			if (q_type) {
			  if (src_type)
			    testval2=s_mass_tK_in_mem[i][bc][src_num];
			  else
			    testval2=s_mass_tpi_in_mem[i][bc];
			} else {
			  testval2=0;
			}

			if (!testval2) { //Only need to reload if not already in memory.
			  
			  VRB.Result("","main()","Already inverted mass %f, %s boundary conditions, %s%sconfiguration number %d.  Reloading.\n",this_mass,bc_type,mom_string,source_string,seqNum);
			  //Checkpoint
			  if (chkpoints)
			    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
			  
			  if (q_type) {
			    if (src_type) {
			      q_strange_tK[i][bc][src_num] = new QPropW(lattice, &qpropw_arg, &common_arg);
			      q_strange_tK[i][bc][src_num]->ReLoad(qio_filename);
			    } else {
			      q_strange_tpi[i][bc] = new QPropW(lattice, &qpropw_arg, &common_arg);
			      q_strange_tpi[i][bc]->ReLoad(qio_filename);
			    }
			  } else {
			    if (src_type) {
			      q_light_tK[0][bc][src_num] = new QPropW(lattice, &qpropw_arg, &common_arg);
			      q_light_tK[0][bc][src_num]->ReLoad(qio_filename);
			    } else {
			      q_light_tpi[0][bc][mom_num][mom_dir] = new QPropW(lattice, &qpropw_arg, &common_arg);
			      q_light_tpi[0][bc][mom_num][mom_dir]->ReLoad(qio_filename);
			    }
			  }
			  VRB.Result("","main()","Reloaded mass %f, %s boundary conditions, %s%sconfiguration number %d.\n",this_mass,bc_type,mom_string,source_string,seqNum);
			  //Checkpoint
			  if (chkpoints)
			    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

			  //Record that we now have the propagator in memory if
			  //it's a strange quark.
			  if (q_type) {
			    if (src_type)
			      s_mass_tK_in_mem[i][bc][src_num]=1;
			    else
			      s_mass_tpi_in_mem[i][bc]=1;
			  }

			} //mass_in_mem (testval2) if statement
		      } //l_mass_contracted[i_light] if statement
		      
		    } else { //mass_done (testval) if statement. This mass not done, must invert.
		      
		      if (!l_mass_contracted[i_light] || qpropw_arg.save_prop) { //Only invert if we need it.

			//Test to see if the propagator is already in memory.
			//(Can't be true for a light quark).
			int testval2;
			if (q_type) {
			  if (src_type)
			    testval2=s_mass_tK_in_mem[i][bc][src_num];
			  else
			    testval2=s_mass_tpi_in_mem[i][bc];
			} else {
			  testval2=0;
			}
			
			if (!testval2) { //Only need to reload if not already in memory.
			  
			  VRB.Result("","main()","Inverting mass %f, %s boundary conditions, %s%sconfiguration number %d.\n",this_mass,bc_type,mom_string,source_string,seqNum);
			  //Checkpoint
			  if (chkpoints)
			    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
			  
			  sprintf(midprop_contractions_fname,"%s/midprop_contractions.%s",midprop_contractions_dir,prop_name_w_config_no);
			  if( (fp = Fopen((char *)(midprop_contractions_fname),"w")) == NULL ){
			    ERR.FileW("main","main", (char *)(midprop_contractions_fname) );
			  }
			  Fclose(fp);
			  common_arg.results=midprop_contractions_fname;
			  if (q_type) {
			    if (src_type)
			      q_strange_tK[i][bc][src_num] = new QPropWWallSrc(lattice, &qpropw_arg, &common_arg);
			    else
			      q_strange_tpi[i][bc] = new QPropWWallSrc(lattice, &qpropw_arg, &common_arg);
			  } else { //q_type if statement
			    if (src_type)
			      q_light_tK[0][bc][src_num] = new QPropWWallSrc(lattice, &qpropw_arg, &common_arg);
			    else if (mom_num==0) { //src_type if statement
			      q_light_tpi[0][bc][0][0] = new QPropWWallSrc(lattice, &qpropw_arg, &common_arg);
			    } else { //src_type if statement
			      int p[3];
			      if (mom_num==1) {
				p[0]=p[1]=p[2]=0;
				p[mom_dir]=1;
			      } else if (mom_num==2) {
				for (int j=0; j<3; j++)
				  if (j==mom_dir)
				    p[j]=0;
				  else
				    p[j]=1;
			      } else { //mom_num if statement
				p[0]=p[1]=p[2]=1;
			      } //mom_num if statement
			      q_light_tpi[0][bc][mom_num][mom_dir] = new QPropWMomCosTwistSrc(lattice, &qpropw_arg, p, &common_arg);
			    } //src_type if statement
			  } //q_type if statement
			  VRB.Result("","main()","Finished inversion, mass %f, %s boundary conditions, %s%sconfiguration number %d.\n",this_mass,bc_type,mom_string,source_string,seqNum);
			  //Checkpoint
			  if (chkpoints)
			    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
			  
			  if (qpropw_arg.save_prop) {
			    if( (fp = Fopen((char *)(logfile),"a")) == NULL ){
			      ERR.FileA("main","main", (char *)(logfile) );
			    }
			    Fprintf(fp,"%0.4f %d %d %d %d\n",this_mass,bc,mom_num,mom_dir,src_time);
			    Fclose(fp);
			  }

			  //Record that we now have the propagator in memory if
			  //it's a strange quark.
			  if (q_type) {
			    if (src_type)
			      s_mass_tK_in_mem[i][bc][src_num]=1;
			    else
			      s_mass_tpi_in_mem[i][bc]=1;
			  }
			  
			} //mass_in_mem (testval2) if statement
			
		      } //!l_mass_contracted[i_light] || qpropw_arg.save_prop if statement
		      
		    } //mass_done (testval) if statement
		    
		  } //mom_dir for loop
		  
		}  //mom_num for loop
		
	      } //bc for loop
	      
	    } //i for loop
	    
	  } //src_num for loop
	  
	} //src_type for loop
	
      } //q_type for loop
      
      //Do contractions for this light mass if not done.
      if (!l_mass_contracted[i_light]) {
	
	VRB.Result("","main()","Starting threept calculations for configuration number %d, light mass %f.\n",seqNum,l_mass[i_light]);
	
	//The following is what alg_threept should see.
	threept_arg.num_light=1;
	threept_arg.l_mass[0]=l_mass[i_light];
	threept_arg.results=wme_results[i_light];
	threept_arg.results_mres_ZA=wme_results_mres_ZA[i_light];
	threept_arg.results_pipi=wme_results_pipi[i_light];
	
	common_arg.results=threept_arg.results;

	//Set the threept_prop_arg.
	ThreePtPropArg threept_prop_arg;
	for (bc=0; bc<2; bc++) {
	  for (mom_num=0; mom_num<4; mom_num++) {
	    
	    if ( !do_first_mom && mom_num==1 ) continue;
	    if ( !do_second_mom && mom_num==2 ) continue;
	    if ( !do_third_mom && mom_num==3 ) continue;
	    
	    int n_dir;
	    if (mom_num==0 || mom_num==3)
	      n_dir=1;
	    else
	      n_dir=3;
	    
	    for (mom_dir=0; mom_dir<n_dir; mom_dir++)
	      threept_prop_arg.q_light_tpi[0][bc][mom_num][mom_dir]=q_light_tpi[0][bc][mom_num][mom_dir];
	    
	  } //mom_num for loop
	  
	  if (do_kaons_tK)
	    for (int j=0; j<num_tK; j++)
	      threept_prop_arg.q_light_tK[0][bc][j]=q_light_tK[0][bc][j];
	  
	} //bc for loop
	
	for (int i=0; i<num_strange; i++) {
	  if (do_kaon_at_walls)
	    for (bc=0; bc<2; bc++)
	      threept_prop_arg.q_strange_tpi[i][bc]=q_strange_tpi[i][bc];
	  int bc_num=1+do_p_plus_a_kaon;
	  for (bc=0; bc<bc_num; bc++)
	    for (int j=0; j<num_tK; j++)
	      threept_prop_arg.q_strange_tK[i][bc][j]=q_strange_tK[i][bc][j];
	}
	
	//Feed everything to alg_threept to do contractions.
	//Checkpoint
	if (chkpoints)
	  chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
	AlgThreePt alg_threept(lattice,&common_arg,&threept_arg,&threept_prop_arg);
	alg_threept.run();
	VRB.Result("","main()","Finished threept calculations for configuration number %d, light mass %f.\n",seqNum,l_mass[i_light]);
	//Checkpoint
	if (chkpoints)
	  chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
	if( (fp = Fopen((char *)(contractions_done_list),"a")) == NULL ){
	  ERR.FileA("main","main", (char *)(contractions_done_list) );
	}
	Fprintf(fp,"%d %0.4f\n",seqNum,l_mass[i_light]);
	Fclose(fp);

	//Undo the changes to the strange quark propagators made by alg_threept.
	for (int m=0; m<num_strange; m++) {
	  
	  //tpi strange quarks
	  //q_strange_tpi[m][0] is (P+A)/2 and q_strange_tpi[m][1] is (P-A)/2
	  
	  if (do_kaon_at_walls) {
	    q_strange_tpi[m][0]->LinComb(*q_strange_tpi[m][1],1.0,1.0); //q_strange_tpi[m][0] is now P
	    q_strange_tpi[m][1]->LinComb(*q_strange_tpi[m][0],-2.0,1.0); //q_strange_tpi[m][1] is now A
	  }

	  //tK strange quarks
	  //q_strange_tK[m][0][j_tmp] is (P+A)/2 and q_strange_tK[m][1][j_tmp] is A
	  
	  if (do_p_plus_a_kaon) {
	    for (int j_tmp=0; j_tmp<num_tK; j_tmp++) {
	      q_strange_tK[m][0][j_tmp]->LinComb(*q_strange_tK[m][1][j_tmp],2.0,-1.0); //q_strange_tK[m][0][j_tmp] is now P
	    }
	  }

	} //m for loop, strange quark mass
	    
	
      } //!l_mass_contracted[i_light] if statement
      
      //Delete light quark propagators for this light mass
      //(set pointer to null once deleted)
      for (int bc_tmp=0; bc_tmp<2; bc_tmp++) {
	for (int mom_num_tmp=0; mom_num_tmp<4; mom_num_tmp++)
	  for (int mom_dir_tmp=0; mom_dir_tmp<3; mom_dir_tmp++) {
	    delete q_light_tpi[0][bc_tmp][mom_num_tmp][mom_dir_tmp];
	    q_light_tpi[0][bc_tmp][mom_num_tmp][mom_dir_tmp]=0;
	  }
	for (int j_tmp=0; j_tmp<num_tK; j_tmp++) {
	  delete q_light_tK[0][bc_tmp][j_tmp];
	  q_light_tK[0][bc_tmp][j_tmp]=0;
	}
      }

      
    } //i_light for loop

    VRB.Result("","main()","Finished mass loop for configuration number %d.\n",seqNum);

    //Delete remaining propagators to avoid memory leak.
    for (int bc_tmp=0; bc_tmp<2; bc_tmp++) {
      for (int mom_num_tmp=0; mom_num_tmp<4; mom_num_tmp++)
	for (int mom_dir_tmp=0; mom_dir_tmp<3; mom_dir_tmp++) {
	  delete q_light_tpi[0][bc_tmp][mom_num_tmp][mom_dir_tmp];
	  q_light_tpi[0][bc_tmp][mom_num_tmp][mom_dir_tmp]=0;
	}
      for (int j_tmp=0; j_tmp<num_tK; j_tmp++) {
	delete q_light_tK[0][bc_tmp][j_tmp];
	q_light_tK[0][bc_tmp][j_tmp]=0;
      }
    }
    for (int i_tmp=0; i_tmp<num_strange; i_tmp++)
      for (int bc_tmp=0; bc_tmp<2; bc_tmp++) {
	delete q_strange_tpi[i_tmp][bc_tmp];
	q_strange_tpi[i_tmp][bc_tmp]=0;
	for (int j_tmp=0; j_tmp<num_tK; j_tmp++) {
	  delete q_strange_tK[i_tmp][bc_tmp][j_tmp];
	  q_strange_tK[i_tmp][bc_tmp][j_tmp]=0;
	}
      }
    

  } //all_masses_done and all_contractions_done if statement

  if (save_gauge_rot_lats) {
    
    //Gauge rotate and save the lattice if this hasn't been done already
    //First check if the gauge rotated lattice has already been saved
    if( (fp = fopen((char *)(gauge_rotated_lats_saved_log),"r")) == NULL ){
      ERR.FileR("main","main", (char *)(gauge_rotated_lats_saved_log) );
    }
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
      VRB.Result("","main()","Gauge rotating the lattice and saving it.\n");
      //Gauge rotate lattice
      VRB.Result("","main()","Gauge rotating lattice.\n");
      rotate_gauge_explicit(lattice);
      VRB.Result("","main()","Gauge rotated lattice.\n");
      //Shift the lattice back
      VRB.Result("","main()","Shifting lattice back.  Shifting by %d, command GDS.Set(0,0,0,%d).\n",t_src,-t_shift);
      GDS.Set(0,0,0,-t_shift);
      lattice.Shift();
      VRB.Result("","main()","Lattice shifted.\n");
      //Save it
      char rotated_lat_fname[200];
      sprintf(rotated_lat_fname,"%s/ckpoint_lat_gauge_rotated.%d",gauge_rotated_lats_dir,seqNum);
      VRB.Result("","main()","Writing gauge rotated lattice.\n");
      WriteLatticeParallel writeLat;
      QioArg qio_arg;
      qio_arg.init(rotated_lat_fname,evo_arg.io_concurrency,0.01,outformat,INT_AUTOMATIC,1); //Have to do it this way to pass it the io_concurrency from evo_arg
      writeLat.write(lattice, qio_arg);
      VRB.Result("","main()","Wrote gauge rotated lattice.\n");
      VRB.Result("","main()","Finished gauge rotating and saving lattice.\n");
      //Checkpoint
      if (chkpoints)
	chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
      //Record in log that this lattice has been saved
      if( (fp = Fopen((char *)(gauge_rotated_lats_saved_log),"a")) == NULL ){
	ERR.FileA("main","main", (char *)(gauge_rotated_lats_saved_log) );
      }
      Fprintf(fp,"%d\n",seqNum);
      Fclose(fp);
    } //config_saved if statement
    
  } //save_gauge_rot_lats if statement
  
  if( (fp = Fopen((char *)(configs_done_log),"a")) == NULL ){
    ERR.FileA("main","main", (char *)(configs_done_log) );
  }
  Fprintf(fp,"%d\n",seqNum);
  Fclose(fp);


  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
  

  fix_gauge.free();

  VRB.Result("","main()","Finished configuration number %d.\n",seqNum);
  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);


  return 0;

}

//---------------------------------------------------------------------------
//Checksum, machine synchronization and timing checkpoint function
//---------------------------------------------------------------------------
void chkpt(const int num_nodes,int& chkpoint_no,const Float dtime_start,Float& dtime_last,Float& dtime_now)
{
  ERR.HdwCheck("main","chkpt");
  int dummy=0;
  Float test_val;
  Float* ptest_val=&test_val;
  
  for (int ii=0; ii<num_nodes; ii++){
    test_val=0.0;
    if (UniqueID()==ii){
      test_val=1.0;
    }
    glb_sum_five(ptest_val);
    while(test_val==0.0)
      dummy=0; //Does absolutely nothing, here as a placeholder for while
  }
  VRB.Result("","main()","Checkpoint no. %d reached.\n",chkpoint_no++);

  dtime_last=dtime_now;
  dtime_now=dclock();
  
  Float time_tmp=dtime_now-dtime_last;
  int hr_tmp=time_tmp/3600.0;
  int min_tmp=(time_tmp-3600.0*hr_tmp)/60.0;
  Float sec_tmp=time_tmp-3600.0*hr_tmp-60.0*min_tmp;
  VRB.Result("","main()","Time since last checkpoint: %d hours %d minutes %f seconds.\n",hr_tmp,min_tmp,sec_tmp);
  
  time_tmp=dtime_now-dtime_start;
  hr_tmp=time_tmp/3600.0;
  min_tmp=(time_tmp-3600.0*hr_tmp)/60.0;
  sec_tmp=time_tmp-3600.0*hr_tmp-60.0*min_tmp;
  VRB.Result("","main()","Time since beginning: %d hours %d minutes %f seconds.\n",hr_tmp, min_tmp,sec_tmp);
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
    VRB.Result("","main()","Error:: direction should be 0,1,2,3\n");
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
  
  VRB.Result("","main()","dir == %d \n",dir);
  VRB.Result("","main()","slice index == %d %d %d\n",slice_ind[0],slice_ind[1],slice_ind[2]);
  
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
