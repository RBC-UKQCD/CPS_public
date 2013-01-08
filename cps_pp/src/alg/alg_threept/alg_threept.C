#include <config.h>
#include <time.h>
//------------------------------------------------------------------
//
// alg_threept.C
//
// AlgThreePt is derived from Alg and is relevant to  
// three point correlation functions with Wilson-type 
// fermions (i.e. Wilson, Clover, Dwf). 
// The type of fermion is determined by the argument to the 
// constructor but it should not be a non-Wilson-type lattice.
//
// PLEASE NOTE that the propagators passed to this method will
// be altered to change propagators with periodic (P) and 
// antiperiodic (A) boundary conditions into P+A and P-A if
// necessary.  This is done so that extra memory isn't taken
// up storing separate P+A and P-A propagators.
//------------------------------------------------------------------

#include <util/qcdio.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/qioarg.h>
#include <util/time_cps.h>
#include <comms/glb.h>
#include <util/checksum.h>
#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
#include <qcdocos/scu_checksum.h>
#endif
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE
#endif
CPS_END_NAMESPACE
//#include <alg/alg_threept.h>
#include <alg/alg_threept.h>
#include <alg/wilson_matrix.h>
#include <alg/spin_matrix.h>
#include <alg/qpropw.h>
CPS_START_NAMESPACE

// Sums a quantity over the whole lattice; 99 is a magic number.
void lat_sum(Float *float_p, int block) {slice_sum(float_p,block,99);}

// Assigns a pointer of type TYPE to PTR with enough memory allocated
// for SIZE objects of that type. Prints error on fail.
#define SMALLOC(PTR,TYPE,SIZE) PTR = (TYPE*)smalloc(SIZE*sizeof(TYPE)); \
  if (PTR == 0) ERR.Pointer(cname,fname,#PTR); \
  VRB.Smalloc(cname,fname,#PTR,PTR,SIZE);

// Checkpoint function for timing.
void chkpt(const int num_nodes,int& chkpoint_no,const Float dtime_start,Float& dtime_last,Float& dtime_now);

#define FFLUSH(A)
//------------------------------------------------------------------
/*!
  \param latt The lattice on which to compute the three-point function.
  \param c_arg The common argument structure for all algorithms.
  \param arg The parameters specific to this algorithm.
 */
//------------------------------------------------------------------
AlgThreePt::AlgThreePt(Lattice& latt, CommonArg *c_arg, ThreePtArg *arg, ThreePtPropArg *p_arg): 
  Alg(latt,c_arg) {

  cname = "AlgThreePt";
  char *fname = "AlgThreePt(L&,CommonArg*,ThreePtArg*,ThreePtPropArg*)";
  VRB.Func(cname,fname);

  if (arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_threept_arg = arg;
  if (p_arg == 0)
    ERR.Pointer(cname,fname, "p_arg");
  alg_threept_prop_arg = p_arg;

  VRB.Input(cname,fname,
	    "stop_rsd = %g\n",IFloat (alg_threept_arg->cg.stop_rsd));
  VRB.Input(cname,fname,
	    "max_num_iter = %d\n", alg_threept_arg->cg.max_num_iter);
  VRB.Input(cname,fname,
	    "mass = %g\n",IFloat (alg_threept_arg->cg.mass));

  char tag[3][5] = {"SP","VA","TT"};
  for (int i=0;i<3;i++) strcpy(gam[i],tag[i]);
  char tag2[3][5] = {"SS","VV","TT"};
  for (int i=0;i<3;i++) strcpy(gam2[i],tag2[i]);
  char tag3[4][10] ={"TR_","TRTR_","TR_MX_","TRTR_MX_"};
  for (int i=0;i<4;i++) strcpy(tra[i],tag3[i]);

  if (arg->results != 0) {
	if ((fp = Fopen((char*)arg->results, "w")) == NULL) {
	  ERR.FileW(cname,fname,(char*)arg->results);
	}
  }

  Fprintf(fp,"Three-Point Code Version 4.9.20\nMASS= STRANGE LIGHT (PUPIL)\n");
  FFLUSH(fp);

  if (arg->results_mres_ZA != 0) {
    if ((fp_mres_ZA = Fopen((char*)arg->results_mres_ZA, "w")) == NULL) {
      ERR.FileW(cname,fname,(char*)arg->results_mres_ZA);
    }
  }

  Fprintf(fp_mres_ZA,"Three-Point Code Version 4.9.20\nMASS= STRANGE LIGHT (PUPIL)\n");
  FFLUSH(fp_mres_ZA);

  if (arg->results_pipi != 0) {
    if ((fp_pipi = Fopen((char*)arg->results_pipi, "w")) == NULL) {
      ERR.FileW(cname,fname,(char*)arg->results_pipi);
    }
  }

  Fprintf(fp_pipi,"Three-Point Code Version 4.9.20\nMASS= STRANGE LIGHT (PUPIL)\n");
  FFLUSH(fp_pipi);

  
}

// Destructor
//------------------------------------------------------------------
AlgThreePt::~AlgThreePt() {

  char *fname = "~AlgThreePt()";
  VRB.Func(cname,fname);

  Fclose(fp);
  Fclose(fp_mres_ZA);
  Fclose(fp_pipi);

}

// Solves for propagators and combines them into matrix elements
//------------------------------------------------------------------
void AlgThreePt::run() {

  //Initialize Timing
  Float dtime_start=dclock();
  Float dtime_last=dtime_start;
  Float dtime_now=dtime_start;
  
  int chkpoint_no=0;
  int chkpoints = alg_threept_arg->chkpoints;

  const int num_nodes=GJP.Xnodes()*GJP.Ynodes()*GJP.Znodes()*GJP.Tnodes()*GJP.Snodes();

  char *fname = "run()";
  VRB.Func(cname,fname);

  VRB.Flow(cname,fname,"Starting AlgThreePt::run()\n");
  //Checkpoint
  if (chkpoints)
    chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

  int num_light = alg_threept_arg->num_light;
  int num_strange = alg_threept_arg->num_strange;
  int num_heavy = alg_threept_arg->num_heavy;
  int num_tK = alg_threept_arg->num_tK;
  int num_hits = alg_threept_arg->num_hits;
  do_susy = alg_threept_arg->do_susy;
  Float *l_mass = alg_threept_arg->l_mass; // light masses
  Float *s_mass = alg_threept_arg->s_mass; // strange masses
  Float *h_mass = alg_threept_arg->h_mass; // heavy masses (pupils only)
  int *tK = alg_threept_arg->tK; // times to put kaon sources

  int do_zero_mom = alg_threept_arg->do_zero_mom;
  int do_first_mom = alg_threept_arg->do_first_mom;
  int do_second_mom = alg_threept_arg->do_second_mom;
  int do_third_mom = alg_threept_arg->do_third_mom;

  int do_p_plus_a_kaon = alg_threept_arg->do_p_plus_a_kaon;
  int do_kaon_at_walls = alg_threept_arg->do_kaon_at_walls;
  int do_kaons_tK = alg_threept_arg->do_kaons_tK;

  Lattice& lat = AlgLattice();
  QPropWArg qp_arg;
  qp_arg.cg = alg_threept_arg->cg;
  qp_arg.gauge_fix_src = TRUE;
  qp_arg.gauge_fix_snk = TRUE;
  qp_arg.save_prop = FALSE; // don't save propagators inside alg_threept
  qp_arg.x=0;
  qp_arg.y=0;
  qp_arg.z=0;
  qp_arg.store_midprop=0;
  qp_arg.save_ls_prop=0;
  qp_arg.do_half_fermion=0;
  qp_arg.ensemble_label=alg_threept_arg->ensemble_label;
  qp_arg.ensemble_id=alg_threept_arg->ensemble_id;
  qp_arg.seqNum=alg_threept_arg->seqNum;

  //Note that the figure8_spectator function (K->Pi Pi) is used in this version
  //of the code in a way that assumes t_src=0 and t_snk=time_size.
  int t_src = alg_threept_arg->t_src;   // source location
  int t_snk = alg_threept_arg->t_snk;   // sink location
  int t_shift = alg_threept_arg->t_shift;  // amount to shift lattice if generating propagators
  int t_op  = alg_threept_arg->t_op;    // operator (& pupil) location
  int width = alg_threept_arg->width;

  // set spatial boundary conditions periodic
  GJP.Xbc(BND_CND_PRD);
  GJP.Ybc(BND_CND_PRD);
  GJP.Zbc(BND_CND_PRD);

#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
  if(!ScuChecksum::ChecksumsOn())
	ScuChecksum::Initialise(true,true);
#endif

  //NOTE:  Propagators passed to this function will now be altered
  //       to change P and A propagators into P+A and P-A so that
  //       extra memory isn't need to store the P+A and P-A propagators.

  QPropW* q_light_tpi[num_light][2][4][3];
  QPropW* q_light_tK[num_light][2][num_tK];
  QPropW* q_strange_tpi[num_strange][2];
  QPropW* q_strange_tK[num_strange][2][num_tK];
  
  //Make these pointers point to the propagators passed via threept_prop_arg.
  for (int i=0; i<num_light; i++) {
    for (int bc=0; bc<2; bc++) {
      for (int mom_num=0; mom_num<4; mom_num++) {
	
	if ( !do_first_mom && mom_num==1 ) continue;
	if ( !do_second_mom && mom_num==2 ) continue;
	if ( !do_third_mom && mom_num==3 ) continue;
	
	int n_dir;
	if (mom_num==0 || mom_num==3)
	  n_dir=1;
	else
	  n_dir=3;
	
	for (int mom_dir=0; mom_dir<n_dir; mom_dir++) {
	  q_light_tpi[i][bc][mom_num][mom_dir]=alg_threept_prop_arg->q_light_tpi[i][bc][mom_num][mom_dir];
	  q_light_tpi[i][bc][mom_num][mom_dir]->ChangeSourceTime(t_src+(t_snk-t_src)*bc);
	}
	
      } //mom_num for loop

      if (do_kaons_tK)
	for (int j=0; j<num_tK; j++) {
	  q_light_tK[i][bc][j]=alg_threept_prop_arg->q_light_tK[i][bc][j];
	  q_light_tK[i][bc][j]->ChangeSourceTime(tK[j]);
	}

    } //bc for loop
  } //i for loop
  
  for (int i=0; i<num_strange; i++) {
    if (do_kaon_at_walls) {
      for (int bc=0; bc<2; bc++) {
	q_strange_tpi[i][bc]=alg_threept_prop_arg->q_strange_tpi[i][bc];
	q_strange_tpi[i][bc]->ChangeSourceTime(t_src+(t_snk-t_src)*bc);
      }
    }
    int bc_num=1+do_p_plus_a_kaon;
    for (int bc=0; bc<bc_num; bc++) {
      for (int j=0; j<num_tK; j++) {
	q_strange_tK[i][bc][j]=alg_threept_prop_arg->q_strange_tK[i][bc][j];
	q_strange_tK[i][bc][j]->ChangeSourceTime(tK[j]);
      }
    }
  }
  
  //Do P+A and P-A

  //Light quarks
  for (int m=0; m<num_light; m++) {
    
    //tpi light quarks

    //Calculate the contractions needed in the denominators
    //for mres and ZA before we do P+A.
    Fprintf(fp_mres_ZA, "MASS= %e %e\n", l_mass[m],l_mass[m]);
    spectrum(*q_light_tpi[m][0][0][0],*q_light_tpi[m][0][0][0],0);
    spectrum(*q_light_tpi[m][1][0][0],*q_light_tpi[m][1][0][0],1);
    box_spectrum(*q_light_tpi[m][0][0][0],*q_light_tpi[m][0][0][0],0);
    box_spectrum(*q_light_tpi[m][1][0][0],*q_light_tpi[m][1][0][0],1);
    
    
    //Change from P,A to (P+A)/2,(P-A)/2, (loosely 
    //called P+A and P-A above). Right now 
    //q_light_tpi[m][0][mom_num][mom_dir] is P
    //and q_light_tpi[m][1][mom_num][mom_dir] is A.
    
    q_light_tpi[m][0][0][0]->Average(*q_light_tpi[m][1][0][0]); //q_light_tpi[m][0][0][0] is now (P+A)/2
    q_light_tpi[m][1][0][0]->LinComb(*q_light_tpi[m][0][0][0],-1.0,1.0);
    //q_light_tpi[m][1][0][0] is now -A+(P+A)/2 = (P-A)/2
    
    //Now do it for the non-zero momentum
    
    for (int mom_num=1; mom_num<4; mom_num++) {
      
      if (!do_first_mom && mom_num==1) continue;
      if (!do_second_mom && mom_num==2) continue;
      if (!do_third_mom && mom_num==3) continue;
      
      int n_dir;
      if (mom_num==3)
	n_dir=1;
      else
	n_dir=3;

      for (int mom_dir=0; mom_dir<n_dir; mom_dir++) {
	q_light_tpi[m][0][mom_num][mom_dir]->Average(*q_light_tpi[m][1][mom_num][mom_dir]);
	//q_light_tpi[m][0][mom_num][mom_dir] is now (P+A)/2
	q_light_tpi[m][1][mom_num][mom_dir]->LinComb(*q_light_tpi[m][0][mom_num][mom_dir],-1.0,1.0);
	//q_light_tpi[m][1][mom_num][mom_dir] is now -A+(P+A)/2 = (P-A)/2
      } //mom_dir for loop
      
    } //mom_num for loop


    //tK light quarks
    
    if (do_kaons_tK) {
      for (int j=0; j<num_tK; j++) {
	q_light_tK[m][0][j]->Average(*q_light_tK[m][1][j]);  //q_light_tK[m][0][j] is now (P+A)/2
      }
    }
    
  } //m (for loop, light quarks)

  //Light quarks
  for (int m=0; m<num_light; m++) {
  for (int n=0; n<= m; n++) {
    
    //tpi light quarks

    QPropW *q_src_m=q_light_tpi[m][0][0][0];
    QPropW *q_src_n=q_light_tpi[n][0][0][0];
    QPropW *q_snk_m=q_light_tpi[m][1][0][0];
    QPropW *q_snk_n=q_light_tpi[n][1][0][0];

    //Calculate the contractions needed in the denominators
    //for mres and ZA before we do P+A.
    Fprintf(fp, "MASS= %e %e\n", l_mass[n],l_mass[m]);
    spectrum(*q_src_n,*q_src_m,2);
    spectrum(*q_snk_n,*q_snk_m,2);
    box_spectrum(*q_src_n,*q_src_m,2);
    box_spectrum(*q_snk_n,*q_snk_m,2);

    figure8(*q_src_n,*q_src_m,*q_snk_n,*q_snk_m);
    figure8(*q_src_n,*q_src_m,*q_snk_m,*q_snk_m,1);
    figure8(*q_snk_n,*q_snk_m,*q_src_m,*q_src_m,1);
    figure8_vacuum(*q_src_n,*q_src_m,*q_snk_m,*q_snk_m,*q_snk_m);
    figure8_vacuum(*q_snk_n,*q_snk_m,*q_src_m,*q_src_m,*q_src_m);
    figure8_spectator_old(*q_src_n,*q_src_m,*q_snk_m,*q_snk_m,*q_snk_m);
    figure8_spectator_old(*q_snk_n,*q_snk_m,*q_src_m,*q_src_m,*q_src_m);
//    fflush(fp);

  }
  }

  
  //Strange quarks
  for (int m=0; m<num_strange; m++) {
    
    //tpi strange quarks
    
    if (do_kaon_at_walls) {
      //Calculate the contractions needed in the denominators
      //for mres and ZA before we do P+A.
      Fprintf(fp_mres_ZA, "MASS= %e %e\n",s_mass[m],s_mass[m]);
      spectrum(*q_strange_tpi[m][0],*q_strange_tpi[m][0],0);
      spectrum(*q_strange_tpi[m][1],*q_strange_tpi[m][1],1);
      box_spectrum(*q_strange_tpi[m][0],*q_strange_tpi[m][0],0);
      box_spectrum(*q_strange_tpi[m][1],*q_strange_tpi[m][1],1);
      
      
      //Change from P,A to (P+A)/2,(P-A)/2, (loosely 
      //called P+A and P-A above). Right now 
      //q_strange_tpi[m][0] is P and q_strange_tpi[m][1] is A.
      
      q_strange_tpi[m][0]->Average(*q_strange_tpi[m][1]); //q_strange_tpi[m][0] is now (P+A)/2
      q_strange_tpi[m][1]->LinComb(*q_strange_tpi[m][0],-1.0,1.0);
      //q_strange_tpi[m][1] is now -A+(P+A)/2 = (P-A)/2
    }

    //tK strange quarks

    if (do_p_plus_a_kaon) {
      for (int j=0; j<num_tK; j++) {
	q_strange_tK[m][0][j]->Average(*q_strange_tK[m][1][j]); //q_strange_tK[m][0][j] is now (P+A)/2
      }
    }
    
  } //m (for loop, strange quarks)

  FFLUSH(fp_mres_ZA);
  
  
  /*----------------------------------------------------------------------------------
  At this point q_light_tpi[m][0][mom_num][mom_dir] and q_strange_tpi[m][0] are 
  "source" propagators (their sources are at t_src) and 
  q_light_tpi[m][1][mom_num][mom_dir] and q_strange_tpi[m][1] are "sink" propagators 
  (their sources are at t_snk).
  The q_strange_tK[m][0][src_num] have their sources at tK[src_num] and are P+A if
  the do_p_plus_a_kaon flag is turned on and are just P if this flag is turned off.
  Either way, DON'T use the q_strange_tK[m][1][src_num] anymore.
  The q_light_tK[m][0][src_num] have their sources at tK[src_num] and are P+A.
  DON'T use the q_light_tK[m][1][src_num] anymore.
  ----------------------------------------------------------------------------------*/


  //Do Pi Pi correlators for all masses first.
  for (int m=0; m<num_light; m++) {
    
    Fprintf(fp_pipi, "MASS= %e %e\n",l_mass[m],l_mass[m]);
    if (do_zero_mom) {
      pipi(*q_light_tpi[m][0][0][0]);
      pipi(*q_light_tpi[m][1][0][0]);
      
      VRB.Flow(cname,fname,"Did Pi Pi correlator 0 momentum, mass=%f.\n",l_mass[m]);
      //Checkpoint
      if (chkpoints)
	chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
      
    }
    
    for (int mom_num=1; mom_num<4; mom_num++) {
      
      if (!do_first_mom && mom_num==1) continue;
      if (!do_second_mom && mom_num==2) continue;
      if (!do_third_mom && mom_num==3) continue;

      int n_dir;
      if (mom_num==3)
	n_dir=1;
      else
	n_dir=3;
      for (int mom_dir=0; mom_dir<n_dir; mom_dir++) {
	pipi(*q_light_tpi[m][0][0][0],*q_light_tpi[m][0][mom_num][mom_dir],mom_num,mom_dir);
	pipi(*q_light_tpi[m][1][0][0],*q_light_tpi[m][1][mom_num][mom_dir],mom_num,mom_dir);
      } //mom_dir for loop

      VRB.Flow(cname,fname,"Did Pi Pi correlator momnum=%d, mass=%f.\n",mom_num,l_mass[m]);
      //Checkpoint
      if (chkpoints)
	chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

    } //mom_num for loop

    FFLUSH(fp_pipi);

  } //m for loop (light quarks)

  //Pi Pi correlators with heavy pions (strange quarks)
  for (int m=0; m<num_strange; m++) {
    
    if (do_zero_mom) {

      Fprintf(fp_pipi, "MASS= %e %e\n",s_mass[m],s_mass[m]);

      if (do_kaon_at_walls) {
	pipi(*q_strange_tpi[m][0]);
	pipi(*q_strange_tpi[m][1]);
      }

      if (do_p_plus_a_kaon) {
	for (int j=0; j<num_tK; j++) {
	  pipi(*q_strange_tK[m][0][j]);
	}
      }

      VRB.Flow(cname,fname,"Did (heavy) Pi Pi correlator 0 momentum, mass=%f.\n",s_mass[m]);
      //Checkpoint
      if (chkpoints)
	chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

    } //do_zero_mom if statement


    FFLUSH(fp_pipi);

  } //m for loop (strange quarks)



  // Loop over strange AND light masses
  for (int m_l=0; m_l<num_light; m_l++) {

    //----------------------------------------------------------------
    //old version n= strange mass, m= light mass
	for (int m_s=0; m_s<num_strange; m_s++) {

	  //Quantities for degenerate masses - Two light quarks (e.g. pion correlator):

	  if (m_s==0) { //so that the same thing isn't repeated num_strange times

	    Fprintf(fp, "MASS= %e %e\n",l_mass[m_l],l_mass[m_l]);
	    
	    //Zero momentum
	    if (do_zero_mom) {
	      spectrum(*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][0][0]);
	      spectrum(*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][0][0]);
	      box_spectrum(*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][0][0]);
	      box_spectrum(*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][0][0]);

	      if (do_kaons_tK) {
		for (int j=0; j<num_tK; j++) {
		  spectrum(*q_light_tK[m_l][0][j],*q_light_tK[m_l][0][j]);
		  box_spectrum(*q_light_tK[m_l][0][j],*q_light_tK[m_l][0][j]);
		}
	      }

	      VRB.Flow(cname,fname,"Did 0 momentum pion correlators, mass=%f.\n",l_mass[m_l]);
	      //Checkpoint
	      if (chkpoints)
		chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

	    }

	    //Non-zero momentum
	    for (int mom_num=1; mom_num<4; mom_num++) {
	      
	      if (!do_first_mom && mom_num==1) continue;
	      if (!do_second_mom && mom_num==2) continue;
	      if (!do_third_mom && mom_num==3) continue;
	      
	      int n_dir;
	      if (mom_num==3)
		n_dir=1;
	      else
		n_dir=3;

	      for (int mom_dir=0; mom_dir<n_dir; mom_dir++) {
		spectrum(*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][mom_num][mom_dir],mom_num,mom_dir);
		spectrum(*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][mom_num][mom_dir],mom_num,mom_dir);
		box_spectrum(*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][mom_num][mom_dir],mom_num,mom_dir);
		box_spectrum(*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][mom_num][mom_dir],mom_num,mom_dir);
	      } //mom_dir for loop

	      VRB.Flow(cname,fname,"Did pion correlators, momnum=%d, mass=%f.\n",mom_num,l_mass[m_l]);
	      //Checkpoint
	      if (chkpoints)
		chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
	      
	    } //mom_num for loop

	  } //m_s if statement

	  //Quantities for degenerate masses - Two strange quarks (e.g. heavy pion correlator):
	  
	  if (do_zero_mom && m_l==0) { //so that the same thing isn't repeated num_light times
	    Fprintf(fp, "MASS= %e %e\n",s_mass[m_s],s_mass[m_s]);
	    if (do_kaon_at_walls) {
	      spectrum(*q_strange_tpi[m_s][0],*q_strange_tpi[m_s][0]);
	      spectrum(*q_strange_tpi[m_s][1],*q_strange_tpi[m_s][1]);
	      box_spectrum(*q_strange_tpi[m_s][0],*q_strange_tpi[m_s][0]);
	      box_spectrum(*q_strange_tpi[m_s][1],*q_strange_tpi[m_s][1]);
	    }
	    if (do_p_plus_a_kaon) {
	      for (int j=0; j<num_tK; j++) {
		spectrum(*q_strange_tK[m_s][0][j],*q_strange_tK[m_s][0][j]);
		box_spectrum(*q_strange_tK[m_s][0][j],*q_strange_tK[m_s][0][j]);
	      }
	    }

	    VRB.Flow(cname,fname,"Did heavy pion correlators mass=%f.\n",s_mass[m_s]);
	    //Checkpoint
	    if (chkpoints)
	      chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
	    
	  }
	  
	  
	  //Quantities for non-degenerate masses.
	  
	  Fprintf(fp, "MASS= %e %e\n",s_mass[m_s],l_mass[m_l]);

	  //Zero Momentum
	  if (do_zero_mom) {

	    if (do_kaon_at_walls) {
	      spectrum(*q_strange_tpi[m_s][0],*q_light_tpi[m_l][0][0][0]);
	      spectrum(*q_strange_tpi[m_s][1],*q_light_tpi[m_l][1][0][0]);
	      box_spectrum(*q_strange_tpi[m_s][0],*q_light_tpi[m_l][0][0][0]);
	      box_spectrum(*q_strange_tpi[m_s][1],*q_light_tpi[m_l][1][0][0]);

	      VRB.Flow(cname,fname,"Did kaon correlators at walls, m_s=%f, m_l=%f.\n",s_mass[m_s],l_mass[m_l]);
	      //Checkpoint
	      if (chkpoints)
		chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);


	      
	      figure8(*q_strange_tpi[m_s][0],*q_light_tpi[m_l][0][0][0],*q_strange_tpi[m_s][1],*q_light_tpi[m_l][1][0][0]);
	      VRB.Flow(cname,fname,"Did figure8 number 1, 0 momentum, m_s=%f, m_l=%f.\n",s_mass[m_s],l_mass[m_l]);
	      //Checkpoint
	      if (chkpoints)
		chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
	      //----------------------------------------------
	      figure8(*q_strange_tpi[m_s][0],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][0][0],1);
	      VRB.Flow(cname,fname,"Did figure8 number 2, 0 momentum, m_s=%f, m_l=%f.\n",s_mass[m_s],l_mass[m_l]);
	      //Checkpoint
	      if (chkpoints)
		chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
	      //----------------------------------------------
	      figure8(*q_strange_tpi[m_s][1],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][0][0],1);
	      VRB.Flow(cname,fname,"Did figure8 number 3, 0 momentum, m_s=%f, m_l=%f.\n",s_mass[m_s],l_mass[m_l]);
	      //Checkpoint
	      if (chkpoints)
		chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
	      //----------------------------------------------
	      figure8_vacuum(*q_strange_tpi[m_s][0],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][0][0]);
	      VRB.Flow(cname,fname,"Did figure8_vacuum number 1, 0 momentum, m_s=%f, m_l=%f.\n",s_mass[m_s],l_mass[m_l]);
	      //Checkpoint
	      if (chkpoints)
		chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
	      //----------------------------------------------
	      figure8_vacuum(*q_strange_tpi[m_s][1],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][0][0]);
	      VRB.Flow(cname,fname,"Did figure8_vacuum number 2, 0 momentum, m_s=%f, m_l=%f.\n",s_mass[m_s],l_mass[m_l]);
	      //Checkpoint
	      if (chkpoints)
		chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);



	      
	      //K->PiPi using the figure8_spectator function.  
	      //NOTE that the figure8_spectator function has been CHANGED from the original so that it works
	      //for the specific way it is used in this version of the code (see note inside the function).
	      
	      //Kaon and pion at very left (t=t_src) and very right (t=t_snk) walls.
	      //----------------------------------------------
	      figure8_spectator(*q_strange_tpi[m_s][0],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][0][0]);
	      VRB.Flow(cname,fname,"Did figure8_spectator number 1, 0 momentum, m_s=%f, m_l=%f.\n",s_mass[m_s],l_mass[m_l]);
	      //Checkpoint
	      if (chkpoints)
		chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
	      //----------------------------------------------
	      figure8_spectator(*q_strange_tpi[m_s][1],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][0][0]);
	      VRB.Flow(cname,fname,"Did figure8_spectator number 2, 0 momentum, m_s=%f, m_l=%f.\n",s_mass[m_s],l_mass[m_l]);
	      //Checkpoint
	      if (chkpoints)
		chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);


	    } //do_kaon_at_walls if statement

	    //Kaon correlator at tK[j]
	    if (do_kaons_tK) {
	      for (int j=0; j<num_tK; j++) {
		spectrum(*q_strange_tK[m_s][0][j],*q_light_tK[m_l][0][j]);
		box_spectrum(*q_strange_tK[m_s][0][j],*q_light_tK[m_l][0][j]);
	      }
	    }
	    
	    VRB.Flow(cname,fname,"Did kaon correlators at times tK[j], m_s=%f, m_l=%f.\n",s_mass[m_s],l_mass[m_l]);
	    //Checkpoint
	    if (chkpoints)
	      chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
	    
	    //Kaon at tK[j] and pion at left and right (to get both a time separation of
	    //tK[j] and Nt-tK[j] between the kaon and the two pions)
	    for (int j=0; j<num_tK; j++) {
	      
	      //Do contractions with pions at BOTH walls
	      figure8_spectator(*q_strange_tK[m_s][0][j],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][0][0]);
	      figure8_spectator(*q_strange_tK[m_s][0][j],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][0][0]);
	      
	      VRB.Flow(cname,fname,"Did figure8_spectator with 0 momentum, tK=%d, m_s=%f, m_l=%f.\n",tK[j],s_mass[m_s],l_mass[m_l]);
	      //Checkpoint
	      if (chkpoints)
		chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
	      
	    } //j for loop
	    
	  } //do_zero_mom if statement
	  

	  //Non-zero Momentum
	  for (int mom_num=1; mom_num<4; mom_num++) {
	    
	    if (!do_first_mom && mom_num==1) continue;
	    if (!do_second_mom && mom_num==2) continue;
	    if (!do_third_mom && mom_num==3) continue;

	    int n_dir;
	    if (mom_num==3)
	      n_dir=1;
	    else
	      n_dir=3;
	    for (int mom_dir=0; mom_dir<n_dir; mom_dir++) {

	      //Kaon with non-zero momentum

	      if (do_kaon_at_walls) {
		spectrum(*q_strange_tpi[m_s][0],*q_light_tpi[m_l][0][mom_num][mom_dir],mom_num,mom_dir);
		spectrum(*q_strange_tpi[m_s][1],*q_light_tpi[m_l][1][mom_num][mom_dir],mom_num,mom_dir);
		box_spectrum(*q_strange_tpi[m_s][0],*q_light_tpi[m_l][0][mom_num][mom_dir],mom_num,mom_dir);
		box_spectrum(*q_strange_tpi[m_s][1],*q_light_tpi[m_l][1][mom_num][mom_dir],mom_num,mom_dir);
	      }

	      VRB.Flow(cname,fname,"Did kaon correlators at walls, momnum=%d, momdir=%d, m_s=%f, m_l=%f.\n",mom_num,mom_dir,s_mass[m_s],l_mass[m_l]);
	      //Checkpoint
	      if (chkpoints)
		chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);

	      //K->PiPi using figure8_spectator for non-zero momentum pions
	      //q_snk1 and q_snk3 are d, q_snk2 and q_spc are u
	      
	      if (do_kaon_at_walls) {
		//Kaon and pion at very left (t=t_src) and very right (t=t_snk) walls.
		figure8_spectator(*q_strange_tpi[m_s][0],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][mom_num][mom_dir],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][mom_num][mom_dir],mom_num,mom_dir);
		figure8_spectator(*q_strange_tpi[m_s][1],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][0][mom_num][mom_dir],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][mom_num][mom_dir],mom_num,mom_dir);
	      } //do_kaon_at_walls if statement

	      VRB.Flow(cname,fname,"Did figure8_spectator at walls, momnum=%d, momdir=%d, m_s=%f, m_l=%f.\n",mom_num,mom_dir,s_mass[m_s],l_mass[m_l]);
	      //Checkpoint
	      if (chkpoints)
		chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
	      
	      //Kaon at tK[j] and pion at left and right
	      for (int j=0; j<num_tK; j++) {
		
		//Do contractions with pions at BOTH walls
		figure8_spectator(*q_strange_tK[m_s][0][j],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][mom_num][mom_dir],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][mom_num][mom_dir],mom_num,mom_dir);
		figure8_spectator(*q_strange_tK[m_s][0][j],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][mom_num][mom_dir],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][mom_num][mom_dir],mom_num,mom_dir);
		
		VRB.Flow(cname,fname,"Did figure8_spectator with momnum=%d, momdir=%d, tK=%d, m_s=%f, m_l=%f.\n",mom_num,mom_dir,tK[j],s_mass[m_s],l_mass[m_l]);
		//Checkpoint
		if (chkpoints)
		  chkpt(num_nodes,chkpoint_no,dtime_start,dtime_last,dtime_now);
		
	      } //j for loop
	      
	    } //mom_dir for loop

	  } //mom_num for loop
	  
	  FFLUSH(fp);
	} // for strange masses
   //----------------------------------------------------------------

#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
	printf("SCU checksum test\n");
	if (!ScuChecksum::CsumSwap())
	  ERR.Hardware(cname,fname, "SCU Checksum mismatch\n");
#endif

  } // for light masses

  //qp_arg.type = SLAB;
  //qp_arg.rng = GAUSS;
  //qp_arg.slab_width = width;
  qp_arg.t = t_op;
  GJP.Tbc(BND_CND_PRD);
  qp_arg.save_prop = FALSE; // do NOT save pupils
  if (do_kaon_at_walls) {
    for (int hit=0; hit<num_hits; hit++) {
      ERR.NotImplemented("alg_threept","run");
      for (int n=0; n<num_light+num_strange+num_heavy; n++) {
	qp_arg.cg.mass = (n<num_light?l_mass[n]:(n<num_light+num_strange?s_mass[n-num_light]:h_mass[n-num_light-num_strange]));
	QPropW *q_ppl = new QPropW(lat, &qp_arg, common_arg);
	q_ppl->Run(); // manually activate inversion
	
	//----------------------------------------------------------------
	for (int m_l=0; m_l<num_light; m_l++) {
	  for (int m_s=0; m_s<num_strange; m_s++) {
	    Fprintf(fp, "MASS= %e %e %e\n",s_mass[m_s],l_mass[m_l],qp_arg.cg.mass);
	    k_to_vac(*q_strange_tpi[m_s][0],*q_light_tpi[m_l][0][0][0],*q_ppl);
	    k_to_vac(*q_strange_tpi[m_s][1],*q_light_tpi[m_l][1][0][0],*q_ppl);
	    eye(*q_strange_tpi[m_s][0],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][1][0][0],*q_ppl);
	    eye(*q_strange_tpi[m_s][1],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][0][0][0],*q_ppl);
	    eye_vacuum(*q_strange_tpi[m_s][0],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][1][0][0],*q_ppl);
	    eye_vacuum(*q_strange_tpi[m_s][1],*q_light_tpi[m_l][1][0][0],*q_light_tpi[m_l][0][0][0],*q_light_tpi[m_l][0][0][0],*q_ppl);
	    FFLUSH(fp);
	  } // for strange masses
	} // for light masses
	//----------------------------------------------------------------
	delete q_ppl;
      } // for pupil masses
#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
      printf("SCU checksum test\n");
      if (!ScuChecksum::CsumSwap())
	ERR.Hardware(cname,fname, "SCU Checksum mismatch\n");
#endif
    } // for hits
  } //do_kaon_at_walls if statement
  
  
}

// Pi Pi correlator
//------------------------------------------------------------------
void AlgThreePt::pipi(QPropW& q_src) {
  
  char *fname = "pipi()";
  VRB.Func(cname,fname);
  
  int t_src = q_src.SourceTime();
  int t, time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex *tr, *trtr;
  SMALLOC(tr,Rcomplex,time_size); //TR piece of Pi Pi correlator
  SMALLOC(trtr,Rcomplex,time_size); //TRTR piece of Pi Pi correlator

  for (t=0; t<time_size; t++)
    tr[t] = trtr[t] = 0.0;

  for (t=0; t<time_size; t++) {
    WilsonMatrix tmp_q = q_src.WallSinkProp(t);
    WilsonMatrix tmp_qh = tmp_q;
    tmp_qh.hconj();
    tr[t] += Trace(tmp_q * tmp_qh, tmp_q * tmp_qh);
    trtr[t] += Trace(tmp_q, tmp_qh) * Trace(tmp_q, tmp_qh);
  }
  // Print out results
  //----------------------------------------------------------------
  for (t=0; t<time_size; t++)
    Fprintf(fp_pipi, "PiPi_TR %d %d  %.16e %.16e\n", t_src, t, tr[t].real(), tr[t].imag());
  for (t=0; t<time_size; t++)
    Fprintf(fp_pipi, "PiPi_TRTR %d %d  %.16e %.16e\n", t_src, t, trtr[t].real(), trtr[t].imag());

  sfree(tr);
  sfree(trtr);

}

// Pi Pi correlator for non-zero momentum
//------------------------------------------------------------------
void AlgThreePt::pipi(QPropW& q_u, QPropW& q_d_mom, int mom_num, int mom_dir) {
  
  char *fname = "pipi()";
  VRB.Func(cname,fname);
  
  int t_src = q_u.SourceTime();
  int t, time_size=GJP.Tnodes()*GJP.TnodeSites();

  int n_comb=0; //Number of different possible momenta for the two pions in the final state
                //(doing both |p1,p2> and |p2,p1> even though they are the same by Bose symmetry).
                //n_comb=2^mom_num if not doing non-zero total momentum correlators, otherwise
                //n_comb=2^(2*mom_num)=(2^mom_num)^2
  if (mom_num==1)
    n_comb=2;
  else if (mom_num==2)
    n_comb=4;
  else if (mom_num==3)
    n_comb=8;
  else
    ERR.General(cname,fname,"Input parameter mom_num must be 1, 2, or 3, but received mom_num=%d\n",mom_num);

  if (mom_dir<0 || mom_dir>2)
    ERR.General(cname,fname,"Input parameter mom_dir must be 0, 1, or 2, but received mom_dir=%d\n",mom_dir);

  if (alg_threept_arg->do_pipi_non_zero_tot_mom)
    n_comb=n_comb*n_comb;
  
  Rcomplex *tr[n_comb], *trtr[n_comb];
  Rcomplex *tr_cos, *trtr_cos; //Cosine sink
  for (int i=0; i<n_comb; i++) {
    SMALLOC(tr[i],Rcomplex,time_size); //TR piece of Pi Pi correlator
    SMALLOC(trtr[i],Rcomplex,time_size); //TRTR piece of Pi Pi correlator
  }
  SMALLOC(tr_cos,Rcomplex,time_size);
  SMALLOC(trtr_cos,Rcomplex,time_size);
  

  for (t=0; t<time_size; t++) {
    for (int i=0; i<n_comb; i++)
      tr[i][t] = trtr[i][t] = 0.0;
    tr_cos[t] = trtr_cos[t] = 0.0;
  }

  char mom_string[n_comb][100]; //String describing each momentum combination
  char mom_info[100]; //String describing the momentum of the cosine source
  int p[3];
  for (int i=0; i<3; i++)
    p[i]=0;
  if (mom_num==1) {
    p[mom_dir]=1;
  } else if (mom_num==2) {
    for (int mom_dir_tmp=0; mom_dir_tmp<3; mom_dir_tmp++)
      if (mom_dir_tmp!=mom_dir)
	p[mom_dir_tmp]=1;
  } else {
    for (int i=0; i<3; i++)
      p[i]=1;
  }
  sprintf(mom_info,"p=(%d,%d,%d)",p[0],p[1],p[2]);


  for (t=0; t<time_size; t++) {

    WilsonMatrix tmp_qu = q_u.WallSinkProp(t);

    int comb_num=0;
    int q1[3], q2[3], p1[3], p2[3];
    //Loop over all possible values of q1 and q2 (each component can be -1 or 1).
    //If mom_num is 1 or 2 then only execute the loop when the components in 
    //non-twisted directions are equal to -1.  Copy q1, q2 into p1, p2, and set
    //all components in non-twisted directions to 0 in p1 and p2 so that p1 and
    //p2 are the actual momenta of the two pions.
    for (q1[0]=-1; q1[0]<=1; q1[0]+=2)
      for (q1[1]=-1; q1[1]<=1; q1[1]+=2)
	for (q1[2]=-1; q1[2]<=1; q1[2]+=2)
	  for (q2[0]=-1; q2[0]<=1; q2[0]+=2)
	    for (q2[1]=-1; q2[1]<=1; q2[1]+=2)
	      for (q2[2]=-1; q2[2]<=1; q2[2]+=2) {

		for (int i=0; i<3; i++) {
		  p1[i]=q1[i];
		  p2[i]=q2[i];
		}

		int do_this_iter=1;
		if (mom_num==1)
		  for (int mom_dir_tmp=0; mom_dir_tmp<3; mom_dir_tmp++) {
		    if (mom_dir_tmp!=mom_dir)
		      if (q1[mom_dir_tmp]!=-1 || q2[mom_dir_tmp]!=-1) {
			do_this_iter=0;
			continue;
		      } else {
			p1[mom_dir_tmp]=0;
			p2[mom_dir_tmp]=0;
		      }
		  }
		if (mom_num==2)
		  if (q1[mom_dir]!=-1 || q2[mom_dir]!=-1)
		    do_this_iter=0;
		  else {
		    p1[mom_dir]=0;
		    p2[mom_dir]=0;
		  }
		if (!do_this_iter) continue;

		//Don't do non-zero total momentum if the do_pipi_non_zero_tot_mom flag is turned off.
		int p_tot[3];
		for (int i=0; i<3; i++)
		  p_tot[i]=p1[i]+p2[i];
		if (!alg_threept_arg->do_pipi_non_zero_tot_mom && (p_tot[0]!=0 || p_tot[1]!=0 || p_tot[2]!=0)) continue;

		WilsonMatrix tmp_qd_p1_h = q_d_mom.TwistMomSinkProp(t,p1);
		tmp_qd_p1_h.hconj();
		WilsonMatrix tmp_qd_p2_h = q_d_mom.TwistMomSinkProp(t,p2);
		tmp_qd_p2_h.hconj();

		tr[comb_num][t] += Trace(tmp_qu * tmp_qd_p1_h, tmp_qu * tmp_qd_p2_h);
		trtr[comb_num][t] += Trace(tmp_qu, tmp_qd_p1_h) * Trace(tmp_qu, tmp_qd_p2_h);

		//String to describe this momentum combination in the output
		sprintf(mom_string[comb_num],"p1=(%d,%d,%d)_p2=(%d,%d,%d)",p1[0],p1[1],p1[2],p2[0],p2[1],p2[2]);

		comb_num++;

	      } //q2[2] for loop
    //end q1[i], q2[i] for loops

    if (comb_num!=n_comb)
      ERR.General(cname,fname,"Number of momentum combinations done incorrect for some reason.\nDid %d, should have done %d.\nmom_num=%d  mom_dir=%d\n",comb_num,n_comb,mom_num,mom_dir);

    //Cosine sink

    //Use p with all non-zero entries equal to +1 (+/-1 doesn't matter for cosine sink).
    //This p was already calculated above for the mom_info string.
    
    WilsonMatrix tmp_qd_cos_h = q_d_mom.TwistCosSinkProp(t,p);
    tmp_qd_cos_h.hconj();

    tr_cos[t] += Trace(tmp_qu * tmp_qd_cos_h, tmp_qu * tmp_qd_cos_h);
    trtr_cos[t] += Trace(tmp_qu, tmp_qd_cos_h) * Trace(tmp_qu, tmp_qd_cos_h);
    
  } // t for loop


  // Print out results
  //----------------------------------------------------------------
  for (int comb_num=0; comb_num<n_comb; comb_num++) {
    for (t=0; t<time_size; t++)
      Fprintf(fp_pipi,"PiPi_Source_Cos_%s_Sink_%s_TR %d %d %.16e %.16e\n", mom_info, mom_string[comb_num], t_src, t, tr[comb_num][t].real(), tr[comb_num][t].imag());
    for (t=0; t<time_size; t++)
      Fprintf(fp_pipi,"PiPi_Source_Cos_%s_Sink_%s_TRTR %d %d %.16e %.16e\n", mom_info, mom_string[comb_num], t_src, t, trtr[comb_num][t].real(), trtr[comb_num][t].imag());
  }
  for (t=0; t<time_size; t++)
    Fprintf(fp_pipi,"PiPi_Source_Cos_%s_Sink_Cos_%s_TR %d %d %.16e %.16e\n", mom_info, mom_info, t_src, t, tr_cos[t].real(), tr_cos[t].imag());
  for (t=0; t<time_size; t++)
    Fprintf(fp_pipi,"PiPi_Source_Cos_%s_Sink_Cos_%s_TRTR %d %d %.16e %.16e\n", mom_info, mom_info, t_src, t, trtr_cos[t].real(), trtr_cos[t].imag());

  for (int i=0; i<n_comb; i++) {
    sfree(tr[i]);
    sfree(trtr[i]);
  }
  sfree(tr_cos);
  sfree(trtr_cos);

}
  
// box-point two-point functions
//------------------------------------------------------------------
void AlgThreePt::spectrum(QPropW& q1_src, QPropW& q2_src, int bc) { 
  
  char *fname = "spectrum()";
  VRB.Func(cname,fname);
  
  /*-------------------------------------------------------------------
    Note that this gives correlators that are the complex conjugate of
    the conventions used in my notes, i.e. here the light quark
    propagator (q2_src) is Hermitian conjugated whereas in the
    convention of my notes the strange quark propagator (q1_src)
    should be Hermitian conjugated.
  -------------------------------------------------------------------*/

  //bc = 0 : periodic boundary conditions, mres/ZA contractions
  //bc = 1 : antiperiodic boundary conditions, mres/ZA contractions
  //bc = 2 : normal spectrum code (P+A already done)
  if ( bc<0 || bc>2)
    ERR.General(cname,fname,"bc must be 0, 1, or 2, received %d\n",bc);

  int t_src = q1_src.SourceTime();
  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex *ps, *sc, *vc, *ax, *tn, *a0, *v0, *psax;
  SMALLOC(ps,  Rcomplex,time_size); // pseudo
  SMALLOC(sc,  Rcomplex,time_size); // scalar
  SMALLOC(vc,  Rcomplex,time_size); // vector
  SMALLOC(ax,  Rcomplex,time_size); // axial
  SMALLOC(tn,  Rcomplex,time_size); // tensor
  SMALLOC(a0,  Rcomplex,time_size); // axial-0
  SMALLOC(v0,  Rcomplex,time_size); // vector-0
  SMALLOC(psax,Rcomplex,time_size); // pseudo-axial
  WilsonMatrix tmp_src;
 
  int t;
  for (t=0; t<time_size; t++)
	ps[t] = sc[t] = vc[t] = ax[t] = tn[t] = a0[t] = v0[t] = psax[t] = 0.0;

  int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  int vol = (GJP.VolNodeSites()/GJP.TnodeSites());
  for (int i=0; i<GJP.VolNodeSites(); i++) {
    t = i/vol + shift_t;
    tmp_src = q2_src[i];
    tmp_src.hconj();
    ps[t] += Trace(q1_src[i], tmp_src);
    tmp_src = q2_src[i];
    //tmp_src.hconj().gl(-5).gr(-5);
    tmp_src.hconj();
    tmp_src.gl(-5).gr(-5);
    sc[t] += Trace(q1_src[i], tmp_src);
    tmp_src = q2_src[i];
    //tmp_src.hconj().gl(-5).gr(-5).gl(3).gr(3);
    tmp_src.hconj();
    tmp_src.gl(-5).gr(-5).gl(3).gr(3);
    v0[t] += Trace(q1_src[i], tmp_src);
    tmp_src = q2_src[i];
    //tmp_src.hconj().gl(-5).gr(-5).gl(-5).gl(3).gr(3).gr(-5);
    tmp_src.hconj();
    tmp_src.gl(-5).gr(-5).gl(-5).gl(3).gr(3).gr(-5);
    a0[t] += Trace(q1_src[i], tmp_src);
    tmp_src = q2_src[i];
    //tmp_src.hconj().gr(-5).gr(3).gr(-5);
    tmp_src.hconj();
    tmp_src.gr(-5).gr(3).gr(-5);
    psax[t] += Trace(q1_src[i], tmp_src);
  }
  for (int dir=0; dir<3; dir++) 
	for (int i=0; i<GJP.VolNodeSites(); i++) {
	  t = i/vol + shift_t;
	  tmp_src = q2_src[i];
	  //tmp_src.hconj().gl(-5).gr(-5).gl(dir).gr(dir);
	  tmp_src.hconj();
	  tmp_src.gl(-5).gr(-5).gl(dir).gr(dir);
	  vc[t] += Trace(q1_src[i], tmp_src);
	  tmp_src = q2_src[i];
	  //tmp_src.hconj().gl(-5).gr(-5).gl(-5).gl(dir).gr(dir).gr(-5);
	  tmp_src.hconj();
	  tmp_src.gl(-5).gr(-5).gl(-5).gl(dir).gr(dir).gr(-5);
	  ax[t] += Trace(q1_src[i], tmp_src);
	}
  for (int mu=0; mu<3; mu++) 
	for (int nu=mu+1; nu<3; nu++) 
	  for (int i=0; i<GJP.VolNodeSites(); i++) {
		t = i/vol + shift_t;
		tmp_src = q2_src[i];
		//tmp_src.hconj().gl(-5).gr(-5).gl(nu).gl(mu).gr(mu).gr(nu);
		tmp_src.hconj();
		tmp_src.gl(-5).gr(-5).gl(nu).gl(mu).gr(mu).gr(nu);
		tn[t] += Trace(q1_src[i], tmp_src);
	  }
  // Global sums
  for (t=0; t<time_size; t++) {
    lat_sum((Float*)&ps[t], 2);
    lat_sum((Float*)&sc[t], 2);
    lat_sum((Float*)&ax[t], 2);
    lat_sum((Float*)&vc[t], 2);
    lat_sum((Float*)&tn[t], 2);
    lat_sum((Float*)&a0[t], 2);
    lat_sum((Float*)&v0[t], 2);
    lat_sum((Float*)&psax[t], 2);
  }
  // Print out results
  //----------------------------------------------------------------

  FILE *fp_spec;
  if (bc==2)
    fp_spec=fp;
  else
    fp_spec=fp_mres_ZA;

  char src_string[10];
  if (bc==2) {
    sprintf(src_string,"%d",t_src);
  } else {
    if (bc==1 && t_src==alg_threept_arg->t_snk)
      sprintf(src_string,"A");
    else if (bc==0 && t_src==alg_threept_arg->t_src)
      sprintf(src_string,"P");
    else
      ERR.General(cname,fname,"Inconsistent bc and source time.\nt_src = %d , bc = %d\n",t_src,bc);
  }
  
  
  for (t=0; t<time_size; t++)
	Fprintf(fp_spec, "bPS_pPS %s %d  %.16e %.16e\n", src_string, t, ps[t].real(),  ps[t].imag());
  for (t=0; t<time_size; t++)
	Fprintf(fp_spec, "bSC_pSC %s %d  %.16e %.16e\n", src_string, t, sc[t].real(),  sc[t].imag());
  for (t=0; t<time_size; t++) 
	Fprintf(fp_spec, "bVC_pVC %s %d  %.16e %.16e\n", src_string, t, vc[t].real(),  vc[t].imag());
  for (t=0; t<time_size; t++)
	Fprintf(fp_spec, "bAX_pAX %s %d  %.16e %.16e\n", src_string, t, ax[t].real(),  ax[t].imag());
  for (t=0; t<time_size; t++)
	Fprintf(fp_spec, "bTN_pTN %s %d  %.16e %.16e\n", src_string, t, tn[t].real(),  tn[t].imag());
  for (t=0; t<time_size; t++)
	Fprintf(fp_spec, "bA0_pA0 %s %d  %.16e %.16e\n", src_string, t, a0[t].real(),  a0[t].imag());
  for (t=0; t<time_size; t++)
	Fprintf(fp_spec, "bV0_pV0 %s %d  %.16e %.16e\n", src_string, t, v0[t].real(),  v0[t].imag());
  for (t=0; t<time_size; t++)
	Fprintf(fp_spec, "bPS_pA0 %s %d  %.16e %.16e\n", src_string, t, psax[t].real(), psax[t].imag());

  sfree(psax);
  sfree(v0);
  sfree(a0);
  sfree(tn);
  sfree(vc);
  sfree(ax);
  sfree(ps);
  sfree(sc);

}

// box-point two-point functions for non-zero momentum
//------------------------------------------------------------------
void AlgThreePt::spectrum(QPropW& q_s, QPropW& q_u_mom, int mom_num, int mom_dir) { 
  
  char *fname = "spectrum()";
  VRB.Func(cname,fname);
  
  /*-------------------------------------------------------------------
    Note that this gives correlators that are the complex conjugate of
    the conventions used in my notes, i.e. here the light quark
    propagator (q_u_mom) is Hermitian conjugated whereas in the
    convention of my notes the strange quark propagator (q_s)
    should be Hermitian conjugated.
  -------------------------------------------------------------------*/

  int t_src = q_u_mom.SourceTime();
  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  
  int n_mom=0; //Number of different momenta given that each component can be plus or minus.
  if (mom_num==1)
    n_mom=2;
  else if (mom_num==2)
    n_mom=4;
  else if (mom_num==3)
    n_mom=8;
  else
    ERR.General(cname,fname,"Input parameter mom_num must be 1, 2, or 3, but received mom_num=%d\n",mom_num);

  if (mom_dir<0 || mom_dir>2)
    ERR.General(cname,fname,"Input parameter mom_dir must be 0, 1, or 2, but received mom_dir=%d\n",mom_dir);

  Rcomplex *ps[n_mom], *sc[n_mom], *vc[n_mom], *ax[n_mom], *tn[n_mom], *a0[n_mom], *v0[n_mom], *psax[n_mom];
  Rcomplex *ps_cos, *sc_cos, *vc_cos, *ax_cos, *tn_cos, *a0_cos, *v0_cos, *psax_cos; //Cosine sink
  for (int i=0; i<n_mom; i++) {
    SMALLOC(ps[i],  Rcomplex,time_size); // pseudo
    SMALLOC(sc[i],  Rcomplex,time_size); // scalar
    SMALLOC(vc[i],  Rcomplex,time_size); // vector
    SMALLOC(ax[i],  Rcomplex,time_size); // axial
    SMALLOC(tn[i],  Rcomplex,time_size); // tensor
    SMALLOC(a0[i],  Rcomplex,time_size); // axial-0
    SMALLOC(v0[i],  Rcomplex,time_size); // vector-0
    SMALLOC(psax[i],Rcomplex,time_size); // pseudo-axial
  }
  SMALLOC(ps_cos,  Rcomplex,time_size); // pseudo
  SMALLOC(sc_cos,  Rcomplex,time_size); // scalar
  SMALLOC(vc_cos,  Rcomplex,time_size); // vector
  SMALLOC(ax_cos,  Rcomplex,time_size); // axial
  SMALLOC(tn_cos,  Rcomplex,time_size); // tensor
  SMALLOC(a0_cos,  Rcomplex,time_size); // axial-0
  SMALLOC(v0_cos,  Rcomplex,time_size); // vector-0
  SMALLOC(psax_cos,Rcomplex,time_size); // pseudo-axial
  
  WilsonMatrix tmp_src;
 
  int t;
  for (t=0; t<time_size; t++) {
    for (int i=0; i<n_mom; i++)
	ps[i][t] = sc[i][t] = vc[i][t] = ax[i][t] = tn[i][t] = a0[i][t] = v0[i][t] = psax[i][t] = 0.0;
    ps_cos[t] = sc_cos[t] = vc_cos[t] = ax_cos[t] = tn_cos[t] = a0_cos[t] = v0_cos[t] = psax_cos[t] = 0.0;
  }
  
  char mom_string[n_mom][100]; //String describing each momentum at the sink
  char mom_info[100]; //String describing the momentum of the cosine source
  int p[3];
  for (int i=0; i<3; i++)
    p[i]=0;
  if (mom_num==1) {
    p[mom_dir]=1;
  } else if (mom_num==2) {
    for (int mom_dir_tmp=0; mom_dir_tmp<3; mom_dir_tmp++)
      if (mom_dir_tmp!=mom_dir)
	p[mom_dir_tmp]=1;
  } else {
    for (int i=0; i<3; i++)
      p[i]=1;
  }
  sprintf(mom_info,"p=(%d,%d,%d)",p[0],p[1],p[2]);

  Float Pi_const=3.141592654;

  int x, y, z;
  
  int shift_x = GJP.XnodeCoor()*GJP.XnodeSites();
  int shift_y = GJP.YnodeCoor()*GJP.YnodeSites();
  int shift_z = GJP.ZnodeCoor()*GJP.ZnodeSites();
  int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  //Local lattice dimensions:
  int size_x = GJP.XnodeSites();
  int size_y = GJP.YnodeSites();
  int size_z = GJP.ZnodeSites();
  int size_t = GJP.TnodeSites();
  int size_xy = size_x*size_y;
  int vol = (GJP.VolNodeSites()/GJP.TnodeSites()); // =size_x*size_y_size_z
  //Global lattice dimensions
  int Size_X = GJP.Xnodes()*GJP.XnodeSites();
  int Size_Y = GJP.Ynodes()*GJP.YnodeSites();
  int Size_Z = GJP.Znodes()*GJP.ZnodeSites();
  int Size_T = GJP.Tnodes()*GJP.TnodeSites();

  for (int i=0; i<GJP.VolNodeSites(); i++) {
    
    //Global coordinates
    t = i/vol + shift_t;
    z = (i%vol)/size_xy + shift_z;
    y = (i%size_xy)/size_x + shift_y;
    x = i%size_x + shift_x;
    
    int mom_ind=0;
    int q[3], p_snk[3];
    //Loop over all possible values of q (each component can be -1 or 1).
    //If mom_num is 1 or 2 then only execute the loop when the components in 
    //non-twisted directions are equal to -1.  Copy q into p_snk, and set
    //all components in non-twisted directions to 0 in p_snk so that it is
    //the actual momentum of the kaon.
    for (q[0]=-1; q[0]<=1; q[0]+=2)
      for (q[1]=-1; q[1]<=1; q[1]+=2)
	for (q[2]=-1; q[2]<=1; q[2]+=2) {
	  
	  for (int ii=0; ii<3; ii++)
	    p_snk[ii]=q[ii];
	  
	  int do_this_iter=1;
	  if (mom_num==1)
	    for (int mom_dir_tmp=0; mom_dir_tmp<3; mom_dir_tmp++) {
	      if (mom_dir_tmp!=mom_dir)
		if (q[mom_dir_tmp]!=-1) {
		  do_this_iter=0;
		  continue;
		} else {
		  p_snk[mom_dir_tmp]=0;
		}
	    }
	  if (mom_num==2)
	    if (q[mom_dir]!=-1)
	      do_this_iter=0;
	    else {
	      p_snk[mom_dir]=0;
	    }
	  if (!do_this_iter) continue;

	  //Calculate mom_fac=exp(i*p.x)
	  //(Note: not exp(-i*p.x) since the conventions are
	  //the complex conjugate of those in my notes, see
	  //note at beginning of this function).
	  Float pdotx=0.0;
	  pdotx+=((Float) p_snk[0]*x)/((Float) Size_X);
	  pdotx+=((Float) p_snk[1]*y)/((Float) Size_Y);
	  pdotx+=((Float) p_snk[2]*z)/((Float) Size_Z);
	  pdotx*=Pi_const;
	  Rcomplex mom_fac(cos(pdotx),sin(pdotx));
          
	  tmp_src = q_u_mom[i];
	  tmp_src.hconj();
	  ps[mom_ind][t] += mom_fac*Trace(q_s[i], tmp_src);
	  tmp_src = q_u_mom[i];
	  //tmp_src.hconj().gl(-5).gr(-5);
	  tmp_src.hconj();
	  tmp_src.gl(-5).gr(-5);
	  sc[mom_ind][t] += mom_fac*Trace(q_s[i], tmp_src);
	  tmp_src = q_u_mom[i];
	  //tmp_src.hconj().gl(-5).gr(-5).gl(3).gr(3);
	  tmp_src.hconj();
	  tmp_src.gl(-5).gr(-5).gl(3).gr(3);
	  v0[mom_ind][t] += mom_fac*Trace(q_s[i], tmp_src);
	  tmp_src = q_u_mom[i];
	  //tmp_src.hconj().gl(-5).gr(-5).gl(-5).gl(3).gr(3).gr(-5);
	  tmp_src.hconj();
	  tmp_src.gl(-5).gr(-5).gl(-5).gl(3).gr(3).gr(-5);
	  a0[mom_ind][t] += mom_fac*Trace(q_s[i], tmp_src);
	  tmp_src = q_u_mom[i];
	  //tmp_src.hconj().gr(-5).gr(3).gr(-5);
	  tmp_src.hconj();
	  tmp_src.gr(-5).gr(3).gr(-5);
	  psax[mom_ind][t] += mom_fac*Trace(q_s[i], tmp_src);

	  //String to describe this momentum combination in the output
	  sprintf(mom_string[mom_ind],"p=(%d,%d,%d)",p_snk[0],p_snk[1],p_snk[2]);

	  mom_ind++;
	  
	} //q[2] for loop
    //End of q[i] for loops
    
    if (mom_ind!=n_mom)
      ERR.General(cname,fname,"Number of momentum directions done incorrect for some reason.\nDid %d, should have done %d.\nmom_num=%d  mom_dir=%d\n",mom_ind,n_mom,mom_num,mom_dir);
    
    //Cosine (point) sink
    
    Float cos_fac=1.0;
    cos_fac*=cos(Pi_const*((Float) p[0]*x)/((Float) Size_X));
    cos_fac*=cos(Pi_const*((Float) p[1]*y)/((Float) Size_Y));
    cos_fac*=cos(Pi_const*((Float) p[2]*z)/((Float) Size_Z));

    tmp_src = q_u_mom[i];
    tmp_src.hconj();
    ps_cos[t] += cos_fac*Trace(q_s[i], tmp_src);
    tmp_src = q_u_mom[i];
    //tmp_src.hconj().gl(-5).gr(-5);
    tmp_src.hconj();
    tmp_src.gl(-5).gr(-5);
    sc_cos[t] += cos_fac*Trace(q_s[i], tmp_src);
    tmp_src = q_u_mom[i];
    //tmp_src.hconj().gl(-5).gr(-5).gl(3).gr(3);
    tmp_src.hconj();
    tmp_src.gl(-5).gr(-5).gl(3).gr(3);
    v0_cos[t] += cos_fac*Trace(q_s[i], tmp_src);
    tmp_src = q_u_mom[i];
    //tmp_src.hconj().gl(-5).gr(-5).gl(-5).gl(3).gr(3).gr(-5);
    tmp_src.hconj();
    tmp_src.gl(-5).gr(-5).gl(-5).gl(3).gr(3).gr(-5);
    a0_cos[t] += cos_fac*Trace(q_s[i], tmp_src);
    tmp_src = q_u_mom[i];
    //tmp_src.hconj().gr(-5).gr(3).gr(-5);
    tmp_src.hconj();
    tmp_src.gr(-5).gr(3).gr(-5);
    psax_cos[t] += cos_fac*Trace(q_s[i], tmp_src);    

  } //i for loop

  for (int dir=0; dir<3; dir++) 
	for (int i=0; i<GJP.VolNodeSites(); i++) {
    
	  //Global coordinates
	  t = i/vol + shift_t;
	  z = (i%vol)/size_xy + shift_z;
	  y = (i%size_xy)/size_x + shift_y;
	  x = i%size_x + shift_x;
	  
	  int mom_ind=0;
	  int q[3], p_snk[3];
	  //Loop over all possible values of q (each component can be -1 or 1).
	  //If mom_num is 1 or 2 then only execute the loop when the components in 
	  //non-twisted directions are equal to -1.  Copy q into p_snk, and set
	  //all components in non-twisted directions to 0 in p_snk so that it is
	  //the actual momentum of the kaon.
	  for (q[0]=-1; q[0]<=1; q[0]+=2)
	    for (q[1]=-1; q[1]<=1; q[1]+=2)
	      for (q[2]=-1; q[2]<=1; q[2]+=2) {
		
		for (int ii=0; ii<3; ii++)
		  p_snk[ii]=q[ii];
		
		int do_this_iter=1;
		if (mom_num==1)
		  for (int mom_dir_tmp=0; mom_dir_tmp<3; mom_dir_tmp++) {
		    if (mom_dir_tmp!=mom_dir)
		      if (q[mom_dir_tmp]!=-1) {
			do_this_iter=0;
			continue;
		      } else {
			p_snk[mom_dir_tmp]=0;
		      }
		  }
		if (mom_num==2)
		  if (q[mom_dir]!=-1)
		    do_this_iter=0;
		  else {
		    p_snk[mom_dir]=0;
		  }
		if (!do_this_iter) continue;
		
		//Calculate mom_fac=exp(i*p.x)
		//(Note: not exp(-i*p.x) since the conventions are
		//the complex conjugate of those in my notes, see
		//note at beginning of this function).
		Float pdotx=0.0;
		pdotx+=((Float) p_snk[0]*x)/((Float) Size_X);
		pdotx+=((Float) p_snk[1]*y)/((Float) Size_Y);
		pdotx+=((Float) p_snk[2]*z)/((Float) Size_Z);
		pdotx*=Pi_const;
		Rcomplex mom_fac(cos(pdotx),sin(pdotx));
		
		tmp_src = q_u_mom[i];
		//tmp_src.hconj().gl(-5).gr(-5).gl(dir).gr(dir);
		tmp_src.hconj();
		tmp_src.gl(-5).gr(-5).gl(dir).gr(dir);
		vc[mom_ind][t] += mom_fac*Trace(q_s[i], tmp_src);
		tmp_src = q_u_mom[i];
		//tmp_src.hconj().gl(-5).gr(-5).gl(-5).gl(dir).gr(dir).gr(-5);
		tmp_src.hconj();
		tmp_src.gl(-5).gr(-5).gl(-5).gl(dir).gr(dir).gr(-5);
		ax[mom_ind][t] += mom_fac*Trace(q_s[i], tmp_src);

		mom_ind++;
		
	      } //q[2] for loop
	  //End of q[i] for loops
	  
	  if (mom_ind!=n_mom)
	    ERR.General(cname,fname,"Number of momentum directions done incorrect for some reason.\nDid %d, should have done %d.\nmom_num=%d  mom_dir=%d\n",mom_ind,n_mom,mom_num,mom_dir);
	  
	  //Cosine (point) sink
	  
	  Float cos_fac=1.0;
	  cos_fac*=cos(Pi_const*((Float) p[0]*x)/((Float) Size_X));
	  cos_fac*=cos(Pi_const*((Float) p[1]*y)/((Float) Size_Y));
	  cos_fac*=cos(Pi_const*((Float) p[2]*z)/((Float) Size_Z));

	  tmp_src = q_u_mom[i];
	  //tmp_src.hconj().gl(-5).gr(-5).gl(dir).gr(dir);
	  tmp_src.hconj();
	  tmp_src.gl(-5).gr(-5).gl(dir).gr(dir);
	  vc_cos[t] += cos_fac*Trace(q_s[i], tmp_src);
	  tmp_src = q_u_mom[i];
	  //tmp_src.hconj().gl(-5).gr(-5).gl(-5).gl(dir).gr(dir).gr(-5);
	  tmp_src.hconj();
	  tmp_src.gl(-5).gr(-5).gl(-5).gl(dir).gr(dir).gr(-5);
	  ax_cos[t] += cos_fac*Trace(q_s[i], tmp_src);

	} //i for loop

  for (int mu=0; mu<3; mu++) 
	for (int nu=mu+1; nu<3; nu++) 
	  for (int i=0; i<GJP.VolNodeSites(); i++) {

	    //Global coordinates
	    t = i/vol + shift_t;
	    z = (i%vol)/size_xy + shift_z;
	    y = (i%size_xy)/size_x + shift_y;
	    x = i%size_x + shift_x;
	    
	    int mom_ind=0;
	    int q[3], p_snk[3];
	    //Loop over all possible values of q (each component can be -1 or 1).
	    //If mom_num is 1 or 2 then only execute the loop when the components in 
	    //non-twisted directions are equal to -1.  Copy q into p_snk, and set
	    //all components in non-twisted directions to 0 in p_snk so that it is
	    //the actual momentum of the kaon.
	    for (q[0]=-1; q[0]<=1; q[0]+=2)
	      for (q[1]=-1; q[1]<=1; q[1]+=2)
		for (q[2]=-1; q[2]<=1; q[2]+=2) {
		  
		  for (int ii=0; ii<3; ii++)
		    p_snk[ii]=q[ii];
		  
		  int do_this_iter=1;
		  if (mom_num==1)
		    for (int mom_dir_tmp=0; mom_dir_tmp<3; mom_dir_tmp++) {
		      if (mom_dir_tmp!=mom_dir)
			if (q[mom_dir_tmp]!=-1) {
			  do_this_iter=0;
			  continue;
			} else {
			  p_snk[mom_dir_tmp]=0;
			}
		    }
		  if (mom_num==2)
		    if (q[mom_dir]!=-1)
		      do_this_iter=0;
		    else {
		      p_snk[mom_dir]=0;
		    }
		  if (!do_this_iter) continue;
		  
		  //Calculate mom_fac=exp(i*p.x)
		  //(Note: not exp(-i*p.x) since the conventions are
		  //the complex conjugate of those in my notes, see
		  //note at beginning of this function).
		  Float pdotx=0.0;
		  pdotx+=((Float) p_snk[0]*x)/((Float) Size_X);
		  pdotx+=((Float) p_snk[1]*y)/((Float) Size_Y);
		  pdotx+=((Float) p_snk[2]*z)/((Float) Size_Z);
		  pdotx*=Pi_const;
		  Rcomplex mom_fac(cos(pdotx),sin(pdotx));

		  tmp_src = q_u_mom[i];
		  //tmp_src.hconj().gl(-5).gr(-5).gl(nu).gl(mu).gr(mu).gr(nu);
		  tmp_src.hconj();
		  tmp_src.gl(-5).gr(-5).gl(nu).gl(mu).gr(mu).gr(nu);
		  tn[mom_ind][t] += mom_fac*Trace(q_s[i], tmp_src);

		  mom_ind++;
		  
		} //q[2] for loop
	    //End of q[i] for loops
	    
	    if (mom_ind!=n_mom)
	      ERR.General(cname,fname,"Number of momentum directions done incorrect for some reason.\nDid %d, should have done %d.\nmom_num=%d  mom_dir=%d\n",mom_ind,n_mom,mom_num,mom_dir);
	    
	    //Cosine (point) sink
	    
	    Float cos_fac=1.0;
	    cos_fac*=cos(Pi_const*((Float) p[0]*x)/((Float) Size_X));
	    cos_fac*=cos(Pi_const*((Float) p[1]*y)/((Float) Size_Y));
	    cos_fac*=cos(Pi_const*((Float) p[2]*z)/((Float) Size_Z));
	    
	    tmp_src = q_u_mom[i];
	    //tmp_src.hconj().gl(-5).gr(-5).gl(nu).gl(mu).gr(mu).gr(nu);
	    tmp_src.hconj();
	    tmp_src.gl(-5).gr(-5).gl(nu).gl(mu).gr(mu).gr(nu);
	    tn_cos[t] += cos_fac*Trace(q_s[i], tmp_src);
	    
	  } //i for loop
  
  
  // Global sums
  for (t=0; t<time_size; t++) {
    for (int i=0; i<n_mom; i++) {
      lat_sum((Float*)&ps[i][t], 2);
      lat_sum((Float*)&sc[i][t], 2);
      lat_sum((Float*)&ax[i][t], 2);
      lat_sum((Float*)&vc[i][t], 2);
      lat_sum((Float*)&tn[i][t], 2);
      lat_sum((Float*)&a0[i][t], 2);
      lat_sum((Float*)&v0[i][t], 2);
      lat_sum((Float*)&psax[i][t], 2);
    }
    lat_sum((Float*)&ps_cos[t], 2);
    lat_sum((Float*)&sc_cos[t], 2);
    lat_sum((Float*)&ax_cos[t], 2);
    lat_sum((Float*)&vc_cos[t], 2);
    lat_sum((Float*)&tn_cos[t], 2);
    lat_sum((Float*)&a0_cos[t], 2);
    lat_sum((Float*)&v0_cos[t], 2);
    lat_sum((Float*)&psax_cos[t], 2);    
  }

  // Print out results
  //----------------------------------------------------------------

  FILE *fp_spec;
  fp_spec=fp;

  char src_string[10];
  sprintf(src_string,"%d",t_src);
  
  
  for (int mom_ind=0; mom_ind<n_mom; mom_ind++) {
    for (t=0; t<time_size; t++)
      Fprintf(fp_spec, "bPS_pPS_Source_Cos_%s_Sink_%s %s %d  %.16e %.16e\n", mom_info, mom_string[mom_ind], src_string, t, ps[mom_ind][t].real(),  ps[mom_ind][t].imag());
    for (t=0; t<time_size; t++)
      Fprintf(fp_spec, "bSC_pSC_Source_Cos_%s_Sink_%s %s %d  %.16e %.16e\n", mom_info, mom_string[mom_ind], src_string, t, sc[mom_ind][t].real(),  sc[mom_ind][t].imag());
    for (t=0; t<time_size; t++) 
      Fprintf(fp_spec, "bVC_pVC_Source_Cos_%s_Sink_%s %s %d  %.16e %.16e\n", mom_info, mom_string[mom_ind], src_string, t, vc[mom_ind][t].real(),  vc[mom_ind][t].imag());
    for (t=0; t<time_size; t++)
      Fprintf(fp_spec, "bAX_pAX_Source_Cos_%s_Sink_%s %s %d  %.16e %.16e\n", mom_info, mom_string[mom_ind], src_string, t, ax[mom_ind][t].real(),  ax[mom_ind][t].imag());
    for (t=0; t<time_size; t++)
      Fprintf(fp_spec, "bTN_pTN_Source_Cos_%s_Sink_%s %s %d  %.16e %.16e\n", mom_info, mom_string[mom_ind], src_string, t, tn[mom_ind][t].real(),  tn[mom_ind][t].imag());
    for (t=0; t<time_size; t++)
      Fprintf(fp_spec, "bA0_pA0_Source_Cos_%s_Sink_%s %s %d  %.16e %.16e\n", mom_info, mom_string[mom_ind], src_string, t, a0[mom_ind][t].real(),  a0[mom_ind][t].imag());
    for (t=0; t<time_size; t++)
      Fprintf(fp_spec, "bV0_pV0_Source_Cos_%s_Sink_%s %s %d  %.16e %.16e\n", mom_info, mom_string[mom_ind], src_string, t, v0[mom_ind][t].real(),  v0[mom_ind][t].imag());
    for (t=0; t<time_size; t++)
      Fprintf(fp_spec, "bPS_pA0_Source_Cos_%s_Sink_%s %s %d  %.16e %.16e\n", mom_info, mom_string[mom_ind], src_string, t, psax[mom_ind][t].real(), psax[mom_ind][t].imag());
  }
  
  for (t=0; t<time_size; t++)
    Fprintf(fp_spec, "bPS_pPS_Source_Cos_%s_Sink_Cos_%s %s %d  %.16e %.16e\n", mom_info, mom_info, src_string, t, ps_cos[t].real(),  ps_cos[t].imag());
  for (t=0; t<time_size; t++)
    Fprintf(fp_spec, "bSC_pSC_Source_Cos_%s_Sink_Cos_%s %s %d  %.16e %.16e\n", mom_info, mom_info, src_string, t, sc_cos[t].real(),  sc_cos[t].imag());
  for (t=0; t<time_size; t++) 
    Fprintf(fp_spec, "bVC_pVC_Source_Cos_%s_Sink_Cos_%s %s %d  %.16e %.16e\n", mom_info, mom_info, src_string, t, vc_cos[t].real(),  vc_cos[t].imag());
  for (t=0; t<time_size; t++)
    Fprintf(fp_spec, "bAX_pAX_Source_Cos_%s_Sink_Cos_%s %s %d  %.16e %.16e\n", mom_info, mom_info, src_string, t, ax_cos[t].real(),  ax_cos[t].imag());
  for (t=0; t<time_size; t++)
    Fprintf(fp_spec, "bTN_pTN_Source_Cos_%s_Sink_Cos_%s %s %d  %.16e %.16e\n", mom_info, mom_info, src_string, t, tn_cos[t].real(),  tn_cos[t].imag());
  for (t=0; t<time_size; t++)
    Fprintf(fp_spec, "bA0_pA0_Source_Cos_%s_Sink_Cos_%s %s %d  %.16e %.16e\n", mom_info, mom_info, src_string, t, a0_cos[t].real(),  a0_cos[t].imag());
  for (t=0; t<time_size; t++)
    Fprintf(fp_spec, "bV0_pV0_Source_Cos_%s_Sink_Cos_%s %s %d  %.16e %.16e\n", mom_info, mom_info, src_string, t, v0_cos[t].real(),  v0_cos[t].imag());
  for (t=0; t<time_size; t++)
    Fprintf(fp_spec, "bPS_pA0_Source_Cos_%s_Sink_Cos_%s %s %d  %.16e %.16e\n", mom_info, mom_info, src_string, t, psax_cos[t].real(), psax_cos[t].imag());
  

  for (int i=0; i<n_mom; i++) {
    sfree(psax[i]);
    sfree(v0[i]);
    sfree(a0[i]);
    sfree(tn[i]);
    sfree(vc[i]);
    sfree(ax[i]);
    sfree(ps[i]);
    sfree(sc[i]);
  }
  sfree(psax_cos);
  sfree(v0_cos);
  sfree(a0_cos);
  sfree(tn_cos);
  sfree(vc_cos);
  sfree(ax_cos);
  sfree(ps_cos);
  sfree(sc_cos);  
  
}

// box-box two-point functions
//------------------------------------------------------------------
void AlgThreePt::box_spectrum(QPropW& q1_src, QPropW& q2_src, int bc) {

  char *fname = "box_spectrum()";
  VRB.Func(cname,fname);
  
  /*-------------------------------------------------------------------
    Note that this gives correlators that are the complex conjugate of
    the conventions used in my notes, i.e. here the light quark
    propagator (q2_src) is Hermitian conjugated whereas in the
    convention of my notes the strange quark propagator (q1_src)
    should be Hermitian conjugated.
  -------------------------------------------------------------------*/

  //bc = 0 : periodic boundary conditions, mres/ZA contractions
  //bc = 1 : antiperiodic boundary conditions, mres/ZA contractions
  //bc = 2 : normal box_spectrum code (P+A already done)
  if ( bc<0 || bc>2)
    ERR.General(cname,fname,"bc must be 0, 1, or 2, received %d\n",bc);

  int t_src = q1_src.SourceTime();
  int t, time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex *ps, *psax;
  SMALLOC(ps,Rcomplex,time_size);   // pseudo-pseudo
  SMALLOC(psax,Rcomplex,time_size); // pseudo-axial0

  for (t=0; t<time_size; t++)
	ps[t] = psax[t] = 0.0;

  for (t=0; t<time_size; t++) {
	WilsonMatrix tmp_q1 = q1_src.WallSinkProp(t);
	WilsonMatrix tmp_q2 = q2_src.WallSinkProp(t);
    tmp_q2.hconj();                     // pion-antiquark-pion
    ps[t] += Trace(tmp_q1, tmp_q2);
	tmp_q2.gr(-5).gr(3).gr(-5);         // pion-antiquark-axial0
	psax[t] += Trace(tmp_q1, tmp_q2);
  }
  // Print out results
  //----------------------------------------------------------------

  FILE *fp_spec;
  if (bc==2)
    fp_spec=fp;
  else
    fp_spec=fp_mres_ZA;

  char src_string[10];
  if (bc==2) {
    sprintf(src_string,"%d",t_src);
  } else {
    if (bc==1 && t_src==alg_threept_arg->t_snk)
      sprintf(src_string,"A");
    else if (bc==0 && t_src==alg_threept_arg->t_src)
      sprintf(src_string,"P");
    else
      ERR.General(cname,fname,"Inconsistent bc and source time.\nt_src = %d , bc = %d\n",t_src,bc);
  }
    
  for (t=0; t<time_size; t++)
	Fprintf(fp_spec, "bPS_bPS %s %d  %.16e %.16e\n", src_string, t, ps[t].real(), ps[t].imag());
  for (t=0; t<time_size; t++)
	Fprintf(fp_spec, "bPS_bA0 %s %d  %.16e %.16e\n", src_string, t, psax[t].real(), psax[t].imag());

  sfree(psax);
  sfree(ps);

}

// box-box two-point functions for non-zero momentum
//------------------------------------------------------------------
void AlgThreePt::box_spectrum(QPropW& q_s, QPropW& q_u_mom, int mom_num, int mom_dir) {

  char *fname = "box_spectrum()";
  VRB.Func(cname,fname);
  
  /*-------------------------------------------------------------------
    Note that this gives correlators that are the complex conjugate of
    the conventions used in my notes, i.e. here the light quark
    propagator (q_u_mom) is Hermitian conjugated whereas in the
    convention of my notes the strange quark propagator (q_s)
    should be Hermitian conjugated.
  -------------------------------------------------------------------*/


  int t_src = q_u_mom.SourceTime();
  int t, time_size=GJP.Tnodes()*GJP.TnodeSites();

  int n_mom=0; //Number of different momenta given that each component can be plus or minus.
  if (mom_num==1)
    n_mom=2;
  else if (mom_num==2)
    n_mom=4;
  else if (mom_num==3)
    n_mom=8;
  else
    ERR.General(cname,fname,"Input parameter mom_num must be 1, 2, or 3, but received mom_num=%d\n",mom_num);

  if (mom_dir<0 || mom_dir>2)
    ERR.General(cname,fname,"Input parameter mom_dir must be 0, 1, or 2, but received mom_dir=%d\n",mom_dir);

  Rcomplex *ps[n_mom], *psax[n_mom];
  Rcomplex *ps_cos, *psax_cos; //Cosine sink
  for (int i=0; i<n_mom; i++) {
    SMALLOC(ps[i],Rcomplex,time_size);   // pseudo-pseudo
    SMALLOC(psax[i],Rcomplex,time_size); // pseudo-axial0
  }
  SMALLOC(ps_cos,Rcomplex,time_size);
  SMALLOC(psax_cos,Rcomplex,time_size);


  for (t=0; t<time_size; t++) {
    for (int i=0; i<n_mom; i++)
      ps[i][t] = psax[i][t] = 0.0;
    ps_cos[t] = psax_cos[t] = 0.0;
  }

  char mom_string[n_mom][100]; //String describing each momentum at the sink
  char mom_info[100]; //String describing the momentum of the cosine source
  int p[3];
  for (int i=0; i<3; i++)
    p[i]=0;
  if (mom_num==1) {
    p[mom_dir]=1;
  } else if (mom_num==2) {
    for (int mom_dir_tmp=0; mom_dir_tmp<3; mom_dir_tmp++)
      if (mom_dir_tmp!=mom_dir)
	p[mom_dir_tmp]=1;
  } else {
    for (int i=0; i<3; i++)
      p[i]=1;
  }
  sprintf(mom_info,"p=(%d,%d,%d)",p[0],p[1],p[2]);


  for (t=0; t<time_size; t++) {

    WilsonMatrix tmp_qs = q_s.WallSinkProp(t);
    
    int mom_ind=0;
    int q[3], p_snk[3];
    //Loop over all possible values of q (each component can be -1 or 1).
    //If mom_num is 1 or 2 then only execute the loop when the components in 
    //non-twisted directions are equal to -1.  Copy q into p_snk, and set
    //all components in non-twisted directions to 0 in p_snk so that it is
    //the actual momentum of the kaon.
    for (q[0]=-1; q[0]<=1; q[0]+=2)
      for (q[1]=-1; q[1]<=1; q[1]+=2)
	for (q[2]=-1; q[2]<=1; q[2]+=2) {
	  
	  for (int i=0; i<3; i++)
	    p_snk[i]=q[i];
	  
	  int do_this_iter=1;
	  if (mom_num==1)
	    for (int mom_dir_tmp=0; mom_dir_tmp<3; mom_dir_tmp++) {
	      if (mom_dir_tmp!=mom_dir)
		if (q[mom_dir_tmp]!=-1) {
		  do_this_iter=0;
		  continue;
		} else {
		  p_snk[mom_dir_tmp]=0;
		}
	    }
	  if (mom_num==2)
	    if (q[mom_dir]!=-1)
	      do_this_iter=0;
	    else {
	      p_snk[mom_dir]=0;
	    }
	  if (!do_this_iter) continue;
	  
	  WilsonMatrix tmp_qu_p_snk_h = q_u_mom.TwistMomSinkProp(t,p_snk);
	  tmp_qu_p_snk_h.hconj();                     // pion-antiquark-pion
	  ps[mom_ind][t] += Trace(tmp_qs, tmp_qu_p_snk_h);
	  tmp_qu_p_snk_h.gr(-5).gr(3).gr(-5);         // pion-antiquark-axial0
	  psax[mom_ind][t] += Trace(tmp_qs, tmp_qu_p_snk_h);
	  
	  //String to describe this momentum combination in the output
	  sprintf(mom_string[mom_ind],"p=(%d,%d,%d)",p_snk[0],p_snk[1],p_snk[2]);
	  
	  mom_ind++;
	  
	} //q[2] for loop
    //end q[i] for loops
    
    if (mom_ind!=n_mom)
      ERR.General(cname,fname,"Number of momentum directions done incorrect for some reason.\nDid %d, should have done %d.\nmom_num=%d  mom_dir=%d\n",mom_ind,n_mom,mom_num,mom_dir);
    
    //Cosine sink
    
    WilsonMatrix tmp_qu_p_cos_h = q_u_mom.TwistCosSinkProp(t,p);
    tmp_qu_p_cos_h.hconj();                   // pion-antiquark-pion
    ps_cos[t] += Trace(tmp_qs, tmp_qu_p_cos_h);
    tmp_qu_p_cos_h.gr(-5).gr(3).gr(-5);       // pion-antiquark-axial0
    psax_cos[t] += Trace(tmp_qs, tmp_qu_p_cos_h);

  } // t for loop
  
  // Print out results
  //----------------------------------------------------------------

  FILE *fp_spec;
  fp_spec=fp;

  char src_string[10];
  sprintf(src_string,"%d",t_src);
  
  for (int mom_ind=0; mom_ind<n_mom; mom_ind++) {
    for (t=0; t<time_size; t++)
        Fprintf(fp_spec, "bPS_bPS_Source_Cos_%s_Sink_%s %s %d  %.16e %.16e\n", mom_info, mom_string[mom_ind], src_string, t, ps[mom_ind][t].real(), ps[mom_ind][t].imag());
    for (t=0; t<time_size; t++)
	Fprintf(fp_spec, "bPS_bA0_Source_Cos_%s_Sink_%s %s %d  %.16e %.16e\n", mom_info, mom_string[mom_ind], src_string, t, psax[mom_ind][t].real(), psax[mom_ind][t].imag());
  }
  
  for (t=0; t<time_size; t++)
    Fprintf(fp_spec, "bPS_bPS_Source_Cos_%s_Sink_Cos_%s %s %d  %.16e %.16e\n", mom_info, mom_info, src_string, t, ps_cos[t].real(), ps_cos[t].imag());
  for (t=0; t<time_size; t++)
    Fprintf(fp_spec, "bPS_bA0_Source_Cos_%s_Sink_Cos_%s %s %d  %.16e %.16e\n", mom_info, mom_info, src_string, t, psax_cos[t].real(), psax_cos[t].imag());

  for (int i=0; i<n_mom; i++) {
    sfree(ps[i]);
    sfree(psax[i]);
  }
  sfree(ps_cos);
  sfree(psax_cos);

}

// K to vacuum diagrams (labelled I2)
//----------------------------------------------------------------
void AlgThreePt::k_to_vac(QPropW& q_str, QPropW& q_src,
						  QPropW& q_ppl) { 

  char *fname = "k_to_vac()";
  VRB.Func(cname,fname);
  ERR.NotImplemented(cname,fname); //The problem is that the current QPropW
                                   //class doesn't have the SourceWidth and
                                   //Source functions.

  /*
  int gat, trt;
  int t_src = q_str.SourceTime();
  int t_op = q_ppl.SourceTime(), width = q_ppl.SourceWidth();
  int t, time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex *op[3][4], *po[3][4];
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	SMALLOC(op[gat][trt],Rcomplex,time_size);
	SMALLOC(po[gat][trt],Rcomplex,time_size);
  }
  WilsonMatrix tmp_str, tmp_ppl; // strange, pupil
  SpinMatrix spn_Kon, spn_ppl;   // kaon, pupil
  Matrix     col_Kon, col_ppl;

  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++) 
	op[gat][trt][t] = po[gat][trt][t] = 0.0;

  int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  int vol = (GJP.VolNodeSites()/GJP.TnodeSites());
  for (int mu=-1; mu<4; mu++) {
	for (int nu=-1; nu<4; nu++) {
	  if (mu<3 && mu>-1 && nu<3 && nu>-1 && mu<=nu) continue;
	  gat = (nu<0?0:(mu<0?1:2));
	  if (!do_susy && gat!=1) continue;
	  for (int i=0; i<GJP.VolNodeSites(); i++) {
		t = i/vol + shift_t;
		Rcomplex cc = conj(q_ppl.Source(i));
		tmp_str = q_str[i];
		//tmp_str.hconj().gr(-5).gr(mu).gr(nu); // pion-antiquark-gammas
		tmp_str.hconj();
		tmp_str.gr(-5).gr(mu).gr(nu); // pion-antiquark-gammas
		tmp_ppl = cc*q_ppl[i];
		tmp_ppl.gr(mu).gr(nu).gr(-5);	      // (pupil)quark-gammas
		op[gat][TR][t] += Trace(q_src[i] * tmp_str, tmp_ppl);
		op[gat][TRTR][t] += Trace(q_src[i] , tmp_str) * tmp_ppl.Trace();
		spn_Kon = ColorTrace(q_src[i], tmp_str);
		spn_ppl = ColorTrace(tmp_ppl);
		op[gat][TR_MX][t] += Tr(spn_Kon,spn_ppl);
		col_Kon = SpinTrace(q_src[i], tmp_str);
		col_ppl = SpinTrace(tmp_ppl);
		op[gat][TRTR_MX][t] += Tr(col_Kon,col_ppl);
 		tmp_str.gr(-5); // parity swap
		tmp_ppl.gr(-5);
		po[gat][TR][t] += Trace(q_src[i] * tmp_str, tmp_ppl);
		po[gat][TRTR][t] += Trace(q_src[i] , tmp_str) * tmp_ppl.Trace();
		spn_Kon = ColorTrace(q_src[i], tmp_str);
		spn_ppl = ColorTrace(tmp_ppl);
		po[gat][TR_MX][t] += Tr(spn_Kon,spn_ppl);
		col_Kon = SpinTrace(q_src[i], tmp_str);
		col_ppl = SpinTrace(tmp_ppl);
		po[gat][TRTR_MX][t] += Tr(col_Kon,col_ppl);
	  }
	}
  }
  // Global sums
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=t_op;t<t_op+width;t++) {
    lat_sum((Float*)&op[gat][trt][t], 2);
    lat_sum((Float*)&po[gat][trt][t], 2);
  }
  // Print out results
  //----------------------------------------------------------------
  for (gat=0;gat<3;gat++) {
	if (!do_susy && gat!=1) continue;
	for (trt=0;trt<4;trt++) for (t=t_op;t<t_op+width;t++)
	  Fprintf(fp," I2%s%s %d %d %d %.16e %.16e\t%.16e %.16e\n",
			  tra[trt], gam[gat], t_src, t, 0,
			  op[gat][trt][t].real(), op[gat][trt][t].imag(),
			  po[gat][trt][t].real(), po[gat][trt][t].imag());
  }
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	sfree(op[gat][trt]);
	sfree(po[gat][trt]);
  }
  */
}

// figure eight diagrams (labelled F8 (B_K style) or F8l (K->pi style))
//----------------------------------------------------------------
void AlgThreePt::figure8(QPropW& q_str, QPropW& q_src,
						 QPropW& q_snk1, QPropW& q_snk2,
						 int is_light) {

  char *fname = "figure8()";
  VRB.Func(cname,fname);

  int gat, trt;
  int t_src = q_str.SourceTime();
  int t_snk = q_snk1.SourceTime();
  int t, time_size = GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex *oo[3][4], *pp[3][4];
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	SMALLOC(oo[gat][trt],Rcomplex,time_size);
	SMALLOC(pp[gat][trt],Rcomplex,time_size);
  }
  WilsonMatrix tmp_str, tmp_snk; // strange, sink
  SpinMatrix spn_src, spn_snk;   // source, sink
  Matrix     col_src, col_snk;

  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++) 
	oo[gat][trt][t] = pp[gat][trt][t] = 0.0;

  int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  int vol = GJP.VolNodeSites()/GJP.TnodeSites();
  for (int mu=-1; mu<4; mu++) {
	for (int nu=-1; nu<4; nu++) {
	  if (mu<3 && mu>-1 && nu<3 && nu>-1 && mu<=nu) continue;
	  gat = (nu<0?0:(mu<0?1:2));
	  if (!do_susy && gat!=1) continue;
	  for (int i=0; i<GJP.VolNodeSites(); i++) {
		t = i/vol + shift_t;
		tmp_str = q_str[i];
		//tmp_str.hconj().gr(-5).gr(mu).gr(nu);  // pion-antiquark-gammas
		tmp_str.hconj();
		tmp_str.gr(-5).gr(mu).gr(nu);  // pion-antiquark-gammas
		tmp_snk = q_snk1[i];
		//tmp_snk.hconj().gr(-5).gr(mu).gr(nu); // pion-antiquark-gammas
		tmp_snk.hconj();
		tmp_snk.gr(-5).gr(mu).gr(nu); // pion-antiquark-gammas
		oo[gat][TR][t] += Trace(q_snk2[i]*tmp_snk,q_src[i]*tmp_str);
		oo[gat][TRTR][t] += Trace(q_snk2[i],tmp_snk)*Trace(q_src[i],tmp_str);
		spn_src = ColorTrace(q_src[i],tmp_str);
		spn_snk = ColorTrace(q_snk2[i],tmp_snk);
		oo[gat][TR_MX][t] += Tr(spn_snk,spn_src);
		col_src = SpinTrace(q_src[i],tmp_str);
		col_snk = SpinTrace(q_snk2[i],tmp_snk);
		oo[gat][TRTR_MX][t] += Tr(col_snk,col_src);
		tmp_str.gr(-5);		          // parity swap
		tmp_snk.gr(-5);
		pp[gat][TR][t] += Trace(q_snk2[i]*tmp_snk,q_src[i]*tmp_str);
		pp[gat][TRTR][t] += Trace(q_snk2[i],tmp_snk)*Trace(q_src[i],tmp_str);
		spn_src = ColorTrace(q_src[i],tmp_str);
		spn_snk = ColorTrace(q_snk2[i],tmp_snk);
		pp[gat][TR_MX][t] += Tr(spn_snk,spn_src);
		col_src = SpinTrace(q_src[i],tmp_str);
		col_snk = SpinTrace(q_snk2[i],tmp_snk);
		pp[gat][TRTR_MX][t] += Tr(col_snk,col_src);
	  }
    }
  }

  // Global sums and Output the correlators
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++) {
    lat_sum((Float*)&oo[gat][trt][t], 2);
    lat_sum((Float*)&pp[gat][trt][t], 2);
  }
  // Print out results
  //----------------------------------------------------------------
  for (gat=0;gat<3;gat++) {
	if (!do_susy && gat!=1) continue;
	for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++)
	  Fprintf(fp,"F8%s%s%s %d %d %d  %.16e %.16e\t%.16e %.16e\n",
			  (is_light?"l":"\0"), tra[trt], gam2[gat],
			  t_src, t, t_snk,
			  oo[gat][trt][t].real(), oo[gat][trt][t].imag(),
			  pp[gat][trt][t].real(), pp[gat][trt][t].imag());
  }
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	sfree(oo[gat][trt]);
	sfree(pp[gat][trt]);
  }

}

// eye diagrams (labelled I3)
//  also scalar and psibar-psi (labelled VS3)
//----------------------------------------------------------------
void AlgThreePt::eye(QPropW& q_str, QPropW& q_spc,
					 QPropW& q_snk, QPropW& q_ppl) {

  char *fname = "eye()";
  VRB.Func(cname,fname);
  ERR.NotImplemented(cname,fname); //The problem is that the current QPropW
                                   //class doesn't have the SourceWidth and
                                   //Source functions.

  /*
  int gat, trt;
  int t_src = q_str.SourceTime();
  int t_op = q_ppl.SourceTime(), width = q_ppl.SourceWidth();
  int t_snk = q_snk.SourceTime();
  int t, time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex *oo[3][4], *pp[3][4], *pbp;
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	SMALLOC(oo[gat][trt],Rcomplex,time_size);
	SMALLOC(pp[gat][trt],Rcomplex,time_size);
  }
  SMALLOC(pbp,Rcomplex,time_size); // psi-bar psi
  WilsonMatrix tmp_str, tmp_snk, tmp_ppl; // strange, sink, pupil
  WilsonMatrix tmp_spc = (Float)0.0;      // spectator
  SpinMatrix spn_eye, spn_ppl; // eye, pupil
  Matrix     col_eye, col_ppl;

  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++) 
	oo[gat][trt][t] = pp[gat][trt][t] = 0.0;
  for (t=0; t<time_size; t++)
	pbp[t] = 0.0;

  int t_snk_mod, t_displace;
  t_snk_mod = t_snk>=0 ? t_snk%time_size
                           : (time_size+(t_snk%time_size))%time_size;
               // puts t_snk into the domain [0,time_size-1] using mod
               // special treatment is needed for t_snk<0
               // this is needed for example because t_snk=time_size for
               //   P-A
  t_displace = t_snk>=0 ? t_snk/time_size : (-t_snk-1)/time_size+1;
               // this is the number of domains that t_snk is away from
               //   the fundamental domain [0,time_size-1]
	       // it is >= 0
  if ( t_displace%2 == 0 )
    tmp_spc = q_spc.WallSinkProp(t_snk_mod);
  else {
    //If t_displace is odd then use the other propagator,
    //i.e. if q_spc is P+A then use P-A and vice versa.
    tmp_spc = q_snk.WallSinkProp(t_snk_mod);
  }

  int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  int vol = GJP.VolNodeSites()/GJP.TnodeSites();
  for (int mu=-1; mu<4; mu++) {
	for (int nu=-1; nu<4; nu++) {
	  if (mu<3 && mu>-1 && nu<3 && nu>-1 && mu<=nu) continue;
	  gat = (nu<0?0:(mu<0?1:2));
	  if (!do_susy && gat!=1) continue;
	  for (int i=0; i<GJP.VolNodeSites(); i++) {
		t = i/vol + shift_t;
		Rcomplex cc = conj(q_ppl.Source(i));
		tmp_str = q_str[i];
		//tmp_str.hconj().gr(-5);   // pion-antiquark
		tmp_str.hconj();
		tmp_str.gr(-5);   // pion-antiquark
		tmp_snk = q_snk[i];
		tmp_snk.gr(-5);           // quark-pion
		tmp_str.gr(mu).gr(nu);			// insert gammas
		tmp_ppl = cc*q_ppl[i];
		if (mu<0 && nu==0) pbp[t] += tmp_ppl.Trace();
		tmp_ppl.gr(mu).gr(nu);			// (pupil)quark-gammas
		oo[gat][TR][t] += Trace(tmp_snk * tmp_spc, tmp_str * tmp_ppl);
		oo[gat][TRTR][t] += Trace(tmp_snk * tmp_spc, tmp_str) * tmp_ppl.Trace();
		spn_eye = ColorTrace(tmp_snk, tmp_spc, tmp_str);
		spn_ppl = ColorTrace(tmp_ppl);
		oo[gat][TR_MX][t] += Tr(spn_eye,spn_ppl);
		col_eye = SpinTrace(tmp_snk, tmp_spc, tmp_str);
		col_ppl = SpinTrace(tmp_ppl);
		oo[gat][TRTR_MX][t] += Tr(col_eye,col_ppl);
		tmp_str.gr(-5);			// parity swap
		tmp_ppl.gr(-5);
		pp[gat][TR][t] += Trace(tmp_snk * tmp_spc, tmp_str * tmp_ppl);
		pp[gat][TRTR][t] += Trace(tmp_snk * tmp_spc, tmp_str) * tmp_ppl.Trace();
		spn_eye = ColorTrace(tmp_snk, tmp_spc, tmp_str);
		spn_ppl = ColorTrace(tmp_ppl);
		pp[gat][TR_MX][t] += Tr(spn_eye,spn_ppl);
		col_eye = SpinTrace(tmp_snk, tmp_spc, tmp_str);
		col_ppl = SpinTrace(tmp_ppl);
		pp[gat][TRTR_MX][t] += Tr(col_eye,col_ppl);
	  }
	}
  }
  
  // Global sums
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=t_op;t<t_op+width;t++) {
    lat_sum((Float*)&oo[gat][trt][t], 2);
    lat_sum((Float*)&pp[gat][trt][t], 2);
  }
  for (t=t_op; t<t_op+width; t++)
    lat_sum((Float*)&pbp[t], 2);
  // Print out results
  //----------------------------------------------------------------
  for (gat=0;gat<3;gat++) {
	if (!do_susy && gat!=1) continue;
	for (trt=0;trt<4;trt++) for (t=t_op;t<t_op+width;t++)
	  Fprintf(fp,"I3%s%s %d %d %d  %.16e %.16e\t%.16e %.16e\n",
			  tra[trt], gam2[gat], t_src, t, t_snk,
			  oo[gat][trt][t].real(), oo[gat][trt][t].imag(),
			  pp[gat][trt][t].real(), pp[gat][trt][t].imag());
  }
  for (t=t_op; t<t_op+width; t++)
	Fprintf(fp,"PBP %d  %.16e %.16e\t%.16e %.16e\n", t,
			pbp[t].real(), pbp[t].imag());

  sfree(pbp);
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	sfree(oo[gat][trt]);
	sfree(pp[gat][trt]);
  }
  */
}

// figure-eight with a spectator quark (labelled F8s), old version from v5_0_3-wme
//----------------------------------------------------------------
void AlgThreePt::figure8_spectator_old(QPropW& q_str, QPropW& q_spc,
								   QPropW& q_snk1,QPropW& q_snk2,
								   QPropW& q_snk3) {
 
  char *fname = "figure8_spectator()";
  VRB.Func(cname,fname);

  int gat, trt;
  int t_src = q_str.SourceTime();
  int t_snk = q_snk1.SourceTime();
  int t, time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex *op[3][4], *po[3][4], *sd;
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	SMALLOC(op[gat][trt],Rcomplex,time_size);
	SMALLOC(po[gat][trt],Rcomplex,time_size);
  }
  SMALLOC(sd,Rcomplex,time_size); // strange-down insertion
  WilsonMatrix tmp_str, tmp_snk;     // source, sink
  WilsonMatrix tmp_spc = (Float)0.0; // spectator
  SpinMatrix spn_eye, spn_vac; // eye, vacuum bubble
  Matrix     col_eye, col_vac;

  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++) 
	op[gat][trt][t] = po[gat][trt][t] = 0.0;
  for (t=0;t<time_size;t++) 
	sd[t] = 0.0;

/* 
   Commented out, assuming this routine is called only for P+A P-A 
   There could be a need for a different parameter checking routine (CJ)
  if ( alg_threept_arg->explicit_src_snk && (t_snk<0 || t_snk>=time_size) )
    ERR.General(cname,fname,"Source/sink time set to %d, should be between 0 and %d since explicit_src_snk is turned on.\n",t_snk,time_size);
*/
  int t_snk_mod, t_displace;
  t_snk_mod = t_snk>=0 ? t_snk%time_size
                           : (time_size+(t_snk%time_size))%time_size;
               // puts t_snk into the domain [0,time_size-1] using mod
               // special treatment is needed for t_snk<0
               // this is needed for example because t_snk=time_size for
               //   P-A
  t_displace = t_snk>=0 ? t_snk/time_size : (-t_snk-1)/time_size+1;
               // this is the number of domains that t_snk is away from
               //   the fundamental domain [0,time_size-1]
	       // it is >= 0
  if ( t_displace%2 == 0 )
    tmp_spc = q_spc.WallSinkProp(t_snk_mod);
  else {
    //If t_displace is odd then use the other propagator,
    //i.e. if q_spc is P+A then use P-A and vice versa.
    tmp_spc = q_snk1.WallSinkProp(t_snk_mod);
  }
  tmp_spc.gl(-5); // pion-quark

  int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  int vol = GJP.VolNodeSites()/GJP.TnodeSites();
  for (int mu=-1; mu<4; mu++) {
	for (int nu=-1; nu<4; nu++) {
	  if (mu>0 && nu>0 && mu<=nu) continue;
	  gat = (nu<0?0:(mu<0?1:2));
	  if (!do_susy && gat!=1) continue;
	  for (int i=0; i<GJP.VolNodeSites(); i++) {
		t = i/vol + shift_t;
		tmp_str = q_str[i];
		//tmp_str.hconj().gr(-5); // pion-antiquark
		tmp_str.hconj();
		tmp_str.gr(-5); // pion-antiquark
		tmp_snk = q_snk2[i];
		//tmp_snk.hconj().gr(-5); // pion-antiquark
		tmp_snk.hconj();
		tmp_snk.gr(-5); // pion-antiquark
		if (mu<0 && nu==0) sd[t] += Trace(tmp_spc * tmp_str, q_snk3[i]);
		tmp_str.gr(mu).gr(nu);        // gammas
		tmp_snk.gr(mu).gr(nu).gr(-5); // gammas
		op[gat][TR][t] += Trace(q_snk1[i]*tmp_spc*tmp_str,q_snk3[i]*tmp_snk);
		op[gat][TRTR][t] += Trace(q_snk1[i]*tmp_spc,tmp_str)*Trace(q_snk3[i],tmp_snk);
		spn_eye = ColorTrace(q_snk1[i]*tmp_spc*tmp_str);
		spn_vac = ColorTrace(q_snk3[i]*tmp_snk);
		op[gat][TR_MX][t] += Tr(spn_eye,spn_vac);
		col_eye = SpinTrace(q_snk1[i]*tmp_spc*tmp_str);
		col_vac = SpinTrace(q_snk3[i]*tmp_snk);
		op[gat][TRTR_MX][t] += Tr(col_eye,col_vac);
		tmp_str.gr(-5); // parity swap
		tmp_snk.gr(-5);
		po[gat][TR][t] += Trace(q_snk1[i]*tmp_spc*tmp_str,q_snk3[i]*tmp_snk);
		po[gat][TRTR][t] += Trace(q_snk1[i]*tmp_spc,tmp_str)*Trace(q_snk3[i],tmp_snk);
		spn_eye = ColorTrace(q_snk1[i]*tmp_spc*tmp_str);
		spn_vac = ColorTrace(q_snk3[i]*tmp_snk);
		po[gat][TR_MX][t] += Tr(spn_eye,spn_vac);
		col_eye = SpinTrace(q_snk1[i]*tmp_spc*tmp_str);
		col_vac = SpinTrace(q_snk3[i]*tmp_snk);
		po[gat][TRTR_MX][t] += Tr(col_eye,col_vac);
	  }
	}
  }
  
  // Global sums
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++) {
    lat_sum((Float*)&op[gat][trt][t], 2);
    lat_sum((Float*)&po[gat][trt][t], 2);
  }
  for (t=0;t<time_size;t++)
	lat_sum((Float*)&sd[t], 2);
  // Print out results
  //----------------------------------------------------------------
  for (gat=0;gat<3;gat++) {
	if (!do_susy && gat!=1) continue;
	for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++)
	  Fprintf(fp,"F8s%s%s %d %d %d  %.16e %.16e\t%.16e %.16e\n",
			  tra[trt], gam[gat], t_src, t, t_snk,
			  op[gat][trt][t].real(), op[gat][trt][t].imag(),
			  po[gat][trt][t].real(), po[gat][trt][t].imag());
  }
  for (t=0;t<time_size;t++) {
	Fprintf(fp,"SD %d %d %d  %.16e %.16e\n", t_src, t, t_snk,
			sd[t].real(), sd[t].imag());
  }
  sfree(sd);
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	sfree(op[gat][trt]);
	sfree(po[gat][trt]);
  }

}

// figure-eight with a spectator quark (labelled F8s)
//----------------------------------------------------------------
void AlgThreePt::figure8_spectator(QPropW& q_str, QPropW& q_spc,
								   QPropW& q_snk1,QPropW& q_snk2,
								   QPropW& q_snk3,
				                                   int mom_num, int mom_dir) {
  
  //Default values: mom_num=0, mom_dir=0

  /*---------------------------------------------------------------------------------
  NOTE that the figure8_spectator function has been CHANGED from the original.
  Let tK be the time of the kaon and tpi the time of the pion, and let t1, t2, t3,
  t_spc, t_str, be the source times of the propagators q_snk1, q_snk2, q_snk3,
  q_spc, and q_str respectively.  The function now requires that t1=t2=t3=tpi
  and t_str=tK.  In addition only the following three cases are allowed:
  
  (i) tK=0, tpi=time_size, t_spc=time_size
  (ii) tK=time_size, tpi=0, t_spc=time_size
  (iii) 0<tK<time_size, tpi=t_spc
  ---------------------------------------------------------------------------------*/

  char *fname = "figure8_spectator()";
  VRB.Func(cname,fname);

  //Check the source times of the propagators for consistency with the note above.
  int t, time_size=GJP.Tnodes()*GJP.TnodeSites();
  int t_src = q_str.SourceTime(); //tK
  int t_snk = q_snk1.SourceTime(); //tpi
  int t_spc = q_spc.SourceTime();
  int t2 = q_snk2.SourceTime();
  int t3 = q_snk3.SourceTime();
  if (t_snk!=t2 || t2!=t3)
    ERR.General(cname,fname,"q_snk1, q_snk2, and q_snk3 must have the same source times but got %d, %d, %d.\n",t_snk,t2,t3);
  if (t_src==0) {
    if (t_snk!=time_size || t_spc!=time_size)
      ERR.General(cname,fname,"Invalid source times: tK=%d, tpi=%d, tspc=%d.\n",t_src,t_snk,t_spc);
  } else if (t_src==time_size) {
    if (t_snk!=0 || t_spc!=time_size)
      ERR.General(cname,fname,"Invalid source times: tK=%d, tpi=%d, tspc=%d.\n",t_src,t_snk,t_spc);
  } else if ( t_src>0 && t_src<time_size) {
    if (t_spc!=t_snk)
      ERR.General(cname,fname,"Invalid source times: tK=%d, tpi=%d, tspc=%d.\n",t_src,t_snk,t_spc);
  } else {
    ERR.General(cname,fname,"Received tK=%d, but should be between 0 and %d.\n",t_src,time_size);
  }

  //Make the appropriate label for the momentum, and check that
  //the values of mom_num and mom_dir are valid.
  if (mom_dir<0 || mom_dir>2)
    ERR.General(cname,fname,"Input parameter mom_dir must be 0, 1, or 2, but received mom_dir=%d\n",mom_dir);
  char mom_string[100];
  int p[3];
  for (int i=0; i<3; i++)
    p[i]=0;
  if (mom_num==0)
    sprintf(mom_string,"");
  else {
    if (mom_num==1) {
      p[mom_dir]=1;
    } else if (mom_num==2) {
      for (int mom_dir_tmp=0; mom_dir_tmp<3; mom_dir_tmp++)
	if (mom_dir_tmp!=mom_dir)
	  p[mom_dir_tmp]=1;
    } else if (mom_num==3) {
      for (int i=0; i<3; i++)
	p[i]=1;
    } else {
      ERR.General(cname,fname,"Input parameter mom_num must be 0, 1, 2, or 3, but received mom_num=%d\n",mom_num);
    }
    sprintf(mom_string,"_Cos_p=(%d,%d,%d)",p[0],p[1],p[2]);
  }

  int gat, trt;
  Rcomplex *op[3][4], *po[3][4], *sd;
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	SMALLOC(op[gat][trt],Rcomplex,time_size);
	SMALLOC(po[gat][trt],Rcomplex,time_size);
  }
  SMALLOC(sd,Rcomplex,time_size); // strange-down insertion
  WilsonMatrix tmp_str, tmp_snk;     // source, sink
  WilsonMatrix tmp_spc = (Float)0.0; // spectator
  SpinMatrix spn_eye, spn_vac; // eye, vacuum bubble
  Matrix     col_eye, col_vac;

  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++) 
	op[gat][trt][t] = po[gat][trt][t] = 0.0;
  for (t=0;t<time_size;t++) 
	sd[t] = 0.0;

  //Calculate tmp_spc appropriate to the source times of the propagators.
  if (t_src==0 || t_src==time_size) {
    tmp_spc = q_spc.WallSinkProp(0);
    tmp_spc.gl(-5); //pion-quark
  } else {
    tmp_spc = q_spc.WallSinkProp(t_src);
    tmp_spc.hconj();
    tmp_spc.gr(-5); //pion-quark
  }

  int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  int vol = GJP.VolNodeSites()/GJP.TnodeSites();
  for (int mu=-1; mu<4; mu++) {
	for (int nu=-1; nu<4; nu++) {
	  if (mu<3 && mu>-1 && nu<3 && nu>-1 && mu<=nu) continue;
	  gat = (nu<0?0:(mu<0?1:2));
	  if (!do_susy && gat!=1) continue;
	  for (int i=0; i<GJP.VolNodeSites(); i++) {
		t = i/vol + shift_t;
		tmp_str = q_str[i];
		//tmp_str.hconj().gr(-5); // pion-antiquark
		tmp_str.hconj();
		tmp_str.gr(-5); // pion-antiquark
		tmp_snk = q_snk2[i];
		//tmp_snk.hconj().gr(-5); // pion-antiquark
		tmp_snk.hconj();
		tmp_snk.gr(-5); // pion-antiquark
		if (mu<0 && nu==0) sd[t] += Trace(tmp_spc * tmp_str, q_snk3[i]);
		tmp_str.gr(mu).gr(nu);        // gammas
		tmp_snk.gr(mu).gr(nu).gr(-5); // gammas
		op[gat][TR][t] += Trace(q_snk1[i]*tmp_spc*tmp_str,q_snk3[i]*tmp_snk);
		op[gat][TRTR][t] += Trace(q_snk1[i]*tmp_spc,tmp_str)*Trace(q_snk3[i],tmp_snk);
		spn_eye = ColorTrace(q_snk1[i]*tmp_spc*tmp_str);
		spn_vac = ColorTrace(q_snk3[i]*tmp_snk);
		op[gat][TR_MX][t] += Tr(spn_eye,spn_vac);
		col_eye = SpinTrace(q_snk1[i]*tmp_spc*tmp_str);
		col_vac = SpinTrace(q_snk3[i]*tmp_snk);
		op[gat][TRTR_MX][t] += Tr(col_eye,col_vac);
		tmp_str.gr(-5); // parity swap
		tmp_snk.gr(-5);
		po[gat][TR][t] += Trace(q_snk1[i]*tmp_spc*tmp_str,q_snk3[i]*tmp_snk);
		po[gat][TRTR][t] += Trace(q_snk1[i]*tmp_spc,tmp_str)*Trace(q_snk3[i],tmp_snk);
		spn_eye = ColorTrace(q_snk1[i]*tmp_spc*tmp_str);
		spn_vac = ColorTrace(q_snk3[i]*tmp_snk);
		po[gat][TR_MX][t] += Tr(spn_eye,spn_vac);
		col_eye = SpinTrace(q_snk1[i]*tmp_spc*tmp_str);
		col_vac = SpinTrace(q_snk3[i]*tmp_snk);
		po[gat][TRTR_MX][t] += Tr(col_eye,col_vac);
	  }
	}
  }
  
  // Global sums
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++) {
    lat_sum((Float*)&op[gat][trt][t], 2);
    lat_sum((Float*)&po[gat][trt][t], 2);
  }
  for (t=0;t<time_size;t++)
	lat_sum((Float*)&sd[t], 2);
  // Print out results
  //----------------------------------------------------------------
  for (gat=0;gat<3;gat++) {
	if (!do_susy && gat!=1) continue;
	for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++)
	  Fprintf(fp,"F8s%s%s%s %d %d %d  %.16e %.16e\t%.16e %.16e\n",
			  tra[trt], gam[gat], mom_string, t_src, t, t_snk,
			  op[gat][trt][t].real(), op[gat][trt][t].imag(),
			  po[gat][trt][t].real(), po[gat][trt][t].imag());
  }
  for (t=0;t<time_size;t++) {
	Fprintf(fp,"SD%s %d %d %d  %.16e %.16e\n", mom_string, t_src, t, t_snk,
			sd[t].real(), sd[t].imag());
  }
  sfree(sd);
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	sfree(op[gat][trt]);
	sfree(po[gat][trt]);
  }

}

// figure-eight with extra quark pair from vacuum (labelled F8v)
//----------------------------------------------------------------
void AlgThreePt::figure8_vacuum(QPropW& q_str, QPropW& q_src,
								QPropW& q_snk1, QPropW& q_vac,
								QPropW& q_snk2) {
 
  char *fname = "figure8_vacuum()";
  VRB.Func(cname,fname);

  int gat, trt;
  int t_src = q_str.SourceTime();
  int t_snk = q_snk1.SourceTime();
  int t, time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex *op[3][4], *po[3][4];
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	SMALLOC(op[gat][trt],Rcomplex,time_size);
	SMALLOC(po[gat][trt],Rcomplex,time_size);
  }
  WilsonMatrix tmp_str, tmp_snk;     // strange, sink
  WilsonMatrix tmp_vac = (Float)0.0; // quark-pair out of vacuum
  SpinMatrix spn_2pi, spn_Kon; // two pions, kaon
  Matrix     col_2pi, col_Kon;

  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++) 
	op[gat][trt][t] = po[gat][trt][t] = 0.0;

  int t_snk_mod, t_displace;
  t_snk_mod = t_snk>=0 ? t_snk%time_size
                           : (time_size+(t_snk%time_size))%time_size;
               // puts t_snk into the domain [0,time_size-1] using mod
               // special treatment is needed for t_snk<0
               // this is needed for example because t_snk=time_size for
               //   P-A
  t_displace = t_snk>=0 ? t_snk/time_size : (-t_snk-1)/time_size+1;
               // this is the number of domains that t_snk is away from
               //   the fundamental domain [0,time_size-1]
	       // it is >= 0
  if ( t_displace%2 == 0 )
    tmp_vac = q_vac.WallSinkProp(t_snk_mod);
  else {
    //If t_displace is odd then use the other propagator,
    //i.e. if q_spc is P+A then use P-A and vice versa.
    tmp_vac = q_src.WallSinkProp(t_snk_mod);
  }
  tmp_vac.gl(-5); // pion-quark

  int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  int vol = GJP.VolNodeSites()/GJP.TnodeSites();
  for (int mu=-1; mu<4; mu++) {
	for (int nu=-1; nu<4; nu++) {
	  if (mu<3 && mu>-1 && nu<3 && nu>-1 && mu<=nu) continue;
	  gat = (nu<0?0:(mu<0?1:2));
	  if (!do_susy && gat!=1) continue;
	  for (int i=0; i<GJP.VolNodeSites(); i++) {
		t = i/vol + shift_t;
		tmp_str = q_str[i];
		//tmp_str.hconj().gr(-5).gr(mu).gr(nu);	     // pion-antiquark-gammas
		tmp_str.hconj();
		tmp_str.gr(-5).gr(mu).gr(nu);	     // pion-antiquark-gammas
		tmp_snk = q_snk1[i];
		//tmp_snk.hconj().gr(-5).gr(mu).gr(nu).gr(-5); // pion-antiquark-gammas
		tmp_snk.hconj();
		tmp_snk.gr(-5).gr(mu).gr(nu).gr(-5); // pion-antiquark-gammas
		op[gat][TR][t] += Trace(q_src[i]*tmp_str,q_snk2[i]*tmp_vac*tmp_snk);
		op[gat][TRTR][t] += Trace(q_src[i],tmp_str)*Trace(q_snk2[i],tmp_vac*tmp_snk);
		spn_Kon = ColorTrace(q_src[i]*tmp_str);
		spn_2pi = ColorTrace(q_snk2[i]*tmp_vac*tmp_snk);
		op[gat][TR_MX][t] += Tr(spn_Kon,spn_2pi);
		col_Kon = SpinTrace(q_src[i]*tmp_str);
		col_2pi = SpinTrace(q_snk2[i]*tmp_vac*tmp_snk);
		op[gat][TRTR_MX][t] += Tr(col_Kon,col_2pi);
		tmp_str.gr(-5); // parity swap
		tmp_snk.gr(-5);
		po[gat][TR][t] += Trace(q_src[i]*tmp_str,q_snk2[i]*tmp_vac*tmp_snk);
		po[gat][TRTR][t] += Trace(q_src[i],tmp_str)*Trace(q_snk2[i],tmp_vac*tmp_snk);
		spn_Kon = ColorTrace(q_src[i]*tmp_str);
		spn_2pi = ColorTrace(q_snk2[i]*tmp_vac*tmp_snk);
		po[gat][TR_MX][t] += Tr(spn_Kon,spn_2pi);
		col_Kon = SpinTrace(q_src[i]*tmp_str);
		col_2pi = SpinTrace(q_snk2[i]*tmp_vac*tmp_snk);
		po[gat][TRTR_MX][t] += Tr(col_Kon,col_2pi);
	  }
	}
  }

  // Global sums
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++) {
    lat_sum((Float*)&op[gat][trt][t], 2);
    lat_sum((Float*)&po[gat][trt][t], 2);
  }
  // Print out results
  //----------------------------------------------------------------
  for (gat=0;gat<3;gat++) {
	if (!do_susy && gat!=1) continue;
	for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++)
	  Fprintf(fp,"F8v%s%s %d %d %d  %.16e %.16e\t%.16e %.16e\n",
			  tra[trt], gam[gat], t_src, t, t_snk,
			  op[gat][trt][t].real(), op[gat][trt][t].imag(),
			  po[gat][trt][t].real(), po[gat][trt][t].imag());
  }
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	sfree(op[gat][trt]);
	sfree(po[gat][trt]);
  }

}

// eye diagrams with extra quark pair from vacuum (labelled I3v)
//----------------------------------------------------------------
void AlgThreePt::eye_vacuum(QPropW& q_str, QPropW& q_spc,
							QPropW& q_vac, QPropW& q_snk,
							QPropW& q_ppl) {
 
  char *fname = "eye_vacuum()";
  VRB.Func(cname,fname);
  ERR.NotImplemented(cname,fname); //The problem is that the current QPropW
                                   //class doesn't have the SourceWidth and
                                   //Source functions.

  /*
  int gat, trt;
  int t_src = q_str.SourceTime();
  int t_snk = q_snk.SourceTime();
  int t_op = q_ppl.SourceTime(), width = q_ppl.SourceWidth();
  int t, time_size=GJP.Tnodes()*GJP.TnodeSites();
  Rcomplex *op[3][4], *po[3][4];
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	SMALLOC(op[gat][trt],Rcomplex,time_size);
	SMALLOC(po[gat][trt],Rcomplex,time_size);
  }
  WilsonMatrix tmp_str, tmp_ppl;   // strange, pupil
  WilsonMatrix tmp_spc = (Float)0.0;      // spectator
  WilsonMatrix tmp_vac = (Float)0.0;      // quark-pair from vacuum
  SpinMatrix spn_eye, spn_ppl; // eye, pupil
  Matrix     col_eye, col_ppl;

  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=0;t<time_size;t++) 
	op[gat][trt][t] = po[gat][trt][t] = 0.0;

  int t_snk_mod, t_displace;
  t_snk_mod = t_snk>=0 ? t_snk%time_size
                           : (time_size+(t_snk%time_size))%time_size;
               // puts t_snk into the domain [0,time_size-1] using mod
               // special treatment is needed for t_snk<0
               // this is needed for example because t_snk=time_size for
               //   P-A
  t_displace = t_snk>=0 ? t_snk/time_size : (-t_snk-1)/time_size+1;
               // this is the number of domains that t_snk is away from
               //   the fundamental domain [0,time_size-1]
	       // it is >= 0
  if ( t_displace%2 == 0 ) {
    tmp_spc = q_spc.WallSinkProp(t_snk_mod);
    tmp_vac = q_vac.WallSinkProp(t_snk_mod);
  } else {
    //If t_displace is odd then use the other propagator,
    //i.e. if q_spc is P+A then use P-A and vice versa.
    tmp_spc = q_vac.WallSinkProp(t_snk_mod);
    tmp_vac = q_spc.WallSinkProp(t_snk_mod);
  }
  tmp_vac.hconj(); // pion-antiquark-pion

  int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  int vol = GJP.VolNodeSites()/GJP.TnodeSites();
  for (int mu=-1; mu<4; mu++) {
	for (int nu=-1; nu<4; nu++) {
	  if (mu<3 && mu>-1 && nu<3 && nu>-1 && mu<=nu) continue;
	  gat = (nu<0?0:(mu<0?1:2));
	  if (!do_susy && gat!=1) continue;
	  for (int i=0; i<GJP.VolNodeSites(); i++) {
		t = i/vol + shift_t;
		Rcomplex cc = conj(q_ppl.Source(i));
		tmp_str = q_str[i];
		//tmp_str.hconj().gr(-5).gr(mu).gr(nu); // pion-antiquark-gammas
		tmp_str.hconj();
		tmp_str.gr(-5).gr(mu).gr(nu); // pion-antiquark-gammas
		tmp_ppl = cc*q_ppl[i];
		tmp_ppl.gr(mu).gr(nu).gr(-5);		  // (pupil)quark-gammas
		op[gat][TR][t] += Trace(q_snk[i]*tmp_vac*tmp_spc*tmp_str,tmp_ppl);
		op[gat][TRTR][t] += Trace(q_snk[i]*tmp_vac,tmp_spc*tmp_str)*tmp_ppl.Trace();
		spn_eye = ColorTrace(q_snk[i]*tmp_vac*tmp_spc*tmp_str);
		spn_ppl = ColorTrace(tmp_ppl);
		op[gat][TR_MX][t] += Tr(spn_eye,spn_ppl);
		col_eye = SpinTrace(q_snk[i]*tmp_vac*tmp_spc*tmp_str);
		col_ppl = SpinTrace(tmp_ppl);
		op[gat][TRTR_MX][t] += Tr(col_eye,col_ppl);
		tmp_str.gr(-5); // parity swap
		tmp_ppl.gr(-5);
		po[gat][TR][t] += Trace(q_snk[i]*tmp_vac*tmp_spc*tmp_str,tmp_ppl);
		po[gat][TRTR][t] += Trace(q_snk[i]*tmp_vac,tmp_spc*tmp_str)*tmp_ppl.Trace();
		spn_eye = ColorTrace(q_snk[i]*tmp_vac*tmp_spc*tmp_str);
		spn_ppl = ColorTrace(tmp_ppl);
		po[gat][TR_MX][t] += Tr(spn_eye,spn_ppl);
		col_eye = SpinTrace(q_snk[i]*tmp_vac*tmp_spc*tmp_str);
		col_ppl = SpinTrace(tmp_ppl);
		po[gat][TRTR_MX][t] += Tr(col_eye,col_ppl);
	  }
	}
  }
  
  // Global sums
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) for (t=t_op;t<t_op+width;t++) {
    lat_sum((Float*)&op[gat][trt][t], 2);
    lat_sum((Float*)&po[gat][trt][t], 2);
  }
  // Print out results
  //----------------------------------------------------------------
  for (gat=0;gat<3;gat++) {
	if (!do_susy && gat!=1) continue;
	for (trt=0;trt<4;trt++) for (t=t_op;t<t_op+width;t++)
	  Fprintf(fp,"I3v%s%s %d %d %d  %.16e %.16e\t%.16e %.16e\n",
			  tra[trt], gam[gat], t_src, t, t_snk,
			  op[gat][trt][t].real(), op[gat][trt][t].imag(),
			  po[gat][trt][t].real(), po[gat][trt][t].imag());
  }
  for (gat=0;gat<3;gat++) for (trt=0;trt<4;trt++) {
	sfree(op[gat][trt]);
	sfree(po[gat][trt]);
  }
  */
}

//---------------------------------------------------------------------------
//Checksum, machine synchronization and timing checkpoint function
//---------------------------------------------------------------------------
void chkpt(const int num_nodes,int& chkpoint_no,const Float dtime_start,Float& dtime_last,Float& dtime_now)
{
  int dummy=0;
  Float test_val;
  Float* ptest_val=&test_val;
  char *cname = "AlgThreePt";
  char *fname = "chkpt";
  ERR.HdwCheck(cname,fname);
  
  for (int ii=0; ii<num_nodes; ii++){
    test_val=0.0;
    if (UniqueID()==ii){
      test_val=1.0;
    }
    glb_sum_five(ptest_val);
    while(test_val==0.0)
      dummy=0; //Does absolutely nothing, here as a placeholder for while
  }
  VRB.Flow(cname,fname,"Checkpoint no. %d reached.\n",chkpoint_no++);

  dtime_last=dtime_now;
  dtime_now=dclock();
  
  Float time_tmp=dtime_now-dtime_last;
  int hr_tmp=time_tmp/3600.0;
  int min_tmp=(time_tmp-3600.0*hr_tmp)/60.0;
  Float sec_tmp=time_tmp-3600.0*hr_tmp-60.0*min_tmp;
  VRB.Flow(cname,fname,"Time since last checkpoint: %d hours %d minutes %f seconds.\n",hr_tmp,min_tmp,sec_tmp);
  
  time_tmp=dtime_now-dtime_start;
  hr_tmp=time_tmp/3600.0;
  min_tmp=(time_tmp-3600.0*hr_tmp)/60.0;
  sec_tmp=time_tmp-3600.0*hr_tmp-60.0*min_tmp;
  VRB.Flow(cname,fname,"Time since beginning: %d hours %d minutes %f seconds.\n",hr_tmp, min_tmp,sec_tmp);
}
//---------------------------------------------------------------------------

CPS_END_NAMESPACE
