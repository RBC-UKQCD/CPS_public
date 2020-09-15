//CK: test to ensure the 1/2(1 +/- \gamma^5) projection operators for WilsonMatrix, glPL, glPR, grPL and grPR, are working correctly


#include<config.h>

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include<util/lattice.h>
#include<util/random.h>
#include<util/time_cps.h>

#include<alg/alg_hmc.h>
#include<alg/common_arg.h>
#include<alg/hmc_arg.h>
#include<alg/hmd_arg.h>

#include<alg/alg_int.h>
#include<alg/int_arg.h>

#include<alg/no_arg.h>
#include<alg/do_arg.h>
#include<alg/alg_plaq.h>
#include<alg/alg_pbp.h>
#include<alg/pbp_arg.h>
#include<alg/alg_remez.h>
#include<alg/alg_wline.h>
#include<alg/wilson_matrix.h>

#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/qcdio.h>
#include<util/WriteLatticePar.h>
#include<util/ReadLatticePar.h>
#include<util/qioarg.h>
#include<util/command_line.h>
//#include <sys/bgl/bgl_sys_all.h>

#undef USE_SCU_CHECKSUMS
#ifdef USE_SCU_CHECKSUMS
#include <qcdocos/scu_checksum.h>
#endif
//--------------------------------------------------------------

USING_NAMESPACE_CPS
using namespace std;

static inline void conv_i_spincolour(const int &i, int &spn, int &clr){
  //map an index 0..11 to spin and colour indices
  
  int rem = i;
  spn = rem % 4;
  rem/=4;
  clr = rem % 3;
}
  
static Complex& elem(WilsonMatrix &of, const int &i, const int &j){
  int spni, spnj, clri, clrj;
  conv_i_spincolour(i,spni,clri);
  conv_i_spincolour(j,spnj,clrj);

  return of(spni,clri,spnj,clrj);
}

static void compare(WilsonMatrix &a, WilsonMatrix &b, const char* descr){
  printf("Comparing %s...",descr);
  for(int i=0;i<12;i++){
    for(int j=0;j<12;j++){
      Rcomplex diff = elem(a,i,j) - elem(b,i,j);
      if(fabs(diff.real())>1e-10 || fabs(diff.imag())>1e-10 ){
	int spni, spnj, clri, clrj;
	conv_i_spincolour(i,spni,clri);
	conv_i_spincolour(j,spnj,clrj);

	printf("fail on (%d %d,%d %d) : (%.3f,%.3f) vs (%.3f,%.3f)\n", descr, spni, clri, spnj, clrj, elem(a,i,j).real(), elem(a,i,j).imag(), elem(b,i,j).real(), elem(b,i,j).imag() );
	exit(-1);
      }
    }
  }
  printf("pass\n");
}

static void print(WilsonMatrix &a){
  for(int i=0;i<12;i++){
    for(int j=0;j<12;j++){
      printf("(%.3f,%.3f) ",elem(a,i,j).real(), elem(a,i,j).imag() );
    }
    printf("\n");
  }
}


int main(int argc, char *argv[])
{ 
  Start(&argc,&argv);
  CommandLine::is(argc,argv);

  char *cname=argv[0];
  char *fname="main()";
  Float dtime;

  int size[] = {4,4,4,4,4};

  DoArg do_arg;
  do_arg.x_sites = size[0];
  do_arg.y_sites = size[1];
  do_arg.z_sites = size[2];
  do_arg.t_sites = size[3];
  do_arg.s_sites = size[4];
  do_arg.x_node_sites = 0;
  do_arg.y_node_sites = 0;
  do_arg.z_node_sites = 0;
  do_arg.t_node_sites = 0;
  do_arg.s_node_sites = 0;
  do_arg.x_nodes = 0;
  do_arg.y_nodes = 0;
  do_arg.z_nodes = 0;
  do_arg.t_nodes = 0;
  do_arg.s_nodes = 0;
  do_arg.updates = 0;
  do_arg.measurements = 0;
  do_arg.measurefreq = 0;
  do_arg.cg_reprod_freq = 10;
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_conf_load_addr = 0x0;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.start_seed_filename = "../rngs/ckpoint_rng.0";
  do_arg.start_conf_filename = "../configurations/ckpoint_lat.0";
  do_arg.start_conf_alloc_flag = 6;
  do_arg.wfm_alloc_flag = 2;
  do_arg.wfm_send_alloc_flag = 2;
  do_arg.start_seed_value = 83209;
  do_arg.beta =   2.25;
  do_arg.c_1 =   -3.3100000000000002e-01;
  do_arg.u0 =   1.0000000000000000e+00;
  do_arg.dwf_height =   1.8000000000000000e+00;
  do_arg.dwf_a5_inv =   1.0000000000000000e+00;
  do_arg.power_plaq_cutoff =   0.0000000000000000e+00;
  do_arg.power_plaq_exponent = 0;
  do_arg.power_rect_cutoff =   0.0000000000000000e+00;
  do_arg.power_rect_exponent = 0;
  do_arg.verbose_level = -1202;//VERBOSE_DEBUG_LEVEL; //-1202;
  do_arg.checksum_level = 0;
  do_arg.exec_task_list = 0;
  do_arg.xi_bare =   1.0000000000000000e+00;
  do_arg.xi_dir = 3;
  do_arg.xi_v =   1.0000000000000000e+00;
  do_arg.xi_v_xi =   1.0000000000000000e+00;
  do_arg.clover_coeff =   0.0000000000000000e+00;
  do_arg.clover_coeff_xi =   0.0000000000000000e+00;
  do_arg.xi_gfix =   1.0000000000000000e+00;
  do_arg.gfix_chkb = 1;
  do_arg.asqtad_KS =   0.0000000000000000e+00;
  do_arg.asqtad_naik =   0.0000000000000000e+00;
  do_arg.asqtad_3staple =   0.0000000000000000e+00;
  do_arg.asqtad_5staple =   0.0000000000000000e+00;
  do_arg.asqtad_7staple =   0.0000000000000000e+00;
  do_arg.asqtad_lepage =   0.0000000000000000e+00;
  do_arg.p4_KS =   0.0000000000000000e+00;
  do_arg.p4_knight =   0.0000000000000000e+00;
  do_arg.p4_3staple =   0.0000000000000000e+00;
  do_arg.p4_5staple =   0.0000000000000000e+00;
  do_arg.p4_7staple =   0.0000000000000000e+00;
  do_arg.p4_lepage =   0.0000000000000000e+00;

 
  GJP.Initialize(do_arg);
  LRG.Initialize();

  LRG.AssignGenerator(0);

  WilsonMatrix wmat;
  for(int i=0;i<12;i++){
    for(int j=0;j<12;j++){
      elem(wmat,i,j) = Rcomplex( LRG.Urand(), LRG.Urand() );
      printf("(%.3f,%.3f) ",elem(wmat,i,j).real(), elem(wmat,i,j).imag() );
    }
    printf("\n");
  }

  //check g5
  {
    WilsonMatrix a = wmat;
    a.gl(-5);
    //first 6 rows should be +
    //second 6 -
    
    WilsonMatrix b = wmat;
    for(int srow =2;srow<4;srow++)
      for(int crow=0;crow<3;crow++)
	for(int scol=0;scol<4;scol++)
	  for(int ccol=0;ccol<3;ccol++)
	    b(srow,crow,scol,ccol) *= -1.0;
    
    compare(a,b,"g5");
  }

  WilsonMatrix one(0.0);
  for(int i=0;i<12;i++) elem(one,i,i) = Rcomplex(1.0,0.0);

  //printf("One:");
  //print(one);

  { //glPL
    WilsonMatrix a = wmat; a.gl(-5);
    a = 0.5*( wmat - a );

    WilsonMatrix b = wmat; b.glPL();
    compare(a,b,"glPL");
  }
  { //glPR
    WilsonMatrix a = wmat; a.gl(-5);
    a = 0.5*( wmat + a );

    WilsonMatrix b = wmat; b.glPR();
    compare(a,b,"glPR");
  }
  { //grPL
    WilsonMatrix a = wmat; a.gr(-5);
    a = 0.5*( wmat - a );

    WilsonMatrix b = wmat; b.grPL();
    compare(a,b,"grPL");
  }
  { //grPR
    WilsonMatrix a = wmat; a.gr(-5);
    a = 0.5*( wmat + a );

    WilsonMatrix b = wmat; b.grPR();
    compare(a,b,"grPR");
  }

  printf("End of program\n");
  return(0);
}

