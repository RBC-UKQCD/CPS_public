#include <stdio.h>
#include <stdlib.h>	// exit()
#if TARGET==QCDOC
#include <sysfunc_cps.h>
#endif
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <alg/alg_plaq.h>
#include <alg/alg_fourier_prop.h>
#include <alg/do_arg.h>
#include <alg/common_arg.h>
#include <alg/fourierprop_arg.h>
#include <alg/no_arg.h>
#include <alg/fix_gauge_arg.h>
#include <alg/alg_fix_gauge.h>
#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

//#include <util/CommandLine.h>

int main(int argc,char *argv[])
{
  printf("Program starting...\n");

  DoArg do_arg;
  
  do_arg.x_nodes      = SizeX();
  do_arg.y_nodes      = SizeY();
  do_arg.z_nodes      = SizeZ();
  do_arg.t_nodes      = SizeT();
  do_arg.s_nodes      = SizeS();

  do_arg.x_sites = 16;
  do_arg.y_sites = 16;
  do_arg.z_sites = 16;
  do_arg.t_sites = 32;
  do_arg.s_sites = 12;
  
  do_arg.x_node_sites = do_arg.x_sites/do_arg.x_nodes;
  do_arg.y_node_sites = do_arg.y_sites/do_arg.y_nodes;
  do_arg.z_node_sites = do_arg.z_sites/do_arg.z_nodes;
  do_arg.t_node_sites = do_arg.t_sites/do_arg.t_nodes;
  do_arg.s_node_sites = do_arg.s_sites/do_arg.s_nodes;

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;
  
   // type of starting lattice

  do_arg.start_conf_kind = START_CONF_LOAD ; 
  
  do_arg.start_conf_load_addr = 
    (u_long)smalloc(do_arg.x_node_sites * do_arg.y_node_sites *
		    do_arg.z_node_sites * do_arg.t_node_sites * 4 * sizeof(Matrix));

  do_arg.start_conf_filename = argv[1];

  // type of the random number generator
  
  do_arg.start_seed_kind      = START_SEED_FIXED ;
  do_arg.start_seed_value     = 3529;

  //  do_arg.colors               = 3;

  do_arg.beta                 = 6.0; // not used 

  //  do_arg.verbose_level        = CommandLine::arg_as_int  ();

  do_arg.dwf_height           = 1.8;

  GJP.Initialize(do_arg);


  //----------------------------------------------------------------
  // Set verbose level
  //----------------------------------------------------------------

  VRB.Level(VERBOSE_RESULT_LEVEL);

  //----------------------------------------------------------------
  // Initialize argument structures
  //----------------------------------------------------------------

  /*
    Declare lattice
  */
  
  GimprRectFdwf lattice;

  /* 
     Load Lattice
  */

  {
    ReadLatticeParallel rd;
    rd.read(lattice, do_arg.start_conf_filename);

    WriteLatticeParallel wt;
    wt.write(lattice, "checklat.lat");
  }

  /*
    Restore the random number stream
  */

  /* random stream disabled for debugging purpose -- Shu Li 12/19/2005 

  // random number stream ?
  const int rnd_stream(CommandLine::arg_as_int());

  // stream lives here
  unsigned *rand_buffer((unsigned*)CommandLine::arg_as_int()); 

  if ( rnd_stream == 1 )
    {
      RandomGenerator Rand;
      Rand.RestoreSeeds( rand_buffer );
    }
  
  */


  CommonArg common_arg;
  CommonArg common_arg_gf;
  CommonArg common_arg_plaq;

  NoArg plaq_arg;
 
  // get STARTING PLAQ 
  {

    common_arg_plaq.results = (void*)"plaq.dat";
    AlgPlaq plaq(lattice,&common_arg_plaq,&plaq_arg);
    plaq.run();

  }

  //----------------------------------------------------------------
  // fix the gauge
  //----------------------------------------------------------------
  
  
  printf("**********************\nFixing Gauge\n********************\n");
  
  FixGaugeArg fix_arg;
  fix_arg.fix_gauge_kind   = FIX_GAUGE_LANDAU;
  fix_arg.hyperplane_start = 0;
  fix_arg.hyperplane_step  = 0;
  fix_arg.hyperplane_num   = 0;
  fix_arg.stop_cond        = 1e-8;
  fix_arg.max_iter_num     = 0;  // disable GF, for debugging, was 20000 --Shu Li 12/19/05
  

  {
    AlgFixGauge fix_gauge(lattice,&common_arg_gf,&fix_arg);
    fix_gauge.run();
  } 

  
  //----------------------------------------------------------------
  // Fourier arguments                               
  //----------------------------------------------------------------
  
  FourierPropArg fourierprop_arg;
  
  fourierprop_arg.cg.stop_rsd     = 1e-8;
  fourierprop_arg.cg.max_num_iter = 5000;
  fourierprop_arg.x_src           = 0;
  fourierprop_arg.y_src           = 0;
  fourierprop_arg.z_src           = 0;
  fourierprop_arg.t_src           = 0;
  fourierprop_arg.results = (void*)"gffprop_point.dat";
  
  //
  // point-point propagator
  //

  const int min_x = -2;
  const int max_x = 2;
  const int min_y = -2;
  const int max_y = 2;
  const int min_z = 0;
  const int max_z = 2;
  const int min_t = 0;
  const int max_t = 4;

  
  int four_sum(0);
  int p1,p2,p3,p4;
  for (p1=min_x;p1<=max_x;++p1)
    for (p2=min_y;p2<=max_y;++p2)
      for (p3=min_z;p3<=max_z;++p3)
	for (p4=min_t;p4<=max_t;++p4) 
          four_sum++;
  
  fourierprop_arg.plist.size(four_sum);

  four_sum=0;
  for ( p1=min_x; p1<=max_x; ++p1 )
    {
      for ( p2=min_y; p2<=max_y; ++p2 )
        {
          for ( p3=min_z; p3<=max_z; ++p3 )
            {
              for ( p4=min_t; p4<=max_t; ++p4 )
                {
                  fourierprop_arg.plist[four_sum]= FourMom( p1, p2, p3, p4 );
                  four_sum++;
                }
            } 
        }    
    } 
  

  // 
  // volume source propagators
  //
  fourierprop_arg.smom.size(4);
 
  MomentaList *smom_comp 
    = (MomentaList*) smalloc( fourierprop_arg.smom.size() * sizeof(MomentaList));
  if (smom_comp == 0x0) ERR.Pointer(":","main","smom_comp");
  
  fourierprop_arg.smom_comp=smom_comp;
  
  smom_comp[0].size(16);

  fourierprop_arg.smom[0]=FourMom(1,1,1,2);

  smom_comp[0][0] = FourMom(1,1,1,2);
  smom_comp[0][1] = FourMom(-1,1,1,2);
  smom_comp[0][2] = FourMom(1,-1,1,2);
  smom_comp[0][3] = FourMom(1,1,-1,2);
  smom_comp[0][4] = FourMom(1,1,1,-2);
  smom_comp[0][5] = FourMom(-1,-1,1,2);
  smom_comp[0][6] = FourMom(1,-1,-1,2);
  smom_comp[0][7] = FourMom(1,1,-1,-2);
  smom_comp[0][8] = FourMom(-1,1,-1,2);
  smom_comp[0][9] = FourMom(-1,1,1,-2);
  smom_comp[0][10]= FourMom(1,-1,1,-2);
  smom_comp[0][11]= FourMom(-1,-1,-1,2);
  smom_comp[0][12]= FourMom(1,-1,-1,-2);
  smom_comp[0][13]= FourMom(-1,1,-1,-2);
  smom_comp[0][14]= FourMom(-1,-1,1,-2);
  smom_comp[0][15]= FourMom(-1,-1,-1,-2);
  
  smom_comp[3].size(16);

  fourierprop_arg.smom[3]=FourMom(2,2,2,4);

  smom_comp[3][0] = FourMom(2,2,2,4);
  smom_comp[3][1] = FourMom(-2,2,2,4);
  smom_comp[3][2] = FourMom(2,-2,2,4);
  smom_comp[3][3] = FourMom(2,2,-2,4);
  smom_comp[3][4] = FourMom(2,2,2,-4);
  smom_comp[3][5] = FourMom(-2,-2,2,4);
  smom_comp[3][6] = FourMom(2,-2,-2,4);
  smom_comp[3][7] = FourMom(2,2,-2,-4);
  smom_comp[3][8] = FourMom(-2,2,-2,4);
  smom_comp[3][9] = FourMom(-2,2,2,-4);
  smom_comp[3][10]= FourMom(2,-2,2,-4);
  smom_comp[3][11]= FourMom(-2,-2,-2,4);
  smom_comp[3][12]= FourMom(2,-2,-2,-4);
  smom_comp[3][13]= FourMom(-2,2,-2,-4);
  smom_comp[3][14]= FourMom(-2,-2,2,-4);
  smom_comp[3][15]= FourMom(-2,-2,-2,-4);

  smom_comp[1].size(24);

  fourierprop_arg.smom[1]=FourMom(0,2,2,0);

  smom_comp[1][0] = FourMom(0,2,2,0);
  smom_comp[1][1] = FourMom(0,-2,2,0);
  smom_comp[1][2] = FourMom(0,2,-2,0);
  smom_comp[1][3] = FourMom(0,-2,-2,0);
  smom_comp[1][4] = FourMom(2,2,0,0);
  smom_comp[1][5] = FourMom(-2,2,0,0);
  smom_comp[1][6] = FourMom(2,-2,0,0);
  smom_comp[1][7] = FourMom(-2,-2,0,0);
  smom_comp[1][8] = FourMom(2,0,2,0);
  smom_comp[1][9] = FourMom(-2,0,2,0);
  smom_comp[1][10]= FourMom(2,0,-2,0);
  smom_comp[1][11]= FourMom(-2,0,-2,0);
  smom_comp[1][12]= FourMom(0,2,0,4);
  smom_comp[1][13]= FourMom(0,-2,0,4);
  smom_comp[1][14]= FourMom(0,2,0,-4);
  smom_comp[1][15]= FourMom(0,-2,0,-4);
  smom_comp[1][16]= FourMom(2,0,0,4);
  smom_comp[1][17]= FourMom(-2,0,0,4);
  smom_comp[1][18]= FourMom(2,0,0,-4);
  smom_comp[1][19]= FourMom(-2,0,0,-4);
  smom_comp[1][20]= FourMom(0,0,2,4);
  smom_comp[1][21]= FourMom(0,0,-2,4);
  smom_comp[1][22]= FourMom(0,0,2,-4);
  smom_comp[1][23]= FourMom(0,0,-2,-4);
  
  smom_comp[2].size(80);

  fourierprop_arg.smom[2]=FourMom(1,1,2,4);

  smom_comp[2][0] = FourMom(1,1,2,4);
  smom_comp[2][1] = FourMom(-1,1,2,4);
  smom_comp[2][2] = FourMom(1,-1,2,4);
  smom_comp[2][3] = FourMom(1,1,-2,4);
  smom_comp[2][4] = FourMom(1,1,2,-4);
  smom_comp[2][5] = FourMom(-1,-1,2,4);
  smom_comp[2][6] = FourMom(-1,1,-2,4);
  smom_comp[2][7] = FourMom(-1,1,2,-4);
  smom_comp[2][8] = FourMom(1,-1,-2,4);
  smom_comp[2][9] = FourMom(1,-1,2,-4);
  smom_comp[2][10]= FourMom(1,1,-2,-4);
  smom_comp[2][11]= FourMom(-1,-1,-2,4);
  smom_comp[2][12]= FourMom(1,-1,-2,-4);
  smom_comp[2][13]= FourMom(-1,-1,2,-4);
  smom_comp[2][14]= FourMom(-1,1,-2,-4);
  smom_comp[2][15]= FourMom(-1,-1,-2,-4);
 
  smom_comp[2][16]= FourMom(1,2,1,4);
  smom_comp[2][17]= FourMom(-1,2,1,4);
  smom_comp[2][18]= FourMom(1,-2,1,4);
  smom_comp[2][19]= FourMom(1,2,-1,4);
  smom_comp[2][20]= FourMom(1,2,1,-4);
  smom_comp[2][21]= FourMom(-1,-2,1,4);
  smom_comp[2][22]= FourMom(-1,2,-1,4);
  smom_comp[2][23]= FourMom(-1,2,1,-4);
  smom_comp[2][24]= FourMom(1,-2,-1,4);
  smom_comp[2][25]= FourMom(1,-2,1,-4);
  smom_comp[2][26]= FourMom(1,2,-1,-4);
  smom_comp[2][27]= FourMom(-1,-2,-1,4);
  smom_comp[2][28]= FourMom(1,-2,-1,-4);
  smom_comp[2][29]= FourMom(-1,-2,1,-4);
  smom_comp[2][30]= FourMom(-1,2,-1,-4);
  smom_comp[2][31]= FourMom(-1,-2,-1,-4);
 
  smom_comp[2][32]= FourMom(2,1,1,4);
  smom_comp[2][33]= FourMom(-2,1,1,4);
  smom_comp[2][34]= FourMom(2,-1,1,4);
  smom_comp[2][35]= FourMom(2,1,-1,4);
  smom_comp[2][36]= FourMom(2,1,1,-4);
  smom_comp[2][37]= FourMom(-2,-1,1,4);
  smom_comp[2][38]= FourMom(-2,1,-1,4);
  smom_comp[2][39]= FourMom(-2,1,1,-4);
  smom_comp[2][40]= FourMom(2,-1,-1,4);
  smom_comp[2][41]= FourMom(2,-1,1,-4);
  smom_comp[2][42]= FourMom(2,1,-1,-4);
  smom_comp[2][43]= FourMom(-2,-1,-1,4);
  smom_comp[2][44]= FourMom(2,-1,-1,-4);
  smom_comp[2][45]= FourMom(-2,-1,1,-4);
  smom_comp[2][46]= FourMom(-2,1,-1,-4);
  smom_comp[2][47]= FourMom(-2,-1,-1,-4);

  smom_comp[2][48]=FourMom(2,1,2,2);
  smom_comp[2][49]=FourMom(-2,1,2,2);
  smom_comp[2][50]=FourMom(2,-1,2,2);
  smom_comp[2][51]=FourMom(2,1,-2,2);
  smom_comp[2][52]=FourMom(2,1,2,-2);
  smom_comp[2][53]=FourMom(-2,-1,2,2);
  smom_comp[2][54]=FourMom(-2,1,-2,2);
  smom_comp[2][55]=FourMom(-2,1,2,-2);
  smom_comp[2][56]=FourMom(2,-1,-2,2);
  smom_comp[2][57]=FourMom(2,-1,2,-2);
  smom_comp[2][58]=FourMom(2,1,-2,-2);
  smom_comp[2][59]=FourMom(-2,-1,-2,2);
  smom_comp[2][60]=FourMom(2,-1,-2,-2);
  smom_comp[2][61]=FourMom(-2,-1,2,-2);
  smom_comp[2][62]=FourMom(-2,1,-2,-2);
  smom_comp[2][63]=FourMom(-2,-1,-2,-2);

  smom_comp[2][64]=FourMom(1,2,2,2);
  smom_comp[2][65]=FourMom(-1,2,2,2);
  smom_comp[2][66]=FourMom(1,-2,2,2);
  smom_comp[2][67]=FourMom(1,2,-2,2);
  smom_comp[2][68]=FourMom(1,2,2,-2);
  smom_comp[2][69]=FourMom(-1,-2,2,2);
  smom_comp[2][70]=FourMom(-1,2,-2,2);
  smom_comp[2][71]=FourMom(-1,2,2,-2);
  smom_comp[2][72]=FourMom(1,-2,-2,2);
  smom_comp[2][73]=FourMom(1,-2,2,-2);
  smom_comp[2][74]=FourMom(1,2,-2,-2);
  smom_comp[2][75]=FourMom(-1,-2,-2,2);
  smom_comp[2][76]=FourMom(1,-2,-2,-2);
  smom_comp[2][77]=FourMom(-1,-2,2,-2);
  smom_comp[2][78]=FourMom(-1,2,-2,-2);
  smom_comp[2][79]=FourMom(-1,-2,-2,-2);

  fourierprop_arg.results2 = (void*)"gffprop2.dat";

  int num_mass = 5;
  Float *mass = (Float*) smalloc ( num_mass * sizeof(Float) );
  
  mass[0] = 0.01;
  mass[1] = 0.02;
  mass[2] = 0.03;
  mass[3] = 0.04;
  mass[4] = 0.05;

  int   *disc = (int*) smalloc ( num_mass * sizeof(int) );
  
  disc[0] = 0;
  disc[1] = 0;
  disc[2] = 0;
  disc[3] = 1;
  disc[4] = 0;

    
  //--------------------
  // fourier transform
  //--------------------

  printf("**********************\nFourier Transform\n********************\n");


  {
    common_arg.results = (void*)"info.dat";
    int i;
    for( i=0; i< num_mass; i++)
      {
        fourierprop_arg.cg.mass = mass[i];
	
	AlgFourierPropDis fprop(lattice,&common_arg,&fourierprop_arg);
	fprop.set_discon(disc[i]==1);
	fprop.run();

      } //mass
  }

  sfree((void*)do_arg.start_conf_load_addr);
  
  sfree(smom_comp);  
  sfree(mass);
  sfree(disc);

  //----------------------------------------------------------------
  // Preserve final random stream
  //----------------------------------------------------------------

  /* disabled random stream for debugging -- Shu Li 12/19/05

  if ( rnd_stream == 1 ) 
    {
      RandomGenerator rand ;
      rand.StoreSeeds(rand_buffer) ;
    }

  */
 
}








