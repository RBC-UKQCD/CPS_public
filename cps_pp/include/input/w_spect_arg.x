

/*---------------------------------------------------------------------------*/
/* The following specifies output format and filenames*/
/*---------------------------------------------------------------------------*/
/**/
/* BARYON_PAST:   same as in phys_v3.11.1 and older versions, which means */
/*                information partially lost during time slice folding.*/
/*                Print out only T/2 + 1 slices.*/
/*                Fitting ansatz: cosh or sinh*/
/**/
/* BARYON_FOLD:   average [up to a minus sign from bnd cnd for baryons]*/
/*                between timeslice i of up Dirac projection and */
/*                timeslice T - i of down projection.*/
/*                Print out T slices.*/
/*                Fitting ansatz: exp*/
/**/
/* BARYON_RAW:    print out up and down projections seperately.*/
/*                Print out T slices.*/
/*                Fitting ansatz: exp*/
/**/
/**/
/* Data Format in General*/
/*        each row:*/
/*             col 1 - measurement counter*/
/*             col 2 - time slice index*/
/*             col 3 - momentum index [0 for baryons]*/
/*             col 4 - the real part of the hadron correlator*/
/*             col 5 - the real part of the hadron correlator if BARYON_RAW*/
/*        repeat above for all momenta from the same time slice*/
/*             [only relevant to mesons] */
/*        repeat above for T [BARYON_FOLD and BARYON_RAW] or */
/*                         T/2+1 [mesons and BARYON_PAST] slices.*/
/* */
/*---------------------------------------------------------------------------*/



class  WspectOutput 
{
  WbaryonFold      fold;    
  
  char *           cg;                  /* datafile for cg iter & res*/
  char *           cg2;                 /* need and extra one for second quark propagator?*/
  char *           pbp;                 /* datafile for pbp and pbg5p*/
  char *           mid_point;           /* datafile for mid-point sink*/
  char *           a0_p;             /* datafile for <A_0 P>*/

  /* the following is as in phys_v3.11.4.xiaodong, i.e. average over polarisation*/
  char *           a1;
  char *           b1;
  char *           pion;
  char *           pion_prime;
  char *           rho;
  char *           rho_prime;

/* The a0 and a1's must be continguous because of part of the code*/
  char *           a0;                  /* datafile for meson a0*/
  char *           a0_prime;
  char *           a1_x;                /* x, y, z are to be understood*/
  char *           a1_y;                /* as the three directions other*/
  char *           a1_z;                /* than the propagation direction.*/
  char *           b1_x;
  char *           b1_y;
  char *           b1_z;
  char *           rho_x;
  char *           rho_y;
  char *           rho_z;
  char *           rho_x_prime;
  char *           rho_y_prime;
  char *           rho_z_prime;
  char *           nucleon;
  char *           nucleon_prime;
  char *           delta_x;
  char *           delta_y;
  char *           delta_z;
  char *           delta_t;
// char *filetail;/*eg. P, W, B4, JN30E2.50,... */

  memfun WspectOutput() ;
//{ pbp=0; mid_point=0; a0_p=0;} //that jsut sounds wrong, so I won't propagate - petrov.

};  


/*---------------------------------------------------------------------------*/
/* The following specifies how the spectrum is measured.*/
/*---------------------------------------------------------------------------*/
/* One WspectArg needed for each SourceKind*/
/*---------------------------------------------------------------------------*/

class WspectArg {

  /* Propagation direction. [0..3] as [x,y,z,t]*/
  int prop_dir;               

  /* Apply multi-sinking from one calcuation of quark prop. - N/A yet.*/
  /*  int num_sink;             */
  
  /* Apply multi-momenta projection.*/
  int num_mom;
				
  /* Type of source for the calculation of the quark propagator*/
  SourceKind source_kind;       

  /*tail of output file names */
/* char *filetail;eg. P, W, B4, JN30E2.50,... */

  /*--------------------------*/
  /*Paremeter of sources*/
  /*--------------------------*/
  /* Relevant only for BOX_W*/
  /**/
  /* Note: 1. begin[prop_dir] and end[prop_dir] do not matter.*/
  /*       2. In all directions including prop_dir, we check that*/
  /*          0 <= begin[] <= end[] <= total_sites*/
  /*       3. We check that at least in one direction other than prop_dir*/
  /*          0 <= begin[] <  end[] <= total_sites*/
  int src_box_b[4];           
  int src_box_e[4];      


  /* int snk_box_b[4];          - N/A yet. */
  /* int snk_box_e[4];          - N/A yet.*/


  /*added by Thomas and Xiaodong*/
  /*relevant to Jacobi sources*/
  Float g_epsi;
  int g_n;
  int g_center[4]; /*center of the Gaussian source*/
      

  /* rescaling the source vector to get the full*/
  /* range of single precision arithmetic*/
  /* important for heavy quarks where the strong exponential*/
  /* decay of propagators may be out of single precision range*/
  /* default value = 1.0 (no rescaling)*/
  Float rescale_factor;

  /* AOTS: time slices to put source on for quark propagators*/
  /*       time slices are in arithmetic series [start, ...]*/
  /*       num = 0 ==> no inversion will be done.*/
  int aots_num;                    
  int aots_start;                  
  int aots_step;                   


  /*------------------------------*/
  /*Mesurement control flags*/
  /*------------------------------*/
  int baryons_on;                
  int normal_mesons_on;
  int extended_mesons_on;
  int extended_mesonsBE_on; /*extended mesons from local B/E fields*/

  /*related to extended_mesons using derivative operators*/
  /*See enum DEVOpKind for explanation*/
  int extended_mesons_op_groupId; /*use source operator sum only*/
  int extended_mesons_first_dev_on; 
  int extended_mesons_second_sym_dev_on; /*do second derivative src_op*/
  int extended_mesons_second_antisym_dev_on;
  int extended_mesons_second_diag_dev_on;

  /*Fuzzing parameters*/
  /*-------------------------------*/
  int fuzzing_on; /* 1/0(on/off)*/
  int sink_fuzzing_only; /*1(yes): use fuzzed links at sink only. 0(no): use at both sink and source*/
  int fuzzing_level;
  int fuzzing_c_num;/*number of coefficients to run*/
  Float fuzzing_c[MAX_FUZZING_C_NUM]; /*multiplier*/
  int fuzzing_hits; /*cabbobo hits*/


  /*related to extended_mesonsBE*/
  int extended_mesonsBE_op_groupId; 
  int extended_mesonsBE_Elec_on;
  int extended_mesonsBE_Magn_on;

  /*Fuzzing parameters for extended_mesonsBE*/
  int BEfuzzing_on; /* 1/0(on/off)*/
  int BEfuzzing_level;
  int BEfuzzing_c_num;/*number of coefficients to run*/
  Float BEfuzzing_c[MAX_FUZZING_C_NUM]; /*multiplier*/
  int BEfuzzing_hits; /*cabbobo hits*/
memfun WspectArg();
};


