/*  The structure type QPropWArg holds the parameters specific to
    QPropW class   */


enum SourceType {
  POINT      = 0 ,  
  VOLUME     = 1 , 
  WALL       = 2 ,
  BOX        = 3 ,
  RANDVOLUME = 4 , 
  RANDWALL   = 5 , 
  RANDSLAB   = 6 , 
  MESSEQ     = 7 , 
  PROT_U_SEQ = 8 , 
  PROT_D_SEQ = 9 ,
  UNDEF      = 10,
  DERIV      = 11,
  GAUSS_GAUGE_INV = 12,
  EXP             = 13,  // for exponential smearing
  SUM_MOM         = 14,
  FL_FUNC         = 15,
  MOM         = 16,
  BOX_4D     = 17
};

enum RandomType {
  GAUSS  = 0, 
  UONE   = 1, 
  ZTWO   = 2, 
  NORAND = 3 } ;

class QPropWArg {

  //! CG arguments for quark propagator inversion
  CgArg cg;		

  //! filename from which propagator may be loaded
  string file<>;

  //! source location
  int x;
  int y;
  int z;
  int t;


  //! Gauge Fixing flags
  int gauge_fix_src;
  int gauge_fix_snk;

  //! should midpoint propagator be stored?
  int store_midprop;
  //! should propagator be saved to disk?
  int save_prop;
  //! should 5d propagator be saved to disk(1) or memory(2)?
  int save_ls_prop;
  //!  should (1+gamma_t)/2 projected sources be used? (good far baryons)
  int do_half_fermion;
  SourceType SeqSmearSink ;

  //! header information of qio write prop
  string ensemble_label<>;
  string ensemble_id<>;
  int seqNum;
  int StartSrcSpin;
  int EndSrcSpin;
  int StartSrcColor;
  int EndSrcColor;
  memfun QPropWArg();
};

class QPropWGFArg {
  //! Gauge Fixing flags
  int gauge_fix_src;
};

class QPropWPointArg {
  int x;
  int y;
  int z;
  memfun QPropWPointArg();
};

class QPropWBoxArg {
  // box source size
  int box_start;
  int box_end;
  int use_xyz_offset; // if we use QPropW.{x,y,z} for the offset of the box 
  memfun QPropWBoxArg();
};

// 4D box source
class QPropW4DBoxArg {
    // These are NOT start/end pairs!
    int box_start[4];
    int box_size[4];
    // momentum, in units of 2*Pi/L, where L is the global size of the
    // lattice in the concerned direction.
    Float mom[4];
    memfun QPropW4DBoxArg();
};

class QPropWRandArg {
  //! random number generator type
  RandomType rng;
  //! random number seed
  int seed;
  memfun QPropWRandArg();
};

class QPropWSlabArg {
  QPropWRandArg rand_arg;
  //! width of slab for random source
  int slab_width;
  memfun QPropWSlabArg();
};

class QPropWExpArg {

  // for exponential smeared source
  // phi(|x-y|) = exp_A*exp(-exp_B*|x-y|) for |x-y| <= exp_C
  //            = 0                       for |x-y| >  exp_C
  Float exp_A;
  Float exp_B;
  Float exp_C;
  memfun QPropWExpArg();
};

class QPropWGaussArg{
  //Gaussian source params
    int gauss_N ;
    Float gauss_W ;
  //Multi Gaussian source params
    int nt;
    int mt[5];
  // smearing scheme for links in Gaussian smearing by YA
    GaussianKernelLinkSmearType gauss_link_smear_type;
    int gauss_link_smear_N;
    Float gauss_link_smear_coeff;
  memfun QPropWGaussArg();
};
