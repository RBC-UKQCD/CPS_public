#include<config.h>
#include<cstdlib>
#ifdef USE_C11_RNG
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
CPS_START_NAMESPACE
/*!\file
  \brief   Methods for the Random Number Generator classes.

*/
//--------------------------------------------------------------------
//
//
//--------------------------------------------------------------------
//---------------------------------------------------------------
//  This is the routine ran3 from Numerical Recipes in C 
//---------------------------------------------------------------
  CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <util/random.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/latrngio.h>
#include <util/data_shift.h>
#include <util/time_cps.h>
#include <comms/glb.h>
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE 

typedef long int SeedInt;

static const int OFFSET = 23;

int RandomGenerator::MBIG = 1000000000;
IFloat RandomGenerator::FAC = 1.0E-09;	// 1.0/MBIG
const int RandomGenerator::state_size;


uint32_t BOOTSTRAP_MAX = 1000000000;

#if (defined USE_C11_MT)
#define RNG_SEED_RANLUX
#elif  (defined USE_C11_RANLUX)
#define RNG_SEED_RANLUX
#else //sitmo
//#define RNG_SEED_RANLUX
#define RNG_SEED_SKIP
#endif
#undef RNG_WARMUP
#undef BOOTSTRAP
const uint64_t RNG_SKIP = 1099511627776; // 2^40
const uint64_t MAX_RNG = 16777216; // 2^24
//const uint64_t RNG_SKIP = 23;

#if (!defined RNG_SEED_RANLUX)  && (!defined RNG_SEED_SKIP )
#error RNG_SEED_RANLUX or RNG_SEED_SKIP should be defined
#endif

//---------------------------------------------------------------


LatRanGen LRG;
const int LatRanGen::default_seed;
//---------------------------------------------------------
LatRanGen::LatRanGen ()
:urand (0., 1.), grand (0., 1.)
{
  cname = "LatRanGen";
  char *fname = "LatRanGen()";
//  printf("%s::%s Entered\n",cname,fname);
  is_initialized = 0;
  UseParIO = 1;
  io_good = 1;
}

LatRanGen::~LatRanGen ()
{
  char *fname = "~LatRanGen()";
//  printf("%s::%s Entered\n",cname,fname);
  if (is_initialized) {
  }
}


/*! Seeds the RNGs according to the method defined in GJP */
//---------------------------------------------------------
void LatRanGen::Initialize ()
{
  if (is_initialized == 1)
    return;

  const char *fname = "Initialize";
  VRB.Func (cname, fname);

//  stringstream cpsran_dump;
//  cpsran_dump.open("cpsran_dump");

  if (0 != GJP.VolNodeSites () % 16)
    ERR.General (cname, fname,
		 "Must have a multiple of 2^4 lattice sites per node.");

  n_rgen_4d = GJP.VolNodeSites () / 16;


  int default_concur = 0;
#if TARGET==BGQ
  default_concur = 1;
#endif


  is_initialized = 1;

//  cpsran = new CPS_RNG[n_rgen_4d];
  cpsran.resize (n_rgen_4d);
//printf("cpsran\n");
// if(!cpsran) ERR.Pointer(cname, fname, "cpsran"); 

  // Lower bounds of global lattice coordinates on this node
  int x_o[5];
  x_o[0] = GJP.XnodeCoor () * GJP.XnodeSites ();
  x_o[1] = GJP.YnodeCoor () * GJP.YnodeSites ();
  x_o[2] = GJP.ZnodeCoor () * GJP.ZnodeSites ();
  x_o[3] = GJP.TnodeCoor () * GJP.TnodeSites ();
  x_o[4] = GJP.SnodeCoor () * GJP.SnodeSites ();

  // Upper bounds of global lattice coordinates on this node
  int x_f[5];
  x_f[0] = x_o[0] + GJP.XnodeSites () - 1;
  x_f[1] = x_o[1] + GJP.YnodeSites () - 1;
  x_f[2] = x_o[2] + GJP.ZnodeSites () - 1;
  x_f[3] = x_o[3] + GJP.TnodeSites () - 1;
  x_f[4] = x_o[4] + GJP.SnodeSites () - 1;


  hx[0] = GJP.XnodeSites () / 2;
  hx[1] = GJP.YnodeSites () / 2;
  hx[2] = GJP.ZnodeSites () / 2;
  hx[3] = GJP.TnodeSites () / 2;
  hx[4] = GJP.SnodeSites () / 2;

  can[0] = 2;
  can[1] = can[0] * GJP.XnodeSites ();
  can[2] = can[1] * GJP.YnodeSites ();
  can[3] = can[2] * GJP.ZnodeSites ();
  can[4] = can[3] * GJP.TnodeSites ();

  int vx[5];
  vx[0] = hx[0] * GJP.Xnodes ();
  vx[1] = hx[1] * GJP.Ynodes ();
  vx[2] = hx[2] * GJP.Znodes ();
  vx[3] = hx[3] * GJP.Tnodes ();
  vx[4] = hx[4] * GJP.Snodes ();

  int index, index_4d;
  index = index_4d = 0;

#ifdef UNIFORM_SEED_TESTING
  if (1) {
#else
  if (GJP.StartSeedKind () == START_SEED_FIXED_UNIFORM) {
#endif
    int x[5];
    SeedInt start_seed = default_seed;
#ifdef RNG_SEED_RANLUX
    VRB.Result(cname,fname,"PRNG seed generated ranlux48\n");
    std::ranlux48 ran_std (start_seed);
#endif


    for (x[3] = 0; x[3] < GJP.TnodeSites (); x[3] += 2) {
      for (x[2] = 0; x[2] < GJP.ZnodeSites (); x[2] += 2) {
	for (x[1] = 0; x[1] < GJP.YnodeSites (); x[1] += 2) {
	  for (x[0] = 0; x[0] < GJP.XnodeSites (); x[0] += 2) {

#ifdef RNG_SEED_RANLUX
	    start_seed = ran_std ();
	    VRB.RNGSeed (cname, fname, "index_4d=%d start_seed=%ld\n", index_4d,
			 start_seed);
	    cpsran[index_4d++].seed (start_seed);
#else
	    cpsran[index_4d].seed (start_seed);
	    uint64_t nskip = RNG_SKIP * (uint64_t) index_4d;
	    VRB.RNGSeed (cname, fname, "index_4d=%d nskip=%ld\n", index_4d,
			 nskip);
	    cpsran[index_4d].discard (nskip);
	    index_4d++;
#endif

	  }
	}
      }
    }
    return;
  }


  // Sort out what the seed should be depending on the GJP.StartSeedKind()

  int start_seed, base_seed;
  uint32_t start_seed_4d, base_seed_4d;

  switch (GJP.StartSeedKind ()) {
  case START_SEED_FILE:
#if 1
    if (!LatRanGen::Read (GJP.StartSeedFilename (), default_concur)) {
      ERR.General (cname, fname, "Reading file %s", GJP.StartSeedFilename ());
    }
#else
    ERR.General (cname, fname, "START_SEED_FILE not implemented yet\n");
#endif
    return;
    break;
  case START_SEED_INPUT_UNIFORM:
  case START_SEED_INPUT_NODE:
  case START_SEED_INPUT:
    base_seed = GJP.StartSeedValue ();
    break;
  case START_SEED_UNIFORM:
  case START_SEED:
    base_seed = SeedS ();
    break;
  default:
    base_seed = default_seed;
  }

#ifdef RNG_SEED_RANLUX
  VRB.Result(cname,fname,"PRNG seed generated ranlux48\n");
  std::ranlux48 ran_std (base_seed);
// warming up rand48 (maybe not necessary?)
  for (int tmp_i = 0; tmp_i < 100; tmp_i++) {
    SeedInt temp = ran_std ();
//      printf("temp=%d\n",temp);
  };
#endif

#ifdef PARALLEL
  if (GJP.StartSeedKind () == START_SEED_INPUT_NODE) {
    int node =
      GJP.XnodeCoor () + GJP.Xnodes () * (GJP.YnodeCoor () +
					  GJP.Ynodes () * (GJP.ZnodeCoor () +
							   GJP.Znodes () *
							   (GJP.TnodeCoor () +
							    GJP.Tnodes () *
							    GJP.SnodeCoor ())));
    base_seed = base_seed + OFFSET * node;
  }
#endif

  // Seed the hypercube RNGs

  int x[5];
//  int index, index_4d;
  index = index_4d = 0;
  uint32_t rng_offset = 0;
#ifdef RNG_SEED_RANLUX
  std::vector < SeedInt > SeedList;
  std::vector < SeedInt >::iterator SeedI;
  uint32_t rng_count = 0;
#endif
#ifdef RNG_SEED_SKIP
#ifdef USE_C11_SITMO
  VRB.Result(cname,fname,"RNG generated by seed and set_counter\n");
#else
  VRB.Result(cname,fname,"RNG generated by seed and skip=%ld\n",RNG_SKIP);
#endif
#endif

  for (x[3] = x_o[3]; x[3] <= x_f[3]; x[3] += 2) {
    for (x[2] = x_o[2]; x[2] <= x_f[2]; x[2] += 2) {
      for (x[1] = x_o[1]; x[1] <= x_f[1]; x[1] += 2) {
	for (x[0] = x_o[0]; x[0] <= x_f[0]; x[0] += 2) {

	  start_seed = start_seed_4d = base_seed;

	  if (GJP.StartSeedKind () == START_SEED ||
	      GJP.StartSeedKind () == START_SEED_INPUT ||
	      GJP.StartSeedKind () == START_SEED_FIXED) {
	    rng_offset = (x[0] / 2 + vx[0] *
			  (x[1] / 2 + vx[1] * (x[2] / 2 + vx[2] * (x[3] / 2))));
#ifdef RNG_SEED_RANLUX
	    while (rng_count < rng_offset) {
	      long int temp = ran_std ();
	      SeedI = std::find (SeedList.begin (), SeedList.end (), temp);
	      if (SeedI != SeedList.end ())
		ERR.General (cname, fname,
			     "index_4d=%d start_seed=%ld Seed Collision(Birthday problem) \n",
			     index_4d, temp);
	      else
		SeedList.push_back (temp);
	      VRB.RNGSeed (cname, fname, "rng_count=%d start_seed=%ld\n",
			   rng_count, temp);
	      rng_count++;
	    }
	    SeedInt temp = ran_std ();
	    SeedI= find (SeedList.begin (), SeedList.end (), temp);
	    if (SeedI != SeedList.end ())
	      ERR.General (cname, fname,
			   "index_4d=%d start_seed=%ld Seed Collision(Birthday problem) \n",
			   index_4d, temp);
	    else
	      SeedList.push_back (temp);
	    start_seed_4d = temp;
	    VRB.RNGSeed (cname, fname, "index_4d=%d start_seed=%ld\n", index_4d,
			 start_seed_4d);
	    rng_count++;
#endif
	  cpsran[index_4d].seed (start_seed_4d);
#ifdef RNG_SEED_SKIP
#ifdef USE_C11_SITMO
//#if 0 
	  cpsran[index_4d].set_counter (0, 0, (uint64_t) rng_offset, 0, 0);
#else
	  assert(index_4d < MAX_RNG); // can be changed to do it sequentially, circumventing the limit
	  uint64_t nskip = RNG_SKIP * (uint64_t) rng_offset;
	  cpsran[index_4d].discard (nskip);
#endif
#endif
	  }
#if 0
	  {
//      volatile int temp = ((int) (dclock()/1000.0))%1000;
	    volatile int temp = rand () % 1000;
	    start_seed += temp;
	    printf ("%g temp start_seed = %d %d", dclock (), temp, start_seed);
	    temp = rand () % 1000;
	    start_seed_4d += temp;
	    printf ("%g temp start_seed_4d = %d %d\n", dclock (), temp,
		    start_seed_4d);
	  }
#endif

//              VRB.Result(cname,fname,"index_4d=%d rng_count=%d start_seed= %d\n",index_4d,rng_count,start_seed_4d);
//      printf("(%d %d %d %d): index_4d=%d rng_count=%d start_seed_4d=%u\n",x[0],x[1],x[2],x[3],index_4d,rng_count,start_seed_4d);
//      std::cout << "cpsran["<<index_4d<<"]:\n"<<cpsran[index_4d]<<std::endl;

#ifdef RNG_WARMUP
	  {
	    int n_warm = ugran_4d[index_4d].Urand (100, 0);
	    VRB.RNGSeed (cname, fname, "index_4d=%d n_warm=%d\n", index_4d,
			n_warm);
	    while (n_warm > 0) {
	      int temp = ugran_4d[index_4d].Urand (100, 0);
	      n_warm--;
	    }
	  }
#endif
#ifdef BOOTSTRAP
	  {
	    std::uniform_int_distribution <> uniform_dist (0, BOOTSTRAP_MAX);
	    int new_seed = uniform_dist (cpsran[index_4d]);
	    VRB.RNGSeed (cname, fname,
			"index_4d=%d start_seed_4d=%d new_seed=%d\n", index_4d,
			start_seed_4d, new_seed);
	    cpsran[index_4d].seed (new_seed);
	  }
#endif
	  if (0) {
	    std::stringstream cpsran_dump;
	    cpsran_dump << cpsran[index_4d];
	    if (!UniqueID ()) {
	      std::cout << "cpsran[" << index_4d << "]" << std::endl;
	      std::cout << cpsran_dump.str () << std::endl;
	      cpsran_dump.seekg (0, cpsran_dump.beg);
	    }
#if 0
	    for (int i_dump = 0; i_dump < 1000 && !cpsran_dump.eof (); i_dump++) {
	      RNGSTATE dump, dump2;
#if 1
	      cpsran_dump >> dump;
#else
	      char temp_num[50];
	      cpsran_dump.get (temp_num, 50, ' ');
	      sscanf (temp_num, "%d", &dump);
	      cpsran_dump.get (temp_num, 50, ' ');
	      sscanf (temp_num, "%d", &dump2);
#endif
	      std::cout << i_dump << " " << dump << std::endl;
	    }
#endif
	  }
	  index_4d++;
	}
      }
    }
  }
}

//---------------------------------------------------------
/*!
  \pre A RNG must be assigned using ::AssignGenerator.
  \return A uniform random number for this  hypercube.
*/
//----------------------------------------------------------------------
IFloat LatRanGen::Urand (FermionFieldDimension frm_dim)
{
  char *fname = "Urand(FermionFieldDimension)";
  return (urand_lo + (urand_hi - urand_lo) * urand (cpsran[rgen_pos_4d]));
}

IFloat LatRanGen::Urand (Float hi, Float lo, FermionFieldDimension frm_dim)
{
  char *fname = "Urand(FermionFieldDimension)";
  return (lo + (hi - lo) * urand (cpsran[rgen_pos_4d]));
}

//---------------------------------------------------------
/*!
  \pre A RNG must be assigned using ::AssignGenerator.
  \return A gaussian random number for this hypercube.
*/
//----------------------------------------------------------------------
IFloat LatRanGen::Grand (FermionFieldDimension frm_dim)
{
  char *fname = "Grand(FermionFieldDimension)";
//  printf("LatRanGen::Grand():%d\n",rgen_pos);
  return grand_mean + sqrt (grand_sigma) * grand (cpsran[rgen_pos_4d]);
}


//---------------------------------------------------------
/*!
  The parameters are set for the RNGs on all hypercubes.
  \param high the upper bound of the distribution range
  \param low the lower bound of the distribution range
*/
//----------------------------------------------------------------------
void LatRanGen::SetInterval (IFloat high, IFloat low)
{
//  ugran[0].SetInterval(high,low);
//    UniformRandomGenerator::SetInterval(high,low);
  urand_lo = low;
  urand_hi = high;
}

//---------------------------------------------------------
/*!
  The parameters are set for the RNGs on all hypercubes. The mean is zero.
  \param sigma the variance of the gaussian distribution.
*/
//----------------------------------------------------------------------
void LatRanGen::SetSigma (IFloat sigma)
{
//  ugran[0].SetSigma(sigma);
//    GaussianRandomGenerator::SetSigma(sigma);
  grand_sigma = sigma;
}

//---------------------------------------------------------
/*!
  For a given lattice site, this identifies and assigns the corresponding
  hypercube RNG.
  \param x The x coordinate of the lattice site.
  \param y The y coordinate of the lattice site.
  \param z The z coordinate of the lattice site.
  \param t The t coordinate of the lattice site.
  \param s The s coordinate of the lattice site.  
  \post  Subsequent calls to \e e.g. Urand will return results from this
  particular hypercubic RNG.  
 */
//----------------------------------------------------------------------
void LatRanGen::AssignGenerator (int x, int y, int z, int t, int s,
				 const int &field_idx)
{
  x = x % GJP.XnodeSites ();
  y = y % GJP.YnodeSites ();
  z = z % GJP.ZnodeSites ();
  t = t % GJP.TnodeSites ();

  if (GJP.SnodeSites () < 2)
    s = 0;
  else
    s = s % GJP.SnodeSites ();

  x /= 2;
  y /= 2;
  z /= 2;
  t /= 2;
  s /= 2;

  rgen_pos_4d = x + hx[0] * (y + hx[1] * (z + hx[2] * t));
}

//----------------------------------------------------------------------
/*!
  For a given lattice site, this identifies and
  assigns the corresponding hypercube RNG.
  \param x Array holding the lattice site coordinates in order X, Y, Z, T.
  \post  Subsequent calls to \e e.g. Urand will return results from this
  particular hypercubic RNG.
*/
//----------------------------------------------------------------------
void LatRanGen::AssignGenerator (const int *x, const int &field_idx)
{
  AssignGenerator (x[0], x[1], x[2], x[3], 0);
}

//---------------------------------------------------------
/*!
  For a given lattice site, this identifies and
  assigns the corresponding hypercube RNG.
  \param i The canonical index of the lattice site.
  \post  Subsequent calls to \e e.g. Urand will return results from this
  particular hypercubic RNG.
*/
//----------------------------------------------------------------------
void LatRanGen::AssignGenerator (int i, const int &field_idx)
{
  char *fname = "AssignGenerator(i)";
  int x = (i / can[0]) % hx[0];
  int y = (i / can[1]) % hx[1];
  int z = (i / can[2]) % hx[2];
  int t = (i / can[3]) % hx[3];
  int s;
  if (GJP.SnodeSites () < 2)
    s = 0;
  else
    s = (i / can[4]) % hx[4];
  rgen_pos_4d = x + hx[0] * (y + hx[1] * (z + hx[2] * t));
  if (rgen_pos_4d >= n_rgen_4d)
    ERR.General (cname, fname, "rgen_pos_4d(%d)>=n_rgen_4d(%d)", rgen_pos_4d,
		 n_rgen_4d);
//  Fprintf(stdout,"i=%d x = %d %d %d %d %d rgen_pos=%d ",i,x,y,z,t,s,rgen_pos);
}


//--------------------------------------------------------------
/*!
  \return A uniform random number; one per node, the same on each node.
*/
// Lrand will return the same random number on every node by
// performing a global sum over all 2^4 hypercubes, and taking the
// average value
//--------------------------------------------------------------
#if 0
IFloat LatRanGen::Lrand ()
{
  Float cntr = 0.0;
  for (int i = 0; i < n_rgen; i++)
    cntr += (Float) ugran[i].Urand ();

  Float divisor =
    (Float) (n_rgen * GJP.Xnodes () * GJP.Ynodes () * GJP.Znodes ()
	     * GJP.Tnodes ());
  if (GJP.Snodes () > 1)
    divisor *= (Float) GJP.Snodes ();
  cntr /= divisor;
  glb_sum_five (&cntr);
  return (IFloat) cntr;

}
#endif


/*!
  \return The number of unsigned ints that comprise the RNG state vector.
*/
int LatRanGen::StateSize () const
{

  return state_size + 2;

}

/*!
  \return The total number of RNG states on the local lattice.
*/
int LatRanGen::NStates (FermionFieldDimension frm_dim) const
{
  return n_rgen_4d;
}

#ifndef USE_C11_RNG
/*!
  \pre A RNG must be assigned using ::AssignGenerator.
  \a state must point to an array of unsigned ints (at least) as long
  as the value returned by ::StateSize.

  \param state The state to be copied from the RNG.
  \param frm_dim If FIVE_D, the default, refers to the normal RNG. If FOUR_D
  then the special RNG four gauge fields with domain-wall fermions is used.  
*/
void LatRanGen::GetState (RNGSTATE * state, FermionFieldDimension frm_dim) const
{
  ugran_4d[rgen_pos_4d].StoreSeeds (state);
}

/*!
  \pre A RNG must be assigned using ::AssignGenerator.
  \pre \a s must be an array with length given by ::StateSize.
  \param s The state to assigned to the RNG on an assigned site.
*/
void LatRanGen::SetState (const RNGSTATE * s, FermionFieldDimension frm_dim)
{

  // ugran[rgen_pos].RestoreSeeds(s);
  ugran_4d[rgen_pos_4d].RestoreSeeds (s);

}

/*!
  \pre \a s must point to a 2-d array with lengths given by ::NStates and 
  ::StateSize.
  \param s The state to assigned to the RNGs on the entire local lattice.
  \param frm_dim If FIVE_D, the default, refers to the normal RNG. If FOUR_D
  then the special RNG four gauge fields with domain-wall fermions is used.  
*/
void LatRanGen::SetStates (unsigned int **s, FermionFieldDimension frm_dim)
{
  if (frm_dim == FIVE_D)
    for (int h = 0; h < n_rgen; h++)
      ugran[h].RestoreSeeds (s[h]);
  else
    for (int h = 0; h < n_rgen_4d; h++)
      ugran_4d[h].RestoreSeeds (s[h]);

}

/*!
  \pre \a s must point to a 2-d array with lengths given by ::NStates and
  ::StateSize.
  \param s The state to be copied from the RNGs on the entire local lattice. 
  \param frm_dim If FIVE_D, the default, refers to the normal RNG. If FOUR_D
  then the special RNG four gauge fields with domain-wall fermions is used.  
*/
void LatRanGen::GetStates (unsigned int **s, FermionFieldDimension frm_dim) const
{

  if (frm_dim == FIVE_D)
    for (int h = 0; h < n_rgen; h++)
      ugran[h].StoreSeeds (s[h]);
  else
    for (int h = 0; h < n_rgen_4d; h++)
      ugran_4d[h].StoreSeeds (s[h]);

}

void LatRanGen::Shift ()
{
  GDS.Shift (ugran_4d, n_rgen_4d * sizeof (UGrandomGenerator));
}
#endif

void LatRanGen::GetAllStates (RNGSTATE * dump)
{
  for (int h = 0; h < NStates (); h++) {
    std::stringstream ss_dump;
    ss_dump << cpsran[h];
    if (!UniqueID () && h == 0) {
      std::cout << "GetAllState::cpsran[0]= " << ss_dump.str () << std::endl;
      ss_dump.seekg (0);
    }
    for (int i = 0; i < StateSize (); i++) {
      RNGSTATE temp;
      ss_dump >> temp;
//              if (!UniqueID()&& h==0)
//      std::cout <<"GetAllState::cpsran[0]= "<< i <<" "<<temp <<std::endl;
      dump[h * StateSize () + i] = temp;
    }
  }
}

void LatRanGen::SetAllStates (RNGSTATE * dump)
{
  for (int h = 0; h < NStates (); h++) {
    std::stringstream ss_dump;
    for (int i = 0; i < StateSize (); i++) {
      ss_dump << dump[h * StateSize () + i] << " ";
    }
    if (!UniqueID () && h == 0) {
      std::cout << "SetAllState::cpsran[" << h << "]= " << ss_dump.
	str () << std::endl;
      ss_dump.seekg (0);
    }
    ss_dump >> cpsran[h];
  }
}

bool LatRanGen::Read (const char *filename, int concur_io_number)
{
  io_good = false;
  QioArg rd_arg (filename, concur_io_number);
  LatRngRead rd;
  RNGSTATE *cpsran_dump = new RNGSTATE[StateSize () * NStates ()];
  if (parIO ())
    rd.setParallel ();
  else
    rd.setSerial ();
  if (do_log)
    rd.setLogDir (log_dir);
  rd.read (cpsran_dump, rd_arg);
  LRG.SetAllStates (cpsran_dump);
  delete[]cpsran_dump;
  return (io_good = rd.good ());
}


bool LatRanGen::Write (const char *filename, int concur_io_number)
{
  io_good = false;
  IntConv intconv;
  INT_FORMAT intFormat = intconv.testHostFormat (sizeof (RNGSTATE));
  QioArg wt_arg (filename, concur_io_number, intFormat);
  VRB.Result (cname, "Write()", "concur_io_number=%d %s\n", concur_io_number,
	      intconv.name (wt_arg.FileIntFormat));
  RNGSTATE *cpsran_dump = new RNGSTATE[StateSize () * NStates ()];
  LRG.GetAllStates (cpsran_dump);
  LatRngWrite wt;

  if (parIO ())
    wt.setParallel ();
  else
    wt.setSerial ();
  if (do_log)
    wt.setLogDir (log_dir);
  wt.write (cpsran_dump, wt_arg);
  delete[]cpsran_dump;
  return (io_good = wt.good ());
}

CPS_END_NAMESPACE
#endif
