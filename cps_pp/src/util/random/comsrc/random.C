#include<config.h>
#ifndef USE_C11_RNG
CPS_START_NAMESPACE
/*!\file
  \brief   Methods for the Random Number Generator classes.
*/
//---------------------------------------------------------------
//  This is the routine ran3 from Numerical Recipes in C 
//---------------------------------------------------------------
#define FIX_GPBC_RNG_BUG	//CK 2015 fix for GPBC seeding bug
  CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <util/random.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/latrngio.h>
#include <util/data_shift.h>
#include <comms/glb.h>
#include <comms/sysfunc_cps.h>
#include <stdio.h>
#include <string.h>
CPS_START_NAMESPACE
#define BOOTSTRAP		//need to be undef'd for regression testing with Gparity
#define RNG_WARMUP		//need to turn off for non-Gparity regression!
static const int OFFSET = 23;
static const int N_WARMUP = 1000;

int RandomGenerator::MBIG = 1000000000;
IFloat RandomGenerator::FAC = 1.0E-09;	// 1.0/MBIG
const int RandomGenerator::state_size;

//CK

RandomGenerator::RandomGenerator (const RandomGenerator & in):inext (in.inext),
inextp (in.inextp)
{
  memcpy (&ma, &in.ma, state_size * sizeof (int));
}

RandomGenerator & RandomGenerator::operator= (const RandomGenerator & in)
{
  inext = in.inext;
  inextp = in.inextp;
  memcpy (&ma, &in.ma, state_size * sizeof (int));
  return *this;
}

bool RandomGenerator::operator== (const RandomGenerator & in) const
{
  if (inext != in.inext || inextp != in.inextp) {
    printf ("Rand false on ints %d:%d   %d:%d\n", inext, in.inext, inextp,
	    in.inextp);
    return false;
  }

  for (int i = 0; i < state_size; i++) {
    if (ma[i] != in.ma[i]) {
      printf ("Rand false on ma[%d] : %d:%d\n", i, ma[i], in.ma[i]);
      return false;
    }
  }
  return true;
}



void RandomGenerator::store (int *buf)
{
  memcpy (buf, ma, state_size * sizeof (int));
  buf[state_size] = inext;
  buf[state_size + 1] = inextp;
}

//! to load from file
void RandomGenerator::load (int *buf)
{
  char *cname = "RandomGenerator";
  char *fname = "load(int *buf)";

  memcpy (ma, buf, state_size * sizeof (int));
  inext = buf[state_size];
  inextp = buf[state_size + 1];
  if (inext >= state_size || inextp >= state_size)
    ERR.General (cname, fname,
		 "Read error: inext = %d or inextp = %d is greater than state_size = %d\n",
		 inext, inextp);

  //printf("RandomGenerator::load got inext %d, inextp %d\n",inext,inextp);
}

/*!
  This method must be called before the RNG is used
  \param idum The seed.	
 */
void RandomGenerator::Reset (int idum)
{
  int i, k, ii;
  int mk, mj;

  const int MSEED = 161803398;

  //-----------------------------------------------------
  //  Initialize ma[state_size] using the seed idum and the large
  //  number MSEED.
  //-----------------------------------------------------
  mj = MSEED - (idum < 0 ? -idum : idum);
  mj %= MBIG;
  if (mj < 0)
    mj = -mj;			// Added by Roy and Pavlos to protect
  // for the case where idum > MSEED
  ma[54] = mj;


  //-----------------------------------------------------
  //  Initialize the reset of the table in a slightly
  //  random order, with numbers that are not especially
  //  random.
  //-----------------------------------------------------
  mk = 1;
  for (i = 1; i < state_size; i++) {
    ii = (21 * i) % state_size - 1;
    ma[ii] = mk;
    mk = mj - mk;
    if (mk < 0)
      mk += MBIG;
    mj = ma[ii];
  }

  // Randomize them by "warming up the generator"
  for (k = 0; k < 4; k++) {
    for (i = 0; i < state_size; i++) {
      ma[i] -= ma[(i + 31) % state_size];
      if (ma[i] < 0)
	ma[i] += MBIG;
    }
  }

  //-----------------------------------------------------
  //  Prepare indices for our first generated number.
  //  the constant 30 is special.
  //-----------------------------------------------------
  inext = -1;
  inextp = 30;
}




//---------------------------------------------------------------

/*!
  \param to Pointer to a buffer of size at least 57*sizeof(int).
*/

void RandomGenerator::StoreSeeds (unsigned int *to) const
{
  *to++ = (unsigned int) inext;
  *to++ = (unsigned int) inextp;
  for (int i = 0; i < 55; ++i)
    *to++ = ma[i];
}

//CK added for testing
GaussianRandomGenerator::
GaussianRandomGenerator (const GaussianRandomGenerator &
			 in):RandomGenerator (in), iset (in.iset),
gset (in.gset)
{
}

GaussianRandomGenerator & GaussianRandomGenerator::
operator= (const GaussianRandomGenerator & in)
{
  RandomGenerator::operator= (in);
  iset = in.iset;
  gset = in.gset;
  return *this;
}

bool GaussianRandomGenerator::operator== (const GaussianRandomGenerator & in) const
{
  if (iset != in.iset || (iset != 0 && gset != in.gset)) {	//if iset is zero gset is regenerated when a random number is drawn
    printf ("False on Grand ints %d:%d   %d:%d\n", iset, in.iset, gset,
	    in.gset);
    return false;
  }
  return RandomGenerator::operator== (in);
}



UGrandomGenerator::UGrandomGenerator (const UGrandomGenerator & in):GaussianRandomGenerator (in),
UniformRandomGenerator
(in)
{
}

UGrandomGenerator & UGrandomGenerator::operator= (const UGrandomGenerator & in)
{
  GaussianRandomGenerator::operator= (in);
  return *this;
}

bool UGrandomGenerator::operator== (const UGrandomGenerator & in) const
{
  return GaussianRandomGenerator::operator== (in);	//UniformRandomNumberGenerator has no data members
}

void UGrandomGenerator::StoreSeeds (unsigned int *to) const
{
  RandomGenerator::StoreSeeds (to);
}

/*!
  \param from Pointer to a buffer of size at least 57*sizeof(int).
*/

void RandomGenerator::RestoreSeeds (const unsigned int *from)
{
  inext = *from++;
  inextp = *from++;

  for (int i = 0; i < 55; ++i)
    ma[i] = *from++;
}

/*!
  \return The number of unsigned ints that comprise the RNG state vector.
*/
int RandomGenerator::StateSize () const
{

  return state_size + 2;

}

IFloat UniformRandomGenerator::A = -0.5;
IFloat UniformRandomGenerator::B = 0.5;
IFloat GaussianRandomGenerator::sigma2 = 1.0;



LatRanGen LRG;
const int LatRanGen::default_seed;
//---------------------------------------------------------
LatRanGen::LatRanGen ()
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
    delete[]ugran;
    delete[]ugran_4d;
  }
}

//CK added copy constructor, operator= and operator== for testing purposes
LatRanGen::LatRanGen (const LatRanGen & in):n_rgen (in.n_rgen),
rgen_pos (in.rgen_pos), is_initialized (in.is_initialized),
n_rgen_4d (in.n_rgen_4d), rgen_pos_4d (in.rgen_pos_4d), UseParIO (in.UseParIO),
io_good (in.io_good), do_log (in.do_log)
{
  char *fname = "LatRanGen(const LatRanGen &in)";

  if (is_initialized) {
    memcpy (&can, &in.can, 5 * sizeof (int));
    memcpy (&hx, &in.hx, 5 * sizeof (int));
    memcpy (&log_dir, &in.log_dir, 200 * sizeof (char));

    ugran = new UGrandomGenerator[n_rgen];
    if (!ugran)
      ERR.Pointer (cname, fname, "ugran");
    ugran_4d = new UGrandomGenerator[n_rgen_4d];
    if (!ugran_4d)
      ERR.Pointer (cname, fname, "ugran_4d");

    for (int i = 0; i < n_rgen; i++) {
      ugran[i] = in.ugran[i];
    }
    for (int i = 0; i < n_rgen_4d; i++) {
      ugran_4d[i] = in.ugran_4d[i];
    }
    cname = "LatRanGen";
  }
}

LatRanGen & LatRanGen::operator= (const LatRanGen & in)
{
  char *fname = "LatRanGen &LatRanGen::operator=(const LatRanGen &in)";

  if (is_initialized) {
    delete[]ugran;
    delete[]ugran_4d;
  }

  n_rgen = in.n_rgen;
  rgen_pos = in.rgen_pos;
  is_initialized = in.is_initialized;
  n_rgen_4d = in.n_rgen_4d;
  rgen_pos_4d = in.rgen_pos_4d;
  UseParIO = in.UseParIO;
  io_good = in.io_good;
  do_log = in.do_log;

  if (is_initialized) {
    memcpy (&can, &in.can, 5 * sizeof (int));
    memcpy (&hx, &in.hx, 5 * sizeof (int));
    memcpy (&log_dir, &in.log_dir, 200 * sizeof (char));

    ugran = new UGrandomGenerator[n_rgen];
    if (!ugran)
      ERR.Pointer (cname, fname, "ugran");
    ugran_4d = new UGrandomGenerator[n_rgen_4d];
    if (!ugran_4d)
      ERR.Pointer (cname, fname, "ugran_4d");

    for (int i = 0; i < n_rgen; i++) {
      ugran[i] = in.ugran[i];
    }
    for (int i = 0; i < n_rgen_4d; i++) {
      ugran_4d[i] = in.ugran_4d[i];
    }
    cname = "LatRanGen";
  }
}

bool LatRanGen::operator== (const LatRanGen & in) const
{
  if (n_rgen != in.n_rgen || rgen_pos != in.rgen_pos
      || is_initialized != in.is_initialized || n_rgen_4d != in.n_rgen_4d
      || rgen_pos_4d != in.rgen_pos_4d) {
    printf ("False on integers %d:%d  %d:%d  %d:%d  %d:%d  %d:%d\n", n_rgen,
	    in.n_rgen, rgen_pos, in.rgen_pos, is_initialized, in.is_initialized,
	    n_rgen_4d, in.n_rgen_4d, rgen_pos_4d, in.rgen_pos_4d);

    return false;
  }
  //Dont care about IO stuff, only the state of the RNG
  //|| UseParIO != in.UseParIO ||
  //io_good != in.io_good || do_log != in.do_log) return false;
  //if(strcmp(log_dir,in.log_dir) != 0) return false;

  for (int i = 0; i < n_rgen; i++) {
    if (!(ugran[i] == in.ugran[i])) {
      printf ("False on ugran[%d]\n", i);
      return false;
    }
  }
  for (int i = 0; i < n_rgen_4d; i++) {
    if (!(ugran_4d[i] == in.ugran_4d[i])) {
      printf ("False on ugran_4d[%d]\n", i);
      return false;
    }
  }
  return true;
}

//added by CK - turns off is_initialised flag, frees memory and runs initialize again
void LatRanGen::Reinitialize ()
{
  is_initialized = 0;
  delete[]ugran;
  delete[]ugran_4d;
  return Initialize ();
}

/*! Seeds the RNGs according to the method defined in GJP */
//---------------------------------------------------------
void LatRanGen::Initialize ()
{
  if (is_initialized == 1)
    return;

  const char *fname = "Initialize";
  VRB.Func (cname, fname);

  if (0 != GJP.VolNodeSites () % 16)
    ERR.General (cname, fname,
		 "Must have a multiple of 2^4 lattice sites per node.");

  n_rgen = n_rgen_4d = GJP.VolNodeSites () / 16;

  if (GJP.SnodeSites () >= 2)
    n_rgen = GJP.VolNodeSites () * GJP.SnodeSites () / 32;
   VRB.Debug(cname,fname,"n_rgen=%d %d\n",n_rgen,n_rgen_4d);

  //CK: G-parity
  //for 4d RNGs, stack 2 4d lattice volumes on top of each other
  //for 5d RNGs, we have the same idea, only each block is a 2^5 hypercube. We still stack 2 on each s layer.
  if (GJP.Gparity ()) {
    n_rgen *= 2;
    n_rgen_4d *= 2;
  }
  int blocks_per_s_layer;
  if (GJP.SnodeSites () >= 2)
    blocks_per_s_layer = n_rgen / (GJP.SnodeSites () / 2);
  else
    blocks_per_s_layer = 1;

  int default_concur = 0;
#if TARGET==BGQ
  default_concur = 1;
#endif

  //For GPBC
  int nodes_4d = 1;
  for (int i = 0; i < 4; i++)
    nodes_4d *= GJP.Nodes (i);

  int n_rgen_glb = n_rgen * nodes_4d * GJP.Snodes ();	//global number of 5D generators
  int n_rgen_4d_glb = n_rgen_4d * nodes_4d;
   VRB.Debug(cname,fname,"n_rgen_glb=%d %d\n",n_rgen_glb,n_rgen_4d_glb);


  is_initialized = 1;

  ugran = new UGrandomGenerator[n_rgen];
  if (!ugran)
    ERR.Pointer (cname, fname, "ugran");
  ugran_4d = new UGrandomGenerator[n_rgen_4d];
  if (!ugran)
    ERR.Pointer (cname, fname, "ugran_4d");

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

  //G-parity, second stacked set of RNGs, start increment from halfway point of RNG array
  //RNGs for U* field are therefore *NOT THE SAME* as those for the U field. Be aware of this when generating random transforms for the gauge fields
  int stk_index = blocks_per_s_layer / 2;
  int stk_index_4d = n_rgen_4d / 2;

  /*PAB: Implement the Britney and Christina tests correctly */
#ifdef UNIFORM_SEED_TESTING
  if (1) {
#else
  if (GJP.StartSeedKind () == START_SEED_FIXED_UNIFORM) {
#endif
    int x[5];
    int start_seed = default_seed;
    int start_seed_stacked = default_seed + n_rgen / 2 * OFFSET;

    for (x[4] = 0; x[4] < GJP.SnodeSites (); x[4] += 2) {
      for (x[3] = 0; x[3] < GJP.TnodeSites (); x[3] += 2) {
	for (x[2] = 0; x[2] < GJP.ZnodeSites (); x[2] += 2) {
	  for (x[1] = 0; x[1] < GJP.YnodeSites (); x[1] += 2) {
	    for (x[0] = 0; x[0] < GJP.XnodeSites (); x[0] += 2) {

	      if (GJP.Gparity ()) {	//fill in both stacked volumes simultaneously with different seeds
		start_seed += OFFSET;
		start_seed_stacked += OFFSET;
		if (x[4] == 0) {
		  ugran_4d[index_4d++].Reset (start_seed);
		  ugran_4d[stk_index_4d++].Reset (start_seed_stacked);
		}
		ugran[index++].Reset (start_seed);
		ugran[stk_index++].Reset (start_seed_stacked);
	      } else {
		start_seed += OFFSET;
		ugran[index++].Reset (start_seed);
		if (x[4] == 0)
		  ugran_4d[index_4d++].Reset (start_seed);
	      }

	    }
	  }
	}
      }
      if (GJP.Gparity ()) {	//moving onto next s, we need to shift both the stacked and original write index forward one 4d volume's worth
	index += blocks_per_s_layer / 2;
	stk_index += blocks_per_s_layer / 2;
      }
    }
    return;
  }


  // Sort out what the seed should be depending on the GJP.StartSeedKind()

  int start_seed, base_seed;
  int start_seed_4d, base_seed_4d;
  int start_seed_stacked, start_seed_stacked_4d;	//CK: G-parity

  switch (GJP.StartSeedKind ()) {
  case START_SEED_FILE:
    if (!LatRanGen::Read (GJP.StartSeedFilename (), default_concur)) {
      ERR.General (cname, fname, "Reading file %s", GJP.StartSeedFilename ());
    }
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
  UGrandomGenerator rng_seed_4d;
  UGrandomGenerator rng_seed_5d;
  rng_seed_4d.Reset (base_seed);
  rng_seed_5d.Reset (base_seed);
  int rng_count = 0, rng_count_4d = 0;


  int x[5];
  //  int index, index_4d;
  index = index_4d = 0;

  for (x[4] = x_o[4]; x[4] <= x_f[4]; x[4] += 2) {
    for (x[3] = x_o[3]; x[3] <= x_f[3]; x[3] += 2) {
      for (x[2] = x_o[2]; x[2] <= x_f[2]; x[2] += 2) {
	for (x[1] = x_o[1]; x[1] <= x_f[1]; x[1] += 2) {
	  for (x[0] = x_o[0]; x[0] <= x_f[0]; x[0] += 2) {

	    start_seed = start_seed_4d = 0;

	    if (GJP.StartSeedKind () == START_SEED ||
		GJP.StartSeedKind () == START_SEED_INPUT ||
		GJP.StartSeedKind () == START_SEED_FIXED) {
//                    start_seed = base_seed
//                        + OFFSET * (x[0]/2 + vx[0]*
	      start_seed = (x[0] / 2 + vx[0] *
			    (x[1] / 2 + vx[1] *
			     (x[2] / 2 + vx[2] *
			      (x[3] / 2 + vx[3] * (x[4] / 2 + 1)))));
//                    start_seed_4d = base_seed
//                        + OFFSET * (x[0]/2 + vx[0]*
	      start_seed_4d = (x[0] / 2 + vx[0] *
			       (x[1] / 2 + vx[1] *
				(x[2] / 2 + vx[2] * (x[3] / 2))));
	    }
	    if (GJP.Gparity ()) {
	      start_seed = base_seed + OFFSET * start_seed;
	      start_seed_4d = base_seed + OFFSET * start_seed_4d;
	      start_seed_stacked = start_seed + n_rgen_glb / 2 * OFFSET;	//offset seed by 23* global number of flavor 0 RNGs
	      start_seed_stacked_4d =
		start_seed_4d + n_rgen_4d_glb / 2 * OFFSET;

	      if (x[4] == x_o[4]) {
		ugran_4d[index_4d].Reset (start_seed_4d);
		ugran_4d[stk_index_4d].Reset (start_seed_stacked_4d);
#ifdef RNG_WARMUP
		{
		  int n_warm = ugran_4d[index_4d].Urand (N_WARMUP, 0);
		  VRB.RNGSeed (cname, fname, "index_4d=%d n_warm=%d\n", index_4d,
			      n_warm);
		  while (n_warm > 0) {
		    int temp = ugran_4d[index_4d].Urand (100, 0);
		    n_warm--;
		  }
		  n_warm = ugran_4d[stk_index_4d].Urand (N_WARMUP, 0);
		  VRB.RNGSeed (cname, fname, "index_4d=%d n_warm=%d\n",
			      stk_index_4d, n_warm);
		  while (n_warm > 0) {
		    int temp = ugran_4d[stk_index_4d].Urand (100, 0);
		    n_warm--;
		  }
		}
#endif
		index_4d++;
		stk_index_4d++;
	      }
	      ugran[index].Reset (start_seed);
	      ugran[stk_index].Reset (start_seed_stacked);
#ifdef RNG_WARMUP
	      {
		int n_warm = ugran[index].Urand (N_WARMUP, 0);
		  VRB.RNGSeed (cname, fname, "index=%d n_warm=%d\n", index, n_warm);
		while (n_warm > 0) {
		  int temp = ugran[index].Urand (100, 0);
		  n_warm--;
		}
		n_warm = ugran[stk_index].Urand (N_WARMUP, 0);
		while (n_warm > 0) {
		  int temp = ugran[stk_index].Urand (100, 0);
		  n_warm--;
		}
	      }
#endif
	      index++;
	      stk_index++;
	    } else {		//GJP.Gparity
	      start_seed = base_seed + OFFSET * start_seed;
	      VRB.Debug (cname, fname,
			   "%d %d %d %d %d index=%d start_seed= %d\n", x[0],
			   x[1], x[2], x[3], x[4], index, start_seed);
	      ugran[index].Reset (start_seed);
#ifdef BOOTSTRAP
	      {
		while (rng_count < start_seed) {
		  rng_seed_5d.Urand (1, 0);
		  rng_count++;
		}
//              int new_seed = ugran[index].Urand(RandomGenerator::MBIG,0);
		int new_seed = rng_seed_5d.Urand (RandomGenerator::MBIG, 0);
		rng_count++;
		VRB.RNGSeed (cname, fname,
			     "index=%d start_seed=%d new_seed=%d\n", index,
			     start_seed, new_seed);
		ugran[index].Reset (new_seed);
	      }
#endif
#ifdef RNG_WARMUP
	      {
		int n_warm = ugran[index].Urand (N_WARMUP, 0);
		VRB.RNGSeed (cname, fname, "index=%d n_warm=%d\n", index, n_warm);
		while (n_warm > 0) {
		  int temp = ugran[index].Urand (100, 0);
		  n_warm--;
		}
	      }
#endif
	      index++;
	      if (x[4] == x_o[4]) {
		start_seed_4d = base_seed + OFFSET * start_seed_4d;
		VRB.RNGSeed (cname, fname, "index_4d=%d start_seed= %d\n",
			     index_4d, start_seed_4d);
		ugran_4d[index_4d].Reset (start_seed_4d);
#ifdef BOOTSTRAP
		{
		  while (rng_count_4d < start_seed_4d) {
		    rng_seed_4d.Urand (1, 0);
		    rng_count_4d++;
		  }
//              int new_seed = ugran[index_4d].Urand(RandomGenerator::MBIG,0);
		  int new_seed = rng_seed_4d.Urand (RandomGenerator::MBIG, 0);
		  rng_count_4d++;
		  VRB.RNGSeed (cname, fname,
			       "index_4d=%d start_seed_4d=%d new_seed=%d\n",
			       index_4d, start_seed_4d, new_seed);
		  ugran_4d[index_4d].Reset (new_seed);
		}
#endif
#ifdef RNG_WARMUP
		{
		  int n_warm = ugran_4d[index_4d].Urand (N_WARMUP, 0);
		VRB.RNGSeed (cname, fname, "index_4d=%d n_warm= %d\n",
			     index_4d, n_warm);
		  while (n_warm > 0) {
		    int temp = ugran_4d[index_4d].Urand (100, 0);
		    n_warm--;
		  }
		}
#endif
		index_4d++;
	      }

	    }
	  }
	}
      }
    }
    if (GJP.Gparity ()) {	//moving onto next s, we need to shift both the stacked and original write index forward on 4d volume's worth
      // printf("Just finished s-layer between s=%d and %d\n",x[4],x[4]+2);
      // printf("index %d  stk_index %d  (blocks per s layer %d [2 stacked fields])\n",index,stk_index,blocks_per_s_layer);

      index += blocks_per_s_layer / 2;
      stk_index += blocks_per_s_layer / 2;
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
  if (frm_dim == FIVE_D)
    return ugran[rgen_pos].Urand ();
  else
    return ugran_4d[rgen_pos_4d].Urand ();
}

IFloat LatRanGen::Urand (Float hi, Float lo, FermionFieldDimension frm_dim)
{
  char *fname = "Urand(FermionFieldDimension)";
  if (frm_dim == FIVE_D)
    return ugran[rgen_pos].Urand (hi, lo);
  else
    return ugran_4d[rgen_pos_4d].Urand (hi, lo);
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
  if (frm_dim == FIVE_D)
    return ugran[rgen_pos].Grand ();
  else
    return ugran_4d[rgen_pos_4d].Grand ();
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
  UniformRandomGenerator::SetInterval (high, low);
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
  GaussianRandomGenerator::SetSigma (sigma);
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

  int nstacked = 1;
  if (GJP.Gparity ())
    nstacked = 2;

  rgen_pos =
    x + hx[0] * (y +
		 hx[1] * (z +
			  hx[2] * (t + hx[3] * (field_idx + nstacked * s))));
  rgen_pos_4d = x + hx[0] * (y + hx[1] * (z + hx[2] * (t + hx[3] * field_idx)));

  if (field_idx != 0 && !GJP.Gparity ())
    ERR.General (cname, "AssignGenerator(x,y,z,t,s,field_idx)",
		 "Non-zero field index not defined when G-parity not in use\n");
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
  AssignGenerator (x[0], x[1], x[2], x[3], 0, field_idx);
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
  char *fname = "AssignGenerator(i,field_idx)";
  int x = (i / can[0]) % hx[0];
  int y = (i / can[1]) % hx[1];
  int z = (i / can[2]) % hx[2];
  int t = (i / can[3]) % hx[3];
  int s;
  if (GJP.SnodeSites () < 2)
    s = 0;
  else
    s = (i / can[4]) % hx[4];
  int nstacked = 1;
  if (GJP.Gparity ())
    nstacked = 2;

  rgen_pos =
    x + hx[0] * (y +
		 hx[1] * (z +
			  hx[2] * (t + hx[3] * (field_idx + nstacked * s))));
  rgen_pos_4d = x + hx[0] * (y + hx[1] * (z + hx[2] * (t + hx[3] * field_idx)));
  if (rgen_pos >= n_rgen)
    ERR.General (cname, fname, "rgen_pos(%d)>=n_rgen(%d)", rgen_pos, n_rgen);
  if (rgen_pos_4d >= n_rgen_4d)
    ERR.General (cname, fname, "rgen_pos_4d(%d)>=n_rgen_4d(%d)", rgen_pos_4d,
		 n_rgen_4d);
  if (field_idx != 0 && !GJP.Gparity ())
    ERR.General (cname, fname,
		 "Non-zero field index not defined when G-parity not in use\n");
//  Fprintf(stdout,"i=%d x = %d %d %d %d %d rgen_pos=%d ",i,x,y,z,t,s,rgen_pos);
}


//--------------------------------------------------------------
/*!
  \return A uniform random number; one per node, the same on each node.
*/
//--------------------------------------------------------------
IFloat LatRanGen::Lrand (Float hi, Float lo)
{
  Float cntr = 0.0;

  //Get one non-zero float from global RNG coordinate zero
  if (!CoorX () && !CoorY () && !CoorZ () && !CoorT () && !CoorS ())
    cntr = (Float) ugran[0].Urand (hi, lo);	//5d RNG, local coord 0

  glb_sum_five (&cntr);
  return (IFloat) cntr;
}

IFloat LatRanGen::Lrand ()
{
  IFloat hi, lo;
  UniformRandomGenerator::GetInterval (hi, lo);
  return Lrand (hi, lo);
}

/*!
  \return The number of unsigned ints that comprise the RNG state vector.
*/
int LatRanGen::StateSize () const
{
  return ugran[0].StateSize ();
}

/*!
  \pre A RNG must be assigned using ::AssignGenerator.
  \a state must point to an array of unsigned ints (at least) as long
  as the value returned by ::StateSize.

  \param state The state to be copied from the RNG.
  \param frm_dim If FIVE_D, the default, refers to the normal RNG. If FOUR_D
  then the special RNG four gauge fields with domain-wall fermions is used.  
*/
void LatRanGen::GetState (unsigned int *state, FermionFieldDimension frm_dim) const
{
  if (frm_dim == FIVE_D)
    ugran[rgen_pos].StoreSeeds (state);
  else
    ugran_4d[rgen_pos_4d].StoreSeeds (state);
}

/*!
  \pre A RNG must be assigned using ::AssignGenerator.
  \pre \a s must be an array with length given by ::StateSize.
  \param s The state to assigned to the RNG on an assigned site.
*/
void LatRanGen::SetState (const unsigned int *s, FermionFieldDimension frm_dim)
{

  // ugran[rgen_pos].RestoreSeeds(s);
  if (frm_dim == FIVE_D)
    ugran[rgen_pos].RestoreSeeds (s);
  else
    ugran_4d[rgen_pos_4d].RestoreSeeds (s);

}

/*!
  \return The total number of RNG states on the local lattice.
*/
int LatRanGen::NStates (FermionFieldDimension frm_dim) const
{
  if (frm_dim == FOUR_D)
    return n_rgen_4d;
  else
    return n_rgen;
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


bool LatRanGen::Read (const char *filename, int concur_io_number)
{
  io_good = false;
  QioArg rd_arg (filename, concur_io_number);
  LatRngRead rd;
  if (parIO ())
    rd.setParallel ();
  else
    rd.setSerial ();
  if (do_log)
    rd.setLogDir (log_dir);
  rd.read (ugran, ugran_4d, rd_arg);
  return (io_good = rd.good ());
}


bool LatRanGen::Write (const char *filename, int concur_io_number)
{
  io_good = false;
  QioArg wt_arg (filename, concur_io_number);
  VRB.Result (cname, "Write()", "concur_io_number=%d\n", concur_io_number);
  LatRngWrite wt;
  if (parIO ())
    wt.setParallel ();
  else
    wt.setSerial ();
  if (do_log)
    wt.setLogDir (log_dir);
  wt.write (ugran, ugran_4d, wt_arg);
  return (io_good = wt.good ());
}

void LatRanGen::Shift ()
{
  GDS.Shift (ugran, n_rgen * sizeof (UGrandomGenerator));
  GDS.Shift (ugran_4d, n_rgen_4d * sizeof (UGrandomGenerator));
}

CPS_END_NAMESPACE
#endif
