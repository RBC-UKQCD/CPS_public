#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief   Methods for the Random Number Generator classes.

  $Id: random.C,v 1.3 2004-04-30 12:18:00 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-04-30 12:18:00 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/comsrc/random.C,v 1.3 2004-04-30 12:18:00 zs Exp $
//  $Id: random.C,v 1.3 2004-04-30 12:18:00 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: random.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/comsrc/random.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

//---------------------------------------------------------------
//  This is the routine ran3 from Numerical Recipes in C 
//---------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/random.h>
#include <util/gjp.h>
#include <util/error.h>
#include <comms/glb.h>
#ifdef PARALLEL
#include <comms/sysfunc.h>
#else
#include <time.h>
#endif
CPS_START_NAMESPACE



/*!
  This method must be called before the RNG is used
  \param idum The seed.	
 */
void RandomGenerator::Reset(int idum)
{
    int i, k, ii;
    int mk, mj;

    const int MSEED = 161803398;
    
    //-----------------------------------------------------
    //  Initialize ma[state_size] using the seed idum and the large
    //  number MSEED.
    //-----------------------------------------------------
    mj = MSEED - (idum<0 ? -idum : idum);
    mj %= MBIG;
    if(mj < 0) mj = -mj; // Added by Roy and Pavlos to protect
                          // for the case where idum > MSEED
    ma[54] = mj;


    //-----------------------------------------------------
    //  Initialize the reset of the table in a slightly
    //  random order, with numbers that are not especially
    //  random.
    //-----------------------------------------------------
    mk = 1;
    for( i = 1; i < state_size; i++) {
        ii = (21*i)%state_size - 1;
        ma[ii] = mk ;
        mk = mj - mk ;
        if ( mk < 0 ) mk += MBIG ;
        mj = ma[ii] ;
    }

    // Randomize them by "warming up the generator"
    for( k = 0; k < 4; k++) {
        for( i = 0; i < state_size; i++) {
	    ma[i] -= ma[(i+31)%state_size] ;
	    if ( ma[i] < 0 ) ma[i] += MBIG ;
        }
    }

    //-----------------------------------------------------
    //  Prepare indices for our first generated number.
    //  the constant 30 is special.
    //-----------------------------------------------------
    inext = -1 ;
    inextp = 30 ;
}




//---------------------------------------------------------------

/*!
  \param to Pointer to a buffer of size at least 57*sizeof(int).
*/

void RandomGenerator::StoreSeeds(unsigned *to) const
{
    *to++ = (unsigned )inext;
    *to++ = (unsigned )inextp;
    for(int i = 0; i < 55; ++i) *to++ = ma[i];
}

/*!
  \param from Pointer to a buffer of size at least 57*sizeof(int).
*/

void RandomGenerator::RestoreSeeds(const unsigned *from)
{
    inext = *from++;
    inextp = *from++;

    for(int i = 0; i < 55; ++i) ma[i] = *from++;
}

/*!
  \return The number of unsigned ints that comprise the RNG state vector.
*/
int RandomGenerator::StateSize() const {

    return state_size+2;

}

//--------------------------------------------------------
IFloat UGrandomGenerator::Rand()
{
  ERR.NotImplemented("UGrandomGenerator", "Rand()");
  return 0.0;
}


//---------------------------------------------------------
LatRanGen::LatRanGen()
{
    cname = "LatRanGen";
    is_initialized = 0;
}

LatRanGen::~LatRanGen() {
	delete[] ugran;
}

/*! Seeds the RNGs according to the method defined in GJP */  
//---------------------------------------------------------
void LatRanGen::Initialize()
{
  if(is_initialized == 1) return;

  const char *fname = "Initialize";
  VRB.Func(cname, fname);

  if(0!=GJP.VolNodeSites()%16)
      ERR.General(cname, fname,
		  "Must have a multiple of 2^4 lattice sites per node.");
  
  n_rgen = GJP.VolNodeSites() / 16;

  

  is_initialized = 1;

  ugran = new UGrandomGenerator[n_rgen];
  if(!ugran) ERR.Pointer(cname, fname, "ugran"); 

  // Lower bounds of global lattice coordinates on this node
  int x_o[4];
  x_o[0] = GJP.XnodeCoor() * GJP.XnodeSites();
  x_o[1] = GJP.YnodeCoor() * GJP.YnodeSites();
  x_o[2] = GJP.ZnodeCoor() * GJP.ZnodeSites();
  x_o[3] = GJP.TnodeCoor() * GJP.TnodeSites();

  // Upper bounds of global lattice coordinates on this node
  int x_f[4];
  x_f[0] = x_o[0] + GJP.XnodeSites() - 1;
  x_f[1] = x_o[1] + GJP.YnodeSites() - 1;
  x_f[2] = x_o[2] + GJP.ZnodeSites() - 1;
  x_f[3] = x_o[3] + GJP.TnodeSites() - 1;

  hx[0] = GJP.XnodeSites() / 2;
  hx[1] = GJP.YnodeSites() / 2;
  hx[2] = GJP.ZnodeSites() / 2;
  hx[3] = GJP.TnodeSites() / 2;

  can[0] = 2;
  can[1] = can[0] * GJP.XnodeSites();
  can[2] = can[1] * GJP.YnodeSites();
  can[3] = can[2] * GJP.ZnodeSites();

  int vx[3];
  vx[0] = hx[0] * GJP.Xnodes();
  vx[1] = hx[1] * GJP.Ynodes();
  vx[2] = hx[2] * GJP.Znodes();


  // Sort out what the seed should be depending on the GJP.StartSeedKind()

  int start_seed, base_seed;

  if(GJP.StartSeedKind()==START_SEED_INPUT_UNIFORM||
     GJP.StartSeedKind()==START_SEED_INPUT_NODE||
     GJP.StartSeedKind()==START_SEED_INPUT)
      base_seed = GJP.StartSeedValue();

  else if(GJP.StartSeedKind()==START_SEED_UNIFORM||
	  GJP.StartSeedKind()==START_SEED)
#ifdef PARALLEL
      base_seed = SeedS();
#else
      base_seed = int (time(NULL));
#endif

  else if(GJP.StartSeedKind()==START_SEED_FIXED_UNIFORM||
	  GJP.StartSeedKind()==START_SEED_FIXED)
      base_seed = default_seed;

  else
      ERR.General(cname, fname,"Unknown StartSeedType %d\n",
		  int(GJP.StartSeedKind()));
    
#ifdef PARALLEL
  if(GJP.StartSeedKind()==START_SEED_INPUT_NODE){
      int node  = GJP.XnodeCoor()
	  + GJP.YnodeCoor() * GJP.Xnodes()
	  + GJP.ZnodeCoor() * GJP.Xnodes() * GJP.Ynodes()
	  + GJP.TnodeCoor() * GJP.Xnodes() * GJP.Ynodes() * GJP.Znodes();
      base_seed = base_seed + 23 * node;
  }
#endif

  // Seed the hypercube RNGs

  int x[4];
  int index = 0;
  
  for(x[3] = x_o[3]; x[3] <= x_f[3]; x[3]+=2) {
      for(x[2] = x_o[2]; x[2] <= x_f[2]; x[2]+=2) {
	  for(x[1] = x_o[1]; x[1] <= x_f[1]; x[1]+=2) {
	      for(x[0] = x_o[0]; x[0] <= x_f[0]; x[0]+=2) {

		  start_seed = base_seed;
		  
		  if(GJP.StartSeedKind()==START_SEED||
		     GJP.StartSeedKind()==START_SEED_INPUT||
		     GJP.StartSeedKind()==START_SEED_FIXED)
		      start_seed = base_seed
			  + 23 * (x[0]/2 + vx[0]*
				  (x[1]/2 + vx[1]*
				   (x[2]/2 + vx[2]*(x[3]/2))) );
		      
		  ugran[index++].Reset(start_seed);
		  
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
IFloat LatRanGen::Urand()
{
  return ugran[rgen_pos].Urand();
}

//---------------------------------------------------------
/*!
  \pre A RNG must be assigned using ::AssignGenerator.
  \return A gaussian random number for this hypercube.
*/
//----------------------------------------------------------------------
IFloat LatRanGen::Grand()
{
  return ugran[rgen_pos].Grand();
}


//---------------------------------------------------------
/*!
  The parameters are set for the RNGs on all hypercubes.
  \param high the upper bound of the distribution range
  \param lower the lower bound of the distribution range
*/
//----------------------------------------------------------------------
void LatRanGen::SetInterval(IFloat high, IFloat low)
{
  for(int i=0; i<n_rgen; i++) ugran[i].SetInterval(high, low);
  
}

//---------------------------------------------------------
/*!
  The parameters are set for the RNGs on all hypercubes. The mean is zero.
  \param sigma the variance of the gaussian distribution.
*/
//----------------------------------------------------------------------
void LatRanGen::SetSigma(IFloat sigma)
{
  for(int i=0; i<n_rgen; i++) ugran[i].SetSigma(sigma);
}

//---------------------------------------------------------
/*!
  For a given lattice site, this identifies and assigns the corresponding
  hypercube RNG.
  \param x The x coordinate of the lattice site.
  \param y The y coordinate of the lattice site.
  \param z The z coordinate of the lattice site.
  \param t The t coordinate of the lattice site.
  \post  Subsequent calls to \e e.g. Urand will return results from this
  particular hypercubic RNG.  
 */
//----------------------------------------------------------------------
void LatRanGen::AssignGenerator(int x, int y, int z, int t)
{
  x = x % GJP.XnodeSites();
  y = y % GJP.YnodeSites();
  z = z % GJP.ZnodeSites();
  t = t % GJP.TnodeSites();

  x/=2;
  y/=2;
  z/=2;
  t/=2;

  rgen_pos = x + hx[0] * (y + hx[1] * (z + hx[2] * t));
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
void LatRanGen::AssignGenerator(const int * x)
{
    AssignGenerator(x[0], x[1], x[2], x[3]);
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
void LatRanGen ::AssignGenerator(int i)
{
  int x = (i/can[0]) % hx[0];
  int y = (i/can[1]) % hx[1];
  int z = (i/can[2]) % hx[2];
  int t = (i/can[3]) % hx[3];

  rgen_pos = x + hx[0] * (y + hx[1] * (z + hx[2] * t));
}


//--------------------------------------------------------------
/*!
  \return A uniform random number; one per node, the same on each node.
*/
// Lrand will return the same random number on every node by
// performing a global sum over all 2^4 hypercubes, and taking the
// average value
//--------------------------------------------------------------
IFloat LatRanGen::Lrand()
{
  Float cntr = 0.0;
  for(int i = 0; i < n_rgen; i++) cntr += (Float) ugran[i].Urand();
  
  Float divisor = (Float) (n_rgen*GJP.Xnodes()*GJP.Ynodes()*GJP.Znodes()
			   *GJP.Tnodes());
  cntr/=divisor;
  glb_sum(&cntr);
  return  (IFloat) cntr;

}


/*!
  \return The number of unsigned ints that comprise the RNG state vector.
*/
int LatRanGen::StateSize() const{

    return ugran[0].StateSize();    
    
}

/*!
  \pre A RNG must be assigned using ::AssignGenerator.
  \a state must point to an array of unsigned ints (at least) as long
  as the value returned by ::StateSize.

  \param state The state to be copied from the RNG.
*/
void LatRanGen::GetState(unsigned *state) const{

    ugran[rgen_pos].StoreSeeds(state);        
    
}

/*!
  \pre A RNG must be assigned using ::AssignGenerator.
  \pre \a s must be an array with length given by ::StateSize.
  \param s The state to assigned to the RNG on an assigned site.
*/
void LatRanGen::SetState(const unsigned* s){

    ugran[rgen_pos].RestoreSeeds(s);

}

/*!
  \return The total number of RNG states on the local lattice.
*/
int LatRanGen::NStates() const{

    return n_rgen;

}

/*!
  \pre \a s must point to a 2-d array with lengths given by ::NStates and 
  ::StateSize.
  \param s The state to assigned to the RNGs on the entire local lattice.
*/
void LatRanGen::SetStates(unsigned **s){

    for(int h=0; h<n_rgen; h++) ugran[h].RestoreSeeds(s[h]);

}
	
/*!

  \pre \a s must point to a 2-d array with lengths given by ::NStates and
  ::StateSize.
  \param s The state to be copied from the RNGs on the entire local lattice. 
*/
void LatRanGen::GetStates(unsigned **s) const{

    for(int h=0; h<n_rgen; h++) ugran[h].StoreSeeds(s[h]);

}

CPS_END_NAMESPACE
