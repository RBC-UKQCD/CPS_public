#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief   Methods for the Random Number Generator classes.

  $Id: random.C,v 1.34 2012-05-15 05:50:09 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-05-15 05:50:09 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/comsrc/random.C,v 1.34 2012-05-15 05:50:09 chulwoo Exp $
//  $Id: random.C,v 1.34 2012-05-15 05:50:09 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: random.C,v $
//  $Revision: 1.34 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/comsrc/random.C,v $
//  $State: Exp $
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
#include <comms/glb.h>
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE

int  RandomGenerator::MBIG  = 1000000000;
IFloat  RandomGenerator::FAC = 1.0E-09;			// 1.0/MBIG
const int RandomGenerator::state_size;

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

void RandomGenerator::StoreSeeds(unsigned int *to) const
{
    *to++ = (unsigned int)inext;
    *to++ = (unsigned int)inextp;
    for(int i = 0; i < 55; ++i) *to++ = ma[i];
}

void UGrandomGenerator::StoreSeeds(unsigned int *to) const
{
    RandomGenerator::StoreSeeds(to);
}

/*!
  \param from Pointer to a buffer of size at least 57*sizeof(int).
*/

void RandomGenerator::RestoreSeeds(const unsigned int *from)
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

IFloat UniformRandomGenerator::A = -0.5;
IFloat UniformRandomGenerator::B = 0.5;
IFloat GaussianRandomGenerator::sigma2 = 1.0;



LatRanGen LRG;
const int LatRanGen::default_seed;
//---------------------------------------------------------
LatRanGen::LatRanGen()
{
    cname = "LatRanGen";
    char *fname = "LatRanGen()";
//  printf("%s::%s Entered\n",cname,fname);
    is_initialized = 0;
    UseParIO = 1;
    io_good = 1;
}

LatRanGen::~LatRanGen() {
    char *fname = "~LatRanGen()";
//  printf("%s::%s Entered\n",cname,fname);
  if(is_initialized){
	delete[] ugran;
	delete[] ugran_4d;
  }
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
  
  n_rgen = n_rgen_4d = GJP.VolNodeSites()/16;

  if (GJP.SnodeSites()>=2)
    n_rgen = GJP.VolNodeSites()*GJP.SnodeSites() / 32;
//  VRB.Flow(cname,fname,"n_rgen=%d\n",n_rgen);

  int default_concur=0;
#if TARGET==BGQ
	default_concur=1;
#endif
  

  is_initialized = 1;

  ugran = new UGrandomGenerator[n_rgen];
  if(!ugran) ERR.Pointer(cname, fname, "ugran"); 
  ugran_4d = new UGrandomGenerator[n_rgen_4d];
  if(!ugran) ERR.Pointer(cname, fname, "ugran_4d"); 

  // Lower bounds of global lattice coordinates on this node
  int x_o[5];
  x_o[0] = GJP.XnodeCoor() * GJP.XnodeSites();
  x_o[1] = GJP.YnodeCoor() * GJP.YnodeSites();
  x_o[2] = GJP.ZnodeCoor() * GJP.ZnodeSites();
  x_o[3] = GJP.TnodeCoor() * GJP.TnodeSites();
  x_o[4] = GJP.SnodeCoor() * GJP.SnodeSites();

  // Upper bounds of global lattice coordinates on this node
  int x_f[5];
  x_f[0] = x_o[0] + GJP.XnodeSites() - 1;
  x_f[1] = x_o[1] + GJP.YnodeSites() - 1;
  x_f[2] = x_o[2] + GJP.ZnodeSites() - 1;
  x_f[3] = x_o[3] + GJP.TnodeSites() - 1;
  x_f[4] = x_o[4] + GJP.SnodeSites() - 1;


  hx[0] = GJP.XnodeSites() / 2;
  hx[1] = GJP.YnodeSites() / 2;
  hx[2] = GJP.ZnodeSites() / 2;
  hx[3] = GJP.TnodeSites() / 2;
  hx[4] = GJP.SnodeSites() / 2;

  can[0] = 2;
  can[1] = can[0] * GJP.XnodeSites();
  can[2] = can[1] * GJP.YnodeSites();
  can[3] = can[2] * GJP.ZnodeSites();
  can[4] = can[3] * GJP.TnodeSites();

  int vx[5];
  vx[0] = hx[0] * GJP.Xnodes();
  vx[1] = hx[1] * GJP.Ynodes();
  vx[2] = hx[2] * GJP.Znodes();
  vx[3] = hx[3] * GJP.Tnodes();
  vx[4] = hx[4] * GJP.Snodes();

  int index, index_4d;
  index = index_4d = 0;

  /*PAB: Implement the Britney and Christina tests correctly*/
#ifdef UNIFORM_SEED_TESTING
  if (1) {
#else
  if ( GJP.StartSeedKind() ==  START_SEED_FIXED_UNIFORM ) {
#endif
    int x[5];
    int start_seed = default_seed;

    for(x[4] = 0; x[4] < GJP.SnodeSites(); x[4]+=2) {
    for(x[3] = 0; x[3] < GJP.TnodeSites(); x[3]+=2) {
    for(x[2] = 0; x[2] < GJP.ZnodeSites(); x[2]+=2) {
    for(x[1] = 0; x[1] < GJP.YnodeSites(); x[1]+=2) {
    for(x[0] = 0; x[0] < GJP.XnodeSites(); x[0]+=2) {

      start_seed += 23;
      ugran[index++].Reset(start_seed);
      if ( x[4] == 0 ) ugran_4d[index_4d++].Reset(start_seed);

    }
    }
    }
    }
    }
    return;
  }



  // Sort out what the seed should be depending on the GJP.StartSeedKind()

  int start_seed, base_seed;
  int start_seed_4d, base_seed_4d;

  switch(GJP.StartSeedKind()){
  case START_SEED_FILE:
	if ( !LatRanGen::Read(GJP.StartSeedFilename(),default_concur) ) {
	      ERR.General(cname, fname,
		  "Reading file %s",GJP.StartSeedFilename());
	} 
	return;
	break;
  case START_SEED_INPUT_UNIFORM:
  case START_SEED_INPUT_NODE:
  case START_SEED_INPUT:
      base_seed = GJP.StartSeedValue();
      break;
  case START_SEED_UNIFORM:
  case START_SEED:
      base_seed = SeedS();
      break;
  default:
      base_seed = default_seed;
  }
    
#ifdef PARALLEL
  if(GJP.StartSeedKind()==START_SEED_INPUT_NODE){
      int node  = 
	GJP.XnodeCoor() + GJP.Xnodes()* (
	GJP.YnodeCoor() + GJP.Ynodes()* (
	GJP.ZnodeCoor() + GJP.Znodes()* (
	GJP.TnodeCoor() + GJP.Tnodes()*  GJP.SnodeCoor() )));
      base_seed = base_seed + 23 * node;
  }
#endif

  // Seed the hypercube RNGs

  int x[5];
//  int index, index_4d;
  index = index_4d = 0;
  
for(x[4] = x_o[4]; x[4] <= x_f[4]; x[4]+=2) {
  for(x[3] = x_o[3]; x[3] <= x_f[3]; x[3]+=2) {
      for(x[2] = x_o[2]; x[2] <= x_f[2]; x[2]+=2) {
	  for(x[1] = x_o[1]; x[1] <= x_f[1]; x[1]+=2) {
	      for(x[0] = x_o[0]; x[0] <= x_f[0]; x[0]+=2) {

		  start_seed = start_seed_4d = base_seed;
		  
		  if(GJP.StartSeedKind()==START_SEED||
		     GJP.StartSeedKind()==START_SEED_INPUT||
		     GJP.StartSeedKind()==START_SEED_FIXED){
		      start_seed = base_seed
			  + 23 * (x[0]/2 + vx[0]*
				 (x[1]/2 + vx[1]*
				 (x[2]/2 + vx[2]*
				 (x[3]/2 + vx[3]*(x[4]/2+1) ))));
		      start_seed_4d = base_seed
			  + 23 * (x[0]/2 + vx[0]*
				 (x[1]/2 + vx[1]*
				 (x[2]/2 + vx[2]*(x[3]/2) )));
		  }
//		  Fprintf(stderr,"%d %d %d %d %d",x[0],x[1],x[2],x[3],x[4]);
		  VRB.Debug(cname,fname,"index=%d start_seed= %d\n",index,start_seed);
		  ugran[index++].Reset(start_seed);
		  if(x[4]==x_o[4]){
		  	VRB.Debug(cname,fname,"index_4d=%d start_seed= %d\n",index_4d,start_seed_4d);
 			ugran_4d[index_4d++].Reset(start_seed_4d);
		  }
	      }
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
IFloat LatRanGen::Urand(FermionFieldDimension frm_dim)
{
  char *fname = "Urand(FermionFieldDimension)";
  if (frm_dim == FIVE_D)
    return ugran[rgen_pos].Urand();
  else
    return ugran_4d[rgen_pos_4d].Urand();
}

IFloat LatRanGen::Urand(Float hi, Float lo, FermionFieldDimension frm_dim)
{
  char *fname = "Urand(FermionFieldDimension)";
  if (frm_dim == FIVE_D)
    return ugran[rgen_pos].Urand(hi,lo);
  else
    return ugran_4d[rgen_pos_4d].Urand(hi,lo);
}

//---------------------------------------------------------
/*!
  \pre A RNG must be assigned using ::AssignGenerator.
  \return A gaussian random number for this hypercube.
*/
//----------------------------------------------------------------------
IFloat LatRanGen::Grand(FermionFieldDimension frm_dim)
{
  char *fname = "Grand(FermionFieldDimension)";
//  printf("LatRanGen::Grand():%d\n",rgen_pos);
  if (frm_dim == FIVE_D)
    return ugran[rgen_pos].Grand();
  else
    return ugran_4d[rgen_pos_4d].Grand();
}


//---------------------------------------------------------
/*!
  The parameters are set for the RNGs on all hypercubes.
  \param high the upper bound of the distribution range
  \param low the lower bound of the distribution range
*/
//----------------------------------------------------------------------
void LatRanGen::SetInterval(IFloat high, IFloat low)
{
//  ugran[0].SetInterval(high,low);
    UniformRandomGenerator::SetInterval(high,low);
}

//---------------------------------------------------------
/*!
  The parameters are set for the RNGs on all hypercubes. The mean is zero.
  \param sigma the variance of the gaussian distribution.
*/
//----------------------------------------------------------------------
void LatRanGen::SetSigma(IFloat sigma)
{
//  ugran[0].SetSigma(sigma);
    GaussianRandomGenerator::SetSigma(sigma);
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
void LatRanGen::AssignGenerator(int x, int y, int z, int t,int s)
{
  x = x % GJP.XnodeSites();
  y = y % GJP.YnodeSites();
  z = z % GJP.ZnodeSites();
  t = t % GJP.TnodeSites();

  if (GJP.SnodeSites()<2) s = 0;
  else s = s % GJP.SnodeSites();

  x/=2;
  y/=2;
  z/=2;
  t/=2;
  s/=2;

  rgen_pos = x + hx[0] * (y + hx[1] * (z + hx[2] * (t +hx[3] * s)));
  rgen_pos_4d = x + hx[0] * (y + hx[1] * (z + hx[2] * t ));
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
    AssignGenerator(x[0], x[1], x[2], x[3], 0);
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
  char *fname = "AssignGenerator(i)";
  int x = (i/can[0]) % hx[0];
  int y = (i/can[1]) % hx[1];
  int z = (i/can[2]) % hx[2];
  int t = (i/can[3]) % hx[3];
  int s;
  if (GJP.SnodeSites()<2) s = 0;
  else  s = (i/can[4]) % hx[4];
  rgen_pos = x + hx[0] * (y + hx[1] * (z + hx[2] * (t + hx[3] * s)));
  rgen_pos_4d = x + hx[0] * (y + hx[1] * (z + hx[2] * t ));
  if (rgen_pos >=  n_rgen)
      ERR.General(cname, fname, "rgen_pos(%d)>=n_rgen(%d)",rgen_pos,n_rgen);
  if (rgen_pos_4d >=  n_rgen_4d)
      ERR.General(cname, fname, "rgen_pos_4d(%d)>=n_rgen_4d(%d)",rgen_pos_4d,n_rgen_4d);
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
IFloat LatRanGen::Lrand()
{
  Float cntr = 0.0;
  for(int i = 0; i < n_rgen; i++) cntr += (Float) ugran[i].Urand();
  
  Float divisor = (Float) (n_rgen*GJP.Xnodes()*GJP.Ynodes()*GJP.Znodes()
			   *GJP.Tnodes());
  if(GJP.Snodes()>1) divisor *= (Float) GJP.Snodes();
  cntr/=divisor;
  glb_sum_five(&cntr);
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
  \param frm_dim If FIVE_D, the default, refers to the normal RNG. If FOUR_D
  then the special RNG four gauge fields with domain-wall fermions is used.  
*/
void LatRanGen::GetState(unsigned int *state,
			 FermionFieldDimension frm_dim) const{
    if (frm_dim == FIVE_D)
    ugran[rgen_pos].StoreSeeds(state);        
    else
    ugran_4d[rgen_pos_4d].StoreSeeds(state);        
}

/*!
  \pre A RNG must be assigned using ::AssignGenerator.
  \pre \a s must be an array with length given by ::StateSize.
  \param s The state to assigned to the RNG on an assigned site.
*/
void LatRanGen::SetState(const unsigned int * s,
			 FermionFieldDimension frm_dim) {

   // ugran[rgen_pos].RestoreSeeds(s);
    if (frm_dim == FIVE_D)
    ugran[rgen_pos].RestoreSeeds(s);        
    else
    ugran_4d[rgen_pos_4d].RestoreSeeds(s);        

}

/*!
  \return The total number of RNG states on the local lattice.
*/
int LatRanGen::NStates(FermionFieldDimension frm_dim) const{
  if (frm_dim == FOUR_D) return n_rgen_4d;
  else return n_rgen;
}

/*!
  \pre \a s must point to a 2-d array with lengths given by ::NStates and 
  ::StateSize.
  \param s The state to assigned to the RNGs on the entire local lattice.
  \param frm_dim If FIVE_D, the default, refers to the normal RNG. If FOUR_D
  then the special RNG four gauge fields with domain-wall fermions is used.  
*/
void LatRanGen::SetStates(unsigned int **s, FermionFieldDimension frm_dim) {
    if (frm_dim == FIVE_D)
    for(int h=0; h<n_rgen; h++) ugran[h].RestoreSeeds(s[h]);
    else
    for(int h=0; h<n_rgen_4d; h++) ugran_4d[h].RestoreSeeds(s[h]);

}
	
/*!
  \pre \a s must point to a 2-d array with lengths given by ::NStates and
  ::StateSize.
  \param s The state to be copied from the RNGs on the entire local lattice. 
  \param frm_dim If FIVE_D, the default, refers to the normal RNG. If FOUR_D
  then the special RNG four gauge fields with domain-wall fermions is used.  
*/
void LatRanGen::GetStates(unsigned int **s,
			  FermionFieldDimension frm_dim) const {

    if (frm_dim == FIVE_D)
    for(int h=0; h<n_rgen; h++) ugran[h].StoreSeeds(s[h]);
    else
    for(int h=0; h<n_rgen_4d; h++) ugran_4d[h].StoreSeeds(s[h]);

}


bool LatRanGen::Read(const char * filename, int concur_io_number) {
  io_good = false;
  QioArg rd_arg(filename, concur_io_number);
  LatRngRead  rd;
  if(parIO()) rd.setParallel();
  else        rd.setSerial();
  if(do_log) rd.setLogDir(log_dir);
  rd.read(ugran,ugran_4d,rd_arg);
  return (io_good = rd.good());
}


bool LatRanGen::Write(const char * filename, int concur_io_number) {
  io_good = false;
  QioArg wt_arg(filename, concur_io_number);
  VRB.Result(cname,"Write()","concur_io_number=%d\n",concur_io_number);
  LatRngWrite  wt;
  if(parIO()) wt.setParallel();
  else        wt.setSerial();
  if(do_log) wt.setLogDir(log_dir);
  wt.write(ugran,ugran_4d,wt_arg);
  return (io_good=wt.good());
}

void LatRanGen::Shift()
{
   GDS.Shift(ugran, n_rgen*sizeof(UGrandomGenerator));
   GDS.Shift(ugran_4d, n_rgen_4d*sizeof(UGrandomGenerator));
}
CPS_END_NAMESPACE
