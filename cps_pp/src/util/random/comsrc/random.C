#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/comsrc/random.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: random.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.5  2002/12/04 17:16:27  zs
//  Merged the new 2^4 RNG into the code.
//  This new RNG is implemented in the LatRanGen class.
//  The following algorithm and utility classes are affected:
//
//  AlgEig                  Fdwf
//  AlgGheatBath            Fstag
//  AlgHmd                  GlobalJobParameter
//  AlgNoise                Lattice
//  AlgPbp                  Matrix
//  AlgThreept              RandomGenerator
//                          Vector
//
//  Revision 1.4  2001/08/16 10:50:38  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:35  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:10  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: random.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/random/comsrc/random.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//  random.C
//---------------------------------------------------------------
//  This is the routine from Numerical Recipes in C PP.283 ran3
//---------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/random.h>
#include<util/gjp.h>
#include<util/error.h>
#include<comms/glb.h>
CPS_START_NAMESPACE


// Static variables
int LatRanGen::n_rgen;
int LatRanGen::rgen_pos;
int LatRanGen::can[4];
int LatRanGen::hx[4];
int LatRanGen::is_initialized = 0;
UGrandomGenerator *LatRanGen::ugran;



const int MBIG  = 1000000000;
const int MSEED = 161803398;
const int MZ    = 0;




void RandomGenerator::Reset(int idum)
{
    int i, k, ii;
    int mk, mj;
 

    //-----------------------------------------------------
    //  Initialize ma[55] using the seed idum and the large
    //  number MSEED.
    //-----------------------------------------------------
    mj = MSEED - (idum<0 ? -idum : idum);
    mj %= MBIG;
    if(mj < MZ) mj = -mj; // Added by Roy and Pavlos to protect
                          // for the case where idum > MSEED
    ma[54] = mj;


    //-----------------------------------------------------
    //  Initialize the reset of the table in a slightly
    //  random order, with numbers that are not especially
    //  random.
    //-----------------------------------------------------
    mk = 1;
    for( i = 1; i < 55; i++) {
        ii = (21*i)%55 - 1;
        ma[ii] = mk ;
        mk = mj - mk ;
        if ( mk < MZ ) mk += MBIG ;
        mj = ma[ii] ;
    }

    // Randomize them by "warming up the generator"
    for( k = 0; k < 4; k++) {
        for( i = 0; i < 55; i++) {
	    ma[i] -= ma[(i+31)%55] ;
	    if ( ma[i] < MZ ) ma[i] += MBIG ;
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
    	// to[] anf from[] are buffers of size 55+4

void RandomGenerator::StoreSeeds(unsigned *to)
{
    *to++ = (unsigned )inext;
    *to++ = (unsigned )inextp;
    for(int i = 0; i < 55; ++i) {
        *to++ = ma[i];
    }
}


void RandomGenerator::RestoreSeeds(unsigned *from)
{
    inext = *from++;
    inextp = *from++;

    for(int i = 0; i < 55; ++i) {
        ma[i] = *from++;
    }
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
  Initialize();
}

//---------------------------------------------------------
void LatRanGen::Initialize()
{
  if(is_initialized == 1) {
    return;
  }

  VRB.Func("LatRanGen", "Initialize()");

  n_rgen = GJP.VolNodeSites() / 16;
  if(n_rgen == 0) {
    return;
  }

  is_initialized = 1;

  ugran = 0;
  ugran = new UGrandomGenerator[n_rgen];
  if(ugran == 0) {
    ERR.Pointer("LatRanGen","Initialize()", "ugran");
    VRB.Pmalloc("LatRanGen", "Initialize()",  "ugran", ugran,
                                        sizeof(UGrandomGenerator)*n_rgen);
  }

  int x_o[4];
  x_o[0] = GJP.XnodeCoor() * GJP.XnodeSites();
  x_o[1] = GJP.YnodeCoor() * GJP.YnodeSites();
  x_o[2] = GJP.ZnodeCoor() * GJP.ZnodeSites();
  x_o[3] = GJP.TnodeCoor() * GJP.TnodeSites();

  int x_f[4];
  x_f[0] = x_o[0] + GJP.XnodeSites() - 1;
  x_f[1] = x_o[1] + GJP.YnodeSites() - 1;
  x_f[2] = x_o[2] + GJP.ZnodeSites() - 1;
  x_f[3] = x_o[3] + GJP.TnodeSites() - 1;

  hx[0] = (GJP.XnodeSites() / 2);
  hx[1] = (GJP.YnodeSites() / 2);
  hx[2] = (GJP.ZnodeSites() / 2);
  hx[3] = (GJP.TnodeSites() / 2);

  can[0] = 2;
  can[1] = can[0] * GJP.XnodeSites();
  can[2] = can[1] * GJP.YnodeSites();
  can[3] = can[2] * GJP.ZnodeSites();

  int index = 0;
  int x[4], vx[3];
  vx[0] = hx[0] * GJP.Xnodes();
  vx[1] = hx[1] * GJP.Ynodes();
  vx[2] = hx[2] * GJP.Znodes();

  int start_seed;
  int base_seed = GJP.StartSeedValue();
  for(x[3] = x_o[3]; x[3] <= x_f[3]; x[3]+=2) {
    for(x[2] = x_o[2]; x[2] <= x_f[2]; x[2]+=2) {
      for(x[1] = x_o[1]; x[1] <= x_f[1]; x[1]+=2) {
        for(x[0] = x_o[0]; x[0] <= x_f[0]; x[0]+=2) {
          if( (GJP.StartSeedKind() == START_SEED_FIXED_UNIFORM) ||
              (GJP.StartSeedKind() == START_SEED_UNIFORM)       ||
              (GJP.StartSeedKind() == START_SEED_INPUT_UNIFORM)) {
            start_seed = base_seed;
          }
          else if((GJP.StartSeedKind() == START_SEED_FIXED) ||
                  (GJP.StartSeedKind() == START_SEED)       ||
                  (GJP.StartSeedKind() == START_SEED_INPUT) ) {
            start_seed = base_seed + 23 * (int(x[0]/2) + vx[0]*(int(x[1]/2) +
                 vx[1]*(int(x[2]/2) + vx[2]*(int(x[3]/2)))) );
          }
          else if(GJP.StartSeedKind() == START_SEED_INPUT_NODE) {
            start_seed = base_seed;
            // This option will not create platform independent results
          }
          ugran[index++].Reset(start_seed);
        }
      }
    }
  }
}

//---------------------------------------------------------
IFloat LatRanGen
::Urand(void)
{
  return ugran[rgen_pos].Urand();
}

//---------------------------------------------------------
IFloat LatRanGen
::Grand(void)
{
  return ugran[rgen_pos].Grand();
}


//---------------------------------------------------------
void LatRanGen
::SetInterval(IFloat high, IFloat low)
{
  for(int i=0; i<n_rgen; i++) {
    ugran[i].SetInterval(high, low);
  }
}

//---------------------------------------------------------
void LatRanGen
::SetSigma(IFloat sigma)
{
  for(int i=0; i<n_rgen; i++) {
    ugran[i].SetSigma(sigma);
  }
}

//---------------------------------------------------------
void LatRanGen
::AssignGenerator(int x, int y, int z, int t)
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

void LatRanGen
::AssignGenerator(int * x)
{
  x[0] = x[0] % GJP.XnodeSites();
  x[1] = x[1] % GJP.YnodeSites();
  x[2] = x[2] % GJP.ZnodeSites();
  x[3] = x[3] % GJP.TnodeSites();

  rgen_pos = int(x[0]/2) + hx[0] * (int(x[1]/2) + hx[1] * (
                int(x[2]/2) + hx[2] * int(x[3]/2)));
}

//---------------------------------------------------------
void LatRanGen
::AssignGenerator(int i)
{
  int x = (i/can[0]) % hx[0];
  int y = (i/can[1]) % hx[1];
  int z = (i/can[2]) % hx[2];
  int t = (i/can[3]) % hx[3];

  rgen_pos = x + hx[0] * (y + hx[1] * (z + hx[2] * t));
}


//--------------------------------------------------------------
// Lrand will return the same random number on every node by
// performing a global sum over all 2^4 hypercubes, and taking the
// average value
//--------------------------------------------------------------
IFloat LatRanGen
::Lrand(void)
{
  Float cntr = 0.0;
  for(int i = 0; i < n_rgen; i++) {
    cntr+= (Float) ugran[i].Urand();
  }
  Float divisor = (Float) (n_rgen*GJP.Xnodes()*GJP.Ynodes()*GJP.Znodes()
                                        *GJP.Tnodes());
  cntr/=divisor;
  glb_sum(&cntr);
  IFloat result = (IFloat) cntr;
  return result;
}

CPS_END_NAMESPACE
