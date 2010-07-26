#include <config.h>
#ifdef PROFILE
#include <stdio.h>
#include <math.h>
#endif
#include <util/gjp.h>
#include <util/time_cps.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/vector.h>
#include <util/pt.h>
#if TARGET==QCDOC
#include <qcdocos.h>
#include <qalloc.h>
#include <ppc_lib.h>
#endif

#define CPP
#undef CPP

CPS_START_NAMESPACE

// ------------------------------------------------------------
// Xiao-Yong Jin    05/07/08
//
// This file is a stripped version of file
// 
//     src/util/dirac_op/d_op_p4/qcdoc/p4_dirac.C
//
// with all knights killed.
//
// The following comments were in the original file.  They are here
// for reference pruposes.
// ------------------------------------------------------------

// --------------------------------------------------------------------------
// Michael Cheng 02/14/05
//
// This file implements the methods that allocate and deallocate memory
// for the P4 dirac operator, as well as implementing the derivative
// term of the fermion matrix.
//
// For this operator, Dslash is similar to Dslash for the Asqtad action.
// However, instead of a straight, three-link Naik term, the P4 action
// uses a bent, knight's move term in addition to the smeared one-link
// part of the operator.
//
// Explicity, in the free field case, the bent link terms look like:
//
// Dslash(Psi(x)) = Sum_mu[Sum_(nu != mu)[
//                  Psi(x+mu+2*nu) + Psi(x+mu-2*nu) -
//                  Psi(x-mu+2*nu) - Psi(x-mu-2*nu)]]
//
// Of course, in the interacting case, the field must be parallel transported
// by multiplication with the appropriate gauge links.
//
// Define the following Parallel transports:
//
// P_mu(x) = U_mu(x);  P_nu_mu(x) = U_nu(x) P_mu(x+nu);  etc...
//
// In terms of these parallel transporers, the bent link term in Dslash
// becomes:
//
// Dslash(Psi(x)) = c12*Sum_mu[Sum_(nu != mu)[
//                  (P_mu_nu_nu(x) + P_nu_nu_mu(x))Psi(x+mu+2*nu) +
//                  (P_mu_-nu_-nu(x) + P_-nu_-nu_mu(x))Psi(x+mu-2*nu) -
//                  (P_-mu_nu_nu(x) + P_nu_nu_-mu(x))Psi(x-mu+2*nu) -
//                  (P_-mu_-nu_-nu(x) + P_-nu_-nu_-mu(x))Psi(x-mu-2*nu)
//
// In this implementation, we rely heavily on the ParTranStaggered_cb class,
// which automatically implements these parallel transports.
//
// We also use a "gathering" technique, where some like terms are gathered
// and summed before all parallel transports are done.
//--------------------------------------------------------------------------
// Christian Schmidt 11/18/05
//
// The derivative of the fermion Matrix with respect to the chemical 
// potential was added.
//--------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
//Assembly routines for recombination
//-------------------------------------------------------------------------

//Recombination when all sites are recombined with the same sign
extern "C" void stag_recom(int sites, Float *src, Float **res);

//Recombination when sites are recombined with alternating sign
extern "C" void stag_recom_n(int sites, Float *src, Float **res);

//Final recombination routine
extern "C" void stag_dsum(int sites, Float * src, Float *f_out, double * one, double * negone);
//-------------------------------------------------------------------------

enum{VECT_LEN=6, VECT_LEN2=8, MATRIX_SIZE=18, SITE_LEN=72, NUM_DIR=8, N=4};

static int size[4];
static int coord[4];
static int vol;

//Holds the initial field
static Vector * fermion[8];

//Holds smeared one link part of Dslash
static IFloat * onelink;

static ParTransStaggered_cb * pt;

//Initializes some important constants, called by P4 lattice class
extern "C"
void stag_dirac_init(const void * gauge_u)
{
  //Using the convention consistent for the staggered parallel transporter
  //that {0,1,2,3} = {X,Y,Z,T}
  size[0] = GJP.XnodeSites();
  size[1] = GJP.YnodeSites();
  size[2] = GJP.ZnodeSites();
  size[3] = GJP.TnodeSites();
  vol = size[0]*size[1]*size[2]*size[3];
}

//Allocates memory for various memory blocks needed for calculation
//Also copies the smeared gauge links
extern "C"
void stag_dirac_init_g (Lattice& lat)
{
  char *cname = "DiracOpStag";
  char *fname = "stag_dirac_init_g()";

  if (vol < 4096) {
      onelink = (IFloat *) fmalloc (NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
  } else {
      onelink = (IFloat *) smalloc (NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
  }

  if(onelink == 0)
    ERR.Pointer(cname,fname, "one_link");
}

// Instantiate the base parallel transport class
extern "C"
void stag_dirac_init_with_lat (Lattice& lat)
{
  pt = new ParTransStaggered_cb(lat);
}

extern "C"
void stag_destroy_dirac_buf()
{
}

//Free memory
extern "C"
void stag_destroy_dirac_buf_g()
{
  if(vol<4096)
    {
      ffree(onelink);
    }
  else
    {
      sfree(onelink);
    }
  delete pt;
}

//-----------------------------------------------------------------------
//This method calculates Dslash on an input field.
//
//f_out - The result of applying Dslash to f_in
//f_in - The input vector field
//cb - 0 if the input field is located on even sites, 1 if odd
//dag - not used in this implementation

#undef PROFILE
extern "C"
void stag_dirac(Vector *f_out, Vector *f_in, int cb, int dag)
{
  if (dag != 0)
    {
      ERR.NotImplemented ("DiracOpStag",
                          "Dslash calls stag_dirac w/ dag != 0");
    }
  //Define the eight different directions.
  //Let mu = dir[i]/2;
  //If dir[i] is even, then this transport in the negative mu direction
  //If dir[i] is odd, then we have transport in the positive mu direction
  char *fname = "stag_dirac()";
  int dir[8] = {0,1,2,3,4,5,6,7};

  IFloat * fp0,*fp2,*fp4,*fp8,*fp10;
  int k,n,mu;

  //Copy the input fields into fermion
  for(int mu = 0; mu < NUM_DIR; mu++)
    fermion[mu] = f_in;

  int nflops=0;
  Float dtime;
#ifdef PROFILE
  dtime = -dclock();
  ParTrans::PTflops = 0;
#endif

  //Calculate the one link term
  pt->run(8,onelink,fermion,dir,(ChkbType) cb,1);

#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT",ParTrans::PTflops,dtime);
  ParTrans::PTflops=0;
  dtime =-dclock();
#endif


  bzero((char*)f_out,(vol/2)*VECT_LEN*sizeof(IFloat));
  #ifdef CPP
   for(n = 0; n < vol/2; n++)
    {
      fp2 = onelink + n*NUM_DIR*VECT_LEN2;
      fp0 = (IFloat *)(f_out) + n*VECT_LEN;

      for(mu = 0; mu < 8; mu+=2)
	{
	  fp8 = fp2+mu*VECT_LEN2;
	  fp10 = fp2+(mu+1)*VECT_LEN2;
	  for(k = 0; k < VECT_LEN; k++)
	    {
	      fp4 = fp0 + k;
	      *(fp4) += *(fp8 + k);
	      *(fp4) -= *(fp10 + k);
	    }
	}
    }
   #else
   double one = 1.0;
   double negone = -1.0;
   stag_dsum(vol/2, (Float *)onelink, (Float *) f_out, &one, &negone);
   #endif


#ifdef PROFILE
  dtime +=dclock();
  nflops = (vol/2)*8*6;
  print_flops(fname,"dsum",nflops,dtime);
  dtime =-dclock();
#endif

  // FIXME: Number of flops!  -- Xiao-Yong Jin
   //Parallel Transport flops:
   //7 parallel transports * (264*vol) flops/transport = 1848*vol
   //Linear combination flops:
   //60*vol + 60*vol + 120*vol = 240*vol
   //Total flop count = 2088*vol

  // FIXME: This value bought from ../qcdoc/dirac.C  -- Xiao-Yong Jin
   DiracOp::CGflops += 285*vol;
}

#undef CPP

CPS_END_NAMESPACE
