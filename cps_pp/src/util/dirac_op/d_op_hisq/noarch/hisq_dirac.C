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
#include <util/qcdio.h>
#if TARGET==QCDOC
#include <qcdocos.h>
#include <qalloc.h>
#include <ppc_lib.h>
#endif

#define CPP

CPS_START_NAMESPACE

// T.Blum copied from modified p4->asqtad for hisq
// 7/23/2015
//--------------------------------------------------------------------------
// Michael Cheng  02/14/05
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

#ifndef CPP
//Recombination when all sites are recombined with the same sign
extern "C" void hisq_recom(int sites, Float *src, Float **res);

//Recombination when sites are recombined with alternating sign
extern "C" void hisq_recom_n(int sites, Float *src, Float **res);

//Final recombination routine
extern "C" void hisq_dsum(int sites, Float * src, Float *f_out, double * one, double * negone);
#endif
//-------------------------------------------------------------------------

enum{VECT_LEN=6, VECT_LEN2=6, MATRIX_SIZE=18, SITE_LEN=72, NUM_DIR=8, N=4};

static int size[4];
static int coord[4];
static int vol;

//Holds the initial field
static Vector * fermion[8];

//Temporary fields used for calculation
static Vector * tmp_frm[8];
static Vector * tmp_frm2[8];

//Holds smeared one link part of Dslash
static Vector * smeared_onelink[8];

//Holds the smeared gauge fields
static IFloat * smeared_gauge;

//One link and bent link coefficients
static IFloat c_onelink;

static ParTransStaggered_cb * pt;

static int LexVector(int * x);
static int LexGauge(int * x, int dir);
static int LexGauge_block(int *x, int dir);
static int LexGauge_block_cb(int *x, int dir);
static void cpy(IFloat *dest, IFloat *src,int len);
static void cpy(IFloat *dest, IFloat *src);
static void dagcpy(IFloat *dest, IFloat *src);

//Sets a pointer to the lattice class
static Fhisq *lat_pt;
void set_pt(Fhisq *lat)
{
  lat_pt = lat;
}


//Initializes some important constants, called by Hisq lattice class
extern "C"
void hisq_dirac_init(const void * gauge_u)
{
  //Using the convention consistent for the staggered parallel transporter
  //that {0,1,2,3} = {X,Y,Z,T}
  size[0] = GJP.XnodeSites();
  size[1] = GJP.YnodeSites();
  size[2] = GJP.ZnodeSites();
  size[3] = GJP.TnodeSites();
  vol = size[0]*size[1]*size[2]*size[3];
  c_onelink = 1.0;
}

//Allocates memory for various memory blocks needed for calculation
//Also copies the smeared gauge links
extern "C"
void hisq_dirac_init_g()
{
  char *cname = "DiracOpHisq";
  char *fname = "hisq_dirac_init_g()";

  if (vol < 4096) {
    for(int i = 0; i < 8; i++) {
      smeared_onelink[i] = (Vector *) fmalloc (vol/2 * VECT_LEN * sizeof(IFloat));
      tmp_frm[i] = (Vector *) fmalloc (vol/2*VECT_LEN*sizeof(IFloat));
      tmp_frm2[i] = (Vector *) fmalloc (vol/2*VECT_LEN*sizeof(IFloat));
    }
  } else {
    for(int i = 0; i < 8; i++) {
      smeared_onelink[i] = (Vector *) smalloc (vol/2 * VECT_LEN * sizeof(IFloat));
      tmp_frm[i] = (Vector *) smalloc (vol/2*VECT_LEN*sizeof(IFloat));
      tmp_frm2[i] = (Vector *) smalloc (vol/2*VECT_LEN*sizeof(IFloat));
    }
  }
  smeared_gauge = (IFloat *) fmalloc(vol*SITE_LEN*sizeof(IFloat));

  if(smeared_gauge == 0)
    ERR.Pointer(cname,fname, "smeared_gauge");
  if(smeared_onelink == 0)
    ERR.Pointer(cname,fname, "smeared_onelink");

  //Copy the smeared links into smeared_gauge
  //The smeared links are calculated by the Smear() method
  //in the lattice class.
  lat_pt->Smear();
  Matrix *Fat = lat_pt->Fields(0);
  IFloat *fp0;
  IFloat *fp1;
  for(coord[3] = 0; coord[3] < size[3]; coord[3]++)
    for(coord[2] = 0; coord[2] < size[2]; coord[2]++)
      for(coord[1] = 0; coord[1] < size[1]; coord[1]++)
	for(coord[0] = 0; coord[0] < size[0]; coord[0]++)
	  for(int i = 0; i < N; i++)
	    {
	      fp0 = smeared_gauge + LexGauge_block_cb(coord, i)*MATRIX_SIZE;
	      fp1 = (IFloat *)(Fat + i*vol + LexVector(coord));
	      dagcpy(fp0,fp1);
	    }
  //Instantiate the base parallel transport class
  pt = new ParTransStaggered_cb(*(lat_pt));
#if 0
  for(int nu=0;nu<4;nu++){
    Matrix* link = (Matrix*)smeared_gauge + nu;
    for(int i=0; i<3; i++)
	for(int j=0; j<3; j++)
	  printf("link %d nu %d  %d %d  %g %g\n",0,nu,i,j,link[nu](i,j).real(),link[nu](i,j).imag());
  }
#endif
}

extern "C"
void hisq_destroy_dirac_buf()
{}

//Free memory
extern "C"
void hisq_destroy_dirac_buf_g()
{
  if(vol<4096)
    {
      ffree(smeared_gauge);
      for(int i = 0; i < 8; i++)
	{
          ffree(smeared_onelink[i]);
	  ffree(tmp_frm[i]);
	  ffree(tmp_frm2[i]);
	}
    }
  else
    {
      ffree(smeared_gauge);
      for(int i = 0; i < 8; i++)
	{
          ffree(smeared_onelink[i]);
	  ffree(tmp_frm[i]);
	  ffree(tmp_frm2[i]);
	}
    }
  delete pt;
}

static void dagcpy(IFloat * dest, IFloat *src)
{
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      {
	*(dest + 6*i +2*j) = *(src + 6*j + 2*i);
	*(dest + 6*i +2*j+1) = -*(src + 6*j + 2*i+1);
      }
}
static void cpy(IFloat * dest, IFloat *src)
{
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      {
	*(dest + 6*i +2*j) = *(src + 6*i + 2*j);
	*(dest + 6*i +2*j+1) = -*(src + 6*i + 2*j+1);
      }
}

//Return index for CANONICAL ordering of fields
static int LexVector(int * x)
{
  return x[0]+size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3]));
}

//Return index for CANONICAL or STAG ordering of gauge links
static int LexGauge(int * x, int dir)
{
  return N*LexVector(x) + dir;
}

//Return index for gauge links block-ordered with sites ordered
//like txyz
static int LexGauge_block(int *x, int dir)
{
  return dir*vol+x[3]+size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2]));
}

static int LexGauge_block_cb(int *x, int dir)
{
  int result = (x[3]+size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2])))/2;
  return dir*vol+result + ((x[0]+x[1]+x[2]+x[3])%2)*vol/2;
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
void hisq_dirac(Vector *f_out, Vector *f_in, int cb, int dag)
{
  //Define the eight different directions.
  //Let mu = dir[i]/2;
  //If dir[i] is even, then this transport in the negative mu direction
  //If dir[i] is odd, then we have transport in the positive mu direction
  char *fname = "hisq_dirac()";
  int dir[8] = {0,1,2,3,4,5,6,7};

  IFloat * fp0,*fp1,*fp2,*fp3,*fhisq,*fp5,*fp6,*fp7,*fp8,*fp9,*fp10;
  Vector * tmp_frm_p[8];
  Float * res_p[4];
  int k,n,mu,nu = 0;

  //Copy the input fields into fermion
  for(int mu = 0; mu < NUM_DIR; mu++)
    fermion[mu] = f_in;

  int nflops=0;
  Float dtime;
#ifdef PROFILE
  dtime = -dclock();
  ParTrans::PTflops = 0;
#endif

  //Calculate the smeared one link term

  //no pad
  pt->run(8,smeared_onelink,fermion,dir,(ChkbType) cb,smeared_gauge);



  //sum the smeared one link terms and place them in fout

  bzero((char*)f_out,(vol/2)*VECT_LEN*sizeof(IFloat));
#ifdef CPP
   for(n = 0; n < vol/2; n++)
    {

      fp0 = (IFloat *)(f_out) + n*VECT_LEN;

      for(mu = 0; mu < 8; mu+=2)
	{

          fp8 = (IFloat*)smeared_onelink[mu] + n*VECT_LEN2;
          fp10 = (IFloat*)smeared_onelink[mu+1] + n*VECT_LEN2;

	  for(k = 0; k < VECT_LEN; k++)
	    {
	      fhisq = fp0 + k;
	      *(fhisq) += *(fp8 + k);
	      *(fhisq) -= *(fp10 + k);
	    }
	}
    }
#else
   double one = 1.0;
   double negone = -1.0;
   hisq_dsum(vol/2, (Float *)smeared_onelink, (Float *) f_out, &one, &negone);
#endif

#ifdef PROFILE
  dtime +=dclock();
  nflops = (vol/2)*(2*8*2+8)*6;
  print_flops(fname,"dsum",nflops,dtime);
  dtime =-dclock();
#endif
  
   DiracOp::CGflops += 2088*vol;
}


//-----------------------------------------------------------------------
//This method calculates dMdmu on an input field.
//
//f_out - The result of applying dMdmu to f_in
//f_in - The input vector field
//cb - 0 if the input field is located on even sites, 1 if odd
//dag - not used in this implementation
//order - the order of the derivative with respect to mu

#undef PROFILE
#define CPP
extern "C"
void hisq_dMdmu(Vector *f_out, Vector *f_in, int cb, int dag, int order)
{
  char *fname = "hisq_dMdmu()";

  ERR.General("DiracOpHisq",fname,"Not Implemented\n");

}

#undef CPP

CPS_END_NAMESPACE
