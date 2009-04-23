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

CPS_START_NAMESPACE

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
extern "C" void p4_recom(int sites, Float *src, Float **res);

//Recombination when sites are recombined with alternating sign
extern "C" void p4_recom_n(int sites, Float *src, Float **res);

//Final recombination routine
extern "C" void p4_dsum(int sites, Float * src, Float *f_out, double * one, double * negone);
#endif
//-------------------------------------------------------------------------

enum{VECT_LEN=6, VECT_LEN2=8, MATRIX_SIZE=18, SITE_LEN=72, NUM_DIR=8, N=4};

static int size[4];
static int coord[4];
static int vol;

//Holds the initial field
static Vector * fermion[8];

//Temporary fields used for calculation
static Vector * tmp_frm[8];
static Vector * tmp_frm2[8];

//Holds fields of the form P_mu_nu_nu(x)
static IFloat * knight_onetwo;

//Holds fields of the form P_nu_nu_mu(x)
static IFloat * knight_twoone;

//Holds smeared one link part of Dslash
static IFloat * smeared_onelink;

//Holds the smeared gauge fields
static IFloat * smeared_gauge;

//One link and bent link coefficients
static IFloat c_onelink;
static IFloat c_knight;

static ParTransStaggered_cb * pt;

static int LexVector(int * x);
static int LexGauge(int * x, int dir);
static int LexGauge_block(int *x, int dir);
static int LexGauge_block_cb(int *x, int dir);
static void cpy(IFloat *dest, IFloat *src,int len);
static void dagcpy(IFloat *dest, IFloat *src);

//Sets a pointer to the lattice class
static Fp4 *lat_pt;
void set_pt(Fp4 *lat)
{
  lat_pt = lat;
}


//Initializes some important constants, called by P4 lattice class
extern "C"
void p4_dirac_init(const void * gauge_u)
{
  //Using the convention consistent for the staggered parallel transporter
  //that {0,1,2,3} = {X,Y,Z,T}
  size[0] = GJP.XnodeSites();
  size[1] = GJP.YnodeSites();
  size[2] = GJP.ZnodeSites();
  size[3] = GJP.TnodeSites();
  vol = size[0]*size[1]*size[2]*size[3];
  c_onelink = GJP.p4_KS_coeff();
  c_knight = GJP.p4_knight_coeff();
}

//Allocates memory for various memory blocks needed for calculation
//Also copies the smeared gauge links
extern "C"
void p4_dirac_init_g()
{
  char *cname = "DiracOpP4";
  char *fname = "p4_dirac_init_g()";

  if (vol < 4096) {
      knight_onetwo = (IFloat *) fmalloc (NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      knight_twoone = (IFloat *) fmalloc (NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      smeared_onelink = (IFloat *) fmalloc (NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      for(int i = 0; i < 8; i++) {
	  tmp_frm[i] = (Vector *) fmalloc (vol/2*VECT_LEN*sizeof(IFloat));
	  tmp_frm2[i] = (Vector *) fmalloc (vol/2*VECT_LEN*sizeof(IFloat));
      }
  } else {
      knight_onetwo = (IFloat *) smalloc (NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      knight_twoone = (IFloat *) smalloc (NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      smeared_onelink = (IFloat *) smalloc (NUM_DIR * vol/2 * VECT_LEN2 * sizeof(IFloat));
      for(int i = 0; i < 8; i++) {
	  tmp_frm[i] = (Vector *) smalloc (vol/2*VECT_LEN*sizeof(IFloat));
	  tmp_frm2[i] = (Vector *) smalloc (vol/2*VECT_LEN*sizeof(IFloat));
      }
  }
  smeared_gauge = (IFloat *) fmalloc(vol*SITE_LEN*sizeof(IFloat));

  if(smeared_gauge == 0)
    ERR.Pointer(cname,fname, "smeared_gauge");
  if(smeared_onelink == 0)
    ERR.Pointer(cname,fname, "smeared_one_link");
  if(knight_onetwo == 0) 
    ERR.Pointer(cname,fname, "knight_onetwo");
  if(knight_twoone == 0) 
    ERR.Pointer(cname,fname, "knight_twoone");

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
  //ParTransStaggered_cb partrans(*(lat_pt));
  pt = new ParTransStaggered_cb(*(lat_pt));	    
}

extern "C"
void p4_destroy_dirac_buf()
{
}

//Free memory
extern "C"
void p4_destroy_dirac_buf_g()
{
  if(vol<4096)
    {
      ffree(knight_onetwo);
      ffree(knight_twoone);
      ffree(smeared_onelink);
      ffree(smeared_gauge);
      for(int i = 0; i < 8; i++)
	{
	  ffree(tmp_frm[i]);
	  ffree(tmp_frm2[i]);
	}
    }
  else
    {
      ffree(knight_onetwo);
      ffree(knight_twoone);
      ffree(smeared_onelink);
      ffree(smeared_gauge);
      for(int i = 0; i < 8; i++)
	{
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
void p4_dirac(Vector *f_out, Vector *f_in, int cb, int dag)
{
  //Define the eight different directions.
  //Let mu = dir[i]/2;
  //If dir[i] is even, then this transport in the negative mu direction
  //If dir[i] is odd, then we have transport in the positive mu direction
  char *fname = "p4_dirac()";
  int dir[8] = {0,1,2,3,4,5,6,7};

  IFloat * fp0,*fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8,*fp9,*fp10;
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

  //First step, gather one link and two link parallel transports from
  //all eight directions

  //Calculation of P_mu(x); result-> knight_onetwo
  pt->run(8,knight_onetwo,fermion,dir,(ChkbType) cb,1);

#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT",ParTrans::PTflops,dtime);
  ParTrans::PTflops = 0;
  dtime = -dclock();
#endif

  //Calculation of P_nu_nu(x); result-> knight_twoone
  pt->run(8,tmp_frm,fermion,dir,(ChkbType) cb);

#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT",ParTrans::PTflops,dtime);
  ParTrans::PTflops = 0;
  dtime = -dclock();
#endif
  pt->run(8,knight_twoone,tmp_frm,dir,(ChkbType) (1-cb),1);

#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT",ParTrans::PTflops,dtime);
  ParTrans::PTflops = 0;
  dtime = -dclock();
#endif

  //Do appropriate linear combinations for the initial one link term
  //
  //Terms of the form P_nu_nu_mu(x) or P_nu_nu_-mu(x) can be combined
  //locally after the first transport in the mu direction before transporting
  //twice in the nu direction

  #ifdef CPP
  for(nu = 0; nu < NUM_DIR; nu+=2)
    {
      bzero((char *)tmp_frm[nu],(vol/2)*VECT_LEN*sizeof(IFloat));
      //bzero((char *)tmp_frm[nu+1],(vol/2)*VECT_LEN*sizeof(IFloat));
      for(n = 0; n < vol/2; n++)
	{
	  fp2 = knight_onetwo + n*NUM_DIR*VECT_LEN2;
          fp1 = (IFloat *)(tmp_frm[nu+1]) + n*VECT_LEN;
	  fp0 = (IFloat *)(tmp_frm[nu]) + n*VECT_LEN;

	  for(mu = 0; mu < NUM_DIR; mu+=2)
	    {
	      if(mu/2 != nu/2)
		{
		  fp3 = fp2+mu*VECT_LEN2;
		  fp4 = fp2+(mu+1)*VECT_LEN2;
		  for(k = 0; k < VECT_LEN; k++)
		    {
		      *(fp0+k) += *(fp3 + k);
		      *(fp0+k) -= *(fp4 + k);
		    }
		}
	    }

	  //moveMem(fp1,fp0,VECT_LEN*sizeof(IFloat));
	  //for(k = 0; k < VECT_LEN; k++)
	  //  *(fp1+k) = *(fp0+k);

	}
	  tmp_frm_p[nu] = tmp_frm[nu];
	  tmp_frm_p[nu+1] = tmp_frm[nu];
    }
  #else
  for(nu = 0; nu < 4; nu++)
    {
//      bzero((char *)tmp_frm[nu],(vol/2)*VECT_LEN*sizeof(IFloat));
      res_p[nu] = (Float *)tmp_frm[nu];
    }
  p4_recom_n(vol/2, (Float *)knight_onetwo, res_p);
  for(nu = 0; nu < 4; nu++)
    {
      tmp_frm_p[2*nu] = tmp_frm[nu];
      tmp_frm_p[2*nu+1] = tmp_frm[nu];
    }
  #endif

#ifdef PROFILE
  dtime +=dclock();
  nflops = 4*(vol/2)*5*6;
  print_flops(fname,"recom",nflops,dtime);
  dtime = -dclock();
#endif

  //Transport the resulting combinations in the nu direction
  pt->run(8,tmp_frm2,tmp_frm_p,dir,(ChkbType) (1-cb));

#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT2",ParTrans::PTflops,dtime);
  ParTrans::PTflops = 0;
  dtime = -dclock();
#endif
  pt->run(8,knight_onetwo,tmp_frm2,dir,(ChkbType) cb, 1);

#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT2",ParTrans::PTflops,dtime);
  ParTrans::PTflops=0;
  dtime =-dclock();
#endif

  //Do linear combination for the initial two-link terms
  //
  //Terms of the form P_mu_nu_nu(x) or P_mu_-nu_-nu(x) can be combined
  //locally after the first transport in the nu direction before transporting
  //mu direction

  #ifdef CPP
  for(mu = 0; mu < NUM_DIR; mu+=2)
    {
      bzero((char *)tmp_frm[mu],(vol/2)*VECT_LEN*sizeof(IFloat));
      //bzero((char *)tmp_frm[mu+1],(vol/2)*VECT_LEN*sizeof(IFloat));
      for(n = 0; n < vol/2; n++)
	{
	  fp2 = knight_twoone + n*NUM_DIR*VECT_LEN2;
	  fp1 = (IFloat *)(tmp_frm[mu+1])+n*VECT_LEN;
	  fp0 = (IFloat *)(tmp_frm[mu]) + n*VECT_LEN;

	  for(nu = 0; nu < NUM_DIR; nu++)
	    if(mu/2 != nu/2)
	      {
		fp3 = fp2+nu*VECT_LEN2;
		for(k = 0; k < VECT_LEN; k++)
		  *(fp0+k) += *(fp3 + k);
	      }

	  //moveMem(fp1,fp0,VECT_LEN*sizeof(IFloat));
	  //for(k = 0; k < VECT_LEN; k++)
	  //  *(fp1+k) = *(fp0+k);

      }
	  tmp_frm_p[mu] = tmp_frm[mu];
	  tmp_frm_p[mu+1] = tmp_frm[mu];
    }
  #else
  for(nu = 0; nu < 4; nu++)
    {
 //     bzero((char *)tmp_frm[nu],(vol/2)*VECT_LEN*sizeof(IFloat));
      res_p[nu] = (Float *)tmp_frm[nu];
    }
  p4_recom(vol/2, (Float *)knight_twoone, res_p);
  for(nu = 0; nu < 4; nu++)
    {
      tmp_frm_p[2*nu] = tmp_frm[nu];
      tmp_frm_p[2*nu+1] = tmp_frm[nu];
    }
  #endif

#ifdef PROFILE
  dtime +=dclock();
  nflops = 4*(vol/2)*5*6;
  print_flops(fname,"recom2",nflops,dtime);
  dtime =-dclock();
#endif

  //Transport these by one link
  pt->run(8,knight_twoone,tmp_frm_p,dir,(ChkbType) cb, 1);

#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT3",ParTrans::PTflops,dtime);
  ParTrans::PTflops=0;
  dtime =-dclock();
#endif

  //Calculate the smeared one link term
  pt->run(8,smeared_onelink,fermion,dir,(ChkbType) cb,1,smeared_gauge);

  //Sum up contributions from the knight's move terms
  //and the smeared one link term and place them in fout


#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT3",ParTrans::PTflops,dtime);
  ParTrans::PTflops=0;
  dtime =-dclock();
#endif


  bzero((char*)f_out,(vol/2)*VECT_LEN*sizeof(IFloat));
  #ifdef CPP
   for(n = 0; n < vol/2; n++)
    {
      fp3 = knight_onetwo + n*NUM_DIR*VECT_LEN2;
      fp2 = smeared_onelink + n*NUM_DIR*VECT_LEN2;
      fp1 = knight_twoone + n*NUM_DIR*VECT_LEN2;
      fp0 = (IFloat *)(f_out) + n*VECT_LEN;

      for(mu = 0; mu < 8; mu+=2)
	{
	  fp5 = fp3+mu*VECT_LEN2;
	  fp6 = fp3+(mu+1)*VECT_LEN2;
	  fp7 = fp1+mu*VECT_LEN2;
	  fp8 = fp2+mu*VECT_LEN2;
	  fp9 = fp1+(mu+1)*VECT_LEN2;
	  fp10 = fp2+(mu+1)*VECT_LEN2;
	  for(k = 0; k < VECT_LEN; k++)
	    {
	      fp4 = fp0 + k;
	      *(fp4) += *(fp5 + k)*c_knight;
	      *(fp4) += *(fp6 + k)*c_knight;
	      *(fp4) += *(fp7 + k)*c_knight;
	      *(fp4) += *(fp8 + k);
	      *(fp4) -= *(fp9 + k)*c_knight;
	      *(fp4) -= *(fp10 + k);
	    }
	}
    }
   #else
   double cknight = (double) c_knight;
   double negcknight = -cknight;
   double one = 1.0;
   double negone = -1.0;
   p4_dsum(vol/2, (Float *)knight_onetwo, (Float *) f_out, &cknight, &cknight);
   p4_dsum(vol/2, (Float *)knight_twoone, (Float *) f_out, &cknight, &negcknight);
   p4_dsum(vol/2, (Float *)smeared_onelink, (Float *) f_out, &one, &negone);
   #endif


#ifdef PROFILE
  dtime +=dclock();
  nflops = (vol/2)*(2*8*2+8)*6;
  print_flops(fname,"dsum",nflops,dtime);
  dtime =-dclock();
#endif
  
#if 0
   #ifdef PROFILE
   dtime += dclock();
   nflops += ParTrans::PTflops;
   #endif
#endif

   //Parallel Transport flops:
   //7 parallel transports * (264*vol) flops/transport = 1848*vol
   //Linear combination flops:
   //60*vol + 60*vol + 120*vol = 240*vol
   //Total flop count = 2088*vol

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
void p4_dMdmu(Vector *f_out, Vector *f_in, int cb, int dag, int order)
{
  //Define the eight different directions.
  //Let mu = dir[i]/2;
  //If dir[i] is even, then this transport in the negative mu direction
  //If dir[i] is odd, then we have transport in the positive mu direction
  char *fname = "p4_dMdmu()";
  int dir[8] = {0,1,2,3,4,5,6,7};

  IFloat * fp0,*fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8,*fp9,*fp10;
  Vector * tmp_frm_p[8];
  Float * res_p[4];
  int k,n,mu,nu = 0;

  double c_one_up   = 1.0;
  double c_one_down = 1.0;
  double c_two_up   = 1.0;
  double c_two_down = 1.0;

  for(k=0; k<order; k++)
    {
      c_one_down *= -1.0;
      c_two_up   *=  2.0;
      c_two_down *= -2.0;
    }

  //Copy the input fields into fermion
  for(int mu = 0; mu < NUM_DIR; mu++)
    fermion[mu] = f_in;

  int nflops=0;
  Float dtime;
#ifdef PROFILE
  dtime = -dclock();
  ParTrans::PTflops = 0;
#endif

  //First step, gather one link and two link parallel transports from
  //all eight directions

  //Calculation of P_mu(x); result-> knight_onetwo
  pt->run(8,knight_onetwo,fermion,dir,(ChkbType) cb,1);

#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT",ParTrans::PTflops,dtime);
  ParTrans::PTflops = 0;
  dtime = -dclock();
#endif

  //Calculation of P_nu_nu(x); result-> knight_twoone
  pt->run(8,tmp_frm,fermion,dir,(ChkbType) cb);

#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT",ParTrans::PTflops,dtime);
  ParTrans::PTflops = 0;
  dtime = -dclock();
#endif
  pt->run(8,knight_twoone,tmp_frm,dir,(ChkbType) (1-cb),1);

#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT",ParTrans::PTflops,dtime);
  ParTrans::PTflops = 0;
  dtime = -dclock();
#endif

  //Do appropriate linear combinations for the initial one link term
  //
  //Terms of the form P_nu_nu_mu(x) or P_nu_nu_-mu(x) can be combined
  //locally after the first transport in the mu direction before transporting
  //twice in the nu direction

  #ifdef CPP
  for(nu = 0; nu < NUM_DIR; nu+=2)
    {
      bzero((char *)tmp_frm[nu],(vol/2)*VECT_LEN*sizeof(IFloat));
      //bzero((char *)tmp_frm[nu+1],(vol/2)*VECT_LEN*sizeof(IFloat));
      for(n = 0; n < vol/2; n++)
	{
	  fp2 = knight_onetwo + n*NUM_DIR*VECT_LEN2;
	  fp0 = (IFloat *)(tmp_frm[nu]) + n*VECT_LEN;

	  if ( nu < 6 )
	    {
	      fp3 = fp2+6*VECT_LEN2;
	      fp4 = fp2+7*VECT_LEN2;
	      for(k = 0; k < VECT_LEN; k++)
		{
		  *(fp0+k) += *(fp3 + k);
		  *(fp0+k) -= *(fp4 + k) * c_one_down;
		}
	      
	    }
	  else
	    {
	      for(mu = 0; mu < 6; mu+=2)
		{
		  fp3 = fp2+mu*VECT_LEN2;
		  fp4 = fp2+(mu+1)*VECT_LEN2;
		  for(k = 0; k < VECT_LEN; k++)
		    {
		      *(fp0+k) += *(fp3 + k);
		      *(fp0+k) -= *(fp4 + k);
		    }
		}
	    }
	}
      tmp_frm_p[nu] = tmp_frm[nu];
      tmp_frm_p[nu+1] = tmp_frm[nu];
    }
  #else
  // This need to be done !!
  // =======================
  //  for(nu = 0; nu < 4; nu++)
  //  {
  //    bzero((char *)tmp_frm[nu],(vol/2)*VECT_LEN*sizeof(IFloat));
  //    res_p[nu] = (Float *)tmp_frm[nu];
  //  }
  // p4_recom_n(vol/2, (Float *)knight_onetwo, res_p);
  // for(nu = 0; nu < 4; nu++)
  //  {
  //    tmp_frm_p[2*nu] = tmp_frm[nu];
  //    tmp_frm_p[2*nu+1] = tmp_frm[nu];
  //  }
  #endif

#ifdef PROFILE
  dtime +=dclock();
  nflops = 4*(vol/2)*5*6;
  print_flops(fname,"recom",nflops,dtime);
  dtime = -dclock();
#endif

  //Transport the resulting combinations in the nu direction
  pt->run(8,tmp_frm2,tmp_frm_p,dir,(ChkbType) (1-cb));

#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT2",ParTrans::PTflops,dtime);
  ParTrans::PTflops = 0;
  dtime = -dclock();
#endif
  pt->run(8,knight_onetwo,tmp_frm2,dir,(ChkbType) cb, 1);

#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT2",ParTrans::PTflops,dtime);
  ParTrans::PTflops=0;
  dtime =-dclock();
#endif

  //Do linear combination for the initial two-link terms
  //
  //Terms of the form P_mu_nu_nu(x) or P_mu_-nu_-nu(x) can be combined
  //locally after the first transport in the nu direction before transporting
  //mu direction

  #ifdef CPP
  for(mu = 0; mu < NUM_DIR; mu+=2)
    {
      bzero((char *)tmp_frm[mu],(vol/2)*VECT_LEN*sizeof(IFloat));
      //bzero((char *)tmp_frm[mu+1],(vol/2)*VECT_LEN*sizeof(IFloat));
      for(n = 0; n < vol/2; n++)
	{
	  fp2 = knight_twoone + n*NUM_DIR*VECT_LEN2;
	  fp0 = (IFloat *)(tmp_frm[mu]) + n*VECT_LEN;

	  if ( mu < 6 ) 
	    {
	      fp3 = fp2+6*VECT_LEN2;
	      fp4 = fp2+7*VECT_LEN2;
	      for(k = 0; k < VECT_LEN; k++)
		{
		  *(fp0+k) += *(fp3 + k) * c_two_up;
		  *(fp0+k) += *(fp4 + k) * c_two_down;
		}
	    }
	  else
	    {
	      for(nu = 0; nu < 6; nu+=2)
		{
		  fp3 = fp2+nu*VECT_LEN2;
		  fp4 = fp2+(nu+1)*VECT_LEN2;
		  for(k = 0; k < VECT_LEN; k++)
		    {
		      *(fp0+k) += *(fp3 + k);
		      *(fp0+k) += *(fp4 + k);
		    }
		}
	    }
	}
      tmp_frm_p[mu] = tmp_frm[mu];
      tmp_frm_p[mu+1] = tmp_frm[mu];
    }

  #else
  // This needs to be done !!
  // ========================
  // for(nu = 0; nu < 4; nu++)
  //   {
  //     bzero((char *)tmp_frm[nu],(vol/2)*VECT_LEN*sizeof(IFloat));
  //     res_p[nu] = (Float *)tmp_frm[nu];
  //   }
  // p4_recom(vol/2, (Float *)knight_twoone, res_p);
  // for(nu = 0; nu < 4; nu++)
  //   {
  //     tmp_frm_p[2*nu] = tmp_frm[nu];
  //     tmp_frm_p[2*nu+1] = tmp_frm[nu];
  //   }
  #endif

#ifdef PROFILE
  dtime +=dclock();
  nflops = 4*(vol/2)*5*6;
  print_flops(fname,"recom2",nflops,dtime);
  dtime =-dclock();
#endif

  //Transport these by one link
  pt->run(8,knight_twoone,tmp_frm_p,dir,(ChkbType) cb, 1);

#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT3",ParTrans::PTflops,dtime);
  ParTrans::PTflops=0;
  dtime =-dclock();
#endif

  //Calculate the smeared one link term
  pt->run(8,smeared_onelink,fermion,dir,(ChkbType) cb,1,smeared_gauge);

  //Sum up contributions from the knight's move terms
  //and the smeared one link term and place them in fout


#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"PT3",ParTrans::PTflops,dtime);
  ParTrans::PTflops=0;
  dtime =-dclock();
#endif


  bzero((char*)f_out,(vol/2)*VECT_LEN*sizeof(IFloat));
  #ifdef CPP

   for(n = 0; n < vol/2; n++)
    {
      fp3 = knight_onetwo   + n*NUM_DIR*VECT_LEN2;
      fp2 = smeared_onelink + n*NUM_DIR*VECT_LEN2;
      fp1 = knight_twoone   + n*NUM_DIR*VECT_LEN2;
      fp0 = (IFloat *)(f_out) + n*VECT_LEN;

      for(mu=0; mu<6; mu+=2)
	{
	  fp5  = fp3 +  mu    * VECT_LEN2;
	  fp6  = fp3 + (mu+1) * VECT_LEN2;
	  fp7  = fp1 +  mu    * VECT_LEN2;
	  fp8  = fp1 + (mu+1) * VECT_LEN2;

	  for(k = 0; k < VECT_LEN; k++)
	    {
	      fp4 = fp0 + k;
	      *(fp4) += *( fp5 + k) * c_knight;
	      *(fp4) += *( fp6 + k) * c_knight;
	      *(fp4) += *( fp7 + k) * c_knight;
	      *(fp4) -= *( fp8 + k) * c_knight;
	    }
	}

      fp5  = fp3 + 6 * VECT_LEN2;
      fp6  = fp3 + 7 * VECT_LEN2;
      fp7  = fp1 + 6 * VECT_LEN2;
      fp8  = fp1 + 7 * VECT_LEN2;
      fp9  = fp2 + 6 * VECT_LEN2;
      fp10 = fp2 + 7 * VECT_LEN2;
      
      for(k = 0; k < VECT_LEN; k++)
	{
	  fp4 = fp0 + k;
	  *(fp4) += *( fp5 + k) * c_knight * c_two_up;
	  *(fp4) += *( fp6 + k) * c_knight * c_two_down;
	  *(fp4) += *( fp7 + k) * c_knight * c_one_up;
	  *(fp4) -= *( fp8 + k) * c_knight * c_one_down;
	  *(fp4) += *( fp9 + k) * c_one_up;
	  *(fp4) -= *(fp10 + k) * c_one_down;
	}
    }
   #else
   // This needs to be done !!
   // ========================
   // double cknight = (double) c_knight;
   // double negcknight = -cknight;
   // double one = 1.0;
   // double negone = -1.0;
   // p4_dsum(vol/2, (Float *)knight_onetwo, (Float *) f_out, &cknight, &cknight);
   // p4_dsum(vol/2, (Float *)knight_twoone, (Float *) f_out, &cknight, &negcknight);
   // p4_dsum(vol/2, (Float *)smeared_onelink, (Float *) f_out, &one, &negone);
   #endif


#ifdef PROFILE
  dtime +=dclock();
  nflops = (vol/2)*(2*8*2+8)*6;
  print_flops(fname,"dsum",nflops,dtime);
  dtime =-dclock();
#endif
  
#if 0
   #ifdef PROFILE
   dtime += dclock();
   nflops += ParTrans::PTflops;
   #endif
#endif

   //Parallel Transport flops:
   //7 parallel transports * (264*vol) flops/transport = 1848*vol
   //Linear combination flops:
   //60*vol + 60*vol + 120*vol = 240*vol
   //Total flop count = 2088*vol

   DiracOp::CGflops += 2088*vol;
}

#undef CPP

CPS_END_NAMESPACE
