#include<config.h>
CPS_START_NAMESPACE
// alg_mom.C
//
// AlgMom calculates the phase factor for each
// lattice site given a number of momenta and the source parameters
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdlib.h>	// exit()
#include <util/qcdio.h>
#include <alg/alg_mom.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/gjp.h>
#include <alg/common_arg.h>
#include <util/smalloc.h>
#include <util/vector.h>           // Tr() and ReTr()
#include <util/verbose.h> 
#include <util/error.h>
#include <comms/glb.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgMom::AlgMom( CommonArg *c_arg, MomArg *arg)
{
  cname = "AlgMom";
  char *fname = "AlgMom(CommonArg*, MomArg*)";
  VRB.Func(cname,fname);

  PI = 3.14159265358979323846; 

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_mom_arg = arg;

  no_of_mom = alg_mom_arg->no_of_momenta;
  deg       = alg_mom_arg->deg;

  // allocate space for momentum factor
  int mom_fact_size = no_of_mom*GJP.VolNodeSites();

  mom_fact = (Complex *) smalloc(mom_fact_size * sizeof(Complex));
  if(mom_fact == 0)
    ERR.Pointer(cname,fname, "mom_fact");
  VRB.Smalloc(cname,fname, "mom_fact", mom_fact, mom_fact_size * sizeof(Complex));

  // local lattice extent
  nx[0] = GJP.XnodeSites();
  nx[1] = GJP.YnodeSites();
  nx[2] = GJP.ZnodeSites();
  nx[3] = GJP.TnodeSites();

  // global lattice extent
  glb_L[0] = GJP.XnodeSites()*GJP.Xnodes();
  glb_L[1] = GJP.YnodeSites()*GJP.Ynodes();
  glb_L[2] = GJP.ZnodeSites()*GJP.Znodes();
  glb_L[3] = GJP.TnodeSites()*GJP.Tnodes();

  // global source location
  glb_sour_center[0] = (alg_mom_arg->src_end[0]+alg_mom_arg->src_begin[0])/2;
  glb_sour_center[1] = (alg_mom_arg->src_end[1]+alg_mom_arg->src_begin[1])/2;
  glb_sour_center[2] = (alg_mom_arg->src_end[2]+alg_mom_arg->src_begin[2])/2;
  glb_sour_center[3] = (alg_mom_arg->src_end[3]+alg_mom_arg->src_begin[3])/2;

  // printf(" source_a = %d %d %d %d \n",src_begin[0],src_begin[1],src_begin[2],src_begin[3]);
  // printf(" source_b = %d %d %d %d \n",src_end[0],src_end[1],src_end[2],src_end[3]);
  // printf(" center   = %d %d %d %d \n",glb_sour_center[0],glb_sour_center[1],glb_sour_center[2],glb_sour_center[3]);

  // propagation direction
  dir = alg_mom_arg->dir;
  // other directions: i,j,k
  i = (alg_mom_arg->dir + 1)%4;
  j = (alg_mom_arg->dir + 2)%4;
  k = (alg_mom_arg->dir + 3)%4;

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgMom::~AlgMom() {
  char *fname = "~AlgMom()";
  VRB.Func(cname,fname);

  // de-allocate space for momentum factor
  VRB.Sfree(cname,fname, "mom_fact", mom_fact);
  sfree(mom_fact);
}





//------------------------------------------------------------------
// run()
//
//------------------------------------------------------------------
void AlgMom::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  int s[4];          // local  coordinates
  int glb[4];        // global coordinates
  int n[4];          // momentum in lattice units e.g. n=(1,0,0,0)
  n[dir]=0;          // momentum in propagation direction always zero 
  Float mom[4];      // momentum in physical units: p_mu = 2*PI*n_mu / L_mu  


  int imom;
  for (imom=0; imom < no_of_mom; imom++){

    // first define momentum in lattice units
    if (deg) {
      switch(imom) {
      case 0: n[i]=0; n[j]=0; n[k]=0;break;
	
      case 1: n[i]=1; n[j]=0; n[k]=0;break;
      case 2: n[i]=1; n[j]=1; n[k]=0;break;
      case 3: n[i]=1; n[j]=1; n[k]=1;break;

      case 4: n[i]=2; n[j]=0; n[k]=0;break;
      case 5: n[i]=2; n[j]=1; n[k]=0;break;
      case 6: n[i]=2; n[j]=1; n[k]=1;break;

      case 7: n[i]=2; n[j]=2; n[k]=0;break;
      case 8: n[i]=2; n[j]=2; n[k]=1;break;
      case 9: n[i]=2; n[j]=2; n[k]=2;break;
      default: n[i]=0; n[j]=0; n[k]=0;break;  // zero-momentum for default

      } // end switch(imom) 

    } else {

      // if (!deg) use different way of indexing momenta !

      switch(imom) {
      case 0: n[i]=0; n[j]=0; n[k]=0;break;
	
	// (1,0,0) and permutations
      case 1: n[i]=1; n[j]=0; n[k]=0;break;
      case 2: n[i]=0; n[j]=1; n[k]=0;break;
      case 3: n[i]=0; n[j]=0; n[k]=1;break;
	

	// (1,1,0) and permutations
      case 4: n[i]=1; n[j]=1; n[k]=0;break;
      case 5: n[i]=1; n[j]=0; n[k]=1;break;
      case 6: n[i]=0; n[j]=1; n[k]=1;break;
		

	// (1,1,1)
      case 7: n[i]=1; n[j]=1; n[k]=1;break;


	// (2,0,0) and permutations
      case 8:  n[i]=2; n[j]=0; n[k]=0;break;
      case 9:  n[i]=0; n[j]=2; n[k]=0;break;
      case 10: n[i]=0; n[j]=0; n[k]=2;break;

	// (2,1,0) and permutations
      case 11: n[i]=2; n[j]=1; n[k]=0;break;
      case 12: n[i]=2; n[j]=0; n[k]=1;break;
      case 13: n[i]=1; n[j]=2; n[k]=0;break;
      case 14: n[i]=1; n[j]=0; n[k]=2;break;
      case 15: n[i]=0; n[j]=2; n[k]=1;break;
      case 16: n[i]=0; n[j]=1; n[k]=2;break;
	
	// (2,1,1) and permutations
      case 17:  n[i]=2; n[j]=1; n[k]=1;break;
      case 18:  n[i]=1; n[j]=2; n[k]=1;break;
      case 19:  n[i]=1; n[j]=1; n[k]=2;break;
	
	
	// (2,2,0) and permutations
      case 20: n[i]=2; n[j]=2; n[k]=0;break;
      case 21: n[i]=2; n[j]=0; n[k]=2;break;
      case 22: n[i]=0; n[j]=2; n[k]=2;break;
		
		
	// (2,2,1) and permutations
      case 23: n[i]=2; n[j]=2; n[k]=1;break;
      case 24: n[i]=2; n[j]=1; n[k]=2;break;
      case 25: n[i]=1; n[j]=2; n[k]=2;break;
		

	// (2,2,2)
      case 26: n[i]=2; n[j]=2; n[k]=2;break;


	// (3,0,0)
      case 27: n[i]=3; n[j]=0; n[k]=0;break;
      case 28: n[i]=0; n[j]=3; n[k]=0;break;
      case 29: n[i]=0; n[j]=0; n[k]=3;break;
	

      default: n[i]=0; n[j]=0; n[k]=0;break;  // zero-momentum for default
	
      } // end switch(imom)

    } // end if (deg) {  } else { }  

    // all momenta n are defined now 

    // convert lattice momenta --> physical: p[mu] = 2 Pi n[mu] / L[mu] 
    for (int mu=0; mu<4; mu++) mom[mu] = 2.*PI*n[mu]/glb_L[mu];
    // printf("mom = %e %e %e %e \n",mom[0],mom[1],mom[2],mom[3]);

    // start loop over the local lattice volume
    for (s[dir] = 0; s[dir] < nx[dir]; s[dir]++) {
      for (s[i] = 0; s[i] < nx[i]; s[i]++) {
	for (s[j] = 0; s[j] < nx[j]; s[j]++) {
	  for (s[k] = 0; s[k] < nx[k]; s[k]++) {
	  

	    glb[0] = s[0] + GJP.XnodeSites()*GJP.XnodeCoor();
	    glb[1] = s[1] + GJP.YnodeSites()*GJP.YnodeCoor();
	    glb[2] = s[2] + GJP.ZnodeSites()*GJP.ZnodeCoor();
	    glb[3] = s[3] + GJP.TnodeSites()*GJP.TnodeCoor();

	    int x[4];
	    // global shifts from origin
	    for (int mu=0; mu<4; mu++) x[mu] = glb[mu]-glb_sour_center[mu];
  
	    // offset for mom_fact(imom,s)
	    int offset = imom + s[0]*no_of_mom + s[1]*no_of_mom*nx[0] + s[2]*no_of_mom*nx[0]*nx[1] + s[3]*no_of_mom*nx[0]*nx[1]*nx[2];

	    // mom-factor exp [ -i*alpha ] = cos(alpha) - i sin(alpha)
	    Float alpha[3]; 
	    alpha[0] = x[i]*mom[i] + x[j]*mom[j] + x[k]*mom[k];


	    // printf("s = %d %d %d %d offset = %d alpha = %e \n",s[0],s[1],s[2],s[3],offset,alpha[0]);
	    if (deg) {      
	      // to be used with SYMMETRIC LATTICES ONLY
	      // i.e. calculate also the equivalent 
	      // alpha[1] = x*py + y*pz + z*px
	      // alpha[2] = x*pz + y*px + z*py

	      alpha[1] = x[i]*mom[j] + x[j]*mom[k] + x[k]*mom[i];
	      alpha[2] = x[i]*mom[k] + x[j]*mom[i] + x[k]*mom[j];

	      // sum and normalise with respect to degenrate number of momenta
	      mom_fact[offset]=Complex( ( cos(alpha[0])+cos(alpha[1])+cos(alpha[2]))/3.0, (-sin(alpha[0])-sin(alpha[1])-sin(alpha[2]))/3.0 );
//	      mom_fact[offset].real( ( cos(alpha[0])+cos(alpha[1])+cos(alpha[2]))/3.0 );
//	      mom_fact[offset].imag( (-sin(alpha[0])-sin(alpha[1])-sin(alpha[2]))/3.0 );

	    } else {
	      // take only the first momentum choice
      
	      mom_fact[offset]=Complex(  cos(alpha[0]) , -sin(alpha[0]) );
	    } // end if (deg) { } else { }

	  } // s[3]
	} // s[2]
      } // s[1]
    } // s[0]
  } // imom
  
  // finished determination of mom_fact[offset] ========

}


//------------------------------------------------------------------
// fact(int imom, int *s)
// returns the complex phase factor for momentum "imom" at site "s"
//------------------------------------------------------------------
Complex AlgMom::fact(int imom, int *s) 
{
  char *fname = "fact()";
  VRB.Func(cname,fname);

  int offset = imom + s[0]*no_of_mom + s[1]*no_of_mom*nx[0] + s[2]*no_of_mom*nx[0]*nx[1] + s[3]*no_of_mom*nx[0]*nx[1]*nx[2];

  return mom_fact[offset];


}

CPS_END_NAMESPACE
