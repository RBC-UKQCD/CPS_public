#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Definitiion of Mom class methods.

  $Id: mom.C,v 1.5 2012-08-15 03:45:46 chulwoo Exp $ 
*/
// mom.C
//
// Mom calculates the phase factor for each
// lattice site given a number of momenta and the source parameters
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdlib.h>	// exit()
#include <util/qcdio.h>
#include <util/mom.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/verbose.h> 
#include <util/error.h>
CPS_START_NAMESPACE

Mom MOM;

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
Mom::Mom()
{
  cname = "Mom";
  char *fname = "Mom(MomArg*)";
  VRB.Func(cname,fname);
}
//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Mom::~Mom() {
  char *fname = "~Mom()";
  VRB.Func(cname,fname);

  // de-allocate space for momentum factor
  VRB.Sfree(cname,fname, "mom_fact", mom_fact);
  sfree(mom_fact);
}



//------------------------------------------------------------------
//Initialisation
/*!
  \param arg The set of run-time parameters for the momentum computation.
  \todo Why is this stuff not in the constructor?
*/
//------------------------------------------------------------------
void Mom::Init(MomArg *arg)
{
  char *fname = "Init()";
  VRB.Func(cname,fname);

  PI = 3.14159265358979323846; 

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  mom_arg = arg;

  no_of_mom = mom_arg->no_of_momenta;
  deg       = mom_arg->deg;
  max_p1    = mom_arg->max_p1;
  max_p2    = mom_arg->max_p2;
  max_p3    = mom_arg->max_p3;
  // allocate space for momentum factor table
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
  glb_sour_center[0] = (mom_arg->src_end[0] + mom_arg->src_begin[0])/2;
  glb_sour_center[1] = (mom_arg->src_end[1] + mom_arg->src_begin[1])/2;
  glb_sour_center[2] = (mom_arg->src_end[2] + mom_arg->src_begin[2])/2;
  glb_sour_center[3] = (mom_arg->src_end[3] + mom_arg->src_begin[3])/2;

  // propagation direction
  dir = mom_arg->dir;
  // other directions: i,j,k
  i = (mom_arg->dir + 1)%4;
  j = (mom_arg->dir + 2)%4;
  k = (mom_arg->dir + 3)%4;

}



//------------------------------------------------------------------
/*!
				
This method actually computes all the phases and stores them internally.
\post A table of momenta and their corresponding numbers by which they
can be identified is written to a
file called \c mom_table.log.
*/
//------------------------------------------------------------------
void Mom::run() 
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  int s[4];          // local  coordinates
  int glb[4];        // global coordinates
  int n[4];          // momentum in lattice units e.g. n=(1,0,0,0)
  n[dir]=0;          // momentum in propagation direction always zero 
  Float mom[4];      // momentum in physical units: p_mu = 2*PI*n_mu / L_mu  

  if (deg && (  (max_p1 < max_p2) || (max_p2 < max_p3) ) ) {
    ERR.General(cname, fname, "mom_deg requires ordered momenta %d %d %d \n", max_p1,max_p2,max_p3) ;
  }

  /*
  FILE *fp;
  if( (fp = Fopen("mom_table.log", "a")) == NULL ) {
    ERR.FileA(cname, fname, "mom_table.log");
  }
  Fprintf(fp,"# phase factors for dir = %d deg-flag = %d \n",dir,deg);
  Fprintf(fp,"# index   n[%d]  n[%d]  n[%d] \n",i,j,k);
  */

  int imom=0;
  for (int p1=0; p1<= max_p1; p1++){
    for (int p2=0; p2<= max_p2; p2++){
      if (p2>p1 && deg) continue;
      for (int p3=0; p3<= max_p3; p3++){
	if (p3>p2 && deg) continue;

	n[i]=p1; n[j]=p2; n[k]=p3;  // momentum in natural units
	
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

		} else {
		  // take only the first momentum choice
      
		  mom_fact[offset]=Complex(  cos(alpha[0]) , -sin(alpha[0]) );
		} // end if (deg) { } else { }

	      } // s[i]
	    } // s[j]
	  } // s[k]
	} // s[dir]
	// Fprintf(fp,"%d         %d      %d      %d \n",imom, p1,p2,p3);
	imom++; // calculate next momentum

      } // p3
    } // p2
  } // p1

  if (no_of_mom != imom)
    ERR.General(cname, fname, "conflict between imom=%d and no_of_mom = %d",imom,no_of_mom);
  
  /*
    Fprintf(fp,"total number of momenta = %d \n",no_of_mom);
    Fprintf(fp,"=================================== \n");
    Fclose(fp);
  */
  // finished determination of mom_fact[offset] ========
}


//------------------------------------------------------------------
// fact(int imom, int *s)
// returns the complex phase factor for momentum "imom" at site "s"
/*!
  \param imom The number referring to the required momentum (a table of
  momenta and their corresponding numbers is written to a file called
  \c mom_table.log by ::run.
  \param s The coordinates of the lattice site at which the phase is
  required.
  \return The complex phase.
*/	  
//------------------------------------------------------------------
Complex Mom::fact(int imom, int *s) 
{
  char *fname = "fact()";
  VRB.Func(cname,fname);

  int offset = imom + s[0]*no_of_mom + s[1]*no_of_mom*nx[0] + s[2]*no_of_mom*nx[0]*nx[1] + s[3]*no_of_mom*nx[0]*nx[1]*nx[2];

  return mom_fact[offset];


}

CPS_END_NAMESPACE
