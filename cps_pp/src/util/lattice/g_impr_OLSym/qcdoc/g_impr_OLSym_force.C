//Force term for the Symanzik-improved gauge action
//Takes advantage of lattice-wide operations to speed HMC evolution

#include <config.h>
#include <math.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pt.h>
#include <util/time.h>
CPS_START_NAMESPACE

//-----------------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float step_size):
// It evolves the canonical momentum mom by step_size
// using the pure gauge force.
//-----------------------------------------------------------------------------
#define PROFILE
void GimprOLSym::EvolveMomGforce(Matrix *mom, Float step_size)
{
  char *fname = "EvolveMomGforce(M*,F)";  //Name of our function
  VRB.Func(cname,fname);                  //Sets name of function


  static Matrix mt0;
  static Matrix *mp0 = &mt0;

  //Sets some kind of flop counter
#ifdef PROFILE
  Float time = -dclock();
  ForceFlops = 0;
  ParTrans::PTflops=0;
#endif

  static int vol = GJP.VolNodeSites();  //Local lattice volume
  const int N = 4;                      //Num of dimensions
  Float tmp_plaq = minus_beta_over_3;
  Float tmp_rect = minus_beta_over_3*rect_coeff;
  Float tmp_cube = minus_beta_over_3*cube_coeff;

  //Pointer to block of unit matrices for each lattice site
  Matrix *Unit = (Matrix *) fmalloc(vol*sizeof(Matrix));

  //Set all these matrices to the identity
  for(int i = 0; i < vol;i++)
    Unit[i] = 1.;

  //Temporary matrices for use in calculation of the staples
  Matrix *tmp1[N];
  Matrix *tmp2[N];
  Matrix *tmp3[N];
  Matrix *tmp4[N];

  for(int i = 0;i<N;i++)
    {
      tmp1[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
      tmp2[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
      tmp3[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
      tmp4[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
      //Set all bytes in tmp3
      bzero(tmp3[i],vol*sizeof(Matrix));
    }

  //Holds the sum of staples associated with each lattice site
  Matrix *result[4];
  for(int i = 0;i<4;i++){
    result[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
  }

  //Array of four pointers to a block of unit matrices
  Matrix *Units[4];
  for(int i = 0; i < N;i++)
    Units[i] = Unit;

  int mu,nu,rho;
  {
    int dirs_p[] = {0,2,4,6,0,2,4};   //Positive directions
    int dirs_m[] = {1,3,5,7,1,3,5};   //Negative directions

    //Instantiate parallel transporter
    ParTransGauge pt(*this);

    //Take a vector field V(x).  Suppose we wish to parallel transport
    //this vector field from the point x to the point x-mu.  That is,
    //we wish to find a new field V'(x-mu) in terms of V(x) such that
    //the new field V'(x-mu) changes appropriately under gauge transformation.
    //
    //The parallel transport that satisfies this property is
    //V'(x-mu) = U_mu(x-mu)V(x).  This combination transforms like
    //a vector field.
    //
    //When using the parallel transport class, one must specify an intial
    //field, a direction, and a final field.  However, the direction
    //that one must specify is the direction of the link matrix, not the
    //direction of the parallel transport.  Thus, to calculate V' as above,
    //use:
    //
    //pt.run(1,V',V,dir_Plus_mu)
    //
    //The new vector field will be indexed according to its new position.
    //
    //The staples for the rectangle and chair terms can be seen as
    //modifications of the standard plaquette.
    //
    //Consider the standard plaquette staple.  Instead of the normal three
    //links that constitute a staple, remove a single link and replace
    //it with some superposition of the original link and its possible
    //staples.  These extra terms are equivalent to the rectangle and
    //the chair terms.
    //
    //In order to minimize the number of matrix multiplications, this
    //philosophy will be used to calculate the staple, instead
    //of calculating the 6+18+72 different terms per link sequentially

      for(nu = 1;nu<4;nu++){
	//First, consider replacing the link U_nu(x)_dagger

	//Plaquette term
	pt.run(N,tmp1,Units,dirs_m+nu);
	for(int i = 0; i<N;i++)
	  {
	    vecTimesEquFloat((IFloat *)tmp1[i],tmp_plaq,vol*18);
	    moveMem(result[i],tmp1[i],vol*sizeof(Matrix));
	  }

	//Rectangle term
	pt.run(N,tmp1,Units,dirs_p);
	pt.run(N,tmp2,tmp1,dirs_m+nu);
	pt.run(N,tmp1,tmp2,dirs_m);
	for(int i = 0; i<N;i++)
	  vaxpy3_m(result[i],&tmp_rect,tmp1[i],result[i],vol*3);

	//Chair terms
	for(rho = 2;rho<4;rho++)
	  {
	    pt.run(N,tmp1,Units,dirs_p+rho);
	    pt.run(N,tmp2,tmp1,dirs_m+nu);
	    pt.run(N,tmp1,tmp2,dirs_m+rho);
	    for(int i = 0; i<N;i++)
	      vaxpy3_m(result[i],&tmp_cube,tmp1[i],result[i],vol*3);
	    pt.run(N,tmp1,Units,dirs_m+rho);
	    pt.run(N,tmp2,tmp1,dirs_m+nu);
	    pt.run(N,tmp1,tmp2,dirs_p+rho);
	    for(int i = 0; i<N;i++)
	      vaxpy3_m(result[i],&tmp_cube,tmp1[i],result[i],vol*3);
	  }

	//Now calculate the rest of the plaquette for this linear combination
	pt.run(N,tmp2,result,dirs_m);
	pt.run(N,result,tmp2,dirs_p+nu);

	//Add this result to the sum of the staples
	for(int i = 0; i < N; i++)
	  vecAddEquVec((IFloat *)tmp3[i],(IFloat *)result[i],vol*18);

	//Next, consider replacing the link U_mu(x+nu)_dagger
	pt.run(N,tmp4,Units,dirs_m+nu);
	//Plaquette term already calculated

	//Rectangle term
	pt.run(N,tmp1,tmp4,dirs_m+nu);
	pt.run(N,tmp2,tmp1,dirs_m);
	pt.run(N,tmp1,tmp2,dirs_p+nu);
	for(int i = 0; i<N;i++)
	  {
	    vecTimesEquFloat((IFloat *)tmp1[i],tmp_rect,vol*18);
	    moveMem(result[i],tmp1[i],vol*sizeof(Matrix));
	  }


	//Chair terms
	for(rho = 2;rho<4;rho++)
	  {
	    pt.run(N,tmp1,tmp4,dirs_p+rho);
	    pt.run(N,tmp2,tmp1,dirs_m);
	    pt.run(N,tmp1,tmp2,dirs_m+rho);
	    for(int i = 0; i<N;i++)
	      vaxpy3_m(result[i],&tmp_cube,tmp1[i],result[i],vol*3);
	    pt.run(N,tmp1,tmp4,dirs_m+rho);
	    pt.run(N,tmp2,tmp1,dirs_m);
	    pt.run(N,tmp1,tmp2,dirs_p+rho);
	    for(int i = 0; i<N;i++)
	      vaxpy3_m(result[i],&tmp_cube,tmp1[i],result[i],vol*3);
	  }

	//Now calculate the rest of the plaquette for this linear combination
	pt.run(N,tmp2,result,dirs_p+nu);

	//Add this result to the sum of the staples
	for(int i = 0; i < N; i++)
	  vecAddEquVec((IFloat *)tmp3[i],(IFloat *)tmp2[i],vol*18);

	//Next, consider replacing the link U_nu(x+mu)
	pt.run(N,tmp2,Units,dirs_m+nu);
	pt.run(N,tmp4,tmp2,dirs_m);
	//Plaquette term already calculated

	//Rectangle term
	pt.run(N,tmp1,tmp4,dirs_m);
	pt.run(N,tmp2,tmp1,dirs_p+nu);
	pt.run(N,tmp1,tmp2,dirs_p);
	for(int i = 0; i<N;i++)
	  {
	    vecTimesEquFloat((IFloat *)tmp1[i],tmp_rect,vol*18);
	    moveMem(result[i],tmp1[i],vol*sizeof(Matrix));
	  }

	//Chair terms
	for(rho = 2;rho<4;rho++)
	  {
	    pt.run(N,tmp1,tmp4,dirs_p+rho);
	    pt.run(N,tmp2,tmp1,dirs_p+nu);
	    pt.run(N,tmp1,tmp2,dirs_m+rho);
	    for(int i = 0; i<N;i++)
	      vaxpy3_m(result[i],&tmp_cube,tmp1[i],result[i],vol*3);
	    pt.run(N,tmp1,tmp4,dirs_m+rho);
	    pt.run(N,tmp2,tmp1,dirs_p+nu);
	    pt.run(N,tmp1,tmp2,dirs_p+rho);
	    for(int i = 0; i<N;i++)
	      vaxpy3_m(result[i],&tmp_cube,tmp1[i],result[i],vol*3);
	  }

	//Add this result to the sum of the staples
	for(int i = 0; i < N; i++)
	  vecAddEquVec((IFloat *)tmp3[i],(IFloat *)result[i],vol*18);

	//-------------------------------------------------------------------------
	//Now, for the plaquette in the negative nu direction
	//First, consider replacing the link U_nu(x-nu)

	//Plaquette term
	pt.run(N,tmp1,Units,dirs_p+nu);
	for(int i = 0; i<N;i++)
	  {
	    vecTimesEquFloat((IFloat *)tmp1[i],tmp_plaq,vol*18);
	    moveMem(result[i],tmp1[i],vol*sizeof(Matrix));
	  }

	//Rectangle term
	pt.run(N,tmp1,Units,dirs_p);
	pt.run(N,tmp2,tmp1,dirs_p+nu);
	pt.run(N,tmp1,tmp2,dirs_m);
	for(int i = 0; i<N;i++)
	  vaxpy3_m(result[i],&tmp_rect,tmp1[i],result[i],vol*3);

	//Chair terms
	for(rho = 2;rho<4;rho++)
	  {
	    pt.run(N,tmp1,Units,dirs_p+rho);
	    pt.run(N,tmp2,tmp1,dirs_p+nu);
	    pt.run(N,tmp1,tmp2,dirs_m+rho);
	    for(int i = 0; i<N;i++)
	      vaxpy3_m(result[i],&tmp_cube,tmp1[i],result[i],vol*3);
	    pt.run(N,tmp1,Units,dirs_m+rho);
	    pt.run(N,tmp2,tmp1,dirs_p+nu);
	    pt.run(N,tmp1,tmp2,dirs_p+rho);
	    for(int i = 0; i<N;i++)
	      vaxpy3_m(result[i],&tmp_cube,tmp1[i],result[i],vol*3);
	  }

	//Now calculate the rest of the plaquette for this linear combination
	pt.run(N,tmp2,result,dirs_m);
	pt.run(N,result,tmp2,dirs_m+nu);

	//Add this result to the sum of the staples
	for(int i = 0; i < N; i++)
	  vecAddEquVec((IFloat *)tmp3[i],(IFloat *)result[i],vol*18);

	//Next, consider replacing the link U_mu(x-nu)_dagger
	pt.run(N,tmp4,Units,dirs_p+nu);
	//Plaquette term already calculated

	//Rectangle term
	pt.run(N,tmp1,tmp4,dirs_p+nu);
	pt.run(N,tmp2,tmp1,dirs_m);
	pt.run(N,tmp1,tmp2,dirs_m+nu);
	for(int i = 0; i<N;i++)
	  {
	    vecTimesEquFloat((IFloat *)tmp1[i],tmp_rect,vol*18);
	    moveMem(result[i],tmp1[i],vol*sizeof(Matrix));
	  }

	//Chair terms
	for(rho = 2;rho<4;rho++)
	  {
	    pt.run(N,tmp1,tmp4,dirs_p+rho);
	    pt.run(N,tmp2,tmp1,dirs_m);
	    pt.run(N,tmp1,tmp2,dirs_m+rho);
	    for(int i = 0; i<N;i++)
	      vaxpy3_m(result[i],&tmp_cube,tmp1[i],result[i],vol*3);
	    pt.run(N,tmp1,tmp4,dirs_m+rho);
	    pt.run(N,tmp2,tmp1,dirs_m);
	    pt.run(N,tmp1,tmp2,dirs_p+rho);
	    for(int i = 0; i<N;i++)
	      vaxpy3_m(result[i],&tmp_cube,tmp1[i],result[i],vol*3);
	  }

	//Now calculate the rest of the plaquette for this linear combination
	pt.run(N,tmp2,result,dirs_m+nu);

	//Add this result to the sum of the staples
	for(int i = 0; i < N; i++)
	  vecAddEquVec((IFloat *)tmp3[i],(IFloat *)tmp2[i],vol*18);

	//Next, consider replacing the link U_nu(x+mu-nu)_dagger
	pt.run(N,tmp2,Units,dirs_p+nu);
	pt.run(N,tmp4,tmp2,dirs_m);
	//Plaquette term already calculated


	//Rectangle term
	pt.run(N,tmp1,tmp4,dirs_m);
	pt.run(N,tmp2,tmp1,dirs_m+nu);
	pt.run(N,tmp1,tmp2,dirs_p);
	for(int i = 0; i<N;i++)
	  {
	    vecTimesEquFloat((IFloat *)tmp1[i],tmp_rect,vol*18);
	    moveMem(result[i],tmp1[i],vol*sizeof(Matrix));
	  }

	//Chair terms
	for(rho = 2;rho<4;rho++)
	  {
	    pt.run(N,tmp1,tmp4,dirs_p+rho);
	    pt.run(N,tmp2,tmp1,dirs_m+nu);
	    pt.run(N,tmp1,tmp2,dirs_m+rho);
	    for(int i = 0; i<N;i++)
	      vaxpy3_m(result[i],&tmp_cube,tmp1[i],result[i],vol*3);
	    pt.run(N,tmp1,tmp4,dirs_m+rho);
	    pt.run(N,tmp2,tmp1,dirs_m+nu);
	    pt.run(N,tmp1,tmp2,dirs_p+rho);
	    for(int i = 0; i<N;i++)
	      vaxpy3_m(result[i],&tmp_cube,tmp1[i],result[i],vol*3);
	  }

	//Add this result to the sum of the staples
	for(int i = 0; i < N; i++)
	  vecAddEquVec((IFloat *)tmp3[i],(IFloat *)result[i],vol*18);

#if 0

	//First calculate the staple in the positive nu direction
        pt.run(N,tmp1,Units,dirs_m+nu);
  	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_p+nu);

	//tmp2 contains the sum of the staples for a given link
	for(int i = 0; i<N;i++)
	vaxpy3_m(tmp2[i],&tmp_plaq,tmp1[i],tmp2[i],vol*3);

	//Calculating one rectangular staple
	pt.run(N,tmp1,Units,dirs_p);
	pt.run(N,result,tmp1,dirs_m+nu);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_p+nu);

	for(int i = 0; i<N;i++)
	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);

	//Calculating another rectangular staple;
	pt.run(N,tmp1,Units,dirs_m+nu);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_p+nu);
	pt.run(N,tmp1,result,dirs_p);

	for(int i = 0; i<N;i++)
	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);

	//Calculating another rectangular staple;
	pt.run(N,tmp1,Units,dirs_m+nu);
	pt.run(N,result,tmp1,dirs_m+nu);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_p+nu);
	pt.run(N,tmp1,result,dirs_p+nu);

	for(int i = 0; i<N;i++)
	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);

	//Calculating the staple in the negative nu direction
	pt.run(N,tmp1,Units,dirs_p+nu);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_m+nu);

	//Add this result into tmp2
	for(int i = 0; i<N;i++)
	vaxpy3_m(tmp2[i],&tmp_plaq,tmp1[i],tmp2[i],vol*3);

	//Calculating one rectangular staple
	pt.run(N,tmp1,Units,dirs_p);
	pt.run(N,result,tmp1,dirs_p+nu);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_m+nu);

	for(int i = 0; i<N;i++)
	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);

	//Calculating another rectangular staple;
	pt.run(N,tmp1,Units,dirs_p+nu);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_m+nu);
	pt.run(N,tmp1,result,dirs_p);

	for(int i = 0; i<N;i++)
	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);

	//Calculating another rectangular staple;
	pt.run(N,tmp1,Units,dirs_p+nu);
	pt.run(N,result,tmp1,dirs_p+nu);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_m+nu);
	pt.run(N,tmp1,result,dirs_m+nu);

	for(int i = 0; i<N;i++)
	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);
#endif

	//Count the flops?
        ForceFlops +=vol*288*N;
      }
      //Multiply on the left by our original link matrix to get force term
      pt.run(N,result,tmp3,dirs_p);
  }

      Matrix mp1;
      for(mu = 0; mu<4;mu++)
	{
	  Matrix *mtmp = result[mu];
	  // Takes TrLessAntiHerm part to get the force
	  for(int i = 0; i<vol;i++)
	    {
	      #if 1
	      mp1.Dagger((IFloat *)mtmp);
	      mtmp->TrLessAntiHermMatrix(mp1);
	      #else
	      mtmp->TrLessAntiHermMatrix();
	      #endif
	      mtmp++;
	    }
	}
      ForceFlops += vol*60;

  int x[4];
  
  for(x[0] = 0; x[0] < GJP.XnodeSites(); ++x[0]) {
    for(x[1] = 0; x[1] < GJP.YnodeSites(); ++x[1]) {
      for(x[2] = 0; x[2] < GJP.ZnodeSites(); ++x[2]) {
	for(x[3] = 0; x[3] < GJP.TnodeSites(); ++x[3]) {
	  
	  //Calculates the array index offset for the gauge links
	  //located at lattice point x
	  int uoff = GsiteOffset(x);
	  
	  for (int mu = 0; mu < 4; ++mu) {
//	    GforceSite(*mp0, x, mu);   //Old force calculation
	    
	    IFloat *ihp = (IFloat *)(mom+uoff+mu);  //The gauge momentum
	    IFloat *dotp = (IFloat *)mp0;
	    IFloat *dotp2 = (IFloat *) (result[mu]+(uoff/4));
	    fTimesV1PlusV2(ihp, step_size, dotp2, ihp, 18);  //Update the gauge momentum
	  }
	}
      }
    }
  }
  ForceFlops += vol*144;
#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops+ParTrans::PTflops,time);
#endif

  //Free some memory
  ffree(Unit);
  for(int i = 0;i<N;i++){
  ffree(tmp1[i]);
  ffree(tmp2[i]);
  ffree(tmp3[i]);
  ffree(tmp4[i]);
  }
  for(int i = 0;i<4;i++) 
  ffree(result[i]);
}
CPS_END_NAMESPACE
