//Hybrid Monte Carlo force term for an improved gauge action with plaquette
//and rectangle terms.  Takes advantage of lattice-wide operations to speed
//evolution of HMC trajectory.

#include <config.h>
#include <math.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pt.h>
#include <util/time_cps.h>
CPS_START_NAMESPACE

//? I don't know what this does
#define PROFILE
//------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float dt):
// It evolves the canonical momentum mom by dt
// using the pure gauge force.
//------------------------------------------------------------------

ForceArg GimprRect::EvolveMomGforce(Matrix *mom, Float dt)
{
  char *fname = "EvolveMomGforce(M*,F)mod";  //Name of our function
  VRB.Func(cname,fname);                     //Sets name of the function

  Float L1=0.0;
  Float L2=0.0;
  Float Linf=0.0;

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
  Float tmp_plaq = plaq_coeff;          //1-8*c_1
  Float tmp_rect = rect_coeff;          //c_1

  //Pointer to block of unit matrices for each lattice site
  Matrix *Unit = (Matrix *) fmalloc(vol*sizeof(Matrix));

  //Set all these matrices to the identity
  for(int i = 0; i < vol;i++)
    Unit[i] = 1.;

  //Temporary matrices for use in calculation of the staples
  Matrix *tmp1[N];
  Matrix *tmp2[N];

  for(int i = 0;i<N;i++)
    {
      tmp1[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
      tmp2[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
      //Set all bytes in tmp2 to zero
      bzero((char*)tmp2[i],vol*sizeof(Matrix));
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

  int mu,nu;
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
      for(nu = 1;nu<4;nu++){

	//First calculate the staple in the positive nu direction
        pt.run(N,tmp1,Units,dirs_m+nu);
  	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_p+nu);

	//tmp2 contains the sum of the staples for a given link
	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_plaq,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_plaq,tmp1[i],tmp2[i],vol);

	//Calculating one rectangular staple
	pt.run(N,tmp1,Units,dirs_p);
	pt.run(N,result,tmp1,dirs_m+nu);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_p+nu);

	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_rect,tmp1[i],tmp2[i],vol);

	//Calculating another rectangular staple;
	pt.run(N,tmp1,Units,dirs_m+nu);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_p+nu);
	pt.run(N,tmp1,result,dirs_p);

	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_rect,tmp1[i],tmp2[i],vol);

	//Calculating another rectangular staple;
	pt.run(N,tmp1,Units,dirs_m+nu);
	pt.run(N,result,tmp1,dirs_m+nu);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_p+nu);
	pt.run(N,tmp1,result,dirs_p+nu);

	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_rect,tmp1[i],tmp2[i],vol);

	//Calculating the staple in the negative nu direction
	pt.run(N,tmp1,Units,dirs_p+nu);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_m+nu);

	//Add this result into tmp2
	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_plaq,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_plaq,tmp1[i],tmp2[i],vol);

	//Calculating one rectangular staple
	pt.run(N,tmp1,Units,dirs_p);
	pt.run(N,result,tmp1,dirs_p+nu);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_m+nu);

	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_rect,tmp1[i],tmp2[i],vol);

	//Calculating another rectangular staple;
	pt.run(N,tmp1,Units,dirs_p+nu);
	pt.run(N,result,tmp1,dirs_m);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_m+nu);
	pt.run(N,tmp1,result,dirs_p);

	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_rect,tmp1[i],tmp2[i],vol);

	//Calculating another rectangular staple;
	pt.run(N,tmp1,Units,dirs_p+nu);
	pt.run(N,result,tmp1,dirs_p+nu);
	pt.run(N,tmp1,result,dirs_m);
	pt.run(N,result,tmp1,dirs_m+nu);
	pt.run(N,tmp1,result,dirs_m+nu);

	for(int i = 0; i<N;i++)
//	vaxpy3_m(tmp2[i],&tmp_rect,tmp1[i],tmp2[i],vol*3);
	tmp2[i]->FTimesV1PlusV2(tmp_rect,tmp1[i],tmp2[i],vol);

	//Count the flops?
        ForceFlops +=vol*288*N;
      }
      //Multiply on the left by our original link matrix to get force term
      pt.run(N,result,tmp2,dirs_p);
  }

#if 1
#pragma omp parallel for default(shared) private(mu) reduction(+:L1,L2)
  for(int index=0;index<4*vol;index++){
    Matrix mp1;
    int i = index%vol;
    mu = index/vol;
    Matrix *mtmp = (result[mu]+i);
    mp1.Dagger((IFloat *)mtmp);
    mtmp->TrLessAntiHermMatrix(mp1);
    IFloat *ihp = (IFloat *)(mom+i*4+mu);  //The gauge momentum
    IFloat *dotp = (IFloat *)mp0;
    IFloat *dotp2 = (IFloat *) (result[mu]+(i));
    fTimesV1PlusV2Single(ihp, dt, dotp2, ihp, 18);  //Update the gauge momentum
    Float norm = ((Matrix*)dotp2)->norm();
    Float tmp = sqrt(norm);
    L1 += tmp;
    L2 += norm;
  }
#else

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
	    fTimesV1PlusV2(ihp, dt, dotp2, ihp, 18);  //Update the gauge momentum
	    Float norm = ((Matrix*)dotp2)->norm();
	    Float tmp = sqrt(norm);
	    L1 += tmp;
	    L2 += norm;
	    Linf = (tmp>Linf ? tmp : Linf);
	  }
	}
      }
    }
  }
}
#endif

#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops+ParTrans::PTflops,time);
#endif

  //Free some memory
  ffree(Unit);
  for(int i = 0;i<N;i++){
  ffree(tmp1[i]);
  ffree(tmp2[i]);
  }
  for(int i = 0;i<4;i++) 
  ffree(result[i]);

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  return ForceArg(dt*L1, dt*sqrt(L2), dt*Linf);

}
CPS_END_NAMESPACE
