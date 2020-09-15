#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <util/qcdio.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <comms/scu.h>
#include <comms/glb.h>

#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/time_cps.h>
#include <alg/do_arg.h>
#include <alg/no_arg.h>
#include <alg/common_arg.h>
#include <alg/hmd_arg.h>
#include <alg/alg_plaq.h>
#include <alg/alg_rnd_gauge.h>
#include <alg/threept_arg.h>
#include <alg/threept_prop_arg.h>
#include <alg/alg_threept.h>
#include <util/smalloc.h>
#include <alg/alg_int.h>

#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

#include <util/command_line.h>
#include <sstream>

#include<unistd.h>
#include<config.h>

#include <alg/qpropw.h>
#include <alg/qpropw_arg.h>

#include <alg/alg_fix_gauge.h>
#include <alg/fix_gauge_arg.h>

#include <util/data_shift.h>

#include <alg/lanc_arg.h>
#include <alg/prop_attribute_arg.h>
#include <alg/gparity_contract_arg.h>
#include <alg/propmanager.h>
#include <alg/alg_gparitycontract.h>

#include <util/gparity_singletodouble.h>
#include <util/enum_func.h>

//some piece of **** defines these elsewhere, so the bfm header gets screwed up
#undef ND
#undef SPINOR_SIZE
#undef HALF_SPINOR_SIZE
#undef GAUGE_SIZE
#undef Nmu
#undef Ncb
#undef NMinusPlus
#undef Minus
#undef Plus
#undef DaggerYes
#undef DaggerNo
#undef SingleToDouble
#undef DoubleToSingle
#undef Odd
#undef Even

#include <alg/int_arg.h>
#include <util/lattice/bfm_evo.h>
#include <util/lattice/bfm_eigcg.h>
#include <util/lattice/fbfm.h>
#include <util/wilson.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/enum_func.h>
#include <util/sproj_tr.h>
#include <util/time_cps.h>
#include <util/lattice/fforce_wilson_type.h>
#include <alg/alg_eig.h>


#include<omp.h>

#ifdef HAVE_BFM
#include <chroma.h>
#endif

using namespace std;
USING_NAMESPACE_CPS

void setup_double_latt(Lattice &double_latt, Matrix* orig_gfield, bool gparity_X, bool gparity_Y){
  //orig latt ( U_0 U_1 ) ( U_2 U_3 ) ( U_4 U_5 ) ( U_6 U_7 )
  //double tatt ( U_0 U_1 U_2 U_3 ) ( U_4 U_5 U_6 U_7 ) ( U_0* U_1* U_2* U_3* ) ( U_4* U_5* U_6* U_7* )

  Matrix *dbl_gfield = double_latt.GaugeField();

  if(!UniqueID()){ printf("Setting up 1f lattice.\n"); fflush(stdout); }
  SingleToDoubleLattice lattdoubler(gparity_X,gparity_Y,orig_gfield,double_latt);
  lattdoubler.Run();
  if(!UniqueID()){ printf("Finished setting up 1f lattice\n"); fflush(stdout); }
}
void setup_double_rng(bool gparity_X, bool gparity_Y){
  //orig 4D rng 2 stacked 4D volumes
  //orig ([R_0 R_1][R'_0 R'_1])([R_2 R_3][R'_2 R'_3])([R_4 R_5][R'_4 R'_5])([R_6 R_7][R'_6 R'_7])
  //double (R_0 R_1 R_2 R_3)(R_4 R_5 R_6 R_7)(R'_0 R'_1 R'_2 R'_3)(R'_4 R'_5 R'_6 R'_7)
  
  //orig 5D rng 2 stacked 4D volumes per ls/2 slice (ls/2 as only one RNG per 2^4 block)

  SingleToDouble4dRNG fourDsetup(gparity_X,gparity_Y);
  SingleToDouble5dRNG fiveDsetup(gparity_X,gparity_Y);
  
  LRG.Reinitialize(); //reset the LRG and prepare for doubled lattice form
  
  if(!UniqueID()){ printf("Setting up 1f 4D RNG\n"); fflush(stdout); }
  fourDsetup.Run();      
  if(!UniqueID()){ printf("Setting up 1f 5D RNG\n"); fflush(stdout); }
  fiveDsetup.Run();    
}
void setup_double_matrixfield(Matrix* double_mat, Matrix* orig_mat, int nmat_per_site, bool gparity_X, bool gparity_Y){
  if(!UniqueID()){ printf("Setting up 1f matrix field.\n"); fflush(stdout); }
  SingleToDoubleMatrixField doubler(gparity_X,gparity_Y,nmat_per_site,orig_mat,double_mat);
  doubler.Run();
  if(!UniqueID()){ printf("Finished setting up 1f matrixfield\n"); fflush(stdout); }
}
void setup_double_4d_vector(Vector *double_vect, Vector* orig_vect, bool gparity_X, bool gparity_Y){
  if(!UniqueID()){ printf("Setting up 1f vector field.\n"); fflush(stdout); }
  SingleToDouble4dVectorField doubler(gparity_X, gparity_Y, orig_vect, double_vect, CANONICAL);
  doubler.Run();
  if(!UniqueID()){ printf("Finished setting up 1f vector field\n"); fflush(stdout); }
}
void setup_double_4d_odd_wilson_vector(Vector *double_vect, Vector* orig_vect, bool gparity_X, bool gparity_Y){
  if(!UniqueID()){ printf("Setting up 1f Wilson-ordered odd-cb vector field.\n"); fflush(stdout); }
  SingleToDouble4dVectorField doubler(gparity_X, gparity_Y, orig_vect, double_vect, WILSON, 1);
  doubler.Run();
  if(!UniqueID()){ printf("Finished setting up 1f Wilson-ordered odd-cb vector field\n"); fflush(stdout); }
}  



void GaugeTransformU(Matrix *gtrans, Lattice &lat);

Float* rand_4d_canonical_fermion(Lattice &lat){
  long f_size = (long)24 * GJP.VolNodeSites();
  if(GJP.Gparity()) f_size*=2;
  Float *v1 = (Float *)pmalloc(sizeof(Float) * f_size);
  printf("Making random gaussian 4d vector\n");
  lat.RandGaussVector((Vector*)v1, 0.5, 2, CANONICAL, FOUR_D);

  if(GJP.Gparity1fX() && GJP.Gparity1fY()){
    if(!UniqueID()){ printf("Putting minus sign on fermion source in UR quadrant\n"); fflush(stdout); }
    //make source on upper-right quadrant negative (RNGs should be correct)
    for(int t=0;t<GJP.TnodeSites();t++){
      for(int z=0;z<GJP.ZnodeSites();z++){
	for(int y=0;y<GJP.YnodeSites();y++){
	  for(int x=0;x<GJP.XnodeSites();x++){
	    int gx = x+GJP.XnodeCoor()*GJP.XnodeSites();
	    int gy = y+GJP.YnodeCoor()*GJP.YnodeSites();

	    if(gx>=GJP.Xnodes()*GJP.XnodeSites()/2 && gy>=GJP.Ynodes()*GJP.YnodeSites()/2){
	      int pos[5] = {x,y,z,t};
	      int f_off = x + GJP.XnodeSites()*(y + GJP.YnodeSites()*(z + GJP.ZnodeSites()*t));
	      f_off *= 24;

	      for(int i=0;i<24;i++) v1[f_off+i] *=-1;
	    }
	  }
	}
      }
    }
  }

  printf("Finished making random gaussian vector\n");
  return v1;
}

static void oldTm_EvolveMomFforce_nogp(Matrix *mom, Vector *chi, 
				       Float mass, Float epsilon, Float dt, Lattice *lattice)
{
//------------------------------------------------------------------
// allocate space for two CANONICAL fermion fields.
//------------------------------------------------------------------
  Matrix *gauge = lattice->GaugeField();
  const char* cname ="";
  const char* fname = "oldTm_EvolveMomFforce_nogp(..)";
  
  size_t f_size = lattice->FsiteSize() * GJP.VolNodeSites() ;
  if(GJP.Gparity()) f_size *=2;

  char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)smalloc(cname, fname, str_v1, f_size*sizeof(Float)) ;

  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)smalloc(cname, fname, str_v2, f_size*sizeof(Float)) ;

//------------------------------------------------------------------
// allocate space for two CANONICAL fermion field on a site.
//------------------------------------------------------------------

  char *str_site_v1 = "site_v1";
  Float *site_v1 = (Float *)smalloc(cname, fname, str_site_v1, lattice->FsiteSize()*sizeof(Float));

  char *str_site_v2 = "site_v2";
  Float *site_v2 = (Float *)smalloc(cname, fname, str_site_v2, lattice->FsiteSize()*sizeof(Float));

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;


  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;
    cg_arg.epsilon = epsilon;

    DiracOpWilsonTm wilson(*lattice, v1, v2, &cg_arg, CNV_FRM_YES) ;
    // chi <- (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
    wilson.CalcHmdForceVecs(chi) ;
  }
#if 0
  VRB.Result(cname,fname,"Being skipped for debugging!");
#else

  int x, y, z, t, lx, ly, lz, lt ;

  lx = GJP.XnodeSites() ;
  ly = GJP.YnodeSites() ;
  lz = GJP.ZnodeSites() ;
  lt = GJP.TnodeSites() ;

//------------------------------------------------------------------
// start by summing first over direction (mu) and then over site
// to allow SCU transfers to happen face-by-face in the outermost
// loop.
//------------------------------------------------------------------

  int mu ;

  Matrix tmp, f ;

  for (mu=0; mu<4; mu++) {
    for (t=0; t<lt; t++)
    for (z=0; z<lz; z++)
    for (y=0; y<ly; y++)
    for (x=0; x<lx; x++) {
      int gauge_offset = x+lx*(y+ly*(z+lz*t)) ;
      int vec_offset = lattice->FsiteSize()*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;

      Float *v1_plus_mu ;
      Float *v2_plus_mu ;
      int vec_plus_mu_offset = lattice->FsiteSize() ;

      Float coeff = -2.0 * dt ;

      switch (mu) {
        case 0 :
          vec_plus_mu_offset *= (x+1)%lx+lx*(y+ly*(z+lz*t)) ;
          if ((x+1) == lx) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;                        
            v2_plus_mu = site_v2 ;                        
            if (GJP.XnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 1 :
          vec_plus_mu_offset *= x+lx*((y+1)%ly+ly*(z+lz*t)) ;
          if ((y+1) == ly) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;                        
            v2_plus_mu = site_v2 ;                        
            if (GJP.YnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 2 :
          vec_plus_mu_offset *= x+lx*(y+ly*((z+1)%lz+lz*t)) ;
          if ((z+1) == lz) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            if (GJP.ZnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 3 :
          vec_plus_mu_offset *= x+lx*(y+ly*(z+lz*((t+1)%lt))) ;
          if ((t+1) == lt) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            if (GJP.TnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
      } // end switch mu

      GnoneFwilsonTm* tmlat = dynamic_cast<GnoneFwilsonTm*>(lattice);
      tmlat->sproj_tr[mu](   (IFloat *)&tmp,
                      (IFloat *)v1_plus_mu,
                      (IFloat *)v2+vec_offset, 1, 0, 0);

      tmlat->sproj_tr[mu+4]( (IFloat *)&f,
                      (IFloat *)v2_plus_mu,
                      (IFloat *)v1+vec_offset, 1, 0, 0);

      tmp += f ;

      f.DotMEqual(*(gauge+gauge_offset), tmp) ;

      tmp.Dagger(f) ;

      f.TrLessAntiHermMatrix(tmp) ;

      f *= coeff ;

      *(mom+gauge_offset) += f ;
      Float norm = f.norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }
  }
#endif

//------------------------------------------------------------------
// deallocate space for two CANONICAL fermion fields on a site.
//------------------------------------------------------------------
  VRB.Sfree(cname, fname, str_site_v2, site_v2) ;
  sfree(site_v2) ;

  VRB.Sfree(cname, fname, str_site_v1, site_v1) ;
  sfree(site_v1) ;

//------------------------------------------------------------------
// deallocate space for two CANONICAL fermion fields.
//------------------------------------------------------------------
  VRB.Sfree(cname, fname, str_v2, v2) ;
  sfree(v2) ;

  VRB.Sfree(cname, fname, str_v1, v1) ;
  sfree(v1) ;
}

//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *frm, Vector *frm,
//                 Float mass, Float epsilon, Float dt):
// It evolves the canonical momentum mom by dt
// using the boson-component of the alg_quotient force.
// Now sets epsilon in cg_arg from new input parameter
// chi <- (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
// phi <- M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
//------------------------------------------------------------------
ForceArg oldTm_EvolveMomFforce_boson_nogp(Matrix *mom, Vector *chi, Vector *eta,
					     Float mass, Float epsilon, Float dt, Lattice *lattice) {
  const char* cname ="";
  const char* fname = "oldTm_EvolveMomFforce_boson_nogp(..)";

  if (lattice->Colors() != 3)       { ERR.General(cname,fname,"Wrong nbr of colors.") ; }
  if (lattice->SpinComponents()!=4) { ERR.General(cname,fname,"Wrong nbr of spin comp.") ;}
  if (mom == 0)            { ERR.Pointer(cname,fname,"mom") ; }
  if (eta == 0)            { ERR.Pointer(cname,fname,"eta") ; }
  if (GJP.Gparity())       { ERR.General(cname,fname,"Not implemented for G-parity boundary conditions") ; }

//------------------------------------------------------------------
// allocate space for two CANONICAL fermion fields.
// these are all full fermion vector sizes ( i.e. *not* preconditioned )
//------------------------------------------------------------------
  Matrix *gauge = lattice->GaugeField();
 
  size_t f_size        ( lattice->FsiteSize() * GJP.VolNodeSites() );

  char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v1 == 0) ERR.Pointer(cname, fname, str_v1) ;
  VRB.Smalloc(cname, fname, str_v1, v1, f_size*sizeof(Float)) ;

  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v2 == 0) ERR.Pointer(cname, fname, str_v2) ;
  VRB.Smalloc(cname, fname, str_v2, v2, f_size*sizeof(Float)) ;

//------------------------------------------------------------------
// two fermion vectors at a single position
//    - these will be used to store off-node
//      field components
//------------------------------------------------------------------

  char *str_site_v1 = "site_v1";
  Float *site_v1 = (Float *)smalloc(lattice->FsiteSize()*sizeof(Float));
  if (site_v1 == 0) ERR.Pointer(cname, fname, str_site_v1) ;
  VRB.Smalloc(cname, fname, str_site_v1, site_v1,
    lattice->FsiteSize()*sizeof(Float)) ;

  char *str_site_v2 = "site_v2";
  Float *site_v2 = (Float *)smalloc(lattice->FsiteSize()*sizeof(Float));
  if (site_v2 == 0) ERR.Pointer(cname, fname, str_site_v2) ;
  VRB.Smalloc(cname, fname, str_site_v2, site_v2,
    lattice->FsiteSize()*sizeof(Float)) ;
  
  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;


  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;
    cg_arg.epsilon = epsilon;

    DiracOpWilsonTm wilson(*lattice, v1, v2, &cg_arg, CNV_FRM_YES) ;

//~~
//~~ fermion version:  	wilson.CalcHmdForceVecs(chi)
//~~ boson version:  	wilson.CalcBsnForceVecs(chi, eta)
//~~
    // chi <- (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
    // phi <- M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
    wilson.CalcBsnForceVecs(chi, eta) ;
  }

#if 0
  VRB.Result(cname,fname,"Being skipped for debugging!");
#else

  // evolve the momenta by the fermion force
  int mu, x, y, z, t;
 
  const int lx(GJP.XnodeSites());
  const int ly(GJP.YnodeSites());
  const int lz(GJP.ZnodeSites());
  const int lt(GJP.TnodeSites());

//------------------------------------------------------------------
// start by summing first over direction (mu) and then over site
// to allow SCU transfers to happen face-by-face in the outermost
// loop.
//------------------------------------------------------------------

// (1-gamma_\mu) Tr_s[v1(x+\mu) v2^{\dagger}(x)] +          
// 		(1+gamma_\mu) Tr_s [v2(x+\mu) v1^{\dagger}(x)]
  for (mu=0; mu<4; mu++) {
    for (t=0; t<lt; t++)
      for (z=0; z<lz; z++)
        for (y=0; y<ly; y++)
          for (x=0; x<lx; x++) {
            // position offset
            int gauge_offset = x+lx*(y+ly*(z+lz*t)) ;
            
            // offset for vector field at this point
            int vec_offset = lattice->FsiteSize()*gauge_offset ;

            // offset for link in mu direction from this point
            gauge_offset = mu+4*gauge_offset ;

            Float *v1_plus_mu ;
            Float *v2_plus_mu ;
            int vec_plus_mu_offset = lattice->FsiteSize() ;

            // sign of coeff (look at momenta update)
            Float coeff = -2.0 * dt ;

            switch (mu) {
              case 0 :
                // next position in mu direction
                vec_plus_mu_offset *= (x+1)%lx+lx*(y+ly*(z+lz*t)) ;
                if ((x+1) == lx) {
                   // off-node
                   // fill site_v1 and site_v2 with v1 and v2 data
                   // from x=0 on next node, need loop because
                   // data is not contiguous in memory 
                   getPlusData( (IFloat *)site_v1,
                              (IFloat *)v1+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
                   getPlusData( (IFloat *)site_v2,
                              (IFloat *)v2+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;

                   v1_plus_mu = site_v1 ;                        
                   v2_plus_mu = site_v2 ;                        

                   // GJP.XnodeBc() gives the forward boundary
                   // condition only (so this should work).
                   if (GJP.XnodeBc()==BND_CND_APRD) coeff = -coeff ;
                } else {
                   //
                   //  on - node: just add offset to v1 and v2
                   // (they are now 1 forward in the mu direction )
                   //
                   v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
                   v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
                }
                break ;

        // Repeat for the other directions
        case 1 :
          vec_plus_mu_offset *= x+lx*((y+1)%ly+ly*(z+lz*t)) ;
          if ((y+1) == ly) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;                        
            v2_plus_mu = site_v2 ;                        
            if (GJP.YnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 2 :
          vec_plus_mu_offset *= x+lx*(y+ly*((z+1)%lz+lz*t)) ;
          if ((z+1) == lz) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            if (GJP.ZnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 3 :
          vec_plus_mu_offset *= x+lx*(y+ly*(z+lz*((t+1)%lt))) ;
          if ((t+1) == lt) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, lattice->FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            if (GJP.TnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
      } // end (the evil) mu switch 

      Matrix tmp_mat1, tmp_mat2;  

      GnoneFwilsonTm *tmlat = dynamic_cast<GnoneFwilsonTm *>(lattice);

// ( 1 - gamma_\mu ) Tr_s [ v1(x+\mu) v2^{\dagger}(x) ]           
      tmlat->sproj_tr[mu](   (IFloat *)&tmp_mat1,   	// output color matrix
                      (IFloat *)v1_plus_mu,		// row vector, NOT conjugated
                      (IFloat *)v2+vec_offset, 	// col vector, IS conjugated
                      1, 0, 0);				// 1 block, 0 strides

// (1 + gamma_\mu)  Tr_s [ v2(x+\mu) v1^{\dagger}(x) ]
      tmlat->sproj_tr[mu+4]( (IFloat *)&tmp_mat2,		// output color matrix
                      (IFloat *)v2_plus_mu,		// row vector, NOT conjugated
                      (IFloat *)v1+vec_offset, 	// col vector, IS conjugated
                      1, 0, 0);				// 1 block, 0 strides

      // exactly what this sounds like
      tmp_mat1 += tmp_mat2 ;
            
      // multiply sum by the link in the \mu direction
      tmp_mat2.DotMEqual(*(gauge+gauge_offset), tmp_mat1) ;
      
      // take tracless antihermitian piece
      // TrLessAntiHermMatrix need to be passed
      // the dagger of the matrix in question
      tmp_mat1.Dagger(tmp_mat2) ;
      tmp_mat2.TrLessAntiHermMatrix(tmp_mat1) ;

      tmp_mat2 *= coeff ;
            
//~~
//~~ fermion version:  	(mom+gauge_offset) += f
//~~ boson version:  	(mom+gauge_offset) -= f
//~~
      *(mom+gauge_offset) -= tmp_mat2 ;

	 Float norm = tmp_mat2.norm();
	 Float tmp = sqrt(norm);
	 L1 += tmp;
	 L2 += norm;
	 Linf = (tmp>Linf ? tmp : Linf);
	 
    }
  }

//------------------------------------------------------------------
// deallocate space for two CANONICAL fermion fields on a site.
//------------------------------------------------------------------
  VRB.Sfree(cname, fname, str_site_v2, site_v2) ;
  sfree(site_v2) ;

  VRB.Sfree(cname, fname, str_site_v1, site_v1) ;
  sfree(site_v1) ;

//------------------------------------------------------------------
// deallocate space for two CANONICAL fermion fields.
//------------------------------------------------------------------
  VRB.Sfree(cname, fname, str_v2, v2) ;
  sfree(v2) ;

  VRB.Sfree(cname, fname, str_v1, v1) ;
  sfree(v1) ;

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();
#endif

  VRB.FuncEnd(cname,fname);
  return ForceArg(L1, sqrt(L2), Linf);
}

static Float* convertCanonicalSingleChkb(Float* cbferm, const int &parity, Lattice &lat){
  long f_size = GJP.VolNodeSites()*24;
  long cb_size = f_size/2;

  Float* canonical = (Float*)pmalloc(f_size * sizeof(Float));
  for(int i=0;i<f_size;i++) canonical[i] = 0.0;
  
  long start = 0;
  if(parity==0) start = cb_size;
  
  int j=0;
  for(int i=start; i<start+cb_size;i++, j++) canonical[i] = cbferm[j];

  lat.Fconvert((Vector*)canonical, CANONICAL, WILSON);

  return canonical;
}


static void generateEigArg(EigArg &eig_arg, const Float &mass, const Float &epsilon) {
  const char *fname = "generateEigArg()";
  const char *cname = "";
  
  EigenDescr eigen;
  eigen.eigen_measure = EIGEN_MEASURE_YES;
  eigen.stop_rsd = 0.00002;
  eigen.max_num_iter = 10000;
  eigen.eig_lo_stem = "eig_low";
  eigen.eig_hi_stem = "eig_hi";

  eig_arg.Cv_fact = 0.0;
  eig_arg.N_eig = 2;
  eig_arg.mass = mass;
  eig_arg.epsilon = epsilon;
  eig_arg.fname = "eigenpooh";

  eig_arg.RitzMatOper = MATPCDAG_MATPC;
  
  int n_masses = 1;

  //!< Setup AlgEig parameters if necessary
  eig_arg.pattern_kind = ARRAY;
  eig_arg.Mass.Mass_len = n_masses;
  eig_arg.Mass.Mass_val = 
    (Float*) smalloc(n_masses*sizeof(Float),"Mass_val", fname, cname);

  //CK: added for twisted mass fermions
  eig_arg.Epsilon.Epsilon_len = n_masses;
  eig_arg.Epsilon.Epsilon_val = 
    (Float*) smalloc(n_masses*sizeof(Float),"Epsilon_val", fname, cname);

  eig_arg.Mass.Mass_val[0] = mass;
  eig_arg.Epsilon.Epsilon_val[0] = epsilon;

  eig_arg.Kalk_Sim = 0;
  eig_arg.MaxCG = eigen.max_num_iter;
  eig_arg.RsdR_a = eigen.stop_rsd;
  eig_arg.RsdR_r = eigen.stop_rsd;
  eig_arg.Rsdlam = eigen.stop_rsd;
  eig_arg.Cv_fact =   0.0;
  eig_arg.N_min = 0;
  eig_arg.N_max = 5000;
  eig_arg.N_KS_max = 0;
  eig_arg.n_renorm = 100;
  eig_arg.ProjApsiP = 0;
  eig_arg.print_hsum = 0;
  eig_arg.hsum_dir = 0;
  eig_arg.ncorr = 0;
}

// class ActionQuotientArg {

//   memfun void resize(u_int nmass);

//   // Mass parameter here is dummy
//   ActionBilinearArg bi_arg;

//   //!< Quotient parameters
//   QuotientDescr quotients<>;

// };


// class QuotientDescr {
//   Float bsn_mass;
//   //! ~~epsilon parameter for twisted mass wilson fermions
//   Float bsn_mass_epsilon;
//   Float frm_mass;
//   //! ~~epsilon parameter for twisted mass wilson fermions
//   Float frm_mass_epsilon;
//   int   chrono;
//   Float stop_rsd_hb;
//   //! The multiplier we used for additional CG inversion in force
//   //! gradient integrator.
//   Float stop_rsd_fg_mult;
//   Float stop_rsd_md;
//   Float stop_rsd_mc;
// } ;



void setupQuoArg(ActionQuotientArg &into, const int &ndet, Float* bsn_masses, Float* bsn_epsilon, Float* frm_masses, Float* frm_epsilon, FclassType fermion_type = F_CLASS_WILSON_TM){
  //bi_arg
  into.bi_arg.fermion = fermion_type;
  into.bi_arg.bilinears.bilinears_len = ndet;
  into.bi_arg.bilinears.bilinears_val = new BilinearDescr[ndet];
  for(int i=0;i<ndet;i++){
    into.bi_arg.bilinears.bilinears_val[i].mass = 0.0;
    into.bi_arg.bilinears.bilinears_val[i].max_num_iter = 5000;
  }
  into.bi_arg.action_arg.force_measure = FORCE_MEASURE_YES;
  into.bi_arg.action_arg.force_label = "Quotient";
  
  //quotients
  into.quotients.quotients_len = ndet;
  into.quotients.quotients_val = new QuotientDescr[ndet];
  for(int i=0;i<ndet;i++){
    QuotientDescr &quo = into.quotients.quotients_val[i];
    quo.bsn_mass = bsn_masses[i];
    quo.bsn_mass_epsilon = bsn_epsilon[i];

    quo.frm_mass = frm_masses[i];
    quo.frm_mass_epsilon = frm_epsilon[i];

    quo.chrono =0;
    quo.stop_rsd_hb =   1.0000000000000000e-10;
    quo.stop_rsd_fg_mult =   1.0000000000000000e+00;
    quo.stop_rsd_md =   1.0000000000000000e-08;
    quo.stop_rsd_mc =   1.0000000000000000e-10;
  }
}


//CK: This is set up to correctly normalize the eigenvalue bounds between Fbfm and CPS to take into account the differing normalization
//of the preconditioned matrix. For given input bounds, the CPS and BFM forces will agree when this normalization is applied.
#define DO_EVAL_NORM_BFM

void setupRatQuoArg(ActionRationalQuotientArg &into, const int &ndet, Float* bsn_masses, Float* bsn_epsilon, Float* frm_masses, Float* frm_epsilon, int *pwr_num, int *pwr_den, FclassType fermion_type = F_CLASS_WILSON_TM){
  //bi_arg
  into.bi_arg.fermion = fermion_type;
  into.bi_arg.bilinears.bilinears_len = ndet;
  into.bi_arg.bilinears.bilinears_val = new BilinearDescr[ndet];
  for(int i=0;i<ndet;i++){
    into.bi_arg.bilinears.bilinears_val[i].mass = 0.0;
    into.bi_arg.bilinears.bilinears_val[i].max_num_iter = 5000;
  }
  into.bi_arg.action_arg.force_measure = FORCE_MEASURE_YES;
  into.bi_arg.action_arg.force_label = "RationalQuotient";
  
  into.spread = 0.0;
  into.remez_generate= 0;
  if(fermion_type == F_CLASS_WILSON_TM) into.rat_poles_file = "rqpoles.wilsonTm.vml";
  else into.rat_poles_file = "rqpoles.fbfm.vml";

  //bsn_mass
  into.bsn_mass.bsn_mass_len = ndet;
  into.bsn_mass.bsn_mass_val = new Float[ndet];
  for(int i=0;i<ndet;i++){
    into.bsn_mass.bsn_mass_val[i] = bsn_masses[i];
  }
  into.bsn_mass_epsilon.bsn_mass_epsilon_len = ndet;
  into.bsn_mass_epsilon.bsn_mass_epsilon_val = new Float[ndet];
  for(int i=0;i<ndet;i++){
    into.bsn_mass_epsilon.bsn_mass_epsilon_val[i] = bsn_epsilon[i];
  }

  //frm_mass
  into.frm_mass.frm_mass_len = ndet;
  into.frm_mass.frm_mass_val = new Float[ndet];
  for(int i=0;i<ndet;i++){
    into.frm_mass.frm_mass_val[i] = frm_masses[i];
  }
  into.frm_mass_epsilon.frm_mass_epsilon_len = ndet;
  into.frm_mass_epsilon.frm_mass_epsilon_val = new Float[ndet];
  for(int i=0;i<ndet;i++){
    into.frm_mass_epsilon.frm_mass_epsilon_val[i] = frm_epsilon[i];
  }

  //bosons
  into.bosons.bosons_len = ndet;
  into.bosons.bosons_val = new RationalDescr[ndet];
  for(int i=0;i<ndet;i++){
    Float kappa = 1.0/2.0/sqrt( (bsn_masses[i] + 4.0)*(bsn_masses[i] + 4.0) + bsn_epsilon[i] * bsn_epsilon[i] );
    Float fbfm_factor = 4*kappa*kappa*(4+bsn_masses[i]); //CPS eigenvalues are fbfm_factor * BFM eigenvalues due to normalization differences

    into.bosons.bosons_val[i].field_type = BOSON;
    into.bosons.bosons_val[i].power_num = pwr_num[i];
    into.bosons.bosons_val[i].power_den = pwr_den[i];
    into.bosons.bosons_val[i].precision = 40;
    into.bosons.bosons_val[i].stop_rsd_fg_mult = 1.0;
    
    ApproxDescr *approx = &into.bosons.bosons_val[i].md_approx;
    approx->approx_type = RATIONAL_APPROX_POWER;
    approx->bounds_type = RATIONAL_BOUNDS_MANUAL;
    approx->lambda_low =   1.0100000000000001e-04;
    approx->lambda_high =   1.000000000000000e+01;
    approx->stop_rsd.stop_rsd_len = 15;
    approx->stop_rsd.stop_rsd_val = new Float[15];
    for(int j=0;j<15;j++) approx->stop_rsd.stop_rsd_val[j] = 1.0e-7;

#ifdef DO_EVAL_NORM_BFM
    if(fermion_type != F_CLASS_WILSON_TM){
      approx->lambda_low /= fbfm_factor;
      approx->lambda_high /= fbfm_factor;
    }
#endif    

    approx = &into.bosons.bosons_val[i].mc_approx;
    approx->approx_type = RATIONAL_APPROX_POWER;
    approx->bounds_type = RATIONAL_BOUNDS_MANUAL;
    approx->lambda_low =   1.0200000000000002e-04;
    approx->lambda_high =   1.000000000000000e+01;
    approx->stop_rsd.stop_rsd_len = 15;
    approx->stop_rsd.stop_rsd_val = new Float[15];
    for(int j=0;j<15;j++) approx->stop_rsd.stop_rsd_val[j]= 1.0e-7; 

#ifdef DO_EVAL_NORM_BFM
    if(fermion_type != F_CLASS_WILSON_TM){
      approx->lambda_low /= fbfm_factor;
      approx->lambda_high /= fbfm_factor;
    }
#endif

    into.bosons.bosons_val[i].stag_bsn_mass = 0.0;
  }
  //fermions
  into.fermions.fermions_len = ndet;
  into.fermions.fermions_val = new RationalDescr[ndet];
  for(int i=0;i<ndet;i++){
    Float kappa = 1.0/2.0/sqrt( (frm_masses[i] + 4.0)*(frm_masses[i] + 4.0) + frm_epsilon[i] * frm_epsilon[i] );
    Float fbfm_factor = 4*kappa*kappa*(4+frm_masses[i]); //CPS eigenvalues are fbfm_factor * BFM eigenvalues due to normalization differences

    into.fermions.fermions_val[i].field_type = FERMION;
    into.fermions.fermions_val[i].power_num = pwr_num[i];
    into.fermions.fermions_val[i].power_den = pwr_den[i];
    into.fermions.fermions_val[i].precision = 40;
    into.fermions.fermions_val[i].stop_rsd_fg_mult = 1.0;
    
    ApproxDescr *approx = &into.fermions.fermions_val[i].md_approx;
    approx->approx_type = RATIONAL_APPROX_POWER;
    approx->bounds_type = RATIONAL_BOUNDS_MANUAL;
    approx->lambda_low =   1.0300000000000003e-04;
    approx->lambda_high =   1.000000000000000e+01;
    approx->stop_rsd.stop_rsd_len = 15;
    approx->stop_rsd.stop_rsd_val = new Float[15];
    for(int j=0;j<15;j++) approx->stop_rsd.stop_rsd_val[j] =   1.0e-7;

#ifdef DO_EVAL_NORM_BFM
    if(fermion_type != F_CLASS_WILSON_TM){
      approx->lambda_low /= fbfm_factor;
      approx->lambda_high /= fbfm_factor;
    }
#endif

    approx = &into.fermions.fermions_val[i].mc_approx;
    approx->approx_type = RATIONAL_APPROX_POWER;
    approx->bounds_type = RATIONAL_BOUNDS_MANUAL;
    approx->lambda_low =   1.0400000000000004e-04;
    approx->lambda_high =   1.0000000000000000e+01;
    approx->stop_rsd.stop_rsd_len = 15;
    approx->stop_rsd.stop_rsd_val = new Float[15];
    for(int j=0;j<15;j++) approx->stop_rsd.stop_rsd_val[j]= 1.0000000000000000e-7; 

#ifdef DO_EVAL_NORM_BFM
    if(fermion_type != F_CLASS_WILSON_TM){
      approx->lambda_low /= fbfm_factor;
      approx->lambda_high /= fbfm_factor;
    }
#endif

    into.fermions.fermions_val[i].stag_bsn_mass = 0.0;
  }

  into.eigen.eigen_measure = EIGEN_MEASURE_YES;
  into.eigen.stop_rsd = 0.00002;
  into.eigen.max_num_iter = 10000;
  into.eigen.eig_lo_stem = "eig_low";
  into.eigen.eig_hi_stem = "eig_hi";
}

//extern void g5theta(Vector *in, int vol, IFloat ctheta, IFloat stheta);



extern "C" {
  inline void cTimesC(IFloat *a, IFloat re, IFloat im);
  void g5theta(Vector *in, int vol, IFloat ctheta, IFloat stheta);
}
inline void cTimesC_2(IFloat *a, IFloat re, IFloat im)
{
        // a points to real part
	IFloat t;               
	t = (*a);                         // save real part
	*a = re * (*a) - im * *(a+1);     // real part
	*(a+1)   = re * *(a+1) + im * t;  // imag part
};


void CPS_CalcHmdForceVecs(Vector *chi, Vector *f_out, Vector *f_in, DiracOpWilsonTm &d_op)
{
  const char *fname = "CalcHmdForceVecs(V*) [CPS version]" ;
  const char *cname = "";
  VRB.Func(cname,fname) ;

  if (f_out == 0)
    ERR.Pointer(cname, fname, "f_out") ;

  if (f_in == 0)
    ERR.Pointer(cname, fname, "f_in") ;

//------------------------------------------------------------------
// f_out stores (chi,rho) (the v of eqn B12): 
// rho = gamma_5(-theta) Dslash chi
// f_in stores (psi,sigma) (the w of eqn B11):
// psi = gamma_5(-theta) MatPc chi
// sigma = gamma_5(-theta) Dslash psi
//------------------------------------------------------------------

  Float ctheta = d_op.get_ctheta();
  Float stheta = d_op.get_stheta();
  Float kappa = d_op.get_kappa();

  Vector *chi_new, *rho, *psi, *sigma ;

  int vol =  GJP.VolNodeSites()/2;
  size_t f_size_cb = 12 * GJP.VolNodeSites() ;
  if(GJP.Gparity()){ 
    vol*=2;
    f_size_cb *= 2; //Layout is   |   odd   |   even  |
                    //            | f0 | f1 | f0 | f1 |
                    //where for each checkerboard, each flavour field occupies one half-volume
  }

  chi_new = f_out ;
  rho = (Vector *)((Float *)f_out + f_size_cb) ;
  psi = f_in ;
  sigma = (Vector *)((Float *)f_in + f_size_cb) ;

  chi_new->CopyVec(chi, f_size_cb) ;

  d_op.MatPc(psi,chi) ;
  g5theta(psi, vol, ctheta, stheta);

  psi->VecTimesEquFloat(-kappa*kappa,f_size_cb) ;

  d_op.Dslash(rho, chi, CHKB_ODD, DAG_NO) ;
  g5theta(rho, vol, ctheta, -stheta);

  d_op.Dslash(sigma, psi, CHKB_ODD, DAG_YES) ;
  g5theta(sigma, vol, ctheta, stheta);

  return ;
}

int InvCgShift_CPS(Vector *out, 
		   Vector *in, 
		   Float src_norm_sq, 
		   Float *true_res,
		   Float *shift,
		   DiracOp &dop,
		   Lattice& lat, CgArg &cg_arg){
  int itr;                       // Current number of CG iterations
  int max_itr;                       // Max number of CG iterations
  Float stp_cnd;                   // Stop if residual^2 <= stp_cnd
  Float res_norm_sq_prv;          // The previous step |residual|^2
  Float res_norm_sq_cur;           // The current step |residual|^2 
  Float a;
  Float b;
  Float d;
  int i, ic, icb;
  const char *fname = "InvCgShift(V*,V*,F,F*) [Duplicate of d_op_base/noarch version]";
  const char* cname = "";

  IFloat *temp;


// Print out input parameters
//------------------------------------------------------------------
  VRB.Input(cname,fname,
	    "stop_rsd = %e\n",IFloat(cg_arg.stop_rsd));
  VRB.Input(cname,fname,
	    "max_num_iter = %d\n",cg_arg.max_num_iter);
  VRB.Input(cname,fname,
	    "mass = %e\n",IFloat(cg_arg.mass));
  VRB.Input(cname,fname,
	    "src_norm_sq = %e\n",IFloat(src_norm_sq));


//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------

// Set the source vector pointer
//------------------------------------------------------------------
  Vector *src = in;

// Set the solution vector pointer
//------------------------------------------------------------------
  Vector *sol = out;

// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------

  size_t f_size_cb;

  if(lat.Fclass() == F_CLASS_CLOVER) {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  } else {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
  }
    
  if(GJP.Gparity()) f_size_cb*=2;

  Vector *res = (Vector *) smalloc(f_size_cb * sizeof(Float));
  Vector *dir = (Vector *) smalloc(f_size_cb * sizeof(Float));
  Vector *mmp = (Vector *) smalloc(f_size_cb * sizeof(Float));

// If src_norm_sq is not provided calculate it
//------------------------------------------------------------------
  if(src_norm_sq == 0){
    src_norm_sq = src->NormSqNode(f_size_cb); //CK: in G-parity situation we want the norm^2 of the whole 2-flavour double-wrapped source
    glb_sum(&src_norm_sq);
  }
  VRB.Flow(cname,fname,"src_norm_sq=%e\n",src_norm_sq);

// Calculate stopping condition
//------------------------------------------------------------------
  stp_cnd = src_norm_sq * cg_arg.stop_rsd * cg_arg.stop_rsd;
  VRB.Flow(cname,fname, 
	   "stp_cnd =%e\n", IFloat(stp_cnd));

// Make IFloat pointers out of Vector pointers
//------------------------------------------------------------------
  IFloat *f_sol = (IFloat *) sol; 
  IFloat *f_dir = (IFloat *) dir; 
  IFloat *f_res = (IFloat *) res; 
  IFloat *f_mmp = (IFloat *) mmp; 

//------------------------------------------------------------------
// Initial step:
// res = src - MatPcDagMatPc * sol
// dir = res
// if( |res|^2 <= stp_cnd ){ 
//   n_count = 0
//   free memory
//   return
// }
//------------------------------------------------------------------
  Float *in_f =  (Float *) sol;
  // Mmp = MatPcDagMatPc * sol
  dop.MatPcDagMatPc(mmp, sol);
  if (shift){
    mmp -> FTimesV1PlusV2(*shift,sol,mmp, f_size_cb);
  }
  //print_vec( mmp, "mmp");

  // res = src
  res->CopyVec(src, f_size_cb);
  //print_vec( res, "res");

  // res -= mmp
  res->VecMinusEquVec(mmp, f_size_cb);
  //print_vec( res, "res");

  // dir = res
  dir->CopyVec(res, f_size_cb);  
  //print_vec( dir, "dir");

  // res_norm_sq_cur = res * res
  res_norm_sq_cur = res->NormSqNode(f_size_cb);
  //printf("res_norm_sq_cur=%e\n",res_norm_sq_cur);
  glb_sum(&res_norm_sq_cur);

  // if( |res|^2 <= stp_cnd ) we are done
  VRB.Flow(cname,fname,
  	   "|res[0]|^2 = %e\n", IFloat(res_norm_sq_cur));
  itr = 0;
  max_itr = 9999;
  if(res_norm_sq_cur <= stp_cnd) max_itr = 0;
  //printf("max_itr=%d\n",max_itr);


//------------------------------------------------------------------
// Loop over CG iterations
//------------------------------------------------------------------

  for(i=0; i < max_itr; i++){
    itr = itr + 1;
    res_norm_sq_prv = res_norm_sq_cur;

    // mmp = MatPcDagMatPc * dir
    // d = <dir, MatPcDagMatPc*dir>

    dop.MatPcDagMatPc(mmp, dir, &d);
    if (shift){
      mmp -> FTimesV1PlusV2(*shift,dir,mmp, f_size_cb);
      Float dir_sq = dir -> NormSqNode(f_size_cb);
      glb_sum(&dir_sq);
      d += (*shift) * dir_sq;
    }
    //printf("d=%e\n",d);
    //print_vec( mmp, "mmp");
  
    glb_sum(&d);
    VRB.Flow(cname,fname, "d = %e\n", IFloat(d));

    // If d = 0 we are done
    if(d == 0.0) {
      VRB.Warn(cname,fname,"d(%e) = 0.0!!\n",d);
      //	exit(5);
      break;
      //??? or should we give a warning or error? Yes we should, really.
    }

    a = res_norm_sq_prv / d;
    VRB.Flow(cname,fname, "a = %e\n", IFloat(a));

    // Set circular buffer
    //    setCbufCntrlReg(4, CBUF_MODE4);

    // sol = a * dir + sol;
    sol->FTimesV1PlusV2(a, dir, sol, f_size_cb);
    //print_vec( sol, "sol");

    // res = - a * (MatPcDagMatPc * dir) + res;
    res->FTimesV1PlusV2(-a, mmp, res, f_size_cb);
    //print_vec( res, "res");

    // res_norm_sq_cur = res * res
    res_norm_sq_cur = res->NormSqNode(f_size_cb);
    glb_sum(&res_norm_sq_cur);

    // if( |res|^2 <= stp_cnd ) we are done
    VRB.Flow(cname,fname,
	     "|res[%d]|^2 = %e\n", itr, IFloat(res_norm_sq_cur));
    if(res_norm_sq_cur <= stp_cnd) break;

    b = res_norm_sq_cur / res_norm_sq_prv;
    VRB.Flow(cname,fname, "b = %e\n", IFloat(b));

    // dir = b * dir + res;
    dir->FTimesV1PlusV2(b, dir, res, f_size_cb);
    //print_vec( dir, "dir");
  }

  // It has not reached stp_cnd: Issue a warning
  if(itr == cg_arg.max_num_iter - 1){
    VRB.Warn(cname,fname,
	      "CG reached max iterations = %d. |res|^2 = %e\n",
	     itr+1, IFloat(res_norm_sq_cur) );
  }

//------------------------------------------------------------------
// Done. Finish up and return
//------------------------------------------------------------------
  // Calculate and set true residual: 
  // true_res = |src - MatPcDagMatPc * sol| / |src|
  dop.MatPcDagMatPc(mmp, sol);
  if (shift){
    mmp -> FTimesV1PlusV2(*shift,sol,mmp, f_size_cb);
  }
  res->CopyVec(src, f_size_cb);
  res->VecMinusEquVec(mmp, f_size_cb);
  res_norm_sq_cur = res->NormSqNode(f_size_cb);
  glb_sum(&res_norm_sq_cur);
  Float tmp = res_norm_sq_cur / src_norm_sq;
  tmp = sqrt(tmp);
  if(true_res != 0){
    *true_res = tmp;
  }
  VRB.Result(cname,fname,
	     "True |res| / |src| = %e, iter = %d\n", IFloat(tmp), itr+1);

  // Free memory
  VRB.Sfree(cname,fname, "mmp", mmp);
  sfree(mmp);
  VRB.Sfree(cname,fname, "dir", dir);
  sfree(dir);
  VRB.Debug("b ============\n");
  VRB.Sfree(cname,fname, "res", res);
  sfree(res);

  VRB.Debug("a ============\n");

  // Return number of iterations
  return itr+1;
}

int InvCg_CPS(Vector *out, 
	      Vector *in, 
	      Float src_norm_sq, 
	      Float *true_res,
	      DiracOp &dop,
	      Lattice& lat, CgArg &cg_arg){
  return InvCgShift_CPS(out,in,src_norm_sq,true_res,NULL,dop,lat,cg_arg);
}



int MInvCG_CPS(Vector **psi, Vector *chi, Float chi_norm, Float *mass, 
	       int Nmass, int isz, Float *RsdCG,
	       MultiShiftSolveType type, Float *alpha, DiracOp &dop, Lattice &lat, CgArg &cg_arg)
{
  const char *fname = "MInvCG(V*,V**,...) [Duplicate of d_op_base/noarch version]";
  const char *cname = "";
  VRB.Func(cname,fname);
    
  if( (lat.Fclass() != F_CLASS_DWF) && (GJP.Snodes()>1) )
    ERR.General(cname,fname,"Fermion class type inconsistent with spread-out S dimension\n");


// Print out input parameters
//------------------------------------------------------------------
  VRB.Result(cname,fname,
	    "number of shifts = %d\n",Nmass);
  VRB.Result(cname,fname,
	    "smallest shift stop_rsd = %e\n",IFloat(RsdCG[0]));
  VRB.Result(cname,fname,
	    "max_num_iter = %d\n",cg_arg.max_num_iter);
  VRB.Result(cname,fname,
	    "mass = %e\n",IFloat(cg_arg.mass));
  VRB.Result(cname,fname,
	    "src_norm_sq = %e\n",IFloat(chi_norm));

//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------

  int iz, k, s;
  size_t f_size;

  if(lat.Fclass() == F_CLASS_CLOVER)
    f_size = lat.FsiteSize()*GJP.VolNodeSites() / 2;
  else
    f_size = lat.FsiteSize()*GJP.VolNodeSites() / (lat.FchkbEvl()+1);

  if(GJP.Gparity()) f_size*=2;

  Vector *r = (Vector*)smalloc(f_size*sizeof(Float),cname,fname,"r");
  
  Vector *Ap = (Vector*)smalloc(f_size*sizeof(Float),cname,fname,"Ap");
  
  Vector **p = (Vector**)smalloc(Nmass*sizeof(Vector*),cname,fname,"p");

  for (s=0; s<Nmass; s++) 
    *(p+s) = (Vector*)smalloc(f_size * sizeof(Float),
			      cname,fname, "p[i]");
    
  int convP;
  int *convsP = (int*)smalloc(Nmass*sizeof(int));
  int *converged = (int*)smalloc(Nmass*sizeof(int));
  
  Float a=0, as, b, bp, b_tmp;
  Float *bs = (Float*)smalloc(Nmass * sizeof(Float));
  Float **z = (Float**)smalloc(2 * sizeof(Float*));
  for (s=0; s<2; s++) *(z+s) = (Float*)smalloc(Nmass * sizeof(Float));
  Float css, ztmp;
  
  Float c, cs, d, cp;
  Float *dot=0;

  Float *rsd_sq = (Float*)smalloc(Nmass*sizeof(Float)); 
  Float *rsdcg_sq = (Float*)smalloc(Nmass*sizeof(Float)); 

  if (type == MULTI) {
    for (int i=0; i<Nmass; i++)
      psi[i] -> VecZero(f_size);
  }
  
  if (type == SINGLE || type == MULTI) {
    r-> CopyVec(chi,f_size);
    cp = chi_norm;
  } else if (type == GENERAL) {
    dop.MatPcDagMatPc(r,psi[0]);
    r -> FTimesV1MinusV2(1.0,chi,r,f_size);
    cp = r -> NormSqGlbSum(f_size);
  }

  for (s=0; s<Nmass; s++) p[s] -> CopyVec(r,f_size);

  for (s=0; s<Nmass; s++) {
    rsdcg_sq[s] = RsdCG[s]*RsdCG[s];
    rsd_sq[s] = cp*rsdcg_sq[s];
    converged[s] = 0;
  }

  /*  d = <p, A.p>  */
  if (mass[0] > 0) {
    dop.MatPcDagMatPc(Ap,p[0]); //Ap = Mpc^dag Mpc p
    Ap -> FTimesV1PlusV2(mass[0],p[0],Ap, f_size); //Ap = mass[0]*p + Ap
    d = p[0] -> ReDotProductGlbSum(Ap, f_size); //d = p.Ap =  p. (Mpc^dag Mpc + mass[0])p
  } else { //CK: Why do we do something different for negative mass??
    dop.MatPcDagMatPc(Ap,p[0],&d);
    glb_sum(&d);
  }
  IFloat *Ap_tmp = (IFloat *)Ap;
  VRB.Flow(cname,fname,"Ap= %e pAp =%e\n",*Ap_tmp,d);

  b = -cp/d;

  VRB.Flow(cname,fname,"b = -cp/d = -%e/%e =%e\n",cp,d,b);

  z[0][0] = 1.0;
  z[1][0] = 1.0;
  bs[0] = b;
  iz = 1;
  
  for (s=0; s<Nmass; s++) {
    if (s==0) continue;
    z[1-iz][s] = 1.0;
    z[iz][s] = 1.0 / ( 1.0 - b*(mass[s] - mass[0]) );
    bs[s] = b*z[iz][s];
  }

  // r[1] += b[0] A.p[0]
  r -> FTimesV1PlusV2(b,Ap,r,f_size);
  // c = |r[1]|^2
  c = r -> NormSqGlbSum(f_size);
  VRB.Flow(cname,fname,"|r[1]|^2 =%e\n", c);

  // Psi[1] -= b[0] p[0] =- b[0] chi;
  if (type == SINGLE) {
    for (s=0; s<Nmass; s++) {
      b_tmp = bs[s] * alpha[s];
      psi[0] -> FTimesV1PlusV2(-b_tmp,chi,psi[0],f_size);
    }
  } else {
    for (s=0; s<Nmass; s++) 
      psi[s]-> FTimesV1PlusV2(-bs[s],chi,psi[s],f_size);  
  }

  // Check the convergance of the first solution
  for (s=0; s<Nmass; s++) convsP[s] = 0;
  
  convP = (c < rsd_sq[0]) ? 1 : 0;
  
  // a[k+1] = |r[k]**2/ |r[k-1]|**2
  a = c/cp;
  
  // for k=1 until MaxCG do
  // if |psi[k+1] - psi[k]| <= RsdCG |psi[k+1]| then return
  for (k=1; k<=cg_arg.max_num_iter && !convP; k++) {
    // a[k+1] = |r[k]**2/ |r[k-1]|**2
    a = c/cp;
    VRB.Flow(cname,fname,"a =%e, |r[%d]]^2 = %e\n",a,k,c);

    // p[k+1] = r[k+1] + a[k+1] p[k]
    //   Compute the shifted as
    //   ps[k+1] = zs[k+1] r[k+1] + a[k+1] ps[k]
    for (s=0; s<Nmass; s++) {
      if (convsP[s]) continue;
      if (s==0) {
	p[s] -> FTimesV1PlusV2(a,p[s],r,f_size);
      } else {
	as = a * (z[iz][s] * bs[s]) / (z[1-iz][s] * b);
	p[s] -> VecTimesEquFloat(as,f_size);	
	p[s] -> FTimesV1PlusV2(z[iz][s],r,p[s],f_size);	
      }
    }
    
    // cp = |r[k]**2
    cp = c;
    
    // b[k] = |r[k]**2 / <p[k], Ap[k]>
    if (mass[0] > 0) {
      dop.MatPcDagMatPc(Ap,p[0],dot);
      Ap -> FTimesV1PlusV2(mass[0],p[0],Ap, f_size);
      d = p[0] -> ReDotProductGlbSum(Ap, f_size);
    } else {
      dop.MatPcDagMatPc(Ap,p[0],&d);
      glb_sum(&d);
    }

    bp = b;
    b = -cp/d;
    
    //Compute the shifted bs and z
    bs[0] = b;
    iz = 1 - iz;
    for (s=0; s<Nmass; s++) {
      if (s==0 || convsP[s]) continue;      
      ztmp = z[1-iz][s]*z[iz][s]*bp / 
	( b*a*(z[iz][s]-z[1-iz][s]) + z[iz][s]*bp*(1-b*(mass[s] - mass[0])));
      bs[s] = b*ztmp / z[1-iz][s];
      z[iz][s] = ztmp;
    }
    
    // r[k+1] += b[k] A.p[k]
    r -> FTimesV1PlusV2(b,Ap,r,f_size);
    // c = |r[k]|**2
    c = r-> NormSqGlbSum(f_size);

    // p[k+1] = r[k+1] + a[k+1] p[k]
    //   Compute the shifted as
    //   ps[k+1] = zs[k+1] r[k+1] + a[k+1] ps[k]
    // Psi[k+1] -= b[k] p[k]

    if (type == SINGLE)
      for (s=0; s<Nmass; s++) {
	if (convsP[s]){ 
	Float *tmp_p = (Float *)psi[0];
	VRB.Result(cname,fname,"bs[%d]=%g psi[%d]=%g\n",s,bs[s],s,*tmp_p);
 	continue;}
	psi[0]->FTimesV1PlusV2(-bs[s]*alpha[s],p[s],psi[0],f_size);
      }
    else
      for (s=0; s<Nmass; s++) {
	if (convsP[s]){ 
 	continue;}
	psi[s]->FTimesV1PlusV2(-bs[s],p[s],psi[s],f_size);
      }
    
    // if |psi[k+1] -psi[k]| <= rsdCG |psi[k+1]| then return
    // or if |r[k+1]| <= RsdCG |chi| then return
    for (s=0; s<Nmass; s++) {
      if (convsP[s]) continue;
      //check norm of shifted residuals
      css = c * z[iz][s] * z[iz][s];
      convsP[s] = (css < rsd_sq[s]) ? 1 : 0;
      if (convsP[s]){
	RsdCG[s] = css;
	converged[s] = k;
      }
    }    
    
    convP = convsP[0];
    // if zero solution has converged, exit unless other solutions have not
    if (convP) for (s=0; s<Nmass; s++) if (!convsP[s]) convP = 0;
  }
  
  if (k >= cg_arg.max_num_iter)   // It has not reached stp_cnd: Issue a warning
    VRB.Warn(cname,fname,"CG reached max iterations = %d. |res|^2 = %e\n",k, css);

  for (s=Nmass-1; s>=0; s--) {
    if (convsP[s]) {
      VRB.Result(cname,fname,"%d shift converged, iter = %d, res^2 = %e\n",
		 s+isz,converged[s],RsdCG[s]);
      RsdCG[s] = sqrt(RsdCG[s]);
    } else {
      RsdCG[s] = c*z[iz][s]*z[iz][s];
      VRB.Result(cname,fname,
		 "%d shift did not converge, iter = %d, res^2 = %e\n",
		 s+isz,k,RsdCG[s]);
      RsdCG[s] = sqrt(RsdCG[s]);
    }
  }

  // free arrays and vectors
  for (s=0; s<Nmass; s++)
    sfree(*(p+s), cname, fname, "p[s]");
    
  sfree(p, cname, fname, "p");
  sfree(Ap, cname, fname, "Ap");
  sfree(r, cname, fname, "r");
  sfree(bs);
  sfree(*(z+1));
  sfree(*z);
  sfree(z);
  sfree(converged);
  sfree(convsP);
  sfree(rsdcg_sq);
  sfree(rsd_sq);
  
  return k;
}


int CPS_FeigSolv(Vector **f_eigenv, Float *lambda,
		 Float *chirality, int *valid_eig,
		 Float **hsum,
		 EigArg *eig_arg, 
		 CnvFrmType cnv_frm,
		 Lattice &lat)
{
  const char* cname = "";
  const char *fname = "FeigSolv(EigArg*,V*,F*,CnvFrmType) [CPS version]";
  VRB.Func(cname,fname);

  //CK: The regular CPS version

  int iter;

  CgArg cg_arg;
  cg_arg.mass = eig_arg->mass;
  cg_arg.epsilon = eig_arg->epsilon; //CK: passes down epsilon parameter for twisted mass Wilson fermions. Irrelevant here but this same function is used in FwilsonTm
  cg_arg.RitzMatOper = eig_arg->RitzMatOper;
  int N_eig = eig_arg->N_eig;
  int i;

  //=========================
  // convert fermion field
  //=========================

  if(cnv_frm == CNV_FRM_YES) //Fixed by CK to allow for single checkerboard input vectors as used in AlgActionRational. Previously it would always convert from CANONICAL to WILSON
    for(i=0; i < N_eig; ++i)  lat.Fconvert(f_eigenv[i], WILSON, CANONICAL);

  //------------------------------------------------------------------
  //  we want both the eigenvalues of D_{hermitian} and
  //  D^{+}D.  To not change the arguments passed to RitzEig,
  //  we pass a float pointer which points to 2 * N_eig values
  //  and return both lambda and lambda^2 from RitzEig
  //------------------------------------------------------------------

  Float * lambda2 = (Float * ) smalloc (N_eig*2*sizeof(Float));
  if ( lambda2 == 0 ) ERR.Pointer(cname,fname, "lambda2");
  
  {
    DiracOpWilsonTm wilson(lat, (Vector*) 0 , (Vector*) 0, &cg_arg, CNV_FRM_NO);
    iter = wilson.RitzEig(f_eigenv, lambda2, valid_eig, eig_arg);
  }

  if(cnv_frm == CNV_FRM_YES) 
    for(i=0; i < N_eig; ++i) lat.Fconvert(f_eigenv[i], CANONICAL, WILSON);


  /*
    the call to RitzEig returns a negative number if either the KS or CG maxes
    out, we wish to cope with this in alg_eig, so "pass it up". Clean up the
    storage order first in case we still want to use the eigenvectors as a
    guess.
  */
  if ( iter < 0 ) { return iter ; }


  // Compute chirality
  int Ncb = NumChkb(cg_arg.RitzMatOper);
  size_t f_size = (GJP.VolNodeSites() * lat.FsiteSize()) * Ncb / 2; //CK: fixed
  if(GJP.Gparity()) f_size *= 2;

  Vector* v1 = (Vector *)smalloc(f_size*sizeof(Float));
  if (v1 == 0)
    ERR.Pointer(cname, fname, "v1");
  VRB.Smalloc(cname, fname, "v1", v1, f_size*sizeof(Float));

  int nspinvect = GJP.VolNodeSites() * Ncb/2;
  if(GJP.Gparity()) nspinvect *= 2;

  for(i=0; i < N_eig; ++i)
  {
    lat.Gamma5(v1, f_eigenv[i], nspinvect);
    chirality[i] = f_eigenv[i]->ReDotProductGlbSum4D(v1, f_size);
  }

  VRB.Sfree(cname, fname, "v1", v1);
  sfree(v1);


  // rescale wilson eigenvalues to the convention  m + Dslash(U)
  Float factor = 4.0 + eig_arg->mass;
    
  
  FILE* fp=Fopen(eig_arg->fname,"a");
  for(i=0; i<N_eig; ++i)
    {
      lambda2[i] *= factor;	 		 //rescale eigenvalue
      lambda2[N_eig + i] *= ( factor * factor ); //rescale squared evalue
      lambda[i]=lambda2[i];                      //copy back
      
      //print out eigenvalue, eigenvalue^2, chirality 
      Fprintf(fp,"%d %g %g %g %d\n",i,
              (float)lambda2[i],
              (float)lambda2[N_eig + i],
	      (float)chirality[i],valid_eig[i]);
    }
  Fclose(fp);
  sfree(lambda2); 


  // Slice-sum the eigenvector density to make a 1D vector
  if (eig_arg->print_hsum){
    for(i=0; i < N_eig; ++i){
      //CK: The vector needs to be in canonical ordering. Thus if CNV_FRM_NO we need to convert to CANONICAL. If Ncb==1 we have
      //    only the odd part, so we will need to fill the even part with zeroes prior to conversion
      Vector* tosum = f_eigenv[i];
      if(cnv_frm == CNV_FRM_NO){
	//Create a temp copy of the eigenvector
	int alloc_size = f_size; if(Ncb==1) alloc_size *= 2;
	Float* full = (Float *)smalloc(alloc_size*sizeof(Float));
	for(int j=0;j<f_size;j++) full[j] = ((Float*)f_eigenv[i])[j];

	//Fill in even part with zero for Ncb==1
	if(Ncb==1) for(int j=f_size;j<alloc_size;j++) full[j] = 0; //zero even part

	//Convert
	if(cnv_frm == CNV_FRM_NO) lat.Fconvert((Vector*)full, CANONICAL, WILSON);
	tosum = (Vector*)full;
      }
      tosum->NormSqArraySliceSum(hsum[i], lat.FsiteSize(), eig_arg->hsum_dir);

      if(cnv_frm == CNV_FRM_NO) sfree(tosum);

    }
    
  }

  // The remaining part in QCDSP version are all about "downloading
  // eigenvectors", supposedly not applicable here.

  // Return the number of iterations
  return iter;
}


class AlgEig_CPSeigsolver : public Alg
{
 private:
    char *cname;

    EigArg *alg_eig_arg;
        // The argument structure for the eig algorithm
 
    int Ncb;       
        // Number of checkerboards for fermion field (1 or 2)

    Vector **eigenv;
        // The eigenvectors (initial and final)

    Float *lambda;
        // The eigenvalues (final)

    Float *chirality;
        // The chirality of the eigenvalues (final)

    int *valid_eig;
        // Whether the eigenvalues are valid or not (final)

    int n_masses;
    // The number of masses in the loop
 public:
    AlgEig_CPSeigsolver(Lattice & latt, CommonArg *c_arg, EigArg *arg);

    virtual ~AlgEig_CPSeigsolver();

    void run(void);
    void run(Float **lambda);

    //Added by CK. Call these only *after* running
    Vector ** getEigenVectors();
};

static void convertSorderToCanonical(Vector* vect, Fbfm* fbfm){
  Fermion_t handle[2] = { fbfm->bd.allocFermion(), fbfm->bd.allocFermion() };
  fbfm->bd.cps_impexFermion_s((Float*)vect,handle,1);
  fbfm->bd.cps_impexFermion((Float*)vect, handle, 0);
  fbfm->bd.freeFermion(handle[0]);
  fbfm->bd.freeFermion(handle[1]);
}
// static void convertFullCbFermToCanonical(Vector* vect, Fbfm* fbfm){
//   size_t f_size_cb = GJP.VolNodeSites()*24/2; if(GJP.Gparity()) f_size_cb*=2;
//   Fermion_t handle[2] = { fbfm->bd.allocFermion(), fbfm->bf.allocFermion() };
//   fbfm->bd.cps_impexcbFermion((Float*)vect,handle[0],1,1);
//   fbfm->bd.cps_impexcbFermion((Float*)vect + f_size_cb,handle[1],1,0);
//   fbfm->bd.cps_impexFermion((Float*)vect, handle[2], 0);
//   fbfm->bd.freeFermion(handle[0]);
//   fbfm->bd.freeFermion(handle[1]);
// }


//Note: calculation of HMD force vecs, v2 is the one that contains Deo^dag and Mprec chi
//In BFM version, v1 is the one that contains Deo^dag and Mprec chi
ForceArg CPS_EvolveMomFforce(Matrix *mom, Vector *v1, Vector* v2, 
			     Float mass, Float epsilon, Float dt, FwilsonTm* latt)
{
  const char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
  const char *cname = "";
  Matrix *gauge = latt->GaugeField() ;

  size_t f_size = latt->FsiteSize() * GJP.VolNodeSites() ;

  Float *site_v1 = (Float *)pmalloc(latt->FsiteSize()*sizeof(Float));
  Float *site_v2 = (Float *)pmalloc(latt->FsiteSize()*sizeof(Float));

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;

  int x, y, z, t, lx, ly, lz, lt ;

  lx = GJP.XnodeSites() ;
  ly = GJP.YnodeSites() ;
  lz = GJP.ZnodeSites() ;
  lt = GJP.TnodeSites() ;

  int lattsz[4] = {lx,ly,lz,lt};

  Matrix tmp, f ;

  for (int mu=0; mu<4; mu++) {
    for (t=0; t<lt; t++)
    for (z=0; z<lz; z++)
    for (y=0; y<ly; y++)
    for (x=0; x<lx; x++) {
      int gauge_offset = x+lx*(y+ly*(z+lz*t)) ;
      int vec_offset = latt->FsiteSize()*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;

      Float *v1_plus_mu ;
      Float *v2_plus_mu ;
      int vec_plus_mu_offset = latt->FsiteSize() ;

      Float coeff = -2.0 * dt ;
      
      /* CK: Replaced nasty 300000 line switch statement with the following:*/
      {
	int pos_p_mu[] = {x,y,z,t};
	pos_p_mu[mu] = (pos_p_mu[mu]+1) % lattsz[mu];
	vec_plus_mu_offset *= pos_p_mu[0] + lx*(pos_p_mu[1] + ly*(pos_p_mu[2] + lz*pos_p_mu[3]));
      }      
      int pos[4] = {x,y,z,t};

      if ((pos[mu]+1) == lattsz[mu]) {
	getPlusData( (IFloat *)site_v1,
		     (IFloat *)v1+vec_plus_mu_offset, latt->FsiteSize(), mu) ;
	getPlusData( (IFloat *)site_v2,
		     (IFloat *)v2+vec_plus_mu_offset, latt->FsiteSize(), mu) ;
	v1_plus_mu = site_v1 ;                        
	v2_plus_mu = site_v2 ;                        
	if (GJP.NodeBc(mu)==BND_CND_APRD) coeff = -coeff ;
      } else {
	v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
	v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
      }

      latt->sproj_tr[mu](   (IFloat *)&tmp,
                      (IFloat *)v1_plus_mu,
                      (IFloat *)v2+vec_offset, 1, 0, 0);

      latt->sproj_tr[mu+4]( (IFloat *)&f,
                      (IFloat *)v2_plus_mu,
                      (IFloat *)v1+vec_offset, 1, 0, 0);

      tmp += f ;

      f.DotMEqual(*(gauge+gauge_offset), tmp) ;

      tmp.Dagger(f) ;
      f.TrLessAntiHermMatrix(tmp) ;
      f *= coeff ;

      *(mom+gauge_offset) += f ;

      Float norm = f.norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }
  }


  pfree(site_v2) ;
  pfree(site_v1) ;

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();
  VRB.FuncEnd(cname,fname);

  return ForceArg(L1, sqrt(L2), Linf);
}


static int no_gparity_test(GnoneFwilsonTm* lattice){
  int mom_size = 18*4*GJP.VolNodeSites();

  Float* mom_test1 = (Float *)pmalloc( sizeof(Float) * mom_size);
  for(int i=0;i<mom_size;i++) mom_test1[i] = 0.0;

  Float* mom_test2 = (Float *)pmalloc( sizeof(Float) * mom_size);
  for(int i=0;i<mom_size;i++) mom_test2[i] = 0.0;

  Float* v1 = rand_4d_canonical_fermion(*lattice);

#if 0
  {
    size_t f_size = 24*GJP.VolNodeSites();
    Vector* f_in_cps = (Vector*)pmalloc(f_size * sizeof(Float));
    Vector* f_out_cps = (Vector*)pmalloc(f_size * sizeof(Float));
    Vector* f_in_bfm = (Vector*)pmalloc(f_size * sizeof(Float));
    Vector* f_out_bfm = (Vector*)pmalloc(f_size * sizeof(Float));
    
    CgArg cg_arg;
    cg_arg.mass = -1.8;
    cg_arg.epsilon = 0.5;
    
    DiracOpWilsonTm dop(*lattice,f_out_bfm,f_in_bfm,&cg_arg,CNV_FRM_NO);

    double time_bfm = -dclock();
    dop.CalcHmdForceVecs((Vector*)v1);
    time_bfm += dclock();

    double time_cps = -dclock();
    CPS_CalcHmdForceVecs((Vector*)v1, f_out_cps, f_in_cps, dop);
    time_cps += dclock();

    printf("MD force vec calc times: CPS %.9e  BFM %.9e\n",time_cps,time_bfm);
    printf("Kappa : %e\n",dop.get_kappa());

    bool fail = false;

    for(int i=0;i<f_size;i++){
      if(fabs(((Float*)f_in_cps)[i]-((Float*)f_in_bfm)[i])>1e-9){
	printf("f_in comparison fail %d: %.9e %.9e, ratio %.9e\n",i, ((Float*)f_in_cps)[i], ((Float*)f_in_bfm)[i], ((Float*)f_in_cps)[i]/((Float*)f_in_bfm)[i]);
	fail = true;
      }else printf("f_in comparison pass %d: %.9e %.9e, ratio %.9e\n",i, ((Float*)f_in_cps)[i], ((Float*)f_in_bfm)[i], ((Float*)f_in_cps)[i]/((Float*)f_in_bfm)[i]);
      if(i >= f_size/2 && fabs(((Float*)f_out_cps)[i]-((Float*)f_out_bfm)[i])>1e-9){
	printf("f_out comparison fail %d: %.9e %.9e, ratio %.9e\n",i, ((Float*)f_out_cps)[i], ((Float*)f_out_bfm)[i], ((Float*)f_out_cps)[i]/((Float*)f_out_bfm)[i]);
	fail = true;
      }else if(i >= f_size/2) printf("f_out comparison pass %d: %.9e %.9e, ratio %.9e\n",i, ((Float*)f_out_cps)[i], ((Float*)f_out_bfm)[i], ((Float*)f_out_cps)[i]/((Float*)f_out_bfm)[i]);
    }
    if(fail){ printf("Failed MD force vec test\n"); exit(-1); }
    else printf("Passed MD force vec test\n");
  }
#endif



  oldTm_EvolveMomFforce_nogp((Matrix*)mom_test1, (Vector*)v1, 0.5, 0.3, 0.1, lattice);  
  lattice->EvolveMomFforce((Matrix*)mom_test2, (Vector*)v1, 0.5, 0.3, 0.1); 

  bool fail = false;
  for(int i=0;i<mom_size;i++){
    int rem = i;
    int midx = rem % 18; rem/=18;
    int mu = rem % 4; rem/=4;

    int x[4];
    for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

    if(fabs(mom_test2[i] - mom_test1[i])>1e-08){
      printf("Fail EvolveMomFforce improved code test fail midx %d mu %d (%d %d %d %d): %f %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test2[i],mom_test1[i]);
      fail=true;
    }//else printf("Pass EvolveMomFforce improved code test midx %d mu %d (%d %d %d %d): %f %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test2[i],mom_test1[i]);
  }
  if(fail){ printf("Failed EvolveMomFforce improved code test\n"); exit(-1); }
  else printf("Passed EvolveMomFforce improved code test\n");

  {
    //DO a convert test
    lattice->Fconvert((Vector*)v1, WILSON, CANONICAL);
    lattice->Fconvert((Vector*)v1, CANONICAL, WILSON);
    printf("Got through convert test\n"); fflush(stdout);
  }
  

  //I had to add FeigSolve and RHMC_EvolveMomFforce, as well as modify AlgActionRationalQuotient for twisted fermions, as these had not been done. Test them here
  //Also added a BFM version which we also test here.
  {
    //Test eigenvalues
    CommonArg c_arg;
    c_arg.set_filename("pooh");
    c_arg.set_label("pooh");
   
    Float mass = -0.6;
    Float epsilon = 0.5;

    EigArg e_arg;
    
    generateEigArg(e_arg, mass, epsilon); //operator is MatPcDagMatPc
    e_arg.Rsdlam = 1e-6;
    e_arg.RsdR_a = e_arg.Rsdlam ;
    e_arg.RsdR_r = e_arg.Rsdlam ;
    //e_arg.N_eig = 1;

    CgArg cg_arg;
    cg_arg.mass = mass;
    cg_arg.epsilon = epsilon;


    Float **evalues = (Float**) pmalloc(2*sizeof(Float*)); //indexed by eigenvalue index, then mass
    for(int i=0;i<2;i++) evalues[i] = (Float*)pmalloc(sizeof(Float));

    LatRanGen LRGbak = LRG;

    AlgEig_CPSeigsolver algeig(*lattice,&c_arg, &e_arg);
    algeig.run(evalues);
    
    LRG = LRGbak;

    Vector **eigenvectors = algeig.getEigenVectors();

#ifdef USE_BFM_TM
    {
      //Check that BFM version gets the same result
      //Can only do one eigenvalue here

      Float **evalues_bfm = (Float**) pmalloc(sizeof(Float*)); //indexed by eigenvalue index, then mass
      for(int i=0;i<1;i++) evalues_bfm[i] = (Float*)pmalloc(sizeof(Float));

      e_arg.N_eig = 1;

      AlgEig algeig_bfm(*lattice,&c_arg, &e_arg);
      algeig_bfm.run(evalues_bfm);

      LRG = LRGbak;

      Vector **eigenvectors_bfm = algeig_bfm.getEigenVectors();
      size_t f_size = GJP.VolNodeSites()*24;
      
      bool fail = false;

      for(long i=0;i<f_size;i++){
	int rem = i;
	int m = rem % 24; rem/=24;
	int pos[4];
	for(int j=0;j<4;j++){ pos[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem/=GJP.NodeSites(j); }

	if( fabs( ((Float*)eigenvectors[0])[i] - ((Float*)eigenvectors_bfm[0])[i]) > 5e-05 ){ 
	  printf("**Eigenvector CPS/BFM test fail, %d %d %d %d, m=%d : %f %f, ratio %f\n",pos[0],pos[1],pos[2],pos[3],
		 m,((Float*)eigenvectors[0])[i],((Float*)eigenvectors_bfm[0])[i], 
		 ((Float*)eigenvectors[0])[i]/((Float*)eigenvectors_bfm[0])[i]); 
	  fail=true; 
	}
      }
      if(fail){
	printf("Eigenvector CPS/BFM test failed\n"); //exit(-1);
      }else printf("Eigenvector CPS/BFM test passed\n");


      fail = false;

      for(int i=0;i<1;i++){
	if( fabs(evalues[i][0] - evalues_bfm[i][0])>5e-05 ){
	  printf("Failed BFM/CPS Eig test %d: %f %f, ratio %f\n",i,evalues[i][0],evalues_bfm[i][0],evalues[i][0]/evalues_bfm[i][0] );
	  fail = true;
	}
      }
      if(fail){
	printf("Failed BFM/CPS Eig test");
	exit(-1);
      }else printf("Passed BFM/CPS Eig test");

    }
#endif


    for(int ev = 0; ev < 2; ev++){
      //lattice->Fconvert(eigenvectors[ev], WILSON, CANONICAL);
    
      size_t f_size = GJP.VolNodeSites()*24/2; //only on odd checkerboard
      Vector* MdagM_v = (Vector*)pmalloc(f_size * sizeof(Float));
      Vector* lambda_v = (Vector*)pmalloc(f_size * sizeof(Float));

      //Apply MdagM to eigenvector
      {
	DiracOpWilsonTm wilson(*lattice, (Vector*) 0 , (Vector*) 0, &cg_arg, CNV_FRM_NO);
	wilson.MatPcDagMatPc(MdagM_v, eigenvectors[ev],0); 
      }      

      //Multiply eigenvector by lambda
      //Note, eigenvalues are multiplied internally by a normalisation factor of 4.0 + mass, which we remove for this comparison
      Float invnorm = 1.0/(4.0+cg_arg.mass);
      for(int ii=0;ii<f_size;ii++) ((Float*)lambda_v)[ii] = evalues[ev][0] * invnorm * ((Float*)eigenvectors[ev])[ii];

      //Note, these vectors are single-checkerboard WILSON ordered. Convert to canonical
      Float* MdagM_v_f = convertCanonicalSingleChkb((Float*)MdagM_v,1,*lattice);
      Float* lambda_v_f = convertCanonicalSingleChkb((Float*)lambda_v,1,*lattice);

      Float eps = 5e-5;

      //Compare
      fail = false;
      for(long i=0;i<f_size*2;i++){
	int rem = i;
	int m = rem % 24; rem/=24;
	int pos[4];
	for(int j=0;j<4;j++){ pos[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem/=GJP.NodeSites(j); }

	if( fabs(MdagM_v_f[i] - lambda_v_f[i]) > eps ){ printf("**Eigen test %d fail, %d %d %d %d, m=%d : %f %f, ratio %f\n",ev,pos[0],pos[1],pos[2],pos[3],
								m,MdagM_v_f[i],lambda_v_f[i], MdagM_v_f[i]/lambda_v_f[i]); fail=true; 
	}
      }
      if(fail){
	printf("Eigen test %d failed\n",ev); exit(-1);
      }else printf("Eigen test %d passed\n",ev);

      pfree(MdagM_v);
      pfree(lambda_v);
      pfree(MdagM_v_f);
      pfree(lambda_v_f);
      printf("Finished deleting memory for test\n"); fflush(stdout);
    }

    for(int i=0;i<2;i++) pfree(evalues[i]);
    pfree(evalues);

  }

  //See if we can coax Fbfm into doing twisted mass fermions
  {
    Float mass = 0.5;
    Float epsilon = 0.5;

    Float* v2 = rand_4d_canonical_fermion(*lattice);
    delete lattice;
    Fbfm::current_arg_idx = 0;
    Fbfm::bfm_args[0].solver = WilsonTM;
    Fbfm::bfm_args[0].node_latt[0]  = QDP::Layout::subgridLattSize()[0];
    Fbfm::bfm_args[0].node_latt[1]  = QDP::Layout::subgridLattSize()[1];
    Fbfm::bfm_args[0].node_latt[2]  = QDP::Layout::subgridLattSize()[2];
    Fbfm::bfm_args[0].node_latt[3]  = QDP::Layout::subgridLattSize()[3];
    Fbfm::bfm_args[0].verbose=1;
    Fbfm::bfm_args[0].reproduce=0;

#if TARGET == BGQ
    omp_set_num_threads(64);
#else 
    omp_set_num_threads(1);
#endif

#if TARGET == BGQ
    bfmarg::Threads(64);
#else
    bfmarg::Threads(1);
#endif

    bfmarg::Reproduce(0);
    bfmarg::ReproduceChecksum(0);
    bfmarg::ReproduceMasterCheck(0);
    bfmarg::Verbose(1);
    bfmarg::onepluskappanorm = 0;

    multi1d<int> procs = QDP::Layout::logicalSize();

    printf("%d dim machine\n\t", procs.size());
    for(int mu=0;mu<4;mu++){
      printf("%d ", procs[mu]);
      if ( procs[mu]>1 ) Fbfm::bfm_args[0].local_comm[mu] = 0;
      else Fbfm::bfm_args[0].local_comm[mu] = 1;
    }
    printf("\nLocal comm = ");
    for(int mu=0;mu<4;mu++){
      printf("%d ", Fbfm::bfm_args[0].local_comm[mu]);
    }
    printf("\n");

    Fbfm::bfm_args[0].precon_5d = 0;
    Fbfm::bfm_args[0].Ls   = 1;
    Fbfm::bfm_args[0].M5   = 0.0;
    Fbfm::bfm_args[0].mass = toDouble(mass);
    Fbfm::bfm_args[0].twistedmass = toDouble(epsilon);
    Fbfm::bfm_args[0].Csw  = 0.0;
    Fbfm::bfm_args[0].max_iter = 10000;
    Fbfm::bfm_args[0].residual = 1e-08;

    size_t f_size = 24*GJP.VolNodeSites();
    

    //Test 1: check we get the same HMD force vectors (up to the appropriate normalization difference)

    //CPS version:  
    //v1 = f_field_out becomes (phi1, g5theta(ctheta,-stheta)D_eo phi1)
    //v2 = f_field_in  becomes ( -kappa^2 g5theta(ctheta,stheta)phi2, -kappa^2 g5theta(ctheta,stheta) Deo^dag g5theta(ctheta,stheta)phi2 )

    //kappa = 1/[2 sqrt( (m+4)^2 + eps^2 )]
    //ctheta = 2 (m+4) kappa
    //stheta = 2 eps kappa
    //g5theta(ctheta,stheta) = ctheta + i stheta g5

    //g5theta(ctheta,stheta)^2 = g5theta(ctheta,-stheta)^2 = ctheta^2 + stheta^2 = 4 kappa^2 ( (m+4)^2 + eps^2 ) = 4 kappa^2 * 1/4 1/kappa^2 = 1

    //BFM version:
    //v1 = (phi1,  MeeInv^dag Meo^dag phi1)
    //v2 = (Boo phi2, Bee MeeInv Meo phi2)

 // * WilsonTM : Mee = Moo = 4+m + i tm g5
 // *            Meeinv    = (4+m - i tm g5)/ ( (4+m)^2 + tm^2 ) 
 // *                      =    (4+m)
 // *                        -----------------   .  ( 1 -i tm/(4+m) g5 )
 // *                        ( (4+m)^2 + tm^2 )
 // *
    //'tm' is same as epsilon
    //Hence Meeinv = 1/(2 kappa) (ctheta - i stheta g5) * 4 kappa^2 = 2 kappa g5theta(ctheta,-stheta)
    //I think Boo, Bee should be 1. According to bfm comments, Meo = -1/2 Deo
    //BFM version is thus
    //v1 = (phi1,  -kappa g5theta(ctheta,stheta) Deo^dag phi1)
    //v2 = (phi2,  -kappa g5theta(ctheta,-stheta) Deo phi2)

    //Comparing to CPS its clear we need to swap v1 and v2 and phi1 and phi2 in the input such that the BFM version gives
    //v1 = (phi1,  -kappa g5theta(ctheta,-stheta) Deo phi1)
    //v2 = (phi2,  -kappa g5theta(ctheta,stheta) Deo^dag phi2)
    
    //we also need to multiply phi2 by g5theta(ctheta,stheta) prior to placing into bfm version

    //then we have
    //v1_odd_CPS = v1_odd_bfm
    //v1_even_CPS = -1/kappa v1_even_BFM
    //v2_odd_CPS = -kappa^2 v2_odd_BFM
    //v2_even_CPS = kappa v2_even_BFM

    Float kappa;
    Float ctheta;
    Float stheta;

    Vector* forcevec1_cps = (Vector*)pmalloc(f_size * sizeof(Float));
    Vector* forcevec2_cps = (Vector*)pmalloc(f_size * sizeof(Float));
    {
      GnoneFwilsonTm *latt_cps = new GnoneFwilsonTm;
      CgArg cg_arg;
      cg_arg.mass = mass;
      cg_arg.epsilon = epsilon;
      {
	DiracOpWilsonTm wilson(*latt_cps, forcevec1_cps, forcevec2_cps, &cg_arg, CNV_FRM_YES);
	wilson.CalcBsnForceVecs((Vector*)v1,(Vector*)v2);
	//v1 and v2 are automatically converted to canonical by dirac op destructor
	kappa = wilson.get_kappa();
	ctheta = wilson.get_ctheta();
	stheta = wilson.get_stheta();

      }
      delete latt_cps;
    }
    printf("WilsonTM parameters: mass %f, epsilon %f, kappa %.9e ctheta %.9e stheta %.9e\n",mass,epsilon,kappa,ctheta,stheta);

    Vector* forcevec1_fbfm = (Vector*)pmalloc(f_size * sizeof(Float));
    Vector* forcevec2_fbfm = (Vector*)pmalloc(f_size * sizeof(Float));
    {
      Vector* v2_times_g5theta = (Vector*)pmalloc(f_size/2 * sizeof(Float));
      memcpy((void*)v2_times_g5theta,(void*)v2,f_size/2*sizeof(Float));
      g5theta((Vector*)v2_times_g5theta,GJP.VolNodeSites()/2,ctheta,stheta);

      GnoneFbfm *latt_fbfm = new GnoneFbfm;
      //technically input vectors should be odd-checkerboard rather than canonical ordered, but as they are random it doesnt matter
      dynamic_cast<Fbfm*>(latt_fbfm)->CalcHmdForceVecsBilinear((Float*)forcevec2_fbfm, (Float*)forcevec1_fbfm, 
							       v2_times_g5theta,(Vector*)v1,mass,epsilon);
      //fermions come out in s-order, convert to canonical
      convertSorderToCanonical(forcevec1_fbfm, dynamic_cast<Fbfm*>(latt_fbfm));
      convertSorderToCanonical(forcevec2_fbfm, dynamic_cast<Fbfm*>(latt_fbfm));
      delete latt_fbfm;
    }

    for(int i=0;i<f_size/6;i++){
      int rem = i;
      int sidx = rem%4; rem/=4;
      int pos[4];
      for(int j=0;j<4;j++){ pos[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem/=GJP.NodeSites(j); }
      
      if( (pos[0]+pos[1]+pos[2]+pos[3]) % 2 == 0){ //even
	forcevec1_fbfm[i] *= -1.0/kappa;
	forcevec2_fbfm[i] *= kappa;
      }else{ //odd
	forcevec2_fbfm[i] *= -kappa*kappa;
      }
    }



    for(long i=0;i<f_size;i++){
      int rem = i;
      int m = rem % 24; rem/=24;
      int pos[4];
      for(int j=0;j<4;j++){ pos[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem/=GJP.NodeSites(j); }

      if( fabs( ((Float*)forcevec1_cps)[i] - ((Float*)forcevec1_fbfm)[i]) > 1e-08 ){ 
	printf("**HMD force vec 1 CPS/FBFM test fail, %d %d %d %d, m=%d : %.9e %.9e, ratio %f\n",pos[0],pos[1],pos[2],pos[3],
	       m,((Float*)forcevec1_cps)[i],((Float*)forcevec1_fbfm)[i], 
	       ((Float*)forcevec1_cps)[i]/((Float*)forcevec1_fbfm)[i]);
	fail=true; 
      }

      if( fabs( ((Float*)forcevec2_cps)[i] - ((Float*)forcevec2_fbfm)[i]) > 1e-08 ){ 
	printf("**HMD force vec 2 CPS/FBFM test fail, %d %d %d %d, m=%d : %.9e %.9e, ratio %f\n",pos[0],pos[1],pos[2],pos[3],
	       m,((Float*)forcevec2_cps)[i],((Float*)forcevec2_fbfm)[i], 
	       ((Float*)forcevec2_cps)[i]/((Float*)forcevec2_fbfm)[i]);
	fail=true; 
      }
    }
    if(fail){
      printf("HMD force vec CPS/FBFM test failed\n"); exit(-1);
    }else printf("HMD force vec CPS/FBFM test passed\n");


    //BFM preconditioned matrix:  chi = Mprec psi = [ Moo -              1/4                      Doe        MeeInv            Deo ] psi
    //CPS preconditioned matrix:  chi = Mprec psi = [  1  -      kappa^2 g5theta(ctheta, -stheta) Doe g5theta(ctheta, -stheta) Deo ] psi

    //Moo = 4+m + i tm g5 = 1/(2 kappa) g5theta(ctheta,stheta)
    //MeeInv = 2 kappa g5theta(ctheta,-stheta)

    //BFM preconditioned matrix:  chi = Mprec psi = [ 1/(2 kappa) g5theta(ctheta,stheta) - kappa/2 Doe g5theta(ctheta,-stheta) Deo ] psi
    //                                            = [1/(2 kappa)] g5theta(ctheta,stheta) (1 - kappa^2 g5theta(ctheta, -stheta) Doe g5theta(ctheta, -stheta) Deo ] psi
    //                                            = [1/(2 kappa)] g5theta(ctheta,stheta) Mprec_CPS psi

    //(I think this extra factor of g5theta(ctheta,stheta) makes up for the lack of this factor in the HMD force vector calculation)
    
    //Confirm this normalization
    
    Vector* chi_CPS;
    {
      GnoneFwilsonTm *latt_cps = new GnoneFwilsonTm;
      CgArg cg_arg;
      cg_arg.mass = mass;
      cg_arg.epsilon = epsilon;
      {
	Vector* tmp_chi_CPS = (Vector*)pmalloc(f_size/2 * sizeof(Float)); //odd checkerboard only

	DiracOpWilsonTm wilson(*latt_cps, (Vector*)NULL, (Vector*)NULL, &cg_arg, CNV_FRM_NO);
	wilson.MatPc(tmp_chi_CPS, (Vector*)v1);
	chi_CPS = (Vector*)convertCanonicalSingleChkb((Float*)tmp_chi_CPS, 1, *latt_cps);
	
	pfree(tmp_chi_CPS);
      }
      delete latt_cps;
    }

    Vector* chi_BFM;
    {
      Vector* tmp_chi_BFM = (Vector*)pmalloc(f_size/2 * sizeof(Float));
      GnoneFbfm *latt_fbfm = new GnoneFbfm;
      dynamic_cast<Fbfm*>(latt_fbfm)->MatPc(tmp_chi_BFM, (Vector*)v1, mass, epsilon, DAG_NO);
      //Fbfm doesn't have an implementation of Fconvert, thus delete the lattice and make a 
      delete latt_fbfm;
      GnoneFwilsonTm *latt_cps = new GnoneFwilsonTm;
      g5theta(tmp_chi_BFM,GJP.VolNodeSites()/2,ctheta,-stheta); //part 1 of the normalization

      chi_BFM = (Vector*)convertCanonicalSingleChkb((Float*)tmp_chi_BFM, 1, *latt_cps);
      pfree(tmp_chi_BFM);
      delete latt_cps;
    }

    for(long i=0;i<f_size;i++){
      int rem = i;
      int m = rem % 24; rem/=24;
      int pos[4];
      for(int j=0;j<4;j++){ pos[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem/=GJP.NodeSites(j); }

      Float cps = ((Float*)chi_CPS)[i];
      Float bfm = ((Float*)chi_BFM)[i]*2*kappa;


      if( fabs(cps-bfm) > 1e-08 ){ 
	printf("**Mpc norm CPS/FBFM test fail, %d %d %d %d, m=%d : %.9e %.9e, ratio %f\n",pos[0],pos[1],pos[2],pos[3],
	       m,cps,bfm,cps/bfm);
	fail=true; 
      }
    }
    if(fail){
      printf("Mpc norm CPS/FBFM test failed\n"); exit(-1);
    }else printf("Mpc norm CPS/FBFM test passed\n");

    //Check normalization also for Mpc^dag
    //Mprec_BFM = [1/(2 kappa)] g5theta(ctheta,stheta) Mprec_CPS
    //thus
    //Mprec_BFM^dag = [1/(2 kappa)]  Mprec_CPS^dag g5theta(ctheta,-stheta)
    //i.e. need to multiply input vector


    {
      GnoneFwilsonTm *latt_cps = new GnoneFwilsonTm;
      CgArg cg_arg;
      cg_arg.mass = mass;
      cg_arg.epsilon = epsilon;
      {
	Vector* tmp_chi_CPS = (Vector*)pmalloc(f_size/2 * sizeof(Float)); //odd checkerboard only

	DiracOpWilsonTm wilson(*latt_cps, (Vector*)NULL, (Vector*)NULL, &cg_arg, CNV_FRM_NO);
	wilson.MatPcDag(tmp_chi_CPS, (Vector*)v1);
	chi_CPS = (Vector*)convertCanonicalSingleChkb((Float*)tmp_chi_CPS, 1, *latt_cps);
	
	pfree(tmp_chi_CPS);
      }
      delete latt_cps;
    }

    {
      Vector* tmp_chi_BFM = (Vector*)pmalloc(f_size/2 * sizeof(Float));
      Vector* tmp_v1 = (Vector*)pmalloc(f_size/2 * sizeof(Float));
      for(int i=0;i<f_size/2;i++) ((Float*)tmp_v1)[i] = v1[i];
      g5theta(tmp_v1,GJP.VolNodeSites()/2,ctheta,stheta);

      GnoneFbfm *latt_fbfm = new GnoneFbfm;
      dynamic_cast<Fbfm*>(latt_fbfm)->MatPc(tmp_chi_BFM, tmp_v1, mass, epsilon, DAG_YES);

      delete latt_fbfm;
      GnoneFwilsonTm *latt_cps = new GnoneFwilsonTm;
      chi_BFM = (Vector*)convertCanonicalSingleChkb((Float*)tmp_chi_BFM, 1, *latt_cps);
      pfree(tmp_chi_BFM);
      delete latt_cps;
    }

    for(long i=0;i<f_size;i++){
      int rem = i;
      int m = rem % 24; rem/=24;
      int pos[4];
      for(int j=0;j<4;j++){ pos[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem/=GJP.NodeSites(j); }

      Float cps = ((Float*)chi_CPS)[i];
      Float bfm = ((Float*)chi_BFM)[i]*2*kappa;


      if( fabs(cps-bfm) > 1e-08 ){ 
	printf("**Mpc norm CPS/FBFM test dag version fail, %d %d %d %d, m=%d : %.9e %.9e, ratio %f\n",pos[0],pos[1],pos[2],pos[3],
	       m,cps,bfm,cps/bfm);
	fail=true; 
      }
    }
    if(fail){
      printf("Mpc norm CPS/FBFM test dag version failed\n"); exit(-1);
    }else printf("Mpc norm CPS/FBFM test dag version passed\n");







    //Given 2 random HMD force vectors, does fbfm::EvolveMomFforce and FwilsonTm give the same result? 
    //FwilsonTm has the calculation of the Hmd force vectors baked in, so I made a new version without that step

    //Note: WilsonTm calculation of HMD force vecs, v2 is the one that contains Deo^dag and Mprec chi
    //In BFM version, v1 is the one that contains Deo^dag and Mprec chi
    //hence we need to swap order of v1 and v2 in arguments to BFM version for comparison

    Float* mom_cps_tm = (Float *)pmalloc( sizeof(Float) * mom_size);
    for(int i=0;i<mom_size;i++) mom_cps_tm[i] = 0.0;

    {
      GnoneFwilsonTm *latt_cps = new GnoneFwilsonTm;
      CPS_EvolveMomFforce((Matrix*)mom_cps_tm, (Vector*)v1, (Vector*)v2, 
			  mass, epsilon, 0.1, dynamic_cast<FwilsonTm*>(latt_cps));
      delete latt_cps;
    }

    Float* mom_bfm_tm = (Float *)pmalloc( sizeof(Float) * mom_size);
    for(int i=0;i<mom_size;i++) mom_bfm_tm[i] = 0.0;

    {
      GnoneFbfm *latt_fbfm = new GnoneFbfm;
      //argument order reversed as above
      FforceWilsonType cal_force((Matrix*)mom_bfm_tm, latt_fbfm->GaugeField(),
				 v2, v1, Fbfm::bfm_args[0].Ls, 0.1);
      ForceArg ret = cal_force.run();

      delete latt_fbfm;
    }
    //The force in the CPS version appears to be exactly 2x larger than the BFM version, why is this? Its explicitly baked into the code, in all CPS versions of EvolveMomFforce (cf. f_wilson.C:568 or f_dwf_base_force.C:189)
    Float fact_cps = 0.5;

    for(int i=0;i<mom_size;i++){
      int rem = i;
      int midx = rem % 18; rem/=18;
      int mu = rem % 4; rem/=4;

      int x[4];
      for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

      if(fabs(fact_cps*mom_cps_tm[i] - mom_bfm_tm[i])>1e-08){
	printf("Fail EvolveMomFforce CPS/BFM code test fail midx %d mu %d (%d %d %d %d): %.9e %.9e, ratio %.9e\n",midx,mu,x[0],x[1],x[2],x[3],fact_cps*mom_cps_tm[i],mom_bfm_tm[i],fact_cps*mom_cps_tm[i]/mom_bfm_tm[i] );
	fail=true;
      }
    }
    if(fail){ printf("Failed EvolveMomFforce CPS/BFM code test\n"); exit(-1); }
    else printf("Passed EvolveMomFforce CPS/BFM code test\n");

    //OK, we know now how to make the forces agree and also the HMD force vectors. Can we combine this knowledge? A certain amount of the normalization differences should cancel out,
    //so we just try to find out what we have to do to the initial random gaussian vector to get the two to agree
    
    for(int i=0;i<mom_size;i++) mom_cps_tm[i] = 0.0;

    {
      GnoneFwilsonTm *latt_cps = new GnoneFwilsonTm;
      latt_cps->EvolveMomFforce((Matrix*)mom_cps_tm, (Vector*) v1, mass, epsilon, 0.1); 
      delete latt_cps;
    }

    for(int i=0;i<mom_size;i++) mom_bfm_tm[i] = 0.0;

    {
      Float *v1_temp = (Float*)pmalloc(f_size*sizeof(Float));
      for(int i=0;i<f_size;i++) v1_temp[i] = 2*kappa*v1[i];

      GnoneFbfm *latt_fbfm = new GnoneFbfm;
      dynamic_cast<Fbfm*>(latt_fbfm)->EvolveMomFforce((Matrix*)mom_bfm_tm, (Vector*) v1_temp, mass, epsilon, 0.1);
      delete latt_fbfm;
      pfree(v1_temp);
    }
    //OK, it seems that the differences in normalization mostly cancel apart from an overall factor of 4*kappa^2 (i.e. P_CPS = 4*kappa^2 P_bfm). 
    //We can fix by taking the phi vector used for the BFM part and multiplying by 2 kappa.  Fixed up the normalization above.

    fact_cps = 1.0;

    for(int i=0;i<mom_size;i++){
      int rem = i;
      int midx = rem % 18; rem/=18;
      int mu = rem % 4; rem/=4;

      int x[4];
      for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

      if(fabs(fact_cps*mom_cps_tm[i] - mom_bfm_tm[i])>1e-08){
	printf("Fail EvolveMomFforce CPS/BFM code test 2 fail midx %d mu %d (%d %d %d %d): %.9e %.9e, ratio %.9e\n",midx,mu,x[0],x[1],x[2],x[3],fact_cps*mom_cps_tm[i],mom_bfm_tm[i],fact_cps*mom_cps_tm[i]/mom_bfm_tm[i] );
	fail=true;
      }
    }
    if(fail){ printf("Failed EvolveMomFforce CPS/BFM code test 2\n"); exit(-1); }
    else printf("Passed EvolveMomFforce CPS/BFM code test 2\n");

    //The normalization difference should go away if we use the proper heatbath, which includes a MatPc^dag factor on the input. Let's test this
    
    for(int i=0;i<mom_size;i++) mom_cps_tm[i] = 0.0;

    {
      GnoneFwilsonTm *latt_cps = new GnoneFwilsonTm;
      Float *phi = (Float*)pmalloc(f_size/2*sizeof(Float));
      latt_cps->SetPhi((Vector*)phi,(Vector*)v1,(Vector*)NULL,mass,epsilon,DAG_YES);
      latt_cps->EvolveMomFforce((Matrix*)mom_cps_tm, (Vector*) phi, mass, epsilon, 0.1); 
      delete latt_cps;
      pfree(phi);
    }

    for(int i=0;i<mom_size;i++) mom_bfm_tm[i] = 0.0;

    {
      GnoneFbfm *latt_fbfm = new GnoneFbfm;
      Float *phi = (Float*)pmalloc(f_size/2*sizeof(Float));

      Float *v1_temp = (Float*)pmalloc(f_size*sizeof(Float));
      for(int i=0;i<f_size;i++) v1_temp[i] = 4*kappa*kappa*v1[i];
      g5theta((Vector*)v1_temp,GJP.VolNodeSites()/2,ctheta,stheta); 

      dynamic_cast<Fbfm*>(latt_fbfm)->SetPhi((Vector*)phi,(Vector*)v1_temp,(Vector*)NULL,mass,epsilon,DAG_YES);
      dynamic_cast<Fbfm*>(latt_fbfm)->EvolveMomFforce((Matrix*)mom_bfm_tm, (Vector*) phi, mass, epsilon, 0.1);
      delete latt_fbfm;
      pfree(phi);
      pfree(v1_temp);
    }
    //Hmm, that just makes it worse! Now there is not just an overall factor but something else.

    //So we know that EvolveMomFforce agree up to P_CPS = 4*kappa^2 P_bfm, which requires the 'phi' vector to be multiplied by 2 kappa before giving to the BFM EvolveMomFforce
    //We therefore need to make the heatbath vector,  phi_BFM = Mpc_BFM^dag v1   to be equal to 2*kappa* phi_CPS =  2*kappa Mpc_CPS^dag v1
    //We know Mpc_BFM^dag = [1/(2 kappa)]  Mpc_CPS^dag g5theta(ctheta,-stheta)
    //hence phi_BFM = [1/(2 kappa)] Mpc_CPS^dag  g5theta(ctheta,-stheta) v1.  We just need to multiply v1 by 4*kappa^2*g5theta(ctheta,stheta) then. I added this above.

    for(int i=0;i<mom_size;i++){
      int rem = i;
      int midx = rem % 18; rem/=18;
      int mu = rem % 4; rem/=4;

      int x[4];
      for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

      if(fabs(fact_cps*mom_cps_tm[i] - mom_bfm_tm[i])>1e-08){
	printf("Fail EvolveMomFforce CPS/BFM code test 3 fail midx %d mu %d (%d %d %d %d): %.9e %.9e, ratio %.9e\n",midx,mu,x[0],x[1],x[2],x[3],fact_cps*mom_cps_tm[i],mom_bfm_tm[i],fact_cps*mom_cps_tm[i]/mom_bfm_tm[i] );
	fail=true;
      }
    }
    if(fail){ printf("Failed EvolveMomFforce CPS/BFM code test 3\n"); exit(-1); }
    else printf("Passed EvolveMomFforce CPS/BFM code test 3\n"); 

    //So in conclusion, the random gaussian vector passed into the BFM evolution code needs to be 4*kappa^2 larger (g5theta is unitary, so it doesn't matter really)

    //Test FmatEvlInv.   out = (Mpc^dag Mpc)^{-1} in.    
    //BFM:  out = (Mpc_BFM^dag Mpc_BFM)^{-1} in =  (  [1/(4 kappa^2)]  Mpc_CPS^dag g5theta(ctheta,-stheta) g5theta(ctheta,stheta) Mpc_CPS  )^{-1} in
    //                                          = 4 kappa^2 (Mpc_CPS^dag Mpc_CPS  )^{-1} in
    Vector* test_inv_CPS;
    {
      GnoneFwilsonTm *latt_cps = new GnoneFwilsonTm;

      CgArg cg_arg;
      cg_arg.mass = mass;
      cg_arg.epsilon = epsilon;
      cg_arg.max_num_iter = 10000;
      cg_arg.stop_rsd = 1e-10;
      cg_arg.true_rsd = 1e-10;
      cg_arg.Inverter = CG;

      Vector* tmp_inv = (Vector*)pmalloc(f_size*sizeof(Float));
      tmp_inv -> VecZero(f_size);

      Float pooh;
      latt_cps->FmatEvlInv(tmp_inv, (Vector*)v1, &cg_arg, &pooh, CNV_FRM_NO);
     
      test_inv_CPS = (Vector*)convertCanonicalSingleChkb((Float*)tmp_inv, 1, *latt_cps);
      pfree(tmp_inv);
      delete latt_cps;
    }

    Vector* test_inv_BFM;
    {
      GnoneFbfm *latt_fbfm = new GnoneFbfm;
	
      CgArg cg_arg;
      cg_arg.mass = mass;
      cg_arg.epsilon = epsilon;
      cg_arg.max_num_iter = 10000;
      cg_arg.stop_rsd = 1e-10;
      cg_arg.true_rsd = 1e-10;
      cg_arg.Inverter = CG;
      
      Vector* tmp_inv = (Vector*)pmalloc(f_size*sizeof(Float));
      tmp_inv -> VecZero(f_size);
      
      Float pooh;
      latt_fbfm->FmatEvlInv(tmp_inv, (Vector*)v1, &cg_arg, &pooh, CNV_FRM_NO);
      
      delete latt_fbfm;

      GnoneFwilsonTm *latt_cps = new GnoneFwilsonTm;
      test_inv_BFM = (Vector*)convertCanonicalSingleChkb((Float*)tmp_inv, 1, *latt_cps);
      delete latt_cps;

      pfree(tmp_inv);
    }
    Float fac_bfm = 1.0/(4*kappa*kappa);

    for(long i=0;i<f_size;i++){
      int rem = i;
      int m = rem % 24; rem/=24;
      int pos[4];
      for(int j=0;j<4;j++){ pos[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem/=GJP.NodeSites(j); }

      Float cps = ((Float*)test_inv_CPS)[i];
      Float bfm = ((Float*)test_inv_BFM)[i]*fac_bfm;

      if( fabs(cps-bfm) > 1e-08 ){ 
	printf("**EvlInv CPS/FBFM test fail, %d %d %d %d, m=%d : %.9e %.9e, ratio %f\n",pos[0],pos[1],pos[2],pos[3],
	       m,cps,bfm,cps/bfm);
	fail=true; 
      }
    }
    if(fail){
      printf("EvlInv CPS/FBFM test failed\n"); exit(-1);
    }else printf("EvlInv CPS/FBFM test passed\n");

    //Show that by correctly normalizing the input vectors we can get the quotient to agree

    //Heatbath:    phi = Mpc_bose (Mpc_bose^dag Mpc_bose)^{-1} Mpc_ferm^dag v1 = (Mpc_bose^dag)^{-1} Mpc_ferm^dag v1
    //BFM:  phi_BFM = ([1/(2 kappa_b)]  Mpc_b_CPS^dag g5theta(ctheta_b,-stheta_b) )^{-1} [1/(2 kappa_f)]  Mpc_f_CPS^dag g5theta(ctheta_f,-stheta_f) v1
    //              = kappa_b/kappa_f g5theta(ctheta_b,stheta_b) (Mpc_b_CPS^dag)^{-1} Mpc_f_CPS^dag g5theta(ctheta_f,-stheta_f) v1

    Vector* phi_quo_CPS;
    //choose some random numbers for the boson mass and twist
    Float eps_b = 0.9;
    Float mass_b = -1.8;
    
    {
      Vector* tmp_phi = (Vector*)pmalloc(f_size/2*sizeof(Float));

      GnoneFwilsonTm *latt_cps = new GnoneFwilsonTm;
      latt_cps->SetPhi(tmp_phi, (Vector*)v1, (Vector*)NULL, mass, epsilon, DAG_YES); //tmp_phi = Mpc_ferm^dag v1
	
      CgArg cg_arg;
      cg_arg.mass = mass_b;
      cg_arg.epsilon = eps_b;
      cg_arg.max_num_iter = 10000;
      cg_arg.stop_rsd = 1e-10;
      cg_arg.true_rsd = 1e-10;
      cg_arg.Inverter = CG;

      Vector* tmp_inv = (Vector*)pmalloc(f_size*sizeof(Float));
      tmp_inv -> VecZero(f_size);
      Float pooh;
      latt_cps->FmatEvlInv(tmp_inv, tmp_phi, &cg_arg, &pooh, CNV_FRM_NO); //tmp_inv = (Mpc_bose^dag Mpc_bose)^{-1} tmp_phi
      
      latt_cps->SetPhi(tmp_phi, tmp_inv, (Vector*)NULL, mass_b, eps_b, DAG_NO); //tmp_phi = Mpc_bose tmp_inv
      phi_quo_CPS = (Vector*)convertCanonicalSingleChkb((Float*)tmp_phi, 1, *latt_cps);
      pfree(tmp_inv);
      pfree(tmp_phi);
      delete latt_cps;
    }

    Vector* phi_quo_BFM;
    {
      Vector* tmp_phi = (Vector*)pmalloc(f_size/2*sizeof(Float));

      Float kappa_ferm = 1.0/2.0/sqrt( (mass+4.0)*(mass+4.0) + epsilon*epsilon );
      Float kappa_bose = 1.0/2.0/sqrt( (mass_b+4.0)*(mass_b+4.0) + eps_b*eps_b );

      Float ctheta_f = 2*kappa_ferm*(mass+4.0);
      Float stheta_f = 2*kappa_ferm*epsilon;

      Float ctheta_b = 2*kappa_bose*(mass_b+4.0);
      Float stheta_b = 2*kappa_bose*eps_b;


      Float *v1_temp = (Float*)pmalloc(f_size*sizeof(Float));
      for(int i=0;i<f_size;i++) v1_temp[i] = kappa_ferm/kappa_bose*v1[i];
      g5theta((Vector*)v1_temp,GJP.VolNodeSites()/2,ctheta_f,stheta_f); //to cancel the extra factor of g5theta(ctheta_f,-stheta_f)
      

      GnoneFbfm *latt_fbfm = new GnoneFbfm;
      latt_fbfm->SetPhi(tmp_phi, (Vector*)v1_temp, (Vector*)NULL, mass, epsilon, DAG_YES);
	
      CgArg cg_arg;
      cg_arg.mass = mass_b;
      cg_arg.epsilon = eps_b;
      cg_arg.max_num_iter = 10000;
      cg_arg.stop_rsd = 1e-10;
      cg_arg.true_rsd = 1e-10;
      cg_arg.Inverter = CG;
      
      Vector* tmp_inv = (Vector*)pmalloc(f_size*sizeof(Float));
      tmp_inv -> VecZero(f_size);
      
      Float pooh;
      latt_fbfm->FmatEvlInv(tmp_inv, tmp_phi, &cg_arg, &pooh, CNV_FRM_NO);
      
      latt_fbfm->SetPhi(tmp_phi, tmp_inv, (Vector*)NULL, mass_b, eps_b, DAG_NO);

      g5theta(tmp_phi,GJP.VolNodeSites()/2,ctheta_b,-stheta_b);  //cancel extra factor of g5theta(ctheta_b,stheta_b)


      delete latt_fbfm;

      GnoneFwilsonTm *latt_cps = new GnoneFwilsonTm;
      phi_quo_BFM = (Vector*)convertCanonicalSingleChkb((Float*)tmp_phi, 1, *latt_cps);
      delete latt_cps;

      pfree(tmp_inv);
      pfree(tmp_phi);
    }

    for(long i=0;i<f_size;i++){
      int rem = i;
      int m = rem % 24; rem/=24;
      int pos[4];
      for(int j=0;j<4;j++){ pos[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem/=GJP.NodeSites(j); }

      Float cps = ((Float*)phi_quo_CPS)[i];
      Float bfm = ((Float*)phi_quo_BFM)[i];

      if( fabs(cps-bfm) > 1e-08 ){ 
	printf("**Quo phi CPS/FBFM test fail, %d %d %d %d, m=%d : %.9e %.9e, ratio %f\n",pos[0],pos[1],pos[2],pos[3],
	       m,cps,bfm,cps/bfm);
	fail=true; 
      }
    }
    if(fail){
      printf("Quo phi CPS/FBFM test failed\n"); exit(-1);
    }else printf("Quo phi CPS/FBFM test passed\n");


    //OK, lets try putting it together

    {
      ActionQuotientArg quo_arg;
      setupQuoArg(quo_arg, 1, &mass_b, &eps_b, &mass, &epsilon, F_CLASS_BFM);
      
      Float* mom_test3 = (Float *)pmalloc( sizeof(Float) * mom_size);
      Float* mom_test4 = (Float *)pmalloc( sizeof(Float) * mom_size);

      AlgMomentum mom;
      Float* m1 = (Float*)mom.getMom();

      printf("Doing Quotient BFM version\n");
      LatRanGen LRGbak = LRG;
      {
	for(int i=0;i<mom_size;i++) m1[i] = 0.0;
    
	AlgActionQuotient quo(mom,quo_arg);
	LRG = LRGbak;
	LRG.AssignGenerator(0);
	Float grand_5d = LRG.Grand(FIVE_D);
	Float grand_4d = LRG.Grand(FOUR_D);
	printf("Random numbers: %f %f\n",grand_5d,grand_4d);

	LatRanGen LRGbak2 = LRG;
	quo.heatbath(); //just use heatbath method to setup other variables, we manually fixup the input RGVs to match the CPS version
	LRG = LRGbak2;

	{
	  Vector **phi = quo.getPhi();
	  
	  Lattice &lat = LatticeFactory::Create(F_CLASS_BFM, G_CLASS_NONE);
	  
	  Vector *tmp1 = (Vector*)pmalloc(f_size/2*sizeof(Float));
	  Vector *tmp2 = (Vector*)pmalloc(f_size*sizeof(Float));

	  lat.RandGaussVector(tmp1, 0.5, 1);

	  //Twist RGV by g5theta
	  Float kappa_ferm = 1.0/2.0/sqrt( (mass+4.0)*(mass+4.0) + epsilon*epsilon );
	  Float kappa_bose = 1.0/2.0/sqrt( (mass_b+4.0)*(mass_b+4.0) + eps_b*eps_b );

	  Float ctheta_f = 2*kappa_ferm*(mass+4.0);
	  Float stheta_f = 2*kappa_ferm*epsilon;

	  Float ctheta_b = 2*kappa_bose*(mass_b+4.0);
	  Float stheta_b = 2*kappa_bose*eps_b;
	  g5theta(tmp1,GJP.VolNodeSites()/2,ctheta_f,stheta_f); 
	  //proceed as usual

	  lat.SetPhi(phi[0], tmp1, (Vector*)NULL, mass, epsilon, DAG_YES);

	  CgArg cg_arg;
	  cg_arg.mass = mass_b;
	  cg_arg.epsilon = eps_b;
	  cg_arg.max_num_iter = 10000;
	  cg_arg.stop_rsd = 1e-10;
	  cg_arg.true_rsd = 1e-10;
	  cg_arg.Inverter = CG;

	  tmp2->VecZero(f_size);
	  Float pooh;
	  lat.FmatEvlInv(tmp2, phi[0], &cg_arg, &pooh, CNV_FRM_NO);

	  lat.SetPhi(phi[0], tmp2, (Vector*)NULL, mass_b, eps_b, DAG_NO);
	  LatticeFactory::Destroy();

	  pfree(tmp1); pfree(tmp2);
        }
	
	printf("Quo bfm initial energy: %f\n",quo.energy());

	quo.evolve(0.1, 1);
    
	for(int i=0;i<mom_size;i++) mom_test3[i] = m1[i];
      }
      LRG = LRGbak;

      ActionQuotientArg quo_arg_2;
      setupQuoArg(quo_arg_2, 1, &mass_b, &eps_b, &mass, &epsilon, F_CLASS_WILSON_TM);

      printf("Doing Quotient CPS version\n");
      {
	for(int i=0;i<mom_size;i++) m1[i] = 0.0;
    
	AlgActionQuotient quo(mom,quo_arg_2);
	LRG = LRGbak;
	LRG.AssignGenerator(0);
	Float grand_5d = LRG.Grand(FIVE_D);
	Float grand_4d = LRG.Grand(FOUR_D);
	printf("Random numbers: %f %f\n",grand_5d,grand_4d);

	quo.heatbath();
	Vector **phi = quo.getPhi();
	printf("Quo cps initial energy: %f\n",quo.energy());

	quo.evolve(0.1, 1);
    
	for(int i=0;i<mom_size;i++) mom_test4[i] = m1[i];
      }
      LRG = LRGbak;

      //lattice = new GnoneFwilsonTm;

      fail = false;
      for(int i=0;i<mom_size;i++){
	int rem = i;
	int midx = rem % 18; rem/=18;
	int mu = rem % 4; rem/=4;

	int x[4];
	for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

	if(fabs(mom_test4[i] - mom_test3[i])>1e-08){
	  printf("**Fbfm Quotient evol test fail midx %d mu %d (%d %d %d %d): %f %f, ratio %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test4[i],mom_test3[i],mom_test4[i]/mom_test3[i]);
	  fail=true;
	}//else printf("Pass Quotient evol test midx %d mu %d (%d %d %d %d): %f %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test2[i],mom_test1[i]);
      }
      if(fail){ printf("Failed Quotient evol test\n"); exit(-1); }
      else printf("Passed Quotient evol test\n");
    }


// 	  if(0){
// 	    //as we are multiplying by the square root of a matrix, we need to multiply by the square root of g5theta as well
// 	    //g5 is in chiral basis, so matrix is diagonal - this is thus quite easy. Represent  theta_c + i theta_s   as  amplitude * e^{i angle}
// 	    Float angle = atan2(stheta_f,ctheta_f); Float amplitude = stheta_f / sin(angle);
// 	    Float root[2] = { sqrt(amplitude)*cos( angle/2.0 ),  sqrt(amplitude)*sin( angle/2.0 ) };

// #define PSIA(f,r,c,s,n)     (f)+(r+2*(c+3*(s+4*(n))))

// 	    IFloat* phiif = reinterpret_cast<IFloat*>(phi[i]);
// 	    for(int n=0;n<GJP.VolNodeSites()/2;n++){
// 	      for(int c=0;c<3;c++){
// 		cTimesC_2(PSIA(phiif,0,c,0,n), root[0], root[1]);
// 		cTimesC_2(PSIA(phiif,0,c,1,n), root[0], root[1]);
// 		cTimesC_2(PSIA(phiif,0,c,2,n), root[0], -root[1]);
// 		cTimesC_2(PSIA(phiif,0,c,3,n), root[0], -root[1]);
// 	      }	
// 	    }
// 	  }


    //OK, now try it for the rational quotient. First try to get heatbaths to agree
    //Heatbath:    phi = (Mpc_bose^dag Mpc_bose)^{- pow} (Mpc_ferm^dag Mpc_ferm)^pow v1
    //BFM:  phi_BFM = ([1/(4 kappa_b^2)]  Mpc_b_CPS^dag Mpc_b_CPS  )^{-pow}  ( [1/(4 kappa_f^2)]  Mpc_f_CPS^dag Mpc_f_CPS)^pow  v1
    //              = (kappa_b^2/kappa_f^2)^pow (Mpc_b_CPS^dag Mpc_b_CPS)^{-pow} (Mpc_f_CPS^dag Mpc_f_CPS)^pow v1

    //Note, pow = 1/2 * pow_in,   where pow_in is the power specified in the input arguments; cf. AlgActionRational::generateApprox line 542

    //After testing, I have determined that the normalization is correct as it is. However the two versions only agree once you tighten the precision of the
    //rational approximation. This is because there is a different normalization: the precision of the approximation will be different. To make them the same
    //we would need to modify the eigenvalue bounds of the approximation such that they match after normalization

    //OK: after fixing the eigenvalue bounds we get excellent agreement: cf. normalization in setupRatQuoArg

    int n_masses = 1;
    ActionRationalQuotientArg rat_quo_arg;
    Float bsn_masses[] = {-1.8};
    Float bsn_epsilon[] = {0.9};
    Float frm_masses[] = {0.5};
    Float frm_epsilon[] = {0.5};
    int pwr_num[] = {1};
    int pwr_den[] = {2};

    setupRatQuoArg(rat_quo_arg, 1, &bsn_masses[0], &bsn_epsilon[0], &frm_masses[0], &frm_epsilon[0], &pwr_num[0], &pwr_den[0], F_CLASS_BFM);

    Float* mom_test3 = (Float *)pmalloc( sizeof(Float) * mom_size);
    Float* mom_test4 = (Float *)pmalloc( sizeof(Float) * mom_size);

    Vector* cps_phi[n_masses];
    Vector* bfm_phi[n_masses];

    AlgMomentum mom;
    Float* m1 = (Float*)mom.getMom();

    printf("Doing BFM version\n");
    LatRanGen LRGbak = LRG;
    {
      for(int i=0;i<mom_size;i++) m1[i] = 0.0;
    
      AlgActionRationalQuotient rat_quo(mom,rat_quo_arg,0);
      LRG = LRGbak;
      LRG.AssignGenerator(0);
      Float grand_5d = LRG.Grand(FIVE_D);
      Float grand_4d = LRG.Grand(FOUR_D);
      printf("Random numbers: %f %f\n",grand_5d,grand_4d);

      LatRanGen LRGbak2 = LRG;
      rat_quo.heatbath();
      LRG = LRGbak2;

      Vector **phi = rat_quo.getPhi();

      {
	GnoneFwilsonTm *lat = new GnoneFwilsonTm;
	for(int i=0;i<n_masses;i++){
	  Float kappa_ferm = 1.0/2.0/sqrt( (frm_masses[i]+4.0)*(frm_masses[i]+4.0) + frm_epsilon[i]*frm_epsilon[i] );
	  Float kappa_bose = 1.0/2.0/sqrt( (bsn_masses[i]+4.0)*(bsn_masses[i]+4.0) + bsn_epsilon[i]*bsn_epsilon[i] );

	  Float ctheta_f = 2*kappa_ferm*(frm_masses[i]+4.0);
	  Float stheta_f = 2*kappa_ferm*frm_epsilon[i];

	  Float ctheta_b = 2*kappa_bose*(bsn_masses[i]+4.0);
	  Float stheta_b = 2*kappa_bose*bsn_epsilon[i];

	  bfm_phi[i] = (Vector*)convertCanonicalSingleChkb((Float*)phi[i], 1, *lat);
	  for(int j=0;j<f_size;j++) ((Float*)bfm_phi[i])[j] = pow(kappa_ferm/kappa_bose, 2*Float(pwr_num[i])/Float(2*pwr_den[i]) )* ((Float*)bfm_phi[i])[j];
	}
	delete lat;
      }

      printf("RatQuo bfm initial energy: %f\n",rat_quo.energy());

      rat_quo.evolve(0.1, 1);
    
      for(int i=0;i<mom_size;i++) mom_test3[i] = m1[i];
    }     
    LRG = LRGbak;

    ActionRationalQuotientArg rat_quo_arg_2;
    setupRatQuoArg(rat_quo_arg_2, 1, &bsn_masses[0], &bsn_epsilon[0], &frm_masses[0], &frm_epsilon[0], &pwr_num[0], &pwr_den[0], F_CLASS_WILSON_TM);

    printf("Doing CPS version\n");
    {
      for(int i=0;i<mom_size;i++) m1[i] = 0.0;
    
      AlgActionRationalQuotient rat_quo(mom,rat_quo_arg_2,0);
      LRG = LRGbak;
      LRG.AssignGenerator(0);
      Float grand_5d = LRG.Grand(FIVE_D);
      Float grand_4d = LRG.Grand(FOUR_D);
      printf("Random numbers: %f %f\n",grand_5d,grand_4d);

      rat_quo.heatbath();
      Vector **phi = rat_quo.getPhi();
      
      {
	GnoneFwilsonTm *lat = new GnoneFwilsonTm;
	for(int i=0;i<n_masses;i++) cps_phi[i] = (Vector*)convertCanonicalSingleChkb((Float*)phi[i], 1, *lat);
	delete lat;
      }

      printf("RatQuo cps initial energy: %f\n",rat_quo.energy());

      rat_quo.evolve(0.1, 1);
    
      for(int i=0;i<mom_size;i++) mom_test4[i] = m1[i];
    }
    LRG = LRGbak;

    //Compare phi

    lattice = new GnoneFwilsonTm;

    for(int midx=0;midx<n_masses;midx++){
      
      for(long i=0;i<f_size;i++){
	int rem = i;
	int m = rem % 24; rem/=24;
	int pos[4];
	for(int j=0;j<4;j++){ pos[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem/=GJP.NodeSites(j); }
	
	Float cps = ((Float*)(cps_phi[midx]))[i];
	Float bfm = ((Float*)(bfm_phi[midx]))[i];
	
	if( fabs(cps-bfm) > 1e-08 ){ 
	  printf("**RatQuo phi CPS/FBFM test fail, %d %d %d %d, m=%d : %.9e %.9e, ratio %.9f\n",pos[0],pos[1],pos[2],pos[3],
		 m,cps,bfm,cps/bfm);
	  fail=true; 
	}
      }
      if(fail){
	printf("RatQuo phi CPS/FBFM test failed\n"); //exit(-1);
      }else printf("RatQuo phi CPS/FBFM test passed\n");
    }

    fail = false;
    for(int i=0;i<mom_size;i++){
      int rem = i;
      int midx = rem % 18; rem/=18;
      int mu = rem % 4; rem/=4;

      int x[4];
      for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

      if(fabs(mom_test4[i] - mom_test3[i])>1e-08){
	printf("**Fbfm twisted mass test fail midx %d mu %d (%d %d %d %d): %.9f %.9f, ratio %.9f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test4[i],mom_test3[i],mom_test4[i]/mom_test3[i]);
	fail=true;
      }else printf("Pass twisted mass test midx %d mu %d (%d %d %d %d): %.9f %.9f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test2[i],mom_test1[i]);
    }
    if(fail){ printf("Failed Fbfm twisted mass test\n"); exit(-1); }
    else printf("Passed Fbfm twisted mass test\n");
  }


  printf("All tests passed\n");
  return 0;
}









#if 1
int main(int argc,char *argv[])
{
  Start(&argc,&argv); //initialises QMP

  CommandLine::is(argc,argv);

  bool gparity_X(false);
  bool gparity_Y(false);

  int arg0 = CommandLine::arg_as_int(0);
  printf("Arg0 is %d\n",arg0);
  if(arg0==0){
    gparity_X=true;
    printf("Doing G-parity HMC test in X direction\n");
  }else if(arg0==1){
    printf("Doing G-parity HMC test in X and Y directions\n");
    gparity_X = true;
    gparity_Y = true;
  }else{
    printf("Doing No-G-parity twisted mass test\n");
  }

  bool save_config(false);
  bool load_config(false);
  bool load_lrg(false);
  bool save_lrg(false);
  char *load_config_file;
  char *save_config_file;
  char *save_lrg_file;
  char *load_lrg_file;
  bool gauge_fix(false);
  bool verbose(false);
  bool unit_gauge(false);

  int size[] = {2,2,2,2,2};

  int i=2;
  while(i<argc){
    char* cmd = argv[i];  
    if( strncmp(cmd,"-save_config",15) == 0){
      if(i==argc-1){ printf("-save_config requires an argument\n"); exit(-1); }
      save_config=true;
      save_config_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-load_config",15) == 0){
      if(i==argc-1){ printf("-save_config requires an argument\n"); exit(-1); }
      load_config=true;
      load_config_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-latt",10) == 0){
      if(i>argc-6){
	printf("Did not specify enough arguments for 'latt' (require 5 dimensions)\n"); exit(-1);
      }
      size[0] = CommandLine::arg_as_int(i); //CommandLine ignores zeroth input arg (i.e. executable name)
      size[1] = CommandLine::arg_as_int(i+1);
      size[2] = CommandLine::arg_as_int(i+2);
      size[3] = CommandLine::arg_as_int(i+3);
      size[4] = CommandLine::arg_as_int(i+4);
      i+=6;
    }else if( strncmp(cmd,"-load_lrg",15) == 0){
      if(i==argc-1){ printf("-load_lrg requires an argument\n"); exit(-1); }
      load_lrg=true;
      load_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-save_lrg",15) == 0){
      if(i==argc-1){ printf("-save_lrg requires an argument\n"); exit(-1); }
      save_lrg=true;
      save_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-gauge_fix",15) == 0){
      gauge_fix=true;
      i++;   
    }else if( strncmp(cmd,"-verbose",15) == 0){
      verbose=true;
      i++;
    }else if( strncmp(cmd,"-unit_gauge",15) == 0){
      unit_gauge=true;
      i++;
    }else{
      if(UniqueID()==0) printf("Unrecognised argument: %s\n",cmd);
      exit(-1);
    }
  }
  

  printf("Lattice size is %d %d %d %d\n",size[0],size[1],size[2],size[3],size[4]);

  DoArg do_arg;
  do_arg.x_sites = size[0];
  do_arg.y_sites = size[1];
  do_arg.z_sites = size[2];
  do_arg.t_sites = size[3];
  do_arg.s_sites = size[4];
  do_arg.x_node_sites = 0;
  do_arg.y_node_sites = 0;
  do_arg.z_node_sites = 0;
  do_arg.t_node_sites = 0;
  do_arg.s_node_sites = 0;
  do_arg.x_nodes = 0;
  do_arg.y_nodes = 0;
  do_arg.z_nodes = 0;
  do_arg.t_nodes = 0;
  do_arg.s_nodes = 0;
  do_arg.updates = 0;
  do_arg.measurements = 0;
  do_arg.measurefreq = 0;
  do_arg.cg_reprod_freq = 10;
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_conf_load_addr = 0x0;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.start_seed_filename = "../rngs/ckpoint_rng.0";
  do_arg.start_conf_filename = "../configurations/ckpoint_lat.0";
  do_arg.start_conf_alloc_flag = 6;
  do_arg.wfm_alloc_flag = 2;
  do_arg.wfm_send_alloc_flag = 2;
  do_arg.start_seed_value = 83209;
  do_arg.beta =   2.25;
  do_arg.c_1 =   -3.3100000000000002e-01;
  do_arg.u0 =   1.0000000000000000e+00;
  do_arg.dwf_height =   1.8000000000000000e+00;
  do_arg.dwf_a5_inv =   1.0000000000000000e+00;
  do_arg.power_plaq_cutoff =   0.0000000000000000e+00;
  do_arg.power_plaq_exponent = 0;
  do_arg.power_rect_cutoff =   0.0000000000000000e+00;
  do_arg.power_rect_exponent = 0;
  do_arg.verbose_level = -1202; //VERBOSE_DEBUG_LEVEL; //-1202;
  do_arg.checksum_level = 0;
  do_arg.exec_task_list = 0;
  do_arg.xi_bare =   1.0000000000000000e+00;
  do_arg.xi_dir = 3;
  do_arg.xi_v =   1.0000000000000000e+00;
  do_arg.xi_v_xi =   1.0000000000000000e+00;
  do_arg.clover_coeff =   0.0000000000000000e+00;
  do_arg.clover_coeff_xi =   0.0000000000000000e+00;
  do_arg.xi_gfix =   1.0000000000000000e+00;
  do_arg.gfix_chkb = 1;
  do_arg.asqtad_KS =   0.0000000000000000e+00;
  do_arg.asqtad_naik =   0.0000000000000000e+00;
  do_arg.asqtad_3staple =   0.0000000000000000e+00;
  do_arg.asqtad_5staple =   0.0000000000000000e+00;
  do_arg.asqtad_7staple =   0.0000000000000000e+00;
  do_arg.asqtad_lepage =   0.0000000000000000e+00;
  do_arg.p4_KS =   0.0000000000000000e+00;
  do_arg.p4_knight =   0.0000000000000000e+00;
  do_arg.p4_3staple =   0.0000000000000000e+00;
  do_arg.p4_5staple =   0.0000000000000000e+00;
  do_arg.p4_7staple =   0.0000000000000000e+00;
  do_arg.p4_lepage =   0.0000000000000000e+00;

  if(verbose) do_arg.verbose_level = VERBOSE_DEBUG_LEVEL;

  if(gparity_X) do_arg.x_bc = BND_CND_GPARITY;
  if(gparity_Y) do_arg.y_bc = BND_CND_GPARITY;

  GJP.Initialize(do_arg);

  LRG.Initialize(); //usually initialised when lattice generated, but I pre-init here so I can load the state from file

#ifdef HAVE_BFM
  cps_qdp_init(&argc,&argv);
  Chroma::initialize(&argc,&argv);
#endif
  
  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }
  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }
  
  GnoneFwilsonTm* lattice = new GnoneFwilsonTm;
					       
  if(!load_config){
    printf("Creating gauge field\n");
    if(!unit_gauge) lattice->SetGfieldDisOrd();
    else lattice->SetGfieldOrd();
  }else{
    ReadLatticeParallel readLat;
    if(UniqueID()==0) printf("Reading: %s (NERSC-format)\n",load_config_file);
    readLat.read(*lattice,load_config_file);
    if(UniqueID()==0) printf("Config read.\n");
  }

  if(save_config){
    if(UniqueID()==0) printf("Saving config to %s\n",save_config_file);

    QioArg wt_arg(save_config_file,0.001);
    
    wt_arg.ConcurIONumber=32;
    WriteLatticeParallel wl;
    wl.setHeader("disord_id","disord_label",0);
    wl.write(*lattice,wt_arg);
    
    if(!wl.good()) ERR.General("main","()","Failed write lattice %s",save_config_file);

    if(UniqueID()==0) printf("Config written.\n");
  }

  if(gauge_fix){
    lattice->FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
    lattice->FixGauge(1e-06,2000);
    if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
  }
  cps_qdp_init(&argc,&argv);

  if(!gparity_X && !gparity_Y) return no_gparity_test(lattice);
  
  LatRanGen LRGbak(LRG);
  Float* v1 = rand_4d_canonical_fermion(*lattice);
  Float* v2 = rand_4d_canonical_fermion(*lattice);
  LRG = LRGbak;

  { //do a simple test of the CANONICAL->WILSON and WILSON->CANONICAL reordering code
    printf("Starting conversion test\n"); fflush(stdout);

    long f_size = (long)24 * GJP.VolNodeSites()*2;  // m + 24*( x + lx*(y+ly*(z+lz*(t+lt*flav))))
    Float *v1_wilson_and_back = (Float *)pmalloc(sizeof(Float) * f_size);
    moveFloat(v1_wilson_and_back, v1, f_size);
    lattice->Fconvert((Vector*)v1_wilson_and_back, WILSON, CANONICAL);
    lattice->Fconvert((Vector*)v1_wilson_and_back, CANONICAL, WILSON);
    bool fail = false;
    for(long i=0;i<f_size;i++){
      int rem = i;
      int m = rem % 24; rem/=24;
      int pos[4];
      for(int j=0;j<4;j++){ pos[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem/=GJP.NodeSites(j); }
      int flav = rem;
      
      if( fabs(v1_wilson_and_back[i] - v1[i]) > 1e-12 ){ printf("**Conversion test fail, %d %d %d %d, flav %d, m=%d : %f %f\n",pos[0],pos[1],pos[2],pos[3],flav,m,v1_wilson_and_back[i],v1[i]); fail=true; }
      //else printf("Conversion test pass, %d %d %d %d, flav %d, m=%d : %f %f\n",pos[0],pos[1],pos[2],pos[3],flav,m,v1_wilson_and_back[i],v1[i]);
    }
    pfree(v1_wilson_and_back);
    if(fail){
      printf("Conversion test failed\n");
      exit(-1);
    }else{ printf("Conversion test passed\n"); fflush(stdout); }
  }

  int mom_size = 18*4*GJP.VolNodeSites()*2;

  Float* mom_test1 = (Float *)pmalloc( sizeof(Float) * mom_size);
  for(int i=0;i<mom_size;i++) mom_test1[i] = 0.0;

  Float* mom_test2 = (Float *)pmalloc( sizeof(Float) * mom_size);
  for(int i=0;i<mom_size;i++) mom_test2[i] = 0.0;

  //Test both fermion and boson forces
  lattice->Fconvert((Vector*)v1, WILSON, CANONICAL);
  lattice->Fconvert((Vector*)v2, WILSON, CANONICAL);
  lattice->EvolveMomFforce((Matrix*)mom_test1, (Vector*)v1, 0.5, 0.3, 0.1); 
  lattice->EvolveMomFforce((Matrix*)mom_test2, (Vector*)v1, (Vector*)v2, 0.5, 0.3, 0.1); 
  lattice->Fconvert((Vector*)v1, CANONICAL, WILSON);
  lattice->Fconvert((Vector*)v2, CANONICAL, WILSON);


  Float *invcg_test_2fout;
  {
    //Test InvCG BFM version vs regular CPS version with G-parity
    //Use v1 in Wilson form
    lattice->Fconvert((Vector*)v1, WILSON, CANONICAL);
  
    CgArg cg_arg;
    cg_arg.mass = -1.8;
    cg_arg.epsilon = 0.5;
    cg_arg.max_num_iter = 10000;
    cg_arg.stop_rsd = 1e-10;
    cg_arg.true_rsd = 1e-10;
    cg_arg.Inverter = CG;
  
    DiracOpWilsonTm dop(*lattice, (Vector*)0, (Vector*)0, &cg_arg, CNV_FRM_NO);

    size_t f_size_cb =  GJP.VolNodeSites() * 24/2 * 2; //extra factor of 2 for G-parity second flavour
    invcg_test_2fout = (Float*)pmalloc(f_size_cb*sizeof(Float));
    
    Float* out_1 = invcg_test_2fout;
    Float* out_2 = (Float*)pmalloc(f_size_cb*sizeof(Float));

    Float true_rsd;
    InvCg_CPS((Vector*)out_2, (Vector*)v1, 0.0, &true_rsd, dop, *lattice, cg_arg);  //regular CPS version
    dop.InvCg((Vector*)out_1, (Vector*)v1, 0.0, &true_rsd); //bfm version
    
    bool fail(false);
    for(int i=0;i<f_size_cb;i++){
      if( fabs(out_1[i]-out_2[i])>1e-08 ){
	printf("InvCGtest fail %d: %f %f, ratio %f\n",i,out_1[i],out_2[i],out_1[i]/out_2[i]);
	fail=true;
      }
    }
    if(fail){
      printf("Failed InvCg test\n"); exit(-1);
    }else printf("Passed InvCg test\n");
    

    pfree(out_2);
    lattice->Fconvert((Vector*)v1, CANONICAL, WILSON);
  }

  Vector** minvcg_test_2fout;

  int nmass =3;
  Float mass[3] = {0.1, 0.2, 0.3};
  {
    //Test MInvCG version vs regular CPS version with G-parity
    lattice->Fconvert((Vector*)v1, WILSON, CANONICAL);
    Float targ_resid = 1e-12;

    CgArg cg_arg;
    cg_arg.mass = -1.8;
    cg_arg.epsilon = 0.5;
    cg_arg.max_num_iter = 10000;
    cg_arg.stop_rsd = targ_resid;
    cg_arg.true_rsd = targ_resid;
    cg_arg.Inverter = CG;
  
    DiracOpWilsonTm dop(*lattice, (Vector*)0, (Vector*)0, &cg_arg, CNV_FRM_NO);
    
    size_t f_size_cb =  GJP.VolNodeSites() * 24/2 * 2;
    Float src_norm_sq = ((Vector*)v1)->NormSqNode(f_size_cb);
    glb_sum(&src_norm_sq);
  
    Vector **out_1 = (Vector**)pmalloc(3 * sizeof(Vector*));
    Vector **out_2 = (Vector**)pmalloc(3 * sizeof(Vector*));
  
    //Float mass[3] = {-1.8, -1.0, -0.4};
    for(int j=0;j<3;j++){
      out_1[j] = (Vector*)pmalloc(f_size_cb*sizeof(Float));
      out_2[j] = (Vector*)pmalloc(f_size_cb*sizeof(Float));
      
      for(int i=0;i<f_size_cb;i++){
	((Float*)out_1[j])[i] = 0.0; ((Float*)out_2[j])[i] = 0.0;
      }
    }
    
    Float true_rsd[3] = {targ_resid, targ_resid, targ_resid};
    
    MInvCG_CPS((Vector**)out_2, (Vector*)v1, src_norm_sq, &mass[0], nmass, 0, &true_rsd[0], MULTI, NULL, dop, *lattice, cg_arg); //cps version
    
    for(int i=0;i<nmass;i++) true_rsd[i] = targ_resid;
    
    dop.MInvCG((Vector**)out_1, (Vector*)v1, src_norm_sq, &mass[0], nmass, 0, &true_rsd[0], MULTI, NULL); //bfm version
  
    minvcg_test_2fout = (Vector**)out_1;

    bool fail=false;
    for(int j=0;j<nmass;j++){
      Float* o1f = (Float*)out_1[j];
      Float* o2f = (Float*)out_2[j];

      for(int i=0;i<f_size_cb;i++){
	if( fabs(o1f[i]-o2f[i])>1e-08 ){
	  printf("MInvCGtest %d fail %d: %f %f, ratio %f\n",j,i,o1f[i],o2f[i],o1f[i]/o2f[i]);
	  fail=true;
	}
      }

      if(!fail) printf("Passed MInvCg test for shift %d\n",j);
    }
    if(fail){
      printf("Failed MInvCg test\n"); exit(-1);
    }else printf("Passed MInvCg test\n");
    
    for(int j=0;j<3;j++){
      //pfree(out_1[j]);
      pfree(out_2[j]);
    }
    //pfree(out_1);
    pfree(out_2);
    lattice->Fconvert((Vector*)v1, CANONICAL, WILSON);
  }

  Vector* v1_wilsonord;
  {
    long f_size = (long)24 * GJP.VolNodeSites()*2; 
    v1_wilsonord = (Vector*)pmalloc(f_size*sizeof(Float));
    for(int i=0;i<f_size;i++){
      ((Float*)v1_wilsonord)[i] = ((Float*)v1)[i];
    }
    lattice->Fconvert((Vector*)v1_wilsonord, WILSON, CANONICAL);
  }
    
  //Calculate eigenvalue of MdagM and compare with 1f version
  Float evalue_2f;
  Float *eigenvector_2f; 
  {
    CommonArg c_arg;
    c_arg.set_filename("pooh");
    c_arg.set_label("pooh");
   
    Float mass = -1.8;

    EigArg e_arg;
    generateEigArg(e_arg, mass, 0.5); //operator is MatPcDagMatPc
    e_arg.Rsdlam = 2e-6;
    e_arg.RsdR_a = e_arg.Rsdlam ;
    e_arg.RsdR_r = e_arg.Rsdlam ;
    e_arg.N_eig = 1;

    CgArg cg_arg;
    cg_arg.mass = mass;
    cg_arg.epsilon = 0.5;

    Float **evalues = (Float**) pmalloc(sizeof(Float*)); //indexed by eigenvalue index, then mass
    for(int i=0;i<1;i++) evalues[i] = (Float*)pmalloc(sizeof(Float));

    LatRanGen LRGbak = LRG;
    AlgEig algeig(*lattice,&c_arg, &e_arg);
    algeig.run(evalues);    
    LRG = LRGbak;

    evalue_2f = evalues[0][0];
    
    size_t f_size_cb = GJP.VolNodeSites()*24/2 * 2;
    eigenvector_2f = (Float*)pmalloc(f_size_cb * sizeof(Float) );

    Float* ev = (Float*)(algeig.getEigenVectors()[0]);
    for(int i=0;i<f_size_cb;i++) eigenvector_2f[i] = ev[i];
  }

  if(UniqueID()==0){ printf("Starting double lattice section\n"); fflush(stdout); }
  
  //Backup lattice
  int array_size = 2*lattice->GsiteSize() * GJP.VolNodeSites() * sizeof(Float);
  Matrix *orig_lattice = (Matrix *) pmalloc(array_size);
  memcpy((void*)orig_lattice, (void*)lattice->GaugeField(), array_size);

  lattice->FreeGauge(); //free memory and reset
  delete lattice; //lattice objects are singleton (scope_lock)
  
  //setup 1f model. Upon calling GJP.Initialize the lattice size will be doubled in the appropriate directions
  //and the boundary condition set to APRD
  if(gparity_X) do_arg.gparity_1f_X = 1;
  if(gparity_Y) do_arg.gparity_1f_Y = 1;

  GJP.Initialize(do_arg);

  if(GJP.Gparity()){ printf("Que?\n"); exit(-1); }
  if(UniqueID()==0) printf("Doubled lattice : %d %d %d %d\n", GJP.XnodeSites()*GJP.Xnodes(),GJP.YnodeSites()*GJP.Ynodes(),
			   GJP.ZnodeSites()*GJP.Znodes(),GJP.TnodeSites()*GJP.Tnodes());
  
  GnoneFwilsonTm doubled_lattice;
  setup_double_latt(doubled_lattice,orig_lattice,gparity_X,gparity_Y);
  setup_double_rng(gparity_X,gparity_Y);
 
  if(gauge_fix){
    doubled_lattice.FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
    doubled_lattice.FixGauge(1e-06,2000);
    if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
  }
 
#ifdef HAVE_BFM
  {
    QDP::multi1d<int> nrow(Nd);  
    for(int i = 0;i<Nd;i++) nrow[i] = GJP.Sites(i);
    QDP::Layout::setLattSize(nrow);
    QDP::Layout::create();
  }
#endif

  long f_size = (long)24 * GJP.VolNodeSites();

  //convert the random CANONICAL vectors from the 2f setup to 1f form
  Float *v1_dbl = (Float *)pmalloc(sizeof(Float) * f_size);
  setup_double_4d_vector((Vector*)v1_dbl,(Vector*)v1, gparity_X,gparity_Y);
  pfree(v1);

  Float *v2_dbl = (Float *)pmalloc(sizeof(Float) * f_size);
  setup_double_4d_vector((Vector*)v2_dbl,(Vector*)v2, gparity_X,gparity_Y);
  pfree(v2);

  {
    //Test conversion
    LatRanGen LRGbak(LRG);
    Float* v1_test = rand_4d_canonical_fermion(doubled_lattice);
    LRG = LRGbak;

    bool fail = false;
    for(int i=0;i<f_size;i++) if( fabs(v1_dbl[i] - v1_test[i])>1e-12 ){ printf("Conversion test fail %d: %f, %f\n",i,v1_dbl[i],v1_test[i]); fail=true; }
    if(fail){
      printf("Conversion test failed\n"); exit(-1);
    }else printf("Conversion test passed\n");
  }


  Float *mom_test1_dbl = (Float *) pmalloc(GJP.VolNodeSites()*18*4*sizeof(Float));
  setup_double_matrixfield((Matrix*)mom_test1_dbl, (Matrix*)mom_test1, 4, gparity_X, gparity_Y);
  pfree(mom_test1);

  Float *mom_test2_dbl = (Float *) pmalloc(GJP.VolNodeSites()*18*4*sizeof(Float));
  setup_double_matrixfield((Matrix*)mom_test2_dbl, (Matrix*)mom_test2, 4, gparity_X, gparity_Y);
  pfree(mom_test2);

  //Double up the Wilson ordered odd-checkerboard solutions of the invcg and minvcg tests
  //First test the doubling by converting v1_dbl to Wilson ordering and comparing to v1_wilsonord doubled
  Float *v1_wilsonord_dbl = (Float *)pmalloc(sizeof(Float) * f_size/2);
  for(int i=0;i<f_size/2;i++) v1_wilsonord_dbl[i] = 0.0;
  
  setup_double_4d_odd_wilson_vector((Vector*)v1_wilsonord_dbl, (Vector*)v1_wilsonord, gparity_X, gparity_Y);
  pfree(v1_wilsonord);

  {
    doubled_lattice.Fconvert((Vector*)v1_dbl, WILSON, CANONICAL);
    long f_size_cb = (long)24 * GJP.VolNodeSites()/2;

    bool fail(false);
    for(int i=0;i<f_size_cb;i++){
      int rem = i; 
      int m = rem %24; rem/=24;
      rem *= 2;
      int pos[4];
      for(int j=0;j<4;j++){
	pos[j] = rem % GJP.NodeSites(j);  rem/= GJP.NodeSites(j);
      }

      //int offset = 24* (index + cboff*parity)/2;

      if(fabs( ((Float*)v1_dbl)[i] - ((Float*)v1_wilsonord_dbl)[i] )>1e-08){
	printf("Odd-CB Wilson ord converion test fail %d,  pos (%d %d %d %d), off %d: %f %f\n",i,pos[0],pos[1],pos[2],pos[3],m,((Float*)v1_dbl)[i],((Float*)v1_wilsonord_dbl)[i] ); fail=true;
      }
    }
    if(fail){ printf("Failed Odd-CB Wilson ord converion test\n"); exit(-1); }
    else printf("Passed Odd-CB Wilson ord converion test\n");

    doubled_lattice.Fconvert((Vector*)v1_dbl,CANONICAL,WILSON);
  }
  pfree(v1_wilsonord_dbl);
  
  //OK, that worked. Now convert solutions of the invcg and minvcg tests
  Float *invcg_test_2fout_dbl = (Float *)pmalloc(sizeof(Float) * f_size/2);
  setup_double_4d_odd_wilson_vector((Vector*)invcg_test_2fout_dbl, (Vector*)invcg_test_2fout, gparity_X, gparity_Y);
  pfree(invcg_test_2fout);

  Vector** minvcg_test_2fout_dbl = (Vector**)pmalloc(3 * sizeof(Vector*));
  for(int j=0;j<3;j++){
    minvcg_test_2fout_dbl[j] = (Vector *)pmalloc(sizeof(Float) * f_size/2);
    setup_double_4d_odd_wilson_vector(minvcg_test_2fout_dbl[j], minvcg_test_2fout[j], gparity_X, gparity_Y);
    pfree(minvcg_test_2fout[j]);
  }
  pfree(minvcg_test_2fout);

  //Also convert 2f eigenvector

  Float *eigenvector_2f_dbl = (Float *)pmalloc(sizeof(Float) * f_size/2);
  for(int i=0;i<f_size/2;i++) eigenvector_2f_dbl[i] = 0.0;
  
  setup_double_4d_odd_wilson_vector((Vector*)eigenvector_2f_dbl, (Vector*)eigenvector_2f, gparity_X, gparity_Y);
  pfree(eigenvector_2f);


  //Do the same as before but in 1f setup

  mom_size = 18*4*GJP.VolNodeSites();

  Float* mom_test1_dbllat = (Float *)pmalloc( sizeof(Float) * mom_size);
  for(int i=0;i<mom_size;i++) mom_test1_dbllat[i] = 0.0;

  Float* mom_test2_dbllat = (Float *)pmalloc( sizeof(Float) * mom_size);
  for(int i=0;i<mom_size;i++) mom_test2_dbllat[i] = 0.0;

  doubled_lattice.Fconvert((Vector*)v1_dbl, WILSON, CANONICAL);
  doubled_lattice.Fconvert((Vector*)v2_dbl, WILSON, CANONICAL);
  doubled_lattice.EvolveMomFforce((Matrix*)mom_test1_dbllat, (Vector*)v1_dbl, 0.5, 0.3, 0.1); 
  doubled_lattice.EvolveMomFforce((Matrix*)mom_test2_dbllat, (Vector*)v1_dbl, (Vector*)v2_dbl, 0.5, 0.3, 0.1); 
  doubled_lattice.Fconvert((Vector*)v1_dbl, CANONICAL, WILSON);
  doubled_lattice.Fconvert((Vector*)v2_dbl, CANONICAL, WILSON);

  //Compare
  bool fail(false);
    
  //compare mom
  for(int i=0;i<mom_size;i++){
    int rem = i;
    int midx = rem % 18; rem/=18;
    int mu = rem % 4; rem/=4;
    int x[4];
    for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

    //printf("XnodeSites is %d, XnodeSites()/2 is %d, this site x coord is %d\n",GJP.XnodeSites(), GJP.XnodeSites()/2 , x[0]);

    if(GJP.Xnodes() == 1 && (x[0] == GJP.XnodeSites()/2-1 || x[0] == GJP.XnodeSites()-1)) continue; //these were surface sites before we doubled the lattice

    if(fabs(mom_test1_dbllat[i] - mom_test1_dbl[i])>1e-08){
      printf("Fermion force fail midx %d mu %d (%d %d %d %d): %f %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test1_dbllat[i],mom_test1_dbl[i]);
      fail=true;
    }//else  printf("Mom test 1 pass midx %d mu %d (%d %d %d %d): %f %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test1_dbllat[i],mom_test1_dbl[i]);
  }
  if(fail){ printf("Failed fermion force test\n"); exit(-1); }
  else printf("Passed fermion force test\n");

  for(int i=0;i<mom_size;i++){
    int rem = i;
    int midx = rem % 18; rem/=18;
    int mu = rem % 4; rem/=4;
    int x[4];
    for(int j=0;j<4;j++){ x[j] = rem % GJP.NodeSites(j) + GJP.NodeSites(j)*GJP.NodeCoor(j); rem /=  GJP.NodeSites(j); }

    //printf("XnodeSites is %d, XnodeSites()/2 is %d, this site x coord is %d\n",GJP.XnodeSites(), GJP.XnodeSites()/2 , x[0]);

    if(GJP.Xnodes() == 1 && (x[0] == GJP.XnodeSites()/2-1 || x[0] == GJP.XnodeSites()-1)) continue; //these were surface sites before we doubled the lattice

    if(fabs(mom_test2_dbllat[i] - mom_test2_dbl[i])>1e-08){
      printf("Boson force fail midx %d mu %d (%d %d %d %d): %f %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test1_dbllat[i],mom_test1_dbl[i]);
      fail=true;
    }//else  printf("Mom test 1 pass midx %d mu %d (%d %d %d %d): %f %f\n",midx,mu,x[0],x[1],x[2],x[3],mom_test1_dbllat[i],mom_test1_dbl[i]);
  }
  if(fail){ printf("Failed boson force test\n"); exit(-1); }
  else printf("Passed boson force test\n");


  //Compare InvCG results 1f/2f
  //Float *invcg_test_2fout;
  {
    doubled_lattice.Fconvert((Vector*)v1_dbl, WILSON, CANONICAL);
  
    CgArg cg_arg;
    cg_arg.mass = -1.8;
    cg_arg.epsilon = 0.5;
    cg_arg.max_num_iter = 10000;
    cg_arg.stop_rsd = 1e-10;
    cg_arg.true_rsd = 1e-10;
    cg_arg.Inverter = CG;
  
    DiracOpWilsonTm dop(doubled_lattice, (Vector*)0, (Vector*)0, &cg_arg, CNV_FRM_NO);

    size_t f_size_cb =  GJP.VolNodeSites() * 24/2;
    invcg_test_2fout = (Float*)pmalloc(f_size_cb*sizeof(Float));
    
    Float* out_1 = (Float*)pmalloc(f_size_cb*sizeof(Float));

    Float true_rsd;
    dop.InvCg((Vector*)out_1, (Vector*)v1_dbl, 0.0, &true_rsd); //bfm version
    
    bool fail(false);
    for(int i=0;i<f_size_cb;i++){
      if( fabs(out_1[i]-invcg_test_2fout_dbl[i])>1e-08 ){
	printf("InvCGtest 1f/2f comparison fail %d: %f %f, ratio %f\n",i,out_1[i],invcg_test_2fout_dbl[i],out_1[i]/invcg_test_2fout_dbl[i]);
	fail=true;
      }
    }
    if(fail){
      printf("Failed 1f/2f comparison InvCg test\n"); exit(-1);
    }else printf("Passed 1f/2f comparison InvCg test\n");
    
    pfree(out_1);
    doubled_lattice.Fconvert((Vector*)v1_dbl, CANONICAL, WILSON);
  }

  //Compare MInvCG results 1f/2f
  {
    doubled_lattice.Fconvert((Vector*)v1_dbl, WILSON, CANONICAL);
    Float targ_resid = 1e-12;

    CgArg cg_arg;
    cg_arg.mass = -1.8;
    cg_arg.epsilon = 0.5;
    cg_arg.max_num_iter = 10000;
    cg_arg.stop_rsd = targ_resid;
    cg_arg.true_rsd = targ_resid;
    cg_arg.Inverter = CG;
  
    DiracOpWilsonTm dop(doubled_lattice, (Vector*)0, (Vector*)0, &cg_arg, CNV_FRM_NO);
    
    size_t f_size_cb =  GJP.VolNodeSites() * 24/2;
    Float src_norm_sq = ((Vector*)v1_dbl)->NormSqNode(f_size_cb);
    glb_sum(&src_norm_sq);
  
    Vector **out_1 = (Vector**)pmalloc(3 * sizeof(Vector*));
  
    //Float mass[3] = {-1.8, -1.0, -0.4};
    for(int j=0;j<3;j++){
      out_1[j] = (Vector*)pmalloc(f_size_cb*sizeof(Float));
      for(int i=0;i<f_size_cb;i++){
	((Float*)out_1[j])[i] = 0.0;
      }
    }
    
    Float true_rsd[3] = {targ_resid, targ_resid, targ_resid};
    
    dop.MInvCG((Vector**)out_1, (Vector*)v1_dbl, src_norm_sq, &mass[0], nmass, 0, &true_rsd[0], MULTI, NULL); //bfm version
  
    bool fail=false;
    for(int j=0;j<nmass;j++){
      Float* o1f = (Float*)out_1[j];
      Float* o2f = (Float*)minvcg_test_2fout_dbl[j];

      for(int i=0;i<f_size_cb;i++){
	if( fabs(o1f[i]-o2f[i])>1e-08 ){
	  printf("MInvCG double latt test %d fail %d: %f %f, ratio %f\n",j,i,o1f[i],o2f[i],o1f[i]/o2f[i]);
	  fail=true;
	}
      }

      if(!fail) printf("Passed MInvCg double latt test for shift %d\n",j);
    }
    if(fail){
      printf("Failed MInvCg double latt test\n"); exit(-1);
    }else printf("Passed MInvCg double latt test\n");
    
    for(int j=0;j<3;j++){
      pfree(out_1[j]);
    }
    pfree(out_1);
    doubled_lattice.Fconvert((Vector*)v1_dbl, CANONICAL, WILSON);
  }

  {
    CommonArg c_arg;
    c_arg.set_filename("pooh");
    c_arg.set_label("pooh");
   
    Float mass = -1.8;

    EigArg e_arg;
    generateEigArg(e_arg, mass, 0.5); //operator is MatPcDagMatPc
    e_arg.Rsdlam = 2e-6;
    e_arg.RsdR_a = e_arg.Rsdlam ;
    e_arg.RsdR_r = e_arg.Rsdlam ;
    e_arg.N_eig = 1;

    CgArg cg_arg;
    cg_arg.mass = mass;
    cg_arg.epsilon = 0.5;

    Float **evalues = (Float**) pmalloc(sizeof(Float*)); //indexed by eigenvalue index, then mass
    for(int i=0;i<1;i++) evalues[i] = (Float*)pmalloc(sizeof(Float));

    LatRanGen LRGbak = LRG;
    AlgEig algeig(doubled_lattice,&c_arg, &e_arg);
    algeig.run(evalues);
    LRG = LRGbak;

    Float evalue_1f = evalues[0][0];
    
    size_t f_size_cb = GJP.VolNodeSites()*24/2;
    Float* eigenvector_1f = (Float*)(algeig.getEigenVectors()[0]);

    //There is a normalization factor of sqrt(2) between the 2f and 1f versions in the 2-directions case 
    //arising from the volume doubling in the second direction. In the 1-direction case the number of spin-color vectors is identical
    //they are just ordered differently. For the 2-directions case the number of spin-color vectors is twice as large (although the
    //extra sites are not independent).
    if(gparity_X && gparity_Y) for(int i=0;i<f_size_cb;i++) eigenvector_1f[i] *= sqrt(2.0);

    bool fail(false);
    for(int i=0;i<f_size_cb;i++){
      int rem = i; 
      int m = rem %24; rem/=24;
      rem *= 2;
      int pos[4];
      for(int j=0;j<4;j++){
	pos[j] = rem % GJP.NodeSites(j);  rem/= GJP.NodeSites(j);
      }
      //pos[0] is special, as rem % Lx gives x[0] - x[0]%2    which maps  x=1 -> 0  x=3 -> 2
      if( (pos[0]+pos[1]+pos[2]+pos[3])%2 == 0 ) pos[0]+=1;

 // ( x[0] - x[0]%2 + GJP.XnodeSites() *
 //           ( x[1] + GJP.YnodeSites() *
 //             (x[2] + GJP.ZnodeSites() * x[3])) +
 //           ((x[0]+x[1]+x[2]+x[3]+1)%2)*GJP.VolNodeSites()
 //         )>>1 ;

      if(fabs( eigenvector_2f_dbl[i] - eigenvector_1f[i] )>1e-08){
	printf("Eigenvector 2f/1f test fail,  pos (%d %d %d %d), off %d: %f %f\n",pos[0],pos[1],pos[2],pos[3],m,eigenvector_2f_dbl[i], eigenvector_1f[i] ); fail=true;
      }
    }
    if(fail){ printf("Failed eigenvector 2f/1f test\n"); exit(-1); }
    else printf("Passed eigenvector 2f/1f test\n");

    if(fabs( evalue_2f - evalue_1f )>1e-08){
      printf("Failed eigenvalue 2f/1f test, %f %f\n",evalue_2f, evalue_1f); exit(-1); 
    }else printf("Passed eigenvalue 2f/1f test\n");
 
    pfree(evalues[0]);
    pfree(evalues);
    pfree(eigenvector_2f_dbl);
  }

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}




void GaugeTransformU(Matrix *gtrans, Lattice &lat){
  Matrix recv_buf;
  Matrix tmp;
  //apply the gauge transformation to U
  int nflav = 1;
  if(GJP.Gparity()) nflav = 2;

  for(int flav=0;flav<nflav;flav++){
    for(int t=0;t<GJP.TnodeSites();t++){
      for(int z=0;z<GJP.ZnodeSites();z++){
	for(int y=0;y<GJP.YnodeSites();y++){
	  for(int x=0;x<GJP.XnodeSites();x++){
	    int pos[4] = {x,y,z,t};
	    int v_x_off = x + GJP.XnodeSites()*(y+GJP.YnodeSites()*(z+GJP.ZnodeSites()*t)) + flav*GJP.VolNodeSites();
	    Matrix &v_x = *(gtrans + v_x_off);

	    for(int mu=0;mu<4;mu++){
	      int u_x_off = lat.GsiteOffset(pos) + mu + flav*4*GJP.VolNodeSites();
	      Matrix &u_x = *(lat.GaugeField() + u_x_off);

	      //get V_x+mu
	      int posp[4] = {x,y,z,t};
	      posp[mu] = (posp[mu]+1)%GJP.NodeSites(mu);

	      Matrix *v_xpmu_ptr = gtrans + posp[0] + GJP.XnodeSites()*(posp[1]+GJP.YnodeSites()*(posp[2]+GJP.ZnodeSites()*posp[3])) + flav*GJP.VolNodeSites();
	      if(pos[mu] == GJP.NodeSites(mu)-1){
		//if node is on the left wall, send the opposite flavour 
		if(GJP.Bc(mu) == BND_CND_GPARITY && GJP.NodeCoor(mu) == 0){
		  if(flav == 1)
		    v_xpmu_ptr-= GJP.VolNodeSites();
		  else
		    v_xpmu_ptr+= GJP.VolNodeSites();		  
		}

		//doesnt need to be fast!
		getPlusData((double *)&recv_buf, (double *)v_xpmu_ptr, 18, mu);
		v_xpmu_ptr = &recv_buf; 
	      }

	      //dagger/transpose it
	      Matrix vdag_xpmu;
	      vdag_xpmu.Dagger(*v_xpmu_ptr);

	      //gauge transform link
	      tmp.DotMEqual(v_x,u_x);
	      u_x.DotMEqual(tmp,vdag_xpmu);
	    }
	  }
	}
      }
    }

  }

}





















//------------------------------------------------------------------
// Constructor 
/*!
  \param latt The lattice on which to compute the condensate.
  \param c_arg The common argument structure for all algorithms.
  \param arg The algorithm parameters.
 */
//------------------------------------------------------------------
AlgEig_CPSeigsolver::AlgEig_CPSeigsolver(Lattice& latt, 
	       CommonArg *c_arg,
	       EigArg *arg) : 
	       Alg(latt, c_arg) 
{
  cname = "AlgEig_CPSeigsolver";
  char *fname = "AlgEig_CPSeigsolver(L&,CommonArg*,EigArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_eig_arg = arg;

  Ncb = NumChkb(alg_eig_arg->RitzMatOper);

  // Set the node size of the full (non-checkerboarded) fermion field
  // NOTE: at this point we must know on what lattice size the operator 
  // will act.
  //----------------------------------------------------------------
  size_t f_size = GJP.VolNodeSites() * latt.FsiteSize() * Ncb / 2;
  if(GJP.Gparity()) f_size*=2;

  VRB.Flow(cname,fname,"f_size=%d\n",0);
//  size_t f_size = GJP.VolNodeSites() * Ncb / 2;
//  exit(1);
  int N_eig = alg_eig_arg->N_eig;

  // Allocate memory for the eigenvectors and eigenvalues
  //----------------------------------------------------------------
  eigenv = (Vector **) smalloc (cname,fname, "eigenv", N_eig * sizeof(Vector *));
  
  for(int n = 0; n < N_eig; ++n)
  {
    eigenv[n] = (Vector *) smalloc(cname,fname, "eigenv[n]", (f_size)* sizeof(Float));
  }

  lambda = (Float *) smalloc(cname,fname, "lambda", 2*N_eig * sizeof(Float));

  chirality = (Float *) smalloc(cname,fname,"chirality", N_eig * sizeof(Float));

  valid_eig = (int *) smalloc(cname,fname,"valid_eig",N_eig * sizeof(int));

  // Print out input parameters
  //----------------------------------------------------------------
  VRB.Input(cname,fname,
	    "N_eig = %d\n",int(N_eig));
  VRB.Input(cname,fname,
	    "MaxCG = %d\n",alg_eig_arg->MaxCG);
  VRB.Input(cname,fname,
	    "Mass_init = %g\n",IFloat(alg_eig_arg->Mass_init));
  VRB.Input(cname,fname,
	    "Mass_final = %g\n",IFloat(alg_eig_arg->Mass_final));
  VRB.Input(cname,fname,
	    "Mass_step = %g\n",IFloat(alg_eig_arg->Mass_step));

  //CK: add check for twisted mass fermions
  if(latt.Fclass() == F_CLASS_WILSON_TM && alg_eig_arg->pattern_kind != ARRAY) ERR.General(cname,fname,"Not implemented for pattern_kind other than ARRAY");

  // Calculate n_masses if necessary
  VRB.Flow(cname,fname,"alg_eig_arg->pattern_kind=%d\n",alg_eig_arg->pattern_kind);
  switch( alg_eig_arg->pattern_kind ) {
  case ARRAY: 
    n_masses = alg_eig_arg->Mass.Mass_len;
    break;
  case LOG:
    n_masses = (int) ((log(alg_eig_arg->Mass_final - alg_eig_arg->Mass_init)
		    / log(alg_eig_arg->Mass_step)) + 1.000001);
    break;
  case LIN:   
    n_masses = (int) (fabs((alg_eig_arg->Mass_final 
      - alg_eig_arg->Mass_init)/alg_eig_arg->Mass_step) + 1.000001); 
    break;
  default: 
    ERR.General(cname, fname,
		"alg_eig_arg->pattern_kind = %d is unrecognized\n", 
		alg_eig_arg->pattern_kind);
    break;
  }  
  VRB.FuncEnd(cname,fname);

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgEig_CPSeigsolver::~AlgEig_CPSeigsolver() {
  char *fname = "~AlgEig_CPSeigsolver()";
  VRB.Func(cname,fname);

  // Free memory
  //----------------------------------------------------------------
  sfree(cname,fname, "valid_eig", valid_eig);

  sfree(cname,fname, "chirality", chirality);

  sfree(cname,fname, "lambda", lambda);

  for(int n = alg_eig_arg->N_eig - 1; n >= 0; --n)
  {
    sfree(cname,fname, "eigenv[n] ",eigenv[n]);
  }

  sfree(cname,fname, "eigenv", eigenv);

  //???
}

//!< Overloaded for backwards compatibility
void AlgEig_CPSeigsolver::run()
{
  run((Float**)0);
}


void gamma_5(IFloat *v_out, IFloat *v_in, int num_sites) 
{
  IFloat *p_out = v_out ;
  IFloat *p_in = v_in ;

  int half_site_size = 12 ;

  for (int site=0; site<num_sites; ++site) {
    int comp;
    for (comp=0; comp<half_site_size; ++comp) {
      *p_out++ = *p_in++ ;
    }
    for (comp=0; comp<half_site_size; ++comp) {
      *p_out++ = -*p_in++ ;
    }
  }
}


//------------------------------------------------------------------
//! Performs the computation.
/*!
  \post The results are written to the file specified in the common_arg
  structure.
*/
//------------------------------------------------------------------
void AlgEig_CPSeigsolver::run(Float **evalues)
{

  Float time = -dclock();
  int iter=0;
  EigArg *eig_arg;
  char *fname = "run()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer eig_arg
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();
  eig_arg = alg_eig_arg;

  const int N_eig = eig_arg->N_eig;
  size_t f_size = GJP.VolNodeSites() * lat.FsiteSize() * Ncb / 2;
  if(GJP.Gparity()) f_size*=2;
  Float **hsum;
  int hsum_len = 0;
  int n;

  if (eig_arg->print_hsum)
  {
    switch(eig_arg->hsum_dir)
    {
    case 0:
      hsum_len = GJP.Xnodes()*GJP.XnodeSites();
      break;
    case 1:
      hsum_len = GJP.Ynodes()*GJP.YnodeSites();
      break;
    case 2:
      hsum_len = GJP.Znodes()*GJP.ZnodeSites();
      break;
    case 3:
      hsum_len = GJP.Tnodes()*GJP.TnodeSites();
      break;
    case 4:
      if (lat.Fclass() == F_CLASS_DWF) 
        hsum_len = GJP.Snodes()*GJP.SnodeSites();
      else
        ERR.General(cname,fname,"Invalid direction\n");
      break;
     default:
      ERR.General(cname,fname,"Invalid direction\n");
    }
    
    if(GJP.Bc(eig_arg->hsum_dir) == BND_CND_GPARITY) hsum_len*=2; //stack hsum for second flavour after first
    hsum = (Float **) smalloc(cname,fname,"hsum",N_eig * sizeof(Float*)); // surely Float* ?
  
    for(n = 0; n < N_eig; ++n)
    {
      hsum[n] = (Float *) smalloc(cname,fname,"hsum[n]",hsum_len * sizeof(Float));
    }
  }
  else
  {
    hsum = (Float **) 0;
  }

  // allocate memory to store the eigenvectors
  Vector** eig_store=0;
  if(eig_arg->ncorr) 
  {
    eig_store = (Vector**)smalloc(cname,fname, "eig_store",N_eig * sizeof(Vector*));
    for(n=0;n<N_eig;++n)
    {
      eig_store[n] = (Vector*) smalloc(cname,fname, "eig_store",f_size * sizeof(Float));
    }
  }


  // Initialize eigenvectors to gaussian
  // and compute eigenvectors
  //----------------------------------------------------------------

  int sign_dm=1;

  // Initialize the cg_arg mass, with the first mass we
  // want to compute for:
  switch( eig_arg->pattern_kind ) {
  case ARRAY: 
    eig_arg->mass = eig_arg->Mass.Mass_val[0]; 
    eig_arg->epsilon = eig_arg->Epsilon.Epsilon_val[0]; //Added by CK for twisted mass fermions
    break;
  case LIN:   
    // Loop over mass values
    sign_dm = (eig_arg->Mass_step < 0.0) ? -1 : 1;
    eig_arg->mass = eig_arg->Mass_init; 
    if (sign_dm*eig_arg->Mass_init > sign_dm*eig_arg->Mass_final)
      ERR.General(cname,fname,"initial and final mass not valid\n");
    break;
  case LOG:   
    eig_arg->mass = eig_arg->Mass_init; 
    break;
  default: 
    ERR.General(cname, fname,
		"eig_arg->pattern_kind = %d is unrecognized\n", 
		eig_arg->pattern_kind);
    break;
  }

  // Loop over masses
  for(int m=0; m<n_masses; m++){
    
    //for(Float mass = eig_arg->Mass_init; 
    //sign_dm*mass <= sign_dm*eig_arg->Mass_final;
    //mass += eig_arg->Mass_step){

    //eig_arg->mass = mass;
    // count the number of masses 
    int count(0);


    // store eigenvectors from previous mass
    if(eig_arg->ncorr && ( count > 0 ) ){
      for ( n=0; n<N_eig; n++ )
      {
        eig_store[n]->CopyVec(eigenv[n],f_size); 
      }
    }


    // random guess every time; do *not* put in the old solution
    for(n = 0; n<N_eig; ++n)
    {
      lat.RandGaussVector(eigenv[n], 0.5, Ncb);

      if(GJP.Gparity1fX() && GJP.Gparity1fY()){
	if(Ncb!=1) ERR.General(cname,fname,"G-parity 1f XY Only set up for odd-even preconditioned fermion vectors\n");
	if(!UniqueID()){ printf("Putting minus sign on fermion source in UR quadrant\n"); fflush(stdout); }
	//make source on upper-right quadrant negative (RNGs should be correct)
	for(int s=0;s<GJP.SnodeSites();s++){
	  for(int t=0;t<GJP.TnodeSites();t++){
	    for(int z=0;z<GJP.ZnodeSites();z++){
	      for(int y=0;y<GJP.YnodeSites();y++){
		for(int x=0;x<GJP.XnodeSites();x++){
		  if( (x+y+z+t+s)%2 == 0) continue; //ferm vect is odd parity only

		  int gx = x+GJP.XnodeCoor()*GJP.XnodeSites();
		  int gy = y+GJP.YnodeCoor()*GJP.YnodeSites();

		  if(gx>=GJP.Xnodes()*GJP.XnodeSites()/2 && gy>=GJP.Ynodes()*GJP.YnodeSites()/2){
		    int pos[5] = {x,y,z,t,s};
		    int f_off = lat.FsiteOffsetChkb(pos) * lat.SpinComponents();

		    for(int spn=0;spn<lat.SpinComponents();spn++) *(eigenv[n]+f_off+spn) *=-1;
		  }
		}
	      }
	    }
	  }
	}
      }


    }

    //==============================================
    // print out mass here in case we want to
    // do any output from within the FeigSolv call
    //==============================================
			      
    if (eig_arg->fname!=0)
    {
      FILE* filep;
      filep=Fopen(eig_arg->fname,"a");
      Fprintf(filep,"mass = %g\n",(Float)eig_arg->mass);
      Fclose(filep); // close file
    }
         
    VRB.Result(cname,fname, "mass = %g\n", (Float)eig_arg->mass);

    // Solve for eigenvectors and eigenvalues.
    // Use eigenv as initial guess. Lambda is not used initially.

    if(Ncb==2)
      iter = CPS_FeigSolv(eigenv, lambda, chirality, valid_eig, 
			  hsum, eig_arg, CNV_FRM_YES,lat);
    else if(Ncb==1)
      iter = CPS_FeigSolv(eigenv, lambda, chirality, valid_eig, 
			  hsum, eig_arg, CNV_FRM_NO,lat);

    //------------------------------------------------------------
    // Solve for eigenvectors and eigenvalues.
    // Lambda is not used initially.
    // This call will return a negative value for the iteration number
    // if either the solver maxes out on one of the limits.

   
    //!< Copy over eigenvalues to return them
    if (evalues != 0) {
      for (int eig=0; eig<eig_arg->N_eig; eig++) {
	evalues[eig][m] = lambda[eig];
      }
    }

    if ( iter < 0 )
    {
      FILE* filep;
      filep=Fopen(eig_arg->fname,"a");
      if ( iter == -1 )
      {
        Fprintf(filep, "maxed out in ritz\n");
      }
      else
      {
        Fprintf(filep, "maxed out for KS steps\n");
      }
      Fclose(filep);
    }
    else
    {   
      //----------------------------------------
      //  spatial correlations of eigen vectors
      //---------------------------------------
      
      if ( eig_arg->fname != 0x0 )
	{
	  
	  FILE* filep;
	  filep=Fopen(eig_arg->fname,"a");
              
	  int i_eig,j_eig;
              
	  // GM5Correlation
              
	  //tmp vector v1
	  Vector* v1 = (Vector *)smalloc(cname, fname, "v1",f_size*sizeof(Float));
              
	  for(i_eig=0;i_eig<N_eig;i_eig++){
	    for(j_eig=i_eig;j_eig<N_eig;j_eig++){
	      gamma_5((Float*)v1, 
		      (Float*)eigenv[j_eig], 
		      f_size/24);
	      cps::Complex cr = eigenv[i_eig]->CompDotProductGlbSum(v1,f_size);
	      Fprintf(filep,"GM5CORR: %d %d %g %g\n",
		      i_eig, j_eig, (Float)cr.real(), (Float)cr.imag());
	    }
	  }
	  sfree(cname,fname, "v1", v1);
              
	  // Correlation with previous eigen vector
	  if(count >0 && eig_arg->ncorr ){
	    for(i_eig=0;i_eig<N_eig;i_eig++){
	      for(j_eig=0;j_eig<N_eig;j_eig++){
		cps::Complex cr = eig_store[i_eig]
		  ->CompDotProductGlbSum(eigenv[j_eig],f_size);
		Fprintf(filep,"NeibCORR: %d %d %g %g\n",
			i_eig, j_eig, (Float)cr.real(), (Float)cr.imag());
	      }
	    }
	  }


	  //------------------------------------------
	  // Print out number of iterations  and hsum
	  //------------------------------------------
	  
	  int i;
	  Fprintf(filep, "  iter = %d\n", iter);
	  if (eig_arg->print_hsum)
	    {
	      for(n = 0; n < eig_arg->N_eig; ++n)
		{
		  for(i = 0; i < hsum_len; ++i)
		    {
		      Fprintf(filep, "  hsum[%d][%d] = %g\n",n,i,
			      (Float)hsum[n][i]);
		    }
		}
	    }
	  Fclose(filep);
	} // output
    } // solver worked

    // If there is another mass loop iteration ahead, we should
    // set the eig_arg->mass to it's next desired value
    if( m < n_masses - 1 ) {
      switch( eig_arg->pattern_kind ) {
      case ARRAY: 
	eig_arg->mass = eig_arg->Mass.Mass_val[m+1];
	eig_arg->epsilon = eig_arg->Epsilon.Epsilon_val[m+1]; //CK: Added for twisted mass fermions
	break;
      case LIN:   
	eig_arg->mass += eig_arg->Mass_step; 
	break;
      case LOG:   
	eig_arg->mass *= eig_arg->Mass_step; 
	break;
      }
    }
    count++;
  } // mass loop

  //==============================
  // deallocate eigenvalue store
  //==============================
  
  if ( eig_arg->ncorr )
    {
      for(n = 0; n < N_eig; ++n)
        {
          sfree(cname,fname,"eig_store[n]",eig_store[n]);
	  if (eig_arg->print_hsum) sfree(cname,fname,"hsum[n]",hsum[n]);
        }
      sfree(cname,fname,"eig_store",eig_store);
      if (eig_arg->print_hsum) sfree(cname,fname,"hsum",hsum);
    }
  time +=dclock();
  print_flops(cname,fname,0,time);

}

Vector ** AlgEig_CPSeigsolver::getEigenVectors(){
  return eigenv;
}











#endif
