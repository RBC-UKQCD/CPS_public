#include<config.h>
#include <util/qcdio.h>
#include <math.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/dwf.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/random.h>
#include <util/error.h>
#include <util/time.h>
#include <comms/scu.h> // GRF
#include <comms/glb.h>

CPS_START_NAMESPACE

#define PROFILE

//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *chi, Float mass, 
//                 Float step_size):
// It evolves the canonical momentum mom by step_size
// using the fermion force.
//------------------------------------------------------------------
void FdwfBase::EvolveMomFforce(Matrix *mom, Vector *chi, 
			   Float mass, Float step_size){
  char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
  VRB.Func(cname,fname);

  Matrix *gauge = GaugeField() ;

  if (Colors() != 3)
    ERR.General(cname,fname,"Wrong nbr of colors.") ;
 
  if (SpinComponents() != 4)
    ERR.General(cname,fname,"Wrong nbr of spin comp.") ;
 
  if (mom == 0)
    ERR.Pointer(cname,fname,"mom") ;
 
  if (chi == 0)
    ERR.Pointer(cname,fname,"chi") ;
 
  //----------------------------------------------------------------
  // allocate space for two CANONICAL fermion fields
  //----------------------------------------------------------------

  int f_size = FsiteSize() * GJP.VolNodeSites() ;
  int f_site_size_4d = 2 * Colors() * SpinComponents();
  int f_size_4d = f_site_size_4d * GJP.VolNodeSites() ;
 
  char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v1 == 0) ERR.Pointer(cname, fname, str_v1) ;
  VRB.Smalloc(cname, fname, str_v1, v1, f_size*sizeof(Float)) ;

  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v2 == 0) ERR.Pointer(cname, fname, str_v2) ;
  VRB.Smalloc(cname, fname, str_v2, v2, f_size*sizeof(Float)) ;

  //----------------------------------------------------------------
  // allocate buffer space for two fermion fields that are assoc
  // with only one 4-D site.
  //----------------------------------------------------------------

  char *str_site_v1 = "site_v1" ;
  Float *site_v1 = (Float *)smalloc(FsiteSize()*sizeof(Float)) ;
  if (site_v1 == 0) ERR.Pointer(cname, fname, str_site_v1) ;
  VRB.Smalloc(cname, fname, str_site_v1, site_v1, FsiteSize()*sizeof(Float)) ;

  char *str_site_v2 = "site_v2" ;
  Float *site_v2 = (Float *)smalloc(FsiteSize()*sizeof(Float)) ;
  if (site_v2 == 0) ERR.Pointer(cname, fname, str_site_v2) ;
  VRB.Smalloc(cname, fname, str_site_v2, site_v2, FsiteSize()*sizeof(Float)) ;


  //----------------------------------------------------------------
  // Calculate v1, v2. Both v1, v2 must be in CANONICAL order after
  // the calculation.
  //----------------------------------------------------------------  

  VRB.Clock(cname, fname, "Before calc force vecs.\n") ;

  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;

    DiracOpDwf dwf(*this, v1, v2, &cg_arg, CNV_FRM_YES) ;
    dwf.CalcHmdForceVecs(chi) ;
  }
#ifdef PROFILE
  Float time = -dclock();
  ForceFlops=0;
#endif
  int mu, x, y, z, t, s, lx, ly, lz, lt, ls ;
 
  lx = GJP.XnodeSites() ;
  ly = GJP.YnodeSites() ;
  lz = GJP.ZnodeSites() ;
  lt = GJP.TnodeSites() ;
  ls = GJP.SnodeSites() ;
  int vol=lx*ly*lz*lt;
  Matrix tmp_mat1, tmp_mat2 ;
 
//------------------------------------------------------------------
// start by summing first over direction (mu) and then over site
// to allow SCU transfers to happen face-by-face in the outermost
// loop.
//------------------------------------------------------------------

  VRB.Clock(cname, fname, "Before loop over links.\n") ;

  for (mu=0; mu<4; mu++) {
    for (t=0; t<lt; t++)
    for (z=0; z<lz; z++)
    for (y=0; y<ly; y++)
    for (x=0; x<lx; x++) {
	  
      int gauge_offset = x+lx*(y+ly*(z+lz*t)) ;
      int vec_offset = f_site_size_4d*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;

      Float *v1_plus_mu ;
      Float *v2_plus_mu ;
      int vec_plus_mu_stride ;
      int vec_plus_mu_offset = f_site_size_4d ;

      Float coeff = -2.0 * step_size ;

      switch (mu) {
        case 0 :
          vec_plus_mu_offset *= (x+1)%lx+lx*(y+ly*(z+lz*t)) ;
          if ((x+1) == lx) {
            for (s=0; s<ls; s++) {
              getPlusData( (IFloat *)site_v1+s*f_site_size_4d,
                (IFloat *)v1+vec_plus_mu_offset+s*f_size_4d,
                f_site_size_4d, mu) ;
              getPlusData( (IFloat *)site_v2+s*f_site_size_4d,
                (IFloat *)v2+vec_plus_mu_offset+s*f_size_4d,
                f_site_size_4d, mu) ;
            } // end for s
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            vec_plus_mu_stride = 0 ;
            if (GJP.XnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
            vec_plus_mu_stride = f_size_4d - f_site_size_4d ;
          }
          break ;
        case 1 :
          vec_plus_mu_offset *= x+lx*((y+1)%ly+ly*(z+lz*t)) ;
          if ((y+1) == ly) {
            for (s=0; s<ls; s++) {
              getPlusData( (IFloat *)site_v1+s*f_site_size_4d,
                (IFloat *)v1+vec_plus_mu_offset+s*f_size_4d,
                f_site_size_4d, mu) ;
              getPlusData( (IFloat *)site_v2+s*f_site_size_4d,
                (IFloat *)v2+vec_plus_mu_offset+s*f_size_4d,
                f_site_size_4d, mu) ;
            } // end for s
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            vec_plus_mu_stride = 0 ;
            if (GJP.YnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
            vec_plus_mu_stride = f_size_4d - f_site_size_4d ;
          }
          break ;
        case 2 :
          vec_plus_mu_offset *= x+lx*(y+ly*((z+1)%lz+lz*t)) ;
          if ((z+1) == lz) {
            for (s=0; s<ls; s++) {
              getPlusData( (IFloat *)site_v1+s*f_site_size_4d,
                (IFloat *)v1+vec_plus_mu_offset+s*f_size_4d,
                f_site_size_4d, mu) ;
              getPlusData( (IFloat *)site_v2+s*f_site_size_4d,
                (IFloat *)v2+vec_plus_mu_offset+s*f_size_4d,
                f_site_size_4d, mu) ;
            } // end for s
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            vec_plus_mu_stride = 0 ;
            if (GJP.ZnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
            vec_plus_mu_stride = f_size_4d - f_site_size_4d ;
          }
          break ;
        case 3 :
          vec_plus_mu_offset *= x+lx*(y+ly*(z+lz*((t+1)%lt))) ;
          if ((t+1) == lt) {
            for (s=0; s<ls; s++) {
              getPlusData( (IFloat *)site_v1+s*f_site_size_4d,
                (IFloat *)v1+vec_plus_mu_offset+s*f_size_4d,
                f_site_size_4d, mu) ;
              getPlusData( (IFloat *)site_v2+s*f_site_size_4d,
                (IFloat *)v2+vec_plus_mu_offset+s*f_size_4d,
                f_site_size_4d, mu) ;
            } // end for s
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            vec_plus_mu_stride = 0 ;
            if (GJP.TnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
            vec_plus_mu_stride = f_size_4d - f_site_size_4d ;
          }
      } // end switch mu 

      sproj_tr[mu]( (IFloat *)&tmp_mat1,
                    (IFloat *)v1_plus_mu,
                    (IFloat *)v2+vec_offset,
                    ls, vec_plus_mu_stride, f_size_4d-f_site_size_4d) ;

      sproj_tr[mu+4]( (IFloat *)&tmp_mat2,
                      (IFloat *)v2_plus_mu,
                      (IFloat *)v1+vec_offset,
                      ls, vec_plus_mu_stride, f_size_4d-f_site_size_4d) ;

      tmp_mat1 += tmp_mat2 ;

      // If GJP.Snodes > 1 sum up contributions from all s nodes
      if(GJP.Snodes() != 1) {
	  glb_sum_multi_dir((Float *)&tmp_mat1, 4, sizeof(Matrix)/sizeof(Float)) ;
      }

      tmp_mat2.DotMEqual(*(gauge+gauge_offset), tmp_mat1) ;

      tmp_mat1.Dagger(tmp_mat2) ;

      tmp_mat2.TrLessAntiHermMatrix(tmp_mat1) ;

      tmp_mat2 *= coeff ;

      *(mom+gauge_offset) += tmp_mat2 ;

    } // end for x,y,z,t
  } // end for mu

  ForceFlops += (2*9*16*ls + 18+ 198+36+24)*vol*4;	  
 
#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif
 
//------------------------------------------------------------------
// deallocate smalloc'd space
//------------------------------------------------------------------
  VRB.Sfree(cname, fname, str_site_v2, site_v2) ;
  sfree(site_v2) ;
 
  VRB.Sfree(cname, fname, str_site_v1, site_v1) ;
  sfree(site_v1) ;
 
  VRB.Sfree(cname, fname, str_v2, v2) ;
  sfree(v2) ;
 
  VRB.Sfree(cname, fname, str_v1, v1) ;
  sfree(v1) ;
 
  return ;
}
CPS_END_NAMESPACE
