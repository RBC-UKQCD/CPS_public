#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of FdwfBase class.

  $Id: f_dwf_base_force.C,v 1.14 2012-08-31 04:55:08 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_dwf_base/noarch/f_dwf_base_force.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_dwf_base_force.C
//
// (R)HMC force term for FdwfBase
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
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
#include <util/time_cps.h>
#include <comms/scu.h> // GRF
#include <comms/glb.h>
CPS_START_NAMESPACE
#undef PROFILE

// CJ: change start
//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *chi, Float mass, 
//                 Float dt):
// It evolves the canonical momentum mom by dt
// using the fermion force.
//------------------------------------------------------------------
ForceArg FdwfBase::EvolveMomFforce(Matrix *mom, Vector *chi, 
			   Float mass, Float dt){
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

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;

  //----------------------------------------------------------------
  // Calculate v1, v2. Both v1, v2 must be in CANONICAL order after
  // the calculation.
  //----------------------------------------------------------------  

  VRB.Clock(cname, fname, "Before calc force vecs.\n") ;
  VRB.Flow(cname, fname, "Before calc force vecs.\n") ;

  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;

    DiracOpDwf dwf(*this, v1, v2, &cg_arg, CNV_FRM_YES) ;
    dwf.CalcHmdForceVecs(chi) ;
  }
  VRB.Flow(cname, fname, "After calc force vecs.\n") ;
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

  Matrix tmp_mat1, tmp_mat2 ;
 
//------------------------------------------------------------------
// start by summing first over direction (mu) and then over site
// to allow SCU transfers to happen face-by-face in the outermost
// loop.
//------------------------------------------------------------------

  VRB.Clock(cname, fname, "Before loop over links.\n") ;

  for (mu=0; mu<4; mu++){
    for (t=0; t<lt; t++){
    for (z=0; z<lz; z++){
    for (y=0; y<ly; y++){
    for (x=0; x<lx; x++){
      int gauge_offset = x+lx*(y+ly*(z+lz*t)) ;
      int vec_offset = f_site_size_4d*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;

      Float *v1_plus_mu ;
      Float *v2_plus_mu ;
      int vec_plus_mu_stride ;
      int vec_plus_mu_offset = f_site_size_4d ;

      Float coeff = -2.0 * dt ;

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
      if(GJP.Snodes() > 1) {
//        if(!UniqueID())printf("%s::%s:GJP.Snodes()=%d\n",cname,fname,GJP.Snodes()); 
	glb_sum_multi_dir((Float *)&tmp_mat1,4,sizeof(Matrix)/sizeof(IFloat));
      }

      tmp_mat2.DotMEqual(*(gauge+gauge_offset), tmp_mat1) ;

      tmp_mat1.Dagger(tmp_mat2) ;

      tmp_mat2.TrLessAntiHermMatrix(tmp_mat1) ;

      tmp_mat2 *= coeff ;

      *(mom+gauge_offset) += tmp_mat2 ;
      Float norm = tmp_mat2.norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);

    } } } } // end for x,y,z,t
  } // end for mu
  ForceFlops += (2*9*16*ls + 18+ 198+36+24)*lx*ly*lz*lt*4;
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
 
  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  VRB.FuncEnd(cname,fname);
  return ForceArg(L1, sqrt(L2), Linf);

}

//------------------------------------------------------------------
// "Odd" fermion force evolution routine written by Chris Dawson, taken 
// verbatim, so performance will suck on qcdoc.
//------------------------------------------------------------------
ForceArg FdwfBase::EvolveMomFforce( Matrix* mom, // momenta
                               Vector* phi, // odd pseudo-fermion field
                               Vector* eta, // very odd pseudo-fermion field
                               Float  mass, 
                               Float dt )
{
  char *fname = "EvolveMomFforce(M*,V*,V*,F,F)";
  VRB.Func(cname,fname);
  
  Matrix *gauge = GaugeField() ;

  if (Colors() != 3)       { ERR.General(cname,fname,"Wrong nbr of colors.") ; }
  if (SpinComponents()!=4) { ERR.General(cname,fname,"Wrong nbr of spin comp.") ;}
  if (mom == 0)            { ERR.Pointer(cname,fname,"mom") ; }
  if (phi == 0)            { ERR.Pointer(cname,fname,"phi") ; }
   
  // allocate space for two CANONICAL fermion fields

  // these are all full fermion vector sizes ( i.e. *not* preconditioned )

  const int f_size        ( FsiteSize() * GJP.VolNodeSites() );
  const int f_size_cb     ( f_size/2 ) ; // f_size must be multiple of 2
  const int f_site_size_4d( 2 * Colors() * SpinComponents() );
  const int f_size_4d     ( f_site_size_4d * GJP.VolNodeSites()) ;
  
  char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v1 == 0) ERR.Pointer(cname, fname, str_v1) ;
  VRB.Smalloc(cname, fname, str_v1, v1, f_size*sizeof(Float)) ;

  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v2 == 0) ERR.Pointer(cname, fname, str_v2) ;
  VRB.Smalloc(cname, fname, str_v2, v2, f_size*sizeof(Float)) ;

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;

#ifdef PROFILE
  Float time = -dclock();
  ForceFlops=0;
#endif

  //Calculate v1, v2. Both must be in CANONICAL order afterwards
  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;
    
    DiracOpDwf dwf(*this, v1, v2, &cg_arg, CNV_FRM_YES) ;
    Float kappa( 1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight())));

    v2->CopyVec(phi,f_size_cb);

    // rescale the input field. As the second half of the this field
    // will be constructed by acting with the PC dslash on v1, this
    // rescales *one* of the full vectors - giving rise to an overall
    // rescaling of the final answer by exactly -\kappa^2
    
    v2->VecTimesEquFloat(-kappa*kappa,f_size_cb);

    // only need one factor of -\kappa^2, so don't rescale the second
    // full vector (v2)
    v1->CopyVec(eta,f_size_cb);
        
    dwf.Dslash(v2+(f_size_cb/6), v2 , CHKB_ODD, DAG_YES);
    dwf.Dslash(v1+(f_size_cb/6), v1 , CHKB_ODD, DAG_NO);
    
    // v1 and v2 are now the vectors needed to contruct the force term
    // written in ( ODD, EVEN ) ordering. They will be converted back
    // into canonical ordering when the destructor is called.
    
  }  

  // two fermion vectors at a single position
  //    - these will be used to store off-node
  //      field components

 
  char *str_site_v1 = "site_v1" ;
  Float *site_v1 = (Float *)smalloc(FsiteSize()*sizeof(Float)) ;
  if (site_v1 == 0) ERR.Pointer(cname, fname, str_site_v1) ;
  VRB.Smalloc(cname, fname, str_site_v1, site_v1, FsiteSize()*sizeof(Float)) ;

  char *str_site_v2 = "site_v2" ;
  Float *site_v2 = (Float *)smalloc(FsiteSize()*sizeof(Float)) ;
  if (site_v2 == 0) ERR.Pointer(cname, fname, str_site_v2) ;
  VRB.Smalloc(cname, fname, str_site_v2, site_v2, FsiteSize()*sizeof(Float)) ;

  // evolve the momenta by the fermion force
  int mu, x, y, z, t, s;
 
  const int lx(GJP.XnodeSites());
  const int ly(GJP.YnodeSites());
  const int lz(GJP.ZnodeSites());
  const int lt(GJP.TnodeSites());
  const int ls(GJP.SnodeSites());
  
  // start by summing first over direction (mu) and then over site to
  // allow SCU transfers to happen face-by-face in the outermost loop.

  VRB.Clock(cname, fname, "Before loop over links.\n") ;

  for (mu=0; mu<4; mu++) {
    for (t=0; t<lt; t++){
      for (z=0; z<lz; z++){
        for (y=0; y<ly; y++){
          for (x=0; x<lx; x++) {
            // position offset
            int gauge_offset = x+lx*(y+ly*(z+lz*t));
            
            // offset for vector field at this point
            // (4d only, no fifth dimension)
            int vec_offset = f_site_size_4d*gauge_offset ;
            
            // offset for link in mu direction from this point
            gauge_offset = mu+4*gauge_offset ; 
            
            Float *v1_plus_mu=NULL ;
            Float *v2_plus_mu=NULL ;
            int vec_plus_mu_stride=0 ;
            int vec_plus_mu_offset = f_site_size_4d ;
            
            // sign of coeff (look at momenta update)
            Float coeff = -2.0 * dt ;
            
            switch (mu) 
              {
              case 0 :
                // next position in mu direction
                vec_plus_mu_offset *= (x+1)%lx+lx*(y+ly*(z+lz*t)) ;
                // vec_plus_mu_offset now the correct
                // offset for a fermion field at this point
                // in the lattice 
                if ((x+1) == lx) 
                  {
                    // off-node
                    for (s=0; s<ls; s++) 
                      {
                        // fill site_v1 and site_v2 with v1 and v2 data
                        // from x=0 on next node, need loop because
                        // data is not contiguous in memory 
                        getPlusData( (Float *)site_v1+s*f_site_size_4d,
                                     (Float *)v1+vec_plus_mu_offset+s*f_size_4d,
                                     f_site_size_4d, mu) ;
                        getPlusData( (Float*)site_v2+s*f_site_size_4d,
                                     (Float*)v2+vec_plus_mu_offset+s*f_size_4d,
                                     f_site_size_4d, mu) ;
                      } // end for s
                    
                    v1_plus_mu = site_v1   ;  
                    v2_plus_mu = site_v2   ;  
                    vec_plus_mu_stride = 0 ;  // field now contiguous
                    
                    // GJP.XnodeBc() gives the forward boundary
                    // condition only (so this should work).
                    if (GJP.XnodeBc()==BND_CND_APRD) coeff = -coeff ;
                  } 
                else 
                  {
                    // on - node
                    //
                    // just add offset to v1 and v2
                    // (they are now 1 forward in the mu direction )
                    //
                    v1_plus_mu = (Float*)v1+vec_plus_mu_offset ;
                    v2_plus_mu = (Float*)v2+vec_plus_mu_offset ;
                    vec_plus_mu_stride = f_size_4d - f_site_size_4d ; // explained below
                  }
                break ;
                // Repeat for the other directions
              case 1 :
                vec_plus_mu_offset *= x+lx*((y+1)%ly+ly*(z+lz*t)) ;
                if ((y+1) == ly) 
                  {
                    for (s=0; s<ls; s++) 
                      {
                        getPlusData( (Float*)site_v1+s*f_site_size_4d,
                                     (Float*)v1+vec_plus_mu_offset+s*f_size_4d,
                                     f_site_size_4d, mu) ;
                        getPlusData( (Float*)site_v2+s*f_site_size_4d,
                                     (Float*)v2+vec_plus_mu_offset+s*f_size_4d,
                                     f_site_size_4d, mu) ;
                      }
                    v1_plus_mu = site_v1 ;
                    v2_plus_mu = site_v2 ;
                    vec_plus_mu_stride = 0 ;
                    if (GJP.YnodeBc()==BND_CND_APRD) coeff = -coeff ;
                  } 
                else 
                  {
                    v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
                    v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
                    vec_plus_mu_stride = f_size_4d - f_site_size_4d ;
                  }
                break ;
              case 2 :
                vec_plus_mu_offset *= x+lx*(y+ly*((z+1)%lz+lz*t)) ;
                if ((z+1) == lz) 
                  {
                    for (s=0; s<ls; s++) 
                      {
                        getPlusData( (Float*)site_v1+s*f_site_size_4d,
                                     (Float*)v1+vec_plus_mu_offset+s*f_size_4d,
                                     f_site_size_4d, mu) ;
                        getPlusData( (Float*)site_v2+s*f_site_size_4d,
                                     (Float*)v2+vec_plus_mu_offset+s*f_size_4d,
                                     f_site_size_4d, mu) ;
                      }
                    v1_plus_mu = site_v1 ;
                    v2_plus_mu = site_v2 ;
                    vec_plus_mu_stride = 0 ;
                    if (GJP.ZnodeBc()==BND_CND_APRD) coeff = -coeff ;
                  } 
                else 
                  {
                    v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
                    v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
                    vec_plus_mu_stride = f_size_4d - f_site_size_4d ;
                  }
                break ;
              case 3 :
                vec_plus_mu_offset *= x+lx*(y+ly*(z+lz*((t+1)%lt))) ;
                if ((t+1) == lt) 
                  {
                    for (s=0; s<ls; s++) 
                      {
                        getPlusData( (Float*)site_v1+s*f_site_size_4d,
                                     (Float*)v1+vec_plus_mu_offset+s*f_size_4d,
                                     f_site_size_4d, mu) ;
                        getPlusData( (Float*)site_v2+s*f_site_size_4d,
                                     (Float*)v2+vec_plus_mu_offset+s*f_size_4d,
                                     f_site_size_4d, mu) ;
                      } 
                    v1_plus_mu = site_v1 ;
                    v2_plus_mu = site_v2 ;
                    vec_plus_mu_stride = 0 ;
                    if (GJP.TnodeBc()==BND_CND_APRD) coeff = -coeff ;
                  } 
                else 
                  {
                    v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
                    v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
                    vec_plus_mu_stride = f_size_4d - f_site_size_4d ;
                  }
              } // end (the evil) mu switch 


            Matrix tmp_mat1, tmp_mat2;  

            // the non-zero stride pattern is due to domain wall
            // fermions ( summing up *ls* different sproj's )
            //
            // f_size_4d-f_site_size_4d is the number of floats
            // between the end of one spinor at s and the start of the 
            // spinor at s+1 
            // 
            // vec_plus_mu_stride is the same, except when
            // this is off boundary, in that case the info
            // is copied into a contiguous block in the above code
            // and vec_plus_mu_stride set to zero
            
            // ( 1 - gamma_\mu ) Tr_s [ v1(x+\mu) v2^{\dagger}(x) ]
            
            sproj_tr[mu]( (Float *)&tmp_mat1,
                          (Float *)v1_plus_mu,
                          (Float *)v2+vec_offset,
                          ls, vec_plus_mu_stride, f_size_4d-f_site_size_4d) ;
            
            // (1 + gamma_\mu)  Tr_s [ v2(x+\mu) v1^{\dagger}(x) ]
            sproj_tr[mu+4]( (Float *)&tmp_mat2,
                            (Float *)v2_plus_mu,
                            (Float *)v1+vec_offset,
                            ls, vec_plus_mu_stride, f_size_4d-f_site_size_4d) ;
            
            
            // exactly what is sounds like
            tmp_mat1 += tmp_mat2 ;
            
            if(GJP.Snodes() != 1) {
              for (s=0; s<(sizeof(Matrix)/sizeof(Float)); ++s) {
                glb_sum_dir((Float *)&tmp_mat1 + s, 4) ;
              }
            }
            
            // multiply sum by the link in the \mu direction
            tmp_mat2.DotMEqual(*(gauge+gauge_offset), tmp_mat1) ;
            
            // take tracless antihermitian piece
            // TrLessAntiHermMatrix need to be passed
            // the dagger of the matrix in question
            tmp_mat1.Dagger(tmp_mat2) ;
            tmp_mat2.TrLessAntiHermMatrix(tmp_mat1) ;

            tmp_mat2 *= coeff ;
            
            // note the minus sign.
            *(mom+gauge_offset) -= tmp_mat2 ;
	    Float norm = tmp_mat2.norm();
	    Float tmp = sqrt(norm);
	    L1 += tmp;
	    L2 += norm;
	    Linf = (tmp>Linf ? tmp : Linf);
	    
          } // end for x
        } // end for y
      } // end for z
    } // end for t
  } // end for mu
  ForceFlops += (2*9*16*ls + 18+ 198+36+24)*lx*ly*lz*lt*4;
#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif
  
  // deallocate smalloc'd space

  VRB.Sfree(cname, fname, str_site_v2, site_v2) ;
  sfree(site_v2) ;
 
  VRB.Sfree(cname, fname, str_site_v1, site_v1) ;
  sfree(site_v1) ;
 
  VRB.Sfree(cname, fname, str_v2, v2) ;
  sfree(v2) ;
  
  VRB.Sfree(cname, fname, str_v1, v1) ;
  sfree(v1) ;

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  VRB.FuncEnd(cname,fname);
  return ForceArg(L1, sqrt(L2), Linf);

}
// CJ: change end

CPS_END_NAMESPACE
