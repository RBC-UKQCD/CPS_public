#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of FdwfBase class.

*/
//--------------------------------------------------------------------
//
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
  if(GJP.Gparity()) return EvolveMomFforceGparity(mom,chi,mass,dt);

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

  size_t f_size = FsiteSize() * GJP.VolNodeSites() ;
  int f_site_size_4d = 2 * Colors() * SpinComponents();
  size_t f_size_4d = f_site_size_4d * GJP.VolNodeSites() ;
 
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

  Matrix *gparity_1f_mombuf;
  if(GJP.Gparity1fX()){
    gparity_1f_mombuf = (Matrix *)fmalloc(cname,fname,"gparity_1f_mombuf",4 * GJP.VolNodeSites() * sizeof(Matrix) ) ;
    for(int i=0;i<4*GJP.VolNodeSites();i++) gparity_1f_mombuf[i].ZeroMatrix();
  }

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
    //v1 becomes (1+D)\chi and v2 becomes (D^\dagger-\kappa^2 M)\chi
    // M is the odd-even preconditioned fermion matrix connecting odd to
    // odd parity sites and \e D is the hopping term connecting odd to
    // even parity sites. Recall that \a chi is defined on odd sites only.
    // The new vectors are in odd-even order.
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
      int vec_offset = f_site_size_4d*gauge_offset ; //offset of current site
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


      // These functions compute 
      // \f$ f_{ij} = Tr_{spins}[ (1 \pm \gamma_\mu) v_i w^\dagger_j ] \f$
      // for spin-colour vectors \a v and \a w where \a i and \a j are
      // colour indices.

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

      tmp_mat2.DotMEqual(*(gauge+gauge_offset), tmp_mat1) ; //F = U dS/sU

      tmp_mat1.Dagger(tmp_mat2) ; //why the dagger when we just make an antiherm matrix out of it below?

      tmp_mat2.TrLessAntiHermMatrix(tmp_mat1) ;

      tmp_mat2 *= coeff ;

      *(mom+gauge_offset) += tmp_mat2 ;

      if(GJP.Gparity1fX()) *(gparity_1f_mombuf+gauge_offset) = tmp_mat2;
      
      Float norm = tmp_mat2.norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);

    } } } } // end for x,y,z,t
  } // end for mu
  ForceFlops += (2*9*16*ls + 18+ 198+36+24)*lx*ly*lz*lt*4;


  //G-parity 1f add \delta p from other 'flavour'
  //Providing you set up the 1f prop source correctly (minus sign on UR quadrant) this code also works for
  //1f G-parity in both X and Y directions
  if(GJP.Gparity1fX()){
    int momsz = GsiteSize() * GJP.VolNodeSites();
    Matrix *buf2 = (Matrix *)fmalloc(cname,fname,"buf2",momsz * sizeof(Float) ) ;

    //Communicate \delta p from first half onto second half and vice versa
    Matrix *data_buf = gparity_1f_mombuf;
    Matrix *send_buf = data_buf;
    Matrix *recv_buf = buf2;

    int vol = GJP.VolNodeSites();
    int size[4] = {GJP.XnodeSites(),GJP.YnodeSites(),GJP.ZnodeSites(),GJP.TnodeSites()};

    if(GJP.Xnodes()>1){
      //pass between nodes
      for(int i=0;i<GJP.Xnodes()/2;i++){
	getMinusData((Float *)recv_buf, (Float *)send_buf, momsz , 0);
	data_buf = recv_buf;
	recv_buf = send_buf;
	send_buf = data_buf;
      }
    }else{
      //shift field by xsites/2
#ifdef USE_OMP
#pragma omp parallel for default(shared) private(mu)
      for(long i=0;i<vol*4;i++){
#else
	for(long i=0;i<vol*4;i++){
#endif
	  //i = mu + 4*(x + Lx*(y+Ly*(z+Lz*t) ) )
	  int mu = i%4;
	  int x = (i/4) % GJP.XnodeSites();
	  int pos_rem = i/4/GJP.XnodeSites(); //(y+Ly*(z+Lz*t)

	  int x_from = (x + GJP.XnodeSites()/2) % GJP.XnodeSites();
	  int i_from = mu + 4*(x_from + GJP.XnodeSites()*pos_rem);
	  buf2[i] = gparity_1f_mombuf[i_from];
	}
	data_buf = buf2;
      }

#ifdef USE_OMP
#pragma omp parallel for default(shared) private(mu)
      for (long i=0;i<vol*4;i++){
  	int pos[4];
	long rest=i;
	mu = rest%4; rest = rest/4;
	for(int j =0; j<4;j++){
	  pos[j]= rest%size[j]; rest = rest/size[j];
	}
#else
	int pos[4];
	for (mu=0; mu<4; mu++)
	  for (pos[3]=0; pos[3]<size[3]; pos[3]++)
	    for (pos[2]=0; pos[2]<size[2]; pos[2]++)
	      for (pos[1]=0; pos[1]<size[1]; pos[1]++)
		for (pos[0]=0; pos[0]<size[0]; pos[0]++){
#endif
		  int gauge_offset = mu + 4*(pos[0]+size[0]*(pos[1]+size[1]*(pos[2]+size[2]*pos[3])));
      
		  //complex conjugate the \delta p from the other flavour
		  Float *m = (Float*) &data_buf[gauge_offset];
		  for(int c=1;c<18;c+=2) m[c] *= -1;
      
		  //add it to the momentum at this site
		  *(mom+gauge_offset) += data_buf[gauge_offset];
		}
        ffree(buf2);
      }

    
      if(GJP.Gparity1fX()) ffree(gparity_1f_mombuf,cname,fname,"gparity_1f_mombuf");


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

#if 0
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

  const size_t f_size        ( FsiteSize() * GJP.VolNodeSites() );
  const size_t f_size_cb     ( f_size/2 ) ; // f_size must be multiple of 2
  const int f_site_size_4d( 2 * Colors() * SpinComponents() );
  const size_t f_size_4d     ( f_site_size_4d * GJP.VolNodeSites()) ;
  
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
#endif
// CJ: change end
//------------------------------------------------------------------
// EvolveMomFforceGparity(Matrix *mom, Vector *chi, Float mass, 
//                 Float dt):
// It evolves the canonical momentum mom by dt for G-parity Bcs
// using the fermion force.
//------------------------------------------------------------------
ForceArg FdwfBase::EvolveMomFforceGparity(Matrix *mom, Vector *chi, 
			   Float mass, Float dt){
  char *fname = "EvolveMomFforceGparity(M*,V*,F,F,F)";
  VRB.Func(cname,fname);
  
  if(!GJP.Gparity()) ERR.General(cname,fname,"Function executed when G-parity BCs not in use\n");

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

  size_t f_size = FsiteSize() * GJP.VolNodeSites() ;
  int f_site_size_4d = 2 * Colors() * SpinComponents();
  size_t f_size_4d = f_site_size_4d * GJP.VolNodeSites() ;

  //CK: Need space for both d and C\bar u^T fields stacked
  //Two 4d volumes are stacked on each Ls such that Ls can be split over multiple nodes
  size_t f_size_alloc = f_size *2 *sizeof(Float);
  int f_single4dsite_alloc = FsiteSize()*sizeof(Float)*2;


  char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)smalloc(f_size_alloc) ;
  if (v1 == 0) ERR.Pointer(cname, fname, str_v1) ;
  VRB.Smalloc(cname, fname, str_v1, v1, f_size_alloc) ;

  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)smalloc(f_size_alloc) ;
  if (v2 == 0) ERR.Pointer(cname, fname, str_v2) ;
  VRB.Smalloc(cname, fname, str_v2, v2, f_size_alloc) ;

  //----------------------------------------------------------------
  // allocate buffer space for two fermion fields that are assoc
  // with only one 4-D site.
  //----------------------------------------------------------------

  char *str_site_v1 = "site_v1" ;
  Float *site_v1 = (Float *)smalloc(f_single4dsite_alloc) ;
  if (site_v1 == 0) ERR.Pointer(cname, fname, str_site_v1) ;
  VRB.Smalloc(cname, fname, str_site_v1, site_v1, f_single4dsite_alloc) ;

  char *str_site_v2 = "site_v2" ;
  Float *site_v2 = (Float *)smalloc(f_single4dsite_alloc) ;
  if (site_v2 == 0) ERR.Pointer(cname, fname, str_site_v2) ;
  VRB.Smalloc(cname, fname, str_site_v2, site_v2, f_single4dsite_alloc) ;

  //CK: G-parity site vectors on each Ls are 2 * flavours 4 spin * 3 clr * 2 re/im
  //Sticking with same layout as 5d guys, stack 2 site vectors on each Ls
  //Ls=0 [ d ][ CubarT ] Ls=1 [ d ][ CubarT ] .....

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
    //Net result is:  (odd, even)
    // f_in = ( -kappa^2 Mprec chi, -kappa^2 D_eo^dag Mprec chi )    f_out = ( chi, D_eo chi )
    // kappa = 1/[2(5-M5)]
    // These are converted into CANONICAL format when DiracOpDWF is destroyed
    // here we have f_out = v1,   f_in = v2
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

  Matrix tmp_mat1, tmp_mat2, tmp_mat3;
 
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

      int pos[] = {x,y,z,t};
      int lattsz[] = {lx,ly,lz,lt};    

      Float *v1_plus_mu ;
      Float *v2_plus_mu ;

      Float *v1_plus_mu_CubarT;
      Float *v2_plus_mu_CubarT;      

      int vec_plus_mu_stride ;
      int vec_plus_mu_offset = f_site_size_4d ;

      Float d_coeff = -2.0 * dt ; //coefficient for first G-parity field
      Float CubarT_coeff = d_coeff; //coefficient for second G-parity field

      /* CK: Replaced nasty 300000 line switch statement with the following:*/
      {
	int pos_p_mu[] = {x,y,z,t};
	pos_p_mu[mu] = (pos_p_mu[mu]+1) % lattsz[mu];
	vec_plus_mu_offset *= pos_p_mu[0] + lx*(pos_p_mu[1] + ly*(pos_p_mu[2] + lz*pos_p_mu[3]));
      }      

      if ((pos[mu]+1) == lattsz[mu]) {
	if(GJP.Bc(mu) == BND_CND_GPARITY && GJP.NodeCoor(mu)==GJP.Nodes(mu)-1){
	  //at global boundary
	  // d <- +CubarT,  CubarT <- -d
	  CubarT_coeff = -CubarT_coeff;
	}else if(GJP.NodeBc(mu)==BND_CND_APRD){
	  // d <- -d,   CubarT <- -CubarT
	  d_coeff = -d_coeff;
	  CubarT_coeff = -CubarT_coeff;
	}

	for (s=0; s<ls; s++) {
	  //if on global lattice boundary we must perform the G-parity flavor swap
	  int d_site_vect_offset = s*f_site_size_4d*2; //where d-fields are stored in site vect
	  int CubarT_site_vect_offset = d_site_vect_offset + f_site_size_4d; //same for CubarT field (offset by 1 4d vect from d)
	      
	  //normally d <- d,  CubarT <- CubarT
	  int d_comms_to_offset = d_site_vect_offset; //place in site vect where comm'd d-field will be stored
	  int CubarT_comms_to_offset = CubarT_site_vect_offset; //....

	  if(GJP.Bc(mu) == BND_CND_GPARITY && GJP.NodeCoor(mu)==GJP.Nodes(mu)-1){
	    //at global boundary
	    // d <- +CubarT,  CubarT <- -d
	    CubarT_comms_to_offset = d_site_vect_offset;
	    d_comms_to_offset = CubarT_site_vect_offset;
	  }

	  //communicate d field
	  getPlusData( (IFloat *)site_v1+d_comms_to_offset,
		       (IFloat *)v1+vec_plus_mu_offset+s*f_size_4d*2,
		       f_site_size_4d, mu) ;
	  getPlusData( (IFloat *)site_v2+d_comms_to_offset,
		       (IFloat *)v2+vec_plus_mu_offset+s*f_size_4d*2,
		       f_site_size_4d, mu) ;
	  //communicate CubarT field
	  getPlusData( (IFloat *)site_v1+CubarT_comms_to_offset,
		       (IFloat *)v1 + vec_plus_mu_offset + f_size_4d + s*f_size_4d*2,
		       f_site_size_4d, mu) ;
	  getPlusData( (IFloat *)site_v2+CubarT_comms_to_offset,
		       (IFloat *)v2 + vec_plus_mu_offset + f_size_4d + s*f_size_4d*2,
		       f_site_size_4d, mu) ;
	} // end for s
	v1_plus_mu = site_v1 ;
	v2_plus_mu = site_v2 ;

	v1_plus_mu_CubarT = site_v1 + f_site_size_4d;
	v2_plus_mu_CubarT = site_v2 + f_site_size_4d;

	vec_plus_mu_stride = f_site_size_4d ; //skip over other stacked field
      } else {
	v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
	v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;

	v1_plus_mu_CubarT = (Float *)v1+vec_plus_mu_offset + f_size_4d ;
	v2_plus_mu_CubarT = (Float *)v2+vec_plus_mu_offset + f_size_4d ;

	vec_plus_mu_stride = 2*f_size_4d - f_site_size_4d ;
      }
      
      //Calculate contributions to force at site
      //from d field

      // These functions compute 
      // \f$ f_{ij} = Tr_{spins}[ (1 \pm \gamma_\mu) v_i w^\dagger_j ] \f$
      // for spin-colour vectors \a v and \a w where \a i and \a j are
      // colour indices.
      //v1 becomes (1+D)\chi and v2 becomes (D^\dagger-\kappa^2 M)\chi  (M o->o, D o->e)

      // f_{ij} = Tr_{spins}[ (1 \pm \gamma_\mu) v_i w^\dagger_j ]

      //tmp_mat1_{ij} = tr_s ( P_+^mu v1(x+mu)_i v2^\dagger(x)_j 
      sproj_tr[mu]( (IFloat *)&tmp_mat1,
                    (IFloat *)v1_plus_mu,
                    (IFloat *)v2+vec_offset,
                    ls, vec_plus_mu_stride, 2*f_size_4d-f_site_size_4d) ;

      //tmp_mat2_{ij} = tr_s ( P_-^mu v2(x+mu)_i v1^\dagger(x)_j 
      sproj_tr[mu+4]( (IFloat *)&tmp_mat2,
                      (IFloat *)v2_plus_mu,
                      (IFloat *)v1+vec_offset,
                      ls, vec_plus_mu_stride, 2*f_size_4d-f_site_size_4d) ;

      tmp_mat1 += tmp_mat2 ;
      tmp_mat1 *= d_coeff; //moved from final stage

      //for flavour 2 field we must flip the projection operator used for the force term
      //and swap the coordinates of the force vectors

      //tmp_mat2_{ij} = tr_s ( P_-^mu v1(x)_i v2^\dagger(x+mu)_j 
      sproj_tr[mu+4]( (IFloat *)&tmp_mat2,
		      (IFloat *)v1+vec_offset+f_size_4d,
		      (IFloat *)v2_plus_mu_CubarT,
		      ls, 2*f_size_4d-f_site_size_4d, vec_plus_mu_stride) ;

      //tmp_mat3_{ij} = tr_s ( P_+^mu v2(x)_i v1^\dagger(x+mu)_j 
      sproj_tr[mu]( (IFloat *)&tmp_mat3,
		    (IFloat *)v2+vec_offset+f_size_4d,
		    (IFloat *)v1_plus_mu_CubarT,
		    ls, 2*f_size_4d-f_site_size_4d, vec_plus_mu_stride) ;
      
      tmp_mat2 += tmp_mat3;

      tmp_mat3.Trans(tmp_mat2);//must transpose outer product on colour index
      tmp_mat3 *= CubarT_coeff;//different coefficient for the two flavours
      tmp_mat1 += tmp_mat3 ;

      // If GJP.Snodes > 1 sum up contributions from all s nodes
      if(GJP.Snodes() > 1) {
	glb_sum_multi_dir((Float *)&tmp_mat1,4,sizeof(Matrix)/sizeof(IFloat));
      }

      tmp_mat2.DotMEqual(*(gauge+gauge_offset), tmp_mat1) ;
      tmp_mat1.Dagger(tmp_mat2) ;
      tmp_mat2.TrLessAntiHermMatrix(tmp_mat1) ;
      

      *(mom+gauge_offset) += tmp_mat2 ;
      Float norm = tmp_mat2.norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);

      //set force for U* link
      Matrix* mom_u = mom+gauge_offset;
      Matrix* mom_ustar = mom+gauge_offset+GJP.VolNodeSites()*4;
      mom_ustar->Conj((IFloat*)mom_u);


    } } } } // end for x,y,z,t
  } // end for mu
  ForceFlops += (2*9*16*ls + 18+ 198+36+24)*lx*ly*lz*lt*4; //almost certainly wrong!
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

//CK: This code is #if'd out in the f_dwf_base.C in the root directory in the latest CPS but a replacement implementation
//    is only provided for the qmp type, not noarch. This fixes the problem
#if 1
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
  if(GJP.Gparity()) { ERR.General(cname,fname,"Not implemented for G-parity boundary conditions") ; }
  // allocate space for two CANONICAL fermion fields

  // these are all full fermion vector sizes ( i.e. *not* preconditioned )

  const size_t f_size        ( FsiteSize() * GJP.VolNodeSites() );
  const size_t f_size_cb     ( f_size/2 ) ; // f_size must be multiple of 2
  const int f_site_size_4d( 2 * Colors() * SpinComponents() );
  const size_t f_size_4d     ( f_site_size_4d * GJP.VolNodeSites()) ;
  
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
#endif


CPS_END_NAMESPACE

