#include <config.h>
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
#include <comms/scu.h>
#include <comms/glb.h>
CPS_START_NAMESPACE

FimprDwf::FimprDwf(){
}

FimprDwf::~FimprDwf(){
}


//------------------------------------------------------------------
// int FmatEvlInv(Vector *f_out, Vector *f_in, 
//                CgArg *cg_arg, 
//                Float *true_res,
//		  CnvFrmType cnv_frm = CNV_FRM_YES):
// It calculates f_out where A * f_out = f_in and
// A is the preconditioned fermion matrix that appears
// in the HMC evolution (even/odd  preconditioning 
// of [Dirac^dag Dirac]. The inversion is done
// with the conjugate gradient. cg_arg is the structure
// that contains all the control parameters, f_in is the
// fermion field source vector, f_out should be set to be
// the initial guess and on return is the solution.
// f_in and f_out are defined on a checkerboard.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int FimprDwf::FmatEvlInv(Vector *f_out, Vector *f_in, 
		     CgArg *cg_arg, 
		     Float *true_res,
		     CnvFrmType cnv_frm)
{
  int iter;
  char *fname = "FmatEvlInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  DiracOpImprDwf dwf(*this, f_out, f_in, cg_arg, cnv_frm);
  
  iter = dwf.InvCg(true_res);

  // Return the number of iterations
  return iter;
}

//------------------------------------------------------------------
// SetPhi(Vector *phi, Vector *frm1, Vector *frm2, Float mass):
// It sets the pseudofermion field phi from frm1, frm2.
// Note that frm2 is not used.
//------------------------------------------------------------------
void FimprDwf::SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		  Float mass){
  char *fname = "SetPhi(V*,V*,V*,F)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = mass;

  if (phi == 0)
    ERR.Pointer(cname,fname,"phi") ;

  if (frm1 == 0)
    ERR.Pointer(cname,fname,"frm1") ;

  DiracOpImprDwf dwf(*this, frm1, 0, &cg_arg, CNV_FRM_NO) ;

  dwf.MatPcDag(phi, frm1) ;

  return ;
}

//------------------------------------------------------------------
// Float BhamiltonNode(Vector *boson, Float mass):
// The boson Hamiltonian of the node sublattice.
//------------------------------------------------------------------
Float FimprDwf::BhamiltonNode(Vector *boson, Float mass){
  char *fname = "BhamiltonNode(V*,F)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = mass;

  if (boson == 0)
    ERR.Pointer(cname,fname,"boson");

  int f_size = GJP.VolNodeSites() * FsiteSize() / 2 ;

  Vector *bsn_tmp = (Vector *)
    smalloc(f_size*sizeof(Float));

  char *str_tmp = "bsn_tmp" ;

  if (bsn_tmp == 0)
    ERR.Pointer(cname,fname,str_tmp) ;

  VRB.Smalloc(cname,fname,str_tmp,bsn_tmp,f_size*sizeof(Float));

  DiracOpImprDwf dwf(*this, boson, bsn_tmp, &cg_arg, CNV_FRM_NO) ;

  dwf.MatPc(bsn_tmp,boson);

  Float ret_val = bsn_tmp->NormSqNode(f_size) ;

  VRB.Sfree(cname,fname,str_tmp,bsn_tmp);

  sfree(bsn_tmp) ;

  // Sum accross s nodes in case Snodes() != 1
  glb_sum_dir(&ret_val, 4) ;

  return ret_val ;
}

//------------------------------------------------------------------
// int FmatInv(Vector *f_out, Vector *f_in, 
//             CgArg *cg_arg, 
//             Float *true_res,
//             CnvFrmType cnv_frm = CNV_FRM_YES,
//             PreserveType prs_f_in = PRESERVE_YES):
// It calculates f_out where A * f_out = f_in and
// A is the fermion matrix (Dirac operator). The inversion
// is done with the conjugate gradient. cg_arg is the 
// structure that contains all the control parameters, f_in 
// is the fermion field source vector, f_out should be set 
// to be the initial guess and on return is the solution.
// f_in and f_out are defined on the whole lattice.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// cnv_frm is used to specify if f_in should be converted 
// from canonical to fermion order and f_out from fermion 
// to canonical. 
// prs_f_in is used to specify if the source
// f_in should be preserved or not. If not the memory usage
// is less by half the size of a fermion vector.
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int FimprDwf::FmatInv(Vector *f_out, Vector *f_in, 
		  CgArg *cg_arg, 
		  Float *true_res,
		  CnvFrmType cnv_frm,
		  PreserveType prs_f_in)
{
  int iter;
  char *fname = "FmatInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  DiracOpImprDwf dwf(*this, f_out, f_in, cg_arg, cnv_frm);
    
  iter = dwf.MatInv(true_res, prs_f_in);

  // Return the number of iterations
  return iter;
}


//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *chi, Float mass, 
//                 Float step_size):
// It evolves the canonical momentum mom by step_size
// using the fermion force.
//------------------------------------------------------------------
void FimprDwf::EvolveMomFforce(Matrix *mom, Vector *chi, 
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

    DiracOpImprDwf dwf(*this, v1, v2, &cg_arg, CNV_FRM_YES) ;
    dwf.CalcHmdForceVecs(chi) ;
  }

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
#if 0
    
    const int N = 4;  // N can be 1, 2 or 4.
	Fconvert(tmp,CANONICAL,STAG);

    enum{plus, minus, n_sign};
    int size;

    // Vectors for which we must allocate memory
    Vector *Pnu[n_sign][N];
    Vector *P3[n_sign][n_sign][N];
    Vector *Prhonu[n_sign][n_sign][N];
    Vector *P5[n_sign][n_sign][n_sign][N];
    Vector *P7[n_sign][n_sign][n_sign][n_sign][N];
    Vector *Psigma7[n_sign][n_sign][n_sign][n_sign][N];
    
    ParTransDwf pt(*this);
    size = GJP.VolNodeSites()*sizeof(Vector);

    for(int w=0; w<N; w++){
	for(int i=0; i<n_sign; i++){
	    Pnu[i][w] = (Vector*)pmalloc(size);
	    if(!Pnu[i]) ERR.Pointer(cname, fname, "Pnu");
	    for(int j=0; j<n_sign; j++){
		P3[i][j][w] = (Vector*)pmalloc(size);
		if(!P3[i][j][w]) ERR.Pointer(cname, fname, "P3");

		Prhonu[i][j][w] = (Vector*)pmalloc(size);
		if(!Prhonu[i][j][w]) ERR.Pointer(cname, fname, "Prhonu");

		for(int k=0; k<n_sign; k++){
		    P5[i][j][k][w] = (Vector*)pmalloc(size);
		    if(!P5[i][j][k][w]) ERR.Pointer(cname, fname, "P5");

  		    for(int l=0; l<n_sign; l++){
	 		P7[i][j][k][l][w] = (Vector*)pmalloc(size);
			if(!P7[i][j][k][l][w]) ERR.Pointer(cname, fname, "P7");

			Psigma7[i][j][k][l][w] = (Vector*)pmalloc(size);
			if(!Psigma7[i][j][k][l][w]) ERR.Pointer(cname, fname, "P7sigma");
		    }
		}
	    }
	}
    }

    // These vectors can be overlapped with previously allocated memory
    
    Vector *Psigmarhonu[n_sign][n_sign][n_sign][N];
    Vector *Prho5[n_sign][n_sign][n_sign][N];
    Vector *Pnunu[n_sign][N];
    Vector *Pnu5[n_sign][n_sign][N];
    Vector *Pnu3[n_sign][n_sign][N];
    Vector *Pnununu[N];

    // If the vector arrays were suitably dynamically allocated we wouldn't
    // have to do all this loopy stuff.
    
    for(int w=0; w<N; w++){
	Pnununu[w] = Prhonu[0][0][w];
	for(int i=0; i<n_sign; i++){
	    Pnunu[i][w] = Psigma7[0][0][0][i][w];;
	    for(int j=0; j<n_sign; j++){
		Pnu5[i][j][w] = P7[0][0][i][j][w];
		Pnu3[i][j][w] = P7[0][0][i][j][w];
		for(int k=0; k<n_sign; k++){
		    Prho5[i][j][k][w] = Psigma7[0][i][j][k][w];
		    Psigmarhonu[i][j][k][w] = Psigma7[0][i][j][k][w];
		}
	    }
	}
    }

    Vector *X = tmp;	

    // Array in which to accumulate the force term
    // This must be initialised to zero 

    Matrix *force[4];
    for(int i=0; i<4; i++){
	force[i] = (Matrix*)pmalloc(GJP.VolNodeSites()*sizeof(Matrix));         
        VRB.Pmalloc(cname, fname,"force[i]", force[i], GJP.VolNodeSites()*sizeof(Matrix));
	if(!force[i]) ERR.Pointer(cname, fname, "force[i]");
    }
   
    // input/output arrays for the parallel transport routines
    Vector *vin[n_sign*N], *vout[n_sign*N];
    int dir[n_sign*N];
	
    int mu[N], nu[N], rho[N], sigma[N];   // Sets of directions
    int w;                                // The direction index 0...N-1
    int ms, ns, rs, ss;                   // Sign of direction
    bool done[4] = {false,false,false,false};  // Flags to tell us which 
                                               // nu directions we have done.

	    
    for (int m=0; m<4; m+=N){                     	    // Loop over mu
	for(w=0; w<N; w++) mu[w] = (m+w)%4; 

	for (int n=m+1; n<m+4; n++){                        // Loop over nu
	    for(w=0; w<N; w++) nu[w] = (n+w)%4;

	    // Pnu = U_nu X

	    for(int i=0; i<N; i++){
		vin[i] = vin[i+N] = X;
		dir[n_sign*i] = (n_sign*nu[i]);        // nu_i
		dir[n_sign*i+1] = (n_sign*nu[i]+1);    // -nu_i
		vout[n_sign*i] = Pnu[minus][i];
		vout[n_sign*i+1] = Pnu[plus][i];
	    }
//	    parallel_transport(pt,vin, dir, n_sign*N, vout);
	    pt.run(n_sign*N, vout, vin, dir);

	    // P3 = U_mu Pnu
	    // ms is the nu sign index, ms is the mu sign index,
	    // w is the direction index
	    for(int i=0; i<N; i++){
		dir[n_sign*i] = (n_sign*mu[i]);        // mu_i
		dir[n_sign*i+1] = (n_sign*mu[i]+1);    // -mu_i
	    }
	    for(ns=0; ns<n_sign; ns++){               // ns is the sign of nu
		for(int i=0; i<N; i++){
		    vin[n_sign*i] = vin[n_sign*i+1] = Pnu[ns][i];
		    vout[n_sign*i] = P3[plus][ns][i];
		    vout[n_sign*i+1] = P3[minus][ns][i];
		}
//		parallel_transport(pt,vin, dir, n_sign*N, vout);
	    	pt.run(n_sign*N, vout, vin, dir);
	    }
	    
	    for(w=0; w<N; w++)
		for(ns=0; ns<n_sign; ns++)
		    force_product_sum(P3[plus][ns][w], Pnu[ns][w],
				       -GJP.staple3_coeff(),
				       force[mu[w]]);

	    for(int r=n+1; r<n+4; r++){                     // Loop over rho
		bool nextr = false;
		for(w=0; w<N; w++){
		    rho[w] = (r+w)%4;		
		    if(rho[w]==mu[w]){
			nextr = true;
			break;
		    }
		}
		if(nextr) continue;

		for(w=0; w<N; w++){                         // sigma
		    for(int s=rho[w]+1; s<rho[w]+4; s++){
			sigma[w] = s%4;
			if(sigma[w]!=mu[w] && sigma[w]!=nu[w]) break;
		    }
		}

		// Prhonu = U_rho Pnu 

		for(int i=0; i<N; i++){
		    dir[n_sign*i] = (n_sign*rho[i]);        
		    dir[n_sign*i+1] = (n_sign*rho[i]+1);    
		}
		for(ns=0; ns<n_sign; ns++){
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = vin[n_sign*i+1] = Pnu[ns][i];
			vout[n_sign*i] = Prhonu[ns][minus][i];
			vout[n_sign*i+1] = Prhonu[ns][plus][i];
		    }
//			parallel_transport(pt,vin, dir, n_sign*N, vout);
			pt.run(n_sign*N, vout, vin, dir);
		}

		// P5 = U_mu Prhonu

		for(int i=0; i<N; i++){
		    dir[n_sign*i] = (n_sign*mu[i]);        
		    dir[n_sign*i+1] = (n_sign*mu[i]+1);    
		}
		for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) {
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = vin[n_sign*i+1] = Prhonu[ns][rs][i];
			vout[n_sign*i] = P5[plus][ns][rs][i];
			vout[n_sign*i+1] = P5[minus][ns][rs][i];
		    }
//		    parallel_transport(pt,vin, dir, n_sign*N, vout);
	    	pt.run(n_sign*N, vout, vin, dir);
		}

		// F_mu += P5 Prhonu^dagger
		      
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++)
			force_product_sum(P5[plus][ns][rs][w],
					  Prhonu[ns][rs][w],
					  GJP.staple5_coeff(),
					  force[mu[w]]);

		// Psigmarhonu = U_sigma P_rhonu
		
		for(int i=0; i<N; i++){
		    dir[n_sign*i] = (n_sign*sigma[i]);        
		    dir[n_sign*i+1] = (n_sign*sigma[i]+1);    
		}
		for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++){
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = vin[n_sign*i+1] = Prhonu[ns][rs][i];
			vout[n_sign*i] = Psigmarhonu[ns][rs][minus][i];
			vout[n_sign*i+1] = Psigmarhonu[ns][rs][plus][i];
		    }
//			parallel_transport(pt,vin, dir, n_sign*N, vout);
	    	pt.run(n_sign*N, vout, vin, dir);
		}

		// P7 = U_mu P_sigmarhonu
		for(int i=0; i<N; i++){
		    dir[n_sign*i] = (n_sign*mu[i]);        
		    dir[n_sign*i+1] = (n_sign*mu[i]+1);    
		}
		for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) for(ss=0; ss<n_sign; ss++){
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = vin[n_sign*i+1] = Psigmarhonu[ns][rs][ss][i];
			vout[n_sign*i] = P7[plus][ns][rs][ss][i];
			vout[n_sign*i+1] = P7[minus][ns][rs][ss][i];
		    }
//		    parallel_transport(pt,vin, dir, n_sign*N, vout);
	    	pt.run(n_sign*N, vout, vin, dir);
		}

		// F_mu -= P7 Psigmarhonu^\dagger
		
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) for(ss=0; ss<n_sign; ss++)
			force_product_sum(P7[plus][ns][rs][ss][w],
					   Psigmarhonu[ns][rs][ss][w],
					   -GJP.staple7_coeff(),
					   force[mu[w]]);

		// F_sigma += P7 Psigmarhonu^\dagger
		// N.B. this is the same as one of the previous products.
		
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
			force_product_sum(P7[plus][ns][rs][minus][w],
					  Psigmarhonu[ns][rs][minus][w],
					  GJP.staple7_coeff(),
					  force[sigma[w]]);

		// F_sigma += Psigmarhonu P7^\dagger
		
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
			force_product_sum(Psigmarhonu[ns][rs][minus][w],
					  P7[minus][ns][rs][minus][w],
					  GJP.staple7_coeff(),
					  force[sigma[w]]);

		// Psigma7 = U_sigma P7 
		for(int i=0; i<N; i++){
		    dir[n_sign*i] = (n_sign*sigma[i]);        
		    dir[n_sign*i+1] = (n_sign*sigma[i]+1);    
		}
		for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++){
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = P7[ms][ns][rs][plus][i];
			vin[n_sign*i+1] = P7[ms][ns][rs][minus][i];
			vout[n_sign*i] = Psigma7[ms][ns][rs][plus][i];
			vout[n_sign*i+1] = Psigma7[ms][ns][rs][minus][i];
		    }
//		    parallel_transport(pt,vin, dir, n_sign*N, vout);
	    	pt.run(n_sign*N, vout, vin, dir);
		}

		// F_sigma += Fsigma7 Frhonu^\dagger

	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
			force_product_sum(Psigma7[plus][ns][rs][plus][w],
					  Prhonu[ns][rs][w],
					  GJP.staple7_coeff(),
					  force[sigma[w]]);

		// F_sigma += Frhonu Fsigma7^\dagger

	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
			force_product_sum(Prhonu[ns][rs][w],
					  Psigma7[minus][ns][rs][plus][w],
					  GJP.staple7_coeff(),
					  force[sigma[w]]);

		// P5 += c_7/c_5 Psigma7

		if(GJP.staple5_coeff()!=0.0)
		    for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) for(ss=0; ss<n_sign; ss++) for(w=0; w<N; w++)
			P5[ms][ns][rs][w]->FTimesV1PlusV2(GJP.staple7_coeff()/GJP.staple5_coeff(), Psigma7[ms][ns][rs][ss][w], P5[ms][ns][rs][w], GJP.VolNodeSites()*VECT_LEN);

		// F_rho -= P5 Prhonu^\dagger
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++)
			force_product_sum(P5[plus][ns][minus][w],
					   Prhonu[ns][minus][w],
					   -GJP.staple5_coeff(),
					   force[rho[w]]);

		// F_rho -= Prhonu P5^\dagger
		    
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++)
			force_product_sum(Prhonu[ns][minus][w],
					   P5[minus][ns][minus][w],
					   -GJP.staple5_coeff(),
					   force[rho[w]]);

		// Prho5 = U_rho P5

		for(int i=0; i<N; i++){
		    dir[n_sign*i] = (n_sign*rho[i]);        
		    dir[n_sign*i+1] = (n_sign*rho[i]+1);    
		}
		for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++){
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = P5[ms][ns][plus][i];
			vin[n_sign*i+1] = P5[ms][ns][minus][i];
			vout[n_sign*i] = Prho5[ms][ns][plus][i];
			vout[n_sign*i+1] = Prho5[ms][ns][minus][i];
		    }
//		    parallel_transport(pt,vin, dir, n_sign*N, vout);
	    	pt.run(n_sign*N, vout, vin, dir);
		}

		// F_rho -= Prho5 Pnu^\dagger
		
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++)
			force_product_sum(Prho5[plus][ns][plus][w],
					   Pnu[ns][w],
					   -GJP.staple5_coeff(),
					   force[rho[w]]);

		// F_rho -= Pnu Prho5^\dagger
		
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++)
			force_product_sum(Pnu[ns][w],
					   Prho5[minus][ns][plus][w],
					   -GJP.staple5_coeff(),
					   force[rho[w]]);
		
		// P3 += c_5/c_3 Prho5

		if(GJP.staple3_coeff()!=0.0)		
		    for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) for(w=0; w<N; w++)
			P3[ms][ns][w]->FTimesV1PlusV2(GJP.staple5_coeff()/GJP.staple3_coeff(), Prho5[ms][ns][rs][w], P3[ms][ns][w], GJP.VolNodeSites()*VECT_LEN);

	    } // rho+sigma loop

	    // Pnunu = U_nu Pnu

	    for(int i=0; i<N; i++){
		dir[n_sign*i] = (n_sign*nu[i]);        
		dir[n_sign*i+1] = (n_sign*nu[i]+1);    
	    }
	    for(int i=0; i<N; i++){
		vin[n_sign*i] = Pnu[minus][i];
		vin[n_sign*i+1] = Pnu[plus][i];
		vout[n_sign*i] = Pnunu[minus][i];
		vout[n_sign*i+1] = Pnunu[plus][i];
	    }
//	    parallel_transport(pt,vin, dir, n_sign*N, vout);
	    pt.run(n_sign*N, vout, vin, dir);

	    // P5 = U_mu Pnunu

	    for(int i=0; i<N; i++){
		dir[n_sign*i] = (n_sign*mu[i]);        
		dir[n_sign*i+1] = (n_sign*mu[i]+1);    
	    }
	    for(ns=0; ns<n_sign; ns++){
		for(int i=0; i<N; i++){
		    vin[n_sign*i] = Pnunu[ns][i];
		    vin[n_sign*i+1] = Pnunu[ns][i];
		    vout[n_sign*i] = P5[plus][ns][0][i];
		    vout[n_sign*i+1] = P5[minus][ns][0][i];
		}
//		parallel_transport(pt,vin, dir, n_sign*N, vout);
	    pt.run(n_sign*N, vout, vin, dir);
	    }

	    // F_mu += P5 Pnunu^\dagger

	    for(w=0; w<N; w++)
		for(ns=0; ns<n_sign; ns++)
		    force_product_sum(P5[plus][ns][0][w],
				      Pnunu[ns][w],
				      GJP.Lepage_coeff(),
				      force[mu[w]]);

	    // F_nu -= P5 Pnunu^\dagger
	    // N.B. this is the same as one of the previous products
	    
	    for(w=0; w<N; w++)
		force_product_sum(P5[plus][minus][0][w],
				   Pnunu[minus][w],
				   -GJP.Lepage_coeff(),
				   force[nu[w]]);
	    
	    // F_nu -= Pnunu P5^\dagger
	    
	    for(w=0; w<N; w++)
		force_product_sum(Pnunu[minus][w],
				   P5[minus][minus][0][w],
				   -GJP.Lepage_coeff(),
				   force[nu[w]]);

	    // Pnu5 = U_nu P5

	    for(int i=0; i<N; i++){
		dir[n_sign*i] = (n_sign*nu[i]);        
		dir[n_sign*i+1] = (n_sign*nu[i]+1);    
	    }
	    for(ms=0; ms<n_sign; ms++){
		for(int i=0; i<N; i++){
		    vin[n_sign*i] =   P5[ms][plus][0][i]; 
		    vin[n_sign*i+1] = P5[ms][minus][0][i];
		    vout[n_sign*i] =   Pnu5[ms][plus][i];
		    vout[n_sign*i+1] = Pnu5[ms][minus][i];
		}
//		parallel_transport(pt,vin, dir, n_sign*N, vout);
	    pt.run(n_sign*N, vout, vin, dir);
	    }

	    // F_nu -= Pnu5 Pnu^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnu5[plus][plus][w],
				   Pnu[plus][w],
				   -GJP.Lepage_coeff(),
				   force[nu[w]]);

	    // F_nu -= Pnu Pnu5^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnu[plus][w],
				   Pnu5[minus][plus][w],
				   -GJP.Lepage_coeff(),
				   force[nu[w]]);

	    // P3 += c_L/c_3 Pnu5

	    if(GJP.staple3_coeff()!=0.0)
		for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++) for(w=0; w<N; w++)
		    P3[ms][ns][w]->FTimesV1PlusV2(GJP.Lepage_coeff()/GJP.staple3_coeff(), Pnu5[ms][ns][w], P3[ms][ns][w], GJP.VolNodeSites()*VECT_LEN);

	    // F_nu += P3 Pnu^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(P3[plus][minus][w],
				  Pnu[minus][w],
				  GJP.staple3_coeff(),
				  force[nu[w]]);

	    // F_nu +=  Pnu P3^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnu[minus][w],
				  P3[minus][minus][w],
				  GJP.staple3_coeff(),
				  force[nu[w]]);
	    
	    // Pnu3 = U_nu P3

	    for(int i=0; i<N; i++)
		dir[i] = (n_sign*nu[i]);        
	    for(ms=0; ms<n_sign; ms++){
		for(int i=0; i<N; i++){
		    vin[i] = P3[ms][plus][i]; 
		    vout[i] = Pnu3[ms][plus][i];
		}
		parallel_transport(pt,vin, dir, N, vout);
	    }

	    // F_nu += Pnu3 X^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnu3[plus][plus][w], X,
				  GJP.staple3_coeff(),
				  force[nu[w]]);

	    // F_nu += X Pnu3^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(X, Pnu3[minus][plus][w], 
				  GJP.staple3_coeff(),
				  force[nu[w]]);

	    // This stuff is to be done once only for each value of nu[w].
	    // Look for  N nu's that haven't been done before.

    	    bool nextn = false;
	    for(w=0; w<N; w++)
		if(done[nu[w]]){
		    nextn = true;
		    break;
		}
	    if(nextn) continue;
	    for(w=0; w<N; w++) done[nu[w]] = true;

	    // Got N new nu's, so do some stuff...
	    
	    // F_nu += Pnu X^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnu[minus][w], X,
				  GJP.KS_coeff(),
				  force[nu[w]]);

	} // nu loop
    } // mu loop


    // Now that we have computed the force, we can update the momenta

    Matrix mf, mfd; 
    for(int s=0; s<GJP.VolNodeSites(); s++){
	Matrix *ip = mom+4*s;
	for (int mu=0; mu<4; mu++){
	    mf = force[mu][s];
	    mfd.Dagger((IFloat*)&mf);
	    mf.TrLessAntiHermMatrix(mfd);
	    mf *= 0.5;				// for MILC compatibility
	    fTimesV1PlusV2((IFloat*)(ip+mu), dt, (IFloat*)&mf,
			   (IFloat*)(ip+mu)+BANK4_BASE, MATRIX_SIZE);
	}
    }
    
    for(int w=0; w<N; w++){
	for(int i=0; i<n_sign; i++){
	    pfree(Pnu[i][w]);
	    for(int j=0; j<n_sign; j++){
		pfree(P3[i][j][w]);
		pfree(Prhonu[i][j][w]);
		for(int k=0; k<n_sign; k++){
		    pfree(P5[i][j][k][w]);
  		    for(int l=0; l<n_sign; l++){
	 		pfree(P7[i][j][k][l][w]);
			pfree(Psigma7[i][j][k][l][w]);
		    }
		}
	    }
	}
    }
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

void FimprDwf::ForceProductSum(const Vector *v, const Vector *w,
		IFloat coeff, Matrix *f){
}

CPS_END_NAMESPACE
