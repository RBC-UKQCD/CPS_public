#include<config.h>
#include<stdio.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Fasqtad::EvolveMomFforce.

  $Id: Fforce.C,v 1.1 2004-01-28 06:06:22 cwj Exp $
*/
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/pt.h>
#include <util/gjp.h>
#include <stdlib.h>
#include <comms/scu_enum.h>
#include <comms/nga_reg.h>
CPS_START_NAMESPACE

const char * par_trans_filename = CWDPREFIX("par_trans");



// The parity of the lattice site n

ChkbType Fasqtad::parity(const int *n){

    int d, p;
    
    for(p=0, d=0; d<4; d++) p += n[d];
    return  p%2?CHKB_ODD:CHKB_EVEN;
    
}


// The outer product of v and w, multiplied by the parity factor and
// coeffient, and added to the force f.
// N.B. The force term is multiplied by -1 on ODD parity sites. This the
// the MILC convention and is here to test against the MILC code.
// This should eventually be changed to EVEN for the CPS. 

void Fasqtad::force_product_sum(const Vector *v, const Vector *w,
				    IFloat coeff, Matrix *f){

    Matrix m0;
    int n, s[4];
    
    for(n=0, s[3]=0; s[3]<GJP.TnodeSites(); s[3]++)
	for(s[2]=0; s[2]<GJP.ZnodeSites(); s[2]++)
	    for(s[1]=0; s[1]<GJP.YnodeSites(); s[1]++)
		for(s[0]=0; s[0]<GJP.XnodeSites(); s[0]++, n++){
		    m0.Cross2(v[n], w[n]);
		    if(parity(s)==CHKB_ODD) m0 *= -1.0;    // MILC stylee!
		    m0 *= coeff;
		    f[n] += m0;
		}

}



// N.B. No optimising provision is made if any of the asqtad coefficients
// are zero.

void Fasqtad::EvolveMomFforce(Matrix *mom, Vector *frm,
				  Float mass, Float dt){

    char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
    VRB.Func(cname,fname);
    int pos[] = {0,0,0,0};
    int f_size = GJP.VolNodeSites()*FsiteSize()*sizeof(Float)/2;

    Vector *tmp = (Vector *)smalloc(2*f_size);
    Vector *tmp_o = &(tmp[GJP.VolNodeSites()/2]);
    moveMem(tmp,frm,f_size);
#if 1
    moveMem(tmp_o,f_tmp,f_size);
#else
    {
    CgArg cg_arg;
    cg_arg.mass=0.0;
    DiracOpAsqtad asqtad(*this,tmp_o,frm,&cg_arg,CNV_FRM_NO);
      asqtad.Dslash(tmp_o,frm,CHKB_EVEN,DAG_NO);
    }
#endif
    Convert(STAG);

    
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
    
    ParTransAsqtad pt(*this);
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
	    parallel_transport(pt,vin, dir, n_sign*N, vout);

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
		parallel_transport(pt,vin, dir, n_sign*N, vout);
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
		    parallel_transport(pt,vin, dir, n_sign*N, vout);
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
		    parallel_transport(pt,vin, dir, n_sign*N, vout);
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
		    parallel_transport(pt,vin, dir, n_sign*N, vout);
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
		    parallel_transport(pt,vin, dir, n_sign*N, vout);
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
		    parallel_transport(pt,vin, dir, n_sign*N, vout);
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
		    parallel_transport(pt,vin, dir, n_sign*N, vout);
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
	    parallel_transport(pt,vin, dir, n_sign*N, vout);

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
		parallel_transport(pt,vin, dir, n_sign*N, vout);
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
		parallel_transport(pt,vin, dir, n_sign*N, vout);
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

	    // F_nu += Pnu Pnunu^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnunu[minus][w], Pnu[plus][w],
				  -GJP.Naik_coeff(),
				  force[nu[w]]);

	    // F_nu += Pnunu Pnu^\dagger
		
	    for(w=0; w<N; w++)
		force_product_sum(Pnu[minus][w], Pnunu[plus][w],
				  GJP.Naik_coeff(),
				  force[nu[w]]);

	    // Pnununu = U_nu Pnunu

	    for(int i=0; i<N; i++){
		dir[i] = (n_sign*nu[i]);        
		vin[i] = Pnunu[minus][i]; 
		vout[i] = Pnununu[i];
	    }
	    parallel_transport(pt,vin, dir, N, vout);
	    
	    // F_nu += Pnununu X^\dagger
		
	    for(w=0; w<N; w++)
		force_product_sum(Pnununu[w], X,
				  GJP.Naik_coeff(),
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
	
   for(int i=0;i<4;i++) pfree(force[i]);
 	sfree(tmp);
//	Fconvert(frm,STAG,CANONICAL);
}



// Dummy routines for testing - will be removed.

int staggered_phase(const int *n, int mu){

//    int d, p;

    // MILC phases

    switch(mu){
    case 3:
	return 1;
    case 0:
	return n[3]%2?-1:1;
    case 1:
	return (n[0]+n[3])%2?-1:1;
    case 2:
	return (n[0]+n[1]+n[3])%2?-1:1;
    }

//     for(p=0, d=0; d<mu; d++) p += n[d];
//     return p%2?-1:1;
    
}

static int called=0;
void Fasqtad::parallel_transport(ParTransAsqtad &pt,Vector **vin, int *dir,
				     int N, Vector **vout){

    Matrix u, ud;
    int n, s[4], spm[4];
    int pos[] = {0,0,0,0};
FILE *fp;
	called++;
	char filename[200];
#if 0
	sprintf(filename,"%s%d",par_trans_filename,called);
	fp = fopen(filename,"w");
	printf("filename = %s called=%d\n",filename,called);
#else
	fp = stdout;
#endif
//	fprintf(fp,"called=%d\n",called);

#if 0
    for(int d=0; d<N; d++){

	bool backwards = false;
	int mu = (int)dir[d];
	if(mu%2) backwards = true;
        mu /= 2;
	
	for(n=0, spm[3]=s[3]=0; s[3]<GJP.TnodeSites(); spm[3]++, s[3]++)
	    for(spm[2]=s[2]=0; s[2]<GJP.ZnodeSites(); spm[2]++, s[2]++)
		for(spm[1]=s[1]=0; s[1]<GJP.YnodeSites(); spm[1]++, s[1]++)
		    for(spm[0]=s[0]=0; s[0]<GJP.XnodeSites(); spm[0]++, s[0]++, n++) {
			if(backwards){
			    --s[mu];
			    spm[mu] = (spm[mu]-1+GJP.NodeSites(mu))%GJP.NodeSites(mu);
			}else spm[mu] = (spm[mu]+1)%GJP.NodeSites(mu);
			Matrix u = *(const_cast<Matrix*>(GetLink(s, mu)));

//			if(staggered_phase(s, mu)==-1) u *= -1.0;
			
			if(backwards){
			    ud.Dagger(u);
			    u = ud;	
			    ++s[mu];
			}
			
			uDotXEqual((IFloat*)&vout[d][n], (IFloat*)&u, (IFloat*)&vin[d][FsiteOffset(spm)]);
			spm[mu] = s[mu];
#if 0
	IFloat * vout_p = (IFloat *)(vin[d] + n);
	fprintf(fp,"vin[%d][%d]=",d,n);
	for(int k=0;k<6;k++) fprintf(fp,"  %0.8e",vout_p[k]);
	fprintf(fp,"\n");
	vout_p = (IFloat *)(vout[d] + n);
	fprintf(fp,"vout[%d][%d]=",d,n);
	for(int k=0;k<6;k++) fprintf(fp,"  %0.8e",vout_p[k]);
	fprintf(fp,"\n");
#endif
		    }
    }
#endif

#if 0
    int vol = GJP.TnodeSites()*
    		 GJP.XnodeSites()*
    		 GJP.YnodeSites()*
    		 GJP.ZnodeSites();
    Vector *vout2[N];
    for(int i=0;i<N;i++){
	vout2[i]= (Vector *)smalloc(vol*sizeof(Vector));
    }
#endif
    pt.run(N,vout,vin,dir);
#if 0
    for(int d=0; d<N; d++)
	for(n=0, spm[3]=s[3]=0; s[3]<GJP.TnodeSites(); spm[3]++, s[3]++)
	for(spm[2]=s[2]=0; s[2]<GJP.ZnodeSites(); spm[2]++, s[2]++)
	for(spm[1]=s[1]=0; s[1]<GJP.YnodeSites(); spm[1]++, s[1]++)
	for(spm[0]=s[0]=0; s[0]<GJP.XnodeSites(); spm[0]++, s[0]++, n++) {
		IFloat * vout_p = (IFloat *)(vin[d] + n);
		fprintf(fp,"vin[%d][%d]=",d,n);
		for(int k=0;k<6;k++) fprintf(fp,"  %0.8e",vout_p[k]);
		fprintf(fp,"\n");
		vout_p = (IFloat *)(vout[d] + n);
		IFloat * vout_p2 = (IFloat *)(vout[d] + n);
		fprintf(fp,"vout[%d][%d]=",d,n);
		for(int k=0;k<6;k++) fprintf(fp,"  %0.8e",vout_p2[k]);
		fprintf(fp,"\n");
////		fprintf(fp,"(vout2-vout)[%d][%d]=",d,n);
////		for(int k=0;k<6;k++) fprintf(fp,"  %0.8e",vout_p2[k]-vout_p[k]);
////		fprintf(fp,"\n");
	}
#endif
//   fclose(fp);
#if 0
    for(int i=0;i<N;i++){
	sfree(vout2[i]);
    }
#endif
}


CPS_END_NAMESPACE
