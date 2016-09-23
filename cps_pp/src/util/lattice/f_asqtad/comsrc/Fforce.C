#include<config.h>

/*!\file
  \brief  Implementation of Fasqtad::EvolveMomFforce.

  $Id: Fforce.C,v 1.10 2008/02/08 18:35:07 chulwoo Exp $
*/
//--------------------------------------------------------------------


#include <util/lattice.h>
#include <util/pt.h>
#include <util/gjp.h>
#include <util/time_cps.h>
#include <util/amalloc.h>
#include <stdio.h>
#include <stdlib.h>

CPS_START_NAMESPACE



#define PROFILE

// N.B. No optimising provision is made if any of the asqtad coefficients
// are zero.

ForceArg Fasqtad::EvolveMomFforce(Matrix *mom, Vector *frm, Float mass, Float dt){

    char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
    VRB.Func(cname,fname);

#ifdef PROFILE    
    ParTrans::PTflops = 0;
    ForceFlops = 0;
    Float dtime;
#endif

    size_t size;

    size = GJP.VolNodeSites()/2*FsiteSize()*sizeof(Float);

    Vector *X = (Vector *)smalloc(2*size);
    Vector *X_e = X;                             // even sites
    Vector *X_o = X+GJP.VolNodeSites()/2;  // odd sites

    // The argument frm should have the CG solution.
    // The FstagTypes protected pointer f_tmp should contain Dslash frm

    #define TESTING
    #undef TESTING
    #ifndef TESTING
    moveMem(X_e, frm, size);
#undef DEBUGGING
#ifdef DEBUGGING
    Vector * f_tmp_save = f_tmp;      // FstagTypes destructor needs this.
    f_tmp = frm+GJP.VolNodeSites()/2; // Debugging only.
#endif
    moveMem(X_o, f_tmp, size);
    #else
    moveMem(X, frm, 2*size);
    #endif
    Fconvert(X, CANONICAL, STAG);

    Convert(STAG);  // Puts staggered phases into gauge field.
    
    const int VAXPY_UNROLL = 6; // For optimised linear algebra functions.
    const int vax_len = GJP.VolNodeSites()*VECT_LEN/VAXPY_UNROLL;
    
    // The number of directions to do simultaneously; N can be 1, 2 or 4.
    int N = 4;  
    if (GJP.VolNodeSites()>512) N = 1;       // This is to fit the vectors 
    else if (GJP.VolNodeSites()>256) N = 2;  // into EDRAM on QCDOC.

    VRB.Flow(cname,fname,"N=%d\n",N);

    enum{plus=0, minus=1, n_sign=2};

    // Array in which to accumulate the force term:
    // this must be initialised to zero 

    Matrix **force = (Matrix**)amalloc(fmalloc, sizeof(Matrix), 2, 4, GJP.VolNodeSites());
    for(int i=0; i<4; i++)
	for(int s=0; s<GJP.VolNodeSites(); s++) force[i][s].ZeroMatrix();


    // Vector arrays for which we must allocate memory

    Vector ***Pnu = (Vector***)amalloc(fmalloc, sizeof(Vector), 3, n_sign, N, GJP.VolNodeSites());
    
    Vector ****P3 = (Vector****)amalloc(fmalloc, sizeof(Vector), 4, n_sign, n_sign, N, GJP.VolNodeSites());

    Vector ****Prhonu = (Vector****)amalloc(fmalloc, sizeof(Vector), 4, n_sign, n_sign, N, GJP.VolNodeSites());

    Vector *****P5 = (Vector*****)amalloc(fmalloc, sizeof(Vector), 5, n_sign, n_sign, n_sign, N, GJP.VolNodeSites());

    Vector ******P7 = (Vector******)amalloc(fmalloc, sizeof(Vector), 6, n_sign, n_sign, n_sign, n_sign, N, GJP.VolNodeSites());

    Vector ******Psigma7 = (Vector******)amalloc(fmalloc, sizeof(Vector), 6, n_sign, n_sign, n_sign, n_sign, N, GJP.VolNodeSites());


    
    // These vectors can be overlapped with previously allocated memory
    
    Vector **Pnununu = Prhonu[0][0];
    Vector ***Pnunu = Psigma7[0][0][0];;
    Vector ****Pnu5 = P7[0][0];
    Vector ****Pnu3 = P7[0][0];
    Vector *****Prho5 = Psigma7[0];
    Vector *****Psigmarhonu = Psigma7[0];


    // input/output arrays for the parallel transport routines
    Vector *vin[n_sign*4], *vout[n_sign*4];
    int dir[n_sign*4];
	
    int mu[4], nu[4], rho[4], sigma[4];   // Sets of directions
    int w;                                // The direction index 0...N-1
    int ms, ns, rs, ss;                   // Sign of direction
    bool done[4] = {false,false,false,false};  // Flags to tell us which 
                                               // nu directions we have done.
    ParTransAsqtad parallel_transport(*this);
	    
    dtime = -dclock();
    for (int m=0; m<4; m+=N){                     	    // Loop over mu
	for(w=0; w<N; w++) mu[w] = (m+w)%4; 

	for (int n=m+1; n<m+4; n++){                        // Loop over nu
	    for(w=0; w<N; w++) nu[w] = (n+w)%4;

	    // Pnu = U_nu X

	    for(int i=0; i<N; i++){
		vin[i] = vin[i+N] = X;
		dir[n_sign*i] = n_sign*nu[i]+plus;        // nu_i
		dir[n_sign*i+1] = n_sign*nu[i]+minus;    // -nu_i
		vout[n_sign*i] = Pnu[minus][i];
		vout[n_sign*i+1] = Pnu[plus][i];
	    }
	    parallel_transport.run(n_sign*N, vout, vin, dir);

	    // P3 = U_mu Pnu
	    // ms is the nu sign index, ms is the mu sign index,
	    // w is the direction index
	    for(int i=0; i<N; i++){
		dir[n_sign*i] = n_sign*mu[i]+plus;        // mu_i
		dir[n_sign*i+1] = n_sign*mu[i]+minus;    // -mu_i
	    }
	    for(ns=0; ns<n_sign; ns++){               // ns is the sign of nu
		for(int i=0; i<N; i++){
		    vin[n_sign*i] = vin[n_sign*i+1] = Pnu[ns][i];
		    vout[n_sign*i] = P3[plus][ns][i];
		    vout[n_sign*i+1] = P3[minus][ns][i];
		}
		parallel_transport.run(n_sign*N, vout, vin, dir);
	    }
	    
	    for(w=0; w<N; w++)
		for(ns=0; ns<n_sign; ns++){
		    force_product_sum(P3[plus][ns][w], Pnu[ns][w],
				      GJP.staple3_coeff(),
				      force[mu[w]]);
		}

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
		    dir[n_sign*i] = n_sign*rho[i]+plus;        
		    dir[n_sign*i+1] = n_sign*rho[i]+minus;    
		}
		for(ns=0; ns<n_sign; ns++){
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = vin[n_sign*i+1] = Pnu[ns][i];
			vout[n_sign*i] = Prhonu[ns][minus][i];
			vout[n_sign*i+1] = Prhonu[ns][plus][i];
		    }
		    parallel_transport.run(n_sign*N, vout, vin, dir);
		}

		// P5 = U_mu Prhonu

		for(int i=0; i<N; i++){
		    dir[n_sign*i] = n_sign*mu[i]+plus;        
		    dir[n_sign*i+1] = n_sign*mu[i]+minus;    
		}
		for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) {
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = vin[n_sign*i+1] = Prhonu[ns][rs][i];
			vout[n_sign*i] = P5[plus][ns][rs][i];
			vout[n_sign*i+1] = P5[minus][ns][rs][i];
		    }
		    parallel_transport.run(n_sign*N, vout, vin, dir);
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
		    parallel_transport.run(n_sign*N, vout, vin, dir);
		}

		// P7 = U_mu P_sigmarhonu
		for(int i=0; i<N; i++){
		    dir[n_sign*i] = n_sign*mu[i]+plus;        
		    dir[n_sign*i+1] = n_sign*mu[i]+minus;    
		}
		for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) for(ss=0; ss<n_sign; ss++){
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = vin[n_sign*i+1] = Psigmarhonu[ns][rs][ss][i];
			vout[n_sign*i] = P7[plus][ns][rs][ss][i];
			vout[n_sign*i+1] = P7[minus][ns][rs][ss][i];
		    }
		    parallel_transport.run(n_sign*N, vout, vin, dir);
		}

		// F_mu -= P7 Psigmarhonu^\dagger
		
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) for(ss=0; ss<n_sign; ss++)
			force_product_sum(P7[plus][ns][rs][ss][w],
					  Psigmarhonu[ns][rs][ss][w],
					  GJP.staple7_coeff(),
					  force[mu[w]]);

		// F_sigma += P7 Psigmarhonu^\dagger
		// N.B. this is the same as one of the previous products.
		
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
			force_product_sum(P7[plus][ns][rs][minus][w],
					  Psigmarhonu[ns][rs][minus][w],
					  -GJP.staple7_coeff(),
					  force[sigma[w]]);

		// F_sigma += Psigmarhonu P7^\dagger
		
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
			force_product_sum(Psigmarhonu[ns][rs][minus][w],
					  P7[minus][ns][rs][minus][w],
					  -GJP.staple7_coeff(),
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
		    parallel_transport.run(n_sign*N, vout, vin, dir);
		}

		// F_sigma += Fsigma7 Frhonu^\dagger

	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
			force_product_sum(Psigma7[plus][ns][rs][plus][w],
					  Prhonu[ns][rs][w],
					  -GJP.staple7_coeff(),
					  force[sigma[w]]);

		// F_sigma += Frhonu Fsigma7^\dagger

	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
			force_product_sum(Prhonu[ns][rs][w],
					  Psigma7[minus][ns][rs][plus][w],
					  -GJP.staple7_coeff(),
					  force[sigma[w]]);

		// P5 += c_7/c_5 Psigma7

		if(GJP.staple5_coeff()!=0.0){
		    Float c75 = -GJP.staple7_coeff()/GJP.staple5_coeff();
		    for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) for(ss=0; ss<n_sign; ss++) for(w=0; w<N; w++)
#if TARGET==QCDOC
			
			vaxpy3(P5[ms][ns][rs][w],&c75, Psigma7[ms][ns][rs][ss][w], P5[ms][ns][rs][w], vax_len);
#else
			
			P5[ms][ns][rs][w]->FTimesV1PlusV2(c75, Psigma7[ms][ns][rs][ss][w], P5[ms][ns][rs][w], GJP.VolNodeSites()*VECT_LEN);
#endif
#ifdef PROFILE
			ForceFlops += 2*GJP.VolNodeSites()*VECT_LEN*N*n_sign*n_sign*n_sign*n_sign;
#endif
		}
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
		    dir[n_sign*i] = n_sign*rho[i]+plus;        
		    dir[n_sign*i+1] = n_sign*rho[i]+minus;    
		}
		for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++){
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = P5[ms][ns][plus][i];
			vin[n_sign*i+1] = P5[ms][ns][minus][i];
			vout[n_sign*i] = Prho5[ms][ns][plus][i];
			vout[n_sign*i+1] = Prho5[ms][ns][minus][i];
		    }
		    parallel_transport.run(n_sign*N, vout, vin, dir);
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

		if(GJP.staple3_coeff()!=0.0){
		    Float c53 = -GJP.staple5_coeff()/GJP.staple3_coeff();
		    for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) for(w=0; w<N; w++)
#if TARGET==QCDOC
			vaxpy3(P3[ms][ns][w],&c53,Prho5[ms][ns][rs][w], P3[ms][ns][w], vax_len);
#else
			P3[ms][ns][w]->FTimesV1PlusV2(c53, Prho5[ms][ns][rs][w], P3[ms][ns][w], GJP.VolNodeSites()*VECT_LEN);
#endif
#ifdef PROFILE
		ForceFlops += 2*GJP.VolNodeSites()*VECT_LEN*N*n_sign*n_sign*n_sign;
#endif
		}
		
	    } // rho+sigma loop

	    // Pnunu = U_nu Pnu

	    for(int i=0; i<N; i++){
		dir[n_sign*i] = n_sign*nu[i]+plus;        
		dir[n_sign*i+1] = n_sign*nu[i]+minus;    
	    }
	    for(int i=0; i<N; i++){
		vin[n_sign*i] = Pnu[minus][i];
		vin[n_sign*i+1] = Pnu[plus][i];
		vout[n_sign*i] = Pnunu[minus][i];
		vout[n_sign*i+1] = Pnunu[plus][i];
	    }
	    parallel_transport.run(n_sign*N, vout, vin, dir);

	    // P5 = U_mu Pnunu

	    for(int i=0; i<N; i++){
		dir[n_sign*i] = n_sign*mu[i]+plus;        
		dir[n_sign*i+1] = n_sign*mu[i]+minus;    
	    }
	    for(ns=0; ns<n_sign; ns++){
		for(int i=0; i<N; i++){
		    vin[n_sign*i] = Pnunu[ns][i];
		    vin[n_sign*i+1] = Pnunu[ns][i];
		    vout[n_sign*i] = P5[plus][ns][0][i];
		    vout[n_sign*i+1] = P5[minus][ns][0][i];
		}
		parallel_transport.run(n_sign*N, vout, vin, dir);
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
		dir[n_sign*i] = n_sign*nu[i]+plus;        
		dir[n_sign*i+1] = n_sign*nu[i]+minus;    
	    }
	    for(ms=0; ms<n_sign; ms++){
		for(int i=0; i<N; i++){
		    vin[n_sign*i] =   P5[ms][plus][0][i]; 
		    vin[n_sign*i+1] = P5[ms][minus][0][i];
		    vout[n_sign*i] =   Pnu5[ms][plus][i];
		    vout[n_sign*i+1] = Pnu5[ms][minus][i];
		}
		parallel_transport.run(n_sign*N, vout, vin, dir);
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

	    if(GJP.staple3_coeff()!=0.0){
		Float cl3 = -GJP.Lepage_coeff()/GJP.staple3_coeff();
		for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++) for(w=0; w<N; w++)
#if TARGET==QCDOC
		    vaxpy3(P3[ms][ns][w], &cl3, Pnu5[ms][ns][w], P3[ms][ns][w], vax_len);
#else
		    
		    P3[ms][ns][w]->FTimesV1PlusV2(cl3, Pnu5[ms][ns][w], P3[ms][ns][w], GJP.VolNodeSites()*VECT_LEN);
#endif
#ifdef PROFILE
		ForceFlops += 2*GJP.VolNodeSites()*VECT_LEN*N*n_sign*n_sign;
#endif
	    }
	    // F_nu += P3 Pnu^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(P3[plus][minus][w],
				  Pnu[minus][w],
				  -GJP.staple3_coeff(),
				  force[nu[w]]);

	    // F_nu +=  Pnu P3^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnu[minus][w],
				  P3[minus][minus][w],
				  -GJP.staple3_coeff(),
				  force[nu[w]]);
	    
	    // Pnu3 = U_nu P3

	    for(int i=0; i<N; i++)
		dir[i] = n_sign*nu[i]+plus;        
	    for(ms=0; ms<n_sign; ms++){
		for(int i=0; i<N; i++){
		    vin[i] = P3[ms][plus][i]; 
		    vout[i] = Pnu3[ms][plus][i];
		}
		parallel_transport.run(N, vout, vin, dir);
	    }

	    // F_nu += Pnu3 X^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnu3[plus][plus][w], X,
				  -GJP.staple3_coeff(),
				  force[nu[w]]);

	    // F_nu += X Pnu3^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(X, Pnu3[minus][plus][w], 
				  -GJP.staple3_coeff(),
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

	    // F_nu += Pnunu Pnu^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnunu[minus][w], Pnu[plus][w],
				  -GJP.Naik_coeff(),
				  force[nu[w]]);

	    // F_nu += Pnu Pnunu^\dagger
	    
	    for(w=0; w<N; w++)
		force_product_sum(Pnu[minus][w], Pnunu[plus][w],
				  GJP.Naik_coeff(),
				  force[nu[w]]);

	    // Pnununu = U_nu Pnunu

	    for(int i=0; i<N; i++){
		dir[i] = n_sign*nu[i]+plus;        
		vin[i] = Pnunu[minus][i]; 
		vout[i] = Pnununu[i];
	    }
	    parallel_transport.run(N, vout, vin, dir);
	    
	    // F_nu += Pnununu X^\dagger
		
	    for(w=0; w<N; w++)
		force_product_sum(Pnununu[w], X,
				  GJP.Naik_coeff(),
				  force[nu[w]]);
		
	    

	} // nu loop
    } // mu loop


    // Now that we have computed the force, we can update the momenta

#ifdef PROFILE
    dtime += dclock();
    print_flops(cname, fname, ParTrans::PTflops+ForceFlops, dtime);
#endif

    ForceArg Fdt = update_momenta(force, dt, mom);
    /*printf("HELLO\n");
    IFloat temp;
    IFloat *tmp_pointer;
    for(int i = 0; i < 4; i++)
      {
	tmp_pointer = (IFloat *)force[i];
	for(int j = 0; j < GJP.VolNodeSites(); j++)
	  for(int k = 0;k < 18; k++)
	    {
	      temp = *(tmp_pointer + 18*j+k);
	      if(temp > 1e-8)
		printf("HMC force[%d][%d][%d] = %e\n",i,j,k,temp);
	    }
	    }*/
    // Tidy up

    sfree(Pnu);
    sfree(P3);
    sfree(Prhonu);
    sfree(P5);
    sfree(P7);
    sfree(Psigma7);
    sfree(X);
    sfree(force);

    Convert(CANONICAL);

#ifdef DEBUGGING
    f_tmp = f_tmp_save;
#endif

    return Fdt;
    
}


CPS_END_NAMESPACE
