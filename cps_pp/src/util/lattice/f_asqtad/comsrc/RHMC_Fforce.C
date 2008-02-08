#include<config.h>
#include<math.h>
#include<util/qcdio.h>

/*!\file
  \brief  Implementation of Fasqtad::RHMC_EvolveMomFforce.

  $Id: RHMC_Fforce.C,v 1.14 2008-02-08 18:35:07 chulwoo Exp $
*/
//--------------------------------------------------------------------


#include <util/lattice.h>
#include <util/pt.h>
#include <util/gjp.h>
#include <util/amalloc.h>
#include <comms/glb.h>
#include <util/time_cps.h>

CPS_START_NAMESPACE

  
/*!
  The momentum is evolved for a single molecular dynamics timestep
  using the force from the RHMC pseudofermion action.
  \param mom The momentum.
  \param sol The solutions from the inversions against the pseudofermion
  fields.
  \param degree The degree of the rational approximation.
  \param alpha The factors with which the solution vectors are weighted.
  \param mass Not used.
  \param dt The molecular dynamics timestep used in the numerical integration.
  \param sol_d The solution vectors multiplied by the D-slash matrix.
   
  \post \a mom is assigned the value of the momentum after the molecular
  dynamics update.
*/
// N.B. No optimising provision is made if any of the asqtad coefficients
// are zero.

#define PROFILE

ForceArg Fasqtad::RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
				       int isz, Float *alpha, Float mass, 
				       Float dt, Vector **sol_d, 
				       ForceMeasure force_measure){

  char *fname = "RHMC_EvolveMomFforce";
  VRB.Func(cname,fname);

  // If we want to measure the force then we must use the slow calculation
  if (force_measure == FORCE_MEASURE_YES) {
    Float L1 = 0.0;
    Float L2 = 0.0;
    Float Linf = 0.0;
      
    int g_size = GJP.VolNodeSites() * GsiteSize();
      
    Matrix *mom_tmp = 
      (Matrix*)smalloc(g_size*sizeof(Float),cname, fname, "mom_tmp");
    ((Vector*)mom_tmp) -> VecZero(g_size);

    char force_label[100];

    for (int i=0; i<degree; i++) {
      f_tmp -> CopyVec(sol_d[i], e_vsize);
      ForceArg Fdt = EvolveMomFforce(mom_tmp,sol[i],mass,dt*alpha[i]);
      sprintf(force_label, "Rational, mass = %e, pole = %d:", mass, i+isz);
      Fdt.print(dt, force_label);
    }

    for (int i=0; i<g_size/18; i++) {
      Float norm = (mom_tmp+i)->norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }

    glb_sum(&L1);
    glb_sum(&L2);
    glb_max(&Linf);

    L1 /= 4.0*GJP.VolSites();
    L2 /= 4.0*GJP.VolSites();
      
    fTimesV1PlusV2((IFloat*)mom, 1.0, (IFloat*)mom_tmp, (IFloat*)mom, g_size);
      
    sfree(mom_tmp, cname, fname, "mom_tmp");
      
    return ForceArg(L1, sqrt(L2), Linf);

  }

#ifdef PROFILE
  Float dtime;
  ParTrans::PTflops=0;
  ForceFlops=0;
#endif

  // The number of directions to do simultaneously - N can be 1, 2 or 4.
  const int N = 1;  
  enum{plus, minus, n_sign};

  const int vol = GJP.VolNodeSites();
  const int f_size = vol*FsiteSize()/2;
    
  for (int i=0; i<degree; i++) {
    sol[i]->VecTimesEquFloat(alpha[i], f_size);
  }

  // X_odd = D X_even    
  for (int p=0; p<degree; p++) {
    Fconvert(sol[p], CANONICAL, STAG);
  }
  Convert(STAG);  // Puts staggered phases into gauge field.

  // These are work Matrices used in force_product_sum routines
  Matrix *mtmp = (Matrix*)fmalloc(vol*sizeof(Matrix), cname, fname, "mtmp");

  /*
  // Matrix fields for which we must allocate memory
  Matrix ***Pnu = (Matrix***)amalloc(fmalloc, sizeof(Matrix), 3, n_sign, N, vol);
  Matrix ****P3 = (Matrix****)amalloc(fmalloc, sizeof(Matrix), 4, n_sign, n_sign, N, vol);
  Matrix ****Prhonu = (Matrix****)amalloc(fmalloc, sizeof(Matrix), 4, n_sign, n_sign, N, vol);
  Matrix *****P5 = (Matrix*****)amalloc(fmalloc, sizeof(Matrix), 5, n_sign, n_sign, n_sign, N, vol);
  Matrix ******P7 = (Matrix******)amalloc(fmalloc, sizeof(Matrix), 6, n_sign, n_sign, n_sign, n_sign,N,vol);
  Matrix ******Psigma7 = (Matrix******)amalloc(fmalloc,sizeof(Matrix),6,n_sign,n_sign,n_sign,n_sign,N,vol);
  Matrix **L = (Matrix**)amalloc(fmalloc, sizeof(Matrix), 2, N, vol);
  Matrix ***Lnu = (Matrix***)amalloc(fmalloc, sizeof(Matrix), 3, n_sign, N, vol);
  Matrix ****Lrhonu = (Matrix****)amalloc(fmalloc, sizeof(Matrix), 4, n_sign, n_sign, N, vol);
  Matrix ****Lmusigmarhonu = (Matrix****)amalloc(fmalloc, sizeof(Matrix), 4, n_sign, n_sign, N, vol);

  */

  Matrix ***Pnu = (Matrix***)smalloc(n_sign*sizeof(Matrix**), cname, fname, "Pnu");
  Matrix ****P3 = (Matrix****)smalloc(n_sign*sizeof(Matrix***), cname, fname, "P3");
  Matrix ****Prhonu = (Matrix****)smalloc(n_sign*sizeof(Matrix***), cname, fname, "Prhonu");
  Matrix *****P5 = (Matrix*****)smalloc(n_sign*sizeof(Matrix****), cname, fname, "P5");
  Matrix ******P7 = (Matrix******)smalloc(n_sign*sizeof(Matrix*****), cname, fname, "P7");
  Matrix ******Psigma7 = 
    (Matrix******)smalloc(n_sign*sizeof(Matrix*****), cname, fname, "Psigma7");
  Matrix ***Lnu = (Matrix***)smalloc(n_sign*sizeof(Matrix**), cname, fname, "Lnu");
  Matrix ****Lrhonu = 
    (Matrix****)smalloc(n_sign*sizeof(Matrix***), cname, fname, "Lrhonu");
  Matrix ****Lmusigmarhonu = 
    (Matrix****)smalloc(n_sign*sizeof(Matrix***), cname, fname, "Lmusigmarhonu");

  Matrix **L = (Matrix**)smalloc(N*sizeof(Matrix*), cname, fname, "L");
  for (int i=0; i<N; i++)  L[i] = (Matrix*)smalloc(vol*sizeof(Matrix), cname, fname, "L[i]");

  for (int i=0; i<n_sign; i++) {
    Pnu[i] = (Matrix**)smalloc(N*sizeof(Matrix*), cname, fname, "Pnu[i]");
    Lnu[i] = (Matrix**)smalloc(N*sizeof(Matrix*), cname, fname, "Lnu[i]");
    for (int j=0; j<N; j++) {
      Pnu[i][j] = (Matrix*)smalloc(vol*sizeof(Matrix), cname, fname, "Pnu[i][j]");
      Lnu[i][j] = (Matrix*)smalloc(vol*sizeof(Matrix), cname, fname, "Lnu[i][j]");
    }
    P3[i] = (Matrix***)smalloc(n_sign*sizeof(Matrix**), cname, fname, "P3[i]");
    Prhonu[i] = (Matrix***)smalloc(n_sign*sizeof(Matrix**), cname, fname, "Prhonu[i]");
    Lrhonu[i] = (Matrix***)smalloc(n_sign*sizeof(Matrix**), cname, fname, "Lrhonu[i]");
    Lmusigmarhonu[i] = (Matrix***)smalloc(n_sign*sizeof(Matrix**), 
					  cname, fname, "Lmusigmarhonu[i]");
    P5[i] = (Matrix****)smalloc(n_sign*sizeof(Matrix***), cname, fname, "P5[i]");
    P7[i] = (Matrix*****)smalloc(n_sign*sizeof(Matrix****), cname, fname, "P7[i]");
    Psigma7[i] = (Matrix*****)smalloc(n_sign*sizeof(Matrix****), cname, fname, "Psigma7[i]");
    for (int j=0; j<n_sign; j++) {
      P3[i][j] = (Matrix**)smalloc(N*sizeof(Matrix*), cname, fname, "P3[i][j]");
      Prhonu[i][j] = (Matrix**)smalloc(N*sizeof(Matrix*), cname, fname, "Prhonu[i][j]");
      Lrhonu[i][j] = (Matrix**)smalloc(N*sizeof(Matrix*), cname, fname, "Lrhonu[i][j]");
      Lmusigmarhonu[i][j] = (Matrix**)smalloc(N*sizeof(Matrix*), 
					      cname, fname, "Lmusigmarhonu[i][j]");
      for (int k=0; k<N; k++) {
	P3[i][j][k] = (Matrix*)smalloc(vol*sizeof(Matrix), cname, fname, "P3[i][j][k]");
	Prhonu[i][j][k] = (Matrix*)smalloc(vol*sizeof(Matrix), "Prhonu[i][j][k]");
	Lrhonu[i][j][k] = (Matrix*)smalloc(vol*sizeof(Matrix), "Lrhonu[i][j][k]");
	Lmusigmarhonu[i][j][k] = 
	  (Matrix*)smalloc(vol*sizeof(Matrix), "Lmusigmarhonu[i][j][k]");
      }
      P5[i][j] = (Matrix***)smalloc(n_sign*sizeof(Matrix**), cname, fname, "P5[i][j]");
      P7[i][j] = (Matrix****)smalloc(n_sign*sizeof(Matrix***), cname, fname, "P7[i][j]");
      Psigma7[i][j] = (Matrix****)
	smalloc(n_sign*sizeof(Matrix***), cname, fname, "Psigma7[i][j]");
      for (int k=0; k<n_sign; k++) {
	P5[i][j][k] = (Matrix**)smalloc(N*sizeof(Matrix*), cname, fname, "P5[i][j][k]");
	for (int l=0; l<N; l++) {
	  P5[i][j][k][l] = (Matrix*)smalloc(vol*sizeof(Matrix), cname, 
					    fname, "P5[i][j][k][l]");
	}
	P7[i][j][k] = (Matrix***)smalloc(n_sign*sizeof(Matrix**), cname, 
					 fname, "P7[i][j][k]");
	Psigma7[i][j][k] = (Matrix***)smalloc(n_sign*sizeof(Matrix**), cname, 
					      fname, "Psigma7[i][j][k]");
	for (int l=0; l<n_sign; l++) {
	  P7[i][j][k][l] = (Matrix**)smalloc(N*sizeof(Matrix*), cname, 
					     fname, "P7[i][j][k][l]");
	  Psigma7[i][j][k][l] = (Matrix**)smalloc(N*sizeof(Matrix*), cname, fname, 
						  "P7[i][j][k][l]");
	  for (int m=0; m<N; m++) {
	    P7[i][j][k][l][m] = (Matrix*)smalloc(vol*sizeof(Matrix), cname, fname, 
						 "P7[i][j][k][l][m]");
	    Psigma7[i][j][k][l][m] = (Matrix*)smalloc(vol*sizeof(Matrix), cname, fname, 
						      "Psigma7[i][j][k][l][m]");
	  }
	}
      }
    }
  }


  // Array in which to accumulate the force term      
  //Matrix **force = (Matrix**)amalloc(fmalloc, sizeof(Matrix), 2, 4, vol);

  Matrix **force = (Matrix**)smalloc(4*sizeof(Matrix*), cname, fname, "force");
  for (int i=0; i<4; i++) {
    force[i] = (Matrix*) fmalloc(vol*sizeof(Matrix), cname, fname, "force[i]");
  }

  ParTransAsqtad parallel_transport(*this);

  // These fields can be overlapped with previously allocated memory
    
  Matrix **Pnununu = Prhonu[0][0];
  Matrix ***Pnunu = Psigma7[0][0][0];
  Matrix ****Pnu5 = P7[0][0];
  Matrix ****Pnu3 = P7[0][0];
  Matrix *****Prho5 = P7[0];
  Matrix *****Psigmarhonu = Psigma7[0];

  Matrix *****Lsigmarhonu = Psigma7[1];
  Matrix ****Lmurhonu = P7[0][0];    
  Matrix ***Lmunu = P7[1][0][0];
  Matrix ***Lnunu = P7[0][0][1];
  Matrix **Lmununu = Psigma7[1][0][0][0];
  Matrix **Lmu = Psigma7[1][0][0][0];	
  Matrix **L3 = Psigma7[1][0][0][0];	
  Matrix **L3nu = Psigma7[1][0][0][1];
  Matrix **L3nunu = Psigma7[1][0][0][0];	


  for(int i=0; i<4; i++)
    for(int s=0; s<vol; s++) force[i][s].ZeroMatrix();

  // input/output arrays for the parallel transport routines
  Matrix *vin[n_sign*N], *vout[n_sign*N];
  int dir[n_sign*N];

	
  int mu[N], nu[N], rho[N], sigma[N];   // Sets of directions
  int w;                                // The direction index 0...N-1
  int ms, ns, rs, ss;                   // Sign of direction
  bool done[4] = {false,false,false,false};  // Flags to tell us which 

#ifdef PROFILE
  dtime = -dclock();
#endif  

  // nu directions we have done.
  for (int m=0; m<4; m+=N){                     	    // Loop over mu
    for(w=0; w<N; w++) mu[w] = (m+w)%4;

    // L = sum X[x]X[x+mu]^dagger
    parallel_transport.vvpd(sol, degree, mu, N, 1, L);

    // F_mu += U L^dagger
    for(w=0; w<N; w++)
      force_product_sum(L[w], mu[w], GJP.KS_coeff(), force[mu[w]]);

    for (int n=m+1; n<m+4; n++){                        // Loop over nu

      for(w=0; w<N; w++) nu[w] = (n+w)%4;
	
      // Pnu_[+/-nu] = U_[-/+nu] 
	
      for(int i=0; i<N; i++){ 
	dir[n_sign*i] = n_sign*nu[i]+plus;        // mu_i
	dir[n_sign*i+1] = n_sign*nu[i]+minus;    // -mu_i
	vout[n_sign*i] = Pnu[minus][i];
	vout[n_sign*i+1] = Pnu[plus][i];
      }

      parallel_transport.shift_link(vout, dir, n_sign*N);

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

      // Lnu[x] = L[x-/+nu]
      for(int i=0; i<N; i++){
	dir[n_sign*i] = n_sign*nu[i]+plus;       // +nu_i
	dir[n_sign*i+1] = n_sign*nu[i]+minus;    // -nu_i
	  
	vin[n_sign*i] = vin[n_sign*i+1] = L[i];
	vout[n_sign*i] = Lnu[minus][i];
	vout[n_sign*i+1] = Lnu[plus][i];
      }

      // L = sum X[x]X[x+mu]^dagger
      parallel_transport.shift_field(vin, dir, n_sign*N, 1, vout);

      // F += P3 (Pnu Lnu)^dagger
      for(w=0; w<N; w++)
	for(ns=0; ns<n_sign; ns++)
	  force_product_sum(P3[plus][ns][w], Pnu[ns][w], Lnu[ns][w],
			    GJP.staple3_coeff(), force[mu[w]], mtmp);


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
	  
	// Lrhonu[x] = Lnu[x+/-rho]
	  
	for(int i=0; i<N; i++){
	  dir[n_sign*i] = n_sign*rho[i]+plus;        
	  dir[n_sign*i+1] = n_sign*rho[i]+minus;
	}
	for(ns=0; ns<n_sign; ns++){
	  for(int i=0; i<N; i++){
	    vin[n_sign*i] = vin[n_sign*i+1] = Lnu[ns][i];
	    vout[n_sign*i] = Lrhonu[ns][minus][i];
	    vout[n_sign*i+1] = Lrhonu[ns][plus][i];
	  }
	  parallel_transport.shift_field(vin, dir, n_sign*N, 1, vout);
	}
	  
	// F_mu += P5 (Prhonu Lrhonu)^dagger
	  
	for(w=0; w<N; w++)
	  for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++)
	    force_product_sum(P5[plus][ns][rs][w], Prhonu[ns][rs][w],
			      Lrhonu[ns][rs][w], GJP.staple5_coeff(),
			      force[mu[w]], mtmp);
	  
	// Psigmarhonu = U_sigma P_rhonu
	  
	for(int i=0; i<N; i++){
	  dir[n_sign*i] = n_sign*sigma[i]+plus;        
	  dir[n_sign*i+1] = n_sign*sigma[i]+minus;    
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
	  
	// Lsigmarhonu = Lrhonu[x+/-sigma]
	  
	for(int i=0; i<N; i++){
	  dir[n_sign*i] = n_sign*sigma[i]+plus;        
	  dir[n_sign*i+1] = n_sign*sigma[i]+minus;    
	}   
	for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++){
	  for(int i=0; i<N; i++){
	    vin[n_sign*i] = vin[n_sign*i+1] = Lrhonu[ns][rs][i];
	    vout[n_sign*i] = Lsigmarhonu[ns][rs][minus][i];
	    vout[n_sign*i+1] = Lsigmarhonu[ns][rs][plus][i];
	  }
	  parallel_transport.shift_field(vin, dir, n_sign*N, 1, vout);
	}
	  
	  
	// F_mu -= P7 (Psigmarhonu Lsigmarhonu)^\dagger
	  
	for(w=0; w<N; w++)
	  for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) for(ss=0; ss<n_sign; ss++)
	    force_product_sum(P7[plus][ns][rs][ss][w], Psigmarhonu[ns][rs][ss][w],
			      Lsigmarhonu[ns][rs][ss][w], GJP.staple7_coeff(),
			      force[mu[w]], mtmp);
	
	// F_sigma += P7 (Psigmarhonu Lsigmarhonu)^\dagger
	// N.B. this is the same as one of the previous products.
	  
	for(w=0; w<N; w++)
	  for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
	    force_product_sum(P7[plus][ns][rs][minus][w], Psigmarhonu[ns][rs][minus][w],
			      Lsigmarhonu[ns][rs][minus][w], -GJP.staple7_coeff(),
			      force[sigma[w]], mtmp);

	// Lmusigmarhonu = Lsigmarhonu[x-mu]
	  
	for(int i=0; i<N; i++) dir[i] = n_sign*mu[i]+minus;
	for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++){
	  for(int i=0; i<N; i++){
	    vin[i] = Lsigmarhonu[ns][rs][minus][i];
	    vout[i] = Lmusigmarhonu[ns][rs][i];
	  }
	  parallel_transport.shift_field(vin, dir, N, 1, vout);
	}
	  
	  
	// F_sigma += Psigmarhonu (P7 Lmusigmarhonu)^\dagger
	  
	for(w=0; w<N; w++)
	  for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
	    force_product_sum(Psigmarhonu[ns][rs][minus][w], P7[minus][ns][rs][minus][w],
			      Lmusigmarhonu[ns][rs][w], -GJP.staple7_coeff(),
			      force[sigma[w]], mtmp);

	// Psigma7 = U_sigma P7 
	for(int i=0; i<N; i++){
	  dir[n_sign*i] = n_sign*sigma[i]+plus;        
	  dir[n_sign*i+1] = n_sign*sigma[i]+minus;    
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
	  
	// F_sigma += Fsigma7 (Frhonu Lrhonu)^\dagger
	  
	for(w=0; w<N; w++)
	  for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
	    force_product_sum(Psigma7[plus][ns][rs][plus][w],	Prhonu[ns][rs][w],
			      Lrhonu[ns][rs][w], -GJP.staple7_coeff(),
			      force[sigma[w]], mtmp);

	// Lmurhonu = Lrhonu[x-mu]
		  
	for(int i=0; i<N; i++) dir[i] = n_sign*mu[i]+minus;
	for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++){
	  for(int i=0; i<N; i++){
	    vin[i] = Lrhonu[ns][rs][i];
	    vout[i] = Lmurhonu[ns][rs][i];
	  }
	  parallel_transport.shift_field(vin, dir, N, 1, vout);
	}
	  
	// F_sigma += Frhonu (Fsigma7 Lmurhonu)^\dagger
	  
	for(w=0; w<N; w++)
	  for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
	    force_product_sum(Prhonu[ns][rs][w], Psigma7[minus][ns][rs][plus][w],
			      Lmurhonu[ns][rs][w], -GJP.staple7_coeff(),
			      force[sigma[w]], mtmp);
	  
	if(GJP.staple5_coeff()!=0.0) {
	  Float tmp = -GJP.staple7_coeff()/GJP.staple5_coeff();
	  for(ms=0; ms<n_sign; ms++) 
	    for(ns=0; ns<n_sign; ns++) 
	      for(rs=0; rs<n_sign; rs++) 
		for(ss=0; ss<n_sign; ss++) 
		  for(w=0; w<N; w++)
		    {
#if TARGET == QCDOC
		      vaxpy3_m(P5[ms][ns][rs][w],
			       &tmp,
			       Psigma7[ms][ns][rs][ss][w],
			       P5[ms][ns][rs][w],
			       vol*MATRIX_SIZE/6);
#else
		      fTimesV1PlusV2((IFloat*)P5[ms][ns][rs][w],
				     tmp,
				     (IFloat*)Psigma7[ms][ns][rs][ss][w],
				     (IFloat*)P5[ms][ns][rs][w],
				     vol*MATRIX_SIZE);
#endif
#ifdef PROFILE
		      ForceFlops += 2*vol*MATRIX_SIZE;
#endif
		    }
	}
	  
	// F_rho -= P5 (Prhonu Lrhonu)^\dagger
	for(w=0; w<N; w++)
	  for(ns=0; ns<n_sign; ns++)
	    force_product_sum(P5[plus][ns][minus][w], Prhonu[ns][minus][w],
			      Lrhonu[ns][minus][w], -GJP.staple5_coeff(),
			      force[rho[w]], mtmp);

	// F_rho -= Prhonu (P5 Lmurhono)^\dagger
	  
	for(w=0; w<N; w++)
	  for(ns=0; ns<n_sign; ns++)
	    force_product_sum(Prhonu[ns][minus][w], P5[minus][ns][minus][w],
			      Lmurhonu[ns][minus][w], -GJP.staple5_coeff(),
			      force[rho[w]], mtmp);

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
	  
	// F_rho -= Prho5 (Pnu Lnu)^\dagger
	  
	for(w=0; w<N; w++)
	  for(ns=0; ns<n_sign; ns++)
	    force_product_sum(Prho5[plus][ns][plus][w], Pnu[ns][w],
			      Lnu[ns][w], -GJP.staple5_coeff(),
			      force[rho[w]], mtmp);
	// Lmunu = Lnu[x-mu]
	  
	for(int i=0; i<N; i++) dir[i] = n_sign*mu[i]+minus;
	for(ns=0; ns<n_sign; ns++){
	  for(int i=0; i<N; i++){
	    vin[i] = Lnu[ns][i];
	    vout[i] = Lmunu[ns][i];
	  }
	  parallel_transport.shift_field(vin, dir, N, 1, vout);
	}
	  
	  
	// F_rho -= Pnu (Prho5 Lmunu)^\dagger
	  
	for(w=0; w<N; w++)
	  for(ns=0; ns<n_sign; ns++)
	    force_product_sum(Pnu[ns][w], Prho5[minus][ns][plus][w],
			      Lmunu[ns][w], -GJP.staple5_coeff(),
			      force[rho[w]], mtmp);

	// P3 += c_5/c_3 Prho5
	  
	if(GJP.staple3_coeff()!=0.0) {
	  Float tmp = -GJP.staple5_coeff()/GJP.staple3_coeff();
	  for(ms=0; ms<n_sign; ms++) 
	    for(ns=0; ns<n_sign; ns++) 
	      for(rs=0; rs<n_sign; rs++) 
		for(w=0; w<N; w++)
		  {
#if TARGET == QCDOC
		    vaxpy3_m(P3[ms][ns][w],
			     &tmp,
			     Prho5[ms][ns][rs][w],
			     P3[ms][ns][w],
			     vol*MATRIX_SIZE/6);
#else
		    fTimesV1PlusV2((IFloat*)P3[ms][ns][w],
				   tmp,
				   (IFloat*)Prho5[ms][ns][rs][w],
				   (IFloat*)P3[ms][ns][w],
				   vol*MATRIX_SIZE);
#endif
#ifdef PROFILE
		    ForceFlops += 2*vol*MATRIX_SIZE;
#endif
		  }
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
	  vin[n_sign*i] = vin[n_sign*i+1] = Pnunu[ns][i];
	  vout[n_sign*i] = P5[plus][ns][0][i];
	  vout[n_sign*i+1] = P5[minus][ns][0][i];
	}
	parallel_transport.run(n_sign*N, vout, vin, dir);
      }
	
      //Lnunu = Lnu[x+/-nu]
	
      for(int i=0; i<N; i++){
	dir[n_sign*i] = n_sign*nu[i]+minus;        
	dir[n_sign*i+1] = n_sign*nu[i]+plus;    
      }
      for(int i=0; i<N; i++){
	vin[n_sign*i] = Lnu[plus][i];
	vin[n_sign*i+1] = Lnu[minus][i];
	vout[n_sign*i] = Lnunu[plus][i];
	vout[n_sign*i+1] = Lnunu[minus][i];
      }

      parallel_transport.shift_field(vin, dir, n_sign*N, 1, vout);
	
      // F_mu += P5 (Pnunu Lnunu)^\dagger
	
      for(w=0; w<N; w++)
	for(ns=0; ns<n_sign; ns++)
	  force_product_sum(P5[plus][ns][0][w], Pnunu[ns][w], Lnunu[ns][w],
			    GJP.Lepage_coeff(), force[mu[w]], mtmp);

      // F_nu -= P5 (Pnunu Lnunu)^\dagger
      // N.B. this is the same as one of the previous products
	
      for(w=0; w<N; w++)
	force_product_sum(P5[plus][minus][0][w], Pnunu[minus][w],
			  Lnunu[minus][w], -GJP.Lepage_coeff(),
			  force[nu[w]], mtmp);

      // Lmununu = Lnunu[x-mu]
	
      for(int i=0; i<N; i++) dir[i] = n_sign*mu[i]+minus;
      for(int i=0; i<N; i++){
	vin[i] = Lnunu[minus][i];
	vout[i] = Lmununu[i];
      }

      parallel_transport.shift_field(vin, dir, N, 1, vout);
	
      // F_nu -= Pnunu (P5 Lmununu)^\dagger
	
      for(w=0; w<N; w++)
	force_product_sum(Pnunu[minus][w], P5[minus][minus][0][w],
			  Lmununu[w], -GJP.Lepage_coeff(),
			  force[nu[w]], mtmp);

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
	
      // F_nu -= Pnu5 (Pnu Lnu)^\dagger
	
      for(w=0; w<N; w++)
	force_product_sum(Pnu5[plus][plus][w], Pnu[plus][w], Lnu[plus][w],
			  -GJP.Lepage_coeff(), force[nu[w]], mtmp);

      // F_nu -= Pnu (Pnu5 Lmunu)^\dagger
	
      for(w=0; w<N; w++)
	force_product_sum(Pnu[plus][w], Pnu5[minus][plus][w],
			  Lmunu[plus][w], -GJP.Lepage_coeff(),
			  force[nu[w]], mtmp);

      // P3 += c_L/c_3 Pnu5
	
      if(GJP.staple3_coeff()!=0.0) {
	Float tmp = -GJP.Lepage_coeff()/GJP.staple3_coeff();
	for(ms=0; ms<n_sign; ms++) 
	  for(ns=0; ns<n_sign; ns++) 
	    for(w=0; w<N; w++)
	      {
#if TARGET == QCDOC
		vaxpy3_m(P3[ms][ns][w],
			 &tmp,
			 Pnu5[ms][ns][w],
			 P3[ms][ns][w],
			 vol*MATRIX_SIZE/6);
#else
		fTimesV1PlusV2((IFloat*)P3[ms][ns][w],
			       tmp,
			       (IFloat*)Pnu5[ms][ns][w],
			       (IFloat*)P3[ms][ns][w],
			       vol*MATRIX_SIZE);
#endif
#ifdef PROFILE
		ForceFlops += 2*vol*MATRIX_SIZE;
#endif			    
	      }
      }

      // F_nu += P3 (Pnu Lnu)^\dagger
	
      for(w=0; w<N; w++)
	force_product_sum(P3[plus][minus][w], Pnu[minus][w], Lnu[minus][w], 
			  -GJP.staple3_coeff(), force[nu[w]], mtmp);

      // F_nu +=  Pnu (P3 Lmunu)^\dagger
	
      for(w=0; w<N; w++)
	force_product_sum(Pnu[minus][w], P3[minus][minus][w], Lmunu[minus][w],
			  -GJP.staple3_coeff(), force[nu[w]], mtmp);

      // Pnu3 = U_nu P3
	
      for(int i=0; i<N; i++) dir[i] = n_sign*nu[i]+plus;        
      for(ms=0; ms<n_sign; ms++){
	for(int i=0; i<N; i++){
	  vin[i] = P3[ms][plus][i]; 
	  vout[i] = Pnu3[ms][plus][i];
	}
	parallel_transport.run(N, vout, vin, dir);
      }
	
      // F_nu += Pnu3 L^\dagger
	
      for(w=0; w<N; w++)
	force_product_sum(Pnu3[plus][plus][w], L[w], -GJP.staple3_coeff(),
			  force[nu[w]]);

      // Lmu = L[x-mu]
	
      for(int i=0; i<N; i++) dir[i] = n_sign*mu[i]+minus;
      for(int i=0; i<N; i++){
	vin[i] = L[i];
	vout[i] = Lmu[i];
      }

      parallel_transport.shift_field(vin, dir, N, 1, vout);
	
      // F_nu += (Pnu3 Lmu)^\dagger
	
      for(w=0; w<N; w++)
	force_product_d_sum(Pnu3[minus][plus][w], Lmu[w], -GJP.staple3_coeff(),
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

#if 1
      // L3 = sum X[x]X[x+3]^dagger
      parallel_transport.vvpd(sol, degree, nu, N, 3, L3);
#else
      Vector ***vx = (Vector ***)smalloc(sizeof(Vector**),cname,fname,"vx");
      vx[0] = (Vector**)smalloc(n_sign*N*sizeof(Vector*),cname,fname,"vx[0]");
      Matrix **mout = (Matrix**)smalloc(n_sign*N*sizeof(Matrix*),cname,fname,"mout");
      for(int deg = 0; deg < degree; deg++)
	{
	  for(int i = 0; i < N; i++)
	    {
	      dir[i] = n_sign*nu[i]+plus;
	      vx[0][i] = sol[deg];
	    }
	  if(deg == 0)
	    parallel_transport.vvpd(sol,vx,1,dir,N,3,L3,1);
	  else
	    parallel_transport.vvpd(sol+deg,vx,1,dir,N,3,L3,0);
	}
      sfree(vx[0], cname, fname, "vx[0]");
      sfree(vx, cname, fname, "vx");
      sfree(mout, cname, fname, "mout");
#endif

      // Pnununu = U_nu Pnunu

      for(int i=0; i<N; i++){
	dir[i] = n_sign*nu[i]+plus;        
	vin[i] = Pnunu[minus][i]; 
	vout[i] = Pnununu[i];
      }
      parallel_transport.run(N, vout, vin, dir);
	
      // F_nu += Pnununu L3^\dagger
	
      for(w=0; w<N; w++)
	{
	  force_product_sum(Pnununu[w], L3[w], GJP.Naik_coeff(), force[nu[w]]);
	  /*printf("mu[%d] = %d\n",w,mu[w]);
	    IFloat temp;
	    IFloat *tmp_pointer = (IFloat *)Pnu[minus][w];
	    for(int i = 0; i < 4; i++)
	    for(int j = 0; j < GJP.VolNodeSites(); j++)
	    for(int k = 0;k < 18; k++)
	    {
	    temp = *(tmp_pointer + 72*j+18*i+k);
	    if(temp > 1e-8)
	    printf("RHMC Pmumu[plus][plus][%d][%d][%d][%d] = %e\n",w,i,j,k,temp);
	    }*/
	}

      // L3nu = L3[x-nu]
	
      for(int i=0; i<N; i++){
	dir[i] = n_sign*nu[i]+minus;        
	vin[i] = L3[i]; 
	vout[i] = L3nu[i];
      }
      parallel_transport.shift_field(vin, dir, N, 1, vout);
	
      // F_nu += Pnunu (Pnu L3nu)^\dagger
	
      for(w=0; w<N; w++)
	force_product_sum(Pnunu[minus][w], Pnu[plus][w], L3nu[w],
			  -GJP.Naik_coeff(), force[nu[w]], mtmp);

      // L3nunu = L3nu[x-nu]
	
      for(int i=0; i<N; i++){
	dir[i] = n_sign*nu[i]+minus;        
	vin[i] = L3nu[i]; 
	vout[i] = L3nunu[i];
      }
      parallel_transport.shift_field(vin, dir, N, 1, vout);
	
      // F_nu += Pnu (Pnunu L3nunu)^\dagger
	
      for(w=0; w<N; w++)
	force_product_sum(Pnu[minus][w], Pnunu[plus][w], L3nunu[w],
			  GJP.Naik_coeff(), force[nu[w]], mtmp);

    } // nu loop
  } // mu loop

#ifdef PROFILE
  dtime += dclock();
  print_flops(cname, fname, ParTrans::PTflops+ForceFlops, dtime);
#endif
    
  // Now that we have computed the force, we can update the momenta
  update_momenta(force, dt, mom);

  // Free allocated memory
    

  for (int i=0; i<N; i++) sfree(L[i], cname, fname, "L[i]");
  sfree(L, cname, fname, "L");
    
  for (int i=0; i<n_sign; i++) {
    for (int j=0; j<n_sign; j++) {
      for (int k=0; k<N; k++) {
	sfree(P3[i][j][k], cname, fname, "P3[i][j][k]");
	sfree(Prhonu[i][j][k], cname, fname, "Prhonu[i][j][k]");
	sfree(Lrhonu[i][j][k], cname, fname, "Lrhonu[i][j][k]");
	sfree(Lmusigmarhonu[i][j][k], cname, fname, "Lmusigmarhonu[i][j][k]");
      }
      for (int k=0; k<n_sign; k++) {
	for (int l=0; l<n_sign; l++) {
	  for (int m=0; m<N; m++) {
	    sfree(P7[i][j][k][l][m], cname, fname, "P7[i][j][k][l][m]");
	    sfree(Psigma7[i][j][k][l][m], cname, fname, "Psigma7[i][j][k][l][m]");
	  }
	  sfree(P7[i][j][k][l], cname, fname, "P7[i][j][k][l]");
	  sfree(Psigma7[i][j][k][l], cname, fname, "Psigma7[i][j][k][l]");
	}
	for (int l=0; l<N; l++) sfree(P5[i][j][k][l], cname, fname, "P5[i][j][k][l]");
	sfree(P7[i][j][k], cname, fname, "P7[i][j][k]");
	sfree(Psigma7[i][j][k], cname, fname, "Psigma[i][j][k]");
	sfree(P5[i][j][k], cname, fname, "P5[i][j][k]");
      }
      sfree(P3[i][j], cname, fname, "P3[i][j]");
      sfree(Prhonu[i][j], cname, fname, "Prhonu[i][j]");
      sfree(Lrhonu[i][j], cname, fname, "Lrhonu[i][j]");
      sfree(Lmusigmarhonu[i][j], cname, fname, "Lmusigmarhonu[i][j]");
      sfree(P5[i][j], cname, fname, "P5[i][j]");
      sfree(P7[i][j], cname, fname, "P7[i][j]");
      sfree(Psigma7[i][j], cname, fname, "Psigma[i][j]");
    }
    for (int j=0; j<N; j++) {
      sfree(Pnu[i][j], cname, fname, "Pnu[i][j]");
      sfree(Lnu[i][j], cname, fname, "Lnu[i][j]");
    }
    sfree(Pnu[i], cname, fname, "Pnu[i]");
    sfree(Lnu[i], cname, fname, "Lnu[i]");
    sfree(P3[i], cname, fname, "P3[i]");
    sfree(Prhonu[i], cname, fname, "Prhonu[i]");
    sfree(P5[i], cname, fname, "P5[i]");
    sfree(P7[i], cname, fname, "P7[i]");
    sfree(Psigma7[i], cname, fname, "Psigma[i]");
    sfree(Lrhonu[i], cname, fname, "Lrhonu[i]");
    sfree(Lmusigmarhonu[i], cname, fname, "Lmusigmarhonu[i]");
  }
  sfree(Pnu, cname, fname, "Pnu");
  sfree(Lnu, cname, fname, "Lnu");
  sfree(P3, cname, fname, "P3");
  sfree(Prhonu, cname, fname, "Prhonu");
  sfree(P5, cname, fname, "P5");
  sfree(P7, cname, fname, "P7");
  sfree(Psigma7, cname, fname, "Psigma");
  sfree(Lrhonu, cname, fname, "Lrhonu");
  sfree(Lmusigmarhonu, cname, fname, "Lmusigmarhonu");
    
  ffree(mtmp, cname, fname, "mtmp");
  for (int i=0; i<4; i++) ffree(force[i], cname, fname, "force[i]");
  ffree(force, cname, fname, "force");
    
  /*
    ffree(Pnu);
    ffree(P3);
    ffree(Prhonu);
    ffree(P5);
    ffree(P7);
    ffree(Psigma7);
    ffree(L);
    ffree(Lnu);
    ffree(Lrhonu);
    ffree(Lmusigmarhonu);
  */
    
  Convert(CANONICAL);
  for (int p=0; p<degree; p++) Fconvert(sol[p], STAG, CANONICAL);

  return ForceArg(0.0,0.0,0.0);
}

CPS_END_NAMESPACE
