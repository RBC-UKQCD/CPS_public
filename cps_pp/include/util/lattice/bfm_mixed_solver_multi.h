// -*- mode:c++; c-basic-offset: 4 -*-
#ifndef INCLUDED_BFM_MIXED_SOLVER_MULTI_H
#define INCLUDED_BFM_MIXED_SOLVER_MULTI_H

#include <util/lattice/bfm_evo.h>
#include <util/lattice/bfm_mixed_solver.h>
#include <alg/enum_int.h>

#include <util/error.h>

#include <omp.h>

#include <cstdio>
#include <cstdlib>
#include <unistd.h>

//disabling. Seems to have been manually copied from C. Kelly
#if 0
namespace mixed_cg {


    // Convert between single/double precision bfm fermions
    // CK: And do it quickly! The original takes as long as an Mprec!
    template<typename Float_out, typename Float_in>
    inline void threaded_convFermion_fast(Fermion_t out, Fermion_t in,
					  bfm_evo<Float_out> &bfm_out,
					  bfm_evo<Float_in> &bfm_in)
    {
        // Simple copy, this shouldn't be called, right?
        if(sizeof(Float_out) == sizeof(Float_in)) {
            return bfm_out.copy(out, in);
        }

#if 0
	const static int nspinco = 12;
	
	int me,thrlen,throff;
	int work = (bfm_out.gparity? 2:1)*bfm_out.cbLs*bfm_out.simd_cbvol*bfm_out.nsimd*nspinco*2;
	bfm_out.thread_work(work,me,thrlen,throff);
	
	Float_in *x = (Float_in *) in + throff;
	Float_out *y = (Float_out *) out + throff;
	for(int s=0;s<thrlen;++s) y[s] = x[s];

	bfm_out.thread_barrier();
#else
	//Use bfm-3.2 (imported) precisionChange method
	if(sizeof(Float_out) == sizeof(double)) bfm_out.precisionChange(in,out,SingleToDouble,0);
	else bfm_out.precisionChange(in,out,DoubleToSingle,0);
#endif
    }
    


    inline double sigma_sum_recurse(const int &i, const int &start, const int &N, double shifts[], const int &nprod_remaining){
	double out = 0.0;

	for(int j=start;j<N;j++){
	    if(j==i) continue; //skip i
	    double toadd = shifts[j];
	    if(nprod_remaining>1){
		toadd *= sigma_sum_recurse(i,j+1,N,shifts,nprod_remaining-1);
	    }
	    out += toadd;
	}
	return out;
    }

    inline double sigma_prod(const int &i, const int &n, const int &N, double shifts[]){
	if(n==N-1) return 1.0;
  
	int prod_size = N-1-n;
	return sigma_sum_recurse(i,0,N,shifts,prod_size);
    }

    /*CK: Implementation of multi-shift with guesses from 
      J.~C.~Osborn,
      ``Initial guesses for multi-shift solvers,''
      PoS LATTICE {\bf 2008} (2008) 029
      [arXiv:0810.1081 [hep-lat]].
    */

    template <typename Float>
    int CGNE_prec_MdagM_multi_shift_with_guesses(Fermion_t psi[], 
						 Fermion_t src,
						 Fermion_t guesses[],
						 double    mass[],
						 double    alpha[],
						 int       nshift,
						 const double mresidual_in[],
						 int single,
						 bfm_evo<Float> &bfm)
    {
	//'single' appears to sum all the solutions into psi[0]
	int me = bfm.thread_barrier();

	//Threads take a local copy of mresidual_in[] so wno worries about race conditions
	double mresidual[nshift]; for(int i=0;i<nshift;++i) mresidual[i] = mresidual_in[i];

	//Form combination of guesses that results in a new source y in the common Krylov space
	// w = \sum_i  c_i guess_i      where   c_i = \Prod_{j!=i} 1/(shift[j]-shift[i])
	Fermion_t w = bfm.threadedAllocFermion();
	bfm.set_zero(w);
	for(int i = 0; i < nshift; i++){
	    double c_i = 1.0;
	    for(int j=0;j<nshift;j++){
		if(j==i) continue;
		c_i *= 1.0/(mass[j] - mass[i]);
	    }
	    bfm.axpy(w, guesses[i], w, c_i);
	    if(bfm.isBoss() && (!me)) printf("CGNE_prec_MdagM_multi_shift_with_guesses: c_%d = %e\n",i,c_i);
	}
	if(bfm.isBoss() && (!me)) printf("CGNE_prec_MdagM_multi_shift_with_guesses: w norm = %e\n",bfm.norm(w));

	//Form y_i = Prod_{j!=i} (MMdag + shift[j]) w     these are all in the same Krylov space
	Fermion_t y[nshift];
	for(int i=0;i<nshift;i++){
	    y[i] = bfm.threadedAllocFermion(mem_fast);
	    bfm.set_zero(y[i]);
	}

	//We can multiply out. Let MMdag = B
	//y_i = Prod_{j!=i, 0<=j<=N-1} (MMdag + shift[j]) w 
	//    = B^{N-1}w 
	//               + (sum_{j!=i, j<=N-1} shift[j])B^{N-2}w 
	//                                                      +  sum_{j!=i, j<=N-2} shift[j] * ( sum_{k!=i,k>j, k<=N-1} shift[k] ) B^{N-3}w 
	//                                                                                                                                   +... 
	//The coefficient of a term B^n is the sum of all the unique products of (N-1)-n elements of the set of shifts excluding shift[i]: { shift[0], shift[1],.. excl shift[i].... shift[N-2], shift[N-1] }
	//We take a product of size zero to have value 1.0 (as we would with factorials). Indices all start from 0 and end at N-1.
	//For example
	//n=N-1 has coefficent 1.0, which is the product of (N-1)-(N-1)=0 shifts
	//n=N-2 is the sum of all products comprising 1 shift:  sum_{j!=i, j<=N-1} shift[j], which
	//n=N-3 is the sum of all products of 2 shifts: shift[0]*(shift[1] + shift[2] + ... shift[N-1]) + shift[1]*(shift[2] + shift[3] + ... shift[N-1]) + ... + shift[N-2]*shift[N-1]   
	//                                               =  sum_{j!=i, j<=N-2} shift[j] * ( sum_{k!=i,k>j;k<=N-1} shift[k] )
	//etc...
	//n=0  is the product of all shifts bar i

	//Note, we contain the running product of MMdag in w: w_n+1 = MMdag w_n    and definining  w_0 = w
	Fermion_t tmp = bfm.threadedAllocFermion(mem_fast); 
	Fermion_t tmp2 = bfm.threadedAllocFermion(mem_fast); 

	for(int n=0;n<nshift;n++){ //last term   MMdag^{nshift-1}
	    //At this w = MMdag^n w
	    for(int i=0;i<nshift;i++){
		double coeff_i = sigma_prod(i,n,nshift,mass);
		if(bfm.isBoss() && (!me)) printf("CGNE_prec_MdagM_multi_shift_with_guesses: n=%d i=%d sum_prod=%e \n",n,i,coeff_i);
		bfm.axpy(y[i],w,y[i],coeff_i);  // y[i] += coeff[i] * MMdag^n w
	    }
	    bfm.Mprec(w,tmp,tmp2,DaggerNo); //tmp = M w
	    bfm.Mprec(tmp,w,tmp2,DaggerYes); //w = Mdag tmp
	}

	Fermion_t r = bfm.threadedAllocFermion(mem_fast); 
	//Also need r = src - Prod_j (MMdag + shift[j]) w 
	//            = src - (MMdag + shift[i]) Prod_{j!=i, j<=N-1} (MMdag + shift[j]) w    for any i
	//Use i=N-1-1:
	//            = src - (MMdag + shift[N-1]) Prod_{j!=i, j<=N-1} (MMdag + shift[j]) w 
	//            = src - (MMdag + shift[N-1]) y[N-1]
	bfm.Mprec(y[nshift-1],tmp,tmp2,DaggerNo);
	bfm.Mprec(tmp,r,tmp2,DaggerYes); //after this, r = MMdag y[N-1]
	bfm.axpy(r,y[nshift-1],r,mass[nshift-1]); // r = (MMdag + shift[N-1])y[N-1]
	bfm.axpy(r,r,src,-1.0); // r = src - (MMdag + shift[N-1])y[N-1]

	//run standard multi-shift with r as the source
	//To achieve the desired residual on the solution we need to modify the residuals to reflect the difference in source norm
	//In the multi-mass solve we are trying to emulate:
	//|resid|^2 = |orig src|^2 * (orig mresidual)^2
	//Keeping this fixed requires so we want
	//(new mresidual) = sqrt(  |orig src|^2 * (orig mresidual)^2 / |new src|^2 )
	Float orig_src_norm2 = bfm.norm(src);
	Float new_src_norm2 = bfm.norm(r);

	for(int i=0;i<nshift;i++) mresidual[i] *= sqrt(orig_src_norm2/new_src_norm2);

	int iter = bfm.CGNE_prec_MdagM_multi_shift(psi, r, mass, alpha, nshift, mresidual, 0); 

	//Add the solutions to the initial guess vectors y
	for(int n=0;n<nshift;n++) bfm.axpy(psi[n],psi[n],y[n],1.0);

	// Check answers 
	if ( bfm.isBoss() && (!me) ) printf("bfm::CGNE_prec_MdagM_multi_shift_with_guesses: Checking solutions\n");

	for(int s=0; s < nshift; s++) { 
	    bfm.Mprec(psi[s],tmp,tmp2,DaggerNo);
	    bfm.Mprec(tmp,w,tmp2,DaggerYes); //reuse w, now w=MMdag psi[s]
	    bfm.axpy(w,psi[s],w,mass[s]); // w += mass[s]*psi[s]
	    bfm.axpy(r,w,src,-1); // r = src - (MMdag+mass[s])*psi[s]
	    double rn = bfm.norm(r);
	    double cn = bfm.norm(src);
	    if ( bfm.isBoss() && !me ) {
		printf("bfm::CGNE_prec_MdagM_multi_shift_with_guesses: shift[%d] true residual %le \n",s,sqrt(rn/cn));
	    }
	}

	if ( single ) {
	    for(int s=1; s < nshift; s++) { 
		bfm.axpy(psi[0],psi[s],psi[0],1.0);
	    }      
	}

	bfm.threadedFreeFermion(tmp);
	bfm.threadedFreeFermion(tmp2);
	bfm.threadedFreeFermion(w);
	bfm.threadedFreeFermion(r);

	for(int i=0;i<nshift;i++){
	    bfm.threadedFreeFermion(y[i]);
	}  
	return iter;
    }


    inline double sigma_sum_recurse_2(const int &i, const int &k, const int &start, const int &N, double shifts[], const int &nprod_remaining){
	double out = 0.0;

	for(int j=start;j<N;j++){
	    if(j==i || j==k) continue; //skip i,k
	    double toadd = shifts[j];
	    if(nprod_remaining>1){
		toadd *= sigma_sum_recurse_2(i,k,j+1,N,shifts,nprod_remaining-1);
	    }
	    out += toadd;
	}
	return out;
    }

    inline double sigma_prod_2(const int &i, const int &k, const int &n, const int &N, double shifts[]){
	if(n==N-2) return 1.0;
  
	int prod_size = N-2-n;
	return sigma_sum_recurse_2(i,k,0,N,shifts,prod_size);
    }

    //CK: Multi-mass with multiple sources, using the method in the Osborn paper cited above
    template <typename Float>
    int CGNE_prec_MdagM_multi_shift_multi_src(Fermion_t psi[], 
					      Fermion_t src[],
					      double    mass[],
					      double    alpha[],
					      int       nshift,
					      const double mresidual_in[],
					      int single,
					      bfm_evo<Float> &bfm)
    {
	//'single' appears to sum all the solutions into psi[0]
	int me = bfm.thread_barrier();

	//Not know if mresidual_in is shared or not. To be safe , threads each take a local copy of the residuals that they can 
	//all simultaneously modify without worrying about race conditions
	double mresidual[nshift];
	for(int i=0;i<nshift;++i) mresidual[i] = mresidual_in[i];

	//Here we need to form  y_i = \sum_{k\neq i} [ \prod_{j\neq i,k} (A + \sigma_j)/(\sigma_j-\sigma_k) ] (b_i - b_k)/(\sigma_i-\sigma_k)
	//                          = \sum_{k\neq i} [ \prod_{j\neq i,k} (A + \sigma_j) ] \prod_{j\neq k} 1/(\sigma_j-\sigma_k) (b_i - b_k)
	//                          = \sum_{k\neq i} [ \prod_{j\neq i,k} (A + \sigma_j) ] c_k (b_i - b_k)
	//where c_k = \prod_{j\neq k} 1/(\sigma_j-\sigma_k) is calculated beforehand

	Fermion_t y[nshift];
	Fermion_t Apowm_src[nshift];
	double c[nshift];
	for(int i=0;i<nshift;i++){
	    y[i] = bfm.threadedAllocFermion(mem_fast);
	    Apowm_src[i] = bfm.threadedAllocFermion(mem_fast);
	    bfm.copy(Apowm_src[i],src[i]);
	    bfm.set_zero(y[i]);

	    //Calculate coefficients c_i
	    c[i] = 1.0;
	    for(int j=0;j<nshift;j++){
		if(j==i) continue;
		c[i] *= 1.0/(mass[j]-mass[i]);
	    }
	}

	//I don't think we can avoid doing ~nshift^2 matrix multiplications, but I believe we can avoid storing nshift^2 vectors
	Fermion_t tmp = bfm.threadedAllocFermion(mem_fast); 
	Fermion_t tmp2 = bfm.threadedAllocFermion(mem_fast); 

	//y_i = \sum_{k\neq i} [ \prod_{j\neq i,k} (A + \sigma_j) ] c_k (b_i - b_k)
	//    = \sum_{k\neq i} c_k {   [ \prod_{j\neq i,k} (A + \sigma_j) ] b_i   -    [ \prod_{j\neq i,k} (A + \sigma_j) ] b_k  }

	//These objects are the fundamental elements  \prod_{j\neq i,k} (A + \sigma_j) ] b_i
	//Multiplying out                                                                     = A^{N-2}   +   A^{N-3}\sum_{j\neq i,k} \sigma_j  + A^{N-4} \sum_{j\neq i,k} \sigma_j \sum_{l>j, l\neq i,k} \sigma_l + ..... 
                                      
	for(int m=0;m<nshift-1;m++){ //last term   MMdag^{nshift-2}
	    for(int i=0;i<nshift;i++){
		for(int k=0;k<nshift;k++){ 	//\sum_{k\neq i}
		    if(k==i) continue;
		    double coeff = c[k]*sigma_prod_2(i,k,m,nshift,mass);
		    bfm.axpy(y[i], Apowm_src[i], y[i], coeff  );
		    bfm.axpy(y[i], Apowm_src[k], y[i], -coeff );
		}
	    }

	    for(int n=0;n<nshift;n++){
		bfm.Mprec(Apowm_src[n],tmp,tmp2,DaggerNo); //tmp = M w
		bfm.Mprec(tmp,Apowm_src[n],tmp2,DaggerYes); //w = Mdag tmp
	    }
	}

	Fermion_t r = bfm.threadedAllocFermion(mem_fast); 

	//Also need r = src_i -  (A+\sigma_i)y_i  for any i    
	bfm.Mprec(y[nshift-1],tmp,tmp2,DaggerNo);
	bfm.Mprec(tmp,r,tmp2,DaggerYes);
	bfm.axpy(r,y[nshift-1],r,mass[nshift-1]);
	bfm.axpy(r,r,src[nshift-1],-1.0);

	//run standard multi-shift with r as the source
	//To achieve the desired residual on the solution we need to modify the residuals to reflect the difference in source norm
	//In the multi-mass solve we are trying to emulate:
	//|resid|^2 = |orig src|^2 * (orig mresidual)^2
	//Keeping this fixed requires so we want
	//(new mresidual) = sqrt(  |orig src|^2 * (orig mresidual)^2 / |new src|^2 )
	Float new_src_norm2 = bfm.norm(r);

	for(int i=0;i<nshift;i++){
	    Float orig_src_norm2 = bfm.norm(src[i]);
	    if ( bfm.isBoss() && !me ) printf("bfm::CGNE_prec_MdagM_multi_src: input src[%d] norm2 %le\n",i,orig_src_norm2);
	    mresidual[i] *= sqrt(orig_src_norm2/new_src_norm2);
	}
	if ( bfm.isBoss() && !me ) printf("bfm::CGNE_prec_MdagM_multi_src: r norm2 %le\n",new_src_norm2);
 
	int iter = bfm.CGNE_prec_MdagM_multi_shift(psi, r, mass, alpha, nshift, mresidual, 0); 

	//Add the solutions to the initial guess vectors y
	for(int n=0;n<nshift;n++) bfm.axpy(psi[n],psi[n],y[n],1.0);

	// Check answers 
	if ( bfm.isBoss() && (!me) ) printf("bfm::CGNE_prec_MdagM_multi_src: Checking solutions\n");

	for(int s=0; s < nshift; s++) { 
	    bfm.Mprec(psi[s],tmp,tmp2,DaggerNo);
	    bfm.Mprec(tmp,Apowm_src[0],tmp2,DaggerYes); //reuse Apowm_src[0], now w=MMdag psi[s]
	    bfm.axpy(Apowm_src[0],psi[s],Apowm_src[0],mass[s]); // Apowm_src[0] += mass[s]*psi[s]
	    bfm.axpy(r,Apowm_src[0],src[s],-1); // r = src - (MMdag+mass[s])*psi[s]
	    double rn = bfm.norm(r);
	    double cn = bfm.norm(src[s]);
	    if ( bfm.isBoss() && !me ) {
		printf("bfm::CGNE_prec_MdagM_multi_src: shift[%d] true residual %le \n",s,sqrt(rn/cn));
	    }
	}

	if ( single ) {
	    for(int s=1; s < nshift; s++) { 
		bfm.axpy(psi[0],psi[s],psi[0],1.0);
	    }      
	}

	bfm.threadedFreeFermion(tmp);
	bfm.threadedFreeFermion(tmp2);
	bfm.threadedFreeFermion(r);

	for(int i=0;i<nshift;i++){
	    bfm.threadedFreeFermion(Apowm_src[i]);
	    bfm.threadedFreeFermion(y[i]);
	}  
	return iter;
    }



    template<typename Float>
    inline void MdagMplusShift(Fermion_t in, Fermion_t out, const double &shift, Fermion_t tmp1, Fermion_t tmp2, bfm_evo<Float> &bfm){
	bfm.Mprec(in, tmp1, tmp2, DaggerNo);
	bfm.Mprec(tmp1, out, tmp2, DaggerYes); 
	bfm.axpy(out, in, out, shift);
    }

    //CK: mixed precision multi-mass using multiple single precision restarted inner loop.
    //    Does not work very well because the rediduals get very large such that the required stopping conditions
    //    are less than single precision accuracy

    // Both sol and src are double precision fermions. Single precision
    // solver is only used internally.
    // mass, alpha and nshift as usual
    // dresidual are the target residuals
    // fresidual are the initial single precision residuals, which are dynamically modified during the solve
    //
    // Things to be set before using this function:
    //
    // double precision solver mass, max iteration
    // number.
    //
    // single precision solver mass, max iteration
    // number.
    //
    // Gauge field must be initialized for both double and single prec
    // solvers.
    //
    // the communication subsystem must be ready for bfm_d to use (due to
    // the way bfmcommspi is written, one must reinitialize the
    // communication object when switching between single and double
    // precisions).
    //
    // max_cycle: the maximum number of restarts will be performed.
    inline int threaded_cg_mixed_restarted_multi_shift_MdagM(Fermion_t sol[], 
							     Fermion_t src,
							     double    mass[],
							     double    alpha[],
							     int       nshift,
							     const double dresidual[],
							     const double fresidual_in[],
							     int single,
							     bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f,
							     int max_cycle){
	int me = bfm_d.thread_barrier();

	double fresidual[nshift]; for(int i=0;i<nshift;++i) fresidual[i] = fresidual_in[i]; //local thread copy that can be modified freely

	Fermion_t src_d = bfm_d.threadedAllocFermion();
	Fermion_t tv1_d = bfm_d.threadedAllocFermion();
	Fermion_t tv2_d = bfm_d.threadedAllocFermion();

	double src_norm = bfm_d.norm(src);
	//Source and solution locations for single precision input and output
	Fermion_t sol_f[nshift];
	Fermion_t src_f[nshift];
	int finished[nshift];
	double stop[nshift];
	for(int n=0;n<nshift;n++){
	    sol_f[n] = bfm_f.threadedAllocFermion();
	    src_f[n] = bfm_f.threadedAllocFermion();
	    finished[n] = 0;
	    stop[n] = src_norm * dresidual[n] * dresidual[n];
	}

	//Do an initial single precision solve with regular multi-mass solver using input residuals
	if(bfm_f.isBoss() && !me) printf("threaded_cg_mixed_restarted_multi_shift_MdagM: Doing initial single precision solve\n");
	{
	    threaded_convFermion(src_f[0], src, bfm_f, bfm_d);
	    switch_comm(bfm_f, bfm_d);
	    bfm_f.CGNE_prec_MdagM_multi_shift(sol_f,src_f[0],mass,alpha,nshift,fresidual,0);
	    for(int n=0;n<nshift;n++) threaded_convFermion(sol[n], sol_f[n], bfm_d, bfm_f);
	    switch_comm(bfm_d, bfm_f);
	}

	//Perform restarted multi-mass until the double prec residual meets the target
	if(bfm_f.isBoss() && !me) printf("threaded_cg_mixed_restarted_multi_shift_MdagM: Starting main iteration loop\n");
    
	int iter = 0;
	for(int i = 0; i < max_cycle; ++i) {
	    // compute double precision rsd and also new RHS vector for each shift
	    int fin_count = 0;
	    for(int n=0;n<nshift;n++){
		if(!finished[n]){
		    MdagMplusShift<double>(sol[n], tv1_d, mass[n], src_d, tv2_d, bfm_d); // tv1_d = (MdagM + mass[n]) * sol[n]
		    double norm = bfm_d.axpy_norm(src_d, tv1_d, src, -1.);

		    //ad hoc stopping condition from cg_mixed implementation
		    if(norm < 100. * stop[n]) finished[n] = 1;

		    if(bfm_f.isBoss() && !me) printf("threaded_cg_mixed_restarted_multi_shift_MdagM: iter = %d shift = %d rsd = %17.10e(d) stop = %17.10e(d) finished = %d\n", i,n, norm, stop[n],finished[n]);

		    while(norm * fresidual[n] * fresidual[n] < stop[n]) fresidual[n] *= 2;

		    threaded_convFermion(src_f[n], src_d, bfm_f, bfm_d);
		}
		fin_count += finished[n];
	    }
	    if(fin_count == nshift) break; //stop when all have finished
	  
	    switch_comm(bfm_f, bfm_d);

	    iter += CGNE_prec_MdagM_multi_shift_multi_src(sol_f, src_f, mass, alpha, nshift, fresidual, 0, bfm_f);

	    switch_comm(bfm_d, bfm_f);
	    for(int n=0;n<nshift;n++){
		threaded_convFermion(tv1_d, sol_f[n], bfm_d, bfm_f);
		bfm_d.axpy(sol[n], tv1_d, sol[n], 1.);
	    }
	}

	bfm_d.threadedFreeFermion(src_d);
	bfm_d.threadedFreeFermion(tv1_d);
	bfm_d.threadedFreeFermion(tv2_d);

	for(int i=0;i<nshift;i++){
	    bfm_f.threadedFreeFermion(sol_f[i]);
	    bfm_f.threadedFreeFermion(src_f[i]);
	}

	if(bfm_f.isBoss() && !me) printf("threaded_cg_mixed_restarted_multi_shift_MdagM: Running double precision multi-mass using single precision version as guess");
	iter += CGNE_prec_MdagM_multi_shift_with_guesses(sol,src,sol,mass,alpha,nshift,dresidual,single,bfm_d);
	return iter;
    }



    //CK: mixed precision multi-mass using single precision solve as guess for double precision
    // Both sol and src are double precision fermions. Single precision
    // solver is only used internally.
    // mass, alpha and nshift as usual
    // dresidual are the target residuals
    // fresidual are the initial single precision residuals, which are dynamically modified during the solve
    //
    // Things to be set before using this function:
    //
    // double precision solver mass, max iteration
    // number.
    //
    // single precision solver mass, max iteration
    // number.
    //
    // Gauge field must be initialized for both double and single prec
    // solvers.
    //
    // the communication subsystem must be ready for bfm_d to use (due to
    // the way bfmcommspi is written, one must reinitialize the
    // communication object when switching between single and double
    // precisions).
    //
    // max_cycle: the maximum number of restarts will be performed.
    inline int threaded_cg_mixed_single_prec_as_guess_multi_shift_MdagM(Fermion_t sol[], 
									Fermion_t src,
									double    mass[],
									double    alpha[],
									int       nshift,
									double dresidual[],
									double fresidual[],
									int single,
									bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f){
	int me = bfm_d.thread_barrier();

	//Source and solution locations for single precision input and output
	Fermion_t sol_f[nshift];
	for(int n=0;n<nshift;n++) sol_f[n] = bfm_f.threadedAllocFermion();
	Fermion_t src_f = bfm_f.threadedAllocFermion();

	int iter = 0;

	//Do an initial single precision solve with regular multi-mass solver using input residuals
	if(bfm_f.isBoss() && !me) printf("threaded_cg_mixed_single_prec_as_guess_multi_shift_MdagM: Doing single precision solve\n");
	{
	    threaded_convFermion(src_f, src, bfm_f, bfm_d);
	    switch_comm(bfm_f, bfm_d);
	    iter += bfm_f.CGNE_prec_MdagM_multi_shift(sol_f,src_f,mass,alpha,nshift,fresidual,0);
	    for(int n=0;n<nshift;n++) threaded_convFermion(sol[n], sol_f[n], bfm_d, bfm_f);
	    switch_comm(bfm_d, bfm_f);
	}

	bfm_f.threadedFreeFermion(src_f);
	for(int n=0;n<nshift;n++) bfm_f.threadedFreeFermion(sol_f[n]);
    
	if(bfm_f.isBoss() && !me) printf("threaded_cg_mixed_single_prec_as_guess_multi_shift_MdagM: Running double precision multi-mass using single precision version as guess");
	iter += CGNE_prec_MdagM_multi_shift_with_guesses(sol,src,sol,mass,alpha,nshift,dresidual,single,bfm_d);
	return iter;
    }


    //CK: "Single Shift Inverter" : Modified version of int bfm::CGNE_prec that solves (MdagM + shift) psi = src
    template<typename Float>
    int threaded_CGNE_MdagM_plus_shift(Fermion_t psi, Fermion_t src, Float shift, bfm_evo<Float> &bfm)
    {
	//Standard CG algorithm from BFM:
	//(Use subscript to label iteration)

	//r_1 = MMdag psi - src,   p_1 = MMdag psi - src,   c_1 = |r_1|^2 = |p_1|^2
	//Iteration:
	//d_k = |M p_k|^2 = p_k^dag M^dag M p_k
	//a_k = c_k / d_k
	//r_k+1 = MMdag p_k - a_k r_k
	//c_k+1 = |r_k+1|^2
	//b_k = c_k+1 / c_k
	//psi = a_k p_k + psi
	//p_k+1 = b_k p_k + r_k+1
    
	//Note: norm(vec) is actually |vec|^2

	//Shift modified version should look similar:
	//r_1 = (MMdag+shift) psi - src,   p_1 = (MMdag+shift) psi - src,   c_1 = |r_1|^2 = |p_1|^2
	//Iteration:
	//d_k = p_k^dag M^dag M p_k + shift * p_k^dag p_k 
	//a_k = c_k / d_k
	//r_k+1 = (MMdag+shift) p_k - a_k r_k
	//c_k+1 = |r_k+1|^2
	//b_k = c_k+1 / c_k
	//psi = a_k p_k + psi
	//p_k+1 = b_k p_k + r_k+1

	int me = bfm.thread_barrier();

	int verbose = bfm.verbose;

	double f;
	double cp,c,a,d,b;
	double residual = bfm.residual;
	int max_iter = bfm.max_iter;

	if ( bfm.isBoss() && (!me) ) { 
	    bfm.InverterEnter();
	}

	Fermion_t p   = bfm.threadedAllocFermion(mem_fast); 
	Fermion_t tmp = bfm.threadedAllocFermion(mem_fast); 
	Fermion_t mp  = bfm.threadedAllocFermion(mem_fast); 
	Fermion_t mmp = bfm.threadedAllocFermion(mem_fast); 
	Fermion_t r   = bfm.threadedAllocFermion(mem_fast); 

	//Initial residual computation & set up
	double guess = bfm.norm(psi);
	d = bfm.Mprec(psi,mp,tmp,DaggerNo);
	bfm.Mprec(mp,mmp,tmp,DaggerYes);
	b = bfm.axpy_norm(mmp,psi,mmp,shift); //MMdag psi + shift*psi
	cp= bfm.axpy_norm (r, mmp, src,-1.0);
	a = bfm.axpy_norm (p, mmp, src,-1.0);

	//a = bfm.norm(p);
	//cp= bfm.norm(r);
	//r_1 = (MMdag+shift) psi - src,   p_1 = (MMdag+shift) psi - src,   c_1 = |r_1|^2 = |p_1|^2

	Float ssq =  bfm.norm(src);
	if ( verbose && bfm.isBoss() && !me ) {
	    printf("mixed_cg::CGNE_MdagM_plus_shift gues %le \n",guess);
	    printf("mixed_cg::CGNE_MdagM_plus_shift src  %le \n",ssq);
	    printf("mixed_cg::CGNE_MdagM_plus_shift  Mp  %le \n",d);
	    printf("mixed_cg::CGNE_MdagM_plus_shift  (MMdag + shift)p %le \n",b);
	    printf("mixed_cg::CGNE_MdagM_plus_shift   r  %le \n",cp);
	    printf("mixed_cg::CGNE_MdagM_plus_shift   p  %le \n",a);
	}
	Float rsq =  residual* residual*ssq;

	//Check if guess is really REALLY good :)
	if ( cp <= rsq ) {
	    if ( verbose && bfm.isBoss() && !me ) {
		printf("mixed_cg::CGNE_MdagM_plus_shift k=0 converged - suspiciously nice guess %le %le\n",cp,rsq);
	    }
	    bfm.threadedFreeFermion(tmp);
	    bfm.threadedFreeFermion(p);
	    bfm.threadedFreeFermion(mp);
	    bfm.threadedFreeFermion(mmp);
	    bfm.threadedFreeFermion(r);
	    if ( bfm.isBoss() && (!me) ) { 
		bfm.InverterExit();
	    }
	    return 0;
	}
	if ( verbose && bfm.isBoss() && !me ) 
	    printf("mixed_cg::CGNE_MdagM_plus_shift k=0 residual %le rsq %le\n",cp,rsq);

	if ( bfm.isBoss() && !me ) {
	    if ( bfm.watchfile ) {
		printf("mixed_cg::CGNE_MdagM_plus_shift watching file \"%s\"\n",bfm.watchfile);
	    }    
	}
	struct timeval start,stop;
	if(bfm.isBoss() && !me)	gettimeofday(&start,NULL);

	for (int k=1;k<=max_iter;k++){
	    bfm.iter=k;
	    uint64_t t_iter_1=GetTimeBase();

	    c=cp; 

	    uint64_t t_mprec_1=GetTimeBase();

	    //d_k = p_k^dag M^dag M p_k + shift * p_k^dag p_k 
	    d = bfm.Mprec(p,mp,tmp,0,1);

	    double norm_p = bfm.norm(p);
	    d += shift * norm_p;

	    uint64_t t_mprec_2=GetTimeBase();
	    a = c/d;

	    uint64_t t_mprec_3=GetTimeBase();

	    bfm.Mprec(mp,mmp,tmp,1); 
	    bfm.axpy(mmp,p,mmp,shift); // mmp = MMdag p + shift * p

	    uint64_t t_mprec_4=GetTimeBase();

	    uint64_t tr1=GetTimeBase();

	    cp = bfm.axpy_norm(r,mmp,r,-a); //r_k+1 = (MMdag+shift) p_k - a_k r_k

	    b = cp/c;
	    uint64_t tr2=GetTimeBase();

	    uint64_t tpsi1=GetTimeBase();
	    bfm.axpy(psi,p,psi,a);
	    uint64_t tpsi2=GetTimeBase();
	    // New (conjugate/M-orthogonal) search direction
	    uint64_t tp1=GetTimeBase();
	    bfm.axpy(p,p,r,b);
	    uint64_t tp2=GetTimeBase();

	    uint64_t t_iter_2=GetTimeBase();

	    // verbose nonsense
	    if ( (bfm.iter==bfm.time_report_iter) && bfm.isBoss() && (!me)  && verbose) {
      
		int lx = bfm.node_latt[0];
		int ly = bfm.node_latt[1];
		int lz = bfm.node_latt[2];
		int lt = bfm.node_latt[3];

		int cb4dsites = (lx*ly*lz*lt)/2;
		printf("fermionCacheFootprint: %ld \n",7*bfm.axpyBytes()/3);
		printf("gauge  CacheFootprint: %ld \n",2*18*8*cb4dsites*2);
		printf("fermionVecBytes      : %ld \n",bfm.axpyBytes()/3);
		printf("axpyBytes            : %ld \n",bfm.axpyBytes());
		printf("axpy      (soln)     : %ld cyc %le MB/s\n",(tpsi2-tpsi1),(double)bfm.axpyBytes()*1600./(tpsi2-tpsi1));
		printf("axpy_norm (residual) : %ld cyc %le MB/s\n",(tr2-tr1),(double)bfm.axpyBytes()*1600./(tr2-tr1));
		printf("axpy      (search)   : %ld cyc %le MB/s\n",(tp2-tp1),(double)bfm.axpyBytes()*1600./(tp2-tp1));
		printf("Iter time            : %ld cyc\n",t_iter_2-t_iter_1);
		printf("linalg time          : %ld cyc\n",t_iter_2-t_iter_1-(t_mprec_2-t_mprec_1)-(t_mprec_4-t_mprec_3));
		printf("Mprec time           : %ld cyc\n",t_mprec_2-t_mprec_1);
		printf("Mprec time           : %ld cyc\n",t_mprec_4-t_mprec_3);
		fflush(stdout);
	    }


	    if ( ((k%100 == 0) && (verbose!=0)) || (verbose > 10)){
		if ( bfm.isBoss() && !me ) {
		    printf("mixed_cg::CGNE_MdagM_plus_shift: k=%d r^2=%le %le %lx\n",k,cp,sqrt(cp/ssq),&bfm);
		}
	    }

	    // Stopping condition
	    if ( cp <= rsq ) { 
		//I did not update the flops count so I have commented them out
		struct timeval diff;

		if ( bfm.isBoss() && !me ){
		    gettimeofday(&stop,NULL);
		    timersub(&stop,&start,&diff);
		}

		if ( bfm.isBoss() && !me ) printf("mixed_cg::CGNE_MdagM_plus_shift converged in %d iterations\n",k);
		if ( bfm.isBoss() && !me ) printf("mixed_cg::CGNE_MdagM_plus_shift converged in %d.%6.6d s\n",diff.tv_sec,diff.tv_usec);


		//double flops = mprecFlops()*2.0 + 2.0*axpyNormFlops() + axpyFlops()*2.0;
		//flops = flops * k;

		//double t = diff.tv_sec*1.0E6 + diff.tv_usec;
		// if ( isBoss()&& !me ) 
		//   printf("mixed_cg::CGNE_MdagM_plus_shift: %d mprec flops/site\n",mprecFlopsPerSite());
		// if ( isBoss()&& !me ) printf("mixed_cg::CGNE_MdagM_plus_shift: %le flops\n",flops);
		// if ( isBoss()&& !me ) printf("mixed_cg::CGNE_MdagM_plus_shift: %le mflops per node\n",flops/t);

		if ( bfm.isBoss() && !me ){ printf("mixed_cg::CGNE_MdagM_plus_shift calculating true resid. V0\n"); fflush(stdout); } //DEBUG
		bfm.Mprec(psi,mp,tmp,0);
		bfm.Mprec(mp,mmp,tmp,1); 
		bfm.axpy(mmp,psi,mmp,shift);

		double resid = bfm.axpy_norm(tmp,src,mmp,-1.0);

		double src_norm = bfm.norm(src);
		double true_residual = sqrt(resid/src_norm);
		if ( bfm.isBoss() && !me ) 
		    printf("mixed_cg::CGNE_MdagM_plus_shift: true residual is %le \n",true_residual);

		if ( bfm.isBoss() && !me ){ printf("mixed_cg::CGNE_MdagM_plus_shift cleaning up\n"); fflush(stdout); } //DEBUG

		bfm.threadedFreeFermion(tmp);
		bfm.threadedFreeFermion(p);
		bfm.threadedFreeFermion(mp);
		bfm.threadedFreeFermion(mmp);
		bfm.threadedFreeFermion(r);
#ifdef LIST_ENGINE
		if ( bfm.list_engine) bfm.L1P_PatternUnconfigure();
#endif
		if ( bfm.isBoss() && (!me) ) { 
		    bfm.InverterExit();
		}
		return k;
	    }

	}
	if ( bfm.isBoss() && !me ) printf("mixed_cg::CGNE_MdagM_plus_shift: CG not converged \n");
	bfm.threadedFreeFermion(tmp);
	bfm.threadedFreeFermion(p);
	bfm.threadedFreeFermion(mp);
	bfm.threadedFreeFermion(mmp);
	bfm.threadedFreeFermion(r);
#ifdef LIST_ENGINE
	if (bfm.list_engine ) bfm.L1P_PatternUnconfigure();
#endif
	if ( bfm.isBoss() && (!me) ) { 
	    bfm.InverterExit();
	}

	return -1;
    }

    //CK: Single precision solve followed by defect correction loop using single shift solver independently
    //    for each shift
    // Both sol and src are double precision fermions. Single precision
    // solver is only used internally.
    //
    // Things to be set before using this function:
    //
    // double precision solver mass, stopping condition, max iteration
    // number.
    //
    // single precision solver mass, stopping condition, max iteration
    // number.
    //
    // Gauge field must be initialized for both double and single prec
    // solvers.
    //
    // the communication subsystem must be ready for bfm_d to use (due to
    // the way bfmcommspi is written, one must reinitialize the
    // communication object when switching between single and double
    // precisions).
    //
    // max_cycle: the maximum number of restarts will be performed.

    //fresidual are the residuals used for the initial single precision solve

    //min_fp_resid: the smallest value for the residual of the initial single-precision multi-mass solve
    inline int threaded_cg_mixed_defect_correction_multi_shift_MdagM(Fermion_t sol[], Fermion_t src,
								     double  mass[], double  alpha[], 
								     bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f, int nshift,
								     double mresidual[], double fresidual[], int single, int max_cycle){
	int me = bfm_d.thread_barrier();
	double frsd = bfm_f.residual; //save original residual for later restoration

	//First we perform the multi-mass inversion using the single-precision solver
	Fermion_t src_f = bfm_f.threadedAllocFermion();
	Fermion_t sol_f[nshift];
	for(int i=0;i<nshift;i++) sol_f[i] = bfm_f.threadedAllocFermion();

	threaded_convFermion(src_f, src, bfm_f, bfm_d);
	switch_comm(bfm_f, bfm_d);

	int single_prec_iter = bfm_f.CGNE_prec_MdagM_multi_shift(sol_f, src_f, mass, alpha, nshift, fresidual, 0);
	
	if(bfm_f.isBoss() && !me) {
	    printf("threaded_cg_mixed_defect_correction_multi_shift_MdagM: single-prec multi-shift iter = %d\n",single_prec_iter);
	}

	//Now we loop through the shifted solutions and do defect-correction on each individually
	switch_comm(bfm_d, bfm_f);
	for(int i=0;i<nshift;i++) threaded_convFermion(sol[i], sol_f[i], bfm_d, bfm_f);
	
	double src_norm = bfm_d.norm(src);

	Fermion_t tv1_d = bfm_d.threadedAllocFermion();
	Fermion_t tv2_d = bfm_d.threadedAllocFermion();
	Fermion_t src_d = bfm_d.threadedAllocFermion();

	int iter = 0;
	for(int shift=0;shift<nshift;shift++){
	    double stop = src_norm * mresidual[shift]*mresidual[shift];
	    bfm_f.thread_barrier(); //make sure no threads have yet to write to bfm_f.residual from previous loop cycle
	    bfm_f.residual = mresidual[shift];

	    for(int i = 0; i < max_cycle; ++i) {
		// compute double precision rsd and also new RHS vector.
		bfm_d.Mprec(sol[shift],   tv1_d, src_d, 0, 0); //here src_d is just used as a temp storage
		bfm_d.Mprec(tv1_d, tv2_d, src_d, 1, 0); // tv2_d = MdagM * sol
		bfm_d.axpy(tv2_d,sol[shift],tv2_d,mass[shift]); //tv2_d = (MdagM + shift)* sol

		double norm = bfm_d.axpy_norm(src_d, tv2_d, src, -1.);

		// Hantao's ad hoc stopping condition
		if(norm < 100. * stop) break;

		if(!me)
		    while(norm * bfm_f.residual * bfm_f.residual < stop) bfm_f.residual *= 2;
		//bfm_f.thread_barrier(); //Not needed because there is a barrier at start of next call

		if(bfm_f.isBoss() && !me) {
		    printf("threaded_cg_mixed_defect_correction_multi_shift_MdagM: shift = %d, defect correction cycle = %d rsd = %17.10e(d) stop = %17.10e(d)  [True resid %17.10e(d), next single prec target resid %17.10e]\n",
			   shift,i, norm, stop, sqrt(norm/src_norm), bfm_f.residual );
		}

		//We need to invert MdagM + shift, for which we cannot use the regular inverter. Use my optimised single-shift inverter
		//Could also use the multi-shift with a single shift, but we can avoid some overhead by using my optimised version
		threaded_convFermion(src_f, src_d, bfm_f, bfm_d);
		switch_comm(bfm_f, bfm_d);

		bfm_f.set_zero(sol_f[shift]);
		iter += threaded_CGNE_MdagM_plus_shift<float>(sol_f[shift],src_f,mass[shift],bfm_f);

		switch_comm(bfm_d, bfm_f);
		threaded_convFermion(tv1_d, sol_f[shift], bfm_d, bfm_f);

		bfm_d.axpy(sol[shift], tv1_d, sol[shift], 1.);
	    }
	    bfm_f.residual = frsd; //restore original single precision residual at end of each step
	}

	bfm_d.threadedFreeFermion(src_d);
	bfm_d.threadedFreeFermion(tv1_d);
	bfm_d.threadedFreeFermion(tv2_d);

	for(int i=0;i<nshift;i++) bfm_f.threadedFreeFermion(sol_f[i]);
	bfm_f.threadedFreeFermion(src_f);

	for(int shift=0;shift<nshift;shift++){
	    if(bfm_d.isBoss() && !me) printf("threaded_cg_mixed_defect_correction_multi_shift_MdagM: doing final inversion for shift %d using corrected solution as guess\n",shift);

	    double restore_resid = bfm_d.residual;
	    bfm_d.thread_barrier(); //make sure all threads get the same value before we change it

	    bfm_d.residual = mresidual[shift];
	    iter += threaded_CGNE_MdagM_plus_shift<double>(sol[shift],src,mass[shift],bfm_d);
	    bfm_d.residual = restore_resid;
	    //bfm_d.thread_barrier(); //Not needed because barrier in next call

	    double sol_norm = bfm_d.norm(sol[shift]);      
	    if(bfm_d.isBoss() && !me) printf("threaded_cg_mixed_defect_correction_multi_shift_MdagM: final sol[%d] norm = %17.10e\n", shift,sol_norm);
	}

	if ( single ) {
	    for(int s=1; s < nshift; s++) { 
		bfm_d.axpy(sol[0],sol[s],sol[0],1.0);
	    }      
	}
	return iter;
    }
    

    //CK 2014: The version below performs the multi-shift with the matrix multiplication in single precision. The residual is stored in single precision, but the search directions and solution are
    //stored in double precision. Every update_freq iterations the residual is corrected in double precision. 
    
    //Note that the final double precision residuals may not be as good as desired, so you may want to perform defect correction on each pole afterwards. I have added a version that does this extra step below.

    inline int threaded_cg_mixed_multi_shift_MdagM_sp_relup_dp(Fermion_t psi[], 
							       Fermion_t src,
							       double    mass[],
							       double    alpha[],
							       int       nshift,
							       double mresidual[],
							       int single,
							       bfm_evo<float> &bfm_f,
							       bfm_evo<double> &bfm_d,
							       int update_freq = 100,
							       int report_freq = -1)
    {
	//NOTE: Assumes bfm_d comms are active
	//update_freq is the frequency at which the reliable update step is performed
	//report_freq prints the double precision true residual when k % report_freq = 0. Use -1 to disable

	int me = bfm_d.thread_barrier();
 
	double    bs [nshift];
	double    rsq[nshift];
	double    z[nshift][2];
	int       converged[nshift];

	const int       primary =0;

	//Primary shift fields CG iteration
	double a,b,c,d;
	double cp,bp; //prev

	//Single precision fields
	Fermion_t r = bfm_f.threadedAllocFermion(mem_slow); //residual vector, single precision  
	Fermion_t tmp = bfm_f.threadedAllocFermion(mem_fast); 
	Fermion_t p = bfm_f.threadedAllocFermion(mem_fast); 
	Fermion_t mp  = bfm_f.threadedAllocFermion(mem_fast); 
	Fermion_t mmp = bfm_f.threadedAllocFermion(mem_fast); 
	Fermion_t src_f = bfm_f.threadedAllocFermion(mem_slow);	
	mixed_cg::threaded_convFermion_fast(src_f, src, bfm_f, bfm_d);

	//Double precision fields
	Fermion_t p_d = bfm_d.threadedAllocFermion(mem_fast); //search direction, double precision
	Fermion_t tmp_d = bfm_d.threadedAllocFermion(mem_fast); 
	Fermion_t mp_d  = bfm_d.threadedAllocFermion(mem_fast); 
	Fermion_t mmp_d = bfm_d.threadedAllocFermion(mem_fast); 
	Fermion_t ps_d [nshift]; // search directions (double precision)

	for(int i=0;i<nshift;i++){
	    ps_d[i] = bfm_d.threadedAllocFermion(mem_slow);
	    converged[i]=0;
	}

#define DEALLOCATE_ALL							\
	bfm_f.threadedFreeFermion(r);					\
	bfm_f.threadedFreeFermion(tmp);					\
	bfm_f.threadedFreeFermion(p);					\
	bfm_f.threadedFreeFermion(mp);					\
	bfm_f.threadedFreeFermion(mmp);					\
	bfm_f.threadedFreeFermion(src_f);				\
	bfm_d.threadedFreeFermion(p_d);					\
	bfm_d.threadedFreeFermion(tmp_d);				\
	bfm_d.threadedFreeFermion(mp_d);				\
	bfm_d.threadedFreeFermion(mmp_d);				\
	for(int s=0;s<nshift;s++) bfm_d.threadedFreeFermion(ps_d[s])

	// Check lightest mass
	for(int s=0;s<nshift;s++){
	    if ( mass[s] < mass[primary] ) {
		printf("First shift not lightest - oops\n");
		exit(-1);
	    }
	}

	cp = bfm_d.norm(src);
	for(int s=0;s<nshift;s++){
	    rsq[s] = cp * mresidual[s] * mresidual[s];
	    bfm_d.copy(ps_d[s],src);
	}
	// r and p for primary
	bfm_f.copy(r,src_f); //residual vector in single prec
	bfm_d.copy(p_d,src);

	double rn = cp; //norm of src = p_d
  
	mixed_cg::switch_comm(bfm_f, bfm_d);
	mixed_cg::threaded_convFermion_fast(p, p_d, bfm_f, bfm_d);

	d= bfm_f.Mprec(p,mp,tmp,DaggerNo,1); //mp = Mpc p,  what is the 'norm' of?? I think its |Mpc p|^2
	bfm_f.Mprec(mp,mmp,tmp,DaggerYes); //mmp = Mpc^dag mp = Mpc^dag Mpc p
	bfm_f.axpy(mmp,p,mmp,mass[0]); //mmp = p*mass[0]+mmp
  
	d += rn*mass[0];
	b = -cp /d;
	if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: b = -cp/d = -%le/%le = %le\n",cp,d,b);

	// Set up the various shift variables
	int       iz=0;

	z[0][1-iz] = 1.0;
	z[0][iz]   = 1.0;
	bs[0]      = b;
	for(int s=1;s<nshift;s++){
	    z[s][1-iz] = 1.0;
	    z[s][iz]   = 1.0/( 1.0 - b*(mass[s]-mass[0]));
	    bs[s]      = b*z[s][iz]; // Sign relative to Mike - FIXME
	}

	c= bfm_f.axpy_norm(r,mmp,r,b);
	if ( bfm_f.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: k=0 residual %le \n",c);

	for(int s=0;s<nshift;s++) {
	    bfm_d.axpby(psi[s],src,src,0.,-bs[s]*alpha[s]); //initialize double prec solutions
	}
  
	// Iteration loop
	for (int k=1;k<=bfm_f.max_iter;k++){
	    a = c /cp;

#define CK_BAGEL_OPTIMISE

#ifndef CK_BAGEL_OPTIMISE
	    mixed_cg::threaded_convFermion_fast(tmp_d, r, bfm_d, bfm_f); //store double prec copy of r in tmp_d

	    bfm_d.axpy(p_d,p_d,tmp_d,a);

	    for(int s=0;s<nshift;s++){
		if ( ! converged[s] ) {
		    if (s==0){
			bfm_d.axpy(ps_d[s],ps_d[s],tmp_d,a);
		    } else{
			double as =a *z[s][iz]*bs[s] /(z[s][1-iz]*b);
			bfm_d.axpby(ps_d[s],tmp_d,ps_d[s],z[s][iz],as); //ps_d[s] = z[s][iz]*tmp_d + as*ps_d[s]
		    }
		}
	    }
#else
	    //Note, I moved the update of the search vectors to further down so it can potentially be combined with the solution vector update
	    
	    double as_uc[nshift+1], z_uc[nshift+1];
	    Fermion_t ps_d_unconv[nshift+1];
	    int nunconv=0;
	    for(int s=0;s<nshift;s++)
		if ( ! converged[s] ) {
		    ps_d_unconv[nunconv] = ps_d[s];
		    z_uc[nunconv] = z[s][iz];

		    if(s==0) as_uc[nunconv] = a;
		    else as_uc[nunconv] = a *z[s][iz]*bs[s] /(z[s][1-iz]*b);

		    ++nunconv;
		}
	    bfm_d.axpy_sy(p_d,p_d,r,a); //r vector in single precision
#endif

	    cp=c;
    
	    mixed_cg::threaded_convFermion_fast(p, p_d, bfm_f, bfm_d);

	    d= bfm_f.Mprec(p,mp,tmp,DaggerNo,1); 
	    bfm_f.Mprec(mp,mmp,tmp,DaggerYes);
	    bfm_f.axpy(mmp,p,mmp,mass[0]);

	    double rn = bfm_f.norm(p);

	    d += rn*mass[0];
	    bp=b;
	    b=-cp/d;
    
	    // Toggle the recurrence history
	    bs[0] = b;
	    iz = 1-iz;

	    for(int s=1;s<nshift;s++){
		if(!converged[s]){
		    double z0 = z[s][1-iz];
		    double z1 = z[s][iz];
		    z[s][iz] = z0*z1*bp
			/ (b*a*(z1-z0) + z1*bp*(1- (mass[s]-mass[0])*b)); 
		    bs[s] = b*z[s][iz]/z0; // NB sign  rel to Mike
		}
	    }
#define CK_BAGEL_OPTIMISE_COMBINE_PSI_PS


#ifndef CK_BAGEL_OPTIMISE_COMBINE_PSI_PS

#  ifdef CK_BAGEL_OPTIMISE //Update the search vectors here rather than above
	    bfm_d.axpby_multi_reusey(ps_d_unconv,ps_d_unconv,r,as_uc,z_uc,nunconv,1);
#  endif
	    for(int s=0;s<nshift;s++){
		int ss = s;
		if(!converged[s])
		    bfm_d.axpy(psi[ss],ps_d[s],psi[ss],-bs[s]*alpha[s]);
	    }
#else
	    //CK_BAGEL_OPTIMISE_COMBINE_PSI_PS, combine the above steps
	    double c_uc[nunconv];
	    Fermion_t psi_d_unconv[nunconv];
	    int off=0;
	    for(int s=0;s<nshift;s++)
		if(!converged[s]){
		    c_uc[off] = -bs[s]*alpha[s];
		    psi_d_unconv[off++] = psi[s];
		}

	    bfm_d.cgmulti_update_srch_sol(psi_d_unconv,ps_d_unconv, r, as_uc, z_uc , c_uc ,nunconv,1);

#endif

	    //Reliable update
	    if(k % update_freq == 0){
		double c_sp = bfm_f.axpy_norm(r,mmp,r,b);

		//Replace r with true residual
		mixed_cg::switch_comm(bfm_d, bfm_f);

		bfm_d.Mprec(psi[0],mp_d,tmp_d,0,1);
		bfm_d.Mprec(mp_d,mmp_d,tmp_d,1);
		bfm_d.axpy(mmp_d,psi[0],mmp_d,mass[0]);

		c = bfm_d.axpy_norm(tmp_d,mmp_d,src,-1.0);

		if( bfm_d.isBoss() && !me) printf("bfmbase::CGNE_prec_multi: reliable update iter %d, replaced |r|^2 = %.12le with |r|^2 = %.12le\n",k,c_sp,c);

		mixed_cg::threaded_convFermion_fast(r, tmp_d, bfm_f, bfm_d);
		mixed_cg::switch_comm(bfm_f, bfm_d);
	    }else{
		c= bfm_f.axpy_norm(r,mmp,r,b);
	    }

	    // Convergence checks
	    int all_converged = 1;
	    if ( ((k%100)==0) && bfm_f.isBoss() && (!me) ) printf("bfmbase::CGNE_prec_multi: k=%d c=%g, shift in current dir for lightest pole %.12e\n",k,c,-bs[0]*alpha[0]);
	    for(int s=0;s<nshift;s++){
		if (!converged[s]){
		    double css  = c * z[s][iz]* z[s][iz];	
		    if(css<rsq[s]) converged[s]=1;
		    else all_converged=0;
		    if(bfm_f.isBoss() && (!me) && converged[s]) printf("bfmbase::CGNE_prec_multi: Shift %d converged on iter %d: test cur %g, targ %g   [Stated true resid %g].\n",s,k,css,rsq[s],sqrt(css/rsq[s])*mresidual[s]);
		    else if ( ((k%100)==0) && bfm_f.isBoss() && (!me) ) printf("bfmbase::CGNE_prec_multi: Shift %d convergence test cur %g, targ %g   [Stated true resid %g].\n",s,css,rsq[s],sqrt(css/rsq[s])*mresidual[s]);	
		}
	    }
	    if(converged[0] && !all_converged){
		if(bfm_f.isBoss() && !me) printf("bfmbase::CGNE_prec_multi: WARNING, shift[0] has converged but not all higher mass poles have. Algorithm ending here!\n");
		all_converged = 1;
	    }

	    if ( all_converged ){
		if ( bfm_f.isBoss() && (!me) )printf("bfmbase::CGNE_prec_multi: k=%d All shifts have converged\n",k);
		if ( bfm_f.isBoss() && (!me) )printf("bfmbase::CGNE_prec_multi: k=%d Checking solutions\n",k);

		// Check answers
		mixed_cg::switch_comm(bfm_d, bfm_f);

		for(int s=0; s < nshift; s++) {
		    //Convert solution to double precision
		    bfm_d.Mprec(psi[s],mp_d,tmp_d,DaggerNo);
		    bfm_d.Mprec(mp_d,mmp_d,tmp_d,DaggerYes);
		    bfm_d.axpy(tmp_d,psi[s],mmp_d,mass[s]);
		    bfm_d.axpy(mp_d,tmp_d,src,-1);
		    double rn = bfm_d.norm(mp_d);
		    double cn = bfm_d.norm(src);
		    if ( bfm_d.isBoss() && !me ) {
			printf("double prec final: shift[%d] true residual %.12le \n",s,sqrt(rn/cn));
		    }
		}

		if ( single ) {
		    for(int s=1; s < nshift; s++) { 
			bfm_d.axpy(psi[0],psi[s],psi[0],1.0);
		    }      
		}

		DEALLOCATE_ALL;

		return k;

	    }else if(report_freq != -1 && k % report_freq == 0){
		mixed_cg::switch_comm(bfm_d, bfm_f);
		for(int s=0; s < nshift; s++) {
		    double css  = c * z[s][iz]* z[s][iz];
		    bfm_d.Mprec(psi[s],mp_d,tmp_d,DaggerNo);
		    bfm_d.Mprec(mp_d,mmp_d,tmp_d,DaggerYes);
		    bfm_d.axpy(tmp_d,psi[s],mmp_d,mass[s]);
		    bfm_d.axpy(mp_d,tmp_d,src,-1);
		    double rn = bfm_d.norm(mp_d);
		    double cn = bfm_d.norm(src);
		    if ( bfm_d.isBoss() && !me ) {
			printf("iter %d, double prec: shift[%d] true residual %.12le, running true residual %.12le [converged = %d]\n",k,s,sqrt(rn/cn),sqrt(css/rsq[s])*mresidual[s],converged[s]);
		    }
		}
		mixed_cg::switch_comm(bfm_f, bfm_d);
	    }
	}

	mixed_cg::switch_comm(bfm_d, bfm_f);
	if ( bfm_d.isBoss() && !me ) printf("bfmbase::CGNE_prec_multi: CG not converged \n");

	DEALLOCATE_ALL;

	return -1;
    }
#undef DEALLOCATE_ALL

    //This version has the following steps:
    //1) Single precision multi-mass solve with reliable update and double precision shift vectors
    //2) Single precision restarted CG with defect correction loop over poles
    //3) Double precision restarted CG with defect correction loop over poles

    inline int threaded_cg_mixed_multi_shift_MdagM_sp_relup_dp_defect_correction(Fermion_t psi[], 
										 Fermion_t src,
										 double    mass[],
										 double    alpha[],
										 int       nshift,
										 double mresidual[],
										 int single,
										 bfm_evo<float> &bfm_f,
										 bfm_evo<double> &bfm_d,
										 int update_freq = 100,
										 int report_freq = -1,
										 int max_cycle = 10){

	int me = bfm_d.thread_barrier();
	double frsd = bfm_f.residual; //save original residual for later restoration

	if (bfm_d.isBoss() && !me) printf("threaded_cg_mixed_multi_shift_MdagM_sp_relup_dp_defect_correction: hello!!!!! nshift = %d\n", nshift);

	// FIXME: Greg test: tighten residuals for all but lightest pole by a factor of ten during single precision solve
	//for (int i = 1; i < nshift; i++) mresidual[i] /= 10.0;

	struct timeval tstart, tstop, tdiff;
	gettimeofday(&tstart,NULL);
	int iter_multi = threaded_cg_mixed_multi_shift_MdagM_sp_relup_dp(psi,src,mass,alpha,nshift,mresidual,0,bfm_f,bfm_d,update_freq,report_freq);
	gettimeofday(&tstop,NULL);
	timersub(&tstop,&tstart,&tdiff);

	// FIXME: Greg test: tighten residuals for all but lightest pole by a factor of ten during single precision solve
	//for (int i = 1; i < nshift; i++) mresidual[i] *= 10.0;

	if (bfm_d.isBoss() && !me) printf("threaded_cg_mixed_multi_shift_MdagM_sp_relup_dp_defect_correction: Initial multi-shift iter = %d, time %d.%6.6d s\n", iter_multi, tdiff.tv_sec, tdiff.tv_usec);
	
	gettimeofday(&tstart,NULL);

	Fermion_t src_f = bfm_f.threadedAllocFermion();
	Fermion_t sol_f = bfm_f.threadedAllocFermion();

	Fermion_t src_d = bfm_d.threadedAllocFermion();

	Fermion_t tv1_d = bfm_d.threadedAllocFermion(mem_fast);
	Fermion_t tv2_d = bfm_d.threadedAllocFermion(mem_fast);

	double src_norm = bfm_d.norm(src);

	int iter = 0;
	for(int shift=0;shift<nshift;shift++){
	    double stop = src_norm * mresidual[shift]*mresidual[shift];
	    bfm_f.thread_barrier(); //ensure all thread writes to bfm_f.residual from previous iteration have completed
	    bfm_f.residual = mresidual[shift];

	    for(int i = 0; i < max_cycle; ++i) {
		// compute double precision rsd and also new RHS vector.
		bfm_d.Mprec(psi[shift],   tv1_d, src_d, 0, 0); //here src_d is just used as a temp storage
		bfm_d.Mprec(tv1_d, tv2_d, src_d, 1, 0); // tv2_d = MdagM * sol
		bfm_d.axpy(tv2_d,psi[shift],tv2_d,mass[shift]); //tv2_d = (MdagM + shift)* sol

		double norm = bfm_d.axpy_norm(src_d, tv2_d, src, -1.);

		// Hantao's ad hoc stopping condition
		if(norm < 100. * stop){
		    if(bfm_f.isBoss() && !me) {
			printf("threaded_cg_mixed_multi_shift_MdagM_sp_relup_dp_defect_correction: shift = %d needs no correction: rsd = %17.10e(d) stop = %17.10e(d)  [True resid %17.10e(d)]\n",
			       shift, norm, stop, sqrt(norm/src_norm));
		    }
		    break;
		}

		if(!me)
		    while(norm * bfm_f.residual * bfm_f.residual < stop) bfm_f.residual *= 2;
		//bfm_f.thread_barrier(); //No need, next call has a barrier

		threaded_convFermion_fast(src_f, src_d, bfm_f, bfm_d);
		switch_comm(bfm_f, bfm_d);

		if(bfm_f.isBoss() && !me) {
		    printf("threaded_cg_mixed_multi_shift_MdagM_sp_relup_dp_defect_correction: shift = %d, defect correction cycle = %d rsd = %17.10e(d) stop = %17.10e(d)  [True resid %17.10e(d), next single prec target resid %17.10e]\n",
			   shift,i, norm, stop, sqrt(norm/src_norm), bfm_f.residual );
		}

		bfm_f.set_zero(sol_f);
		iter += threaded_CGNE_MdagM_plus_shift<float>(sol_f,src_f,mass[shift],bfm_f);

		switch_comm(bfm_d, bfm_f);
		threaded_convFermion_fast(tv1_d, sol_f, bfm_d, bfm_f);

		bfm_d.axpy(psi[shift], tv1_d, psi[shift], 1.);
	    }
	    bfm_f.residual = frsd; //restore original single precision residual at end of each step
	}

	gettimeofday(&tstop,NULL);
	timersub(&tstop,&tstart,&tdiff);
	if(bfm_d.isBoss() && !me) printf("threaded_cg_mixed_multi_shift_MdagM_sp_relup_dp_defect_correction: defect correction time %d.%6.6d s\n", tdiff.tv_sec,tdiff.tv_usec);

	gettimeofday(&tstart,NULL);

	for(int shift=0;shift<nshift;shift++){
	    if(bfm_d.isBoss() && !me){ printf("threaded_cg_mixed_multi_shift_MdagM_sp_relup_dp_defect_correction: doing final inversion for shift %d using corrected solution as guess\n",shift); fflush(stdout); }
	    bfm_d.thread_barrier(); //ensure writes to bfm_d.residual from previous iteration have completed
	    double restore_resid = bfm_d.residual;
	    bfm_d.thread_barrier(); //ensure all threads have the same value
	    bfm_d.residual = mresidual[shift];

	    iter += threaded_CGNE_MdagM_plus_shift<double>(psi[shift],src,mass[shift],bfm_d);
	    bfm_d.residual = restore_resid;
	}

	if ( single ) {
	    for(int s=1; s < nshift; s++) { 
		bfm_d.axpy(psi[0],psi[s],psi[0],1.0);
	    }      
	}

	bfm_d.threadedFreeFermion(src_d);
	bfm_f.threadedFreeFermion(src_f);

	bfm_d.threadedFreeFermion(tv1_d);
	bfm_d.threadedFreeFermion(tv2_d);

	bfm_f.threadedFreeFermion(sol_f);

	gettimeofday(&tstop,NULL);
	timersub(&tstop,&tstart,&tdiff);
	if(bfm_d.isBoss() && !me) printf("threaded_cg_mixed_multi_shift_MdagM_sp_relup_dp_defect_correction: finishing up time %d.%6.6d s\n", tdiff.tv_sec,tdiff.tv_usec);

	return iter;
    }
};
#endif


#endif
