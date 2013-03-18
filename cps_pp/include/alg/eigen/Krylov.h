#ifndef KRYLOV_H
#define KRYLOV_H

#include <bfm_qdp.h>

#include <omp.h>
#include "Deflate.h"
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <time.h>

#include <util/lattice/bfm_eigcg.h>
#include <alg/lanc_arg.h>
#include <util/time_cps.h>

//Remove dwfa_setup
//Remove save5d, prectype, M5, mass
//Remove lx, ly, lz, lt, Ls
//Remove u
//Remove  int Arnoldi_Factor(int start, int end);
//Fix the quit failure problem

//enum EigenType {D, DDAG, G5D, DDAGD}; 
using cps::EigenType; 
using cps::D; 
using cps::DDAG; 
using cps::G5D; 
using cps::DDAGD; 
using cps::dclock;

/** Different ways to sort **/
template <class T> inline bool cf_large(T i, T j) { return (abs(i) < abs(j)); }
template <class T> inline bool cf_small(T i, T j) { return (abs(i) > abs(j)); }

/** Find the largest/smallest eigenvalues **/
enum ordertype {small, large};	
/** Doing arnoldi or lanczos ? **/
enum krylovtype {arn, lan};


/** A sorting algorithm so bad it makes small children cry.
    keeps track of which eval goes with which evec		**/
template <class T> 
void inline EigenSort(vector<T> &evals, vector<vector<T> > &evecs, ordertype comp)
{
    int N = evals.size();
    vector<T> unsorted(N); 
    unsorted = evals;

    if(comp == small) 
        sort(evals.begin(), evals.end(), cf_small<T>);
    if(comp == large) 
        sort(evals.begin(), evals.end(), cf_large<T>);

    vector<vector<T> > sortedvecs(N);
    vector<int> derange(N);

    for(int j = 0; j < N; ++j)
	{
            double min = abs(evals[j] - unsorted[0]); 
            derange[j] = 0;
            for(int i = 1; i < N; ++i)
		{
                    double tmin = abs(evals[j] - unsorted[i]);
                    if(tmin < min)
			{
                            derange[j] = i; 
                            min = tmin;
			}	
		}
	}

    for(int i = 0; i < N; ++i)
        sortedvecs[i] = evecs[derange[i]];
    evecs = sortedvecs;
}

template <class T> 
void EigenSort(multi1d<Complex> &vals, multi1d<LatticeFermion> &evecs)
{
    int N = vals.size();
    vector<T> evals(N); 
    for(int i = 0; i < N; ++i)
        evals[i] = toDouble( real(vals[i]) );
    vector<T> unsorted(N); unsorted = evals;

    sort(evals.begin(), evals.end(), cf_small<T>);

    multi1d<LatticeFermion> sortedvecs(N);
    vector<int> derange(N);

    for(int j = 0; j < N; ++j)
	{
            double min = abs(evals[j] - unsorted[0]); 
            derange[j] = 0;
            for(int i = 1; i < N; ++i)
		{
                    double tmin = abs(evals[j] - unsorted[i]);
                    if(tmin < min)
			{
                            derange[j] = i; 
                            min = tmin;
			}	
		}
	}

    for(int i = 0; i < N; ++i)
        sortedvecs[i] = evecs[derange[i]];
    for(int i = 0; i < N; ++i)
        evecs[i] = sortedvecs[i];

    multi1d<Complex> sorted(N);
    for(int i = 0; i < N; ++i)
        sorted[i] = vals[derange[i]];
    for(int i = 0; i < N; ++i)
        vals[i] = sorted[N - i - 1];
}


template <class T> 
void EigenSort(multi1d<Complex> &vals, multi1d<multi1d<LatticeFermion> > &evecs)
{

    int N = vals.size();
    int Ls = evecs[0].size();

    vector<T> evals(N); for(int i=0;i<N;i++){evals[i] = toDouble( real(vals[i]) );}
    vector<T> unsorted(N); unsorted = evals;

    sort(evals.begin(),evals.end(),cf_small<T>);

    multi1d<multi1d<LatticeFermion> > sortedvecs(N); for(int j=0;j<N;j++){sortedvecs[j].resize(Ls);}
    vector<int> derange(N);

    for(int j=0;j<N;j++){
        double min = abs(evals[j] - unsorted[0]); derange[j] = 0;
        for(int i=1;i<N;i++){
            double tmin = abs(evals[j] - unsorted[i]);
            if(tmin < min){derange[j] = i; min = tmin;}	
        }
    }

    for(int i=0;i<N;i++){
        for(int s = 0;s<Ls;s++){
            sortedvecs[i][s] = evecs[derange[i]][s];}}
    for(int i=0;i<N;i++){
        for(int s = 0;s<Ls;s++){evecs[i][s] = sortedvecs[i][s];}}

    multi1d<Complex> sorted(N);
    for(int i=0;i<N;i++){sorted[i] = vals[derange[i]];}
    for(int i=0;i<N;i++){vals[i] = sorted[N-i-1];}

}

/** An abstract iterative Krylov eigensolver.  **/
template <class S> 
class Krylov : public bfm_internal<S> 
{
public:

    int M;  ///Dimension M Krylov space
    int K;  ///Want K converged vectors
    int get; //Actually number of eigen vectors you will get
    int P;  ///Num shifts

    int con; ///number converged so far
    int maxits; ///maxiterations
    bool stop;	///converged
    double conv;  ///Convergence
    double convQR;	///convergence of intermediate QR solves
    ordertype comp; ///Type of evals wanted
    Matrix<S> H;	///Hessenberg matrix
    krylovtype kr; ///Lanczos or Arnoldi
    int mvprod; ///Counts matrix vector products
    int innerprod; ///Counts inner products
    bool lock; ///Use locking transofrmation or not
    int iterations;///counts its
    bool cont;	///Continuation of pervious calc
    bool refine;	///Improve previous calc
    EigenType D;	///The operator type

    double rho;   ///orthagonality requirement for Lanczos
    int ref;      ///orthagonality requirement for Lanczos

    int prec; ///preconditioned = 1 or not = 0

    int Tn;	///Order of Chebyshev polynomial
    bool sh; 	///Shifting or not
    double ch_rad; ///Spectral radius
    double ch_off; ///Spectral offset (ie. find eigenvalues of magnitude less than this)
    double ch_shift; ///Shift the peak

    multi1d<S> shift_extra;

    bfm_evo<S> &dop;

    ///Basis and residual
    multi1d<bfm_fermion> bq;
    multi1d<bfm_fermion> bf;

    ///Evals and evecs of Hessenberg/Tridiagonal matrix
    vector<S> evals;
    vector<vector<S> > evecs;
    //eigen vaules 
    multi1d<S> bl;

    Krylov(bfm_evo<S> &dwf);
    Krylov(bfm_evo<S> &dwf, cps::LancArg &lanc_arg);
    virtual ~Krylov();
    void Run();
    void ObliqueProj(bfm_fermion b, bfm_fermion x);

    //time stat
    double dslash_time;
    double tot_time;
    double shift_time;

protected:
    ///Useful temporaries
    //compact fermion
    Fermion_t tmp1[2], tmp2[2]; 
    //fermion, not compact
    Fermion_t tmp;
    Fermion_t invec[2], outvec[2];
    bool initialized;

    void Krylov_init();
    virtual void init() = 0;
    void Check();
    void Reset();
    void ImplicitRestart(	int TM, vector<S> &evals,  vector<vector<S> > &evecs);
    void TestConv(int SS, vector<S> &tevals, vector<vector<S> > &tevecs);
    void True_Residual( void );
    int Lanczos_Factor(int start, int end);
    void G5R(bfm_fermion result, bfm_fermion input);
    void Abort(int ff, vector<S> &evals,  vector<vector<S> > &evecs);
    void Shift(Matrix<S> &Q, vector<S> shifts);

    virtual void herm_mult(bfm_fermion input, bfm_fermion &result) =0;
    virtual void dwf_multiply(bfm_fermion input, bfm_fermion &result) =0;
    void Gram(bfm_fermion r, multi1d<bfm_fermion> &q, int N);
    void times_real(multi1d<bfm_fermion> &q, Matrix<S> &Q, int N);
    void times(multi1d<bfm_fermion> &q, Matrix<S> &Q, int N);

    void axpy(bfm_fermion r, double a, bfm_fermion x, bfm_fermion y);
    double axpy_norm(bfm_fermion r, double a, bfm_fermion x, bfm_fermion y);
    void axpby(bfm_fermion r, double a, bfm_fermion x, double b, bfm_fermion y);
    double axpby_norm(bfm_fermion r, double a, bfm_fermion x, double b, bfm_fermion y);

    void bfm_to_qdp(Fermion_t in[2], multi1d<LatticeFermion> &out);
    void bfm_to_qdp(Fermion_t in[2], LatticeFermion &out);
    void qdp_to_bfm(multi1d<LatticeFermion> &in, Fermion_t out[2]);
    void qdp_to_bfm(LatticeFermion in, Fermion_t out[2]);
    void free_fermion(bfm_fermion x);
    void zero_fermion(bfm_fermion x);
    void init_fermion(bfm_fermion x);
    void equate(bfm_fermion x, bfm_fermion y);
    void scale(bfm_fermion x, S s); 
    void scale(bfm_fermion x, complex<double> s); 
    double norm(bfm_fermion x);
    double norm2(bfm_fermion x);
    void gamma5(bfm_fermion x); 
    double innerProduct_real(bfm_fermion x, bfm_fermion y);
    complex<double> innerProduct(bfm_fermion x, bfm_fermion y);
    double Chebyshev(double x);
    void Chebyshev(bfm_fermion input, bfm_fermion &v1);
    double qpoly(double alpha, double beta, double x);
    void qpoly(double beta, double alpha, bfm_fermion input, bfm_fermion &result);
    void shift_multiply(bfm_fermion input, bfm_fermion &result);
    void multiply(bfm_fermion input, bfm_fermion &result);
};

/// Constructor
template <class S> 
Krylov<S>::Krylov(bfm_evo<S> &dwf): dop(dwf)
{
    con = 0;
    conv = 1e-10;	//TODO residual underestimates true error, converge more accurately, usually fine.
    convQR = 1e-14;
    maxits = 10000;
    comp = small;
    stop = false;
    srand48(123456789);
    mvprod = 0;
    innerprod = 0;
    lock = true;
    iterations = 0;
    cont = false;
    refine = false;
    initialized = false;
    rho = 1.0;
    ref=0;
    sh = false;

    shift_extra.resize(0);

    //origninally in dwfa_setup
    //dwfa.residual = conv/100.0;  ///Do PV inverse to this accuracy
    //origninally in dwfa_setup end

}

/// Constructor
template <class S> 
Krylov<S>::Krylov(bfm_evo<S> &dwf, cps::LancArg &lanc_arg): dop(dwf)
{
    con = 0;
    conv = 1e-10;	//TODO residual underestimates true error, converge more accurately, usually fine.
    convQR = 1e-14;
    maxits = 10000;
    comp = small;
    stop = false;
    srand48(123456789);
    mvprod = 0;
    innerprod = 0;
    iterations = 0;
    cont = false;
    refine = false;
    initialized = false;
    rho = 1.0;
    ref=0;

    shift_extra.resize(0);

    conv = lanc_arg.stop_rsd;
    convQR = lanc_arg.qr_rsd;
    M = lanc_arg.N_use;
    K = lanc_arg.N_get;	
    get = lanc_arg.N_true_get;
    Tn = lanc_arg.ch_ord;
    ch_off = lanc_arg.ch_beta;
    ch_rad = lanc_arg.ch_alpha;
    sh = lanc_arg.ch_sh;
    ch_shift = lanc_arg.ch_mu;
    prec = lanc_arg.precon;
    lock = lanc_arg.lock;
    D = lanc_arg.EigenOper;	  	
    maxits = lanc_arg.maxits;

    dop.mass = lanc_arg.mass;
    //for mobius
    dop.GeneralisedFiveDimEnd();
    dop.GeneralisedFiveDimInit();
    //origninally in dwfa_setup
    //dwfa.residual = conv/100.0;  ///Do PV inverse to this accuracy
    //origninally in dwfa_setup end
	
    dslash_time = 0.0;
    tot_time = 0.0;
    shift_time = 0.0;
}

///Basic initialization
template <class S> void Krylov<S>::Krylov_init(){

    this->P = this->M-this->K;
    this->H.resize(this->M); this->H.Fill(0.0);
    this->evals.resize(this->M);
    this->evecs.resize(this->M);
    if(this->sh) this->ch_rad += this->ch_shift;
    this->conv = this->conv/100.0;
    (this->Tn > 0) ? this->comp = large : this->comp = small;
}

///Destructor
template <class S> 
Krylov<S>::~Krylov()
{
    //release memory for eigenvectors
    for(int i = 0; i < get; i++)
        this->free_fermion(bq[i]);
}

///Run the Eigensolver
template <class S> 
void Krylov<S>::Run()
{
    tot_time -= dclock();
    init();
    int ff = this->con;

    QDPIO::cout << "Krylov: Run()"<<endl;

    if(this->kr == lan) 
	{
            QDPIO::cout << "Krylov: Lanczos_Factor() ff = "<< ff << " " << this->M << endl;
            ff = Lanczos_Factor(ff, this->M);
	}

    if(ff < this->M)
	{
            QDPIO::cout << "Krylov: aborting ff "<<ff <<" "<<this->M<<endl;
            Abort(ff, this->evals, this->evecs);
	}

    int itcount = 0;
    bool stop = false;

    for(int it = 0; it < this->maxits && (this->con < this->get); ++it)
	{
            this->iterations++;
            itcount++;	
            QDPIO::cout << "Krylov: Iteration --> " << it << endl;
            int lock_num = this->lock ? this->con : 0;
            std::vector<S> tevals(this->M - lock_num );
            std::vector<std::vector<S> > tevecs(this->M - lock_num);

            //chekc true residual
            // True_Residual();
            //check residual of polynominal 
            TestConv(this->M, tevals, tevecs);
            if(this->con >= this->get)
                break;
            ImplicitRestart(ff, tevals,tevecs);
	}

    if(this->kr == lan)
        SymmEigensystem<S>(this->H, this->evals, this->evecs, this->convQR);
    Check();
	
    tot_time += dclock();

    QDPIO::cout << "Done : "
		<< "Total time = "<< tot_time << "sec, "
		<< "Dlash time = "<< dslash_time << "sec, "
		<< "Shift time = "<< shift_time << "sec. "
		<<endl;
}

///H - shift I = QR; H = Q* H Q
template <class S> void Krylov<S>::Shift(Matrix<S> &Q, vector<S> shifts)
{
    shift_time -= dclock();

    int P = shifts.size();
    int M = Q.dim;
    Q.Unity();
    int lock_num = this->lock ? this->con : 0;

    for(int i=0;i<P;i++){

        S x, y, z;
        vector<S> ck(3), v(3);

        x = H(lock_num+0,lock_num+0)-shifts[i];
        y = H(lock_num+1,lock_num+0);
        ck[0] = x; ck[1] = y; ck[2] = 0; 

        normalize(ck);	///Normalization cancels in PHP anyway
        S beta;

        Householder_vector<S>(ck, 0, 2, v, beta);
        Householder_mult<S>(H,v,beta,0,lock_num+0,lock_num+2,0);
        Householder_mult<S>(H,v,beta,0,lock_num+0,lock_num+2,1);
        ///Accumulate eigenvector
        Householder_mult<S>(Q,v,beta,0,lock_num+0,lock_num+2,1);
        int sw = 0;

        for(int k=lock_num+0;k<M-2;k++){

            x = H(k+1,k); 
            y = H(k+2,k); 
            z = (S)0.0;
            if(k+3 <= M-1){
                z = H(k+3,k);
            }else{
                sw = 1; v[2] = 0.0;
            }

            ck[0] = x; ck[1] = y; ck[2] = z;

            normalize(ck);

            Householder_vector<S>(ck, 0, 2-sw, v, beta);
            Householder_mult<S>(H,v, beta,0,k+1,k+3-sw,0);
            Householder_mult<S>(H,v, beta,0,k+1,k+3-sw,1);
            ///Accumulate eigenvector
            Householder_mult<S>(Q,v, beta,0,k+1,k+3-sw,1);
        }
    }
	
    shift_time += dclock();
}

///Converge on invariant subspace and quit "elegantly"
template <class S> void Krylov<S>::Abort(int ff, vector<S> &evals,  vector<vector<S> > &evecs)
{
    QDPIO::cout << "Check for converged eigenvectors before aborting" << endl;
    QDPIO::cout << "Subspace dimension " << ff << " < " << this->M << endl; 

    int lock_num = this->lock ? this->con : 0;
    std::vector<S> tevals(ff-lock_num);
    std::vector<std::vector<S> > tevecs(ff-lock_num);

    TestConv(ff, tevals, tevecs);
    Matrix<S> AH(ff); 
    AH = this->H.GetSubMtx(0,ff,0,ff);

    this->H.resize(ff);
    this->H = AH;

    evals.resize(ff);
    evecs.resize(ff);
    this->K = this->con;
    this->M = ff;
}

///result = gamma5 R_5 input.
template <class S> 
void Krylov<S>::G5R(bfm_fermion result, bfm_fermion input)
{
    int Ls = dop.Ls;
    int x[4];
    int Nspinco = 12;
    int cb0, cb1, site, bidx, rb_idx;
    S sgn;
    // PAB
    // Very inefficient
    //for(int s=0;s<this->Ls;s++){
    for(int s=0;s<Ls;s++){
        for ( x[3]=0; x[3]<dop.node_latt[3];x[3]++ ) { 
            for ( x[2]=0; x[2]<dop.node_latt[2];x[2]++ ) { 
                for ( x[1]=0; x[1]<dop.node_latt[1];x[1]++ ) { 
                    for ( x[0]=0; x[0]<dop.node_latt[0];x[0]++ ) { 

                        site = x[0]+x[1]+x[2]+x[3];   
                        if(dop.precon_5d == 1){
                            cb0 = ((site+s)&0x1);
                            cb1 = ((site+Ls-1-s)&0x1);}
                        else{
                            cb0 = ((site)&0x1);
                            cb1 = ((site)&0x1);}

                        for ( int co=0;co<Nspinco;co++ ) { 
                            for ( int reim=0;reim<2;reim++ ) { 

                                bidx = dop.bagel_idx5d(x,Ls-1-s,reim,co,Nspinco,1);
                                rb_idx = dop.bagel_idx5d(x,s,reim,co,Nspinco,1);	///TODO Must be an easier way...
                                sgn = 1.0;
                                if(co > 5) sgn = -1.0;


                                if(this->prec == 0){
                                    S * forward = (S *)input[cb1];
                                    S * backward = (S *)result[cb0];
                                    backward[rb_idx] = sgn * forward[bidx]; 
                                }
                                else{
                                    S * forward = (S *)input[1];
                                    S * backward = (S *)result[1];
                                    backward[rb_idx] = sgn * forward[bidx]; 
                                }

                            }}
                    }}}}
    }

}

template <class S> int Krylov<S>::Lanczos_Factor(int start, int end){

    S beta;  
    complex<S> alpha;

    QDPIO::cout<<"Lanczos_Factor start/end " <<start <<"/"<<end<<endl;
    ///Starting from scratch, bq[0] contains a random vector and |bq[0]| = 1
    if(start == 0){
        start++;
        Chebyshev(this->bq[0],this->bf[0]);						//bf = A bq[0]
        alpha = this->innerProduct(this->bq[0],this->bf[0]);this->innerprod++;	//alpha =  bq[0]^dag A bq[0]
        this->axpy(this->bf[0] ,-1.0*real(alpha), this->bq[0], this->bf[0]);	//bf =  A bq[0] - alpha bq[0]
        this->H( real(alpha),0,0);							//bq[0] . bf =  bq[0]^dag A bq[0] - alpha |bq[0]|^2 = 0
    }
    int re = 0;
    start--;

    beta = this->axpy_norm(this->bf[0],0.0,this->bf[0],this->bf[0]);		//|bf|^2
    S sqbt = sqrt(beta);

    if( this->cont){beta = 0;sqbt = 0;}	
    for(int j=start+1;j<end;j++){
        QDPIO::cerr << "Factor j " << j+1 << endl;
        if(this->cont){
            this->axpy(this->bq[j] , 0.0 , this->bf[0], this->bf[0]);this->cont = false;
        }
        else{
            this->axpby(this->bq[j] , (1.0/sqbt) , this->bf[0], 0.0, this->bf[0]);	//bq = bf/|bf|
            this->H( sqbt,j,j-1); this->H( sqbt,j-1,j);
        }

        Chebyshev(this->bq[j],this->bf[0]);									//bf = A bq[j]
        this->axpy(this->bf[0], -1.0*sqbt, this->bq[j-1], this->bf[0]);					//bf = A bq[j] - beta bq[j-1]
        alpha = this->innerProduct(this->bq[j],this->bf[0]); this->innerprod++;				//alpha = bq[j]^dag A bq[j]
        S fnorm = this->axpy_norm(this->bf[0] , -1.0*real(alpha), this->bq[j], this->bf[0]);		//bf = A bq[j] - beta bq[j-1] - alpha bq[j]

        re = 0;
        S bck = sqrt( real( conj(alpha)*alpha ) + beta );
        beta = fnorm;
        sqbt = sqrt(beta);

        ///Iterative refinement of orthogonality V = [ bq[0]  bq[1]  ...  bq[M] ]
        while( re == this->ref || (sqbt < this->rho * bck && re < 5) ){

            //bex = V^dag bf
            vector<complex<S> > bex(j+1);
            for(int k=0;k<j+1;k++){bex[k] = this->innerProduct(this->bq[k],this->bf[0]); this->innerprod++;}

            this->zero_fermion(this->tmp2);
            //tmp2 = V s
            for(int l=0;l<j+1;l++){
                S nrm = norm2(this->bq[l]);
                this->axpy(this->tmp1,0.0,this->bq[l],this->bq[l]); this->scale(this->tmp1,bex[l]); 	//tmp1 = V[j] bex[j]
                this->axpy(this->tmp2,1.0,this->tmp2,this->tmp1);					//tmp2 += V[j] bex[j]
            }

            //bf = bf - V V^dag bf.   Subtracting off any component in span { V[j] } 
            S btc = this->axpy_norm(this->bf[0],-1.0,this->tmp2,this->bf[0]);
            alpha = alpha + bex[j];	      sqbt = sqrt(real(btc));	      
            S nmbex = 0;for(int k=0;k<j+1;k++){nmbex = nmbex + real( conj(bex[k])*bex[k]  );}
            bck = sqrt( nmbex );
            re++;
        }
        this->H( real(alpha) ,j,j);
        if(re > 1) QDPIO::cout << "orthagonality refined " << re << " times" << endl;
    }
    return end;
}

/**Print out the true residuals**/
template <class S> void Krylov<S>::True_Residual( void ){

    std::vector<S> mevals(this->M);
    std::vector<std::vector<S> > mevecs(this->M);
    Wilkinson<S>(this->H, mevals, mevecs, this->convQR);
    EigenSort(mevals,mevecs,this->comp);

    bfm_fermion bS; this->init_fermion(bS); 
    bfm_fermion bv; this->init_fermion(bv);
    bfm_fermion bvv; this->init_fermion(bvv);

    this->con = 0;

    for(int i=this->M-1; i>(this->M-1-this->get) && i>=0; i--){

        this->axpby(bS,0.0,bS,0.0,bS);

        double bdiff;
        S eval_i;
        if(this->kr == lan){

            for(int j=0;j<this->bq.size();j++){
                this->axpy(bS,real(mevecs[i][j]), this->bq[j], bS);
            }

            this->herm_mult(bS,bv);
            eval_i = this->innerProduct_real(bS,bv);	///Rayleigh quotient
            double bS_norm = sqrt( this->axpy_norm(bS,0,bS,bS) );
            eval_i = eval_i/bS_norm/bS_norm;
            this->axpby(bS,0,bS,1.0/bS_norm,bS);
            double bv_norm = sign( real(eval_i) ) * sqrt( this->axpy_norm(bv,0,bv,bv) );
            this->axpby(bv,0,bv,1.0/bv_norm,bv);

            bdiff = this->axpy_norm(bvv, -1.0 , bS, bv);
        }

        bdiff = sqrt(bdiff);
        if(bdiff < 10*conv)
            {
                QDPIO::cout <<"true residual " << this->M-1-i << " = " <<  bdiff << " of "
                            << "  :  (" << real(eval_i) << "," << imag(eval_i) << ")" << endl;
                this->con++;
            }
        else break;

    }
    this->free_fermion(bS); 
    this->free_fermion(bv);
    this->free_fermion(bvv);

}

template <class S> 
void Krylov<S>::TestConv(int SS, vector<S> &tevals, vector<vector<S> > &tevecs)
{
    QDPIO::cout << "Converged " << this->con << " so far." << endl;
    int lock_num = this->lock ? this->con : 0;
    ///Active Factorization
    Matrix<S> AH(SS - lock_num );
    AH = this->H.GetSubMtx(lock_num, SS, lock_num, SS);

    //if(this->kr == arn){ QReigensystem<S>(AH, tevals, tevecs, this->convQR);}
    if(this->kr == lan)
        Wilkinson<S>(AH, tevals, tevecs, this->convQR);
    EigenSort(tevals, tevecs, this->comp);

    Real resid_nrm;
    //if(this->kr == arn) resid_nrm =  this->axpy_norm(bf[SS-1],0.0,bf[SS-1],bf[SS-1]) ;
    if(this->kr == lan) 
        resid_nrm =  this->axpy_norm(bf[0], 0., bf[0], bf[0]) ;

    if(!this->lock) this->con = 0;

    for(int i = SS - lock_num - 1; i >= SS - this->K && i >= 0; --i)
        {
            double diff = 0;
            diff = abs(tevecs[i][this->M - 1 - lock_num]) * toDouble(resid_nrm); 

            QDPIO::cout << "residual estimate " << SS-1-i << " " << diff << " of (" << real(tevals[i]) 
                        << " , " << imag(tevals[i]) << ")" << endl;

            if(diff < this->conv)
                {
                    Matrix<S> Q(this->M);

                    if(this->lock)
			{
                            bool herm = true; 
                            //if(this->kr == arn){ herm = false;}  
                            Deflate<S>::Lock(this->H, Q, tevals[i], this->con, this->convQR, SS, herm);
                            //if(this->kr == arn){
                            //  this->times(this->bq,Q,this->bq.size());
                            //  this->scale(this->bf[this->M-1] , Q(this->M-1,this->M-1) );
                            //}
                            if(this->kr == lan)
				{
                                    this->times_real(this->bq, Q, this->bq.size());
                                    this->axpby(this->bf[0], real(Q(this->M - 1, this->M - 1)), this->bf[0], 0., this->bf[0]);
				}
                            lock_num++;
			}
                    this->con++;
                    QDPIO::cerr << " converged on eval " << this->con << " of " << this->get << endl;
                }
            else break;
        }

    QDPIO::cout << "Got " << this->con << " so far " << endl;	
}

template <class S> 
void Krylov<S>::ImplicitRestart(int TM, vector<S> &evals,  vector<vector<S> > &evecs)
{
    /// Sort by smallest imaginary part or whatever
    EigenSort(evals, evecs, this->comp);
    ///Assign shifts
    if(this->K - this->con < 4) this->P = (this->M - this->K-1); //one
    //if(this->K - this->con == 1) this->P = (this->M - this->K - ( (this->M -this->K)/2 ) ); //two
    //if(this->K - this->con < (this->M -this->K)/2 ) this->P = (this->M - this->K - ( (this->M -this->K)/2 ) ); //three Best?
    //if(this->K - this->con == 1) this->P = (this->M - this->K-1); //one

    std::vector<S> shifts(this->P + this->shift_extra.size());
    for(int k = 0; k < this->P; ++k)
        shifts[k] = evals[k]; 

    /// Shift to form a new H and q
    Matrix<S> Q(TM);
    Q.Unity();
    Shift(Q, shifts);

    int ff = this->K;
    /// Shifted H defines a new K step Arnoldi factorization
    S  beta = this->H(ff, ff-1); 
    S  sig = Q(TM - 1, ff - 1);
    QDPIO::cout << "beta = " << real(beta) << " sig = " << real(sig) << endl;
    if(this->kr == lan) 
	{
            QDPIO::cout << "TM = " << TM << endl;
            /// q -> q Q
            QDPIO::cout << this->norm2(this->bq[0]) << " -- before" << endl;

            this->times_real(this->bq, Q, TM);
            QDPIO::cout << this->norm2(this->bq[0]) << " -- after " << ff << endl;
            this->axpby(bf[0], real(beta), this->bq[ff], real(sig), this->bf[0]);
            /// Do the rest of the factorization
            ff = Lanczos_Factor(ff, this->M);
	}

    if(ff < this->M)
        Abort(ff, evals, evecs);
}

template <class S> void Krylov<S>::multiply(bfm_fermion input, bfm_fermion &result){
    EigenType tmp = this->D;
    if(this->D == G5D){ this->D = DDAGD;}
    herm_mult(input, result);
    this->D = tmp;
}
template <class S> void Krylov<S>::shift_multiply(bfm_fermion input, bfm_fermion &result){
    bfm_fermion v0; this->init_fermion(v0);
    EigenType tmp = this->D;

    herm_mult(input, v0);
    this->axpy(v0, -this->ch_shift, input, v0);
    herm_mult(v0, result);
    this->axpy(result, -this->ch_shift, v0, result);

    this->free_fermion(v0);
    this->D = tmp;
}
///q(Q)= 2*( Q^2/s^2 ) - ( 1+r ) I
template <class S> void Krylov<S>::qpoly(double beta, double alpha, bfm_fermion input, bfm_fermion &result){

    if(this->sh){
        this->shift_multiply(input,result);
    }
    else{
        this->multiply(input,result);
    }
    alpha = alpha*alpha;
    beta = beta*beta;
    double amb = (alpha - beta);
    double apb = (alpha + beta);
    this->axpby(result, 2.0/amb, result, -apb/amb, input);
}
template <class S> double Krylov<S>::qpoly(double alpha, double beta, double x){
    return (2.0*x*x - (alpha + beta)) / (beta - alpha);
}

/// n order Chebyshev
template <class S> void Krylov<S>::Chebyshev(bfm_fermion input, bfm_fermion &v1){

    //if(this->kr == arn) this->dwf_multiply(input,v1); 
    if(Tn == 0){
        if(this->kr == lan) this->herm_mult(input,v1); 
        return;
    }
    bfm_fermion v0; this->init_fermion(v0); 
    this->axpy(v0,0.0,input, input);

    this->qpoly(ch_rad,ch_off,v0,v1);

    if(Tn == 1){this->free_fermion(v0); return;}
    bfm_fermion tmp; this->init_fermion(tmp);

    for(int i=1;i<Tn;i++){
        this->qpoly(ch_rad,ch_off,v1,tmp);
        this->axpby(tmp,2.0,tmp,-1.0,v0);
        this->axpy(v0,0.0,v1,v1);
        this->axpy(v1,0.0,tmp,tmp);
    }

    this->free_fermion(tmp);
    this->free_fermion(v0);
}
/// n order Chebyshev 
template <class S> double Krylov<S>::Chebyshev(double x){

    double v0 = 1;
    if(Tn == 0){return x;}
    double v1 = this->qpoly(ch_off, ch_rad, x);
    double xbar = v1;
    if(Tn == 1){return v1;}
    double tmp;

    for(int i=1;i<Tn;i++){
        tmp = 2.0*xbar*v1-1.0*v0;
        v0=v1;
        v1=tmp;
    }

    return v1;
}

template <class S> 
complex<double> Krylov<S>::innerProduct(bfm_fermion x, bfm_fermion y)
{
    complex<double> dot;
#pragma omp parallel 
    {
        complex<double> mydot = 0;
        for(int cb = this->prec; cb < 2; cb++)
            {
                mydot += dop.inner(x[cb], y[cb]);
            }
        dot = mydot;
    }
    return dot;
}

template <class S> 
double Krylov<S>::innerProduct_real(bfm_fermion x, bfm_fermion y)
{
    double nrm;
#pragma omp parallel 
    {
        double mynrm = 0;
        for(int cb = this->prec; cb < 2; cb++)
            {
                mynrm += dop.inner_real(x[cb], y[cb]);
            }
        nrm = mynrm;
    }
    return nrm;
}

template <class S> void Krylov<S>::gamma5(bfm_fermion x)
{ 
#pragma omp parallel 
    {
        for(int cb = this->prec; cb < 2; cb++)
            {
                dop.chiralproj_axpby(x[cb], x[cb], -2.0, -1, 1);
            }
    }
}

template <class S> 
double Krylov<S>::norm2(bfm_fermion x)
{
    double nrm;
#pragma omp parallel 
    {
        double mynrm = 0;
        for(int cb = this->prec; cb < 2; cb++)
            {
                mynrm += dop.norm(x[cb]);
            }
        nrm = mynrm;
    }
    return nrm;
}

template <class S> 
double Krylov<S>::norm(bfm_fermion x)
{
    return sqrt(norm2(x));
}

template <class S> 
void Krylov<S>::scale(bfm_fermion x, complex<double> s)
{ 
#pragma omp parallel 
    {
        for(int cb = this->prec; cb < 2; cb++)
            { 
                dop.scale(x[cb], real(s), imag(s));
            }
    }
}

template <class S> 
void Krylov<S>::scale(bfm_fermion x, S s)
{ 
#pragma omp parallel 
    {
        for(int cb = this->prec; cb < 2; cb++)
            { 
                dop.scale(x[cb], s);
            }
    }
}

template <class S> 
void Krylov<S>::equate(bfm_fermion x, bfm_fermion y)
{
#pragma omp parallel 
    {
        for(int cb = this->prec; cb < 2; cb++)
            { 
                dop.equate(x[cb], y[cb]);  
            }
    }
}

template <class S> 
void Krylov<S>::init_fermion(bfm_fermion x)
{
#pragma omp parallel 
    {
        for(int cb = this->prec; cb < 2; cb++)
            {
                x[cb] = dop.threadedAllocCompactFermion();
            } 
    }
}
template <class S> 
void Krylov<S>::zero_fermion(bfm_fermion x)
{
#pragma omp parallel 
    {
        for(int cb = this->prec; cb < 2; cb++)
            {
                dop.set_zero(x[cb]);
            } 
    }
}

template <class S> 
void Krylov<S>::free_fermion(bfm_fermion x)
{
#pragma omp parallel 
    {
        for(int cb = this->prec; cb < 2; cb++)
            {
                dop.threadedFreeFermion(x[cb]);
            } 
    }
}

template <class S> void Krylov<S>::qdp_to_bfm(LatticeFermion in, Fermion_t out[2]){
    for(int cb=this->prec;cb<2;cb++){
        dop.importFermion(in,out[cb],cb);
    }
}
template <class S> void Krylov<S>::qdp_to_bfm(multi1d<LatticeFermion> &in, Fermion_t out[2]){
    for(int cb=this->prec;cb<2;cb++){
        dop.importFermion(in,out[cb],cb);
    }
}
template <class S> void Krylov<S>::bfm_to_qdp(Fermion_t in[2], LatticeFermion &out){
    for(int cb=this->prec;cb<2;cb++){
        dop.exportFermion(out,in[cb],cb);
    }
}
template <class S> void Krylov<S>::bfm_to_qdp(Fermion_t in[2], multi1d<LatticeFermion> &out){
    for(int cb=this->prec;cb<2;cb++){
        dop.exportFermion(out,in[cb],cb);
    }
}

template <class S> 
double Krylov<S>::axpby_norm(bfm_fermion r, double a, bfm_fermion x, double b, bfm_fermion y)
{
    double sum = 0;
#pragma omp parallel 
    {
        double mysum = 0;
        for(int cb = this->prec; cb < 2; cb++)
            {
                mysum += dop.axpby_norm(r[cb], x[cb], y[cb], a, b);
            }
        sum = mysum;
    }
    return sum;
}

template <class S> 
void Krylov<S>::axpby(bfm_fermion r, double a, bfm_fermion x, double b, bfm_fermion y)
{
#pragma omp parallel 
    {
        for(int cb = this->prec; cb < 2; cb++)
            {
                dop.axpby(r[cb], x[cb], y[cb], a, b);
            }
    }
}

template <class S> 
double Krylov<S>::axpy_norm(bfm_fermion r, double a, bfm_fermion x, bfm_fermion y)
{
    double sum;
#pragma omp parallel 
    {
        double mysum = 0;
        for(int cb = this->prec; cb < 2; cb++)
            {
                mysum += dop.axpy_norm(r[cb], x[cb], y[cb], a);
            }
        sum = mysum;
    }
    return sum;
}

template <class S> 
void Krylov<S>::axpy(bfm_fermion r, double a, bfm_fermion x, bfm_fermion y)
{
#pragma omp parallel 
    {
        for(int cb = this->prec; cb <2 ; cb++)
            {
                dop.axpy(r[cb], x[cb], y[cb], a);
            }
    }
}
/// q -> q Q
template <class S> void Krylov<S>::times(multi1d<bfm_fermion> &q, Matrix<S> &Q, int N){

    multi1d<bfm_fermion> SS( N ); 
    for(int i=0;i<N ;i++){
        this->init_fermion(SS[i]);
    }
    bfm_fermion tmp; 
    this->init_fermion(tmp); 


    for(int j=0;j<N ;j++){
        this->axpby(SS[j],0.0,SS[j],0.0,SS[j]);
        for(int k=0;k<N ;k++){
            this->axpy(tmp,0.0,q[k],q[k]);
            this->scale(tmp,Q(k,j));
            this->axpy(SS[j],1.0,SS[j],tmp);
        }
    }

    for(int j=0;j<N ;j++){
        axpy(q[j],0.0,SS[j],SS[j]);
    }

    this->free_fermion(tmp);
    for(int j=0;j<N ;j++){
        this->free_fermion(SS[j]);
    }
}
/// Gram Schmidting.
template <class S> 
void Krylov<S>::Gram(bfm_fermion r, multi1d<bfm_fermion> &q, int N)
{
    bfm_fermion tmp; 
    this->init_fermion(tmp); 

    {
        for(int i=0;i<N;i++){
            double norm = this->axpy_norm(tmp, 0.0, q[i], q[i]);
            complex<double> sub = ( this->innerProduct(q[i],r)/norm );
            this->scale(tmp,sub);
            this->axpy(r , -1.0 , tmp, r);
        }

        double nrm = this->norm(r);
        this->axpby(r,  1.0/nrm,r, 0.0,r );
    }
    this->free_fermion(tmp);
}

template <class S> 
void Krylov<S>::ObliqueProj(bfm_fermion x, bfm_fermion b)
{
    QDPIO::cout << "sol norm 1 cpt = " << this->norm(x) << endl;
    QDPIO::cout << "source norm 1 cpt = " << this->norm(b) << endl;

    std::vector<complex<S> > Wb(this->K);
    for(int k=0;k<this->K; k++){
        Wb[k] = this->innerProduct(this->bq[k],b);

    }
    multi1d<bfm_fermion> AW(this->K);
    for(int k=0;k<this->K; k++){
        this->init_fermion(AW[k]);
        this->herm_mult(this->bq[k],AW[k]);
    }

    Matrix<complex<S> > WAW(this->K);
    for(int j=0;j<this->K; j++){
        for(int k=0;k<this->K; k++){
            WAW(this->innerProduct(this->bq[j],AW[k]),j,k);
        }}

    WAW.Out();

    std::vector<complex<S> > WAWWb(this->K);
    print(Wb);
    CG_Matrix(WAW,Wb,WAWWb);
    print(WAWWb);
    QDPIO::cout << "norm x " << this->norm(x) << endl;
    this->axpby(x,0.0,x,0.0,x);
    bfm_fermion tmpt;
    this->init_fermion(tmpt);
    for(int k=0;k<this->K; k++){
        this->axpy(tmpt,0,this->bq[k],this->bq[k]);
        this->scale(tmpt,WAWWb[k]);
        this->axpy(x,1.0,x,tmpt);
        QDPIO::cout << "norm x " << this->norm(x) << endl;
    }
    this->free_fermion(tmpt);
    QDPIO::cout << "norm x " << this->norm(x) << endl;
    for(int k=0;k<this->K; k++){
        this->free_fermion(AW[k]);
    }

}

//copy from Qi's EigCG code
//QZ is saved in column major format. 
//perform V = V*QZ;
template<class Float> void eigcg_vec_mult(Float** V, const int m, double *QZ, const int n, const int f_size_cb, const int nthread, const int me)
{
    //implementation 4
    //reorder QZ to 4x4 blocks
    Float **Vptr = new Float*[m];
    assert(f_size_cb%nthread==0);
    int each = f_size_cb/nthread;
    for(int i=0;i<m;i++)Vptr[i] = V[i] + me*each;

    double *aux = new double[m*n];
    double *pt=aux;
    const int BS = 4;
    const int m_cblocks = m / BS + (m%BS ? 1 : 0);
    const int n_cblocks = n / BS + (n%BS ? 1 : 0);
    for(int c=0;c<n_cblocks;++c)// row direction
	{
            const int i=c*BS;
            for(int r=0;r<m_cblocks;++r) //column direction
		{
                    const int j=r*BS;
                    for(int ci=i;(ci<i+BS) && (ci<n); ci++)
                        for(int rj=j;(rj<j+BS)&&(rj<m);rj++)
                            *(pt++)=QZ[ci*m+rj];
		}
	}
    const int BBSS=24;
    int xlen=BBSS*n;
    double *x = new double[BBSS*n];
    for(int block=0;block<each/BBSS;block++)
	{
            memset(x,0,xlen*sizeof(double));
            matrix_dgemm(BBSS,n,m,Vptr,aux,x);
	}
    delete [] aux;
    delete [] x;
    delete [] Vptr;
}

/// q -> q Q
template <class Float> void Krylov<Float>::times_real(multi1d<bfm_fermion> &q, Matrix<Float> &Q, int N)
{

    Float **V;
    V = new Float* [N];

    double *QZ;
    QZ = new double [N*N];
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            QZ[i*N + j] = Q(j, i);

    int f_size_cb = dop.cbLs*24*dop.node_cbvol;

    for(int cb = this->prec; cb < 2; cb++)
	{
            for(int i = 0; i < N; i++)
                V[i] = (Float*)(q[i][cb]);

#pragma omp parallel
            {
                int me = dop.thread_barrier();
                eigcg_vec_mult(V, N, QZ, N, f_size_cb, dop.nthread, me);
            }
	}
    delete [] V;
    delete [] QZ;
	
}

///Check
template <class S> void Krylov<S>::Check(){
    /*vector<multi1d<LatticeFermion> > goodvec(0);*/
    vector<S> goodval(this->get);
    EigenSort(this->evals,this->evecs,this->comp);
    bfm_fermion bv; 
    this->init_fermion(bv);  
    bfm_fermion bvv; 
    this->init_fermion(bvv);

    int NM = this->M;
    int Nget = this->get;
    S **V;
    V = new S* [NM];

    double *QZ;
    QZ = new double [NM*NM];
    for(int i = 0; i < NM; i++)
        for(int j = 0; j < NM; j++)
            QZ[i*NM+j] = this->evecs[i][j];

    int f_size_cb = 24*dop.cbLs*dop.node_cbvol;

    for(int cb = this->prec; cb < 2; cb++)
	{
            for(int i = 0; i < NM; i++)
                V[i] = (S*)(this->bq[i][cb]);
#pragma omp parallel
            {
                int me = dop.thread_barrier();
                eigcg_vec_mult(V, NM, QZ, NM, f_size_cb, dop.nthread, me);
            }
	}
    delete [] V;
    delete [] QZ;

    for(int i=NM-1; i > NM-Nget-1 && i>=0;i--)
	{
            bfm_fermion bS;
            bS[0] = this->bq[i][0];
            bS[1] = this->bq[i][1];

            int ii = NM-1-i;
            double bdiff ;
            if(this->kr == lan)
		{
                    this->herm_mult(bS,bv);
                    goodval[ii] = this->innerProduct_real(bS,bv);	///Rayleigh quotient
                    double bS_norm = sqrt( this->axpy_norm(bS,0,bS,bS) );
                    goodval[ii] = goodval[ii]/bS_norm/bS_norm;
                    this->axpby(bS,0,bS,1.0/bS_norm,bS);
                    double bv_norm = sign( real(goodval[ii]) ) * sqrt( this->axpy_norm(bv,0,bv,bv) );
                    this->axpby(bv,0,bv,1.0/bv_norm,bv);	
                    bdiff = this->axpy_norm(bvv, -1.0 , bS, bv);
		}

            bdiff = sqrt(bdiff);
            string cv= " not converged";  
            if(bdiff < 10*this->conv ){ 
                cv = " converged";
            }else{
                QDPIO::cout << "WARNING LARGE RESIDUAL = " << bdiff 
                            << " on eval " << ii << endl;
            }
            QDPIO::cout <<"diff = " <<  bdiff << " < " << 10*this->conv 
			<< "  :  (" << real(goodval[ii]) << "," << imag(goodval[ii]) 
			<< ") eval # " << ii << cv << endl;
	}

    this->evals.resize(Nget);
    this->bl.resize(Nget);

    for(int i=0;i<Nget;i++)
	{
            this->evals[i] = goodval[i]; 
            this->bl[i] = goodval[i];
            QDPIO::cout<<i<<":"<<evals[i]<<endl;
	}

    for(int i = 0; i < Nget; i++)
        this->axpy(this->bq[i], 0.0, bv, this->bq[i+NM-Nget]);

    for(int i = 0; i < Nget/2; i++)
	{
            this->axpy(bv, 0.0, bv, this->bq[i]);
            this->axpy(this->bq[i], 0.0, bv, this->bq[Nget-1-i]);
            this->axpy(this->bq[Nget-1-i], 0.0, bv, bv);
	}

    this->free_fermion(bv);
    this->free_fermion(bvv);

    //release overhead memory, keep only the eigenvectors
    for(int i = Nget; i < NM; i++)
        free_fermion(bq[i]);
#pragma omp parallel 
    {
        dop.threadedFreeFermion(tmp);
    }
    free_fermion(tmp1);
    free_fermion(tmp2);
    free_fermion(invec);
    free_fermion(outvec);
    for(int i = 0; i < bf.size(); i++)
        free_fermion(bf[i]);

    QDPIO::cout << "total number of vector inner products " << this->innerprod << endl; 
    QDPIO::cout << "total number of matrix vector products " << this->mvprod << endl;
}
#endif

