#ifndef BFM_KRY5D_H
#define BFM_KRY5D_H

#include "Krylov.h"
#include <cstdlib>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <complex>

namespace BFM_Krylov{

template <class S> 
class Krylov_5d : public Krylov<S>
{
  public:
    Krylov_5d(bfm_evo<S> &dwf) : Krylov<S>(dwf){}
    Krylov_5d(bfm_evo<S> &dwf, cps::LancArg &lanc_arg) : Krylov<S>(dwf, lanc_arg) {}
    virtual ~Krylov_5d() {}
    void eigs_to_qdp(QDP::multi1d<QDP::multi1d<QDP::LatticeFermion> > &con_vecs, QDP::multi1d<QDP::Complex> &con_vals);

    //Copied from Daiqian. Change stored eigenvectors 'bq' to single precision
    void toSingle();
  protected:
    void init();

    virtual void dwf_multiply(bfm_fermion input, bfm_fermion &result);
    void dwf_multiply(bfm_fermion input, bfm_fermion &result, bool dag);
    virtual void herm_mult(bfm_fermion input, bfm_fermion &result);
};

template <class S> 
void Krylov_5d<S>::toSingle() // Change bq to single precision
{
	int words = 24 * this->dop.node_cbvol * this->dop.cbLs * (1 
#ifdef BFM_GPARITY
+ this->dop.gparity
#endif
);

	for(int i = 0; i < this->get; i++) {
		bfm_fermion bq_tmp;
		for(int cb = this->prec; cb < 2; cb++) {
			bq_tmp[cb] = (float *)bfm_alloc(words * sizeof(float),mem_slow);
			if(bq_tmp[cb] == 0){ printf("bfmbase::allocFermion\n"); fflush(stdout); exit(-1);}
			for(int j = 0; j < words; j++) {
				((float*)bq_tmp[cb])[j] = ((S*)(this->bq[i][cb]))[j];
				//QDPIO::cout<<"i="<<i<<" "<<"j="<<j<<" "<<((float*)bq_tmp[cb])[j]<<" "<<((S*)this->bq[i][cb])[j]<<endl; // test passed
			}
		}
		this->free_fermion(this->bq[i]);
		for(int cb = this->prec; cb < 2; cb++) {
			this->bq[i][cb] = bq_tmp[cb];
		}
	}
}


template <class S> 
void Krylov_5d<S>::init()
{
  if(this->initialized)
    return;

  QDPIO::cout << "initializing" << std::endl;
  this->Krylov_init();
  if(this->D == G5D && this->dop.solver == DWF && this->prec == 1)
  {
    QDPIO::cerr << "5d preconditioned gamma5 D operator is not Hermitian. " << std::endl;
    exit(0);
  }

  this->init_fermion(this->tmp1); 
  this->init_fermion(this->tmp2);

#pragma omp parallel 
  {
    this->tmp = this->dop.threadedAllocFermion();
    for(int cb = this->prec; cb < 2; cb++)
      {			
	this->invec[cb] = this->dop.threadedAllocFermion();
	this->outvec[cb] = this->dop.threadedAllocFermion();
      }
  }
  

  this->bq.resize(this->M); 
  for(int i = 0; i < this->M; ++i)
    this->init_fermion(this->bq[i]);

  int sizebf;
  if(this->kr == lan)
    sizebf = 1;
  else{
    QDPIO::cerr << "bf vector unknown size for thie krylovtype\n";
    exit(-1);
  }

  this->bf.resize(sizebf); 

#ifdef BFM_GPARITY
  if(this->dop.gparity) QDPIO::cout << "Krylov init for 2f G-parity\n";
  else 
#endif
  QDPIO::cout << "Krylov init for standard fermions\n";

  for(int i = 0; i < sizebf; ++i){
    this->init_fermion(this->bf[i]);
    this->zero_fermion(this->bf[i]);
  }
  Fermion of = zero;
  SpinVector tspin = zero;
  
  //CK: seems to generate a source vector with 1.0 on every spin-color index
  QDP::Complex cno = cmplx(Real(1), Real(0));
  for(int spn = 0; spn < Nd; spn++) 
    pokeSpin(tspin, cno, spn);
  //  for(int col = 0; col < Nd; col++)  //shouldn't this be 3, not 4??
  //  pokeColor(of, tspin, col);
  for(int col = 0; col < 3; col++)
    pokeColor(of, tspin, col);

#ifdef BFM_GPARITY
  if(!this->dop.gparity)
#endif
    {
    multi1d<LatticeFermion> st(this->dop.Ls);
    //CK: place the 4d source vector on both walls of the 5th dimension and set other s-slices to 0
    st[0] = of; 
    st[this->dop.Ls - 1] = of; 
    for(int k = 1; k < this->dop.Ls - 1; ++k)
      st[k] = zero;
    this->qdp_to_bfm(st, this->bq[0]);
  }
#ifdef BFM_GPARITY
else{
    //CK: for G-parity, on each s-slice, 2 4d fields are expected: s=0 |f0 f1| s=1 |f0 f1| ....
    //As the source is 1.0 on every spin-color index, it makes sense to make it 1.0 on both flavour indices too
    printf("Making multi1d<LatticeFermion> of size %d\n", 2*this->dop.Ls);
    multi1d<LatticeFermion> st(2*this->dop.Ls);
    //CK: place the 4d source vector on both walls of the 5th dimension and set other s-slices to 0
    st[0] = of; st[1] = of;
    st[2*this->dop.Ls - 1] = of;  st[2*this->dop.Ls - 2] = of; 

    for(int k = 2; k < 2*this->dop.Ls - 2; ++k)
      st[k] = zero;
    this->qdp_to_bfm(st, this->bq[0]);
  }
#endif
  double hnorm = this->norm(this->bq[0]);
  printf("Starting vector norm %g\n",hnorm);
  
  this->axpby(this->bq[0], 0, this->bq[0], 1.0/hnorm, this->bq[0]);
  this->initialized = true;
}

template <class S> 
void Krylov_5d<S>::eigs_to_qdp(QDP::multi1d<QDP::multi1d<QDP::LatticeFermion> > &con_vecs, QDP::multi1d<QDP::Complex> &con_vals)
{
	int Nev = this->evals.size();
	con_vecs.resize(Nev); for(int i = 0; i<Nev; i++){con_vecs[i].resize(this->dop.Ls);}
	con_vals.resize(Nev);
	for(int i = 0; i<Nev; i++){
		con_vals[i] = this->evals[i];
		this->bfm_to_qdp(this->bq[i], con_vecs[i]);
	}
}

template <class S> 
void Krylov_5d<S>::dwf_multiply(bfm_fermion input, bfm_fermion &result)
{
  bool dag = false;
  if(this->D == DDAG)
    dag = true;
    
#pragma omp parallel 
  {
		if(this->prec == 1)
			this->dop.CompactMprec(input[1], result[1], this->invec[1], this->outvec[1], this->tmp, dag, 0) ;
		else 
			this->dop.CompactMunprec(input, result, this->invec, this->outvec, this->tmp, dag) ;
  }
}

template <class S> 
void Krylov_5d<S>::dwf_multiply(bfm_fermion input, bfm_fermion &result, bool dag)
{
#pragma omp parallel 
  {
		if(this->prec == 1)
			this->dop.CompactMprec(input[1], result[1], this->invec[1], this->outvec[1], this->tmp, dag, 0) ;
		else 
			this->dop.CompactMunprec(input, result, this->invec, this->outvec, this->tmp, dag) ;
  }
}

template <class S> 
void Krylov_5d<S>::herm_mult(bfm_fermion input, bfm_fermion &result)
{
	this->dslash_time -= dclock();
  if(this->D == G5D)
  {
    dwf_multiply(input, result, false);  
    this->mvprod++;  //D
    this->axpy(this->tmp1 , 0.0, result, result);
    this->G5R(result, this->tmp1);
  } 
  else if(this->D == DDAGD)
  {
    dwf_multiply(input, this->tmp1, false);
    dwf_multiply(this->tmp1, result, true);    
    this->mvprod += 2;
  }
  else{
    QDPIO::cerr << "Krylov_5d<S>::herm_mult : Don't know what to do with that operator" << std::endl;
    exit(-1); //Added by CK: why would we continue??
  }

  this->dslash_time += dclock();
  

}

template <class T> 
class Lanczos_5d : public Krylov_5d<T>
{
  public:
    Lanczos_5d(bfm_evo<T> &dwf) : Krylov_5d<T>(dwf){this->kr = lan;}
    Lanczos_5d(bfm_evo<T> &dwf, cps::LancArg &lanc_arg) : Krylov_5d<T>(dwf, lanc_arg) {this->kr = lan;} //Jasper
    virtual ~Lanczos_5d(){}

    void do_init(){ this->init(); } //CK: testing workaround for fact that init is protected

    //Jasper
    //void Initial_State(multi1d< multi1d<LatticeFermion> > Eigs, multi1d<Complex> val)
    //{
    //  this->init();
    //  if(this->cont)
    //  {
    //    QDPIO::cout << "Lanczos5d: continue " << Eigs.size() << endl;
    //    for(int i=0;i<Eigs.size();i++){
    //      this->qdp_to_bfm(Eigs[i],this->bq[i]); 
    //      this->H( this->Chebyshev( toDouble(real(val[i]) )) , i, i);	
    //    }
    //    std::vector<double> tevals(Eigs.size());
    //    std::vector<std::vector<double> > tevecs(Eigs.size());
    //    this->TestConv(Eigs.size(), tevals, tevecs);

    //    multi1d<LatticeFermion> tmp(this->Ls);
    //    for(int s = 0; s< this->Ls; s++){gaussian(tmp[s]);}
    //    this->qdp_to_bfm(tmp, this->bf[0]);
    //    this->Gram(this->bf[0],this->bq,Eigs.size());
    //  } 
    //  else if(this->refine)
    //  {
    //    QDPIO::cout << "Lanczos5d: refining " << endl;
    //    multi1d<LatticeFermion> gg(this->Ls);
    //    for(int s = 0; s< this->Ls; s++) gg[s] = Eigs[0][s];

    //    for(int i=1;i<Eigs.size();i++){
    //      for(int s = 0; s< this->Ls; s++) gg[s] = gg[s] + Eigs[i][s]; 
    //    }

    //    Real nm = 0.0; for(int s = 0; s< this->Ls; s++) nm += norm2(gg[s]);
    //    nm = sqrt(nm); nm = 1.0/nm;
    //    for(int s = 0; s< this->Ls; s++) gg[s] = nm*gg[s];
    //    this->qdp_to_bfm(gg,this->bf[0]); 
    //  }
    //}
};

}

#endif

