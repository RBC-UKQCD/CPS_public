#ifndef KRY5D_H
#define KRY5D_H

#include "Krylov.h"
#include <cstdlib>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <complex>

template <class S> 
class Krylov_5d : public Krylov<S>
{
  public:
    Krylov_5d(bfm_evo<S> &dwf) : Krylov<S>(dwf){}
    Krylov_5d(bfm_evo<S> &dwf, cps::LancArg &lanc_arg) : Krylov<S>(dwf, lanc_arg) {}
    virtual ~Krylov_5d() {}
    void eigs_to_qdp(multi1d<multi1d<LatticeFermion> > &con_vecs, multi1d<Complex> &con_vals);

  protected:
    void init();
    virtual void dwf_multiply(bfm_fermion input, bfm_fermion &result);
    void dwf_multiply(bfm_fermion input, bfm_fermion &result, bool dag);
    virtual void herm_mult(bfm_fermion input, bfm_fermion &result);
};

template <class S> 
void Krylov_5d<S>::init()
{
  if(this->initialized)
    return;

  QDPIO::cout << "initializing" << endl;

  this->Krylov_init();
  if(this->D == G5D && this->dop.solver == DWF && this->prec == 1)
  {
    QDPIO::cerr << "5d preconditioned gamma5 D operator is not Hermitian. " << endl;
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
  this->bf.resize(sizebf); 
  for(int i = 0; i < sizebf; ++i)
  {
    this->init_fermion(this->bf[i]);
    this->zero_fermion(this->bf[i]);
  }
  multi1d<LatticeFermion> st(this->dop.Ls);
  Fermion of = zero;
  SpinVector tspin = zero;
  Complex cno = cmplx(Real(1), Real(0));
  for(int spn = 0; spn < Nd; spn++) 
    pokeSpin(tspin, cno, spn);
  for(int col = 0; col < Nd; col++) 
    pokeColor(of, tspin, col);

  st[0] = of; 
  st[this->dop.Ls - 1] = of; 
  for(int k = 1; k < this->dop.Ls - 1; ++k)
    st[k] = zero;
  this->qdp_to_bfm(st, this->bq[0]);
  double hnorm = this->norm(this->bq[0]);
  this->axpby(this->bq[0], 0, this->bq[0], 1.0/hnorm, this->bq[0]);
  this->initialized = true;
}

template <class S> 
void Krylov_5d<S>::eigs_to_qdp(multi1d<multi1d<LatticeFermion> > &con_vecs, multi1d<Complex> &con_vals)
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
    G5R(result, this->tmp1);
  } 
  else if(this->D == DDAGD)
  {
    dwf_multiply(input, this->tmp1, false);
    dwf_multiply(this->tmp1, result, true);    
    this->mvprod += 2;
  }
  else
    QDPIO::cerr << "Don't know what to do with that operator" << endl;
	this->dslash_time += dclock();
}

template <class T> 
class Lanczos_5d : public Krylov_5d<T>
{
  public:
    Lanczos_5d(bfm_evo<T> &dwf) : Krylov_5d<T>(dwf){this->kr = lan;}
    Lanczos_5d(bfm_evo<T> &dwf, cps::LancArg &lanc_arg) : Krylov_5d<T>(dwf, lanc_arg) {this->kr = lan;} //Jasper
    virtual ~Lanczos_5d(){}
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
#endif

