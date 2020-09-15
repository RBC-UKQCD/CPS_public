#ifndef BFM_KRY4D_H
#define BFM_KRY4D_H

#include "Krylov.h"
#include <cstdlib>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <complex>
#include <bfm_qdp_g5d.h>

namespace BFM_Krylov{


//#include "dwf_4d.h"
/**
A 4d eigensolver.
**/

template <class S> class Krylov_4d : public Krylov<S>{

  public:

  bfmDwf4d dop4d;
  bfm_qdp<S>  dop5d;

  int total_inner;
  int fivecount;

  LatticeFermion tmp_lf;
  SpinMatrix ident, Gamma5, Pp, Pm;
  
  Krylov_4d(bfm_qdp<S> &dwf) : Krylov<S>(dwf){this->prec = 0; this->prectype = 1;}
  ~Krylov_4d(){dop4d.end();dop5d.end();}
  void eigs_to_qdp( multi1d<LatticeFermion> &con_vecs, multi1d<Complex> &con_vals );

  protected:
  void init();
  void four_to_five(Fermion_t input[2], Fermion_t result5d[2]);
  void five_to_four(Fermion_t input5d[2], Fermion_t result[2], int sgn);
  void dwf_multiply(bfm_fermion input, bfm_fermion &result, bool dagger, int sgn);
  virtual void dwf_multiply(bfm_fermion input, bfm_fermion &result);
  virtual void herm_mult(Fermion_t input[2], bfm_fermion &result);

};


 template <class S> void Krylov_4d<S>::eigs_to_qdp( multi1d<LatticeFermion> &con_vecs, multi1d<Complex> &con_vals ){

	int Nev = this->evals.size();
	con_vecs.resize(Nev); 
	con_vals.resize(Nev);
	for(int i = 0; i<Nev; i++){
		con_vals[i] = this->evals[i];
		this->bfm_to_qdp(this->bq[i], con_vecs[i]);
	}
 }

  template <class S> void Krylov_4d<S>::init(){

    if(!this->initialized){

    this->Krylov_init();
    QDPIO::cout << "Krylov4d: Initialising effective 5d DWF." << endl;
    //this->dwfa.solver = DWFrb4d;
    this->dwfa.precon_5d = this->prectype;
    this->dwfa.Ls = this->Ls;
    this->dwfa.residual = this->conv;
    dop5d.init(this->dwfa);
#pragma omp parallel 
  {
    omp_set_num_threads(dop5d.nthread);
#pragma omp for 
    for(int i=0;i<dop5d.nthread;i++) {
    this->tmp = dop5d.threadedAllocFermion();
}}    
    QDPIO::cout << "Krylov4d: Initialising 4d DWF for linalg" << endl;
    //this->dwfa.solver = HtCayleyTanh;
    this->dwfa.precon_5d = 0;
    this->dwfa.Ls = 1;
    this->dop.init(this->dwfa);
    
    total_inner = 0;
    ///4d initialize
    QDPIO::cout << "Krylov4d: Initialising effective 4d DWF object." << endl;
    //dop4d.init(toDouble(this->mq),toDouble(this->M5),(int)this->Ls,this->u, this->dwfa.residual);
    dop4d.init(this->dwfa.solver,toDouble(this->mq),toDouble(this->M5),(int)this->Ls,this->u, this->dwfa.residual, this->dwfa.max_iter);

    this->bq.resize(this->M); for(int m=0;m<this->M;m++){this->init_fermion(this->bq[m]);}
    this->init_fermion(this->tmp1); 
    this->init_fermion(this->tmp2);
    
    //if(!this->cont && !this->refine){

      QDPIO::cout << "Krylov4d: Starting from scratch." << endl;
      LatticeFermion st;
      //gaussian(st);
      Fermion of = zero;
      SpinVector tspin = zero;
      Complex cno = cmplx(Real(1),Real(0));
      for(int spn=0;spn<Nd;spn++) pokeSpin(tspin, cno, spn);
      for(int col=0;col<Nd;col++) pokeColor(of, tspin, col);
      st = of; 

      this->qdp_to_bfm(st,this->bq[0]);  
      double norm = sqrt( this->axpy_norm(this->bq[0],0.0,this->bq[0],this->bq[0]) );
      this->axpby(this->bq[0], 1.0/norm, this->bq[0], 0.0, this->bq[0]);
    //}
      	int sizebf;
	if(this->kr == lan){sizebf = 1;}
	//else if(this->kr == arn){sizebf = this->M;}
	this->bf.resize(sizebf); 
	for(int i=0; i< sizebf; i++){
        	this->init_fermion(this->bf[i]);
		this->zero_fermion(this->bf[i]);
	}
	this->initialized=true;	
        fivecount = 0;
    }
  }

  template <class S> void Krylov_4d<S>::four_to_five(Fermion_t input[2], Fermion_t result5d[2]){

#pragma omp parallel 
  {
    omp_set_num_threads(dop5d.nthread);
#pragma omp for 
    for(int i=0;i<dop5d.nthread;i++) {
      for(int cb=0;cb<2;cb++){
	dop5d.axpby(result5d[cb],result5d[cb],result5d[cb],0.0,0.0);
      }
    }
}
    int x[4];
    int Nspinco=12;
    int site, cb4, cb5;
    int sp, b5idx, bidx;
    S *b5;
    S *b;

   // QDPIO::cout << "four_to_five" << endl;


    for ( x[3]=0; x[3]<this->dop.node_latt[3];x[3]++ ) { 
    for ( x[2]=0; x[2]<this->dop.node_latt[2];x[2]++ ) {   
    for ( x[1]=0; x[1]<this->dop.node_latt[1];x[1]++ ) { 
    for ( x[0]=0; x[0]<this->dop.node_latt[0];x[0]++ ) { 
	    
      sp = 0;
      site = x[0]+x[1]+x[2]+x[3];
      cb4 = ((site)&0x1);
      cb5 = ((site+sp)&0x1);

      for ( int co=Nspinco/2;co<Nspinco;co++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 
	  b5 = (S *)result5d[cb5];
	  b = (S *)input[cb4];
	  b5idx = dop5d.bagel_idx5d(x,sp,reim,co,Nspinco,1);
	  bidx = this->dop.bagel_idx(x,reim,co,Nspinco,1);
	  b5[b5idx] = b[bidx]; 
      }}

      sp = this->Ls - 1;
      cb4 = ((site)&0x1);
      cb5 = ((site+sp)&0x1);
      
      for ( int co=0;co<Nspinco/2;co++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 
        b5 = (S *)result5d[cb5];
	b = (S *)input[cb4];
        b5idx = dop5d.bagel_idx5d(x,sp,reim,co,Nspinco,1);
        bidx = this->dop.bagel_idx(x,reim,co,Nspinco,1);
	b5[b5idx] = b[bidx]; 
      }}

    }}}}

  }

  template <class S> void Krylov_4d<S>::five_to_four(Fermion_t input5d[2], Fermion_t result[2], int sgn){

    int x[4];
    int Nspinco=12;
    int site, cb4, cb5;
    int sp, b5idx, bidx;
    S * b5;
    S * b;

    this->axpby(this->tmp1,0.0,this->tmp1,0.0,this->tmp1);
    this->axpby(this->tmp2,0.0,this->tmp2,0.0,this->tmp2);

    for ( x[3]=0; x[3]<this->dop.node_latt[3];x[3]++ ) { 
    for ( x[2]=0; x[2]<this->dop.node_latt[2];x[2]++ ) { 
    for ( x[1]=0; x[1]<this->dop.node_latt[1];x[1]++ ) { 
    for ( x[0]=0; x[0]<this->dop.node_latt[0];x[0]++ ) { 

      sp = 0;
      site = x[0]+x[1]+x[2]+x[3];
      cb4 = ((site)&0x1);
      cb5 = ((site+sp)&0x1);

      for ( int co=Nspinco/2;co<Nspinco;co++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 
        S * b5 = (S *)input5d[cb5];
	S * b = (S *)result[cb4];
        int b5idx = dop5d.bagel_idx5d(x,sp,reim,co,Nspinco,1);
        int bidx = this->dop.bagel_idx(x,reim,co,Nspinco,1);
	b[bidx] = sgn*b5[b5idx]; 
      }}

      sp = this->Ls - 1;
      cb4 = ((site)&0x1);
      cb5 = ((site+sp)&0x1);

      for ( int co=0;co<Nspinco/2;co++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 
        S * b5 = (S *)input5d[cb5];
	S * b = (S *)result[cb4];
        int b5idx = dop5d.bagel_idx5d(x,sp,reim,co,Nspinco,1);
        int bidx = this->dop.bagel_idx(x,reim,co,Nspinco,1);
	b[bidx] = b5[b5idx]; 
      }}

    }}}}

  }


  template <class S> void Krylov_4d<S>::dwf_multiply(bfm_fermion input, bfm_fermion &result, bool dagger, int sgn){

/*    Fermion_t five[2], five_res[2];
#pragma omp parallel 
  {
    omp_set_num_threads(dop5d.nthread);
#pragma omp for 
    for(int i=0;i<dop5d.nthread;i++) {
    for(int cb=0;cb<2;cb++){
      five[cb]     = dop5d.threadedAllocFermion(); 
      five_res[cb] = dop5d.threadedAllocFermion(); 
    }
}
}
    four_to_five(input, five);

    dop4d.Munprec(five,five_res,dagger); 
    this->mvprod++;

    five_to_four(five_res, result, sgn);
#pragma omp parallel 
  {
    omp_set_num_threads(dop5d.nthread);
#pragma omp for 
    for(int i=0;i<dop5d.nthread;i++) {
    for(int cb=0;cb<2;cb++){
      dop5d.threadedFreeFermion(five_res[cb]); 
      dop5d.threadedFreeFermion(five[cb]);
    }
}}*/


LatticeFermion in, res;
this->bfm_to_qdp(input,in);
dop4d.Munprec(in,res,dagger); 
if(sgn==-1){
	Gamma G5(15);
	SpinMatrix ident = Real(1.0);
	SpinMatrix Gamma5 = (G5 * ident);
	res = Gamma5 * res;
}
this->qdp_to_bfm(res,result);



  }

  template <class S> void Krylov_4d<S>::dwf_multiply(bfm_fermion input, bfm_fermion &result){

      bool dagger = false;
      if(this->D == DDAG){ dagger = true;}
      dwf_multiply(input, result, dagger, 1);  
      
    }


  template <class S> void Krylov_4d<S>::herm_mult(Fermion_t input[2], bfm_fermion &result){
    if(this->D == G5D){
      dwf_multiply(input, result, false, -1);
    } 
    else if(this->D == DDAGD){
      bfm_fermion Dq; this->init_fermion(Dq);
      dwf_multiply(input, Dq, false, 1); 
      dwf_multiply(Dq, result, true, 1);
      this->free_fermion(Dq);
    }
    else{
      QDPIO::cerr << "Don't know what to do with that operator" << endl; exit(1);
    }
  }


template <class T> class Lanczos_4d : public Krylov_4d<T>{
public:

Lanczos_4d(bfm_qdp<T> &dwf) : Krylov_4d<T>(dwf){this->kr = lan;}
void Initial_State(multi1d<LatticeFermion> Eigs, multi1d<Complex> val){
		
		this->init();
		if(this->cont){
  			QDPIO::cout << "Lanczos4d: continue " << this->bq.size() << endl;
			for(int i=0;i<Eigs.size();i++){
				this->qdp_to_bfm(Eigs[i],this->bq[i]); 
				this->H( this->Chebyshev( toDouble(real(val[i]) )) , i, i);	
			}
			std::vector<double> tevals(Eigs.size());
        		std::vector<std::vector<double> > tevecs(Eigs.size());
			this->TestConv(Eigs.size(), tevals, tevecs);
    
			LatticeFermion tmp; 
			gaussian(tmp);
			this->qdp_to_bfm(tmp,this->bf[0]); 
			this->Gram(this->bf[0],this->bq,Eigs.size());
		} 		
		else if(this->refine){
  			QDPIO::cout << "Lanczos4d: refining " << endl;
			LatticeFermion gg;
			gg = Eigs[0];

			for(int i=1;i<Eigs.size();i++){
				gg = gg + Eigs[i]; 
			}

			Real norm = 1.0/sqrt( norm2(gg) );
			gg = norm*gg;
			this->qdp_to_bfm(gg,this->bq[0]); 
		}
}

};





/**
An Arnoldi eigensolver.
**/
/*template <class T> class Arnoldi_4d : public Krylov_4d<complex<T>, T, complex<double> >{
public:



Arnoldi_4d(bfm_qdp<T> dwf, EigenDescr Ed, eigentype c, multi1d<LatticeColorMatrix> UU, multi1d<LatticeFermion> Eigs,
multi1d<Complex> val) : Krylov_4d<complex<T>, T, complex<double> >(dwf, Ed, c, UU, Eigs, val){

  	this->kr = arn;
  	this->bf.resize(this->M); for(int i=0;i<this->M;i++)this->init_fermion(this->bf[i]);

}


};*/

}

#endif
