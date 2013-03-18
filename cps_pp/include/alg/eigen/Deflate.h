#ifndef DEFLATE_H
#define DEFLATE_H

/**********************************************
*
*  TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*  This is a massive lag point. Figure out how to make it
*  exponentially faster or don't do it.
*
************************************************/


#include "Francis.h"
#include <cstdlib>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <complex>
#include <time.h>

/** Return smallest number number such that 1 + epsilon > 1.
 stolen from Tony Kennedy*/
static double findUlp() {
	double epsilon = 1., ope;
	do { epsilon *= 1./2.; ope = 1. + epsilon; } while (ope > 1.);
	return 2. * epsilon;
}
/** Unit of least precision (ULP), smallest number such that 1 + ulp > 1*/
const double ulp = findUlp();
/** Tolerance for floating point computations*/
const double eps = ulp * (double) 3.;

/**Deflation functions live here. Locking and (one day) purging**/
template <class T> class Deflate{
public:

/**
	There is some matrix Q such that for any vector y
	Q.e_1 = y and Q is unitary.
**/
static T orthQ(Matrix<T> &Q, vector<T> y){

int N = y.size();	//Matrix Size
Q.Fill(0.0);
T tau;

for(int i=0;i<N;i++){Q(y[i],i,0);}
T sig = conj(y[0])*y[0];
T tau0 = abs(sqrt(sig));

for(int j=1;j<N;j++){
sig += conj(y[j])*y[j]; 
tau = abs(sqrt(sig) ); 	


if(abs(tau0) > 0.0){

	T gam = conj( (y[j]/tau)/tau0 );
	for(int k=0;k<=j-1;k++){  Q(-gam*y[k],k,j); }
	Q(tau0/tau,j,j);
}
else{
	Q(1.0,j-1,j);
}
tau0 = tau;
}

return tau;
}


/**
	There is some matrix Q such that for any vector y
	Q.e_1 = y and Q is unitary.
**/
static double orthQM(Matrix<T> &Q, Matrix<T> &H, vector<T> y){

int N = y.size();	//Matrix Size
Q.Fill(0.0);
T tau;

int ii=0; while(real(y[ii]) == 0 && ii < N){y[ii] = eps/(double)N; ii++;}
  
T sig = conj(y[0])*y[0];
T tau0 = abs(y[0]);

for(int j=1;j<N;j++){ 
  sig += conj(y[j])*y[j];
  tau = sqrt(sig) ; 

if(abs(tau0) > 0.0){	
  T gam = -1.0*conj( ( y[j]/tau)/tau0 );
  T rho =  tau0/tau;
  for(int k=0;k<=j-1;k++){  Q(gam*y[k],k,j); } 
  Q(rho,j,j);
	if(j < (N-1) && abs(tau) < 0.05){
	  T taup = sqrt(sig + conj(y[j+1])*y[j+1]);
	  T alpha = 0.0; vector<T> tmp(j+1); for(int k=0;k<j+1;k++){tmp[k] = Q(k,j);}
			      Matrix<T> Hj(j+1); Hj = H.GetSubMtx(0, j+1, 0, j+1);
			      vector<T> tmp2(j+1); tmp2 = Hj*tmp; for(int k=0;k<j+1;k++){alpha+=conj(y[k])*tmp2[k];} //TODO
	  T delta = y[j+1]*H(j,j+1)*rho;
	  T sig0 = -1.0;
	      if(real( sign(alpha)*sign(delta) ) < 0.0){sig0 = 1.0;}
	      if(abs(delta + alpha) > eps*abs(taup) ){
		  if(abs(rho) < sqrt(eps)/100.0){
		    T qjj = Q(j,j); Q(qjj*sig0*(eps*taup + abs(alpha))/abs(delta),j,j);
		  }
		  else{
			if(abs(alpha) > abs(delta)){
			   T phi = sig0*(eps*taup + abs(delta))/abs(alpha);
			   for(int k=0;k<j+1;k++){y[k] = y[k]*phi;}
			   sig = sig*phi*phi; tau = tau * abs(phi);
			}
			else{
			   T psi = sig0*(eps*taup + abs(alpha))/abs(delta);
			   for(int k=j+1;k<N;k++){y[k] = y[k]*psi;}
			}
		  }  
	      }
	}
}
else{
	Q(1.0,j-1,j);
}
tau0 = tau;
}
T nm = real( norm(y) );
for(int i=0;i<N;i++){Q(y[i]/nm,i,0);}
return 0;
}

/** Left shift on a vector**/
static void SL(vector<T> &y){

T front = y[0];
for(int i=0;i<y.size()-1;i++){y[i] = y[i+1];}
y[y.size()-1] = front; 

}
/** Left shift on a matrix**/
static void SL(Matrix<T> &H){

int N = H.dim;

for(int r=0;r<N;r++){
	T front = H(r,0);
	for(int i=0;i<N-1;i++){H(  H(r,i+1)  ,r,i);}
	H(front , r , N-1); 
}

}

/**
	There is some matrix Q such that for any vector y
	Q.e_k = y and Q is unitary.
**/
static T orthU(Matrix<T> &Q, vector<T> y){

T tau = orthQ(Q,y);
SL(Q);
return tau;

}


/**
	Wind up with a matrix with the first con rows untouched

say con = 2
	Q is such that Qdag H Q has {x, x, val, 0, 0, 0, 0, ...} as 1st colum
	and the matrix is upper hessenberg
	and with f and Q appropriately modidied with Q is the arnoldi factorization

**/

static void Lock(Matrix<T> &H, 	///Hess mtx	
	Matrix<T> &Q, 	///Lock Transform
	T val, 		///value to be locked
	int con, 	///number already locked
	double small,
	int dfg,
	bool herm){	

	//ForceTridiagonal(H);

	int M = H.dim;
	vector<T> vec(M-con);
	Matrix<T> AH(M-con);
	AH = H.GetSubMtx(con, M, con, M);
	Matrix<T> QQ(M-con);
	Q.Unity(); QQ.Unity();
	
	vector<T> evals(M-con);
	vector<vector<T> > evecs(M-con);


	if(herm){Wilkinson<T>(AH, evals, evecs, small);}
	else{QReigensystem<T>(AH, evals, evecs, small);}

	int k=0;
	double cold = abs( val - evals[k]); 
	for(int i=1;i<M-con;i++){
		double cnew = abs( val - evals[i]);
		if( cnew < cold ){k = i; cold = cnew;}
		
	}
	vec = evecs[k];

	complex<double> tau;
	orthQ(QQ,vec);
	//orthQM(QQ,AH,vec);

	AH = QQ.Hermitian()*AH;
	AH = AH*QQ;


	for(int i=con;i<M;i++){
	for(int j=con;j<M;j++){
		Q( QQ(i-con, j-con) , i,j);
		H( AH(i-con, j-con) , i,j);
	}
	}




	for(int j = M-1; j>con+2; j--){

	Matrix<T> U(j-1-con);
	vector<T> z(j-1-con); T nm = norm(z); for(int k = con+0;k<j-1;k++){z[k-con] = conj( H(j,k+1) );}
	normalize(z);

	double tmp = 0;
	for(int i=0;i<z.size()-1;i++){tmp = tmp + abs(z[i]);}
	if(tmp < small/( (double)z.size()-1.0) ){ continue;}	

	tau = orthU(U,z);

	vector<vector<T> > Hb(j-1-con);	
	for(int k=con+1;k<j;k++){Hb[k-1-con].resize(M);}
	
	for(int a = 0;a<M;a++){
	for(int b = 0;b<j-1-con;b++){T sum = 0;
	for(int c = 0;c<j-1-con;c++){sum += H.A[a][con+1+c]*U.A[c][b];}//sum += H(a,con+1+c)*U(c,b);}
	Hb[b][a] = sum;}}
	
	for(int k=con+1;k<j;k++){
	for(int l=0;l<M;l++){H.A[l][k] = Hb[k-1-con][l];}}//H(Hb[k-1-con][l] , l,k);}}

	vector<vector<T> > Qb(j-1-con);	
	for(int k=con+1;k<j;k++){Qb[k-1-con].resize(M);} 
	
	for(int a = 0;a<M;a++){
	for(int b = 0;b<j-1-con;b++){T sum = 0;
	for(int c = 0;c<j-1-con;c++){sum += Q.A[a][con+1+c]*U.A[c][b];}//sum += Q(a,con+1+c)*U(c,b);}
	Qb[b][a] = sum;}}
	
	for(int k=con+1;k<j;k++){
	for(int l=0;l<M;l++){Q.A[l][k] = Qb[k-1-con][l];}}//Q(Qb[k-1-con][l] , l,k);}}

	vector<vector<T> > Hc(M);	
	for(int k=0;k<M;k++){Hc[k].resize(j-1-con);}
	
	
	for(int a = 0;a<j-1-con;a++){
	for(int b = 0;b<M;b++){T sum = 0;
	for(int c = 0;c<j-1-con;c++){sum += conj( U.A[c][a] )*H.A[con+1+c][b];}//sum += conj( U(c,a) )*H(con+1+c,b);}
	Hc[b][a] = sum;}}
	
	for(int k=0;k<M;k++){
	for(int l=con+1;l<j;l++){H.A[l][k] = Hc[k][l-1-con];}}//H(Hc[k][l-1-con] , l,k);}}


	}


	
}


/**
	Wind up with a matrix with the first con rows untouched

say con = 2
	Q is such that Qdag H Q has {x, x, val, 0, 0, 0, 0, ...} as 1st colum
	and the matrix is upper hessenberg
	and with f and Q appropriately modidied with Q is the arnoldi factorization

**/
/*
static void Lock(Matrix<T> &H, 	///Hess mtx	
	Matrix<T> &Q, 	///Lock Transform
	T val, 		///value to be locked
	int con, 	///number already locked
	double small,
	int dfg,
	bool herm){	

	int M = H.dim;
	vector<T> vec(M-con);
	Matrix<T> AH(M-con);
	

	for(int i=0;i<M-con;i++){
	for(int j=0;j<M-con;j++){
		AH( H(i+con, j+con), i, j);
	}
	}

	vector<T> evals(M-con);
	vector<vector<T> > evecs(M-con);


	if(herm){Wilkinson<T>(AH, evals, evecs, small);}
	else{Eigensystem<T>(AH, evals, evecs, small);}

	int k=0;
	double cold = abs( val - evals[k]); 
	for(int i=1;i<M-con;i++){
		double cnew = abs( val - evals[i]);
		if( cnew < cold ){k = i; cold = cnew;}
		
	}

	
	vec = evecs[k];





	Matrix<T> QQ(M-con);
	Q.Unity(); QQ.Unity();
	complex<T> tau = orthQ(QQ,vec);

	for(int i=con;i<M;i++){
	for(int j=con;j<M;j++){
		Q( QQ(i-con, j-con) , i,j);
	}
	}

	H = Q.Hermitian()*H;
	H = H*Q;


	for(int j = M-1; j>con+2; j--){

	Matrix<T> U(j-1-con);
	vector<T> z(j-1-con); for(int k = con+0;k<j-1;k++){z[k-con] = conj( H(j,k+1) );}
	normalize(z);
	tau = orthU(U,z);

double tmp = 0;
for(int i=0;i<z.size()-1;i++){tmp = tmp + abs(z[i]);}
if(tmp < small/( (double)z.size()-1.0) ){ U.Unity();}	

	vector<vector<T> > Hb(j-1-con);	
	for(int k=con+1;k<j;k++){Hb[k-1-con].resize(M);
	for(int l=0;l<M;l++){Hb[k-1-con][l] = H(l,k); }} 
 	Hb = times(Hb,U,j-1-con);	 
	for(int k=con+1;k<j;k++){
	for(int l=0;l<M;l++){H(Hb[k-1-con][l] , l,k);}}

	vector<vector<T> > Qb(j-1-con);	
	for(int k=con+1;k<j;k++){Qb[k-1-con].resize(M);
	for(int l=0;l<M;l++){Qb[k-1-con][l] = Q(l,k); }} 
 	Qb = times(Qb,U,j-1-con);	 
	for(int k=con+1;k<j;k++){
	for(int l=0;l<M;l++){Q(Qb[k-1-con][l] , l,k);}}


	
	vector<vector<T> > Hc(M);	
	for(int k=0;k<M;k++){Hc[k].resize(j-1-con);
	for(int l=con+1;l<j;l++){Hc[k][l-1-con] = H(l,k);}}
	U = U.Hermitian();
	Hc = times(U,Hc,j-1-con);  
	for(int k=0;k<M;k++){
	for(int l=con+1;l<j;l++){H(Hc[k][l-1-con] , l,k);}}

	}



}*/

};

#endif

