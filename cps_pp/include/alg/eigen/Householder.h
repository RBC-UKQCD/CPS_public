#ifndef HOUSEHOLDER_H
#define HOUSEHOLDER_H

#include <cstdlib>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <complex>
#include <algorithm>
#include "RandomMatrix.h"


using namespace std;

/** Comparison function for finding the max element in a vector **/
template <class T> bool cf(T i, T j) { return abs(i) < abs(j); }


/** 
	Calculate a real Givens angle 
					**/
template <class T> inline void Givens_calc(T y, T z, T &c, T &s){
	double mz = (double)abs(z);
	if(mz==0.0){
		c = 1; s = 0;
	}
	if(mz >= (double)abs(y)){
		T t = -y/z;
		s = (T)1.0 / sqrt ((T)1.0 + t * t);
		c = s * t;
	}
	else{
		T t = -z/y;
		c = (T)1.0 / sqrt ((T)1.0 + t * t);
		s = c * t;
	}
}

template <class T> inline void Givens_mult(Matrix<T> &A,  int i, int k, T c, T s, int dir){

int q = A.dim;

if(dir == 0){
for(int j=0;j<q;j++){
 T nu = A.A[i][j];
 T w  = A.A[k][j];
 A.A[i][j] = (c*nu + s*w);
 A.A[k][j] = (-s*nu + c*w);
}
}

if(dir == 1){
for(int j=0;j<q;j++){
 T nu = A.A[j][i];
 T w  = A.A[j][k];
 A.A[j][i] = (c*nu - s*w);
 A.A[j][k] = (s*nu + c*w);
}
}

}

/**
	from input = x;
	Compute the complex Householder vector, v, such that
		P = (I - b v transpose(v) )
		b = 2/v.v
		    
		P | x |    | x | k = 0
		  | x |    | 0 | 
		  | x | =  | 0 |
		  | x |    | 0 | j = 3
		  | x |	   | x |

	These are the "Unreduced" Householder vectors.

**/

template <class T> inline void Householder_vector(vector<T> input, int k, int j, vector<T> &v, T &beta){


int N = input.size();
T m = *max_element(input.begin() + k, input.begin() + j + 1, cf<T> );

if(abs(m) > 0.0){
T alpha = 0;

	for(int i=k; i<j+1; i++){
		v[i] = input[i]/m;
	alpha = alpha + v[i]*conj(v[i]);
	}
	alpha = sqrt(alpha);
	beta = (T)1.0/(alpha*(alpha + abs(v[k]) ));
	if(abs(v[k]) > 0.0){v[k] = v[k] + (v[k]/abs(v[k]))*alpha;}
	else{v[k] = -alpha;}
}
else{for(int i=k; i<j+1; i++){v[i] = 0.0;} }

}

/**
	from input = x;
	Compute the complex Householder vector, v, such that
		P = (I - b v transpose(v) )
		b = 2/v.v
		    
		Px = alpha*e_dir

	These are the "Unreduced" Householder vectors.

**/

template <class T> inline void Householder_vector(vector<T> input, int k, int j, int dir, vector<T> &v, T &beta){


int N = input.size();
T m = *max_element(input.begin() + k, input.begin() + j + 1, cf);

if(abs(m) > 0.0){
T alpha = 0;

	for(int i=k; i<j+1; i++){
		v[i] = input[i]/m;
	alpha = alpha + v[i]*conj(v[i]);
	}
	alpha = sqrt(alpha);
	beta = 1.0/(alpha*(alpha + abs(v[dir]) ));
	if(abs(v[dir]) > 0.0){v[dir] = v[dir] + (v[dir]/abs(v[dir]))*alpha;}
	else{v[dir] = -alpha;}
}
else{for(int i=k; i<j+1; i++){v[i] = 0.0;} }

}


/**
	Compute the product PA if trans = 0
			    AP if trans = 1
		P = (I - b v transpose(v) )
		b = 2/v.v
		start at element l of matrix A
		v is of length j - k + 1 of v are nonzero

**/


template <class T> inline void Householder_mult(Matrix<T> &A , vector<T> v, T beta, int l, int k, int j, int trans){

if(abs(beta) > 0.0){
int N = A.dim;
Matrix<T> tmp(N,0); 
T s;
for(int p=l; p<N; p++){
	s = 0;
	if(trans==0){for(int i=k;i<j+1;i++){ s = s + conj(v[i-k])*A(i,p);}}
	else{for(int i=k;i<j+1;i++){ s = s + v[i-k]*A(p,i);}}
	s = beta*s;
	if(trans==0){for(int i=k;i<j+1;i++){ tmp(tmp(i,p) - s*v[i-k], i,p);}}
	else{for(int i=k;i<j+1;i++){ tmp(tmp(p,i) - s*conj(v[i-k]), p,i);}}
}
for(int p=l; p<N; p++){
	if(trans==0){for(int i=k;i<j+1;i++){ A(A(i,p) + tmp(i,p), i,p);}}
	else{for(int i=k;i<j+1;i++){ A(A(p,i) + tmp(p,i), p,i);}}
}

}

}


/**
	Compute the product PA if trans = 0
			    AP if trans = 1
		P = (I - b v transpose(v) )
		b = 2/v.v
		start at element l of matrix A
		v is of length j - k + 1 of v are nonzero
		A is tridiagonal
**/


template <class T> inline void Householder_mult_tri(Matrix<T> &A , vector<T> v, T beta, int l, int M, int k, int j, int trans){

if(abs(beta) > 0.0){
int N = A.dim;
Matrix<T> tmp(N,0); 
T s;
for(int p=l; p<M; p++){
	s = 0;
	if(trans==0){for(int i=k;i<j+1;i++){ s = s + conj(v[i-k])*A(i,p);}}
	else{for(int i=k;i<j+1;i++){ s = s + v[i-k]*A(p,i);}}
	s = beta*s;
	if(trans==0){for(int i=k;i<j+1;i++){ tmp(tmp(i,p) - s*v[i-k], i,p);}}
	else{for(int i=k;i<j+1;i++){ tmp(tmp(p,i) - s*conj(v[i-k]), p,i);}}
}
for(int p=l; p<M; p++){
	if(trans==0){for(int i=k;i<j+1;i++){ A(A(i,p) + tmp(i,p), i,p);}}
	else{for(int i=k;i<j+1;i++){ A(A(p,i) + tmp(p,i), p,i);}}
}

}

}

#endif
