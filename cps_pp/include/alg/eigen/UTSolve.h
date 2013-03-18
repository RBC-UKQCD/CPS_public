#ifndef UTSOLVE_H
#define UTSOLVE_H

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


/** 
Find eigenvectors of an upper quasi triangular matrix U (N x N)

eg.

    x  x  x  x  x  x  x  x  x
    0  x  x  x  x  x  x  x  x
    0  0  x  x  x  x  x  x  x
    0  0  0  r  r  x  x  x  x
    0  0  0  r  r  x  x  x  x
    0  0  0  0  0  0  x  x  x
    0  0  0  0  0  0  0  x  x
    0  0  0  0  0  0  0  0  x

It can handle the r r block if you tell it where to look with trows, here trows[3] = 1 rest are zero.
                  r r
I am aware that this is stupid.


The Eigenvalues are specified in uvals (N).
Gives the vector of vectors uvecs (N x N) of corresponding eigenvectors.

So many special cases...
**/
template <class T> void UTSymmEigenvectors(Matrix<T > &U, vector<int> trows, vector<T > uvals, vector<vector<T > > &uvecs){

int N = U.dim;
uvecs.resize(N);
T a,b,c,d,L1,L2,Norm,lab;

for(int eno = 0;eno<N;eno++){

if(trows[eno] == 1){
	uvecs[eno].resize(N); 
	uvecs[eno+1].resize(N); 
	for(int j=0;j<N;j++){uvecs[eno][j] = 0; uvecs[eno+1][j] = 0;}

		a = U(eno,eno);   	b = U(eno,eno+1);
		c = U(eno+1,eno);     	d = U(eno+1,eno+1);
		L1 = uvals[eno];
		L2 = uvals[eno+1];

		/**
		| a  b | |x|  =  L|x|	(eno)
		| c  d | |y|	  |y|	(eno+1)

	      solutions: 
			y = 1; x = b / (L-a)
	      or	x = 1; y = c / (L-d)

	      **/

	lab = (L1-a) / b;
	Norm = 1.0/sqrt(1.0 + lab*lab ); 
	uvecs[eno][eno] = Norm*lab; 
	uvecs[eno][eno+1] = Norm;
	
	lab = (L2-a) / b;
	Norm = sqrt(1.0 + lab*lab ); 
	uvecs[eno+1][eno] = Norm*lab; 
	uvecs[eno+1][eno+1] = Norm;
	  
	eno++;
}
else{
    uvecs[eno].resize(N); for(int j=0;j<N;j++){uvecs[eno][j] = 0;}
    uvecs[eno][eno] = 1;
}
  
}

}

template <class T> void UTSymmEigensystem(Matrix<T > &U, vector<int> trows, vector<T > uvals, vector<vector<T > > &uvecs){

int N = U.dim;
uvecs.resize(N);
uvals.resize(N);
T a,b,c,d,L1,L2,apd,amd,bc, Norm, lab;

for(int eno = 0;eno<N;eno++){

if(trows[eno] == 1){
	uvecs[eno].resize(N); 
	uvecs[eno+1].resize(N); 
	for(int j=0;j<N;j++){uvecs[eno][j] = 0; uvecs[eno+1][j] = 0;}

		a = U(eno,eno);   	b = U(eno,eno+1);
		c = U(eno+1,eno);     	d = U(eno+1,eno+1);

		apd = a+d;
		amd = a-d;
		bc = (T)4.0*b*c;
		L1 = (T)0.5*( apd + sqrt(amd*amd+bc) ); 
		L2 = (T)0.5*( apd - sqrt(amd*amd+bc) );
		uvals[eno] = L1;
		uvals[eno+1] = L2;
		/**
		| a  b | |x|  =  L|x|	(eno)
		| c  d | |y|	  |y|	(eno+1)

	      solutions: 
			y = 1; x = b / (L-a)
	      or	x = 1; y = c / (L-d)

	      **/

	lab = (L1-a) / b;
	Norm = 1.0/sqrt(1.0 + lab*lab ); 
	uvecs[eno][eno] = Norm*lab; 
	uvecs[eno][eno+1] = Norm;
	
	lab = (L2-a) / b;
	 Norm = sqrt(1.0 + lab*lab ); 
	uvecs[eno+1][eno] = Norm*lab; 
	uvecs[eno+1][eno+1] = Norm;
	  
	eno++;
}
else{
    uvals[eno] = U(eno,eno);
    uvecs[eno].resize(N); for(int j=0;j<N;j++){uvecs[eno][j] = 0;}
    uvecs[eno][eno] = 1;
}
  
}

}

template <class T> void UTeigenvectors(Matrix<T > &U, vector<int> trows, vector<T > uvals, vector<vector<T > > &uvecs){

int N = uvals.size();
uvecs.resize(N);


for(int i=0; i<N; i++){
uvecs[i].resize(N);for(int j=0;j<N;j++){uvecs[i][j] = 0;}
uvecs[i][i] = 1;
}



T a,b,c,d,x1,x2,L1,L2;
/// Modified backsubstitution
for(int eno = 0;eno<N;eno++){
  //cerr << "utsolve " << eno << endl;
	int tag = 0;
///If find a 2 x 2 block compute the last two components seperately
	if(trows[eno] == 1){
		//for(int i=0;i<N;i++){uvecs[eno-1][i] = 0;} 
		uvecs[eno-1][eno-1] = 0;
		
		a = U(eno-1,eno-1);   b = U(eno-1,eno);
		c = U(eno,eno-1);     d = U(eno,eno);


L1 = uvals[eno-1];
L2 = uvals[eno];

/**
| a  b | |x|  =  L|x|	(eno-1)
| c  d | |y|	  |y|	(eno)

solutions: 
	y = 1; x = b / (L-a)
or	x = 1; y = c / (L-d)

avoid dividing by small number choose whichever (L-?) is larger
**/
if( abs(L1-a) > abs(L1-d) ){
	uvecs[eno-1][eno-1] = b/(L1-a); 
	uvecs[eno-1][eno] = 1;
}
else{
	uvecs[eno-1][eno-1] = 1;
	uvecs[eno-1][eno] = c/(L1-d); 
}

if( abs(L2-a) > abs(L2-d) ){
	uvecs[eno][eno-1] = b/(L2-a); 
	uvecs[eno][eno] = 1;
}
else{
	uvecs[eno][eno-1] = 1;
	uvecs[eno][eno] = c/(L2-d); 
}

tag = 2;	///precomputed last two components
		
	}
	for(int i=(eno-tag);i>=0;i--){
		///Simple eval. if tag = 2 do two evecs at once
		if(trows[i] == 0){

		for(int j=i+1;j<=eno;j++){
			uvecs[eno][i] = uvecs[eno][i] - U(i,j)*uvecs[eno][j];
		if(tag==2){uvecs[eno-1][i] = uvecs[eno-1][i] - U(i,j)*uvecs[eno-1][j];}			
			}
		if(i<eno){
			if( abs(uvals[i] - uvals[eno]) != 0.0){ uvecs[eno][i] = uvecs[eno][i]/( uvals[i] - uvals[eno]); } 
			else{uvecs[eno][i] = 0.0;}
		if(tag==2){
			if( abs(uvals[i] - uvals[eno-1]) != 0.0){ uvecs[eno-1][i] = uvecs[eno-1][i]/( uvals[i] - uvals[eno-1]); } 
			else{uvecs[eno-1][i] = 0.0;}
			}
		}
		}
		///Double eval. if tag = 2 do two evecs at once
		else{
			T y1=0, y2=0, y3=0, y4=0;
			for(int j=i+1;j<=eno;j++){
				y1 = y1 - U(i-1,j)*uvecs[eno][j];
				y2 = y2 - U(i,j)*uvecs[eno][j];
				if(tag==2){y3 = y3 - U(i-1,j)*uvecs[eno-1][j];
				y4 = y4 - U(i,j)*uvecs[eno-1][j];}
			}
				a = U(i-1,i-1)-uvals[eno]; b = U(i-1,i);
				c = U(i,i-1); d = U(i,i)-uvals[eno];

				x1 = (-d*y1+b*y2)/(b*c-a*d);
				x2 = (-c*y1+a*y2)/(a*d-b*c);
				uvecs[eno][i] = x2;
				uvecs[eno][i-1] = x1;
			if(tag==2){
				a = U(i-1,i-1)-uvals[eno-1]; b = U(i-1,i);
				c = U(i,i-1); d = U(i,i)-uvals[eno-1];

				x1 = (-d*y3+b*y4)/(b*c-a*d);
				x2 = (-c*y3+a*y4)/(a*d-b*c);

				uvecs[eno-1][i] = x2;
				uvecs[eno-1][i-1] = x1;
			}
			i--;
		}

		}
}


/// Normalize eigenvectors
	for(int eno=0;eno<N;eno++) normalize(uvecs[eno]);

}

/** 
Find eigenvalues of an upper quasi triangular matrix U (N x N)

eg.

    x  x  x  x  x  x  x  x  x 
    0  x  x  x  x  x  x  x  x 
    0  0  x  x  x  x  x  x  x 
    0  0  0  r  r  x  x  x  x 
    0  0  0  r  r  x  x  x  x 
    0  0  0  0  0  x  x  x  x 
    0  0  0  0  0  0  x  x  x 
    0  0  0  0  0  0  0  x  x 
    0  0  0  0  0  0  0  0  x 

It can handle the r r block if you tell it where to look with trows, here trows[3] = 1 rest are zero.
                  r r
I am aware that this is stupid.

The Eigenvalues are specified in uvals (N).
**/

template <class T> void UTeigenvalues(Matrix<T > &U, vector<int> trows, vector<T > &uvals){

int N = U.dim;
uvals.resize(N); 
for(int eno = 0;eno<N;eno++){
	if(trows[eno] == 1){
		T apd,amd,bc;
		apd = U(eno-1,eno-1) + U(eno,eno);
		amd = U(eno-1,eno-1) - U(eno,eno);
		bc = (T)4.0*U(eno-1,eno)*U(eno,eno-1);
		uvals[eno-1]   = (T)0.5*( apd + sqrt(amd*amd+bc) ); 
		uvals[eno] = (T)0.5*( apd - sqrt(amd*amd+bc) );
	}
	else{
		uvals[eno] = U(eno,eno);
	}	
}
}


template <class T> void UTeigensystem(Matrix<T> &U, vector<int> trows, vector<T> &uvals, vector<vector<T> > &uvecs){

UTeigenvalues(U, trows, uvals);
UTeigenvectors(U, trows, uvals, uvecs);

}

/**
	Determine the structure of the Reduced Matrix U
e.g.
	x  x  x  x  x
	x  x  x  x  x
	0  0  x  x  x
	0  0  0  x  x
	0  0  0  x  x

	here trows = {1 , 0 , 0,  1,  0}

**/

template <class T> void UTrows(Matrix<T> &U, vector<int> &trows){

fill(trows,0);
T nrm = abs( U.Norm() ); 
int M = U.dim;
for(int i=1;i<M;i++){ 
	if( (T)( abs(nrm) + abs( U(i, i-1) ) ) > (T)abs(nrm) ){
		trows[i] = 1;
	}  
}

}

#endif








