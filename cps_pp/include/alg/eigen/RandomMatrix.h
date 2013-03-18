#ifndef RANDOMMATRIX_H
#define RANDOMMATRIX_H

#include <cstdlib>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <complex>
#include "Matrix.h"

/**Random real vector**/

template <class T> void RandomReal(vector<T> &A, int seed){

srand48(seed);for(int i=0;i<A.size();i++){A[i] = drand48();}

}

/**Random complex vector**/
template <class T> void RandomComplex(vector<T> &A, int seed){

srand48(seed);for(int i=0;i<A.size();i++){A[i] = T( drand48(),drand48() );}

}

/**Random unit vector**/
template <class T> void RandomComplexUnit(vector<T> &A, int seed){

srand48(seed);for(int i=0;i<A.size();i++){A[i] = T( drand48(),drand48() );}
T norm = norm(A); for(int i=0;i<A.size();i++){A[i] /= norm;}

}

/**Random unit vector**/
template <class T> void RandomRealUnit(vector<T> &A, int seed){

srand48(seed);for(int i=0;i<A.size();i++){A[i] = drand48();}
T norm = norm(A); for(int i=0;i<A.size();i++){A[i] /= norm;}

}

/**Random dense vector(vector(**/
template <class T> void RandomReal(vector<vector<T> > &A, int seed){

int M = A.size();
srand48(seed);
for(int i=0;i<M;i++){
for(int j=0;j<M;j++){ A[i][j] = drand48(); }}

}

/**Random dense vector(vector(**/
template <class T> void RandomComplex(vector<vector<T> > &A, int seed){

int M = A.size();
srand48(seed);
for(int i=0;i<M;i++){
for(int j=0;j<M;j++){ A[i][j] = T( drand48(),drand48() ); }}

}

/**Random dense matrix**/

template <class T> void RandomComplex(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
for(int i=0;i<M;i++){
for(int j=0;j<M;j++){A( T( drand48(),drand48() )  ,  i , j);}}

}

/**Random real matrix**/

template <class T> void RandomReal(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
for(int i=0;i<M;i++){
for(int j=0;j<M;j++){A( drand48()  ,  i , j);}}

}

/**Random unitary matrix
 
Adag A = A Adag = I

**/

template <class T> void RandomUnitary(Matrix<T> &A, int seed){

int M = A.dim;
vector< vector<T> > q(M);

srand48(seed);
for(int i=0;i<M;i++){q[i].resize(M);
for(int j=0;j<M;j++){q[i][j] = T( drand48(),drand48() );}}

for(int j=0;j<M;j++){
for(int i=0;i<j;i++){
q[j] = sub(q[j]   ,   times( inner(q[i],q[j])/norm2(q[i]), q[i] )    );
}
normalize(q[j]);
}
Matrix<T> AT(q); 
A = AT;

}

/**Random Hermitian matrix
 
Adag = A

**/
template <class T> void RandomHermitian(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
for(int i=0;i<M;i++){
for(int j=i;j<M;j++){

if(j != i){A( T( drand48(),drand48() )  ,  i , j);}
else{ A( T( drand48(), 0.0 )  ,  i , j); }

}}

for(int i=0;i<M;i++){
for(int j=i;j<M;j++){A( conj( A(i,j) )  , j, i);}}

}

/**Random Symmetric matrix
 
Atranspose = A

**/
template <class T> void RandomSymmetric(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
for(int i=0;i<M;i++){
for(int j=i;j<M;j++){A( T( drand48(),drand48() )  ,  i , j);}}

for(int i=0;i<M;i++){
for(int j=i;j<M;j++){A( A(i,j)  , j, i);}}

}

/**Random Real Symmetric matrix
 
Atranspose = A

**/
template <class T> void RandomRealSymmetric(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
for(int i=0;i<M;i++){
for(int j=i;j<M;j++){A( drand48() ,  i , j);}}

for(int i=0;i<M;i++){
for(int j=i;j<M;j++){A( A(i,j)  , j, i);}}

}

/**Random diagonal matrix
 
x  0  0  0
0  x  0  0
0  0  x  0
0  0  0  x

x is complex
**/
template <class T> void RandomComplexDiagonal(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
A.Fill(0.0);
for(int j=0;j<M;j++){A( T( drand48(),drand48() )  ,  j , j);}

}

/**Random diagonal matrix
 
x  0  0  0
0  x  0  0
0  0  x  0
0  0  0  x

x is real
**/
template <class T> void RandomRealDiagonal(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
A.Fill(0.0);
for(int j=0;j<M;j++){A( T( drand48(),0.0 )  ,  j , j);}

}

/**Random unitary diagonal matrix
 
x  0  0  0
0  x  0  0
0  0  x  0
0  0  0  x

x is e^it

so Adag A = A Adag = I
**/
template <class T> void RandomUnitaryDiagonal(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
A.Fill(0.0);
for(int j=0;j<M;j++){
T tmp = T( -1.0 + 2.0*drand48(), -1.0 + 2.0*drand48() );
A( tmp/abs(tmp)  ,  j , j);}

}

/**Random Upper Hessenberg matrix
 
x  x  x  x  x
x  x  x  x  x
0  x  x  x  x
0  0  x  x  x
0  0  0  x  x

**/
template <class T> void RandomHessenberg(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
A.Fill(0.0);
for(int i=0;i<M;i++){
for(int j=i;j<M;j++){A( T( drand48(),drand48() )  ,  i , j);}}
for(int i=1;i<M;i++){A( T( drand48(),drand48() ) ,  i, i-1);}

}

template <class T> void RandomRealHessenberg(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
A.Fill(0.0);
for(int i=0;i<M;i++){
for(int j=i;j<M;j++){A(  drand48() ,  i , j);}}
for(int i=1;i<M;i++){A(  drand48() ,  i, i-1);}

}

/**Random Tridiagonal matrix
 
x  x  0  0  0
x  x  x  0  0
0  x  x  x  0
0  0  x  x  x
0  0  0  x  x

**/
template <class T> void RandomTridiagonal(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
A.Fill(0.0);

for(int i=0;i<M;i++){A( T( drand48(),drand48() ) ,  i, i);}
for(int i=1;i<M;i++){A( T( drand48(),drand48() ) ,  i, i-1);}
for(int i=0;i<M-1;i++){A( T( drand48(),drand48() ) ,  i, i+1);}

}


/**Random Symmetric Tridiagonal matrix
 
x  x  0  0  0
x  x  x  0  0
0  x  x  x  0
0  0  x  x  x
0  0  0  x  x

**/
template <class T> void RandomSymmTridiagonal(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
A.Fill(0.0);
for(int i=0;i<M;i++){A( T( drand48(),drand48() ) ,  i, i);}

for(int i=0;i<M-1;i++){
	T off = T( drand48(),drand48() );
	 A( off ,  i+1, i);
	 A( off ,  i, i+1);
}

}

/**Random Hermitian Tridiagonal matrix
 
x  x  0  0  0
x  x  x  0  0
0  x  x  x  0
0  0  x  x  x
0  0  0  x  x

**/
template <class T> void RandomHermTridiagonal(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
A.Fill(0.0);
for(int i=0;i<M;i++){A( T( drand48(),0 ) ,  i, i);}

for(int i=0;i<M-1;i++){
	T off = T( drand48(),drand48() );
	 A( off ,  i+1, i);
	 A( conj(off) ,  i, i+1);
}

}

/**Random Hermitian Tridiagonal matrix
 
x  x  0  0  0
x  x  x  0  0
0  x  x  x  0
0  0  x  x  x
0  0  0  x  x

**/
template <class T> void RandomRealSymmTridiagonal(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
A.Fill(0.0);
for(int i=0;i<M;i++){A( drand48() ,  i, i);}

for(int i=0;i<M-1;i++){
	T off = drand48();
	 A( off ,  i+1, i);
	 A( off ,  i, i+1);
}

}

/**Random Upper Triangular matrix
 
x  x  x  x  x
0  x  x  x  x
0  0  x  x  x
0  0  0  x  x
0  0  0  0  x

**/
template <class T> void RandomTriangular(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
A.Fill(0.0);
for(int i=0;i<M;i++){
for(int j=i;j<M;j++){A( T( drand48(),drand48() )  ,  i , j);}}

}

/**Random Reduced matrix
 
    x  x  x  x  x  x  x  x  x
    0  x  x  x  x  x  x  x  x
    0  0  x  x  x  x  x  x  x
    0  0  0  x  x  x  x  x  x
    0  0  0  x  x  x  x  x  x
    0  0  0  0  0  0  x  x  x
    0  0  0  0  0  0  0  x  x
    0  0  0  0  0  0  0  0  x

2 by 2 blocks at random on the diagonal
**/
template <class T> void RandomReduced(Matrix<T> &A, int seed){

int M = A.dim;
srand48(seed);
A.Fill(0.0);
for(int i=0;i<M;i++){
for(int j=i;j<M;j++){A( T( drand48(),drand48() )  ,  i , j);}}

for(int i=1;i<M;i++){ if(drand48() > 0.5){A( T( drand48(),drand48() ) ,  i, i-1); i++;}  }


}
#endif
/*
int main(void){

int M = 5;

Matrix<complex<double> > A(M) , B(M), C(M);

RandomUnitary(A,1234567);
B = A.Hermitian();
A = B*A;
cout << "Unitary " << endl;
A.Out();

RandomHermitian(A,1234567);
B = A.Hermitian();
A = A - B;
cout << "Hermitian " << endl;
A.Out();

RandomSymmetric(A,1234567);
B = A.Transpose();
A = A - B;
cout << "Symmetric " << endl;
A.Out();

RandomComplexDiagonal(A,1234567);
cout << "Complex Diagonal " << endl;
A.Out();

RandomRealDiagonal(A,1234567);
cout << "Real Diagonal " << endl;
A.Out();

RandomUnitaryDiagonal(A,1234567);
B = A.Hermitian();
A = B*A;
cout << "Unitary Diagonal" << endl;
A.Out();

RandomHessenberg(A,1234567);
cout << "Hessenberg" << endl;
A.Out();

RandomTriangular(A,1234567);
cout << "Triangular" << endl;
A.Out();

RandomReduced(A,1234567);
cout << "Reduced" << endl;
A.Out();

return 0;
}
*/
