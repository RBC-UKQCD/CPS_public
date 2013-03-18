#ifndef MATRIX_H
#define MATRIX_H

#include <qdp.h>
#include <qdp_util.h>

#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>

#include <bfm.h>
#include <bfm_qdp.h>




#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <complex>
#include <typeinfo>

#define NUM_BFM_THREADS (64)
using namespace std;

typedef Fermion_t bfm_fermion [2];
/** print a fermion **/
/*void print_fermion(Fermion in){

  for(int s=0;s<Nd;s++){
  for(int c=0;c<Nc;c++){
  
	SpinVector tmp1 = peekColor(in,c);
	Complex C = peekSpin(tmp1,s);
	QDPIO::cout << C << "  ";
	
  }
  QDPIO::cout << endl;
  }
  QDPIO::cout << endl;

}

void print_fermion2(Fermion in){
  QDPIO::cout << endl;
  for(int s=0;s<Nd;s++){
  for(int c=0;c<Nc;c++){
  
	SpinVector tmp1 = peekColor(in,c);
	Complex C = peekSpin(tmp1,s);
	QDPIO::cout << "(" << s << "," << c << ") = " << C << endl;
	
  }
  }
  QDPIO::cout << endl;

}*/


/** 
	Do complex stuff to floats/doubles 
					**/
//float conj(float z){return z;}
static float real(float z){return z;}
static float imag(float z){return 0;}

//double conj(double z){return z;}
static double real(double z){return z;}
static double imag(double z){return 0;}


/** Chroma complex to C++ complex(double) **/
template <class T> complex<T> toComplex(Complex z){return complex<T>(toDouble(real(z) ) , toDouble(imag(z) ) ); }
template <class T> T toComplex(Real z){return (T)toDouble(z); }

/** C++ complex(double) to Chroma complex **/
template <class T> Complex toComplex(T z){return cmplx(Real(real(z) ) , Real(imag(z) ) ); }


/** Sign function **/
template <class T> T sign(T p){return ( p/abs(p) );}


/**
*****************************************
*					*
*	Complex Vector operations	*
*					*
*****************************************
**/
/**Conj of a vector **/
template <class T> vector<T> conj(vector<T> p){
	vector<T> q(p.size());
	for(int i=0;i<p.size();i++){q[i] = conj(p[i]);}
	return q;
}

/** Norm of a vector**/
template <class T> T norm(vector<T> p){
	T sum = 0;
	for(int i=0;i<p.size();i++){sum = sum + p[i]*conj(p[i]);}
	return abs(sqrt(sum));
}

/** Norm squared of a vector **/
template <class T> T norm2(vector<T> p){
	T sum = 0;
	for(int i=0;i<p.size();i++){sum = sum + p[i]*conj(p[i]);}
	return abs((sum));
}

/** Sum elements of a vector **/
template <class T> T trace(vector<T> p){
	T sum = 0;
	for(int i=0;i<p.size();i++){sum = sum + p[i];}
	return sum;
}

/** Fill a vector with constant c **/
template <class T> void fill(vector<T> &p, T c){
	for(int i=0;i<p.size();i++){p[i] = c;}
}

/** Normalize a vector **/
template <class T> void normalize(vector<T> &p){
	T m = norm(p);
	if( abs(m) > 0.0) for(int i=0;i<p.size();i++){p[i] /= m;}
}

/** print components of a vector in mathematica readable (list) form **/
template <class T> void print(vector<T> p, int w){
	cout << fixed << setprecision(w) << setw(w) << "{ ";
	for(int i=0;i<p.size();i++){
		i==p.size()-1 ? cout << real(p[i]) << " +I*" << imag(p[i]) << " }; ":cout << real(p[i]) << " +I*" << imag(p[i]) << " , ";
	}
	cout << endl;
}

template <class T> void print(vector<T> p){
	print(p,3);
}

/** print components of a vector of vectors in mathematica readable (matrix) form **/
template <class T> void print(vector< vector<T> > Q){
	int N = Q.size();
	cout << fixed << setprecision(14) << setw(14) << "{ ";
	for(int i = 0;i<N;i++){
	int M = Q[i].size();
	cout << "{ ";
	for(int j = 0;j<M;j++){
		j==M-1?cout << real(Q[i][j]) << " + I*" << imag(Q[i][j]) << " ":cout << real(Q[i][j]) << " + I*" << imag(Q[i][j]) << " , ";
	}
	i==N-1?cout << "}":cout << "},"<< endl;
	}
	cout << "};" << endl;
}
/** print components of a vector of vectors in mathematica readable (matrix) form **/
template <class T> void printSimpleReal(vector< vector<T> > Q){
	int N = Q.size();
	cout << fixed << setprecision(3) << scientific;
	for(int i = 0;i<N;i++){
	int M = Q[i].size();
	for(int j = 0;j<M;j++){
		cout << real(Q[i][j]) << " , ";
	}
	cout << endl;
	}
	cout << endl;
}

/** Vector by scalar **/
template <class T, class U> vector<T> times(vector<T> p, U s){
	for(int i=0;i<p.size();i++){p[i] *= s;}
	return p;
}
template <class T, class U> vector<T> times(U s, vector<T> p){
	for(int i=0;i<p.size();i++){p[i] *= s;}
	return p;
}

/** inner product of a and b = conj(a) . b **/
template <class T> T inner(vector<T> a, vector<T> b){
	T m = 0.;
	for(int i=0;i<a.size();i++){m = m + conj(a[i])*b[i];}
	return m;
}

/** sum of a and b = a + b **/
template <class T> vector<T> add(vector<T> a, vector<T> b){
	vector<T> m(a.size());
	for(int i=0;i<a.size();i++){m[i] = a[i] + b[i];}
	return m;
}

/** sum of a and b = a - b **/
template <class T> vector<T> sub(vector<T> a, vector<T> b){
	vector<T> m(a.size());
	for(int i=0;i<a.size();i++){m[i] = a[i] - b[i];}
	return m;
}


/** 
*********************************
*				*
*	Square matricies	*
* 				*
*********************************
**/
/*
	Why use pointers and not vector<vector<T> > ?
	Pointers -> confusing, must pass matricies by reference, freeing memory problems ...
	Vectors  -> So so so so so slow

	Want computations to finish before computer becomes obsolete -> use pointers
*/
	
template <class T> class Matrix {
public:
int dim;		/// dimension
T **A;			/// data

/**Default**/
Matrix(void){
dim=0;
try{
   A = new T*[dim];
}
catch (bad_alloc& ba){
    cerr << "bad_alloc caught: " << ba.what() << endl;
}
}

/** Construct new matrix unspecified entries**/
Matrix(int d){
dim = d; 
try{
   A = new T*[dim];
}
  catch (bad_alloc& ba)
  {
    cerr << "bad_alloc caught: " << ba.what() << endl;
  }
for(int i=0;i<dim;i++){
try{
A[i] = new T [dim];
}
  catch (bad_alloc& ba)
  {
    cerr << "bad_alloc caught: " << ba.what() << endl;
  }
}
}

/** Construct new matrix entries all = val**/
Matrix(int d, T val){ 
dim = d; 
try{
A = new T*[dim];
}
  catch (bad_alloc& ba)
  {
    cerr << "bad_alloc caught: " << ba.what() << endl;
  }
for(int i=0;i<dim;i++){
try{
A[i] = new T [dim];
}
  catch (bad_alloc& ba)
  {
    cerr << "bad_alloc caught: " << ba.what() << endl;
  }
}
Fill(val);
}

/** Destroy a matrix **/
~Matrix(){
	for(int i=0;i<dim;++i){ delete[] A[i];}
	delete[] A;
}

/** Resize a matrix destroying everything in it **/
void resize(int N){
	for(int i=0;i<dim;++i){ delete[] A[i];}
	delete[] A;
dim = N; 
try{
   A = new T*[dim];
}
  catch (bad_alloc& ba)
  {
    cerr << "bad_alloc caught: " << ba.what() << endl;
  }
for(int i=0;i<dim;i++){
try{
A[i] = new T [dim];
}
  catch (bad_alloc& ba)
  {
    cerr << "bad_alloc caught: " << ba.what() << endl;
  }
}
}

/** initialize a matrix with a vector(vector(
my convection (for reasons that are not obvious here but made writing some loops easier)
**/
Matrix(vector< vector<T> > B){ 
dim = B.size(); 
try{
   A = new T*[dim];
}
  catch (bad_alloc& ba)
  {
    cerr << "bad_alloc caught: " << ba.what() << endl;
  }
for(int i=0;i<dim;i++){
try{
A[i] = new T [dim];
}
  catch (bad_alloc& ba)
  {
    cerr << "bad_alloc caught: " << ba.what() << endl;
  }
}
for(int i = 0; i<dim; i++){
	for(int j=0;j<dim;j++){A[i][j] = B[j][i];}
	}
		}

/** Fill matricies with the value x **/
void Fill(T x){
		for(int i=0;i<dim;i++){
			for(int j=0;j<dim;j++){A[i][j] = x;} 
		} 
	}

/** print matrix in Mathematica Readable form : default precision **/
void Out(){
		OutSimpleReal(6);
	}

void OutSimpleReal(int precis){

		for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
		cout << fixed << setw(precis) << setprecision(precis) << real(A[i][j]) << " , ";
		}
		cout << endl;			
		} 
		cout << endl;
	}
	
/** print matrix in Mathematica Readable form : user precision **/
void Out(int precis){
		cout << "{";
		for(int i=0;i<dim;i++){
			cout << "{";
			for(int j=0;j<dim;j++){
j == (dim-1) ? cout << fixed << setprecision(precis)<< setw(precis) <<  real(A[i][j]) << " + I*" << imag(A[i][j]) << " " : cout << fixed << setprecision(precis) << setw(precis) << real(A[i][j]) << " + I*" << imag(A[i][j]) << "\t,"; } 
			i == (dim-1) ? cout << "}};" : cout << "},";
			cout << endl;
		} 
	}

/** print matrix in Mathematica Readable form : scientific **/
void OutSci(int precis){
		cout << "{";
		for(int i=0;i<dim;i++){
			cout << "{";
			for(int j=0;j<dim;j++){
j == (dim-1) ? cout << fixed << scientific << setprecision(precis)<< setw(precis) <<  real(A[i][j]) << " + I*" << imag(A[i][j]) << " " : cout << fixed << scientific << setprecision(precis) << setw(precis) << real(A[i][j]) << " + I*" << imag(A[i][j]) << "\t,"; } 
			i == (dim-1) ? cout << "}};" : cout << "},";
			cout << endl;
		} 
	}

/** print matrix in c form **/
void UseOut(){
		
		for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
cout << fixed << scientific <<  setprecision(14) << setw(14) << "H(complex<double> ( " <<real(A[i][j]) << "," << imag(A[i][j]) << ") , " << i << " , " << j << " );\n"; } 
			cout << endl;
		} 
	}

/** print matrix in c form **/
void UseOutR(){
		
		for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
cout << fixed << scientific <<  setprecision(14) << setw(14) << "H(" <<real(A[i][j]) <<  " , " << i << " , " << j << " );\n"; } 
			cout << endl;
		} 
	}

/** Transpose of a matrix **/
Matrix<T> Transpose(){
	Matrix C(dim);
	for(int i=0;i<dim;i++){
			for(int j=0;j<dim;j++){C.A[i][j] = A[j][i];} 
			} 
	return C;
}

/** Set Matrix to unit matrix **/
void Unity(){
	for(int i=0;i<dim;i++){
	for(int j=0;j<dim;j++){
i==j?A[i][j] = 1:A[i][j] = 0;
			} 
			} 
}


/** Add C * I to matrix **/
void PlusUnit(T c){
	for(int i=0;i<dim;i++){A[i][i] = A[i][i] + c;} 
}

/** return the Hermitian conjugate of matrix **/
Matrix<T> Hermitian(){
	Matrix C(dim);
	for(int i=0;i<dim;i++){
			for(int j=0;j<dim;j++){C.A[i][j] = conj(A[j][i]);} 
			} 
	return C;
}

/** return diagonal entries as a vector **/
vector<T> diag(){
	vector<T> d(dim);
	for(int i=0;i<dim;i++){d[i] = A[i][i];}
	return d;
}

/** Left multiply by a vector **/
vector<T> LeftMult(vector<T> &B){
	vector<T> C(dim);
	for(int j=0;j<dim;j++){
			T sum = 0.0;
			for(int i=0;i<dim;i++){
				sum += B[i] * A[i][j];}
				C[j] =  sum;
			}
	return C; 
}

/** return 1/diagonal entries as a vector **/
vector<T> inv_diag(){
	vector<T> d(dim);
	for(int i=0;i<dim;i++){d[i] = 1.0/A[i][i];}
	return d;
}
/** sizeof **/
inline int size(Matrix<T> &B){
	return dim;
	} 
/** Matrix Addition **/
inline Matrix<T> operator+ (Matrix<T> &B){
	Matrix C(dim);
	for(int i=0;i<dim;i++){
			for(int j=0;j<dim;j++){C.A[i][j] = A[i][j] +  B.A[i][j];} 
			} 
	return C;
	} 

/** Matrix Subtraction **/
inline Matrix<T> operator- (Matrix<T> &B){
	Matrix C(dim);
	for(int i=0;i<dim;i++){
			for(int j=0;j<dim;j++){C.A[i][j] = A[i][j] -  B.A[i][j];} 
			} 
	return C;
	} 

/** Matrix scalar multiplication **/
inline Matrix<T> operator* (T c){
	Matrix C(dim);
	for(int i=0;i<dim;i++){
			for(int j=0;j<dim;j++){C.A[i][j] = A[i][j]*c;} 
			} 
	return C;
	} 

/** Matrix Matrix multiplication **/
inline Matrix<T> operator* (Matrix<T> &B){
	Matrix C(dim);
	for(int i=0;i<dim;i++){
			for(int j=0;j<dim;j++){
				T sum = 0.0;
				for(int k=0;k<dim;k++){sum += A[i][k]*B.A[k][j];}
				C(sum,i,j);
				}
			}
	return C; 
	} 

/** Matrix vector multiplication **/
inline vector<T> operator* (vector<T> &B){
	vector<T> C(dim);
	for(int i=0;i<dim;i++){
			T sum = 0.0;
			for(int j=0;j<dim;j++){
				sum += A[i][j]*B[j];}
				C[i] =  sum;
			}
	return C; 
	} 

/** Matrix equivalence **/
inline Matrix<T>& operator=(const Matrix<T>& B){
	dim = B.dim;
	for(int i=0;i<dim;i++){
			for(int j=0;j<dim;j++){A[i][j] = B.A[i][j];} 
			}
	return (*this); 
	}

/** Some version of Matrix norm **/
inline T Norm(){

T norm = 0;
for(int i=0;i<dim;i++){
for(int j=0;j<dim;j++){
	norm += abs(A[i][j]);
}}
return norm;

}

/** Some version of Matrix norm **/
inline T LargestDiag(){

T ld = abs(A[0][0]);
for(int i=1;i<dim;i++){
  
  T cf = abs(A[i][i]);
  if(abs(cf) > abs(ld) ){ld = cf;}
  
}

return ld;

}

/** Look for entries on the leading subdiagonal that are smaller than 'small' **/
template <class U> 
int Chop_subdiag(T norm, int offset, U small)
{
  for(int l = dim - 1 - offset; l >= 1; l--) 
  {             		
    if((U)abs(A[l][l - 1]) < (U)small) 
    {
      A[l][l-1]=(U)0.0;
      return l;
    }
  }
  return 0;
}

/** Look for entries on the leading subdiagonal that are smaller than 'small' **/
template <class U> 
int Chop_symm_subdiag(T norm, int offset, U small)
{
  for(int l = dim - 1 - offset; l >= 1; l--) 
  {             		
    if((U)abs(A[l][l - 1]) < (U)small) 
    {
      A[l][l - 1] = (U)0.0;
      A[l - 1][l] = (U)0.0;
      return l;
    }
  }
  return 0;
}

/**Assign a submatrix to a larger one**/
void AssignSubMtx(int row_st, int row_end, int col_st, int col_end, Matrix<T> &S){
	
	Matrix<T> H(S.dim); H = S;
	for(int i = row_st; i<row_end; i++){
		for(int j = col_st; j<col_end; j++){
			A[i][j] = H(i - row_st, j - col_st);
		}
	}

}
/**Get a square submatrix**/
Matrix<T> GetSubMtx(int row_st, int row_end, int col_st, int col_end){

	Matrix<T> H(row_end - row_st);
	if( (row_end - row_st) == (col_end - col_st) ){	

	for(int i = row_st; i<row_end; i++){
		for(int j = col_st; j<col_end; j++){
			H( A[i][j] , i - row_st, j - col_st);
		}
	}
	}
	else{
		cerr << "Square matricies only use vector for rectangular matricies" << endl;
		exit(1);
	}
	return H;

}
/**Assign a submatrix to a larger one NB remember vector vectors are transposes of the matricies they represent**/
void AssignSubMtx(int row_st, int row_end, int col_st, int col_end, vector< vector<T> > &S){
	
	for(int i = row_st; i<row_end; i++){
		for(int j = col_st; j<col_end; j++){
			A[i][j] = S[j - col_st][i - row_st];
		}
	}

}
/**Get a square submatrix**/
vector< vector<T> > GetSubMtxVec(int row_st, int row_end, int col_st, int col_end){

	vector< vector<T> > H(col_end - col_st);
	

	for(int i = col_st; i<col_end; i++){
		H[i - col_st].resize(row_end - row_st);
		for(int j = row_st; j<row_end; j++){
			H[i - col_st][j - row_st] = A[j][i];
		}
	}
	
	return H;
	
}
/** Assign value x to position i, j **/
inline void operator() (T x, int i, int j){
//if(i >= dim || j >= dim){cout << "Array Index Bound - (" << i << "," << j << ") >= (" << dim << "," << dim << ")" << endl;
//			exit(1);}
			A[i][j] = x;
	} 
/** Assign value x to position i, j **/
inline void operator() (complex<T> x, int i, int j){
//if(i >= dim || j >= dim){cout << "Array Index Bound - (" << i << "," << j << ") >= (" << dim << "," << dim << ")" << endl;
//			exit(1);}
			A[i][j] = real(x);
	} 
/** Return value at position i, j **/
inline T operator() (int i, int j){
//if(i >= dim || j >= dim){cout << "Array Index Bound - (" << i << "," << j << ") >= (" << dim << "," << dim << ")" << endl;
	//		exit(1);}
			return A[i][j];
	}
/** compute b_i A_ij b_j **/
T proj(vector<T> B){
T C = 0;
	for(int i=0;i<dim;i++){
			T sum = 0.0;
			for(int j=0;j<dim;j++){
				sum += A[i][j]*B[j];}
				C +=  B[i]*sum;
			}
	return C; 
	}
};

/**
*************************************************************

Matrix vector products

For some stupid reason if I have M = vector(vector())
   col   row	
M  [i]  [j] 

The functions below transpose M before multiplication

*************************************************************
**/



/** 
	 vector *  matrix
	(N x M) * (M x M) 
										inner index			**/
template <class T> vector< vector<T> > times(vector< vector<T> > Q, Matrix<T> &R, int M){
int N = Q[0].size();

	vector<vector<T > > TMP(M); T sum;
	for(int i = 0;i<N;i++){
	for(int j = 0;j<M;j++){
		TMP[j].resize(N); sum = 0;
	for(int k = 0;k<M;k++){
		sum += Q[k][i]*R(k,j);
	}
	TMP[j][i] = sum;
	}}

return TMP;

}

/** 
	 matrix *  vector
	(M x M) * (M x N) 
										inner index			**/
template <class T> vector< vector<T> > times(Matrix<T> &R, vector< vector<T> > Q, int M){
int N = Q.size();

	vector<vector<T > > TMP(N); T sum;
	for(int i = 0;i<M;i++){
	for(int j = 0;j<N;j++){
		TMP[j].resize(M); sum = 0;
	for(int k = 0;k<M;k++){
		sum += R(i,k)*Q[j][k];
	}
	TMP[j][i] = sum;
	}}

return TMP;

}

/** 
	 vector *  matrix
	(N) * (N x N) 
										inner index			**/
template <class T> vector<T> times(vector<T> Q, Matrix<T> &R){
int N = Q.size();

	vector<T> TMP(N); T sum;
	for(int i = 0;i<N;i++){
		sum = 0;
	for(int k = 0;k<N;k++){
		sum += Q[k]*R(k,i);
	}
	TMP[i] = sum;
	}

return TMP;

}


/** 
	 vector *  vector
	(N x M) * (M) 
										inner index			**/
template <class T> vector<T> times(vector< vector<T> > Q, vector<T> R){
int N = Q[0].size();
int M = R.size();

vector<T> S(N);T sum = 0;
for(int i=0;i<N;i++){
	sum = 0;
	for(int j=0;j<M;j++){
		sum = sum + Q[j][i]*R[j];
	}
	S[i] = sum;
}

return S;

}

/** 
	 vector *  vector
	(N x M) * (M x P) 
										inner index			**/
template <class T> vector<vector<T> > times(vector< vector<T> > Q, vector< vector<T> > R){
int N = Q[0].size();
int M = Q.size();
int P = R.size();

vector< vector<T> > S(P);T sum = 0;
for(int i=0;i<P;i++){
	S[i].resize(N);
	for(int j=0;j<N;j++){
		sum = 0;
	for(int k=0;k<M;k++){
		sum = sum + Q[k][j]*R[i][k];
	}
		S[i][j] = sum;
	}
}

return S;

}


template <class T> void CG_Matrix(Matrix<T> &A, vector<T> source, vector<T> &solution){

	int N = A.dim;
	
  	vector<T> r(N), d(N), q(N), b(N), md(N);

	T delta_new, delta_old, delta_0, alpha, beta;
	int i = 0;
	Matrix<T> AH(N); AH = A.Hermitian();
	//this->D = DDAG_5D;
	//this->dwf_multiply(source, b);
	//this->axpy(solution, 0.0, b, b );
	b = AH*source;
	solution = b;


	r = A*solution;
	r = AH*r;
	//this->D = DDAGD_5D;
	//this->herm_mult(solution, r);

	//delta_new = this->axpy_norm(r, -1.0, r, b );
	r = sub(b,r); delta_new = norm2(r);
	delta_0 = delta_new;
	//this->axpy(d, 0.0, r, r );
	d = r;
	int imax = 10000;
	
	
	while(i < imax && abs(sqrt(delta_new)) > (10e-8) ){
	
	/*this->D = D_5D;
	this->dwf_multiply(d, md);*/
	md = A*d;
	//alpha = this->axpy_norm(md,0,md,md);
	alpha = norm2(md);

	/*this->D = DDAG_5D;
	this->dwf_multiply(md, q);*/
	q = AH*md;
	alpha = delta_new/alpha;
	//this->axpy(solution, alpha, d, solution);
	solution = add(solution, times(alpha,d) );
	//beta = this->axpy_norm(r, -1.0*alpha, q, r);
	r = add(r, times(q,-1.0*alpha) ); beta = norm2(r);

	delta_old = delta_new;
	delta_new = beta;
	beta = delta_new/delta_old;
	//this->axpy(d, beta, d, r);
	d = add(r, times(beta,d) );
	i++;
	}

  QDPIO::cout << "Matrix::Unpreconditioned CG converged in " << i << " iterations" << endl;
  /*this->D = D_5D;
  dwf_multiply(solution,r); */
  r = A*solution;
  double diff = abs(norm( sub(r, source) ) );
  QDPIO::cout << "Matrix::True residual = " << (diff) << endl;


}



/// q -> q Q
template <class T> void times(multi1d<LatticeFermion> &q, Matrix<T> &Q){
	multi1d<LatticeFermion> S( q.size() );
	for(int j=0;j<q.size();j++){
		S[j] = zero;
		for(int k=0;k<q.size();k++){
			S[j] = S[j] + q[k] * toComplex( Q(k,j) ); 
		}
	}
	for(int j=0;j<q.size();j++){q[j] = S[j];}
}



/// q -> q Q
template <class T> void times(multi1d<LatticeFermion> &q, Matrix<T> &Q, int N){
	multi1d<LatticeFermion> S( N );
	for(int j=0;j<N;j++){
		S[j] = zero;
		for(int k=0;k<N;k++){
			S[j] = S[j] + q[k] * toComplex( Q(k,j) ); 
		}
	}
	for(int j=0;j<N;j++){q[j] = S[j];}
}




#endif
