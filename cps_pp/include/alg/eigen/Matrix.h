#ifndef BFM_MATRIX_H
#define BFM_MATRIX_H

#include <qdp.h>
#include <qdp_util.h>

//#include <chroma.h>
//#include <actions/ferm/invert/syssolver_linop_cg_array.h>

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
namespace BFM_Krylov{

  typedef Fermion_t bfm_fermion [2];

  /** 
      Do complex stuff to floats/doubles 
  **/
  inline float conj(float z){return z;}
  static float real(float z){return z;}
  static float imag(float z){return 0;}

  inline double conj(double z){return z;}
  static double real(double z){return z;}
  static double imag(double z){return 0;}

  /** Chroma complex to C++ complex(double) **/
  template <class T> std::complex<T> toComplex(QDP::Complex z){return std::complex<T>(toDouble(real(z) ) , toDouble(imag(z) ) ); }
  template <class T> T toComplex(Real z){return (T)toDouble(z); }

  /** C++ complex(double) to Chroma complex **/
  template <class T> QDP::Complex toComplex(T z){return cmplx(Real(real(z) ) , Real(imag(z) ) ); }


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
  template <class T> std::vector<T> conj(const std::vector<T> &p){
    std::vector<T> q(p.size());
    for(int i=0;i<p.size();i++){q[i] = conj(p[i]);}
    return q;
  }

  /** Norm of a vector**/
  template <class T> T norm(const std::vector<T> &p){
    T sum = 0;
    for(int i=0;i<p.size();i++){sum = sum + p[i]*conj(p[i]);}
    return abs(sqrt(sum));
  }

  /** Norm squared of a vector **/
  template <class T> T norm2(const std::vector<T> &p){
    T sum = 0;
    for(int i=0;i<p.size();i++){sum = sum + p[i]*conj(p[i]);}
    return abs((sum));
  }

  /** Sum elements of a vector **/
  template <class T> T trace(const std::vector<T> &p){
    T sum = 0;
    for(int i=0;i<p.size();i++){sum = sum + p[i];}
    return sum;
  }

  /** Fill a vector with constant c **/
  template <class T> void fill(std::vector<T> &p, const T c){
    for(int i=0;i<p.size();i++){p[i] = c;}
  }

  /** Normalize a vector **/
  template <class T> void normalize(std::vector<T> &p){
    T m = norm(p);
    if( abs(m) > 0.0) for(int i=0;i<p.size();i++){p[i] /= m;}
  }

  /** print components of a vector in mathematica readable (list) form **/
  template <class T> void print(const std::vector<T> &p, const int w){
    using namespace std;
    cout << fixed << setprecision(w) << setw(w) << "{ ";
    for(int i=0;i<p.size();i++){
      i==p.size()-1 ? cout << real(p[i]) << " +I*" << imag(p[i]) << " }; ":cout << real(p[i]) << " +I*" << imag(p[i]) << " , ";
    }
    cout << endl;
  }

  template <class T> void print(const std::vector<T> &p){
    print(p,3);
  }

  /** print components of a vector of vectors in mathematica readable (matrix) form **/
  template <class T> void print(const std::vector< std::vector<T> > &Q){
    using namespace std;
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
  template <class T> void printSimpleReal(const std::vector< std::vector<T> > &Q){
    using namespace std;
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
  template <class T, class U> std::vector<T> times(std::vector<T> p, const U s){
    for(int i=0;i<p.size();i++){p[i] *= s;}
    return p;
  }
  template <class T, class U> std::vector<T> times(const U s, std::vector<T> p){
    for(int i=0;i<p.size();i++){p[i] *= s;}
    return p;
  }

  /** inner product of a and b = conj(a) . b **/
  template <class T> T inner(const std::vector<T> &a, const std::vector<T> &b){
    T m = 0.;
    for(int i=0;i<a.size();i++){m += std::conj(a[i])*b[i];}
    return m;
  }

  /** sum of a and b = a + b **/
  template <class T> std::vector<T> add(const std::vector<T> &a, const std::vector<T> &b){
    std::vector<T> m(a.size());
    for(int i=0;i<a.size();i++){m[i] = a[i] + b[i];}
    return m;
  }

  /** sum of a and b = a - b **/
  template <class T> std::vector<T> sub(const std::vector<T> &a, const std::vector<T> &b){
    std::vector<T> m(a.size());
    for(int i=0;i<a.size();i++){m[i] = a[i] - b[i];}
    return m;
  }


  /** 
*********************************
*				*
*	Square matrices	*
* 				*
*********************************
**/

  //CK faster version
  template <class T> class Matrix {
    T* v;

  public:
    int dim;		/// dimension

    void alloc_matrix(){
      v = (T*)malloc(dim*dim*sizeof(T));
    }      
    void free_matrix(){
      if(v != NULL) free(v);
      v = NULL;
    }

    /**Default**/
    Matrix(void): v(NULL){
      dim=0;
    }

    /** Construct new matrix unspecified entries**/
    Matrix(const int d){
      dim = d; 
      alloc_matrix();
    }

    Matrix(const Matrix<T> &r): dim(r.dim){
      alloc_matrix();
      memcpy(v, r.v, dim*dim*sizeof(T));
    }


    /** Construct new matrix entries all = val**/
    Matrix(const int d, const T val){ 
      dim = d; 
      alloc_matrix();
      if(val == 0.0) memset(v, 0, dim*dim*sizeof(T));
      for(int i=0;i<dim*dim;i++) v[i] = val;
    }

    /** Destroy a matrix **/
    ~Matrix(){
      free_matrix();
    }

    /** Resize a matrix destroying everything in it **/
    void resize(const int N){
      if(dim == N) return;
      if(v != NULL) free(v);
      dim = N; 
      alloc_matrix();
    }

    /** initialize a matrix with a vector(vector(
	my convection (for reasons that are not obvious here but made writing some loops easier)
    **/
    Matrix(const std::vector< std::vector<T> > &B){ 
      dim = B.size();
      alloc_matrix();
      for(int i = 0; i<dim; i++){
	for(int j=0;j<dim;j++){v[i+dim*j] = B[j][i];}
      }
    }

    /** Assign value x to position i, j **/
    inline void operator() (const T x, const int i, const int j){
      v[i+dim*j] = x;
    } 
    /** Assign real value of x to position i, j **/
    inline void operator() (const std::complex<T> x, const int i, const int j){
      v[i+dim*j] = real(x);
    } 
    /** Return value at position i, j **/
    inline T operator() (const int i, const int j) const{
      return v[i+dim*j];
    }
    inline T& operator() (const int i, const int j){
      return v[i+dim*j];
    }



    /** Fill matricies with the value x **/
    void Fill(const T x){
      for(int i=0;i<dim*dim;i++) v[i] = x;
    }

    /** print matrix in Mathematica Readable form : default precision **/
    void Out() const{
      OutSimpleReal(6);
    }

    void OutSimpleReal(const int precis) const{
      using namespace std;
      for(int i=0;i<dim;i++){
	for(int j=0;j<dim;j++){
	  cout << fixed << setw(precis) << setprecision(precis) << real(v[i+dim*j]) << " , ";
	}
	cout << endl;			
      } 
      cout << endl;
    }
	
    /** print matrix in Mathematica Readable form : user precision **/
    void Out(const int precis) const{
      using namespace std;
      cout << "{";
      for(int i=0;i<dim;i++){
	cout << "{";
	for(int j=0;j<dim;j++){
	  j == (dim-1) ? 
	    cout << fixed << setprecision(precis)<< setw(precis) <<  real(v[i+dim*j]) << " + I*" << imag(v[i+dim*j]) << " " : 
	    cout << fixed << setprecision(precis) << setw(precis) << real(v[i+dim*j]) << " + I*" << imag(v[i+dim*j]) << "\t,"; } 
	i == (dim-1) ? cout << "}};" : cout << "},";
	cout << endl;
      } 
    }

    /** print matrix in Mathematica Readable form : scientific **/
    void OutSci(const int precis) const{
      using namespace std;
      cout << "{";
      for(int i=0;i<dim;i++){
	cout << "{";
	for(int j=0;j<dim;j++){
	  j == (dim-1) ? 
	    cout << fixed << scientific << setprecision(precis)<< setw(precis) <<  real(v[i+dim*j]) << " + I*" << imag(v[i+dim*j]) << " " : 
	    cout << fixed << scientific << setprecision(precis) << setw(precis) << real(v[i+dim*j]) << " + I*" << imag(v[i+dim*j]) << "\t,"; } 
	i == (dim-1) ? cout << "}};" : cout << "},";
	cout << endl;
      } 
    }

    /** print matrix in c form **/
    void UseOut() const{
      using namespace std;
      for(int i=0;i<dim;i++){
	for(int j=0;j<dim;j++){
	  cout << fixed << scientific <<  setprecision(14) << setw(14) << "H(complex<double> ( " <<real(v[i+dim*j]) << "," << imag(v[i+dim*j]) << ") , " << i << " , " << j << " );\n"; } 
	cout << endl;
      } 
    }

    /** print matrix in c form **/
    void UseOutR() const{
      using namespace std;
      for(int i=0;i<dim;i++){
	for(int j=0;j<dim;j++){
	  cout << fixed << scientific <<  setprecision(14) << setw(14) << "H(" <<real(v[i+dim*j]) <<  " , " << i << " , " << j << " );\n"; } 
	cout << endl;
      } 
    }

    /** Transpose of a matrix **/
    Matrix<T> Transpose(){
      Matrix C(dim);
      for(int i=0;i<dim;i++){
	for(int j=0;j<dim;j++){C.v[i+dim*j] = v[j+dim*i];} 
      } 
      return C;
    }

    /** Set Matrix to unit matrix **/
    void Unity(){
      memset(v,0, dim*dim*sizeof(T));
      for(int i=0;i<dim;i++) v[i+dim*i] = 1.;
    }


    /** Add C * I to matrix **/
    void PlusUnit(const T c){
      for(int i=0;i<dim;i++){v[i+dim*i] += c;} 
    }

    /** return the Hermitian conjugate of matrix **/
    Matrix<T> Hermitian() const{
      Matrix C(dim);
      for(int i=0;i<dim;i++){
	for(int j=0;j<dim;j++)
	  C.v[i+dim*j] = conj(v[j+dim*i]);
      } 
      return C;
    }

    /** return diagonal entries as a vector **/
    std::vector<T> diag() const{
      std::vector<T> d(dim);
      for(int i=0;i<dim;i++){d[i] = v[i+dim*i];}
      return d;
    }

    /** Left multiply by a vector **/
    std::vector<T> LeftMult(const std::vector<T> &B) const{
      std::vector<T> C(dim);
      for(int j=0;j<dim;j++){
	T sum = 0.0;
	for(int i=0;i<dim;i++){
	  sum += B[i] * v[i+dim*j];}
	C[j] =  sum;
      }
      return C; 
    }

    /** return 1/diagonal entries as a vector **/
    std::vector<T> inv_diag() const{
      std::vector<T> d(dim);
      for(int i=0;i<dim;i++){d[i] = 1.0/v[i+dim*i];}
      return d;
    }
    /** sizeof **/
    inline int size() const{ 
      return dim;
    } 
    /** Matrix Addition **/
    inline Matrix<T> operator+ (const Matrix<T> &B) const{
      Matrix C(dim);
      for(int i=0;i<dim;i++){
	for(int j=0;j<dim;j++){C(i,j) = (*this)(i,j) +  B(i,j);} 
      } 
      return C;
    } 

    /** Matrix Subtraction **/
    inline Matrix<T> operator- (const Matrix<T> &B) const{
      Matrix C(dim);
      for(int i=0;i<dim;i++){
	for(int j=0;j<dim;j++){C(i,j) = (*this)(i,j) -  B(i,j);} 
      } 
      return C;
    } 

    /** Matrix scalar multiplication **/
    inline Matrix<T> operator* (const T c) const{
      Matrix C(dim);
      for(int i=0;i<dim*dim;i++) C.v[i] = v[i]*c;
      return C;
    } 

    /** Matrix Matrix multiplication **/
    inline Matrix<T> operator* (const Matrix<T> &B) const{
      Matrix C(dim);
      for(int i=0;i<dim;i++){
	for(int j=0;j<dim;j++){
	  C(i,j) = 0.0;
	  for(int k=0;k<dim;k++){ C(i,j) += (*this)(i,k)*B(k,j);}
	}
      }
      return C; 
    } 

    /** Matrix vector multiplication **/
    inline std::vector<T> operator* (const std::vector<T> &B){
      std::vector<T> C(dim);
      for(int i=0;i<dim;i++){
	C[i] = 0.0;
	for(int j=0;j<dim;j++){
	  C[i] += (*this)(i,j)*B[j];}
      }
      return C; 
    } 

    /** Matrix equivalence **/
    inline Matrix<T>& operator=(const Matrix<T>& B){
      if(dim != B.dim){
	free_matrix();
	dim = B.dim;
	alloc_matrix();
      }
      memcpy(v, B.v, dim*dim*sizeof(T));
      return (*this); 
    }

    /** Some version of Matrix norm **/
    inline T Norm() const{ 
      T norm = 0;
      for(int i=0;i<dim*dim;i++) norm += abs(v[i]);
      return norm;
    }

    /** Some version of Matrix norm **/
    inline T LargestDiag() const{

      T ld = abs((*this)(0,0));
      for(int i=1;i<dim;i++){
  
	T cf = abs((*this)(i,i));
	if(abs(cf) > abs(ld) ){ld = cf;}
  
      }
      return ld;
    }

    /** Look for entries on the leading subdiagonal that are smaller than 'small' **/
    template <class U> 
    int Chop_subdiag(const T norm, const int offset, const U small){
      for(int l = dim - 1 - offset; l >= 1; l--){             		
	if((U)abs( (*this)(l,l - 1)) < (U)small){
	  (*this)(l,l-1)=(U)0.0;
	  return l;
	}
      }
      return 0;
    }

    /** Look for entries on the leading subdiagonal that are smaller than 'small' **/
    template <class U> 
    int Chop_symm_subdiag(const T norm, const int offset, const U small){
      for(int l = dim - 1 - offset; l >= 1; l--){             		
	if((U)abs( (*this)(l,l - 1) ) < (U)small){
	  (*this)(l,l - 1) = (U)0.0;
	  (*this)(l - 1,l) = (U)0.0;
	  return l;
	}
      }
      return 0;
    }

    /**Assign a submatrix to a larger one**/
    void AssignSubMtx(const int row_st, const int row_end, const int col_st, const int col_end, const Matrix<T> &S){
      for(int i = row_st; i<row_end; i++){
	for(int j = col_st; j<col_end; j++){
	  v[i+dim*j] = S(i - row_st, j - col_st);
	}
      }

    }
    /**Get a square submatrix**/
    Matrix<T> GetSubMtx(const int row_st, const int row_end, const int col_st, const int col_end) const{
      Matrix<T> H(row_end - row_st);
      if( (row_end - row_st) == (col_end - col_st) ){	

	for(int i = row_st; i<row_end; i++){
	  for(int j = col_st; j<col_end; j++){
	    H( v[i+dim*j] , i - row_st, j - col_st);
	  }
	}
      }
      else{
	std::cerr << "Square matrices only use vector for rectangular matrices" << std::endl;
	exit(1);
      }
      return H;
    }
    /**Assign a submatrix to a larger one NB remember vector vectors are transposes of the matricies they represent**/
    void AssignSubMtx(const int row_st, const int row_end, const int col_st, const int col_end, const std::vector< std::vector<T> > &S){	
      for(int i = row_st; i<row_end; i++){
	for(int j = col_st; j<col_end; j++){
	  v[i+dim*j] = S[j - col_st][i - row_st];
	}
      }
    }
    /**Get a square submatrix**/
    std::vector< std::vector<T> > GetSubMtxVec(const int row_st, const int row_end, const int col_st, const int col_end) const{
      std::vector< std::vector<T> > H(col_end - col_st);
	
      for(int i = col_st; i<col_end; i++){
	H[i - col_st].resize(row_end - row_st);
	for(int j = row_st; j<row_end; j++){
	  H[i - col_st][j - row_st] = v[j+dim*i];
	}
      }	
      return H;	
    }

    /** compute b_i A_ij b_j **/
    T proj(const std::vector<T> &B){
      T C = 0;
      for(int i=0;i<dim;i++){
	T sum = 0.0;
	for(int j=0;j<dim;j++){
	  sum += v[i+dim*j]*B[j];}
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
  template <class T> std::vector< std::vector<T> > times(const std::vector< std::vector<T> > &Q, const Matrix<T> &R, const int M){
    int N = Q[0].size();

    std::vector<std::vector<T > > TMP(M); T sum;
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
  template <class T> std::vector< std::vector<T> > times(const Matrix<T> &R, const std::vector< std::vector<T> > &Q, const int M){
    int N = Q.size();

    std::vector<std::vector<T > > TMP(N); T sum;
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
  template <class T> std::vector<T> times(const std::vector<T> &Q, const Matrix<T> &R){
    int N = Q.size();

    std::vector<T> TMP(N); T sum;
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
  template <class T> std::vector<T> times(const std::vector< std::vector<T> > &Q, const std::vector<T> &R){
    int N = Q[0].size();
    int M = R.size();

    std::vector<T> S(N);T sum = 0;
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
  template <class T> std::vector<std::vector<T> > times(const std::vector< std::vector<T> > &Q, const std::vector< std::vector<T> > &R){
    int N = Q[0].size();
    int M = Q.size();
    int P = R.size();

    std::vector< std::vector<T> > S(P);T sum = 0;
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

  template <class T> void Gram(std::vector<std::vector<T> > &q){

    int M = q.size();
    for(int j=0;j<M;j++){
      for(int i=0;i<j;i++){
	q[j] = sub( q[j]   ,  times(  ( inner(q[i],q[j])/norm2(q[i]) ) , q[i])  );
      }
      q[j] = times( q[j] , ( 1.0/sqrt( norm2(q[j]) ) ) );
    }

  }



  template <class T> void CG_Matrix(const Matrix<T> &A, const std::vector<T> &source, std::vector<T> &solution){

    int N = A.dim;
	
    std::vector<T> r(N), d(N), q(N), b(N), md(N);

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

    QDPIO::cout << "Matrix::Unpreconditioned CG converged in " << i << " iterations" << std::endl;
    /*this->D = D_5D;
      dwf_multiply(solution,r); */
    r = A*solution;
    double diff = abs(norm( sub(r, source) ) );
    QDPIO::cout << "Matrix::True residual = " << (diff) << std::endl;


  }


  static void Gram(LatticeFermion &r, multi1d<LatticeFermion> &q, int N){
    for(int i=0;i<N;i++){
      r = r   -  ( innerProduct(q[i],r)/innerProduct(q[i],q[i]) )*q[i];
    }
    r = r*( 1.0/sqrt( norm2(r) ) );
  }

  static void Gram(LatticeFermion &r, multi1d<LatticeFermion> &q){
    int N = q.size();
    for(int i=0;i<N;i++){
      r = r   -  ( innerProduct(q[i],r)/innerProduct(q[i],q[i]) )*q[i];
    }
    r = r*( 1.0/sqrt( norm2(r) ) );
  }

  static void Gram(multi1d<LatticeFermion> &q){
    int M = q.size();
    for(int j=0;j<M;j++){
      for(int i=0;i<j;i++){
	q[j] = q[j]   -  ( innerProduct(q[i],q[j])/innerProduct(q[i],q[i]) )*q[i];
      }
      q[j] = q[j]*( 1.0/sqrt( norm2(q[j]) ) );
    }
  }

  static void Gram(std::vector<LatticeFermion> &q){
    int M = q.size();
    for(int j=0;j<M;j++){
      for(int i=0;i<j;i++){
	q[j] = q[j]   -  ( innerProduct(q[i],q[j])/innerProduct(q[i],q[i]) )*q[i];
      }
      q[j] = q[j]*( 1.0/sqrt( norm2(q[j]) ) );
    }
  }


  static void Gram(multi1d<LatticeFermion> &q, const int M){
    for(int j=0;j<M;j++){
      for(int i=0;i<j;i++){
	q[j] = q[j]   -  ( innerProduct(q[i],q[j])/innerProduct(q[i],q[i]) )*q[i];
      }
      q[j] = q[j]*( 1.0/sqrt( norm2(q[j]) ) );
    }
  }

  static void Gram(std::vector<LatticeFermion> &q, const int M){
    for(int j=0;j<M;j++){
      for(int i=0;i<j;i++){
	q[j] = q[j]   -  ( innerProduct(q[i],q[j])/innerProduct(q[i],q[i]) )*q[i];
      }
      q[j] = q[j]*( 1.0/sqrt( norm2(q[j]) ) );
    } 
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
  template <class T> void times(multi1d<LatticeFermion> &q, const Matrix<T> &Q, const int N){
    multi1d<LatticeFermion> S( N );
    for(int j=0;j<N;j++){
      S[j] = zero;
      for(int k=0;k<N;k++){
	S[j] = S[j] + q[k] * toComplex( Q(k,j) ); 
      }
    }
    for(int j=0;j<N;j++){q[j] = S[j];}
  }


}

#endif
