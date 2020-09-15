#ifndef BFM_FRANCIS_H
#define BFM_FRANCIS_H

#include <cstdlib>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <complex>
#include <algorithm>
#include "UTSolve.h"
#include "Householder.h"
#include <util/time_cps.h>

//#define USE_LAPACK
//#define USE_EIGEN

#ifdef USE_LAPACK
#warning "Using LAPACK in Lanczos"
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>
#endif
#ifdef USE_EIGEN
#warning "Using EIGEN package in Lanczos"
#include <Eigen/Dense>
#endif

namespace BFM_Krylov{

  // template <class T> int SymmEigensystem(Matrix<T > &Ain, std::vector<T> &evals, std::vector<std::vector<T> > &evecs, double small);
  // template <class T> int Eigensystem(Matrix<T > &Ain, std::vector<T> &evals, std::vector<std::vector<T> > &evecs, double small);

  /**

     Find the eigenvalues of an upper hessenberg matrix using the Francis QR algorithm.

     H = 		
     x  x  x  x  x  x  x  x  x 
     x  x  x  x  x  x  x  x  x 
     0  x  x  x  x  x  x  x  x 
     0  0  x  x  x  x  x  x  x 
     0  0  0  x  x  x  x  x  x 
     0  0  0  0  x  x  x  x  x 
     0  0  0  0  0  x  x  x  x 
     0  0  0  0  0  0  x  x  x 
     0  0  0  0  0  0  0  x  x 

     Factorization is P T P^H where T is upper triangular (mod cc blocks) and P is orthagonal/unitary. 

  **/
  template <class T> int QReigensystem(Matrix<T> &Hin, std::vector<T> &evals, std::vector<std::vector<T> > &evecs, double small){
    Matrix<T > H(Hin.dim); H = Hin; ///I don't want to modify the input but matricies must be passed by reference
    for(int i=0;i<evecs.size();i++){evals[i] = 0;for(int j=0;j<evecs[i].size();j++){evecs[i][j] = 0;}}

    int N = H.dim;
    int M = N;
    T s,t,x=0,y=0,z=0;
    T u,d;
    T apd,amd,bc;
    std::vector<T > p(N,0);
    T nrm = H.Norm();		///Matrix Norm
    int n, m;
    int e = 0;
    int it = 0;
    int tot_it = 0;
    int l = 0;
    int r = 0;
    Matrix<T > P(N); P.Unity();
    std::vector<int> trows(N,0);

    /// Check if the matrix is really hessenberg, if not abort
    double sth = 0;
    for(int j=0;j<N;j++){for(int i=j+2;i<N;i++){sth = abs(H(i,j));
	if(sth > small){
	  std::cout << "Non hessenberg H = " << sth << " > " << small << std::endl; 
	  exit(1);
	} 
      }}

    /** Check for convergence 

	x  x  x  x  x
	0  x  x  x  x
	0  0  x  x  x
	0  0  x  x  x
	0  0  0  0  x

	for this matrix l = 4
    **/
    do{
      do{

	l = H.Chop_subdiag(nrm,e,small); 

	r = 0;		///May have converged on more than one eval

	///Single eval
	if(l == N-1){
	  evals[e] = H(l,l); 
	  N--; e++; r++; it = 0;
	}
	///Double eval
	if(l == N-2){
	  trows[l+1] = 1;		///Needed for UTSolve
	  apd = H(l,l) + H(l+1,l+1);
	  amd = H(l,l) - H(l+1,l+1);
	  bc =  (T)4.0*H(l+1,l)*H(l,l+1);

	  evals[e] = (T)0.5*( apd + sqrt(amd*amd + bc) );
	  evals[e+1] = (T)0.5*( apd - sqrt(amd*amd + bc) ); 
	  N-=2; e+=2; r++; it = 0;
	}
      }while(r>0);

      if(N ==0){break;}

      std::vector<T > ck(3), v(3);
      for(int m = N-3; m >= l; m--){
	///Starting vector essentially random shift. 
	if(it%10 == 0 && N >= 3 && it > 0){
	  s = (T)1.618033989*( abs( H(N-1,N-2) ) + abs( H(N-2,N-3) ) );
	  t = (T)0.618033989*( abs( H(N-1,N-2) ) + abs( H(N-2,N-3) ) );
	  x = H(m,m)*H(m,m) + H(m,m+1)*H(m+1,m) - s*H(m,m) + t;
	  y = H(m+1,m)*(H(m,m) + H(m+1,m+1) - s);
	  z = H(m+1,m)*H(m+2,m+1);
	}
	///Starting vector implicit Q theorem
	else{
	  s = (H(N-2,N-2) + H(N-1,N-1));
	  t = (H(N-2,N-2)*H(N-1,N-1) - H(N-2,N-1)*H(N-1,N-2));
	  x = H(m,m)*H(m,m) + H(m,m+1)*H(m+1,m) - s*H(m,m) + t;
	  y = H(m+1,m)*(H(m,m) + H(m+1,m+1) - s);
	  z = H(m+1,m)*H(m+2,m+1);
	}

	ck[0] = x; ck[1] = y; ck[2] = z;
	if(m == l) break;

	/** Some stupid thing from numerical recipies, seems to work**/
	u=abs(H(m,m-1))*(abs(y)+abs(z));
	d=abs(x)*(abs(H(m-1,m-1))+abs(H(m,m))+abs(H(m+1,m+1)));
	if ((T)abs(u+d) == (T)abs(d) ){l = m; break;}
	//if (u < small){l = m; break;}
      }

      if(it > 100000){
	std::cout << "QReigensystem: bugger it got stuck after 100000 iterations" << std::endl;
	std::cout << "got " << e << " evals " << l << " " << N << std::endl;
	exit(1);
      }

      normalize(ck);		///Normalization cancels in PHP anyway
      T beta;

      Householder_vector<T >(ck, 0, 2, v, beta);
      Householder_mult<T >(H,v,beta,0,l,l+2,0);
      Householder_mult<T >(H,v,beta,0,l,l+2,1);
      ///Accumulate eigenvector
      Householder_mult<T >(P,v,beta,0,l,l+2,1);

      int sw = 0;			///Are we on the last row?

      for(int k=l;k<N-2;k++){
	x = H(k+1,k); 
	y = H(k+2,k); 
	z = (T)0.0;
	if(k+3 <= N-1){z = H(k+3,k);}
	else{sw = 1; v[2] = (T)0.0;}

	ck[0] = x; ck[1] = y; ck[2] = z;

	normalize(ck);

	Householder_vector<T >(ck, 0, 2-sw, v, beta);
	Householder_mult<T >(H,v, beta,0,k+1,k+3-sw,0);
	Householder_mult<T >(H,v, beta,0,k+1,k+3-sw,1);
	///Accumulate eigenvector
	Householder_mult<T >(P,v, beta,0,k+1,k+3-sw,1);
      }

      it++; 
      tot_it++;
    }while(N > 1);

    N = evals.size();

    ///Annoying - UT solves in reverse order;
    std::vector<T> tmp(N);for(int i=0;i<N;i++){tmp[i] = evals[N-i-1];} evals = tmp; 
    UTeigenvectors(H, trows, evals, evecs);

    for(int i=0;i<evals.size();i++){evecs[i] = P*evecs[i]; normalize(evecs[i]);}

    return tot_it;
  }



  /**
     turn a matrix A =  	
     x  x  x  x  x 
     x  x  x  x  x
     x  x  x  x  x
     x  x  x  x  x
     x  x  x  x  x	
     into
     x  x  x  x  x
     x  x  x  x  x
     0  x  x  x  x
     0  0  x  x  x
     0  0  0  x  x
     with householder rotations
     Slow.
  */

  template <class T> void Hess(Matrix<T > &A, Matrix<T> &Q, int start){
    int N = A.dim;	//Matrix Size

    std::vector<T > p(N,0);
    for(int k=start;k<N-2;k++){
      std::vector<T > ck(N-k-1), v(N-k-1);
      for(int i=k+1;i<N;i++){ck[i-k-1] = A(i,k);}	///kth column

      normalize(ck);		///Normalization cancels in PHP anyway
      T beta;
      Householder_vector<T >(ck, 0, ck.size()-1, v, beta);	///Householder vector

      Householder_mult<T>(A,v,beta,start,k+1,N-1,0);	///A -> PA
      Householder_mult<T >(A,v,beta,start,k+1,N-1,1);  ///PA -> PAP^H 
      ///Accumulate eigenvector
      Householder_mult<T >(Q,v,beta,start,k+1,N-1,1);  ///Q -> QP^H
    }
  }

  ///Tridiagonalize a matrix
  template <class T> void Tri(Matrix<T > &A, Matrix<T> &Q, int start){
    int N = A.dim;	//Matrix Size
    Hess(A,Q,start);
  }

  ///Tridiagonalize a matrix
  template <class T> 
  void ForceTridiagonal(Matrix<T > &A){
    int N = A.dim;	//Matrix Size
    for(int l=0;l<N-2;l++){
      for(int k=l+2;k<N;k++){
	A(0,l,k);
	A(0,k,l);
      }
    }
  }





  /** 
      Find the eigenvalues of an upper Hessenberg matrix using the Wilkinson QR algorithm.

      H = 		
      x  x  0  0  0  0  
      x  x  x  0  0  0  
      0  x  x  x  0  0  
      0  0  x  x  x  0  
      0  0  0  x  x  x  
      0  0  0  0  x  x  

      Factorization is P T P^H where T is upper triangular (mod cc blocks) and P is orthagonal/unitary.  **/

  //Rudy's original version with modifications passed down the ages.....

  template <class T> 
  int basic_Wilkinson(Matrix<T> &Hin, std::vector<T> &evals, std::vector<std::vector<T> > &evecs, double small)
  {
    return basic_Wilkinson(Hin, evals, evecs, small, small);
  }

  template <class T> 
  int basic_Wilkinson(Matrix<T> &Hin, std::vector<T> &evals, std::vector<std::vector<T> > &evecs, double small, double tol)
  {
    Matrix<T> H(Hin.dim); 
    H = Hin; ///I don't want to modify the input but matricies must be passed by reference (CK: I don't think the author really understands why we usually pass by reference!)

    //Scale a matrix by its "norm"
    //double Hnorm = abs( Hin.LargestDiag() ); H =  H*(1.0/Hnorm);
    double Hnorm = abs(Hin.Norm()); 
    H = H * (1.0 / Hnorm);

    for(int i = 0; i < evecs.size(); ++i){
      evals[i] = 0;
      for(int j = 0; j < evecs[i].size(); ++j)
	evecs[i][j] = 0;
    }

    int N = H.dim;
    int M = N;
    T s, t, x = 0, y = 0, z = 0;
    T u, d;
    T apd, amd, bc;
    std::vector<T> p(N, 0);
    T nrm = H.Norm();		///Matrix Norm
    int n, m;
    int e = 0;
    int it = 0;
    int tot_it = 0;
    int l = 0;
    int r = 0;
    Matrix<T> P(N); 
    P.Unity();
    std::vector<int> trows(N, 0);

    /// Check if the matrix is really symm tridiag
    double sth = 0;
    for(int j = 0; j < N; ++j){
	for(int i = j + 2; i < N; ++i){
	  if(abs(H(i, j)) > tol || abs(H(j, i)) > tol){
	    QDPIO::cout << "Non Tridiagonal H(" << i << ","<< j << ") = |" << Real( real( H(j,i) ) ) << "| > " << tol << std::endl; 
	    QDPIO::cout << "Warning tridiagonalize and call again" << std::endl;
	    exit(1);
	  } 
	}
    }
    do{
      do{
	//Jasper
	//Check if the subdiagonal term is small enough (<small)
	//if true then it is converged.
	//check start from H.dim - e - 1
	//How to deal with more than 2 are converged?
	//What if Chop_symm_subdiag return something int the middle?
	//--------------
	l = H.Chop_symm_subdiag(nrm, e, small); 
	r = 0;		///May have converged on more than one eval

	//Jasper 
	//In this case
	// x  x  0  0  0  0  
	// x  x  x  0  0  0  
	// 0  x  x  x  0  0  
	// 0  0  x  x  x  0  
	// 0  0  0  x  x  0  
	// 0  0  0  0  0  x  <- l
	//--------------
	///Single eval
	if(l == N - 1){ 
	  evals[e] = H(l, l); 
	  N--; 
	  e++; 
	  r++; 
	  it = 0; 
	}
	//Jasper
	// x  x  0  0  0  0  
	// x  x  x  0  0  0  
	// 0  x  x  x  0  0  
	// 0  0  x  x  0  0  
	// 0  0  0  0  x  x  <- l
	// 0  0  0  0  x  x  
	//--------------
	///Double eval
	if(l == N - 2){
	  trows[l + 1] = 1;		///Needed for UTSolve
	  apd = H(l, l) + H(l + 1, l + 1);
	  amd = H(l, l) - H(l + 1, l + 1);
	  bc =  (T) 4.0 * H(l + 1, l) * H(l, l + 1);
	  
	  evals[e] = (T) 0.5 * (apd + sqrt(amd * amd + bc));
	  evals[e + 1] = (T) 0.5 * (apd - sqrt(amd * amd + bc)); 
	  N -= 2; 
	  e += 2; 
	  r++; 
	  it = 0;
	}
      }while(r > 0);

      //Jasper
      //Already converged
      //--------------
      if(N == 0)
	break;

      std::vector<T> ck(2), v(2);
      for(int m = N - 3; m >= l; m--){
	///Starting vector essentially random shift. 
	if(it%10 == 0 && N >= 3 && it > 0){
	  t = abs(H(N - 1, N - 2)) + abs(H(N - 2, N - 3));
	  x = H(m, m) - t;
	  z = H(m + 1, m);
	}
	///Starting vector implicit Q theorem
	else{
	  d = (H(N - 2, N - 2) - H(N - 1, N - 1)) * (T) 0.5;
	  t = H(N - 1, N - 1) - H(N - 1, N - 2) * H(N - 1, N - 2) / (d + sign(d) * sqrt(d * d + H(N - 1, N - 2) * H(N - 1, N - 2)));
	  x = H(m, m) - t;
	  z = H(m + 1, m);
	}
	//Jasper
	//why it is here????
	//-----------------------
	if(m == l)
	  break;

	u = abs(H(m, m - 1)) * (abs(y) + abs(z));
	d = abs(x) * (abs(H(m - 1, m - 1)) + abs(H(m, m)) + abs(H(m + 1, m + 1)));
	if ((T)abs(u + d) == (T)abs(d)){
	  l = m; 
	  break;
	}
      }
      //Jasper
      if(it > 1000000){
	QDPIO::cout << "Wilkinson: bugger it got stuck after 100000 iterations" << std::endl;
	QDPIO::cout << "got " << e << " evals " << l << " " << N << std::endl;
	exit(1);
      }

      T s, c;
      Givens_calc<T>(x, z, c, s);
      Givens_mult<T>(H, l, l + 1, c, -s, 0);
      Givens_mult<T>(H, l, l + 1, c,  s, 1);
      Givens_mult<T>(P, l, l + 1, c,  s, 1);

      for(int k = l; k < N - 2; ++k){
	x = H(k + 1,k); 
	z = H(k + 2,k); 

	Givens_calc<T>(x, z, c, s);
	Givens_mult<T>(H, k + 1, k + 2, c, -s, 0);
	Givens_mult<T>(H, k + 1, k + 2, c,  s, 1);
	Givens_mult<T>(P, k + 1, k + 2, c,  s, 1);
      }
      it++; 
      tot_it++;
    }while(N > 1);

    N = evals.size();

    ///Annoying - UT solves in reverse order;
    std::vector<T> tmp(N);
    for(int i = 0; i < N; ++i)
      tmp[i] = evals[N-i-1];

    evals = tmp; 

    UTeigenvectors(H, trows, evals, evecs);

    for(int i = 0; i < evals.size(); ++i){
      evecs[i] = P * evecs[i]; 
      normalize(evecs[i]); 
      evals[i] = evals[i] * Hnorm;
    }

    return tot_it;
  }

  //An implementation using the Eigen package courtesy of Luchang
#ifdef USE_EIGEN
  template <class T>
  int eigen_Wilkinson(Matrix<T> &Ain, std::vector<T> &evals, std::vector<std::vector<T> > &evecs, double small){
    using namespace Eigen;
    const int size = Ain.dim;
    evals.resize(size);
    evecs.resize(size);
    for (int i = 0; i < size; i++) {
      evecs[i].resize(size);
    }
    MatrixXd A(size, size);
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
	A(i, j) = 0;
      }
    }
    for (int i = 0; i < size; i++) {
      for (int j = std::max(0, i-1); j <= std::min(size-1, i+1); j++) {
	A(i, j) = Ain(i, j);
      }
    }
    SelfAdjointEigenSolver<MatrixXd> eigensolver(A);
    for (int i = 0; i < size; i++) {
      evals[i] = eigensolver.eigenvalues()(i);
    }
    for (int i = 0; i < size; i++) {
      std::vector<T>& vec = evecs[i];
      for (int j = 0; j < size; j++) {
	vec[j] = eigensolver.eigenvectors()(j, i);
      }
    }
    return 0;
  }
#endif

 //An implementation using the LAPACK package courtesy of Luchang
#ifdef USE_LAPACK
  template <class T>
  int lapackcj_Wilkinson(Matrix<T> &AH, vector<T> &tevals, vector<vector<T> > &tevecs, double small) {
    double time = -cps::dclock();
    const int size = AH.dim;
    tevals.resize(size);
    tevecs.resize(size);
    int NN = tevals.size();
    for(int i=0;i<NN;i++) tevecs[i].resize(NN);
    double evals_tmp[NN];
    double evec_tmp[NN][NN];
    memset(evec_tmp[0],0,sizeof(double)*NN*NN);
    double AA[NN][NN];
    //    double ZZ[NN][NN];
    double DD[NN];
    double EE[NN];
    memset(AA, 0, sizeof(double)*NN*NN);
    //    memset(ZZ, 0, sizeof(double)*NN*NN);
    for (int i = 0; i< NN; i++)
      for (int j = i - 1; j <= i + 1; j++)
	if ( j < NN && j >= 0 ) {
	  AA[i][j] = AH(i,j);
	  if (i==j) DD[i] = AA[i][j];
	  if (i==j) evals_tmp[i] = AA[i][j];
	  if (j==(i-1)) EE[j] = AA[i][j];
	  //        if (i<20 && j<20)
	  QDPIO:: cout << "AA["<<i<<"]["<<j<<"]="<<AA[i][j]<<endl;
	}
    int evals_found;
    int lwork = ( (18*NN) > (1+4*NN+NN*NN)? (18*NN):(1+4*NN+NN*NN)) ;
    int liwork =  3+NN*10 ;
    int iwork[liwork];
    double work[lwork];
    int isuppz[2*NN];
    char jobz = 'V'; // calculate evals & evecs
    char range = 'I'; // calculate all evals
    //    char range = 'A'; // calculate all evals
    char uplo = 'U'; // refer to upper half of original matrix
    char compz = 'I'; // Compute eigenvectors of tridiagonal matrix
    int ifail[NN];
    int info;
    // dsyevx(&jobz, &range, &uplo, NN,
    //     AA, NN,
    //     0, 0, 0, 0, // these four are ignored if second parameteris 'A'
    //     0, // tolerance
    //     &evals_found, evals_tmp, evec_tmp, NN,
    //     work, lwork, iwork,
    //     ifail, info);
    int total = QMP_get_number_of_nodes();
    int node = QMP_get_node_number();
    int interval = (NN/total)+1;
    double vl = 0.0, vu = 0.0;
    int il = interval*node+1 , iu = interval*(node+1);
    if (iu > NN)  iu=NN;
    double tol = 0.0;
    if (0)
      {
	memset(evals_tmp,0,sizeof(double)*NN);
	if ( il <= NN){
	  printf("total=%d node=%d il=%d iu=%d\n",total,node,il,iu);
	  LAPACK_dsyevx(&jobz, &range, &uplo, &NN,
			(double*)AA, &NN,
			&vl, &vu, &il, &iu, // these four are ignored if second parameteris 'A'
			&tol, // tolerance
			&evals_found, evals_tmp, (double*)evec_tmp, &NN,
			work, &lwork, iwork,
			ifail, &info);
	  for (int i = iu-1; i>= il-1; i--){
	    printf("node=%d evals_found=%d evals_tmp[%d] = %g\n",node,evals_found, i - (il-1),evals_tmp[i - (il-1)]);
	    evals_tmp[i] = evals_tmp[i - (il-1)];
	    if (il>1) evals_tmp[i-(il-1)]=0.;
	    for (int j = 0; j< NN; j++){
	      evec_tmp[i][j] = evec_tmp[i - (il-1)][j];
	      if (il>1) evec_tmp[i-(il-1)][j]=0.;
	    }
	  }
	}
	{
	  QMP_sum_double_array(evals_tmp,NN);
	  QMP_sum_double_array((double *)evec_tmp,NN*NN);
	}
      } else
      if (1) {
	memset(evals_tmp,0,sizeof(double)*NN);
	if ( il <= NN){
	  printf("total=%d node=%d il=%d iu=%d\n",total,node,il,iu);
	  LAPACK_dstegr(&jobz, &range, &NN,
			(double*)DD, (double*)EE,
			&vl, &vu, &il, &iu, // these four are ignored if second parameteris 'A'
			&tol, // tolerance
			&evals_found, evals_tmp, (double*)evec_tmp, &NN,
			isuppz,
			work, &lwork, iwork, &liwork,
			&info);
	  for (int i = iu-1; i>= il-1; i--){
	    printf("node=%d evals_found=%d evals_tmp[%d] = %g\n",node,evals_found, i - (il-1),evals_tmp[i - (il-1)]);
	    evals_tmp[i] = evals_tmp[i - (il-1)];
	    if (il>1) evals_tmp[i-(il-1)]=0.;
	    for (int j = 0; j< NN; j++){
	      evec_tmp[i][j] = evec_tmp[i - (il-1)][j];
	      if (il>1) evec_tmp[i-(il-1)][j]=0.;
	    }
	  }
	}
	{
	  QMP_sum_double_array(evals_tmp,NN);
	  QMP_sum_double_array((double *)evec_tmp,NN*NN);
	}
      } else
	if(0) {
	  LAPACK_dsteqr(&compz, &NN,
			(double*)evals_tmp, (double*)EE,
			(double*)evec_tmp, &NN,
			work, &info);
	} else
	  {
	    LAPACK_dstedc(&compz, &NN,
			  (double*)evals_tmp, (double*)EE,
			  (double*)evec_tmp, &NN,
			  work, &lwork,
			  iwork, &liwork,
			  &info);
	  }
    for(int i=0;i<NN;i++){
      for(int j=0;j<NN;j++)
	tevecs[i][j]=evec_tmp[i][j];
      tevals[i]=evals_tmp[i];
    }

    cps::print_time("lapackcj_Wilkinson","run time",time+cps::dclock());
  }
#endif



  template <class T>
  int Wilkinson(Matrix<T> &Ain, std::vector<T> &evals, std::vector<std::vector<T> > &evecs, double small){
#ifdef USE_LAPACK
    return lapackcj_Wilkinson(Ain, evals, evecs, small);
#elif USE_EIGEN
    return eigen_Wilkinson(Ain, evals, evecs, small);
#else
    return basic_Wilkinson(Ain, evals, evecs, small);
#endif
  }





  ///Solve a symmetric eigensystem, not necessarily in tridiagonal form
  //Rudy's original implementation
  template <class T> int basic_SymmEigensystem(Matrix<T > &Ain, std::vector<T> &evals, std::vector<std::vector<T> > &evecs, double small){
    int N = Ain.dim;
    Matrix<T > A(N); A = Ain; 
    Matrix<T > Q(N);Q.Unity();
    Tri(A,Q,0);
    int it = basic_Wilkinson<T>(A, evals, evecs, small);
    for(int k=0;k<N;k++){evecs[k] = Q*evecs[k];}
    return it;
  }

  //I guess Luchang wrote this....
#ifdef USE_EIGEN
  template <class T>
  int eigen_SymmEigensystem(Matrix<T> &Ain, std::vector<T> &evals, std::vector<std::vector<T> > &evecs, double small){
    using namespace Eigen;
    const int size = Ain.dim;
    evals.resize(size);
    evecs.resize(size);
    for (int i = 0; i < size; i++) {
      evecs[i].resize(size);
    }
    MatrixXd A(size, size);
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
	A(i, j) = Ain(i, j);
      }
    }
    SelfAdjointEigenSolver<MatrixXd> eigensolver(A);
    for (int i = 0; i < size; i++) {
      evals[i] = eigensolver.eigenvalues()(i);
    }
    for (int i = 0; i < size; i++) {
      std::vector<T>& vec = evecs[i];
      for (int j = 0; j < size; j++) {
	vec[j] = eigensolver.eigenvectors()(j, i);
      }
    }
    return 0;
  }
#endif

  template <class T>
  int SymmEigensystem(Matrix<T> &Ain, std::vector<T> &evals, std::vector<std::vector<T> > &evecs, double small){
#ifdef USE_EIGEN
    return eigen_SymmEigensystem(Ain, evals, evecs, small);
#else
    return basic_SymmEigensystem(Ain, evals, evecs, small);
#endif
  }


  ///Solve a general eigensystem, not necessarily in tridiagonal form
  template <class T> int Eigensystem(Matrix<T > &Ain, std::vector<T> &evals, std::vector<std::vector<T> > &evecs, double small){

    int N = Ain.dim;
    Matrix<T > A(N); A = Ain; 
    Matrix<T > Q(N);Q.Unity();
    Hess(A,Q,0);
    int it = QReigensystem<T>(A, evals, evecs, small);
    for(int k=0;k<N;k++){evecs[k] = Q*evecs[k];}

    return it;
  }


}

#endif



