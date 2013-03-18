#ifndef FRANCIS_H
#define FRANCIS_H

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

using namespace std;

template <class T> int SymmEigensystem(Matrix<T > &Ain, vector<T> &evals, vector<vector<T> > &evecs, double small);
template <class T> int Eigensystem(Matrix<T > &Ain, vector<T> &evals, vector<vector<T> > &evecs, double small);

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
template <class T> int QReigensystem(Matrix<T> &Hin, vector<T> &evals, vector<vector<T> > &evecs, double small){


  Matrix<T > H(Hin.dim); H = Hin; ///I don't want to modify the input but matricies must be passed by reference
  for(int i=0;i<evecs.size();i++){evals[i] = 0;for(int j=0;j<evecs[i].size();j++){evecs[i][j] = 0;}}

  int N = H.dim;
  int M = N;
  T s,t,x=0,y=0,z=0;
  T u,d;
  T apd,amd,bc;
  vector<T > p(N,0);
  T nrm = H.Norm();		///Matrix Norm
  int n, m;
  int e = 0;
  int it = 0;
  int tot_it = 0;
  int l = 0;
  int r = 0;
  Matrix<T > P(N); P.Unity();
  vector<int> trows(N,0);

  /// Check if the matrix is really hessenberg, if not abort
  double sth = 0;
  for(int j=0;j<N;j++){for(int i=j+2;i<N;i++){sth = abs(H(i,j));
    if(sth > small){
      cout << "Non hessenberg H = " << sth << " > " << small << endl; 
      exit(1);
    } 
  }}



  do{

    //cerr << "Francis QR Step N = " << N << endl;


    /** Check for convergence 

      x  x  x  x  x
      0  x  x  x  x
      0  0  x  x  x
      0  0  x  x  x
      0  0  0  0  x

      for this matrix l = 4
     **/
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


    vector<T > ck(3), v(3);
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
      cout << "bugger it got stuck after 100000 iterations" << endl;
      cout << "got " << e << " evals " << l << " " << N << endl;
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
  vector<T> tmp(N);for(int i=0;i<N;i++){tmp[i] = evals[N-i-1];} evals = tmp; 
  UTeigenvectors(H, trows, evals, evecs);

  for(int i=0;i<evals.size();i++){evecs[i] = P*evecs[i]; normalize(evecs[i]);}

  return tot_it;
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
template <class T> 
int Wilkinson(Matrix<T> &Hin, vector<T> &evals, vector<vector<T> > &evecs, double small)
{
  return Wilkinson(Hin, evals, evecs, small, small);
}

template <class T> 
int Wilkinson(Matrix<T> &Hin, vector<T> &evals, vector<vector<T> > &evecs, double small, double tol)
{
  Matrix<T> H(Hin.dim); 
  H = Hin; ///I don't want to modify the input but matricies must be passed by reference

  //Scale a matrix by its "norm"
  //double Hnorm = abs( Hin.LargestDiag() ); H =  H*(1.0/Hnorm);
  double Hnorm = abs(Hin.Norm()); 
  H = H * (1.0 / Hnorm);

  for(int i = 0; i < evecs.size(); ++i)
  {
    evals[i] = 0;
    for(int j = 0; j < evecs[i].size(); ++j)
      evecs[i][j] = 0;
  }

  int N = H.dim;
  int M = N;
  T s, t, x = 0, y = 0, z = 0;
  T u, d;
  T apd, amd, bc;
  vector<T> p(N, 0);
  T nrm = H.Norm();		///Matrix Norm
  int n, m;
  int e = 0;
  int it = 0;
  int tot_it = 0;
  int l = 0;
  int r = 0;
  Matrix<T> P(N); 
  P.Unity();
  vector<int> trows(N, 0);

  /// Check if the matrix is really symm tridiag
  double sth = 0;
  for(int j = 0; j < N; ++j)
  {
    for(int i = j + 2; i < N; ++i)
    {
      if(abs(H(i, j)) > tol || abs(H(j, i)) > tol)
      {
        QDPIO::cout << "Non Tridiagonal H(" << i << ","<< j << ") = |" << Real( real( H(j,i) ) ) << "| > " << tol << endl; 
        QDPIO::cout << "Warning tridiagonalize and call again" << endl;
        exit(1);
        //return;
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
      if(l == N - 1)
      { 
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
      if(l == N - 2)
      {
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

    vector<T> ck(2), v(2);
    for(int m = N - 3; m >= l; m--)
    {
      ///Starting vector essentially random shift. 
      if(it%10 == 0 && N >= 3 && it > 0)
      {
        t = abs(H(N - 1, N - 2)) + abs(H(N - 2, N - 3));
        x = H(m, m) - t;
        z = H(m + 1, m);
      }
      ///Starting vector implicit Q theorem
      else
      {
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
      if ((T)abs(u + d) == (T)abs(d))
      {
        l = m; 
        break;
      }
    }
    //Jasper
    if(it > 1000000)
    {
      QDPIO::cout << "bugger it got stuck after 100000 iterations" << endl;
      QDPIO::cout << "got " << e << " evals " << l << " " << N << endl;
      exit(1);
    }

    T s, c;
    Givens_calc<T>(x, z, c, s);
    Givens_mult<T>(H, l, l + 1, c, -s, 0);
    Givens_mult<T>(H, l, l + 1, c,  s, 1);
    Givens_mult<T>(P, l, l + 1, c,  s, 1);

    for(int k = l; k < N - 2; ++k)
    {
      x = H.A[k + 1][k]; 
      z = H.A[k + 2][k]; 

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
  vector<T> tmp(N);
  for(int i = 0; i < N; ++i)
    tmp[i] = evals[N-i-1];
  evals = tmp; 

  UTeigenvectors(H, trows, evals, evecs);
  //UTSymmEigenvectors(H, trows, evals, evecs);

  for(int i = 0; i < evals.size(); ++i)
  {
    evecs[i] = P * evecs[i]; 
    normalize(evecs[i]); 
    evals[i] = evals[i] * Hnorm;
  }

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


  vector<T > p(N,0);
  for(int k=start;k<N-2;k++){
    //cerr << "hess" << k << endl;

    vector<T > ck(N-k-1), v(N-k-1);
    for(int i=k+1;i<N;i++){ck[i-k-1] = A(i,k);}	///kth column

    normalize(ck);		///Normalization cancels in PHP anyway
    T beta;
    Householder_vector<T >(ck, 0, ck.size()-1, v, beta);	///Householder vector

    Householder_mult<T>(A,v,beta,start,k+1,N-1,0);	///A -> PA
    Householder_mult<T >(A,v,beta,start,k+1,N-1,1);  ///PA -> PAP^H 
    ///Accumulate eigenvector
    Householder_mult<T >(Q,v,beta,start,k+1,N-1,1);  ///Q -> QP^H
  }


  /*for(int l=0;l<N-2;l++){
    for(int k=l+2;k<N;k++){
    A(0,k,l);
    }
    }*/


}

///Tridiagonalize a matrix
template <class T> void Tri(Matrix<T > &A, Matrix<T> &Q, int start){

  int N = A.dim;	//Matrix Size

  Hess(A,Q,start);

  /*for(int l=0;l<N-2;l++){

    for(int k=l+2;k<N;k++){
    A(0,l,k);
    }
    }*/


}

///Tridiagonalize a matrix
template <class T> void ForceTridiagonal(Matrix<T > &A){

  int N = A.dim;	//Matrix Size
  for(int l=0;l<N-2;l++){
    for(int k=l+2;k<N;k++){
      A(0,l,k);
      A(0,k,l);
    }
  }


}

///Solve a symmetric eigensystem, not necessarily in tridiagonal form
template <class T> int SymmEigensystem(Matrix<T > &Ain, vector<T> &evals, vector<vector<T> > &evecs, double small){

  int N = Ain.dim;
  Matrix<T > A(N); A = Ain; 
  Matrix<T > Q(N);Q.Unity();
  Tri(A,Q,0);
  int it = Wilkinson<T>(A, evals, evecs, small);
  for(int k=0;k<N;k++){evecs[k] = Q*evecs[k];}

  return it;
}

///Solve a general eigensystem, not necessarily in tridiagonal form
template <class T> int Eigensystem(Matrix<T > &Ain, vector<T> &evals, vector<vector<T> > &evecs, double small){

  int N = Ain.dim;
  Matrix<T > A(N); A = Ain; 
  Matrix<T > Q(N);Q.Unity();
  Hess(A,Q,0);
  int it = QReigensystem<T>(A, evals, evecs, small);
  for(int k=0;k<N;k++){evecs[k] = Q*evecs[k];}

  return it;
}


#endif



