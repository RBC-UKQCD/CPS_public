//------------------------------------------------------------------
//
// wilson_matrix.h
//
// Header file WilsonMatrix class for 
// Wilson-like quarks.
//
// For now this is specific to three colors. The constructor
// will exit if the number of colors is not equal to three.
//
//------------------------------------------------------------------


#ifndef INCLUDED_WILSONMATRIX_H
#define INCLUDED_WILSONMATRIX_H

#include <math.h>
#include <string.h>
#include <util/rcomplex.h>
#include <util/vector.h>
#include <util/wilson.h>

#include <alg/spin_matrix.h>

CPS_START_NAMESPACE
/* 
   The following structures used to define the quark propagator
   were taken from version 4 of the MILC Code: 
   su3_vector, wilson_vector, color_wilson_vector, wilson_matrix.
*/
typedef struct { Complex c[3]; } su3_vector;
//typedef struct { su3_vector d[4]; } wilson_vector;

/*!
  The WilsonVector class is used to build the WilsonMatrix class.
  It also provides several routines used for the sequential baryon
  source. Finally it provides the rotation routines from the chiral
  to the Dirac basis. The Chiral basis conventions are given in the 
  WilsonMatrix documentation.

  The Dirac basis is

       \f[
       \gamma_0 = \gamma_x =
       \left(\begin{tabular}{rrrr}
	0 & 0 &  0 & -i \\
	0 & 0 & -i &  0 \\
	0 & i &  0 &  0 \\
	i & 0 &  0 &  0 
	\end{tabular}\right),\;\;\;\;
	\gamma_1 = \gamma_y = 
	\left(\begin{tabular}{rrrr}
	0 &  0 &  0 & 1 \\
	0 &  0 & -1 & 0 \\
	0 & -1 &  0 & 0 \\
	1 &  0 &  0 & 0 
	\end{tabular}\right)
	\f]

	\f[
	\gamma_2 = \gamma_z = 
	\left(\begin{tabular}{rrrr}
	0 &  0 & -i & 0 \\
	0 &  0 &  0 & i \\
	i &  0 &  0 & 0 \\
	0 & -i &  0 & 0
	\end{tabular}\right), \;\;\;\;
	\gamma_3 = \gamma_t = 
	\left(\begin{tabular}{rrrr}
        1 &  0  &  0 &  0 \\
        0 &  1  &  0 &  0 \\
        0 &  0  & -1 &  0 \\
        0 &  0  &  0 & -1 
	\end{tabular}\right)
	\f]

	\f[
	\gamma_5 = \left(\begin{tabular}{rrrr}
	0 & 0 & 1 & 0 \\
	0 & 0 & 0 & 1 \\
	1 & 0 & 0 & 0 \\
	0 & 1 & 0 & 0
	\end{tabular}\right)
	\f]
*/

/*! This is equivalent with the wilson_vector */
class WilsonVector 
{
public:
  su3_vector d[4] ;

  WilsonVector() {;}

  ~WilsonVector(){;}

  void Zero(){for(int s(0);s<4;s++)for(int cc(0); cc<3 ; cc++)d[s].c[cc]=0.0;}

  WilsonVector& operator+=(WilsonVector& rhs){
    for(int s(0);s<4;s++)for(int cc(0); cc<3 ; cc++)
      d[s].c[cc]+=rhs.d[s].c[cc] ;
    return *this ;
  }

  WilsonVector& operator-=(WilsonVector& rhs){
    for(int s(0);s<4;s++)for(int cc(0); cc<3 ; cc++)
      d[s].c[cc]-=rhs.d[s].c[cc] ;
    return *this ;
  }

  WilsonVector& conj(){
    for(int s(0);s<4;s++)for(int cc(0); cc<3 ; cc++)
      d[s].c[cc].imag(-d[s].c[cc].imag()) ;
    return *this ;
  }

  WilsonVector& operator*=(Complex& rhs){
    for(int s(0);s<4;s++)for(int cc(0); cc<3 ; cc++)
      d[s].c[cc]*=rhs ;
    return *this ;
  }

  WilsonVector& operator*=(const Float rhs){
    for(int s(0);s<4;s++)for(int cc(0); cc<3 ; cc++)
      d[s].c[cc]*=rhs ;
    return *this ;
  }

  WilsonVector& gamma(int dir) ;
  
  /*!
    Rotate from Dirac to Chiral basis.
    It multiplies the by the rotation matrix 
    that rotates from the Chiral basis to the Dirac.

    \f[
    W_{s',c} = \sum_t V_{s',t} G_{t,c} 
    \f]
    Where 
    \f[
      V = V^\dagger = \frac{\sqrt 2}{2}\left(
      \begin{tabular}{rrrr}
       1 &   0 &   1 &   0 \\
       0 &   1 &   0 &   1 \\
       1 &   0 &  -1 &   0 \\
       0 &   1 &   0 &  -1 
      \end{tabular}\right)
     \f]

     A gamma matrix in the Chiral basis (\f$ G_c \f$) transforms
     to the Dirac basis (\f$ G_{d} \f$) as following
     \f[
     G_d = V^\dagger G_c V 
     \f]

 
  */
  WilsonVector& DiracToChiral() ;
  

  /*!
    Rotate from Chiral to Dirac basis.
    It multiplies by the rotation matrix that rotates 
    from the Dirac basis to the Chiral basis

    Since
    \f[
    V = V^\dagger
    \f]
    
    this routine is identical with DiracToChiral
  */ 
  WilsonVector& ChiralToDirac(){ return DiracToChiral() ; }

  /*!
    Multiplies  by \f$ \frac{1}{2}(1+\gamma_t)\$
    \f[
    V_{s,c}=\sum_{s_1}\left.\frac{1}{2}(1+\gamma_t)\right|_{s,s_1} W_{s_1,c}  
    \f]
  */
  WilsonVector& PParProject();
  

} ;

typedef WilsonVector wilson_vector;
typedef struct { wilson_vector c[3]; } color_wilson_vector;
typedef struct { color_wilson_vector d[4]; } wilson_matrix;


/*!
  The WilsonMatrix class is used to build quark propagators.
     
  The gamma matrices used are in the Chiral basis:

     \f[
     \gamma_0 = \gamma_x =
     \left(\begin{tabular}{rrrr}
       0 &  0  &  0 &  i \\
       0 &  0  &  i &  0 \\
       0 & -i  &  0 &  0 \\
      -i &  0  &  0 &  0 
      \end{tabular}\right),\;\;\;\;
      \gamma_1 = \gamma_y = 
      \left(\begin{tabular}{rrrr}
       0 &  0  &  0 & -1 \\
       0 &  0  &  1 &  0 \\
       0 &  1  &  0 &  0 \\
      -1 &  0  &  0 &  0 
      \end{tabular}\right)
      \f]

      \f[
      \gamma_2 = \gamma_z = 
      \left(\begin{tabular}{rrrr}
       0 &  0  &  i &  0 \\
       0 &  0  &  0 & -i \\
      -i &  0  &  0 &  0 \\
       0 &  i  &  0 &  0 
       \end{tabular}\right), \;\;\;\;
       \gamma_3 = \gamma_t = 
       \left(\begin{tabular}{rrrr}
       0 &  0  &  1 &  0 \\
       0 &  0  &  0 &  1 \\
       1 &  0  &  0 &  0 \\
       0 &  1  &  0 &  0 
       \end{tabular}\right)
       \f]

       \f[
       \gamma_5 = 
       \left(\begin{tabular}{rrrr}
       1 &  0  &  0 &  0 \\
       0 &  1  &  0 &  0 \\
       0 &  0  & -1 &  0 \\
       0 &  0  &  0 & -1 
       \end{tabular}\right)
       \f]
*/
class WilsonMatrix
{
  wilson_matrix p;
  
public:
  

  WilsonMatrix();
  WilsonMatrix(const WilsonMatrix& rhs);
  WilsonMatrix(const wilson_matrix& rhs);
  WilsonMatrix(const Float& rhs);
  WilsonMatrix(const Rcomplex& rhs);
  WilsonMatrix(int source_spin, int source_color, const wilson_vector&);
  WilsonMatrix(int source_spin, int source_color, int sink_spin, 
               int sink_color, const Rcomplex&);

  // Access to elements 

  void Element(int source_spin, int source_color, 
               int sink_spin, int sink_color, const Rcomplex& z);

   /*!
    Return the complex number references by
    s1 - source_spin
    c1 - source_colour
    s2 - sink_spin
    c2 - sink_colour
   */
  Complex& operator()( int s1, int c1, int s2, int c2 ) 
  {
    return p.d[s1].c[c1].d[s2].c[c2];
  }
  
  /*!
    Return the complex number references by
    s1 - source_spin
    c1 - source_colour
    s2 - sink_spin
    c2 - sink_colour
    ( const version )
  */
  Complex  operator()( int s1, int c1, int s2, int c2 ) const
  {
    return p.d[s1].c[c1].d[s2].c[c2];
  }
  
  //! hermitean conjugate the WilsonMatrix
  void hconj();
 
  //! mult the prop by gamma_dir on the left
  WilsonMatrix& gl(int dir); 

  //! mult the prop by gamma_dir on the left
  WilsonMatrix& gr(int dir); 

  //! right mult by sigma_mu_nu
  WilsonMatrix& sigmaR(int mu, int nu); 

  //! make a copy of the hermitean conjugate
  WilsonMatrix conj_cp(); 

  //! get a sol. vector
  wilson_vector& sol(int source_spin, int source_color); 

  void load_vec(int source_spin, int source_color, const wilson_vector&);
  void load_row(int source_spin, int source_color, const wilson_vector&);
  Rcomplex Trace();
  const wilson_matrix& wmat() const; // get p 
  WilsonMatrix& LeftTimesEqual(const WilsonMatrix& rhs);
  
  //! Projects positive parity on the sink
  WilsonMatrix& PParProjectSink()   ;

  //! Projects positive parity on the source
  WilsonMatrix& PParProjectSource() ;

  //! sink chiral to dirac rotation
  WilsonMatrix& SinkChiralToDirac() ;
  
  // operator functions

  WilsonMatrix& operator= (const WilsonMatrix& rhs);
  WilsonMatrix& operator= (const wilson_matrix& rhs);
  WilsonMatrix& operator= (const Float& rhs);
  WilsonMatrix& operator+=(const WilsonMatrix& rhs);
  WilsonMatrix& operator-=(const WilsonMatrix& rhs);
  WilsonMatrix& operator*=(const WilsonMatrix& rhs);
  WilsonMatrix& operator*=(const Float& rhs);
  WilsonMatrix& operator*=(const Rcomplex& rhs);
  
  WilsonMatrix& UMultSource   ( Matrix& U, WilsonMatrix& W) ;
  WilsonMatrix& UdagMultSource( Matrix& U, WilsonMatrix& W) ;
  
  // algebra shortcuts

  /*!
    Logically equivalent to   += fact * x , but 
    does not require a temporary WilsonMtarix
  */
  WilsonMatrix& AddMult( const Rcomplex& fact, const WilsonMatrix& x );
  
  // Baryon things
  WilsonMatrix& ccl(int dir);
  WilsonMatrix& ccr(int dir);
  WilsonMatrix& diq(const WilsonMatrix& rhs);
  WilsonMatrix& joint(const WilsonMatrix& rhs);
  
  friend Rcomplex Trace(const WilsonMatrix& p1, const WilsonMatrix& p2);
  
};


// some proto-types for functions that operate on WilsonMatrices
#ifdef _TARTAN

/*! C = A * B */
extern "C" void wmatMultwmat(IFloat* C, const IFloat* A, const IFloat* B);

/*! C = Tr (A * B) */
extern "C" void Tracewmatwmat(IFloat* C, const IFloat* A, const IFloat* B);
	
#endif

//! times operator
extern WilsonMatrix operator*(const WilsonMatrix& lhs, const WilsonMatrix& rhs);

//! times operator
extern WilsonMatrix operator*(const Float& num, const WilsonMatrix& mat);

//! times operator
extern WilsonMatrix operator*(const WilsonMatrix& mat, const Float& num);

//! times operator
extern WilsonMatrix operator*(const Rcomplex& num, const WilsonMatrix& mat);

//! times operator
extern WilsonMatrix operator*(const WilsonMatrix& mat, const Rcomplex& num);

//! plus operator
extern WilsonMatrix operator+(const WilsonMatrix& lhs, const WilsonMatrix& rhs);

//! minus operator
extern WilsonMatrix operator-(const WilsonMatrix& lhs, const WilsonMatrix& rhs);

//! left multiply by gamma_dir
extern void mult_by_gamma_left( int dir, const wilson_matrix& src, 
				wilson_matrix& dest );

//! right multiply by gamma_dir
extern void mult_by_gamma_right(int dir, const wilson_matrix& src, 
				wilson_matrix& dest );

//! Spin and Color trace of a 2 WilsonMatrices
extern Rcomplex Trace(const WilsonMatrix& p1, const WilsonMatrix& p2);

//! Spin trace of a WilsonMatrix
extern Matrix SpinTrace(const WilsonMatrix& Wmat); 

//! Spin trace of two WilsonMatrices
extern Matrix SpinTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2); 

//! Spin trace of two WilsonMatrices
extern Matrix SpinTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2,
						const WilsonMatrix& Wmat3); 
//! Color trace of a WilsonMatrix
extern SpinMatrix ColorTrace(const WilsonMatrix& Wmat);

//! Color trace of three WilsonMatrices
extern SpinMatrix ColorTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2);

//! Color trace of three WilsonMatrices
extern SpinMatrix ColorTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2,
						const WilsonMatrix& Wmat3);

//! trace of two (color) Matrices
extern Rcomplex Tr(const Matrix& a, const Matrix& b);

//! trace of two SpinMatrices
extern Rcomplex Tr(const SpinMatrix& a, const SpinMatrix& b);

CPS_END_NAMESPACE

#endif


