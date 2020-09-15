// -*- mode:c++;c-basic-offset:4 -*-
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
#include <util/error.h>

#define INLINE_WILSON_MATRIX

#define TIMESPLUSONE(a,b) { b=a; }
#define TIMESMINUSONE(a,b) { b=-a; }
//#define TIMESPLUSI(a,b) { b.real(-a.imag()); b.imag(a.real()); }
//#define TIMESMINUSI(a,b) { b.real(a.imag()); b.imag(-a.real()); }
#define TIMESPLUSI(a,b) { b=Complex(-a.imag(),a.real()); }
#define TIMESMINUSI(a,b) { b=Complex(a.imag(),-a.real()); }

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
  su3_vector d[4];

  WilsonVector() {}

  ~WilsonVector() {}

  void Zero()
    {
        for(int s1=0;s1<4;s1++)
            for(int c1=0;c1<3;c1++)
                d[s1].c[c1] = 0.0;
    }

  WilsonVector& operator+=(WilsonVector& rhs) {
    for(int s1=0;s1<4;s1++) for(int c1=0; c1<3;c1++)
      d[s1].c[c1] += rhs.d[s1].c[c1];
    return *this;
  }

  WilsonVector& operator-=(WilsonVector& rhs) {
    for(int s1=0;s1<4;s1++) for(int c1=0;c1<3;c1++)
      d[s1].c[c1] -= rhs.d[s1].c[c1];
    return *this;
  }

  WilsonVector& conj() {
    for(int s1=0;s1<4;s1++) for(int c1=0;c1<3;c1++)
	  d[s1].c[c1]=Complex(d[s1].c[c1].real(),-d[s1].c[c1].imag());
//      d[s1].c[c1].imag(-d[s1].c[c1].imag());
    return *this;
  }

  WilsonVector& operator*=(Complex& rhs) {
    for(int s1=0;s1<4;s1++) for(int c1=0;c1<3;c1++)
      d[s1].c[c1] *= rhs;
    return *this;
  }

  WilsonVector& operator*=(const Float rhs) {
    for(int s1=0;s1<4;s1++)for(int c1=0;c1<3;c1++)
      d[s1].c[c1] *= rhs;
    return *this;
  }

  WilsonVector& gamma(int dir);
  WilsonVector& LeftTimesEqual(const Matrix& rhs);

  //Added by CK
  WilsonVector & ccl(int dir);

  //Added by CK
  IFloat norm() const{
    Float out(0.0);
    for(int s1=0;s1<4;s1++)
      for(int c1=0;c1<3;c1++)
	out += std::norm(d[s1].c[c1]);
    return out;
  }
  
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
    WilsonVector& ChiralToDirac(){ return DiracToChiral(); }

    /*!
      Multiplies  by \f$ \frac{1}{2}(1+\gamma_t)\$
      \f[
      V_{s,c}=\sum_{s_1}\left.\frac{1}{2}(1+\gamma_t)\right|_{s,s_1} W_{s_1,c}  
      \f]
    */
    WilsonVector& PParProject();

};

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
class WilsonMatrix;

// single precision WilsonMatrix, only a container to save memory.
class WilsonMatrixS
{
public:
    WilsonMatrixS() {}
    WilsonMatrixS(const WilsonMatrix &w);
    const WilsonMatrixS &operator=(const WilsonMatrix &w);

    float a[288];

    // DJM: add accessor
    std::complex<float> operator()(int s1, int c1, int s2, int c2) const {
      return std::complex<float>( a[72*s1+24*c1+6*s2+2*c2], a[72*s1+24*c1+6*s2+2*c2+1] );
    }
};

WilsonMatrix& eq_mult( WilsonMatrix& xmat,
		       const WilsonMatrix& amat,
		       const WilsonMatrix& bmat );

class WilsonMatrix
{
private:
    static const char *cname;
  wilson_matrix p;
  
public:
    WilsonMatrix() {}

    // no need to redefine copy constructor.
  
    WilsonMatrix(const wilson_matrix& rhs) {
        p = rhs;
    }
//    WilsonMatrix(const WilsonMatrix& rhs)

#ifndef INLINE_WILSON_MATRIX
  WilsonMatrix(const Float& rhs);
#else 
  WilsonMatrix(const Float& rhs)
{
	for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
	  for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
	    p.d[s1].c[c1].d[s2].c[c2]=Complex(rhs,0.0);
        }
}
#endif

  WilsonMatrix(const Rcomplex& rhs);

  WilsonMatrix(int sink_spin, int sink_color, const wilson_vector&);

  WilsonMatrix(int sink_spin, int sink_color, int source_spin, 
               int source_color, const Rcomplex&);

    // added by Hantao
    WilsonMatrix(const WilsonMatrixS &ws) {
        for(int i = 0; i < 144; ++i) {
            int j = i;
            int c2 = j % 3; j /= 3;
            int s2 = j % 4; j /= 4;
            int c1 = j % 3; j /= 3;
            int s1 = j % 4;

            p.d[s1].c[c1].d[s2].c[c2] = Complex(ws.a[2*i], ws.a[2*i+1]);
        }
    }

    const WilsonMatrix &operator=(const WilsonMatrixS &ws) {
        for(int i = 0; i < 144; ++i) {
            int j = i;
            int c2 = j % 3; j /= 3;
            int s2 = j % 4; j /= 4;
            int c1 = j % 3; j /= 3;
            int s1 = j % 4;

            p.d[s1].c[c1].d[s2].c[c2] = Complex(ws.a[2*i], ws.a[2*i+1]);
        }
        return *this;
    }

  // Access to elements 

  Rcomplex* ptr() { return reinterpret_cast<Rcomplex*>(p.d); }
  const Rcomplex* ptr() const { 
    return reinterpret_cast<const Rcomplex*>(p.d);
  }

  void Element(int sink_spin, int sink_color, 
               int source_spin, int source_color, const Rcomplex& z);

   /*!
    Return the complex number references by
    s1 - sink_spin
    c1 - sink_colour
    s2 - source_spin
    c2 - source_colour
   */
  Complex& operator()(int s1, int c1, int s2, int c2) 
	{ return p.d[s1].c[c1].d[s2].c[c2]; }
  
    /*!
      Return the complex number references by
      s1 - sink_spin
      c1 - sink_colour
      s2 - source_spin
      c2 - source_colour
      ( const version )
    */
    Complex  operator()(int s1, int c1, int s2, int c2) const
    { return p.d[s1].c[c1].d[s2].c[c2]; }
  
  //CK: Set to unit matrix
  void Unit(){ 
    for(int s1=0;s1<4;s1++)
      for(int c1=0;c1<3;c1++)
	for(int s2=0;s2<4;s2++)
	  for(int c2=0;c2<3;c2++)
	    p.d[s1].c[c1].d[s2].c[c2] = (s1 == s2 && c1 == c2) ? Complex(1.0,0.0) : Complex(0.0,0.0);
  }


  //! hermitean conjugate the WilsonMatrix
  //WilsonMatrix& hconj();
#ifndef INLINE_WILSON_MATRIX
  void hconj();
#else
    void hconj()
    {
        int c1, c2;
        int s1, s2;
        wilson_matrix mat=p;

        for(s1=0;s1<4;s1++)
            for(c1=0;c1<3;c1++)
                for(s2=0;s2<4;s2++)
                    for(c2=0;c2<3;c2++)
                        p.d[s2].c[c2].d[s1].c[c1] = conj(mat.d[s1].c[c1].d[s2].c[c2]);
	
    }
#endif

  void dump(); // print out the prop
 

  //! complex conjugate the WilsonMatrix
  void cconj();

  //! transpose the WilsonMatrix
  void transpose();
  
  // added by Daiqian
  //! transpose the color index of WilsonMatrix
  void transpose_color() {
    wilson_matrix mat=p;
    for(int s2=0;s2<4;s2++)
      for(int c1=0;c1<3;c1++)
	for(int s1=0;s1<4;s1++)
	  for(int c2=0;c2<3;c2++)
	    p.d[s2].c[c2].d[s1].c[c1] = mat.d[s2].c[c1].d[s1].c[c2];
  }


  // added by CK
  //! norm^2 of the WilsonMatrix
  IFloat norm() const{
    IFloat out(0.0);
    for(int s1=0;s1<4;s1++)
      for(int c1=0;c1<3;c1++)
	out += p.d[s1].c[c1].norm();
    return out;
  }

#ifndef INLINE_WILSON_MATRIX
  //! mult the prop by gamma_dir on the left
  WilsonMatrix& gl(int dir); 
  //! mult the prop by gamma_dir*(1-gamma_5) on the left, and return the new matrix
    WilsonMatrix glL(int dir)const;
  //! mult the prop by gamma_dir*(1+gamma_5) on the left, and return the new matrix
    WilsonMatrix glR(int dir)const;

  //! mult the prop by gamma_dir*gamma_5 on the left, and return the new matrix
    WilsonMatrix glA(int dir)const;
  //! glA another version. result = gamma_dir*gamma_5*from
    WilsonMatrix& glA(const WilsonMatrix & from, int dir);
#else
WilsonMatrix& gl(int dir)
{
    int i; /*color*/
    int c2,s2;    /* column indices, color and spin */
    wilson_matrix src=p;

    switch(dir){
    case 0:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    TIMESPLUSI(  src.d[3].c[i].d[s2].c[c2],
                                 p.d[0].c[i].d[s2].c[c2] );
                    TIMESPLUSI(  src.d[2].c[i].d[s2].c[c2],
                                 p.d[1].c[i].d[s2].c[c2] );
                    TIMESMINUSI( src.d[1].c[i].d[s2].c[c2],
                                 p.d[2].c[i].d[s2].c[c2] );
                    TIMESMINUSI( src.d[0].c[i].d[s2].c[c2],
                                 p.d[3].c[i].d[s2].c[c2] );
                }
        break;
    case 1:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    TIMESMINUSONE( src.d[3].c[i].d[s2].c[c2],
                                   p.d[0].c[i].d[s2].c[c2] );
                    TIMESPLUSONE(  src.d[2].c[i].d[s2].c[c2],
                                   p.d[1].c[i].d[s2].c[c2] );
                    TIMESPLUSONE(  src.d[1].c[i].d[s2].c[c2],
                                   p.d[2].c[i].d[s2].c[c2] );
                    TIMESMINUSONE( src.d[0].c[i].d[s2].c[c2],
                                   p.d[3].c[i].d[s2].c[c2] );
                }
        break;
    case 2:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    TIMESPLUSI(  src.d[2].c[i].d[s2].c[c2],
                                 p.d[0].c[i].d[s2].c[c2] );
                    TIMESMINUSI( src.d[3].c[i].d[s2].c[c2],
                                 p.d[1].c[i].d[s2].c[c2] );
                    TIMESMINUSI( src.d[0].c[i].d[s2].c[c2],
                                 p.d[2].c[i].d[s2].c[c2] );
                    TIMESPLUSI(  src.d[1].c[i].d[s2].c[c2],
                                 p.d[3].c[i].d[s2].c[c2] );
                }
	break;
    case 3:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    TIMESPLUSONE( src.d[2].c[i].d[s2].c[c2],
                                  p.d[0].c[i].d[s2].c[c2] );
                    TIMESPLUSONE( src.d[3].c[i].d[s2].c[c2],
                                  p.d[1].c[i].d[s2].c[c2] );
                    TIMESPLUSONE( src.d[0].c[i].d[s2].c[c2],
                                  p.d[2].c[i].d[s2].c[c2] );
                    TIMESPLUSONE( src.d[1].c[i].d[s2].c[c2],
                                  p.d[3].c[i].d[s2].c[c2] );
                }
        break;
    case -5:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    TIMESPLUSONE(  src.d[0].c[i].d[s2].c[c2],
                                   p.d[0].c[i].d[s2].c[c2] );
                    TIMESPLUSONE(  src.d[1].c[i].d[s2].c[c2],
                                   p.d[1].c[i].d[s2].c[c2] );
                    TIMESMINUSONE( src.d[2].c[i].d[s2].c[c2],
                                   p.d[2].c[i].d[s2].c[c2] );
                    TIMESMINUSONE( src.d[3].c[i].d[s2].c[c2],
                                   p.d[3].c[i].d[s2].c[c2] );
                }
        break;
  case PL: // (1-g5)/2
    for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	  p.d[0].c[i].d[s2].c[c2]=0.0;
	  p.d[1].c[i].d[s2].c[c2]=0.0;
	  TIMESPLUSONE( src.d[2].c[i].d[s2].c[c2],
			p.d[2].c[i].d[s2].c[c2] );
	  TIMESPLUSONE( src.d[3].c[i].d[s2].c[c2],
			p.d[3].c[i].d[s2].c[c2] );
        }
        break;
  case PR: // (1+g5)/2
    for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	  TIMESPLUSONE(  src.d[0].c[i].d[s2].c[c2],
			 p.d[0].c[i].d[s2].c[c2] );
	  TIMESPLUSONE(  src.d[1].c[i].d[s2].c[c2],
			 p.d[1].c[i].d[s2].c[c2] );
	  p.d[2].c[i].d[s2].c[c2]=0.0;
	  p.d[3].c[i].d[s2].c[c2]=0.0;
	}
        break;
    default:
//	ERR.General("WilsonMatrix", "gl()", "BAD CALL TO gl()\n");
	break;
    }
    return *this;
}

WilsonMatrix glR(int dir) const
{
  int i; /*color*/
  int c2,s2;    /* column indices, color and spin */
  wilson_matrix result;

  switch(dir){
    case 0:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
			result.d[0].c[i].d[s2].c[c2]=0.0;
			result.d[1].c[i].d[s2].c[c2]=0.0;
            TIMESMINUSI( 2.*p.d[1].c[i].d[s2].c[c2],
                result.d[2].c[i].d[s2].c[c2] );
            TIMESMINUSI( 2.*p.d[0].c[i].d[s2].c[c2],
                result.d[3].c[i].d[s2].c[c2] );
        }
        break;
    case 1:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
			result.d[0].c[i].d[s2].c[c2]=0.0;
			result.d[1].c[i].d[s2].c[c2]=0.0;
            TIMESPLUSONE(  2.*p.d[1].c[i].d[s2].c[c2],
                result.d[2].c[i].d[s2].c[c2] );
            TIMESMINUSONE( 2.*p.d[0].c[i].d[s2].c[c2],
                result.d[3].c[i].d[s2].c[c2] );
        }
        break;
    case 2:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
			result.d[0].c[i].d[s2].c[c2]=0.0;
			result.d[1].c[i].d[s2].c[c2]=0.0;
            TIMESMINUSI( 2.*p.d[0].c[i].d[s2].c[c2],
                result.d[2].c[i].d[s2].c[c2] );
            TIMESPLUSI(  2.*p.d[1].c[i].d[s2].c[c2],
                result.d[3].c[i].d[s2].c[c2] );
        }
	break;
    case 3:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
			result.d[0].c[i].d[s2].c[c2]=0.0;
			result.d[1].c[i].d[s2].c[c2]=0.0;
            TIMESPLUSONE( 2.*p.d[0].c[i].d[s2].c[c2],
                result.d[2].c[i].d[s2].c[c2] );
            TIMESPLUSONE( 2.*p.d[1].c[i].d[s2].c[c2],
                result.d[3].c[i].d[s2].c[c2] );
        }
        break;
    default:
//	ERR.General("wilsonmatrix","glL(int)","BAD CALL TO gl()\n");
	break;
  }
	return WilsonMatrix(result);
}
WilsonMatrix glL(int dir) const
{
  int i; /*color*/
  int c2,s2;    /* column indices, color and spin */
  wilson_matrix result;

  switch(dir){
    case 0:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESPLUSI(  2.*p.d[3].c[i].d[s2].c[c2],
                result.d[0].c[i].d[s2].c[c2] );
            TIMESPLUSI(  2.*p.d[2].c[i].d[s2].c[c2],
                result.d[1].c[i].d[s2].c[c2] );
			result.d[2].c[i].d[s2].c[c2]=0.0;
			result.d[3].c[i].d[s2].c[c2]=0.0;
        }
        break;
    case 1:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESMINUSONE( 2.*p.d[3].c[i].d[s2].c[c2],
                result.d[0].c[i].d[s2].c[c2] );
            TIMESPLUSONE(  2.*p.d[2].c[i].d[s2].c[c2],
                result.d[1].c[i].d[s2].c[c2] );
			result.d[2].c[i].d[s2].c[c2]=0.0;
			result.d[3].c[i].d[s2].c[c2]=0.0;
        }
        break;
    case 2:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESPLUSI(  2.*p.d[2].c[i].d[s2].c[c2],
                result.d[0].c[i].d[s2].c[c2] );
            TIMESMINUSI( 2.*p.d[3].c[i].d[s2].c[c2],
                result.d[1].c[i].d[s2].c[c2] );
			result.d[2].c[i].d[s2].c[c2]=0.0;
			result.d[3].c[i].d[s2].c[c2]=0.0;
        }
	break;
    case 3:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESPLUSONE( 2.*p.d[2].c[i].d[s2].c[c2],
                result.d[0].c[i].d[s2].c[c2] );
            TIMESPLUSONE( 2.*p.d[3].c[i].d[s2].c[c2],
                result.d[1].c[i].d[s2].c[c2] );
			result.d[2].c[i].d[s2].c[c2]=0.0;
			result.d[3].c[i].d[s2].c[c2]=0.0;
        }
        break;
    default:
//	ERR.General("wilsonmatrix","glR(int)","BAD CALL TO gl()\n");
	break;
  }
	return WilsonMatrix(result);
}

//multiply gamma(i)gamma(5) on the left and return a new one
WilsonMatrix glA(int dir) const
{
  int i; /*color*/
  int c2,s2;    /* column indices, color and spin */
  wilson_matrix result;

  switch(dir){
    case 0:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESMINUSI(  p.d[3].c[i].d[s2].c[c2],
                result.d[0].c[i].d[s2].c[c2] );
            TIMESMINUSI(  p.d[2].c[i].d[s2].c[c2],
                result.d[1].c[i].d[s2].c[c2] );
            TIMESMINUSI( p.d[1].c[i].d[s2].c[c2],
                result.d[2].c[i].d[s2].c[c2] );
            TIMESMINUSI( p.d[0].c[i].d[s2].c[c2],
                result.d[3].c[i].d[s2].c[c2] );
        }
        break;
    case 1:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESPLUSONE( p.d[3].c[i].d[s2].c[c2],
                result.d[0].c[i].d[s2].c[c2] );
            TIMESMINUSONE(  p.d[2].c[i].d[s2].c[c2],
                result.d[1].c[i].d[s2].c[c2] );
            TIMESPLUSONE(  p.d[1].c[i].d[s2].c[c2],
                result.d[2].c[i].d[s2].c[c2] );
            TIMESMINUSONE( p.d[0].c[i].d[s2].c[c2],
                result.d[3].c[i].d[s2].c[c2] );
        }
        break;
    case 2:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESMINUSI(  p.d[2].c[i].d[s2].c[c2],
                result.d[0].c[i].d[s2].c[c2] );
            TIMESPLUSI( p.d[3].c[i].d[s2].c[c2],
                result.d[1].c[i].d[s2].c[c2] );
            TIMESMINUSI( p.d[0].c[i].d[s2].c[c2],
                result.d[2].c[i].d[s2].c[c2] );
            TIMESPLUSI(  p.d[1].c[i].d[s2].c[c2],
                result.d[3].c[i].d[s2].c[c2] );
        }
	break;
    case 3:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESMINUSONE( p.d[2].c[i].d[s2].c[c2],
                result.d[0].c[i].d[s2].c[c2] );
            TIMESMINUSONE( p.d[3].c[i].d[s2].c[c2],
                result.d[1].c[i].d[s2].c[c2] );
            TIMESPLUSONE( p.d[0].c[i].d[s2].c[c2],
                result.d[2].c[i].d[s2].c[c2] );
            TIMESPLUSONE( p.d[1].c[i].d[s2].c[c2],
                result.d[3].c[i].d[s2].c[c2] );
        }
        break;
    default:
//	ERR.General(cname, "glA()", "BAD CALL TO glA(dir = %d\n",dir);
		//VRB.Result(cname,fname,"BAD CALL TO glA(int)\n");
	break;
  }
	return WilsonMatrix(result);
}
WilsonMatrix& glA(const WilsonMatrix & from, int dir)
{
  int i; /*color*/
  int c2,s2;    /* column indices, color and spin */
  const wilson_matrix & from_mat=from.wmat();

  switch(dir){
    case 0:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESMINUSI(  from_mat.d[3].c[i].d[s2].c[c2],
                p.d[0].c[i].d[s2].c[c2] );
            TIMESMINUSI(  from_mat.d[2].c[i].d[s2].c[c2],
                p.d[1].c[i].d[s2].c[c2] );
            TIMESMINUSI( from_mat.d[1].c[i].d[s2].c[c2],
                p.d[2].c[i].d[s2].c[c2] );
            TIMESMINUSI( from_mat.d[0].c[i].d[s2].c[c2],
                p.d[3].c[i].d[s2].c[c2] );
        }
        break;
    case 1:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESPLUSONE( from_mat.d[3].c[i].d[s2].c[c2],
                p.d[0].c[i].d[s2].c[c2] );
            TIMESMINUSONE(  from_mat.d[2].c[i].d[s2].c[c2],
                p.d[1].c[i].d[s2].c[c2] );
            TIMESPLUSONE(  from_mat.d[1].c[i].d[s2].c[c2],
                p.d[2].c[i].d[s2].c[c2] );
            TIMESMINUSONE( from_mat.d[0].c[i].d[s2].c[c2],
                p.d[3].c[i].d[s2].c[c2] );
        }
        break;
    case 2:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESMINUSI(  from_mat.d[2].c[i].d[s2].c[c2],
                p.d[0].c[i].d[s2].c[c2] );
            TIMESPLUSI( from_mat.d[3].c[i].d[s2].c[c2],
                p.d[1].c[i].d[s2].c[c2] );
            TIMESMINUSI( from_mat.d[0].c[i].d[s2].c[c2],
                p.d[2].c[i].d[s2].c[c2] );
            TIMESPLUSI(  from_mat.d[1].c[i].d[s2].c[c2],
                p.d[3].c[i].d[s2].c[c2] );
        }
	break;
    case 3:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESMINUSONE( from_mat.d[2].c[i].d[s2].c[c2],
                p.d[0].c[i].d[s2].c[c2] );
            TIMESMINUSONE( from_mat.d[3].c[i].d[s2].c[c2],
                p.d[1].c[i].d[s2].c[c2] );
            TIMESPLUSONE( from_mat.d[0].c[i].d[s2].c[c2],
                p.d[2].c[i].d[s2].c[c2] );
            TIMESPLUSONE( from_mat.d[1].c[i].d[s2].c[c2],
                p.d[3].c[i].d[s2].c[c2] );
        }
        break;
    default:
//	ERR.General(cname, "glA()", "BAD CALL TO glA(from=%p, dir = %d\n", &from,dir);
		//VRB.Result(cname,fname,"BAD CALL TO glA(int)\n");
	break;
  }
	return *this;
}
#endif

  //! glA another version. this -> gamma_dir*gamma_5*this
  inline WilsonMatrix& glAx(int dir){
    WilsonMatrix tmp(*this);
    glA(tmp,dir);
    return *this;
  }

#ifndef INLINE_WILSON_MATRIX
  //! mult the prop by gamma_dir on the left, and return the new matrix
    WilsonMatrix glV(int dir)const;
  //! glV another version. result = gamma_dir*from
    WilsonMatrix& glV(const WilsonMatrix & from, int dir);
//  void glV(const WilsonMatrix & from, int dir);
#else
//multiply gamma(i) on the left and return a new one
WilsonMatrix glV(int dir) const
{
  int i; /*color*/
  int c2,s2;    /* column indices, color and spin */
  wilson_matrix result;

  switch(dir){
    case 0:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESPLUSI(  p.d[3].c[i].d[s2].c[c2],
                result.d[0].c[i].d[s2].c[c2] );
            TIMESPLUSI(  p.d[2].c[i].d[s2].c[c2],
                result.d[1].c[i].d[s2].c[c2] );
            TIMESMINUSI( p.d[1].c[i].d[s2].c[c2],
                result.d[2].c[i].d[s2].c[c2] );
            TIMESMINUSI( p.d[0].c[i].d[s2].c[c2],
                result.d[3].c[i].d[s2].c[c2] );
        }
        break;
    case 1:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESMINUSONE( p.d[3].c[i].d[s2].c[c2],
                result.d[0].c[i].d[s2].c[c2] );
            TIMESPLUSONE(  p.d[2].c[i].d[s2].c[c2],
                result.d[1].c[i].d[s2].c[c2] );
            TIMESPLUSONE(  p.d[1].c[i].d[s2].c[c2],
                result.d[2].c[i].d[s2].c[c2] );
            TIMESMINUSONE( p.d[0].c[i].d[s2].c[c2],
                result.d[3].c[i].d[s2].c[c2] );
        }
        break;
    case 2:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESPLUSI(  p.d[2].c[i].d[s2].c[c2],
                result.d[0].c[i].d[s2].c[c2] );
            TIMESMINUSI( p.d[3].c[i].d[s2].c[c2],
                result.d[1].c[i].d[s2].c[c2] );
            TIMESMINUSI( p.d[0].c[i].d[s2].c[c2],
                result.d[2].c[i].d[s2].c[c2] );
            TIMESPLUSI(  p.d[1].c[i].d[s2].c[c2],
                result.d[3].c[i].d[s2].c[c2] );
        }
	break;
    case 3:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESPLUSONE( p.d[2].c[i].d[s2].c[c2],
                result.d[0].c[i].d[s2].c[c2] );
            TIMESPLUSONE( p.d[3].c[i].d[s2].c[c2],
                result.d[1].c[i].d[s2].c[c2] );
            TIMESPLUSONE( p.d[0].c[i].d[s2].c[c2],
                result.d[2].c[i].d[s2].c[c2] );
            TIMESPLUSONE( p.d[1].c[i].d[s2].c[c2],
                result.d[3].c[i].d[s2].c[c2] );
        }
        break;
    case -5:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESPLUSONE(  p.d[0].c[i].d[s2].c[c2],
                result.d[0].c[i].d[s2].c[c2] );
            TIMESPLUSONE(  p.d[1].c[i].d[s2].c[c2],
                result.d[1].c[i].d[s2].c[c2] );
            TIMESMINUSONE( p.d[2].c[i].d[s2].c[c2],
                result.d[2].c[i].d[s2].c[c2] );
            TIMESMINUSONE( p.d[3].c[i].d[s2].c[c2],
                result.d[3].c[i].d[s2].c[c2] );
        }
        break;
    default:
		//VRB.Result(cname,fname,"BAD CALL TO glV()\n");
//	ERR.General(cname, "glV()", "BAD CALL TO glV(dir = %d\n", dir);
	break;
  }
	return WilsonMatrix(result);
}

    //! glV another version. result = gamma_dir*from
    WilsonMatrix& glV(const WilsonMatrix & from, int dir)
    {
        int i; /*color*/
        int c2,s2;    /* column indices, color and spin */
        const wilson_matrix & from_mat=from.wmat();

        switch(dir){
        case 0:
            for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                        TIMESPLUSI(  from_mat.d[3].c[i].d[s2].c[c2],
                                     p.d[0].c[i].d[s2].c[c2] );
                        TIMESPLUSI(  from_mat.d[2].c[i].d[s2].c[c2],
                                     p.d[1].c[i].d[s2].c[c2] );
                        TIMESMINUSI( from_mat.d[1].c[i].d[s2].c[c2],
                                     p.d[2].c[i].d[s2].c[c2] );
                        TIMESMINUSI( from_mat.d[0].c[i].d[s2].c[c2],
                                     p.d[3].c[i].d[s2].c[c2] );
                    }
            break;
        case 1:
            for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                        TIMESMINUSONE( from_mat.d[3].c[i].d[s2].c[c2],
                                       p.d[0].c[i].d[s2].c[c2] );
                        TIMESPLUSONE(  from_mat.d[2].c[i].d[s2].c[c2],
                                       p.d[1].c[i].d[s2].c[c2] );
                        TIMESPLUSONE(  from_mat.d[1].c[i].d[s2].c[c2],
                                       p.d[2].c[i].d[s2].c[c2] );
                        TIMESMINUSONE( from_mat.d[0].c[i].d[s2].c[c2],
                                       p.d[3].c[i].d[s2].c[c2] );
                    }
            break;
        case 2:
            for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                        TIMESPLUSI(  from_mat.d[2].c[i].d[s2].c[c2],
                                     p.d[0].c[i].d[s2].c[c2] );
                        TIMESMINUSI( from_mat.d[3].c[i].d[s2].c[c2],
                                     p.d[1].c[i].d[s2].c[c2] );
                        TIMESMINUSI( from_mat.d[0].c[i].d[s2].c[c2],
                                     p.d[2].c[i].d[s2].c[c2] );
                        TIMESPLUSI(  from_mat.d[1].c[i].d[s2].c[c2],
                                     p.d[3].c[i].d[s2].c[c2] );
                    }
            break;
        case 3:
            for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                        TIMESPLUSONE( from_mat.d[2].c[i].d[s2].c[c2],
                                      p.d[0].c[i].d[s2].c[c2] );
                        TIMESPLUSONE( from_mat.d[3].c[i].d[s2].c[c2],
                                      p.d[1].c[i].d[s2].c[c2] );
                        TIMESPLUSONE( from_mat.d[0].c[i].d[s2].c[c2],
                                      p.d[2].c[i].d[s2].c[c2] );
                        TIMESPLUSONE( from_mat.d[1].c[i].d[s2].c[c2],
                                      p.d[3].c[i].d[s2].c[c2] );
                    }
            break;
        case -5:
            for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                        TIMESPLUSONE(  from_mat.d[0].c[i].d[s2].c[c2],
                                       p.d[0].c[i].d[s2].c[c2] );
                        TIMESPLUSONE(  from_mat.d[1].c[i].d[s2].c[c2],
                                       p.d[1].c[i].d[s2].c[c2] );
                        TIMESMINUSONE( from_mat.d[2].c[i].d[s2].c[c2],
                                       p.d[2].c[i].d[s2].c[c2] );
                        TIMESMINUSONE( from_mat.d[3].c[i].d[s2].c[c2],
                                       p.d[3].c[i].d[s2].c[c2] );
                    }
            break;
        default:
//            ERR.General(cname, "glV", "BAD CALL TO glV(), dir = %d\n", dir);
            break;
        }
        return *this;
    }

#endif

#ifndef INLINE_WILSON_MATRIX
    //! mult the prop by gamma_dir*gamma_5 on the left
    WilsonMatrix& grA(const WilsonMatrix & from, int dir);
    //! mult the prop by gamma_dir on the left
    WilsonMatrix& grV(const WilsonMatrix & from, int dir);

  //! mult the prop by gamma_dir on the left
  WilsonMatrix& gr(int dir); 

#else
WilsonMatrix& gr(int dir)
{
    int i; /*color*/
    int c1, s1;    /* row indices, color and spin */
    wilson_matrix src=p;

    switch(dir){
    case 0:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESMINUSI( src.d[s1].c[c1].d[3].c[i],
                                 p.d[s1].c[c1].d[0].c[i] );
                    TIMESMINUSI( src.d[s1].c[c1].d[2].c[i],
                                 p.d[s1].c[c1].d[1].c[i] );
                    TIMESPLUSI(  src.d[s1].c[c1].d[1].c[i],
                                 p.d[s1].c[c1].d[2].c[i] );
                    TIMESPLUSI(  src.d[s1].c[c1].d[0].c[i],
                                 p.d[s1].c[c1].d[3].c[i] );
                }
        break;
    case 1:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESMINUSONE( src.d[s1].c[c1].d[3].c[i],
                                   p.d[s1].c[c1].d[0].c[i] );
                    TIMESPLUSONE(  src.d[s1].c[c1].d[2].c[i],
                                   p.d[s1].c[c1].d[1].c[i] );
                    TIMESPLUSONE(  src.d[s1].c[c1].d[1].c[i],
                                   p.d[s1].c[c1].d[2].c[i] );
                    TIMESMINUSONE( src.d[s1].c[c1].d[0].c[i],
                                   p.d[s1].c[c1].d[3].c[i] );
                }
        break;
    case 2:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESMINUSI( src.d[s1].c[c1].d[2].c[i],
                                 p.d[s1].c[c1].d[0].c[i] );
                    TIMESPLUSI(  src.d[s1].c[c1].d[3].c[i],
                                 p.d[s1].c[c1].d[1].c[i] );
                    TIMESPLUSI(  src.d[s1].c[c1].d[0].c[i],
                                 p.d[s1].c[c1].d[2].c[i] );
                    TIMESMINUSI( src.d[s1].c[c1].d[1].c[i],
                                 p.d[s1].c[c1].d[3].c[i] );
                }
        break;
    case 3:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESPLUSONE( src.d[s1].c[c1].d[2].c[i],
                                  p.d[s1].c[c1].d[0].c[i] );
                    TIMESPLUSONE( src.d[s1].c[c1].d[3].c[i],
                                  p.d[s1].c[c1].d[1].c[i] );
                    TIMESPLUSONE( src.d[s1].c[c1].d[0].c[i],
                                  p.d[s1].c[c1].d[2].c[i] );
                    TIMESPLUSONE( src.d[s1].c[c1].d[1].c[i],
                                  p.d[s1].c[c1].d[3].c[i] );
                }
        break;
    case -5:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESPLUSONE(  src.d[s1].c[c1].d[0].c[i],
                                   p.d[s1].c[c1].d[0].c[i] );
                    TIMESPLUSONE(  src.d[s1].c[c1].d[1].c[i],
                                   p.d[s1].c[c1].d[1].c[i] );
                    TIMESMINUSONE( src.d[s1].c[c1].d[2].c[i],
                                   p.d[s1].c[c1].d[2].c[i] );
                    TIMESMINUSONE( src.d[s1].c[c1].d[3].c[i],
                                   p.d[s1].c[c1].d[3].c[i] );
                }
        break;
    default:
//       ERR.General(cname, "gr()", "BAD CALL TO gr()\n");
        break;
    }
    return *this;
}

WilsonMatrix& grV(const WilsonMatrix & from, int dir)
{
    int i; /*color*/
    int c1, s1;    /* row indices, color and spin */
    const wilson_matrix & src=from.wmat();

    switch(dir){
    case 0:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESMINUSI( src.d[s1].c[c1].d[3].c[i],
                                 p.d[s1].c[c1].d[0].c[i] );
                    TIMESMINUSI( src.d[s1].c[c1].d[2].c[i],
                                 p.d[s1].c[c1].d[1].c[i] );
                    TIMESPLUSI(  src.d[s1].c[c1].d[1].c[i],
                                 p.d[s1].c[c1].d[2].c[i] );
                    TIMESPLUSI(  src.d[s1].c[c1].d[0].c[i],
                                 p.d[s1].c[c1].d[3].c[i] );
                }
        break;
    case 1:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESMINUSONE( src.d[s1].c[c1].d[3].c[i],
                                   p.d[s1].c[c1].d[0].c[i] );
                    TIMESPLUSONE(  src.d[s1].c[c1].d[2].c[i],
                                   p.d[s1].c[c1].d[1].c[i] );
                    TIMESPLUSONE(  src.d[s1].c[c1].d[1].c[i],
                                   p.d[s1].c[c1].d[2].c[i] );
                    TIMESMINUSONE( src.d[s1].c[c1].d[0].c[i],
                                   p.d[s1].c[c1].d[3].c[i] );
                }
        break;
    case 2:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESMINUSI( src.d[s1].c[c1].d[2].c[i],
                                 p.d[s1].c[c1].d[0].c[i] );
                    TIMESPLUSI(  src.d[s1].c[c1].d[3].c[i],
                                 p.d[s1].c[c1].d[1].c[i] );
                    TIMESPLUSI(  src.d[s1].c[c1].d[0].c[i],
                                 p.d[s1].c[c1].d[2].c[i] );
                    TIMESMINUSI( src.d[s1].c[c1].d[1].c[i],
                                 p.d[s1].c[c1].d[3].c[i] );
                }
        break;
    case 3:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESPLUSONE( src.d[s1].c[c1].d[2].c[i],
                                  p.d[s1].c[c1].d[0].c[i] );
                    TIMESPLUSONE( src.d[s1].c[c1].d[3].c[i],
                                  p.d[s1].c[c1].d[1].c[i] );
                    TIMESPLUSONE( src.d[s1].c[c1].d[0].c[i],
                                  p.d[s1].c[c1].d[2].c[i] );
                    TIMESPLUSONE( src.d[s1].c[c1].d[1].c[i],
                                  p.d[s1].c[c1].d[3].c[i] );
                }
        break;
    case -5:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESPLUSONE(  src.d[s1].c[c1].d[0].c[i],
                                   p.d[s1].c[c1].d[0].c[i] );
                    TIMESPLUSONE(  src.d[s1].c[c1].d[1].c[i],
                                   p.d[s1].c[c1].d[1].c[i] );
                    TIMESMINUSONE( src.d[s1].c[c1].d[2].c[i],
                                   p.d[s1].c[c1].d[2].c[i] );
                    TIMESMINUSONE( src.d[s1].c[c1].d[3].c[i],
                                   p.d[s1].c[c1].d[3].c[i] );
                }
        break;
    default:
        //VRB.Result(cname,fname,"BAD CALL TO gl()\n");
        break;
    }
    return *this;
}

WilsonMatrix& grA(const WilsonMatrix & from, int dir)
{
    int i; /*color*/
    int c1, s1;    /* row indices, color and spin */
    const wilson_matrix & src=from.wmat();

    switch(dir){
    case 0:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESMINUSI( src.d[s1].c[c1].d[3].c[i],
                                 p.d[s1].c[c1].d[0].c[i] );
                    TIMESMINUSI( src.d[s1].c[c1].d[2].c[i],
                                 p.d[s1].c[c1].d[1].c[i] );
                    TIMESMINUSI(  src.d[s1].c[c1].d[1].c[i],
                                  p.d[s1].c[c1].d[2].c[i] );
                    TIMESMINUSI(  src.d[s1].c[c1].d[0].c[i],
                                  p.d[s1].c[c1].d[3].c[i] );
                }
        break;
    case 1:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESMINUSONE( src.d[s1].c[c1].d[3].c[i],
                                   p.d[s1].c[c1].d[0].c[i] );
                    TIMESPLUSONE(  src.d[s1].c[c1].d[2].c[i],
                                   p.d[s1].c[c1].d[1].c[i] );
                    TIMESMINUSONE(  src.d[s1].c[c1].d[1].c[i],
                                    p.d[s1].c[c1].d[2].c[i] );
                    TIMESPLUSONE( src.d[s1].c[c1].d[0].c[i],
                                  p.d[s1].c[c1].d[3].c[i] );
                }
        break;
    case 2:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESMINUSI( src.d[s1].c[c1].d[2].c[i],
                                 p.d[s1].c[c1].d[0].c[i] );
                    TIMESPLUSI(  src.d[s1].c[c1].d[3].c[i],
                                 p.d[s1].c[c1].d[1].c[i] );
                    TIMESMINUSI(  src.d[s1].c[c1].d[0].c[i],
                                  p.d[s1].c[c1].d[2].c[i] );
                    TIMESPLUSI( src.d[s1].c[c1].d[1].c[i],
                                p.d[s1].c[c1].d[3].c[i] );
                }
        break;
    case 3:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESPLUSONE( src.d[s1].c[c1].d[2].c[i],
                                  p.d[s1].c[c1].d[0].c[i] );
                    TIMESPLUSONE( src.d[s1].c[c1].d[3].c[i],
                                  p.d[s1].c[c1].d[1].c[i] );
                    TIMESMINUSONE( src.d[s1].c[c1].d[0].c[i],
                                   p.d[s1].c[c1].d[2].c[i] );
                    TIMESMINUSONE( src.d[s1].c[c1].d[1].c[i],
                                   p.d[s1].c[c1].d[3].c[i] );
                }
        break;
  case PL: // 1-g5
    for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
	  p.d[s1].c[c1].d[0].c[i]=0.0;
	  p.d[s1].c[c1].d[1].c[i]=0.0;
	  TIMESPLUSONE( src.d[s1].c[c1].d[2].c[i],
			p.d[s1].c[c1].d[2].c[i] );
	  TIMESPLUSONE( src.d[s1].c[c1].d[3].c[i],
			p.d[s1].c[c1].d[3].c[i] );
        }
    break;
  case PR: // 1+g5
    for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
	  TIMESPLUSONE(  src.d[s1].c[c1].d[0].c[i],
			 p.d[s1].c[c1].d[0].c[i] );
	  TIMESPLUSONE(  src.d[s1].c[c1].d[1].c[i],
			 p.d[s1].c[c1].d[1].c[i] );
	  p.d[s1].c[c1].d[2].c[i]=0.0;
	  p.d[s1].c[c1].d[3].c[i]=0.0;
        }
    break;
  default:
    //VRB.Result(cname,fname,"BAD CALL TO gl()\n");
    break;
  }
  return *this;
}
#endif

  //! right mult by sigma_mu_nu
  WilsonMatrix& sigmaR(int mu, int nu); 
  //! make a copy of the hermitean conjugate
  WilsonMatrix conj_cp(); 

  //Added by CK
  //! left multiply by 1/2(1-gamma_5)
  WilsonMatrix& glPL();
  //! left multiply by 1/2(1+gamma_5)
  WilsonMatrix& glPR(); 
  //! right multiply by 1/2(1-gamma_5)
  WilsonMatrix& grPL(); 
  //! right multiply by 1/2(1+gamma_5)
  WilsonMatrix& grPR(); 

  //! get a sol. vector
  wilson_vector& sol(int source_spin, int source_color); 

  void load_vec(int sink_spin, int sink_color, const wilson_vector&);
#if 0
    void load_row(int source_spin, int source_color, const wilson_vector&rhs)
    {
        for(int s1 = 0; s1 < 4; ++s1) {
            for(int c1 = 0; c1 < 3; ++c1) {
                p.d[s1].c[c1].d[source_spin].c[source_color]
                    = rhs.d[s1].c[c1];
            }
        }
    }
#endif  //temporarily commented out
  void save_row(int source_spin, int source_color, wilson_vector&);
  void load_row(int source_spin, int source_color, const wilson_vector&);
  void load_elem(int i, int j, int k, int l, Rcomplex elem);

//  Rcomplex Trace();
#ifndef INLINE_WILSON_MATRIX
  Rcomplex Trace() const;
#else
// trace of WilsonMatrix
Rcomplex Trace() const
{
    int c1;
    int s1;
    Rcomplex tr(0.0,0.0);

    for(s1=0;s1<4;++s1){
        for(c1=0;c1<3;++c1){
	    tr+=p.d[s1].c[c1].d[s1].c[c1];
        }
    }
    return tr;
}
#endif
  const wilson_matrix& wmat() const { return p; }
  WilsonMatrix& LeftTimesEqual(const WilsonMatrix& rhs);
  WilsonMatrix& LeftTimesEqual(const Matrix& rhs);
  
  //! Projects positive parity on the sink
  WilsonMatrix& PParProjectSink();
  WilsonMatrix& NParProjectSink();

    //! Projects positive parity on the source
    WilsonMatrix& PParProjectSource();

    //! sink chiral to dirac rotation
    WilsonMatrix& SinkChiralToDirac();
  
  // operator functions

    // another "equal" operator for WilsonMatrix
    WilsonMatrix& operator= (const wilson_matrix& rhs) {
        p=rhs;
        return *this;
    }

#ifndef INLINE_WILSON_MATRIX
    WilsonMatrix& operator= (const Float& rhs);
#else
// equal member operator for WilsonMatrix
WilsonMatrix& operator=(const Float& rhs)
{

    for(int s1=0;s1<4;++s1){
        for(int c1=0;c1<3;++c1){
	    for(int s2=0;s2<4;++s2){
#if 1
	      for(int c2=0;c2<3;++c2){
		    p.d[s1].c[c1].d[s2].c[c2]=Complex(rhs,0.0);
	      }
#else
                for(int c2=0;c2<3;++c2){
		    p.d[s1].c[c1].d[s2].c[c2].real(rhs);
		    p.d[s1].c[c1].d[s2].c[c2].imag(0.0);
                }
#endif
	    }
        }
    }
    return *this;
} 
#endif

#ifndef INLINE_WILSON_MATRIX
    WilsonMatrix& operator+=(const WilsonMatrix& rhs);
    WilsonMatrix& operator-=(const WilsonMatrix& rhs);
    WilsonMatrix& operator*=(const WilsonMatrix& rhs);
#else
// plus-equal member operator for WilsonMatrix
WilsonMatrix& operator+=(const WilsonMatrix& rhs)
{
    int c1, c2;
    int s1, s2;

    for(s1=0;s1<4;++s1){
        for(c1=0;c1<3;++c1){
	    for(s2=0;s2<4;++s2){
                for(c2=0;c2<3;++c2){
                    p.d[s1].c[c1].d[s2].c[c2] += rhs.p.d[s1].c[c1].d[s2].c[c2];
                }
	    }
        }
    }
    return *this;
}

// minus-equal member operator for WilsonMatrix
WilsonMatrix& operator-=(const WilsonMatrix& rhs)
{
    int c1, c2;
    int s1, s2;

    for(s1=0;s1<4;++s1){
        for(c1=0;c1<3;++c1){
	    for(s2=0;s2<4;++s2){
                for(c2=0;c2<3;++c2){
                    p.d[s1].c[c1].d[s2].c[c2] -= rhs.p.d[s1].c[c1].d[s2].c[c2];
                }
	    }
        }
    }
    return *this;
}
    WilsonMatrix& operator*=(const WilsonMatrix& rhs)
{
    wilson_matrix temp=p;
    eq_mult(*this,temp,rhs);
    return *this;
} 
#endif

#ifndef INLINE_WILSON_MATRIX
    WilsonMatrix& operator*=(const Float& rhs);
    WilsonMatrix& operator*=(const Rcomplex& rhs);
#else
// times-equal member operator for WilsonMatrix
WilsonMatrix& operator*=(const Float& rhs)
{

    for(int s1=0;s1<4;++s1){
        for(int c1=0;c1<3;++c1){
	    for(int s2=0;s2<4;++s2){
                for(int c2=0;c2<3;++c2){
		    p.d[s1].c[c1].d[s2].c[c2] *= rhs;
                }
	    }
        }
    }
    return *this;
} 

// times-equal member operator for WilsonMatrix
WilsonMatrix& operator*=(const Rcomplex& rhs)
{
    for(int s1=0;s1<4;++s1){
        for(int c1=0;c1<3;++c1){
	    for(int s2=0;s2<4;++s2){
                for(int c2=0;c2<3;++c2){
		    p.d[s1].c[c1].d[s2].c[c2] *= rhs;
                }
	    }
        }
    }
    return *this;
}
#endif

    friend WilsonMatrix operator*(const WilsonMatrix &wm, const SpinMatrix &sm);
    friend WilsonMatrix operator*(const SpinMatrix &sm, const WilsonMatrix &wm);
    friend WilsonMatrix operator*(const WilsonMatrix &wm, const Matrix &cm);
    friend WilsonMatrix operator*(const Matrix &cm, const WilsonMatrix &wm);

    WilsonMatrix& UMultSource   ( Matrix& U, WilsonMatrix& W);
    WilsonMatrix& UdagMultSource( Matrix& U, WilsonMatrix& W);
  
    // algebra shortcuts

    /*!
      Logically equivalent to   += fact * x , but 
      does not require a temporary WilsonMtarix
    */
    WilsonMatrix& AddMult( const Rcomplex& fact, const WilsonMatrix& x );
  
  // Baryon things
  //CK 2015: For clarity  A.ccl(1) = C^-1 A  A.ccl(-1) = C A
  //                      A.ccr(1) = A C     A.ccr(-1) = A C^-1
  //Poor conventions I know
  WilsonMatrix& ccl(int dir);
  WilsonMatrix& ccr(int dir);
  WilsonMatrix& diq(const WilsonMatrix& rhs);
  WilsonMatrix& joint(const WilsonMatrix& rhs);
  
    friend Rcomplex Trace(const WilsonMatrix& p1, const WilsonMatrix& p2);
  
    SpinMatrix ColorComponent(int row, int col)const {
        SpinMatrix ret;
        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                ret(i, j) = p.d[i].c[row].d[j].c[col];
            }
        }
        return ret;
    }
  //Added by CK
  void ColorTrace(SpinMatrix &into) const{
    //SpinMatrix mapping is i*4+j
    into[0] = p.d[0].c[0].d[0].c[0] + p.d[0].c[1].d[0].c[1] + p.d[0].c[2].d[0].c[2]; //0,0
    into[1] = p.d[0].c[0].d[1].c[0] + p.d[0].c[1].d[1].c[1] + p.d[0].c[2].d[1].c[2]; //0,1
    into[2] = p.d[0].c[0].d[2].c[0] + p.d[0].c[1].d[2].c[1] + p.d[0].c[2].d[2].c[2]; //0,2
    into[3] = p.d[0].c[0].d[3].c[0] + p.d[0].c[1].d[3].c[1] + p.d[0].c[2].d[3].c[2]; //0,3

    into[4] = p.d[1].c[0].d[0].c[0] + p.d[1].c[1].d[0].c[1] + p.d[1].c[2].d[0].c[2]; //1,0
    into[5] = p.d[1].c[0].d[1].c[0] + p.d[1].c[1].d[1].c[1] + p.d[1].c[2].d[1].c[2]; //1,1
    into[6] = p.d[1].c[0].d[2].c[0] + p.d[1].c[1].d[2].c[1] + p.d[1].c[2].d[2].c[2]; //1,2
    into[7] = p.d[1].c[0].d[3].c[0] + p.d[1].c[1].d[3].c[1] + p.d[1].c[2].d[3].c[2]; //1,3

    into[8] = p.d[2].c[0].d[0].c[0] + p.d[2].c[1].d[0].c[1] + p.d[2].c[2].d[0].c[2]; //2,0
    into[9] = p.d[2].c[0].d[1].c[0] + p.d[2].c[1].d[1].c[1] + p.d[2].c[2].d[1].c[2]; //2,1
    into[10] = p.d[2].c[0].d[2].c[0] + p.d[2].c[1].d[2].c[1] + p.d[2].c[2].d[2].c[2]; //2,2
    into[11] = p.d[2].c[0].d[3].c[0] + p.d[2].c[1].d[3].c[1] + p.d[2].c[2].d[3].c[2]; //2,3

    into[12] = p.d[3].c[0].d[0].c[0] + p.d[3].c[1].d[0].c[1] + p.d[3].c[2].d[0].c[2]; //3,0
    into[13] = p.d[3].c[0].d[1].c[0] + p.d[3].c[1].d[1].c[1] + p.d[3].c[2].d[1].c[2]; //3,1
    into[14] = p.d[3].c[0].d[2].c[0] + p.d[3].c[1].d[2].c[1] + p.d[3].c[2].d[2].c[2]; //3,2
    into[15] = p.d[3].c[0].d[3].c[0] + p.d[3].c[1].d[3].c[1] + p.d[3].c[2].d[3].c[2]; //3,3   
  }


};

// added by Hantao
inline WilsonMatrixS::WilsonMatrixS(const WilsonMatrix &w) {
    for(int i = 0; i < 144; ++i) {
        int j = i;
        int c2 = j % 3; j /= 3;
        int s2 = j % 4; j /= 4;
        int c1 = j % 3; j /= 3;
        int s1 = j % 4;
        
        a[2*i  ] = std::real(w(s1, c1, s2, c2));
        a[2*i+1] = std::imag(w(s1, c1, s2, c2));
    }
}

inline const WilsonMatrixS &WilsonMatrixS::operator=(const WilsonMatrix &w) {
    for(int i = 0; i < 144; ++i) {
        int j = i;
        int c2 = j % 3; j /= 3;
        int s2 = j % 4; j /= 4;
        int c1 = j % 3; j /= 3;
        int s1 = j % 4;
        
        a[2*i  ] = std::real(w(s1, c1, s2, c2));
        a[2*i+1] = std::imag(w(s1, c1, s2, c2));
    }
    return *this;
}


#ifndef INLINE_WILSON_MATRIX
WilsonMatrix& eq_mult( WilsonMatrix& xmat,
		       const WilsonMatrix& amat,
		       const WilsonMatrix& bmat );
#else
inline void cmad( Rcomplex& x, const Rcomplex& y, const Rcomplex& z ) { x += y * z; }
inline void cmeq( Rcomplex& x, const Rcomplex& y, const Rcomplex& z ) { x = y * z; }
inline WilsonMatrix& eq_mult( WilsonMatrix& xmat,
		       const WilsonMatrix& amat,
		       const WilsonMatrix& bmat )
{
  const Rcomplex* a(amat.ptr());
  const Rcomplex* b(bmat.ptr());
  Rcomplex* xoff(xmat.ptr());
  register Rcomplex const *point;
  for (int i1=0;i1<12;++i1)
    {
      point = b;
      register const Rcomplex& aval(*a);
      cmeq(xoff[0] ,aval, point[0]);
      cmeq(xoff[1] ,aval, point[1]);
      cmeq(xoff[2] ,aval, point[2]);
      cmeq(xoff[3] ,aval, point[3]);
      cmeq(xoff[4] ,aval, point[4]);
      cmeq(xoff[5] ,aval, point[5]);
      cmeq(xoff[6] ,aval, point[6]);
      cmeq(xoff[7] ,aval, point[7]);
      cmeq(xoff[8] ,aval, point[8]);
      cmeq(xoff[9] ,aval, point[9]);
      cmeq(xoff[10],aval, point[10]);
      cmeq(xoff[11],aval, point[11]);
      a++;
      point+=12;
      for (int i3=1;i3<12;++i3)
	{
	  register const Rcomplex& aval(*a);
	  cmad(xoff[0] ,aval, point[0]);
	  cmad(xoff[1] ,aval, point[1]);
	  cmad(xoff[2] ,aval, point[2]);
	  cmad(xoff[3] ,aval, point[3]);
	  cmad(xoff[4] ,aval, point[4]);
	  cmad(xoff[5] ,aval, point[5]);
	  cmad(xoff[6] ,aval, point[6]);
	  cmad(xoff[7] ,aval, point[7]);
	  cmad(xoff[8] ,aval, point[8]);
	  cmad(xoff[9] ,aval, point[9]);
	  cmad(xoff[10],aval, point[10]);
	  cmad(xoff[11],aval, point[11]);
 	  a++;
	  point+=12;
	}
      xoff+=12;
    }
  return xmat;
}
#endif

// some proto-types for functions that operate on WilsonMatrices
//#ifdef _TARTAN
#if 0

/*! C = A * B */
extern "C" void wmatMultwmat(IFloat* C, const IFloat* A, const IFloat* B);

/*! C = Tr (A * B) */
extern "C" void Tracewmatwmat(IFloat* C, const IFloat* A, const IFloat* B);
	
#endif

//! times operator
#ifndef INLINE_WILSON_MATRIX
extern WilsonMatrix operator*(const WilsonMatrix& lhs, const WilsonMatrix& rhs);
extern WilsonMatrix operator*(const Float& num, const WilsonMatrix& mat);
extern WilsonMatrix operator*(const WilsonMatrix& mat, const Float& num);
extern WilsonMatrix operator*(const Rcomplex& num, const WilsonMatrix& mat);
extern WilsonMatrix operator*(const WilsonMatrix& mat, const Rcomplex& num);
#else
static inline WilsonMatrix operator*(const WilsonMatrix& lhs, const WilsonMatrix& rhs)
{
    WilsonMatrix result(lhs);
    return result *= rhs;
}

static inline WilsonMatrix operator*(const Float& num, const WilsonMatrix& mat)
{
    WilsonMatrix result(mat);
    return result *= num;
}

static inline WilsonMatrix operator*(const WilsonMatrix& mat, const Float& num)
{
    WilsonMatrix result(mat);
    return result *= num;
}

static inline WilsonMatrix operator*(const Rcomplex& num, const WilsonMatrix& mat)
{
    WilsonMatrix result(mat);
    return result *= num;
}

static inline WilsonMatrix operator*(const WilsonMatrix& mat, const Rcomplex& num)
{
    WilsonMatrix result(mat);
    return result *= num;
}
#endif

extern WilsonMatrix operator*(const WilsonMatrix &wm, const SpinMatrix &sm);
extern WilsonMatrix operator*(const SpinMatrix &sm, const WilsonMatrix &wm);
extern WilsonMatrix operator*(const WilsonMatrix &wm, const Matrix &cm);
extern WilsonMatrix operator*(const Matrix &cm, const WilsonMatrix &wm);

#ifndef INLINE_WILSON_MATRIX
//! plus operator
extern WilsonMatrix operator+(const WilsonMatrix& lhs, const WilsonMatrix& rhs);
//! minus operator
extern WilsonMatrix operator-(const WilsonMatrix& lhs, const WilsonMatrix& rhs);
#else
inline WilsonMatrix operator+(const WilsonMatrix& lhs, const WilsonMatrix& rhs)
{
    WilsonMatrix result(lhs);
    return result += rhs;
}

inline WilsonMatrix operator-(const WilsonMatrix& lhs, const WilsonMatrix& rhs)
{
    WilsonMatrix result(lhs);
    return result -= rhs;
}
#endif

//! left multiply by gamma_dir
extern void mult_by_gamma_left( int dir, const wilson_matrix& src, 
				wilson_matrix& dest );

//! right multiply by gamma_dir
extern void mult_by_gamma_right(int dir, const wilson_matrix& src, 
				wilson_matrix& dest );

inline Rcomplex Trace(const WilsonMatrix& p1){ return p1.Trace(); }

//! Spin and Color trace of a 2 WilsonMatrices
#ifndef INLINE_WILSON_MATRIX
extern Rcomplex Trace(const WilsonMatrix& p1, const WilsonMatrix& p2);
#else
inline Rcomplex Trace( const WilsonMatrix& amat,
                const WilsonMatrix& bmat )
{
    const Rcomplex* a(amat.ptr());
    const Rcomplex* b(bmat.ptr());
    Rcomplex tr(0,0);
    for (int i1(0);i1<12;i1++)
        {
            for (int i2(0);i2<12;i2++)
                {
                    cmad(tr,*a,b[i2*12+i1]);
                    a++;
                }
        }
    return tr;
}
#endif

#ifndef INLINE_WILSON_MATRIX
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
#else
// Spin trace of WilsonMatrix
inline Matrix SpinTrace(const WilsonMatrix& Wmat)
{
    Matrix tr=0.0;

    for(int s1=0;s1<4;++s1){
        for(int c1=0;c1<3;++c1){
	    for(int c2=0;c2<3;++c2){
	        tr(c1,c2)+=Wmat.wmat().d[s1].c[c1].d[s1].c[c2];
	    }
        }
    }
    return tr;
}
 
// Spin trace of two WilsonMatrices
inline Matrix SpinTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2)
{
    Matrix tr=0.0;

    for(int s1=0;s1<4;++s1){
        for(int c1=0;c1<3;++c1){
	    for(int s2=0;s2<4;++s2){
                for(int c2=0;c2<3;++c2){
                    for(int c3=0;c3<3;++c3){
                        tr(c1,c3)+=Wmat.wmat().d[s1].c[c1].d[s2].c[c2]
                            *Wmat2.wmat().d[s2].c[c2].d[s1].c[c3];
                    }
                }
	    }
        }
    }
    return tr;
}
 
// Spin trace of three WilsonMatrices
inline Matrix SpinTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2,
                 const WilsonMatrix& Wmat3)
{
    Matrix tr=0.0;

    for(int c1=0;c1<3;++c1){
        for(int c4=0;c4<3;++c4){
            for(int s1=0;s1<4;++s1){
                for(int s2=0;s2<4;++s2){
                    for(int c2=0;c2<3;++c2){
                        for(int s3=0;s3<4;++s3){
                            for(int c3=0;c3<3;++c3){
                                tr(c1,c4)+=Wmat.wmat().d[s1].c[c1].d[s2].c[c2]
                                    *Wmat2.wmat().d[s2].c[c2].d[s3].c[c3]
                                    *Wmat3.wmat().d[s3].c[c3].d[s1].c[c4];
                            }
                        }
                    }
                }
	    }
        }
    }
    return tr;
}
// Color trace of a WilsonMatrix
inline SpinMatrix ColorTrace(const WilsonMatrix& Wmat)
{
    SpinMatrix tr=(Float)0.0;

    for(int s1=0;s1<4;++s1){
        for(int c1=0;c1<3;++c1){
            for(int s2=0;s2<4;++s2){
                tr(s1,s2)+=Wmat.wmat().d[s1].c[c1].d[s2].c[c1];
            }
        }
    }
    return tr;
}

// Color trace of two WilsonMatrices
inline SpinMatrix ColorTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2)
{
    SpinMatrix tr=(Float)0.0;

    for(int s3=0;s3<4;++s3){
        for(int s1=0;s1<4;++s1){
            for(int c1=0;c1<3;++c1){
                for(int s2=0;s2<4;++s2){
                    for(int c2=0;c2<3;++c2){
                        tr(s1,s3)+=Wmat.wmat().d[s1].c[c1].d[s2].c[c2]
                            *Wmat2.wmat().d[s2].c[c2].d[s3].c[c1];
                    }
                }
            }
        }
    }
    return tr;
}

// Color trace of three WilsonMatrices
inline SpinMatrix ColorTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2,
                      const WilsonMatrix& Wmat3)
{
    SpinMatrix tr=(Float)0.0;

    for(int s1=0;s1<4;++s1){
        for(int s4=0;s4<4;++s4){
            for(int c1=0;c1<3;++c1){
                for(int s2=0;s2<4;++s2){
                    for(int c2=0;c2<3;++c2){
                        for(int s3=0;s3<4;++s3){
                            for(int c3=0;c3<3;++c3){
                                tr(s1,s4)+=Wmat.wmat().d[s1].c[c1].d[s2].c[c2]
                                    *Wmat2.wmat().d[s2].c[c2].d[s3].c[c3]
                                    *Wmat3.wmat().d[s3].c[c3].d[s4].c[c1];
                            }
                        }
                    }
                }
            }
        }
    }
    return tr;
}

// Trace of 2 (color) Matrices.
inline Rcomplex Tr(const Matrix& a, const Matrix& b)
{
    Rcomplex tr(0.0,0.0);

    for(int c1=0;c1<3;++c1){
        for(int c2=0;c2<3;++c2){
            tr+=a(c1,c2)*b(c2,c1);
        }
    }

    return tr;
}

// Trace of 2 SpinMatrices.
inline Rcomplex Tr(const SpinMatrix& a, const SpinMatrix& b)
{
    Rcomplex tr(0.0,0.0);

    for(int s1=0;s1<4;++s1){
        for(int s2=0;s2<4;++s2){
            tr+=a(s1,s2)*b(s2,s1);
        }
    }

    return tr;
}
#endif


CPS_END_NAMESPACE

#endif
