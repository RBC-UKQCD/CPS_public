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

#undef INLINE_WILSON_MATRIX

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
    // WilsonMatrix(const WilsonMatrix& rhs);
  
    WilsonMatrix(const wilson_matrix& rhs) {
        p = rhs;
    }

  WilsonMatrix(const Float& rhs);

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
  //! mult the prop by gamma_dir on the left, and return the new matrix
    WilsonMatrix glV(int dir)const;

#ifndef INLINE_WILSON_MATRIX
  //! glV another version. result = gamma_dir*from
    WilsonMatrix& glV(const WilsonMatrix & from, int dir);
#else

#define TIMESPLUSONE(a,b) { b=a; }
#define TIMESMINUSONE(a,b) { b=-a; }
#define TIMESPLUSI(a,b) { b.real(-a.imag()); b.imag(a.real()); }
#define TIMESMINUSI(a,b) { b.real(a.imag()); b.imag(-a.real()); }

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
            ERR.General(cname, "glV", "BAD CALL TO glV(), dir = %d\n", dir);
            break;
        }
        return *this;
    }

#undef TIMESPLUSONE
#undef TIMESMINUSONE
#undef TIMESPLUSI
#undef TIMESMINUSI

#endif

    //! mult the prop by gamma_dir*gamma_5 on the left
    WilsonMatrix& grA(const WilsonMatrix & from, int dir);
    //! mult the prop by gamma_dir on the left
    WilsonMatrix& grV(const WilsonMatrix & from, int dir);

  //! mult the prop by gamma_dir on the left
  WilsonMatrix& gr(int dir); 

  //! right mult by sigma_mu_nu
  WilsonMatrix& sigmaR(int mu, int nu); 

  //! make a copy of the hermitean conjugate
  WilsonMatrix conj_cp(); 

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

  Rcomplex Trace();
  const wilson_matrix& wmat() const; // get p 
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

  WilsonMatrix& operator= (const Float& rhs);
  WilsonMatrix& operator+=(const WilsonMatrix& rhs);
  WilsonMatrix& operator-=(const WilsonMatrix& rhs);
  WilsonMatrix& operator*=(const WilsonMatrix& rhs);
  WilsonMatrix& operator*=(const Float& rhs);
  WilsonMatrix& operator*=(const Rcomplex& rhs);
  
  WilsonMatrix& UMultSource   ( Matrix& U, WilsonMatrix& W);
  WilsonMatrix& UdagMultSource( Matrix& U, WilsonMatrix& W);
  
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
  
    SpinMatrix ColorComponent(int row, int col)const {
        SpinMatrix ret;
        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                ret(i, j) = p.d[i].c[row].d[j].c[col];
            }
        }
        return ret;
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


WilsonMatrix& eq_mult( WilsonMatrix& xmat,
		       const WilsonMatrix& amat,
		       const WilsonMatrix& bmat );

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
#else
static inline WilsonMatrix operator*(const WilsonMatrix& lhs, const WilsonMatrix& rhs)
{
    WilsonMatrix result(lhs);
    return result *= rhs;
}
#endif

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


