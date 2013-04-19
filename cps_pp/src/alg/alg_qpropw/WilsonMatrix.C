//------------------------------------------------------------------
//
//
// The class functions for WilsonMatrix 
//
// For now this is specific to three colors. The constructor
// will exit if the number of colors is not equal to three.
//
//------------------------------------------------------------------

#include <comms/glb.h>

#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <stdio.h>
#include <alg/wilson_matrix.h>

CPS_START_NAMESPACE

const char *WilsonMatrix::cname = "WilsonMatrix";

//------------------------------------------------------------------
//------------------------------------------------------------------
// The WilsonMatrix class member functions.
//------------------------------------------------------------------
//------------------------------------------------------------------

// WilsonMatrix::WilsonMatrix(const WilsonMatrix& rhs)
// { p=rhs.p; }

// copy complex element
WilsonMatrix::WilsonMatrix(int sink_spin, int sink_color, 
			   int source_spin, int source_color, const Rcomplex& z)

{
    p.d[sink_spin].c[sink_color].d[source_spin].c[source_color]=z;
    return;
}

// copy complex element
void WilsonMatrix::Element(int sink_spin, int sink_color, 
			   int source_spin, int source_color, const Rcomplex& z)
{
    p.d[sink_spin].c[sink_color].d[source_spin].c[source_color]=z;
    return;
}

// copy wilson vector
WilsonMatrix::WilsonMatrix(int sink_spin, int sink_color, const wilson_vector& z)
{
    p.d[sink_spin].c[sink_color]=z;
    return;
}

WilsonMatrix::WilsonMatrix(const Float& rhs)
{

	for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
	  for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
	    p.d[s1].c[c1].d[s2].c[c2]=Complex(rhs,0.0);
        }
}

WilsonMatrix::WilsonMatrix(const Rcomplex& rhs)
{
    int c1;
    int s1;

    for(s1=0;s1<4;++s1){
        for(c1=0;c1<3;++c1){
	    p.d[s1].c[c1].d[s1].c[c1] = rhs;
        }
    }
}


void WilsonMatrix::load_vec(int sink_spin, int sink_color, const wilson_vector& rhs)
{
    p.d[sink_spin].c[sink_color]=rhs;
}

void WilsonMatrix::load_row(int source_spin, int source_color, const wilson_vector& rhs)
{
     int c1;
     int s1;

     for(s1=0;s1<4;++s1){
         for(c1=0;c1<3;++c1){
 	    p.d[s1].c[c1].d[source_spin].c[source_color] = rhs.d[s1].c[c1];
         }
     }
}

void WilsonMatrix::save_row(int source_spin, int source_color, wilson_vector& rhs)
{
        int c1;
        int s1;

	for(s1=0;s1<4;++s1){
	  for(c1=0;c1<3;++c1){
//	    p.d[s1].c[c1].d[source_spin].c[source_color] = rhs.d[s1].c[c1];
		rhs.d[s1].c[c1] =  p.d[s1].c[c1].d[source_spin].c[source_color];
	  }
	}
}
void WilsonMatrix::load_elem(int i, int j, int k, int l, Rcomplex elem) {
  p.d[i].c[j].d[k].c[l] = elem;
}

// return the propagator
const wilson_matrix& WilsonMatrix::wmat() const
{
    return p;
}

// return a wilson vector
wilson_vector& WilsonMatrix::sol(int sink_spin, int sink_color)
{
    return p.d[sink_spin].c[sink_color];
}

/*! 
  Hermitian conjugate of the propagator
  \f[ P_{s_1,c_1;s_2,c_2} = P^*_{s_2,c_2;s_1,c_1}\f]
*/
void WilsonMatrix::hconj()
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

void WilsonMatrix::dump()
{
        int c1, c2;
        int s1, s2;
	wilson_matrix mat=p;

	for(s1=0;s1<4;s1++)
	  for(c1=0;c1<3;c1++)
	    for(s2=0;s2<4;s2++)
	      for(c2=0;c2<3;c2++)
		printf("%d %d %d %d %e %e\n",
		       c1,s1,c2,s2,
		       p.d[s2].c[c2].d[s1].c[c1].real(), 
		       p.d[s2].c[c2].d[s1].c[c1].imag()
		       );
	
}
// return the hermitean conjugate 
WilsonMatrix WilsonMatrix::conj_cp()
{
    int c1, c2;
    int s1, s2;
    WilsonMatrix mat;

    for(s1=0;s1<4;s1++)
        for(c1=0;c1<3;c1++)
	    for(s2=0;s2<4;s2++)
                for(c2=0;c2<3;c2++){
                    mat.p.d[s2].c[c2].d[s1].c[c1] = conj(p.d[s1].c[c1].d[s2].c[c2]);
                }
    return mat;
}

// trace of WilsonMatrix
Rcomplex WilsonMatrix::Trace()
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
 
// plus-equal member operator for WilsonMatrix
WilsonMatrix& WilsonMatrix::operator+=(const WilsonMatrix& rhs)
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
WilsonMatrix& WilsonMatrix::operator-=(const WilsonMatrix& rhs)
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


// optimized times-equal member operator for WilsonMatrix
/****************************************************************
   The code below is temporarily isolated for purposes of merging
   with CPS main branch,  01/10/05, Oleg Loktik
-------------------- Quarantine starts --------------------------

#ifdef PARALLEL
 WilsonMatrix& WilsonMatrix::operator*=(const WilsonMatrix& rhs)
{
	wilson_matrix temp=p;

	wmatMultwmat((IFloat*)&p, 
		     (const IFloat*)&temp, 
		     (const IFloat*)&rhs.p);
	return *this;
}

#else
-------------------- Quarantine ends --------------------------*/


//Use C version 
// times-equal member operator for WilsonMatrix
WilsonMatrix& WilsonMatrix::operator*=(const WilsonMatrix& rhs)
{
    wilson_matrix temp=p;
    eq_mult(*this,temp,rhs);
    return *this;
} 

// Part of quarantine #endif

// Left times-equal member for WilsonMatrix
WilsonMatrix& WilsonMatrix::LeftTimesEqual(const WilsonMatrix& rhs)
{
    wilson_matrix temp=p;
    eq_mult(*this,rhs,temp);
    return *this;
} 

// Left times-equal member for WilsonMatrix
/*WilsonMatrix& WilsonMatrix::LeftTimesEqual(const WilsonMatrix& rhs)
{
        int c1, c2, c3;
        int s1, s2, s3;
        wilson_matrix temp=p;

        for(s1=0;s1<4;++s1){
          for(c1=0;c1<3;++c1){
            for(s2=0;s2<4;++s2){
              for(c2=0;c2<3;++c2){
                p.d[s1].c[c1].d[s2].c[c2]=0.0;
                for(s3=0;s3<4;++s3){
                  for(c3=0;c3<3;++c3){
                    p.d[s1].c[c1].d[s2].c[c2]+=
                        rhs.p.d[s1].c[c1].d[s3].c[c3]*
                        temp.d[s3].c[c3].d[s2].c[c2];
                  }
                }
              }
            }
          }
        }
        return *this;
}
*/
WilsonMatrix& WilsonMatrix::LeftTimesEqual(const Matrix& rhs) {
  wilson_matrix temp = p;

  for(int s1=0;s1<4;s1++) for(int c1=0;c1<3;c1++)
        for(int s2=0;s2<4;s2++) for(int c2=0;c2<3;c2++) {
          p.d[s1].c[c1].d[s2].c[c2] = 0.0;
          for(int c3=0;c3<3;c3++)
            p.d[s1].c[c1].d[s2].c[c2] +=
              rhs(c1,c3) * temp.d[s1].c[c3].d[s2].c[c2];
          }
  return *this;
}

// Multiply WilsonVector by su(3) matrix
WilsonVector& WilsonVector::LeftTimesEqual(const Matrix& rhs) {

  for(int s=0;s<4;s++){
    su3_vector temp = d[s];
    for(int c2=0;c2<3;c2++) {
      d[s].c[c2] = 0.0;
      for(int c1=0;c1<3;c1++)
        d[s].c[c2] += rhs(c2,c1) * temp.c[c1];
    }
  }
  return *this;
}

// times-equal member operator for WilsonMatrix
WilsonMatrix& WilsonMatrix::operator*=(const Float& rhs)
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
WilsonMatrix& WilsonMatrix::operator*=(const Rcomplex& rhs)
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

// equal member operator for WilsonMatrix
WilsonMatrix& WilsonMatrix::operator=(const Float& rhs)
{

    for(int s1=0;s1<4;++s1){
        for(int c1=0;c1<3;++c1){
	    for(int s2=0;s2<4;++s2){
	      for(int c2=0;c2<3;++c2){
		    p.d[s1].c[c1].d[s2].c[c2]=Complex(rhs,0.0);
	      }
	    }
        }
    }
    return *this;
} 

WilsonMatrix& WilsonMatrix::AddMult( const Rcomplex&    fact,
                                     const WilsonMatrix&   x )
{
    int s1,c1,s2,c2;
    for( s1=0; s1<4; ++s1 ){
        for( c1=0; c1<3; ++c1 ){
            for( s2=0; s2<4; ++s2 ){
                for( c2=0; c2<3; ++c2 ){

                    p.d[s1].c[c1].d[s2].c[c2] 
                        += x.p.d[s1].c[c1].d[s2].c[c2]*fact;
                }
            }
        }
    }
    return *this;
}

//
// global functions
//-----------------------------------------------------------------------------

WilsonMatrix operator*(const WilsonMatrix& lhs, const WilsonMatrix& rhs)
{
    WilsonMatrix result(lhs);
    return result *= rhs;
}

WilsonMatrix operator*(const Float& num, const WilsonMatrix& mat)
{
    WilsonMatrix result(mat);
    return result *= num;
}

WilsonMatrix operator*(const WilsonMatrix& mat, const Float& num)
{
    WilsonMatrix result(mat);
    return result *= num;
}

WilsonMatrix operator*(const Rcomplex& num, const WilsonMatrix& mat)
{
    WilsonMatrix result(mat);
    return result *= num;
}

WilsonMatrix operator*(const WilsonMatrix& mat, const Rcomplex& num)
{
    WilsonMatrix result(mat);
    return result *= num;
}

WilsonMatrix operator+(const WilsonMatrix& lhs, const WilsonMatrix& rhs)
{
    WilsonMatrix result(lhs);
    return result += rhs;
}

WilsonMatrix operator-(const WilsonMatrix& lhs, const WilsonMatrix& rhs)
{
    WilsonMatrix result(lhs);
    return result -= rhs;
}

inline void cmad( Rcomplex& x, const Rcomplex& y, const Rcomplex& z )
{
    x += Rcomplex(y.real()*z.real() - y.imag()*z.imag(),
                  y.imag()*z.real() + y.real()*z.imag());
    // x.real()+=y.real()*z.real();
    // x.real()-=y.imag()*z.imag();
    // x.imag()+=y.imag()*z.real();
    // x.imag()+=y.real()*z.imag();
}

Rcomplex Trace( const WilsonMatrix& amat,
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

// Spin trace of WilsonMatrix
Matrix SpinTrace(const WilsonMatrix& Wmat)
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
Matrix SpinTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2)
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
Matrix SpinTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2,
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


#define TIMESPLUSONE(a,b) { b=a; }
#define TIMESMINUSONE(a,b) { b=-a; }
//#define TIMESPLUSI(a,b) { b.real(-a.imag()); b.imag(a.real()); }
//#define TIMESMINUSI(a,b) { b.real(a.imag()); b.imag(-a.real()); }
#define TIMESPLUSI(a,b) { b=Complex(-a.imag(),a.real()); }
#define TIMESMINUSI(a,b) { b=Complex(a.imag(),-a.real()); }

/*!
  Left Multiplication by Dirac gamma's
  \verbatim
  Chiral basis
  gamma(XUP)    gamma(YUP)    gamma(ZUP)    gamma(TUP)    gamma(FIVE)
  0  0  0  i    0  0  0 -1    0  0  i  0    0  0  1  0    1  0  0  0
  0  0  i  0    0  0  1  0    0  0  0 -i    0  0  0  1    0  1  0  0
  0 -i  0  0    0  1  0  0   -i  0  0  0    1  0  0  0    0  0 -1  0
  -i  0  0  0   -1  0  0  0    0  i  0  0    0  1  0  0    0  0  0 -1
  \endverbatim
*/

WilsonMatrix& WilsonMatrix::gl(int dir)
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
	ERR.General("WilsonMatrix", "gl()", "BAD CALL TO gl()\n");
	break;
    }
    return *this;
}

WilsonMatrix WilsonMatrix::glR(int dir)const
{
    int i; /*color*/
    int c2,s2;    /* column indices, color and spin */
    wilson_matrix result;

    switch(dir){
    case 0:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    result.d[0].c[i].d[s2].c[c2]=0.0;
                    result.d[1].c[i].d[s2].c[c2]=0.0;
                    TIMESMINUSI( 2*p.d[1].c[i].d[s2].c[c2],
                                 result.d[2].c[i].d[s2].c[c2] );
                    TIMESMINUSI( 2*p.d[0].c[i].d[s2].c[c2],
                                 result.d[3].c[i].d[s2].c[c2] );
                }
        break;
    case 1:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    result.d[0].c[i].d[s2].c[c2]=0.0;
                    result.d[1].c[i].d[s2].c[c2]=0.0;
                    TIMESPLUSONE(  2*p.d[1].c[i].d[s2].c[c2],
                                   result.d[2].c[i].d[s2].c[c2] );
                    TIMESMINUSONE( 2*p.d[0].c[i].d[s2].c[c2],
                                   result.d[3].c[i].d[s2].c[c2] );
                }
        break;
    case 2:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    result.d[0].c[i].d[s2].c[c2]=0.0;
                    result.d[1].c[i].d[s2].c[c2]=0.0;
                    TIMESMINUSI( 2*p.d[0].c[i].d[s2].c[c2],
                                 result.d[2].c[i].d[s2].c[c2] );
                    TIMESPLUSI(  2*p.d[1].c[i].d[s2].c[c2],
                                 result.d[3].c[i].d[s2].c[c2] );
                }
	break;
    case 3:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    result.d[0].c[i].d[s2].c[c2]=0.0;
                    result.d[1].c[i].d[s2].c[c2]=0.0;
                    TIMESPLUSONE( 2*p.d[0].c[i].d[s2].c[c2],
                                  result.d[2].c[i].d[s2].c[c2] );
                    TIMESPLUSONE( 2*p.d[1].c[i].d[s2].c[c2],
                                  result.d[3].c[i].d[s2].c[c2] );
                }
        break;
    default:
	ERR.General(cname, "glR()", "BAD CALL TO glR()\n");
	break;
    }
    return WilsonMatrix(result);
}

WilsonMatrix WilsonMatrix::glL(int dir)const
{
    int i; /*color*/
    int c2,s2;    /* column indices, color and spin */
    wilson_matrix result;

    switch(dir){
    case 0:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    TIMESPLUSI(  2*p.d[3].c[i].d[s2].c[c2],
                                 result.d[0].c[i].d[s2].c[c2] );
                    TIMESPLUSI(  2*p.d[2].c[i].d[s2].c[c2],
                                 result.d[1].c[i].d[s2].c[c2] );
                    result.d[2].c[i].d[s2].c[c2]=0.0;
                    result.d[3].c[i].d[s2].c[c2]=0.0;
                }
        break;
    case 1:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    TIMESMINUSONE( 2*p.d[3].c[i].d[s2].c[c2],
                                   result.d[0].c[i].d[s2].c[c2] );
                    TIMESPLUSONE(  2*p.d[2].c[i].d[s2].c[c2],
                                   result.d[1].c[i].d[s2].c[c2] );
                    result.d[2].c[i].d[s2].c[c2]=0.0;
                    result.d[3].c[i].d[s2].c[c2]=0.0;
                }
        break;
    case 2:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    TIMESPLUSI(  2*p.d[2].c[i].d[s2].c[c2],
                                 result.d[0].c[i].d[s2].c[c2] );
                    TIMESMINUSI( 2*p.d[3].c[i].d[s2].c[c2],
                                 result.d[1].c[i].d[s2].c[c2] );
                    result.d[2].c[i].d[s2].c[c2]=0.0;
                    result.d[3].c[i].d[s2].c[c2]=0.0;
                }
	break;
    case 3:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    TIMESPLUSONE( 2*p.d[2].c[i].d[s2].c[c2],
                                  result.d[0].c[i].d[s2].c[c2] );
                    TIMESPLUSONE( 2*p.d[3].c[i].d[s2].c[c2],
                                  result.d[1].c[i].d[s2].c[c2] );
                    result.d[2].c[i].d[s2].c[c2]=0.0;
                    result.d[3].c[i].d[s2].c[c2]=0.0;
                }
        break;
    default:
	ERR.General(cname, "glL()", "BAD CALL TO glL()\n");
	break;
    }
    return WilsonMatrix(result);
}

//multiply gamma(i) on the left: result = gamma(i)*from
WilsonMatrix& WilsonMatrix::glV(const WilsonMatrix &from, int dir)
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

//multiply gamma(i) on the left and return a new one
WilsonMatrix WilsonMatrix::glV(int dir)const
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
    default:
        ERR.General(cname, "glV", "BAD CALL TO glV(), dir = %d\n", dir);
	break;
    }
    return WilsonMatrix(result);
}

//multiply gamma(i)gamma(5) on the left: result = gamma(i)*gamma(5)*from
WilsonMatrix& WilsonMatrix::glA(const WilsonMatrix & from, int dir)
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
        ERR.General(cname, "glA", "BAD CALL TO glA(), dir = %d\n", dir);
	break;
    }
    return *this;
}
//multiply gamma(i)gamma(5) on the left and return a new one
WilsonMatrix WilsonMatrix::glA(int dir)const
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
        ERR.General(cname, "glA", "BAD CALL TO glA(), dir = %d\n", dir);
	break;
    }
    return WilsonMatrix(result);
}
/*!
  Right Multiplication by Dirac gamma's

  \verbatim
  Chiral basis
  gamma(XUP)    gamma(YUP)    gamma(ZUP)    gamma(TUP)    gamma(FIVE)
  0  0  0  i    0  0  0 -1    0  0  i  0    0  0  1  0    1  0  0  0
  0  0  i  0    0  0  1  0    0  0  0 -i    0  0  0  1    0  1  0  0
  0 -i  0  0    0  1  0  0   -i  0  0  0    1  0  0  0    0  0 -1  0
  -i  0  0  0   -1  0  0  0    0  i  0  0    0  1  0  0    0  0  0 -1
  \endverbatim
*/
WilsonMatrix& WilsonMatrix::gr(int dir)
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
        ERR.General(cname, "gr()", "BAD CALL TO gr()\n");
        break;
    }
    return *this;
}

WilsonMatrix& WilsonMatrix::grV(const WilsonMatrix & from, int dir)
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

WilsonMatrix& WilsonMatrix::grA(const WilsonMatrix & from, int dir)
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

// Color trace of a WilsonMatrix
SpinMatrix ColorTrace(const WilsonMatrix& Wmat)
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
SpinMatrix ColorTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2)
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
SpinMatrix ColorTrace(const WilsonMatrix& Wmat, const WilsonMatrix& Wmat2,
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
Rcomplex Tr(const Matrix& a, const Matrix& b)
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
Rcomplex Tr(const SpinMatrix& a, const SpinMatrix& b)
{
    Rcomplex tr(0.0,0.0);

    for(int s1=0;s1<4;++s1){
        for(int s2=0;s2<4;++s2){
            tr+=a(s1,s2)*b(s2,s1);
        }
    }

    return tr;
}


// Things for baryons



//---------------------------------------------------------------------------

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// This part is for Baryon 2 and 3 point functions
// Functions : ccl
//             ccr 
//
// by S. Sasaki                          
//------------------------------------------------------------------------
// For spectrum of baryon
//------------------------------------------------------------------------
// Multiplication by Charge Conjucation Matrix
/*
  Chiral basis (gamma_x and gamma_z are minus continuum convention)
  CConj          InvCConj     CConj*gamma(FIVE)
  0  1  0  0      0 -1  0  0      0  1  0  0    
  -1  0  0  0      1  0  0  0     -1  0  0  0    
  0  0  0 -1      0  0  0  1      0  0  0  1    
  0  0  1  0      0  0 -1  0      0  0 -1  0    
*/
#ifndef TIMESPLUSONE
#define TIMESPLUSONE(a,b) { b=a; }
#define TIMESMINUSONE(a,b) { b=-a; }
//#define TIMESPLUSI(a,b) { b.real(-a.imag()); b.imag(a.real()); }
//#define TIMESMINUSI(a,b) { b.real(a.imag()); b.imag(-a.real()); }
#define TIMESPLUSI(a,b) { b=Complex(-a.imag(),a.real()); }
#define TIMESMINUSI(a,b) { b=Complex(a.imag(),-a.real()); }
#endif

// Left mult by Charge conjugation matrix C or gamma5*C
WilsonMatrix& WilsonMatrix::ccl(int dir)
{
    int i; /*color*/
    int c2,s2;    /* column indices, color and spin */
    wilson_matrix src=p;

    switch(dir){
    case 1:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    TIMESPLUSONE(  src.d[3].c[i].d[s2].c[c2],
                                   p.d[2].c[i].d[s2].c[c2] );
                    TIMESMINUSONE( src.d[2].c[i].d[s2].c[c2],
                                   p.d[3].c[i].d[s2].c[c2] );
                    TIMESMINUSONE( src.d[1].c[i].d[s2].c[c2],
                                   p.d[0].c[i].d[s2].c[c2] );
                    TIMESPLUSONE(  src.d[0].c[i].d[s2].c[c2],
                                   p.d[1].c[i].d[s2].c[c2] );        	
                }
        break;
    case -1:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    TIMESMINUSONE(  src.d[3].c[i].d[s2].c[c2],
                                    p.d[2].c[i].d[s2].c[c2] );
                    TIMESPLUSONE( src.d[2].c[i].d[s2].c[c2],
                                  p.d[3].c[i].d[s2].c[c2] );
                    TIMESPLUSONE( src.d[1].c[i].d[s2].c[c2],
                                  p.d[0].c[i].d[s2].c[c2] );
                    TIMESMINUSONE(  src.d[0].c[i].d[s2].c[c2],
                                    p.d[1].c[i].d[s2].c[c2] );        	
                }
        break;
    case 5:
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
                    TIMESMINUSONE( src.d[3].c[i].d[s2].c[c2],
                                   p.d[2].c[i].d[s2].c[c2] );
                    TIMESPLUSONE(  src.d[2].c[i].d[s2].c[c2],
                                   p.d[3].c[i].d[s2].c[c2] );
                    TIMESMINUSONE( src.d[1].c[i].d[s2].c[c2],
                                   p.d[0].c[i].d[s2].c[c2] );
                    TIMESPLUSONE(  src.d[0].c[i].d[s2].c[c2],
                                   p.d[1].c[i].d[s2].c[c2] );
                }
        break;
    default:
	//VRB.Result(cname,fname,"BAD CALL TO ccl()\n");
	break;
    }
    return *this;
}
// Left right by Charge conjugation matrix C or gamma5*C
WilsonMatrix& WilsonMatrix::ccr(int dir)
{
    int i; /*color*/
    int c1,s1;    /* column indices, color and spin */
    wilson_matrix src=p;

    switch(dir){
    case 1:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESPLUSONE(  src.d[s1].c[c1].d[3].c[i],
                                   p.d[s1].c[c1].d[2].c[i] );
                    TIMESMINUSONE( src.d[s1].c[c1].d[2].c[i],
                                   p.d[s1].c[c1].d[3].c[i] );
                    TIMESMINUSONE( src.d[s1].c[c1].d[1].c[i],
                                   p.d[s1].c[c1].d[0].c[i] );
                    TIMESPLUSONE(  src.d[s1].c[c1].d[0].c[i],
                                   p.d[s1].c[c1].d[1].c[i] );
                }
        break;
    case -1:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESMINUSONE(  src.d[s1].c[c1].d[3].c[i],
                                    p.d[s1].c[c1].d[2].c[i] );
                    TIMESPLUSONE( src.d[s1].c[c1].d[2].c[i],
                                  p.d[s1].c[c1].d[3].c[i] );
                    TIMESPLUSONE( src.d[s1].c[c1].d[1].c[i],
                                  p.d[s1].c[c1].d[0].c[i] );
                    TIMESMINUSONE(  src.d[s1].c[c1].d[0].c[i],
                                    p.d[s1].c[c1].d[1].c[i] );
                }
        break;
    case 5:
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
                    TIMESMINUSONE( src.d[s1].c[c1].d[3].c[i],
                                   p.d[s1].c[c1].d[2].c[i] );
                    TIMESPLUSONE(  src.d[s1].c[c1].d[2].c[i],
                                   p.d[s1].c[c1].d[3].c[i] );
                    TIMESMINUSONE( src.d[s1].c[c1].d[1].c[i],
                                   p.d[s1].c[c1].d[0].c[i] );
                    TIMESPLUSONE(  src.d[s1].c[c1].d[0].c[i],
                                   p.d[s1].c[c1].d[1].c[i] );
                }
        break;
    default:
	//VRB.Result(cname,fname,"BAD CALL TO ccr()\n");
	break;
    }
    return *this;
}


//------------------------------------------------------------------------
// Diquark state: 
// D(j,a';i,a) = e(a,b,c) * e(a',b',c') * M_u(k,b;i,b') * M_d(k,c;j,c')
WilsonMatrix& WilsonMatrix::diq(const WilsonMatrix& Wmat)
{
    char *fname = "Diquark(WM&, WM&)";	
    char *cname = "Diquark";	
    VRB.Func(cname, fname);

    int eps[3][3][3];
    int weps;

    for(int i=0; i<3; i++ )
        for(int j=0; j<3; j++ )
            for(int k=0; k<3; k++ )
                {
                    eps[i][j][k]=0;
                }
    eps[0][1][2]=1;eps[1][2][0]=1;eps[2][0][1]=1;
    eps[0][2][1]= -1;eps[1][0][2]= -1;eps[2][1][0]= -1;

    wilson_matrix mat=p; // d-quark
	
    for(int s3=0;s3<4;++s3){
        for(int s1=0;s1<4;++s1){
            for(int a3=0;a3<3;++a3){
                for(int a1=0;a1<3;++a1){

                    p.d[s3].c[a3].d[s1].c[a1] = 0.0;

                    for(int b1=0;b1<3;++b1){
                        for(int b3=0;b3<3;++b3){
                            for(int c1=0;c1<3;++c1){
                                for(int c3=0;c3<3;++c3){

                                    weps = eps[a1][b1][c1]*eps[a3][b3][c3];

                                    for(int s2=0;s2<4;++s2){
                                        if (weps == 1){ 
                                            p.d[s3].c[a3].d[s1].c[a1] 
                                                +=   Wmat.wmat().d[s2].c[b1].d[s1].c[b3]  // u-quark
                                                * mat.d[s2].c[c1].d[s3].c[c3];         // d-quark 
                                        }

                                        if (weps == -1){
                                            p.d[s3].c[a3].d[s1].c[a1] 
                                                += - Wmat.wmat().d[s2].c[b1].d[s1].c[b3] // u-quark
                                                * mat.d[s2].c[c1].d[s3].c[c3];        // d-quark 
                                        }
                                    } 
                                }
                            }
                        }
                    }

                }
            }
        }
    }
    return *this;
}

// As far as I can tell this is identical with the operator*= function
// Kostas
//------------------------------------------------------------------------
/* : 
   M_joint(i,a;j,b) = M_backward(i,a;k,c) * M_forward(k,c;j,b)
   This version is modefied on Nov. 16th, 1999
   NOTE: This function is available for both non_transposed and transposed 
   quark propagator with the proper usage.
*/
WilsonMatrix& WilsonMatrix::joint(const WilsonMatrix& Wmat)
{

    wilson_matrix mat=p;

    for(int s1=0;s1<4;++s1){
        for(int a1=0;a1<3;++a1){
            for(int s2=0;s2<4;++s2){
                for(int a2=0;a2<3;++a2){
                    p.d[s1].c[a1].d[s2].c[a2] = 0.0;
                    for(int s3=0;s3<4;++s3){
                        for(int a3=0;a3<3;++a3){			
                            p.d[s1].c[a1].d[s2].c[a2]
                                += mat.d[s1].c[a1].d[s3].c[a3]
                                *  Wmat.wmat().d[s3].c[a3].d[s2].c[a2];		  
                        } /* a2 sum */
                    } /* s2 sum */
                } /* a2 sum */
            } /* s2 sum */
        } /* a1 sum */
    } /* s1 sum */
	 
    return *this;
}
//------------------------------------------------------------------------

/* 
   Multiply a Wilson vector by a gamma matrix
   usage:  gamma( int dir )

   gamma(0)  X
   0  0  0  i
   0  0  i  0
   0 -i  0  0
   -i  0  0  0

   gamma(1)  Y
   0  0  0 -1
   0  0  1  0
   0  1  0  0
   -1  0  0  0

   gamma(2)  Z
   0  0  i  0
   0  0  0 -i
   -i  0  0  0
   0  i  0  0

   gamma(3)  T
   0  0  1  0
   0  0  0  1
   1  0  0  0
   0  1  0  0

   gamma(-5)  gamma_5
   1  0  0  0
   0  1  0  0
   0  0 -1  0
   0  0  0 -1
*/

WilsonVector& WilsonVector::gamma(int dir){
    int i; /*color*/
 
    WilsonVector src=*this;

    switch(dir){
    case 0: //X
	for(i=0;i<3;i++){
	    TIMESPLUSI(  src.d[3].c[i], d[0].c[i] );
	    TIMESPLUSI(  src.d[2].c[i], d[1].c[i] );
	    TIMESMINUSI( src.d[1].c[i], d[2].c[i] );
	    TIMESMINUSI( src.d[0].c[i], d[3].c[i] );
	}
	break;
    case 1: //Y
	for(i=0;i<3;i++){
	    TIMESMINUSONE( src.d[3].c[i], d[0].c[i] );
	    TIMESPLUSONE(  src.d[2].c[i], d[1].c[i] );
	    TIMESPLUSONE(  src.d[1].c[i], d[2].c[i] );
	    TIMESMINUSONE( src.d[0].c[i], d[3].c[i] );
	}
	break;
    case 2: //Z
	for(i=0;i<3;i++){
	    TIMESPLUSI(  src.d[2].c[i], d[0].c[i] );
	    TIMESMINUSI( src.d[3].c[i], d[1].c[i] );
	    TIMESMINUSI( src.d[0].c[i], d[2].c[i] );
	    TIMESPLUSI(  src.d[1].c[i], d[3].c[i] );
	}
	break;
    case 3: //T
	for(i=0;i<3;i++){
	    TIMESPLUSONE( src.d[2].c[i], d[0].c[i] );
	    TIMESPLUSONE( src.d[3].c[i], d[1].c[i] );
	    TIMESPLUSONE( src.d[0].c[i], d[2].c[i] );
	    TIMESPLUSONE( src.d[1].c[i], d[3].c[i] );
	}
	break;
    case -5: //GAMMA_5
	for(i=0;i<3;i++){
	    TIMESPLUSONE(  src.d[0].c[i], d[0].c[i] );
	    TIMESPLUSONE(  src.d[1].c[i], d[1].c[i] );
	    TIMESMINUSONE( src.d[2].c[i], d[2].c[i] );
	    TIMESMINUSONE( src.d[3].c[i], d[3].c[i] );
	}
	break;
    default:
        ERR.General("WilsonVector","gamma(int)","BAD CALL\n");
    }
    return *this ;
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
WilsonVector& WilsonVector::DiracToChiral()
{
    WilsonVector src(*this) ;
    Float norm(0.5*sqrt(2.0)) ;
    for(int cc(0);cc<3;cc++)
        {
            d[0].c[cc] = src.d[0].c[cc] + src.d[2].c[cc] ; d[0].c[cc] *= norm ;
            d[1].c[cc] = src.d[1].c[cc] + src.d[3].c[cc] ; d[1].c[cc] *= norm ;
            d[2].c[cc] = src.d[0].c[cc] - src.d[2].c[cc] ; d[2].c[cc] *= norm ;
            d[3].c[cc] = src.d[1].c[cc] - src.d[3].c[cc] ; d[3].c[cc] *= norm ;
        }
    return *this ;
}

/*!
  Multiplies  by \f$ \frac{1}{2}(1+\gamma_t)\$
  \f[
  V_{s,c}=\sum_{s_1}\left.\frac{1}{2}(1+\gamma_t)\right|_{s,s_1} W_{s_1,c}  
  \f]
*/
WilsonVector& WilsonVector::PParProject()
{
    WilsonVector tmp(*this);
  
    const Float half(0.5);
  
    tmp.gamma(3) ;
    *this *= half ;
    tmp *= half ;
  
    *this += tmp ; // *this = 0.5 (1+gamma_t) *this = 0.5 *this + 0.5 tmp 
  
    return *this ;
}


/*!
  Multilies the WilsonMatrix source indices with the gauge Matrix U
  \f[
  V^{s,c}_{s',c'} = \sum_{c_1}   W^{s,c}_{s',c_1} U_{c_1,c'}
  \f]
  
*/
WilsonMatrix& WilsonMatrix::UMultSource(Matrix& U, WilsonMatrix& W)
{
    Matrix tmp ;
    tmp.Trans(U) ; //uDotXEqual does \sum_t U_(c,t) X_t
    // we want  \sum_t U_(t,c) X_t Thus transpose
    for(int sink_s=0;sink_s<4;sink_s++)
        for(int sink_c=0;sink_c<3;sink_c++)
            for(int source_s=0;source_s<4;source_s++)
                // uDagDotXEqual needs U in CRAM and W in RAM ???
                uDotXEqual((IFloat *)p.d[sink_s].c[sink_c].d[source_s].c,
                           (const IFloat *)(&tmp),
                           (const IFloat *)W.p.d[sink_s].c[sink_c].d[source_s].c) ;
        
    return *this ;
}


/*!
  Multilies the WilsonMatrix source indices with the dagger 
  of the gauge Matrix U
  
  \f[
  V^{s,c}_{s',c'} = \sum_{c_1} W^{s,c}_{s',c_1}   U^{\dagger}_{c_1,c'}
  \f]

*/
WilsonMatrix& WilsonMatrix::UdagMultSource(Matrix& U, WilsonMatrix& W)
{
    Matrix tmp ;
    tmp.Trans(U) ; // uDagDotXEqual does \sum_t U_(c,t) X_t 
    // we want  \sum_t U_(t,c) X_t Thus transpose
    for(int sink_s=0;sink_s<4;sink_s++)
        for(int sink_c=0;sink_c<3;sink_c++)
            for(int source_s=0;source_s<4;source_s++) 
                // uDagDotXEqual needs U in CRAM and W in RAM ???
                uDagDotXEqual((IFloat *)p.d[sink_s].c[sink_c].d[source_s].c,
                              (const IFloat *)(&tmp),
                              (const IFloat *)W.p.d[sink_s].c[sink_c].d[source_s].c) ;
    return *this ;     
}


/*!
  Multiplies the  sink indices by \f$ \frac{1}{2}(1+\gamma_t)\$
  \f[
  V^{s,c}_{s',c'} = \sum_{s_1} \left.\frac{1}{2}(1+\gamma_t)\right|_{s,s_1}
  W^{s_1,c}_{s',c'}   
  \f]

*/
WilsonMatrix& WilsonMatrix::PParProjectSink()
{
    WilsonMatrix tmp(*this);
    const Float half(0.5) ;

    tmp.gl(3) ;
    *this *= half ;
    tmp *= half ;

    *this += tmp ; // *this = 0.5 (1+gamma_t) *this = 0.5 *this + 0.5 tmp 

    return *this ;
}

WilsonMatrix& WilsonMatrix::NParProjectSink()
{
    WilsonMatrix tmp(*this);
    const Float half(0.5) ;

    tmp.gl(3) ;
    *this *= half ;
    tmp *= half ;

    *this -= tmp ; // *this = 0.5 (1-gamma_t) *this = 0.5 *this + 0.5 tmp 

    return *this ;
}

/*!
  Multiplies the  source indices by \f$ \frac{1}{2}(1+\gamma_t)\$
  \f[
  V^{s,c}_{s',c'} = \sum_{s_1} W^{s,c}_{s_1,c'}   
  \left.\frac{1}{2}(1+\gamma_t)\right|_{s_1,s'}
  \f]

*/
WilsonMatrix& WilsonMatrix::PParProjectSource()
{
    WilsonMatrix tmp(*this);
    const Float half(0.5) ;

    tmp.gr(3) ;
    *this *= half ;
    tmp *= half ;

    *this += tmp ; // *this = *this  0.5 (1+gamma_t) = 0.5 *this + 0.5 tmp 

    return *this ;
}

/*!
  Multiplies the sink indices with the rotation matrix that takes
  you from the Chiral to the Dirac basis.
  \f[
  V^{s,c}_{s',c'} = \sum_{d} W^{s,c}_{d,c'} V^\dagger_{d,s'}
  \f]
*/
WilsonMatrix&  WilsonMatrix::SinkChiralToDirac() 
{
    for(int sp(0) ; sp<4 ; sp++) // sink spin 
        for(int cc(0) ; cc<3; cc++) // sink color
            p.d[sp].c[cc].ChiralToDirac() ;

    return *this ;
}

// Right mult by sigma_mu_nu = 1/2 (gamma_mu gamma_nu - gamma_nu gamma_mu)
WilsonMatrix& WilsonMatrix::sigmaR(int mu, int nu)
{
    wilson_matrix src=p;

    int i; 

    if(nu<mu){int k=nu;nu=mu;mu=k;}

    switch(mu){

    case 0:

        switch(nu){

        case 1:
            for(i=0;i<3;i++)for(int s1=0;s1<4;s1++)for(int j=0;j<3;j++){
                        TIMESMINUSI( src.d[s1].c[j].d[0].c[i],
                                     p.d[s1].c[j].d[0].c[i] );
                        TIMESPLUSI(  src.d[s1].c[j].d[1].c[i],
                                     p.d[s1].c[j].d[1].c[i] );
                        TIMESMINUSI( src.d[s1].c[j].d[2].c[i],
                                     p.d[s1].c[j].d[2].c[i] );
                        TIMESPLUSI(  src.d[s1].c[j].d[3].c[i],
                                     p.d[s1].c[j].d[3].c[i] );	
                    }
            break;
        case 2:
            for(i=0;i<3;i++)for(int s1=0;s1<4;s1++)for(int j=0;j<3;j++){
                        TIMESPLUSONE(  src.d[s1].c[j].d[1].c[i],
                                       p.d[s1].c[j].d[0].c[i] );
                        TIMESMINUSONE( src.d[s1].c[j].d[0].c[i],
                                       p.d[s1].c[j].d[1].c[i] );
                        TIMESPLUSONE(  src.d[s1].c[j].d[3].c[i],
                                       p.d[s1].c[j].d[2].c[i] );
                        TIMESMINUSONE( src.d[s1].c[j].d[2].c[i],
                                       p.d[s1].c[j].d[3].c[i] );	
                    }
            break;
        case 3:
            for(i=0;i<3;i++)for(int s1=0;s1<4;s1++)for(int j=0;j<3;j++){
                        TIMESPLUSI( src.d[s1].c[j].d[1].c[i],
                                    p.d[s1].c[j].d[0].c[i] );
                        TIMESPLUSI( src.d[s1].c[j].d[0].c[i],
                                    p.d[s1].c[j].d[1].c[i] );
                        TIMESMINUSI( src.d[s1].c[j].d[3].c[i],
                                     p.d[s1].c[j].d[2].c[i] );
                        TIMESMINUSI( src.d[s1].c[j].d[2].c[i],
                                     p.d[s1].c[j].d[3].c[i] );	
                    }
            break;
        default:
            ERR.General("","sigmaR(int,int)","BAD CALL\n");
            break;
        }
        break;

    case 1:

        switch(nu){

        case 2:

            for(i=0;i<3;i++)for(int s1=0;s1<4;s1++)for(int j=0;j<3;j++){
                        TIMESMINUSI( src.d[s1].c[j].d[1].c[i],
                                     p.d[s1].c[j].d[0].c[i] );
                        TIMESMINUSI( src.d[s1].c[j].d[0].c[i],
                                     p.d[s1].c[j].d[1].c[i] );
                        TIMESMINUSI( src.d[s1].c[j].d[3].c[i],
                                     p.d[s1].c[j].d[2].c[i] );
                        TIMESMINUSI( src.d[s1].c[j].d[2].c[i],
                                     p.d[s1].c[j].d[3].c[i] );
                    }	
            break;
        case 3:
            for(i=0;i<3;i++)for(int s1=0;s1<4;s1++)for(int j=0;j<3;j++){
                        TIMESPLUSONE(  src.d[s1].c[j].d[1].c[i],
                                       p.d[s1].c[j].d[0].c[i] );
                        TIMESMINUSONE( src.d[s1].c[j].d[0].c[i],
                                       p.d[s1].c[j].d[1].c[i] );
                        TIMESMINUSONE( src.d[s1].c[j].d[3].c[i],
                                       p.d[s1].c[j].d[2].c[i] );
                        TIMESPLUSONE(  src.d[s1].c[j].d[2].c[i],
                                       p.d[s1].c[j].d[3].c[i] );	
                    }
            break;
        default:
            ERR.General("","sigmaR(int,int)","BAD CALL\n");
            break;
        }
        break;
    case 2:

        switch(nu){

        case 3:
            for(i=0;i<3;i++)for(int s1=0;s1<4;s1++)for(int j=0;j<3;j++){
                        TIMESPLUSI(  src.d[s1].c[j].d[0].c[i],
                                     p.d[s1].c[j].d[0].c[i] );
                        TIMESMINUSI( src.d[s1].c[j].d[1].c[i],
                                     p.d[s1].c[j].d[1].c[i] );
                        TIMESMINUSI( src.d[s1].c[j].d[2].c[i],
                                     p.d[s1].c[j].d[2].c[i] );
                        TIMESPLUSI(  src.d[s1].c[j].d[3].c[i],
                                     p.d[s1].c[j].d[3].c[i] );	
                    }
            break;
        default:
            ERR.General("","sigmaR(int,int)","BAD CALL\n");
            break;
        }
        break;
    default:
        ERR.General("","sigmaR(int,int)","BAD CALL\n");
        break;
    }
    return *this;
}

CPS_END_NAMESPACE
