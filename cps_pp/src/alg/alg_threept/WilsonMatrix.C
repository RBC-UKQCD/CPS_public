#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:40 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_threept/WilsonMatrix.C,v 1.7 2004-08-18 11:57:40 zs Exp $
//  $Id: WilsonMatrix.C,v 1.7 2004-08-18 11:57:40 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: WilsonMatrix.C,v $
//  $Revision: 1.7 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_threept/WilsonMatrix.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
//
// The class functions for WilsonMatrix 
//
// For now this is specific to three colors. The constructor
// will exit if the number of colors is not equal to three.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <alg/wilson_matrix.h>
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/glb.h>
#include <comms/sysfunc.h>
CPS_START_NAMESPACE
#endif

//------------------------------------------------------------------
//------------------------------------------------------------------
// The WilsonMatrix class member functions.
//------------------------------------------------------------------
//------------------------------------------------------------------


WilsonMatrix::WilsonMatrix()
{}

WilsonMatrix::WilsonMatrix(const wilson_matrix& rhs)
{ p=rhs; };

WilsonMatrix::WilsonMatrix(const WilsonMatrix& rhs)
{ p=rhs.p; };

// copy complex element
WilsonMatrix::WilsonMatrix(int source_spin, int source_color, int sink_spin, int sink_color, const Rcomplex& z)
{
	p.d[source_spin].c[source_color].d[sink_spin].c[sink_color]=z;
	return;
}

// copy wilson vector
WilsonMatrix::WilsonMatrix(int src_spin, int src_color, const wilson_vector& z)
{
	p.d[src_spin].c[src_color]=z;
	return;
}

WilsonMatrix::WilsonMatrix(const Float& rhs)
{

	for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
	  for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
	    p.d[s1].c[c1].d[s2].c[c2].real(rhs);
	    p.d[s1].c[c1].d[s2].c[c2].imag(0.0);
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


void WilsonMatrix::load_vec(int spin, int color, const wilson_vector& rhs)
{
         p.d[spin].c[color]=rhs;
}

// return the propagator
const wilson_matrix& WilsonMatrix::wmat() const
{
	return p;
}

// return a wilson vector
wilson_vector& WilsonMatrix::sol(int src_spin, int src_color)
{
	return p.d[src_spin].c[src_color];
}

// conjugate the propagator
void WilsonMatrix::hconj()
{
        int c1, c2;
        int s1, s2;
	wilson_matrix mat=p;

	for(s1=0;s1<4;s1++)
	  for(c1=0;c1<3;c1++)
	    for(s2=0;s2<4;s2++)
	      for(c2=0;c2<3;c2++){
	        p.d[s2].c[c2].d[s1].c[c1] = conj(mat.d[s1].c[c1].d[s2].c[c2]);
	      }
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

// "equal" operator for WilsonMatrix
WilsonMatrix& WilsonMatrix::operator=(const WilsonMatrix& rhs)
{
	p=rhs.p;
	return *this;
}

// another "equal" operator for WilsonMatrix
WilsonMatrix& WilsonMatrix::operator=(const wilson_matrix& rhs)
{
	p=rhs;
	return *this;
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
 
// times-equal member operator for WilsonMatrix
WilsonMatrix& WilsonMatrix::operator*=(const WilsonMatrix& rhs)
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
		    	temp.d[s1].c[c1].d[s3].c[c3] * rhs.p.d[s3].c[c3].d[s2].c[c2];
	          }
	        }
	      }
	    }
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
		    p.d[s1].c[c1].d[s2].c[c2].real(rhs);
		    p.d[s1].c[c1].d[s2].c[c2].imag(0.0);
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


// trace of product of two WilsonMatrices
Rcomplex Trace(const WilsonMatrix& p1, const WilsonMatrix& p2)
{
        int c1, c2;
        int s1, s2;
	Rcomplex trace(0.0,0.0);

	for(s1=0;s1<4;++s1){
	  for(c1=0;c1<3;++c1){
	    for(s2=0;s2<4;++s2){
	      for(c2=0;c2<3;++c2){
		trace+=p1.p.d[s1].c[c1].d[s2].c[c2] * 
			p2.p.d[s2].c[c2].d[s1].c[c1];
	      }
	    }
	  }
	}
	return trace;
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

// Multiplication by Dirac gamma's
/*
 Chiral basis
 gamma(XUP)    gamma(YUP)    gamma(ZUP)    gamma(TUP)    gamma(FIVE)
 0  0  0  i    0  0  0 -1    0  0  i  0    0  0  1  0    1  0  0  0
 0  0  i  0    0  0  1  0    0  0  0 -i    0  0  0  1    0  1  0  0
 0 -i  0  0    0  1  0  0   -i  0  0  0    1  0  0  0    0  0 -1  0
-i  0  0  0   -1  0  0  0    0  i  0  0    0  1  0  0    0  0  0 -1
*/

#define TIMESPLUSONE(a,b) { b=a; }
#define TIMESMINUSONE(a,b) { b=-a; }
#define TIMESPLUSI(a,b) { b.real(-a.imag()); b.imag(a.real()); }
#define TIMESMINUSI(a,b) { b.real(a.imag()); b.imag(-a.real()); }

// Left mult by gamma_dir
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
    default:
	//VRB.Result(cname,fname,"BAD CALL TO gl()\n");
	break;
  }
	return *this;
}

// Right mult by gamma_dir
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
        //VRB.Result(cname,fname,"BAD CALL TO gl()\n");
        break;
  }
        return *this;
}

// Color trace of a WilsonMatrix
SpinMatrix ColorTrace(const WilsonMatrix& Wmat)
{
        SpinMatrix tr=0.0;

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
        SpinMatrix tr=0.0;

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
        SpinMatrix tr=0.0;

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


CPS_END_NAMESPACE
