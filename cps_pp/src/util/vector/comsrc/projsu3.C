/*
 * File: projsu3.C
 * Author: Thomas Dumitrescu
 * Last Modified: 06/27/2006
 * ------------------------------------------------------------------
 * Implements Matrix::ProjSU3 which projects an arbitrary nonsingular
 * Matrix onto SU(3) by polar decomposition.
 */

/* YA changed:  6/21/2007
   1) double -> Float
      still assuming Float is typedeffed as double,
      for DOUBLE_EQUAL and TEN_TOTHE_SIGFIG to be appropriate
      precision is defiend in include/util/wfm_config.h ?
      This program is anyway for the backward compatibility
      check of Olegs code. An SU(3) projection is already
      in AlgSmear which will be used in production run.
   2) deleted a space between "\" and NL in #define ROTATE,
      to prevent the warning from compiler. It looks it was
      correctly compiled as authour intended, at least on linux gcc.
   3) #define DOUBLE_EQUAL(a, b) can exec zero div.
*/

//#include <float.h>
#include <math.h>
#include <util/vector.h>


#define MAX_JACOBI_SWEEPS 50
#define TEN_TOTHE_SIGFIG 1.0e15      // should equal 10^(<float.h>:DBL_DIG)


#define ROTATE(a, i, j, k, l) g = a[i][j]; h = a[k][l]; a[i][j] = g - s*(h + g*tau); \
                              a[k][l] = h + s*(g - h*tau);
//#define DOUBLE_EQUAL(a, b) ((((a) == 0.0) && ((b) == 0.0)) || (fabs((a)/(b) - 1.0) < 1.0e-8))
#define DOUBLE_EQUAL(a, b) ((((a) == 0.0) && ((b) == 0.0)) || \
 (fabs(a-b) < fabs(b)*1.0e-8))


CPS_START_NAMESPACE


int diagonalize(Float a[][2*COLORS], Float d[], Float v[][2*COLORS]);


/*
 * projects a general 3 by 3 Matrix M onto SU(3)
 *  ------------------------------------------------------------------
 * this method replaces the matrix M by the properly normalized
 * unitary part U of the polar decomposition M = U*H where H is a
 * positive hermitian matrix and the decomposition is unique if M is
 * nonsingular. The principal brach of the cube root is used in the
 * normalization factor (det(U))^(1/3). If M is singular the method
 * terminates without altering M and returns 1. If M is nonsingular
 * the method returns 0 upon sucessful termination. If the Jacobi
 * matrix diagonalization used to obtain H^(-1) from M does not reach
 * machine precision in MAX_JACOBI_SWEEPS a warning message is printed
 * to stdout.
 */ 

int Matrix::ProjSU3(void)
{
  int i, j, k, l, r;
  Float radius, theta;
  Float det[2];
  Float d[2*COLORS];
  Float a[2*COLORS][2*COLORS];
  Float v[2*COLORS][2*COLORS];
  Complex norm;
  Matrix h, temp;
 
  (*this).Det(det);
  if (DOUBLE_EQUAL(det[0], 0.0) && DOUBLE_EQUAL(det[1], 0.0))
    return (1);
  // determine H^2
  temp.Dagger(*this);
  h.DotMEqual(temp, *this);
  // define augmented real problem to diagonalize H^2
  for (i = 0; i < COLORS; i++) {
    for (j = 0; j < COLORS; j++) {
      a[i][j] = h(i,j).real();
      a[i + COLORS][j] = h(i,j).imag();
      a[i][j + COLORS] = - a[i + COLORS][j];
      a[i + COLORS][j + COLORS] = a[i][j];
    }
  }
  if ((i = diagonalize(a, d, v)) == -1) {
    printf("WARNING: cps::Matrix::ProjSU3 did not reach desired precision in ");
    printf("MAX_JACOBI_SWEEPS = %d steps\n", MAX_JACOBI_SWEEPS);
  }
  // extract COLORS compelx eigenvectors from 2*COLORS real eigenvectors
  l = 0;
  for (k = 0; k < COLORS; k++)  {
    for (r = l + 1; r < 2*COLORS; r++) {
      if (!DOUBLE_EQUAL(d[l],d[r]))
	break; 
    }
    // tripple degeneracy
    if ((r - l) == 2*COLORS) {
      temp = 1.0;
      break;
    // double degeneracy
    } else if ((r - l) == (2*COLORS - 2)) {
      for (i = 0; i < COLORS; i++) {
	temp(i, k)=Complex(v[i][l], v[COLORS + i][l]);
      }
      k++;
      i = 0;
      if (l == 0) {
	for (i = 0; i < COLORS; i++) {
	  temp(i, k + 1)=Complex(v[i][r], v[COLORS + i][r]);
	  i = k + 1;
	}
      }
      // construct a third orthonormal eigenvector by taking the cross product
      temp(0,k) = conj(temp(1, k - 1)*temp(2, i) - temp(2, k - 1)*temp(1, i));
      temp(1,k) = conj(temp(2, k - 1)*temp(0, i) - temp(0, k - 1)*temp(2, i));
      temp(2,k) = conj(temp(0, k - 1)*temp(1, i) - temp(1, k - 1)*temp(0, i));
      break;
    // simple eigenvalue
    } else {
      for (i = 0; i < COLORS; i++) {
	temp(i, k)=Complex(v[i][l], v[COLORS + i][l]);
      }
    }
    l = r;
  }
  // determine H^(-1)
  h = 0.0;
  for (k = 0; k < COLORS; k++) {
    d[2*k] = sqrt(d[2*k]);
    for (i = 0; i < COLORS; i++) {
      for (j = 0; j < COLORS; j++) {
	h(i,j) += temp(i,k)*conj(temp(j,k))/d[2*k];
      }
    }
  }
  // determine U and normalize to unit determinant
  temp.DotMEqual(*this, h);
  *this = temp;
  (*this).Det(det);
  radius = pow(det[0]*det[0] + det[1]*det[1], -1.0/6.0);
   theta = -atan2(det[1], det[0])/3.0;
   norm=Complex(radius*cos(theta), radius*sin(theta));
   h = norm;
   temp.DotMEqual(h, *this);
   *this = temp;
  return (0);
}     
	

//------------------ local function implementations ------------------ 


/* diagonalize a real symmetric 6 by 6 matrix (jacobi method)
 * ----------------------------------------------------------
 * a = input matrix, d = output vector of sorted eigenvalues,
 * v = output column matrix of normalized sorted eigenvectors
 * ----------------------------------------------------------
 * returns the number of jacobi sweeps necessary to
 * diagonalize a to machine precision; returns -1 if
 * number of jacobi sweeps exceeds MAX_JACOBI_SWEEPS.
 * For a careful description of this algorithm, see
 * "Numerical Recipies in C: The Art of Scientific
 * Computing".
 */  

int diagonalize(Float a[][2*COLORS], Float d[], Float v[][2*COLORS])
{
  int i, j, k, iq, ip;
  Float tresh, theta, tau, t, sm, s, h, g, c;
  Float b[2*COLORS];
  Float z[2*COLORS];

  for (ip = 0; ip < 2*COLORS; ip++) {
    for (iq = 0; iq < 2*COLORS; iq++) {
      v[ip][iq] = 0.0;
    }
    v[ip][ip] = 1.0;
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }
  for (i = 0; i < MAX_JACOBI_SWEEPS; i++) {
    sm = 0.0;
    for (ip = 0; ip < 2*COLORS - 1; ip++) {
      for (iq = ip + 1; iq < 2*COLORS; iq++) {
	sm += fabs(a[ip][iq]);
      }
    }
    if (sm == 0.0)
      // terminate if off-diagonal elements are 0.0 to machine precision
      break;
    if (i < 4) {
      tresh = 0.2*sm/(4*COLORS*COLORS);
    } else {
      tresh = 0.0;
    }
    for (ip = 0; ip < 2*COLORS - 1; ip++) {
      for (iq = ip + 1; iq < 2*COLORS; iq++) {
	// g controls numerical precision
	g = TEN_TOTHE_SIGFIG*100.0*fabs(a[ip][iq]);
	if ((i > 4) && ((fabs(d[ip]) + g) == fabs(d[ip])) 
	    && ((fabs(d[iq]) + g) == fabs(d[iq]))) {
	  a[ip][iq] = 0.0;
	} else if (fabs(a[ip][iq]) > tresh) {
	  h = d[iq] - d[ip];
	  if ((fabs(h) + g) == fabs(h)) {
	    t = (a[ip][iq])/h;
	  } else {
	    theta = 0.5*h/(a[ip][iq]);
	    t = 1.0/(fabs(theta) + sqrt(1.0 + theta*theta));
	    if (theta < 0.0)
	      t = -t;
	  }
	  c = 1.0/sqrt(1 + t*t);
	  s = t*c;
	  tau = s/(1.0 + c);
	  h = t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq] = 0.0;
	  for (j = 0; j < ip; j++) {
	    ROTATE(a, j, ip, j, iq);
	  }
	  for (j = ip + 1; j < iq; j++) {
	    ROTATE(a, ip, j, j, iq)
	  }
	  for (j = iq + 1; j < 2*COLORS; j++) {
	    ROTATE(a, ip, j, iq, j)
	  }
	  for (j = 0; j < 2*COLORS; j++) {
	    ROTATE(v, j, ip, j, iq)
	  }
	}
      }
    }
    for (ip = 0; ip < 2*COLORS; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  for (i = 0; i < 2*COLORS; i++) {
    g = d[(k = i)];
    for (j = i + 1; j < 2*COLORS; j++) {
      if (d[j] >= g)
	g = d[(k = j)];
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = g;
      for (j = 0; j < 2*COLORS; j++) {
	g = v[j][i];
        v[j][i] = v[j][k];
	v[j][k] = g;
      }
    }
  } 
  return ((i < MAX_JACOBI_SWEEPS) ? (i) : (-1));
}


CPS_END_NAMESPACE
