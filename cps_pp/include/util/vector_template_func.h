#ifndef VECTOR_TEMPLATE_FUNC_H
#define VECTOR_TEMPLATE_FUNC_H
#include<config.h>
#include <string.h>		/* memcpy */
#include <util/vector.h>
#include <util/time_cps.h>
CPS_START_NAMESPACE

/*!
  If the vectors are real, this function computes the scalar product.
  \param a A vector
  \param b Another vector
  \param len The size of the vectors.
  \return The scalar product of the vectors.
 */
template<typename AFloat,typename BFloat>
IFloat dotProduct(const AFloat *a, const BFloat *b, int len)
{
    IFloat sum = 0.0;
    for(int i = 0; i < len; ++i) {
    	sum += *a++ * *b++;
    }
    return sum;
}

/*!
  \param a The vector to be multiplied
  \param b The real scalar
  \param len The length of the vectors.
  \post \a a is the multiplied vector.
 */
template<typename AFloat>
void vecTimesEquFloat(AFloat *a, IFloat b, int len)
{
    for(int i = 0; i < len; ++i) {
    	*a++ *= b;
    }
}

/*! Multiplication by a real scalar
  \param u The input vector
  \param fb The real scalar
  \param len The number of real numbers in the vectors.
  \post This vector is multiplied by \a fb
*/
template<typename AFloat, typename BFloat>
void vecEqualsVecTimesEquFloat(AFloat *a, BFloat *b, Float c, int len)
{
  for (int i=0; i<len; i++) {
    *a++ = c * *b++;
  }

}

/*!
  \param a A vector to be added to.
  \param b Another vector.
  \param len The length of the vectors.
  \post \a a is the sum vector.
 */
template<typename AFloat, typename BFloat>
void vecAddEquVec(AFloat *a, const BFloat *b, int len)
{
    for(int i = 0; i < len; ++i) {
    	*a++ += *b++;
    }
}

/*!
  \param a A vector to be subtracted from.
  \param b Another vector.
  \param len The length of the vectors.
  \post \a a is the difference vector.
 */
template<typename AFloat, typename BFloat>
void vecMinusEquVec(AFloat *a, const BFloat *b, int len)
{
    for(int i = 0; i < len; ++i) {
    	*a++ -= *b++;
    }
}

/*!
  \param a The resulting vector
  \param b A real scalar
  \param c A vector
  \param d A vector
  \param len The length of the vectors.
 */
template<typename AFloat, typename CFloat, typename DFloat>
void fTimesV1PlusV2(AFloat *a, IFloat b, const CFloat *c,
	const DFloat *d, int len)
{
    for(int i = 0; i < len; ++i) {
    	*a++ = b * *c++ + *d++;
    }
}

/*!
  \param a The resulting vector
  \param b A real scalar
  \param c A vector
  \param d A vector
  \param len The length of the vectors.
 */
template<typename AFloat, typename CFloat, typename DFloat>
void fTimesV1MinusV2(AFloat *a, IFloat b, const CFloat *c,
	const DFloat *d, int len)
{
    for(int i = 0; i < len; ++i) {
    	*a++ = b * *c++ - *d++;
    }
}

/*! The 3x3 complex matrix is assumed to be stored in a linear form
  where the real part of the (i,j) element is at array position [6i+2j]
  and the imaginary part of the (i,j) element is at array position [6i+2j+1].
  \param a The resulting matrix
  \param b A real scalar factor
  \param c A matrix
  \param n This must be 18 in order for this function to do anything meaningful
*/
template<typename AFloat, typename CFloat>
void oneMinusfTimesMatrix(AFloat *a, IFloat b, const CFloat *c, int n)
{
    IFloat *p = a;
    for(int i = 0; i < n; ++i) {
        *p++ = -b * *c++;
    }
    *a += 1.0;    *(a+8) += 1.0;    *(a+16) += 1.0;
}

/*!
  \param a The resulting negated vector
  \param b A vector
  \param len The length of the vectors.
 */
template<typename AFloat>
void vecNegative(AFloat *a, const IFloat *b, int len)
{
    for(int i = 0; i < len; ++i) {
        *a++ = -*b++;
    }
}

/*!
  \param c_r The real part of the scalar product of the vectors.
  \param c_i The imaginary part of the scalar product of the vectors.  
  \param a A complex vector
  \param b Another complex vector
  \param len The length of the vectors
 */
template<typename AFloat, typename BFloat>
void compDotProduct(IFloat *c_r, IFloat *c_i, 
		    const AFloat *a, const BFloat *b, int len)
{
    *c_r = *c_i = 0.0;
    for(int i = 0; i < len; i += 2, a += 2, b += 2) 
    {
      *c_r += *a * *b     + *(a+1) * *(b+1);   // real part
      *c_i += *a * *(b+1) - *(a+1) * *b;       // imag part
    }
}

/*!
  \param re The real part of the scalar complex factor
  \param im The imaginary part of the scalar complex factor 
  \param a The resulting vector
  \param c A vector
  \param d Another vector
  \param len The length of the vectors
 */
template<typename AFloat, typename CFloat, typename DFloat>
void cTimesV1PlusV2(AFloat *a, IFloat re, IFloat im, const CFloat *c,
	const DFloat *d, int len)
{
    for(int i = 0; i < len; i += 2, c += 2) 
    {
      Float c_re = *c; Float c_im = *(c+1);
      *a++ = re * c_re     - im * c_im+ *d++;   // real part
      *a++ = re * c_im + im * c_re   + *d++;   // imag part
    }
}

/*!
  \param re The real part of the scalar complex factor
  \param im The imaginary part of the scalar complex factor 
  \param a The resulting vector
  \param c A vector
  \param d Another vector
  \param len The length of the vectors
 */
template<typename AFloat, typename CFloat, typename DFloat>
void cTimesV1MinusV2(AFloat *a, IFloat re, IFloat im, const CFloat *c,
	const DFloat *d, int len)
{
    for(int i = 0; i < len; i += 2, c += 2) 
    {
      Float c_re = *c; Float c_im = *(c+1);
      *a++ = re * c_re     - im * c_im - *d++;   // real part
      *a++ = re * c_im + im * c_re   - *d++;   // imag part
    }
}

/*! Assign vector to zero.
  \param len The number of real numbers in the vectors.
  \post This vector has the value 0.
*/
template<typename AFloat>
void vecZero(AFloat *a, int len) {

  for (int i=0; i<len; i++)
    *a++ = 0.0;

}

CPS_END_NAMESPACE
#endif
