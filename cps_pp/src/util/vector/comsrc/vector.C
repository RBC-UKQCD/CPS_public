#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/comsrc/vector.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: vector.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.8  2002/12/04 17:16:27  zs
//  Merged the new 2^4 RNG into the code.
//  This new RNG is implemented in the LatRanGen class.
//  The following algorithm and utility classes are affected:
//
//  AlgEig                  Fdwf
//  AlgGheatBath            Fstag
//  AlgHmd                  GlobalJobParameter
//  AlgNoise                Lattice
//  AlgPbp                  Matrix
//  AlgThreept              RandomGenerator
//                          Vector
//
//  Revision 1.7  2002/03/11 22:27:13  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.4.2.1  2002/03/08 16:36:52  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.4  2001/08/16 10:50:39  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:40  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:11  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: vector.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/comsrc/vector.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// vector.C
//
// The vector class functions.
//
// For now this is specific to three colors. The constructor
// will exit if the number of colors is not equal to three.
//
// This file contains the definitions of the Vector and Matrix 
// class functions.
//
// Float is defined in vector.h.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/vector.h>
#include<util/random.h>
#include<util/gjp.h>
#include<comms/glb.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
//------------------------------------------------------------------
// The Matrix class.
// For now Matrix is a class of general 3x3 complex matrices.
// If the number of colors is not 3 the constructor will exit.
//------------------------------------------------------------------
//------------------------------------------------------------------

Matrix::Matrix()
{}

//------------------------------------------------------------------
Matrix::Matrix(IFloat c) 
{ *this = c; }


//------------------------------------------------------------------
Matrix::Matrix(const Complex& c) 
{ *this = c; }


//------------------------------------------------------------------
Matrix::Matrix(const Matrix& m) 
{ *this = m; }


//------------------------------------------------------------------
Matrix& Matrix::operator = (IFloat c)
{
  this -> ZeroMatrix();
  u[0] = u[8] = u[16] = c;
  return *this;
}

//------------------------------------------------------------------
Matrix& Matrix::operator = (const Complex& c)
{
  this -> ZeroMatrix();
  ((Complex*)u)[0] = ((Complex*)u)[4] = ((Complex*)u)[8] = c;
  return *this;
}

//------------------------------------------------------------------
// 0, 1, 2
// 3, 4, 5
// 6, 7, 8
void Matrix::Trans(const IFloat *m) 
{
  Complex *dst = (Complex *)u;
  const Complex *src = (const Complex *)m;
  dst[0] = src[0]; dst[1] = src[3]; dst[2] = src[6];
  dst[3] = src[1]; dst[4] = src[4]; dst[5] = src[7];
  dst[6] = src[2]; dst[7] = src[5]; dst[8] = src[8];
}

//------------------------------------------------------------------
void Matrix::UnitMatrix(void)
{
  IFloat *p = (IFloat *)u;

  for(int i = 0; i < 18; ++i) {
    *p++ = 0.0;
  }
  p = (IFloat *)u;
  //u[0].re = 1.0; u[4].re = 1.0; u[8].re = 1.0;
  *p = 1.0;
  *(p+8) = 1.0;
  *(p+16) = 1.0;
}

//------------------------------------------------------------------
void Matrix::ZeroMatrix(void)
{
  IFloat *p = (IFloat *)u;

  for(int i = 0; i < 18; ++i) {
    *p++ = 0.0;
  }
}


//------------------------------------------------------------------
/*
 *  calculate |U~dag U - I|^2
 */
IFloat Matrix::ErrorSU3() const
{
  Matrix tmp1, tmp2;

  tmp1.Dagger((const IFloat *)this);	// tmp1 = U~dag
  tmp2.DotMEqual(tmp1, *this); 
  
  tmp2 -= 1.0;	// tmp2 = U~dag U - I
  
  const IFloat *p = (const IFloat *)&tmp2;
  return dotProduct(p, p, COLORS*COLORS*2);
}


//------------------------------------------------------------------
Complex& Matrix::operator()(int i, int j)
{ return ((Complex*)u)[i*COLORS+j]; }


//------------------------------------------------------------------
const Complex& Matrix::operator()(int i, int j) const
{ return ((Complex*)u)[i*COLORS+j]; }


//------------------------------------------------------------------
Complex Matrix::Tr() const
{ return ((Complex*)u)[0] + ((Complex*)u)[4] + ((Complex*)u)[8]; }

//------------------------------------------------------------------
Complex Matrix::Char6() const
{
  Complex tmp(reChar6((IFloat *)u), imChar6((IFloat *)u)) ;
  return tmp ;
}

//------------------------------------------------------------------
Complex Matrix::Char8() const
{
  Complex tmp(reChar8((IFloat *)u)) ;
  return tmp ;
}

//------------------------------------------------------------------
Complex Matrix::Char10() const
{
  Complex tmp(reChar10((IFloat *)u), imChar10((IFloat *)u)) ;
  return tmp ;
}


//------------------------------------------------------------------
//------------------------------------------------------------------
// The Vector class.
// For now Vector is a class of general 3 component column
// vectors. If the number of colors is not 3 the constructor
// will exit.
// Vector is not automatically normalized when constructed.
//------------------------------------------------------------------
//------------------------------------------------------------------
Vector::Vector()
{}


//------------------------------------------------------------------
// Returns the dot product (v*, v) summed over the whole lattice.
// len is the number of real numbers in the array.
//------------------------------------------------------------------
Float Vector::NormSqGlbSum(int len)
{
  IFloat sum = dotProduct((IFloat *)&v, (IFloat *)&v, len);
  glb_sum((Float *)&sum);
  return Float(sum);
}


//------------------------------------------------------------------
// Returns the real part of the dot product (v, b) 
// summed over all nodes.
// len is the number of real numbers in the array.
//------------------------------------------------------------------
Float Vector::ReDotProductGlbSum(const Vector *b, int len)
{
  IFloat sum = dotProduct((IFloat *)&v, (IFloat *)b, len);
  glb_sum((Float *)&sum);
  return Float(sum);
}

//------------------------------------------------------------------
// Returns the complex part of the dot product (v, b) 
// summed over all nodes.
// len is the number of real numbers in the array.
//------------------------------------------------------------------
Complex Vector::CompDotProductGlbSum(const Vector *b, int len)
{
  IFloat c_r, c_i;
  compDotProduct(&c_r, &c_i, (IFloat *)&v, (IFloat *)b, len);
  glb_sum((Float *)&c_r);
  glb_sum((Float *)&c_i);
  return Complex(c_r,c_i);
}


//------------------------------------------------------------------
// Slice sum norm square a vector perp. to direction dir
void Vector::NormSqArraySliceSum(Float *f_out, const int size, const int dir)
{
  char *cname = "Vector";
  char *fname = "NormSqArraySliceSum";
  Float *v1 = (Float *)smalloc(GJP.VolNodeSites()*sizeof(Float));
  if (v1 == 0)
    ERR.Pointer(cname, fname, "v1");
  VRB.Smalloc(cname, fname, "v1", v1, GJP.VolNodeSites()*sizeof(Float));

  IFloat *vv = (IFloat *)&v;
  for(int j=0; j < GJP.VolNodeSites(); j++, vv += size)
    v1[j] = Float(dotProduct(vv, vv, size));

  SliceArraySum(f_out, v1, dir);

  VRB.Sfree(cname, fname, "v1", v1);
  sfree(v1);
}

//------------------------------------------------------------------
// Slice sum a lattice of Floats to form an array of length L[dir]
// lattice size of virtual grid in direction dir
void Vector::SliceArraySum(Float *sum, const Float *f_in, const int dir)
{
  char *cname = "Vector";
  char *fname = "SliceArraySum";

  // Slice-sum the vector density to make a 1D vector
  if (dir < 0 || dir >= 4)
    ERR.General(cname,fname,"direction out of bounds");

  int nxx[4] = {GJP.Xnodes()*GJP.XnodeSites(),GJP.Ynodes()*GJP.YnodeSites(),
		GJP.Znodes()*GJP.ZnodeSites(),GJP.Tnodes()*GJP.TnodeSites()};
  int nx[4] = {GJP.XnodeSites(),GJP.YnodeSites(),
	       GJP.ZnodeSites(),GJP.TnodeSites()};
  int coord[4] = {GJP.XnodeCoor(),GJP.YnodeCoor(),
		  GJP.ZnodeCoor(),GJP.TnodeCoor()};
  int len2 = 1;
  for(int i=0; i < dir; ++i)
    len2 *= nx[i];

  int len1 = len2 * nx[dir];
  int plane_offset = nx[dir]*coord[dir];

  int j;
  for(j=0; j < nxx[dir]; ++j)
    sum[j] = 0.0;

  for(j=0; j < GJP.VolNodeSites(); ++j)
  {
    int hypsec = j % len1;
    int plane  = hypsec / len2;
    sum[plane+plane_offset] += f_in[j];
  }

  for(j=0; j < nxx[dir]; ++j)
  {
    glb_sum(&(sum[j]));
  }
}


//------------------------------------------------------------------
// 5D Slice sum a lattice of Floats to form an array of length L[dir]
// lattice size of virtual grid in direction dir
void Vector::SliceArraySumFive(Float *sum, const Float *f_in, const int dir)
{
  char *cname = "Vector";
  char *fname = "SliceArraySumFive";

  // Slice-sum the vector density to make a 1D vector
  if (dir < 0 || dir >= 5)
    ERR.General(cname,fname,"direction out of bounds");

  int nxx[5] = {GJP.Xnodes()*GJP.XnodeSites(),GJP.Ynodes()*GJP.YnodeSites(),
		GJP.Znodes()*GJP.ZnodeSites(),GJP.Tnodes()*GJP.TnodeSites(),
                GJP.Snodes()*GJP.SnodeSites()};
  int nx[5] = {GJP.XnodeSites(),GJP.YnodeSites(),
	       GJP.ZnodeSites(),GJP.TnodeSites(),
               GJP.SnodeSites()};
  int coord[5] = {GJP.XnodeCoor(),GJP.YnodeCoor(),
		  GJP.ZnodeCoor(),GJP.TnodeCoor(),
                  GJP.SnodeCoor()};
  int len2 = 1;
  for(int i=0; i < dir; ++i)
    len2 *= nx[i];

  int len1 = len2 * nx[dir];
  int plane_offset = nx[dir]*coord[dir];

  int j;
  for(j=0; j < nxx[dir]; ++j)
    sum[j] = 0.0;

  for(j=0; j < GJP.VolNodeSites()*GJP.SnodeSites(); ++j)
  {
    int hypsec = j % len1;
    int plane  = hypsec / len2;
    sum[plane+plane_offset] += f_in[j];
  }

  for(j=0; j < nxx[dir]; ++j)
  {
    glb_sum(&(sum[j]));
  }
}


//
//  y  -=  U x
//
void uDotXMinus(IFloat* y, const IFloat* u, const IFloat* x)
{
    *y    -= *u      * *x     - *(u+1)  * *(x+1) + *(u+2)  * *(x+2)
	     - *(u+3)  * *(x+3) + *(u+4)  * *(x+4) - *(u+5)  * *(x+5);
    *(y+1)-= *u      * *(x+1) + *(u+1)  * *x     + *(u+2)  * *(x+3)
	     + *(u+3)  * *(x+2) + *(u+4)  * *(x+5) + *(u+5)  * *(x+4);
    *(y+2)-= *(u+6)  * *x     - *(u+7)  * *(x+1) + *(u+8)  * *(x+2)
	     - *(u+9)  * *(x+3) + *(u+10) * *(x+4) - *(u+11) * *(x+5);
    *(y+3)-= *(u+6)  * *(x+1) + *(u+7)  * *x     + *(u+8)  * *(x+3)
	     + *(u+9)  * *(x+2) + *(u+10) * *(x+5) + *(u+11) * *(x+4);
    *(y+4)-= *(u+12) * *x     - *(u+13) * *(x+1) + *(u+14) * *(x+2)
	     - *(u+15) * *(x+3) + *(u+16) * *(x+4) - *(u+17) * *(x+5);
    *(y+5)-= *(u+12) * *(x+1) + *(u+13) * *x     + *(u+14) * *(x+3)
	     + *(u+15) * *(x+2) + *(u+16) * *(x+5) + *(u+17) * *(x+4);
}

//
//  y  +=  U~dag x
 //
void uDagDotXPlus(IFloat* y, const IFloat* u, const IFloat* x)
{
    *y    += *u      * *x     + *(u+1)  * *(x+1) + *(u+6)  * *(x+2)
	     + *(u+7)  * *(x+3) + *(u+12) * *(x+4) + *(u+13) * *(x+5);
    *(y+1)+= *u      * *(x+1) - *(u+1)  * *x     + *(u+6)  * *(x+3)
	     - *(u+7)  * *(x+2) + *(u+12) * *(x+5) - *(u+13) * *(x+4);
    *(y+2)+= *(u+2)  * *x     + *(u+3)  * *(x+1) + *(u+8)  * *(x+2)
	     + *(u+9)  * *(x+3) + *(u+14) * *(x+4) + *(u+15) * *(x+5);
    *(y+3)+= *(u+2)  * *(x+1) - *(u+3)  * *x     + *(u+8)  * *(x+3)
	     - *(u+9)  * *(x+2) + *(u+14) * *(x+5) - *(u+15) * *(x+4);
    *(y+4)+= *(u+4)  * *x     + *(u+5)  * *(x+1) + *(u+10) * *(x+2)
	     + *(u+11) * *(x+3) + *(u+16) * *(x+4) + *(u+17) * *(x+5);
    *(y+5)+= *(u+4)  * *(x+1) - *(u+5)  * *x     + *(u+10) * *(x+3)
	     - *(u+11) * *(x+2) + *(u+16) * *(x+5) - *(u+17) * *(x+4);
}

//
//  y   =  U~dag x
//
void uDagDotXEqual(IFloat* y, const IFloat* u, const IFloat* x)
{
    *y     =  *u      * *x     + *(u+1)  * *(x+1) + *(u+6)  * *(x+2)
	    + *(u+7)  * *(x+3) + *(u+12) * *(x+4) + *(u+13) * *(x+5);
    *(y+1) =  *u      * *(x+1) - *(u+1)  * *x     + *(u+6)  * *(x+3)
	    - *(u+7)  * *(x+2) + *(u+12) * *(x+5) - *(u+13) * *(x+4);
    *(y+2) =  *(u+2)  * *x     + *(u+3)  * *(x+1) + *(u+8)  * *(x+2)
	    + *(u+9)  * *(x+3) + *(u+14) * *(x+4) + *(u+15) * *(x+5);
    *(y+3) =  *(u+2)  * *(x+1) - *(u+3)  * *x     + *(u+8)  * *(x+3)
	    - *(u+9)  * *(x+2) + *(u+14) * *(x+5) - *(u+15) * *(x+4);
    *(y+4) =  *(u+4)  * *x     + *(u+5)  * *(x+1) + *(u+10) * *(x+2)
	    + *(u+11) * *(x+3) + *(u+16) * *(x+4) + *(u+17) * *(x+5);
    *(y+5) =  *(u+4)  * *(x+1) - *(u+5)  * *x     + *(u+10) * *(x+3)
	    - *(u+11) * *(x+2) + *(u+16) * *(x+5) - *(u+17) * *(x+4);
}


CPS_END_NAMESPACE
