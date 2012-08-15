#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-08-15 03:45:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_w_spect/w_momenta.C,v 1.8 2012-08-15 03:45:46 chulwoo Exp $
//  $Id: w_momenta.C,v 1.8 2012-08-15 03:45:46 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: w_momenta.C,v $
//  $Revision: 1.8 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_w_spect/w_momenta.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


 
CPS_END_NAMESPACE
#include <alg/w_all.h>
#include <util/error.h>
#include <util/verbose.h>              // VRB
#include <util/smalloc.h>              // smalloc, sfree
CPS_START_NAMESPACE


CPS_END_NAMESPACE
#include <stdlib.h>                          
CPS_START_NAMESPACE

char * WspectMomenta::d_class_name = "WspectMomenta";


//---------------------------------------------------------------------------
//  math lib for the calculation of non-zero spatial momenta
//---------------------------------------------------------------------------
CPS_END_NAMESPACE
#include <math.h>
CPS_START_NAMESPACE

static IFloat COS(int n, int m)   {                // cos(2 PI  n / m)
  Float pi = 3.14159265358979323846;              // ^^^^^^^^^^^^^^^^
  int sign = 1;                                   // 1 for pos, -1 for neg

  if (n < 0)        {n = -n;}                     // =>  n>=0
  if (m < 0)        {m = -m;}                     // =>  m>=0
  while (n >= m)    {n -= m;}                     // => [0, 2 PI) angle

  n *= 2;                                         // angle = n/m PI  then
  if (n > m)        {n = 2*m - n;}                // => [0, PI]

  n *= 2;                                         // angle = n/m PI/2  then
  if (n > m)        {n = 2*m - n; sign = -1;}     // => [0, PI/2]

  IFloat answer = sign * (n == 0 ? 1 : (n == m ? 0 :
				       cos((double)(n * pi/(2*m)))));
  return answer;  
}


static IFloat SIN(int n, int m)   {                // sin(2 PI  n / m), 
  Float pi = 3.14159265358979323846;              // ^^^^^^^^^^^^^^^^
  int sign = 1;                                   // 1 for pos, -1 for neg

  if (m < 0)        {m = -m; n = -n;}             // =>  m>=0
  while (n < 0)     {n += m;}                     // =>  n>=0
  while (n >= m)    {n -= m;}                     // => [0, 2 PI) angle

  n *= 2;                                         // angle = n/m PI  then
  if (n > m)        {n = 2*m - n; sign = -1;}     // => [0, PI]

  n *= 2;                                         // angle = n/m PI/2  then
  if (n > m)        {n = 2*m - n;}                // => [0, PI/2]

  IFloat answer = sign * (n == 0 ? 0.0 : (n == m ? 1:
					 sin((double)(n * pi /(2*m)))));
  return answer;
}



//---------------------------------------------------------------------------
// WspectMomenta::CTOR
//---------------------------------------------------------------------------
WspectMomenta::WspectMomenta(const WspectHyperRectangle & whr, 
			     const int center2[LORENTZs],
			     int num)
  : d_num_non_zero(num),
    d_data_p(0),
    d_tdir(whr.dir())
{
  // validate the num of non-zero mementa to be among [1, 6]
  //-------------------------------------------------------------------------
  if (d_num_non_zero < 1 || d_num_non_zero > MAX_NUM) {
    d_num_non_zero = 0;
    return;
  }

  // set d_glb_min, d_glb_max, d_size, d_L
  //-------------------------------------------------------------------------
  {
    int dir[LORENTZs-1];                      // three spatial directions
    
    int l, s;                                 // Lorentz, space index
    for (l = 0, s = 0; l < LORENTZs; ++l) {
      if (l != d_tdir) {
	dir[s++] = l;
      }
    }

    d_size = COMPLEXs * d_num_non_zero * sizeof(IFloat);
    const int * low = whr.glbMin();
    const int * upp = whr.glbMax();    
    for (s = 0; s < LORENTZs-1; ++s) {
      for (s = 0; s < LORENTZs-1; ++s) {
        l = dir[s];
        d_glb_min2[s] = 2 * low[l] - center2[l];
        d_glb_max2[s] = 2 * upp[l] - center2[l];
        d_size *= upp[l] - low[l] + 1;
      }
      d_L2 = glb_sites[dir[0]] * 2;
    }
  }
  

  // allocate memory dynamically
  //-------------------------------------------------------------------------
  {
    d_data_p = (Complex *)smalloc(d_size);
    if (!d_data_p)
      ERR.Pointer(d_class_name, ctor_str, empty_str);
    VRB.Smalloc(d_class_name, ctor_str, empty_str, d_data_p, d_size);
  }
   

  // loop over the spatial coordinates  - the order of loop DOES matter.
  //-------------------------------------------------------------------------
  Complex *cp = d_data_p;    

  for (int z = d_glb_min2[2]; z <= d_glb_max2[2]; z+=2) {
    for (int y = d_glb_min2[1]; y <= d_glb_max2[1]; y+=2) {
      for (int x = d_glb_min2[0]; x <= d_glb_max2[0]; x+=2) {
	switch (d_num_non_zero)  {    // no break between cases!
	case MAX_NUM:                 // spatial [2, 2, 2] 
          cp[5]=Complex(COS(2*(x+y+z), d_L2), SIN(2*(x+y+z), d_L2));
	case 5:                       // spatial [1, 1, 1]
          cp[4]=Complex(COS(x+y+z, d_L2), SIN(x+y+z, d_L2));
	case 4:                       // spatial [0, 2, 2] + permutations
          cp[3]=Complex(
			(COS(2*(x+y),d_L2) + COS(2*(y+z),d_L2) + COS(2*(x+z),d_L2))/3.0,
          (SIN(2*(x+y),d_L2) + SIN(2*(y+z),d_L2) + SIN(2*(x+z),d_L2))/3.0);
	case 3:                       // spatial [0, 1, 1] + permutations
          cp[2]=Complex((COS(x+y, d_L2) + COS(y+z, d_L2) + COS(x+z, d_L2))/3.0,
          (SIN(x+y, d_L2) + SIN(y+z, d_L2) + SIN(x+z, d_L2))/3.0);
	case 2:                       // spatial [0, 0, 2] + permutations
          cp[1]=Complex((COS(2*x, d_L2) + COS(2*y, d_L2) + COS(2*z, d_L2))/3.0,
          (SIN(2*x, d_L2) + SIN(2*y, d_L2) + SIN(2*z, d_L2))/3.0);
	case 1:                       // spatial [0, 0, 1] + permutations
          cp[0]=Complex((COS(x, d_L2) + COS(y, d_L2) + COS(z, d_L2))/3.0,
          (SIN(x, d_L2) + SIN(y, d_L2) + SIN(z, d_L2))/3.0);
	}
	cp += d_num_non_zero;
      }
    }
  } 
}






//---------------------------------------------------------------------------
// WspectMomenta::DTOR
//---------------------------------------------------------------------------
WspectMomenta::~WspectMomenta()
{
  VRB.Func(d_class_name, dtor_str);
  VRB.Sfree(d_class_name, dtor_str, empty_str, d_data_p);
  sfree(d_data_p);  
}



//---------------------------------------------------------------------------
// WspectMomenta::operator[](const int lcl_site[]) 
//---------------------------------------------------------------------------
const Complex *
WspectMomenta::operator[](const int lcl_site[LORENTZs]) const 
{
  const Complex * answer = d_data_p;
  if (d_data_p)
    answer += siteOffset(lcl_site, d_tdir) * d_num_non_zero;  
  return answer;  
}




//---------------------------------------------------------------------------
// WspectMomenta::dumpData(..)
//---------------------------------------------------------------------------
void
WspectMomenta::dumpData() const
{
  if (!d_data_p)
    return;
  
  // loop over the spatial coordinates  - the order of loop DOES matter.
  //-------------------------------------------------------------------------
  Complex *cp = d_data_p;    

  for (int z = d_glb_min2[2]; z <= d_glb_max2[2]; z+=2) {
    for (int y = d_glb_min2[1]; y <= d_glb_max2[1]; y+=2) {
      for (int x = d_glb_min2[0]; x <= d_glb_max2[0]; x+=2) {
        printf("momenta at space[%d][%d][%d]:\n", z/2, y/2, x/2);

	switch (d_num_non_zero)  {    // no break between cases!
	case MAX_NUM:   
          printf("\t 222 [%e, %e]\n", cp[5].real(), cp[5].imag());
        case 5:
          printf("\t 111 [%e, %e]\n", cp[4].real(), cp[4].imag());
        case 4:
          printf("\t 022 [%e, %e]\n", cp[3].real(), cp[3].imag());
        case 3:
          printf("\t 011 [%e, %e]\n", cp[2].real(), cp[2].imag());
        case 2:
          printf("\t 002 [%e, %e]\n", cp[1].real(), cp[1].imag());
        case 1:
          printf("\t 001 [%e, %e]\n", cp[0].real(), cp[0].imag());
        }
        cp += d_num_non_zero;
      }
    }
  }
}



CPS_END_NAMESPACE
