#include<config.h>
CPS_START_NAMESPACE
#ifndef INCLUDED_W_MOMENTA
#define INCLUDED_W_MOMENTA

//---------------------------------------------------------------------------
// class WspectMomenta
//     all the momenta here are up to Lorentz transformations.
//     e.g.:  momentum 001 means average over 001, 010 and 100.
// Assumption:
//     same [global] length in all spatial directions.
// Storage order:
//     if prop_dir = t,  Complex[Z][Y][X][momenta]
//                   z,  Complex[T][Y][X][momenta]
//                   y,  Complex[T][Z][X][momenta]
//                   x,  Complex[T][Z][Y][momenta]
//---------------------------------------------------------------------------
class WspectMomenta  : public WspectGinfo 
{
public:
  // only support 1-6 num of nonzero momenta.
  enum {MAX_NUM = MAX_NUM_MOMENTA - 1};  


  // CTOR               -- holds momenta 001, 002, 011, 022, 111, 222 
  //-------------------------------------------------------------------------
  WspectMomenta(const WspectHyperRectangle & whr,
		const int center2[LORENTZs], 
		int num_non_zero_momenta = WspectMomenta::MAX_NUM); 
		

  // DTOR
  //-------------------------------------------------------------------------
  ~WspectMomenta();


  // ACCESSOR operator   -- return momenta 001, 002, 011, 022, 111, 222 ...
  //-------------------------------------------------------------------------
  const Complex * operator[](const int lcl_site[LORENTZs]) const;


  // ACCESSOR operator   -- check whether there is any data
  //-------------------------------------------------------------------------
  operator const void*() const               {return d_data_p;}


  // ACCESSOR function   -- return number of non-zero momenta              
  //-------------------------------------------------------------------------
  int numNonZeroMom() const                    {return d_num_non_zero;}


  // FOR THE PURPOSE OF DEBUGGING 
  //-------------------------------------------------------------------------
  void dumpData() const;  

private:
  static char * d_class_name;  
  Complex *     d_data_p;               // the spacial momenta table
  int           d_size;                 // size of d_data_p in unit of IFloat  
  int           d_num_non_zero;  

  //  the min and max global space coordinates on this node up to a 
  //  translation [where the quark propagator source centers]  [twice]
  //  d_L:  the spatial length over the entire machine. [ twice ]
  //  d_tdir:  the propagation direction
  int           d_glb_min2[LORENTZs-1];  
  int           d_glb_max2[LORENTZs-1];
  int           d_L2;                    
  int           d_tdir;  
};

#endif // ! _INCLUDED_W_MOMENTA

CPS_END_NAMESPACE
