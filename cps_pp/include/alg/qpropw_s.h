//------------------------------------------------------------------
//qpropw_s.h : single precision qpropw contrctured from qpropw.

#ifndef INCLUDED_QPROPW_S_H
#define INCLUDED_QPROPW_S_H

#include <alg/wilson_matrix.h>
#include <alg/qpropw.h>

CPS_START_NAMESPACE

class QPropWS{

protected:
  //! pointer to 4d prop

	//Changed by bzy, use double precision for the random source (cost tiny memory)
	//! random src pointer if exist
	Float* rsrc;
  
  //! The class name
  char *cname; 
  
private:
  float* prop;
friend class QPropWRand;
public:

	//constructor from QPropW
  QPropWS(QPropW& rhs, bool rand=false); 
  //! copy constructor
  //QPropWS(const QPropWS& rhs); 

  //! Comunicate Wilson Matrices...
  //WilsonMatrix& GetMatrix(const int *, WilsonMatrix&) const;

  virtual Complex rand_src(int i) const;

  /*! This is a better name for the WallWallProp */
  WilsonMatrix WallSinkProp(int t_sink);
  
  //QPropW& operator=(const QPropW& rhs);

  /*! Returns the prop */
  WilsonMatrix operator[](int i) const;
  
  // DESTRUCTORS
  virtual ~QPropWS();
};
  
CPS_END_NAMESPACE

#endif
