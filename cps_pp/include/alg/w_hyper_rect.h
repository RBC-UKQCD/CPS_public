#include<config.h>
CPS_START_NAMESPACE

#ifndef INCLUDED_W_HYPER_RECT
#define INCLUDED_W_HYPER_RECT
//---------------------------------------------------------------------------
// class WspectHyperRectangle  
//     to be used as the on-node part of a hyperplane.
//---------------------------------------------------------------------------
class WspectHyperRectangle : public WspectGinfo 
{
public:
  // CTOR
  WspectHyperRectangle(int hyper_rectangle_direction, 
		       int hyper_rectangle_global_index);  
  // DTOR
  ~WspectHyperRectangle()     {}
  
  // ACCESSORs
  const int * lclMin()   const       {return d_lcl_min;}
  const int * lclMax()   const       {return d_lcl_max;}
  const int * glbMin()   const       {return d_glb_min;}
  const int * glbMax()   const       {return d_glb_max;}
        int   onNode()   const       {return d_is_on_node;}
        int   dir()      const       {return d_dir;}
        int   glbCoord() const       {return d_glb_min[d_dir];}
        int   lclCoord() const       {return d_lcl_min[d_dir];}

private:
  // not implemented
  WspectHyperRectangle(const WspectHyperRectangle &);
  WspectHyperRectangle& operator=(const WspectHyperRectangle &);

  // static data member
  static char *d_class_name;

  // non-static data members
  int          d_dir;  

  int          d_is_on_node;

  int          d_lcl_min[LORENTZs];
  int          d_lcl_max[LORENTZs];

  int          d_glb_min[LORENTZs];
  int          d_glb_max[LORENTZs];
};

#endif // ! _INCLUDED_W_HYPER_RECT

CPS_END_NAMESPACE
