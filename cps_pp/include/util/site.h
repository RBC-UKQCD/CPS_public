//------------------------------------------------------------------
//
// site.h
//
// Header file for the  Site class 
//
// This class provides utilities to help you
// move arround the Lattice
//
// February 2001
//
// Kostas Orginos
//
//------------------------------------------------------------------

#ifndef INCLUDED_SITE_H
#define INCLUDED_SITE_H

#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>

CPS_START_NAMESPACE
class Site 
{
private:

  int s        ; // site index
  int x    [4] ; // site vector x:0 y:1 z:2 t:3 
  int size [4] ; // dimensions of local volume
  int vol  [5] ;
  int shift[4] ; // the shift need to convert to physical coordinates 

  bool _looping; // looping over node ?

private:
  
  /*!
    compute x[0 -> 3] given s
  */
  void computeVec()
  {
    int ss,i ;
    x[0] = s%size[0] ;
    for(i=1;i<4;i++){
      ss = s/vol[i] ;
      x[i] = ss%size[i] ;
    }
  }
  
  /*!
    computes s given x[0->3] 
  */
  void computeSite()
  {
    s = x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3])) ;
  }
  
  
  void setVolumes()
  {
      vol[0] = 1 ;
      vol[1] = vol[0]*size[0] ;
      vol[2] = vol[1]*size[1] ;
      vol[3] = vol[2]*size[2] ;
      vol[4] = vol[3]*size[3] ;
    }
  
  /*!
    computes the shifts to take local co-ords into
    global co-ords
  */
  void computeShifts()
    {
      shift[0] = GJP.XnodeCoor()*GJP.XnodeSites() ;
      shift[1] = GJP.YnodeCoor()*GJP.YnodeSites() ;
      shift[2] = GJP.ZnodeCoor()*GJP.ZnodeSites() ;
      shift[3] = GJP.TnodeCoor()*GJP.TnodeSites() ;
    }

public:
  
  Site();
  Site( int site );
  Site( int x0, int x1, int x2, int x3 );

  ~Site() {;}

  /*
    grab the values of everything of interest
    for the current site
  */
  
  int Index() const { return s ; }
  
  int X() const { return x[0] ; }
  int Y() const { return x[1] ; }
  int Z() const { return x[2] ; }
  int T() const { return x[3] ; }

  int physX() const { return x[0] + shift[0] ; }
  int physY() const { return x[1] + shift[1] ; }
  int physZ() const { return x[2] + shift[2] ; }
  int physT() const { return x[3] + shift[3] ; }

  /*!
    ugly, but useful for interfacing with the
    link buffer routines
  */
  int* pos()  { return x; } 

  int Coor    ( const int d )  const { return x[d]          ; }
  int physCoor( const int d )  const { return x[d]+shift[d] ; }

  /*
    node information
  */

  int Size    ( const int d )  const { return size[d]       ; }
  int Vol     ( const int d )  const { return vol[d]        ; }

  int nodeBc(const int mu) const { return GJP.NodeBc(mu) ; }

  /*!
    Returns the index of the neighbor in the plus mu direction
   */
  int plusIndex(const int mu) const 
    {
      const int Xmu (x[mu]+1) ;
      int indx(s+vol[mu]);
      if(Xmu==size[mu])
	indx -= size[mu]*vol[mu] ;
      return indx ;
    }

  /*!
    Returns the index of the neighbor in the minus mu direction
   */
  int minusIndex(const int mu) const 
    {
      const int Xmu(x[mu]-1) ;
      int indx(s-vol[mu]);
      if(Xmu<0)
	indx += size[mu]*vol[mu] ;
      return indx ;
    }
  
  /*
    for looping over a node
  */
  
  void Begin()
  {
    s=0 ;
    x[0]=x[1]=x[2]=x[3] = 0 ; 
  }

  void Begin(int mu)
  {
    x[mu]=0 ; 
    computeSite() ;
  }
  
  bool End()              const { return (s<vol[4])       ; }
  bool End(const int mu ) const { return (x[mu]<size[mu]) ; }
  
  bool resetEnd(const int mu )
  { 
    bool E(x[mu]<size[mu]) ;
    if(E) 
      return E;
    else
      Begin(mu) ;
    return E ; 
  }

  void nextSite();

  void nextSite( const int mu )
  { 
    x[mu]++ ;
    s+=vol[mu] ;
  }

  void nextSiteExcept(const int mu);


  /*! 
     another way of looping over a node
     syntax:
     
     site x;
     while ( x.LoopsOverNode() )
     {
        // do stuff
     }

     the site index will be zeroed before the start of the
     loop and be zero after the loop finishes
  */

  bool LoopsOverNode();

} ;

#endif //!INCLUDED_SITE_H
CPS_END_NAMESPACE
