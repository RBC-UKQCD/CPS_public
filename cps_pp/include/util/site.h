#ifndef INCLUDED_SITE_H
#define INCLUDED_SITE_H
#include <config.h>
#include <util/gjp.h>
CPS_START_NAMESPACE

/*!
  Class for looping over a single node. Also provides conversion between local
  and global coordinates.
  
  It keeps track of the x,y,z, and t positions separately (as oppposed to just
  the index). If you really need high performance, and whatever you're doing
  at each site is very quick, then you can probably get much better
  performance by hard-coding the loop.
*/
class Site 
{
private:
  //! site index
  int s        ; 
  //! site vector x:0 y:1 z:2 t:3 
  int x    [4] ;
  //! dimensions of local volume
  int size [4] ; 

  int vol  [5] ;
  //! the shift need to convert to physical coordinates 
  int shift[4] ; 
  //! looping over node ?
  bool _looping; 

private:
  
  /*!
    compute x[0 -> 3] given s
  */
  void computeVec();
  
  /*!
    computes s given x[0->3] 
  */
  void computeSite()
  {
    s = x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3])) ;
  }
  
  /*! 
    computes size of the 1-,2-,3- etc.. volumes 
  */
  void setVolumes();

  
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

  int Coor    ( int d )  const { return x[d]          ; }
  int physCoor( int d )  const { return x[d]+shift[d] ; }

  /*
    node information
  */

  int Size    ( int d )  const { return size[d]       ; }
  int Vol     ( int d )  const { return vol[d]        ; }

  int nodeBc(const int mu) const { return GJP.NodeBc(mu) ; }

  /*!
    Returns the index of the neighbor in the plus mu direction
   */
  int plusIndex( int mu ) const;


  /*!
    Returns the index of the neighbor in the minus mu direction
   */
  int minusIndex( int mu ) const;

  
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
  bool End(int mu ) const { return (x[mu]<size[mu]) ; }
  
  bool resetEnd( int mu )
  { 
    bool E(x[mu]<size[mu]) ;
    if(E) 
      return E;
    else
      Begin(mu) ;
    return E ; 
  }

  void nextSite();

  void nextSite( int mu )
  { 
    x[mu]++ ;
    s+=vol[mu] ;
  }

  void nextSiteExcept( int mu );


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

CPS_END_NAMESPACE
#endif //!INCLUDED_SITE_H





