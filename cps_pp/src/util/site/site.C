#include <config.h>
#include <util/site.h>

CPS_START_NAMESPACE

Site::Site():
  s(0),
  _looping(false)
{
  x[0]=x[1]=x[2]=x[3] = 0 ;
  size[0] =  GJP.XnodeSites() ;
  size[1] =  GJP.YnodeSites() ;
  size[2] =  GJP.ZnodeSites() ;
  size[3] =  GJP.TnodeSites() ;
  setVolumes() ;  
  computeShifts() ;   
}

Site::Site(int site):
  s(site),
  _looping(false)
{
  size[0] =  GJP.XnodeSites() ;
  size[1] =  GJP.YnodeSites() ;
  size[2] =  GJP.ZnodeSites() ;
  size[3] =  GJP.TnodeSites() ;
  setVolumes() ;
  computeVec() ;
  computeShifts() ;   
}

Site::Site(int x0,int x1,int x2, int x3):
  _looping(false)
{
  x[0] = x0 ;
  x[1] = x1 ;
  x[2] = x2 ;
  x[3] = x3 ;
  size[0] =  GJP.XnodeSites() ;
  size[1] =  GJP.YnodeSites() ;
  size[2] =  GJP.ZnodeSites() ;
  size[3] =  GJP.TnodeSites() ;
  setVolumes() ;
  computeSite() ;
  computeShifts() ;   
}

void Site::nextSite(void)
{
  int nu(0) ;
  
  s++ ;
  x[nu]++ ;
  while(x[nu]>=size[nu])
    {
      x[nu]=0 ;
      if(nu<3)
        {
          nu++ ;
          x[nu]++ ;
        }
    }
}

void Site::nextSiteExcept(const int mu)
{ 
  int nu(0) ;
  
  
  if(mu==0)
    {
      nu++ ;
      s+=vol[nu] ;
      x[nu]++ ;
    }
  else
    {
      s++ ;
      x[nu]++ ;
    }
  while(x[nu]>=size[nu])
    {
      x[nu]=0 ;
      if(nu<3)
        {
          nu++ ;
          if(nu==mu)
            {
              s-- ;
              nu++ ;
              s+=vol[nu] ;
              if(nu>3) nu-- ;
            }
          x[nu]++ ;
        }
    }
}


bool Site::LoopsOverNode()
{
  if (_looping)
    {
      int nu(0) ;
  
      s++ ;
      if ( s >= vol[4] )
        {
          // should never be >, check for this?
          _looping = false;
          Begin(); // overkill
          return _looping;
        }
      x[nu]++ ;
      while(x[nu]>=size[nu])
        {
          x[nu]=0 ;
          if(nu<3)
            {
              nu++ ;
              x[nu]++ ;
            }
        }
    }
  else
    {
      Begin();
      _looping = true;
    }
  return _looping;
}

CPS_END_NAMESPACE
