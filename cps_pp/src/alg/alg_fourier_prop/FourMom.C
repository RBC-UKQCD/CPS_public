#include <alg/FourMom.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>

CPS_START_NAMESPACE

//--------------------------------
// Constructors
//--------------------------------

FourMom::FourMom()
{
  int i;
  for (i=0;i<4;++i){ val[i]=0; }
}

FourMom::FourMom(const int a[])
{
  int i;
  for (i=0;i<4;++i) { val[i]=a[i]; }
}

FourMom::FourMom(int a,int b,int c,int d)
{
  val[0]=a;
  val[1]=b;
  val[2]=c;
  val[3]=d;
}

FourMom::FourMom(const FourMom& a)
{
  int i;
  for (i=0;i<4;++i) { val[i]=a.val[i]; }
}

//-----------------------------------
// operators 
//-----------------------------------

FourMom& FourMom::operator=(const FourMom& a)
{
  int i;
  for (i=0;i<4;++i) { val[i]=a.val[i]; }
  return (*this);
}

FourMom& FourMom::operator+=(const FourMom& a)
{
  int i;
  for (i=0;i<4;++i) { val[i]+=a.val[i]; }
  return (*this);
}


FourMom& FourMom::operator-=(const FourMom& a)
{
  int i;
  for (i=0;i<4;++i) { val[i]-=a.val[i]; }
  return (*this);
}


FourMom& FourMom::operator*=(int a)
{
  int i;
  for (i=0;i<4;++i) { val[i]*=a; }
  return (*this);
}



// MomentaList is not needed because I think the STL vector<> makes a
// better container  -- Sam, 11/28/05

MomentaList::MomentaList():
  vals ( 0x0 ),
  _size( 0 ),
  cname( "MomentaList" )
{;}


MomentaList::MomentaList(int size):
  vals (0x0 ),
  _size(size)
{
  cname="MomentaList";
  alloc(size);
}


void MomentaList::alloc( int size )
{
  char* fname="alloc(int)";
  _size=size;
  dealloc();
  vals=(FourMom*) smalloc(_size*sizeof(FourMom));
  if ( vals==0x0) ERR.Pointer(cname,fname,"vals");
  VRB.Smalloc(cname,fname,"vals",vals,_size*sizeof(FourMom));
}

void MomentaList::dealloc()
{
  char* fname="dealloc()";
  if ( vals != 0x0)
    {
      VRB.Sfree(cname, fname, "vals", vals);
      sfree(vals);
    }
}

CPS_END_NAMESPACE


