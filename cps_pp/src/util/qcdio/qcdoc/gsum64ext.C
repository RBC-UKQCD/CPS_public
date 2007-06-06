//PAB
//
#include <qalloc.h>
#include <config.h>
#include <util/error.h>
#include <util/gsum64ext.h>

CPS_START_NAMESPACE

void Gsum64Ext::Init(const SCUAxis *axis_p, int Nd)
{
  int ndouble;
  int i;
  int dim;
  int MyCoord;
  int PeerCoord;
  SCUAxis axs;

  if ( GsumArrayFW == NULL  ) { 
      GsumArrayFW = (double *)qalloc(QNONCACHE,sizeof(double) * MAX_NODE_DIM);
      if(!GsumArrayFW) ERR.Pointer("Gsum64Ext", "Init", "GsumArrayFW");
      GsumArrayBW = (double *)qalloc(QNONCACHE,sizeof(double) * MAX_NODE_DIM);
      if(!GsumArrayBW) ERR.Pointer("Gsum64Ext", "Init", "GsumArrayBW");      
  }

  // Fill out machine size info
  Gnodes[SCU_X]=SizeX();
  Gnodes[SCU_Y]=SizeY();
  Gnodes[SCU_Z]=SizeZ();
  Gnodes[SCU_T]=SizeT();
  Gnodes[SCU_S]=SizeS();
  Gnodes[SCU_W]=SizeW();

  // Fill out machine location info
  GCoord[SCU_X]=CoorX();
  GCoord[SCU_Y]=CoorY();
  GCoord[SCU_Z]=CoorZ();
  GCoord[SCU_T]=CoorT();
  GCoord[SCU_S]=CoorS();
  GCoord[SCU_W]=CoorW();

  /*
   * Anything not mentioned in the callers list is
   * _NOT_ summed over.
   */
  for ( dim = 0 ; dim< 6; dim++ ) {
    ReduceDims[axis[dim]] = 0;
  }

  /*
   * Bum args check - revert to default & print message
   */

  if ( Nd > 6 ) { 
    printf("Gsum64Ext::Init() bad args");
    Nd = 4; 
    axis_p = axis;
  }
 
  for ( dim=0 ; dim < Nd; dim ++ ) { 

    axs = axis_p[dim];

    if ( ! nstat_p->SCULocal[axs] ) { 

      ReduceDims[axs] = 1;
      MyCoord = GCoord[axs];

      Lookup[axs][MyCoord] = &MyVal;

      /*
       *   FW = n/2         BW= n/2-1
       *    
       *L               ->
       *F -> ->  ->  -> 
       *B                  <- <- <-
       *  0  1   2   3  4  5  6  7 
       *
       *  Have one local word and 1 send wire per process
       */
      ndouble = Gnodes[axs]/2;
      if (ndouble > 0 ) { 
	PassForward[axs].Init(GsumArrayFW,
			      ndouble,
			      1,
			      rdir[axs],
			      &sdir[axs],
			      1,
			      PassThru::ProcessA);
	
	for ( i=1 ; i <= ndouble ; i ++ ) {
	  PeerCoord = (MyCoord + Gnodes[axs] -i )%Gnodes[axs];
	  Lookup[axs][PeerCoord] = &GsumArrayFW[i-1];
	}
	SndForward[axs] = 1;
      } else {
	SndForward[axs] = 0;
      }

      ndouble = Gnodes[axs]/2 - 1;
      if ( ndouble > 0 ) { 
	PassBackward[axs].Init(GsumArrayBW,
			       ndouble,
			       1,
			       rdir[axs+6],
			       &sdir[axs+6],
			       1,
			       PassThru::ProcessB);
	
	for ( i=1 ; i <= ndouble ; i ++ ) {
	  PeerCoord = (MyCoord + Gnodes[axs] +i )%Gnodes[axs];
	  Lookup[axs][PeerCoord] = &GsumArrayBW[i-1];
	}
	SndBackward[axs] = 1;
      } else { 
	SndBackward[axs] = 0;
      }
    } else {

      ReduceDims[axs] = 0;

    }
  }
}

//double Gsum64Ext::ReduceAll(double value,GsumReduceType type) 
void Gsum64Ext::ReduceAll(GsumReduceType type)
{
  SCUAxis axs;
  SCUAxis axsnext;
  int DoPrep;
  int dim;
  
  //  MyVal = value; /*Copy value into our reduction variable*/

  axsnext = axis[0];
  if ( ReduceDims[axsnext] ) {
    Prepare(axsnext);
  }

  //  unsigned long word;

 /*Is there a syscall for node_stat_p->PhysMachDim?*/
  for ( dim = 0; dim < Ndim; dim ++ ) { 

    axs = axis[dim];

    //    printf("Gsum64Ext: Summing down AppAxis %d",axs);

    DoPrep = 0;
    if ( dim + 1< Ndim ) { /* 
			    * This overlaps the synchronise for consecutive dimensions
			    * Giving a speed up to x,y,z or t,x,y,z or t,x,y,z,s summs
			    */
      axsnext = axis[dim+1];
      if ( ReduceDims[axsnext] ) DoPrep = 1;
    }

    /*If we want a sum over this dimension*/
    if ( ReduceDims[axs] ) { 
      Comm(axs);
    }

    /*If we want a sum over the next dimension prepare it*/
    if ( DoPrep ) {
      Prepare(axsnext);
    }

    if ( ReduceDims[axs] ) { 
      Reduce(axs,type);
    }

  }

  //  return ( MyVal );
  return;
}
void Gsum64Ext::Reduce(SCUAxis axis,GsumReduceType type)
{
  if(! nstat_p->SCULocal[axis]) { // if this dimension non-trivial
    switch(gsumValType) {
    case VAL_DOUBLE:
      {
	TypeSafeReducer<double> reducer;
	reducer.Reduce(Lookup[axis],Gnodes[axis],type);
	reducer.convert(&MyVal,&reducer.result);
	break;
      }
    case VAL_INT:
      {
	TypeSafeReducer<int> reducer;
	reducer.Reduce(Lookup[axis],Gnodes[axis],type);
	reducer.convert(&MyVal,&reducer.result);
	break;
      }
    case VAL_UINT:
      {
	TypeSafeReducer<unsigned int> reducer;
	reducer.Reduce(Lookup[axis],Gnodes[axis],type);
	reducer.convert(&MyVal,&reducer.result);
	break;
      }
    }
  }
  return;

  /*
  double axis_sum;
  double axis_min;
  double axis_max;
  int node;
  int first = 1;

  switch ( type ) {

  case SumReduce:

    if ( ! nstat_p->SCULocal[axis] ) {  // If this dimension is non-trivial
      axis_sum = 0;
      for ( node=0; node<Gnodes[axis]; node++ ) {
	axis_sum += *Lookup[axis][node];
      }
      MyVal = axis_sum;
    }
    break;

  case MinReduce:

    if ( ! nstat_p->SCULocal[axis] ) {  // If this dimension is non-trivial
      for ( node=0; node<Gnodes[axis]; node++ ) {
	if ( first ) {
	  axis_min = *Lookup[axis][node];
	  first = 0;
	} else if ( axis_min > *Lookup[axis][node] ) { 
	  axis_min = *Lookup[axis][node];
	}
      }
      MyVal = axis_min;
    }
    break;

  case MaxReduce:

    if ( ! nstat_p->SCULocal[axis] ) {   // If this dimension is non-trivial
      for ( node=0; node<Gnodes[axis]; node++ ) {
	if ( first ) {
	  axis_max = *Lookup[axis][node];
	  first = 0;
	} else if ( axis_max < *Lookup[axis][node] ) { 
	  axis_max = *Lookup[axis][node];
	}
      }
      MyVal = axis_max;
    }

    break;
  }
  return;
*/
}

void Gsum64Ext::Prepare(SCUAxis axis)
{
  if ( SndForward[axis] )   PassForward [axis].Prepare();
  if ( SndBackward[axis] )  PassBackward[axis].Prepare();

  return;
}
void Gsum64Ext::Comm(SCUAxis axis)
{
  if ( SndForward[axis] )   PassForward [axis].Start(MyVal);
  if ( SndBackward[axis] )  PassBackward[axis].Start(MyVal);

  if ( SndForward[axis] )   PassForward [axis].Complete();
  if ( SndBackward[axis] )  PassBackward[axis].Complete();

  return;
}

// double reducer
double Gsum64Ext::Sum(double value) 
{
  gsumValType = VAL_DOUBLE;
  TypeSafeReducer<double>  rd;
  rd.convert(&MyVal,&value);

  ReduceAll(SumReduce);

  double result;
  rd.convert(&result, &MyVal);
  return result;
}
double Gsum64Ext::Min(double value) 
{
  gsumValType = VAL_DOUBLE;
  TypeSafeReducer<double>  rd;
  rd.convert(&MyVal,&value);

  ReduceAll(MinReduce);

  double result;
  rd.convert(&result, &MyVal);
  return result;
}
double Gsum64Ext::Max(double value) 
{
  gsumValType = VAL_DOUBLE;
  TypeSafeReducer<double>  rd;
  rd.convert(&MyVal,&value);

  ReduceAll(MaxReduce);

  double result;
  rd.convert(&result, &MyVal);
  return result;
}


// int reducer
int Gsum64Ext::Sum(int value) 
{
  gsumValType = VAL_INT;
  TypeSafeReducer<int>  rd;
  rd.convert(&MyVal,&value);

  ReduceAll(SumReduce);

  int result;
  rd.convert(&result, &MyVal);
  return result;
}
int Gsum64Ext::Min(int value) 
{
  gsumValType = VAL_INT;
  TypeSafeReducer<int>  rd;
  rd.convert(&MyVal,&value);

  ReduceAll(MinReduce);

  int result;
  rd.convert(&result, &MyVal);
  return result;
}
int Gsum64Ext::Max(int value) 
{
  gsumValType = VAL_INT;
  TypeSafeReducer<int>  rd;
  rd.convert(&MyVal,&value);

  ReduceAll(MaxReduce);

  int result;
  rd.convert(&result, &MyVal);
  return result;
}


// uint reducer
unsigned int Gsum64Ext::Sum(unsigned int value) 
{
  gsumValType = VAL_UINT;
  TypeSafeReducer<unsigned int>  rd;
  rd.convert(&MyVal,&value);

  ReduceAll(SumReduce);

  unsigned int result;
  rd.convert(&result, &MyVal);
  return result;
}
unsigned int Gsum64Ext::Min(unsigned int value) 
{
  gsumValType = VAL_UINT;
  TypeSafeReducer<unsigned int>  rd;
  rd.convert(&MyVal,&value);

  ReduceAll(MinReduce);

  unsigned int result;
  rd.convert(&result, &MyVal);
  return result;
}
unsigned int Gsum64Ext::Max(unsigned int value) 
{
  gsumValType = VAL_UINT;
  TypeSafeReducer<unsigned int>  rd;
  rd.convert(&MyVal,&value);

  ReduceAll(MaxReduce);

  unsigned int result;
  rd.convert(&result, &MyVal);
  return result;
}




const SCUDir  Gsum64Ext::sdir[] = { SCU_TP, SCU_XP, SCU_YP, SCU_ZP, SCU_SP, SCU_WP,
				 SCU_TM, SCU_XM, SCU_YM, SCU_ZM, SCU_SM, SCU_WM } ;

const SCUDir  Gsum64Ext::rdir[] = { SCU_TM, SCU_XM, SCU_YM, SCU_ZM, SCU_SM, SCU_WM,
				 SCU_TP, SCU_XP, SCU_YP, SCU_ZP, SCU_SP, SCU_WP } ;

const SCUAxis Gsum64Ext::axis[] = { SCU_T,  SCU_X,  SCU_Y,  SCU_Z,  SCU_S,  SCU_W  } ;
 
double Gsum64Ext::MyVal;
double *Gsum64Ext::GsumArrayFW;
double *Gsum64Ext::GsumArrayBW;

int Gsum64Ext::Gnodes[Ndim];
int Gsum64Ext::GCoord[Ndim];

CPS_END_NAMESPACE

