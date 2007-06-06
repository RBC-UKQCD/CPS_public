// PAB

#ifndef _GSUM64_EXT_H_
#define _GSUM64_EXT_H_

#include <config.h>

#include <qcdocos.h>
#include <qcdocos/gsum64.h>

CPS_START_NAMESPACE

enum GsumReduceType { SumReduce, MinReduce,  MaxReduce };


template <class Tp> 
class TypeSafeReducer {   // always assume  sizeof(Tp) < sizeof(double) = 8
 public:
  Tp  result;

  inline void convert(double * tgt,  Tp * src)      {  *(Tp*)tgt = *src;       }
  inline void convert(Tp * tgt, const double * src) {  *tgt = *(const Tp*)src; }

  inline Tp Reduce(double ** buf, int size, GsumReduceType type) {
    switch(type) {
    case SumReduce:
      Sum(buf,size);
      break;
    case MinReduce:
      Min(buf,size);
      break;
    case MaxReduce:
      Max(buf,size);
      break;
    }
    return result;
  }

  inline Tp Sum(double ** buf, int size) {
    Tp res, dat;
    convert(&res, buf[0]);
    for(int i=1;i<size;i++)  {
      convert(&dat, buf[i]);
      res += dat;
    }
    return (result = res);
  }

  inline Tp Max(double ** buf, int size) {
    Tp res, dat;
    convert(&res,buf[0]);
    for(int i=1;i<size;i++)  {
      convert(&dat, buf[i]);
      res = (res>dat?res:dat);
    }
    return (result = res);
  }

  inline Tp Min(double ** buf, int size) {
    Tp res, dat;
    convert(&res,buf[0]);
    for(int i=1;i<size;i++)  {
      convert(&dat, buf[i]);
      res = (res<dat?res:dat);
    }
    return (result = res);
  }
};



class Gsum64Ext{
  // this class sums int and unsigned int, with the buffer as type double *
  // it requires that Passthru class always send raw 64 bits from the buffer 

protected:

  //  enum ReduceType { SumReduce, MinReduce,  MaxReduce };

  static const int Ndim = 6;
  static const int Ndir = 12;
  static const int MAX_NODE_DIM=32;

  static int Gnodes[Ndim];
  static int GCoord[Ndim];

  /*
   * Constant definitions
   */

  static const SCUDir  sdir[];

  static const SCUDir  rdir[];

  static const SCUAxis axis[];


  /*
   * No two gsums can be in progress concurrently, so can share these 
   */
  static double MyVal;
  static double *GsumArrayFW;
  static double *GsumArrayBW;

  /*
   * Potentially will want concurrent 1d,2d,3d, 4d, 5d Gsum instances however, so 
   * want to make lookup local to this Gsum64Ext instance.
   */
  double *Lookup[Ndim][MAX_NODE_DIM]; /*Look table for getting right sum order*/

  int ReduceDims[Ndim];
  int SndForward[Ndim];
  int SndBackward[Ndim];


  PassThru PassForward [Ndim];
  PassThru PassBackward[Ndim];

  void Prepare(SCUAxis axis);
  void Comm   (SCUAxis axis);

  //  double ReduceAll(double,ReduceType);
  void ReduceAll(GsumReduceType);

  void Reduce (SCUAxis axis,GsumReduceType);

  enum VAL_TYPE { VAL_DOUBLE, VAL_INT, VAL_UINT }  gsumValType;

public:

  Gsum64Ext() { Init() ; };

  Gsum64Ext(const SCUAxis *axp,int nd) { Init(axp,nd) ; };

  ~Gsum64Ext() {};

  void Init(const SCUAxis *axis_p = axis ,int Nd = 6);

  double Sum(double);
  double Min(double);
  double Max(double);

  int Sum(int);
  int Min(int);
  int Max(int);

  unsigned int Sum(unsigned int);
  unsigned int Min(unsigned int);
  unsigned int Max(unsigned int);

};
CPS_END_NAMESPACE

#endif
