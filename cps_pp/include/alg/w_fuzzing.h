#include<config.h>
CPS_START_NAMESPACE
/*w_fuzzing.h
 *
 *class WspectFuzzing : public WspectGinfo
 */

#ifndef INCLUDED_WSPECT_FUZZING_H
#define INCLUDED_WSPECT_FUZZING_H

CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/vector.h>        // Float, Complex
CPS_START_NAMESPACE

//forward declarations
class Lattice;
class WspectArg;

//------------------------------------------------------------------
//class WspectFuzzing:public WspectGinfo  declaration
//------------------------------------------------------------------

class WspectFuzzing : public WspectGinfo {
  //---------------
  //private members
  //---------------  
private:
  static char *d_class_name;
  
  Matrix *fuzzedlink_p;  //storage [site][mu][c1][c2][re,img],  mu=0,1,2
  Matrix *fuzzedlink_tmp_p; //temporary buffer
  int d_size;  //lclVolume*3*3*3*2 IFloats

  //fuzzing algorithm parameters
  int fuzzing_on;
  int fuzzing_level;
  Float fuzzing_c; //multiplier
  int fuzzing_hits; //cabbobo hits

  //constants
  const Lattice &d_lat;  
  const int prop_dir;
  const int sliceonly;/*flag indicating whether the buffer store fuzzedlinks
		       *on a single slice(==1) or all slices(==0)
		       */
  
  int slice; //current slice location

  //-----------------
  //public functions
  //------------------
public:  

  WspectFuzzing(Lattice &lat, WspectArg & warg, int singleslice, int fuzzing_c_index); //allocate memory
  ~WspectFuzzing();  //dealocate memory

  Float getFuzzingC() const{return fuzzing_c;}

  const Matrix *GetLink(const int *site, int dir)const ;

  void run(); //run fuzzing on all slices(sliceonly must be zero)
  void run(int wall); //run fuzzing on a single wall
  //for debugging
  void display();
  void displaySU3(Matrix *m);

  //-----------------
  //private functions
  //-------------------
private:
  //stape_s, multiply, projSU3, cabbibo
  Matrix staple_s(Matrix *gf_p,int site[4],int mu);
  void fuzz(Matrix *FU, const Matrix *U, int lclWall, Float c, int flevel, int hit_iter);
  void cabbibo(int hit, int maxhit, Matrix &X, Matrix &T);

};

#endif




CPS_END_NAMESPACE
