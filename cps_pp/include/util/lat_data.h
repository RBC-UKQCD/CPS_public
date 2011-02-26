#include <config.h>
#include <stdio.h>

#ifndef LAT_DATA_H
#define LAT_DATA_H
#include <util/vector.h>
#include <util/data_types.h>

CPS_START_NAMESPACE
enum { NEW = 0, INITTED};
enum LatDataAlloc { DEFAULT = 0,
       FAST = 1};
class LatData{
  private:
    static char *cname;
    int status;
  protected:
//    static LatDataAlloc DEFAULT_FLAG;
    int size; //number of IFlots per site
    int vol;  //number of sites
    IFloat *data; 
    StrOrdType str_ord;
    int Check ( const LatData &lat_data);
  public:
    LatData(const LatData &lat);
    LatData();
#if 0
    LatData(LatDataAlloc flags, int size, int vol){
      status = NEW;
      Init(flags,size,vol);
    }
    LatData(int size, int vol){
      status = NEW;
      Init(DEFAULT,size,vol);
    }
#endif
    void  Init(LatDataAlloc flags, int size, int vol);
    void  Init(int size, int vol){
      Init(DEFAULT,size,vol);
    }
    int Vol() { return vol;}
    LatData &operator=(const LatData &lat);
    ~LatData();
    const IFloat *Field(int pos=0, int n=0)
      { return (data+pos*size+n);}
    const int Size(){return size*vol;}

    friend void glb_sum_multi_dir(LatData &dat, int dir);

};

class LatVector: public LatData{
  private:
    int n_vec;
    int vec_size;
  public: 
    LatVector(LatDataAlloc flag, int n_vec , int vol )
      { Init(flag,n_vec,vol); }
    LatVector(int n_vec=1 , int vol =0)
      { Init(DEFAULT,n_vec,vol); }
    void Init(LatDataAlloc flag, int n_vec , int vol );
    ~LatVector();
    Vector *Vec(int pos=0, int vec_row=0);
    Float *Field(int pos=0, int vec_row=0,int n=0); 

    void FTimesV1PlusV2(const Float &fb, LatVector *c,
                        LatVector *d){
        fTimesV1PlusV2(data,Float(fb), c->Field(),d->Field(),size*vol);
    }
};

class LatMatrix: virtual public LatData{
  private:
    static char *cname;
    int n_mat;
    int mat_size;
    int Check ( const LatMatrix &lat_data);
  public: 
    LatMatrix();
    LatMatrix(const LatMatrix &lat);
    LatMatrix(LatDataAlloc flag, int n_vec = 1 , int vol = 0 ){
//      printf("LatMatrix::LatMatrix(f,i,i)\n");
      Init(flag,n_vec,vol);
    }
    LatMatrix(int n_vec, int vol = 0){
//      printf("LatMatrix::LatMatrix(i,i)\n");
      Init(DEFAULT,n_vec,vol);
    }
    ~LatMatrix();
    void Init(LatDataAlloc flag, int n_vec , int vol );
    Matrix *Mat(int pos=0, int mat_row=0); 
    Float *Field(int pos=0, int mat_row=0, int n=0) ; 
    LatMatrix &operator=(const LatMatrix &lat);
    void operator= (IFloat c); 
    void operator+= (LatMatrix & lat_mat); 
    void MulDag (LatMatrix &c);
    void TrLessAntiHermMatrix();
    Float norm();
};
CPS_END_NAMESPACE
#endif
