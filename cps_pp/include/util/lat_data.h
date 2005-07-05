#include <config.h>
#include <stdio.h>

#ifndef LAT_DATA_H
#define LAT_DATA_H
#include <util/vector.h>
#include <util/data_types.h>

CPS_START_NAMESPACE
enum { NEW = 0, INITTED};
class LatData{
  private:
    char *cname;
    int status;
  protected:
    static int DEFAULT_FLAG;
    int size; //number of IFlots per site
    int vol;  //number of sites
  public:
    IFloat *data; 
    LatData(const LatData &lat);
    LatData(){status = NEW;};
    LatData(int flags, int size, int vol){
      Init(flags,size,vol);
    }
    LatData(int size, int vol){
      Init(DEFAULT_FLAG,size,vol);
    }
    void  Init(int flags, int size, int vol);
    void  Init(int size, int vol){
      Init(DEFAULT_FLAG,size,vol);
    }
    LatData &operator=(const LatData &lat);
    ~LatData();
    const IFloat *Field(int pos=0, int n=0)
      { return (data+pos*size+n);}
    const int Size(){return size*vol;}

    friend void glb_sum_multi_dir(LatData &dat, int dir);

};

class LatVector: public LatData{
  private:
    int vec_size;
  public: 
    LatVector(int flag, int n_vec , int vol )
      { Init(flag,n_vec,vol); }
    LatVector(int n_vec , int vol =0)
      { Init(DEFAULT_FLAG,n_vec,vol); }
    void Init(int flag, int n_vec , int vol );
    ~LatVector();
    Vector *Vec(int pos=0, int vec_row=0);
    Float *Field(int pos=0, int vec_row=0,int n=0); 

    void FTimesV1PlusV2(const Float &fb, LatVector *c,
                        LatVector *d){
        fTimesV1PlusV2(data,Float(fb), c->Field(),d->Field(),size*vol);
    }
};

class LatMatrix: public LatData{
  private:
    int mat_size;
  public: 
    LatMatrix(int flag, int n_vec , int vol );
    LatMatrix(int n_vec = 1, int vol = 0){
      LatMatrix(DEFAULT_FLAG,n_vec,vol);
    }
    ~LatMatrix();
    Matrix *Mat(int pos=0, int mat_row=0); 
    Float *Field(int pos=0, int mat_row=0, int n=0) ; 
};
CPS_END_NAMESPACE
#endif
