#include <config.h>

#ifndef LAT_DATA_H
#define LAT_DATA_H
#include <util/vector.h>

CPS_START_NAMESPACE
enum { NEW = 0, INITTED};
class LatData{
  private:
    char *cname;
    int status;
    static int DEFAULT_FLAG;
  protected:
    int size; //number of IFlots per site
    int vol;  //number of sites
    IFloat *data; 
  public:
    LatData(){cname = "LatData"; status = NEW; data=NULL;};
    inline LatData(int flags, int size, int vol){
      Init(flags,size,vol);
    }
    inline LatData(int size, int vol){
      LatData(DEFAULT_FLAG,size,vol);
    }
    int Init(int flags, int size, int vol);
    inline int Init(int size, int vol){
      Init(DEFAULT_FLAG,size,vol);
    }
    ~LatData();
    IFloat *Field(int pos=0, int n=0);
    int Size(){return size*vol;}
};

class LatVector: public LatData{
  private:
    int vec_size;
  public: 
    LatVector(int n_vec = 1, int vol = 0);
    ~LatVector(){};
    Vector *Field(int pos, int vec_row); 
};

class LatMatrix: public LatData{
  private:
    int mat_size;
  public: 
    LatMatrix(int n_vec = 1, int vol = 0);
    ~LatMatrix();
    Matrix *Field(int pos, int mat_row); 
};
CPS_END_NAMESPACE
#endif
