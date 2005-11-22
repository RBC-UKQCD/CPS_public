#ifndef DATA_SHIFT_H
#define DATA_SHIFT_H
#include <config.h>
CPS_START_NAMESPACE

class GlobalDataShift{
  private:
    static char *cname;
    int shifts[5];
    int origin[5];
    void *addr;
    long long data_len;
    void *temp_buf;
    void Shift(int i, int n_disp);
  public:
  GlobalDataShift();
  ~GlobalDataShift(){};
  void Set(int x, int y, int z, int t, int s=0);
  void SetOrigin(int x, int y, int z, int t, int s=0);
  int Origin(int dir)
    {return origin[dir];}
  void Shift(void *addr, long long len);
};

extern GlobalDataShift GDS;
CPS_END_NAMESPACE
#endif
