#include <config.h>
#include <util/qcdio.h>
#include <comms/sysfunc.h>
#include <util/sproj_tr.h>
#include <util/vector.h>
#include <util/time.h>

USING_NAMESPACE_CPS

int main(int argc, char *argv[])
{

  int num,stride,stride2;
  if(argc<4){
    printf("Usage: %s num stride stride\n",argv[0]);
    exit(-1);
  }
  sscanf(argv[1],"%d",&num);
  printf("num=%d\n",num);
  sscanf(argv[2],"%d",&stride);
  printf("stride=%d\n",stride);
  sscanf(argv[3],"%d",&stride2);
  printf("stride2=%d\n",stride2);
  IFloat *final = (IFloat *)qalloc(QFAST,sizeof(Matrix));
//  IFloat final[18];
  IFloat *v =(IFloat *)qalloc(QFAST,sizeof(Vector)*4*num*(stride+1));
  IFloat *w =(IFloat *)qalloc(QFAST,sizeof(Vector)*4*num*(stride2+1));
  for(int i = 0; i<(6*4*num*(stride+1)) ;i++) v[i] = (IFloat)i;
  for(int i = 0; i<(6*4*num*(stride2+1)) ;i++) w[i] = (IFloat)i;
  for(int i = 0; i<18 ;i++) final[i] = (IFloat)i;
//  sprojTrXp(final,v,w,num,stride,stride2);
  Float time = -dclock();
  sprojTrXp(final,v,w,num,stride,stride2);
  time +=dclock();
  print_flops(num*168,time);
  time = -dclock();
  sprojTrXm(final,v,w,num,stride,stride2);
  time +=dclock();
  print_flops(num*168,time);
  qfree(final);
  qfree(v);
  qfree(w);
}
