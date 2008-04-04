#include <config.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/data_shift.h>
CPS_START_NAMESPACE

GlobalDataShift GDS;

char *GlobalDataShift::cname = "GlobalDataShift";

GlobalDataShift::GlobalDataShift(){
  for(int i =0;i<5;i++) shifts[i]=0;
  for(int i =0;i<5;i++) origin[i]=0;
}

void GlobalDataShift::Set(int x, int y, int z, int t, int s)
{
  char *fname = "Set(x,y,z,t,s)";
  shifts[0]=x;
  shifts[1]=y;
  shifts[2]=z;
  shifts[3]=t;
  shifts[4]=s;
  for(int i =0;i<5;i++){
//    if( (GJP.Bc(i)!=BND_CND_PRD) && (shifts[i]!=0) )
//      ERR.General(cname,fname,
//        "Bc(%d) cannot be shifted due to the boundary condition\n",i);
    shifts[i] = shifts[i]%GJP.Nodes(i);
//    printf("shifts[%d]=%d\n",i,shifts[i]);
    VRB.Flow(cname,fname,"shifts[%d]=%d\n",i,shifts[i]);
  }
}

void GlobalDataShift::SetOrigin(int x, int y, int z, int t, int s)
{
  char *fname = "SetOrigin(x,y,z,t,s)";
  origin[0]=x;
  origin[1]=y;
  origin[2]=z;
  origin[3]=t;
  origin[4]=s;
  for(int i =0;i<5;i++){
    if( (GJP.Bc(i)!=BND_CND_PRD) && (origin[i]!=0) )
      ERR.General(cname,fname,
        "Bc(%d) cannot be shifted due to the boundary condition\n",i);
    origin[i] = origin[i]%GJP.Nodes(i);
//    printf("origin[%d]=%d\n",i,origin[i]);
    VRB.Flow(cname,fname,"origin[%d]=%d\n",i,origin[i]);
  }
}

void GlobalDataShift::Shift(void *a, long long len){
  char *fname = "Shift(a,i)";
  VRB.Func(cname,fname);
  if (len%8)
    ERR.General(cname,fname,"buffer length should be a multiple of 8\n");
  addr = a; data_len = len;
  temp_buf = smalloc(cname,fname,"temp_buf",len);

  for(int i =0;i<5;i++){
    Shift(i, shifts[i]);
  }
  sfree(cname,fname,"temp_buf",temp_buf);
  temp_buf=NULL;
}
CPS_END_NAMESPACE
