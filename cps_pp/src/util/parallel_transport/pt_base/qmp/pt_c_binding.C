#include <config.h>
#ifdef USE_QMP
#include <util/pt.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/error.h>
//#include "pt_int.h"

CPS_START_NAMESPACE
static PT StaticPT;

void pt_init(Lattice &Lat){
  PTArg pt_arg;
//Size of local volume in all four directions
  pt_arg.size[0] = GJP.XnodeSites();
  pt_arg.size[1] = GJP.YnodeSites();
  pt_arg.size[2] = GJP.ZnodeSites();
  pt_arg.size[3] = GJP.TnodeSites();
  //-------------------------------------------------
  //Added by Michael C.
  for(int i = 0; i < 4;i++)
    if(GJP.Nodes(i) == 1)
      pt_arg.local[i] = 1;
    else
      pt_arg.local[i] = 0;
  //-------------------------------------------------
  int temp = 0;
  for(int i = 0;i<4;i++){
    temp += GJP.NodeSites(i)*GJP.NodeCoor(i);
  }
  if (temp%2==0)
    pt_arg.evenodd = PT_EVEN;
  else
    pt_arg.evenodd = PT_ODD;
  pt_arg.gauge_field_addr = (IFloat *)Lat.GaugeField();

  StrOrdType str_ord = Lat.StrOrd();
  printf("str_ord=%d\n",str_ord);

  switch(str_ord){
    case CANONICAL:
      pt_arg.g_str_ord = PT_XYZT;
      pt_arg.v_str_ord = PT_XYZT;
      pt_arg.v_str_ord_cb = PT_TXYZ;
      pt_arg.g_conj = 0;
      break;
    case WILSON:
      pt_arg.g_str_ord = PT_XYZT_CB_O;
      pt_arg.v_str_ord = PT_XYZT_CB_O;
      pt_arg.v_str_ord_cb = PT_TXYZ;
      pt_arg.g_conj = 0;
      break;
    case STAG:
      pt_arg.g_str_ord = PT_XYZT;
      pt_arg.v_str_ord = PT_XYZT;
      pt_arg.v_str_ord_cb = PT_TXYZ;
      pt_arg.g_conj = 1;
      break;
    default:
      break;
  }
  pt_arg.prec = sizeof(IFloat);
  StaticPT.init(&pt_arg);
}

void pt_init_g(){
  StaticPT.init_g();
}

void pt_delete(){
  StaticPT.delete_buf();
}

void pt_delete_g(){
  StaticPT.delete_g_buf();
}

void pt_1vec(int n, IFloat **mout, IFloat **min, int const *dir){
  StaticPT.vec(n,mout,min,dir);
  ParTrans::PTflops +=33*n*PT::vol;
}

void pt_mat(int n, IFloat **mout, IFloat **min, int const *dir){
  StaticPT.mat(n,(matrix**)mout,(matrix**)min,dir);
  ParTrans::PTflops +=198*n*PT::vol;
}

#if 1
void pt_vvpd(IFloat **vect, int n_vect, const int *dir,
             int n_dir, int hop, IFloat **sum){
  StaticPT.vvpd(vect, n_vect, dir, n_dir, hop, sum);
  ParTrans::PTflops +=90*n_vect*n_dir*PT::vol;
}

void pt_vvpd(IFloat **vect2, IFloat ***vect, int n_vect, const int *dir,
             int n_dir, int hop, IFloat **sum, int overwrite){
  StaticPT.vvpd(vect2, vect, n_vect, dir, n_dir, hop, sum, overwrite);
  ParTrans::PTflops +=90*n_vect*n_dir*PT::vol;
}
#endif

void pt_shift_link(IFloat **u, const int *dir, int n_dir){
  StaticPT.shift_link(u,dir,n_dir);
}

void pt_shift_field(IFloat **v, const int *dir, int n_dir,
                    int hop, IFloat **u){
  StaticPT.shift_field(v,dir,n_dir,hop,u);
}

void pt_shift_field_vec(IFloat **v, const int *dir, int n_dir,
                    int hop, IFloat **u){
  StaticPT.shift_field_vec(v,dir,n_dir,hop,u);
}

void pt_mat_cb(int n, Float **mout, Float **min, const int *dir, ChkbType cb)
{
  StaticPT.mat_cb(n,mout,min,dir,cb);
}

void pt_mat_cb(int n, Float **mout, Float **min, const int *dir, ChkbType cb, Float * new_gauge_field)
{
  StaticPT.mat_cb(n,mout,min,dir,cb,new_gauge_field);
}

void pt_1vec_cb(int n, Float **vout, Float **vin, const int *dir, ChkbType cb)
{
  StaticPT.vec_cb(n,vout,vin,dir,cb);
  ParTrans::PTflops +=33*n*PT::vol;
}

void pt_1vec_cb(int n, Float **vout, Float **vin, const int *dir, ChkbType cb, Float * new_gauge_field)
{
  StaticPT.vec_cb(n,vout,vin,dir,cb,new_gauge_field);
  ParTrans::PTflops +=33*n*PT::vol;
}

void pt_1vec_cb(int n, Float *vout, Float **vin, const int *dir, ChkbType cb, int pad)
{
  StaticPT.vec_cb(n,vout,vin,dir,cb,pad);
  ParTrans::PTflops +=33*n*PT::vol;
}

void pt_1vec_cb(int n, Float *vout, Float **vin, const int *dir, ChkbType cb, int pad, Float * new_gauge_field)
{
  StaticPT.vec_cb(n,vout,vin,dir,cb,pad,new_gauge_field);
  ParTrans::PTflops +=33*n*PT::vol;
}
CPS_END_NAMESPACE
#endif
