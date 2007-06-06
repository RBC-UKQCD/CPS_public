#include <config.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/asqtad_int.h>
#include <util/asqtad.h>

//AsqD asqd;
AsqD *asqd_p;
CPS_START_NAMESPACE
static Fasqtad *lat_pt;

extern "C" void asqtad_dirac_init(Fasqtad *lat){
  AsqDArg arg;
  lat_pt = lat;
  arg.size[0] = GJP.TnodeSites();
  arg.size[1] = GJP.XnodeSites();
  arg.size[2] = GJP.YnodeSites();
  arg.size[3] = GJP.ZnodeSites();
  int vol = 1;
  for(int i =0;i<4;i++) vol *= arg.size[i];
  arg.NP[0] = GJP.Tnodes();
  arg.NP[1] = GJP.Xnodes();
  arg.NP[2] = GJP.Ynodes();
  arg.NP[3] = GJP.Znodes();
  arg.coor[0] = GJP.TnodeCoor();
  arg.coor[1] = GJP.XnodeCoor();
  arg.coor[2] = GJP.YnodeCoor();
  arg.coor[3] = GJP.ZnodeCoor();
  arg.c1 = GJP.KS_coeff();
  arg.c2 = GJP.Naik_coeff();
  arg.c3 = GJP.staple3_coeff();
  arg.c5 = GJP.staple5_coeff();
  arg.c7 = GJP.staple7_coeff();
  arg.c6 = GJP.Lepage_coeff();
  Float *tmp_fat = (Float *)lat_pt->Fields(0);
  Float *tmp_naik = (Float *)lat_pt->Fields(1);
  Float *tmp_naik_m = (Float *)lat_pt->Fields(2);
  for(int i =0;i<4;i++){
    arg.Fat[i] = tmp_fat+vol*18*i;
    arg.Naik[i] = tmp_naik+vol*18*i;
    arg.NaikM[i] = tmp_naik_m+vol*18*i;
  }
//  arg.Fat = (IFloat *)lat_pt->Fields(0);
//  arg.Naik = (IFloat *)lat_pt->Fields(1);
//  arg.NaikM = (IFloat *)lat_pt->Fields(2);
//  arg.NaikM = NULL;
  asqd_p = new AsqD;
  asqd_p->init(&arg);
}

void asqtad_dirac_init_g(IFloat *frm_tmp){
  lat_pt->Smear();
  asqd_p->init_g(frm_tmp);
}

void asqtad_destroy_dirac_buf(){
  asqd_p->destroy_buf();
  delete asqd_p;
}

void asqtad_destroy_dirac_buf_g(){
  asqd_p->destroy_buf_g();
}
CPS_END_NAMESPACE
