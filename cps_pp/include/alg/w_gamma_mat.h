#include<config.h>
CPS_START_NAMESPACE
/*w_gamma_mat.h
 *define 4-D dirac matrices
 *to reduce memory usage
 *Total memory used is about 4x4x2x10=0x140
 */

#ifndef W_GAMMA_MAT_H
#define W_GAMMA_MAT_H

CPS_END_NAMESPACE
#include <util/data_types.h>
CPS_START_NAMESPACE

//1->x, 2->y, 3->z, 4->t
//gamma_i[D1][D2][COMPLEXES]

enum WGammaMatrix{
  WUNIT=0,
  WGAM_1=1,
  WGAM_2=2,
  WGAM_3=4,
  WGAM_4=8,
  WGAM_5=15,
  WGAM_1_2=3, //gam_x*gam_y
  WGAM_1_3=5,
  WGAM_1_4=9,
  WGAM_1_5=14,
  WGAM_2_3=6,
  WGAM_2_4=10,
  WGAM_5_2=13,
  WGAM_3_4=12,
  WGAM_3_5=11,
  WGAM_5_4=7
  //...
};


//definition of 16 Gamma Matrices in w_gamma_mat.C

extern IFloat WGamma[16][4][4][2];


#endif

/*
 *do sign*Trace(gam1*G1*gam2*transpose(G2))
 *sign=1 or -1
 *using data in d_zero_mom_proj[D1x][D1y][D2x][D2y]
 *(momproj in short)
 *Explicitly: 
 *result=Sum_d1d2d3d4(gamm1[d1][d2]mompropj[d2][d3][d1][d4]gam2[d3][d4])
 */

/*?????if prop_dir!=3, DEV_I's actual polirization changes
 *when preparing gamm matrix,must take care of this
 *???New function IFloat *getGamMat(int pol_dir1, int pol_dir2=-1)
 *                      which get Gamm(pol_dir1) if pol_dir=-1
 *                                Gamm(pol_dir1)*Gam(pol_dir2) if pol_dir=0(X),1(Y),2(Z),3(T)
 *use pol=WspectQuarkdev::devDirToPolDir(i,d_prop_dir) to translate first
 */

//WspectExtendedMesons::traceDirac(IFloat* gam1,IFloat* gam2 ,int mesonId,int lclw){
//

//?? write the same one for WspectMesons too???
//result_p point to the correlator data entry(Complex number)
//sign is the sign of contribution to correlator of this trace term




CPS_END_NAMESPACE
