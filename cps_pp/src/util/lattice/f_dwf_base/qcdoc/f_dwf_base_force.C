#include<config.h>
#include<qalloc.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of FdwfBase class.

  $Id: f_dwf_base_force.C,v 1.12 2008-02-08 18:35:08 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_dwf_base/qcdoc/f_dwf_base_force.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_dwf_base.C
//
// FdwfBase is derived from FwilsonTypes and is relevant to
// domain wall fermions, moved from old Fdwf class
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <math.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/dwf.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/random.h>
#include <util/error.h>
#include <util/time_cps.h>
#include <comms/scu.h> // GRF
#include <comms/glb.h>
CPS_START_NAMESPACE


static int offset (int *size, int *pos, int mu = -1){
  if (mu>3) printf("FdwfBase::offset: Error!\n");
  if (mu>-1) pos[mu] = (pos[mu]+1)%size[mu];
  int result = 
  pos[0] + size[0] *( pos[1] + size[1] *( pos[2] + size[2] *( pos[3])));
  if (mu>-1) pos[mu] = (pos[mu]-1+size[mu])%size[mu];
  return result;
}

#undef PROFILE
// CJ: change start
//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *chi, Float mass, 
//                 Float dt):
// It evolves the canonical momentum mom by dt
// using the fermion force.
//------------------------------------------------------------------
ForceArg FdwfBase::EvolveMomFforce(Matrix *mom, Vector *chi, 
			   Float mass, Float dt){
  char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
  VRB.Func(cname,fname);
  Matrix *gauge = GaugeField() ;

  if (Colors() != 3)
    ERR.General(cname,fname,"Wrong nbr of colors.") ;
 
  if (SpinComponents() != 4)
    ERR.General(cname,fname,"Wrong nbr of spin comp.") ;
 
  if (mom == 0)
    ERR.Pointer(cname,fname,"mom") ;
 
  if (chi == 0)
    ERR.Pointer(cname,fname,"chi") ;
 
  //----------------------------------------------------------------
  // allocate space for two CANONICAL fermion fields
  //----------------------------------------------------------------

  int f_size = FsiteSize() * GJP.VolNodeSites() ;
  int f_site_size_4d = 2 * Colors() * SpinComponents();
  int f_size_4d = f_site_size_4d * GJP.VolNodeSites() ;
 
  char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)fmalloc(cname,fname,str_v1,f_size*sizeof(Float));
//  if (v1 == 0) ERR.Pointer(cname, fname, str_v1) ;
//  VRB.Smalloc(cname, fname, str_v1, v1, f_size*sizeof(Float)) ;

  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)fmalloc(cname,fname,str_v2,f_size*sizeof(Float)) ;
//  if (v2 == 0) ERR.Pointer(cname, fname, str_v2) ;
//  VRB.Smalloc(cname, fname, str_v2, v2, f_size*sizeof(Float)) ;

//  LatMatrix MomDiff(QFAST,4);
//  Matrix *mom_diff = MomDiff.Mat();

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;

  //----------------------------------------------------------------
  // Calculate v1, v2. Both v1, v2 must be in CANONICAL order after
  // the calculation.
  //----------------------------------------------------------------  

  VRB.Clock(cname, fname, "Before calc force vecs.\n") ;

#ifdef PROFILE
  Float time = -dclock();
  ForceFlops=0;
#endif
  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;

    DiracOpDwf dwf(*this, v1, v2, &cg_arg, CNV_FRM_NO) ;
    dwf.CalcHmdForceVecs(chi) ;
    Fconvert(v1,CANONICAL,WILSON);
    Fconvert(v2,CANONICAL,WILSON);
  }
#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif

  int mu, x, y, z, t, s, lx, ly, lz, lt, ls ;
  int size[5],surf[4],blklen[4],stride[4],numblk[4];
 
  size[0] = lx = GJP.XnodeSites() ;
  size[1] = ly = GJP.YnodeSites() ;
  size[2] = lz = GJP.ZnodeSites() ;
  size[3] = lt = GJP.TnodeSites() ;
  size[4] = ls = GJP.SnodeSites() ;

  blklen[0] = sizeof(Float)*FsiteSize()/size[4];
  numblk[0] = GJP.VolNodeSites()/size[0]*size[4];
  stride[0] = blklen[0] * (size[0]-1);
  for (int i =1;i<4;i++){
    blklen[i] = blklen[i-1] * size[i-1];
	numblk[i] = numblk[i-1] / size[i];
    stride[i] = blklen[i] * (size[i]-1);
  }
  for (int i =0;i<4;i++){
    surf[i] = GJP.VolNodeSites()/size[i];
  }

  //----------------------------------------------------------------
  // allocate buffer space for two fermion fields that are assoc
  // with only one 4-D site.
  //----------------------------------------------------------------

  unsigned long flag = QFAST;
  if (ls<4) flag|= QNONCACHE;
  char *str_site_v1 = "site_v1" ;
  char *str_site_v2 = "site_v2" ;

  Float *v1_buf[4],*v2_buf[4];
  int pos[4];
  Float *v1_buf_pos[4];
  Float *v2_buf_pos[4];
  for(int i =0;i<4;i++) {
    v1_buf[i]=(Float *)fmalloc(cname,fname,"v1_buf",surf[i]*FsiteSize()*sizeof(Float)) ;
    v2_buf[i]=(Float *)fmalloc(cname,fname,"v2_buf",surf[i]*FsiteSize()*sizeof(Float)) ;
  }

//  Matrix tmp_mat1, tmp_mat2 ;
  Matrix *tmp_mat1,*tmp_mat2;
  tmp_mat1 = (Matrix *)fmalloc(cname,fname,"tmp_mat1",sizeof(Matrix)*2);
  tmp_mat2 = tmp_mat1+1;

  SCUDirArgIR Send[4];
  SCUDirArgIR Recv[4];
  SCUDMAInst *dma[2];
  for(int i = 0;i<2;i++) dma[i] = new SCUDMAInst;
  void *addr[2];
  int f_bytes = sizeof(Float)*f_site_size_4d;
  int st_bytes = sizeof(Float)*f_size_4d - f_bytes;
  for (mu=0; mu<4; mu++){
    if ( GJP.Nodes(mu) >1){
      dma[0]->Init(v1_buf[mu],surf[mu]*ls*f_bytes);
      dma[1]->Init(v2_buf[mu],surf[mu]*ls*f_bytes);
      Recv[mu].Init(gjp_scu_dir[2*mu],SCU_REC,dma,2);
      dma[0]->Init(v1,blklen[mu],numblk[mu],stride[mu]);
      dma[1]->Init(v2,blklen[mu],numblk[mu],stride[mu]);
      Send[mu].Init(gjp_scu_dir[2*mu+1],SCU_SEND,dma,2);
    }
  }
  for(int i = 0;i<2;i++) delete dma[i];

  VRB.Clock(cname, fname, "Before loop over links.\n") ;

#ifdef PROFILE
  time = -dclock();
  ForceFlops=0;
#endif

  sys_cacheflush(0);
  for (mu=0; mu<4; mu++){
    if ( GJP.Nodes(mu) >1){
      Recv[mu].StartTrans();
      Send[mu].StartTrans();
    }
  }
  for (mu=0; mu<4; mu++){
    for (pos[3]=0; pos[3]<size[3]; pos[3]++){
    for (pos[2]=0; pos[2]<size[2]; pos[2]++){
    for (pos[1]=0; pos[1]<size[1]; pos[1]++){
    for (pos[0]=0; pos[0]<size[0]; pos[0]++){
      int gauge_offset = offset(size,pos);
      int vec_offset = f_site_size_4d*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;

//      printf("%d %d %d %d %d\n",mu,pos[0],pos[1],pos[2],pos[3]);
      Float *v1_plus_mu ;
      Float *v2_plus_mu ;
      int vec_plus_mu_stride ;
      int vec_plus_mu_offset = f_site_size_4d ;

      Float coeff = -2.0 * dt ;

          vec_plus_mu_offset *= offset(size,pos,mu);
      if ((GJP.Nodes(mu)>1)&&((pos[mu]+1) == size[mu]) ) {
      } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
            vec_plus_mu_stride = f_size_4d - f_site_size_4d ;
            if ((pos[mu]+1) == size[mu])
            if (GJP.NodeBc(mu)==BND_CND_APRD) coeff = -coeff ;

//     Float time2 = -dclock();
      sproj_tr[mu]( (IFloat *)tmp_mat1,
                    (IFloat *)v1_plus_mu,
                    (IFloat *)v2+vec_offset,
                    ls, vec_plus_mu_stride, f_size_4d-f_site_size_4d) ;

      sproj_tr[mu+4]( (IFloat *)tmp_mat2,
                      (IFloat *)v2_plus_mu,
                      (IFloat *)v1+vec_offset,
                      ls, vec_plus_mu_stride, f_size_4d-f_site_size_4d) ;
 //    time2 += dclock();
 //    print_flops("","sproj",2*9*16*ls,time2);

      *tmp_mat1 += *tmp_mat2 ;

      // If GJP.Snodes > 1 sum up contributions from all s nodes
      if(GJP.Snodes() > 1) {
	  glb_sum_multi_dir((Float *)tmp_mat1, 4, sizeof(Matrix)/sizeof(IFloat) ) ;
      }

      tmp_mat2->DotMEqual(*(gauge+gauge_offset), *tmp_mat1) ;

      tmp_mat1->Dagger(*tmp_mat2) ;

      tmp_mat2->TrLessAntiHermMatrix(*tmp_mat1) ;

      *tmp_mat2 *= coeff ;

      *(mom+gauge_offset) += *tmp_mat2 ;

      Float norm = tmp_mat2->norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }

    } } } } // end for x,y,z,t
  } // end for mu
  for (mu=0; mu<4; mu++){
    if ( GJP.Nodes(mu) >1){
      Recv[mu].TransComplete();
      Send[mu].TransComplete();
    }
  }
 
//------------------------------------------------------------------
// start by summing first over direction (mu) and then over site
// to allow SCU transfers to happen face-by-face in the outermost
// loop.
//------------------------------------------------------------------

  for(int i = 0;i<4;i++){
    v1_buf_pos[i] = v1_buf[i];
    v2_buf_pos[i] = v2_buf[i];
  }

  for (mu=0; mu<4; mu++){
    for (pos[3]=0; pos[3]<size[3]; pos[3]++){
    for (pos[2]=0; pos[2]<size[2]; pos[2]++){
    for (pos[1]=0; pos[1]<size[1]; pos[1]++){
    for (pos[0]=0; pos[0]<size[0]; pos[0]++){
      int gauge_offset = offset(size,pos);
      int vec_offset = f_site_size_4d*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;

      Float *v1_plus_mu ;
      Float *v2_plus_mu ;
      int vec_plus_mu_stride ;
      int vec_plus_mu_offset = f_site_size_4d ;

      Float coeff = -2.0 * dt ;

//      printf("%d %d %d %d %d\n",mu,pos[0],pos[1],pos[2],pos[3]);
          vec_plus_mu_offset *= offset(size,pos,mu);
      if ((GJP.Nodes(mu)>1)&&((pos[mu]+1) == size[mu]) ) {
            v1_plus_mu = v1_buf_pos[mu] ;
            v2_plus_mu = v2_buf_pos[mu] ;
 			v1_buf_pos[mu] += f_site_size_4d;
 			v2_buf_pos[mu] += f_site_size_4d;
            vec_plus_mu_stride = (surf[mu] -1)*f_site_size_4d ;
            if (GJP.NodeBc(mu)==BND_CND_APRD) coeff = -coeff ;

//     Float time2 = -dclock();
      sproj_tr[mu]( (IFloat *)tmp_mat1,
                    (IFloat *)v1_plus_mu,
                    (IFloat *)v2+vec_offset,
                    ls, vec_plus_mu_stride, f_size_4d-f_site_size_4d) ;

      sproj_tr[mu+4]( (IFloat *)tmp_mat2,
                      (IFloat *)v2_plus_mu,
                      (IFloat *)v1+vec_offset,
                      ls, vec_plus_mu_stride, f_size_4d-f_site_size_4d) ;
 //    time2 += dclock();
 //    print_flops("","sproj",2*9*16*ls,time2);

      *tmp_mat1 += *tmp_mat2 ;

      // If GJP.Snodes > 1 sum up contributions from all s nodes
      if(GJP.Snodes() >1 ) {
	  glb_sum_multi_dir((Float *)tmp_mat1, 4, sizeof(Matrix)/sizeof(IFloat) ) ;
      }

      tmp_mat2->DotMEqual(*(gauge+gauge_offset), *tmp_mat1) ;

      tmp_mat1->Dagger(*tmp_mat2) ;

      tmp_mat2->TrLessAntiHermMatrix(*tmp_mat1) ;

      *tmp_mat2 *= coeff ;

      *(mom+gauge_offset) += *tmp_mat2 ;
      Float norm = tmp_mat2->norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }

    } } } } // end for x,y,z,t
  } // end for mu

  ForceFlops += (2*9*16*ls + 18+ 198+36+24)*lx*ly*lz*lt*4;
#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif
 
//------------------------------------------------------------------
// deallocate smalloc'd space
//------------------------------------------------------------------

  for(int i =0;i<4;i++) {
    ffree(v1_buf[i],cname,fname,"v1_buf");
    ffree(v2_buf[i],cname,fname,"v2_buf");
  }
 
  VRB.Sfree(cname, fname, str_v2, v2) ;
  ffree(v2) ;
 
  VRB.Sfree(cname, fname, str_v1, v1) ;
  ffree(v1) ;

  ffree(tmp_mat1);
 
  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  return ForceArg(L1, sqrt(L2), Linf);
}
// CJ: change end

CPS_END_NAMESPACE
