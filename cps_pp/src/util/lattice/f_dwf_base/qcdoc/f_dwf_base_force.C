#include<config.h>
#include<qalloc.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of FdwfBase class.

  $Id: f_dwf_base_force.C,v 1.2 2005-02-12 00:54:40 chulwoo Exp $
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
#include <util/time.h>
#include <comms/scu.h> // GRF
#include <comms/glb.h>
CPS_START_NAMESPACE

#define PROFILE

static int offset (int *size, int *pos, int mu = -1){
  if (mu>3) printf("FdwfBase::offset: Error!\n");
  if (mu>-1) pos[mu] = (pos[mu]+1)%size[mu];
  int result = 
  pos[0] + size[0] *( pos[1] + size[1] *( pos[2] + size[2] *( pos[3])));
  if (mu>-1) pos[mu] = (pos[mu]-1+size[mu])%size[mu];
  return result;
}

// CJ: change start
//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *chi, Float mass, 
//                 Float step_size):
// It evolves the canonical momentum mom by step_size
// using the fermion force.
//------------------------------------------------------------------
void FdwfBase::EvolveMomFforce(Matrix *mom, Vector *chi, 
			   Float mass, Float step_size){
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
  Vector *v1 = (Vector *)fmalloc(f_size*sizeof(Float)) ;
  if (v1 == 0) ERR.Pointer(cname, fname, str_v1) ;
  VRB.Smalloc(cname, fname, str_v1, v1, f_size*sizeof(Float)) ;

  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)fmalloc(f_size*sizeof(Float)) ;
  if (v2 == 0) ERR.Pointer(cname, fname, str_v2) ;
  VRB.Smalloc(cname, fname, str_v2, v2, f_size*sizeof(Float)) ;

  //----------------------------------------------------------------
  // Calculate v1, v2. Both v1, v2 must be in CANONICAL order after
  // the calculation.
  //----------------------------------------------------------------  

  VRB.Clock(cname, fname, "Before calc force vecs.\n") ;

  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;

    DiracOpDwf dwf(*this, v1, v2, &cg_arg, CNV_FRM_YES) ;
    dwf.CalcHmdForceVecs(chi) ;
  }

  int mu, x, y, z, t, s, lx, ly, lz, lt, ls ;
  int size[5],surf[4];
 
  size[0] = lx = GJP.XnodeSites() ;
  size[1] = ly = GJP.YnodeSites() ;
  size[2] = lz = GJP.ZnodeSites() ;
  size[3] = lt = GJP.TnodeSites() ;
  size[4] = ls = GJP.SnodeSites() ;

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
  for(int i =0;i<4;i++) {
    v1_buf[i]=(Float *)smalloc(cname,fname,"v1_buf",surf[i]*FsiteSize()*sizeof(Float)) ;
    v2_buf[i]=(Float *)smalloc(cname,fname,"v2_buf",surf[i]*FsiteSize()*sizeof(Float)) ;
  }

    Float *site_v1=(Float *)qalloc(QFAST|QNONCACHE,FsiteSize()*sizeof(Float)) ;
    if (site_v1 == 0) ERR.Pointer(cname, fname, str_site_v1) ;
    VRB.Smalloc(cname, fname, str_site_v1, site_v1, FsiteSize()*sizeof(Float)) ;
    Float *site_v2=(Float *)qalloc(QFAST|QNONCACHE,FsiteSize()*sizeof(Float)) ;
    if (site_v2 == 0) ERR.Pointer(cname, fname, str_site_v2) ;
    VRB.Smalloc(cname, fname, str_site_v2, site_v2, FsiteSize()*sizeof(Float)) ;


  Matrix tmp_mat1, tmp_mat2 ;

  SCUDirArgIR Send[4];
  SCUDirArgIR Recv[4];
  SCUDMAInst *dma[2];
  for(int i = 0;i<2;i++) dma[i] = new SCUDMAInst;
  void *addr[2];
  int f_bytes = sizeof(Float)*f_site_size_4d;
  int st_bytes = sizeof(Float)*f_size_4d - f_bytes;
  for (mu=0; mu<4; mu++){
    if ( GJP.Nodes(mu) >1){
      dma[0]->Init(site_v1,ls*f_bytes);
      dma[1]->Init(site_v2,ls*f_bytes);
      Recv[mu].Init(gjp_scu_dir[2*mu],SCU_REC,dma,2);
      dma[0]->Init(v1,f_bytes,ls,st_bytes);
      dma[1]->Init(v2,f_bytes,ls,st_bytes);
      Send[mu].Init(gjp_scu_dir[2*mu+1],SCU_SEND,dma,2);
    }
  }
  for(int i = 0;i<2;i++) delete dma[i];

  VRB.Clock(cname, fname, "Before loop over links.\n") ;

#ifdef PROFILE
  Float time = -dclock();
  ForceFlops=0;
#endif
 
//------------------------------------------------------------------
// start by summing first over direction (mu) and then over site
// to allow SCU transfers to happen face-by-face in the outermost
// loop.
//------------------------------------------------------------------
  sys_cacheflush(0);

  int pos[4];
  for (mu=0; mu<4; mu++){
    for (pos[3]=0; pos[3]<size[3]; pos[3]++){
    for (pos[2]=0; pos[2]<size[2]; pos[2]++){
    for (pos[1]=0; pos[1]<size[1]; pos[1]++){
    for (pos[0]=0; pos[0]<size[0]; pos[0]++){
//      int gauge_offset = x+lx*(y+ly*(z+lz*t)) ;
      int gauge_offset = offset(size,pos);
      int vec_offset = f_site_size_4d*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;

      Float *v1_plus_mu ;
      Float *v2_plus_mu ;
      int vec_plus_mu_stride ;
      int vec_plus_mu_offset = f_site_size_4d ;

      Float coeff = -2.0 * step_size ;

          vec_plus_mu_offset *= offset(size,pos,mu);
          if ((GJP.Nodes(mu)>1)&&((pos[mu]+1) == size[mu]) ) {
#if 1
          addr[0] = (IFloat *)v1+vec_plus_mu_offset;
          addr[1] = (IFloat *)v2+vec_plus_mu_offset;
          Send[mu].Addr(addr,2);
	  Recv[mu].StartTrans();
	  Send[mu].StartTrans();
	  Recv[mu].TransComplete();
	  Send[mu].TransComplete();
#else
            for (s=0; s<ls; s++) {
              getPlusData( (IFloat *)site_v1+s*f_site_size_4d,
                (IFloat *)v1+vec_plus_mu_offset+s*f_size_4d,
                f_site_size_4d, mu) ;
              getPlusData( (IFloat *)site_v2+s*f_site_size_4d,
                (IFloat *)v2+vec_plus_mu_offset+s*f_size_4d,
                f_site_size_4d, mu) ;
            } // end for s
#endif
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            vec_plus_mu_stride = 0 ;
//            if (GJP.XnodeBc()==BND_CND_APRD) coeff = -coeff ;
            if (GJP.NodeBc(mu)==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
            vec_plus_mu_stride = f_size_4d - f_site_size_4d ;
          }

//     Float time2 = -dclock();
      sproj_tr[mu]( (IFloat *)&tmp_mat1,
                    (IFloat *)v1_plus_mu,
                    (IFloat *)v2+vec_offset,
                    ls, vec_plus_mu_stride, f_size_4d-f_site_size_4d) ;

      sproj_tr[mu+4]( (IFloat *)&tmp_mat2,
                      (IFloat *)v2_plus_mu,
                      (IFloat *)v1+vec_offset,
                      ls, vec_plus_mu_stride, f_size_4d-f_site_size_4d) ;
 //    time2 += dclock();
 //    print_flops("","sproj",2*9*16*ls,time2);

      tmp_mat1 += tmp_mat2 ;

      // If GJP.Snodes > 1 sum up contributions from all s nodes
      if(GJP.Snodes() != 1) {
	  glb_sum_multi_dir((Float *)&tmp_mat1, 4, sizeof(Matrix)/sizeof(IFloat) ) ;
      }

      tmp_mat2.DotMEqual(*(gauge+gauge_offset), tmp_mat1) ;

      tmp_mat1.Dagger(tmp_mat2) ;

      tmp_mat2.TrLessAntiHermMatrix(tmp_mat1) ;

      tmp_mat2 *= coeff ;

      *(mom+gauge_offset) += tmp_mat2 ;

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
  VRB.Sfree(cname, fname, str_site_v2, site_v2) ;
  qfree(site_v2) ;
 
  VRB.Sfree(cname, fname, str_site_v1, site_v1) ;
  qfree(site_v1) ;

  for(int i =0;i<4;i++) {
    sfree(cname,fname,"v1_buf",v1_buf[i]);
    sfree(cname,fname,"v2_buf",v2_buf[i]);
  }
 
  VRB.Sfree(cname, fname, str_v2, v2) ;
  ffree(v2) ;
 
  VRB.Sfree(cname, fname, str_v1, v1) ;
  ffree(v1) ;
 
  return ;
}
// CJ: change end

CPS_END_NAMESPACE
