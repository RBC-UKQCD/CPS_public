#include<config.h>
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
 
  char *str_v0 = "v[0]" ;
  Vector *v[2];
  v[0]= (Vector *)smalloc (cname, fname, str_v0, f_size*sizeof(Float)) ;

  char *str_v1 = "v[1]" ;
  v[1] = (Vector *)smalloc (cname, fname, str_v1, f_size*sizeof(Float)) ;


  //----------------------------------------------------------------
  // Calculate v[0], v[1]. Both v[0], v[1] must be in CANONICAL order after
  // the calculation.
  //----------------------------------------------------------------  

  VRB.Clock(cname, fname, "Before calc force vecs.\n") ;

  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;

    DiracOpDwf dwf(*this, v[0], v[1], &cg_arg, CNV_FRM_YES) ;
    dwf.CalcHmdForceVecs(chi) ;
  }

#ifdef PROFILE
  Float time = -dclock();
  ForceFlops=0;
#endif
  Matrix tmp_mat1, tmp_mat2 ;
  VRB.Clock(cname, fname, "Before loop over links.\n") ;

  int x[4];
  int xpm[4];
  int mu,i;
  int offset_4d_x;
  int offset_4d_xpm;
  int vol=GJP.VolNodeSites();
  int ls=node_sites[4];
  int v_stride=f_site_size_4d*(vol-1);
  int v_offset,vpm_offset;
  IFloat *v_p,*w_p;
  Matrix m;
  Matrix **out;
  BndCndType bc[4]={GJP.XnodeBc(),GJP.YnodeBc(),GJP.ZnodeBc(),GJP.TnodeBc()};
  Float coeff0 = -2 * step_size;
  Float coeff;

  //-----------------------------------------------------------
  //allocate memory for color matrix arrays
  //-----------------------------------------------------------
  out = (Matrix **)smalloc (cname,fname,"out",4*sizeof(int));

  for(i=0;i<4;i++){
    out[i]=(Matrix *)smalloc(cname,fname,"out[i]",18*vol*sizeof(Float));
  }

  for(mu=0;mu<4;mu++) 
    for(int n=0;n<vol;n++) out[mu][n].ZeroMatrix();  
  
  //------------------------------------------
  //singlenode version
  //------------------------------------------
  
  for(mu=0;mu<4;mu++)
    for(xpm[3]=x[3]=0;x[3]<node_sites[3];xpm[3]++,x[3]++)
      for(xpm[2]=x[2]=0;x[2]<node_sites[2];xpm[2]++,x[2]++)
	for(xpm[1]=x[1]=0;x[1]<node_sites[1];xpm[1]++,x[1]++)
	  for(xpm[0]=x[0]=0;x[0]<node_sites[0];xpm[0]++,x[0]++){
	    xpm[mu] = (x[mu]+1)%node_sites[mu];
	    offset_4d_xpm = xpm[0]+node_sites[0]*(xpm[1]+node_sites[1]*(xpm[2]+
									node_sites[2]*xpm[3]));
	    offset_4d_x   = x[0]+node_sites[0]*(x[1]+node_sites[1]*(x[2]+
								    node_sites[2]*x[3]));
	    v_offset = 4*offset_4d_x;
	    vpm_offset = 4*offset_4d_xpm;
	    
	    v_p = (IFloat *)(v[0]+vpm_offset);
	    w_p = (IFloat *)(v[1]+v_offset);
	    
	    sproj_tr[mu]((IFloat *)&m,v_p,w_p,ls,v_stride,v_stride);
	    out[mu][offset_4d_x] += m;
	    
	    v_p = (IFloat *)(v[0]+v_offset);
	    w_p = (IFloat *)(v[1]+vpm_offset);
	    
	    sproj_tr[mu+4]((IFloat *)&m,w_p,v_p,ls,v_stride,v_stride);
	    out[mu][offset_4d_x] += m;
	    
	  }
  
  
  //--------------------------------------------
  //multiply the above matrix out[mu][x] by gauge field
  //and update mom
  //-------------------------------------------
  int gauge_offset;
  for(int mu=0;mu<4;mu++)
    
    for(x[3]=0;x[3]<node_sites[3];x[3]++)
      for(x[2]=0;x[2]<node_sites[2];x[2]++)
	for(x[1]=0;x[1]<node_sites[1];x[1]++)
	  for(x[0]=0;x[0]<node_sites[0];x[0]++){
	    
	    coeff=coeff0;
	    if(x[mu]==node_sites[mu]-1&&bc[mu]==BND_CND_APRD) coeff = -coeff;

	    offset_4d_x  = x[0]+node_sites[0]*(x[1]+node_sites[1]*(x[2]+
								     node_sites[2]*x[3]));
	    
	    //	    if(GJP.Snodes()>1) 
	    //   glb_sum_multi_dir((Float *)(out[mu]+offset_4d_x),4,sizeof(Matrix)/sizeof(IFloat));
	    gauge_offset = 4*offset_4d_x+mu;
	    tmp_mat2.DotMEqual(*(gauge+gauge_offset),out[mu][offset_4d_x]);
	    
	    tmp_mat1.Dagger(tmp_mat2);
	    tmp_mat2.TrLessAntiHermMatrix(tmp_mat1);

	    tmp_mat2 *= coeff;
	    
	    *(mom+gauge_offset) += tmp_mat2;

	  }
 ForceFlops += (2*9*16*ls + 18+ 198+36+24)*vol*4;	  
 
#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif
 
//------------------------------------------------------------------
// deallocate smalloc'd space
//------------------------------------------------------------------

  for(i=0;i<4;i++){
    VRB.Sfree(cname,fname,"out[i]",out[i]);
    sfree(out[i]);    
  }
  
  VRB.Sfree(cname,fname,"out",out);
  sfree(out);
  
  VRB.Sfree(cname, fname, str_v1, v[1]) ;
  sfree(v[1]) ;
  
  VRB.Sfree(cname, fname, str_v0, v[0]) ;
  sfree(v[0]) ;
	      
  return ;
  
}
CPS_END_NAMESPACE
