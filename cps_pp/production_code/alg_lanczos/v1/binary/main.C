
#include <conf.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <alg/alg_plaq.h>
#include <alg/qpropw.h>
#include <alg/alg_lanczos.h>
#include <alg/alg_meas.h>
#include <alg/meas_arg.h>   
#include <alg/meson.h>
#include <alg/alg_meas.h>   
#include <alg/nuc2pt.h>
#include <alg/alg_fix_gauge.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/command_line.h>
#include <util/ReadLatticePar.h>
#include <util/qcdio.h>
#include <util/site.h>

#include <alg/matrixpolynomial_arg.h>

#include <util/eigen_container.h>

#include <util/time_cps.h>

USING_NAMESPACE_CPS
using namespace std;

void CalMesons(QPropW &q, char *out_dir, int traj, char *src_str="POINT");
void CalNucleons(QPropW &q, char *out_dir, int traj, char *src_str="POINT");
void ReadGaugeField(const MeasArg &meas_arg);
void WriteGaugeField(const MeasArg &meas_arg);


#if TARGET == BGL
extern "C"
int omp_get_thread_num(){return 1;}
#endif

//------------------------------------------
// should go to eigen_container.C later, someday 
//------------------------------------------
// needed to declare globally
std::vector<EigenCache*> cps::EigenCacheList(0);

//Search contents that match to arguments, return 0 if not found
EigenCache* cps::EigenCacheListSearch( char* fname_root_bc, int neig )
{
  EigenCache* ecache=0;

  for(int i=0; i< EigenCacheList.size(); ++i  )
    if( EigenCacheList[i]-> is_cached( fname_root_bc, neig ) )
      ecache = EigenCacheList[i];

  return ecache;
}


// Cleanup list EigenCache, it also destroies contents pointed by the elements.
void cps::EigenCacheListCleanup( )
{
  for(size_t i=0;i< EigenCacheList. size(); ++i){
    EigenCacheList[i]-> dealloc(); 
  }
  EigenCacheList.clear() ;
}

//----------------------------------------------------

int main(int argc,char *argv[])
{

  DoArg do_arg;
  NoArg no_arg;
  CommonArg common_arg;
  LanczosArg eig_arg;
  QPropWArg qpropw_arg;
  QPropWBoxArg qpropw_box_arg;
  MatrixPolynomialArg cheby_arg;
  MeasArg meas_arg;
  
  char *cname = argv[0] ;
  char *fname = "main()" ;
  //char* filename;
  //filename = (char*) smalloc( 128*sizeof(char) );
  
  Start(&argc, &argv);
  print_asctime("grand start");
  
  //------------------------------
  //set log directory
  //------------------------------
  //char log_dir[255];
  //sprintf(&log_dir[0],"IOLOG");
  
  if ( argc!=8) { 
    if(!UniqueID())printf("(exe) work-directory LMA_SHIFT[0] LMA_SHIFT[1] LMA_SHIFT[2] LMA_SHIFT[3] do_cg do_eigv_read (-qmp-geom ....)\n");
    exit(-1);
  }
  
  // load defaults
  
  // move to working sub-directory (first arg on command line after exec.)
  // should be relative to /host.../username
  chdir(argv[1]);

  int LMA_SHIFT[4];
  for(int i=0;i<4;++i) LMA_SHIFT[i] = atoi(argv[2+i]);

  int do_cg = atoi(argv[6]);
  
  int do_eigenv_read= atoi(argv[7]);
  
  if(!meas_arg.Decode("meas_arg.vml", "meas_arg"))
    { VRB.Flow(cname, fname, "Can't open meas_arg.vml!\n"); exit(-1);}
  
  if ( !do_arg.Decode("do_arg.vml","do_arg") ) 
    { 
      printf("Decoding of do_arg failed\n"); exit(-1);
    }

  if ( !eig_arg.Decode("lanczos_arg.vml","lanczos_arg") ) 
    { 
      printf("Decoding of lanczos_arg failed\n"); exit(-1);
    }

  if ( !cheby_arg.Decode("cheby_arg.vml","cheby_arg") ) 
    { 
      printf("Decoding of cheby_arg failed\n"); exit(-1);
    }
  eig_arg.matpoly_arg = (Float*)&cheby_arg; // ugly thing to workaround vml 
  
  
  VRB.Level(do_arg.verbose_level);
  GJP.Initialize(do_arg);
  

  int do_write_lattice(0);
  if(do_write_lattice)
  {
    meas_arg. TrajCur = meas_arg.TrajStart;
    WriteGaugeField(meas_arg);
    End();
    return 0;
  }

  
  for(meas_arg.TrajCur=meas_arg.TrajStart;
      meas_arg.TrajCur<=meas_arg.TrajLessThanLimit;
      meas_arg.TrajCur+=meas_arg.TrajIncrement) {
    

    int traj = meas_arg. TrajCur;
    ReadGaugeField(meas_arg);
    GimprRectFdwf lattice;
    
    common_arg.set_filename("plaq.dat");
    AlgPlaq plaq(lattice,&common_arg,&no_arg);
    plaq.run();
  

    /*
     *  Compute eigen vectors and values
     *
     * Set nk_lanczos_vector to 0 in the qpropw.vml file, 
     * if the eigenvectors are already computed and saved.
     */
    
    if( eig_arg.nk_lanczos_vectors > 0) {
      print_asctime("eigenvector comp start");
      char fname[1024];
      snprintf(fname, 1024,"%s/eig4dee.mass%g.traj%04d",
	       eig_arg.file, eig_arg.mass, traj );
      eig_arg.file = fname;

      char resname[1024];
      snprintf(resname,1024,"%s.%04d", eig_arg.results, traj);
      eig_arg.results = resname;
      
      
      int init_flag = 0; // 0 : wall start, 1 :  compression, 2 : compress & decompress, 3: refinement using Ritz
      int ncompress = 10; // number of evecs in the compressed linear combination
      
      //char* comp_file="/scratch1/izubuchi/cps_eigen/DW_b2.13_16x32_ms0.032-mu0.01/eig4dee.mass0.01.traj2200";
      
      char* comp_file;
      
      if( init_flag == 1){
	eig_arg. nk_lanczos_vectors = ncompress;
	eig_arg. np_lanczos_vectors = ncompress;
      }
      AlgLanczos  eig(lattice, &common_arg, &eig_arg);
      eig.run( init_flag, ncompress, comp_file );
    }
    print_asctime("eigenvector comp ends");
    
    
    /*
     *  First read in the eigenvectors and cache them
     *
     *  This is not absolutely necessary, but it's good to have checked
     *  make sure the consistency between eigenv in the file, ensemble loaded, and parameters set
     *
  *
  *  if you use both periodic and anti-periodic boundary conditions you should check both.
  *

   *  Dear Chulwoo:
        I couldn't do this automatically due to a unique scope rule that dirac operator has
      since sometimes I am calling nev_check from dirac (such as inv_lowmodeapprox.C) and that
      I don't know how to use the already instanciated dirac operator.
   */

  qpropw_arg.Decode("qpropw_arg.vml","qpropw_arg");
  qpropw_arg.Encode("qpropw_arg.dat","qpropw_arg");


  
  if(do_eigenv_read)
  {
    print_asctime("eigenv_read start");
    if( GJP.Snodes() != 1) ERR.NotImplemented(cname,fname,"currently only doing I/O for local Ls");
     
    int neig=qpropw_arg.cg.neig;

    const int n_fields =  GJP.SnodeSites();  //   *nk ; 
    const int f_size_per_site = lattice.FsiteSize() / GJP.SnodeSites()   / 2;
    
    // root name of eigen value/vector to be read 
    char efname[1024];
    snprintf(efname, 1024,"%s/eig4dee.mass%g.traj%04d.bc%d%d%d%d",
	     qpropw_arg.cg.fname_eigen, qpropw_arg.cg.mass, traj,
	     GJP.Bc(0),GJP.Bc(1),GJP.Bc(2),GJP.Bc(3));


    // search for eigen cache
    EigenCache* ecache;
    if( (ecache = EigenCacheListSearch( efname, neig ) ) == 0 ){
      ecache = new EigenCache();  
      
    EigenCacheList. push_back( ecache );
  }

  EigenContainer eigcon( lattice, efname, neig, f_size_per_site, n_fields, ecache );


    Float* eval = eigcon.load_eval();
    
    for(int iev=0; iev< neig;++iev){
      Vector* evec= eigcon.nev_load( iev );

      if(iev % 20 == 0){  // Let's check every 20
      //if(0){  // Let's not check

	Float res_chk, eval_chk;
	Float mass = qpropw_arg.cg.mass;  	
	

	eigcon. nev_check( evec, mass, &res_chk, &eval_chk );

	Float ev = *(eval+iev);
	if( fabs(res_chk) > 1e-5 ) 
	  //VRB.Warn
	  ERR.General(cname,fname,"nev_check index %d eigval %g mass %e res %e > 1e-5",
			   iev, ev, mass, res_chk);

	if( fabs(eval_chk- ev) > 1e-5 ) 
	  // VRB.Warn
	  ERR.General(cname,fname,"nev_check index %d mass %e eval_chk %e eval %e, abs_err %e > 1e-5", 
		      iev, mass, eval_chk, ev, fabs(eval_chk-ev) );

      }
    }
    print_asctime("eigenv_read ends");
  }
  

  
  /*
   *  Measurement example
   */

  
  char qprop_dir[256];
  snprintf(qprop_dir,256, "%s",qpropw_arg.file);  //really just the directory, so save it, and modify qpropw_arg.file below
  
  char eigen_dir[256];
  snprintf(eigen_dir,256, "%s", qpropw_arg.cg.fname_eigen);

  /*
   *  CG computation
      Either the full CG or CG with the low-mode deflation
   */

  //int  do_cg = 1; 
  if( do_cg )
  { 
    print_asctime("full cg starts");
    // set solver 
    //qpropw_arg. cg. Inverter = CG;
    qpropw_arg. cg. Inverter = CG_LOWMODE_DEFL;

    // root name of eigen value/vector to be read 
    char efname[1024];
    snprintf(efname, 1024,"%s/eig4dee.mass%g.traj%04d",
	     eigen_dir, qpropw_arg.cg.mass, traj );
    qpropw_arg.cg.fname_eigen = efname;

    // name for propagators
    char qfname[1024];
    snprintf(qfname,1024,
	     "%sprop.point-src_%02d_%02d_%02d_%02d.mass%g.traj%04d",
	     qprop_dir,
	     qpropw_arg.x, qpropw_arg.y, qpropw_arg.z, qpropw_arg.t,
	     qpropw_arg.cg.mass, traj);  
    qpropw_arg.file = qfname;
    
    CommonArg c_arg_qprop;      
    c_arg_qprop.set_filename( "qpropw_fullcg.dat" );

    QPropWPointSrc prop(lattice, &qpropw_arg, &c_arg_qprop);
    prop.Run();
    
    print_asctime("full prop.Run ends");


    char src_str[1024];
    sprintf(src_str,"point_mass%g_fullcg_src%02d_%02d_%02d_%02d", 
	    qpropw_arg.cg.mass,
	    qpropw_arg.x, qpropw_arg.y, qpropw_arg.z, qpropw_arg.t );
    if(!UniqueID()) printf("doing %s\n", src_str);
    CalMesons(prop, ".", traj, src_str);
    CalNucleons(prop, ".", traj, src_str);
    print_asctime("full hadrons end");
  }

  /*
   *  Now do the LMA(Low Mode Approximation)
   */

  // number of shifts in each direction
  int *NUM_SHFT = LMA_SHIFT;
  int nshft = NUM_SHFT[0]*NUM_SHFT[1]*NUM_SHFT[2]*NUM_SHFT[3];

  
  // original source location
  int orig_src[4]= { qpropw_arg.x,qpropw_arg.y,qpropw_arg.z,qpropw_arg.t };

  int sft[4];
  for(int isft_=0; isft_< nshft;isft_++){
    print_asctime("LMA %d / %d starts",isft_,nshft);

    // reset the source location
    qpropw_arg. x =  orig_src[0];
    qpropw_arg. y =  orig_src[1];
    qpropw_arg. z =  orig_src[2];
    qpropw_arg. t =  orig_src[3];


    // compute the shift steps
    int isft = isft_;
    sft[0] = isft % NUM_SHFT[0]; isft /= NUM_SHFT[0];
    sft[1] = isft % NUM_SHFT[1]; isft /= NUM_SHFT[1];
    sft[2] = isft % NUM_SHFT[2]; isft /= NUM_SHFT[2];
    sft[3] = isft % NUM_SHFT[3]; isft /= NUM_SHFT[3];
    
    // do shifts
    qpropw_arg. x += (GJP.Sites(0)/ NUM_SHFT[0]) *sft[0]; 
    qpropw_arg. y += (GJP.Sites(1)/ NUM_SHFT[1]) *sft[1];
    qpropw_arg. z += (GJP.Sites(2)/ NUM_SHFT[2]) *sft[2];
    qpropw_arg. t += (GJP.Sites(3)/ NUM_SHFT[3]) *sft[3];
    
    // make sure this is within the nodes
    qpropw_arg. x %= GJP.Sites(0);
    qpropw_arg. y %= GJP.Sites(1);
    qpropw_arg. z %= GJP.Sites(2);
    qpropw_arg. t %= GJP.Sites(3);

    // root name of eigen value/vector to be read 
    char efname[1024];
    snprintf(efname, 1024,"%s/eig4dee.mass%g.traj%04d",
	     eigen_dir, qpropw_arg.cg.mass, traj );
    qpropw_arg.cg.fname_eigen = efname;
    
    //-----------------------------
    // Low Mode Approximation part
    //-----------------------------
    {
      // set inverter to be  LMA
      qpropw_arg. cg. Inverter = LOWMODEAPPROX;

      // name for propagators to save (if you want to save)
      char qfname[1024];
      snprintf(qfname,1024, "%sprop_mass%g_lma%d_src_%02d_%02d_%02d_%02d",
	       qprop_dir,
	       qpropw_arg.cg.mass,  
	       qpropw_arg.cg.neig,
	       qpropw_arg.x, qpropw_arg.y,qpropw_arg.z, qpropw_arg.t);
      
      qpropw_arg.file = qfname;

      CommonArg c_arg_qprop;      
      c_arg_qprop.set_filename( "qpropw_lma.dat" );
      qpropw_arg. save_prop = 0; // turned off
      
      QPropWPointSrc prop(lattice, &qpropw_arg, &c_arg_qprop);
      prop.Run();
      print_asctime("LMA %d / %d prop.Run() ends",isft_,nshft);
      
      char src_str[1024];
      sprintf(src_str,"point_mass%g_lma%d_src%02d_%02d_%02d_%02d", 
	      qpropw_arg.cg.mass,
	      qpropw_arg.cg.neig,
	      qpropw_arg.x, qpropw_arg.y, qpropw_arg.z, qpropw_arg.t );
      
      if(!UniqueID()) printf("doing LMA %s\n", src_str);
      CalMesons(prop, ".", traj, src_str);
      CalNucleons(prop, ".", traj, src_str);
      print_asctime("LMA %d / %d hadrons end",isft_,nshft);
    }

    //----------------------------------
    // The All Mode Approximation part
    //----------------------------------
    if (qpropw_arg. cg. ama_stop_rsd > 0.0 )
    {
      print_asctime("AMA %d / %d starts",isft,nshft);
      
      // set inverter to be  CG_LOWMODE_DEFLA
      qpropw_arg. cg. Inverter = CG_LOWMODE_DEFL;
      // and relax the stopping criteria in ama_stop_rsd
      Float stop_rsd_save = qpropw_arg. cg. stop_rsd;
      qpropw_arg. cg. stop_rsd = qpropw_arg. cg. ama_stop_rsd;
      
      // name for propagators to save (if you want to save)
      char qfname[1024];
      snprintf(qfname,1024, "%sprop_mass%g_ama%d_src_%02d_%02d_%02d_%02d",
	       qprop_dir,
	       qpropw_arg.cg.mass,  
	       qpropw_arg.cg.neig,
	       qpropw_arg.x, qpropw_arg.y,qpropw_arg.z, qpropw_arg.t);
      
      qpropw_arg.file = qfname;

      CommonArg c_arg_qprop;      
      c_arg_qprop.set_filename( "qpropw_ama.dat" );
      qpropw_arg. save_prop = 0; // turned off
      
      QPropWPointSrc prop(lattice, &qpropw_arg, &c_arg_qprop);
      prop.Run();
      print_asctime("AMA %d / %d prop.Run() ends",isft,nshft);
    
      // recover the original residue
      qpropw_arg. cg. stop_rsd = stop_rsd_save;
      
      char src_str[1024];
      sprintf(src_str,"point_mass%g_ama%d_src%02d_%02d_%02d_%02d", 
	      qpropw_arg.cg.mass,
	      qpropw_arg.cg.neig,
	      qpropw_arg.x, qpropw_arg.y, qpropw_arg.z, qpropw_arg.t );
      
      if(!UniqueID()) printf("doing AMA %s\n", src_str);
      CalMesons(prop, ".", traj, src_str);
      CalNucleons(prop, ".", traj, src_str);
    print_asctime("AMA %d / %d hadrons end",isft_,nshft);
    }


  }//loop(isft)

  EigenCacheListCleanup();
  }// trajectory loop
  

  print_asctime("grand finale, Yeah!");
  End();
  
  return 0;
} 


void CalNucleons(QPropW &q, char *out_dir, int traj, char *src_str){
  Nuc2pt nuc(NUC_G5C, POINT);
  nuc.calcNucleon( q );

  char file[256];
  sprintf(file, "%s/nucleon_%s.dat.%d", out_dir,src_str,traj);
  FILE *fp = Fopen(&file[0],"w");

  {
    Nuc2pt nuc(NUC_G5C, POINT);
    nuc.calcNucleon( q );
    nuc.Print(fp);
  }

  Fclose(fp);
  
  
}

void CalMesons(QPropW &q, char *out_dir, int traj, char *src_str){
 
  char file[256];
  sprintf(file, "%s/meson_%s.dat.%d", out_dir,src_str,traj);
  FILE *fp = Fopen(&file[0],"w");
 
  int mu, nu;
  //vector and scalar
  for ( mu = 0; mu < 4; mu++ ){
    char str[256];
    sprintf(str, "GAM_%d", mu+1);
    Meson mes(mu, str);
    mes.Zero();
    mes.setMass(q.Mass(),q.Mass());
    mes.calcMeson(q,q);
    if(!UniqueID())
      mes.Print(fp);
  }

  //pseudoscalar
  {
    mu = -5;
    char str[256];
    sprintf(str, "GAM_5");
    Meson mes(mu, str);
    mes.Zero();
    mes.setMass(q.Mass(),q.Mass());
    mes.calcMeson(q,q);
    if(!UniqueID())
    mes.Print(fp);
  }
  //tensor
  for ( mu = 0; mu < 4; mu++ ){
    for ( nu = mu+1; nu < 4; nu++ ){
      char str[256];
      sprintf(str, "GAM_%d%d",mu+1,nu+1);
      Meson mes1(mu, nu, str);
      mes1.Zero();
      mes1.setMass(q.Mass(),q.Mass());
      mes1.calcMeson(q,q);
      if(!UniqueID())mes1.Print(fp);

      sprintf(str, "GAM_%d%d",nu+1,mu+1);
      Meson mes2(nu, mu, str);
      mes2.Zero();
      mes2.setMass(q.Mass(),q.Mass());
      mes2.calcMeson(q,q);
      if(!UniqueID()) mes2.Print(fp);
    }
  }
  //axialvector
  nu = -5;
  for ( mu = 0; mu < 4; mu++ ){
    char str[256];
    sprintf(str, "GAM_%d5", mu+1);
    Meson mes(mu, nu, str);
    mes.Zero();
    mes.setMass(q.Mass(),q.Mass());
    mes.calcMeson(q,q);
    if(!UniqueID())
    mes.Print(fp);
  }

  // FIXME:  PUT IN THE  A0  SCALAR MESON, or WE WILL MISS A LOT !!
  
  Fclose(fp);

}


void ReadGaugeField(const MeasArg &meas_arg)
{
	char *cname = "main";
	char *fname = "ReadGaugeField";

	GnoneFnone lat;
	char lat_file[256];
	sprintf(lat_file,"%s.%d",meas_arg.GaugeStem,meas_arg.TrajCur);
	VRB.Result(cname, fname, lat_file);
  QioArg rd_arg(lat_file,0.001);//chkprec=0.001
	
	rd_arg.ConcurIONumber=meas_arg.IOconcurrency;
	
	ReadLatticeParallel rl;
	rl.read(lat,rd_arg);
	if(!rl.good()) ERR.General(cname,fname,"Failed read lattice %s",lat_file);	
}

void WriteGaugeField(const MeasArg &meas_arg)
{
	char *cname = "main";
	char *fname = "WriteGaugeField";

	GnoneFnone lat;
	char lat_file[256];
	sprintf(lat_file,"%s.%d",meas_arg.GaugeStem,meas_arg.TrajCur);
	VRB.Result(cname, fname, lat_file);
	QioArg wt_arg(lat_file,0.001);//chkprec=0.001
	
	wt_arg.ConcurIONumber=meas_arg.IOconcurrency;
	
	WriteLatticeParallel wl;
	wl.write(lat,wt_arg);
	if(!wl.good()) ERR.General(cname,fname,"Failed write lattice %s",lat_file);	
}

