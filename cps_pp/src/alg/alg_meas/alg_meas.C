#include<config.h>

#ifdef USE_BFM
#include <util/lattice/fbfm.h>
#endif

#include <alg/no_arg.h>
#include <alg/common_arg.h>
#include <alg/pbp_arg.h>
#include <alg/fix_gauge_arg.h>
#include <alg/w_spect_arg.h>
#include <alg/eig_arg.h>
#include <alg/pot_arg.h>

#include <alg/alg_plaq.h>
#include <alg/alg_pbp.h>
#include <alg/alg_fix_gauge.h>
#include <alg/alg_w_spect.h>
#include <alg/alg_eig.h>
#include <alg/alg_pot.h>

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief PAB... Definitions of the AlgMeas class methods.
  
  $Id: alg_meas.C,v 1.11 2013-04-05 17:46:30 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:46:30 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_meas/alg_meas.C,v 1.11 2013-04-05 17:46:30 chulwoo Exp $
//  $Id: alg_meas.C,v 1.11 2013-04-05 17:46:30 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_meas.C,v $
//  $Revision: 1.11 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_meas/alg_meas.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_meas.C
//
// AlgMeas is a measurement class derived from Alg and it measures 
// whatever the heck you tell it to measure in the MeasTask list.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <alg/alg_meas.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/glb.h>
#include <util/qcdio.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
/*!
  \param latt The lattice on which to compute the plaquette.
  \param c_arg The common argument structure for all algorithms.
  \param arg A dummy argument structure.
 */
//------------------------------------------------------------------
static WspectOutput woe;

AlgMeas::AlgMeas(CommonArg *c_arg,
	         MeasArg *arg) : 
  Alg(LatticeFactory::Create(arg->Fermion,arg->Gluon), c_arg) 
{

  cname = "AlgMeas";
  char *fname = "AlgMeas(L&,CommonArg*,MeasArg*)";
  VRB.Func(cname,fname);

  alg_meas_arg = arg;

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgMeas::~AlgMeas() {
  char *fname = "~AlgMeas()";
  VRB.Func(cname,fname);
  LatticeFactory::Destroy();
}


//------------------------------------------------------------------
//! Performs the computation.
/*!
  Loops over tasks and then calls RunTask.
*/
//------------------------------------------------------------------
void AlgMeas::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  // Get the Lattice
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // Loop over tasks
  for ( int task=0;
	task<alg_meas_arg->TaskList.TaskList_len;
	task++){
    RunTask(&alg_meas_arg->TaskList.TaskList_val[task]);
  }

}

/*Assume:
 * i)  the lattice is loaded
 * ii) the RNGs are loaded
 */
void AlgMeas::RunTask(MeasTask *task)
{
  char output_file[256];
  char output_tmp[256];
  char meta_file[256];
  char cwd[256];

  sprintf(output_file,"%s.%d",task->OutputFilestem,alg_meas_arg->TrajCur);
  CommonArg ca;
  ca.set_filename(output_file);

  char *output_directory;

    /*
     * For wilson spectrum the output filenames were mangled by someone
     * so we change directory to the output file stem + Traj (creating) 
     * and let all the stupidly named files go there.
     * Take the cwd to get back.
     */
  if ( task->Measurement == MeasAlgWspect ) { 
    strcpy(output_tmp,output_file);
    output_directory = output_tmp;
  } else { 
    TruncateFile(output_file);
    strcpy(output_tmp,output_file);
    output_directory = Dirname(output_tmp);
    sprintf(meta_file,"%s.arg.vml.%d",output_file,alg_meas_arg->TrajCur);
  }

  /* 
   * Outputs the DoArg.TrajCur, and MeasArg.TrajCur and MeasTask.TrajCur for posterity
   */

  switch( task->Measurement ) {
    
  case MeasAlgPlaq:
    {
      NoArg na; 
      AlgPlaq plaq(AlgLattice(),&ca,&na);
      plaq.run();
      Document(output_directory,task);
    }
    break;
  case MeasAlgPbp:
    {

      PbpArg pa; 
      if ( ! pa.Decode(task->ArgFilename,"dummy") ) exit(-1);

      AlgPbp pbp(AlgLattice(),&ca,&pa);
      pbp.run();
      if ( UniqueID() == 0 )
	pa.Encode(meta_file,"PbpArg");
      Document(output_directory,task);
    }
    break;
  case MeasAlgWspect:
    {

      WspectArg wa; 
      CgArg cg; 
      if ( ! wa.Decode(task->ArgFilename,"dummy") ) exit(-1);
      if ( ! woe.Decode(wa.WspectOutputFile,"dummy") ) exit(-1);
      if ( ! cg.Decode(wa.CgArgFile,"dummy") ) exit (-1);
      ca.results = (void *) &woe;
      AlgWspect ws(AlgLattice(),&ca,&wa,&cg,1);

      /*
       * Wspect doesnt respect output names so use a nifty little "chdir" and back again to work
       * around this.
       */
      mkdir(output_directory,0777);
      chdir(output_directory);
      TruncateWspectFiles();
      Document("./",task);
      sprintf(meta_file,"%s/w_spect_arg.%d",output_directory,alg_meas_arg->TrajCur);
      ws.run();
      if ( UniqueID() == 0 ) {
	wa.Encode(meta_file,"WspectArg");
	sprintf(meta_file,"%s/cg_arg.%d",output_directory,alg_meas_arg->TrajCur);
	cg.Encode(meta_file,"CgArg");
	sprintf(meta_file,"%s/w_spect_output.%d",output_directory,alg_meas_arg->TrajCur);
	woe.Encode(meta_file,"WspectOutput");
      }
      chdir(alg_meas_arg->WorkDirectory);
    }
    break;
  case MeasAlgEig:
    {
      EigArg ea; 
      if ( ! ea.Decode(task->ArgFilename,"dummy") ) exit(-1);
      
      AlgEig eig(AlgLattice(),&ca,&ea);
      eig.run();
      if ( UniqueID() == 0 )
	ea.Encode(meta_file,"EigArg");
      Document(output_directory,task);
    }
    break;
  case MeasAlgPot:
    {
      PotArg po; 
      if (! po.Decode(task->ArgFilename,"dummy") ) exit (-1);

      AlgPot poo(AlgLattice(),&ca,&po);
      poo.run();
      if ( UniqueID() == 0 )
	po.Encode(meta_file,"PotArg");
      Document(output_directory,task);
    }
    break;
  case MeasAlgFixGauge:
    {
      FixGaugeArg fga; 
      if ( ! fga.Decode(task->ArgFilename,"dummy") ) exit(-1);
      AlgFixGauge fg(AlgLattice(),&ca,&fga);
      fg.run();
      if ( UniqueID() == 0 )
	fga.Encode(meta_file,"FixGaugeArg");
      Document(output_directory,task);
    }
    break;
  case MeasAlgFixGaugeFree:
    AlgLattice().FixGaugeFree();
    break;
  default: 
    break;
  }

}
void AlgMeas::Document(char * output_directory,MeasTask *task)
{
  char filename[256];
  sprintf(filename,"%s/meas_arg.%d", output_directory, alg_meas_arg->TrajCur);

  alg_meas_arg->Encode(filename,"meas_arg");

  sprintf(filename,"%s/meas_task.%d", output_directory, alg_meas_arg->TrajCur);
  task->Encode(filename,"meas_task");

  sprintf(filename,"%s/do_arg.%d", output_directory, alg_meas_arg->TrajCur);
  GJP.GetDoArg()->Encode(filename,"do_arg");

  return;
}

Lattice *LatticeFactory::lat_p;

void LatticeFactory::Destroy(void)
{
  delete lat_p;
}

Lattice & LatticeFactory::Create(FclassType fermion,GclassType gluon)
{
  /* BFM VALENCE ANALYSIS */
#ifdef USE_BFM
  if ( (fermion == F_CLASS_BFM) && (gluon == G_CLASS_NONE ) ) {
    lat_p = new GnoneFbfm;
    return *lat_p; 
  }
  if ( (fermion == F_CLASS_BFM) && (gluon == G_CLASS_IMPR_RECT ) ) {
    lat_p = new GimprRectFbfm;
    return *lat_p;
  }
#endif

  /* DOMAIN WALL VALENCE ANALYSIS */
  if ( (fermion == F_CLASS_DWF) && (gluon == G_CLASS_NONE ) ) {
    lat_p = new GnoneFdwf ;
    return *lat_p; 
  }
  if ( (fermion == F_CLASS_DWF) && (gluon == G_CLASS_WILSON ) ) {
    lat_p = new GwilsonFdwf;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_DWF) && (gluon == G_CLASS_IMPR_RECT ) ) {
    lat_p = new GimprRectFdwf ;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_MOBIUS) && (gluon == G_CLASS_NONE ) ) {
    lat_p = new GnoneFmobius ;
    return *lat_p; 
  }
  if ( (fermion == F_CLASS_MOBIUS) && (gluon == G_CLASS_WILSON ) ) {
    lat_p = new GwilsonFmobius;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_MOBIUS) && (gluon == G_CLASS_IMPR_RECT ) ) {
    lat_p = new GimprRectFmobius ;
    return *lat_p;
  }

  /* WILSON VALENCE ANALYSIS */
  if ( (fermion == F_CLASS_NAIVE) && (gluon == G_CLASS_NONE ) ) {
    lat_p = new GnoneFnaive ;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_WILSON) && (gluon == G_CLASS_NONE ) ) {
    lat_p = new GnoneFwilson ;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_WILSON) && (gluon == G_CLASS_WILSON ) ) {
    lat_p = new GwilsonFwilson;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_WILSON) && (gluon == G_CLASS_IMPR_RECT ) ) {
    lat_p =  new GimprRectFwilson;
    return *lat_p;
  }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// added for twisted-mass wilson fermions
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /* twisted-mass WILSON VALENCE ANALYSIS */
  if ( (fermion == F_CLASS_WILSON_TM) && (gluon == G_CLASS_NONE ) ) {
    lat_p = new GnoneFwilsonTm ;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_WILSON_TM) && (gluon == G_CLASS_WILSON ) ) {
    lat_p = new GwilsonFwilsonTm;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_WILSON_TM) && (gluon == G_CLASS_IMPR_RECT ) ) {
    lat_p =  new GimprRectFwilsonTm;
    return *lat_p;
  }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /* CLOVER VALENCE ANALYSIS */
  if ( (fermion == F_CLASS_CLOVER) && (gluon == G_CLASS_NONE ) ) {
    lat_p = new GnoneFclover;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_CLOVER) && (gluon == G_CLASS_WILSON ) ) {
    lat_p = new GwilsonFclover;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_CLOVER) && (gluon == G_CLASS_IMPR_RECT ) ) {
    lat_p = new GimprRectFclover;
    return *lat_p;
  }

  /* STAGGERED VALENCE ANALYSIS */
  if ( (fermion == F_CLASS_STAG) && (gluon == G_CLASS_NONE ) ) {
    lat_p = new GnoneFstag;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_STAG) && (gluon == G_CLASS_WILSON ) ) {
    lat_p = new GwilsonFstag;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_STAG) && (gluon == G_CLASS_IMPR_RECT ) ) {
    lat_p = new GimprRectFstag;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_STAG) && (gluon == G_CLASS_IMPR_OLSYM ) ) {
    lat_p = new GimprOLSymFstag;
    return *lat_p;
  }

  /* ASQTAD STAGGERED VALENCE ANALYSIS */
  if ( (fermion == F_CLASS_ASQTAD) && (gluon == G_CLASS_NONE ) ) {
    lat_p = new GnoneFasqtad;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_ASQTAD) && (gluon == G_CLASS_WILSON ) ) {
    lat_p = new GwilsonFasqtad;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_ASQTAD) && (gluon == G_CLASS_IMPR_RECT ) ) {
    lat_p = new GimprRectFasqtad;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_ASQTAD) && (gluon == G_CLASS_IMPR_OLSYM ) ) {
    lat_p = new GimprOLSymFasqtad;
    return *lat_p;
  }

  // Only enabled for QCDOC target currently
  /* P4 STAGGERED VALENCE ANALYSIS */
  if ( (fermion == F_CLASS_P4) && (gluon == G_CLASS_NONE ) ) {
    lat_p = new GnoneFp4;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_P4) && (gluon == G_CLASS_WILSON ) ) {
    lat_p = new GwilsonFp4;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_P4) && (gluon == G_CLASS_IMPR_RECT ) ) {
    lat_p = new GimprRectFp4;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_P4) && (gluon == G_CLASS_IMPR_OLSYM ) ) {
    lat_p = new GimprOLSymFp4;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_P4) && (gluon == G_CLASS_TADPOLE_RECT ) ) {
    lat_p = new GtadpoleRectFp4;
    return *lat_p;
  }

  /* F_NONE VALENCE ANALYSIS */
  if ( (fermion == F_CLASS_NONE) && (gluon == G_CLASS_NONE ) ) {
    lat_p = new GnoneFnone;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_NONE) && (gluon == G_CLASS_WILSON ) ) {
    lat_p = new GwilsonFnone;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_NONE) && (gluon == G_CLASS_IMPR_RECT ) ) {
    lat_p = new GimprRectFnone;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_NONE) && (gluon == G_CLASS_IMPR_OLSYM ) ) {
    lat_p = new GimprOLSymFnone;
    return *lat_p;
  }
  if ( (fermion == F_CLASS_NONE) && (gluon == G_CLASS_TADPOLE_RECT ) ) {
    lat_p = new GtadpoleRectFnone;
    return *lat_p;
  }

  ERR.General("LatticeFactory","Create()",
	      "Lattice type (fermion = %d, gluon = %d) not defined\n",
	      fermion, gluon);

}
char *AlgMeas::Dirname (char *path)
{
  static const char dot[] = ".";
  char *last_slash;
 
  /* Find last '/'.  */
  last_slash = path != NULL ? strrchr (path, '/') : NULL;
 
  if (last_slash == path)
    /* The last slash is the first character in the string.  We have to
       return "/".  */
    ++last_slash;
  else if (last_slash != NULL && last_slash[1] == '\0')
    /* The '/' is the last character, we have to look further.  */
    last_slash = (char *) memchr (path, last_slash - path, '/');
 
  if (last_slash != NULL)
    /* Terminate the path.  */
    last_slash[0] = '\0';
  else
    /* This assignment is ill-designed but the XPG specs require to
       return a string containing "." in any case no directory part is
       found and so a static and constant string is required.  */
    path = (char *) dot;
 
  return path;
}
void AlgMeas::TruncateWspectFiles(void)
{
  TruncateFile(woe.cg);
  TruncateFile(woe.cg2);
  TruncateFile(woe.pbp);
  TruncateFile(woe.mid_point);
  TruncateFile(woe.a0_p);
  TruncateFile(woe.a1);
  TruncateFile(woe.b1);
  TruncateFile(woe.pion);
  TruncateFile(woe.pion_prime);
  TruncateFile(woe.rho);
#if 0
  TruncateFile(woe.a0);
  TruncateFile(woe.a0_prime);
  TruncateFile(woe.a1_x);
  TruncateFile(woe.a1_y);
  TruncateFile(woe.a1_z);
  TruncateFile(woe.b1_x);
  TruncateFile(woe.b1_y);
  TruncateFile(woe.b1_z);
  TruncateFile(woe.rho_x);
  TruncateFile(woe.rho_y);
  TruncateFile(woe.rho_z);
  TruncateFile(woe.rho_x_prime);
  TruncateFile(woe.rho_y_prime);
  TruncateFile(woe.rho_z_prime);
#else
  TruncateFile(woe.meson_name00);
  TruncateFile(woe.meson_name01);
  TruncateFile(woe.meson_name02);
  TruncateFile(woe.meson_name03);
  TruncateFile(woe.meson_name04);
  TruncateFile(woe.meson_name05);
  TruncateFile(woe.meson_name06);
  TruncateFile(woe.meson_name07);
  TruncateFile(woe.meson_name08);
  TruncateFile(woe.meson_name09);
  TruncateFile(woe.meson_name10);
  TruncateFile(woe.meson_name11);
  TruncateFile(woe.meson_name12);
  TruncateFile(woe.meson_name13);
  TruncateFile(woe.meson_name14);
  TruncateFile(woe.meson_name15);
#endif
  TruncateFile(woe.nucleon);
  TruncateFile(woe.nucleon_prime);
  TruncateFile(woe.delta_x);
  TruncateFile(woe.delta_y);
  TruncateFile(woe.delta_z);
  TruncateFile(woe.delta_t);
  /*
   * Hack to work around a hack.... Grr...
   */
  char *meson_names[] = 
    {
      "a0.dat" , 
      "rho_x.dat" , 
      "rho_y.dat"    ,
      "b1_z.dat"     ,
      "rho_z.dat"    ,
      "b1_y.dat",
      "b1_x.dat",
      "pion_prime.dat",
      "a0_prime.dat",
      "rho_x_prime.dat",
      "rho_y_prime.dat",
      "a1_z.dat" ,
      "rho_z_prime.dat" , 
      "a1_y.dat" , 
      "a1_x.dat" , 
      "pion.dat" 
    };
  for ( int i=0;i<16;i++ ) { 
    TruncateFile(meson_names[i]);
  }
}
void AlgMeas::TruncateFile(char *foo)
{
  if ( UniqueID() == 0 ) {
    FILE *fp = fopen(foo,"w");
    if ( fp ) fclose(fp);
  }
  sync();
}


CPS_END_NAMESPACE
