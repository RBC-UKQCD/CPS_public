#include <config.h>
#include <util/gjp.h>
#include <util/error.h>
#ifdef USE_QMP
#include <qmp.h>
#endif

#ifdef USE_GRID
#include <util/lattice/fgrid.h>
#endif

#ifdef USE_QUDA
#include <invert_quda.h>
#endif

//QMP definitions for CPS::Start() and CPS::End()

//extern "C" void _mcleanup();
CPS_START_NAMESPACE
#ifdef USE_QMP
void Start(){
  ERR.General("cps","Start()","Start(&argc,&argv should be used with QMP");
//  QMPSCU::init_qmp();
}

void Start(int * argc, char *** argv) {
//  printf("Start(%d %p)\n",*argc,*argv);
  //Initialize QMP
#ifdef USE_GRID
//  Grid::Grid_init(argc,argv);
//  Fgrid::grid_initted=true;
#endif
  QMPSCU::init_qmp(argc, argv);
  GJP.setArg(argc,argv);
#ifdef USE_QUDA
  initQuda(-1000);
#endif
}
void End(){
  // printf("End()\n");
#ifdef USE_QUDA
  endQuda();
#endif
// generates a lot of error messages. Turning off for now.
//  QMPSCU::destroy_qmp();
  // printf("destroy_qmp()\n");
//  _mcleanup();
  // printf("End()\n");
}
#else
void Start(int * argc, char *** argv) {
  GJP.setArg(argc,argv);
#ifdef USE_QUDA
  initQuda(QudaParam.device);
#endif
}
void End(){
#ifdef USE_QUDA
  endQuda();
#endif
}
#endif
CPS_END_NAMESPACE
