#include <config.h>
#include <util/gjp.h>
#include <util/error.h>
#ifdef USE_QMP
#include <qmp.h>
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
  QMPSCU::init_qmp(argc, argv);
  GJP.setArg(argc,argv);
#ifdef USE_QUDA
  initQuda(-1000);
#endif
}
void End(){
  QMPSCU::destroy_qmp();
//  _mcleanup();
#ifdef USE_QUDA
  endQuda();
#endif
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
