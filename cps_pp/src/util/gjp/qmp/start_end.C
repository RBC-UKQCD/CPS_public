#include <config.h>
#include <util/gjp.h>
#include <qmp.h>

//QMP definitions for CPS::Start() and CPS::End()

extern "C" void _mcleanup();
CPS_START_NAMESPACE
void Start(){
  QMPSCU::init_qmp();
}

void Start(int * argc, char *** argv) {
//  printf("Start(%d %p)\n",*argc,*argv);
  //Initialize QMP
  QMPSCU::init_qmp(argc, argv);
}
void End(){
  QMPSCU::destroy_qmp();
//  _mcleanup();
}
CPS_END_NAMESPACE
