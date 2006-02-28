#include<config.h>
#include<util/qcdio.h>
#include<util/error.h>
#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
#include <qcdocos/scu_checksum.h>
#endif
CPS_START_NAMESPACE
void Error::HdwCheck(char *cname, char *fname){
#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
  if(!ScuChecksum::ChecksumsOn())
  ScuChecksum::Initialise(true,true);
  if ( ! ScuChecksum::CsumSwap() )
    ERR.Hardware(cname,fname, "SCU Checksum mismatch\n");
#else
  printf("This version of QOS does not have SCU checksum class\n");
#endif
}

CPS_END_NAMESPACE
