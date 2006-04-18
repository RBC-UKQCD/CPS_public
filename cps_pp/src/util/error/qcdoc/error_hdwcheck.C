#include<config.h>
#include<util/qcdio.h>
#include<util/error.h>
CPS_START_NAMESPACE
void Error::HdwCheck(const char *cname, const char *fname){
#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
  if(!ScuChecksum::ChecksumsOn())
  ScuChecksum::Initialise(true,true);
  if ( ! ScuChecksum::CsumSwap() ){
    fprintf(stderr,"%s::%s: SCU Checksum mismatch\n",cname,fname);
    exit(exit_value[hardware]);
  }
#else
  printf("This version of QOS does not have SCU checksum class\n");
#endif
}

CPS_END_NAMESPACE
