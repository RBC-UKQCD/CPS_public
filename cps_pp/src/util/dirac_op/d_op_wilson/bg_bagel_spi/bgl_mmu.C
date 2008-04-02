#include <stdio.h>
#include "wfm.h"
#include "wfm_internal.h"

#if 1

void wfm::mmu_print(void){}
void wfm::mmu_optimise(void){}

#else 

#include "wfm_ppclib.h"

#define pid  0x30
#define mmucr  0x3b2
#define ccr0 0x3b3


void wfm::mmu_print(void)
{    
  if ( isBoss() ) { 
    printf("MSR  : %8.8x\n",mfmsr());
    printf("PID  : %8.8x\n",mfspr(pid));
    printf("CCR0 : %8.8x\n",mfspr(ccr0));
    printf("MMUCR: %8.8x\n",mfspr(mmucr));

    for(int t=0;t<64;t++){
	printf("tlb[%d] :",t);
	printf(" %8.8x ",tlbre(t,0));
	printf(" %8.8x ",tlbre(t,1));
	printf(" %8.8x ",tlbre(t,2));
	printf("\n");
    }
  }
}
void wfm::mmu_optimise(void)
{
    unsigned long word;
    // set:
    // swoa in mmucr
    mtspr(mmucr,mfspr(mmucr)|0x01000000);
    // gdcbt in ccr0
    word=mfspr(ccr0)|0x6000;// guarantee touches
    word=word&(~0x00200000); 
    mtspr(ccr0,word);
    // map a n/c tlb covering entire dram
    // switch to user mode?
}
#endif
