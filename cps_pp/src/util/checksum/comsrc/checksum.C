#include <config.h>

CPS_START_NAMESPACE

//------------------------------------------------------------------
//checksum.C
//
//Checksum class is checksuming evolution-related quantities, such as
//gauge configuration, momentum, CG source/solution vector, 
//CG intermediate vectors, global sums, etc.
//
//Meifeng Lin, April, 2005
//-------------------------------------------------------------------


CPS_END_NAMESPACE
#include <util/checksum.h>
//#include <qcdocos/gint.h>
#include <stdio.h>
#include <string.h> //strcpy
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/gjp.h>
#include <util/time_cps.h>

CPS_START_NAMESPACE

CheckSum CSM;

//intialize global and static variables
const int MAX_CSUM_LENGTH = 200000;
const char *csum_id[] = {"CSUM_EVL_LAT","CSUM_EVL_MOM","CSUM_GLB_LOC","CSUM_GLB_SUM","CSUM_EVL_SRC","CSUM_EVL_SOL","CSUM_MMP_SUM","CSUM_EVL_MMP" };
const char *comment = "evolution step";
int CheckSum::traj_cnt = 0;
int CheckSum::counter = 0;
unsigned long CheckSum::sum_sum[CSUM_ALL];

#define PROFILE

CheckSum::CheckSum(){
  cname = "CheckSum";
  char *fname = "CheckSum()";
//  VRB.Func(cname,fname);
  is_initialized = 0;
  csum = NULL;
 }


CheckSum::~CheckSum(){
  char *fname = "~CheckSum()";
//  printf("%s::%s Entered\n",cname,fname);
  VRB.Func(cname,fname);
  if(is_initialized){
    sfree(cname,fname,"csum",csum);
    is_initialized = 0;
  }
  VRB.FuncEnd(cname,fname);
}

void CheckSum::Initialize(){
  char *fname = "Initialize()";
  VRB.Func(cname,fname);
  Initialize(MAX_CSUM_LENGTH);
}

void CheckSum::Initialize(int len){
  char *fname = "Initialize(int)";
  VRB.Func(cname,fname);

  if( is_initialized ) return;
  filename = "csum.log";
  length = len;
  
  csum = (unsigned long *)smalloc(cname,fname,"csum",length*sizeof(unsigned long));

  for(int i = 0; i < CSUM_ALL; i++){
    csum_switch[i] = 0;
    sum_sum[i] = 0x0;
  }
  is_initialized = 1;
}


void CheckSum::Activate(CsumType csum_type){
  char *fname = "Activate(CsumType)";
  VRB.Func(cname,fname);
  VRB.Flow(cname,fname,"csum_type=%d",csum_type);

  if ( csum_type >= 0 && csum_type < CSUM_ALL) csum_switch[csum_type] = 1;
  else if (csum_type == CSUM_ALL) 
    for( int i = 0; i < CSUM_ALL; i++) csum_switch[i] = 1;
}

void CheckSum::Deactivate(CsumType csum_type){
  char *fname = "Deactivate(CsumType)";
  VRB.Func(cname,fname);

  if ( csum_type >= 0 && csum_type < CSUM_ALL) csum_switch[csum_type] = 0;
  else if (csum_type == CSUM_ALL) 
    for( int i = 0; i < CSUM_ALL; i++) csum_switch[i] = 0;
}


void CheckSum::Activate(int level){
  char *fname = "Activate(int)";
  VRB.Func(cname,fname);
  VRB.Flow(cname,fname,"level=%d",level);
  int activate = 0;
  if ( level >= 0 && level <= CSUM_ALL) activate = level;
  else if ( level > CSUM_ALL ) activate = CSUM_ALL;
  for ( int i = 0; i < activate; i++) csum_switch[i] = 1;
}

void CheckSum::SaveCsum(CsumType csum_type,unsigned long sum){
  char *fname = "SaveCsum(CsumType, unsigned long)";
  VRB.Func(cname,fname);
  VRB.Flow(cname,fname,"csum_type=%d sum=%p",csum_type,sum);

  if ( csum_type >= CSUM_ALL || csum_type < 0 ) 
    ERR.General(cname,fname,"Unknown checksum type\n");

  if ( csum_switch[csum_type] )
    SaveCsum(csum_id[csum_type],sum);
  VRB.FuncEnd(cname,fname);
}


void CheckSum::SaveCsum(const char *csum_char,unsigned long sum){
  char *fname = "SaveCsum(const char*, unsigned long)";
  VRB.Func(cname,fname);
    unsigned long *ptr = csum + counter;
  
  if ( csum != NULL ){
#if 0
    *ptr = (unsigned long)csum_char;
    *(ptr + 1) = sum;
    counter+=2;
    if( counter >= length ) counter = 0;
#else
//    printf("Current counter = %d; pointer = %p\n",counter,ptr);
#endif
  }
  VRB.FuncEnd(cname,fname);
}

void CheckSum::SaveComment(unsigned long num){
  SaveCsum(comment,num);
}

void CheckSum::AccumulateCsum(CsumType csum_type,unsigned long sum){
  sum_sum[csum_type] ^= sum;

}

void CheckSum::Clear(CsumType csum_type){
  sum_sum[csum_type] = 0x0;
}

void CheckSum::SaveCsumSum(CsumType csum_type){
  SaveCsum(csum_type,sum_sum[csum_type]);
}

void CheckSum::SetFilename(char *new_name){
  char *fname = "SetFilename(char *)";
  VRB.Func(cname,fname);
  strcpy(filename,new_name);
}

void CheckSum::Print(){
  char *fname = "print()";
  VRB.Func(cname,fname);

  char logfile[256];
  sprintf(&logfile[0],"%d.%s",UniqueID(),filename);
#if 1
//  printf("Current counter = %d; traj = %d\n",counter,traj_cnt);
#else
  FILE *fp = fopen(logfile,"a");
  Print(fp);
  fclose(fp);
#endif
}

void CheckSum::Print(FILE *fp){
  char *fname = "print(FILE *)";
  VRB.Func(cname,fname);
  VRB.Result(cname,fname,"Start unloading checksums\n");
#ifdef PROFILE
  Float time = -dclock();
#endif 
  if ( fp == NULL ) ERR.General(cname,fname,"File pointer not initialized\n");
  unsigned long step_cnt=0;
  for( int i = 0; i < counter; i+=2){
    if(strcmp(comment,(char *)csum[i]) == 0) { //check if it is a comment string
      step_cnt = csum[i+1]; 
      if ( step_cnt == 0 ) traj_cnt++; //if start of the trajectory, increase the trajectory counter 
    }
    else{
      fprintf(fp,"%d\t%d\t%s\t%#x\n",traj_cnt,step_cnt,(char *)csum[i],csum[i+1]);
    }
  }

  fflush(fp);
  counter = 0;
  //  traj_cnt++;
  // Gint::SynchMachine();
#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,0,time);
#endif

}

CPS_END_NAMESPACE
