#include<config.h>
#ifdef USE_QMP
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of local_checksum and global_checksum outine.

*/
//--------------------------------------------------------------------
//  CVS keywords
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qmp/scu/check_sum.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
//  local_checksum
//
//  Calculating checksum on a set of data on one node
//
//  global_checksum
//
//  Calculating checksum on each node and compare against each other
//--------------------------------------------------------------


CPS_END_NAMESPACE
#include <util/data_types.h>
#include<comms/glb.h>
#include <comms/glb_sum_internal.h>
//#include<comms/scu.h>
//#include<util/gjp.h>
//#include<comms/double64.h>
//#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE

unsigned int local_checksum(Float * float_p, int len) {
  unsigned int csum = 0;
  unsigned int * uint_p = (unsigned int *) float_p;
  for(int i=0;i < len * sizeof(Float) / sizeof(unsigned int); i++) {
//    csum += *uint_p;
    csum ^= *uint_p;
    uint_p ++;
  }
  return csum;
}

unsigned int global_checksum(Float * float_p, int len) {
  static int initted = 0;
  unsigned int locsum;

#ifdef UNIFORM_SEED_TESTING
  locsum = test_checksum(float_p,len);
#else
  locsum = local_checksum(float_p, len);
#if 1
  glb_sum_internal2(&locsum,4,0); // 0 for XOR
#else
  if (!initted)
  if (sizeof(unsigned int) != sizeof(unsigned long))
    ERR.General("","global_checksum",
    "sizeof(unsigned int)(%d) != sizeof(unsigned long)(%d)\n",
    sizeof(unsigned int),sizeof(unsigned long));
  QMP_xor_ulong ((unsigned long*)&locsum);
#endif
#endif
  initted=1;
  return locsum;
}

unsigned int test_checksum(Float * float_p, int len) {
  static int gcsum_count = 1;
  static int peRank = 0;
  if (gcsum_count==1){
    peRank = QMP_get_node_number();
//    fprintf(stderr,"peRank=%d\n",peRank);
  }
  unsigned int locsum = local_checksum(float_p, len);
  unsigned int bosssum;

//  if(UniqueID() == 0) bosssum = locsum;
  if(CoorX()==0 &&
  CoorY()==0 &&
  CoorZ()==0 &&
  CoorT()==0 ) bosssum = locsum;
  else                bosssum = 0;
  glb_sum_internal2(&bosssum,4,0); // 0 for XOR
  if(bosssum != locsum) {
      fprintf( stderr,
               "GCheckSum %d : Node %d : Oops I did it again: me (%d,%d,%d,%d,%d) %u != boss %u\n",
               gcsum_count,
	       UniqueID(),
               CoorX(),
               CoorY(),
               CoorZ(),
               CoorT(),
               CoorS(),
               locsum,  bosssum );
      exit(-30);
  }
#if 0
  if(gcsum_count>2 && peRank==0){
    fprintf(stderr,"peRank=%d exiting\n",peRank);
    exit(-24);
  }
#endif
  gcsum_count++;

  return bosssum;
}

CPS_END_NAMESPACE
#endif
