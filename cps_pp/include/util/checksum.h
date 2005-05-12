#include <config.h>
#include <stdio.h>

CPS_START_NAMESPACE

//------------------------------------------------------------
//
// checksum.h
//
// Hearder file for the Checksum class
//
// It checksums various quantities and put them into memory.
// They are printed out to log files for each node upon user request.
// The checksums are organized on a first-come-first-served basis.
// The fields of the output are following:
// <Traj. No.> <step No.> <checksum identifier> <32bit checksum>
// The pointers to checksum identifiers are stored in front of each checksum. 
// CSUM_EVL_LAT -- local gauge field checksum
// CSUM_EVL_MOM -- local momentum checksum
// CSUM_EVL_SRC -- local CG source vector checksum
// CSUM_EVL_SOL -- local CG solution vector checksum
// CSUM_EVL_MMP -- local intermediate solution vector checksum
// CSUM_MMP_SUM -- checksum of all MMPs at each step
// CSUM_GLB_SUM -- running checksum of all global sums at each CG step
// CSUM_GLB_LOC -- running checksum of all local values to be globalsummed 
//                 at each CG step
//
// Meifeng Lin, April 2005
//-------------------------------------------------------------------

#ifndef INCLUDED_CHECKSUM_H
#define INCLUDED_CHECKSUM_H

enum CsumType{
  CSUM_EVL_LAT,
  CSUM_EVL_MOM,
  CSUM_GLB_LOC,
  CSUM_GLB_SUM,
  CSUM_EVL_SRC,
  CSUM_EVL_SOL,
  CSUM_MMP_SUM,
  CSUM_EVL_MMP,
  CSUM_ALL
};


class CheckSum{
 private:
  char *cname;
  unsigned long *csum;
  char *filename;
  int csum_switch[CSUM_ALL];
  int length;
  int is_initialized;
  static int traj_cnt;
  static int counter;
  static unsigned long sum_sum[CSUM_ALL];
  void SaveCsum(const char *csum_char, unsigned long sum);
 
 public:
  CheckSum();
  ~CheckSum();

  //Initialize array using default length
  void Initialize();

  //Initialize array using user-defined length
  void Initialize(int buf_len);

  //Activate a specific checksum
  void Activate(CsumType csum_type);

  //Activate a number of checksums
  void Activate(int level);

  //Deactivate a specific checksum
  void Deactivate(CsumType csum_type);

  //Overwrite the default filename. This is to be used with Print().
  void SetFilename(char *new_name);

  //Default print method: checksums are *appended* to (UniqueID()).csum.log.
  //in host
  void Print();

  //Print to file or terminal(stdout,stderr,etc.) as defined
  void Print(FILE *fp);

  //Save checksum into the array and increase the counter
  void SaveCsum(CsumType csum_type,unsigned long sum);

  //Save the step counter into the array
  void SaveComment(unsigned long num);

  //This is for buffering the running checksums.
  //Each CsumType has a separate temp variable sum_sum.
  //This function XOR all the numbers entered to sum_sum.
  void AccumulateCsum(CsumType csum_type,unsigned long sum);

  //This function zeros sum_sum for the particular checksum type.
  void Clear(CsumType csum_type);

  //This puts the accumulated checksum into the array and increases 
  //the counter as appropriate.
  void SaveCsumSum(CsumType);
};


extern CheckSum CSM;

#endif
CPS_END_NAMESPACE

