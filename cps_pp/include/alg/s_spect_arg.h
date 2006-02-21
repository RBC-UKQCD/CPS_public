#include<config.h>
CPS_START_NAMESPACE
//  s_spect_arg.h

#ifndef INCLUDED_S_SPECT_ARG_H
#define INCLUDED_S_SPECT_ARG_H
CPS_END_NAMESPACE
#include <util/vector.h>	// Float
#include <alg/cg_arg.h>
CPS_START_NAMESPACE

enum { POINT = 0, WALLZ, WALL2Z };	// StagQuarkSrc type
enum { LOCAL, NONLOCAL };		// StagQuarkArg sln
enum { HDM_X = 0, HDM_Y, HDM_Z, HDM_T };// dir

enum { MAX_LENGTH = 40 };

struct StagQuarkSrc {
   int type;
   int origin[4];
   int end[4];
   int dir;
};

struct StagQuarkArg {
  int qid;
  CgArg cg;
  StagQuarkSrc src;
  int sln;
};

struct StagMesonArg {
  int qid0;
  int qid1;
  int dir;
  Float *meson_buf;
};

struct StagMomMesonArg {
  int qid0;
  int qid1;
  int dir;
  int no_of_momenta;
  Float *meson_buf;
};

struct StagNucleonArg {
  int qid0;
  int qid1;
  int qid2;
  int dir;
  Float *nucleon_buf;
};

struct StagNonLocalArg {
  int qid0;
  int qid1;
  int qid2;
  int dir;
  Float *nlocal_buf;
};

struct NLStagMesonArg {
  int qid0[8];
  int dir;
  Float *nlocal_buf;
};


#endif

CPS_END_NAMESPACE
