//VML file for staggered spectroscopy

enum StagQuarkSrcType{ 
	S_QUARK_POINT = 0, 
	WALLZ = 1, 
	WALL2Z = 2
};	/* StagQuarkSrc type*/
enum StagQuarkLocalType { 
	LOCAL = 0, 
	NONLOCAL =1 
};		/* StagQuarkArg sln*/

enum StagQuarkDir { 
	HDM_X = 0, 
	HDM_Y = 1, 
	HDM_Z = 2, 
	HDM_T = 3 
};/* dir*/

//enum { MAX_LENGTH = 40 };

class StagQuarkSrc {
  StagQuarkSrcType type;
  int origin[4];
  int end[4];
  StagQuarkDir dir;
};

class StagQuarkArg {
  int qid;
  CgArg cg;
  StagQuarkSrc src;
  StagQuarkLocalType sln;
};
  
class StagMesonArg {
  int qid0;
  int qid1;
  StagQuarkDir dir;
  unsigned long meson_buf;
};

class StagMomMesonArg {
  int qid0;
  int qid1;
  StagQuarkDir dir;
  int no_of_momenta;
  unsigned long meson_buf;
};

class StagNucleonArg {
  int qid0;
  int qid1;
  int qid2;
  StagQuarkDir dir;
  unsigned long nucleon_buf;
};

class StagNonLocalArg {
  int qid0;
  int qid1;
  int qid2;
  StagQuarkDir dir;
  unsigned long nlocal_buf;
};

class NLStagMesonArg {
  int qid0[8];
  StagQuarkDir dir;
  unsigned long nlocal_buf;
};


/*struct StagQuarkSrc {
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
};*/

//#endif

//CPS_END_NAMESPACE
