/* */

enum VP_TYPE { 
     CONSERVED_LOCAL			= 0,
     CONSERVED_LOCAL_LOWMODE		= 1,
     CONSERVED_LOCAL_PTSRC_TEST         = 2,		
     CONSERVED_LOCAL_LOOP_PTSRC		= 3,
     CONSERVED_LOCAL_TWISTED		= 4
};

class MuonArg {

u_long u1lat;
int source_time;
int oper_time_start;
int oper_time_end;
int operator_gamma;
int x[3];
int n_source;
int source_inc;
Float loop_mass;
Float line_mass;
Float charge;
int GaugeFix;
 int NPROJ;
 Float EPS; 
 int XMOM;
 int YMOM;
 int ZMOM;
 int Nmom;
 int MaxMomSq;
 int NConfs;
 int NHITS;
 string DIRVP<>; 
 string DIRML<>;
 string VPTAG<>;
 int DO_MUONLINE; 
 int DO_VACPOL; 
 int tINC;
 int ptStart;
 int ptINC;
 int conf;
 VP_TYPE vp_kind;
 string EIGTAG<>;
 memfun MuonArg();
};

