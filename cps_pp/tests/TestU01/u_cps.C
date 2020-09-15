#include "unif01.h"
#include "gdef.h"
extern "C"{
#include "util.h"
}
#include <stdlib.h>
#include <config.h>
#include <util/gjp.h>
#include <util/latheader.h>

USING_NAMESPACE_CPS
static int CPS_initted=0;
typedef struct{
  int size[4];
  unsigned long seed;
  int Ntake;
  int Ndrop;
} CPS_param;

typedef struct{
  int pos[4];
  int i_take;
} CPS_state;

const unsigned long TWO_32 = 4294967296;

extern "C" {

double CPS_GetU01 (void *p, void *s){
	CPS_param *param = (CPS_param *)p;
	CPS_state *state = (CPS_state *)s;

	int *pos = state->pos;
	LRG.AssignGenerator(pos);
	double temp = LRG.Urand(0.,1.);
	state->i_take++;
	if (state->i_take >= param->Ntake){
		for(int i=0;i<param->Ndrop;i++) LRG.Urand();
		pos[0] +=2;
		for(int dir=0;dir<3;dir++)
		if (pos[dir] >= param->size[dir]){
			pos[dir]=0; pos[dir+1] += 2;
		}
		if (pos[3] >= param->size[3]) pos[3] =0;
		state->i_take = 0;
	}
	return temp;

}

unsigned long CPS_GetBits (void *p, void *s){
	double temp = CPS_GetU01(p,s)*TWO_32;
	
	return (unsigned long) temp;
}

void CPS_Write (void *s){
	CPS_state *state = (CPS_state *)s;
//	printf("CPS_Write() Not implemented!\n");
	int *pos = state->pos;
	printf("pos: %d %d %d %d i_take: %d\n",pos[0],pos[1],pos[2],pos[3],state->i_take);
}
unif01_Gen * Create_CPS (int *argc, char ***argv, char *doarg_name,  int take, int drop){

	
	unif01_Gen *gen;
	CPS_param *param;
	CPS_state *state;

	if (CPS_initted) {printf("CPS already initted!\n");exit(-10);}
	Start(argc, argv);
	DoArg do_arg;
	if (! do_arg.Decode(doarg_name,"do_arg"))
	ERR.General("","Create_CPS","Bad do_arg.vml\n");
	GJP.Initialize(do_arg);
	LRG.Initialize();
	
	gen = (unif01_Gen *) util_Malloc(sizeof(unif01_Gen));
	param = (CPS_param*) util_Malloc(sizeof(CPS_param));
	state = (CPS_state *) util_Malloc(sizeof(CPS_state));
	
	for(int i=0;i<4;i++) param->size[i] = GJP.NodeSites(i);
	param->Ntake = take;
	param->Ndrop = drop;
	for(int i=0;i<4;i++) state->pos[i] = 0;
	state->i_take=0;

	printf ("Ndrop=%d Ntake=%d\n", param->Ndrop, param->Ntake);

	CPS_initted=1;
	gen->param = param;
	gen->state = state;
//	char *rng_name =new char[(LatRngHeader::rng_datatype).size()+1];
 //	strncpy(rng_name,(LatRngHeader::rng_datatype).c_str(),(LatRngHeader::rng_datatype).size());
	gen->name  = (char *)cps::LatRngHeader::RNGString;
	return gen;
}

unif01_Gen * unif_CreateCPS (int *argc, char ***argv, char *doarg_name,  int take, int drop){

	unif01_Gen *gen;
	gen = Create_CPS(argc,argv,doarg_name,take,drop);
	gen->GetBits = &CPS_GetBits;
	gen->GetU01 = &CPS_GetU01;
	gen->Write = &CPS_Write;
	return gen;
}


void unif_DeleteCPS( unif01_Gen *gen){
	if (!gen) return;
	CPS_param *param;
	CPS_state *state;
	state = (CPS_state *) gen->state;
	util_Free(state);
	util_Free(gen->param);
	util_Free(gen);
}



}
