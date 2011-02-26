#ifndef ALGTWOPION_H
#define ALGTWOPION_H
#include<config.h>
#include<util/vector.h> //A class implementing a general 3 component complex vector
#include<util/lattice.h>
#include<util/rcomplex.h>
#include<alg/alg_base.h>
#include<util/momentum.h> //use ThreeMom class to control momentum
#include<alg/common_arg.h>
#include<alg/qpropw_arg.h>
#include<alg/qpropw.h>
#include<alg/eigcg_arg.h>
using namespace cps;
class AlgPion : public Alg
{
private:
	char* cname;
	QPropWArg *lqprop_arg;
	EigCGArg *eigcg_arg;
	char *filestub;
	void writeCorr(Rcomplex *corr, char *filename);
public:
	AlgPion(Lattice & latt, CommonArg* comm_arg, QPropWArg *lqprop_arg, EigCGArg *eigcg_arg);
	virtual ~AlgPion();

	void runpion();
};
#endif
