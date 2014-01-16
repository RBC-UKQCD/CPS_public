#ifndef HDCGCONTROLLER_H
#define HDCGCONTROLLER_H

#include<vector>
#include<complex>
#include<stdio.h>
#include<stdlib.h>

#include<bfm_qdp.h>
#include<BfmMultiGrid.h>
#include<util/time_cps.h>
#include<comms/sysfunc_cps.h>
//#include<comms/sysfunc_qmp.h>

#include <bfm_hdcg_wrapper.h>
class HDCGInstance{
	public:
	static BfmMultiGridParams  Params;
	static HDCG_wrapper  * _instance ;
	static HDCG_wrapper *getInstance(){return _instance;} 
	static HDCG_wrapper *setInstance(HDCG_wrapper *_new){_instance = _new;} 
	static void free(){
		if (_instance){ 
			_instance->HDCG_end();
			delete _instance;
		}
		_instance=NULL;
	}
	static void SetDefault(){
		if (!CPS_NAMESPACE::UniqueID()) printf("HDCGInstance()\n");
// Partially tuned for 64^3, from PAB's message
#if 0
		Params.NumberSubspace=40;
		Params.SubspaceRationalLs=8;
		Params.SubspaceRationalLo=0.001;
		Params.SubspaceRationalResidual=1e-5;
		Params.SubspaceRationalRefine=1;
		Params.SubspaceRationalRefineResidual=1e-5;
		Params.SubspaceRationalRefineLo=0.001;
		Params.SubspaceSurfaceDepth=12;
		Params.PreconditionerKrylovResidual=1e-5;
		Params.PreconditionerKrylovIterMax=8;
		Params.PreconditionerKrylovShift=1.0;
		Params.LittleDopSolverResidualInner=5e-3;
		Params.LittleDopSolverResidualVstart=1e-4;
		Params.LittleDopSolverResidualSubspace=1.0e-7;
		Params.LittleDopSolverIterMax = 10000;
		Params.LittleDopSolver = LittleDopSolverCG;
#else  
// 48^3
		Params.NumberSubspace=44;
		Params.SubspaceRationalMass=0.00078;
		Params.SubspaceRationalLs=12; // Ls/2
		Params.SubspaceRationalLo=0.001;
		Params.SubspaceRationalResidual=1e-5;
		Params.SubspaceRationalRefine=1;
		Params.SubspaceRationalRefineResidual=1e-4;
		Params.SubspaceRationalRefineLo=0.003;
		Params.SubspaceSurfaceDepth=24;
		Params.PreconditionerKrylovResidual=1e-5;
		Params.PreconditionerKrylovIterMax=6;
		Params.PreconditionerKrylovShift=1.0;
		Params.LittleDopSolverResidualInner=5e-3;
		Params.LittleDopSolverResidualVstart=1e-5;
		Params.LittleDopSolverResidualSubspace=1.0e-7;
		Params.LittleDopSolverIterMax = 10000;
//		Params.LittleDopSolver = LittleDopSolverCG;
#endif
		if (!CPS_NAMESPACE::UniqueID()) printf("HDCGInstance::Params defaults set\n");
	}
//	~HDCGInstance(){}
};

#if 0
template <class Float>
class HDCGController
{
	private:
#if 0
		const int nev; //typical nev=8;
		const int m; //typical m=24; should be bigger than 2*nev
		const int max_def_len;//the size of the deflation space that we try to accumulate
		const double max_eig_cut; //throw away those fake low modes
		const bool always_restart; //it is better not to usually
#endif

		const int Ls;
		const int NumberSubspace;
		int block[5];
		int quadrant[4];

//		Float* V; //the space for eigen vectors
//		std::vector<Float*> U;	  //at some point, we want to use float for U no matter what.
		//To save the usage of memory, 
		//it is fine to use float since the low modes are not 
		//accurately any way and it has not effect on the final result
		

//		const int vec_len; //for convenient use between cps and bfm, this is as an input parameter. 
		//better be infered from bfm_qdp class, it is the length of Fermion_t
	public:
//		int def_len; //keep the record of the number of low modes we achieved
//		std::complex<double> *H;
//		std::vector<double> restart;
		

	private:
		static HDCGController* _instance;
		static BfmMultiGrid<Float>* _bfm_mg;
		static char *cname;
		HDCGController(HDCGController<Float> &){};
		HDCGController<Float> & operator=(HDCGController const &){};
		HDCGController(int _Ls, int _NS,int _block[],int _quad[]):
			Ls(_Ls),NumberSubspace(_NS)
			{
				for(int i=0;i<5;i++) block[i]=_block[i];
				for(int i=0;i<4;i++) quadrant[i]=_quad[i];
			}
#if 0
		HDCGController(int _nev, int _m, int _max_def_len, double _max_eig_cut, std::vector<double> &_restart, bool _ar, int cg_vec_len):
			nev(_nev),m(_m),max_def_len(_max_def_len),max_eig_cut(_max_eig_cut),restart(_restart),vec_len(cg_vec_len),always_restart(_ar)
		{
			if(m<=2*nev)
			{
				printf("m should be larger than 2*nev !\n");
				exit(-1);
			}
			def_len = 0;

			V = new Float[vec_len*m];
			if(V==NULL){
				printf("fail to malloc space for pointer V in eigcgcontroller\n");
				exit(-1);
			}

			U.resize(max_def_len,NULL);
			for(int i=0;i<max_def_len;i++)
			{
				U[i] = new Float[vec_len];
				if(U[i]==NULL){
					printf("fail to malloc space for pointer U[%d] \n",i);
					exit(-1);
				}
			}
			H=new std::complex<double>[max_def_len*max_def_len];
			//The max length. In the first few uses, only def_len*def_len is used.
			if(H==NULL)
			{
					printf("fail to malloc space for pointer H\n");
					exit(-1);
			}
			if(restart.size()>2)
			{
				//check if it is decreasing order!
				for(unsigned int i=1;i<restart.size();i++)
				{
					if(restart[i]>=restart[i-1])
					{
						printf("restart should in decreasing order!\n");
						exit(-1);
					}
				}
			}
		}
#endif

	public:
		static void free()
		{
			freeHDCG();
			delete _instance;
			_instance = NULL;
		}
		static HDCGController<Float>* getInstance()
		{
//			if(NULL == _instance)
			if(0)
			{
				printf("Need to setup/initialize the HDCGController instance first!\n");
				exit(-1);
				//_instance = new HDCGController<Float>();
			}
			return _instance;
		}
		static HDCGController<Float>* setInstance(int _Ls, int _NS,int _block[],int _quad[])
		{
			if(NULL == _instance)
			{
				_instance = new HDCGController<Float>(_Ls, _NS,_block,_quad);
			}
			else
			{
				printf("HDCGController should be initialized once only! or should be deleted before another use for a different case!\n");
				exit(-1);
			}
			return _instance;
		}
		~HDCGController()
		{
			printf(" ~HDCGController() called!!" );
//			if(H!=NULL) delete [] H; H=NULL;
//			if(V!=NULL) delete [] V; V=NULL;
//			for(int i=0;i<max_def_len;i++){if(U[i]!=NULL)delete [] U[i];U[i]=NULL;}
		}
#if 0
		int get_nev() const{return nev;}
		int get_m() const{return m;}
		int get_max_def_len() const{return max_def_len;}
		int get_vec_len() const{return vec_len;}
		double get_max_eig_cut() const{return max_eig_cut;}
		Float* getU(int i){return U.at(i);}
		Float* getV(int i){return V+i*vec_len;}
		bool is_always_restart() const{return always_restart;}
#endif

#if 1
		static BfmMultiGrid<Float>* getHDCG()
		{
//			if(NULL == _bfm_mg)
//			{
//				printf("Need to setup/initialize BfmMultiGrid first!\n");
//				exit(-1);
//				_instance = new HDCGController<Float>();
//			}
			return _bfm_mg;
		}
#endif

#if  1
		void   setHDCG(bfm_qdp<double> &dop, bfm_qdp<float> &dop_sp)
		{
			if(NULL == _bfm_mg)
			{
				Float t1 = -CPS_NAMESPACE::dclock();
				_bfm_mg = new BfmMultiGrid<double>(Ls,NumberSubspace,block,quadrant,&dop,&dop_sp);
//Initialization sequence
				_bfm_mg->RelaxSubspace(&dop);
				_bfm_mg->ComputeLittleMatrixColored();
				_bfm_mg->LdopDeflationBasisInit(8);
				_bfm_mg->LdopDeflationBasisInit(16);
				_bfm_mg->LdopDeflationBasisInit(32);
				_bfm_mg->LdopDeflationBasisInit(64);
				_bfm_mg->LdopDeflationBasisInit(128);
				_bfm_mg->LdopDeflationBasisDiagonalise(128);
				_bfm_mg->SinglePrecSubspace();
				
				t1 +=CPS_NAMESPACE::dclock();
				CPS_NAMESPACE::print_flops(cname,"setHDCG()",0,t1);

			}
			else
			{
				printf("Should be initialized once only! or should be deleted before another use for a different case!\n");
				exit(-1);
			}
//			return _bfm_mg;
		}
		void   freeHDCG(){
			delete _bfm_mg;
			_bfm_mg = NULL;
		}
#endif
};

template<class Float> HDCGController<Float>* HDCGController<Float>::_instance =  NULL;
template<class Float> BfmMultiGrid<Float>* HDCGController<Float>::_bfm_mg =  NULL;
template<class Float> char* HDCGController<Float>::cname = "HDCGController";


#endif

#endif
