#ifndef EIGCGCONTROLLER_H
#define EIGCGCONTROLLER_H

//by Qi Liu 2012
#include<vector>
#include<complex>
#include<stdio.h>
#include<stdlib.h>

template <class Float>
class EigCGController
{
	private:
		const int nev; //typical nev=8;
		const int m; //typical m=24; should be bigger than 2*nev
		const int max_def_len;//the size of the deflation space that we try to accumulate
		const double max_eig_cut; //throw away those fake low modes
		const bool always_restart; //it is better not to usually

		Float* V; //the space for eigen vectors
		std::vector<Float*> U;	  //at some point, we want to use float for U no matter what.
		//To save the usage of memory, 
		//it is fine to use float since the low modes are not 
		//accurately any way and it has not effect on the final result
		

		const int vec_len; //for convenient use between cps and bfm, this is as an input parameter. 
		//better be infered from bfm_qdp class, it is the length of Fermion_t
	public:
		int def_len; //keep the record of the number of low modes we achieved
		std::complex<double> *H;
		std::vector<double> restart;

	private:
		static EigCGController* _instance;
		EigCGController(EigCGController<Float> &){};
		EigCGController<Float> & operator=(EigCGController const &){};
		EigCGController(int _nev, int _m, int _max_def_len, double _max_eig_cut, std::vector<double> &_restart, bool _ar, int cg_vec_len):
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

	public:
		static void free()
		{
			delete _instance;
			_instance = NULL;
		}
		static EigCGController<Float>* getInstance()
		{
			if(NULL == _instance)
			{
				printf("Need to setup/initialize the instance first!\n");
				exit(-1);
				//_instance = new EigCGController<Float>();
			}
			return _instance;
		}
		static EigCGController<Float>* setInstance(int nev, int m, int max_def_len, double max_eig_cut, std::vector<double> &restart, bool always_restart, int vec_len)
		{
			if(NULL == _instance)
			{
				_instance = new EigCGController<Float>(nev, m, max_def_len, max_eig_cut, restart, always_restart, vec_len);
			}
			else
			{
				printf("Should be initialized once only! or should be deleted before another use for a different case!\n");
				exit(-1);
			}
		}
		~EigCGController()
		{
			if(H!=NULL) delete [] H; H=NULL;
			if(V!=NULL) delete [] V; V=NULL;
			for(int i=0;i<max_def_len;i++){if(U[i]!=NULL)delete [] U[i];U[i]=NULL;}
		}

		int get_nev() const{return nev;}
		int get_m() const{return m;}
		int get_max_def_len() const{return max_def_len;}
		int get_vec_len() const{return vec_len;}
		double get_max_eig_cut() const{return max_eig_cut;}
		Float* getU(int i){return U.at(i);}
		Float* getV(int i){return V+i*vec_len;}
		bool is_always_restart() const{return always_restart;}
};

template<class Float> EigCGController<Float>* EigCGController<Float>::_instance =  NULL;

#endif
