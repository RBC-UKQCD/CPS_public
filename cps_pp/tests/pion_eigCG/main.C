//my class
#include"AlgPion.h"
//#include<alg/qpropw_arg.h>
//#include<alg/qpropw.h>


//argument class
#include<alg/eigcg_arg.h>
#include<alg/common_arg.h> //used to save the information of a class
#include<alg/do_arg.h>     //globle parameters
#include<alg/meas_arg.h>   //measurement regulation: how offen,lattic..
#include<alg/cg_arg.h>     //Conjugate Gradient Method (CG,BICGSTAB)
#include<util/qioarg.h>    //argument class for lattice read and write
//algrithm class
#include<alg/alg_meas.h>   //do the measurement??

//general purpose
#include<config.h>		   //???
#include<util/lattice.h>   //???so we included all the lattice classes?? like GwilsonFdwf
#include<util/lattice/lattice_types.h>
#include<util/gjp.h>       //Global job parameter. What does it do?
#include<util/verbose.h>   //control the out put verbosity
#include<util/error.h>     //dealing with error output
#include<util/qcdio.h>     //what for??
#include<util/qioarg.h>    //qioarg argument file
#include<util/ReadLatticePar.h>
						   //read the lattice in. parallel

//c++ classes
#include<ctime>
#include<iostream>
#include<sstream>
#include<unistd.h> //use the chdir(char *) function
#include<sys/stat.h>

//namespace include
using namespace cps;

bool testLatticeFit(DoArg &do_arg);
void ReadGaugeField(const MeasArg &meas_arg);
void ReadRNG(const MeasArg &meas_arg);

int main(int argc,char *argv[])
{
//CPviolation vmlDirectory starttraj endtraj

	if(argc!=2){cout<<"CPviolation vmlDirectory "<<endl;exit(0);}

	char *cname=argv[0];
	char *fname="main(int,char**)";

//-----------------------------------------
//--Argument files initialize--------------
	chdir(argv[1]);
	CommonArg common_arg("","");
	DoArg do_arg;
	MeasArg meas_arg;
	QPropWArg lqpropw_arg;
	EigCGArg eigcg_arg;
	
	if(!do_arg.Decode("do_arg.vml","do_arg")){std::cout<<"Can't open do_arg.vml!"<<std::endl;exit(1);}
	if(!meas_arg.Decode("meas_arg.vml","meas_arg")){std::cout<<"Can't open meas_arg!"<<std::endl;exit(1);}
	if(!testLatticeFit(do_arg)){std::cerr<<"Lattice does not fit on this computer!"<<std::endl;exit(1);}	
	if(!lqpropw_arg.Decode("lqpropw_arg.vml","lqpropw_arg")){std::cout<<"Can't open lqpropw_arg.vm!"<<std::endl;exit(1);}
	if(!eigcg_arg.Decode("eigcg_arg.vml","eigcg_arg")){std::cout<<"Can't read eigcg_arg correctly!"<<std::endl;exit(1);}

//----Initialization end----------------------

	//begin to work
	chdir(meas_arg.WorkDirectory);
	char newdir[100];
	sprintf(newdir,"pion_mass_%1.2f",lqpropw_arg.cg.mass);
	mkdir(newdir,0775);

	Start();
	do_arg.start_seed_value=(unsigned)time(NULL);
	GJP.Initialize(do_arg);
	LRG.Initialize();	

	{
		GnoneFnone lat; //get a lattice instantiation
	}
	/***********************************************
	main loop for each configuration
	***********************************************/
	std::cout<<"main loop begin"<<std::endl;
	for(int conf=meas_arg.TrajStart; conf<=meas_arg.TrajLessThanLimit;conf+=meas_arg.TrajIncrement)
	{
		std::cout<<"configuration :"<<conf<<std::endl;
		meas_arg.TrajCur = conf;
		ReadGaugeField(meas_arg);
		//ReadRNG(meas_arg);
		{
			char comm_arg_filename[200];
			sprintf(comm_arg_filename,"./pion_mass_%1.2f/traj_%d",lqpropw_arg.cg.mass,conf);
			common_arg.set_filename(comm_arg_filename);

			GimprRectFdwf lat;
			AlgPion pion(lat,&common_arg,&lqpropw_arg, &eigcg_arg);
			pion.runpion();
		}
	}
	std::cout<<"The End!"<<endl;
	/*End of the main configuration loop*/
	End();
	return 0;
}

bool testLatticeFit(DoArg & do_arg)
{
  do_arg.x_nodes = SizeX();
  do_arg.y_nodes = SizeY();
  do_arg.z_nodes = SizeZ();
  do_arg.t_nodes = SizeT();
  do_arg.s_nodes = SizeS();

  do_arg.x_node_sites = do_arg.x_sites/do_arg.x_nodes;
  do_arg.y_node_sites = do_arg.y_sites/do_arg.y_nodes;
  do_arg.z_node_sites = do_arg.z_sites/do_arg.z_nodes;
  do_arg.t_node_sites = do_arg.t_sites/do_arg.t_nodes;
  do_arg.s_node_sites = do_arg.s_sites/do_arg.s_nodes;

  if (do_arg.x_sites!=do_arg.x_node_sites*do_arg.x_nodes) return false; 
  if (do_arg.y_sites!=do_arg.y_node_sites*do_arg.y_nodes) return false;
  if (do_arg.z_sites!=do_arg.z_node_sites*do_arg.z_nodes) return false;
  if (do_arg.t_sites!=do_arg.t_node_sites*do_arg.t_nodes) return false;
  if (do_arg.s_sites!=do_arg.s_node_sites*do_arg.s_nodes) return false;

	return true;
}

void ReadGaugeField(const MeasArg &meas_arg)
{
	char *cname = "main";
	char *fname = "ReadGaugeField";

	GnoneFnone lat;
	std::stringstream lat_file;
	lat_file<<meas_arg.GaugeStem<<'.'<<meas_arg.TrajCur;
	QioArg rd_arg(lat_file.str().c_str(),0.001);
		//why do we have to check precision? what is for?
	rd_arg.ConcurIONumber=meas_arg.IOconcurrency;
	
	ReadLatticeParallel rl;
	rl.read(lat,rd_arg);
	if(!rl.good())ERR.General(cname,fname,"Failed read lattice %s",lat_file.str().c_str());	
}

void ReadRNG(const MeasArg &meas_arg)
{
	char *cname = "main";
	char *fname = "ReadRNG(MeasArg&)";
	std::stringstream rng_file;
	rng_file<<meas_arg.RNGStem<<'.'<<meas_arg.TrajCur;
	if ( !LRG.Read(rng_file.str().c_str()) ) 
		ERR.General(cname,fname,"Failed RNG file %s",rng_file.str().c_str());
}
