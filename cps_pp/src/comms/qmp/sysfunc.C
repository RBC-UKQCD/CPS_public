#include<config.h>
#ifdef USE_QMP
/*----------------------------------------------------------*/
/*!\file
  \brief  Definitions for the QMP implementation of the QCDSP/QCDOC communications layer.
  
*/
/*----------------------------------------------------------------------
  The Sysfunc Comms Interface: sysfunc.C

  The QMP implementation of the QCDOC/QCDSP SCU comms-layer.

  M. Cheng michaelc@phys.columbia.edu
  -----------------------------------------------------------

----------------------------------------------------------*/

#include <comms/sysfunc_qmp.h>
//#include <comms/scu_dir_arg.h>
#include <util/qcdio.h>
#include <util/error.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <qmp.h>
#if TARGET == BGL
#include <sys/bgl/bgl_sys_all.h>
#endif
#if TARGET == BGP
void spi_init();
#endif
CPS_START_NAMESPACE

/*!\namespace cps
  \brief Main namespace for all CPS classes, functions <em>etc</em>.
*/
/*!\namespace cps::QMPSCU
  \brief Namespace for the QMP emulations of the SCU.
*/



//Useful global variables
namespace QMPSCU {
  
  static bool initialized = false;
  
  //! Number of grid dimensions.
  static int NDIM = 5;
 
  static int peRank;        /*!< Rank/identify of this  process */
  static int peNum;          /*!< Total number of processors */
#ifdef UNIFORM_SEED_NO_COMMS
  static const int pePos[4] = {0,0,0,0};
  static const int peGrid[4] = {1,1,1,1};
#else
  static const int* pePos;  /*!< Position of this process in the grid.*/ 
  static const int* peGrid; /*!< Number of processors in each direction */
#endif

  //Clean up resources used by QMP
  void destroy_qmp() {
    QMP_finalize_msg_passing();
  }

  //Initialize QMP with null command line
  void init_qmp() {
    ERR.General("","init_qmp()","default arguments no loger supported. Call init_qmp(&argc,&argv) via (CPS_NAMESPACE)::Start(&argc,&argv)");
  }

//Initialize QMP
//Get Allocated machine size, and declare logical machine

void init_qmp(int * argc, char ***argv) {

#if 0
  printf("init_qmp(%d %p)\n",*argc,*argv);
  for(int i = 0; i<*argc;i++){
    printf("argv[%d](before)=%s\n",i,(*argv)[i]); 
  }
#endif

#if 0
   spi_init();
#endif
  
    QMP_thread_level_t prv;
#ifndef UNIFORM_SEED_NO_COMMS
    QMP_status_t init_status = QMP_init_msg_passing(argc, argv, QMP_THREAD_SINGLE, &prv);
    if (init_status) printf("QMP_init_msg_passing returned %d\n",init_status);
    peRank = QMP_get_node_number();
    peNum = QMP_get_number_of_nodes();
    if(!peRank)printf("QMP_init_msg_passing returned %d\n",init_status);

    if (init_status != QMP_SUCCESS) {
      QMP_error("%s\n",QMP_error_string(init_status));
    }

    // check QMP thread level
    // Added by Hantao
    if(peRank == 0) {
        switch(prv) {
        case QMP_THREAD_SINGLE:
            printf("QMP thread level = QMP_THREAD_SINGLE\n");
            break;
        case QMP_THREAD_FUNNELED:
            printf("QMP thread level = QMP_THREAD_FUNNELED\n");
            break;
        case QMP_THREAD_SERIALIZED:
            printf("QMP thread level = QMP_THREAD_SERIALIZED\n");
            break;
        case QMP_THREAD_MULTIPLE:
            printf("QMP thread level = QMP_THREAD_MULTIPLE\n");
            break;
        default:
            printf("QMP thread level = no idea what this is, boom!\n");
        }
    }

    //Check to make sure that this machine is a GRID machine
    //Exit if not GRID machine
    QMP_ictype qmp_type = QMP_get_msg_passing_type();

    //Get information about the allocated machine
    peNum = QMP_get_number_of_nodes();
    NDIM = QMP_get_allocated_number_of_dimensions();
    peGrid = QMP_get_allocated_dimensions();
    pePos = QMP_get_allocated_coordinates();

    if(peRank==0){
      for(int i = 0; i<*argc;i++){
        printf("argv[%d])(after)=%s\n",i,(*argv)[i]); 
      }
    }
#else
    QMP_status_t init_status = QMP_SUCCESS;
    peRank=0;
    peNum=1;
    NDIM=4;
#endif

//#if (TARGET == BGL) || (TARGET == BGP)
  if (NDIM>5){
    peNum = 1;
    for(int i = 0;i<5;i++)
	peNum *= peGrid[i];
    peRank = peRank % peNum;
  }
  int if_print=1;
  for(int i = 0;i<NDIM;i++)
  if (pePos[i]>=2) if_print=0;

  if (if_print){
      printf("Rank=%d Num=%d NDIM=%d\n",peRank,peNum,NDIM);
      printf("dim:");
      for(int i = 0;i<NDIM;i++)
        printf(" %d",peGrid[i]);
      printf("\n");
      printf("pos:");
      for(int i = 0;i<NDIM;i++)
        printf(" %d",pePos[i]);
      printf("\n");

#if 0
    int rc;
    BGLPersonality pers;
    rts_get_personality(&pers, sizeof(pers));
    printf("from personality: %d %d %d %d\n",pers.xCoord,pers.yCoord,pers.zCoord,rts_get_processor_id());
#endif
  }


//     printf("from personality:\n");

#if 0
    if ( (qmp_type!= QMP_GRID) && (qmp_type !=QMP_MESH)  ) {
      QMP_error("CPS on QMP only implemented for GRID or MESH, not (%d) machines\n",qmp_type);
    }
#endif

//     printf("QMP_declare_logical_topology(peGrid, NDIM)\n");
#ifndef UNIFORM_SEED_NO_COMMS
    //Declare the logical topology (Redundant for GRID machines)
    if (QMP_declare_logical_topology(peGrid, NDIM) != QMP_SUCCESS) {
      QMP_error("Node %d: Failed to declare logical topology\n",peRank);
      exit(-4);
    }
#endif
    initialized = true;
  printf("Rank=%d init_qmp() done\n",peRank);
    
  }
    
}//End namespace QMPSCU {


/*-------------------------------------------------------------------------*/
/* Definitions of the actual emulated SCU system functions                 */
/*-------------------------------------------------------------------------*/

int UniqueID() {return QMPSCU::peRank;}

int CoorX() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::pePos[0];}
int CoorY() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::pePos[1];}
int CoorZ() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::pePos[2];}
int CoorT() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} 
               if (QMPSCU::NDIM>3) return QMPSCU::pePos[3];
               else return 0; }
int CoorS() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} 
               if (QMPSCU::NDIM>4) return QMPSCU::pePos[4];
               else return 0; }
int CoorW(){return 0;}

int SizeX() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::peGrid[0];}
int SizeY() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::peGrid[1];}
int SizeZ() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::peGrid[2];}
int SizeT() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} 
               if (QMPSCU::NDIM>3) return QMPSCU::peGrid[3];
               else return 1; }
int SizeS() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} 
               if (QMPSCU::NDIM>4) return QMPSCU::peGrid[4];
               else return 1; }
int SizeW() {return 1;}

int NumNodes() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::peNum;}

//----------------------------------------------------------------
/*
  The seed can be different for each node and can
  change every time the machine is reset.

  \note The behaviour of the MPI and serial implementations
  may differ from that of the QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
static const unsigned int SERIAL_SEED = 112319;


unsigned int Seed(){return SERIAL_SEED;}   //!< Gets a RNG seed.
//----------------------------------------------------------------
/*
  The seed is the same for each node (spatially fixed, hence the S), but
  can change every time the machine is reset.

  \note The behaviour of the MPI and serial implementations may differ
  from that of the QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int SeedS(){return SERIAL_SEED;}  //!< Gets a RNG seed.
//----------------------------------------------------------------
/*
  SeedT can be different for each node, but is fixed in time (the T), so it is
  unchanged by a reset.

  \note The behaviour of the MPI and serial implementations may
  differ from that of the QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int SeedT(){return SERIAL_SEED;}  //!< Gets a RNG seed.
//----------------------------------------------------------------
/*
  SeedST is the same for each node (spatially fixed, hence the S), and the
  same after every reset (fixed time, hence T).

  \note The behaviour of the MPI and serial implementations may differ
  from that of the QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int SeedST(){return SERIAL_SEED;} //!< Gets a RNG seed.

//----------------------------------------------------------------
/*
  This function blocks further code execution until all
  nodes in the machine have begun executing the code in the sync()
  routine.
  \return 0
*/
//----------------------------------------------------------------
#ifndef HAVE_SYNC
#ifdef UNIFORM_SEED_NO_COMMS
unsigned int sync(){return 1;}
#else
unsigned int sync() {
QMP_status_t sync_status = QMP_barrier(); 
if (sync_status != QMP_SUCCESS) {
      QMP_error("Error in QMP sync:%s\n", QMP_error_string(sync_status));
      return 0;
}
return 1;
}
#endif
#endif

//----------------------------------------------------------------
/*
  On QCDSP this function returns the explicit wire
  number (0 - 7) of the physics direction given by \a dir. In the MPI
  version this returns the internal direction from the cartesian
  communicator which corresponds to the given physics direction.
  \param dir The physics (lattice) direction.
  \return The number used by the comms layer to represents that direction.

  Possibly.
*/
/* In this implementation, this just returns the integer value
  associated with the direction from the SCUDir enum */
//int SCURemap( SCUDir dir ) {
//    return (int)dir;
//}


CPS_END_NAMESPACE
#endif
