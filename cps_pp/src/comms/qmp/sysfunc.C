#include<config.h>
/*----------------------------------------------------------*/
/*!\file
  \brief  Definitions for the QMP implementation of the QCDSP/QCDOC communications layer.
  
*/
/*----------------------------------------------------------------------
/* The Sysfunc Comms Interface: sysfunc.C

  The QMP implementation of the QCDOC/QCDSP SCU comms-layer.

  M. Cheng michaelc@phys.columbia.edu
  -----------------------------------------------------------

/*----------------------------------------------------------*/

#include <comms/sysfunc_qmp.h>
//#include <comms/scu_dir_arg.h>
#include <util/qcdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <qmp.h>
CPS_START_NAMESPACE

/*!\namespace cps
  \brief Main namespace for all CPS classes, functions <em>etc</em>.
*/
/*!\namespace cps::QMPSCU
  \brief Namespace for the QMP emulations of the SCU.
*/

//Useful global variables
namespace QMPSCU {
  
  bool initialized = false;
  
  //! Number of grid dimensions.
  static int NDIM = 5;
 
  static int peRank;        /*!< Rank/identify of this  process */
  static int peNum;          /*!< Total number of processors */
  static const int* pePos;  /*!< Position of this process in the grid.*/ 
  static const int* peGrid; /*!< Number of processors in each direction */

  //Clean up resources used by QMP
  void destroy_qmp() {
    QMP_finalize_msg_passing();
  }

  //Initialize QMP with null command line
#if 1
  void init_qmp() {
    int argc=2;
    char *argv[2];
    argv[0] = "-qmp-geom";
    argv[1] = "native";
    init_qmp(&argc, (char ***)&argv);
  }
#endif

  //Initialize QMP
  //Get Allocated machine size, and declare logical machine
  void init_qmp(int * argc, char ***argv) {
//    printf("init_qmp(%d %p)\n",*argc,*argv);
    QMP_thread_level_t prv;
    QMP_status_t init_status = QMP_init_msg_passing(argc, argv, QMP_THREAD_SINGLE, &prv);
    if (init_status != QMP_SUCCESS) {
      QMP_error("%s\n",QMP_error_string(init_status));
    }

    //Check to make sure that this machine is a GRID machine
    //Exit if not GRID machine
    QMP_ictype qmp_type = QMP_get_msg_passing_type();

    //Get information about the allocated machine
    peRank = QMP_get_node_number();
    peNum = QMP_get_number_of_nodes();
    NDIM = QMP_get_allocated_number_of_dimensions();
    peGrid = QMP_get_allocated_dimensions();
    pePos = QMP_get_allocated_coordinates();

    if(peRank==0){
#if 1
      for(int i = 0; i<*argc;i++){
        printf("argv[%d]=%s\n",i,(*argv)[i]); 
      }
#endif
      printf("Rank=%d Num=%d NDIM=%d\n",peRank,peNum,NDIM);
      for(int i = 0;i<NDIM;i++){
        printf("dim %d: %d %d\n",i, peGrid[i],pePos[i]);
      }
    }

#if 1
    if ( (qmp_type!= QMP_GRID) && (qmp_type !=QMP_MESH)  ) {
      QMP_error("CPS on QMP only implemented for GRID or MESH, not (%d) machines\n",qmp_type);
    }
#endif

    //Declare the logical topology (Redundant for GRID machines)
    if (QMP_declare_logical_topology(peGrid, NDIM) != QMP_SUCCESS) {
      QMP_error("Node %d: Failed to declare logical topology\n",peRank);
      exit(-4);
    }
    initialized = true;
    
  }
    
}//End namespace


/*-------------------------------------------------------------------------*/
/* Definitions of the actual emulated SCU system functions                 */
/*-------------------------------------------------------------------------*/

int UniqueID() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::peRank;}

int CoorT() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::pePos[3];}
int CoorX() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::pePos[0];}
int CoorY() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::pePos[1];}
int CoorZ() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::pePos[2];}
int CoorS(){return 0;}

int SizeT() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::peGrid[3];}
int SizeX() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::peGrid[0];}
int SizeY() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::peGrid[1];}
int SizeZ() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::peGrid[2];}
int SizeS() {return 1;}

int NumNodes() {if (!QMPSCU::initialized) {QMPSCU::init_qmp();} return QMPSCU::peNum;}

//----------------------------------------------------------------
/*
  The seed can be different for each node and can
  change every time the machine is reset.

  \note The behaviour of the MPI and serial implementations
  may differ from that of the QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int Seed(){}
//----------------------------------------------------------------
/*
  The seed is the same for each node (spatially fixed, hence the S), but
  can change every time the machine is reset.

  \note The behaviour of the MPI and serial implementations may differ
  from that of the QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int SeedS(){}
//----------------------------------------------------------------
/*
  SeedT can be different for each node, but is fixed in time (the T), so it is
  unchanged by a reset.

  \note The behaviour of the MPI and serial implementations may
  differ from that of the QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int SeedT(){}
//----------------------------------------------------------------
/*
  SeedST is the same for each node (spatially fixed, hence the S), and the
  same after every reset (fixed time, hence T).

  \note The behaviour of the MPI and serial implementations may differ
  from that of the QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int SeedST(){}

//----------------------------------------------------------------
/*
  This function blocks further code execution until all
  nodes in the machine have begun executing the code in the sync()
  routine.
  \return 0
*/
//----------------------------------------------------------------
unsigned int sync() {
QMP_status_t sync_status = QMP_barrier(); 
if (sync_status != QMP_SUCCESS) {
      QMP_error("Error in QMP sync:%s\n", QMP_error_string(sync_status));
}
return 0;
}

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
int SCURemap( SCUDir dir ) {
    return (int)dir;
}


CPS_END_NAMESPACE
