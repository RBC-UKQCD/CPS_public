#ifndef __QIOGENERAL__
#define __QIOGENERAL__


#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>


#include <util/gjp.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>

#include <comms/scu.h>

#include <alg/do_arg.h>

#include <qio.h>
#include <qioxml.h>

#include <qmp.h>

#define EES_ADDON
#if TARGET == BGQ
#undef EES_ADDON //Chroma stuck with older QIO for now
#endif

#ifdef EES_ADDON
// supported only since qio v2.3.4  (not v2.2.X)
#define SPINCOLORCHECK
#define SPINCOLORINSERT
#endif // EES_ADDON


// some definitions for QIO/QMP
// and DEBUG-options...

#undef  DEBUG_PARANOID

#undef DEBUG_NodeIndex
#undef DEBUG_GetCoords
#undef DEBUG_GetField
#undef DEBUG_PutField
#undef DEBUG_NodeNumber
#undef DEBUG_ReadField
#undef DEBUG_openInput
#undef DEBUG_GlobalWrite
#undef DEBUG_GlobalRead
#undef DEBUG_Init
#undef DEBUG_commTchunk
#undef DEBUG_GetProp
#undef DEBUG_PutProp
#undef DEBUG_ReadPropagator
#undef DEBUG_WritePropagator
#undef DEBUG_ReadSpinColor

#ifdef DEBUG_PARANOID
 #define DEBUG_NodeIndex
 #define DEBUG_GetCoords
 #define DEBUG_GetField
 #define DEBUG_PutField
 #define DEBUG_NodeNumber
 #define DEBUG_ReadField
 #define DEBUG_openInput
 #define DEBUG_GlobalWrite
 #define DEBUG_GlobalRead
 #define DEBUG_Init
 #define DEBUG_commTchunk
 #define DEBUG_GetProp
 #define DEBUG_PutProp
 #define DEBUG_ReadPropagator
 #define DEBUG_WritePropagator
 #define DEBUG_ReadSpinColor
#endif //DEBUG_PARANOID

#undef DEBUG_PAIRRECORD
#undef DEBUG_PAIRRECORD_STD
#ifdef RUN_PAIRRECORD_CHECK
 #define DEBUG_PAIRRECORD
 #undef DEBUG_PAIRRECORD_STD
#endif

#ifdef RUN_PAIRRECORD_STD_CHECK
 #define DEBUG_PAIRRECORD
 #define DEBUG_PAIRRECORD_STD
#endif

#define DO_recordInfo
#define DO_readDimSize


#undef PRINT_checksums
// some bug...???



#define QIO_RW_DIMENSION 4

/* default OUTPUT-format (default input is QIO_UNKNOWN) */
/* one of QIO_UNKNOWN, QIO_SINGLEFILE, QIO_PARTFILE, QIO_MULTIFILE */ 
#if ( TARGET == QCDOC ) || (TARGET == BGP)
 #define QIO_VOLFMT QIO_PARTFILE
#else
 #define QIO_VOLFMT QIO_SINGLEFILE
#endif



/* one of QIO_SERIAL, QIO_PARALLEL */
/* for safety added for QCDOC */

//#if (TARGET == QCDOC) || (TARGET == BGP)
#if 1
  #define QIO_SERPAR QIO_PARALLEL
#else
  #define QIO_SERPAR QIO_SERIAL
#endif

// Uncomment the following to enforce the PARTFILE, which may or may not
// improve the speed of I/O on cluster
//#define QIO_VOLFMT QIO_PARTFILE

// The following three lines FORCE PARTFILE with io-node being smaller than
// total nodes. For RICC and FNAL. 

//#define USE_QIO_SPARSE_PARTFILE

#ifdef USE_QIO_SPARSE_PARTFILE
#define QIO_VOLFMT QIO_PARTFILE
#define QIO_SPARSE_PARTFILE
#define QIO_SERPAR QIO_PARALLEL
// Number of nodes, per which one io-node is designated.
// Set it to zero if you want all nodes to be io-node.

// For FNAL DS
#define QIO_SPARSE_PARTFILE_NODES 32

// For RICC
//#define QIO_SPARSE_PARTFILE_NODES 8

#endif


/* one of QIO_ILDGNO, ??ILDGLAT?? */
 #define QIO_ILDGSTYLE QIO_ILDGLAT 

/* one of QIO_VERB_DEBUG, QIO_VERB_REG, QIO_VERB_MED, QIO_VERB_LOW, QIO_VERB_OFF */
 #define QIO_VERB_LEVEL QIO_VERB_LOW

/* debug-level for qmp */
 #define QMP_VERB_LEVEL 0 


#define TOLERANCE 1e-5 //for Plaq, LinkTr check

/* xml-file header for gauge-config */
#define QIO_XML_FILE_GAUGE "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"

#define MAX_HEADER_LINE 255 //maximal length for char-string ensemble-label,-id
#define QIO_INFO_STRING_MAX 1023 // maximal length (bytes) for info-field

/* xml-file header for propagator */
#define QIO_XML_FILE_PROP "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"

/* spin and color range for propagator-output in QIO */
#define QIO_PROP_SPIN_MAX 4
#define QIO_PROP_COLOR_MAX 3



/***********************************************************************************************/


CPS_START_NAMESPACE
using namespace std;

//! source types
enum QIO_PROP_SOURCE_TYPES {QIO_UNKNOWN_SOURCE=0, QIO_SCALAR_SOURCE, QIO_FULL_SOURCE};

//! initialize everything needed for QIO
#ifndef USE_QIO
//#include <util/qio_dummy.h>
class qio_init {

 private:

  char *cname;
//  QMP_bool_t qmp_run;

 public:

  qio_init(int argc, char* argv[]):cname("qio_init")
    {

      const char * fname = "qio_init()";
      ERR.NotImplemented(cname,fname);
      
    }



    virtual ~qio_init()
      {
	const char * fname = "~qio_init()";
	ERR.NotImplemented(cname,fname);
      }

  void qio_setLayout()
      {
	const char * fname = "qio_setLayout()";
	ERR.NotImplemented(cname,fname);
      }

  void qio_setFilesystem()
      {
	const char * fname = "qio_setFilesystem()";
	ERR.NotImplemented(cname,fname);
      }
  

};

// general functions for reading global data


//void qio_putGlobal(char *buf, size_t index, int count, void *arg);

//void qio_putGlobalSingle(char *buf, size_t index, int count, void *arg);


// writing global data


//void qio_getGlobal(char *buf, size_t index, int count, void *arg);

//void qio_getGlobalSingle(char *buf, size_t index, int count, void *arg);

#else

//! initialize everything needed for QIO
class qio_init {

 private:

  char *cname;
  QMP_bool_t qmp_run;

 public:

  qio_init(int argc, char* argv[]):cname("qio_init")
    {

      #ifdef DEBUG_Init
      printf("qio_init is up\n");
      #endif //DEBUG_Init

      const char * fname = "qio_init()";
      VRB.Func(cname,fname);

      // check for 5th Dim

      
      //if(GJP.Snodes() > 1)
      //{
      //  ERR.General(cname,fname,"QIO-WARNING: QIO (so far) only tested without 5th dim. parallelized. BE CAREFUL!!!\n");
      //exit(-1);
      //}


      
      if(GJP.Snodes() > 1)
	{
	  // check if local-T-size can be splitted by S-nodes
	  int Tchunk=GJP.NodeSites(3)/GJP.Snodes();

	  if( Tchunk == 0 )
	    {
	      ERR.General(cname,fname,"QIO-ERROR: local-T-size too small for number of nodes in S-direction!!!!\n");
	      exit(-1);
	    }

	  if( (Tchunk*GJP.Snodes()) != GJP.NodeSites(3) )
	    {
	      ERR.General(cname,fname,"QIO-ERROR: number of nodes in S-direction is not a divider of locat-T-size!!!\n");
	      exit(-1);
	    }
	}



      // start QMP
      
      qmp_run = QMP_is_initialized();
      QMP_status_t qmp_init;  
      QMP_thread_level_t qmp_level;


      if( !qmp_run )
	{
	  VRB.Result(cname,fname,"QIO starts QMP...\n");

	  qmp_init = QMP_init_msg_passing(&argc, &argv,QMP_THREAD_SINGLE , &qmp_level); 
	  
	  QMP_verbose(QMP_VERB_LEVEL); //this sometimes magically fixes a QMP/QIO-bug...


	  

	}
      else
	{
	  VRB.Result(cname,fname,"QMP was already up !?!?!?\n");

	}
      QIO_verbose(QIO_VERB_LEVEL);


      VRB.FuncEnd(cname,fname);
      
    }



    virtual ~qio_init()
      {
	#ifdef DEBUG_Init
	printf("shut down qio_init\n");
	#endif //DEBUG_Init

	const char * fname = "~qio_init()";
	VRB.Func(cname,fname);

//	QMP_bool_t qmp_run(QMP_is_initialized());
	
	// finish QMP
	
	if(!qmp_run)
	  {
	    VRB.Result(cname,fname," ...finish QMP\n");
	    QMP_finalize_msg_passing();
	  }
	else
	  {
	    VRB.Result(cname,fname," no QMP up to be finished ?!?!?\n");
	  }

	VRB.FuncEnd(cname,fname);

      }


  QIO_Layout layout;
  QIO_Filesystem fs; 
  QIO_Filesystem* pointer_fs; //so that we could sustitue NULL for default case
  
  void qio_setLayout();

  void qio_setFilesystem();  
  

};

// general functions for reading global data


void qio_putGlobal(char *buf, size_t index, int count, void *arg);

void qio_putGlobalSingle(char *buf, size_t index, int count, void *arg);


// writing global data


void qio_getGlobal(char *buf, size_t index, int count, void *arg);

void qio_getGlobalSingle(char *buf, size_t index, int count, void *arg);


#endif // USE_QIO


CPS_END_NAMESPACE
#endif // __QIOGENERAL__
