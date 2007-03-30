#ifndef __QIOGENERAL__
#define __QIOGENERAL__

#include <config.h>
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

#include <alg/do_arg.h>

#include <qio.h>

#include <util/qio_xmlinfo.h>

extern "C"
{
#include <util/qio_xml_miss.h>
}

#include <qmp.h>


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
#endif //DEBUG_PARANOID

#define DO_recordInfo
#define DO_readDimSize


#undef PRINT_checksums
// some bug...???



#define QIO_RW_DIMENSION 4


/* default OUTPUT-format (default input is QIO_UNKNOWN) */
/* one of QIO_UNKNOWN, QIO_SINGLEFILE, QIO_PARTFILE, QIO_MULTIFILE */ 
 #define QIO_VOLFMT QIO_PARTFILE

/* one of QIO_SERIAL, QIO_PARALLEL */
 #define QIO_SERPAR QIO_SERIAL

/* one of QIO_ILDGNO, ??ILDGLAT?? */
 #define QIO_ILDGSTYLE QIO_ILDGLAT 

/* one of QIO_VERB_DEBUG, QIO_VERB_REG, QIO_VERB_MED, QIO_VERB_LOW, QIO_VERB_OFF */
 #define QIO_VERB_LEVEL QIO_VERB_LOW

/* debug-level for qmp */
 #define QMP_VERB_LEVEL 0 


#define TOLERANCE 1e-5 //for Plaq, LinkTr check

/* xml-file header for gauge-config */
#define QIO_XML_FILE_GAUGE "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Dummy QCDML</title>"

#define MAX_HEADER_LINE 255 //maximal length for char-string ensemble-label,-id

/***********************************************************************************************/


CPS_START_NAMESPACE
using namespace std;


class qio_init {

 private:

  char *cname;

 public:

  qio_init(int argc, char* argv[]):cname("qio_init")
    {

      #ifdef DEBUG_Init
      printf("qio_init is up\n");
      #endif //DEBUG_Init

      const char * fname = "qio_init()";
      VRB.Func(cname,fname);

      // check for 5th Dim

      
      if(GJP.Snodes() > 1)
	{
	  ERR.General(cname,fname,"QIO-ERROR: QIO (so far) only works without 5th dim. parallelized, SORRY!\n");
	  exit(-1);
	}


      // start QMP
      
      QMP_bool_t qmp_run(QMP_is_initialized());
      QMP_status_t qmp_init;  
      QMP_thread_level_t qmp_level;


      if( !qmp_run )
	{
	  VRB.Flow(cname,fname,"QIO starts QMP...\n");
	  //VRB.Result(cname,fname,"QIO starts QMP...\n");

	  qmp_init = QMP_init_msg_passing(&argc, &argv,QMP_THREAD_SINGLE , &qmp_level); 
	  
	  QMP_verbose(QMP_VERB_LEVEL); //this sometimes magically fixes a QMP/QIO-bug...


	  

	}
      else
	{
	  VRB.Flow(cname,fname,"QMP was already up !?!?!?\n");
	  //VRB.Result(cname,fname,"QMP was already up !?!?!?\n");

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

	QMP_bool_t qmp_run(QMP_is_initialized());
	
	// finish QMP
	
	if(qmp_run)
	  {
	    VRB.Flow(cname,fname," ...finish QMP\n");
	    //VRB.Result(cname,fname," ...finish QMP\n");
	    QMP_finalize_msg_passing();
	  }
	else
	  {
	    VRB.Flow(cname,fname," no QMP up to be finished ?!?!?\n");
	    //VRB.Result(cname,fname," no QMP up to be finished ?!?!?\n");
	  }

	VRB.FuncEnd(cname,fname);

      }


    QIO_Layout layout;

    void qio_setLayout();


};


CPS_END_NAMESPACE
#endif // __QIOGENERAL__
