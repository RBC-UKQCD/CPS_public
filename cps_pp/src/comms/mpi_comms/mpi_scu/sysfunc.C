#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*!\file
  \brief  Definitions for the MPI implementation of the QCDSP/QCDOC communications layer.
  
  $Id: sysfunc.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
/*----------------------------------------------------------------------
/* The Sysfunc Comms Interface: sysfunc.C

  The MPI implementation of the QCDSP SCU comms-layer.

  A.N.Jackson: ajackson@epcc.ed.ac.uk                      
  -----------------------------------------------------------
  CVS keywords
 
  $Author: zs $
  $Date: 2003-07-24 16:53:54 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi_comms/mpi_scu/sysfunc.C,v 1.2 2003-07-24 16:53:54 zs Exp $
  $Id: sysfunc.C,v 1.2 2003-07-24 16:53:54 zs Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: sysfunc.C,v $
  $Revision: 1.2 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi_comms/mpi_scu/sysfunc.C,v $
  $State: Exp $  */
/*----------------------------------------------------------*/

CPS_END_NAMESPACE
#include<comms/sysfunc.h>
CPS_START_NAMESPACE

/* Allow the MPI stuff to be switched out, thus avoiding compiler
   errors (for the time being). */
#ifdef INCLUDE_MPI_SCU

#ifndef INCLUDED_SYSFUNC_C
#define INCLUDED_SYSFUNC_C

//extern "C" {

/*! Size of buffers for storing strings */
#define STRING_MAX_LEN  10000              

/*! Boolean false is zero */
#define FALSE           0                  
/*! Boolean true is non-zero */
#define TRUE           -1                  

/*! The total number of MPI int types */
#define NUM_INT_TYPES   4                  
/*! The total number of MPI IFloat types */
#define NUM_FLOAT_TYPES 3                  

/*----------------------------------------------------------*/
/* Static global storage for communication handles etc:     */
/*----------------------------------------------------------*/

extern int commsMPI_init = FALSE; 
/*!<
  The initial value of the global comms-layer
  initialization flag is FALSE, indicating that the communications system has
  not been set up.
*/

extern int commsMPI_Datasize = COMMS_DATASIZE;
/*!<
  Definition of the default size (in bytes) for the basic elements (ints and
  IFloats) to be communicated .
*/ 

static int          commsMPI_peRank;        /*!< Rank/identify of this  process */
static int          commsMPI_peNum;         /*!< Total number of processors */
static int          commsMPI_peGrid[NDIM];  /*!< Run-time defined: no. of  
					     processors in each direction.*/
static int          commsMPI_pePos[NDIM];   /*!< Position of this process in the grid */ 
static int          commsMPI_root;          /*!< Specify the root processor by rank */ 
static MPI_Comm     commsMPI_comm;          /*!< MPI communicator for the NDIM 
					     cartesian PE topology. */
static char         commsMPI_logFileName[STRING_MAX_LEN]; /*!< String containing the 
							   name of the file to be 
							   used for verbose output. */
static FILE         *commsMPI_logFile;      /*!< Pointer to the logfile output stream. */
static unsigned int commsMPI_RNGseed;       /*!< Seed for the RNG */
static char         commsMPI_seedFileName[STRING_MAX_LEN]; /*!< Filename of the file 
							    that holds the list of 
							    RNG seeds. */
MPIRequestManager   *commsMPI_ReqMan;       /*!< Pointer to the MPI request manager */
static int          commsMPI_nnList[2*NDIM];/*!< Look-up table of NN PEs, 
					     indexed by SCUDIR. */
static MPI_Datatype *commsMPI_mpi_dt;        /*!< Current MPI datatype */


/*-------------------------------------------------------------------------*/
/* The actual subroutine definitions:                                      */
/*-------------------------------------------------------------------------*/
int UniqueID() { if( !commsMPI_init ) SCUCommsInit(); return commsMPI_peRank+1; }

int CoorT() { if( !commsMPI_init ) SCUCommsInit(); return commsMPI_pePos[SCU_T]; }
int CoorX() { if( !commsMPI_init ) SCUCommsInit(); return commsMPI_pePos[SCU_X]; }
int CoorY() { if( !commsMPI_init ) SCUCommsInit(); return commsMPI_pePos[SCU_Y]; }
int CoorZ() { if( !commsMPI_init ) SCUCommsInit(); return commsMPI_pePos[SCU_Z]; }

int SizeT() { if( !commsMPI_init ) SCUCommsInit(); return commsMPI_peGrid[SCU_T]; }
int SizeX() { if( !commsMPI_init ) SCUCommsInit(); return commsMPI_peGrid[SCU_X]; }
int SizeY() { if( !commsMPI_init ) SCUCommsInit(); return commsMPI_peGrid[SCU_Y]; }
int SizeZ() { if( !commsMPI_init ) SCUCommsInit(); return commsMPI_peGrid[SCU_Z]; }

int NumNodes() { if( !commsMPI_init ) SCUCommsInit(); return commsMPI_peNum; }

//----------------------------------------------------------------
/*!
  The seed is (in general) different for each node and is (in general)
  changed every time the machine is reset.

  The behaviour of the MPI version may differ from that of the like-named 
  QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int Seed(){ return SCUReadSeedFile(); }
//----------------------------------------------------------------
/*!
  The seed  is the same for each node (spatially fixed, hence the S), but
  changes every time the machine is reset.

  The behaviour of the MPI version may differ from that of the like-named
  QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int SeedS(){ if( !commsMPI_init ) SCUCommsInit(); return commsMPI_RNGseed; }
//----------------------------------------------------------------
/*!
  SeedT is different for each node, but is fixed in time (the T), so it is
  unchanged by a reset.

  The behaviour of the MPI version may differ from that of the like-named 
  QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int SeedT(){ return Seed(); }
//----------------------------------------------------------------
/*!
  SeedST is the same for each node (spatially fixed, hence the S), and the
  same after every reset (fixed time, hence T).

  The behaviour of the MPI version may differ from that of the like-named 
  QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int SeedST(){ return SeedS(); }

//----------------------------------------------------------------
/*!
  This function blocks further code execution until all
  nodes in the machine have begun executing the code in the sync()
  routine.
  \return 0
*/
//----------------------------------------------------------------
unsigned int sync() {
    if( !commsMPI_init ) SCUCommsInit(); 
    MPI_Barrier( commsMPI_comm );
    return 0;
}

//----------------------------------------------------------------
/*!
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


/*! SCUTrans (multiple, overloaded):
  N.B. For each direction, send and recieve requests must be in the
  same order!  */
/*!
  Performs the communication specified in \a arg.
  \param arg The object fully specifying what data to send to (or receive from)
  where.
*/
//----------------------------------------------------------------
void SCUTrans( SCUDirArg * arg ) {
    SCUTrans( &arg, 1 );
}

//----------------------------------------------------------------
/*!
  Performs the multiple communications specified in \a arg.
  \param arg A pointer to an array of objects, each fully specifying what
  data to send to (or receive from) where.
*/
//----------------------------------------------------------------
void SCUTrans( SCUDirArg ** arg, int n ) {
    int i;

    for( i=0; i<n; i++ ) {
	SCUTrans_mpi( arg[i]->Addr(), 
		      arg[i]->Datatype(), 
		      arg[i]->CommDir(), 
		      arg[i]->CommType() );

    }
}

//----------------------------------------------------------------
/*!
  This function does a number of transfers in the same direction of data
  with the same block length, stride and number
  of blocks but with different addresses. These addresses are specified as
  offsets to the base address.
  \param arg The object containing the information about the structure of the
  data, its base address, \e etc.
  \param offset The array of offsets from the base address, defining addresses
  to/from where data will be sent/received.
  \param n The number of data transfers.
*/
//----------------------------------------------------------------
void SCUTrans( SCUDirArg * arg, unsigned int * offset, int n ) {
    int i;

    for( i=0; i<n; i++ ) {
	SCUTrans_mpi( (void*)((unsigned int)arg->Addr() + offset[i]), 
		      arg->Datatype(), 
		      arg->CommDir(), 
		      arg->CommType() );

    }
}

//----------------------------------------------------------------
/*!
  The block length, stride and number of blocks involved in the data
  transfer are given to to the underlying communications layer,
  but no transfers are done.
  \param arg The object containing the information about the structure of the
  data.
*/
//----------------------------------------------------------------
void SCUSetDMA( SCUDirArg * arg ) { SCUSetDMA( &arg, 1 ); }

//----------------------------------------------------------------
/*!
  The block length, stride and number of blocks involved in a number of data
  transfer are given to to the underlying communications layer,
  but no transfers are done.
  \param arg A pointer ot an array of objects containing the information about
  the structure of the data.
  \param n The number of sets of data to be transferred.
//----------------------------------------------------------------*/
void SCUSetDMA( SCUDirArg ** arg, int n ) {
    int i;

    // remove old settings if necessary:
    if( commsMPI_mpi_dt != NULL ) {
	delete[] commsMPI_mpi_dt;
    }

    // Remember the set of datatypes:
    commsMPI_mpi_dt = new MPI_Datatype[n];
    for( i=0; i<n; i++ ) 
      commsMPI_mpi_dt[i] = arg[i]->Datatype();

    return;
}

//----------------------------------------------------------------
/*!
  Performs the communication specified by its arguments.

  \pre The transfer must have been set up using ::SCUTransAddr

  \param arg The object specifiying the base address of the data, the
  direction of the transfer and whether to send of receive the data. 
*/
//----------------------------------------------------------------
void SCUTransAddr( SCUDirArg * arg ) { SCUTransAddr( &arg, 1 ); }

//----------------------------------------------------------------
/*!
  Performs the communications specified by its arguments.

  \pre The transfers must have been set up using ::SCUTransAddr

  \param arg A pointer to an array of objects specifiying, for each transfer,
  the base address of the data, the direction of the transfer and whether to
  send of receive the data.
  \param n The number of sets of data to be transferred.
*/  
//----------------------------------------------------------------
void SCUTransAddr( SCUDirArg ** arg, int n ) {
    int i;

    for( i=0; i<n; i++ ) {
	SCUTrans_mpi( arg[i]->Addr(), 
		      commsMPI_mpi_dt[i], 
		      arg[i]->CommDir(), 
		      arg[i]->CommType() );

    }
}


//----------------------------------------------------------------
/*!
  This function returns only when all pending communications operations
  have completed.
*/
//----------------------------------------------------------------
void SCUTransComplete( void ) {
    int numreq;
    MPI_Status *stat;

    // Check how many requests are pending:
    numreq = commsMPI_ReqMan->NumReq();
    
    if( numreq > 0 ) {
#ifdef VERBOSE
    qfprintf_all(commsMPI_logFile,"SCUTransComplete: numreq = %i\n", numreq );	
#endif
	/* Grab the requires status array */
	stat = new MPI_Status[numreq];
	/* WAIT for all the previously initiated communications to finish */
	MPI_Waitall( numreq, commsMPI_ReqMan->ReqArray(), stat );
	/* Throw away the status */
	delete[] stat;
#ifdef VERBOSE
    qfprintf_all(commsMPI_logFile,"SCUTransComplete: finished.\n", numreq );	
#endif
    }

    // Clear the request manager:
    commsMPI_ReqMan->Clear();

}

//----------------------------------------------------------------
/*!
  This function finds the parameters relevant for the parallel
  decomposition of the lattice, sets up the communications layer
  and defines the grid topology.

  It also  defines a root node (which is useful for IO)
  and opens logfiles (if required).

  \ingroup mpicomms
*/
//----------------------------------------------------------------
void SCUCommsInit( void ) {
    int   idirn, dummy_argc, root_check, *root_array, ir;
    int   grid_periodicity[NDIM];  /* Array used to specify periodic BCs */
    int   pe_reorder = FALSE;      /* Flag to (dis)allow PE reordering for the cart-comm */
    char  **dummy_argv, *PEnum;    /* Dummy argv and PE number string */
    int   ordPEnum;                /* No. of orders of magnitude (base10) of number of PEs */

    //! If we have already been initialized, don't try to do it twice:
    if( commsMPI_init ) return;

    //! Set-up the default values for the comms parameters:
    strcpy(commsMPI_logFileName,"comlog");
    commsMPI_RNGseed = 1;
    strcpy(commsMPI_seedFileName,"rng.dat");

    // Parses the communications parameters:
    MPIParseCommsParam();

    /* Initialise MPI, using dummy argc and argv arguments */
    MPI_Init(&dummy_argc, &dummy_argv);
    commsMPI_init = TRUE;

    /* Look up processor number */
    MPI_Comm_rank( MPI_COMM_WORLD, &commsMPI_peRank );

    /* Look up number of processors */
    MPI_Comm_size( MPI_COMM_WORLD, &commsMPI_peNum );

#ifdef VERBOSE
    /* Initialise the log-file, which may actually be stdout or stderr */
    if( strcmp(commsMPI_logFileName,"stderr") == 0 ) {
	commsMPI_logFile = stderr;
    } else if( strcmp(commsMPI_logFileName,"stdout") == 0 ) {
	commsMPI_logFile = stdout;
    } else {
	/* Create a logfile name with the PE number as a suffix,
	   such that we have logfile.01, logfile.02 ... logfile.15 etc */
	ordPEnum = 1 + (int)log10(commsMPI_peNum);
	PEnum = (char*)malloc( ( ordPEnum + 1 ) * sizeof(char) );
	strcat(commsMPI_logFileName,".");
	/* Create filenumber based on PE number + leading zeros */
	sprintf(PEnum,"%i", ((int)(exp(((double)ordPEnum)*log(10.0))))+commsMPI_peRank);
	strcat(commsMPI_logFileName,&PEnum[1]); /* Skip the leading 1 */
	free( PEnum );
	/* Open the logfile associated with this PE */
	commsMPI_logFile = fopen(commsMPI_logFileName,"w");
    }
    /* Inform user that initialization has completed */
    qfprintf_all(commsMPI_logFile,"SCUCommsInit:  Initializing...\n");
#endif

    /* Define the (cartesian, periodic) topology of the MPI communicator */
    for( idirn = 0; idirn < NDIM; idirn++ ) grid_periodicity[idirn] = TRUE;
    MPI_Cart_create( MPI_COMM_WORLD,   /* Original communicator */
		     NDIM,             /* No. dimensions */
		     commsMPI_peGrid,   /* No. PEs in each direction */
		     grid_periodicity, /* Periodicity of PE grid in each direction */
		     pe_reorder,       /* True/false flag to allow PE re-ordering */
		     &commsMPI_comm     /* The new, cartesian, communicator */
		     );

    /* Look up processor number again, in case the new communicator has
     renumbered the processor ranks */
    MPI_Comm_rank( MPI_COMM_WORLD, &commsMPI_peRank );

    /* Look-up processor position */
    MPI_Cart_coords( commsMPI_comm, commsMPI_peRank, NDIM, commsMPI_pePos );
    
    /* identify the root processor as that which lies at pos[i]=0 forall i */
    /* calculate the identifier on every process */
    root_check = 0;
    for( idirn = 0; idirn < NDIM; idirn++ ) root_check+=commsMPI_pePos[idirn];
    /* Gather the values of root_check from every PE onto every PE */
    root_array = (int*)malloc(commsMPI_peNum*sizeof(int));
    if( root_array == NULL ) 
      SCURaiseError("malloc failed for root_array in SCUCommsInit!");
    MPI_Allgather( &root_check, /* Pointer to number to be gathered */
		   1,           /* i.e. gathering a single item */
		   MPI_INT,     /* which is a standard C integer */
		   root_array,  /* Pointer to the array which will recieve the data */
		   1,           /* One thing from each PE */
		   MPI_INT,     /* and that thing is an int. */
		   commsMPI_comm /* Using the cartesian communicator */
		   );
    /* Every PE goes through the list and identifies the root PE */
    for( ir = 0; ir < commsMPI_peNum; ir++ ) {
	if( root_array[ir] == 0 ) commsMPI_root = ir;
    }
    /* Free the memory associated with the root-checking array */
    free(root_array);

#ifdef VERBOSE
    /* Log that the initialization has completed and give this PEs rank */
    qfprintf_all(commsMPI_logFile,"SCUCommsInit:  Initialization complete [PE=%i of %i, ROOT_PE=%i].\n",commsMPI_peRank, commsMPI_peNum, commsMPI_root );
#endif

    /* Initialise the MPI Request handler */
    commsMPI_ReqMan = new MPIRequestManager();

    /* Initialise the table of NNs, indexed by SCUDir */
    int dir_index, nnMinus;
    for( int idim = 0; idim < NDIM; idim++ ) {
	for( int idir = -1; idir < +3; idir+=2 ) {
	    if( idim == 0 && idir == +1 ) dir_index = SCU_TP;
	    if( idim == 0 && idir == -1 ) dir_index = SCU_TM;
	    if( idim == 1 && idir == +1 ) dir_index = SCU_XP;
	    if( idim == 1 && idir == -1 ) dir_index = SCU_XM;
	    if( idim == 2 && idir == +1 ) dir_index = SCU_YP;
	    if( idim == 2 && idir == -1 ) dir_index = SCU_YM;
	    if( idim == 3 && idir == +1 ) dir_index = SCU_ZP;
	    if( idim == 3 && idir == -1 ) dir_index = SCU_ZM;
	    MPI_Cart_shift( commsMPI_comm,      /* Using the cartesian communicator */
			idim,              /* Do this dimension */
			idir,                /* Look up first neighbour (on +ve side) */
			&nnMinus,/* Store identity of neighbour PEs (-ve direction)*/
			&(commsMPI_nnList[dir_index])     
			    /* Store identity of neighbour PEs (+ve direction) */
			);
	}
    }


}


/*!   \ingroup mpicomms */

void SCUCommsFinalize( void ) {
    if( commsMPI_init ) MPI_Finalize();
}


/* The global summation */
/*! \ingroup mpicomms collectivecomms */

void SCUGlobalSum(Type_tag t,   /*!< In: Type of data being summed */
		  size_t tsize, /*!< In: Size of the data type */
		  int n,        /*!< In: Number of values to sum */
		  void *ivec,   /*!< In: Vector of input values */
		  void *ovec    /*!< Out: Vector of output values */
		  ) {
    MPI_Datatype mpitype; /*!< This will hold the MPI_Datatype for type (t + size) */

    if( !commsMPI_init ) SCUCommsInit(); 

#ifdef VERBOSE
    qfprintf_all(commsMPI_logFile,"COMM_gsum: Performing a global summation.\n");
#endif

    /* Check args make sense */
    if( n <= 0 )
      SCURaiseError("error in COMM_gsum: no. of values to sum is <= 0!");
    if( ivec == NULL )
      SCURaiseError("error in COMM_gsum: input vector points to NULL!");
    if( ovec == NULL )
      SCURaiseError("error in COMM_gsum: output vector points to NULL!");

    /* Map the requested type onto an MPI_Datatype */
    mpitype = SCUMPITypeConv( t, tsize );

    /* Invoke the relevent MPI call, so that all processors get the global sum */
    MPI_Allreduce(ivec,         /*!< Array containing data to be summed */
	          ovec,         /*!< Array to receive the summations */
		  n,            /*!< Number of items in the array */
		  mpitype,      /*!< MPI datatype corresponding to Type_tag */
		  MPI_SUM,      /*!< Do a global sum operation */
		  commsMPI_comm  /*!< Use the cartesian communicator */
		  );

}

/* SCU-layer error handler:
   Should map onto the ERR class for the QCDSP code. */
/*!
  Prints an error message to \c stdout and causes the program to exit
  immediately with the value \a EXIT_FAILURE.
  \param errstr The messsage.

  \ingroup mpicomms  
*/

void SCURaiseError( char* errstr ) {
    /* Report the error: */
    printf("error: %s\n", errstr);

    /* Finish with MPI if it has been initialised: */
    if( commsMPI_init ) MPI_Finalize();

    exit(EXIT_FAILURE);
}


// Extra error wrapper to deal with string literals. */
/*!
  Prints an error message to \c stdout and causes the program to exit
  immediately with the value \a EXIT_FAILURE.
  \param errstr The messsage.

  \ingroup mpicomms
*/
void SCURaiseError( const char* errstring ) { 
    SCURaiseError( const_cast<char*>(errstring) ); 
}


/*-------------------------------------------------------------------------*/
/*                   Implementation-specific subroutines:                  */
/*              If this were a class, these would be private.              */
/*-------------------------------------------------------------------------*/

//----------------------------------------------------------------
/*!
  The lowest level MPI comms subroutine, on which all other comms calls
  are based.  
*/
//----------------------------------------------------------------
void SCUTrans_mpi( void* addr, MPI_Datatype mpi_dt, SCUDir dir, SCUXR sendrx ) {
    int j, nnPE;
    MPI_Request request;

    // Determine the NN in the given direction:
    nnPE = commsMPI_nnList[dir];

    // Initiate the send or recieve:
    if( sendrx == SCU_SEND ) {
	MPI_Issend( addr,            /* base-address of the data */
		    1,               /* Number of items to send, one datatype */
		    mpi_dt,          /* MPI datatype to send */
		    nnPE,            /* ID of destination PE */
		    dir,             /* Message-tag based on dirn */
		    commsMPI_comm,   /* The communicator */
		    &request         /* RETURNS, the request handle */
		    );
    } else {
	MPI_Irecv( addr,                /* base-address of the data */
		   1,                   /* Number of items to recieve, one struct */
		   mpi_dt,              /* MPI datatype to recv */
		   nnPE,                /* ID of source PE */
		   dir-((dir%2)*2-1),   /* Tag based on dirn */
		   commsMPI_comm,       /* The communicator */
		   &request             /* RETURNS, the request handle */
		   );
    }
    // Add the new request to the req. handler:
    commsMPI_ReqMan->AddRequest(request);

    return;
}

//----------------------------------------------------------------
/*!
  Looks up and parses the run-time user parameters specified via COMMS_ENVVAR.
  i.e. Lots of messy string handling et cetera.
 
  \todo Currently, all PEs open the file and look up the required 
  information.  It would perhaps be quicker to get one PE to look 
  in the file and then distribute the information.
  In fact, if only one node is capable of I/O, this would be 
  neccessary, so it should be done.
*/                                                                     
/* ----------------------------------------------------------------- */
void MPIParseCommsParam(void) {
    enum { NULL_READ, GRID_READ, LOGF_READ, SEED_READ, SEEDFILE_READ};
    char  default_filename[STRING_MAX_LEN] = COMMS_DEFFILE;
    char  *comm_def, f_line[STRING_MAX_LEN], *def_token, *tok_pos;
    int   i, def_flag, read_state, idirn, io_state;
    FILE  *fp;

    /* NULL the file pointer in case this routine fails */
    commsMPI_logFile = NULL;

    /* First, get the environment variable that tells us where to look */
    comm_def = getenv(COMMS_ENVVAR);
    if( comm_def == NULL || strlen(comm_def) == 0 ) {
	/* IF not environment variable, use the default filename */
	comm_def = default_filename;
    }

    /* Determine if the string holds the definition info, and set flag to 1 if this is the case */
    def_flag = 0;
    for( i=0; i<strlen(comm_def) && def_flag == 0; i++ ) if( comm_def[i] == '{' ) def_flag = 1;

    /* If flag is zero, assume it is a file pointer, and attempt to load the file into comm_def */
    if( def_flag == 0 ) {
	fp = fopen(comm_def,"r");
	if( fp == NULL ) SCURaiseError("Could not open comms definition file!");
	/* Lookup the file size and define a suitably sized buffer */
	fseek(fp, 0, SEEK_END);
	comm_def = (char*)malloc( ftell(fp) * sizeof(char) );
	fseek(fp, 0, SEEK_SET);
	/* Read the file line-by-line, and put the whole thing into comm_def: */
	if( comm_def == NULL ) SCURaiseError("malloc failed for file buffer in COMM_userParamLookup!");
	strcpy(comm_def,"");
	while( fscanf(fp,"%[^\n]\n",f_line) != EOF ) strcat(comm_def,f_line);
	fclose(fp);
    }

    /* Set initial (null/invalid) values for the user parameters: */
    for( idirn = 0; idirn < NDIM; idirn++ ) commsMPI_peGrid[idirn] = -1;
    strcpy(commsMPI_logFileName,"!!empty!!");

    /* Now attempt to decipher the definition string held in comm_def */
    /* This is done by breaking the string down into a stream of tokens */
    read_state = NULL_READ;
    tok_pos = comm_def;
    while( def_token = (char*)MPICommsStringTokenizer(comm_def, 
							 "{} =,\n\t", 
							 &tok_pos) ) {
		
	/* Look up the number of processors in each direction. */
	/* If we find the `grid' token, change into GRID_READ mode */
	if( strcmp(def_token,"GRID") == 0 ) {
	    read_state = GRID_READ;
	    idirn = 0;
	} else if( read_state == GRID_READ ) { 
	    /* After `grid', read NDIM*ints into the peGrid array */
	    if( idirn < NDIM ) {
		io_state = sscanf(def_token,"%i",&commsMPI_peGrid[idirn]);
		idirn++;
		if( idirn == NDIM ) read_state = NULL_READ; /* Have we finished? */
	    }
	}

	/* Get the name of the logfile for verbose output */
	if( strcmp(def_token,"LOGFILE") == 0 ) {
	    read_state = LOGF_READ;
	} else if( read_state == LOGF_READ ) {
	    /* Grab the filename token */
	    strcpy(commsMPI_logFileName,def_token);
	    read_state = NULL_READ; /* and finish */
	}

	/* Get the specified RNG seed */
	if( strcmp(def_token,"SEED") == 0 ) {
	    read_state = SEED_READ;
	} else if( read_state == SEED_READ ) {
	    io_state = sscanf(def_token,"%i", &commsMPI_RNGseed );
	    read_state = NULL_READ;
	}

	/* Get the specified RNG seeds filename*/
	if( strcmp(def_token,"SEEDFILE") == 0 ) {
	    read_state = SEEDFILE_READ;
	} else if( read_state == SEEDFILE_READ ) {
	    io_state = sscanf(def_token,"%s", &commsMPI_seedFileName );
	    read_state = NULL_READ;
	}

    }

    /* Free the memory associated with the def buffer if required */
    if( def_flag == 0 ) free(comm_def);

    /* If any necessary parameters have not been specified properly, exit */
    /* Checking the processor-element grid specification: */
    for( idirn = 0; idirn < NDIM; idirn++ ) {
	if( commsMPI_peGrid[idirn] < 0 ) {
	    SCURaiseError("Processor array dimensions have not been specified correctly!");
	}
    }

    /* If any optional parameters have not been specified, use defaults */
    /* Checking the log-file definition, defaults to stderr */
    if( strcmp(commsMPI_logFileName,"!!empty!!") == 0 ) {
	strcpy(commsMPI_logFileName,"stderr");
    }

}

//----------------------------------------------------------------
/*!
  String tokenizer, coded here to ensure portability:
*/
//----------------------------------------------------------------

char *MPICommsStringTokenizer(char* str, const char* delim, char** tok_pos ) {
    char *tokenstr, *substr;
    int i, tokenstate, toki, isgap, idel, tok_find;
    
    substr = *tok_pos;

    if( substr[0] != '0' ) {
	// Not at the end of the string, so find the next token:
	tokenstr = (char*)malloc( strlen(substr) );
	tokenstate = 0; toki = 0;
	for( i=0; i<=strlen(substr); i++ ) {
	    // Determine if the current character is one of the delimiters;
	    idel = 0; tok_find = 0;
	    while( idel < strlen(delim)+1 ) { //<The end-of-string 0 is a delimiter.
		if( substr[i] == delim[idel] ) {
		    tok_find = idel+1;
		    idel = strlen(delim)+1;
		}
		idel++;
	    }
	    if( tok_find == 0 ) {
		isgap = 0;
	    } else {
		isgap = 1;
	    }

	    if( i == 0 && isgap == 0 ) tokenstate = 1;
	    if( tokenstate == 0 && isgap == 0 ) {
		    tokenstate = 1;
		    // A token has begun:
		    tokenstr[toki] = substr[i];
		    toki++;
	    } else if( tokenstate == 1 ) {
		if( isgap == 0 ) {
		    // A token continues:
		    tokenstr[toki] = substr[i];
		    toki++;
		} else {
		    // We have found a token, so return it:
		    *tok_pos = &(substr[i]);
		    tokenstr[toki] = 0;
		    return( tokenstr );
		}
	    }
	    //printf("tokenizer: %i %i %i\n", tokenstate, isgap, toki );
	}
    }

    // If the code gets to here we are at the end of the string:
    return(NULL);
}

//----------------------------------------------------------------
/*!
  On-the-fly type+size -> MPI_Datatype conversion.  There are
 probably somewhat quicker ways of doing this via look-up tables, but
 for now, the look-up has been left explicit.
*/
//----------------------------------------------------------------

MPI_Datatype SCUMPITypeConv( Type_tag t, size_t tsize ) {
    int i, mpisize;
    MPI_Datatype int_types[NUM_INT_TYPES] = { MPI_CHAR, MPI_SHORT, MPI_INT, MPI_LONG };
    MPI_Datatype IFloat_types[NUM_FLOAT_TYPES] = { MPI_FLOAT, MPI_DOUBLE, MPI_LONG_DOUBLE };
    char err_str[STRING_MAX_LEN];

    /* Check arguments make sense */
    if( t != TYPE_IFloat && t != TYPE_int )
      SCURaiseError("error in COMM_typeconv: unknown data type!");
    if( tsize <= 0 )
      SCURaiseError("error in COMM_typeconv: size of data type is <= 0!");

    /* Go through int OR IFloat types */
    /* int types */
    if( t == TYPE_int ) {
	for( i=0; i < NUM_INT_TYPES; i++ ) {
	    MPI_Type_size( int_types[i], &mpisize );  /* Get size of this MPI type */
	    if( tsize == mpisize ) {                  /* If we have a match... */
		return( int_types[i] );               /* Return the matched type. */
	    }
	}
	/* if this executes, no suitable type has not been found, so raise an error */
	sprintf(err_str,"no suitable %i-byte int type among MPI primitive types",tsize);
	SCURaiseError(err_str);

    /* IFloat types */
    } else if( t == TYPE_IFloat ) {
	for( i=0; i < NUM_FLOAT_TYPES; i++ ) {
	    MPI_Type_size( IFloat_types[i], &mpisize );  /* Get size of this MPI type */
	    if( tsize == mpisize ) {                    /* If we have a match... */
		return( IFloat_types[i] );               /* Return the matched type. */
	    }
	}
	/* if this executes, no suitable type has not been found, so raise an error */
	sprintf(err_str,"no suitable %i-byte IFloat type among MPI primitive types",tsize);
	SCURaiseError(err_str);

    }

    /* This statement should never execute, however just to check, and keep 
       the compiler happy, we shall say: */
    SCURaiseError("running unrunnable section of SCUMPITypeConv!  Possible memory access problem?");
    return(( MPI_Datatype) NULL );

}

//----------------------------------------------------------------
/*!
  Reads a seed for every PE from a file specified during intialisation.
*/
//----------------------------------------------------------------

unsigned int SCUReadSeedFile( void ) {
    FILE *seedfp = NULL;
    int i, io_state, n = commsMPI_peNum;
    unsigned int *iseed, seed; 

#ifdef VERBOSE
    qfprintf_all(commsMPI_logFile,"SCUReadSeedFile: Opening seed file %s.\n", commsMPI_seedFileName);
#endif

    //! Check we have actually been initialised:
    if( !commsMPI_init ) SCUCommsInit();

    //! Create the seeds buffer:
    iseed = new unsigned int[n+1];
    
    //! Open the file:
    seedfp = fopen(commsMPI_seedFileName, "r" );
    if( seedfp == NULL ) 
      SCURaiseError("could not open seed file!\n");
    
    //! Read in n seeds:
    i = 0; while( i < n && fscanf(seedfp,"%u",&(iseed[i])) != EOF ) i++;
    
    //! Close the file:
    fclose(seedfp); seedfp = NULL;
    
    //! Die if the file ended before all seeds had been read in:
    if( i < n )
      SCURaiseError("not enough seeds have been supplied in the seed file!");
    //! XXX, EXTREME warning.  This killed one thread and then hung.
    
    //! Get the seed which belongs to this PE:
    seed = iseed[commsMPI_peRank];

    //! Delete the seeds buffer:
    delete [] iseed;

    //! Return this PE's seed:
    return seed;
}
    

//}// End of extern "C".
#endif

#endif /* INCLUDE_MPI_SCU */
CPS_END_NAMESPACE
