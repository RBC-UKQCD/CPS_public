#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*! The SCUDirArg Class: scu_dir_arg.C

  The MPI implementation of the comms-layer data structures. These
  objects will define and commit a new MPI_Datatype based on the user
  spec, which will be accessed by the sysfunc.C layer.

  A.N.Jackson: ajackson@epcc.ed.ac.uk                       
  -----------------------------------------------------------
  CVS keywords
 
  $Author: mcneile $
  $Date: 2003-06-22 13:34:47 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi_comms/mpi_scu/scu_dir_arg.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
  $Id: scu_dir_arg.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: scu_dir_arg.C,v $
  $Revision: 1.1.1.1 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi_comms/mpi_scu/scu_dir_arg.C,v $
  $State: Exp $  */
/*----------------------------------------------------------*/

CPS_END_NAMESPACE
#include<config.h>
CPS_START_NAMESPACE

/* Allow the MPI stuff to be switched out, thus avoiding compiler
   errors (for the time being). */
#ifdef INCLUDE_MPI_SCU

#ifndef INCLUDED_SCU_DIR_ARG
#define INCLUDED_SCU_DIR_ARG

CPS_END_NAMESPACE
#include<comms/comms/scu_dir_arg.h>
#include<comms/sysfunc.h>
CPS_START_NAMESPACE

//----------------------------------
SCUDirArg::SCUDirArg() {
    // Initialise with default (null) values;
    SCUDirArgInitValues();
}

SCUDirArg::SCUDirArg(void* addr, SCUDir dir, SCUXR sendxr, int blklen, 
		     int numblk, int stride) {
    // Full initialise:
    SCUDirArgInitValues();
    Init( addr, dir, sendxr, blklen, numblk, stride );
    return;
}

SCUDirArg::~SCUDirArg() {
    // Free any associated MPI datatype:
    if( mpi_dt_ != MPI_DATATYPE_NULL ) {
    	MPI_Type_free( &mpi_dt_ ); 
    }
    return;
}

//----------------------------------
void SCUDirArg::SCUDirArgInitValues(void) {
    addr_ = NULL;
    dir_ = SCU_NoDir;
    sendrx_ = SCU_NoXR;
    blklen_ = -1;
    numblk_ = -1;
    stride_ = -1;
    mpi_dt_ = MPI_DATATYPE_NULL;
    mpi_datasize_ = COMMS_DATASIZE;
    return;
}


//----------------------------------
void SCUDirArg::Init(void* addr, SCUDir dir, SCUXR sendxr, int blklen, 
		     int numblk, int stride ) {
    // Copy arguments into the instance variables, and check them.
    SCUDirArgInitValues();
    addr_ = addr;
    dir_ = dir;
    sendrx_ = sendxr;
    blklen_ = blklen;
    numblk_ = numblk;
    stride_ = stride;

    // Attempt to create the MPI datatype;
    Create_Datatype();
    return;
}

//----------------------------------
void* SCUDirArg::Addr() {
    return addr_;
}

void* SCUDirArg::Addr(void* addr) {
    void* oldaddr;

    oldaddr = addr_;
    addr_ = addr;

    return oldaddr;
}

//----------------------------------
int SCUDirArg::Blklen() {
    return blklen_;
}

int SCUDirArg::Blklen( int blklen ) {
    int oldblen;

    oldblen = blklen_;
    blklen_ = blklen;
    Create_Datatype();
    return oldblen;
}

//----------------------------------
int SCUDirArg::Numblk() {
    return numblk_;
}

int SCUDirArg::Numblk( int numblk ) {
    int oldnblk;

    oldnblk = numblk_;
    numblk_ = numblk;
    Create_Datatype();
    return oldnblk;
}

//----------------------------------
int SCUDirArg::Stride() {
    return stride_;
}

int SCUDirArg::Stride( int stride ) {
    int oldstride;

    oldstride = stride_;
    stride_ = stride;
    Create_Datatype();
    return oldstride;
}

//----------------------------------
void SCUDirArg::Reload(void* addr, int blklen, int numblk = 1, int stride = 1) {

    addr_ = addr;
    blklen_ = blklen;
    numblk_ = numblk;
    stride_ = stride;
    Create_Datatype();
    return;
}

//----------------------------------
MPI_Datatype SCUDirArg::Datatype() {
    return mpi_dt_;
}

//----------------------------------
SCUDir SCUDirArg::CommDir() {
    return dir_;
}

//----------------------------------
SCUXR SCUDirArg::CommType() {
    return sendrx_;
}

//----------------------------------
void SCUDirArg::SetDataSize( int mpi_datasize ) {
    //! Change the size (in bytes) of the fundamental data element:
    mpi_datasize_ = mpi_datasize;
    //printf("SetDataSize()\n");
    /*! If there is an mpi datatype defined for this class,
      redefine it to use the new data size: */
    Create_Datatype();
    return;
}

//----------------------------------
void SCUDirArg::Create_Datatype() {
    MPI_Datatype mpitype;

    // Initalize the MPI layer, if necessary:
    if( !commsMPI_init ) SCUCommsInit();

    // Check all the args.

    // Dispose of the previous datatype, if there is one.
    if( mpi_dt_ != MPI_DATATYPE_NULL ) {
	MPI_Type_free( &mpi_dt_ );
    }


    // Define the (new) datatype.
    /* Translate the Type_tag+size_t into an MPI type */
    mpitype = SCUMPITypeConv( TYPE_IFloat, mpi_datasize_ );
    

    if( stride_ <= 1 && numblk_ <= 1 ) {
	// Contiguous:
	MPI_Type_contiguous( blklen_, mpitype, &mpi_dt_ );
    } else {
	// Block-strided, converting the QCDSP stride into an MPI stride:
	MPI_Type_vector( numblk_, blklen_, stride_+blklen_-1, mpitype, &mpi_dt_ );
    }

    // Also commit the type:
    MPI_Type_commit( &mpi_dt_ );
}


#endif

#endif /* INCLUDE_MPI_SCU */
CPS_END_NAMESPACE
