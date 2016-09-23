#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*!\file
  \brief  Definition of the SCUDirArg class.

  $Id: scu_dir_arg.C,v 1.5 2008/02/08 18:35:06 chulwoo Exp $
*/
/*----------------------------------------------------------*/
/* The SCUDirArg Class: scu_dir_arg.C

  The MPI implementation of the comms-layer data structures. These
  objects will define and commit a new MPI_Datatype based on the user
  spec, which will be accessed by the sysfunc.C layer.

  A.N.Jackson: ajackson@epcc.ed.ac.uk                       
  -----------------------------------------------------------
  CVS keywords
 
  $Author: chulwoo $
  $Date: 2008/02/08 18:35:06 $
  $Header: /space/cvs/cps/cps++/src/comms/mpi/scu/scu_dir_arg.C,v 1.5 2008/02/08 18:35:06 chulwoo Exp $
  $Id: scu_dir_arg.C,v 1.5 2008/02/08 18:35:06 chulwoo Exp $
  $Name: v5_0_16_hantao_io_test_v7 $
  $Locker:  $
  $RCSfile: scu_dir_arg.C,v $
  $Revision: 1.5 $
  $Source: /space/cvs/cps/cps++/src/comms/mpi/scu/scu_dir_arg.C,v $
  $State: Exp $  */
/*----------------------------------------------------------*/



#ifndef INCLUDED_SCU_DIR_ARG
#define INCLUDED_SCU_DIR_ARG

CPS_END_NAMESPACE
#include<comms/scu_dir_arg.h>
#include<comms/sysfunc_cps.h>
CPS_START_NAMESPACE

//----------------------------------
/*!
  The default constructor initialises an instance of this class as follows
  - base address of data to be sent = NULL
  - direction to send data = ::SCU_NoDir (the null direction)
  - send or receive flag = ::SCU_NoXR (neither)
  - length of each block = -1
  - number of blocks = -1
  - stride = -1
  - MPI datatype = MPI_DATATYPE_NULL
  - size of each basic data element = #COMMS_DATASIZE
  .
  \post A data structure with these characteristics is created.  
*/
//----------------------------------
SCUDirArg::SCUDirArg() {
    // Initialise with default (null) values;
    SCUDirArgInitValues();
}

/*!
  \param addr The base address of the data to be sent/received.
  \param dir The direction in which to communicate.
  \param sendxr The flag indicating whether to send or receive this data.
  \param blklen  The length in floating point or integer numbers of each block of data.
  \param numblk The number of blocks (default = 1).
  \param stride The stride (the interval from the start of one block to the next) in floating point or integer numbers (default = 1).

  \post A data structure with the desired characteristics is created.
 */

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
/*!
  \param addr The base address of the data to be sent/received.
  \param dir The direction in which to communicate.
  \param sendxr The flag indicating whether to send or receive this data.
  \param blklen  The length of each block of data in floating point or integer numbers .
  \param numblk The number of blocks.
  \param stride The stride (the interval from the start of one block to the next) in floating point or integer numbers 

  \post A data structure with the desired characteristics is created.
 */
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
/*!
  \return The base address of the data to be communicated.
*/
//----------------------------------
void* SCUDirArg::Addr() {
    return addr_;
}

/*!
  \param addr The new base address of the data to be communicated.
  \return The previous base address
*/

void* SCUDirArg::Addr(void* addr) {
    void* oldaddr;

    oldaddr = addr_;
    addr_ = addr;

    return oldaddr;
}

//----------------------------------
/*!
  \return The block-length in floating point or integer numbers .
*/
//----------------------------------
int SCUDirArg::Blklen() {
    return blklen_;
}

/*!
  \param The new block length in floating point or integer numbers .
  \return The previous block-length.
 */
int SCUDirArg::Blklen( int blklen ) {
    int oldblen;

    oldblen = blklen_;
    blklen_ = blklen;
    Create_Datatype();
    return oldblen;
}

//----------------------------------
/*!
  \return The number of blocks.
*/
//----------------------------------
int SCUDirArg::Numblk() {
    return numblk_;
}

/*!
  \param The new number of blocks.
  \return The previous number of blocks.
*/
int SCUDirArg::Numblk( int numblk ) {
    int oldnblk;

    oldnblk = numblk_;
    numblk_ = numblk;
    Create_Datatype();
    return oldnblk;
}

//----------------------------------
/*!
  \return The stride in floating point or integer numbers.
*/
//----------------------------------
int SCUDirArg::Stride() {
    return stride_;
}

/*!
  \param The new stride in floating point or integer numbers.
  \return The previous stride.
*/
int SCUDirArg::Stride( int stride ) {
    int oldstride;

    oldstride = stride_;
    stride_ = stride;
    Create_Datatype();
    return oldstride;
}

//----------------------------------
/*!
  \param addr The base address of the data to be sent/received.
  \param dir The direction in which to communicate.
  \param sendxr The flag indicating whether to send or receive this data.
  \param blklen The length of each block of data in floating point or integer numbers.
  \param numblk The number of blocks (default = 1).
  \param stride The stride (the interval from the start of one block to the next) (default = 1) in floating point or integer numbers.
*/
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
/*!
  \return The datatype.
*/
//----------------------------------
MPI_Datatype SCUDirArg::Datatype() {
    return mpi_dt_;
}

//----------------------------------
/*!
  \return The direction flag.
*/
//----------------------------------
SCUDir SCUDirArg::CommDir() {
    return dir_;
}

//----------------------------------
/*!
  \return The flag indicating whether this data is sent or received.
*/
//----------------------------------
SCUXR SCUDirArg::CommType() {
    return sendrx_;
}

//----------------------------------
/*!
  If there is an mpi datatype defined for this class,
  it is redefined to use the new data size.
  \param mpi_datasize The new data size in bytes.
*/
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
    if( !MPISCU::Is_Initialised ) MPISCU::CommsInit();

    // Check all the args.

    // Dispose of the previous datatype, if there is one.
    if( mpi_dt_ != MPI_DATATYPE_NULL ) {
	MPI_Type_free( &mpi_dt_ );
    }


    // Define the (new) datatype.
    /* Translate the Type_tag+size_t into an MPI type */
    mpitype = MPISCU::MPITypeConv( TYPE_IFloat, mpi_datasize_ );
    

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


CPS_END_NAMESPACE
