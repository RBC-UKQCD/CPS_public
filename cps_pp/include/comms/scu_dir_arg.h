#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*! The SCUDirArg Class: scu_dir_arg.h

  Declarations for the MPI implementation of the comms-layer data
  structures.

  A.N.Jackson: ajackson@epcc.ed.ac.uk                       
  -----------------------------------------------------------
  CVS keywords
 
  $Author: mcneile $
  $Date: 2003-06-22 13:34:52 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/scu_dir_arg.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
  $Id: scu_dir_arg.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: scu_dir_arg.h,v $
  $Revision: 1.1.1.1 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/scu_dir_arg.h,v $
  $State: Exp $  */
/*----------------------------------------------------------*/

CPS_END_NAMESPACE
#include<config.h>
CPS_START_NAMESPACE

/* Allow the MPI stuff to be switched out, thus avoiding compiler
   errors (for the time being). */
#ifdef INCLUDE_MPI_SCU

#ifndef INCLUDED_SCU_DIR_ARG_H
#define INCLUDED_SCU_DIR_ARG_H

CPS_END_NAMESPACE
#include <mpi.h>
#include<config.h>
#include<comms/comms/scu_enum.h>
CPS_START_NAMESPACE

/*! enum to handle boolean logic in a neated manner: */
enum { FALSE, TRUE };

/*! Typedefs used to implement portable types               */
typedef enum {
  TYPE_IFloat,                    /*!< identifies IFloat type */
  TYPE_int                       /*!< identifies integer type */
} Type_tag;
#define NUM_INT_TYPES   4                  /*!< The total number of MPI int types */
#define NUM_FLOAT_TYPES 3                  /*!< The total number of MPI IFloat types */

/*! Default size of the basic communications element (in bytes) */
#ifndef COMMS_DATASIZE
#define COMMS_DATASIZE 4
#endif

/*----------------------------------------------------------*/
class SCUDirArg 
{
  private:

    /*----------------------------------
       Private data:
    ----------------------------------*/
    void* addr_;             //!< Base-address of the data to be sent/recieved.
    SCUDir dir_;             //!< Direction in which to communicate.
    SCUXR sendrx_;           //!< Send/receive flag.
    int blklen_;             //!< Length of each block (bytes or things?)
    int numblk_;             //!< Number of blocks.
    int stride_;             //!< Gap from start of one block to the next.
    int bsize_;              //!< Size, in bytes, of the basic data element.
    MPI_Datatype mpi_dt_;    //!< MPI datatype for this data structure.
    int mpi_datasize_;       //!< The basic data-element size.


    /*--------------------------------
       Private methods:
    ----------------------------------*/
    //! This puts initial values into the instance variables:
    void SCUDirArgInitValues(void);
    //! This checks the parameters and creates the datatype.
    void Create_Datatype(void);

  public:
    /*--------------------------------
       Public interface:
    ----------------------------------*/

    //! Default constructor:
    SCUDirArg();
    //! Parameterized constructor:
    SCUDirArg( void* addr, SCUDir dir, SCUXR sendxr, 
	       int blklen, int numblk = 1, int stride = 1);
    //! Destructor:
    ~SCUDirArg();

    //! Initialise, sets up a previously `empty' SCUDirArg:
    void Init( void* addr, SCUDir dir, SCUXR sendxr, 
	       int blklen, int numblk = 1, int stride = 1);

    //! Get the base-address:
    void * Addr();
    //! Set the base-address, returns the previous base-address:
    void * Addr(void* addr);

    //! Get the block-length:
    int Blklen();
    //! Set the block-length, returns the previous block-length:
    int Blklen( int blklen );

    //! Get the number of blocks:
    int Numblk();
    //! Set the number of blocks, returns the previous value:
    int Numblk( int numblk );

    //! Get the stride:
    int Stride();
    //! Set the stride, returns the previous stride value:
    int Stride( int stride );

    //! Set the address, blocklength, #blocks and stride:
    void Reload( void* addr, int blklen, int numblk, int stride );

    //! Get the MPI Datatype created by this SCUDirArg:
    MPI_Datatype Datatype( void );

    //! Get the direction:
    SCUDir CommDir( void );

    //! Get the send/recieve flag of this SCUDirArg:
    SCUXR CommType( void );

    //! Allows the size of the basic data-element (IFloat or int) to be set (in bytes):
    void SetDataSize( int mpi_datasize );

};

#endif

#endif /* INCLUDE_MPI_SCU */
CPS_END_NAMESPACE
