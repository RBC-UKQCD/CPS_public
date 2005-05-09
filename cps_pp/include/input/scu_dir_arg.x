#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*!\file
  \brief  Declaration of the SCUDirArg class.
  
  $Id: scu_dir_arg.x,v 1.3 2005-05-09 07:14:56 chulwoo Exp $
*/
/*----------------------------------------------------------*/
/* The SCUDirArg Class: scu_dir_arg.h

  Declarations for the MPI implementation of the comms-layer data
  structures.

  A.N.Jackson: ajackson@epcc.ed.ac.uk                       
  -----------------------------------------------------------
  CVS keywords
 
  $Author: chulwoo $
  $Date: 2005-05-09 07:14:56 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/input/scu_dir_arg.x,v 1.3 2005-05-09 07:14:56 chulwoo Exp $*/
  $Id: scu_dir_arg.x,v 1.3 2005-05-09 07:14:56 chulwoo Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: scu_dir_arg.x,v $
  $Revision: 1.3 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/input/scu_dir_arg.x,v $*/
  $State: Exp $  */
/*----------------------------------------------------------*/

CPS_END_NAMESPACE
CPS_START_NAMESPACE


#ifndef INCLUDED_SCU_DIR_ARG_H
#define INCLUDED_SCU_DIR_ARG_H

CPS_END_NAMESPACE
#include <mpi.h>
#include <comms/scu_enum.h>
CPS_START_NAMESPACE

/**! Labels for Boolean logic values. */
enum { FALSE, TRUE };



/*! Flags to identify  portable data types               */
typedef enum {
  TYPE_IFloat,                   /*!< Identifies ::IFloat type */
  TYPE_int                       /*!< Identifies integer type */
} Type_tag;

#define NUM_INT_TYPES   4        /*!< The total number of MPI integer types */
#define NUM_FLOAT_TYPES 3        /*!< The total number of MPI IFloat types */

/*! Default size of the basic communications element (in bytes) */
#ifndef COMMS_DATASIZE
#define COMMS_DATASIZE 4
#endif


/*----------------------------------------------------------*/
/*! A class describing the data to be communicated,*/
/*!
   The SCUDirArg class is part of the MPI communications layer. It
   describes the data to be communicated. It is assumed to be of a regularly
   strided form, so it is described in terms of its base address, the number
   of blocks its block length and the stride between each block. The latter
   two are in units of the actual numbers (not bytes).

   Also described here is whether this data should be sent or received and
   the direction in which the data is to communicated.

   At present this code is hardwired to assume that the data is always
   of type IFloat.

   \ingroup comms    
*/
/*----------------------------------------------------------------------*/
class SCUDirArg 
{
  private:

    /*----------------------------------
       Private data:
    ----------------------------------*/
    void* addr_;             /*!< Base-address of the data to be sent/recieved.*/
    SCUDir dir_;             /*!< Direction in which to communicate.*/
    SCUXR sendrx_;           /*!< Send/receive flag.*/
    int blklen_;             /*!< Length of each block (bytes or things?)*/
    int numblk_;             /*!< Number of blocks.*/
    int stride_;             /*!< Gap from start of one block to the next.*/
    int bsize_;              /*!< Size, in bytes, of the basic data element.*/
    MPI_Datatype mpi_dt_;    /*!< MPI datatype for this data structure.*/
    int mpi_datasize_;       /*!< The basic data-element size.*/


    /*--------------------------------
       Private methods:
    ----------------------------------*/
    /*! This puts initial values into the instance variables:*/
    void SCUDirArgInitValues(void);
    /*! This checks the parameters and creates the datatype.*/
    void Create_Datatype(void);

  public:
    /*--------------------------------
       Public interface:
    ----------------------------------*/

    /*! Default constructor*/
    SCUDirArg();
    /*! Parameterized constructor*/
    SCUDirArg( void* addr, SCUDir dir, SCUXR sendxr, 
	       int blklen, int numblk = 1, int stride = 1);
    /*! Destructor:*/
    ~SCUDirArg();

    /*! Initialise (or re-initialise) the datatype.*/
    void Init( void* addr, SCUDir dir, SCUXR sendxr, 
	       int blklen, int numblk = 1, int stride = 1);

    /*! Get the base-address.*/
    void * Addr();
    /*! Set the base-address.*/
    void * Addr(void* addr);

    /*! Get the block-length.*/
    int Blklen();
    /*! Set the block-length.*/
    int Blklen( int blklen );

    /*! Get the number of blocks.*/
    int Numblk();
    /*! Set the number of blocks.*/
    int Numblk( int numblk );

    /*! Get the stride.*/
    int Stride();
    /*! Set the stride.*/
    int Stride( int stride );

    /*! Reset the data structure parameters parameters.*/
    void Reload( void* addr, int blklen, int numblk, int stride );

    /*! Get the MPI Datatype created by this SCUDirArg.*/
    MPI_Datatype Datatype( void );

    /*! Get the direction.*/
    SCUDir CommDir( void );

    /*! Get the send/receive flag of this SCUDirArg.*/
    SCUXR CommType( void );

    /*! Change the size of the fundamental data element.*/
    void SetDataSize( int mpi_datasize );

};


#endif



CPS_END_NAMESPACE
