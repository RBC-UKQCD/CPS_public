/* begin_generated_IBM_copyright_prolog                             */
/*                                                                  */
/* This is an automatically generated copyright prolog.             */
/* After initializing,  DO NOT MODIFY OR MOVE                       */
/*  --------------------------------------------------------------- */
/*                                                                  */
/* (C) Copyright IBM Corp.  2007, 2007                              */
/* IBM CPL License                                                  */
/*                                                                  */
/*  --------------------------------------------------------------- */
/*                                                                  */
/* end_generated_IBM_copyright_prolog                               */

#ifndef _UPC_H_ // Prevent multiple inclusion
#define _UPC_H_

/*!
 * \file spi/UPC.h
 *
 * \brief Master include file for all BGP UPC System Programming Intefaces.
 */

#include <common/namespace.h>
#include <stdint.h>

#include <spi/UPC_Events.h>

__BEGIN_DECLS

/*
 *  BGP UPC Error Codes
 */
typedef enum BGP_UPC_RC {
  BGP_UPC_SUCCESS = 0,
  BGP_UPC_IS_ACTIVE = -1,
  BGP_UPC_INCOMPATIBLE_COUNTER_FOR_CURRENT_COUNTER_MODE = -2,
  BGP_UPC_BUFFER_PTR_IS_NULL = -3,
  BGP_UPC_BUFFER_SIZE_TOO_SMALL = -4,
  BGP_UPC_INVALID_READ_TYPE = -5,
  BGP_UPC_COUNTER_VALUE_TO_SET_NOT_LESS_THAN_CURRENT_THRESHOLD_VALUE = -6,
  BGP_UPC_COUNTER_THRESHOLD_VALUE_TO_SET_MUST_BE_GREATER_THAN_CURRENT_COUNTER_VALUE = -7,
  BGP_UPC_INVALID_EVENT_EDGE = -8,
  BGP_UPC_INVALID_EVENT_ID = -9,
  BGP_UPC_CANNOT_SPECIFY_ELAPSED_TIME_EVENT = -10,
  BGP_UPC_EVENT_ID_NOT_ACTIVE = -11,
  BGP_UPC_EVENT_ID_CURRENTLY_ACTIVE = -12,
  BGP_UPC_INVALID_USER_MODE = -13
} BGP_UPC_RC_t;

/*
 *  Type definitions, structures, values to pass on BGP UPC interfaces
 */
/*
 *  BGP_UPC_Initialize_Counter_Config first parameter
 *
 *  BGP UPC has 4 basic modes:
 *    modes 0, 1 are the main modes
 *    mode 0 gives info on processors 0 and 1
 *    mode 1 gives info on processors 2 and 3
 *    modes 2, 3 are mostly hold events useful to vhdl designers
 *
 *  BGP_UPC_MODE_DEFAULT is a special mode that is primarily
 *  used at job initialization.  The compute card number is
 *  the component card number from the universal component identifier.
 *  The following rules describe default mode:
 *    Regardless of SMP/DUAL/VN mode:
 *      - even numbered compute cards use BGP_UPC_MODE_0
 *      - odd numbered compute cards use BGP_UPC_MODE_1
 *
 *  NOTE:  When a BG/P job is first loaded onto the system, the system
 *         automatically initializes each UPC unit with essentially
 *         the following API call:
 *
 *         BGP_UPC_Initialize_Counter_Config(BGP_UPC_MODE_DEFAULT, BGP_UPC_CFG_EDGE_DEFAULT)
 *
 *         All the counters are initialized, but are not started.  Therefore, if you want to
 *         take advantage of the default settings, your BG/P code can simply do the
 *         following for each application node:
 *
 *         int main(int argc, char * argv[]) {
 *           .
 *           .
 *           BGP_UPC_Initialize();  // Every application node must first do this
 *                                  // call to issue any other UPC API calls...
 *                                  // Probably best done early in your program.
 *                                  // NOTE:  This initialize call in your code
 *                                  //        is NOT initializing the UPC counters.
 *                                  //        It is initializing global UPC
 *                                  //        variables so that each node
 *                                  //        can issue other UPC API calls
 *                                  //        in a coordinated fashion.
 *           .
 *           .
 *           BGP_UPC_Start(0);      // Start the counters with no reset.
 *                                  // This uses the default initialization of
 *                                  // the UPC units and the counters have already
 *                                  // been reset (passing 0 on call).
 *           Run_the_Code();        // Run the program
 *           BGP_UPC_Stop();        // Stop the counters
 *           .
 *           .
 *           // Process the counters
 *              You could print them...
 *                BGP_UPC_Print_Counter_Values();
 *              -or- read them into a buffer to process...
 *                BGP_UPC_Read_Counter_Values(pBuffer, pBufSize, BGP_UPC_READ_EXCLUSIVE)
 *              NOTE:  Processing will depend on whether the mode is
 *                     SMP/DUAL/VN, as one UPC unit could be counting
 *                     for more than one application node.
 *                     If you are running on more then one node
 *                     (i.e., multiple ranks), then the UPC counters
 *                     for each node will have been initialized as documented above.
 *                     Your code *may* need to concerned with the fact that
 *                     one UPC unit was counting for more than one node.
 *                     If you are only running on one application node
 *                     (i.e., one active rank) for the purpose of simply
 *                     analyzing the performance of your code, then this is of
 *                     no concern and the values you print/read/process
 *                     are the values you want.
 *           }
 */
typedef enum BGP_UPC_Mode {
  BGP_UPC_MODE_MIXED = -2,
  BGP_UPC_MODE_DEFAULT = -1,
  BGP_UPC_MODE_0 = 0,
  BGP_UPC_MODE_1,
  BGP_UPC_MODE_2,
  BGP_UPC_MODE_3,
} BGP_UPC_Mode_t;

/*
 *  BGP_UPC_Initialize_Counter_Config second parameter
 *
 *  BGP UPC has 4 basic edge modes:
 *    High level sensitive
 *    Low to high edge sensitive (Rise)
 *    Hige to low edge sensitive (Fall)
 *    Low level sensitive
 *
 *  The defualt edge (BGP_UPC_CFG_EDGE_DEFAULT) has the
 *  system choose the 'best' edge for a given event.
 */
typedef enum BGP_UPC_Event_Edge {
  BGP_UPC_CFG_LEVEL_HIGH =   0x00000000,
  BGP_UPC_CFG_EDGE_DEFAULT = 0x00000001,
  BGP_UPC_CFG_EDGE_RISE =    0x00000004,
  BGP_UPC_CFG_EDGE_FALL =    0x00000008,
  BGP_UPC_CFG_LEVEL_LOW =    0x0000000C
} BGP_UPC_Event_Edge_t;

/*
 *  BGP_UPC_Read_Counter_Values third parameter
 *
 *  The 2 standard read types for read counter values, shared and exclusive...
 *    BGP_UPC_READ_SHARED reads the counters without stopping
 *      the UPC unit.  Unless the user knows the the UPC unit
 *      is not active, there is no guaranteed atomicity across the
 *      entire read operation.
 *    BGP_UPC_READ_EXCLUSIVE reads the counters with the UPC unit
 *      stopped.  If the read operation had to stop the UPC
 *      unit, the read operation will re-start the UPC unit
 *      after the read option has completed.
 */
typedef enum BGP_UPC_Read_Type {
  BGP_UPC_READ_SHARED     = 0,
  BGP_UPC_READ_EXCLUSIVE  = 1
} BGP_UPC_Read_Type_t;

/*
 * Return structure from BGP_UPC_Read_Counter_Config
 */
typedef struct BGP_UPC_Read_Counter_Config_Struct {
  int32_t rank;                       // Rank
  int32_t core;                       // Core
  int32_t upc_number;                 // UPC Number
  int32_t number_processes_per_upc;   // Number of processes per UPC unit
  BGP_UPC_Event_Id_t event_id;        // Event Id (echoed back)
  BGP_UPC_Event_Edge_t event_edge;    // Event edge
  BGP_UPC_Mode_t mode;                // Mode
  uint32_t reserved_1;                // Reserved for alignment
  char location[24];                  // UPC Location
  uint32_t reserved_2;                // Reserved for alignment
  uint32_t reserved_3;                // Reserved for alignment
  int64_t threshold_value;            // Threshold Value
} BGP_UPC_Read_Counter_Config_Struct_t;

/*
 * Return structure from BGP_UPC_Read_Counters
 */
typedef struct BGP_UPC_Read_Counters_Struct {
  int32_t rank;                       // Rank
  int32_t core;                       // Core
  int32_t upc_number;                 // UPC Number
  int32_t number_processes_per_upc;   // Number of processes per UPC unit
  BGP_UPC_Mode_t mode;                // User mode
  int32_t number_of_counters;         // Number of counter values returned
  char location[24];                  // Location
  int64_t elapsed_time;               // Elapsed time
  uint32_t reserved_1;                // Reserved for alignment
  uint32_t reserved_2;                // Reserved for alignment
  int64_t counter[];                  // Counter values
} BGP_UPC_Read_Counters_Struct_t;


/**
 * BGP_UPC_Check_Active
 *
 *  Check for an active UPC unit
 *
 *  USE:
 *    Returns whether the UPC unit is started.
 *    NOTE:  This is a 'point in time' read of the UPC
 *           active status, as the UPC is not allocated.
 *
 *  INPUT:
 *    None
 *
 *  ERRORS:
 *    None
 *
 *  RETURNS:
 *    0 - UPC unit(s) are stopped
 *    1 - UPC unit(s) are started
 */
extern int32_t BGP_UPC_Check_Active(void);

/**
 * BGP_UPC_Check_Active_Event
 *
 *  Check for an active monitored event
 *
 *  USE:
 *    Returns whether an event is currently being moniitored.
 *
 *  INPUT:
 *    pEventId - Event id
 *
 *  ERRORS:
 *    Invalid event id
 *
 *  RETURNS:
 *    0 - Event not being monitored
 *    1 - Event is being monitored
 *        or error
 */
extern int32_t BGP_UPC_Check_Active_Event(const BGP_UPC_Event_Id_t pEventId);

/**
 * BGP_UPC_Get_Counter_Mode
 *
 *  Returns the current counter mode
 *
 *  USE:
 *    Returns the current counter mode for the UPC
 *    unit(s).
 *    NOTE:  If counter mode 0 is active, then
 *           the currently defined events are
 *           for user modes 0 and 1.  If counter
 *           mode 1 is active, then the currently
 *           defined events are for user modes
 *           2 and 3.  Only events with a user mode
 *           of 0 or 1 can be added with a counter
 *           mode of 0;  only events with a user
 *           mode of 2 and 3 can be added with a
 *           counter mode of 1.
 *    NOTE:  To effectively change the counter mode,
 *           BGP_UPC_Initialize_Counter_Config() must
 *           be run with an appropriate user mode value.
 *    NOTE:  This is a 'point in time' read of the UPC
 *           counter mode value, as the UPC is not allocated.
 *
 *  INPUT:
 *    None
 *
 *  ERRORS:
 *    None
 *
 *  RETURNS:
 *    Counter mode value.
 */
extern int32_t BGP_UPC_Get_Counter_Mode(void);

/**
 * BGP_UPC_Get_Counter_State_GenNum
 *
 *  Get the counter state generation number.
 *
 *  USE:
 *    The counter state generation number is returned.
 *    Each time the counter configuration for the UPC
 *    unit is altered, the counter state generation number
 *    is incremented.  The generation number is returned
 *    on each interface that alters the counter
 *    configuration.  This generation number can be
 *    interrogated between two processing points to
 *    determine if the counter configuration has been
 *    altered.
 *
 *    Currently, the following constitutes a configuration
 *    change:
 *      - a new monitored event.
 *        NOTE:  This includes starting monitoring for an
 *               event that is already being monitored
 *      - a new counter value specified for a monitored event.
 *
 *    If the generation number has changed between
 *    two processing points, then the UPC_Read_Counter_Config
 *    interface can be used to retrieve the current
 *    configuration for a given counter.
 *
 *  INPUT:
 *    None
 *
 *  ERRORS:
 *    None
 *
 *  RETURNS:
 *    Counter state generation number
 */
extern int32_t BGP_UPC_Get_Counter_State_GenNum(void);

/**
 * BGP_UPC_Get_Counter_Threshold_Value
 *
 *  Get the counter threshold value for a given event id.
 *
 *  USE:
 *    Retrieves the threshold value for a given event id.
 *
 *  INPUT:
 *    pEventId - Event id
 *
 *  ERRORS:
 *    Invalid event id
 *    Event id that is not currently active
 *
 *  RETURNS:
 *    Counter threshold value, or error.  If a threshold value is not
 *    in effect for the counter, zero is returned.
 */
extern int64_t BGP_UPC_Get_Counter_Threshold_Value(const BGP_UPC_Event_Id_t pEventId);

/**
 * BGP_UPC_Get_Start_Stop_GenNum
 *
 *  Get the start/stop generation number.
 *
 *  USE:
 *    The start/stop generation number is returned.
 *    Each time the UPC unit is started and/or stopped,
 *    the start/stop generation number is incremented.
 *    This value can be interrogated to determine if
 *    the UPC unit has been started/stopped by another
 *    thread.  (e.g., Save away the generation number
 *    after starting the UPC and later interrogate the
 *    generation number.  If it is the same value, it
 *    is assured that no other thread has stopped and/or
 *    started the UPC since this program last started it.)
 *
 *  INPUT:
 *    None
 *
 *  ERRORS:
 *    None
 *
 *  RETURNS:
 *    Start/stop generation number
 */
extern int32_t BGP_UPC_Get_Start_Stop_GenNum(void);

/**
 * BGP_UPC_Get_Event_Description
 *
 *  Get the event's description in text format.
 *
 *  USE:
 *    The event's description is copied into the
 *    provided buffer, up to the specified length.
 *
 *  INPUT:
 *    pEventId  - Event id
 *    pLength   - Maximum length to copy
 *    pBuffer   - Receiving buffer
 *
 *  ERRORS:
 *    Invalid event id
 *
 *  RETURNS:
 *    Return code or error
 */
extern BGP_UPC_RC_t BGP_UPC_Get_Event_Description(const BGP_UPC_Event_Id_t pEventId, const uint32_t pLength, char* pBuffer);

/**
 * BGP_UPC_Get_Event_Name
 *
 *  Get the event's name (Mnemonic-like).
 *
 *  USE:
 *    The event's short name is copied into the
 *    provided buffer, up to the specified length.
 *
 *  INPUT:
 *    pEventId  - Event id
 *    pLength   - Maximum length to copy
 *    pBuffer   - Receiving buffer
 *
 *  ERRORS:
 *    Invalid event id
 *
 *  RETURNS:
 *    Return code or error
 */
extern BGP_UPC_RC_t BGP_UPC_Get_Event_Name(const BGP_UPC_Event_Id_t pEventId, const uint32_t pLength, char* pBuffer);

/**
 * BGP_UPC_Initialize
 *
 *  Initialize the BGP UPC environment for a given core.
 *
 *  USE:
 *    The UPC environment will be initialized so that other calls
 *    can be made using the UPC SPI interfaces.
 *    Either:
 *      All application nodes must first invoke this initialize
 *      routine before any other UPC SPI interface is called
 *      from any of the nodes
 *        -or-
 *      No application node should invoke this initialize routine.
 *
 *  INPUT:
 *    None
 *
 *  ERRORS:
 *    None
 *
 *  RETURNS:
 *    Nothing
 */
extern void BGP_UPC_Initialize(void);

/**
 * BGP_UPC_Initialize_Counter_Config
 *
 *  Initialize the BGP UPC counter configuration for a given core.
 *
 *  USE:
 *    The UPC will be initialized to monitor the events as defined by
 *    the input user mode.  All counter values are set to zero and the
 *    threshold for each counter is set to inactive.
 *
 *  INPUT:
 *    pMode      - all counters will be configured as defined by this
 *                 user mode
 *    pEventEdge - all counters will be configured as defined using this
 *                 event edge
 *
 *  ERRORS:
 *    Invalid user mode
 *    Invalid event edge
 *    UPC unit is currently active
 *
 *  RETURNS:
 *    Counter state generation number or error
 */
extern int32_t BGP_UPC_Initialize_Counter_Config(const BGP_UPC_Mode_t pUserMode, const BGP_UPC_Event_Edge_t pEventEdge);

/**
 * BGP_UPC_Monitor_Event
 *
 *   Start monitoring for a given event id.
 *
 *  USE:
 *    Set a counter event to be monitored.  The
 *    counter value is set to zero.
 *
 *  INPUT:
 *    pEventId    - Event id
 *    pEventEdge  - Event edge
 *
 *  ERRORS:
 *    Invalid event id
 *    Invalid event edge
 *    UPC unit is currently active
 *    Event id already being monitored
 *    Event id not from correct counter mode
 *
 *  RETURNS:
 *    Counter state generation number or error
 */
extern int32_t BGP_UPC_Monitor_Event(const BGP_UPC_Event_Id_t pEventId, const BGP_UPC_Event_Edge_t pEventEdge);

/**
 * BGP_UPC_Print_Config
 *
 *   Print the UPC configuration.
 *
 *  USE:
 *    Prints the UPC configuration.
 *
 *  INPUT:
 *    None
 *
 *  ERRORS:
 *    None
 *
 *  RETURNS:
 *    Nothing
 */
extern void BGP_UPC_Print_Config(void);

/**
 * BGP_UPC_Print_Counter_Value
 *
 *  Print total count for the specified counter.
 *
 *  USE:
 *    Prints total count for the specified counter.
 *
 *    NOTE:  If a READ_EXCLUSIVE type is being performed and
 *           the UPC unit had to be stopped to perform that
 *           read, the start/stop generation number will be
 *           incremented by one.  If the UPC unit did not
 *           have to be stopped, the start/stop generation
 *           number is not incremented.
 *
 *  INPUT:
 *    pReadType - Read type, shared or exclusive
 *
 *  ERRORS:
 *    Event id that is not currently active
 *    Invalid read type
 *
 *  RETURNS:
 *    Return code or error
 */
extern BGP_UPC_RC_t BGP_UPC_Print_Counter_Value(const BGP_UPC_Event_Id_t pEventId, const int32_t pReadType);

/**
 * BGP_UPC_Print_Counter_Values
 *
 *  Print total counts for all counters which have non-zero total counts.
 *
 *  USE:
 *    Prints total counts for all counters which have non-zero total counts.
 *
 *    NOTE:  If a READ_EXCLUSIVE type is being performed and
 *           the UPC unit had to be stopped to perform that
 *           read, the start/stop generation number will be
 *           incremented by one.  If the UPC unit did not
 *           have to be stopped, the start/stop generation
 *           number is not incremented.
 *
 *  INPUT:
 *    pReadType - Read type, shared or exclusive
 *
 *  ERRORS:
 *    Invalid read type
 *
 *  RETURNS:
 *    Return code or error
 */
extern BGP_UPC_RC_t BGP_UPC_Print_Counter_Values(const int32_t pReadType);

/**
 * BGP_UPC_Read_Counter_Config
 *
 *   Read the counter configuration.
 *
 *  USE:
 *    The current counter configuration is output to the provided buffer
 *    space.  The buffer will only be filled up through the specifed buffer
 *    size.
 *
 *    See structure BGP_UPC_Read_Config_Struct for
 *    the output mapping.
 *
 *  INPUT:
 *    pEventId  - Event id
 *    pBuffer   - Pointer to buffer
 *    pBufSize  - Size of buffer
 *
 *  ERRORS:
 *    Invalid event id
 *    Event id that is not currently active
 *    Buffer pointer not set
 *    Buffer size not equal to or greater than minimum value
 *
 *  RETURNS:
 *    Return code or error
 */
extern BGP_UPC_RC_t BGP_UPC_Read_Counter_Config(const BGP_UPC_Event_Id_t pEventId, void* pBuffer, const uint32_t pBufSize);

/**
 * BGP_UPC_Read_Counter
 *
 *   Read a counter value.
 *   NOTE:  Exactly the same as BGP_UPC_Read_Counter_Value
 *
 *  USE:
 *    Read a counter value for a given event id.
 *
 *    NOTE:  If a READ_EXCLUSIVE type is being performed and
 *           the UPC unit had to be stopped to perform that
 *           read, the start/stop generation number will be
 *           incremented by one.  If the UPC unit did not
 *           have to be stopped, the start/stop generation
 *           number is not incremented.
 *
 *  INPUT:
 *    pEventId  - Event id
 *    pReadType - Read type, shared or exclusive
 *
 *  ERRORS:
 *    Invalid event id
 *    Event id that is not currently active
 *    Invalid read type
 *
 *  RETURNS:
 *    Returns counter value or error
 */
extern int64_t BGP_UPC_Read_Counter(const BGP_UPC_Event_Id_t pEventId, const int32_t pReadType);

/**
 * BGP_UPC_Read_Counters
 *
 *   Read all counter values, and only the counter values.
 *
 *  USE:
 *    Read all current counter values into the specified
 *    buffer.  The buffer will only be filled up through
 *    the specified buffer length.
 *
 *    An array of int64_t values is written to the output buffer.
 *
 *    NOTE:  If a READ_EXCLUSIVE type is being performed and
 *           the UPC unit had to be stopped to perform that
 *           read, the start/stop generation number will be
 *           incremented by one.  If the UPC unit did not
 *           have to be stopped, the start/stop generation
 *           number is not incremented.
 *
 *           If space is provided to retrieve values for
 *           reserved counters, those values are returned
 *           as zero.
 *
 *    NOTE:  Constants BGP_UPC_MINIMUM_LENGTH_READ_COUNTERS_ONLY
 *           and BGP_UPC_MAXIMUM_LENGTH_READ_COUNTERS_ONLY
 *           can be used as the size of buffer.  The former allows
 *           for just enough room for a single counter value.
 *           The latter allows room forall potential counters.
 *           Room for all potential counters is (256 * 8).
 *
 *  INPUT:
 *    pBuffer   - Pointer to buffer
 *    pBufSize  - Size of buffer
 *    pReadType - Read type, shared or exclusive
 *
 *  ERRORS:
 *    Buffer pointer not set
 *    Buffer size not equal to or greater than minimum value
 *    Invalid read type
 *
 *  RETURNS:
 *    Start/stop generation number or error
 */
extern int32_t BGP_UPC_Read_Counters(void* pBuffer, const uint32_t pBufSize, const int32_t pReadType);
/**
 * BGP_UPC_Read_Counter_Value
 *
 *   Read a counter value.
 *
 *  USE:
 *    Read a counter value for a given event id.
 *
 *    NOTE:  If a READ_EXCLUSIVE type is being performed and
 *           the UPC unit had to be stopped to perform that
 *           read, the start/stop generation number will be
 *           incremented by one.  If the UPC unit did not
 *           have to be stopped, the start/stop generation
 *           number is not incremented.
 *
 *  INPUT:
 *    pEventId  - Event id
 *    pReadType - Read type, shared or exclusive
 *
 *  ERRORS:
 *    Invalid event id
 *    Event id that is not currently active
 *    Invalid read type
 *
 *  RETURNS:
 *    Counter value or error
 */
extern int64_t BGP_UPC_Read_Counter_Value(const BGP_UPC_Event_Id_t pEventId, const int32_t pReadType);

/**
 * BGP_UPC_Read_Counter_Values
 *
 *   Read all counter values.
 *
 *  USE:
 *    Read all current counter values into the specified
 *    buffer.  The buffer will only be filled up through
 *    the specified buffer length.
 *
 *    See structure BGP_UPC_Read_Counters_Struct for
 *    the output mapping.
 *
 *    NOTE:  If a READ_EXCLUSIVE type is being performed and
 *           the UPC unit had to be stopped to perform that
 *           read, the start/stop generation number will be
 *           incremented by one.  If the UPC unit did not
 *           have to be stopped, the start/stop generation
 *           number is not incremented.
 *
 *           If space is provided to retrieve values for
 *           reserved counters, those values are returned
 *           as zero.
 *
 *    NOTE:  Constants BGP_UPC_MINIMUM_LENGTH_READ_COUNTERS_STRUCTURE
 *           and BGP_UPC_MAXIMUM_LENGTH_READ_COUNTERS_STRUCTURE
 *           can be used as the size of buffer.  The former allows
 *           for just enough room for the structure's header.
 *           The latter allows for the header plus room for
 *           all potential counters.  Room for all potential
 *           counters is (256 * 8).
 *
 *  INPUT:
 *    pBuffer   - Pointer to buffer
 *    pBufSize  - Size of buffer
 *    pReadType - Read type, shared or exclusive
 *
 *  ERRORS:
 *    Buffer pointer not set
 *    Buffer size not equal to or greater than minimum value
 *    Invalid read type
 *
 *  RETURNS:
 *    Start/stop generation number or error
 */
extern int32_t BGP_UPC_Read_Counter_Values(void* pBuffer, const uint32_t pBufSize, const int32_t pReadType);

/**
 * BGP_UPC_Set_Counter_Value
 *
 *   Set an individual counter value.
 *
 *  USE:
 *    Sets a counter value for a given event id.
 *
 *  INPUT:
 *    pEventId  - Event id
 *    pValue    - Counter value
 *
 *  ERRORS:
 *    Invalid event id
 *    Event id that is not currently active
 *    UPC unit is currently active
 *
 *  RETURNS:
 *    Counter state generation number or error
 */
extern int32_t BGP_UPC_Set_Counter_Value(const BGP_UPC_Event_Id_t pEventId, const int64_t pValue);

/**
 * BGP_UPC_Set_Counter_Threshold_Value
 *
 *   Set an individual counter threshold value.
 *
 *  USE:
 *    Set an counter threshold value for a given
 *    event id.  If a threshold value of zero is passed,
 *    the threshold is set to inactive for the
 *    counter.
 *
 *  ERRORS:
 *    Invalid event id
 *    Event id that is not currently active
 *    UPC unit is currently active
 *
 *  INPUT:
 *    pEventId  - Event id
 *    pValue    - Threshold value
 *
 *  RETURNS:
 *    Return code or error
 */
extern BGP_UPC_RC_t BGP_UPC_Set_Counter_Threshold_Value(const BGP_UPC_Event_Id_t pEventId, const int64_t pValue);

/**
 * BGP_UPC_Start
 *
 *  Start the UPC unit with the current configuration.
 *
 *  USE:
 *    Start the UPC unit with the current configuration.
 *    All future prints/reads for counter values will be
 *    calculated relative to this new starting point.
 *
 *    If the UPC unit is currently active, this operation
 *    is effectively a no-op (i.e., the relative point from
 *    which counter values are reported is NOT altered.)
 *    It is an error case if the UPC unit is already
 *    active and a request is made to start the UPC and
 *    reset the counter values.
 *
 *  INPUT:
 *    pResetCounters - Indicates if counter values are to be
 *                     reset.
 *                       0             - no reset
 *                       anything else - reset
 *
 *  ERRORS:
 *    Resetting the counters is not allowed if the UPC unit is
 *      already active
 *
 *  RETURNS:
 *    Returns the start/stop generation number or error
 */
extern int32_t BGP_UPC_Start(const int32_t pResetCounters);

/**
 * BGP_UPC_Stop
 *
 *  Stops the UPC unit.
 *
 *  USE:
 *    Stops the UPC unit.
 *
 *    If the UPC unit is currently not active, this operation
 *    is effectively a no-op.
 *
 *  ERROR:
 *    None
 *
 *  INPUT:
 *    None
 *
 *  RETURNS:
 *    Returns the start/stop generation number
 */
extern int32_t BGP_UPC_Stop(void);

/**
 * BGP_UPC_Zero_Counter_Value
 *
 *  Zeroes an individual counter value.
 *
 *  USE:
 *    Set the counter value for a given event id to zero.
 *
 *  ERROR:
 *    Event id out of range
 *    Event id is not active
 *
 *  INPUT:
 *    pEventId  - Event id
 *
 *  RETURNS:
 *    Counter state generation number
 */
extern int32_t BGP_UPC_Zero_Counter_Value(const BGP_UPC_Event_Id_t pEventId);

/**
 * BGP_UPC_Zero_Counter_Values
 *
 *  Zeroes all counter values.
 *
 *  USE:
 *    Sets all counter values to zero.
 *
 *  ERROR:
 *    None
 *
 *  INPUT:
 *    None
 *
 *  OUTPUT:
 *    Counter state generation number
 */
extern int32_t BGP_UPC_Zero_Counter_Values(void);

__END_DECLS

#endif
