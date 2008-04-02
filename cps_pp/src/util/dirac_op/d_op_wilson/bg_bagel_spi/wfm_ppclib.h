/* File : ppcLib.h 
   Purpose : A collection of the C routines, macros and structs 
             provided by ppcLib.

   Author: bj
           PAB, add dcbi,dcbt,dcbf,icbi,icbt,iccci,dccci Aug 2003
   Date  : 25 - Sept. 2001 

   Comments: 
*/

#ifndef __ppc_lib_h__
#define __ppc_lib_h__

#ifdef __cplusplus
extern "C" {
#endif
/* Macros for setting and retrieving special purpose registers */

#define stringify(s)    tostring(s)
#define tostring(s)     #s

  /* -----                                      ----- */
  /* -----      MSR Accessor functions (macros) ----- */
  /* -----                                      ----- */

 /*! Input a word from a MSR 
   *
   * This function inputs a word (unsigned long) from a 
   * special purpose register 
   *
   * \param the SPR number
   * 
   * \returns the value of the register 
   */

#define mfmsr()         ({unsigned int rval; \
                        asm volatile("mfmsr %0" : "=r" (rval)); rval;})

  /*! Output a word into the MSR 
   * 
   * This function outputs a word (unsigned long) to a 
   * special purpose register.
   *
   * \param SPR number
   * \param the word to be output
   */
#define mtmsr(v)        asm volatile("mtmsr %0" : : "r" (v))


  /* -----                                      ----- */
  /* -----      SPR Accessor functions(macros)  ----- */
  /* -----                                      ----- */
  
  /*! Output a word into a special purpose register 
   * 
   * This function outputs a word (unsigned long) to a 
   * special purpose register.
   *
   * \param SPR number
   * \param the word to be output
   */
#define mtspr(rn, v)    asm volatile("mtspr " stringify(rn) ",%0" : : "r" (v))  

  /*! Input a word from a special purpose register 
   *
   * This function inputs a word (unsigned long) from a 
   * special purpose register 
   *
   * \param the SPR number
   * 
   * \returns the value of the register 
   */
#define mfspr(rn)       ({unsigned int rval; \
                        asm volatile("mfspr %0," stringify(rn) \
                                     : "=r" (rval)); rval;})


  /* -----                                     ----- */
  /* -----      DCR Accessor functions         ----- */
  /* -----                                     ----- */
  /*! Output an unsigned long into a DCR Register 
   *  
   * This function outputs an unsigned long to a DCR Register
   * \param DCR Register Number 
   * \param The value to place in the DCR Register
   *
   * \returns nothing
   */
#define mtdcr(rn, v)    asm volatile("mtdcr " stringify(rn) ",%0" : : "r" (v))

  /*! Input from a DCR Register 
   * 
   * This function gets a value from a DCR register
   * \param DCR Register Number 
   * 
   * \returns the contents of the DCR Ring
   */
#define mfdcr(rn)       ({unsigned int rval; \
                        asm volatile("mfdcr %0," stringify(rn) \
                                     : "=r" (rval)); rval;})

  /* -----                                     ----- */
  /* ----- Memory Mapped IO accessor functions ----- */
  /* -----                                     ----- */


  /*! Input a byte from a memory mapped register 
   *
   * This function reads a byte from a memory mapped register.
   * \param address an unsigned long value containing the address
   * of the device
   *
   * \returns an unsigned char containing the byte read 
   */
  unsigned char   in8(unsigned long address);

  /*! Output a byte to a memory mapped register 
   *
   * This function writes a byte to a memory mapped register.
   * \param address an unsigned long value containing the address
   * of the device
   *
   * \param value is an unsigned character containing the byte
   * to be written to the device 
   */
  void            out8(unsigned long address, unsigned char value);
  
  /*! Input a 16bit unsigned short from a memory mapped register 
   *
   * This function reads a 16 bit long value (unsigned short)
   * from a memory mapped register.
   * \param address an unsigned long value containing the address
   * of the device
   *
   * \returns an unsigned short containing the the value read.
   */
  unsigned short  in16(unsigned long address);
  
  /*! Output a 16bit unsigned short to a memory mapped register 
   *
   * This function writes a 16 bit long value (unsigned short)
   * to memory mapped register.
   * \param address an unsigned long value containing the address
   * of the device
   * \param value an unsigned short containing the value to be written
   */
  void            out16(unsigned long address, unsigned short value);
  
  /*! Input a 32bit word from a memory mapped register 
   *
   * This function reads a 32 bit long value (unsigned long)
   * from a memory mapped register.
   * 
   * \param address an unsigned long value containing the address
   * of the device
   *
   * \returns an unsigned long containing the the value read.
   */
  unsigned long   in32(unsigned long address);
  
  /*! Output a 32bit word to a memory mapped register 
   *
   * This function writes a 32 bit long value (unsigned long)
   * to a memory mapped register.
   * 
   * \param address an unsigned long value containing the address
   * of the device
   * \param value an unsigned long value containing the value 
   *        to be written 
   */
  void            out32(unsigned long address, unsigned long value);
  
  /*! Input a 16bit unsigned short from a memory mapped register 
   * and reverse byte order 
   *
   * This function reads a 16 bit long value (unsigned short)
   * from a memory mapped register and returns it with its
   * byte order reversed.
   *
   * \param address an unsigned long value containing the address
   * of the device
   *
   * \returns an unsigned short containing the the value read
   *         after reversing its byte order
   */
  unsigned short  in16r(unsigned int);
  
  /*! Output a 16bit unsigned short to a memory mapped register 
   * and reverse byte order
   *
   * This function writes a 16 bit long value (unsigned short)
   * to memory mapped register, after reversing the byte order
   * of the value to be written.
   *
   * \param address an unsigned long value containing the address
   * of the device
   * \param value an unsigned short containing the value to be written
   * (in original byte order)
   */
  void            out16r(unsigned int, unsigned short value);
  
  /*! Input a 32bit word from a memory mapped register and reverse byte
   *  order
   *
   * This function reads a 32 bit long value (unsigned long)
   * from a memory mapped register, and reverses its byte 
   * order. 
   * 
   * \param address an unsigned long value containing the address
   * of the device
   *
   * \returns an unsigned long containing the the value read after
   *  having its bytes reversed.
   */
  unsigned long   in32r(unsigned int);
  
  /*! Output a 32bit word to a memory mapped register and reverse
   * the byte order.
   *
   * This function writes a 32 bit long value (unsigned long)
   * to a memory mapped register after reversing its byte order.
   * 
   * \param address an unsigned long value containing the address
   * of the device
   * \param value an unsigned long value containing the value 
   *        to be written (original byte order)
   */
  void            out32r(unsigned int, unsigned long value);


  /*! \brief Time Base Structure
   *
   *  This function defines the structure, that can be used
   *  to query the processor timebas registers 
   */
  typedef struct {                       
    unsigned long   tb_tbu;		/*!< Time-Base Upper                 */
    unsigned long   tb_tbl;		/*!< Time-Base Lower                 */
  } tb_t;

  /*!  Execute Illegal Instruction
   *
   * This is a debugging function that should produce an illegal 
   * opcode.
   */
  void		ppcAbend(void);

  /*! And MSR with value
   *
   * This function performs a Boolean AND on the BITS in the
   * MSR (machine state register) and the supplied word and
   * The value of the MSR bits before the AND is returned
   * 
   * \params value an unsigned long to AND the MSR with
   * \returns an unsigned long containing the MSR bits before the AND.
   */
  unsigned long	ppcAndMsr(unsigned long value);


  /*! Halt Processor 
   * 
   * This function will cause the processor to go into
   * an infinite loop thereby preventing it from performing
   * any more useful computation. This is to simulate
   * a processor HALT instruction
   */
  void		ppcHalt(void);


  /*! Instruction Cache Synchronize 
   *
   * This function will cause the execution of an isync
   * instruction to synchronize the I-Cache 
   */
  void		ppcIsync(void);


  /*! Read the Time Base
   *
   * This function will read the time base registers
   * TBU and TBL.
   *
   * \param clock_data a pointer to a tb_t. The results
   * are written into this struct
   */
void		ppcMftb(tb_t *clock_data);

  /*! Write the Time Base
   * 
   * This function will write the time base registers
   * TBU and TBL (but only in Supervisor mode)
   *
   * \param clock_data a pointer pointing to a tb_t.
   * The contents of the struct will be written to the register
   */
  void		ppcMttb(tb_t *clock_data);


  /*! Or With MSR
   * 
   * This function will OR the bits in the input parameter
   * word with the contents of the MSR, leaving the results
   * in the MSR. The contents of the MSR from before the OR
   * are returned 
   *
   * \param valua an unsigned long containing the bits to OR with
   */
  unsigned long	ppcOrMsr(unsigned long value);
  
  /* XLC provides mechanisms for doing sync and eieio inline  */
  /* in the interest of portability, these are function calls  */

  /*! Processor Synchronise 
   *
   * This function causes the processor to execute a sync instruction.
   *
   */
  void		ppcSync(void);

  /*! Enforce in-order execution of IO
   *
   * This function causes the processor to execute an EIEIO instruction
   *
   */
  void 		ppcEieio(void);

  /****************************************************************
   * Aug 2003 additions, Pab
   ****************************************************************/
  /*! DCBT load cache line
   */
#define dcbt(A) asm volatile ( "dcbt %0 , %1" : : "r" (0)  , "r" (A) )

  /*! DCBI invalidate cache line
   */
#define dcbi(A) asm volatile ( "dcbi %0 , %1" : : "r" (0)  , "r" (A) )

  /*! DCBF flush cache line
   */
#define dcbf(A) asm volatile ( "dcbf %0 , %1" : : "r" (0)  , "r" (A) )

  /*! dccci invalidate Dcache 
   */
#define dccci(A) asm volatile ( "dccci %0 , %1" : : "r" (0)  , "r" (0) )


  /*! ICBT load Icache line
   */
#define icbt(A) asm volatile ( "icbt %0 , %1" : : "r" (0)  , "r" (A) )

  /*! ICBI invalidate Icache line
   */
#define icbi(A) asm volatile ( "icbt %0 , %1" : : "r" (0)  , "r" (A) )

  /*! iccci invalidate Icache 
   */
#define iccci(A) asm volatile ( "iccci %0 , %1" : : "r" (0)  , "r" (0) )

  /*! Write TLB
   */
#define tlbwe(E,W,N) asm volatile ( "tlbwe %0 , %1 , " stringify(N)  : : "r" (W)  , "r" (E) )

  /*! Read TLB
   */
#define tlbre(E,N)       ({unsigned int rval; \
                        asm volatile("tlbre %0, %1, " stringify(N) \
                                     : "=r" (rval) : "r" (E)); rval;})

#ifdef __cplusplus
}
#endif


#endif
