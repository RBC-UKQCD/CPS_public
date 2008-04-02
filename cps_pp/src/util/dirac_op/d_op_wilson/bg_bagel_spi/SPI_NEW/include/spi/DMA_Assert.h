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

#ifndef __DMA_ASSERT_H_ /* Prevent multiple inclusion */
#define __DMA_ASSERT_H_

/*!
 * \file spi/DMA_Assert.h
 *
 * \brief DMA SPI Assert Macros
 *
 * Two sets of assert macros are provided:
 * - Kernel Asserts
 * - User-mode Asserts
 *
 * When DMA SPIs are used within the kernel, a special assert routine is called
 * that does NOT abort.  It just prints the assertion and the location and
 * continues.
 *
 * When DMA SPIs are used within user-mode code, the normal assert routine is
 * called, which prints the assertion and location and aborts.
 *
 * Several levels of asserts are provided, and #define variables control which
 * levels are activated.  The following assert macros are available:
 *
 * SPI_abort  - Always active and always issues assert(0).
 *              Primarily used for unimplemented code paths.
 *              Not available in the kernel.
 * SPI_assert - Active by default, or when ASSERT_PROD is defined.
 *              Meant to flag user errors.
 * SPI_assert_debug - Active by default.  Meant to flag coding
 *                    errors before shipping.
 *
 * The following #defines control which level of asserts are compiled into
 * the code.  Only one of ASSERT_ABORT, ASSERT_PROD (or nothing) should
 * be specified.
 * - ASSERT_ABORT means that the "abort" level is the only level
 *   of asserts that is active.  Other levels are turned off.
 * - ASSERT_PROD means that "abort" and "assert" levels are active.
 *   "assert_debug" is turned off.
 * - Not specifying ASSERT_ABORT or ASSER_PROD means that all
 *   levels of asserts ("abort", "assert", "assert_debug") are
 *   active.
 */

#include <common/namespace.h>


__BEGIN_DECLS


#include <stdio.h>

/* ============================================================ */

#ifdef __CNK__

/*!
 * \brief Production-level Kernel Assert.
 *
 * This production level of assert will be active during normal production
 * code execution.
 *
 * When in the kernel, just do a printf, but don't exit.
 */
#define SPI_assert(x)       DMA_KernelAssert(x)

/*!
 * \brief Debug-level Kernel Assert.
 *
 * This debug level of assert will only be active during in-house debugging.
 *
 * When in the kernel, just do a printf, but don't exit.
 */
#define SPI_assert_debug(x) DMA_KernelAssert(x)

#ifdef NDEBUG

/*!
 * \brief No Debug Kernel Assert Internal Macro
 *
 * This macro is used internally for when asserts are turned off via the NDEBUG
 * flag.  It does nothing.
 */
#define DMA_KernelAssert( __assert_test ) ((void)0)

/* ============================================================ */

#else /* not NDEBUG */

/*!
 * \brief Kernel Assert Internal Function
 *
 * This function is called when the kernel determines that it needs to assert.
 * It prints the assertion that failed and the code location, but does not
 * abort.  The kernel should continue executing.
 *
 * \param[in]  Pointer to the assertion string that failed the test
 * \param[in]  Pointer to the name of the source file that coded the assert
 * \param[in]  Line number within the source file that coded the assert
 */
extern inline void __DMA_KernelAssert( const char *__assertion,
				       const char *__file,
				       int __line )
{
   printf("Assertion Failed: %s, file %s, line %d.\n",
	  __assertion,
	  __file,
	  __line );
}


/*!
 * \brief Kernel Assert Internal Macro
 *
 * This macro is used internally when asserts are turned on (the NDEBUG flag
 * is not specified).  It tests the assertion.  If the assertion is true, it
 * does nothing.  If the assertion is false, it invokes the __DMA_KernelAssert
 * internal function to print out the assert information.
 *
 * \param[in]  Pointer to the assertion string that failed the test
 * \param[in]  Pointer to the name of the source file that coded the assert
 * \param[in]  Line number within the source file that coded the assert
 */
#define DMA_KernelAssert( __assert_test ) \
                           ((__assert_test) ? ((void)0) : \
                           __DMA_KernelAssert( #__assert_test, __FILE__, __LINE__ ))


#endif /* NDEBUG */

/* ============================================================ */

#else /* not __CNK__ */

#include <assert.h>

#ifdef ASSERT_ABORT

/*!
 * \brief Abort-level Abort Assert
 *
 * This macro is defined when the ASSERT_ABORT level of asserts is active.
 *
 * This macro will assert(0).
 *
 */
#define SPI_abort()         assert(0)

/*!
 * \brief Abort-level Production Assert
 *
 * This macro is defined when the ASSERT_ABORT level of asserts is active.
 * This macro will not assert.  It will simply execute the assert test, but
 * because abort-level-only asserts are active, it will not assert.
 *
 */
#define SPI_assert(x)

/*!
 * \brief Abort-level Debug Assert
 *
 * This macro is defined when the ASSERT_ABORT level of asserts is active.
 * This macro will not assert.  It will simply execute the assert test, but
 * because abort-level-only asserts are active, it will not assert.
 *
 */
#define SPI_assert_debug(x)

/* ============================================================ */

#else /* Not ASSERT_ABORT */

#ifdef ASSERT_PROD

/*!
 * \brief Production-level Abort Assert
 *
 * This macro is defined when the ASSERT_PROD level of asserts is active.
 *
 * This macro will assert(0).
 *
 */
#define SPI_abort()         assert(0)

/*!
 * \brief Production-level Production Assert
 *
 * This macro is defined when the ASSERT_PROD level of asserts is active.
 *
 * This macro invokes the standard assert() function with the specified
 * assert test.
 */
#define SPI_assert(x)       assert(x)

/*!
 * \brief Production-level Debug Assert
 *
 * This macro is defined when the ASSERT_PROD level of asserts is active.
 *
 * This macro will not assert.  It will simply execute the assert test, but
 * because production-level-only asserts are active, it will not assert.
 */
#define SPI_assert_debug(x)

/* ============================================================ */

#else /* Not ASSERT_PROD */

/*!
 * \brief Debug-level Abort Assert
 *
 * This macro is defined when all levels of asserts are desired (neither the
 * ASSERT_ABORT nor ASSERT_PROD level of asserts is active.  This is the
 * default).
 *
 * This macro will assert(0).
 *
 */
#define SPI_abort()         assert(0)

/*!
 * \brief Debug-level Production Assert
 *
 * This macro is defined when all levels of asserts are desired (neither the
 * ASSERT_ABORT nor ASSERT_PROD level of asserts is active.  This is the
 * default).
 *
 * This macro invokes the standard assert() function with the specified
 * assert test.
 */
#define SPI_assert(x)       assert(x)

/*!
 * \brief Debug-level Debug Assert
 *
 * This macro is defined when all levels of asserts are desired (neither the
 * ASSERT_ABORT nor ASSERT_PROD level of asserts is active.  This is the
 * default).
 *
 * This macro invokes the standard assert() function with the specified
 * assert test.
 */
#define SPI_assert_debug(x) assert(x)

#endif

#endif

#endif /* __CNK__ */


__END_DECLS


#endif
