#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of rfloat operator methods.

  $Id: rfloat_rnd_op.C,v 1.4 2004-08-18 11:58:08 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:08 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rfloat/noarch/rfloat_rnd_op.C,v 1.4 2004-08-18 11:58:08 zs Exp $
//  $Id: rfloat_rnd_op.C,v 1.4 2004-08-18 11:58:08 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rfloat/noarch/rfloat_rnd_op.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// rfloat_rnd_op.C
//
// These overloaded operators perform the actual rounding.
// This file should not be compiled for use on the DSP.
// This file is provided for consistency when compiling
// for a machine that does rounding. If one needs to compile
// for the DSP the rfload_rnd_op.asm file should be used
// instead. The code in the .asm file does the actual rounding.
//
// More overloaded operators have been directly defined
// in rfloat.h but they all use the operators below.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/rfloat.h>
CPS_START_NAMESPACE


rfloat& rfloat::operator+=(IFloat a)
{x += a;  return *this;}

rfloat& rfloat::operator-=(IFloat a)
{x -= a;  return *this;}

rfloat& rfloat::operator*=(IFloat a)
{x *= a;  return *this;}

rfloat& rfloat::operator/=(IFloat a)
{x /= a;  return *this;}

rfloat operator-(const rfloat& a)
{return -a.x;}


CPS_END_NAMESPACE
