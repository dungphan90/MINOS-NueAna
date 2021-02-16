////////////////////////////////////////////////////////////////////////////
// $Id: Systematic.h,v 1.3 2008/03/06 00:14:42 boehm Exp $
//
// Systematic
//
// Systematic defines Systematic_t which is an enumeration of the
// possible nue systematic effects
//
// Author: C. Smith
//
////////////////////////////////////////////////////////////////////////////

#ifndef SYSTEMATIC_H
#define SYSTEMATIC_H

#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

namespace Systematic {

  typedef enum ESystematic {
    kNorm      = 0,
    kEMCalib   = 1,
    kHadCalib  = 2,
    kRelCalib  = 3,
    kMA_QE     = 4,
    kMA_RES    = 5,
    kKNO       = 6,
    kTrkPlane  = 7,
    kPIDShift  = 8,
    kSKZP      = 9,
    kOscProb   = 10,
    kShwDev    = 11,
    kTauProd   = 12,
    kPIDSkew   = 13,
    kTrkLike   = 14,
    kNCScale   = 15,
    kCCShwE    = 16,
    kUnknown   = 17
  } Systematic_t;

  // no ctor or dtor's - this class consists of only static members

  // Translation enum to/from character strings
  const Char_t*             AsString(Systematic_t syst);
  Systematic::Systematic_t  StringToEnum(const Char_t* chars);
  Double_t                  GetDefaultValue(Systematic_t syst);

}

#endif // SYSTEMATIC_H
