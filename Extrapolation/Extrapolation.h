////////////////////////////////////////////////////////////////////////////
// $Id: Extrapolation.h,v 1.1 2006/11/11 16:07:58 cbs Exp $
//
// Extrapolation
//
// Extrapolation defines Extrapolation_t which is an enumeration of the
// possible nue selection criteria
//
// Author: C. Smith
//
////////////////////////////////////////////////////////////////////////////

#ifndef EXTRAPOLATION_T
#define EXTRAPOLATION_T

#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

namespace Extrapolation {

  typedef enum EExtrapolation {
    kNone     = 0,
    kNorm     = 1,
    kFN       = 2,
    kMM       = 3,
    kFDCC     = 4, 
    kNDFit    = 5,
    kUnknown  = 6
  } Extrapolation_t;

  // no ctor or dtor's - this class consists of only static members

  // Translation enum to/from character strings
  const Char_t*           AsString(Extrapolation_t selection);
  Extrapolation::Extrapolation_t  StringToEnum(const Char_t* chars);
  
}

#endif // EXTRAPOLATION_H
