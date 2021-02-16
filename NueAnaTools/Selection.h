////////////////////////////////////////////////////////////////////////////
// $Id: Selection.h,v 1.15 2017/02/27 18:13:05 wingmc Exp $
//
// Selection
//
// Selection defines Selection_t which is an enumeration of the
// possible nue selection criteria
//
// Author: C. Smith
//
////////////////////////////////////////////////////////////////////////////

#ifndef SELECTION_H
#define SELECTION_H

#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

namespace Selection {

  typedef enum ESelection {
      kNone     = 0,
      kDataQual = 1,
      kFid      = 2,
      kBasic    = 3,
      kPre      = 4,
      kCuts     = 5,
      kANN6     = 6,
      kANN30    = 7,
      kSSPID    = 8, 
      kMCNN     = 9,
      kBDT      = 10,
      kKNue     = 11, 
      kMDA      = 12,
      kANN2PE   = 13,
      kANN2PE_DAIKON04 = 14,
      kANN14_DAIKON04 = 15,
      kCC       = 16,
      kParticlePID = 17,
      kLEMNNBAR = 18,
      kLEMBAR = 19,
      kLEM4 = 20,
      kANN4FHC = 21,
      kANN4RHC = 22,
      kLEMAmby = 23,
      kLEMLSND = 24,
      kLEMAmbyE50N491 = 25,
      kLEMAmbyE50N591 = 26,
      kLEMAmbyE50N691 = 27,
      kLEMAmbyE50N5111 = 28,
      kLEMAmbyE50N6111 = 29,
      kLEMAmbyE50S491 = 30,
      kLEMAmbyE50S591 = 31,
      kLEMAmbyE50S691 = 32,
      kLEMAmbyE50S5111 = 33,
      kLEMAmbyE50S6111 = 34,
      kLEMAmbySAE50S491 = 35,     
      kUnknown  = 36
  } Selection_t;

  // no ctor or dtor's - this class consists of only static members

  // Translation enum to/from character strings
  const Char_t*           AsString(Selection_t selection);
  Selection::Selection_t  StringToEnum(const Char_t* chars);

}

#endif // SELECTION_H
