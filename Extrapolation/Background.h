////////////////////////////////////////////////////////////////////////////
// $Id: Background.h,v 1.3 2014/08/27 17:57:57 wingmc Exp $
//
// Background
//
// Background defines Background_t which is an enumeration of the
// possible nue signal and background
//
// Author: C. Smith
//
////////////////////////////////////////////////////////////////////////////

#ifndef BACKGROUND_H
#define BACKGROUND_H

#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

namespace Background {

  typedef enum EBackground {
    kNueCC    = 0,
    kNC       = 1,
    kNuMuCC   = 2,
    kBNueCC   = 3, 
    kNuTauCC  = 4,
    kPiBNueCC = 5,
    kKaBNueCC = 6,
    kSelCC    = 7,
    kNuMuNC   = 8,
    kBNueNC   = 9,
    kUnknown  = 10
  } Background_t;

  // no ctor or dtor's - this class consists of only static members

  // Translation enum to/from character strings
  const Char_t*             AsString(Background_t background);
  Background::Background_t  StringToEnum(const Char_t* chars);
  Background::Background_t  TranslateFromMC(Int_t iaction,Int_t inu,
					    Int_t inunoosc,Int_t parentid=0);
  Background::Background_t  TranslateFromNueClass(int nueClass);

}

#endif // BACKGROUND_H
