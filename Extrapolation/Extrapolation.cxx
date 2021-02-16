// Body for Extrapolation namespace so that CINT recognizes its existence

#include "NueAna/Extrapolation/Extrapolation.h"
#include "TString.h"

//_____________________________________________________________________________
const Char_t* Extrapolation::AsString(Extrapolation_t extrapolation)
{
  switch (extrapolation) {
  case kNone:      return "None";      break;
  case kNorm:      return "Norm";      break;
  case kFN:        return "F/N";       break;
  case kMM:        return "Matrix";    break;
  case kFDCC:      return "FDCCFlux";  break;
  case kNDFit:     return "NDFit";     break;
  case kUnknown:   return "Unknown";   break;
  default:         return "?Unknown?"; break;
  }
}

//_____________________________________________________________________________
Extrapolation::Extrapolation_t Extrapolation::StringToEnum(const Char_t* chars)
{
  TString theString(chars);
  if(theString.Contains("None"))     return kNone;
  if(theString.Contains("Norm"))     return kNorm;
  if(theString.Contains("F/N"))      return kFN;
  if(theString.Contains("Matrix"))   return kMM;
  if(theString.Contains("FDCCFlux")) return kFDCC;
  if(theString.Contains("NDFit"))    return kNDFit;
  
  return kUnknown;
}
