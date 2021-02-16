// Body for Systematic namespace so that CINT recognizes its existence

#include "NueAna/Extrapolation/Systematic.h"
#include "TString.h"

//_____________________________________________________________________________
const Char_t* Systematic::AsString(Systematic_t syst)
{
   switch (syst) {
   case kNorm     : return "Norm";      break;
   case kEMCalib  : return "EMCalib";   break;
   case kHadCalib : return "HadCalib";  break;
   case kRelCalib : return "RelCalib";  break;
   case kMA_QE    : return "MA_QE";     break;
   case kMA_RES   : return "MA_RES";    break;
   case kKNO      : return "KNO";       break;
   case kTrkPlane : return "TrkPlane";  break;
   case kPIDShift : return "PIDShift";  break;
   case kSKZP     : return "SKZP";      break;
   case kOscProb  : return "OscProb";   break;
   case kShwDev   : return "ShwDev";    break;
   case kTauProd  : return "TauProd";   break;
   case kPIDSkew  : return "PIDSkew";   break;
   case kTrkLike  : return "TrkLike";   break;
   case kNCScale  : return "NCScale";   break;
   case kCCShwE   : return "CCShwEnergy"; break;
   case kUnknown  : return "Unknown";   break;
   default        : return "?Unknown?"; break;
   }
}

//_____________________________________________________________________________
Systematic::Systematic_t Systematic::StringToEnum(const Char_t* chars)
{
  TString theString(chars);
  if(theString.Contains("Norm"))      return kNorm;
  if(theString.Contains("EMCalib"))   return kEMCalib;
  if(theString.Contains("HadCalib"))  return kHadCalib;  
  if(theString.Contains("RelCalib"))  return kRelCalib;
  if(theString.Contains("MA_QE"))     return kMA_QE;
  if(theString.Contains("MA_RES"))    return kMA_RES;
  if(theString.Contains("KNO"))       return kKNO;
  if(theString.Contains("TrkPlane"))  return kTrkPlane;
  if(theString.Contains("TrkLike"))  return kTrkLike;
  if(theString.Contains("PIDShift"))  return kPIDShift;
  if(theString.Contains("SKZP"))      return kSKZP;
  if(theString.Contains("OscProb"))   return kOscProb;
  if(theString.Contains("ShwDev"))    return kShwDev;
  if(theString.Contains("TauProd"))   return kTauProd;
  if(theString.Contains("PIDSkew"))   return kPIDSkew;
  if(theString.Contains("CCShwEnergy")) return kCCShwE;
  if(theString.Contains("NCScale"))   return kNCScale;
  
  return kUnknown;
}

Double_t Systematic::GetDefaultValue(Systematic_t syst)
{
   switch (syst) {
   case kNorm     : return 1; break;
   case kEMCalib  : return 0; break;
   case kHadCalib : return 0; break;
   case kRelCalib : return 0; break;
   case kMA_QE    : return 1; break;
   case kMA_RES   : return 1; break;
   case kKNO      : return 1; break;
   case kTrkPlane : return 0; break;
   case kTrkLike  : return 0; break;
   case kPIDShift : return 0; break;
   case kSKZP     : return 0; break;
   case kOscProb  : return 1; break;
   case kShwDev   : return 0; break;
   case kTauProd  : return 0; break;
   case kPIDSkew  : return 0; break;
   case kNCScale  : return 1; break;
   case kCCShwE   : return 0; break;
   case kUnknown  : return 0; break;
   default        : return 0; break;
   }
}
