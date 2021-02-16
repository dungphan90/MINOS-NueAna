// Body for Background namespace so that CINT recognizes its existence

#include "NueAna/Extrapolation/Background.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "TString.h"

//_____________________________________________________________________________
const Char_t* Background::AsString(Background_t background)
{
   switch (background) {
   case kNueCC:     return "NueCC";      break;
   case kNC:        return "NC";         break;
   case kNuMuCC:    return "NuMuCC";     break;
   case kBNueCC:    return "BNueCC";     break;
   case kNuTauCC:   return "NuTauCC";    break;
   case kPiBNueCC:  return "PiBNueCC";   break;
   case kKaBNueCC:  return "KaBNueCC";   break;    
   case kSelCC:     return "SelCC";      break; 
   case kNuMuNC:    return "NuMuNC";     break;
   case kBNueNC:    return "BNueNC";     break;
   case kUnknown:   return "Unknown";    break;
   default:         return "?Unknown?";  break;
   }
}

//_____________________________________________________________________________
Background::Background_t Background::StringToEnum(const Char_t* chars)
{
  TString theString(chars);
  if(theString.Contains("NueCC"))    return kNueCC;
  if(theString.Contains("NC"))       return kNC;
  if(theString.Contains("NuMuCC"))   return kNuMuCC;
  if(theString.Contains("BNueCC"))   return kBNueCC;
  if(theString.Contains("NuTauCC"))  return kNuTauCC;
  if(theString.Contains("PiBNueCC")) return kPiBNueCC;
  if(theString.Contains("KaBNueCC")) return kKaBNueCC;
  if(theString.Contains("SelCC"))    return kSelCC;
  if(theString.Contains("NuMuNC"))   return kNuMuNC;
  if(theString.Contains("BNueNC"))   return kBNueNC;
  
  return kUnknown;
}

Background::Background_t Background::TranslateFromMC(Int_t iaction,Int_t inu,
						     Int_t inunoosc,Int_t parentid)
{
  if(iaction==0) return Background::kNC;
  if(iaction==0 && TMath::Abs(inu)==14 && TMath::Abs(inunoosc)==14){
     return Background::kNuMuNC;
  }
  if(iaction==0 && TMath::Abs(inu)==12 && TMath::Abs(inunoosc)==12){
     return Background::kBNueNC;
  }
  if(iaction==1) {
    if(TMath::Abs(inu)==12) {
      if(TMath::Abs(inunoosc)==12) {
	if(parentid==0) return Background::kBNueCC;
	else if(TMath::Abs(parentid)==211) return Background::kPiBNueCC;
	else if(TMath::Abs(parentid)==130 || 
		TMath::Abs(parentid)==321 ||
		TMath::Abs(parentid)==310) return Background::kKaBNueCC;
      }
      else if(TMath::Abs(inunoosc)==14) {
	return Background::kNueCC;
      }
    }
    else if(TMath::Abs(inu)==14) return Background::kNuMuCC;
    else if(TMath::Abs(inu)==16) return Background::kNuTauCC;
  }
  return Background::kUnknown;
}

Background::Background_t Background::TranslateFromNueClass(Int_t nueClass)
{
  if(nueClass == ClassType::NC)    return Background::kNC;
  if(nueClass == ClassType::numu)  return Background::kNuMuCC;
  if(nueClass == ClassType::nue)   return Background::kNueCC;
  if(nueClass == ClassType::nutau) return Background::kNuTauCC;
  if(nueClass == ClassType::bnue)  return Background::kBNueCC;

  return Background::kUnknown;
}

