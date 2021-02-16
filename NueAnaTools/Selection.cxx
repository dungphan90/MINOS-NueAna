// Body for Selection namespace so that CINT recognizes its existence

#include "NueAna/NueAnaTools/Selection.h"
#include "TString.h"

//_____________________________________________________________________________
const Char_t* Selection::AsString(Selection_t selection)
{
   switch (selection) {
   case kNone:      return "None";       break;
   case kDataQual:  return "DataQual";   break;
   case kFid:       return "Fid";        break;
   case kBasic:     return "Basic";      break;
   case kPre:       return "Presel";     break;
   case kCuts:      return "Cuts";       break;
   case kANN6:      return "ANN6";       break;
   case kANN30:     return "ANN30";      break;
   case kSSPID:     return "SSPID";      break;
   case kMDA:       return "MDA";        break;
   case kBDT:       return "BDT";        break;
   case kKNue:      return "KNue";       break;
   case kMCNN:      return "MCNN";       break;
   case kLEMNNBAR:  return "LEMNNBAR";   break;
   case kLEMBAR:    return "LEMBAR";     break;
   case kLEM4:      return "LEM4";       break;
   case kLEMAmby:   return "LEMAmby";    break;
   case kLEMLSND:   return "LEMLSND";    break;
   case kANN2PE:    return "ANN2PE";     break;
   case kANN2PE_DAIKON04:    return "ANN2PE_DAIKON04";     break;
   case kANN14_DAIKON04:    return "ANN14_DAIKON04";     break;
   case kANN4FHC:   return "ANN4FHC";    break;
   case kANN4RHC:   return "ANN4RHC";    break;
   case kCC:        return "CC";         break;
   case kUnknown:   return "Unknown";    break;
   case kParticlePID: return "ParticlePID"; break;

   case kLEMAmbyE50N491:   return "LEMAmbyE50N491";    break;
   case kLEMAmbyE50N591:   return "LEMAmbyE50N591";    break;
   case kLEMAmbyE50N691:   return "LEMAmbyE50N691";    break;
   case kLEMAmbyE50N5111:   return "LEMAmbyE50N5111";    break;
   case kLEMAmbyE50N6111:   return "LEMAmbyE50N6111";    break;
   case kLEMAmbyE50S491:   return "LEMAmbyE50S491";    break;
   case kLEMAmbyE50S591:   return "LEMAmbyE50S591";    break;
   case kLEMAmbyE50S691:   return "LEMAmbyE50S691";    break;
   case kLEMAmbyE50S5111:   return "LEMAmbyE50S5111";    break;
   case kLEMAmbyE50S6111:   return "LEMAmbyE50S6111";    break;
   case kLEMAmbySAE50S491:   return "LEMAmbyE50S491";    break;

   default:         return "?Unknown?";  break;
   }
}

//_____________________________________________________________________________
Selection::Selection_t Selection::StringToEnum(const Char_t* chars)
{
  TString theString(chars);
  if(theString.Contains("None"))  return kNone;
  if(theString.Contains("DataQual")) return kDataQual;
  if(theString.Contains("Basic")) return kBasic;
  if(theString.Contains("Fid"))   return kFid;
  if(theString.Contains("Presel"))  return kPre;  
  if(theString.Contains("Cuts"))  return kCuts;
  if(theString.Contains("ANN6"))   return kANN6;
  if(theString.Contains("ANN30"))   return kANN30;
  if(theString.Contains("SSPID")) return kSSPID;
  if(theString.Contains("MDA"))   return kMDA;
  if(theString.Contains("BDT"))   return kBDT;
  if(theString.Contains("KNue"))  return kKNue;
  if(theString.Contains("MCNN"))  return kMCNN;
  if(theString.Contains("LEMNNBAR")) return kLEMNNBAR;
  if(theString.Contains("LEMBAR")) return kLEMBAR;
  if(theString.Contains("LEM4")) return kLEM4;
  if(theString.Contains("LEMAmby")) return kLEMAmby;
  if(theString.Contains("LEMLSND")) return kLEMLSND;
  if(theString.Contains("ANN2PE"))   return kANN2PE;
  if(theString.Contains("ANN2PE_DAIKON04"))   return kANN2PE_DAIKON04;
  if(theString.Contains("ANN14_DAIKON04"))   return kANN14_DAIKON04;
  if(theString.Contains("ANN4FHC"))   return kANN4FHC;
  if(theString.Contains("ANN4RHC"))   return kANN4RHC;
  if(theString.Contains("ParticlePID"))  return kParticlePID;
  if(theString.Contains("CC"))    return kCC;

  if(theString.Contains("LEMAmbyE50N491")) return kLEMAmbyE50N491;
  if(theString.Contains("LEMAmbyE50N591")) return kLEMAmbyE50N591;
  if(theString.Contains("LEMAmbyE50N691")) return kLEMAmbyE50N691;
  if(theString.Contains("LEMAmbyE50N5111")) return kLEMAmbyE50N5111;
  if(theString.Contains("LEMAmbyE50N6111")) return kLEMAmbyE50N6111;
  if(theString.Contains("LEMAmbyE50S491")) return kLEMAmbyE50S491;
  if(theString.Contains("LEMAmbyE50S591")) return kLEMAmbyE50S591;
  if(theString.Contains("LEMAmbyE50S691")) return kLEMAmbyE50S691;
  if(theString.Contains("LEMAmbyE50S5111")) return kLEMAmbyE50S5111;
  if(theString.Contains("LEMAmbyE50S6111")) return kLEMAmbyE50S6111;
  if(theString.Contains("LEMAmbySAE50S491")) return kLEMAmbySAE50S491;
  
  return kUnknown;
}
