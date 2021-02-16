#ifndef NUEEXTRAPOLATION_H
#define NUEEXTRAPOLATION_H
#include "NueAna/Extrapolation/NueBackground.h"
#include "NueAna/Extrapolation/Extrapolation.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "Conventions/Detector.h"

class NueExtrapolation
{

 public:
  NueExtrapolation();
  virtual ~NueExtrapolation();
  
  TH1* GetSpectrum(Detector::Detector_t,Double_t, NueBackground*);  
  void Clear() {fBg = NULL;}
  
  Extrapolation::Extrapolation_t GetExtrapMethod() {return fExtrapMethod;}

  void SetNueBackground(NueBackground *bg) { fBg = bg; }

 protected:
  
  NueBackground *fBg;
  
  void SetHelperFile(std::string);
  TFile *fHelperFile;

  TH1 *None(const char*);
  Extrapolation::Extrapolation_t fExtrapMethod;
  
  virtual TH1 *Extrapolate(const char*,Bool_t);

};
#endif //NUEEXTRAPOLATION_H
