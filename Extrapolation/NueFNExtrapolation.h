#ifndef NUEFNEXTRAPOLATION_H
#define NUEFNEXTRAPOLATION_H
#include "NueAna/Extrapolation/NueBackground.h"
#include "NueAna/Extrapolation/Extrapolation.h"
#include "NueAna/Extrapolation/NueExtrapolation.h"
#include "NueAna/Extrapolation/NueSystematic.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "Conventions/Detector.h"

class NueFNExtrapolation : public NueExtrapolation
{
  
 public:
  NueFNExtrapolation();
  NueFNExtrapolation(std::string);
  virtual ~NueFNExtrapolation();

  void UseTreeEntry(Int_t treeEntry=0) {fTreeEntry = treeEntry;}
  Int_t GetMaxTreeEntry() {if(fHelperTree) return fHelperTree->GetEntries(); return 0;}
  TH1 *GetRatio(string);
  NueSystematic *GetCurrentSystematic();

 protected:
  
  TH1 *Extrapolate(const char*,Bool_t);

 private:

  TH1D  *fHelperRatio;
  TTree *fHelperTree;
  Int_t fTreeEntry;

};
#endif //NUEFNEXTRAPOLATION_H
