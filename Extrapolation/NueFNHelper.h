#ifndef NUEFNHELPER_H
#define NUEFNHELPER_H
#include "NueAna/NueAnaTools/Selection.h"
#include "NueAna/Extrapolation/Background.h"
#include "NueAna/Extrapolation/NueSystematic.h"
#include "NueAna/Extrapolation/NueExtrapHelper.h"
#include "NueAna/Extrapolation/FNHists.h"
#include "NueAna/NueRecord.h"
#include "Conventions/Detector.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TChain.h"

class NueFNHelper : public NueExtrapHelper
{

 public:

  NueFNHelper(Int_t nx,Double_t lx,Double_t ux,
	      Int_t ny=0,Double_t ly=0,Double_t uy=0);
  NueFNHelper(Int_t nx,  Double_t *xbins,
	      Int_t ny=0,Double_t *ybins=0);
  virtual ~NueFNHelper();

  void MakeHelpers(Selection::Selection_t);
  void WriteFile(std::string tag);
  
  void AddNueSystematic(NueSystematic *nueSys);

 private:
  
  std::map<NueSystematic*,std::map<Background::Background_t,FNHists*> > fFarNearEnergyHists;

};
#endif //NUEFNHELPER_H
