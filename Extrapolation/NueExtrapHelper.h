#ifndef NUEEXTRAPHELPER_H
#define NUEEXTRAPHELPER_H
#include "NueAna/NueAnaTools/Selection.h"
#include "NueAna/Extrapolation/Background.h"
#include "NueAna/Extrapolation/NueSystematic.h"
#include "NueAna/NueRecord.h"
#include "Conventions/Detector.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TChain.h"
#include "NueAna/NueMini.h"

class NueExtrapHelper
{

 public:

  NueExtrapHelper();
  NueExtrapHelper(Int_t nx,Double_t lx,Double_t ux,
		  Int_t ny=0,Double_t ly=0,Double_t uy=0);
  NueExtrapHelper(Int_t nx,  Double_t *xbins,
		  Int_t ny=0,Double_t *ybins=0);
  virtual ~NueExtrapHelper();
  
  void SetChains(TChain*,TChain*,Double_t,Double_t);

  virtual void AddNueSystematic(NueSystematic *nueSys);

  virtual void MakeHelpers(Selection::Selection_t);
  virtual void WriteFile(std::string);
  
 protected:
  
  void     Init();
  Bool_t   PassBasicCuts();
  Bool_t   PassCuts(NueRecord *nr, Selection::Selection_t sel);
  Bool_t   PassCuts(Selection::Selection_t);
  void     SetUpNueAnaChain(TChain*);
  void     SetUpNueMiniChain(TChain*);
  Double_t GetNueEnergy(NueRecord*, Selection::Selection_t);
  Double_t GetNueEnergy(Selection::Selection_t);

  Int_t     fNXBins;
  Double_t *fXBins;
  Int_t     fNYBins;
  Double_t *fYBins;
  
  TChain   *fNearChain;
  TChain   *fFarChain;
  Double_t  fNearPOT;
  Double_t  fFarPOT;
  
  Selection::Selection_t fCurSel;

  //Relevant variables from AnaNue files
  NueRecord *fRecord;
  NueMini   *fMini;
};
#endif //NUEEXTRAPHELPER_H
