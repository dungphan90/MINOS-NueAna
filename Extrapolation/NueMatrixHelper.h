#ifndef NUEMATRIXHELPER_H
#define NUEMATRIXHELPER_H
#include <string>
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "NueAna/Extrapolation/MatrixHists.h"
#include "NueAna/Extrapolation/NueSystematic.h"
#include "NueAna/Extrapolation/Background.h"
#include "NueAna/Extrapolation/NueExtrapHelper.h"

class NueMatrixHelper : public NueExtrapHelper
{

 public:
  NueMatrixHelper(Int_t,Double_t,Double_t);
  NueMatrixHelper(Int_t,Double_t*);
  virtual ~NueMatrixHelper();
  
  void MakeANANUEPlots(Selection::Selection_t);
  //void MakeSNTPPlots(TChain*,TChain*);
  //void MakeFLUXPlots(TChain*,TGraph*,TH1F *mikehist=0);
  //void ReadXSecPlots(TFile*);
  void WriteFile(std::string);
  
  virtual void AddNueSystematic(NueSystematic*);

 private:

  std::map<NueSystematic*,std::map<Background::Background_t,MatrixHists*> > fMatrixHists;

};
#endif //NUEMATRIXHELPER_H
