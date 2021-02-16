#ifndef NUEBACKGROUND_H
#define NUEBACKGROUND_H
#include <string>
#include <ostream>
#include "TH1.h"
#include "Conventions/Detector.h"
#include "NueAna/NueAnaTools/Selection.h"
#include "NueAna/Extrapolation/Background.h"

class NueBackground
{

 public:

  NueBackground();
  NueBackground(std::string,TH1*,Detector::Detector_t,
		Background::Background_t,Selection::Selection_t,Double_t);
  ~NueBackground() { delete fHist;}
  
  const TH1 *GetHist() {return fHist;}
  const char *GetName() {return fName.c_str();}
  Detector::Detector_t GetDetector() {return fDet;}
  Background::Background_t GetBackground() {return fBg;}
  Selection::Selection_t GetSelection() {return fSel;}
  Double_t GetPOT() {return fPOT;}
  void Print(std::ostream &) const;
  void Draw(Float_t pot=0,Int_t col=1,const char *opt="");

 private:
  
  std::string fName;
  TH1 *fHist; //note this can be a pointer to TH1F/TH1D/TH2F/etc.
  Detector::Detector_t fDet;
  Background::Background_t fBg;
  Selection::Selection_t fSel;
  Double_t fPOT;

};
#endif //NUEBACKGROUND_H
