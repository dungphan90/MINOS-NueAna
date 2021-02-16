#include "NueAna/Extrapolation/NueExtrapolation.h"
#include "TTree.h"
#include <stdlib.h>
#include <string>
#include <iostream>

using namespace std;

NueExtrapolation::NueExtrapolation() :
  fBg(0),
  fHelperFile(0),
  fExtrapMethod(Extrapolation::kNone)
{
}

NueExtrapolation::~NueExtrapolation()
{
  if(fHelperFile) fHelperFile->Close();
}

void NueExtrapolation::SetHelperFile(std::string filename) 
{
  fHelperFile = new TFile(filename.c_str(),"READ");
}

TH1* NueExtrapolation::GetSpectrum(Detector::Detector_t det,
				   Double_t pot,
				   NueBackground* bg)
{
  this->SetNueBackground(bg);
  TH1 *hist = NULL;
  char histName[256];
  sprintf(histName,"%s_Prediction_%s_%s_%s_%.2fe20POT",
	  Detector::AsString(det),
	  Selection::AsString(fBg->GetSelection()),
	  Background::AsString(fBg->GetBackground()),
	  Extrapolation::AsString(fExtrapMethod),pot*1e-20);

  //first check whether this is already the right detector:
  if(det == fBg->GetDetector()) hist = this->None(histName);
  else {
    //are we going near->far or far->near?
    Bool_t extrapNtoF = true;
    if(det == Detector::kNear) extrapNtoF = false;
    hist = this->Extrapolate(histName,extrapNtoF); 
  }
  hist->Scale(pot/fBg->GetPOT());
  return hist;
}

TH1 *NueExtrapolation::None(const char *histName)
{
  return (TH1*) fBg->GetHist()->Clone(histName);
}

TH1* NueExtrapolation::Extrapolate(const char *histName,Bool_t /*NtoF*/)
{
  return this->None(histName);
}
