#include "NueBackground.h"

NueBackground::NueBackground()
{
  fName = "";
  fHist = NULL;
  fDet = Detector::kUnknown;
  fBg = Background::kUnknown;
  fSel = Selection::kUnknown;
  fPOT = 0;  
}

NueBackground::NueBackground(std::string name,TH1* hist,
			     Detector::Detector_t det,
			     Background::Background_t bg,
			     Selection::Selection_t sel, Double_t pot)
{
  fName = name;
  fHist = (TH1*) hist->Clone(name.c_str());
  fHist->SetDirectory(0);
  fDet = det;
  fBg = bg;
  fSel = sel;
  fPOT = pot;
}

void NueBackground::Print(std::ostream &out) const
{
  out << fName << std::endl;
  out << Background::AsString(Background::EBackground(fBg)) << std::endl; 
  out << Detector::AsString(Detector::EDetector(fDet)) << std::endl;
  out << Selection::AsString(Selection::ESelection(fSel)) << std::endl;
  out << fPOT << " POT" << std::endl;
  out << "Histogram integral: " << fHist->Integral() << std::endl;
}

void NueBackground::Draw(Float_t pot,Int_t col,const char *opt)
{
  fHist->SetLineColor(col);
  if(pot>0) {
    fHist->Scale(pot/fPOT);
    fPOT = pot;
  }
  fHist->Draw(opt);  
}
