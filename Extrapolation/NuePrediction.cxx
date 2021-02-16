#include "NueAna/Extrapolation/NuePrediction.h"
#include "TLegend.h"
#include "TH2.h"
#include <iostream>

NuePrediction::NuePrediction()
{
  fFinalPred = NULL;
}

NuePrediction::~NuePrediction()
{
  this->Clear();
  fExtrapMap.clear();
}

void NuePrediction::ClearExtrap() 
{
  fExtrapMap.clear();
}

void NuePrediction::Clear() 
{
  fBgVec.clear();
  fExtrapVec.clear();
  this->ClearPred();
}

void NuePrediction::ClearPred()
{
  delete fFinalPred;
  fFinalPred = NULL;
  std::vector<TH1*>::iterator beg = fCurrentPred.begin();
  std::vector<TH1*>::iterator end = fCurrentPred.end();
  while(beg!=end) {
    delete (*beg);
    beg++;
  }
  fCurrentPred.clear();
}

Bool_t NuePrediction::AddExtrapolation(NueExtrapolation *ex,
				       Extrapolation::Extrapolation_t extrap)
{
  fExtrapMap[extrap] = ex;
  return true;
}

Bool_t NuePrediction::AddBackground(NueBackground *bg,
				    Extrapolation::Extrapolation_t extrap)
{
  fBgVec.push_back(bg);
  fExtrapVec.push_back(extrap);
  return true;
}

TH1D* NuePrediction::GetPrediction(Detector::Detector_t det,Double_t pot)
{
  if(fBgVec.begin()==fBgVec.end()) {
    return NULL;
    std::cerr << "No NueBackground objects added!" << std::endl;
  }
  std::vector<NueBackground*>::iterator bgBeg = fBgVec.begin();
  std::vector<NueBackground*>::iterator bgEnd = fBgVec.end();
  std::vector<Extrapolation::Extrapolation_t>::iterator extrapBeg = fExtrapVec.begin();
  std::vector<Extrapolation::Extrapolation_t>::iterator extrapEnd = fExtrapVec.end();
  
  this->ClearPred();
  while(bgBeg!=bgEnd && extrapBeg!=extrapEnd){
    fCurrentPred.push_back(fExtrapMap[(*extrapBeg)]->GetSpectrum(det,pot,(*bgBeg)) );
    bgBeg++;
    extrapBeg++;
  }
  char name[256];
  sprintf(name,"%s_Prediction_%s_%.2f#times10^{20}POT",
	  Detector::AsString(Detector::EDetector(det)),
	  Selection::AsString(Selection::ESelection(fBgVec[0]->GetSelection())),
	  pot*1e-20);
  fFinalPred = new TH1D(name,name,
			fCurrentPred[0]->GetNbinsX(),
			fCurrentPred[0]->GetBinLowEdge(1),
			fCurrentPred[0]->GetBinLowEdge(fCurrentPred[0]->
						       GetNbinsX()+1));
  fFinalPred->Sumw2();

  std::vector<TH1*>::iterator beg = fCurrentPred.begin();
  std::vector<TH1*>::iterator end = fCurrentPred.end();
  while(beg!=end) {
    //determine what type of histogram this is:
    if((*beg)->InheritsFrom("TH2")) {
      TH2 *hist = (TH2*) (*beg);
      TH1D *proj = hist->ProjectionX();
      fFinalPred->Add(proj);
    }
    else {
      TH1D *spec = (TH1D*) (*beg);
      fFinalPred->Add(spec);
    }
    beg++;
  }
  return fFinalPred;
}

void NuePrediction::Draw()
{
  if(!fFinalPred) return;

  fFinalPred->Draw("e");
  Int_t col = 2;

  TLegend *leg = new TLegend(0.2,0.2,0.4,0.4);
  leg->AddEntry(fFinalPred,"Total Prediction","lp");
  char legname[256];

  std::vector<TH1*>::iterator beg = fCurrentPred.begin();
  std::vector<TH1*>::iterator end = fCurrentPred.end();
  std::vector<NueBackground*>::iterator bgBeg = fBgVec.begin();
  std::vector<NueBackground*>::iterator bgEnd = fBgVec.end();
  std::vector<Extrapolation::Extrapolation_t>::iterator extrapBeg = fExtrapVec.begin();
  std::vector<Extrapolation::Extrapolation_t>::iterator extrapEnd = fExtrapVec.end();

  while(beg!=end && bgBeg!=bgEnd && extrapBeg!=extrapEnd) {
    if(col==5) col+=1;
    sprintf(legname,"%s_%s",
	    Background::AsString(Background::
				 EBackground((*bgBeg)->GetBackground())),
	    Extrapolation::AsString(Extrapolation::
				    EExtrapolation((*extrapBeg))));
    if((*beg)->InheritsFrom("TH2")) {
      TH2* hist = (TH2*) (*beg);
      TH1D *proj = hist->ProjectionX();
      proj->SetLineColor(col);
      proj->Draw("histsames");
      leg->AddEntry(proj,legname,"l");
    }
    else {
      (*beg)->SetLineColor(col);
      (*beg)->Draw("histsames");
      leg->AddEntry((*beg),legname,"l");
    }

    beg++;
    bgBeg++;
    extrapBeg++;
    col++;
  }
  leg->SetBorderSize(1);
  leg->Draw();
}
