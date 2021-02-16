#include "TList.h"
#include "TIterator.h"
#include "TObject.h"
#include "NueAna/Extrapolation/NueFNHelper.h"
#include <iostream>

using namespace std;

NueFNHelper::NueFNHelper(Int_t nx,Double_t lx,Double_t ux,
			 Int_t ny,Double_t ly,Double_t uy) :
  NueExtrapHelper(nx,lx,ux,ny,ly,uy)
{
}

NueFNHelper::NueFNHelper(Int_t nx,Double_t *x,
			 Int_t ny,Double_t *y) :
  NueExtrapHelper(nx,x,ny,y)
{
}

NueFNHelper::~NueFNHelper()
{
}

void NueFNHelper::AddNueSystematic(NueSystematic *nueSys)
{
  fFarNearEnergyHists[nueSys];
  Int_t max_bg_index = 0;
  while(strcmp(Background::
	       AsString(Background::EBackground(max_bg_index)),
	       "?Unknown?")!=0) {
    gDirectory->cd("/");
    string fnh_name = string(Background::
			     AsString(Background::EBackground(max_bg_index)));
    fnh_name += "_" + string(nueSys->GetName());
    FNHists *fnh = new FNHists(fnh_name.c_str());    
    fnh->fDirectory->cd();
    fnh->fND_RecoEnergy = new TH1D("ND_RecoEnergy","ND Reco Energy",fNXBins,fXBins);
    fnh->fND_RecoEnergy->Sumw2();
    fnh->fFD_RecoEnergy = new TH1D("FD_RecoEnergy","FD Reco Energy",fNXBins,fXBins);
    fnh->fFD_RecoEnergy->Sumw2();
    fnh->fND_TrueEnergy = new TH1D("ND_TrueEnergy","ND True Energy",fNXBins,fXBins);
    fnh->fND_TrueEnergy->Sumw2();
    fnh->fFD_TrueEnergy = new TH1D("FD_TrueEnergy","FD True Energy",fNXBins,fXBins);
    fnh->fFD_TrueEnergy->Sumw2();
    (fFarNearEnergyHists[nueSys])[Background::EBackground(max_bg_index)] = fnh;
    max_bg_index++;
  }
  gDirectory->cd("/");
}

void NueFNHelper::MakeHelpers(Selection::Selection_t sel) 
{  
  if(!fFarChain && !fNearChain) return;
  fCurSel = sel;
  std::map<NueSystematic*,
    std::map<Background::Background_t,FNHists*> >::iterator mapBeg = 
    fFarNearEnergyHists.begin();
  std::map<NueSystematic*,
    std::map<Background::Background_t,FNHists*> >::iterator mapEnd = 
    fFarNearEnergyHists.end();
  //near:
  Int_t nEntries = fNearChain->GetEntries();
  for(int i=0;i<nEntries;i++){    
    if(i%10000==0) std::cout << "Processed " 
			     << 100*Float_t(i)/Float_t(nEntries) 
			     << "% of Near" << std::endl;
    mapBeg = fFarNearEnergyHists.begin();
    while(mapBeg!=mapEnd) {
      fNearChain->GetEntry(i);
      Double_t totWeight = 1;
      if((mapBeg->first)) totWeight = (mapBeg->first)->UpdateRecord(fRecord,sel);
      if(this->PassCuts(sel)) {	
	Double_t recoEnergy = this->GetNueEnergy(sel);
	Background::Background_t bg = 
	  Background::TranslateFromMC(fRecord->mctrue.interactionType,
				      fRecord->mctrue.nuFlavor,
				      fRecord->mctrue.nonOscNuFlavor,
				      fRecord->fluxinfo.tptype);
	(mapBeg->second)[bg]->fND_RecoEnergy->Fill(recoEnergy,totWeight); 
	(mapBeg->second)[bg]->fND_TrueEnergy->Fill(fRecord->mctrue.nuEnergy,totWeight); 
	Background::Background_t bg2 = 
	  Background::TranslateFromMC(fRecord->mctrue.interactionType,
				      fRecord->mctrue.nuFlavor,
				      fRecord->mctrue.nonOscNuFlavor);
	if(bg2!=bg) {
	  (mapBeg->second)[bg2]->fND_RecoEnergy->Fill(recoEnergy,totWeight); 
	  (mapBeg->second)[bg2]->fND_TrueEnergy->Fill(fRecord->mctrue.nuEnergy,totWeight);
	}
      }
      mapBeg++;
    }
  }

  // far:
  nEntries = fFarChain->GetEntries();
  for(int i=0;i<nEntries;i++){
    if(i%10000==0) std::cout << "Processed " 
			     << 100*Float_t(i)/Float_t(nEntries) 
			     << "% of Far" << std::endl;
    mapBeg = fFarNearEnergyHists.begin();
    while(mapBeg!=mapEnd) {
      fFarChain->GetEntry(i);
      Double_t totWeight = 1;
      if((mapBeg->first)) totWeight = (mapBeg->first)->UpdateRecord(fRecord,sel);
      if(this->PassCuts(sel)) {	
	Double_t recoEnergy = this->GetNueEnergy(sel);
	if(TMath::Abs(fRecord->fluxinfo.ptype)==13 &&
	   fRecord->fluxinfo.tptype==0) fRecord->fluxinfo.tptype = 211;
	Background::Background_t bg = 
	  Background::TranslateFromMC(fRecord->mctrue.interactionType,
				      fRecord->mctrue.nuFlavor,
				      fRecord->mctrue.nonOscNuFlavor,
				      fRecord->fluxinfo.tptype);
	(mapBeg->second)[bg]->fFD_RecoEnergy->Fill(recoEnergy,totWeight); 
	(mapBeg->second)[bg]->fFD_TrueEnergy->Fill(fRecord->mctrue.nuEnergy,totWeight); 
	Background::Background_t bg2 = 
	  Background::TranslateFromMC(fRecord->mctrue.interactionType,
				      fRecord->mctrue.nuFlavor,
				      fRecord->mctrue.nonOscNuFlavor);
	if(bg2!=bg) {
	  (mapBeg->second)[bg2]->fFD_RecoEnergy->Fill(recoEnergy,totWeight);
	  (mapBeg->second)[bg2]->fFD_TrueEnergy->Fill(fRecord->mctrue.nuEnergy,totWeight);
	}
      }
      mapBeg++;
    }
  }
}

void NueFNHelper::WriteFile(std::string tag)
{
  std::map<NueSystematic*,
    std::map<Background::Background_t,FNHists*> >::iterator mapBeg = 
    fFarNearEnergyHists.begin();
  std::map<NueSystematic*,
    std::map<Background::Background_t,FNHists*> >::iterator mapEnd = 
    fFarNearEnergyHists.end();
  if(mapBeg==mapEnd) return;
  
  std::string filename = "EnergySpectraHelper_" + tag + ".root";
  TFile *file = new TFile(filename.c_str(),"RECREATE");
  file->cd();
  
  char selection[256]; 
  sprintf(selection,"%s",Selection::AsString(fCurSel));  
  TTree *tree = new TTree("energytree","energytree");
  tree->Branch("Selection",selection,"Selection/C");
  tree->Branch("nearPOT",&fNearPOT,"nearPOT/D");
  tree->Branch("farPOT",&fFarPOT,"farPOT/D");  
  
  while(mapBeg!=mapEnd){
    std::map<Background::Background_t,FNHists*>::iterator FNbeg = (mapBeg->second).begin();
    std::map<Background::Background_t,FNHists*>::iterator FNend = (mapBeg->second).end();
    while(FNbeg!=FNend) {
      TDirectory *filedir = file->mkdir(FNbeg->second->fDirectory->GetName());
      filedir->cd();
      TList *list = FNbeg->second->fDirectory->GetList();
      TIter iter(list->MakeIterator());
      TObject *ob = 0;
      while((ob = iter())) ob->Write();
      file->cd();
      FNbeg++;
    }
    if((mapBeg->first)) (mapBeg->first)->MakeBranches(tree);
    tree->Fill();
    mapBeg++;
  }

  tree->Write();
  delete file;
}
