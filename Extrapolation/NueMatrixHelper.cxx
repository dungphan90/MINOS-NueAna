#include "NueMatrixHelper.h"
#include "TCut.h"
#include "TGraph.h"
#include "TRandom.h"
#include <vector>
#include <fstream>
#include <iostream>
#include "StandardNtuple/NtpStRecord.h"
#include "MCNtuple/NtpMCTruth.h"
#include "MCNtuple/NtpMCStdHep.h"
#include "MCReweight/NuParent.h"
#include "MCReweight/Zfluk.h"
#include "MCReweight/Zbeam.h"
#include "Registry/Registry.h"
#include "MCReweight/MCEventInfo.h"
#include "MCReweight/MCReweight.h"
#include "MCReweight/NeugenWeightCalculator.h"
//#include "reweight.C"
//#include "mcinfo.C"

using namespace std;

NueMatrixHelper::NueMatrixHelper(Int_t nbins,Double_t lowx,Double_t uppx) :
  NueExtrapHelper(nbins,lowx,uppx)
{
}

NueMatrixHelper::NueMatrixHelper(Int_t nbins,Double_t *bins) :
  NueExtrapHelper(nbins,bins)
{
}

void NueMatrixHelper::AddNueSystematic(NueSystematic *nueSys)
{
  fMatrixHists[nueSys];
  Int_t max_bg_index = 0;
  while(strcmp(Background::
	       AsString(Background::EBackground(max_bg_index)),
	       "?Unknown?")!=0) {
    gDirectory->cd("/");
    string mh_name = string(Background::
			    AsString(Background::EBackground(max_bg_index)));
    mh_name += "_" + string(nueSys->GetName());
    MatrixHists *mh = new MatrixHists(mh_name.c_str());
    mh->fDirectory->cd();
    
    //helper histos:
    mh->fRecoVsTrueEnergy_ND = 
      new TH2D("RecoVsTrueEnergy_ND",
	       "Reco Vs True Energy (NearDet)",
	       fNXBins,fXBins,fNXBins,fXBins);
    mh->fRecoVsTrueEnergy_ND->Sumw2();

    mh->fRecoVsTrueEnergy_FD = 
      new TH2D("RecoVsTrueEnergy_FD",
	       "Reco Vs True Energy (FarDet)",
	       fNXBins,fXBins,fNXBins,fXBins);
    mh->fRecoVsTrueEnergy_FD->Sumw2();

    mh->fEfficiency_ND = 
      new TH1D("Efficiency_ND",
	       "NuMu CC Selection Efficiency with True Energy (NearDet)",
	       fNXBins,fXBins);
    mh->fEfficiency_ND->Sumw2();

    mh->fEfficiency_FD = 
      new TH1D("Efficiency_FD",
	       "NuMu CC Selection Efficiency with True Energy (FarDet)",
	       fNXBins,fXBins);
    mh->fEfficiency_FD->Sumw2();

    mh->fPurity_ND = 
      new TH1D("Purity_ND",
	       "NuMu CC Selection Purity with True Energy (NearDet)",
	       fNXBins,fXBins);
    mh->fPurity_ND->Sumw2();

    mh->fPurity_FD = 
      new TH1D("Purity_FD",
	       "NuMu CC Selection Purity with True Energy (FarDet)",
	       fNXBins,fXBins);
    mh->fPurity_FD->Sumw2();
  
    mh->fFDVsNDMatrix = 
      new TH2D("FDVsNDMatrix",
	       "Number of FD Vs ND Events with True Energy",
	       fNXBins,fXBins,fNXBins,fXBins);
    mh->fFDVsNDMatrix->Sumw2();

    mh->fFDVsNDMatrixRW = 
      new TH2D("FDVsNDMatrixRW",
	       "Number of FD Vs ND Events with True Energy (with Near Reweight)",
	       fNXBins,fXBins,fNXBins,fXBins);
    mh->fFDVsNDMatrixRW->Sumw2();

    mh->fFDVsNDMatrixXSec = 
      new TH2D("FDVsNDMatrixXSec",
	       "Number of FD Vs ND Events with True Energy (with XSec)",
	       fNXBins,fXBins,fNXBins,fXBins);
    mh->fFDVsNDMatrixXSec->Sumw2();
    
    mh->fFDVsNDMatrixXSecRW = 
      new TH2D("FDVsNDMatrixXSecRW",
	       "Number of FD Vs ND Events with True Energy (with XSec + Near Reweight)",
	       fNXBins,fXBins,fNXBins,fXBins);
    mh->fFDVsNDMatrixXSecRW->Sumw2();
    
    mh->fXSec_CC = 
      new TH1D("XSec_CC","NuMu CC XSection with True Energy",
	       fNXBins,fXBins);
    mh->fXSec_CC->Sumw2();
    
    mh->fFracErrOnPred = 
      new TH1D("FracErrOnPred",
	       "Fractional Error on Energy Spectrum with Reco Energy",
	       fNXBins,fXBins);
    mh->fFracErrOnPred->Sumw2();
    
    //Check list histos:
    mh->fRecoEnergyAllEvents_ND = 
      new TH1D("RecoEnergyAllEvents_ND",
	       "Selected Events with Reco Energy (NearDet)",
	       fNXBins,fXBins);
    mh->fRecoEnergyAllEvents_ND->Sumw2();
    
    mh->fRecoEnergyCCOnlyEvents_ND = 
      new TH1D("RecoEnergyCCOnlyEvents_ND",
	       "NuMu CC Selected Events with Reco Energy (NearDet)",
	       fNXBins,fXBins);
    mh->fRecoEnergyCCOnlyEvents_ND->Sumw2();
    
    mh->fTrueEnergyCCOnlyEvents_ND = 
      new TH1D("TrueEnergyCCOnlyEvents_ND",
	       "NuMu CC Selected Events with True Energy (NearDet)",
	       fNXBins,fXBins);
    mh->fTrueEnergyCCOnlyEvents_ND->Sumw2();
    
    mh->fTrueEnergyTrueCCFidEvents_ND = 
      new TH1D("TrueEnergyTrueCCFidEvents_ND",
	       "True Fid NuMu CC Events with True Energy (NearDet)",
	       fNXBins,fXBins);
    mh->fTrueEnergyTrueCCFidEvents_ND->Sumw2();
    
    mh->fTrueEnergyNuFlux_ND = 
      new TH1D("TrueEnergyNuFlux_ND",
	       "Neutrino Flux with True Energy (NearDet)",
	       fNXBins,fXBins);
    mh->fTrueEnergyNuFlux_ND->Sumw2();
    
    mh->fTrueEnergyNuFluxRW_ND = 
      new TH1D("TrueEnergyNuFluxRW_ND",
	       "Neutrino Flux with True Energy (NearDet with Reweighting)",
	       fNXBins,fXBins);
    mh->fTrueEnergyNuFluxRW_ND->Sumw2();
    
    mh->fTrueEnergyNuFlux_FD = 
      new TH1D("TrueEnergyNuFlux_FD",
	       "Neutrino Flux with True Energy (FarDet)",
	       fNXBins,fXBins);
    mh->fTrueEnergyNuFlux_FD->Sumw2();
    
    mh->fTrueEnergyNuFluxRW_FD = 
      new TH1D("TrueEnergyNuFluxRW_FD",
	       "Neutrino Flux with True Energy (FarDet with Reweighting)",
	       fNXBins,fXBins);
    mh->fTrueEnergyNuFluxRW_FD->Sumw2();
    
    mh->fTrueEnergyCCFlux_ND = 
      new TH1D("TrueEnergyCCFlux_ND",
	       "NuMu CC Flux with True Energy (NearDet)",
	       fNXBins,fXBins);
    mh->fTrueEnergyCCFlux_ND->Sumw2();
    
    mh->fTrueEnergyCCFluxRW_ND = 
      new TH1D("TrueEnergyCCFluxRW_ND",
	       "NuMu CC Flux with True Energy with (NearDet with Reweighting)",
	       fNXBins,fXBins);
    mh->fTrueEnergyCCFluxRW_ND->Sumw2();
    
    mh->fTrueEnergyCCFlux_FD = 
      new TH1D("TrueEnergyCCFlux_FD",
	       "NuMu CC Flux with True Energy (FarDet)",
	       fNXBins,fXBins);
    mh->fTrueEnergyCCFlux_FD->Sumw2();
    
    mh->fTrueEnergyCCFluxRW_FD = 
      new TH1D("TrueEnergyCCFluxRW_FD",
	       "NuMu CC Flux with True Energy (FarDet with Reweighing)",
	       fNXBins,fXBins);
    mh->fTrueEnergyCCFluxRW_FD->Sumw2();
    
    mh->fTrueEnergyTrueCCFidEvents_FD = 
      new TH1D("TrueEnergyTrueCCFidEvents_FD",
	       "True Fid NuMu CC Events with True Energy (FarDet)",
	       fNXBins,fXBins);
    mh->fTrueEnergyTrueCCFidEvents_FD->Sumw2();
    
    mh->fTrueEnergyCCOnlyEvents_FD = 
      new TH1D("TrueEnergyCCOnlyEvents_FD",
	       "NuMu CC Selected Events with True Energy (FarDet)",
	       fNXBins,fXBins);
    mh->fTrueEnergyCCOnlyEvents_FD->Sumw2();
    
    mh->fRecoEnergyCCOnlyEvents_FD = 
      new TH1D("RecoEnergyCCOnlyEvents_FD",
	       "NuMu CC Selected Events with Reco Energy (FarDet)",
	       fNXBins,fXBins);
    mh->fRecoEnergyCCOnlyEvents_FD->Sumw2();
  
    mh->fRecoEnergyAllEvents_FD = 
      new TH1D("RecoEnergyAllEvents_FD",
	     "Selected Events with Reco Energy (FarDet)",
	       fNXBins,fXBins);
    mh->fRecoEnergyAllEvents_FD->Sumw2();
    (fMatrixHists[nueSys])[Background::EBackground(max_bg_index)] = mh;
    max_bg_index++;
  }
  gDirectory->cd("/");
}

NueMatrixHelper::~NueMatrixHelper()
{
}

void NueMatrixHelper::MakeANANUEPlots(Selection::Selection_t sel)
{
  if(!fFarChain || !fNearChain) return;  
  fCurSel = sel;

  std::map<NueSystematic*,
    std::map<Background::Background_t,MatrixHists*> >::iterator mapBeg = 
    fMatrixHists.begin();
  std::map<NueSystematic*,
    std::map<Background::Background_t,MatrixHists*> >::iterator mapEnd = 
    fMatrixHists.end();
  if(mapBeg==mapEnd) return;

  //ND:
  Int_t ndEntries = fNearChain->GetEntries();
  for(int i=0;i<ndEntries;i++){
    if(i%100000==0) cout << "Processing event " << i << " / " << ndEntries << endl;
    mapBeg = fMatrixHists.begin();
    while(mapBeg!=mapEnd) {
      fNearChain->GetEvent(i);
      Double_t totWeight = 1;
      if((mapBeg->first)) totWeight = (mapBeg->first)->UpdateRecord(fRecord,sel);
      if(this->PassCuts(sel)) {
	Double_t recoEnergy = this->GetNueEnergy(sel);
	//fill histos:
	std::map<Background::Background_t,MatrixHists*>::iterator mhBeg = mapBeg->second.begin();
	std::map<Background::Background_t,MatrixHists*>::iterator mhEnd = mapBeg->second.end();
	while(mhBeg!=mhEnd){
	  mhBeg->second->fRecoEnergyAllEvents_ND->Fill(recoEnergy,totWeight);
	  mhBeg++;
	}
	Background::Background_t bg = 
	  Background::TranslateFromMC(fRecord->mctrue.interactionType,
				      fRecord->mctrue.nuFlavor,
				      fRecord->mctrue.nonOscNuFlavor,
				      fRecord->fluxinfo.tptype);		
	(mapBeg->second)[bg]->fRecoEnergyCCOnlyEvents_ND->Fill(recoEnergy,totWeight);
	(mapBeg->second)[bg]->fTrueEnergyCCOnlyEvents_ND->Fill(fRecord->mctrue.nuEnergy,totWeight);
	(mapBeg->second)[bg]->fPurity_ND->Fill(recoEnergy,totWeight);
	(mapBeg->second)[bg]->fEfficiency_ND->Fill(fRecord->mctrue.nuEnergy,totWeight);
	(mapBeg->second)[bg]->fRecoVsTrueEnergy_ND->Fill(fRecord->mctrue.nuEnergy,
							 recoEnergy,totWeight);	
	Background::Background_t bg2 = 
	  Background::TranslateFromMC(fRecord->mctrue.interactionType,
				      fRecord->mctrue.nuFlavor,
				      fRecord->mctrue.nonOscNuFlavor);
	if(bg2!=bg) {
	  (mapBeg->second)[bg2]->fRecoEnergyCCOnlyEvents_ND->Fill(recoEnergy,totWeight);
	  (mapBeg->second)[bg2]->fTrueEnergyCCOnlyEvents_ND->Fill(fRecord->mctrue.nuEnergy,totWeight);
	  (mapBeg->second)[bg2]->fPurity_ND->Fill(recoEnergy,totWeight);
	  (mapBeg->second)[bg2]->fEfficiency_ND->Fill(fRecord->mctrue.nuEnergy,totWeight);
	  (mapBeg->second)[bg2]->fRecoVsTrueEnergy_ND->Fill(fRecord->mctrue.nuEnergy,
							    recoEnergy,totWeight);
	}
      }
    }
  }
  
  //FD:
  Int_t fdEntries = fFarChain->GetEntries();
  for(int i=0;i<fdEntries;i++){
    if(i%100000==0) cout << "Processing event " << i << " / " << ndEntries << endl;
    mapBeg = fMatrixHists.begin();
    while(mapBeg!=mapEnd) {
      fFarChain->GetEvent(i);
      Double_t totWeight = 1;
      if((mapBeg->first)) totWeight = (mapBeg->first)->UpdateRecord(fRecord,sel);
      if(this->PassCuts(sel)) {
	Double_t recoEnergy = this->GetNueEnergy(sel);
	//fill histos:
	std::map<Background::Background_t,MatrixHists*>::iterator mhBeg = mapBeg->second.begin();
	std::map<Background::Background_t,MatrixHists*>::iterator mhEnd = mapBeg->second.end();
	while(mhBeg!=mhEnd){
	  mhBeg->second->fRecoEnergyAllEvents_FD->Fill(recoEnergy,totWeight);
	  mhBeg++;
	}
	Background::Background_t bg = 
	  Background::TranslateFromMC(fRecord->mctrue.interactionType,
				      fRecord->mctrue.nuFlavor,
				      fRecord->mctrue.nonOscNuFlavor,
				      fRecord->fluxinfo.tptype);		
	(mapBeg->second)[bg]->fRecoEnergyCCOnlyEvents_FD->Fill(recoEnergy,totWeight);
	(mapBeg->second)[bg]->fTrueEnergyCCOnlyEvents_FD->Fill(fRecord->mctrue.nuEnergy,totWeight);
	(mapBeg->second)[bg]->fPurity_FD->Fill(recoEnergy,totWeight);
	(mapBeg->second)[bg]->fEfficiency_FD->Fill(fRecord->mctrue.nuEnergy,totWeight);
	(mapBeg->second)[bg]->fRecoVsTrueEnergy_FD->Fill(fRecord->mctrue.nuEnergy,
							 recoEnergy,totWeight);	
	Background::Background_t bg2 = 
	  Background::TranslateFromMC(fRecord->mctrue.interactionType,
				      fRecord->mctrue.nuFlavor,
				      fRecord->mctrue.nonOscNuFlavor);
	if(bg2!=bg) {
	  (mapBeg->second)[bg2]->fRecoEnergyCCOnlyEvents_FD->Fill(recoEnergy,totWeight);
	  (mapBeg->second)[bg2]->fTrueEnergyCCOnlyEvents_FD->Fill(fRecord->mctrue.nuEnergy,totWeight);
	  (mapBeg->second)[bg2]->fPurity_FD->Fill(recoEnergy,totWeight);
	  (mapBeg->second)[bg2]->fEfficiency_FD->Fill(fRecord->mctrue.nuEnergy,totWeight);
	  (mapBeg->second)[bg2]->fRecoVsTrueEnergy_FD->Fill(fRecord->mctrue.nuEnergy,
							    recoEnergy,totWeight);
	}
      }
    }
  }
  
  mapBeg = fMatrixHists.begin();
  while(mapBeg!=mapEnd) {
    std::map<Background::Background_t,MatrixHists*>::iterator mhBeg = mapBeg->second.begin();
    std::map<Background::Background_t,MatrixHists*>::iterator mhEnd = mapBeg->second.end();
    while(mhBeg!=mhEnd) {
      for(int i=1;i<=fNXBins;i++){ //loop over true bins
	for(int j=1;j<=fNXBins+1;j++){ //loop over reco bins
	  if(mhBeg->second->fRecoEnergyCCOnlyEvents_ND->GetBinContent(j)>0 && 
	     mhBeg->second->fRecoVsTrueEnergy_ND->GetBinContent(i,j)>0) {
	    Float_t error = ( mhBeg->second->fRecoVsTrueEnergy_ND->GetBinError(i,j) / 
			      mhBeg->second->fRecoVsTrueEnergy_ND->GetBinContent(i,j) );
	    mhBeg->second->fRecoVsTrueEnergy_ND->
	      SetBinContent(i,j,mhBeg->second->fRecoVsTrueEnergy_ND->GetBinContent(i,j)/
			    mhBeg->second->fRecoEnergyCCOnlyEvents_ND->GetBinContent(j));
	    mhBeg->second->fRecoVsTrueEnergy_ND->
	      SetBinError(i,j,error*mhBeg->second->fRecoVsTrueEnergy_ND->GetBinContent(i,j));
	  }
	  else {
	    mhBeg->second->fRecoVsTrueEnergy_ND->SetBinContent(i,j,0);
	    mhBeg->second->fRecoVsTrueEnergy_ND->SetBinError(i,j,0);
	  }
	  if(mhBeg->second->fRecoEnergyCCOnlyEvents_FD->GetBinContent(j)>0 && 
	     mhBeg->second->fRecoVsTrueEnergy_FD->GetBinContent(i,j)>0) {
	    Float_t error = ( mhBeg->second->fRecoVsTrueEnergy_FD->GetBinError(i,j) / 
			      mhBeg->second->fRecoVsTrueEnergy_FD->GetBinContent(i,j) );
	    mhBeg->second->fRecoVsTrueEnergy_FD->
	      SetBinContent(i,j,mhBeg->second->fRecoVsTrueEnergy_FD->GetBinContent(i,j)/
			    mhBeg->second->fRecoEnergyCCOnlyEvents_FD->GetBinContent(j));
	    mhBeg->second->fRecoVsTrueEnergy_FD->
	      SetBinError(i,j,error*mhBeg->second->fRecoVsTrueEnergy_FD->GetBinContent(i,j));
	  }
	  else {
	    mhBeg->second->fRecoVsTrueEnergy_FD->SetBinContent(i,j,0);
	    mhBeg->second->fRecoVsTrueEnergy_FD->SetBinError(i,j,0);
	  }
	}
      }
    }
  }  
}

/*
Bool_t NueMatrixHelper::MakeSNTPPlots(TChain *ndChain,TChain *fdChain)
{
  cout << "In MakeSNTPPlots" << endl;

  Zbeam beamWeighter;
  //skzp1:
  //Double_t beamWeighterPars[6] = {-0.441888,0.0303515,-0.0494623,
  //                              -2.63948e-07,-2.20597,0.0};
  //skzp2:
  //  Double_t beamWeighterPars[6] = {0.74087,0.23129,-0.0084731,
  //-2.014e-07,-2.9149,0.0};

  Zfluk reWeighter;
  reWeighter.UseParameterization("SKZP");
  //skzp1:
  //Double_t reWeightPars[5] = {0.535219,0.485084,1.00619,1.58987,0.870674};
  //skzp2:
  //  Double_t reWeightPars[7] = {0.55227,-4.3194,1.0035,0.49904,
  //0.92193,1.8081,0.86612};
  
  MCReweight *mcr = &MCReweight::Instance();
  NeugenWeightCalculator *n = new NeugenWeightCalculator();
  mcr->AddWeightCalculator(n);
  Registry *rwtconfig = new Registry();
  rwtconfig->UnLockValues();
  rwtconfig->UnLockKeys();
  rwtconfig->Set("neugen:config_name","MODBYRS");
  rwtconfig->Set("neugen:config_no",3);
  rwtconfig->LockValues();
  rwtconfig->LockKeys();
  mcr->ComputeWeight(0,rwtconfig);
  NuParent *np=0;

  if(ndChain && fdChain){
    NtpStRecord *strecord = NULL;
    ndChain->SetBranchAddress("NtpStRecord",&strecord);
    Int_t ndEntries = ndChain->GetEntries();
    cout << "Starting ND" << endl;
    for(int i=0;i<ndEntries;i++) {
      if(i%100000==0) cout << "Processing event " << i << " / " << ndEntries << endl;
      ndChain->GetEvent(i);
      TClonesArray& mcArray = *(strecord->mc);      
      Int_t mcEntries = mcArray.GetEntriesFast();
      for(int j=0;j<mcEntries;j++) {
	NtpMCTruth *ntpTruth = dynamic_cast<NtpMCTruth *>(mcArray[j]);	
	if(ntpTruth->vtxz>=1.0 && ntpTruth->vtxz<=5.0 && 
	   sqrt(pow(ntpTruth->vtxx-1.4885,2) + 
		pow(ntpTruth->vtxy-0.1397,2))<=1.0) {
	  if(ntpTruth->iaction==1 && ntpTruth->inu==14) {
	    Double_t flux_weight = reWeighter.GetWeight(ntpTruth->flux.tptype,
							ntpTruth->flux.tpx,
							ntpTruth->flux.tpy,
							ntpTruth->flux.tpz,reWeightPars);
	    for(int j=0;j<6;j++)
	      flux_weight *= beamWeighter.GetWeight(1,2,j+1,beamWeighterPars[j],
						    ntpTruth->p4neu[3]);
	    if(!doSKZP) flux_weight = 1.;
	    MCEventInfo ei;
	    ei.UseStoredXSec(false); 
	    ei.nuE=ntpTruth->p4neu[3]; ei.nuPx=ntpTruth->p4neu[0]; 
	    ei.nuPy=ntpTruth->p4neu[1]; ei.nuPz=ntpTruth->p4neu[2];
	    ei.tarE=ntpTruth->p4tgt[3]; ei.tarPx=ntpTruth->p4tgt[0]; 
	    ei.tarPy=ntpTruth->p4tgt[1]; ei.tarPz=ntpTruth->p4tgt[2];
	    ei.y=ntpTruth->y; ei.x=ntpTruth->x; ei.q2=ntpTruth->q2; ei.w2=ntpTruth->w2;
	    ei.iaction=ntpTruth->iaction; ei.inu=ntpTruth->inu;
	    ei.iresonance=ntpTruth->iresonance; 
	    ei.initial_state=GetInitialState(strecord,ntpTruth);
	    ei.nucleus=GetNucleus(int(ntpTruth->z),int(ntpTruth->a));
	    ei.had_fs=GetHadronicFinalState(strecord,ntpTruth);
	    Double_t xsec_weight = 1.;
	    if(doMODBYRS3) xsec_weight = mcr->ComputeWeight(&ei,np);

	    fTrueEnergyTrueCCFidEvents_ND->Fill(ntpTruth->p4neu[3],
						flux_weight*xsec_weight);
	  }
	}
      }
    }
 
    fdChain->SetBranchAddress("NtpStRecord",&strecord);
    Int_t fdEntries = fdChain->GetEntries();
    cout << "Starting FD" << endl;
    for(int i=0;i<fdEntries;i++) {
      if(i%100000==0) cout << "Processing event " << i << " / " << fdEntries << endl;
      fdChain->GetEvent(i);
      TClonesArray& mcArray = *(strecord->mc);      
      Int_t mcEntries = mcArray.GetEntriesFast();
      for(int j=0;j<mcEntries;j++) {
	NtpMCTruth *ntpTruth = dynamic_cast<NtpMCTruth *>(mcArray[j]);
	if(ntpTruth->vtxz>=0.5 && ntpTruth->vtxz<=28.0 && 
	   (ntpTruth->vtxz>=16.2 || ntpTruth->vtxz<=14.3) && 
	   (ntpTruth->vtxx*ntpTruth->vtxx + 
	    ntpTruth->vtxy*ntpTruth->vtxy)<=14.0) {
	  if(ntpTruth->iaction==1 && ntpTruth->inu==14) {
	    Double_t flux_weight = reWeighter.GetWeight(ntpTruth->flux.tptype,
							ntpTruth->flux.tpx,
							ntpTruth->flux.tpy,
							ntpTruth->flux.tpz,reWeightPars);
	    for(int j=0;j<6;j++)
	      flux_weight *= beamWeighter.GetWeight(2,2,j+1,beamWeighterPars[j],
						    ntpTruth->p4neu[3]);
	    if(!doSKZP) flux_weight = 1.;
	    MCEventInfo ei;
	    ei.UseStoredXSec(false); 
	    ei.nuE=ntpTruth->p4neu[3]; ei.nuPx=ntpTruth->p4neu[0]; 
	    ei.nuPy=ntpTruth->p4neu[1]; ei.nuPz=ntpTruth->p4neu[2];
	    ei.tarE=ntpTruth->p4tgt[3]; ei.tarPx=ntpTruth->p4tgt[0]; 
	    ei.tarPy=ntpTruth->p4tgt[1]; ei.tarPz=ntpTruth->p4tgt[2];
	    ei.y=ntpTruth->y; ei.x=ntpTruth->x; ei.q2=ntpTruth->q2; ei.w2=ntpTruth->w2;
	    ei.iaction=ntpTruth->iaction; ei.inu=ntpTruth->inu;
	    ei.iresonance=ntpTruth->iresonance; 
	    ei.initial_state=GetInitialState(strecord,ntpTruth);
	    ei.nucleus=GetNucleus(int(ntpTruth->z),int(ntpTruth->a)); 
	    ei.had_fs=GetHadronicFinalState(strecord,ntpTruth);
	    Double_t xsec_weight = 1.;
	    if(doMODBYRS3) xsec_weight = mcr->ComputeWeight(&ei,np);

	    fTrueEnergyTrueCCFidEvents_FD->Fill(ntpTruth->p4neu[3],
						flux_weight*xsec_weight);
	  }
	}
      }
    }
    return true;
  }
  return false;
}

Bool_t NueMatrixHelper::MakeFLUXPlots(TChain *fluxChain,TGraph *xsec_Graph,
				   TH1F *mikehist)
{
  cout << "In MakeFLUXPlots" << endl;

  Zbeam beamWeighter;
  //skzp1:
  //Double_t beamWeighterPars[6] = {-0.441888,0.0303515,-0.0494623,
  //                              -2.63948e-07,-2.20597,0.0};
  //skzp2:
  //  Double_t beamWeighterPars[6] = {0.74087,0.23129,-0.0084731,
  //-2.014e-07,-2.9149,0.0};

  Zfluk reWeighter;
  reWeighter.UseParameterization("SKZP");
  //skzp1:
  //Double_t reWeightPars[5] = {0.535219,0.485084,1.00619,1.58987,0.870674};
  //skzp2:
  //  Double_t reWeightPars[7] = {0.55227,-4.3194,1.0035,0.49904,
  //0.92193,1.8081,0.86612};

  
  if(fluxChain) {
    Int_t          run;
    Int_t          evtno;
    Float_t        Ndxdz;
    Float_t        Ndydz;
    Float_t        Npz;
    Float_t        Nenergy;
    Float_t        Ndxdznea;
    Float_t        Ndydznea;
    Float_t        Nenergyn;
    Float_t        Nwtnear;
    Float_t        Ndxdzfar;
    Float_t        Ndydzfar;
    Float_t        Nenergyf;
    Float_t        Nwtfar;
    Int_t          Norig;
    Int_t          Ndecay;
    Int_t          Ntype;
    Float_t        Vx;
    Float_t        Vy;
    Float_t        Vz;
    Float_t        pdpx;
    Float_t        pdpy;
    Float_t        pdpz;
    Float_t        ppdxdz;
    Float_t        ppdydz;
    Float_t        pppz;
    Float_t        ppenergy;
    Int_t          ppmedium;
    Int_t          ptype;
    Float_t        ppvx;
    Float_t        ppvy;
    Float_t        ppvz;
    Float_t        muparpx;
    Float_t        muparpy;
    Float_t        muparpz;
    Float_t        mupare;
    Float_t        Necm;
    Float_t        Nimpwt;
    Float_t        xpoint;
    Float_t        ypoint;
    Float_t        zpoint;
    Float_t        tvx;
    Float_t        tvy;
    Float_t        tvz;
    Float_t        tpx;
    Float_t        tpy;
    Float_t        tpz;
    Int_t          tptype;
    Int_t          tgen;
    Int_t          tgptype;
    Float_t        tgppx;
    Float_t        tgppy;
    Float_t        tgppz;
    Float_t        tprivx;
    Float_t        tprivy;
    Float_t        tprivz;
    Float_t        beamx;
    Float_t        beamy;
    Float_t        beamz;
    Float_t        beampx;
    Float_t        beampy;
    Float_t        beampz;
    fluxChain->SetBranchAddress("run",&run);
    fluxChain->SetBranchAddress("evtno",&evtno);
    fluxChain->SetBranchAddress("Ndxdz",&Ndxdz);
    fluxChain->SetBranchAddress("Ndydz",&Ndydz);
    fluxChain->SetBranchAddress("Npz",&Npz);
    fluxChain->SetBranchAddress("Nenergy",&Nenergy);
    fluxChain->SetBranchAddress("Ndxdznea",&Ndxdznea);
    fluxChain->SetBranchAddress("Ndydznea",&Ndydznea);
    fluxChain->SetBranchAddress("Nenergyn",&Nenergyn);
    fluxChain->SetBranchAddress("Nwtnear",&Nwtnear);
    fluxChain->SetBranchAddress("Ndxdzfar",&Ndxdzfar);
    fluxChain->SetBranchAddress("Ndydzfar",&Ndydzfar);
    fluxChain->SetBranchAddress("Nenergyf",&Nenergyf);
    fluxChain->SetBranchAddress("Nwtfar",&Nwtfar);
    fluxChain->SetBranchAddress("Norig",&Norig);
    fluxChain->SetBranchAddress("Ndecay",&Ndecay);
    fluxChain->SetBranchAddress("Ntype",&Ntype);
    fluxChain->SetBranchAddress("Vx",&Vx);
    fluxChain->SetBranchAddress("Vy",&Vy);
    fluxChain->SetBranchAddress("Vz",&Vz);
    fluxChain->SetBranchAddress("pdpx",&pdpx);
    fluxChain->SetBranchAddress("pdpy",&pdpy);
    fluxChain->SetBranchAddress("pdpz",&pdpz);
    fluxChain->SetBranchAddress("ppdxdz",&ppdxdz);
    fluxChain->SetBranchAddress("ppdydz",&ppdydz);
    fluxChain->SetBranchAddress("pppz",&pppz);
    fluxChain->SetBranchAddress("ppenergy",&ppenergy);
    fluxChain->SetBranchAddress("ppmedium",&ppmedium);
    fluxChain->SetBranchAddress("ptype",&ptype);
    fluxChain->SetBranchAddress("ppvx",&ppvx);
    fluxChain->SetBranchAddress("ppvy",&ppvy);
    fluxChain->SetBranchAddress("ppvz",&ppvz);
    fluxChain->SetBranchAddress("muparpx",&muparpx);
    fluxChain->SetBranchAddress("muparpy",&muparpy);
    fluxChain->SetBranchAddress("muparpz",&muparpz);
    fluxChain->SetBranchAddress("mupare",&mupare);
    fluxChain->SetBranchAddress("Necm",&Necm);
    fluxChain->SetBranchAddress("Nimpwt",&Nimpwt);
    fluxChain->SetBranchAddress("xpoint",&xpoint);
    fluxChain->SetBranchAddress("ypoint",&ypoint);
    fluxChain->SetBranchAddress("zpoint",&zpoint);
    fluxChain->SetBranchAddress("tvx",&tvx);
    fluxChain->SetBranchAddress("tvy",&tvy);
    fluxChain->SetBranchAddress("tvz",&tvz);
    fluxChain->SetBranchAddress("tpx",&tpx);
    fluxChain->SetBranchAddress("tpy",&tpy);
    fluxChain->SetBranchAddress("tpz",&tpz);
    fluxChain->SetBranchAddress("tptype",&tptype);
    fluxChain->SetBranchAddress("tgen",&tgen);
    fluxChain->SetBranchAddress("tgptype",&tgptype);
    fluxChain->SetBranchAddress("tgppx",&tgppx);
    fluxChain->SetBranchAddress("tgppy",&tgppy);
    fluxChain->SetBranchAddress("tgppz",&tgppz);
    fluxChain->SetBranchAddress("tprivx",&tprivx);
    fluxChain->SetBranchAddress("tprivy",&tprivy);
    fluxChain->SetBranchAddress("tprivz",&tprivz);
    fluxChain->SetBranchAddress("beamx",&beamx);
    fluxChain->SetBranchAddress("beamy",&beamy);
    fluxChain->SetBranchAddress("beamz",&beamz);
    fluxChain->SetBranchAddress("beampx",&beampx);
    fluxChain->SetBranchAddress("beampy",&beampy);
    fluxChain->SetBranchAddress("beampz",&beampz);
    Int_t entries = fluxChain->GetEntries();
    for(int i=0;i<entries;i++){
      if(i%100000==0) cout << "Processing event " << i 
			   << " / " << entries << endl;
      fluxChain->GetEvent(i);
      if(Ntype!=56) continue;
      Double_t flux_weight_near = reWeighter.GetWeight(tptype,tpx,tpy,tpz,
						       reWeightPars);
      Double_t flux_weight_far = reWeighter.GetWeight(tptype,tpx,tpy,tpz,
						      reWeightPars);

      for(int j=0;j<6;j++) {
	flux_weight_near *= beamWeighter.GetWeight(1,2,j+1,beamWeighterPars[j],
						  Nenergyn);
	flux_weight_far *= beamWeighter.GetWeight(2,2,j+1,beamWeighterPars[j],
						  Nenergyf);
      }
      if(!doSKZP) {
	flux_weight_near = 1.;
	flux_weight_far = 1.;
      }

      Double_t xsec_nd = double(xsec_Graph->Eval(Nenergyn,0,""));
      Double_t xsec_fd = double(xsec_Graph->Eval(Nenergyf,0,""));
      if(xsec_nd<=0.0) xsec_nd = 0.0;
      if(xsec_fd<=0.0) xsec_fd = 0.0;

      fTrueEnergyNuFlux_ND->Fill(Nenergyn,Nimpwt*Nwtnear*flux_weight_near);
      fTrueEnergyNuFlux_FD->Fill(Nenergyf,Nimpwt*Nwtfar*flux_weight_far);      
      fFDVsNDMatrix->Fill(Nenergyn,Nenergyf,Nimpwt*Nwtfar*flux_weight_far);

      fTrueEnergyCCFlux_ND->Fill(Nenergyn,Nimpwt*Nwtnear*xsec_nd*flux_weight_near);
      fTrueEnergyCCFlux_FD->Fill(Nenergyf,Nimpwt*Nwtfar*xsec_fd*flux_weight_far);
      fFDVsNDMatrixXSec->Fill(Nenergyn,Nenergyf,Nimpwt*Nwtfar*xsec_fd*flux_weight_far);

      //do reweighting for front face of ND:
      BeamVar beamvar;
      beamvar.brun = run;
      beamvar.bevtno = evtno;
      beamvar.bNdxdz = Ndxdz;
      beamvar.bNdydz = Ndydz;
      beamvar.bNpz = Npz;
      beamvar.bNenergy = Nenergy;
      beamvar.bNdxdznea = Ndxdznea;
      beamvar.bNdydznea = Ndydznea;
      beamvar.bNenergyn = Nenergyn;
      beamvar.bNwtnear = Nwtnear;
      beamvar.bNdxdzfar = Ndxdzfar;
      beamvar.bNdydzfar = Ndydzfar;
      beamvar.bNenergyf = Nenergyf;
      beamvar.bNwtfar = Nwtfar;
      beamvar.bNorig = Norig;
      beamvar.bNdecay = Ndecay;
      beamvar.bNtype = Ntype;
      beamvar.bVx = Vx;
      beamvar.bVy = Vy;
      beamvar.bVz = Vz;
      beamvar.bpdpx = pdpx;
      beamvar.bpdpy = pdpy;
      beamvar.bpdpz = pdpz;
      beamvar.bppdxdz = ppdxdz;
      beamvar.bppdydz = ppdydz;
      beamvar.bpppz = pppz;
      beamvar.bppenergy = ppenergy;
      beamvar.bppmedium = ppmedium;
      beamvar.bptype = ptype;
      beamvar.bppvx = ppvx;
      beamvar.bppvy = ppvy;
      beamvar.bppvz = ppvz;
      beamvar.bmuparpx = muparpx;
      beamvar.bmuparpy = muparpy;
      beamvar.bmuparpz = muparpz;
      beamvar.bmupare = mupare;
      beamvar.bNecm = Necm;
      beamvar.bNimpwt = Nimpwt;
      beamvar.bxpoint = xpoint;
      beamvar.bypoint = ypoint;
      beamvar.bzpoint = zpoint;
      beamvar.btvx = tvx;
      beamvar.btvy = tvy;
      beamvar.btvz = tvz;
      beamvar.btpx = tpx;
      beamvar.btpy = tpy;
      beamvar.btpz = tpz;
      beamvar.btptype = tptype;
      beamvar.btgen = tgen;
      beamvar.btgptype = tgptype;
      beamvar.btgppx = tgppx;
      beamvar.btgppy = tgppy;
      beamvar.btgppz = tgppz;
      beamvar.btprivx = tprivx;
      beamvar.btprivy = tprivy;
      beamvar.btprivz = tprivz;
      beamvar.bbeamx = beamx;
      beamvar.bbeamy = beamy;
      beamvar.bbeamz = beamz;
      beamvar.bbeampx = beampx;
      beamvar.bbeampy = beampy;
      beamvar.bbeampz = beampz;

      // beam center (0,0) in detector coord system:
      // (148.844,13.97)
      // z offset is 104000
      // beamdydz = -0.0581

      //re-use neutrinos a few times to spread out importance weights...
      for(int j=0;j<10;j++){
	//Generate random x,y,z in detector coords:
	Double_t r = 1000;
	Double_t rand_x = 0;
	Double_t rand_y = 0;
	Double_t rand_z = gRandom->Uniform(100.0,500.0);
	while(r>100.) {
	  rand_x = gRandom->Uniform(-100,100);
	  rand_y = gRandom->Uniform(-100,100);
	  r = TMath::Sqrt(rand_x*rand_x + rand_y*rand_y);
	}

	//now translate to beam coords:
	Double_t X = rand_x;
	Double_t Y = (rand_y*TMath::Cos(3.323155*TMath::DegToRad()) + 
		      rand_z*TMath::Sin(3.323155*TMath::DegToRad()));
	Double_t Z = 104000.+(rand_z*TMath::Cos(3.323155*TMath::DegToRad()) - 
			      rand_y*TMath::Sin(3.323155*TMath::DegToRad()));
	ArrAy newvals_nd = nuwte(X,Y,Z,beamvar);

	//get new fd weight and energy:
	ArrAy newvals_fd = nuwte(0,0,73200000.,beamvar);

	//recalc xsec
	xsec_nd = double(xsec_Graph->Eval(newvals_nd.new_ene,0,""));
	if(xsec_nd<=0.0) xsec_nd = 0.0; 

	xsec_fd = double(xsec_Graph->Eval(newvals_fd.new_ene,0,""));
	if(xsec_fd<=0.0) xsec_fd = 0.0; 

	//recalc flux weight
	flux_weight_near = reWeighter.GetWeight(tptype,tpx,tpy,tpz,reWeightPars);
	flux_weight_far = reWeighter.GetWeight(tptype,tpx,tpy,tpz,reWeightPars);
	for(int j=0;j<6;j++) {
	  flux_weight_near *= beamWeighter.GetWeight(1,2,j+1,beamWeighterPars[j],
						     newvals_nd.new_ene);
	  flux_weight_far *= beamWeighter.GetWeight(2,2,j+1,beamWeighterPars[j],
						    newvals_fd.new_ene);
	}
	if(!doSKZP) {
	  flux_weight_near = 1.;
	  flux_weight_far = 1.;
	}

	//fill matrix
	fFDVsNDMatrixRW->Fill(newvals_nd.new_ene,newvals_fd.new_ene,
			      Nimpwt*newvals_fd.new_weight*flux_weight_far);
	fTrueEnergyNuFluxRW_ND->Fill(newvals_nd.new_ene,
				     Nimpwt*newvals_nd.new_weight*flux_weight_near);
	fTrueEnergyNuFluxRW_FD->Fill(newvals_fd.new_ene,
				     Nimpwt*newvals_fd.new_weight*flux_weight_far);
	
	fFDVsNDMatrixXSecRW->Fill(newvals_nd.new_ene,newvals_fd.new_ene,
				  Nimpwt*newvals_fd.new_weight*xsec_fd*flux_weight_far);
	fTrueEnergyCCFluxRW_ND->Fill(newvals_nd.new_ene,
				     Nimpwt*newvals_nd.new_weight*xsec_nd*flux_weight_near);
	fTrueEnergyCCFluxRW_FD->Fill(newvals_fd.new_ene,
				     Nimpwt*newvals_fd.new_weight*xsec_fd*flux_weight_far);
	
      }
    }
    //re-norm the histo:
    for(int i=1;i<=fNBins;i++){
      for(int j=1;j<=fNBins;j++){
	if( fTrueEnergyNuFlux_ND->GetBinContent(i)>0 && 
	    fFDVsNDMatrix->GetBinContent(i,j)>0 ) {
	  Float_t error = TMath::Power(fFDVsNDMatrix->GetBinError(i,j) / 
				       fFDVsNDMatrix->GetBinContent(i,j),2);
	  error = TMath::Sqrt(error);
	  fFDVsNDMatrix->SetBinContent(i,j,
				       fFDVsNDMatrix->GetBinContent(i,j)/
				       fTrueEnergyNuFlux_ND->GetBinContent(i));
	  fFDVsNDMatrix->SetBinError(i,j,error * 
				     fFDVsNDMatrix->GetBinContent(i,j));
	}
	else {
	  fFDVsNDMatrix->SetBinContent(i,j,0.);
	  fFDVsNDMatrix->SetBinError(i,j,0.);
	}

	if( fTrueEnergyNuFluxRW_ND->GetBinContent(i)>0 && 
	    fFDVsNDMatrixRW->GetBinContent(i,j)>0 ) {
	  Float_t error = TMath::Power(fFDVsNDMatrixRW->GetBinError(i,j) / 
				       fFDVsNDMatrixRW->GetBinContent(i,j),2);
	  error = TMath::Sqrt(error);
	  fFDVsNDMatrixRW->SetBinContent(i,j,
					 fFDVsNDMatrixRW->GetBinContent(i,j)/
					 fTrueEnergyNuFluxRW_ND->GetBinContent(i));
	  fFDVsNDMatrixRW->SetBinError(i,j,error * 
				       fFDVsNDMatrixRW->GetBinContent(i,j));
	}
	else {
	  fFDVsNDMatrix->SetBinContent(i,j,0.);
	  fFDVsNDMatrix->SetBinError(i,j,0.);
	}

	if( fTrueEnergyCCFlux_ND->GetBinContent(i)>0 &&
	    fFDVsNDMatrixXSec->GetBinContent(i,j)>0 ) {
	  Float_t error = TMath::Power(fFDVsNDMatrixXSec->GetBinError(i,j) / 
				       fFDVsNDMatrixXSec->GetBinContent(i,j),2);
	  error = TMath::Sqrt(error);
	  fFDVsNDMatrixXSec->SetBinContent(i,j,
					   fFDVsNDMatrixXSec->GetBinContent(i,j)/
					   fTrueEnergyCCFlux_ND->GetBinContent(i));
	  fFDVsNDMatrixXSec->SetBinError(i,j,error * 
					 fFDVsNDMatrixXSec->GetBinContent(i,j));
	}
	else {
	  fFDVsNDMatrixXSec->SetBinContent(i,j,0.);
	  fFDVsNDMatrixXSec->SetBinError(i,j,0.);
	}

	if( fTrueEnergyCCFluxRW_ND->GetBinContent(i)>0 &&
	    fFDVsNDMatrixXSecRW->GetBinContent(i,j)>0 ) {
	  Float_t error = TMath::Power(fFDVsNDMatrixXSecRW->GetBinError(i,j) / 
				       fFDVsNDMatrixXSecRW->GetBinContent(i,j),2);
	  error = TMath::Sqrt(error);
	  fFDVsNDMatrixXSecRW->SetBinContent(i,j,
					     fFDVsNDMatrixXSecRW->GetBinContent(i,j)/
					     fTrueEnergyCCFluxRW_ND->GetBinContent(i));
	  fFDVsNDMatrixXSecRW->SetBinError(i,j,error *
					   fFDVsNDMatrixXSecRW->GetBinContent(i,j));
	}
	else {
	  fFDVsNDMatrixXSecRW->SetBinContent(i,j,0.);
	  fFDVsNDMatrixXSecRW->SetBinError(i,j,0.);
	}
      }
    }
    return true;
  }
  return false;
}

Bool_t NueMatrixHelper::ReadXSecPlots(TFile *mikeFile) 
{
  cout << "In ReadXSecPlots" << endl;
  if(mikeFile && mikeFile->IsOpen() && !mikeFile->IsZombie()) {
    TH1F *mike_xsec_cc = (TH1F*) mikeFile->Get("h_numu_cc_tot");
    if(!mike_xsec_cc) return false;
    for(int i=1;i<=fNBins;i++) {
      Float_t xval = fXSec_CC->GetBinCenter(i);
      Int_t theBin = mike_xsec_cc->FindBin(xval);
      Float_t est = 0;
      if(TMath::Abs(mike_xsec_cc->GetBinCenter(theBin)-xval)<1e-6) 
	est = mike_xsec_cc->GetBinContent(theBin);
      else if(theBin==1 || mike_xsec_cc->GetBinCenter(theBin)>xval) {
	est = mike_xsec_cc->GetBinContent(theBin);
	est += ( mike_xsec_cc->GetBinContent(theBin+1) - 
		 mike_xsec_cc->GetBinContent(theBin) ) * 
	  ( (xval - mike_xsec_cc->GetBinCenter(theBin)) / 
	    (mike_xsec_cc->GetBinCenter(theBin+1) - 
	     mike_xsec_cc->GetBinCenter(theBin)) );
      }
      else if(theBin==mike_xsec_cc->GetNbinsX() || 
	      mike_xsec_cc->GetBinCenter(theBin)<xval){
	est = mike_xsec_cc->GetBinContent(theBin-1);
	est += ( mike_xsec_cc->GetBinContent(theBin) - 
		 mike_xsec_cc->GetBinContent(theBin-1) ) * 
	  ( (xval - mike_xsec_cc->GetBinCenter(theBin-1)) / 
	    (mike_xsec_cc->GetBinCenter(theBin) - 
	     mike_xsec_cc->GetBinCenter(theBin-1)) );
      }
      fXSec_CC->SetBinContent(i,est);
    }
    mikeFile->Close();
    return true;
  }
  return false;
}
*/

void NueMatrixHelper::WriteFile(std::string tag)
{
  std::map<NueSystematic*,
    std::map<Background::Background_t,MatrixHists*> >::iterator mapBeg = 
    fMatrixHists.begin();
  std::map<NueSystematic*,
    std::map<Background::Background_t,MatrixHists*> >::iterator mapEnd = 
    fMatrixHists.end();
  if(mapBeg==mapEnd) return;

  std::string filename = "MatrixHelper_" + tag + ".root";
  TFile *file = new TFile(filename.c_str(),"RECREATE");
  file->cd();

  char selection[256]; 
  sprintf(selection,"%s",Selection::AsString(fCurSel));  
  TTree *tree = new TTree("matrixtree","matrixtree");
  tree->Branch("Selection",selection,"Selection/C");
  tree->Branch("nearPOT",&fNearPOT,"nearPOT/D");
  tree->Branch("farPOT",&fFarPOT,"farPOT/D");

  while(mapBeg!=mapEnd){
    std::map<Background::Background_t,MatrixHists*>::iterator Matbeg = (mapBeg->second).begin();
    std::map<Background::Background_t,MatrixHists*>::iterator Matend = (mapBeg->second).end();
    while(Matbeg!=Matend) {
      Matbeg->second->fDirectory->Write();
      Matbeg++;
    }
    if((mapBeg->first)) (mapBeg->first)->MakeBranches(tree);
    tree->Fill();
    mapBeg++;
  }
  tree->Write();
  delete file;
}
