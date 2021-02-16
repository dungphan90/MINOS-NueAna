#include "NueAna/Extrapolation/NueFNExtrapolation.h"
#include "TTree.h"
#include <stdlib.h>
#include <string>
#include <iostream>

using namespace std;

NueFNExtrapolation::NueFNExtrapolation() :
  NueExtrapolation(),
  fHelperRatio(0),fTreeEntry(-1)
{
  fExtrapMethod = Extrapolation::kFN;
}

NueFNExtrapolation::NueFNExtrapolation(std::string name) :
  NueExtrapolation(),
  fHelperRatio(0),fTreeEntry(-1)
{
  fExtrapMethod = Extrapolation::kFN;
  this->SetHelperFile(name);
  fHelperTree = (TTree*) fHelperFile->Get("energytree");  
}

NueFNExtrapolation::~NueFNExtrapolation()
{
  if(fHelperRatio) delete fHelperRatio;
  if(fHelperFile) fHelperFile->Close();
}

TH1* NueFNExtrapolation::Extrapolate(const char *histName,Bool_t extrapNtoF)
{
  if(!fHelperFile) return NULL;
  //get appropriate F/N ratio histogram:
  this->GetRatio("RecoEnergy");

  //get background histogram
  TH1D* extrapHist = NULL;
  if(fBg->GetHist()->InheritsFrom("TH2")){
    TH2* hist = (TH2*) fBg->GetHist();
    extrapHist = (TH1D*) hist->ProjectionX(histName);
  }
  else extrapHist = (TH1D*) fBg->GetHist()->Clone(histName);  

  //apply ratio:
  for(int i=0;i<extrapHist->GetNbinsX()+1;i++){
    Double_t ratio = fHelperRatio->GetBinContent(i);
    if(!extrapNtoF && ratio!=0) ratio = 1/ratio;
    extrapHist->SetBinContent(i,extrapHist->GetBinContent(i)*ratio);
    extrapHist->SetBinError(i,extrapHist->GetBinError(i)*ratio);
  }
  return extrapHist;
}

NueSystematic *NueFNExtrapolation::GetCurrentSystematic()
{
  if(!fHelperFile || !fHelperTree) return NULL;
  if(fTreeEntry==-1) return NULL;
  fHelperTree->ResetBranchAddresses();

  char sysname[256] = "";
  //Double_t beamWeightPars[6] = {};
//  Double_t hadProdPars[7] = {};
  Double_t theta12=0,theta23=0,theta13=0;
  Double_t deltaMSq23=0,deltaMSq12=0;
  Double_t deltaCP=0;Int_t massHierarchy=0;
  fHelperTree->SetBranchAddress("SysName",sysname);
  fHelperTree->SetBranchAddress("Theta12",&theta12);
  fHelperTree->SetBranchAddress("Theta23",&theta23);
  fHelperTree->SetBranchAddress("Theta13",&theta13);
  fHelperTree->SetBranchAddress("DeltaMSq23",&deltaMSq23);
  fHelperTree->SetBranchAddress("DeltaMSq12",&deltaMSq12);
  fHelperTree->SetBranchAddress("DeltaCP",&deltaCP);
  fHelperTree->SetBranchAddress("MassHierarchy",&massHierarchy);
  fHelperTree->GetEntry(fTreeEntry);

  NueSystematic *nueSys = new NueSystematic(string(sysname));
//  nueSys->SetSKZPParams(beamWeightPars,hadProdPars);
  nueSys->SetOscParams(theta12,theta23,theta13,deltaMSq12,
		       deltaMSq23,deltaCP,massHierarchy);
  
  Int_t max_sys_index = 0;
  while(strcmp(Systematic::AsString(Systematic::ESystematic(max_sys_index)),
	       "?Unknown?")!=0) {
    fHelperTree->ResetBranchAddresses();
    Double_t tempDouble = 0;
    fHelperTree->SetBranchAddress(Systematic::AsString(Systematic::ESystematic(max_sys_index)),
				  &tempDouble);
    fHelperTree->GetEntry(fTreeEntry);
    if(tempDouble>-9998) {
      nueSys->AddSystematic(Systematic::ESystematic(max_sys_index),tempDouble);
    }
    max_sys_index++;
  }
  return nueSys;
}

TH1 *NueFNExtrapolation::GetRatio(string distname) 
{
  if(!fBg) return NULL;
  if(!fHelperFile) return NULL;
  if(fHelperRatio) {delete fHelperRatio; fHelperRatio = NULL;}

  fHelperTree->ResetBranchAddresses();
  Double_t nearPOT  = 0;  fHelperTree->SetBranchAddress("nearPOT",&nearPOT);
  Double_t farPOT   = 0;  fHelperTree->SetBranchAddress("farPOT",&farPOT);
  char sysname[256] = ""; fHelperTree->SetBranchAddress("SysName",sysname);
  fHelperTree->GetEntry(fTreeEntry);
  //nearPOT = 2982*400*2.568e13;  
  //farPOT  = 990.*1.02716e20/3.;

  char fddir[256] = "";
  sprintf(fddir,"%s_%s/FD_%s",Background::AsString(fBg->GetBackground()),
	  sysname,distname.c_str());
  char nddir[256] = "";
  if(fBg->GetBackground()==Background::kNuTauCC || 
     fBg->GetBackground()==Background::kNueCC) {
    //these are appearances in FD, so cannot construct F/N ratio directly
    //instead assume that ND numuCC spectrum must have been passed if this
    //extrapolation method was chosen and so construct ratio using ND numuCC
    //as denominator
    sprintf(nddir,"%s_%s/ND_%s",Background::AsString(Background::kNuMuCC),
	    sysname,distname.c_str());
  }
  else sprintf(nddir,"%s_%s/ND_%s",Background::AsString(fBg->GetBackground()),
	       sysname,distname.c_str());

  TH1 *helperHistFD = (TH1*) fHelperFile->Get(fddir);    
  TH1 *helperHistND = (TH1*) fHelperFile->Get(nddir);    
  fHelperRatio = (TH1D*) helperHistFD->Clone("helperRatio");
  fHelperRatio->SetDirectory(0);
  fHelperRatio->Divide(helperHistND);
  if(nearPOT!=0) fHelperRatio->Scale(nearPOT/farPOT);
  return (TH1*) fHelperRatio;
}
