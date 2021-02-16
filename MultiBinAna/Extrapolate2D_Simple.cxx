#define Extrapolate2D_Simple_C

#include "NueAna/MultiBinAna/Extrapolate2D_Simple.h"

Extrapolate2D_Simple::Extrapolate2D_Simple()
{
  int i;
  
  nReco = 100;
  nRecoCC = 100;
  for(i=0;i<nReco+1;i++)
  {
    RecoEdges.push_back(i*1.);
    RecoCCEdges.push_back(i*1.);
  }
  nTrue = 1200;
  for(i=0;i<nTrue+1;i++)
  {
    TrueEdges.push_back(i*0.1);
  }
  nPID = 20;
  for(i=0;i<nPID+1;i++)
  {
    PIDEdges.push_back(i*0.05);
  }
  
  SetNearPOT();
  SetFarPOT();
  FNforBeamNue = false;
  SetRunPeriod();
  InitializeOscCalc();
  
  FDComponents.clear();
  FDComponents.push_back(Background::kNC);
  FDComponents.push_back(Background::kNuMuCC);
  FDComponents.push_back(Background::kNuTauCC);
  FDComponents.push_back(Background::kNueCC);
  FDComponents.push_back(Background::kBNueCC);
  
  Theta13=0;
  Theta12=0;
  Theta23=0;
  DeltaMSq23=0;
  DeltaMSq12=0;
  DeltaCP=0;
  
  ReadExtrapFromFile = false;
  
  MRE_infile = "NULL";
  ReadMREFromFile = false;
  
  RebinE = false;
  RebinP = false;
  
  PrintResult = false;
  
  WriteOutput = false;
  
  Init = false;//ie haven't called GetPrediction() yet
  
  ReadError = false;
  
  UseInputWithNo3D = false;
  
  Oscillated = false;//haven't called OscillatePrediction() yet
  
  UseSeparateNuNuBar=false;
  
  return;
}
Extrapolate2D_Simple::~Extrapolate2D_Simple()
{
}
void Extrapolate2D_Simple::ReadNDDataFile()
{
  
  //should have file with 2d total nd data and nd mc and nd bnue mc component
  //rebin
  //subtract bnue from data and mc and take ratio
  //scale ndmc cc and nc by ratio
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(NDData_infile.c_str())))
  {
    cout<<"Failure to read ND data file."<<endl;
    ReadError = true;
    return;
  }
  
  double *r = new double[nReco+1];
  double *p = new double[nPID+1];
  int i;
  for(i=0;i<nReco+1;i++)
  {
    r[i] = RecoEdges.at(i);
  }
  for(i=0;i<nPID+1;i++)
  {
    p[i] = PIDEdges.at(i);
  }
  
  TFile *f = new TFile(gSystem->ExpandPathName(NDData_infile.c_str()),"READ");
  
  TH2D *h;
  
  NDData_Total = (TH2D*)f->Get(Form("NDData_Total_%s",PID.c_str()));
  NDMC_Total = (TH2D*)f->Get(Form("NDMC_Total_%s",PID.c_str()));
  NDMC_BNueCC = (TH2D*)f->Get(Form("NDMC_BNueCC_%s",PID.c_str()));
  NDMC_NuMuCC = (TH2D*)f->Get(Form("NDMC_NuMuCC_%s",PID.c_str()));
  NDMC_NC = (TH2D*)f->Get(Form("NDMC_NC_%s",PID.c_str()));
  
  if(NDData_Total->GetNbinsX()!=nPID)
  {
    RebinP = true;
  }
  if(NDData_Total->GetNbinsY()!=nReco)
  {
    RebinE = true;
  }
  
  if(RebinE || RebinP)
  {
    MBH.Rebin2DHist(NDData_Total,nPID,p,nReco,r);
    MBH.Rebin2DHist(NDMC_Total,nPID,p,nReco,r);
    MBH.Rebin2DHist(NDMC_BNueCC,nPID,p,nReco,r);
    MBH.Rebin2DHist(NDMC_NuMuCC,nPID,p,nReco,r);
    MBH.Rebin2DHist(NDMC_NC,nPID,p,nReco,r);
  }
  
  NDData_Ratio = (TH2D*)NDData_Total->Clone("NDData_Ratio");
  NDData_Ratio->Add(NDMC_BNueCC,-1.);
  h = (TH2D*)NDMC_Total->Clone("h");
  h->Add(NDMC_BNueCC,-1.);
  NDData_Ratio->Divide(h);
  
  NDData[Background::kBNueCC] = (TH2D*)NDMC_BNueCC->Clone("NDData_BNueCC");
  NDData[Background::kNC] = (TH2D*)NDMC_NC->Clone("NDData_NC");
  NDData[Background::kNC]->Multiply(NDData_Ratio);
  NDData[Background::kNuMuCC] = (TH2D*)NDMC_NuMuCC->Clone("NDData_NuMuCC");
  NDData[Background::kNuMuCC]->Multiply(NDData_Ratio);
  
  return;
}
void Extrapolate2D_Simple::ReadFNFile()
{
  if(ReadError) return;
  
  FD_RecoVsPID.clear();
  FD_True2Reco_Fid.clear();
  FD_Eff.clear();
  FNRatio.clear();
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(FN_infile.c_str())))
  {
    cout<<"Failure to read F/N file."<<endl;
    ReadError = true;
    return;
  }
  
  TFile *f = new TFile(FN_infile.c_str(),"READ");
  
  double np,fp;
  TTree *tree = (TTree*)f->Get("paramtree");
  tree->SetBranchAddress("nearPOT",&np);
  tree->SetBranchAddress("farPOT",&fp);
  tree->SetBranchAddress("Theta12",&Theta12);
  tree->SetBranchAddress("Theta13",&Theta13);
  tree->SetBranchAddress("Theta23",&Theta23);
  tree->SetBranchAddress("DeltaMSq23",&DeltaMSq23);
  tree->SetBranchAddress("DeltaMSq12",&DeltaMSq12);
  tree->SetBranchAddress("DeltaCP",&DeltaCP);
  tree->GetEntry(0);
  
  osc.SetOscParam(OscPar::kTh12,Theta12);
  osc.SetOscParam(OscPar::kTh13,Theta13);
  osc.SetOscParam(OscPar::kTh23,Theta23);
  osc.SetOscParam(OscPar::kDeltaM23,DeltaMSq23);
  osc.SetOscParam(OscPar::kDeltaM12,DeltaMSq12);
  osc.SetOscParam(OscPar::kDelta,DeltaCP);
  SetNominalOscProb();
  
  TH3D *h3;
  
  h3 = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_NC",PID.c_str()));
  FD_RecoVsPID[Background::kNC] = (TH2D*)h3->Project3D("yx");
  FD_RecoVsPID[Background::kNC]->SetName("FD_RecoVsPID_NC");
  
  h3 = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_NuMuToNuMu",PID.c_str()));
  FD_RecoVsPID[Background::kNuMuCC] = (TH2D*)h3->Project3D("yx");
  FD_RecoVsPID[Background::kNuMuCC]->SetName("FD_RecoVsPID_NuMuCC");
  h3 = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_BNueToNuMu",PID.c_str()));
  FD_RecoVsPID[Background::kNuMuCC]->Add((TH2D*)h3->Project3D("yx"));
  
  h3 = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_BNueToBNue",PID.c_str()));
  FD_RecoVsPID[Background::kBNueCC] = (TH2D*)h3->Project3D("yx");
  FD_RecoVsPID[Background::kBNueCC]->SetName("FD_RecoVsPID_BNueCC");
  
  h3 = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_NuMuToNue",PID.c_str()));
  FD_RecoVsPID[Background::kNueCC] = (TH2D*)h3->Project3D("yx");
  FD_RecoVsPID[Background::kNueCC]->SetName("FD_RecoVsPID_NueCC");
  
  h3 = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_NuMuToNuTau",PID.c_str()));
  FD_RecoVsPID[Background::kNuTauCC] = (TH2D*)h3->Project3D("yx");
  FD_RecoVsPID[Background::kNuTauCC]->SetName("FD_RecoVsPID_NuTauCC");
  h3 = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_BNueToNuTau",PID.c_str()));
  FD_RecoVsPID[Background::kNuTauCC]->Add((TH2D*)h3->Project3D("yx"));
  
  if(f->Read("FDMC_Nu")!=0 && f->Read("FDMC_NuBar")!=0)//if these directories exists, then we can do separate nu and nubars
  {
    UseSeparateNuNuBar = true;
  }
  else
  {
    cout<<"Warning: FDMC_Nu and FDMC_NuBar directories not found, so nus and nubars will be combined."<<endl;
  }
  
  if(UseSeparateNuNuBar)
  {
    FD_TrueVsRecoVsPID[qNuMuToNue] = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_NuMuToNue",PID.c_str()));
    FD_TrueVsRecoVsPID[qNuMuToNue]->SetName(Form("FD_TrueVsRecoVs%s_NuMuToNue",PID.c_str()));
    
    FD_Nu_TrueVsRecoVsPID[qNuMuToNue] = (TH3D*)f->Get(Form("FDMC_Nu/TrueVsRecoVs%s_NuMuToNue",PID.c_str()));
    FD_Nu_TrueVsRecoVsPID[qNuMuToNue]->SetName(Form("FD_Nu_TrueVsRecoVs%s_NuMuToNue",PID.c_str()));
    
    FD_NuBar_TrueVsRecoVsPID[qNuMuToNue] = (TH3D*)f->Get(Form("FDMC_NuBar/TrueVsRecoVs%s_NuMuToNue",PID.c_str()));
    FD_NuBar_TrueVsRecoVsPID[qNuMuToNue]->SetName(Form("FD_NuBar_TrueVsRecoVs%s_NuMuToNue",PID.c_str()));
  }
  
  FD_TrueVsReco_Fid = (TH2D*)f->Get("FDMC/TrueVsReco_Fid");
  FD_TrueVsReco_Fid_NuMuCC = (TH2D*)f->Get("FDMC/TrueVsReco_Fid_NuMuCC");
  FD_TrueVsReco_CClike = (TH2D*)f->Get("FDMC/TrueVsReco_CClike");
  FD_TrueVsReco_CClike_NuMuCC = (TH2D*)f->Get("FDMC/TrueVsReco_CClike_NuMuCC");
  
  FD_TrueVsReco_Fid_NueCC = (TH2D*)f->Get("FDMC/TrueVsReco_Fid_NueCC");
  FD_TrueVsReco_Fid_NuTauCC = (TH2D*)f->Get("FDMC/TrueVsReco_Fid_NuTauCC");
  FD_True_Fid_NueCC = (TH1D*)f->Get("FDMC/True_Fid_NueCC");
  FD_True_Fid_NuTauCC = (TH1D*)f->Get("FDMC/True_Fid_NuTauCC");
  
  h3 = (TH3D*)f->Get(Form("NDMC/TrueVsRecoVs%s_NC",PID.c_str()));
  ND_RecoVsPID[Background::kNC] = (TH2D*)h3->Project3D("yx");
  ND_RecoVsPID[Background::kNC]->SetName("ND_RecoVsPID_NC");
  h3 = (TH3D*)f->Get(Form("NDMC/TrueVsRecoVs%s_NuMuToNuMu",PID.c_str()));
  ND_RecoVsPID[Background::kNuMuCC] = (TH2D*)h3->Project3D("yx");
  ND_RecoVsPID[Background::kNuMuCC]->SetName("ND_RecoVsPID_NuMuCC");
  h3 = (TH3D*)f->Get(Form("NDMC/TrueVsRecoVs%s_BNueToBNue",PID.c_str()));
  ND_RecoVsPID[Background::kBNueCC] = (TH2D*)h3->Project3D("yx");
  ND_RecoVsPID[Background::kBNueCC]->SetName("ND_RecoVsPID_BNueCC");
  
  ND_TrueVsReco_Fid = (TH2D*)f->Get("NDMC/TrueVsReco_Fid");
  ND_TrueVsReco_Fid_NuMuCC = (TH2D*)f->Get("NDMC/TrueVsReco_Fid_NuMuCC");
  ND_TrueVsReco_CClike = (TH2D*)f->Get("NDMC/TrueVsReco_CClike");
  ND_TrueVsReco_CClike_NuMuCC = (TH2D*)f->Get("NDMC/TrueVsReco_CClike_NuMuCC");
  
  NDData_Reco_CClike = (TH1D*)f->Get("NDData/Reco_CClike");
  
  FD_RecoVsPID[Background::kNC]->Scale(nPOTFar/fp);
  FD_RecoVsPID[Background::kNuMuCC]->Scale(nPOTFar/fp);
  FD_RecoVsPID[Background::kBNueCC]->Scale(nPOTFar/fp);
  FD_RecoVsPID[Background::kNueCC]->Scale(nPOTFar/fp);
  FD_RecoVsPID[Background::kNuTauCC]->Scale(nPOTFar/fp);
  
  FD_TrueVsReco_Fid->Scale(nPOTFar/fp);
  FD_TrueVsReco_Fid_NuMuCC->Scale(nPOTFar/fp);
  FD_TrueVsReco_CClike->Scale(nPOTFar/fp);
  FD_TrueVsReco_CClike_NuMuCC->Scale(nPOTFar/fp);
  
  FD_TrueVsReco_Fid_NueCC->Scale(nPOTFar/fp);
  FD_TrueVsReco_Fid_NuTauCC->Scale(nPOTFar/fp);
  FD_True_Fid_NueCC->Scale(nPOTFar/fp);
  FD_True_Fid_NuTauCC->Scale(nPOTFar/fp);
  
  ND_RecoVsPID[Background::kNC]->Scale(nPOTNear/np);
  ND_RecoVsPID[Background::kNuMuCC]->Scale(nPOTNear/np);
  ND_RecoVsPID[Background::kBNueCC]->Scale(nPOTNear/np);
  
  ND_TrueVsReco_Fid->Scale(nPOTNear/np);
  ND_TrueVsReco_Fid_NuMuCC->Scale(nPOTNear/np);
  ND_TrueVsReco_CClike->Scale(nPOTNear/np);
  ND_TrueVsReco_CClike_NuMuCC->Scale(nPOTNear/np);
  
  NDData_Reco_CClike->Scale(nPOTNear/np);
  
  double *r = new double[nReco+1];
  double *p = new double[nPID+1];
  double *t = new double[nTrue+1];
  double *rcc = new double[nRecoCC+1];
  int i;
  for(i=0;i<nReco+1;i++)
  {
    r[i] = RecoEdges.at(i);
  }
  for(i=0;i<nPID+1;i++)
  {
    p[i] = PIDEdges.at(i);
  }
  for(i=0;i<nTrue+1;i++)
  {
    t[i] = TrueEdges.at(i);
  }
  for(i=0;i<nRecoCC+1;i++)
  {
    rcc[i] = RecoCCEdges.at(i);
  }
  
  if(FD_RecoVsPID[Background::kNC]->GetNbinsX()!=nPID)
  {
    RebinP = true;
  }
  if(FD_RecoVsPID[Background::kNC]->GetNbinsY()!=nReco)
  {
    RebinE = true;
  }
  
  if(RebinE || RebinP)
  {
    RebinInputHists();
  }
  
  int ip,ir,it;
  double temp;
  
  TH1D *h1;
  
  FNRatio[Background::kNuMuCC] = (TH2D*)FD_RecoVsPID[Background::kNuMuCC]->Clone("FNRatio_NuMuCC");
  FNRatio[Background::kNC] = (TH2D*)FD_RecoVsPID[Background::kNC]->Clone("FNRatio_NC");
  FNRatio[Background::kBNueCC] = (TH2D*)FD_RecoVsPID[Background::kBNueCC]->Clone("FNRatio_BNueCC");
  
  FNRatio[Background::kNuMuCC]->Divide(ND_RecoVsPID[Background::kNuMuCC]);
  FNRatio[Background::kNC]->Divide(ND_RecoVsPID[Background::kNC]);
  FNRatio[Background::kBNueCC]->Divide(ND_RecoVsPID[Background::kBNueCC]);
  
  ND_Reco_CClike = ND_TrueVsReco_CClike->ProjectionX("ND_Reco_CClike");
  
  FD_Reco_CClike = FD_TrueVsReco_CClike->ProjectionX("FD_Reco_CClike");
  FD_Reco2True_CClike = new TH2D("FD_Reco2True_CClike","",nRecoCC,rcc,nTrue,t);
  for(ir=0;ir<nRecoCC;ir++)
  {
    for(it=0;it<nTrue;it++)
    {
      temp=0;
      if(FD_Reco_CClike->GetBinContent(ir+1)>0)
      {
        temp = FD_TrueVsReco_CClike->GetBinContent(ir+1,it+1)/FD_Reco_CClike->GetBinContent(ir+1);
      }
      FD_Reco2True_CClike->SetBinContent(ir+1,it+1,temp);
    }
  }
  
  FD_Purity_CC = new TH1D("FD_Purity_CC","",nTrue,t);
  FD_Eff_CC = new TH1D("FD_Eff_CC","",nTrue,t);
  TH1D *FD_True_CClike = FD_TrueVsReco_CClike->ProjectionY("FD_True_CClike");
  TH1D *FD_True_CClike_NuMuCC = FD_TrueVsReco_CClike_NuMuCC->ProjectionY("FD_True_CClike_NuMuCC");
  TH1D *FD_True_Fid_NuMuCC = FD_TrueVsReco_Fid_NuMuCC->ProjectionY("FD_True_Fid_NuMuCC");
  FD_Purity_CC->Divide(FD_True_CClike_NuMuCC,FD_True_CClike,1,1);
  FD_Eff_CC->Divide(FD_True_CClike_NuMuCC,FD_True_Fid_NuMuCC,1,1);
  
  FD_True2Reco_Fid[Background::kNueCC] = new TH2D("FD_True2Reco_Fid_NueCC","",nReco,r,nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      temp=0;
      if(FD_True_Fid_NueCC->GetBinContent(it+1)>0)
      {
        temp = FD_TrueVsReco_Fid_NueCC->GetBinContent(ir+1,it+1)/FD_True_Fid_NueCC->GetBinContent(it+1);
      }
      FD_True2Reco_Fid[Background::kNueCC]->SetBinContent(ir+1,it+1,temp);
    }
  }
  
  FD_True2Reco_Fid[Background::kNuTauCC] = new TH2D("FD_True2Reco_Fid_NuTauCC","",nReco,r,nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      temp=0;
      if(FD_True_Fid_NuTauCC->GetBinContent(it+1)>0)
      {
        temp = FD_TrueVsReco_Fid_NuTauCC->GetBinContent(ir+1,it+1)/FD_True_Fid_NuTauCC->GetBinContent(it+1);
      }
      FD_True2Reco_Fid[Background::kNuTauCC]->SetBinContent(ir+1,it+1,temp);
    }
  }
  
  FD_Eff[Background::kNueCC] = new TH2D("FD_Eff_NueCC","",nPID,p,nReco,r);
  h1 = FD_TrueVsReco_Fid_NueCC->ProjectionX();
  for(ip=0;ip<nPID;ip++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      temp=0;
      if(h1->GetBinContent(ir+1)>0)
      {
        temp = FD_RecoVsPID[Background::kNueCC]->GetBinContent(ip+1,ir+1)/h1->GetBinContent(ir+1);
      }
      FD_Eff[Background::kNueCC]->SetBinContent(ip+1,ir+1,temp);
    }
  }
  
  FD_Eff[Background::kNuTauCC] = new TH2D("FD_Eff_NuTauCC","",nPID,p,nReco,r);
  h1 = FD_TrueVsReco_Fid_NuTauCC->ProjectionX();
  for(ip=0;ip<nPID;ip++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      temp=0;
      if(h1->GetBinContent(ir+1)>0)
      {
        temp = FD_RecoVsPID[Background::kNuTauCC]->GetBinContent(ip+1,ir+1)/h1->GetBinContent(ir+1);
      }
      FD_Eff[Background::kNuTauCC]->SetBinContent(ip+1,ir+1,temp);
    }
  }
  
  delete [] r;
  delete [] p;
  delete [] t;
  
  return;
}
void Extrapolate2D_Simple::ReadFNFileNo3D()
{
  if(ReadError) return;
  
  FD_RecoVsPID.clear();
  FD_True2Reco_Fid.clear();
  FD_Eff.clear();
  FNRatio.clear();
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(FN_infile.c_str())))
  {
    cout<<"Failure to read F/N file."<<endl;
    ReadError = true;
    return;
  }
  
  TFile *f = new TFile(FN_infile.c_str(),"READ");
  
  double np,fp;
  TTree *tree = (TTree*)f->Get("paramtree");
  tree->SetBranchAddress("nearPOT",&np);
  tree->SetBranchAddress("farPOT",&fp);
  tree->SetBranchAddress("Theta12",&Theta12);
  tree->SetBranchAddress("Theta13",&Theta13);
  tree->SetBranchAddress("Theta23",&Theta23);
  tree->SetBranchAddress("DeltaMSq23",&DeltaMSq23);
  tree->SetBranchAddress("DeltaMSq12",&DeltaMSq12);
  tree->SetBranchAddress("DeltaCP",&DeltaCP);
  tree->GetEntry(0);
  
  osc.SetOscParam(OscPar::kTh12,Theta12);
  osc.SetOscParam(OscPar::kTh13,Theta13);
  osc.SetOscParam(OscPar::kTh23,Theta23);
  osc.SetOscParam(OscPar::kDeltaM23,DeltaMSq23);
  osc.SetOscParam(OscPar::kDeltaM12,DeltaMSq12);
  osc.SetOscParam(OscPar::kDelta,DeltaCP);
  SetNominalOscProb();
  
  FD_RecoVsPID[Background::kNC] = (TH2D*)f->Get(Form("FDMC/RecoVs%s_NC",PID.c_str()));
  FD_RecoVsPID[Background::kNC]->SetName("FD_RecoVsPID_NC");
  
  FD_RecoVsPID[Background::kNuMuCC] = (TH2D*)f->Get(Form("FDMC/RecoVs%s_NuMuCC",PID.c_str()));
  FD_RecoVsPID[Background::kNuMuCC]->SetName("FD_RecoVsPID_NuMuCC");
  
  FD_RecoVsPID[Background::kBNueCC] = (TH2D*)f->Get(Form("FDMC/RecoVs%s_BNueCC",PID.c_str()));
  FD_RecoVsPID[Background::kBNueCC]->SetName("FD_RecoVsPID_BNueCC");
  
  FD_RecoVsPID[Background::kNueCC] = (TH2D*)f->Get(Form("FDMC/RecoVs%s_NueCC",PID.c_str()));
  FD_RecoVsPID[Background::kNueCC]->SetName("FD_RecoVsPID_NueCC");
  
  FD_RecoVsPID[Background::kNuTauCC] = (TH2D*)f->Get(Form("FDMC/RecoVs%s_NuTauCC",PID.c_str()));
  FD_RecoVsPID[Background::kNuTauCC]->SetName("FD_RecoVsPID_NuTauCC");
  
  FD_TrueVsReco_Fid = (TH2D*)f->Get("FDMC/TrueVsReco_Fid");
  FD_TrueVsReco_Fid_NuMuCC = (TH2D*)f->Get("FDMC/TrueVsReco_Fid_NuMuCC");
  FD_TrueVsReco_CClike = (TH2D*)f->Get("FDMC/TrueVsReco_CClike");
  FD_TrueVsReco_CClike_NuMuCC = (TH2D*)f->Get("FDMC/TrueVsReco_CClike_NuMuCC");
  
  FD_TrueVsReco_Fid_NueCC = (TH2D*)f->Get("FDMC/TrueVsReco_Fid_NueCC");
  FD_TrueVsReco_Fid_NuTauCC = (TH2D*)f->Get("FDMC/TrueVsReco_Fid_NuTauCC");
  FD_True_Fid_NueCC = (TH1D*)f->Get("FDMC/True_Fid_NueCC");
  FD_True_Fid_NuTauCC = (TH1D*)f->Get("FDMC/True_Fid_NuTauCC");
  
  ND_RecoVsPID[Background::kNC] = (TH2D*)f->Get(Form("NDMC/RecoVs%s_NC",PID.c_str()));
  ND_RecoVsPID[Background::kNC]->SetName("ND_RecoVsPID_NC");
  ND_RecoVsPID[Background::kNuMuCC] = (TH2D*)f->Get(Form("NDMC/RecoVs%s_NuMuCC",PID.c_str()));
  ND_RecoVsPID[Background::kNuMuCC]->SetName("ND_RecoVsPID_NuMuCC");
  ND_RecoVsPID[Background::kBNueCC] = (TH2D*)f->Get(Form("NDMC/RecoVs%s_BNueCC",PID.c_str()));
  ND_RecoVsPID[Background::kBNueCC]->SetName("ND_RecoVsPID_BNueCC");
  
  ND_TrueVsReco_Fid = (TH2D*)f->Get("NDMC/TrueVsReco_Fid");
  ND_TrueVsReco_Fid_NuMuCC = (TH2D*)f->Get("NDMC/TrueVsReco_Fid_NuMuCC");
  ND_TrueVsReco_CClike = (TH2D*)f->Get("NDMC/TrueVsReco_CClike");
  ND_TrueVsReco_CClike_NuMuCC = (TH2D*)f->Get("NDMC/TrueVsReco_CClike_NuMuCC");
  
  NDData_Reco_CClike = (TH1D*)f->Get("NDData/Reco_CClike");
  
  FD_RecoVsPID[Background::kNC]->Scale(nPOTFar/fp);
  FD_RecoVsPID[Background::kNuMuCC]->Scale(nPOTFar/fp);
  FD_RecoVsPID[Background::kBNueCC]->Scale(nPOTFar/fp);
  FD_RecoVsPID[Background::kNueCC]->Scale(nPOTFar/fp);
  FD_RecoVsPID[Background::kNuTauCC]->Scale(nPOTFar/fp);
  
  FD_TrueVsReco_Fid->Scale(nPOTFar/fp);
  FD_TrueVsReco_Fid_NuMuCC->Scale(nPOTFar/fp);
  FD_TrueVsReco_CClike->Scale(nPOTFar/fp);
  FD_TrueVsReco_CClike_NuMuCC->Scale(nPOTFar/fp);
  
  FD_TrueVsReco_Fid_NueCC->Scale(nPOTFar/fp);
  FD_TrueVsReco_Fid_NuTauCC->Scale(nPOTFar/fp);
  FD_True_Fid_NueCC->Scale(nPOTFar/fp);
  FD_True_Fid_NuTauCC->Scale(nPOTFar/fp);
  
  ND_RecoVsPID[Background::kNC]->Scale(nPOTNear/np);
  ND_RecoVsPID[Background::kNuMuCC]->Scale(nPOTNear/np);
  ND_RecoVsPID[Background::kBNueCC]->Scale(nPOTNear/np);
  
  ND_TrueVsReco_Fid->Scale(nPOTNear/np);
  ND_TrueVsReco_Fid_NuMuCC->Scale(nPOTNear/np);
  ND_TrueVsReco_CClike->Scale(nPOTNear/np);
  ND_TrueVsReco_CClike_NuMuCC->Scale(nPOTNear/np);
  
  NDData_Reco_CClike->Scale(nPOTNear/np);
  
  double *r = new double[nReco+1];
  double *p = new double[nPID+1];
  double *t = new double[nTrue+1];
  double *rcc = new double[nRecoCC+1];
  int i;
  for(i=0;i<nReco+1;i++)
  {
    r[i] = RecoEdges.at(i);
  }
  for(i=0;i<nPID+1;i++)
  {
    p[i] = PIDEdges.at(i);
  }
  for(i=0;i<nTrue+1;i++)
  {
    t[i] = TrueEdges.at(i);
  }
  for(i=0;i<nRecoCC+1;i++)
  {
    rcc[i] = RecoCCEdges.at(i);
  }
  
  if(FD_RecoVsPID[Background::kNC]->GetNbinsX()!=nPID)
  {
    RebinP = true;
  }
  if(FD_RecoVsPID[Background::kNC]->GetNbinsY()!=nReco)
  {
    RebinE = true;
  }
  
  if(RebinE || RebinP)
  {
    RebinInputHists();
  }
  
  int ip,ir,it;
  double temp;
  
  TH1D *h1;
  
  FNRatio[Background::kNuMuCC] = (TH2D*)FD_RecoVsPID[Background::kNuMuCC]->Clone("FNRatio_NuMuCC");
  FNRatio[Background::kNC] = (TH2D*)FD_RecoVsPID[Background::kNC]->Clone("FNRatio_NC");
  FNRatio[Background::kBNueCC] = (TH2D*)FD_RecoVsPID[Background::kBNueCC]->Clone("FNRatio_BNueCC");
  
  FNRatio[Background::kNuMuCC]->Divide(ND_RecoVsPID[Background::kNuMuCC]);
  FNRatio[Background::kNC]->Divide(ND_RecoVsPID[Background::kNC]);
  FNRatio[Background::kBNueCC]->Divide(ND_RecoVsPID[Background::kBNueCC]);
  
  ND_Reco_CClike = ND_TrueVsReco_CClike->ProjectionX("ND_Reco_CClike");
  
  FD_Reco_CClike = FD_TrueVsReco_CClike->ProjectionX("FD_Reco_CClike");
  FD_Reco2True_CClike = new TH2D("FD_Reco2True_CClike","",nRecoCC,rcc,nTrue,t);
  for(ir=0;ir<nRecoCC;ir++)
  {
    for(it=0;it<nTrue;it++)
    {
      temp=0;
      if(FD_Reco_CClike->GetBinContent(ir+1)>0)
      {
        temp = FD_TrueVsReco_CClike->GetBinContent(ir+1,it+1)/FD_Reco_CClike->GetBinContent(ir+1);
      }
      FD_Reco2True_CClike->SetBinContent(ir+1,it+1,temp);
    }
  }
  
  FD_Purity_CC = new TH1D("FD_Purity_CC","",nTrue,t);
  FD_Eff_CC = new TH1D("FD_Eff_CC","",nTrue,t);
  TH1D *FD_True_CClike = FD_TrueVsReco_CClike->ProjectionY("FD_True_CClike");
  TH1D *FD_True_CClike_NuMuCC = FD_TrueVsReco_CClike_NuMuCC->ProjectionY("FD_True_CClike_NuMuCC");
  TH1D *FD_True_Fid_NuMuCC = FD_TrueVsReco_Fid_NuMuCC->ProjectionY("FD_True_Fid_NuMuCC");
  FD_Purity_CC->Divide(FD_True_CClike_NuMuCC,FD_True_CClike,1,1);
  FD_Eff_CC->Divide(FD_True_CClike_NuMuCC,FD_True_Fid_NuMuCC,1,1);
  
  FD_True2Reco_Fid[Background::kNueCC] = new TH2D("FD_True2Reco_Fid_NueCC","",nReco,r,nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      temp=0;
      if(FD_True_Fid_NueCC->GetBinContent(it+1)>0)
      {
        temp = FD_TrueVsReco_Fid_NueCC->GetBinContent(ir+1,it+1)/FD_True_Fid_NueCC->GetBinContent(it+1);
      }
      FD_True2Reco_Fid[Background::kNueCC]->SetBinContent(ir+1,it+1,temp);
    }
  }
  
  FD_True2Reco_Fid[Background::kNuTauCC] = new TH2D("FD_True2Reco_Fid_NuTauCC","",nReco,r,nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      temp=0;
      if(FD_True_Fid_NuTauCC->GetBinContent(it+1)>0)
      {
        temp = FD_TrueVsReco_Fid_NuTauCC->GetBinContent(ir+1,it+1)/FD_True_Fid_NuTauCC->GetBinContent(it+1);
      }
      FD_True2Reco_Fid[Background::kNuTauCC]->SetBinContent(ir+1,it+1,temp);
    }
  }
  
  FD_Eff[Background::kNueCC] = new TH2D("FD_Eff_NueCC","",nPID,p,nReco,r);
  h1 = FD_TrueVsReco_Fid_NueCC->ProjectionX();
  for(ip=0;ip<nPID;ip++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      temp=0;
      if(h1->GetBinContent(ir+1)>0)
      {
        temp = FD_RecoVsPID[Background::kNueCC]->GetBinContent(ip+1,ir+1)/h1->GetBinContent(ir+1);
      }
      FD_Eff[Background::kNueCC]->SetBinContent(ip+1,ir+1,temp);
    }
  }
  
  FD_Eff[Background::kNuTauCC] = new TH2D("FD_Eff_NuTauCC","",nPID,p,nReco,r);
  h1 = FD_TrueVsReco_Fid_NuTauCC->ProjectionX();
  for(ip=0;ip<nPID;ip++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      temp=0;
      if(h1->GetBinContent(ir+1)>0)
      {
        temp = FD_RecoVsPID[Background::kNuTauCC]->GetBinContent(ip+1,ir+1)/h1->GetBinContent(ir+1);
      }
      FD_Eff[Background::kNuTauCC]->SetBinContent(ip+1,ir+1,temp);
    }
  }
  
  delete [] r;
  delete [] p;
  delete [] t;
  
  return;
}
void Extrapolate2D_Simple::ReadFiles()
{
  ReadNDDataFile();
  if(UseInputWithNo3D)
  {
    ReadFNFileNo3D();
  }
  else
  {
    ReadFNFile();
  }
  ReadXSecFile();
  return;
}
void Extrapolate2D_Simple::FarNearPred(Background::Background_t bg)
{
  if(bg!=Background::kNuMuCC && bg!=Background::kNC && bg!=Background::kBNueCC)
  {
    cout<<"Can only use FarNearPred() for cc numu, nc and beam nues"<<endl;
    return;
  }
  
  Pred[bg]->Multiply(FNRatio[bg],NDData[bg],1,1);  
  
  return;
}
void Extrapolate2D_Simple::AppearancePred(Background::Background_t bg)
{
  if(bg!=Background::kNueCC && bg!=Background::kNuTauCC)
  {
    cout<<"Can only use AppearancePred() for taus and signal nues"<<endl;
    return;
  }
  
  int ip,ir,it;
  double temp,sum;
  double nu,nubar;
  for(ip=0;ip<nPID;ip++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      sum=0;
      
      for(it=0;it<nTrue;it++)
      {
	temp = Pred_CC_Fid->GetBinContent(it+1)*XSecWeight[bg]->GetBinContent(it+1)*FD_True2Reco_Fid[bg]->GetBinContent(ir+1,it+1)*FD_Eff[bg]->GetBinContent(ip+1,ir+1);
        nu=temp;
	nubar=temp;
	
	if(bg==Background::kNueCC)
	{
	  temp*=NominalOscProb[qNuMuToNue]->GetBinContent(it+1);
          Pred_3D[qNuMuToNue]->SetBinContent(ip+1,ir+1,it+1,temp);
	  
	  if(UseSeparateNuNuBar)
	  {
            if(FD_TrueVsRecoVsPID[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)>0)
            {
              nu*=(NominalOscProb[qNuMuToNue]->GetBinContent(it+1)*FD_Nu_TrueVsRecoVsPID[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)/FD_TrueVsRecoVsPID[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1));
            }
            else
            {
              nu=0;
            }
	    Pred_Nu_3D[qNuMuToNue]->SetBinContent(ip+1,ir+1,it+1,nu);
	    
            if(FD_TrueVsRecoVsPID[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)>0)
            {
              nubar*=(NominalOscProb_NuBar[qNuMuToNue]->GetBinContent(it+1)*FD_NuBar_TrueVsRecoVsPID[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)/FD_TrueVsRecoVsPID[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1));
            }
            else
            {
              nubar=0;
            }
	    Pred_NuBar_3D[qNuMuToNue]->SetBinContent(ip+1,ir+1,it+1,nubar);
	  }
	}
	if(bg==Background::kNuTauCC)
	{
	  temp*=NominalOscProb[qNuMuToNuTau]->GetBinContent(it+1);
          Pred_3D[qNuMuToNuTau]->SetBinContent(ip+1,ir+1,it+1,temp);
	}
	
	if(UseSeparateNuNuBar && bg==Background::kNueCC)
	{
	  sum += (nu+nubar);
	}
	else 
	{
	  sum+=temp;
	}
      }
      Pred[bg]->SetBinContent(ip+1,ir+1,sum);
    }
  }
  
  
  return;
}
void Extrapolate2D_Simple::BeamNuePred()
{
  Background::Background_t bg = Background::kBNueCC;
  Pred[bg]->Add(FD_RecoVsPID[bg]);
  
  return;
}
void Extrapolate2D_Simple::OscillatePrediction()
{
  if(!Init)
  {
    cout<<"Can't call OscillatePrediction(): The prediction wasn't properly initialized with GetPrediction()."<<endl;
    return;
  }
  
  if(UseSeparateNuNuBar)
  {
    OscillatePrediction_SepNuNuBar();
  }
  else
  {
    OscillatePrediction_NoSep();
  }
  
  if(PrintResult)
  {
    cout.precision(2);
    cout<<fixed;
    cout<<"sin^2(2Th13) = "<<TMath::Sin(2.*osc.GetOscParam(OscPar::kTh13))*TMath::Sin(2.*osc.GetOscParam(OscPar::kTh13))<<endl;
    cout<<"NC: "<<Pred[Background::kNC]->Integral()<<endl;
    cout<<"NuMuCC: "<<Pred[Background::kNuMuCC]->Integral()<<endl;
    cout<<"BNueCC: "<<Pred[Background::kBNueCC]->Integral()<<endl;
    cout<<"NuTauCC: "<<Pred[Background::kNuTauCC]->Integral()<<endl;
    cout<<"Total: "<<Pred_TotalBkgd->Integral()<<endl;
    cout<<"NueCC: "<<Pred[Background::kNueCC]->Integral()<<endl;
  }
  
  Oscillated=true;
  
  return;
}
void Extrapolate2D_Simple::OscillatePrediction_SepNuNuBar()
{
  Background::Background_t bg = Background::kNueCC;
  Pred[bg]->Reset("ICE");
  
  int ip,ir,it;
  double temp,sum,E;
  for(ip=0;ip<nPID;ip++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      sum=0;
      
      for(it=0;it<nTrue;it++)
      {
	E = Pred_CC_Fid->GetXaxis()->GetBinCenter(it+1);
	
	//nu
        temp=0;
	osc.SetOscParam(OscPar::kNuAntiNu,1);
	if(NominalOscProb[qNuMuToNue]->GetBinContent(it+1)>1e-20)
	{
	  temp = Pred_Nu_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.MuToElec(E)/NominalOscProb[qNuMuToNue]->GetBinContent(it+1);
	}
	sum+=temp;
	
	//nubar
	temp=0;
	osc.SetOscParam(OscPar::kNuAntiNu,-1);
	if(NominalOscProb_NuBar[qNuMuToNue]->GetBinContent(it+1)>1e-20)
	{
	  temp = Pred_NuBar_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.MuToElec(E)/NominalOscProb_NuBar[qNuMuToNue]->GetBinContent(it+1);
	}
        sum+=temp;
      }
      Pred[bg]->SetBinContent(ip+1,ir+1,sum);
    }
  }  
  
  Set1DPredHists();
  
  return;
}
void Extrapolate2D_Simple::OscillatePrediction_NoSep()
{
  Background::Background_t bg = Background::kNueCC;
  Pred[bg]->Reset("ICE");
  
  osc.SetOscParam(OscPar::kNuAntiNu,1);//assuming neutrino probability
  
  int ip,ir,it;
  double temp,sum,E;
  for(ip=0;ip<nPID;ip++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      sum=0;
      
      for(it=0;it<nTrue;it++)
      {
        E = Pred_CC_Fid->GetXaxis()->GetBinCenter(it+1);
        
        temp = Pred_CC_Fid->GetBinContent(it+1)*XSecWeight[bg]->GetBinContent(it+1)*FD_True2Reco_Fid[bg]->GetBinContent(ir+1,it+1)*FD_Eff[bg]->GetBinContent(ip+1,ir+1);
        
        temp*=osc.MuToElec(E);
        
        sum+=temp;
      }
      Pred[bg]->SetBinContent(ip+1,ir+1,sum);
    }
  }  
  
  Set1DPredHists();
  
  return;
}
void Extrapolate2D_Simple::SetupPredHists()
{
  Pred_3D.clear();
  Pred_Nu_3D.clear();
  Pred_NuBar_3D.clear();
  Pred.clear();
  Pred_TotalBkgd = NULL;
  Pred_CC_Fid = NULL;
  Pred_TotalBkgd_VsBinNumber = NULL;
  Pred_Signal_VsBinNumber = NULL;
  
  double *r = new double[nReco+1];
  double *p = new double[nPID+1];
  double *t = new double[nTrue+1];
  int i;
  for(i=0;i<nReco+1;i++)
  {
    r[i] = RecoEdges.at(i);
  }
  for(i=0;i<nPID+1;i++)
  {
    p[i] = PIDEdges.at(i);
  }
  for(i=0;i<nTrue+1;i++)
  {
    t[i] = TrueEdges.at(i);
  }
  
  string name;
  
  //include these because there is enough info to calculate them (even without the 3D histograms) and it is useful in ErrorCalc
  Pred_3D[qNuMuToNue] = new TH3D("Pred_3D_NuMuToNue","",nPID,p,nReco,r,nTrue,t);
  Pred_3D[qNuMuToNuTau] = new TH3D("Pred_3D_NuMuToNuTau","",nPID,p,nReco,r,nTrue,t);
  
  Pred_Nu_3D[qNuMuToNue] = new TH3D("Pred_Nu_3D_NuMuToNue","",nPID,p,nReco,r,nTrue,t);
  Pred_Nu_3D[qNuMuToNuTau] = new TH3D("Pred_Nu_3D_NuMuToNuTau","",nPID,p,nReco,r,nTrue,t);
  
  Pred_NuBar_3D[qNuMuToNue] = new TH3D("Pred_NuBar_3D_NuMuToNue","",nPID,p,nReco,r,nTrue,t);
  Pred_NuBar_3D[qNuMuToNuTau] = new TH3D("Pred_NuBar_3D_NuMuToNuTau","",nPID,p,nReco,r,nTrue,t);
  
  TH2D *h2;
  for(unsigned int j=0;j<FDComponents.size();j++)
  {
    name = "Pred_" + string(Background::AsString(FDComponents[j]));
    h2 = new TH2D(name.c_str(),"",nPID,p,nReco,r);
    Pred[FDComponents[j]] = h2;
  }
  
  Pred_TotalBkgd = new TH2D("Pred_TotalBkgd","",nPID,p,nReco,r);
  
  Pred_CC_Fid = new TH1D("Pred_CC_Fid","",nTrue,t);
  
  int totbin = nReco*nPID;
  Pred_TotalBkgd_VsBinNumber = new TH1D("Pred_TotalBkgd_VsBinNumber","",totbin,-0.5,totbin-0.5);
  Pred_Signal_VsBinNumber = new TH1D("Pred_Signal_VsBinNumber","",totbin,-0.5,totbin-0.5);
  
  delete [] r;
  delete [] p;
  delete [] t;
  
  return;
}
void Extrapolate2D_Simple::RebinInputHists()
{
  double *r = new double[nReco+1];
  double *p = new double[nPID+1];
  double *t = new double[nTrue+1];
  int i;
  for(i=0;i<nReco+1;i++)
  {
    r[i] = RecoEdges.at(i);
  }
  for(i=0;i<nPID+1;i++)
  {
    p[i] = PIDEdges.at(i);
  }
  for(i=0;i<nTrue+1;i++)
  {
    t[i] = TrueEdges.at(i);
  }
  
  unsigned int j;
  Background::Background_t bg;
  for(j=0;j<FDComponents.size();j++)
  {
    bg = FDComponents.at(j);
    MBH.Rebin2DHist(FD_RecoVsPID[bg],nPID,p,nReco,r);
    if(bg!=Background::kNueCC && bg!=Background::kNuTauCC)
    {
      MBH.Rebin2DHist(ND_RecoVsPID[bg],nPID,p,nReco,r);
    }
  }
  MBH.Rebin2DHist(FD_TrueVsReco_Fid_NueCC,nReco,r,0,0);
  MBH.Rebin2DHist(FD_TrueVsReco_Fid_NuTauCC,nReco,r,0,0);
  
  if(UseSeparateNuNuBar)
  {
    MBH.Rebin3DHist(FD_TrueVsRecoVsPID[qNuMuToNue],nPID,p,nReco,r,0,0);
    MBH.Rebin3DHist(FD_Nu_TrueVsRecoVsPID[qNuMuToNue],nPID,p,nReco,r,0,0);
    MBH.Rebin3DHist(FD_NuBar_TrueVsRecoVsPID[qNuMuToNue],nPID,p,nReco,r,0,0);
  }
  
  return;
}
