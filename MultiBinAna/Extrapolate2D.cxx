#define Extrapolate2D_C

#include "NueAna/MultiBinAna/Extrapolate2D.h"

Extrapolate2D::Extrapolate2D() {
  nReco = 100;
  nRecoCC = 100;
  // RecoEnergy bin edges
  for (int i = 0; i < nReco + 1; i++) {
    RecoEdges.push_back(i * 1.);
    RecoCCEdges.push_back(i * 1.);
  }

  nTrue = 1200;
  for (int i = 0; i < nTrue + 1; i++) {
    TrueEdges.push_back(i * 0.1);
  }
  
  nPID = 20;
  for (int i = 0; i < nPID + 1; i++) {
    PIDEdges.push_back(i * 0.05);
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
//   FDComponents.push_back(Background::kNuMuNC);
//   FDComponents.push_back(Background::kBNueNC);
  
  Theta13=0;
  Theta12=0;
  Theta23=0;
  DeltaMSq23=0;
  DeltaMSq12=0;
  DeltaCP=0;
  //
  Eps_ee=0;
  Eps_emu=0;
  Eps_etau=0;
  Eps_mumu=0;
  Eps_mutau=0;
  Eps_tautau=0;
  Delta_emu=0;
  Delta_etau=0;
  Delta_mutau=0;
  //
  Theta14=0;
  Theta24=0;
  Theta34=0;
  Dm41=0;
  Delta14=0;
  Delta24=0;

  ReadExtrapFromFile = false;
  
  MRE_infile = "NULL";
  ReadMREFromFile = false;
  
  RebinE = false;
  RebinP = false;
  
  PrintResult = false;
  
  WriteOutput = false;
  
  Init = false;//ie haven't called GetPrediction() yet
  
  ReadError = false;
  
  Oscillated = false;//haven't called OscillatePrediction() yet
  
  UseSeparateNuNuBar=false;
  
  return;
}

Extrapolate2D::~Extrapolate2D() {
  if(WriteOutput) ExtrapFile->Close();
}

void Extrapolate2D::SetRecoBins(Int_t n, Double_t *e) {
  nReco = n;
  RecoEdges.clear();
  for(int i = 0; i < nReco + 1; i++) {
    RecoEdges.push_back(e[i]);
  }
  
  return;
}

void Extrapolate2D::SetPIDBins(Int_t n, Double_t *e) {
  int i;
  nPID = n;
  PIDEdges.clear();
  for(i=0;i<nPID+1;i++) {
    PIDEdges.push_back(e[i]);
  }
  
  return;
}

void Extrapolate2D::SetTrueBins(Int_t n, Double_t *e) {
  int i;
  nTrue = n;
  TrueEdges.clear();
  for(i=0;i<nTrue+1;i++) {
    TrueEdges.push_back(e[i]);
  }
  
  return;
}

void Extrapolate2D::SetOutputFile(string filename) {
  outFileName = filename;
  WriteOutput = true;
  ExtrapFile = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  return;
}

void Extrapolate2D::SetPID(string pid) {
  PID = pid;
  PID_NDData = pid;
  return;
}

void Extrapolate2D::SetNDDataPID(string pid) {
  PID_NDData = pid;
  return;
}

void Extrapolate2D::SetNDDataFile(string s) {
  NDData_infile = s;
  return;
}

void Extrapolate2D::SetFNFile(string s) {
  FN_infile = s;
  return;
}

void Extrapolate2D::SetMREFile(string s) {
  MRE_infile = s;
  ReadMREFromFile = true;
  return;
}

void Extrapolate2D::SetXSecFile(string s) {
  XSec_infile = s;
  return;
}

void Extrapolate2D::SetReadExtrapFromFile(string s) {
  inFileName = s;
  ReadExtrapFromFile = true;
  return;
}

void Extrapolate2D::InitializeOscCalc() {
  Double_t dm2_12 = 7.59e-5;  // updated 2014
  Double_t dm2_23 = 2.41e-3;
  Double_t par[24] = {0};
  par[OscPar::kL] = 735.0;
  par[OscPar::kTh23] = 3.1415926/4.0;
  par[OscPar::kTh12] = 0.59365;    // Sin2(2Th12) = 0.86
  par[OscPar::kTh13] = 0.147;  // Daya Bay Neutrino 2014
  par[OscPar::kDeltaM23] = dm2_23; //normal heirarchy
  par[OscPar::kDeltaM12] = dm2_12;
  par[OscPar::kDensity] = 2.75; //standard rock density
  par[OscPar::kDelta] = 0;
  par[OscPar::kNuAntiNu] = 1;
  
  // Initialize NSI parameters:
  par[OscPar::kEps_ee] =  0;
  par[OscPar::kEps_emu] =  0;
  par[OscPar::kEps_etau] =  0;
  par[OscPar::kEps_mumu] =  0;
  par[OscPar::kEps_mutau] = 0; 
  par[OscPar::kEps_tautau] = 0; 
  par[OscPar::kDelta_emu] = 0;
  par[OscPar::kDelta_etau] = 0;
  par[OscPar::kDelta_mutau] = 0;
  
  par[OscPar::kTh14] = 0;
  par[OscPar::kTh24] = 0;
  par[OscPar::kTh34] = 0;
  par[OscPar::kDm41] = 0;
  par[OscPar::kDelta14] = 0;
  par[OscPar::kDelta24] = 0;
  osc.SetOscParam(par);
  
  return;
}

void Extrapolate2D::SetOscPar(OscPar::OscPar_t par, double val) {
  osc.SetOscParam(par,val);
  
  return;
}

void Extrapolate2D::SetSinSq2Th13(double val) {
  Double_t th13 = TMath::ASin(TMath::Sqrt(val))/2;
  osc.SetOscParam(OscPar::kTh13,th13);
  
  return;
}

void Extrapolate2D::SetSinSq2Th14(double val) {
  Double_t th14 = TMath::ASin(TMath::Sqrt(val))/2;
  osc.SetOscParam(OscPar::kTh14,th14);
  
  return;
}

void Extrapolate2D::SetSinSqTh14(double val) {
  Double_t th14 = TMath::ASin(TMath::Sqrt(val));
  osc.SetOscParam(OscPar::kTh14,th14);
  
  return;
}

void Extrapolate2D::SetSinSqTh24(double val) {
  Double_t th24 = TMath::ASin(TMath::Sqrt(val));
  osc.SetOscParam(OscPar::kTh24,th24);
  
  return;
}

void Extrapolate2D::SetDeltaCP(double val) {
  osc.SetOscParam(OscPar::kDelta,val);
  
  return;
}

void Extrapolate2D::SetEps_etau(double val) {
  osc.SetOscParam(OscPar::kEps_etau,val);
  
  return;
}

void Extrapolate2D::SetDm41(double val) {
  osc.SetOscParam(OscPar::kDm41,val);
  
  return;
}

void Extrapolate2D::SetDelta_etau(double val) {
  osc.SetOscParam(OscPar::kDelta_etau,val);
  
  return;
}

void Extrapolate2D::InvertMassHierarchy() {
  double temp = -1.*osc.GetOscParam(OscPar::kDeltaM23);
  osc.SetOscParam(OscPar::kDeltaM23,temp);
  return;
}

void Extrapolate2D::SetNominalOscProb() {
  NominalOscProb.clear();
  
  int it;
  double *t = new double[nTrue+1];
  for(it=0;it<nTrue+1;it++) {
    t[it] = TrueEdges.at(it);
  }
  
  TH1D *h;
  double E;
  
  //osc prob doesn't depend on run of course, but giving a unique name avoids a memory leak warning
  
  //for neutrinos
  osc.SetOscParam(OscPar::kNuAntiNu,1);
  h = new TH1D(Form("NominalOscProb_NC_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++) {
    h->SetBinContent(it+1,1);
  }
  NominalOscProb[qNC]=h;
  /*  
  h = new TH1D(Form("NominalOscProb_NuMuNC_%i",RunPeriod),"",nTrue,t);
   for(it=0;it<nTrue;it++)
  {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,(1-osc.OscillateLSND(18,14,E)));}
    else{h->SetBinContent(it+1,0);}
  }
  NominalOscProb[qNuMuNC]=h; 
 
  h = new TH1D(Form("NominalOscProb_BNueNC_%i",RunPeriod),"",nTrue,t);
   for(it=0;it<nTrue;it++)
  {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,(1-osc.OscillateLSND(18,12,E)));}
    else{h->SetBinContent(it+1,0);}
  }
  NominalOscProb[qBNueNC]=h; 
  */
  h = new TH1D(Form("NominalOscProb_NuMuToNuMu_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++) {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,osc.OscillateLSND(14,14,E));}
    else if(OscMethod == 2){h->SetBinContent(it+1,osc.OscillateNSI(14,14,E));}
    else if(OscMethod == 1){h->SetBinContent(it+1,osc.Oscillate(14,14,E));}
    else{h->SetBinContent(it+1,osc.MuToMu(E));}
  }
  NominalOscProb[qNuMuToNuMu]=h;
  
  h = new TH1D(Form("NominalOscProb_BNueToNuMu_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++) {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,osc.OscillateLSND(14,12,E));}
    else if(OscMethod == 2){h->SetBinContent(it+1,osc.OscillateNSI(14,12,E));}
    else if(OscMethod == 1){h->SetBinContent(it+1,osc.Oscillate(14,12,E));}
    else{h->SetBinContent(it+1,osc.ElecToMu(E));}
  }
  NominalOscProb[qBNueToNuMu]=h;
  
  h = new TH1D(Form("NominalOscProb_BNueToBNue_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++) {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,osc.OscillateLSND(12,12,E));}
    else if(OscMethod == 2){h->SetBinContent(it+1,osc.OscillateNSI(12,12,E));}
    else if(OscMethod == 1){h->SetBinContent(it+1,osc.Oscillate(12,12,E));}
    else{h->SetBinContent(it+1,osc.ElecToElec(E));}
  }
  NominalOscProb[qBNueToBNue]=h;
  
  h = new TH1D(Form("NominalOscProb_NuMuToNue_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++) {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,osc.OscillateLSND(12,14,E));}
    else if(OscMethod == 2){h->SetBinContent(it+1,osc.OscillateNSI(12,14,E));}
    else if(OscMethod == 1){h->SetBinContent(it+1,osc.Oscillate(12,14,E));}    
    else{h->SetBinContent(it+1,osc.MuToElec(E));}
  }
  NominalOscProb[qNuMuToNue]=h;
  
  h = new TH1D(Form("NominalOscProb_NuMuToNuTau_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++) {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,osc.OscillateLSND(16,14,E));}
    else if(OscMethod == 2){h->SetBinContent(it+1,osc.OscillateNSI(16,14,E));}
    else if(OscMethod == 1){h->SetBinContent(it+1,osc.Oscillate(16,14,E));}
    else{h->SetBinContent(it+1,osc.MuToTau(E));}
  }
  NominalOscProb[qNuMuToNuTau]=h;
  
  h = new TH1D(Form("NominalOscProb_BNueToNuTau_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,osc.OscillateLSND(16,12,E));}
    else if(OscMethod == 2){h->SetBinContent(it+1,osc.OscillateNSI(16,12,E));}
    else if(OscMethod == 1){h->SetBinContent(it+1,osc.Oscillate(16,12,E));}
    else{h->SetBinContent(it+1,osc.ElecToTau(E));}
  }
  NominalOscProb[qBNueToNuTau]=h;
  
  //for antineutrinos
  osc.SetOscParam(OscPar::kNuAntiNu,-1);
  h = new TH1D(Form("NominalOscProb_NuBar_NC_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    h->SetBinContent(it+1,1);
  }
  NominalOscProb_NuBar[qNC]=h;
  /*
  h = new TH1D(Form("NominalOscProb_NuBar_NuMuNC_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,(1-osc.OscillateLSND(-18,-14,E)));}
    else{h->SetBinContent(it+1,0);}
  }
  NominalOscProb_NuBar[qNuMuNC]=h;

  h = new TH1D(Form("NominalOscProb_NuBar_BNueNC_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,(1-osc.OscillateLSND(-18,-12,E)));}
    else{h->SetBinContent(it+1,0);}
  }
  NominalOscProb_NuBar[qBNueNC]=h;
  */
  h = new TH1D(Form("NominalOscProb_NuBar_NuMuToNuMu_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,osc.OscillateLSND(-14,-14,E));}
    else if(OscMethod == 2){h->SetBinContent(it+1,osc.OscillateNSI(-14,-14,E));}
    else if(OscMethod == 1){h->SetBinContent(it+1,osc.Oscillate(-14,-14,E));}
    else{h->SetBinContent(it+1,osc.MuToMu(E));}
  }
  NominalOscProb_NuBar[qNuMuToNuMu]=h;
  
  h = new TH1D(Form("NominalOscProb_NuBar_BNueToNuMu_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,osc.OscillateLSND(-14,-12,E));}
    else if(OscMethod == 2){h->SetBinContent(it+1,osc.OscillateNSI(-14,-12,E));}
    else if(OscMethod == 1){h->SetBinContent(it+1,osc.Oscillate(-14,-12,E));}
    else{h->SetBinContent(it+1,osc.ElecToMu(E));}
  }
  NominalOscProb_NuBar[qBNueToNuMu]=h;
  
  h = new TH1D(Form("NominalOscProb_NuBar_BNueToBNue_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,osc.OscillateLSND(-12,-12,E));}
    else if(OscMethod == 2){h->SetBinContent(it+1,osc.OscillateNSI(-12,-12,E));}
    else if(OscMethod == 1){h->SetBinContent(it+1,osc.Oscillate(-12,-12,E));}
    else{h->SetBinContent(it+1,osc.ElecToElec(E));}
  }
  NominalOscProb_NuBar[qBNueToBNue]=h;
  
  h = new TH1D(Form("NominalOscProb_NuBar_NuMuToNue_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,osc.OscillateLSND(-12,-14,E));}
    else if(OscMethod == 2){h->SetBinContent(it+1,osc.OscillateNSI(-12,-14,E));}
    else if(OscMethod == 1){h->SetBinContent(it+1,osc.Oscillate(-12,-14,E));} 
    else{h->SetBinContent(it+1,osc.MuToElec(E));}
  }
  NominalOscProb_NuBar[qNuMuToNue]=h;
  
  h = new TH1D(Form("NominalOscProb_NuBar_NuMuToNuTau_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,osc.OscillateLSND(-16,-14,E));}
    else if(OscMethod == 2){h->SetBinContent(it+1,osc.OscillateNSI(-16,-14,E));}
    else if(OscMethod == 1){h->SetBinContent(it+1,osc.Oscillate(-16,-14,E));}
    else{h->SetBinContent(it+1,osc.MuToTau(E));}
  }
  NominalOscProb_NuBar[qNuMuToNuTau]=h;
  
  h = new TH1D(Form("NominalOscProb_NuBar_BNueToNuTau_%i",RunPeriod),"",nTrue,t);
  for(it=0;it<nTrue;it++)
  {
    E = h->GetBinCenter(it+1);
    if(OscMethod == 3){h->SetBinContent(it+1,osc.OscillateLSND(-16,-12,E));}
    else if(OscMethod == 2){h->SetBinContent(it+1,osc.OscillateNSI(-16,-12,E));}
    else if(OscMethod == 1){h->SetBinContent(it+1,osc.Oscillate(-16,-12,E));}
    else{h->SetBinContent(it+1,osc.ElecToTau(E));}
  }
  NominalOscProb_NuBar[qBNueToNuTau]=h;
  
//   cout<<osc.MuToElec(2.0)<<endl;
  osc.SetOscParam(OscPar::kNuAntiNu,1);
//   cout<<osc.MuToElec(2.0)<<endl;
  return;
}
void Extrapolate2D::ReadNDDataFile()
{
  //read in ND decomposition results
  
  if(RunPeriod<0 || RunPeriod>3)
  {
    cout<<"Failure to read NDData file: choose run period 0 (all runs), 1, 2, or 3"<<endl;
    ReadError = true;
    return;
  }
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(NDData_infile.c_str())))
  {
    cout<<"Failure to read ND data file."<<endl;
    ReadError = true;
    return;
  }
  
  TFile *f = new TFile(gSystem->ExpandPathName(NDData_infile.c_str()),"READ");
  
  TH2D *h;
  
  h = (TH2D*)f->Get(Form("NDData_NuMuCC_%s_Run%i",PID_NDData.c_str(),RunPeriod));
  h->Scale(nPOTNear/1e19);//assumes NDData files have been normalized to 1e19
  NDData[Background::kNuMuCC] = h;
  
  h = (TH2D*)f->Get(Form("NDData_NC_%s_Run%i",PID_NDData.c_str(),RunPeriod));
  h->Scale(nPOTNear/1e19);//assumes NDData files have been normalized to 1e19
  NDData[Background::kNC] = h;
  
  h = (TH2D*)f->Get(Form("NDData_BNueCC_%s_Run%i",PID_NDData.c_str(),RunPeriod));
  h->Scale(nPOTNear/1e19);//assumes NDData files have been normalized to 1e19
  NDData[Background::kBNueCC] = h;
  
  if(NDData[Background::kNC]->GetNbinsX()!=nPID || NDData[Background::kNC]->GetNbinsY()!=nReco)
  {
    cout<<"Warning: ND Data input histograms don't have the correct binning."<<endl;
  }
  
  return;
}
void Extrapolate2D::ReadFNFile()
{
  if(ReadError) return;
  
  FD_TrueVsRecoVsPID.clear();
  ND_DataOverMC_RecoVsPID.clear();
  FD_True2Reco_Fid.clear();
  FD_True2Reco_NuBar_Fid.clear();
  FD_Eff.clear();
  FD_Eff_NuBar.clear();
  FNRatio.clear();
  NDMC.clear();
  FDMC.clear();
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(FN_infile.c_str())))
  {
    cout<<"Failure to read F/N file."<<endl;
    ReadError = true;
    return;
  }
  
  TFile *f = new TFile(gSystem->ExpandPathName(FN_infile.c_str()),"READ");
  
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
//   //NSI
   tree->SetBranchAddress("Eps_ee",&Eps_ee);
   tree->SetBranchAddress("Eps_emu",&Eps_emu);
   tree->SetBranchAddress("Eps_etau",&Eps_etau);
   tree->SetBranchAddress("Eps_mumu",&Eps_mumu);
   tree->SetBranchAddress("Eps_mutau",&Eps_mutau);
   tree->SetBranchAddress("Eps_tautau",&Eps_tautau);
// 
   tree->SetBranchAddress("Delta_emu",&Delta_emu);
   tree->SetBranchAddress("Delta_etau",&Delta_etau);
   tree->SetBranchAddress("Delta_mutau",&Delta_mutau);
//   //LSND
   tree->SetBranchAddress("Theta14",&Theta14);
   tree->SetBranchAddress("Theta24",&Theta24);
   tree->SetBranchAddress("Theta34",&Theta34);
   tree->SetBranchAddress("Dm41",&Dm41);
   tree->SetBranchAddress("Delta14",&Delta14);
   tree->SetBranchAddress("Delta24",&Delta24);   
   //
  
  //these parameters not incorporated into extrap file yet
  Eps_ee=0;
  Eps_emu=0;
  Eps_etau=0;
  Eps_mumu=0;
  Eps_mutau=0;
  Eps_tautau=0;
  
  Delta_emu=0;
  Delta_etau=0;
  Delta_mutau=0;
  
  Theta14=0;
  Theta24=0;
  Theta34=0;
  Dm41=0;
  Delta14=0;
  Delta24=0;
  
  tree->GetEntry(0);
  
  osc.SetOscParam(OscPar::kTh12,Theta12);
  osc.SetOscParam(OscPar::kTh13,Theta13);
  osc.SetOscParam(OscPar::kTh23,Theta23);
  osc.SetOscParam(OscPar::kDeltaM23,DeltaMSq23);
  osc.SetOscParam(OscPar::kDeltaM12,DeltaMSq12);
  osc.SetOscParam(OscPar::kDelta,DeltaCP);
  //NSI
  osc.SetOscParam(OscPar::kEps_ee,Eps_ee);
  osc.SetOscParam(OscPar::kEps_emu,Eps_emu);
  osc.SetOscParam(OscPar::kEps_etau,Eps_etau);
  osc.SetOscParam(OscPar::kEps_mumu,Eps_mumu);
  osc.SetOscParam(OscPar::kEps_mutau,Eps_mutau);
  osc.SetOscParam(OscPar::kEps_tautau,Eps_tautau);

  osc.SetOscParam(OscPar::kDelta_emu,Delta_emu);
  osc.SetOscParam(OscPar::kDelta_etau,Delta_etau);
  osc.SetOscParam(OscPar::kDelta_mutau,Delta_mutau);
  //LSND
  osc.SetOscParam(OscPar::kTh14,Theta14);
  osc.SetOscParam(OscPar::kTh24,Theta24);
  osc.SetOscParam(OscPar::kTh34,Theta34);
  osc.SetOscParam(OscPar::kDm41,Dm41);
  osc.SetOscParam(OscPar::kDelta14,Delta14);
  osc.SetOscParam(OscPar::kDelta24,Delta24);
  //
  SetNominalOscProb();
  
  //total
  FD_TrueVsRecoVsPID[qNC] = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_NC",PID.c_str()));
  FD_TrueVsRecoVsPID[qNuMuToNuMu] = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_NuMuToNuMu",PID.c_str()));
  FD_TrueVsRecoVsPID[qBNueToNuMu] = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_BNueToNuMu",PID.c_str()));
  FD_TrueVsRecoVsPID[qBNueToBNue] = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_BNueToBNue",PID.c_str()));
  FD_TrueVsRecoVsPID[qNuMuToNue] = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_NuMuToNue",PID.c_str()));
  FD_TrueVsRecoVsPID[qNuMuToNuTau] = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_NuMuToNuTau",PID.c_str()));
  FD_TrueVsRecoVsPID[qBNueToNuTau] = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_BNueToNuTau",PID.c_str()));
  
  FD_TrueVsRecoVsPID[qNC]->SetName(Form("FD_TrueVsRecoVs%s_NC",PID.c_str()));
  FD_TrueVsRecoVsPID[qNuMuToNuMu]->SetName(Form("FD_TrueVsRecoVs%s_NuMuToNuMu",PID.c_str()));
  FD_TrueVsRecoVsPID[qBNueToNuMu]->SetName(Form("FD_TrueVsRecoVs%s_BNueToNuMu",PID.c_str()));
  FD_TrueVsRecoVsPID[qBNueToBNue]->SetName(Form("FD_TrueVsRecoVs%s_BNueToBNue",PID.c_str()));
  FD_TrueVsRecoVsPID[qNuMuToNue]->SetName(Form("FD_TrueVsRecoVs%s_NuMuToNue",PID.c_str()));
  FD_TrueVsRecoVsPID[qNuMuToNuTau]->SetName(Form("FD_TrueVsRecoVs%s_NuMuToNuTau",PID.c_str()));
  FD_TrueVsRecoVsPID[qBNueToNuTau]->SetName(Form("FD_TrueVsRecoVs%s_BNueToNuTau",PID.c_str()));

  if( OscMethod == 3){
  FD_TrueVsRecoVsPID_NC_NueFrac = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_NC_NoOscNue",PID.c_str()));
  FD_TrueVsRecoVsPID_NC_NueBarFrac = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_NC_NoOscNueBar",PID.c_str()));
  FD_TrueVsRecoVsPID_NC_NuMuFrac = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_NC_NoOscNuMu",PID.c_str()));
  FD_TrueVsRecoVsPID_NC_NuMuBarFrac = (TH3D*)f->Get(Form("FDMC/TrueVsRecoVs%s_NC_NoOscNuMuBar",PID.c_str()));
  
  FD_TrueVsRecoVsPID_NC_NueFrac->SetName(Form("FD_TrueVsRecoVs%s_NC_NueFrac",PID.c_str()));
  FD_TrueVsRecoVsPID_NC_NueBarFrac->SetName(Form("FD_TrueVsRecoVs%s_NC_NueBarFrac",PID.c_str()));
  FD_TrueVsRecoVsPID_NC_NuMuFrac->SetName(Form("FD_TrueVsRecoVs%s_NC_NuMuFrac",PID.c_str()));
  FD_TrueVsRecoVsPID_NC_NuMuBarFrac->SetName(Form("FD_TrueVsRecoVs%s_NC_NuMuBarFrac",PID.c_str()));
  }

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
    //nu only
    FD_Nu_TrueVsRecoVsPID[qNC] = (TH3D*)f->Get(Form("FDMC_Nu/TrueVsRecoVs%s_NC",PID.c_str()));
    FD_Nu_TrueVsRecoVsPID[qNuMuToNuMu] = (TH3D*)f->Get(Form("FDMC_Nu/TrueVsRecoVs%s_NuMuToNuMu",PID.c_str()));
    FD_Nu_TrueVsRecoVsPID[qBNueToNuMu] = (TH3D*)f->Get(Form("FDMC_Nu/TrueVsRecoVs%s_BNueToNuMu",PID.c_str()));
    FD_Nu_TrueVsRecoVsPID[qBNueToBNue] = (TH3D*)f->Get(Form("FDMC_Nu/TrueVsRecoVs%s_BNueToBNue",PID.c_str()));
    FD_Nu_TrueVsRecoVsPID[qNuMuToNue] = (TH3D*)f->Get(Form("FDMC_Nu/TrueVsRecoVs%s_NuMuToNue",PID.c_str()));
    FD_Nu_TrueVsRecoVsPID[qNuMuToNuTau] = (TH3D*)f->Get(Form("FDMC_Nu/TrueVsRecoVs%s_NuMuToNuTau",PID.c_str()));
    FD_Nu_TrueVsRecoVsPID[qBNueToNuTau] = (TH3D*)f->Get(Form("FDMC_Nu/TrueVsRecoVs%s_BNueToNuTau",PID.c_str()));
    
    FD_Nu_TrueVsRecoVsPID[qNC]->SetName(Form("FD_Nu_TrueVsRecoVs%s_NC",PID.c_str()));
    FD_Nu_TrueVsRecoVsPID[qNuMuToNuMu]->SetName(Form("FD_Nu_TrueVsRecoVs%s_NuMuToNuMu",PID.c_str()));
    FD_Nu_TrueVsRecoVsPID[qBNueToNuMu]->SetName(Form("FD_Nu_TrueVsRecoVs%s_BNueToNuMu",PID.c_str()));
    FD_Nu_TrueVsRecoVsPID[qBNueToBNue]->SetName(Form("FD_Nu_TrueVsRecoVs%s_BNueToBNue",PID.c_str()));
    FD_Nu_TrueVsRecoVsPID[qNuMuToNue]->SetName(Form("FD_Nu_TrueVsRecoVs%s_NuMuToNue",PID.c_str()));
    FD_Nu_TrueVsRecoVsPID[qNuMuToNuTau]->SetName(Form("FD_Nu_TrueVsRecoVs%s_NuMuToNuTau",PID.c_str()));
    FD_Nu_TrueVsRecoVsPID[qBNueToNuTau]->SetName(Form("FD_Nu_TrueVsRecoVs%s_BNueToNuTau",PID.c_str()));
    
    //nubar only
    FD_NuBar_TrueVsRecoVsPID[qNC] = (TH3D*)f->Get(Form("FDMC_NuBar/TrueVsRecoVs%s_NC",PID.c_str()));
    FD_NuBar_TrueVsRecoVsPID[qNuMuToNuMu] = (TH3D*)f->Get(Form("FDMC_NuBar/TrueVsRecoVs%s_NuMuToNuMu",PID.c_str()));
    FD_NuBar_TrueVsRecoVsPID[qBNueToNuMu] = (TH3D*)f->Get(Form("FDMC_NuBar/TrueVsRecoVs%s_BNueToNuMu",PID.c_str()));
    FD_NuBar_TrueVsRecoVsPID[qBNueToBNue] = (TH3D*)f->Get(Form("FDMC_NuBar/TrueVsRecoVs%s_BNueToBNue",PID.c_str()));
    FD_NuBar_TrueVsRecoVsPID[qNuMuToNue] = (TH3D*)f->Get(Form("FDMC_NuBar/TrueVsRecoVs%s_NuMuToNue",PID.c_str()));
    FD_NuBar_TrueVsRecoVsPID[qNuMuToNuTau] = (TH3D*)f->Get(Form("FDMC_NuBar/TrueVsRecoVs%s_NuMuToNuTau",PID.c_str()));
    FD_NuBar_TrueVsRecoVsPID[qBNueToNuTau] = (TH3D*)f->Get(Form("FDMC_NuBar/TrueVsRecoVs%s_BNueToNuTau",PID.c_str()));
    
    FD_NuBar_TrueVsRecoVsPID[qNC]->SetName(Form("FD_NuBar_TrueVsRecoVs%s_NC",PID.c_str()));
    FD_NuBar_TrueVsRecoVsPID[qNuMuToNuMu]->SetName(Form("FD_NuBar_TrueVsRecoVs%s_NuMuToNuMu",PID.c_str()));
    FD_NuBar_TrueVsRecoVsPID[qBNueToNuMu]->SetName(Form("FD_NuBar_TrueVsRecoVs%s_BNueToNuMu",PID.c_str()));
    FD_NuBar_TrueVsRecoVsPID[qBNueToBNue]->SetName(Form("FD_NuBar_TrueVsRecoVs%s_BNueToBNue",PID.c_str()));
    FD_NuBar_TrueVsRecoVsPID[qNuMuToNue]->SetName(Form("FD_NuBar_TrueVsRecoVs%s_NuMuToNue",PID.c_str()));
    FD_NuBar_TrueVsRecoVsPID[qNuMuToNuTau]->SetName(Form("FD_NuBar_TrueVsRecoVs%s_NuMuToNuTau",PID.c_str()));
    FD_NuBar_TrueVsRecoVsPID[qBNueToNuTau]->SetName(Form("FD_NuBar_TrueVsRecoVs%s_BNueToNuTau",PID.c_str()));
  }
  
  FD_TrueVsReco_Fid = (TH2D*)f->Get("FDMC/TrueVsReco_Fid");
  FD_TrueVsReco_Fid_NuMuCC = (TH2D*)f->Get("FDMC/TrueVsReco_Fid_NuMuCC");
  FD_TrueVsReco_CClike = (TH2D*)f->Get("FDMC/TrueVsReco_CClike");
  FD_TrueVsReco_CClike_NuMuCC = (TH2D*)f->Get("FDMC/TrueVsReco_CClike_NuMuCC");
  
  FD_TrueVsReco_Fid->SetName("FD_TrueVsReco_Fid");
  FD_TrueVsReco_Fid_NuMuCC->SetName("FD_TrueVsReco_Fid_NuMuCC");
  FD_TrueVsReco_CClike->SetName("FD_TrueVsReco_CClike");
  FD_TrueVsReco_CClike_NuMuCC->SetName("FD_TrueVsReco_CClike_NuMuCC");
  
  if(UseSeparateNuNuBar)
  {
    FD_TrueVsReco_Fid_NuMuBarCC = (TH2D*)f->Get("FDMC/TrueVsReco_Fid_NuMuBarCC");
    FD_TrueVsReco_CClike_Pos = (TH2D*)f->Get("FDMC/TrueVsReco_CClike_Pos");
    FD_TrueVsReco_CClike_Pos_NuMuBarCC = (TH2D*)f->Get("FDMC/TrueVsReco_CClike_Pos_NuMuBarCC");
    FD_TrueVsReco_CClike_Neg = (TH2D*)f->Get("FDMC/TrueVsReco_CClike_Neg");
    FD_TrueVsReco_CClike_Neg_NuMuCC = (TH2D*)f->Get("FDMC/TrueVsReco_CClike_Neg_NuMuCC");
    
    FD_TrueVsReco_Fid_NuMuBarCC->SetName("FD_TrueVsReco_Fid_NuMuBarCC");
    FD_TrueVsReco_CClike_Pos->SetName("FD_TrueVsReco_CClike_Pos");
    FD_TrueVsReco_CClike_Pos_NuMuBarCC->SetName("FD_TrueVsReco_CClike_Pos_NuMuBarCC");
    FD_TrueVsReco_CClike_Neg->SetName("FD_TrueVsReco_CClike_Neg");
    FD_TrueVsReco_CClike_Neg_NuMuCC->SetName("FD_TrueVsReco_CClike_Neg_NuMuCC");
  }
  
  FD_TrueVsReco_Fid_NueCC = (TH2D*)f->Get("FDMC/TrueVsReco_Fid_NueCC");
  FD_TrueVsReco_Fid_NuTauCC = (TH2D*)f->Get("FDMC/TrueVsReco_Fid_NuTauCC");
  FD_True_Fid_NueCC = (TH1D*)f->Get("FDMC/True_Fid_NueCC");
  FD_True_Fid_NuTauCC = (TH1D*)f->Get("FDMC/True_Fid_NuTauCC");
  
  if(UseSeparateNuNuBar)
  {
    FD_TrueVsReco_Fid_NueBarCC = (TH2D*)f->Get("FDMC/TrueVsReco_Fid_NueBarCC");
    FD_TrueVsReco_Fid_NuTauBarCC = (TH2D*)f->Get("FDMC/TrueVsReco_Fid_NuTauBarCC");
    FD_True_Fid_NueBarCC = (TH1D*)f->Get("FDMC/True_Fid_NueBarCC");
    FD_True_Fid_NuTauBarCC = (TH1D*)f->Get("FDMC/True_Fid_NuTauBarCC");
  }
  
  ND_TrueVsRecoVsPID[qNC] = (TH3D*)f->Get(Form("NDMC/TrueVsRecoVs%s_NC",PID.c_str()));
  ND_TrueVsRecoVsPID[qNuMuToNuMu] = (TH3D*)f->Get(Form("NDMC/TrueVsRecoVs%s_NuMuToNuMu",PID.c_str()));
  ND_TrueVsRecoVsPID[qBNueToBNue] = (TH3D*)f->Get(Form("NDMC/TrueVsRecoVs%s_BNueToBNue",PID.c_str()));
  
  ND_TrueVsRecoVsPID[qNC]->SetName(Form("ND_TrueVsRecoVs%s_NC",PID.c_str()));
  ND_TrueVsRecoVsPID[qNuMuToNuMu]->SetName(Form("ND_TrueVsRecoVs%s_NuMuToNuMu",PID.c_str()));
  ND_TrueVsRecoVsPID[qBNueToBNue]->SetName(Form("ND_TrueVsRecoVs%s_BNueToBNue",PID.c_str()));
  
  ND_TrueVsReco_Fid = (TH2D*)f->Get("NDMC/TrueVsReco_Fid");
  ND_TrueVsReco_Fid_NuMuCC = (TH2D*)f->Get("NDMC/TrueVsReco_Fid_NuMuCC");
  ND_TrueVsReco_CClike = (TH2D*)f->Get("NDMC/TrueVsReco_CClike");
  ND_TrueVsReco_CClike_NuMuCC = (TH2D*)f->Get("NDMC/TrueVsReco_CClike_NuMuCC");
  
  ND_TrueVsReco_Fid->SetName("ND_TrueVsReco_Fid");
  ND_TrueVsReco_Fid_NuMuCC->SetName("ND_TrueVsReco_Fid_NuMuCC");
  ND_TrueVsReco_CClike->SetName("ND_TrueVsReco_CClike");
  ND_TrueVsReco_CClike_NuMuCC->SetName("ND_TrueVsReco_CClike_NuMuCC");
  
  if(UseSeparateNuNuBar)
  {
    ND_TrueVsReco_Fid_NuMuBarCC = (TH2D*)f->Get("NDMC/TrueVsReco_Fid_NuMuBarCC");
    ND_TrueVsReco_CClike_Pos = (TH2D*)f->Get("NDMC/TrueVsReco_CClike_Pos");
    ND_TrueVsReco_CClike_Neg = (TH2D*)f->Get("NDMC/TrueVsReco_CClike_Neg");
    ND_TrueVsReco_CClike_Pos_NuMuBarCC = (TH2D*)f->Get("NDMC/TrueVsReco_CClike_Pos_NuMuBarCC");
    ND_TrueVsReco_CClike_Neg_NuMuCC = (TH2D*)f->Get("NDMC/TrueVsReco_CClike_Neg_NuMuCC");
    
    ND_TrueVsReco_Fid_NuMuBarCC->SetName("ND_TrueVsReco_Fid_NuMuBarCC");
    ND_TrueVsReco_CClike_Pos->SetName("ND_TrueVsReco_CClike_Pos");
    ND_TrueVsReco_CClike_Neg->SetName("ND_TrueVsReco_CClike_Neg");
    ND_TrueVsReco_CClike_Pos_NuMuBarCC->SetName("ND_TrueVsReco_CClike_Pos_NuMuBarCC");
    ND_TrueVsReco_CClike_Neg_NuMuCC->SetName("ND_TrueVsReco_CClike_Neg_NuMuCC");
  }
  
  NDData_Reco_CClike = (TH1D*)f->Get("NDData/Reco_CClike");
  if(UseSeparateNuNuBar)
  {
    NDData_Reco_CClike_Pos = (TH1D*)f->Get("NDData/Reco_CClike_Pos");
    NDData_Reco_CClike_Neg = (TH1D*)f->Get("NDData/Reco_CClike_Neg");
  }
  
  FD_TrueVsRecoVsPID[qNC]->Scale(nPOTFar/fp);
  FD_TrueVsRecoVsPID[qNuMuToNuMu]->Scale(nPOTFar/fp);
  FD_TrueVsRecoVsPID[qBNueToNuMu]->Scale(nPOTFar/fp);
  FD_TrueVsRecoVsPID[qBNueToBNue]->Scale(nPOTFar/fp);
  FD_TrueVsRecoVsPID[qNuMuToNue]->Scale(nPOTFar/fp);
  FD_TrueVsRecoVsPID[qNuMuToNuTau]->Scale(nPOTFar/fp);
  FD_TrueVsRecoVsPID[qBNueToNuTau]->Scale(nPOTFar/fp);
  
  FD_TrueVsReco_Fid->Scale(nPOTFar/fp);
  FD_TrueVsReco_Fid_NuMuCC->Scale(nPOTFar/fp);
  FD_TrueVsReco_CClike->Scale(nPOTFar/fp);
  FD_TrueVsReco_CClike_NuMuCC->Scale(nPOTFar/fp);
  
  FD_TrueVsReco_Fid_NueCC->Scale(nPOTFar/fp);
  FD_TrueVsReco_Fid_NuTauCC->Scale(nPOTFar/fp);
  FD_True_Fid_NueCC->Scale(nPOTFar/fp);
  FD_True_Fid_NuTauCC->Scale(nPOTFar/fp);
  
  ND_TrueVsRecoVsPID[qNC]->Scale(nPOTNear/np);
  ND_TrueVsRecoVsPID[qNuMuToNuMu]->Scale(nPOTNear/np);
  ND_TrueVsRecoVsPID[qBNueToBNue]->Scale(nPOTNear/np);
  
  ND_TrueVsReco_Fid->Scale(nPOTNear/np);
  ND_TrueVsReco_Fid_NuMuCC->Scale(nPOTNear/np);
  ND_TrueVsReco_CClike->Scale(nPOTNear/np);
  ND_TrueVsReco_CClike_NuMuCC->Scale(nPOTNear/np);
  
  NDData_Reco_CClike->Scale(nPOTNear/np);
  
  if(UseSeparateNuNuBar)
  {
    FD_Nu_TrueVsRecoVsPID[qNC]->Scale(nPOTFar/fp);
    FD_Nu_TrueVsRecoVsPID[qNuMuToNuMu]->Scale(nPOTFar/fp);
    FD_Nu_TrueVsRecoVsPID[qBNueToNuMu]->Scale(nPOTFar/fp);
    FD_Nu_TrueVsRecoVsPID[qBNueToBNue]->Scale(nPOTFar/fp);
    FD_Nu_TrueVsRecoVsPID[qNuMuToNue]->Scale(nPOTFar/fp);
    FD_Nu_TrueVsRecoVsPID[qNuMuToNuTau]->Scale(nPOTFar/fp);
    FD_Nu_TrueVsRecoVsPID[qBNueToNuTau]->Scale(nPOTFar/fp);
    
    FD_NuBar_TrueVsRecoVsPID[qNC]->Scale(nPOTFar/fp);
    FD_NuBar_TrueVsRecoVsPID[qNuMuToNuMu]->Scale(nPOTFar/fp);
    FD_NuBar_TrueVsRecoVsPID[qBNueToNuMu]->Scale(nPOTFar/fp);
    FD_NuBar_TrueVsRecoVsPID[qBNueToBNue]->Scale(nPOTFar/fp);
    FD_NuBar_TrueVsRecoVsPID[qNuMuToNue]->Scale(nPOTFar/fp);
    FD_NuBar_TrueVsRecoVsPID[qNuMuToNuTau]->Scale(nPOTFar/fp);
    FD_NuBar_TrueVsRecoVsPID[qBNueToNuTau]->Scale(nPOTFar/fp);
    
    FD_TrueVsReco_Fid_NuMuBarCC->Scale(nPOTFar/fp);
    
    FD_TrueVsReco_CClike_Pos->Scale(nPOTFar/fp);
    FD_TrueVsReco_CClike_Neg->Scale(nPOTFar/fp);
    
    FD_TrueVsReco_CClike_Pos_NuMuBarCC->Scale(nPOTFar/fp);
    FD_TrueVsReco_CClike_Neg_NuMuCC->Scale(nPOTFar/fp);
    
    FD_TrueVsReco_Fid_NueBarCC->Scale(nPOTFar/fp);
    FD_TrueVsReco_Fid_NuTauBarCC->Scale(nPOTFar/fp);
    FD_True_Fid_NueBarCC->Scale(nPOTFar/fp);
    FD_True_Fid_NuTauBarCC->Scale(nPOTFar/fp);
    
    ND_TrueVsReco_Fid_NuMuBarCC->Scale(nPOTNear/np);
    
    ND_TrueVsReco_CClike_Pos->Scale(nPOTNear/np);
    ND_TrueVsReco_CClike_Neg->Scale(nPOTNear/np);
    
    ND_TrueVsReco_CClike_Pos_NuMuBarCC->Scale(nPOTNear/np);
    ND_TrueVsReco_CClike_Neg_NuMuCC->Scale(nPOTNear/np);
    
    NDData_Reco_CClike_Pos->Scale(nPOTNear/np);
    NDData_Reco_CClike_Neg->Scale(nPOTNear/np);
  }
  
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
  
  if(FD_TrueVsRecoVsPID[qNC]->GetNbinsX()!=nPID)
  {
    RebinP = true;
  }
  if(FD_TrueVsRecoVsPID[qNC]->GetNbinsY()!=nReco)
  {
    RebinE = true;
  }
  
  if(RebinE || RebinP)
  {
    cout<<"Warning: F/N input hists don't have correct binning. Rebinning them now..."<<endl;
    RebinInputHists();
    cout<<"Hists rebinned!"<<endl;
  }
  
  if ( OscMethod == 3){
  TH3D *NCTotal = (TH3D*)FD_TrueVsRecoVsPID_NC_NueFrac->Clone("NCTotal");
  NCTotal->Add(FD_TrueVsRecoVsPID_NC_NueBarFrac);
  NCTotal->Add(FD_TrueVsRecoVsPID_NC_NuMuFrac);
  NCTotal->Add(FD_TrueVsRecoVsPID_NC_NuMuBarFrac);
  
  FD_TrueVsRecoVsPID_NC_NueFrac->Divide(NCTotal);
  FD_TrueVsRecoVsPID_NC_NueBarFrac->Divide(NCTotal);
  FD_TrueVsRecoVsPID_NC_NuMuFrac->Divide(NCTotal);
  FD_TrueVsRecoVsPID_NC_NuMuBarFrac->Divide(NCTotal);
  }
  
  int ip,ir,it;
  double temp;
  
  TH2D *h2;
  TH1D *h1;
  
  h2 = (TH2D*)FD_TrueVsRecoVsPID[qNuMuToNuMu]->Project3D("yx");
  FDMC[Background::kNuMuCC] = (TH2D*)h2->Clone("FDMC_NuMuCC");
  FNRatio[Background::kNuMuCC] = (TH2D*)h2->Clone("FNRatio_NuMuCC");
  
  h2 = (TH2D*)FD_TrueVsRecoVsPID[qBNueToNuMu]->Project3D("yx");
  FDMC[Background::kNuMuCC]->Add(h2);
  FNRatio[Background::kNuMuCC]->Add(h2);
  
  h2 = (TH2D*)FD_TrueVsRecoVsPID[qNC]->Project3D("yx");
  FDMC[Background::kNC] = (TH2D*)h2->Clone("FDMC_NC");
  FNRatio[Background::kNC] = (TH2D*)h2->Clone("FNRatio_NC");
  
  h2 = (TH2D*)FD_TrueVsRecoVsPID[qBNueToBNue]->Project3D("yx");
  FDMC[Background::kBNueCC] = (TH2D*)h2->Clone("FDMC_BNueCC");
  FNRatio[Background::kBNueCC] = (TH2D*)h2->Clone("FNRatio_BNueCC");
  
  h2 = (TH2D*)ND_TrueVsRecoVsPID[qNuMuToNuMu]->Project3D("yx");
  NDMC[Background::kNuMuCC] = (TH2D*)h2->Clone("NDMC_NuMuCC");
  FNRatio[Background::kNuMuCC]->Divide(h2);
  ND_DataOverMC_RecoVsPID[Background::kNuMuCC] = new TH2D("ND_DataOverMC_RecoVsPID_NuMuCC","",nPID,p,nReco,r);
  for(ip=0;ip<nPID;ip++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      temp=0;
      if(h2->GetBinContent(ip+1,ir+1)>0)
      {
        temp = NDData[Background::kNuMuCC]->GetBinContent(ip+1,ir+1)/h2->GetBinContent(ip+1,ir+1);
      }
      ND_DataOverMC_RecoVsPID[Background::kNuMuCC]->SetBinContent(ip+1,ir+1,temp);
    }
  }
  
  h2 = (TH2D*)ND_TrueVsRecoVsPID[qNC]->Project3D("yx");
  NDMC[Background::kNC] = (TH2D*)h2->Clone("NDMC_NC");
  FNRatio[Background::kNC]->Divide(h2);
  ND_DataOverMC_RecoVsPID[Background::kNC] = new TH2D("ND_DataOverMC_RecoVsPID_NC","",nPID,p,nReco,r);
  for(ip=0;ip<nPID;ip++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      temp=0;
      if(h2->GetBinContent(ip+1,ir+1)>0)
      {
        temp = NDData[Background::kNC]->GetBinContent(ip+1,ir+1)/h2->GetBinContent(ip+1,ir+1);
      }
      ND_DataOverMC_RecoVsPID[Background::kNC]->SetBinContent(ip+1,ir+1,temp);
    }
  }
  
  h2 = (TH2D*)ND_TrueVsRecoVsPID[qBNueToBNue]->Project3D("yx");
  NDMC[Background::kBNueCC] = (TH2D*)h2->Clone("NDMC_BNueCC");
  FNRatio[Background::kBNueCC]->Divide(h2);
  ND_DataOverMC_RecoVsPID[Background::kBNueCC] = new TH2D("ND_DataOverMC_RecoVsPID_BNueCC","",nPID,p,nReco,r);
  for(ip=0;ip<nPID;ip++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      temp=0;
      if(h2->GetBinContent(ip+1,ir+1)>0)
      {
        temp = NDData[Background::kBNueCC]->GetBinContent(ip+1,ir+1)/h2->GetBinContent(ip+1,ir+1);
      }
      ND_DataOverMC_RecoVsPID[Background::kBNueCC]->SetBinContent(ip+1,ir+1,temp);
    }
  }
  
  ND_Reco_CClike = ND_TrueVsReco_CClike->ProjectionX("ND_Reco_CClike");
  if(UseSeparateNuNuBar)
  {
    ND_Reco_CClike_Pos = ND_TrueVsReco_CClike_Pos->ProjectionX("ND_Reco_CClike_Pos");
    ND_Reco_CClike_Neg = ND_TrueVsReco_CClike_Neg->ProjectionX("ND_Reco_CClike_Neg");
  }
  
  FD_Reco_CClike = FD_TrueVsReco_CClike->ProjectionX("FD_Reco_CClike");
  if(UseSeparateNuNuBar)
  {
    FD_Reco_CClike_Pos = FD_TrueVsReco_CClike_Pos->ProjectionX("FD_Reco_CClike_Pos");
    FD_Reco_CClike_Neg = FD_TrueVsReco_CClike_Neg->ProjectionX("FD_Reco_CClike_Neg");
  }
  
  FD_Reco2True_CClike = new TH2D("FD_Reco2True_CClike","",nRecoCC,rcc,nTrue,t);
  if(UseSeparateNuNuBar)
  {
    FD_Reco2True_CClike_Pos = new TH2D("FD_Reco2True_CClike_Pos","",nRecoCC,rcc,nTrue,t);
    FD_Reco2True_CClike_Neg = new TH2D("FD_Reco2True_CClike_Neg","",nRecoCC,rcc,nTrue,t);
  }
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
      
      if(UseSeparateNuNuBar)
      {
        temp=0;
        if(FD_Reco_CClike_Pos->GetBinContent(ir+1)>0)
        {
          temp = FD_TrueVsReco_CClike_Pos->GetBinContent(ir+1,it+1)/FD_Reco_CClike_Pos->GetBinContent(ir+1);
        }
        FD_Reco2True_CClike_Pos->SetBinContent(ir+1,it+1,temp);
        
        temp=0;
        if(FD_Reco_CClike_Neg->GetBinContent(ir+1)>0)
        {
          temp = FD_TrueVsReco_CClike_Neg->GetBinContent(ir+1,it+1)/FD_Reco_CClike_Neg->GetBinContent(ir+1);
        }
        FD_Reco2True_CClike_Neg->SetBinContent(ir+1,it+1,temp);
      }
    }
  }
  
  FD_Purity_CC = new TH1D("FD_Purity_CC","",nTrue,t);
  if(UseSeparateNuNuBar)
  {
    FD_Purity_NuMuCC_Neg = new TH1D("FD_Purity_NuMuCC_Neg","",nTrue,t);
    FD_Purity_NuMuBarCC_Pos = new TH1D("FD_Purity_NuMuBarCC_Pos","",nTrue,t);
  }
  
  FD_Eff_CC = new TH1D("FD_Eff_CC","",nTrue,t);
  if(UseSeparateNuNuBar)
  {
    FD_Eff_NuMuCC_Neg = new TH1D("FD_Eff_NuMuCC_Neg","",nTrue,t);
    FD_Eff_NuMuBarCC_Pos = new TH1D("FD_Eff_NuMuBarCC_Pos","",nTrue,t);
  }
  
  if(UseSeparateNuNuBar)
  {
    TH1D *FD_True_CClike_Neg = FD_TrueVsReco_CClike_Neg->ProjectionY("FD_True_CClike_Neg");
    TH1D *FD_True_CClike_Neg_NuMuCC = FD_TrueVsReco_CClike_Neg_NuMuCC->ProjectionY("FD_True_CClike_Neg_NuMuCC");
    TH1D *FD_True_Fid_NuMuCC = FD_TrueVsReco_Fid_NuMuCC->ProjectionY("FD_True_Fid_NuMuCC");
    FD_Purity_NuMuCC_Neg->Divide(FD_True_CClike_Neg_NuMuCC,FD_True_CClike_Neg,1,1);
    FD_Eff_NuMuCC_Neg->Divide(FD_True_CClike_Neg_NuMuCC,FD_True_Fid_NuMuCC,1,1);
    
    TH1D *FD_True_CClike_Pos = FD_TrueVsReco_CClike_Pos->ProjectionY("FD_True_CClike_Pos");
    TH1D *FD_True_CClike_Pos_NuMuBarCC = FD_TrueVsReco_CClike_Pos_NuMuBarCC->ProjectionY("FD_True_CClike_Pos_NuMuBarCC");
    TH1D *FD_True_Fid_NuMuBarCC = FD_TrueVsReco_Fid_NuMuBarCC->ProjectionY("FD_True_Fid_NuMuBarCC");
    FD_Purity_NuMuBarCC_Pos->Divide(FD_True_CClike_Pos_NuMuBarCC,FD_True_CClike_Pos,1,1);
    FD_Eff_NuMuBarCC_Pos->Divide(FD_True_CClike_Pos_NuMuBarCC,FD_True_Fid_NuMuBarCC,1,1);
  }
  else
  {
    TH1D *FD_True_CClike = FD_TrueVsReco_CClike->ProjectionY("FD_True_CClike");
    TH1D *FD_True_CClike_NuMuCC = FD_TrueVsReco_CClike_NuMuCC->ProjectionY("FD_True_CClike_NuMuCC");
    TH1D *FD_True_Fid_NuMuCC = FD_TrueVsReco_Fid_NuMuCC->ProjectionY("FD_True_Fid_NuMuCC");
    FD_Purity_CC->Divide(FD_True_CClike_NuMuCC,FD_True_CClike,1,1);
    FD_Eff_CC->Divide(FD_True_CClike_NuMuCC,FD_True_Fid_NuMuCC,1,1);
  }
  
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
  
  if(UseSeparateNuNuBar)
  {
    FD_True2Reco_NuBar_Fid[Background::kNueCC] = new TH2D("FD_True2Reco_Fid_NueBarCC","",nReco,r,nTrue,t);
    for(it=0;it<nTrue;it++)
    {
      for(ir=0;ir<nReco;ir++)
      {
        temp=0;
        if(FD_True_Fid_NueBarCC->GetBinContent(it+1)>0)
        {
          temp = FD_TrueVsReco_Fid_NueBarCC->GetBinContent(ir+1,it+1)/FD_True_Fid_NueBarCC->GetBinContent(it+1);
        }
        FD_True2Reco_NuBar_Fid[Background::kNueCC]->SetBinContent(ir+1,it+1,temp);
      }
    }
      
    FD_True2Reco_NuBar_Fid[Background::kNuTauCC] = new TH2D("FD_True2Reco_Fid_NuTauBarCC","",nReco,r,nTrue,t);
    for(it=0;it<nTrue;it++)
    {
      for(ir=0;ir<nReco;ir++)
      {
        temp=0;
        if(FD_True_Fid_NuTauBarCC->GetBinContent(it+1)>0)
        {
          temp = FD_TrueVsReco_Fid_NuTauBarCC->GetBinContent(ir+1,it+1)/FD_True_Fid_NuTauBarCC->GetBinContent(it+1);
        }
       FD_True2Reco_NuBar_Fid[Background::kNuTauCC]->SetBinContent(ir+1,it+1,temp);
      }
    }
  }
  
  if(UseSeparateNuNuBar)
  {
    FD_Eff[Background::kNueCC] = new TH2D("FD_Eff_NueCC","",nPID,p,nReco,r);
    h2 = (TH2D*)FD_Nu_TrueVsRecoVsPID[qNuMuToNue]->Project3D("yx");
    h1 = FD_TrueVsReco_Fid_NueCC->ProjectionX();
    for(ip=0;ip<nPID;ip++)
    {
      for(ir=0;ir<nReco;ir++)
      {
        temp=0;
        if(h1->GetBinContent(ir+1)>0)
        {
          temp = h2->GetBinContent(ip+1,ir+1)/h1->GetBinContent(ir+1);
        }
        FD_Eff[Background::kNueCC]->SetBinContent(ip+1,ir+1,temp);
      }
    }
    if(ReadMREFromFile) FD_Eff[Background::kNueCC]->Multiply(MREEffRatio);
    
    FD_Eff_NuBar[Background::kNueCC] = new TH2D("FD_Eff_NueBarCC","",nPID,p,nReco,r);
    h2 = (TH2D*)FD_NuBar_TrueVsRecoVsPID[qNuMuToNue]->Project3D("yx");
    h1 = FD_TrueVsReco_Fid_NueBarCC->ProjectionX();
    for(ip=0;ip<nPID;ip++)
    {
      for(ir=0;ir<nReco;ir++)
      {
        temp=0;
        if(h1->GetBinContent(ir+1)>0)
        {
          temp = h2->GetBinContent(ip+1,ir+1)/h1->GetBinContent(ir+1);
        }
        FD_Eff_NuBar[Background::kNueCC]->SetBinContent(ip+1,ir+1,temp);
      }
    }
    if(ReadMREFromFile) FD_Eff_NuBar[Background::kNueCC]->Multiply(MREEffRatio);
  }
  else
  {
    FD_Eff[Background::kNueCC] = new TH2D("FD_Eff_NueCC","",nPID,p,nReco,r);
    h2 = (TH2D*)FD_TrueVsRecoVsPID[qNuMuToNue]->Project3D("yx");
    h1 = FD_TrueVsReco_Fid_NueCC->ProjectionX();
    for(ip=0;ip<nPID;ip++)
    {
      for(ir=0;ir<nReco;ir++)
      {
        temp=0;
        if(h1->GetBinContent(ir+1)>0)
        {
          temp = h2->GetBinContent(ip+1,ir+1)/h1->GetBinContent(ir+1);
        }
        FD_Eff[Background::kNueCC]->SetBinContent(ip+1,ir+1,temp);
      }
    }
    if(ReadMREFromFile) FD_Eff[Background::kNueCC]->Multiply(MREEffRatio);
  }
  
  if(UseSeparateNuNuBar)
  {
    FD_Eff[Background::kNuTauCC] = new TH2D("FD_Eff_NuTauCC","",nPID,p,nReco,r);
    h2 = (TH2D*)FD_Nu_TrueVsRecoVsPID[qNuMuToNuTau]->Project3D("yx");
    h2->Add((TH2D*)FD_Nu_TrueVsRecoVsPID[qBNueToNuTau]->Project3D("yx"));
    h1 = FD_TrueVsReco_Fid_NuTauCC->ProjectionX();
    for(ip=0;ip<nPID;ip++)
    {
      for(ir=0;ir<nReco;ir++)
      {
        temp=0;
        if(h1->GetBinContent(ir+1)>0)
        {
          temp = h2->GetBinContent(ip+1,ir+1)/h1->GetBinContent(ir+1);
        }
        FD_Eff[Background::kNuTauCC]->SetBinContent(ip+1,ir+1,temp);
      }
    }
    
    FD_Eff_NuBar[Background::kNuTauCC] = new TH2D("FD_Eff_NuTauBarCC","",nPID,p,nReco,r);
    h2 = (TH2D*)FD_NuBar_TrueVsRecoVsPID[qNuMuToNuTau]->Project3D("yx");
    h2->Add((TH2D*)FD_NuBar_TrueVsRecoVsPID[qBNueToNuTau]->Project3D("yx"));
    h1 = FD_TrueVsReco_Fid_NuTauBarCC->ProjectionX();
    for(ip=0;ip<nPID;ip++)
    {
      for(ir=0;ir<nReco;ir++)
      {
        temp=0;
        if(h1->GetBinContent(ir+1)>0)
        {
          temp = h2->GetBinContent(ip+1,ir+1)/h1->GetBinContent(ir+1);
        }
        FD_Eff_NuBar[Background::kNuTauCC]->SetBinContent(ip+1,ir+1,temp);
      }
    }
  }
  else
  {
    FD_Eff[Background::kNuTauCC] = new TH2D("FD_Eff_NuTauCC","",nPID,p,nReco,r);
    h2 = (TH2D*)FD_TrueVsRecoVsPID[qNuMuToNuTau]->Project3D("yx");
    h2->Add((TH2D*)FD_TrueVsRecoVsPID[qBNueToNuTau]->Project3D("yx"));
    h1 = FD_TrueVsReco_Fid_NuTauCC->ProjectionX();
    for(ip=0;ip<nPID;ip++)
    {
      for(ir=0;ir<nReco;ir++)
      {
        temp=0;
        if(h1->GetBinContent(ir+1)>0)
        {
          temp = h2->GetBinContent(ip+1,ir+1)/h1->GetBinContent(ir+1);
        }
        FD_Eff[Background::kNuTauCC]->SetBinContent(ip+1,ir+1,temp);
      }
    }
  }
  
  delete [] r;
  delete [] p;
  delete [] t;
  
  return;
}
void Extrapolate2D::ReadMREFile()
{
  if(!ReadMREFromFile)
  {
    MREEffRatio = 0;
    return;
  }
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(MRE_infile.c_str())))
  {
    cout<<"Failure to read MRE file."<<endl;
    ReadError = true;
    return;
  }
  
  TFile *f = new TFile(gSystem->ExpandPathName(MRE_infile.c_str()),"READ");
  MREEffRatio = (TH2D*)f->Get(Form("%s/eff_ratio",PID.c_str()));
  
  if(MREEffRatio->GetNbinsX()!=nPID || MREEffRatio->GetNbinsY()!=nReco)
  {
    cout<<"Warning: MRE correction histogram doesn't have the correct binning."<<endl;
  }
  
  return;
}
void Extrapolate2D::ReadXSecFile()
{
  //read in tau/nue cross section file
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(XSec_infile.c_str())))
  {
    cout<<"Failure to read cross section file."<<endl;
    ReadError = true;
    return;
  }
  
  TFile *fxsec = new TFile(gSystem->ExpandPathName(XSec_infile.c_str()),"READ");
  string names[5] = {"tot", "qe", "res", "dis", "coh"};
  
  for(int i = 0; i < 5; i++)
  {
    string id = "h_numu_cc_" + names[i];
    fxsec->GetObject(id.c_str(), fNuMuCCXSec[i]);
    id = "h_nutau_cc_" + names[i];
    fxsec->GetObject(id.c_str(), fNuTauCCXSec[i]);
    id = "h_nue_cc_" + names[i];
    fxsec->GetObject(id.c_str(),fNueCCXSec[i]);
    
    id = "h_numubar_cc_" + names[i];
    fxsec->GetObject(id.c_str(), fNuMuBarCCXSec[i]);
    id = "h_nutaubar_cc_" + names[i];
    fxsec->GetObject(id.c_str(), fNuTauBarCCXSec[i]);
    id = "h_nuebar_cc_" + names[i];
    fxsec->GetObject(id.c_str(),fNueBarCCXSec[i]);
  }
  
  XSecWeight.clear();
  XSecWeight_NuBar.clear();
  
  TH1D *hnue,*hnutau;
  hnue = (TH1D*)fNueCCXSec[0]->Clone("hnue");
  hnue->Divide(fNuMuCCXSec[0]);
  hnutau = (TH1D*)fNuTauCCXSec[0]->Clone("hnutau");
  hnutau->Divide(fNuMuCCXSec[0]);
  
  TH1D *hnuebar,*hnutaubar;
  hnuebar = (TH1D*)fNueBarCCXSec[0]->Clone("hnue");
  hnuebar->Divide(fNuMuBarCCXSec[0]);
  hnutaubar = (TH1D*)fNuTauBarCCXSec[0]->Clone("hnutau");
  hnutaubar->Divide(fNuMuBarCCXSec[0]);
  
  //just in case binning of xsec hists is not the same as the true binning in this class:
  //this works as if the range is different, as long as the bin width is the same
  
  int n = hnue->GetNbinsX();
  int i,k;
  double E;
  double *t = new double[nTrue+1];
  for(i=0;i<nTrue+1;i++)
  {
    t[i] = TrueEdges.at(i);
  }
  
  TH1D *XSecWeight_NueCC = new TH1D("XSecWeight_NueCC","",nTrue,t);
  TH1D *XSecWeight_NuTauCC = new TH1D("XSecWeight_NuTauCC","",nTrue,t);
  TH1D *XSecWeight_NueBarCC = new TH1D("XSecWeight_NueBarCC","",nTrue,t);
  TH1D *XSecWeight_NuTauBarCC = new TH1D("XSecWeight_NuTauBarCC","",nTrue,t);
  
  if((XSecWeight_NueCC->GetXaxis()->GetBinWidth(1)-hnue->GetXaxis()->GetBinWidth(1))>1e-5)
  {
    cout<<"Warning: the set true energy bin width is not the same as that in the input cross section file"<<endl;
  }
  
  for(i=0;i<nTrue;i++)
  {
    E = XSecWeight_NueCC->GetBinCenter(i+1);
    for(k=0;k<n;k++)
    {
      if(E>=hnue->GetXaxis()->GetBinLowEdge(k+1) && E<hnue->GetXaxis()->GetBinUpEdge(k+1))
      {
	XSecWeight_NueCC->SetBinContent(i+1,hnue->GetBinContent(k+1));
	XSecWeight_NuTauCC->SetBinContent(i+1,hnutau->GetBinContent(k+1));
	break;
      }
    }
  }
  
  XSecWeight[Background::kNueCC] = XSecWeight_NueCC;
  XSecWeight[Background::kNuTauCC] = XSecWeight_NuTauCC;
  
  for(i=0;i<nTrue;i++)
  {
    E = XSecWeight_NueBarCC->GetBinCenter(i+1);
    for(k=0;k<n;k++)
    {
      if(E>=hnuebar->GetXaxis()->GetBinLowEdge(k+1) && E<hnuebar->GetXaxis()->GetBinUpEdge(k+1))
      {
	XSecWeight_NueBarCC->SetBinContent(i+1,hnuebar->GetBinContent(k+1));
	XSecWeight_NuTauBarCC->SetBinContent(i+1,hnutaubar->GetBinContent(k+1));
	break;
      }
    }
  }
  
  XSecWeight_NuBar[Background::kNueCC] = XSecWeight_NueBarCC;
  XSecWeight_NuBar[Background::kNuTauCC] = XSecWeight_NuTauBarCC;
  
  return;
}
void Extrapolate2D::ReadFiles()
{
  ReadNDDataFile();
  cout << "ND Data File" << endl;
  ReadMREFile();
  cout << "MRE Data File" << endl;
  ReadFNFile();
  cout << "FN File" << endl;
  ReadXSecFile();
  cout << "Xsec File" << endl;
  
  return;
}
void Extrapolate2D::ReadExtrap()
{
  Pred_3D.clear();
  Pred_Nu_3D.clear();
  Pred_NuBar_3D.clear();
  Pred.clear();
  FNRatio.clear();
  Pred_TotalBkgd = NULL;
  Pred_NueCCSignal = NULL;
  Pred_NueBarCCSignal = NULL;
  Pred_CC_Fid = NULL;
  Pred_NuMuCC_Fid = NULL;
  Pred_NuMuBarCC_Fid = NULL;
  Pred_TotalBkgd_VsBinNumber = NULL;
  Pred_Signal_VsBinNumber = NULL;
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(inFileName.c_str())))
  {
    cout<<"Failure to read extrap input file."<<endl;
    ReadError = true;
    return;
  }
  
  TFile *f = new TFile(gSystem->ExpandPathName(inFileName.c_str()),"READ");
  
  double npot,fpot;
  TTree *tree = (TTree*)f->Get("paramtree");
  tree->SetBranchAddress("nearPOT",&npot);
  tree->SetBranchAddress("farPOT",&fpot);
  tree->SetBranchAddress("Theta12",&Theta12);
  tree->SetBranchAddress("Theta13",&Theta13);
  tree->SetBranchAddress("Theta23",&Theta23);
  tree->SetBranchAddress("DeltaMSq23",&DeltaMSq23);
  tree->SetBranchAddress("DeltaMSq12",&DeltaMSq12);
  tree->SetBranchAddress("DeltaCP",&DeltaCP);
//   //NSI
   tree->SetBranchAddress("Eps_ee",&Eps_ee);
   tree->SetBranchAddress("Eps_emu",&Eps_emu);
   tree->SetBranchAddress("Eps_etau",&Eps_etau);
   tree->SetBranchAddress("Eps_mumu",&Eps_mumu);
   tree->SetBranchAddress("Eps_mutau",&Eps_mutau);
   tree->SetBranchAddress("Eps_tautau",&Eps_tautau);
// 
   tree->SetBranchAddress("Delta_emu",&Delta_emu);
   tree->SetBranchAddress("Delta_etau",&Delta_etau);
   tree->SetBranchAddress("Delta_mutau",&Delta_mutau);
//   //LSND
   tree->SetBranchAddress("Theta14",&Theta14);
   tree->SetBranchAddress("Theta24",&Theta24);
   tree->SetBranchAddress("Theta34",&Theta34);
   tree->SetBranchAddress("Dm41",&Dm41);
   tree->SetBranchAddress("Delta14",&Delta14);
   tree->SetBranchAddress("Delta24",&Delta24);  
  //
  
   /*
  //these parameters not incorporated into extrap file yet
  Eps_ee=0;
  Eps_emu=0;
  Eps_etau=0;
  Eps_mumu=0;
  Eps_mutau=0;
  Eps_tautau=0;
  
  Delta_emu=0;
  Delta_etau=0;
  Delta_mutau=0;
  
  Theta14=0;
  Theta24=0;
  Theta34=0;
  Dm41=0;
  Delta14=0;
  Delta24=0;
*/  

  tree->GetEntry(0);
  
  SetNearPOT(npot);
  SetFarPOT(fpot);
  
  osc.SetOscParam(OscPar::kTh12,Theta12);
  osc.SetOscParam(OscPar::kTh13,Theta13);
  osc.SetOscParam(OscPar::kTh23,Theta23);
  osc.SetOscParam(OscPar::kDeltaM23,DeltaMSq23);
  osc.SetOscParam(OscPar::kDeltaM12,DeltaMSq12);
  osc.SetOscParam(OscPar::kDelta,DeltaCP);
  //NSI
  osc.SetOscParam(OscPar::kEps_ee,Eps_ee);
  osc.SetOscParam(OscPar::kEps_emu,Eps_emu);
  osc.SetOscParam(OscPar::kEps_etau,Eps_etau);
  osc.SetOscParam(OscPar::kEps_mumu,Eps_mumu);
  osc.SetOscParam(OscPar::kEps_mutau,Eps_mutau);
  osc.SetOscParam(OscPar::kEps_tautau,Eps_tautau);

  osc.SetOscParam(OscPar::kDelta_emu,Delta_emu);
  osc.SetOscParam(OscPar::kDelta_etau,Delta_etau);
  osc.SetOscParam(OscPar::kDelta_mutau,Delta_mutau);
  //LSND
  osc.SetOscParam(OscPar::kTh14,Theta14);
  osc.SetOscParam(OscPar::kTh24,Theta24);
  osc.SetOscParam(OscPar::kTh34,Theta34);
  osc.SetOscParam(OscPar::kDm41,Dm41);
  osc.SetOscParam(OscPar::kDelta14,Delta14);
  osc.SetOscParam(OscPar::kDelta24,Delta24);
  //
  SetNominalOscProb();
  
  TH2D *h2;
  Background::Background_t bg;
  
  for(unsigned int i=0;i<FDComponents.size();i++)
  {
    bg = FDComponents[i];
    string st = string(Background::AsString(bg));
    
    h2 = (TH2D*)f->Get(Form("Pred_%s",st.c_str()));
    Pred[bg] = h2;
  }
  
  Pred_3D[qNC] = (TH3D*)f->Get("Pred_3D_NC");
  Pred_3D[qNuMuToNuMu] = (TH3D*)f->Get("Pred_3D_NuMuToNuMu");
  Pred_3D[qBNueToNuMu] = (TH3D*)f->Get("Pred_3D_BNueToNuMu");
  Pred_3D[qBNueToBNue] = (TH3D*)f->Get("Pred_3D_BNueToBNue");
  Pred_3D[qNuMuToNue] = (TH3D*)f->Get("Pred_3D_NuMuToNue");
  Pred_3D[qNuMuToNuTau] = (TH3D*)f->Get("Pred_3D_NuMuToNuTau");
  Pred_3D[qBNueToNuTau] = (TH3D*)f->Get("Pred_3D_BNueToNuTau");
  
  if(f->Read("Pred_Nu_3D_NC")!=0)
  {
    UseSeparateNuNuBar=true;
  }
  else
  {
    cout<<"Warning: Pred_Nu_3D or Pred_NuBar_3D histograms not found in input file, so nus and nubars will be combined."<<endl;
  }
  
  if(UseSeparateNuNuBar)
  {
    Pred_Nu_3D[qNC] = (TH3D*)f->Get("Pred_Nu_3D_NC");
    Pred_Nu_3D[qNuMuToNuMu] = (TH3D*)f->Get("Pred_Nu_3D_NuMuToNuMu");
    Pred_Nu_3D[qBNueToNuMu] = (TH3D*)f->Get("Pred_Nu_3D_BNueToNuMu");
    Pred_Nu_3D[qBNueToBNue] = (TH3D*)f->Get("Pred_Nu_3D_BNueToBNue");
    Pred_Nu_3D[qNuMuToNue] = (TH3D*)f->Get("Pred_Nu_3D_NuMuToNue");
    Pred_Nu_3D[qNuMuToNuTau] = (TH3D*)f->Get("Pred_Nu_3D_NuMuToNuTau");
    Pred_Nu_3D[qBNueToNuTau] = (TH3D*)f->Get("Pred_Nu_3D_BNueToNuTau");
    
    Pred_NuBar_3D[qNC] = (TH3D*)f->Get("Pred_NuBar_3D_NC");
    Pred_NuBar_3D[qNuMuToNuMu] = (TH3D*)f->Get("Pred_NuBar_3D_NuMuToNuMu");
    Pred_NuBar_3D[qBNueToNuMu] = (TH3D*)f->Get("Pred_NuBar_3D_BNueToNuMu");
    Pred_NuBar_3D[qBNueToBNue] = (TH3D*)f->Get("Pred_NuBar_3D_BNueToBNue");
    Pred_NuBar_3D[qNuMuToNue] = (TH3D*)f->Get("Pred_NuBar_3D_NuMuToNue");
    Pred_NuBar_3D[qNuMuToNuTau] = (TH3D*)f->Get("Pred_NuBar_3D_NuMuToNuTau");
    Pred_NuBar_3D[qBNueToNuTau] = (TH3D*)f->Get("Pred_NuBar_3D_BNueToNuTau");
  }
  
  FNRatio[Background::kNC] = (TH2D*)f->Get("FNRatio_NC");
  FNRatio[Background::kNuMuCC] = (TH2D*)f->Get("FNRatio_NuMuCC");
  FNRatio[Background::kBNueCC] = (TH2D*)f->Get("FNRatio_BNueCC");
  
  Pred_TotalBkgd = (TH2D*)f->Get("Pred_TotalBkgd");
  if(UseSeparateNuNuBar)
  {
    Pred_NueCCSignal = (TH2D*)f->Get("Pred_NueCCSignal");
    Pred_NueBarCCSignal = (TH2D*)f->Get("Pred_NueBarCCSignal");
  }
  
  Pred_CC_Fid = (TH1D*)f->Get("Pred_CC_Fid");
  if(UseSeparateNuNuBar)
  {
    Pred_NuMuCC_Fid = (TH1D*)f->Get("Pred_NuMuCC_Fid");
    Pred_NuMuBarCC_Fid = (TH1D*)f->Get("Pred_NuMuBarCC_Fid");
  }
  
  NDData[Background::kNC] = (TH2D*)f->Get(Form("NDData_NC_%s_Run%i",PID_NDData.c_str(),RunPeriod));
  NDData[Background::kNuMuCC] = (TH2D*)f->Get(Form("NDData_NuMuCC_%s_Run%i",PID_NDData.c_str(),RunPeriod));
  NDData[Background::kBNueCC] = (TH2D*)f->Get(Form("NDData_BNueCC_%s_Run%i",PID_NDData.c_str(),RunPeriod));
  NDData_Reco_CClike = (TH1D*)f->Get("Reco_CClike");
  if(UseSeparateNuNuBar)
  {
    NDData_Reco_CClike_Neg = (TH1D*)f->Get("Reco_CClike_Neg");
    NDData_Reco_CClike_Pos = (TH1D*)f->Get("Reco_CClike_Pos");
  }
  
  Int_t j;
  Int_t np = Pred_3D[qNC]->GetNbinsX();
  Int_t nr = Pred_3D[qNC]->GetNbinsY();
  Int_t nt = Pred_3D[qNC]->GetNbinsZ();
  double *p = new double[np+1];
  double *r = new double[nr+1];
  double *t = new double[nt+1];
  for(j=0;j<np+1;j++)
  {
    p[j] = Pred_3D[qNC]->GetXaxis()->GetBinLowEdge(j+1);
  }
  for(j=0;j<nr+1;j++)
  {
    r[j] = Pred_3D[qNC]->GetYaxis()->GetBinLowEdge(j+1);
  }
  for(j=0;j<nt+1;j++)
  {
    t[j] = Pred_3D[qNC]->GetZaxis()->GetBinLowEdge(j+1);
  }
  SetPIDBins(np,p);
  SetRecoBins(nr,r);
  SetTrueBins(nt,t);
  
  int totbin = nReco*nPID;
  Pred_TotalBkgd_VsBinNumber = new TH1D("Pred_TotalBkgd_VsBinNumber","",totbin,-0.5,totbin-0.5);
  Pred_Signal_VsBinNumber = new TH1D("Pred_Signal_VsBinNumber","",totbin,-0.5,totbin-0.5);
  Set1DPredHists();
  
  FD_TrueVsRecoVsPID_NC_NueFrac = (TH3D*)f->Get(Form("FD_TrueVsRecoVs%s_NC_NueFrac",PID.c_str()));
  FD_TrueVsRecoVsPID_NC_NueBarFrac = (TH3D*)f->Get(Form("FD_TrueVsRecoVs%s_NC_NueBarFrac",PID.c_str()));
  FD_TrueVsRecoVsPID_NC_NuMuFrac = (TH3D*)f->Get(Form("FD_TrueVsRecoVs%s_NC_NuMuFrac",PID.c_str()));
  FD_TrueVsRecoVsPID_NC_NuMuBarFrac = (TH3D*)f->Get(Form("FD_TrueVsRecoVs%s_NC_NuMuBarFrac",PID.c_str()));
  
  if(PrintResult)
  {
    cout.precision(4);
    cout<<fixed;
    cout<<"sin^2(2Th13) = "<<TMath::Sin(2.*osc.GetOscParam(OscPar::kTh13))*TMath::Sin(2.*osc.GetOscParam(OscPar::kTh13))<<endl;
    cout<<"NC: "<<Pred[Background::kNC]->Integral()<<endl;
    cout<<"NuMuCC: "<<Pred[Background::kNuMuCC]->Integral()<<endl;
    cout<<"BNueCC: "<<Pred[Background::kBNueCC]->Integral()<<endl;
    cout<<"NuTauCC: "<<Pred[Background::kNuTauCC]->Integral()<<endl;
    cout<<"Total: "<<Pred_TotalBkgd->Integral()<<endl;
    cout<<"NueCC: "<<Pred[Background::kNueCC]->Integral()<<endl;
  }
  
  delete [] p;
  delete [] r;
  delete [] t;
  
  return;
}
void Extrapolate2D::GetPrediction()
{
  if(Init)
  {
//     cout<<"You only need to call GetPrediction() once!"<<endl;
    return;
  }
  
  if(ReadExtrapFromFile)
  {
    ReadExtrap();
    if(ReadError) return;
    
    Init = true;
    
    return;
  }
  
  SetupPredHists();
  cout << "Setup Hists" << endl;
  ReadFiles();
  cout << "Read Files" << endl;
  if(ReadError) return;
  GetNoOscCCPrediction();
  cout << "Got No Osc Pred" << endl;
  
  for(unsigned int i=0;i<FDComponents.size();i++)
  {
    if(FDComponents[i]==Background::kNuMuCC || FDComponents[i]==Background::kNC)
    {
      FarNearPred(FDComponents[i]);
    }
    else if(FDComponents[i]==Background::kNueCC || FDComponents[i]==Background::kNuTauCC)
    {
      AppearancePred(FDComponents[i]);
    }
    else if(FDComponents[i]==Background::kBNueCC)
    {
      BeamNuePred();
    }
    else
    {
      cout<<"Unknown component"<<endl;
    }
    if(FDComponents[i]!=Background::kNueCC)
    {
      Pred_TotalBkgd->Add(Pred[FDComponents[i]]);
    }
  }
  cout << "Survived for loop" << endl;  

  Set1DPredHists();
  cout << "Set 1D Pred Hists" << endl;
  if(PrintResult)
  {
    cout.precision(4);
    cout<<fixed;
    cout<<"sin^2(2Th13) = "<<TMath::Sin(2.*osc.GetOscParam(OscPar::kTh13))*TMath::Sin(2.*osc.GetOscParam(OscPar::kTh13))<<endl;
    cout<<"NC: "<<Pred[Background::kNC]->Integral()<<endl;
    cout<<"NuMuCC: "<<Pred[Background::kNuMuCC]->Integral()<<endl;
    cout<<"BNueCC: "<<Pred[Background::kBNueCC]->Integral()<<endl;
    cout<<"NuTauCC: "<<Pred[Background::kNuTauCC]->Integral()<<endl;
    cout<<"Total: "<<Pred_TotalBkgd->Integral()<<endl;
    cout<<"NueCC: "<<Pred[Background::kNueCC]->Integral()<<endl;
  }
  
  if(WriteOutput) WriteToFile();
  
  Init = true;
  
  return;
}
void Extrapolate2D::FarNearPred(Background::Background_t bg)
{
  if(bg!=Background::kNuMuCC && bg!=Background::kNC && bg!=Background::kBNueCC)
  {
    cout<<"Can only use FarNearPred() for cc numu, nc and beam nues"<<endl;
    return;
  }
  
  int ip,ir,it;
  double temp,sum;
  
  for(ip=0;ip<nPID;ip++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      sum=0;
      for(it=0;it<nTrue;it++)
      {
	if(bg==Background::kNC)
	{
	  temp = FD_TrueVsRecoVsPID[qNC]->GetBinContent(ip+1,ir+1,it+1)*ND_DataOverMC_RecoVsPID[bg]->GetBinContent(ip+1,ir+1);
	  Pred_3D[qNC]->SetBinContent(ip+1,ir+1,it+1,temp);
	  sum+=temp;
	  
	  if(UseSeparateNuNuBar)
	  {
	    temp = FD_Nu_TrueVsRecoVsPID[qNC]->GetBinContent(ip+1,ir+1,it+1)*ND_DataOverMC_RecoVsPID[bg]->GetBinContent(ip+1,ir+1);
	    Pred_Nu_3D[qNC]->SetBinContent(ip+1,ir+1,it+1,temp);
	    
	    temp = FD_NuBar_TrueVsRecoVsPID[qNC]->GetBinContent(ip+1,ir+1,it+1)*ND_DataOverMC_RecoVsPID[bg]->GetBinContent(ip+1,ir+1);
	    Pred_NuBar_3D[qNC]->SetBinContent(ip+1,ir+1,it+1,temp);
	  }
	}
	if(bg==Background::kNuMuCC)
	{
	  temp = FD_TrueVsRecoVsPID[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*ND_DataOverMC_RecoVsPID[bg]->GetBinContent(ip+1,ir+1);
	  Pred_3D[qNuMuToNuMu]->SetBinContent(ip+1,ir+1,it+1,temp);
	  sum+=temp;
	  
	  temp = FD_TrueVsRecoVsPID[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*ND_DataOverMC_RecoVsPID[bg]->GetBinContent(ip+1,ir+1);
	  Pred_3D[qBNueToNuMu]->SetBinContent(ip+1,ir+1,it+1,temp);
	  sum+=temp;
	  
	  if(UseSeparateNuNuBar)
	  {
	    temp = FD_Nu_TrueVsRecoVsPID[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*ND_DataOverMC_RecoVsPID[bg]->GetBinContent(ip+1,ir+1);
	    Pred_Nu_3D[qNuMuToNuMu]->SetBinContent(ip+1,ir+1,it+1,temp);
	    
	    temp = FD_Nu_TrueVsRecoVsPID[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*ND_DataOverMC_RecoVsPID[bg]->GetBinContent(ip+1,ir+1);
	    Pred_Nu_3D[qBNueToNuMu]->SetBinContent(ip+1,ir+1,it+1,temp);
	    
	    temp = FD_NuBar_TrueVsRecoVsPID[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*ND_DataOverMC_RecoVsPID[bg]->GetBinContent(ip+1,ir+1);
	    Pred_NuBar_3D[qNuMuToNuMu]->SetBinContent(ip+1,ir+1,it+1,temp);
	    
	    temp = FD_NuBar_TrueVsRecoVsPID[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*ND_DataOverMC_RecoVsPID[bg]->GetBinContent(ip+1,ir+1);
	    Pred_NuBar_3D[qBNueToNuMu]->SetBinContent(ip+1,ir+1,it+1,temp);
	  }
	}
	if(bg==Background::kBNueCC)
	{
	  temp = FD_TrueVsRecoVsPID[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*ND_DataOverMC_RecoVsPID[bg]->GetBinContent(ip+1,ir+1);
	  Pred_3D[qBNueToBNue]->SetBinContent(ip+1,ir+1,it+1,temp);
	  sum+=temp;
	  
	  if(UseSeparateNuNuBar)
	  {
	    temp = FD_Nu_TrueVsRecoVsPID[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*ND_DataOverMC_RecoVsPID[bg]->GetBinContent(ip+1,ir+1);
	    Pred_Nu_3D[qBNueToBNue]->SetBinContent(ip+1,ir+1,it+1,temp);
	    
	    temp = FD_NuBar_TrueVsRecoVsPID[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*ND_DataOverMC_RecoVsPID[bg]->GetBinContent(ip+1,ir+1);
	    Pred_NuBar_3D[qBNueToBNue]->SetBinContent(ip+1,ir+1,it+1,temp);
	  }
	}
	
      }
      Pred[bg]->SetBinContent(ip+1,ir+1,sum);
    }
  }
  
  return;
}
void Extrapolate2D::AppearancePred(Background::Background_t bg)
{
  if(bg!=Background::kNueCC && bg!=Background::kNuTauCC)
  {
    cout<<"Can only use AppearancePred() for taus and signal nues"<<endl;
    return;
  }
  
  int ip,ir,it;
  double temp,sum,sumnu,sumnubar;
  double nu,nubar;
  for(ip=0;ip<nPID;ip++)
  {
    for(ir=0;ir<nReco;ir++)
    {
      sum=0;
      sumnu=0;
      sumnubar=0;
      
      for(it=0;it<nTrue;it++)
      {
        if(UseSeparateNuNuBar)
        {
          nu = Pred_NuMuCC_Fid->GetBinContent(it+1)*XSecWeight[bg]->GetBinContent(it+1)*FD_True2Reco_Fid[bg]->GetBinContent(ir+1,it+1)*FD_Eff[bg]->GetBinContent(ip+1,ir+1);
          
          nubar = Pred_NuMuBarCC_Fid->GetBinContent(it+1)*XSecWeight_NuBar[bg]->GetBinContent(it+1)*FD_True2Reco_NuBar_Fid[bg]->GetBinContent(ir+1,it+1)*FD_Eff_NuBar[bg]->GetBinContent(ip+1,ir+1);
        }
        else
        {
	  temp = Pred_CC_Fid->GetBinContent(it+1)*XSecWeight[bg]->GetBinContent(it+1)*FD_True2Reco_Fid[bg]->GetBinContent(ir+1,it+1)*FD_Eff[bg]->GetBinContent(ip+1,ir+1);
        }
	
	if(bg==Background::kNueCC)
	{
	  if(UseSeparateNuNuBar)
	  {
            nu*=NominalOscProb[qNuMuToNue]->GetBinContent(it+1);
	    Pred_Nu_3D[qNuMuToNue]->SetBinContent(ip+1,ir+1,it+1,nu);
	    
            nubar*=NominalOscProb_NuBar[qNuMuToNue]->GetBinContent(it+1);
	    Pred_NuBar_3D[qNuMuToNue]->SetBinContent(ip+1,ir+1,it+1,nubar);
	    
	    Pred_3D[qNuMuToNue]->SetBinContent(ip+1,ir+1,it+1,nu+nubar);
	  }
	  else
          {
            temp*=NominalOscProb[qNuMuToNue]->GetBinContent(it+1);
            Pred_3D[qNuMuToNue]->SetBinContent(ip+1,ir+1,it+1,temp);
          }
	}
	if(bg==Background::kNuTauCC)
	{
	  if(UseSeparateNuNuBar)
	  {
            nu*=NominalOscProb[qNuMuToNuTau]->GetBinContent(it+1);
	    Pred_Nu_3D[qNuMuToNuTau]->SetBinContent(ip+1,ir+1,it+1,nu);

            nubar*=NominalOscProb_NuBar[qNuMuToNuTau]->GetBinContent(it+1);
	    Pred_NuBar_3D[qNuMuToNuTau]->SetBinContent(ip+1,ir+1,it+1,nubar);
	    
	    Pred_3D[qNuMuToNuTau]->SetBinContent(ip+1,ir+1,it+1,nu+nubar);
	    
	    nu += FD_Nu_TrueVsRecoVsPID[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1);
	    Pred_Nu_3D[qBNueToNuTau]->SetBinContent(ip+1,ir+1,it+1,FD_Nu_TrueVsRecoVsPID[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1));
	    
	    nubar += FD_NuBar_TrueVsRecoVsPID[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1);
	    Pred_NuBar_3D[qBNueToNuTau]->SetBinContent(ip+1,ir+1,it+1,FD_NuBar_TrueVsRecoVsPID[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1));
	    
	    Pred_3D[qBNueToNuTau]->SetBinContent(ip+1,ir+1,it+1,FD_Nu_TrueVsRecoVsPID[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1)+FD_NuBar_TrueVsRecoVsPID[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1));
	  }
	  else
          {
            temp*=NominalOscProb[qNuMuToNuTau]->GetBinContent(it+1);
            Pred_3D[qNuMuToNuTau]->SetBinContent(ip+1,ir+1,it+1,temp);
            
            temp += FD_TrueVsRecoVsPID[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1);//adds in contributino fron nue->nutau
            Pred_3D[qBNueToNuTau]->SetBinContent(ip+1,ir+1,it+1,FD_TrueVsRecoVsPID[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1));
          }
	  
	}
	
	if(UseSeparateNuNuBar)
	{
	  sum += (nu+nubar);
          sumnu += nu;
          sumnubar += nubar;
	}
	else 
	{
	  sum+=temp;
	}
      }
      Pred[bg]->SetBinContent(ip+1,ir+1,sum);
      if(UseSeparateNuNuBar && bg==Background::kNueCC)
      {
        Pred_NueCCSignal->SetBinContent(ip+1,ir+1,sumnu);
        Pred_NueBarCCSignal->SetBinContent(ip+1,ir+1,sumnubar);
      }
    }
  }
  
  
  return;
}
void Extrapolate2D::BeamNuePred()
{
  Background::Background_t bg = Background::kBNueCC;
  int ip,ir,it;
  double temp,sum;
  
  if(FNforBeamNue)
  {
    cout<<"using F/N method for beam nues (based on ND data)"<<endl;
    FarNearPred(bg);
  }
  else
  {
    cout<<"using far MC method for beam nues"<<endl;
    
    for(ip=0;ip<nPID;ip++)
    {
      for(ir=0;ir<nReco;ir++)
      {
        sum=0;
        for(it=0;it<nTrue;it++)
        {
	  temp = FD_TrueVsRecoVsPID[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1);
	  Pred_3D[qBNueToBNue]->SetBinContent(ip+1,ir+1,it+1,temp);
	  sum+=temp;
	  
	  if(UseSeparateNuNuBar)
	  {
	    Pred_Nu_3D[qBNueToBNue]->SetBinContent(ip+1,ir+1,it+1,temp);
	    Pred_NuBar_3D[qBNueToBNue]->SetBinContent(ip+1,ir+1,it+1,temp);
	  }
        }
        Pred[bg]->SetBinContent(ip+1,ir+1,sum);
      }
    }
  }
  
  return;
}
void Extrapolate2D::WriteToFile()
{
  ExtrapFile->cd();
  
  TTree *paramtree = new TTree("paramtree","paramtree");
  paramtree->Branch("nearPOT",&nPOTNear,"nearPOT/D");
  paramtree->Branch("farPOT",&nPOTFar,"farPOT/D");
  paramtree->Branch("Theta12",&Theta12,"Theta12/D");
  paramtree->Branch("Theta13",&Theta13,"Theta13/D");
  paramtree->Branch("Theta23",&Theta23,"Theta23/D");
  paramtree->Branch("DeltaMSq23",&DeltaMSq23,"DeltaMSq23/D");
  paramtree->Branch("DeltaMSq12",&DeltaMSq12,"DeltaMSq12/D");
  paramtree->Branch("DeltaCP",&DeltaCP,"DeltaCP/D");
  paramtree->Branch("Eps_ee",&Eps_ee);
  paramtree->Branch("Eps_emu",&Eps_emu);
  paramtree->Branch("Eps_etau",&Eps_etau);
  paramtree->Branch("Eps_mumu",&Eps_mumu);
  paramtree->Branch("Eps_mutau",&Eps_mutau);
  paramtree->Branch("Eps_tautau",&Eps_tautau);
  paramtree->Branch("Delta_emu",&Delta_emu);
  paramtree->Branch("Delta_etau",&Delta_etau);
  paramtree->Branch("Delta_mutau",&Delta_mutau);
  paramtree->Branch("Theta14",&Theta14,"Theta14/D");
  paramtree->Branch("Theta24",&Theta24,"Theta24/D");
  paramtree->Branch("Theta34",&Theta34,"Theta34/D");
  paramtree->Branch("Dm41",&Dm41,"Dm41/D");
  paramtree->Branch("Delta14",&Delta14);
  paramtree->Branch("Delta24",&Delta24);


  paramtree->Fill();
  
  for(unsigned int i=0;i<FDComponents.size();i++)
  {
    Pred[FDComponents[i]]->Write();
  }
  
  Pred_3D[qNC]->Write();
  Pred_3D[qNuMuToNuMu]->Write();
  Pred_3D[qBNueToNuMu]->Write();
  Pred_3D[qBNueToBNue]->Write();
  Pred_3D[qNuMuToNue]->Write();
  Pred_3D[qNuMuToNuTau]->Write();
  Pred_3D[qBNueToNuTau]->Write();
  
  if(UseSeparateNuNuBar)
  {
    Pred_Nu_3D[qNC]->Write();
    Pred_Nu_3D[qNuMuToNuMu]->Write();
    Pred_Nu_3D[qBNueToNuMu]->Write();
    Pred_Nu_3D[qBNueToBNue]->Write();
    Pred_Nu_3D[qNuMuToNue]->Write();
    Pred_Nu_3D[qNuMuToNuTau]->Write();
    Pred_Nu_3D[qBNueToNuTau]->Write();
    
    Pred_NuBar_3D[qNC]->Write();
    Pred_NuBar_3D[qNuMuToNuMu]->Write();
    Pred_NuBar_3D[qBNueToNuMu]->Write();
    Pred_NuBar_3D[qBNueToBNue]->Write();
    Pred_NuBar_3D[qNuMuToNue]->Write();
    Pred_NuBar_3D[qNuMuToNuTau]->Write();
    Pred_NuBar_3D[qBNueToNuTau]->Write();
  }
  
  Pred_TotalBkgd->Write();
  if(UseSeparateNuNuBar)
  {
    Pred_NueCCSignal->Write();
    Pred_NueBarCCSignal->Write();
  }
  
  Pred_CC_Fid->Write();
  if(UseSeparateNuNuBar)
  {
    Pred_NuMuCC_Fid->Write();
    Pred_NuMuBarCC_Fid->Write();
  }
  FNRatio[Background::kNC]->Write();
  FNRatio[Background::kNuMuCC]->Write();
  FNRatio[Background::kBNueCC]->Write();
  
  NDData[Background::kNC]->Write(Form("NDData_NC_%s_Run%i",PID_NDData.c_str(),RunPeriod));
  NDData[Background::kNuMuCC]->Write(Form("NDData_NuMuCC_%s_Run%i",PID_NDData.c_str(),RunPeriod));
  NDData[Background::kBNueCC]->Write(Form("NDData_BNueCC_%s_Run%i",PID_NDData.c_str(),RunPeriod));
  NDData_Reco_CClike->Write();
  if(UseSeparateNuNuBar)
  {
    NDData_Reco_CClike_Neg->Write();
    NDData_Reco_CClike_Pos->Write();
  }
  
  if( OscMethod == 3){
  FD_TrueVsRecoVsPID_NC_NueFrac->Write();
  FD_TrueVsRecoVsPID_NC_NueBarFrac->Write();
  FD_TrueVsRecoVsPID_NC_NuMuFrac->Write();
  FD_TrueVsRecoVsPID_NC_NuMuBarFrac->Write();
  }
  paramtree->Write();
  
  //ExtrapFile->Close();
  
  return;
}
void Extrapolate2D::OscillatePrediction()
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
    cout.precision(4);
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
void Extrapolate2D::OscillatePrediction_SepNuNuBar()
{
  //will use Pred_Nu_3D and Pred_NuBar_3D to reweight Pred for oscillation
  
  for(unsigned int i=0;i<FDComponents.size();i++)
  {
    Pred[FDComponents[i]]->Reset();
  }
  Pred_TotalBkgd->Reset();
  Pred_NueCCSignal->Reset();
  Pred_NueBarCCSignal->Reset();
  
  int it,ir,ip;
  double E;
  double temp,sum,sumnu,sumnubar;
  
  for(unsigned int i=0;i<FDComponents.size();i++)
  {
    Background::Background_t bg = FDComponents[i];
    for(ip=0;ip<nPID;ip++)
    {
      for(ir=0;ir<nReco;ir++)
      {
	sum=0;
        sumnu=0;
        sumnubar=0;
        for(it=0;it<nTrue;it++)
        {
	  E = Pred_3D[qNC]->GetZaxis()->GetBinCenter(it+1);
	  
	  temp=0;
	  
	  if(bg==Background::kNuMuCC)
	  {
	    //nu
	    osc.SetOscParam(OscPar::kNuAntiNu,1);
	    if(NominalOscProb[qNuMuToNuMu]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp += (Pred_Nu_3D[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(14,14,E)/NominalOscProb[qNuMuToNuMu]->GetBinContent(it+1));}
	      else if(OscMethod == 2){temp += (Pred_Nu_3D[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(14,14,E)/NominalOscProb[qNuMuToNuMu]->GetBinContent(it+1));}
	      else if(OscMethod == 1){temp += (Pred_Nu_3D[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(14,14,E)/NominalOscProb[qNuMuToNuMu]->GetBinContent(it+1));}
	      else{temp += (Pred_Nu_3D[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.MuToMu(E)/NominalOscProb[qNuMuToNuMu]->GetBinContent(it+1));}
	    }
	    if(NominalOscProb[qBNueToNuMu]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp += (Pred_Nu_3D[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(14,12,E)/NominalOscProb[qBNueToNuMu]->GetBinContent(it+1));}
	      else if(OscMethod == 2){temp += (Pred_Nu_3D[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(14,12,E)/NominalOscProb[qBNueToNuMu]->GetBinContent(it+1));}
	      else if(OscMethod == 1){temp += (Pred_Nu_3D[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(14,12,E)/NominalOscProb[qBNueToNuMu]->GetBinContent(it+1));}
	      else{temp += (Pred_Nu_3D[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.ElecToMu(E)/NominalOscProb[qBNueToNuMu]->GetBinContent(it+1));}
	    }
	    //nubar
	    osc.SetOscParam(OscPar::kNuAntiNu,-1);
	    if(NominalOscProb_NuBar[qNuMuToNuMu]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp += (Pred_NuBar_3D[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(-14,-14,E)/NominalOscProb_NuBar[qNuMuToNuMu]->GetBinContent(it+1));}
	      else if(OscMethod == 2){temp += (Pred_NuBar_3D[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(-14,-14,E)/NominalOscProb_NuBar[qNuMuToNuMu]->GetBinContent(it+1));}
	      else if(OscMethod == 1){temp += (Pred_NuBar_3D[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(-14,-14,E)/NominalOscProb_NuBar[qNuMuToNuMu]->GetBinContent(it+1));}
	      else{temp += (Pred_NuBar_3D[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.MuToMu(E)/NominalOscProb_NuBar[qNuMuToNuMu]->GetBinContent(it+1));}
	    }
	    if(NominalOscProb_NuBar[qBNueToNuMu]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp += (Pred_NuBar_3D[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(-14,-12,E)/NominalOscProb_NuBar[qBNueToNuMu]->GetBinContent(it+1));}
	      else if(OscMethod == 2){temp += (Pred_NuBar_3D[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(-14,-12,E)/NominalOscProb_NuBar[qBNueToNuMu]->GetBinContent(it+1));}
	      else if(OscMethod == 1){temp += (Pred_NuBar_3D[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(-14,-12,E)/NominalOscProb_NuBar[qBNueToNuMu]->GetBinContent(it+1));}
	      else{temp += (Pred_NuBar_3D[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.ElecToMu(E)/NominalOscProb_NuBar[qBNueToNuMu]->GetBinContent(it+1));}
	    }
	  }
	  if(bg==Background::kBNueCC)
	  {
	    //nu
	    osc.SetOscParam(OscPar::kNuAntiNu,1);
	    if(NominalOscProb[qBNueToBNue]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp += (Pred_Nu_3D[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(12,12,E)/NominalOscProb[qBNueToBNue]->GetBinContent(it+1));}
	      else if(OscMethod == 2){temp += (Pred_Nu_3D[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(12,12,E)/NominalOscProb[qBNueToBNue]->GetBinContent(it+1));}
	      else if(OscMethod == 1){temp += (Pred_Nu_3D[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(12,12,E)/NominalOscProb[qBNueToBNue]->GetBinContent(it+1));}
	      else{temp += (Pred_Nu_3D[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*osc.ElecToElec(E)/NominalOscProb[qBNueToBNue]->GetBinContent(it+1));}
	    }
	    //nubar
	    osc.SetOscParam(OscPar::kNuAntiNu,-1);
	    if(NominalOscProb_NuBar[qBNueToBNue]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp += (Pred_NuBar_3D[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(-12,-12,E)/NominalOscProb_NuBar[qBNueToBNue]->GetBinContent(it+1));}
	      else if(OscMethod == 2){temp += (Pred_NuBar_3D[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(-12,-12,E)/NominalOscProb_NuBar[qBNueToBNue]->GetBinContent(it+1));}
	      else if(OscMethod == 1){temp += (Pred_NuBar_3D[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(-12,-12,E)/NominalOscProb_NuBar[qBNueToBNue]->GetBinContent(it+1));}
	      else{temp += (Pred_NuBar_3D[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*osc.ElecToElec(E)/NominalOscProb_NuBar[qBNueToBNue]->GetBinContent(it+1));}
	    }
	  }
	  if(bg==Background::kNueCC)
	  {
	    //nu
	    osc.SetOscParam(OscPar::kNuAntiNu,1);
	    if(NominalOscProb[qNuMuToNue]->GetBinContent(it+1)>1e-20)
	    {
              if(OscMethod == 3){
	      temp += (Pred_Nu_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(12,14,E)/NominalOscProb[qNuMuToNue]->GetBinContent(it+1));
              sumnu += (Pred_Nu_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(12,14,E)/NominalOscProb[qNuMuToNue]->GetBinContent(it+1));
	      }
              else if(OscMethod == 2){
	      temp += (Pred_Nu_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(12,14,E)/NominalOscProb[qNuMuToNue]->GetBinContent(it+1));
              sumnu += (Pred_Nu_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(12,14,E)/NominalOscProb[qNuMuToNue]->GetBinContent(it+1));
	      }
              else if(OscMethod == 1){
	      temp += (Pred_Nu_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(12,14,E)/NominalOscProb[qNuMuToNue]->GetBinContent(it+1));
              sumnu += (Pred_Nu_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(12,14,E)/NominalOscProb[qNuMuToNue]->GetBinContent(it+1));
	      }
              else{
	      temp += (Pred_Nu_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.MuToElec(E)/NominalOscProb[qNuMuToNue]->GetBinContent(it+1));
              sumnu += (Pred_Nu_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.MuToElec(E)/NominalOscProb[qNuMuToNue]->GetBinContent(it+1));
	      }
	    }
	     //nubar
	    osc.SetOscParam(OscPar::kNuAntiNu,-1);
	    if(NominalOscProb_NuBar[qNuMuToNue]->GetBinContent(it+1)>1e-20)
	    {
              if(OscMethod == 3){
	      temp += (Pred_NuBar_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(-12,-14,E)/NominalOscProb_NuBar[qNuMuToNue]->GetBinContent(it+1));
              sumnubar += (Pred_NuBar_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(-12,-14,E)/NominalOscProb_NuBar[qNuMuToNue]->GetBinContent(it+1));
	      }
              else if(OscMethod == 2){
	      temp += (Pred_NuBar_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(-12,-14,E)/NominalOscProb_NuBar[qNuMuToNue]->GetBinContent(it+1));
              sumnubar += (Pred_NuBar_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(-12,-14,E)/NominalOscProb_NuBar[qNuMuToNue]->GetBinContent(it+1));
	      }
              else if(OscMethod == 1){
	      temp += (Pred_NuBar_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(-12,-14,E)/NominalOscProb_NuBar[qNuMuToNue]->GetBinContent(it+1));
              sumnubar += (Pred_NuBar_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(-12,-14,E)/NominalOscProb_NuBar[qNuMuToNue]->GetBinContent(it+1));
	      }
              else{
	      temp += (Pred_NuBar_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.MuToElec(E)/NominalOscProb_NuBar[qNuMuToNue]->GetBinContent(it+1));
              sumnubar += (Pred_NuBar_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.MuToElec(E)/NominalOscProb_NuBar[qNuMuToNue]->GetBinContent(it+1));
	      }
	    }
	  }
	  if(bg==Background::kNuTauCC)
	  {
	    //nu
	    osc.SetOscParam(OscPar::kNuAntiNu,1);
	    if(NominalOscProb[qNuMuToNuTau]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp += Pred_Nu_3D[qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(16,14,E)/NominalOscProb[qNuMuToNuTau]->GetBinContent(it+1);}
	      else if(OscMethod == 2){temp += Pred_Nu_3D[qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(16,14,E)/NominalOscProb[qNuMuToNuTau]->GetBinContent(it+1);}
	      else if(OscMethod == 1){temp += Pred_Nu_3D[qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(16,14,E)/NominalOscProb[qNuMuToNuTau]->GetBinContent(it+1);}
	      else{temp += Pred_Nu_3D[qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.MuToTau(E)/NominalOscProb[qNuMuToNuTau]->GetBinContent(it+1);}
	    }
	    if(NominalOscProb[qBNueToNuTau]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp += Pred_Nu_3D[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(16,12,E)/NominalOscProb[qBNueToNuTau]->GetBinContent(it+1);}
	      else if(OscMethod == 2){temp += Pred_Nu_3D[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(16,12,E)/NominalOscProb[qBNueToNuTau]->GetBinContent(it+1);}
	      else if(OscMethod == 1){temp += Pred_Nu_3D[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(16,12,E)/NominalOscProb[qBNueToNuTau]->GetBinContent(it+1);}
	      else{temp += Pred_Nu_3D[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.ElecToTau(E)/NominalOscProb[qBNueToNuTau]->GetBinContent(it+1);}
	    }
	    //nubar
	    osc.SetOscParam(OscPar::kNuAntiNu,-1);
	    if(NominalOscProb_NuBar[qNuMuToNuTau]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp += Pred_NuBar_3D[qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(-16,-14,E)/NominalOscProb_NuBar[qNuMuToNuTau]->GetBinContent(it+1);}
	      else if(OscMethod == 2){temp += Pred_NuBar_3D[qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(-16,-14,E)/NominalOscProb_NuBar[qNuMuToNuTau]->GetBinContent(it+1);}
	      else if(OscMethod == 1){temp += Pred_NuBar_3D[qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(-16,-14,E)/NominalOscProb_NuBar[qNuMuToNuTau]->GetBinContent(it+1);}
	      else{temp += Pred_NuBar_3D[qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.MuToTau(E)/NominalOscProb_NuBar[qNuMuToNuTau]->GetBinContent(it+1);}
	    }
	    if(NominalOscProb_NuBar[qBNueToNuTau]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp += Pred_NuBar_3D[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(-16,-12,E)/NominalOscProb_NuBar[qBNueToNuTau]->GetBinContent(it+1);}
	      else if(OscMethod == 2){temp += Pred_NuBar_3D[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(-16,-12,E)/NominalOscProb_NuBar[qBNueToNuTau]->GetBinContent(it+1);}
	      else if(OscMethod == 1){temp += Pred_NuBar_3D[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(-16,-12,E)/NominalOscProb_NuBar[qBNueToNuTau]->GetBinContent(it+1);}
	      else{temp += Pred_NuBar_3D[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.ElecToTau(E)/NominalOscProb_NuBar[qBNueToNuTau]->GetBinContent(it+1);}
	    }
	  }
	  if(bg==Background::kNC)
	  {
 	    if(OscMethod == 3)
	    {
	      osc.SetOscParam(OscPar::kNuAntiNu,1);
	      temp = FD_TrueVsRecoVsPID_NC_NueFrac->GetBinContent(ip+1,ir+1,it+1)*(1-osc.OscillateLSND(18,12,E)) + FD_TrueVsRecoVsPID_NC_NuMuFrac->GetBinContent(ip+1,ir+1,it+1)*(1-osc.OscillateLSND(18,14,E));
	      osc.SetOscParam(OscPar::kNuAntiNu,-1);
	      temp += (FD_TrueVsRecoVsPID_NC_NueBarFrac->GetBinContent(ip+1,ir+1,it+1)*(1-osc.OscillateLSND(-18,-12,E)) + FD_TrueVsRecoVsPID_NC_NuMuBarFrac->GetBinContent(ip+1,ir+1,it+1)*(1-osc.OscillateLSND(-18,-14,E)));
	      temp *= Pred_3D[qNC]->GetBinContent(ip+1,ir+1,it+1);
	    }
            else{ temp = Pred_3D[qNC]->GetBinContent(ip+1,ir+1,it+1);}
	  }
	  sum+=temp;
        }
	Pred[bg]->SetBinContent(ip+1,ir+1,sum);
        if(bg==Background::kNueCC)
        {
          Pred_NueCCSignal->SetBinContent(ip+1,ir+1,sumnu);
          Pred_NueBarCCSignal->SetBinContent(ip+1,ir+1,sumnubar);
        }
      }
    }
    if(FDComponents[i]!=Background::kNueCC)
    {
      Pred_TotalBkgd->Add(Pred[bg]);
    }
  }
  
  Set1DPredHists();
  
  return;
}
void Extrapolate2D::OscillatePrediction_NoSep()
{
  //will use Pred_3D to reweight Pred for oscillation
  
  for(unsigned int i=0;i<FDComponents.size();i++)
  {
    Pred[FDComponents[i]]->Reset();
  }
  Pred_TotalBkgd->Reset();
  
  int it,ir,ip;
  double E;
  double temp,sum;
  
  osc.SetOscParam(OscPar::kNuAntiNu,1);//this function assumes neutrino probabilites
  
  for(unsigned int i=0;i<FDComponents.size();i++)
  {
    Background::Background_t bg = FDComponents[i];
    for(ip=0;ip<nPID;ip++)
    {
      for(ir=0;ir<nReco;ir++)
      {
	sum=0;
        for(it=0;it<nTrue;it++)
        {
	  E = Pred_3D[qNC]->GetZaxis()->GetBinCenter(it+1);
	  
	  temp=0;
	  
	  if(bg==Background::kNuMuCC)
	  {
	    if(NominalOscProb[qNuMuToNuMu]->GetBinContent(it+1)>1e-20)
	    { 
	      if(OscMethod == 3){temp += Pred_3D[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(14,14,E)/NominalOscProb[qNuMuToNuMu]->GetBinContent(it+1);}
	      else if(OscMethod == 2){temp += Pred_3D[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(14,14,E)/NominalOscProb[qNuMuToNuMu]->GetBinContent(it+1);}
	      else if(OscMethod == 1){temp += Pred_3D[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(14,14,E)/NominalOscProb[qNuMuToNuMu]->GetBinContent(it+1);}
	      else{temp += Pred_3D[qNuMuToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.MuToMu(E)/NominalOscProb[qNuMuToNuMu]->GetBinContent(it+1);}
	    }
	    if(NominalOscProb[qBNueToNuMu]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp += Pred_3D[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(14,12,E)/NominalOscProb[qBNueToNuMu]->GetBinContent(it+1);}
	      else if(OscMethod == 2){temp += Pred_3D[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(14,12,E)/NominalOscProb[qBNueToNuMu]->GetBinContent(it+1);}
	      else if(OscMethod == 1){temp += Pred_3D[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(14,12,E)/NominalOscProb[qBNueToNuMu]->GetBinContent(it+1);}
	      else{temp += Pred_3D[qBNueToNuMu]->GetBinContent(ip+1,ir+1,it+1)*osc.ElecToMu(E)/NominalOscProb[qBNueToNuMu]->GetBinContent(it+1);}
	    }
	  }
	  if(bg==Background::kBNueCC)
	  {
	    if(NominalOscProb[qBNueToBNue]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp = Pred_3D[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(12,12,E)/NominalOscProb[qBNueToBNue]->GetBinContent(it+1);}
	      else if(OscMethod == 2){temp = Pred_3D[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(12,12,E)/NominalOscProb[qBNueToBNue]->GetBinContent(it+1);}
	      else if(OscMethod == 1){temp = Pred_3D[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(12,12,E)/NominalOscProb[qBNueToBNue]->GetBinContent(it+1);}
	      else{temp = Pred_3D[qBNueToBNue]->GetBinContent(ip+1,ir+1,it+1)*osc.ElecToElec(E)/NominalOscProb[qBNueToBNue]->GetBinContent(it+1);}
	    }
	  }
	  if(bg==Background::kNueCC)
	  {
	    if(NominalOscProb[qNuMuToNue]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp = Pred_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(12,14,E)/NominalOscProb[qNuMuToNue]->GetBinContent(it+1);}
	      else if(OscMethod == 2){temp = Pred_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(12,14,E)/NominalOscProb[qNuMuToNue]->GetBinContent(it+1);}
	      else if(OscMethod == 1){temp = Pred_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(12,14,E)/NominalOscProb[qNuMuToNue]->GetBinContent(it+1);}
	      else{temp = Pred_3D[qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1)*osc.MuToElec(E)/NominalOscProb[qNuMuToNue]->GetBinContent(it+1);}
	    }
	  }
	  if(bg==Background::kNuTauCC)
	  {
	    if(NominalOscProb[qNuMuToNuTau]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp += Pred_3D[qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(16,14,E)/NominalOscProb[qNuMuToNuTau]->GetBinContent(it+1);}
	      else if(OscMethod == 2){temp += Pred_3D[qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(16,14,E)/NominalOscProb[qNuMuToNuTau]->GetBinContent(it+1);}
	      else if(OscMethod == 1){temp += Pred_3D[qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(16,14,E)/NominalOscProb[qNuMuToNuTau]->GetBinContent(it+1);}
	      else{temp += Pred_3D[qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.MuToTau(E)/NominalOscProb[qNuMuToNuTau]->GetBinContent(it+1);}
	    }
	    if(NominalOscProb[qBNueToNuTau]->GetBinContent(it+1)>1e-20)
	    {
	      if(OscMethod == 3){temp += Pred_3D[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateLSND(16,12,E)/NominalOscProb[qBNueToNuTau]->GetBinContent(it+1);}
	      else if(OscMethod == 2){temp += Pred_3D[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.OscillateNSI(16,12,E)/NominalOscProb[qBNueToNuTau]->GetBinContent(it+1);}
	      else if(OscMethod == 1){temp += Pred_3D[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.Oscillate(16,12,E)/NominalOscProb[qBNueToNuTau]->GetBinContent(it+1);}
	      else{temp += Pred_3D[qBNueToNuTau]->GetBinContent(ip+1,ir+1,it+1)*osc.ElecToTau(E)/NominalOscProb[qBNueToNuTau]->GetBinContent(it+1);}
	    }
	  }
	  if(bg==Background::kNC)
	  {
            if(OscMethod == 3)
	    {
	      osc.SetOscParam(OscPar::kNuAntiNu,1);
	      temp = FD_TrueVsRecoVsPID_NC_NueFrac->GetBinContent(ip+1,ir+1,it+1)*(1-osc.OscillateLSND(18,12,E)) + FD_TrueVsRecoVsPID_NC_NuMuFrac->GetBinContent(ip+1,ir+1,it+1)*(1-osc.OscillateLSND(18,14,E));
	      osc.SetOscParam(OscPar::kNuAntiNu,-1);
	      temp += (FD_TrueVsRecoVsPID_NC_NueBarFrac->GetBinContent(ip+1,ir+1,it+1)*(1-osc.OscillateLSND(-18,-12,E)) + FD_TrueVsRecoVsPID_NC_NuMuBarFrac->GetBinContent(ip+1,ir+1,it+1)*(1-osc.OscillateLSND(-18,-14,E)));
	      temp *= Pred_3D[qNC]->GetBinContent(ip+1,ir+1,it+1);
	    }
	    else {temp = Pred_3D[qNC]->GetBinContent(ip+1,ir+1,it+1);}
	  }
	  sum+=temp;
        }
	Pred[bg]->SetBinContent(ip+1,ir+1,sum);
      }
    }
    if(FDComponents[i]!=Background::kNueCC)
    {
      Pred_TotalBkgd->Add(Pred[bg]);
    }
  }
  
  Set1DPredHists();
  
  return;
}
void Extrapolate2D::SetupPredHists()
{
  Pred_3D.clear();
  Pred_Nu_3D.clear();
  Pred_NuBar_3D.clear();
  Pred.clear();
  Pred_TotalBkgd = NULL;
  Pred_NueCCSignal = NULL;
  Pred_NueBarCCSignal = NULL;
  Pred_CC_Fid = NULL;
  Pred_NuMuCC_Fid = NULL;
  Pred_NuMuBarCC_Fid = NULL;
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
  
  TH3D *h3;
  string name;
  
  h3 = new TH3D("Pred_3D_NC","",nPID,p,nReco,r,nTrue,t);
  Pred_3D[qNC] = h3;
  h3 = new TH3D("Pred_3D_NuMuToNuMu","",nPID,p,nReco,r,nTrue,t);
  Pred_3D[qNuMuToNuMu] = h3;
  h3 = new TH3D("Pred_3D_BNueToNuMu","",nPID,p,nReco,r,nTrue,t);
  Pred_3D[qBNueToNuMu] = h3;
  h3 = new TH3D("Pred_3D_BNueToBNue","",nPID,p,nReco,r,nTrue,t);
  Pred_3D[qBNueToBNue] = h3;
  h3 = new TH3D("Pred_3D_NuMuToNue","",nPID,p,nReco,r,nTrue,t);
  Pred_3D[qNuMuToNue] = h3;
  h3 = new TH3D("Pred_3D_NuMuToNuTau","",nPID,p,nReco,r,nTrue,t);
  Pred_3D[qNuMuToNuTau] = h3;
  h3 = new TH3D("Pred_3D_BNueToNuTau","",nPID,p,nReco,r,nTrue,t);
  Pred_3D[qBNueToNuTau] = h3;
  
  //nu only
  h3 = new TH3D("Pred_Nu_3D_NC","",nPID,p,nReco,r,nTrue,t);
  Pred_Nu_3D[qNC] = h3;
  h3 = new TH3D("Pred_Nu_3D_NuMuToNuMu","",nPID,p,nReco,r,nTrue,t);
  Pred_Nu_3D[qNuMuToNuMu] = h3;
  h3 = new TH3D("Pred_Nu_3D_BNueToNuMu","",nPID,p,nReco,r,nTrue,t);
  Pred_Nu_3D[qBNueToNuMu] = h3;
  h3 = new TH3D("Pred_Nu_3D_BNueToBNue","",nPID,p,nReco,r,nTrue,t);
  Pred_Nu_3D[qBNueToBNue] = h3;
  h3 = new TH3D("Pred_Nu_3D_NuMuToNue","",nPID,p,nReco,r,nTrue,t);
  Pred_Nu_3D[qNuMuToNue] = h3;
  h3 = new TH3D("Pred_Nu_3D_NuMuToNuTau","",nPID,p,nReco,r,nTrue,t);
  Pred_Nu_3D[qNuMuToNuTau] = h3;
  h3 = new TH3D("Pred_Nu_3D_BNueToNuTau","",nPID,p,nReco,r,nTrue,t);
  Pred_Nu_3D[qBNueToNuTau] = h3;
  
  //nubar only
  h3 = new TH3D("Pred_NuBar_3D_NC","",nPID,p,nReco,r,nTrue,t);
  Pred_NuBar_3D[qNC] = h3;
  h3 = new TH3D("Pred_NuBar_3D_NuMuToNuMu","",nPID,p,nReco,r,nTrue,t);
  Pred_NuBar_3D[qNuMuToNuMu] = h3;
  h3 = new TH3D("Pred_NuBar_3D_BNueToNuMu","",nPID,p,nReco,r,nTrue,t);
  Pred_NuBar_3D[qBNueToNuMu] = h3;
  h3 = new TH3D("Pred_NuBar_3D_BNueToBNue","",nPID,p,nReco,r,nTrue,t);
  Pred_NuBar_3D[qBNueToBNue] = h3;
  h3 = new TH3D("Pred_NuBar_3D_NuMuToNue","",nPID,p,nReco,r,nTrue,t);
  Pred_NuBar_3D[qNuMuToNue] = h3;
  h3 = new TH3D("Pred_NuBar_3D_NuMuToNuTau","",nPID,p,nReco,r,nTrue,t);
  Pred_NuBar_3D[qNuMuToNuTau] = h3;
  h3 = new TH3D("Pred_NuBar_3D_BNueToNuTau","",nPID,p,nReco,r,nTrue,t);
  Pred_NuBar_3D[qBNueToNuTau] = h3;
  
  TH2D *h2;
  for(unsigned int j=0;j<FDComponents.size();j++)
  {
    name = "Pred_" + string(Background::AsString(FDComponents[j]));
    h2 = new TH2D(name.c_str(),"",nPID,p,nReco,r);
    Pred[FDComponents[j]] = h2;
  }
  
  Pred_TotalBkgd = new TH2D("Pred_TotalBkgd","",nPID,p,nReco,r);
  Pred_NueCCSignal = new TH2D("Pred_NueCCSignal","",nPID,p,nReco,r);
  Pred_NueBarCCSignal = new TH2D("Pred_NueBarCCSignal","",nPID,p,nReco,r);
  
  Pred_CC_Fid = new TH1D("Pred_CC_Fid","",nTrue,t);
  Pred_NuMuCC_Fid = new TH1D("Pred_NuMuCC_Fid","",nTrue,t);
  Pred_NuMuBarCC_Fid = new TH1D("Pred_NuMuBarCC_Fid","",nTrue,t);
  
  int totbin = nReco*nPID;
  Pred_TotalBkgd_VsBinNumber = new TH1D("Pred_TotalBkgd_VsBinNumber","",totbin,-0.5,totbin-0.5);
  Pred_Signal_VsBinNumber = new TH1D("Pred_Signal_VsBinNumber","",totbin,-0.5,totbin-0.5);
  
  delete [] r;
  delete [] p;
  delete [] t;
  
  return;
}
void Extrapolate2D::GetNoOscCCPrediction()
{
  //get FD CC prediction for use in tau/nue prediction
  int ir,it;
  double sum;
  double pos,neg;
  
  if(UseSeparateNuNuBar)
  {
    for(it=0;it<nTrue;it++)
    {
      pos=0.0;
      neg=0.0;
      for(ir=0;ir<nRecoCC;ir++)
      {
        if(ND_Reco_CClike_Pos->GetBinContent(ir+1)>0)
        {
          pos += (NDData_Reco_CClike_Pos->GetBinContent(ir+1)*(FD_Reco_CClike_Pos->GetBinContent(ir+1)/ND_Reco_CClike_Pos->GetBinContent(ir+1))*FD_Reco2True_CClike_Pos->GetBinContent(ir+1,it+1));
        }
        if(ND_Reco_CClike_Neg->GetBinContent(ir+1)>0)
        {
          neg += (NDData_Reco_CClike_Neg->GetBinContent(ir+1)*(FD_Reco_CClike_Neg->GetBinContent(ir+1)/ND_Reco_CClike_Neg->GetBinContent(ir+1))*FD_Reco2True_CClike_Neg->GetBinContent(ir+1,it+1));
        }
      }
      
      if(FD_Eff_NuMuBarCC_Pos->GetBinContent(it+1)>0)
      {
        pos = pos*FD_Purity_NuMuBarCC_Pos->GetBinContent(it+1)/FD_Eff_NuMuBarCC_Pos->GetBinContent(it+1);
      }
      else
      {
        pos=0.0;
      }
      Pred_NuMuBarCC_Fid->SetBinContent(it+1,pos);
      
      if(FD_Eff_NuMuCC_Neg->GetBinContent(it+1)>0)
      {
        neg = neg*FD_Purity_NuMuCC_Neg->GetBinContent(it+1)/FD_Eff_NuMuCC_Neg->GetBinContent(it+1);
      }
      else
      {
        neg=0.0;
      }
      Pred_NuMuCC_Fid->SetBinContent(it+1,neg);
    }
    
    Pred_CC_Fid->Add(Pred_NuMuCC_Fid);
    Pred_CC_Fid->Add(Pred_NuMuBarCC_Fid);
  }
  else
  {
    for(it=0;it<nTrue;it++)
    {
      sum=0.0;
      for(ir=0;ir<nRecoCC;ir++)
      {
        if(ND_Reco_CClike->GetBinContent(ir+1)>0)
        {
          sum += (NDData_Reco_CClike->GetBinContent(ir+1)*(FD_Reco_CClike->GetBinContent(ir+1)/ND_Reco_CClike->GetBinContent(ir+1))*FD_Reco2True_CClike->GetBinContent(ir+1,it+1));
        }
      }
      if(FD_Eff_CC->GetBinContent(it+1)>0)
      {
        sum = sum*FD_Purity_CC->GetBinContent(it+1)/FD_Eff_CC->GetBinContent(it+1);
      }
      else
      {
        sum=0.0;
      }
      Pred_CC_Fid->SetBinContent(it+1,sum);
    }
  }
  
  return;
}
void Extrapolate2D::Set1DPredHists()
{
  int ip,ir,i;
  
  i=0;
  for(ir=0;ir<nReco;ir++)
  {
    for(ip=0;ip<nPID;ip++)
    {
      Pred_TotalBkgd_VsBinNumber->SetBinContent(i+1,Pred_TotalBkgd->GetBinContent(ip+1,ir+1));
      Pred_Signal_VsBinNumber->SetBinContent(i+1,Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1));
      i++;
    }
  }
  
  return;
}
void Extrapolate2D::RebinInputHists()
{
  int nr = nReco;
  int np = nPID;
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
  
  MBH.Rebin2DHist(FD_TrueVsReco_Fid_NueCC,nr,r,0,0);
  MBH.Rebin2DHist(FD_TrueVsReco_Fid_NuTauCC,nr,r,0,0);
  if(UseSeparateNuNuBar)
  {
    MBH.Rebin2DHist(FD_TrueVsReco_Fid_NueBarCC,nr,r,0,0);
    MBH.Rebin2DHist(FD_TrueVsReco_Fid_NuTauBarCC,nr,r,0,0);
  }

  MBH.Rebin3DHist(ND_TrueVsRecoVsPID[qNC],np,p,nr,r,0,0);
  MBH.Rebin3DHist(ND_TrueVsRecoVsPID[qNuMuToNuMu],np,p,nr,r,0,0);
  MBH.Rebin3DHist(ND_TrueVsRecoVsPID[qBNueToBNue],np,p,nr,r,0,0);
  
  MBH.Rebin3DHist(FD_TrueVsRecoVsPID[qNC],np,p,nr,r,0,0);
  MBH.Rebin3DHist(FD_TrueVsRecoVsPID[qNuMuToNuMu],np,p,nr,r,0,0);
  MBH.Rebin3DHist(FD_TrueVsRecoVsPID[qBNueToNuMu],np,p,nr,r,0,0);
  MBH.Rebin3DHist(FD_TrueVsRecoVsPID[qBNueToBNue],np,p,nr,r,0,0);
  MBH.Rebin3DHist(FD_TrueVsRecoVsPID[qNuMuToNue],np,p,nr,r,0,0);
  MBH.Rebin3DHist(FD_TrueVsRecoVsPID[qNuMuToNuTau],np,p,nr,r,0,0);
  MBH.Rebin3DHist(FD_TrueVsRecoVsPID[qBNueToNuTau],np,p,nr,r,0,0);
  
  if( OscMethod == 3){
  MBH.Rebin3DHist(FD_TrueVsRecoVsPID_NC_NueFrac,np,p,nr,r,0,0);
  MBH.Rebin3DHist(FD_TrueVsRecoVsPID_NC_NueBarFrac,np,p,nr,r,0,0);
  MBH.Rebin3DHist(FD_TrueVsRecoVsPID_NC_NuMuFrac,np,p,nr,r,0,0);
  MBH.Rebin3DHist(FD_TrueVsRecoVsPID_NC_NuMuBarFrac,np,p,nr,r,0,0);
  }
  if(UseSeparateNuNuBar)
  {
    MBH.Rebin3DHist(FD_Nu_TrueVsRecoVsPID[qNC],np,p,nr,r,0,0);
    MBH.Rebin3DHist(FD_Nu_TrueVsRecoVsPID[qNuMuToNuMu],np,p,nr,r,0,0);
    MBH.Rebin3DHist(FD_Nu_TrueVsRecoVsPID[qBNueToNuMu],np,p,nr,r,0,0);
    MBH.Rebin3DHist(FD_Nu_TrueVsRecoVsPID[qBNueToBNue],np,p,nr,r,0,0);
    MBH.Rebin3DHist(FD_Nu_TrueVsRecoVsPID[qNuMuToNue],np,p,nr,r,0,0);
    MBH.Rebin3DHist(FD_Nu_TrueVsRecoVsPID[qNuMuToNuTau],np,p,nr,r,0,0);
    MBH.Rebin3DHist(FD_Nu_TrueVsRecoVsPID[qBNueToNuTau],np,p,nr,r,0,0);
    
    MBH.Rebin3DHist(FD_NuBar_TrueVsRecoVsPID[qNC],np,p,nr,r,0,0);
    MBH.Rebin3DHist(FD_NuBar_TrueVsRecoVsPID[qNuMuToNuMu],np,p,nr,r,0,0);
    MBH.Rebin3DHist(FD_NuBar_TrueVsRecoVsPID[qBNueToNuMu],np,p,nr,r,0,0);
    MBH.Rebin3DHist(FD_NuBar_TrueVsRecoVsPID[qBNueToBNue],np,p,nr,r,0,0);
    MBH.Rebin3DHist(FD_NuBar_TrueVsRecoVsPID[qNuMuToNue],np,p,nr,r,0,0);
    MBH.Rebin3DHist(FD_NuBar_TrueVsRecoVsPID[qNuMuToNuTau],np,p,nr,r,0,0);
    MBH.Rebin3DHist(FD_NuBar_TrueVsRecoVsPID[qBNueToNuTau],np,p,nr,r,0,0);
  }
  
  delete [] r;
  delete [] p;
  
  return;
}
