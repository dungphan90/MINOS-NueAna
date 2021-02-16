#include "NueAna/MultiBinAna/ExtrapFileGen.h"

ExtrapFileGen::ExtrapFileGen()
{
  Initialize();
  return;
}
void ExtrapFileGen::Initialize()
{
  SetOutputFile();
  SetNearPOT();
  SetFarPOT();
  SetFileLists();
  
  fout = 0;
  paramtree = 0;
  
  nReco = 0;
  nPID = 0;
  nTrue = 0;
  
  fakedata=false;
  
  return;
}
void ExtrapFileGen::AddPID(string pidname)
{
  PIDlist.push_back(pidname);
  return;
}
void ExtrapFileGen::SetFileLists(string nddata,string ndmc,string fdmcbeam,string fdmcnue,string fdmctau)
{
  FileList_NDData = nddata;
  FileList_NDMC = ndmc;
  FileList_FDMC_Beam = fdmcbeam;
  FileList_FDMC_Nue = fdmcnue;
  FileList_FDMC_Tau = fdmctau;
  return;
}
void ExtrapFileGen::SetOutputFile(std::string name)
{
   outFileName = name;
   
   return;
}
void ExtrapFileGen::WriteToFile()
{
//   double par[100]={0};
//   NueStandard::GetOscParam(par);
//   Theta12 = par[OscPar::kTh12];
//   Theta13 = par[OscPar::kTh13];
//   Theta23 = par[OscPar::kTh23];
//   DeltaMSq23 = par[OscPar::kDeltaM23];
//   DeltaMSq12 = par[OscPar::kDeltaM12];
//   DeltaCP = par[OscPar::kDelta];
  
  fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  paramtree = new TTree("paramtree","paramtree");
  paramtree->Branch("nearPOT",&fNearPOT,"nearPOT/D");
  paramtree->Branch("farPOT",&fFarPOT,"farPOT/D");
  paramtree->Branch("Theta12",&Theta12,"Theta12/D");
  paramtree->Branch("Theta13",&Theta13,"Theta13/D");
  paramtree->Branch("Theta23",&Theta23,"Theta23/D");
  paramtree->Branch("DeltaMSq23",&DeltaMSq23,"DeltaMSq23/D");
  paramtree->Branch("DeltaMSq12",&DeltaMSq12,"DeltaMSq12/D");
  paramtree->Branch("DeltaCP",&DeltaCP,"DeltaCP/D");
  paramtree->Fill();
  
  fout->cd();
  TDirectory *ndmcdir = fout->mkdir("NDMC");
  ndmcdir->cd();
  TList *list = ExHists_NDMC->fDirectory->GetList();
  TIter iter(list->MakeIterator());
  TObject *ob = 0;
  while((ob = iter())) 
  {
    ob->Write();
  }
  
  fout->cd();
  TDirectory *nddatadir = fout->mkdir("NDData");
  nddatadir->cd();
  list = ExHists_NDData->fDirectory->GetList();
  iter = list->MakeIterator();
  while((ob = iter())) 
  {
    ob->Write();
  }
  
  fout->cd();
  TDirectory *fdmcdir = fout->mkdir("FDMC");
  fdmcdir->cd();
  list = ExHists_FDMC->fDirectory->GetList();
  iter = list->MakeIterator();
  while((ob = iter())) 
  {
    ob->Write();
  }
  
  fout->cd();
  TDirectory *fdmcdirnu = fout->mkdir("FDMC_Nu");
  fdmcdirnu->cd();
  list = ExHists_Nu_FDMC->fDirectory->GetList();
  iter = list->MakeIterator();
  while((ob = iter())) 
  {
    ob->Write();
  }
  
  fout->cd();
  TDirectory *fdmcdirnubar = fout->mkdir("FDMC_NuBar");
  fdmcdirnubar->cd();
  list = ExHists_NuBar_FDMC->fDirectory->GetList();
  iter = list->MakeIterator();
  while((ob = iter())) 
  {
    ob->Write();
  }
  
  fout->cd();
  paramtree->Write();
  fout->Close();
}

void ExtrapFileGen::InitializeHistograms()
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
  
  gDirectory->cd("/");
  unsigned int np = PIDlist.size();
  unsigned int j;
  
  ExHists_NDMC = new ExtrapHists2D("ExtrapHists_NDMC");
  ExHists_NDMC->fDirectory->cd();
  for(j=0;j<np;j++)
  {
    ExHists_NDMC->AddPID(PIDlist[j]);
  }
  ExHists_NDMC->SetPIDBins(nPID,p);
  ExHists_NDMC->SetRecoBins(nReco,r);
  ExHists_NDMC->SetTrueBins(nTrue,t);
  ExHists_NDMC->InitializeHists();
  
  ExHists_FDMC = new ExtrapHists2D("ExtrapHists_FDMC");
  ExHists_FDMC->fDirectory->cd();
  for(j=0;j<np;j++)
  {
    ExHists_FDMC->AddPID(PIDlist[j]);
  }
  ExHists_FDMC->SetPIDBins(nPID,p);
  ExHists_FDMC->SetRecoBins(nReco,r);
  ExHists_FDMC->SetTrueBins(nTrue,t);
  ExHists_FDMC->InitializeHists();
  
  ExHists_Nu_FDMC = new ExtrapHists2D("ExtrapHists_Nu_FDMC");
  ExHists_Nu_FDMC->fDirectory->cd();
  for(j=0;j<np;j++)
  {
    ExHists_Nu_FDMC->AddPID(PIDlist[j]);
  }
  ExHists_Nu_FDMC->SetPIDBins(nPID,p);
  ExHists_Nu_FDMC->SetRecoBins(nReco,r);
  ExHists_Nu_FDMC->SetTrueBins(nTrue,t);
  ExHists_Nu_FDMC->InitializeHists();
  
  ExHists_NuBar_FDMC = new ExtrapHists2D("ExtrapHists_NuBar_FDMC");
  ExHists_NuBar_FDMC->fDirectory->cd();
  for(j=0;j<np;j++)
  {
    ExHists_NuBar_FDMC->AddPID(PIDlist[j]);
  }
  ExHists_NuBar_FDMC->SetPIDBins(nPID,p);
  ExHists_NuBar_FDMC->SetRecoBins(nReco,r);
  ExHists_NuBar_FDMC->SetTrueBins(nTrue,t);
  ExHists_NuBar_FDMC->InitializeHists();
  
  ExHists_NDData = new ExtrapHists2D("ExtrapHist_NDData");
  ExHists_NDData->fDirectory->cd();
  ExHists_NDData->SetPIDBins(nPID,p);
  ExHists_NDData->SetRecoBins(nReco,r);
  ExHists_NDData->SetTrueBins(nTrue,t);
  ExHists_NDData->InitializeHists();
  
  return;
}
void ExtrapFileGen::Run()
{
  //first read in one file to get the bin edges
  ifstream flist;
  string line;
  flist.open(gSystem->ExpandPathName(FileList_NDMC.c_str()));
  getline(flist,line);
  TH3D *h3;
  int i;
  if(!line.empty())
  {
    if(gSystem->AccessPathName(gSystem->ExpandPathName(line.c_str())))
      {
	cout<<line<<" doesn't exist.  Quitting..."<<endl;
	return;
      }
      TFile f(line.c_str(),"READ");
      if(f.Read("NDMC")==0)
      {
	cout<<"NDMC directory doesn't exist in "<<line<<".  Quitting..."<<endl;
	return;
      }
      h3 = (TH3D*)f.Get(Form("NDMC/TrueVsRecoVs%s_NC",PIDlist[0].c_str()));
      nPID = h3->GetNbinsX();
      nReco = h3->GetNbinsY();
      nTrue = h3->GetNbinsZ();
      for(i=0;i<nPID+1;i++)
      {
	PIDEdges.push_back(h3->GetXaxis()->GetBinLowEdge(i+1));
      }
      for(i=0;i<nReco+1;i++)
      {
	RecoEdges.push_back(h3->GetYaxis()->GetBinLowEdge(i+1));
      }
      for(i=0;i<nTrue+1;i++)
      {
	TrueEdges.push_back(h3->GetZaxis()->GetBinLowEdge(i+1));
      }
  }
  
  InitializeHistograms();
  
  unsigned int np = PIDlist.size();
  unsigned int j;
  int n=0;
  
  gROOT->cd("/");
  TTree *paramtree_;
  double del,dm2_23,dm2_12,t12,t13,t23,npot;
  
  if(!fakedata)//if using real ND data
  {
  cout<<"NDData:"<<endl;
  double tot_ndata=0;
  ifstream filelist;
  filelist.open(gSystem->ExpandPathName(FileList_NDData.c_str()));
  if(!filelist.good())
  {
    cout<<"Problem reading "<<FileList_NDData<<". Quitting..."<<endl;
    return;
  }
  n=0;
  while(!filelist.eof())
  {
    cout<<n<<endl;
    getline(filelist,line);
    if(!line.empty())
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(line.c_str())))
      {
	cout<<line<<" doesn't exist.  Skipping..."<<endl;
	continue;
      }
      TFile f(line.c_str(),"READ");
      if(f.Read("NDData")==0 || f.Read("paramtree")==0)
      {
	cout<<"NDData directory or paramtree doesn't exist.  Skipping..."<<endl;
	continue;
      }
      
      ExHists_NDData->Reco_CClike->Add((TH1D*)f.Get("NDData/Reco_CClike"));
      ExHists_NDData->Reco_CClike_Pos->Add((TH1D*)f.Get("NDData/Reco_CClike_Pos"));
      ExHists_NDData->Reco_CClike_Neg->Add((TH1D*)f.Get("NDData/Reco_CClike_Neg"));
      
      paramtree_ = (TTree*)f.Get("paramtree");
      paramtree_->SetBranchAddress("POT",&npot);
      paramtree_->GetEntry(0);
      tot_ndata+=npot;
    }
    n++;
  }
  filelist.close();
  
  ExHists_NDData->Reco_CClike->Scale(fNearPOT/tot_ndata);
  ExHists_NDData->Reco_CClike_Pos->Scale(fNearPOT/tot_ndata);
  ExHists_NDData->Reco_CClike_Neg->Scale(fNearPOT/tot_ndata);
  }
  
  cout<<"NDMC:"<<endl;
  double tot_nmc=0;
  ifstream filelist2;
  filelist2.open(gSystem->ExpandPathName(FileList_NDMC.c_str()));
  if(!filelist2.good())
  {
    cout<<"Problem reading "<<FileList_NDMC<<". Quitting..."<<endl;
    return;
  }
  n=0;
  while(!filelist2.eof())
  {
    cout<<n<<endl;
    getline(filelist2,line);
    if(!line.empty())
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(line.c_str())))
      {
	cout<<line<<" doesn't exist.  Skipping..."<<endl;
	continue;
      }
      TFile f(line.c_str(),"READ");
      if(f.Read("NDMC")==0 || f.Read("paramtree")==0)
      {
	cout<<"NDMC directory or paramtree doesn't exist.  Skipping..."<<endl;
	continue;
      }
      
      ExHists_NDMC->TrueVsReco_Fid->Add((TH2D*)f.Get("NDMC/TrueVsReco_Fid"));
      ExHists_NDMC->TrueVsReco_Fid_NuMuCC->Add((TH2D*)f.Get("NDMC/TrueVsReco_Fid_NuMuCC"));
      ExHists_NDMC->TrueVsReco_Fid_NuMuBarCC->Add((TH2D*)f.Get("NDMC/TrueVsReco_Fid_NuMuBarCC"));
      ExHists_NDMC->TrueVsReco_CClike->Add((TH2D*)f.Get("NDMC/TrueVsReco_CClike"));
      ExHists_NDMC->TrueVsReco_CClike_Pos->Add((TH2D*)f.Get("NDMC/TrueVsReco_CClike_Pos"));
      ExHists_NDMC->TrueVsReco_CClike_Neg->Add((TH2D*)f.Get("NDMC/TrueVsReco_CClike_Neg"));
      ExHists_NDMC->TrueVsReco_CClike_NuMuCC->Add((TH2D*)f.Get("NDMC/TrueVsReco_CClike_NuMuCC"));
      ExHists_NDMC->TrueVsReco_CClike_Pos_NuMuBarCC->Add((TH2D*)f.Get("NDMC/TrueVsReco_CClike_Pos_NuMuBarCC"));
      ExHists_NDMC->TrueVsReco_CClike_Neg_NuMuCC->Add((TH2D*)f.Get("NDMC/TrueVsReco_CClike_Neg_NuMuCC"));
      
      for(j=0;j<np;j++)
      {
        ExHists_NDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Add((TH3D*)f.Get(Form("NDMC/TrueVsRecoVs%s_NC",PIDlist[j].c_str())));
        ExHists_NDMC->TrueVsRecoVsPID_NuMuToNuMu[PIDlist[j]]->Add((TH3D*)f.Get(Form("NDMC/TrueVsRecoVs%s_NuMuToNuMu",PIDlist[j].c_str())));
        ExHists_NDMC->TrueVsRecoVsPID_BNueToBNue[PIDlist[j]]->Add((TH3D*)f.Get(Form("NDMC/TrueVsRecoVs%s_BNueToBNue",PIDlist[j].c_str())));
      }
      
      paramtree_ = (TTree*)f.Get("paramtree");
      paramtree_->SetBranchAddress("POT",&npot);
      paramtree_->GetEntry(0);
      tot_nmc+=npot;
    }
    n++;
  }
  filelist2.close();
  
  ExHists_NDMC->TrueVsReco_Fid->Scale(fNearPOT/tot_nmc);
  ExHists_NDMC->TrueVsReco_Fid_NuMuCC->Scale(fNearPOT/tot_nmc);
  ExHists_NDMC->TrueVsReco_Fid_NuMuBarCC->Scale(fNearPOT/tot_nmc);
  ExHists_NDMC->TrueVsReco_CClike->Scale(fNearPOT/tot_nmc);
  ExHists_NDMC->TrueVsReco_CClike_Pos->Scale(fNearPOT/tot_nmc);
  ExHists_NDMC->TrueVsReco_CClike_Neg->Scale(fNearPOT/tot_nmc);
  ExHists_NDMC->TrueVsReco_CClike_NuMuCC->Scale(fNearPOT/tot_nmc);
  ExHists_NDMC->TrueVsReco_CClike_Pos_NuMuBarCC->Scale(fNearPOT/tot_nmc);
  ExHists_NDMC->TrueVsReco_CClike_Neg_NuMuCC->Scale(fNearPOT/tot_nmc);
  
  for(j=0;j<np;j++)
  {
    ExHists_NDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Scale(fNearPOT/tot_nmc);
    ExHists_NDMC->TrueVsRecoVsPID_NuMuToNuMu[PIDlist[j]]->Scale(fNearPOT/tot_nmc);
    ExHists_NDMC->TrueVsRecoVsPID_BNueToBNue[PIDlist[j]]->Scale(fNearPOT/tot_nmc);
  }
  
  if(fakedata)
  {
    ExHists_NDData->Reco_CClike->Add(ExHists_NDMC->TrueVsReco_CClike->ProjectionX());
    ExHists_NDData->Reco_CClike_Pos->Add(ExHists_NDMC->TrueVsReco_CClike_Pos->ProjectionX());
    ExHists_NDData->Reco_CClike_Neg->Add(ExHists_NDMC->TrueVsReco_CClike_Neg->ProjectionX());
  }
  
  gROOT->cd("/");
  vector<TH3D*> TrueVsRecoVsPID_NC_Beam;
  vector<TH3D*> TrueVsRecoVsPID_NC_Nue;
  vector<TH3D*> TrueVsRecoVsPID_NC_Tau;
  for(j=0;j<np;j++)
  {
    TrueVsRecoVsPID_NC_Beam.push_back((TH3D*)ExHists_NDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Clone(Form("NC_Beam_%s",PIDlist[j].c_str())));
    TrueVsRecoVsPID_NC_Nue.push_back((TH3D*)ExHists_NDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Clone(Form("NC_Nue_%s",PIDlist[j].c_str())));
    TrueVsRecoVsPID_NC_Tau.push_back((TH3D*)ExHists_NDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Clone(Form("NC_Tau_%s",PIDlist[j].c_str())));
    
    TrueVsRecoVsPID_NC_Beam[j]->Reset();
    TrueVsRecoVsPID_NC_Nue[j]->Reset();
    TrueVsRecoVsPID_NC_Tau[j]->Reset();
  }
  vector<TH3D*> TrueVsRecoVsPID_Nu_NC_Beam;
  vector<TH3D*> TrueVsRecoVsPID_Nu_NC_Nue;
  vector<TH3D*> TrueVsRecoVsPID_Nu_NC_Tau;
  for(j=0;j<np;j++)
  {
    TrueVsRecoVsPID_Nu_NC_Beam.push_back((TH3D*)ExHists_NDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Clone(Form("Nu_NC_Beam_%s",PIDlist[j].c_str())));
    TrueVsRecoVsPID_Nu_NC_Nue.push_back((TH3D*)ExHists_NDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Clone(Form("Nu_NC_Nue_%s",PIDlist[j].c_str())));
    TrueVsRecoVsPID_Nu_NC_Tau.push_back((TH3D*)ExHists_NDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Clone(Form("Nu_NC_Tau_%s",PIDlist[j].c_str())));
    
    TrueVsRecoVsPID_Nu_NC_Beam[j]->Reset();
    TrueVsRecoVsPID_Nu_NC_Nue[j]->Reset();
    TrueVsRecoVsPID_Nu_NC_Tau[j]->Reset();
  }
  vector<TH3D*> TrueVsRecoVsPID_NuBar_NC_Beam;
  vector<TH3D*> TrueVsRecoVsPID_NuBar_NC_Nue;
  vector<TH3D*> TrueVsRecoVsPID_NuBar_NC_Tau;
  for(j=0;j<np;j++)
  {
    TrueVsRecoVsPID_NuBar_NC_Beam.push_back((TH3D*)ExHists_NDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Clone(Form("NuBar_NC_Beam_%s",PIDlist[j].c_str())));
    TrueVsRecoVsPID_NuBar_NC_Nue.push_back((TH3D*)ExHists_NDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Clone(Form("NuBar_NC_Nue_%s",PIDlist[j].c_str())));
    TrueVsRecoVsPID_NuBar_NC_Tau.push_back((TH3D*)ExHists_NDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Clone(Form("NuBar_NC_Tau_%s",PIDlist[j].c_str())));
    
    TrueVsRecoVsPID_NuBar_NC_Beam[j]->Reset();
    TrueVsRecoVsPID_NuBar_NC_Nue[j]->Reset();
    TrueVsRecoVsPID_NuBar_NC_Tau[j]->Reset();
  }
  
  cout<<"FDMC (Beam):"<<endl;
  double tot_fbmc=0;
  ifstream filelist3;
  filelist3.open(gSystem->ExpandPathName(FileList_FDMC_Beam.c_str()));
  if(!filelist3.good())
  {
    cout<<"Problem reading "<<FileList_FDMC_Beam<<". Quitting..."<<endl;
    return;
  }
  n=0;
  while(!filelist3.eof())
  {
    cout<<n<<endl;
    getline(filelist3,line);
    if(!line.empty())
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(line.c_str())))
      {
	cout<<line<<" doesn't exist.  Skipping..."<<endl;
	continue;
      }
      TFile f(line.c_str(),"READ");
      if(f.Read("FDMC_Beam")==0 || f.Read("paramtree")==0)
      {
	cout<<"FDMC_Beam directory or paramtree doesn't exist.  Skipping..."<<endl;
	continue;
      }
      
      ExHists_FDMC->TrueVsReco_Fid->Add((TH2D*)f.Get("FDMC_Beam/TrueVsReco_Fid"));
      ExHists_FDMC->TrueVsReco_Fid_NuMuCC->Add((TH2D*)f.Get("FDMC_Beam/TrueVsReco_Fid_NuMuCC"));
      ExHists_FDMC->TrueVsReco_Fid_NuMuBarCC->Add((TH2D*)f.Get("FDMC_Beam/TrueVsReco_Fid_NuMuBarCC"));
      ExHists_FDMC->TrueVsReco_CClike->Add((TH2D*)f.Get("FDMC_Beam/TrueVsReco_CClike"));
      ExHists_FDMC->TrueVsReco_CClike_Pos->Add((TH2D*)f.Get("FDMC_Beam/TrueVsReco_CClike_Pos"));
      ExHists_FDMC->TrueVsReco_CClike_Neg->Add((TH2D*)f.Get("FDMC_Beam/TrueVsReco_CClike_Neg"));
      ExHists_FDMC->TrueVsReco_CClike_NuMuCC->Add((TH2D*)f.Get("FDMC_Beam/TrueVsReco_CClike_NuMuCC"));
      ExHists_FDMC->TrueVsReco_CClike_Pos_NuMuBarCC->Add((TH2D*)f.Get("FDMC_Beam/TrueVsReco_CClike_Pos_NuMuBarCC"));
      ExHists_FDMC->TrueVsReco_CClike_Neg_NuMuCC->Add((TH2D*)f.Get("FDMC_Beam/TrueVsReco_CClike_Neg_NuMuCC"));
      
      for(j=0;j<np;j++)
      {
        TrueVsRecoVsPID_NC_Beam[j]->Add((TH3D*)f.Get(Form("FDMC_Beam/TrueVsRecoVs%s_NC",PIDlist[j].c_str())));
        ExHists_FDMC->TrueVsRecoVsPID_NuMuToNuMu[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Beam/TrueVsRecoVs%s_NuMuToNuMu",PIDlist[j].c_str())));
        ExHists_FDMC->TrueVsRecoVsPID_BNueToBNue[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Beam/TrueVsRecoVs%s_BNueToBNue",PIDlist[j].c_str())));
	
	ExHists_FDMC->TrueVsRecoVsPID_NC_NoOscNue[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Beam/TrueVsRecoVs%s_NC_NoOscNue",PIDlist[j].c_str())));
	ExHists_FDMC->TrueVsRecoVsPID_NC_NoOscNueBar[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Beam/TrueVsRecoVs%s_NC_NoOscNueBar",PIDlist[j].c_str())));
	ExHists_FDMC->TrueVsRecoVsPID_NC_NoOscNuMu[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Beam/TrueVsRecoVs%s_NC_NoOscNuMu",PIDlist[j].c_str())));
	ExHists_FDMC->TrueVsRecoVsPID_NC_NoOscNuMuBar[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Beam/TrueVsRecoVs%s_NC_NoOscNuMuBar",PIDlist[j].c_str())));
      }
      
      paramtree_ = (TTree*)f.Get("paramtree");
      paramtree_->SetBranchAddress("POT",&npot);
      if(n==0)
      {
        paramtree_->SetBranchAddress("DeltaCP",&del);
        paramtree_->SetBranchAddress("DeltaMSq12",&dm2_12);
        paramtree_->SetBranchAddress("DeltaMSq23",&dm2_23);
        paramtree_->SetBranchAddress("Theta12",&t12);
        paramtree_->SetBranchAddress("Theta23",&t23);
        paramtree_->SetBranchAddress("Theta13",&t13);
      }
      paramtree_->GetEntry(0);
      tot_fbmc+=npot;
    }
    n++;
  }
  filelist3.close();
  
  cout<<"FDMC Nu (Beam):"<<endl;
  filelist3.clear();
  filelist3.open(gSystem->ExpandPathName(FileList_FDMC_Beam.c_str()));
  n=0;
  while(!filelist3.eof())
  {
    cout<<n<<endl;
    getline(filelist3,line);
    if(!line.empty())
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(line.c_str())))
      {
	cout<<line<<" doesn't exist.  Skipping..."<<endl;
	continue;
      }
      TFile f(line.c_str(),"READ");
      if(f.Read("FDMC_Beam_Nu")==0 || f.Read("paramtree")==0)
      {
	cout<<"FDMC_Beam_Nu directory or paramtree doesn't exist.  Skipping..."<<endl;
	continue;
      }
      
      for(j=0;j<np;j++)
      {
	TrueVsRecoVsPID_Nu_NC_Beam[j]->Add((TH3D*)f.Get(Form("FDMC_Beam_Nu/TrueVsRecoVs%s_NC",PIDlist[j].c_str())));
        ExHists_Nu_FDMC->TrueVsRecoVsPID_NuMuToNuMu[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Beam_Nu/TrueVsRecoVs%s_NuMuToNuMu",PIDlist[j].c_str())));
        ExHists_Nu_FDMC->TrueVsRecoVsPID_BNueToBNue[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Beam_Nu/TrueVsRecoVs%s_BNueToBNue",PIDlist[j].c_str())));
      }
    }
    n++;
  }
  filelist3.close();
  
  cout<<"FDMC NuBar (Beam):"<<endl;
  filelist3.clear();
  filelist3.open(gSystem->ExpandPathName(FileList_FDMC_Beam.c_str()));
  n=0;
  while(!filelist3.eof())
  {
    cout<<n<<endl;
    getline(filelist3,line);
    if(!line.empty())
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(line.c_str())))
      {
	cout<<line<<" doesn't exist.  Skipping..."<<endl;
	continue;
      }
      TFile f(line.c_str(),"READ");
      if(f.Read("FDMC_Beam_NuBar")==0 || f.Read("paramtree")==0)
      {
	cout<<"FDMC_Beam_NuBar directory or paramtree doesn't exist.  Skipping..."<<endl;
	continue;
      }
      
      for(j=0;j<np;j++)
      {
	TrueVsRecoVsPID_NuBar_NC_Beam[j]->Add((TH3D*)f.Get(Form("FDMC_Beam_NuBar/TrueVsRecoVs%s_NC",PIDlist[j].c_str())));
        ExHists_NuBar_FDMC->TrueVsRecoVsPID_NuMuToNuMu[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Beam_NuBar/TrueVsRecoVs%s_NuMuToNuMu",PIDlist[j].c_str())));
        ExHists_NuBar_FDMC->TrueVsRecoVsPID_BNueToBNue[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Beam_NuBar/TrueVsRecoVs%s_BNueToBNue",PIDlist[j].c_str())));
      }
    }
    n++;
  }
  filelist3.close();
  
  ExHists_FDMC->TrueVsReco_Fid->Scale(fFarPOT/tot_fbmc);
  ExHists_FDMC->TrueVsReco_Fid_NuMuCC->Scale(fFarPOT/tot_fbmc);
  ExHists_FDMC->TrueVsReco_Fid_NuMuBarCC->Scale(fFarPOT/tot_fbmc);
  ExHists_FDMC->TrueVsReco_CClike->Scale(fFarPOT/tot_fbmc);
  ExHists_FDMC->TrueVsReco_CClike_Pos->Scale(fFarPOT/tot_fbmc);
  ExHists_FDMC->TrueVsReco_CClike_Neg->Scale(fFarPOT/tot_fbmc);
  ExHists_FDMC->TrueVsReco_CClike_NuMuCC->Scale(fFarPOT/tot_fbmc);
  ExHists_FDMC->TrueVsReco_CClike_Pos_NuMuBarCC->Scale(fFarPOT/tot_fbmc);
  ExHists_FDMC->TrueVsReco_CClike_Neg_NuMuCC->Scale(fFarPOT/tot_fbmc);
  for(j=0;j<np;j++)
  {
    TrueVsRecoVsPID_NC_Beam[j]->Scale(fFarPOT/tot_fbmc);
    ExHists_FDMC->TrueVsRecoVsPID_NuMuToNuMu[PIDlist[j]]->Scale(fFarPOT/tot_fbmc);
    ExHists_FDMC->TrueVsRecoVsPID_BNueToBNue[PIDlist[j]]->Scale(fFarPOT/tot_fbmc);
    
    ExHists_FDMC->TrueVsRecoVsPID_NC_NoOscNue[PIDlist[j]]->Scale(fFarPOT/tot_fbmc);
    ExHists_FDMC->TrueVsRecoVsPID_NC_NoOscNueBar[PIDlist[j]]->Scale(fFarPOT/tot_fbmc);
    ExHists_FDMC->TrueVsRecoVsPID_NC_NoOscNuMu[PIDlist[j]]->Scale(fFarPOT/tot_fbmc);
    ExHists_FDMC->TrueVsRecoVsPID_NC_NoOscNuMuBar[PIDlist[j]]->Scale(fFarPOT/tot_fbmc);
    
    TrueVsRecoVsPID_Nu_NC_Beam[j]->Scale(fFarPOT/tot_fbmc);
    ExHists_Nu_FDMC->TrueVsRecoVsPID_NuMuToNuMu[PIDlist[j]]->Scale(fFarPOT/tot_fbmc);
    ExHists_Nu_FDMC->TrueVsRecoVsPID_BNueToBNue[PIDlist[j]]->Scale(fFarPOT/tot_fbmc);
    
    TrueVsRecoVsPID_NuBar_NC_Beam[j]->Scale(fFarPOT/tot_fbmc);
    ExHists_NuBar_FDMC->TrueVsRecoVsPID_NuMuToNuMu[PIDlist[j]]->Scale(fFarPOT/tot_fbmc);
    ExHists_NuBar_FDMC->TrueVsRecoVsPID_BNueToBNue[PIDlist[j]]->Scale(fFarPOT/tot_fbmc);
  }
  
  cout<<"FDMC (Nue):"<<endl;
  double tot_fnmc=0;
  ifstream filelist4;
  filelist4.open(gSystem->ExpandPathName(FileList_FDMC_Nue.c_str()));
  if(!filelist4.good())
  {
    cout<<"Problem reading "<<FileList_FDMC_Nue<<". Quitting..."<<endl;
    return;
  }
  n=0;
  while(!filelist4.eof())
  {
    cout<<n<<endl;
    getline(filelist4,line);
    if(!line.empty())
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(line.c_str())))
      {
	cout<<line<<" doesn't exist.  Skipping..."<<endl;
	continue;
      }
      TFile f(line.c_str(),"READ");
      if(f.Read("FDMC_Nue")==0 || f.Read("paramtree")==0)
      {
	cout<<"FDMC_Nue directory or paramtree doesn't exist.  Skipping..."<<endl;
	continue;
      }
      
      ExHists_FDMC->TrueVsReco_Fid_NueCC->Add((TH2D*)f.Get("FDMC_Nue/TrueVsReco_Fid_NueCC"));
      ExHists_FDMC->True_Fid_NueCC->Add((TH1D*)f.Get("FDMC_Nue/True_Fid_NueCC"));
      ExHists_FDMC->TrueVsReco_Fid_NueBarCC->Add((TH2D*)f.Get("FDMC_Nue/TrueVsReco_Fid_NueBarCC"));
      ExHists_FDMC->True_Fid_NueBarCC->Add((TH1D*)f.Get("FDMC_Nue/True_Fid_NueBarCC"));
      
      for(j=0;j<np;j++)
      {
        TrueVsRecoVsPID_NC_Nue[j]->Add((TH3D*)f.Get(Form("FDMC_Nue/TrueVsRecoVs%s_NC",PIDlist[j].c_str())));
        ExHists_FDMC->TrueVsRecoVsPID_NuMuToNue[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Nue/TrueVsRecoVs%s_NuMuToNue",PIDlist[j].c_str())));
        ExHists_FDMC->TrueVsRecoVsPID_BNueToNuMu[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Nue/TrueVsRecoVs%s_BNueToNuMu",PIDlist[j].c_str())));
      }
      
      paramtree_ = (TTree*)f.Get("paramtree");
      paramtree_->SetBranchAddress("POT",&npot);
      paramtree_->GetEntry(0);
      tot_fnmc+=npot;
    }
    n++;
  }
  filelist4.close();
  
  cout<<"FDMC Nu (Nue):"<<endl;
  filelist4.clear();
  filelist4.open(gSystem->ExpandPathName(FileList_FDMC_Nue.c_str()));
  n=0;
  while(!filelist4.eof())
  {
    cout<<n<<endl;
    getline(filelist4,line);
    if(!line.empty())
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(line.c_str())))
      {
	cout<<line<<" doesn't exist.  Skipping..."<<endl;
	continue;
      }
      TFile f(line.c_str(),"READ");
      if(f.Read("FDMC_Nue_Nu")==0 || f.Read("paramtree")==0)
      {
	cout<<"FDMC_Nue_Nu directory or paramtree doesn't exist.  Skipping..."<<endl;
	continue;
      }
      
      for(j=0;j<np;j++)
      {
	TrueVsRecoVsPID_Nu_NC_Nue[j]->Add((TH3D*)f.Get(Form("FDMC_Nue_Nu/TrueVsRecoVs%s_NC",PIDlist[j].c_str())));
        ExHists_Nu_FDMC->TrueVsRecoVsPID_NuMuToNue[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Nue_Nu/TrueVsRecoVs%s_NuMuToNue",PIDlist[j].c_str())));
        ExHists_Nu_FDMC->TrueVsRecoVsPID_BNueToNuMu[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Nue_Nu/TrueVsRecoVs%s_BNueToNuMu",PIDlist[j].c_str())));
      }
    }
    n++;
  }
  filelist4.close();
  
  cout<<"FDMC NuBar (Nue):"<<endl;
  filelist4.clear();
  filelist4.open(gSystem->ExpandPathName(FileList_FDMC_Nue.c_str()));
  n=0;
  while(!filelist4.eof())
  {
    cout<<n<<endl;
    getline(filelist4,line);
    if(!line.empty())
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(line.c_str())))
      {
	cout<<line<<" doesn't exist.  Skipping..."<<endl;
	continue;
      }
      TFile f(line.c_str(),"READ");
      if(f.Read("FDMC_Nue_NuBar")==0 || f.Read("paramtree")==0)
      {
	cout<<"FDMC_Nue_NuBar directory or paramtree doesn't exist.  Skipping..."<<endl;
	continue;
      }
      
      for(j=0;j<np;j++)
      {
	TrueVsRecoVsPID_NuBar_NC_Nue[j]->Add((TH3D*)f.Get(Form("FDMC_Nue_NuBar/TrueVsRecoVs%s_NC",PIDlist[j].c_str())));
        ExHists_NuBar_FDMC->TrueVsRecoVsPID_NuMuToNue[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Nue_NuBar/TrueVsRecoVs%s_NuMuToNue",PIDlist[j].c_str())));
        ExHists_NuBar_FDMC->TrueVsRecoVsPID_BNueToNuMu[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Nue_NuBar/TrueVsRecoVs%s_BNueToNuMu",PIDlist[j].c_str())));
      }
    }
    n++;
  }
  filelist4.close();
  
  ExHists_FDMC->TrueVsReco_Fid_NueCC->Scale(fFarPOT/tot_fnmc);
  ExHists_FDMC->True_Fid_NueCC->Scale(fFarPOT/tot_fnmc);
  ExHists_FDMC->TrueVsReco_Fid_NueBarCC->Scale(fFarPOT/tot_fnmc);
  ExHists_FDMC->True_Fid_NueBarCC->Scale(fFarPOT/tot_fnmc);
  for(j=0;j<np;j++)
  {
    TrueVsRecoVsPID_NC_Nue[j]->Scale(fFarPOT/tot_fnmc);
    ExHists_FDMC->TrueVsRecoVsPID_NuMuToNue[PIDlist[j]]->Scale(fFarPOT/tot_fnmc);
    ExHists_FDMC->TrueVsRecoVsPID_BNueToNuMu[PIDlist[j]]->Scale(fFarPOT/tot_fnmc);
    
    TrueVsRecoVsPID_Nu_NC_Nue[j]->Scale(fFarPOT/tot_fnmc);
    ExHists_Nu_FDMC->TrueVsRecoVsPID_NuMuToNue[PIDlist[j]]->Scale(fFarPOT/tot_fnmc);
    ExHists_Nu_FDMC->TrueVsRecoVsPID_BNueToNuMu[PIDlist[j]]->Scale(fFarPOT/tot_fnmc);
    
    TrueVsRecoVsPID_NuBar_NC_Nue[j]->Scale(fFarPOT/tot_fnmc);
    ExHists_NuBar_FDMC->TrueVsRecoVsPID_NuMuToNue[PIDlist[j]]->Scale(fFarPOT/tot_fnmc);
    ExHists_NuBar_FDMC->TrueVsRecoVsPID_BNueToNuMu[PIDlist[j]]->Scale(fFarPOT/tot_fnmc);
  }
  
  cout<<"FDMC (Tau):"<<endl;
  double tot_ftmc=0;
  ifstream filelist5;
  filelist5.open(gSystem->ExpandPathName(FileList_FDMC_Tau.c_str()));
  if(!filelist5.good())
  {
    cout<<"Problem reading "<<FileList_FDMC_Tau<<". Quitting..."<<endl;
    return;
  }
  n=0;
  while(!filelist5.eof())
  {
    cout<<n<<endl;
    getline(filelist5,line);
    if(!line.empty())
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(line.c_str())))
      {
	cout<<line<<" doesn't exist.  Skipping..."<<endl;
	continue;
      }
      TFile f(line.c_str(),"READ");
      if(f.Read("FDMC_Tau")==0 || f.Read("paramtree")==0)
      {
	cout<<"FDMC_Tau directory or paramtree doesn't exist.  Skipping..."<<endl;
	continue;
      }
      
      ExHists_FDMC->TrueVsReco_Fid_NuTauCC->Add((TH2D*)f.Get("FDMC_Tau/TrueVsReco_Fid_NuTauCC"));
      ExHists_FDMC->True_Fid_NuTauCC->Add((TH1D*)f.Get("FDMC_Tau/True_Fid_NuTauCC"));
      ExHists_FDMC->TrueVsReco_Fid_NuTauBarCC->Add((TH2D*)f.Get("FDMC_Tau/TrueVsReco_Fid_NuTauBarCC"));
      ExHists_FDMC->True_Fid_NuTauBarCC->Add((TH1D*)f.Get("FDMC_Tau/True_Fid_NuTauBarCC"));
      
      for(j=0;j<np;j++)
      {
        TrueVsRecoVsPID_NC_Tau[j]->Add((TH3D*)f.Get(Form("FDMC_Tau/TrueVsRecoVs%s_NC",PIDlist[j].c_str())));
        ExHists_FDMC->TrueVsRecoVsPID_NuMuToNuTau[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Tau/TrueVsRecoVs%s_NuMuToNuTau",PIDlist[j].c_str())));
        ExHists_FDMC->TrueVsRecoVsPID_BNueToNuTau[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Tau/TrueVsRecoVs%s_BNueToNuTau",PIDlist[j].c_str())));
      }
      
      paramtree_ = (TTree*)f.Get("paramtree");
      paramtree_->SetBranchAddress("POT",&npot);
      paramtree_->GetEntry(0);
      tot_ftmc+=npot;
    }
    n++;
  }
  filelist5.close();
  
  cout<<"FDMC Nu (Tau):"<<endl;
  filelist5.clear();
  filelist5.open(gSystem->ExpandPathName(FileList_FDMC_Tau.c_str()));
  n=0;
  while(!filelist5.eof())
  {
    cout<<n<<endl;
    getline(filelist5,line);
    if(!line.empty())
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(line.c_str())))
      {
	cout<<line<<" doesn't exist.  Skipping..."<<endl;
	continue;
      }
      TFile f(line.c_str(),"READ");
      if(f.Read("FDMC_Tau_Nu")==0 || f.Read("paramtree")==0)
      {
	cout<<"FDMC_Tau_Nu directory or paramtree doesn't exist.  Skipping..."<<endl;
	continue;
      }
      
      for(j=0;j<np;j++)
      {
	TrueVsRecoVsPID_Nu_NC_Tau[j]->Add((TH3D*)f.Get(Form("FDMC_Tau_Nu/TrueVsRecoVs%s_NC",PIDlist[j].c_str())));
        ExHists_Nu_FDMC->TrueVsRecoVsPID_NuMuToNuTau[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Tau_Nu/TrueVsRecoVs%s_NuMuToNuTau",PIDlist[j].c_str())));
        ExHists_Nu_FDMC->TrueVsRecoVsPID_BNueToNuTau[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Tau_Nu/TrueVsRecoVs%s_BNueToNuTau",PIDlist[j].c_str())));
      }
    }
    n++;
  }
  filelist5.close();
  
  cout<<"FDMC NuBar (Tau):"<<endl;
  filelist5.clear();
  filelist5.open(gSystem->ExpandPathName(FileList_FDMC_Tau.c_str()));
  n=0;
  while(!filelist5.eof())
  {
    cout<<n<<endl;
    getline(filelist5,line);
    if(!line.empty())
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(line.c_str())))
      {
	cout<<line<<" doesn't exist.  Skipping..."<<endl;
	continue;
      }
      TFile f(line.c_str(),"READ");
      if(f.Read("FDMC_Tau_NuBar")==0 || f.Read("paramtree")==0)
      {
	cout<<"FDMC_Tau_NuBar directory or paramtree doesn't exist.  Skipping..."<<endl;
	continue;
      }
      
      for(j=0;j<np;j++)
      {
	TrueVsRecoVsPID_NuBar_NC_Tau[j]->Add((TH3D*)f.Get(Form("FDMC_Tau_NuBar/TrueVsRecoVs%s_NC",PIDlist[j].c_str())));
        ExHists_NuBar_FDMC->TrueVsRecoVsPID_NuMuToNuTau[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Tau_NuBar/TrueVsRecoVs%s_NuMuToNuTau",PIDlist[j].c_str())));
        ExHists_NuBar_FDMC->TrueVsRecoVsPID_BNueToNuTau[PIDlist[j]]->Add((TH3D*)f.Get(Form("FDMC_Tau_NuBar/TrueVsRecoVs%s_BNueToNuTau",PIDlist[j].c_str())));
      }
    }
    n++;
  }
  filelist5.close();
  
  ExHists_FDMC->TrueVsReco_Fid_NuTauCC->Scale(fFarPOT/tot_ftmc);
  ExHists_FDMC->True_Fid_NuTauCC->Scale(fFarPOT/tot_ftmc);
  ExHists_FDMC->TrueVsReco_Fid_NuTauBarCC->Scale(fFarPOT/tot_ftmc);
  ExHists_FDMC->True_Fid_NuTauBarCC->Scale(fFarPOT/tot_ftmc);
  for(j=0;j<np;j++)
  {
    TrueVsRecoVsPID_NC_Tau[j]->Scale(fFarPOT/tot_ftmc);
    ExHists_FDMC->TrueVsRecoVsPID_NuMuToNuTau[PIDlist[j]]->Scale(fFarPOT/tot_ftmc);
    ExHists_FDMC->TrueVsRecoVsPID_BNueToNuTau[PIDlist[j]]->Scale(fFarPOT/tot_ftmc);
    
    TrueVsRecoVsPID_Nu_NC_Tau[j]->Scale(fFarPOT/tot_ftmc);
    ExHists_Nu_FDMC->TrueVsRecoVsPID_NuMuToNuTau[PIDlist[j]]->Scale(fFarPOT/tot_ftmc);
    ExHists_Nu_FDMC->TrueVsRecoVsPID_BNueToNuTau[PIDlist[j]]->Scale(fFarPOT/tot_ftmc);
    
    TrueVsRecoVsPID_NuBar_NC_Tau[j]->Scale(fFarPOT/tot_ftmc);
    ExHists_NuBar_FDMC->TrueVsRecoVsPID_NuMuToNuTau[PIDlist[j]]->Scale(fFarPOT/tot_ftmc);
    ExHists_NuBar_FDMC->TrueVsRecoVsPID_BNueToNuTau[PIDlist[j]]->Scale(fFarPOT/tot_ftmc);
  }
  
  for(j=0;j<np;j++)
  {
    ExHists_FDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Add(TrueVsRecoVsPID_NC_Beam[j]);
    ExHists_FDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Add(TrueVsRecoVsPID_NC_Nue[j]);
    ExHists_FDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Add(TrueVsRecoVsPID_NC_Tau[j]);
    
    ExHists_Nu_FDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Add(TrueVsRecoVsPID_Nu_NC_Beam[j]);
    ExHists_Nu_FDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Add(TrueVsRecoVsPID_Nu_NC_Nue[j]);
    ExHists_Nu_FDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Add(TrueVsRecoVsPID_Nu_NC_Tau[j]);
    
    ExHists_NuBar_FDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Add(TrueVsRecoVsPID_NuBar_NC_Beam[j]);
    ExHists_NuBar_FDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Add(TrueVsRecoVsPID_NuBar_NC_Nue[j]);
    ExHists_NuBar_FDMC->TrueVsRecoVsPID_NC[PIDlist[j]]->Add(TrueVsRecoVsPID_NuBar_NC_Tau[j]);
  }
  
  DeltaCP = del;
  DeltaMSq12 = dm2_12;
  DeltaMSq23 = dm2_23;
  Theta12 = t12;
  Theta13 = t13;
  Theta23 = t23;
  
  WriteToFile();
  
  return;
}
