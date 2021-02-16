#include "NueAna/MultiBinAna/SysFileGen2D.h"

SysFileGen2D::SysFileGen2D()
{
  Initialize();
  return;
}
void SysFileGen2D::Initialize()
{
  SetOutputFile();
  SetNearPOT();
  SetFarPOT();
  SetDefaultOsc();
  SetCutLevel();
  
  fout = 0;
  paramtree = 0;
  
  return;
}
void SysFileGen2D::AddPID(std::string pid)
{
  PIDlist.push_back(pid);
  return;
}

void SysFileGen2D::SetDefaultOsc() 
{
  NueStandard::SetDefaultOscParam();
  return; 
}

void SysFileGen2D::RunSystematics()
{
  for(unsigned int j = 0; j < fSysNames.size(); j++)
  {
    std::cout<<fSysNames[j]<<endl;
    fCurrentSysName = fSysNames[j];
    ResetHistograms();  //Clear out whatever happened last time
    FillHistograms();
    WriteToFile();
  }
}

void SysFileGen2D::FillHistograms()
{
}

void SysFileGen2D::SetOutputFile(std::string name)
{
   outFileName = name;
   
   return;
}

void SysFileGen2D::FinalizeHistograms()
{
  std::map<Background::Background_t, SysHists2D*>::iterator bgiter = fSysHistsMap.begin();
  std::map<Background::Background_t, SysHists2D*>::iterator last  = fSysHistsMap.end();
  
  int ir;
  int npid,nreco;
  double r,p;
  double c1,c2;
  
  while(bgiter != last)
  {
    SysHists2D* sysh = bgiter->second;
    
    sysh->fDirectory->cd();
    
    npid = sysh->GetPIDBins();
    nreco = sysh->GetRecoBins();

    std::map<std::string, TH2D*>::iterator piditer = sysh->ND_RecoVsPID.begin();
    std::map<std::string, TH2D*>::iterator lastpid = sysh->ND_RecoVsPID.end();
    
    while(piditer != lastpid)
    {
      TH2D *h = piditer->second;
      p = h->GetXaxis()->GetBinCenter(npid);
      for(ir=0;ir<nreco;ir++)
      {
        r = h->GetYaxis()->GetBinCenter(ir+1);
        h->Fill(p,r,h->GetBinContent(npid+1,ir+1));//fill last pid bin with overflow content
	h->SetBinContent(npid+1,ir+1,0);//set overflow content to 0
      }
      piditer++;
    }
    
    piditer = sysh->FD_RecoVsPID.begin();
    lastpid = sysh->FD_RecoVsPID.end();
    
    while(piditer != lastpid)
    {
      TH2D *h = piditer->second;
      p = h->GetXaxis()->GetBinCenter(npid);
      for(ir=0;ir<nreco;ir++)
      {
        r = h->GetYaxis()->GetBinCenter(ir+1);
        h->Fill(p,r,h->GetBinContent(npid+1,ir+1));//fill last pid bin with overflow content
	h->SetBinContent(npid+1,ir+1,0);//set overflow content to 0
      }
      piditer++;
    }
    
    piditer = sysh->ND_RecoVsPID_ErrorHist.begin();
    lastpid = sysh->ND_RecoVsPID_ErrorHist.end();
    
    while(piditer != lastpid)
    {
      TH2D *h = piditer->second;
      for(ir=0;ir<nreco;ir++)
      {
        c1 = h->GetBinContent(npid,ir+1);//content of last pid bin
        c2 = h->GetBinContent(npid+1,ir+1);//content of overflow bin
        h->SetBinContent(npid,ir+1,sqrt(c1*c1 + c2*c2));//add content of overflow bin to last pid bin (sum in quadrature because this is an error histogram)
        h->SetBinContent(npid+1,ir+1,0);//set overflow content to 0
      }
      piditer++;
    }
    
    piditer = sysh->FD_RecoVsPID_ErrorHist.begin();
    lastpid = sysh->FD_RecoVsPID_ErrorHist.end();
    
    while(piditer != lastpid)
    {
      TH2D *h = piditer->second;
      
      for(ir=0;ir<nreco;ir++)
      {
        c1 = h->GetBinContent(npid,ir+1);//content of last pid bin
        c2 = h->GetBinContent(npid+1,ir+1);//content of overflow bin
        h->SetBinContent(npid,ir+1,sqrt(c1*c1 + c2*c2));//add content of overflow bin to last pid bin (sum in quadrature because this is an error histogram)
        h->SetBinContent(npid+1,ir+1,0);//set overflow content to 0
      }
      piditer++;
    }
    
    bgiter++;
  }
  
  return;
}
void SysFileGen2D::WriteToFile()
{
  FinalizeHistograms();
  
  char selection[256];
  
  double par[100]={0};
  NueStandard::GetOscParam(par);
  Theta12 = par[OscPar::kTh12];
  Theta13 = par[OscPar::kTh13];
  Theta23 = par[OscPar::kTh23];
  DeltaMSq23 = par[OscPar::kDeltaM23];
  DeltaMSq12 = par[OscPar::kDeltaM12];
  DeltaCP = par[OscPar::kDelta];
  
  if(fout == 0)
  {
    fout = new TFile(outFileName.c_str(),"RECREATE");
    paramtree = new TTree("paramtree","paramtree");
    paramtree->Branch("Selection",selection,"Selection/C");
    paramtree->Branch("nearPOT",&fNearPOT,"nearPOT/D");
    paramtree->Branch("farPOT",&fFarPOT,"farPOT/D");
    paramtree->Branch("Theta12",&Theta12,"Theta12/D");
    paramtree->Branch("Theta13",&Theta13,"Theta13/D");
    paramtree->Branch("Theta23",&Theta23,"Theta23/D");
    paramtree->Branch("DeltaMSq23",&DeltaMSq23,"DeltaMSq23/D");
    paramtree->Branch("DeltaMSq12",&DeltaMSq12,"DeltaMSq12/D");
    paramtree->Branch("DeltaCP",&DeltaCP,"DeltaCP/D");
  }
  
  fout->cd();
  
  sprintf(selection,"%s",Selection::AsString(CutLevel));

  std::map<Background::Background_t,SysHists2D*>::iterator bgiter = fSysHistsMap.begin();
  std::map<Background::Background_t,SysHists2D*>::iterator last = fSysHistsMap.end();

  while(bgiter!=last) {
    std::string sysh_name = bgiter->second->fDirectory->GetName();
    sysh_name += "_" + fCurrentSysName + "_" + std::string(Selection::AsString(CutLevel));

    TDirectory *filedir = fout->mkdir(sysh_name.c_str());
    filedir->cd();
    TList *list = bgiter->second->fDirectory->GetList();
    TIter iter(list->MakeIterator());
    TObject *ob = 0;
    while((ob = iter())) ob->Write();
    fout->cd();
    bgiter++;
  }

  fout->cd();

  paramtree->Fill();

  if(fCurrentSysName == fSysNames[fSysNames.size()-1])
  {
    paramtree->Write();
    fout->Close();
  }
}

void SysFileGen2D::InitializeHistograms()
{

  vector<Background::Background_t> bgs;
  bgs.push_back(Background::kNC);
  bgs.push_back(Background::kNuMuCC);
  bgs.push_back(Background::kNuTauCC);
  bgs.push_back(Background::kNueCC);
  bgs.push_back(Background::kBNueCC);

  for(unsigned int i = 0; i < bgs.size(); i++){
    gDirectory->cd("/");
    std::string sysh_name = std::string(Background::AsString(bgs[i]));
    SysHists2D *sysh = new SysHists2D(sysh_name.c_str());
    sysh->fDirectory->cd();
    for(unsigned int j=0;j<PIDlist.size();j++)
    {
      sysh->AddPID(PIDlist.at(j));
    }
    sysh->InitializeHists();
 
    (fSysHistsMap)[bgs[i]] = sysh;
  }
}
void SysFileGen2D::ResetHistograms()
{

  if(fSysHistsMap.size() == 0) InitializeHistograms();
  
  std::map<Background::Background_t, SysHists2D*>::iterator bgiter = fSysHistsMap.begin();
  std::map<Background::Background_t, SysHists2D*>::iterator last  = fSysHistsMap.end();
  
  while(bgiter != last)
  {
    SysHists2D* sysh = bgiter->second;
    
    sysh->fDirectory->cd();
    
    sysh->ND_TrueVsReco->Reset("ICE");
    sysh->FD_TrueVsReco->Reset("ICE");
    sysh->ND_TrueVsReco_ErrorHist->Reset("ICE");
    sysh->FD_TrueVsReco_ErrorHist->Reset("ICE");
    
    std::map<std::string, TH2D*>::iterator piditer = sysh->ND_RecoVsPID.begin();
    std::map<std::string, TH2D*>::iterator lastpid = sysh->ND_RecoVsPID.end();
    while(piditer != lastpid)
    {
      TH2D *h = piditer->second;
      h->Reset("ICE");
      piditer++;
    }
    piditer = sysh->FD_RecoVsPID.begin();
    lastpid = sysh->FD_RecoVsPID.end();
    while(piditer != lastpid)
    {
      TH2D *h = piditer->second;
      h->Reset("ICE");
      piditer++;
    }
    piditer = sysh->ND_RecoVsPID_ErrorHist.begin();
    lastpid = sysh->ND_RecoVsPID_ErrorHist.end();
    while(piditer != lastpid)
    {
      TH2D *h = piditer->second;
      h->Reset("ICE");
      piditer++;
    }
    piditer = sysh->FD_RecoVsPID_ErrorHist.begin();
    lastpid = sysh->FD_RecoVsPID_ErrorHist.end();
    while(piditer != lastpid)
    {
      TH2D *h = piditer->second;
      h->Reset("ICE");
      piditer++;
    }
    
    bgiter++;
  }
}
double SysFileGen2D::GetNueRecoEnergy(NueRecord *nr)
{
  return nr->srevent.phNueGeV;
}
double SysFileGen2D::GetCCRecoEnergy(NueRecord *nr)
{
  double CCEn = 0.0;
  if(nr->srtrack.phCCGeV > 0) CCEn += nr->srtrack.phCCGeV;
  if(nr->srshower.phCCGeV > 0) CCEn += nr->srshower.phCCGeV;
  
  return CCEn;
}
