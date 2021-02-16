#include "NueAna/MultiBinAna/ExtrapHistsGen.h"

ExtrapHistsGen::ExtrapHistsGen()
{
  Initialize();
  return;
}
void ExtrapHistsGen::Initialize()
{
  SetOutputFile();
  SetFileList();
  SetDataType();
  
  fout = 0;
  paramtree = 0;
  ExHists = 0;
  
  NueStandard::SetDefaultOscParam();
  
  nReco=0;
  nTrue=0;
  nPID=0;
  
  isRHC=false;
  
  return;
}
void ExtrapHistsGen::AddPID(string pidname,Selection::Selection_t pidsel)
{
  PIDlist[pidname] = pidsel;
  return;
}
void ExtrapHistsGen::SetFileList(string list)
{
  FileList = list;
  return;
}
void ExtrapHistsGen::SetDataType(string type)
{
  DataType = type;
  return;
}
void ExtrapHistsGen::SetOutputFile(std::string name)
{
   outFileName = name;
   
   return;
}
void ExtrapHistsGen::SetRecoBins(int n, double *edges)
{
  int i;
  
  nReco = n;
  
  RecoEdges.clear();
  for(i=0;i<nReco+1;i++)
  {
    RecoEdges.push_back(edges[i]);
  }
  
  return;
}
void ExtrapHistsGen::SetTrueBins(int n, double *edges)
{
  int i;
  
  nTrue = n;
  
  TrueEdges.clear();
  for(i=0;i<nTrue+1;i++)
  {
    TrueEdges.push_back(edges[i]);
  }
  
  return;
}
void ExtrapHistsGen::SetPIDBins(int n, double *edges)
{
  int i;
  
  nPID = n;
  
  PIDEdges.clear();
  for(i=0;i<nPID+1;i++)
  {
    PIDEdges.push_back(edges[i]);
  }
  
  return;
}
void ExtrapHistsGen::InitializeHistograms()
{
  gDirectory->cd("/");
  map<string,Selection::Selection_t>::iterator piditer;
  map<string,Selection::Selection_t>::iterator last;
  string p;
  
  int i;
  double *redges=0,*pedges=0,*tedges=0; 
  if(nReco!=0)
  {
    redges = new double[nReco+1];
    for(i=0;i<nReco+1;i++)
    {
      redges[i] = RecoEdges.at(i);
    }
  }
  if(nPID!=0)
  {
    pedges = new double[nPID+1];
    for(i=0;i<nPID+1;i++)
    {
      pedges[i] = PIDEdges.at(i);
    }
  }
  if(nTrue!=0)
  {
    tedges = new double[nTrue+1];
    for(i=0;i<nTrue+1;i++)
    {
      tedges[i] = TrueEdges.at(i);
    }
  }
  
  ExHists = new ExtrapHists2D(Form("ExtrapHists_%s",DataType.c_str()));
  if(nReco!=0) ExHists->SetRecoBins(nReco,redges);
  if(nPID!=0) ExHists->SetPIDBins(nPID,pedges);
  if(nTrue!=0) ExHists->SetTrueBins(nTrue,tedges);
  ExHists->fDirectory->cd();
  piditer = PIDlist.begin();
  last = PIDlist.end();
  while(piditer!=last)
  {
    p = piditer->first;
    ExHists->AddPID(p);
    piditer++;
  }
  ExHists->InitializeHists();
  
  ExHists_Nu = new ExtrapHists2D(Form("ExtrapHists_Nu_%s",DataType.c_str()));
  if(nReco!=0) ExHists_Nu->SetRecoBins(nReco,redges);
  if(nPID!=0) ExHists_Nu->SetPIDBins(nPID,pedges);
  if(nTrue!=0) ExHists_Nu->SetTrueBins(nTrue,tedges);
  ExHists_Nu->fDirectory->cd();
  piditer = PIDlist.begin();
  last = PIDlist.end();
  while(piditer!=last)
  {
    p = piditer->first;
    ExHists_Nu->AddPID(p);
    piditer++;
  }
  ExHists_Nu->InitializeHists();
  
  ExHists_NuBar = new ExtrapHists2D(Form("ExtrapHists_NuBar_%s",DataType.c_str()));
  if(nReco!=0) ExHists_NuBar->SetRecoBins(nReco,redges);
  if(nPID!=0) ExHists_NuBar->SetPIDBins(nPID,pedges);
  if(nTrue!=0) ExHists_NuBar->SetTrueBins(nTrue,tedges);
  ExHists_NuBar->fDirectory->cd();
  piditer = PIDlist.begin();
  last = PIDlist.end();
  while(piditer!=last)
  {
    p = piditer->first;
    ExHists_NuBar->AddPID(p);
    piditer++;
  }
  ExHists_NuBar->InitializeHists();
  
  return;
}
void ExtrapHistsGen::FillHistograms()
{
  if(FileList=="NULL")
  {
    cout<<"No file list. Quitting..."<<endl;
    return;
  }
  
  NueRecord* nr = new NueRecord();
  TChain *chain = new TChain("ana_nue");
  chain->SetBranchAddress("NueRecord",&nr);
  
  NuePOT *np = new NuePOT();
  TChain *pchain = new TChain("pottree");
  pchain->SetBranchAddress("NuePOT",&np);
  
  Int_t i;
  double totpot=0;
  Int_t n;
  ifstream filelist;
  string line;
  
  filelist.open(gSystem->ExpandPathName(FileList.c_str()));
  if(!filelist.good())
  {
    cout<<"Problem reading "<<FileList<<". Quitting..."<<endl;
    return;
  }
  while(!filelist.eof())
  {
    getline(filelist,line);
    if(!line.empty())
    {
      chain->Add(gSystem->ExpandPathName(line.c_str()));
      pchain->Add(gSystem->ExpandPathName(line.c_str()));
    }
  }
  filelist.close();
  
  n = pchain->GetEntries();
  totpot=0;
  for(i=0;i<n;i++)
  {
    pchain->GetEntry(i);
    totpot+=np->pot;
  }
  totpot = totpot*1e12;
  cout<<totpot<<" POT"<<endl;
  
  int nevt = chain->GetEntries();
  cout<<nevt<<" events"<<endl;
  if(nevt==0)
  {
    cout<<"No events.  Quitting..."<<endl;
    return;
  }
  
  Double_t weight=1;
  Double_t oscweight=1;
  bool isMC=false;
  Double_t CCEn = 0;
  Double_t NueEn = 0;
  Detector::Detector_t det=Detector::kUnknown;
  
  map<string,Selection::Selection_t>::iterator piditer;
  map<string,Selection::Selection_t>::iterator last;
  string pidname;
  Selection::Selection_t pidsel;
  
  string dt="NULL";
  chain->GetEvent(0);
  isMC = (nr->GetHeader().GetVldContext().GetSimFlag()==SimFlag::kMC);
  det = nr->GetHeader().GetVldContext().GetDetector();
  if(!isMC && det==Detector::kNear) dt="NDData";
  if(isMC && det==Detector::kNear) dt="NDMC";
  if(isMC && det==Detector::kFar)
  {
    if(nr->mctrue.nuFlavor==-9999 && nevt>1) chain->GetEvent(1);//try a couple of more entries in case the nuFlavor isn't set correctly in the first entry
    if(nr->mctrue.nuFlavor==-9999 && nevt>5) chain->GetEvent(5);
    if((TMath::Abs(nr->mctrue.nuFlavor)==14 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==14) || (TMath::Abs(nr->mctrue.nuFlavor)==12 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==12))//beam MC
    {
      dt="FDMC_Beam";
    }
    else if((TMath::Abs(nr->mctrue.nuFlavor)==14 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==12) || (TMath::Abs(nr->mctrue.nuFlavor)==12 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==14))
    {
      dt="FDMC_Nue";
    }
    else if((TMath::Abs(nr->mctrue.nuFlavor)==16 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==14) || (TMath::Abs(nr->mctrue.nuFlavor)==16 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==12))
    {
      dt="FDMC_Tau";
    }
    else
    {
      cout<<"Unknown data type.  ("<<TMath::Abs(nr->mctrue.nonOscNuFlavor)<<", "<<TMath::Abs(nr->mctrue.nuFlavor)<<")  Quitting..."<<endl;
      return;
    }
  }
  SetDataType(dt);
  
  InitializeHistograms();
  ExHists->SetPOT(totpot);
  ExHists_Nu->SetPOT(totpot);
  ExHists_NuBar->SetPOT(totpot);
  
  for(i=0;i<nevt;i++)
//   for(i=0;i<1000;i++)
  {
    chain->GetEvent(i);
    if(i%1000==0)
    {
      cout<<i<<endl;
    }
    
    isMC = (nr->GetHeader().GetVldContext().GetSimFlag()==SimFlag::kMC);
    det = nr->GetHeader().GetVldContext().GetDetector();
    
    if(!NueStandard::PassesSelection(nr,Selection::kFid)) continue;
    
    NueConvention::NueEnergyCorrection(nr);
    if(isRHC) NueConvention::RHCNueEnergyCorrection(nr);
    
    weight = 1;
    if(isMC) weight = NueStandard::GetMCWeights(nr);
    if(!isMC && det == Detector::kNear) weight = NueStandard::GetNDDataWeights(nr);
    
    oscweight = 1;
    if(det == Detector::kFar && isMC)//oscillation reweighting
    {
      oscweight = NueStandard::GetOscWeight(nr->mctrue.nuFlavor,nr->mctrue.nonOscNuFlavor,nr->mctrue.nuEnergy);
    }
    
    CCEn = GetCCRecoEnergy(nr);
    NueEn = GetNueRecoEnergy(nr);
    
    if(!isMC && det == Detector::kNear)//near data
    {
      if(NueStandard::PassesSelection(nr,Selection::kCC))
      {
        ExHists->Reco_CClike->Fill(CCEn,weight);
	if(nr->srtrack.fitMomentum>0)
	{
	  ExHists->Reco_CClike_Pos->Fill(CCEn,weight);
	}
	else
	{
	  ExHists->Reco_CClike_Neg->Fill(CCEn,weight);
	}
      }
      continue;
    }
    
    if(isMC && det == Detector::kNear)//near MC
    {
      ExHists->TrueVsReco_Fid->Fill(CCEn,nr->mctrue.nuEnergy,weight);
      if(nr->mctrue.interactionType>0 && nr->mctrue.nuFlavor==14)//numuCC
      {
        ExHists->TrueVsReco_Fid_NuMuCC->Fill(CCEn,nr->mctrue.nuEnergy,weight);
      }
      if(nr->mctrue.interactionType>0 && nr->mctrue.nuFlavor==-14)//numubarCC
      {
        ExHists->TrueVsReco_Fid_NuMuBarCC->Fill(CCEn,nr->mctrue.nuEnergy,weight);
      }
    
      if(NueStandard::PassesSelection(nr, Selection::kCC))
      {
        ExHists->TrueVsReco_CClike->Fill(CCEn,nr->mctrue.nuEnergy,weight);
	if(nr->srtrack.fitMomentum>0)
	{
	  ExHists->TrueVsReco_CClike_Pos->Fill(CCEn,nr->mctrue.nuEnergy,weight);
	}
	else
	{
	  ExHists->TrueVsReco_CClike_Neg->Fill(CCEn,nr->mctrue.nuEnergy,weight);
	}
	
        if(nr->mctrue.interactionType>0 && TMath::Abs(nr->mctrue.nuFlavor)==14)//numu/numubar CC
        {
          ExHists->TrueVsReco_CClike_NuMuCC->Fill(CCEn,nr->mctrue.nuEnergy,weight);
        }
        if(nr->srtrack.fitMomentum>0 && nr->mctrue.interactionType>0 && nr->mctrue.nuFlavor==-14 )//pos track, numubarCC
	{
	  ExHists->TrueVsReco_CClike_Pos_NuMuBarCC->Fill(CCEn,nr->mctrue.nuEnergy,weight);
	}
	if(nr->srtrack.fitMomentum<0 && nr->mctrue.interactionType>0 && nr->mctrue.nuFlavor==14 )//neg track, numuCC
	{
	  ExHists->TrueVsReco_CClike_Neg_NuMuCC->Fill(CCEn,nr->mctrue.nuEnergy,weight);
	}
      }
      
      if(NueStandard::PassesSelection(nr, Selection::kPre))
      {
        piditer = PIDlist.begin();
        last = PIDlist.end();
        while(piditer!=last)
        {
          pidname = piditer->first;
          pidsel = piditer->second;
          if(nr->mctrue.interactionType==0)//NC
          {
            ExHists->TrueVsRecoVsPID_NC[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight);
          }
          else//CC
          {
            if(TMath::Abs(nr->mctrue.nuFlavor)==14)//numu
            {
              ExHists->TrueVsRecoVsPID_NuMuToNuMu[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight);
            }
            if(TMath::Abs(nr->mctrue.nuFlavor)==12)//beam nue
            {
              ExHists->TrueVsRecoVsPID_BNueToBNue[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight);
            }
          }
          piditer++;
        }
      }
      
      continue;
    }
    
    if(isMC && det == Detector::kFar)
    {
      if((TMath::Abs(nr->mctrue.nuFlavor)==14 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==14) || (TMath::Abs(nr->mctrue.nuFlavor)==12 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==12))//beam MC only
      {
        ExHists->TrueVsReco_Fid->Fill(CCEn,nr->mctrue.nuEnergy,weight);
        if(nr->mctrue.interactionType>0 && nr->mctrue.nuFlavor==14)//numuCC
        {
          ExHists->TrueVsReco_Fid_NuMuCC->Fill(CCEn,nr->mctrue.nuEnergy,weight);
        }
        if(nr->mctrue.interactionType>0 && nr->mctrue.nuFlavor==-14)//numubarCC
        {
          ExHists->TrueVsReco_Fid_NuMuBarCC->Fill(CCEn,nr->mctrue.nuEnergy,weight);
        }
    
        if(NueStandard::PassesSelection(nr, Selection::kCC))
        {
          ExHists->TrueVsReco_CClike->Fill(CCEn,nr->mctrue.nuEnergy,weight);
	  if(nr->srtrack.fitMomentum>0)
	  {
	    ExHists->TrueVsReco_CClike_Pos->Fill(CCEn,nr->mctrue.nuEnergy,weight);
	  }
	  else
	  {
	    ExHists->TrueVsReco_CClike_Neg->Fill(CCEn,nr->mctrue.nuEnergy,weight);
	  }
	  
          if(nr->mctrue.interactionType>0 && TMath::Abs(nr->mctrue.nuFlavor)==14)//numu/numubar CC
          {
            ExHists->TrueVsReco_CClike_NuMuCC->Fill(CCEn,nr->mctrue.nuEnergy,weight);
          }
          if(nr->srtrack.fitMomentum>0 && nr->mctrue.interactionType>0 && nr->mctrue.nuFlavor==-14 )//pos track, numubarCC
	  {
	    ExHists->TrueVsReco_CClike_Pos_NuMuBarCC->Fill(CCEn,nr->mctrue.nuEnergy,weight);
	  }
	  if(nr->srtrack.fitMomentum<0 && nr->mctrue.interactionType>0 && nr->mctrue.nuFlavor==14 )//neg track, numuCC
	  {
	    ExHists->TrueVsReco_CClike_Neg_NuMuCC->Fill(CCEn,nr->mctrue.nuEnergy,weight);
	  }
        }
      }
      
      if(TMath::Abs(nr->mctrue.nuFlavor)==12 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==14 && nr->mctrue.interactionType>0)//oscillated nueCC only
      {
	if(nr->mctrue.nuFlavor==12)
	{
          ExHists->TrueVsReco_Fid_NueCC->Fill(NueEn,nr->mctrue.nuEnergy,weight*oscweight);
          ExHists->True_Fid_NueCC->Fill(nr->mctrue.nuEnergy,weight*oscweight);
	}
	if(nr->mctrue.nuFlavor==-12)
	{
          ExHists->TrueVsReco_Fid_NueBarCC->Fill(NueEn,nr->mctrue.nuEnergy,weight*oscweight);
          ExHists->True_Fid_NueBarCC->Fill(nr->mctrue.nuEnergy,weight*oscweight);
	}
      }
      
      if(TMath::Abs(nr->mctrue.nuFlavor)==16 && nr->mctrue.interactionType>0)//oscillated nutauCC only
      {
	if(nr->mctrue.nuFlavor==16)
	{
          ExHists->TrueVsReco_Fid_NuTauCC->Fill(NueEn,nr->mctrue.nuEnergy,weight*oscweight);
          ExHists->True_Fid_NuTauCC->Fill(nr->mctrue.nuEnergy,weight*oscweight);
	}
	if(nr->mctrue.nuFlavor==-16)
	{
          ExHists->TrueVsReco_Fid_NuTauBarCC->Fill(NueEn,nr->mctrue.nuEnergy,weight*oscweight);
          ExHists->True_Fid_NuTauBarCC->Fill(nr->mctrue.nuEnergy,weight*oscweight);
	}
      }
      
      if(NueStandard::PassesSelection(nr, Selection::kPre))
      {
        piditer = PIDlist.begin();
        last = PIDlist.end();
        while(piditer!=last)
        {
          pidname = piditer->first;
          pidsel = piditer->second;
          if(nr->mctrue.interactionType==0)//NC
          {
            ExHists->TrueVsRecoVsPID_NC[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
	    
	    if((TMath::Abs(nr->mctrue.nuFlavor)==14 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==14) || (TMath::Abs(nr->mctrue.nuFlavor)==12 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==12))//beam MC only
	    {
	      if(nr->mctrue.nonOscNuFlavor==14)
	      {
	        ExHists->TrueVsRecoVsPID_NC_NoOscNuMu[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight);
	      }
	      if(nr->mctrue.nonOscNuFlavor==-14)
	      {
	        ExHists->TrueVsRecoVsPID_NC_NoOscNuMuBar[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight);
	      }
	      if(nr->mctrue.nonOscNuFlavor==12)
	      {
	        ExHists->TrueVsRecoVsPID_NC_NoOscNue[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight);
	      }
	      if(nr->mctrue.nonOscNuFlavor==-12)
	      {
	        ExHists->TrueVsRecoVsPID_NC_NoOscNueBar[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight);
	      }
	    }
          }
          else//CC
          {
            if(TMath::Abs(nr->mctrue.nuFlavor)==14 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==14)//numu->numu
            {
              ExHists->TrueVsRecoVsPID_NuMuToNuMu[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
	      if(nr->mctrue.nuFlavor>0)  ExHists_Nu->TrueVsRecoVsPID_NuMuToNuMu[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
	      else ExHists_NuBar->TrueVsRecoVsPID_NuMuToNuMu[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
            }
            if(TMath::Abs(nr->mctrue.nuFlavor)==14 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==12)//beam nue->numu
            {
              ExHists->TrueVsRecoVsPID_BNueToNuMu[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
	      if(nr->mctrue.nuFlavor>0) ExHists_Nu->TrueVsRecoVsPID_BNueToNuMu[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
	      else ExHists_NuBar->TrueVsRecoVsPID_BNueToNuMu[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
            }
            if(TMath::Abs(nr->mctrue.nuFlavor)==12 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==12)//beam nue->beam nue
            {
              ExHists->TrueVsRecoVsPID_BNueToBNue[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
	      if(nr->mctrue.nuFlavor>0) ExHists_Nu->TrueVsRecoVsPID_BNueToBNue[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
	      else ExHists_NuBar->TrueVsRecoVsPID_BNueToBNue[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
            }
            if(TMath::Abs(nr->mctrue.nuFlavor)==12 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==14)//numu->nue
            {
              ExHists->TrueVsRecoVsPID_NuMuToNue[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
	      if(nr->mctrue.nuFlavor>0) ExHists_Nu->TrueVsRecoVsPID_NuMuToNue[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
	      else ExHists_NuBar->TrueVsRecoVsPID_NuMuToNue[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
            }
            if(TMath::Abs(nr->mctrue.nuFlavor)==16 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==14)//numu->nutau
            {
              ExHists->TrueVsRecoVsPID_NuMuToNuTau[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
	      if(nr->mctrue.nuFlavor>0) ExHists_Nu->TrueVsRecoVsPID_NuMuToNuTau[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
	      else ExHists_NuBar->TrueVsRecoVsPID_NuMuToNuTau[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
            }
            if(TMath::Abs(nr->mctrue.nuFlavor)==16 && TMath::Abs(nr->mctrue.nonOscNuFlavor)==12)//beam nue->nutau
            {
              ExHists->TrueVsRecoVsPID_BNueToNuTau[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
	      if(nr->mctrue.nuFlavor>0) ExHists_Nu->TrueVsRecoVsPID_BNueToNuTau[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
	      else ExHists_NuBar->TrueVsRecoVsPID_BNueToNuTau[pidname]->Fill(NueStandard::GetPIDValue(nr,pidsel),NueEn,nr->mctrue.nuEnergy,weight*oscweight);
            }
          }
          piditer++;
        }
      }
    }
  }
  
  WriteToFile();
  
  return;
  
}
void ExtrapHistsGen::FinalizeHistograms()
{
  int ir,it;
  int npid,nreco,ntrue;
  double r,p,t;
  
  npid = ExHists->GetPIDBins();
  nreco = ExHists->GetRecoBins();
  ntrue = ExHists->GetTrueBins();
  
  map<std::string, TH3D*>::iterator piditer;
  map<std::string, TH3D*>::iterator lastpid;
  string pidname;
  
  ExHists->fDirectory->cd();
  piditer = ExHists->TrueVsRecoVsPID_NC.begin();
  lastpid = ExHists->TrueVsRecoVsPID_NC.end();
  while(piditer != lastpid)
  {
    pidname = piditer->first;
    p = ExHists->TrueVsRecoVsPID_NC[pidname]->GetXaxis()->GetBinCenter(npid);//value of pid in last bin
    
    for(it=0;it<ntrue;it++)
    {
      for(ir=0;ir<nreco;ir++)
      {
        r = ExHists->TrueVsRecoVsPID_NC[pidname]->GetYaxis()->GetBinCenter(ir+1);
        t = ExHists->TrueVsRecoVsPID_NC[pidname]->GetZaxis()->GetBinCenter(it+1);
        
        ExHists->TrueVsRecoVsPID_NC[pidname]->Fill(p,r,t,ExHists->TrueVsRecoVsPID_NC[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists->TrueVsRecoVsPID_NC[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists->TrueVsRecoVsPID_NuMuToNuMu[pidname]->Fill(p,r,t,ExHists->TrueVsRecoVsPID_NuMuToNuMu[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists->TrueVsRecoVsPID_NuMuToNuMu[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists->TrueVsRecoVsPID_BNueToNuMu[pidname]->Fill(p,r,t,ExHists->TrueVsRecoVsPID_BNueToNuMu[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists->TrueVsRecoVsPID_BNueToNuMu[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists->TrueVsRecoVsPID_BNueToBNue[pidname]->Fill(p,r,t,ExHists->TrueVsRecoVsPID_BNueToBNue[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists->TrueVsRecoVsPID_BNueToBNue[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists-> TrueVsRecoVsPID_NuMuToNue[pidname]->Fill(p,r,t,ExHists->TrueVsRecoVsPID_NuMuToNue[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists->TrueVsRecoVsPID_NuMuToNue[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists->TrueVsRecoVsPID_NuMuToNuTau[pidname]->Fill(p,r,t,ExHists->TrueVsRecoVsPID_NuMuToNuTau[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists->TrueVsRecoVsPID_NuMuToNuTau[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists->TrueVsRecoVsPID_BNueToNuTau[pidname]->Fill(p,r,t,ExHists->TrueVsRecoVsPID_BNueToNuTau[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists->TrueVsRecoVsPID_BNueToNuTau[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
	
	ExHists->TrueVsRecoVsPID_NC_NoOscNue[pidname]->Fill(p,r,t,ExHists->TrueVsRecoVsPID_NC_NoOscNue[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists->TrueVsRecoVsPID_NC_NoOscNue[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
	
	ExHists->TrueVsRecoVsPID_NC_NoOscNuMu[pidname]->Fill(p,r,t,ExHists->TrueVsRecoVsPID_NC_NoOscNuMu[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists->TrueVsRecoVsPID_NC_NoOscNuMu[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
	
	ExHists->TrueVsRecoVsPID_NC_NoOscNueBar[pidname]->Fill(p,r,t,ExHists->TrueVsRecoVsPID_NC_NoOscNueBar[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists->TrueVsRecoVsPID_NC_NoOscNueBar[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
	
	ExHists->TrueVsRecoVsPID_NC_NoOscNuMuBar[pidname]->Fill(p,r,t,ExHists->TrueVsRecoVsPID_NC_NoOscNuMuBar[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists->TrueVsRecoVsPID_NC_NoOscNuMuBar[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
      }
    }
    
    piditer++;
  }
  
  ExHists_Nu->fDirectory->cd();
  piditer = ExHists_Nu->TrueVsRecoVsPID_NC.begin();
  lastpid = ExHists_Nu->TrueVsRecoVsPID_NC.end();
  while(piditer != lastpid)
  {
    pidname = piditer->first;
    p = ExHists_Nu->TrueVsRecoVsPID_NC[pidname]->GetXaxis()->GetBinCenter(npid);//value of pid in last bin
    
    for(it=0;it<ntrue;it++)
    {
      for(ir=0;ir<nreco;ir++)
      {
        r = ExHists_Nu->TrueVsRecoVsPID_NC[pidname]->GetYaxis()->GetBinCenter(ir+1);
        t = ExHists_Nu->TrueVsRecoVsPID_NC[pidname]->GetZaxis()->GetBinCenter(it+1);
        
        ExHists_Nu->TrueVsRecoVsPID_NC[pidname]->Fill(p,r,t,ExHists_Nu->TrueVsRecoVsPID_NC[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_Nu->TrueVsRecoVsPID_NC[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists_Nu->TrueVsRecoVsPID_NuMuToNuMu[pidname]->Fill(p,r,t,ExHists_Nu->TrueVsRecoVsPID_NuMuToNuMu[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_Nu->TrueVsRecoVsPID_NuMuToNuMu[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists_Nu->TrueVsRecoVsPID_BNueToNuMu[pidname]->Fill(p,r,t,ExHists_Nu->TrueVsRecoVsPID_BNueToNuMu[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_Nu->TrueVsRecoVsPID_BNueToNuMu[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists_Nu->TrueVsRecoVsPID_BNueToBNue[pidname]->Fill(p,r,t,ExHists_Nu->TrueVsRecoVsPID_BNueToBNue[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_Nu->TrueVsRecoVsPID_BNueToBNue[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists_Nu-> TrueVsRecoVsPID_NuMuToNue[pidname]->Fill(p,r,t,ExHists_Nu->TrueVsRecoVsPID_NuMuToNue[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_Nu->TrueVsRecoVsPID_NuMuToNue[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists_Nu->TrueVsRecoVsPID_NuMuToNuTau[pidname]->Fill(p,r,t,ExHists_Nu->TrueVsRecoVsPID_NuMuToNuTau[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_Nu->TrueVsRecoVsPID_NuMuToNuTau[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists_Nu->TrueVsRecoVsPID_BNueToNuTau[pidname]->Fill(p,r,t,ExHists_Nu->TrueVsRecoVsPID_BNueToNuTau[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_Nu->TrueVsRecoVsPID_BNueToNuTau[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
      }
    }
    
    piditer++;
  }
  
  ExHists_NuBar->fDirectory->cd();
  piditer = ExHists_NuBar->TrueVsRecoVsPID_NC.begin();
  lastpid = ExHists_NuBar->TrueVsRecoVsPID_NC.end();
  while(piditer != lastpid)
  {
    pidname = piditer->first;
    p = ExHists_NuBar->TrueVsRecoVsPID_NC[pidname]->GetXaxis()->GetBinCenter(npid);//value of pid in last bin
    
    for(it=0;it<ntrue;it++)
    {
      for(ir=0;ir<nreco;ir++)
      {
        r = ExHists_NuBar->TrueVsRecoVsPID_NC[pidname]->GetYaxis()->GetBinCenter(ir+1);
        t = ExHists_NuBar->TrueVsRecoVsPID_NC[pidname]->GetZaxis()->GetBinCenter(it+1);
        
        ExHists_NuBar->TrueVsRecoVsPID_NC[pidname]->Fill(p,r,t,ExHists_NuBar->TrueVsRecoVsPID_NC[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_NuBar->TrueVsRecoVsPID_NC[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists_NuBar->TrueVsRecoVsPID_NuMuToNuMu[pidname]->Fill(p,r,t,ExHists_NuBar->TrueVsRecoVsPID_NuMuToNuMu[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_NuBar->TrueVsRecoVsPID_NuMuToNuMu[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists_NuBar->TrueVsRecoVsPID_BNueToNuMu[pidname]->Fill(p,r,t,ExHists_NuBar->TrueVsRecoVsPID_BNueToNuMu[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_NuBar->TrueVsRecoVsPID_BNueToNuMu[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists_NuBar->TrueVsRecoVsPID_BNueToBNue[pidname]->Fill(p,r,t,ExHists_NuBar->TrueVsRecoVsPID_BNueToBNue[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_NuBar->TrueVsRecoVsPID_BNueToBNue[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists_NuBar-> TrueVsRecoVsPID_NuMuToNue[pidname]->Fill(p,r,t,ExHists_NuBar->TrueVsRecoVsPID_NuMuToNue[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_NuBar->TrueVsRecoVsPID_NuMuToNue[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists_NuBar->TrueVsRecoVsPID_NuMuToNuTau[pidname]->Fill(p,r,t,ExHists_NuBar->TrueVsRecoVsPID_NuMuToNuTau[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_NuBar->TrueVsRecoVsPID_NuMuToNuTau[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
        ExHists_NuBar->TrueVsRecoVsPID_BNueToNuTau[pidname]->Fill(p,r,t,ExHists_NuBar->TrueVsRecoVsPID_BNueToNuTau[pidname]->GetBinContent(npid+1,ir+1,it+1));//fill last pid bin with overflow content
        ExHists_NuBar->TrueVsRecoVsPID_BNueToNuTau[pidname]->SetBinContent(npid+1,ir+1,it+1,0);//set overflow content to 0
        
      }
    }
    
    piditer++;
  }
  
  return;
}
void ExtrapHistsGen::WriteToFile()
{
  FinalizeHistograms();
  
  double par[100]={0};
  NueStandard::GetOscParam(par);
  Theta12 = par[OscPar::kTh12];
  Theta13 = par[OscPar::kTh13];
  Theta23 = par[OscPar::kTh23];
  DeltaMSq23 = par[OscPar::kDeltaM23];
  DeltaMSq12 = par[OscPar::kDeltaM12];
  DeltaCP = par[OscPar::kDelta];
  
  double p = ExHists->GetPOT();
  
  fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  paramtree = new TTree("paramtree","paramtree");
  paramtree->Branch("POT",&p,"POT/D");
  paramtree->Branch("Theta12",&Theta12,"Theta12/D");
  paramtree->Branch("Theta13",&Theta13,"Theta13/D");
  paramtree->Branch("Theta23",&Theta23,"Theta23/D");
  paramtree->Branch("DeltaMSq23",&DeltaMSq23,"DeltaMSq23/D");
  paramtree->Branch("DeltaMSq12",&DeltaMSq12,"DeltaMSq12/D");
  paramtree->Branch("DeltaCP",&DeltaCP,"DeltaCP/D");
  paramtree->Fill();
  
  fout->cd();
  TDirectory *dir = fout->mkdir(DataType.c_str());
  dir->cd();
  TList *list = ExHists->fDirectory->GetList();
  TIter iter(list->MakeIterator());
  TObject *ob = 0;
  while((ob = iter())) 
  {
    ob->Write();
  }
  
  string n;
  
  fout->cd();
  n = DataType + "_Nu";
  TDirectory *dirnu = fout->mkdir(n.c_str());
  dirnu->cd();
  TList *listnu = ExHists_Nu->fDirectory->GetList();
  TIter iternu(listnu->MakeIterator());
  TObject *obnu = 0;
  while((obnu = iternu())) 
  {
    obnu->Write();
  }
  
  fout->cd();
  n = DataType + "_NuBar";
  TDirectory *dirnubar = fout->mkdir(n.c_str());
  dirnubar->cd();
  TList *listnubar = ExHists_NuBar->fDirectory->GetList();
  TIter iternubar(listnubar->MakeIterator());
  TObject *obnubar = 0;
  while((obnubar = iternubar())) 
  {
    obnubar->Write();
  }
  
  fout->cd();
  paramtree->Write();
  fout->Close();
}
double ExtrapHistsGen::GetNueRecoEnergy(NueRecord *nr)
{
  return nr->srevent.phNueGeV;
}
double ExtrapHistsGen::GetCCRecoEnergy(NueRecord *nr)
{
  double CCEn = 0.0;
  if(nr->srtrack.phCCGeV > 0) CCEn += nr->srtrack.phCCGeV;
  if(nr->srshower.phCCGeV > 0) CCEn += nr->srshower.phCCGeV;
  
  return CCEn;
}
