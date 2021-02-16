#include "NueAna/MultiBinAna/ExtrapHists2D.h"

ExtrapHists2D::ExtrapHists2D(std::string name)
{
  fDirectory = new TDirectory(name.c_str(),name.c_str());
  
  int i;
  
  nReco = 100;
  for(i=0;i<nReco+1;i++)
  {
    RecoEdges.push_back(i*1.);
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
}

ExtrapHists2D::~ExtrapHists2D()
{
}
void ExtrapHists2D::SetRecoBins(int n, double *edges)
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
void ExtrapHists2D::SetTrueBins(int n, double *edges)
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
void ExtrapHists2D::SetPIDBins(int n, double *edges)
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
void ExtrapHists2D::SetPOT(double p)
{
  POT = p;
  return;
}
void ExtrapHists2D::AddPID(std::string pid)
{
  PIDlist.push_back(pid);
  
  return;
}
void ExtrapHists2D::InitializeHists()
{
  double *redges = new double[nReco+1];
  double *tedges = new double[nTrue+1];
  double *pedges = new double[nPID+1];
  int i;
  for(i=0;i<nReco+1;i++)
  {
    redges[i] = RecoEdges.at(i);
  }
  for(i=0;i<nTrue+1;i++)
  {
    tedges[i] = TrueEdges.at(i);
  }
  for(i=0;i<nPID+1;i++)
  {
    pedges[i] = PIDEdges.at(i);
  }
  
  unsigned int np = PIDlist.size();
  unsigned int j;
  
  TrueVsReco_Fid = new TH2D("TrueVsReco_Fid","",nReco,redges,nTrue,tedges);
  TrueVsReco_Fid_NuMuCC = new TH2D("TrueVsReco_Fid_NuMuCC","",nReco,redges,nTrue,tedges);
  TrueVsReco_Fid_NuMuBarCC = new TH2D("TrueVsReco_Fid_NuMuBarCC","",nReco,redges,nTrue,tedges);
  TrueVsReco_CClike = new TH2D("TrueVsReco_CClike","",nReco,redges,nTrue,tedges);
  TrueVsReco_CClike_Pos = new TH2D("TrueVsReco_CClike_Pos","",nReco,redges,nTrue,tedges);
  TrueVsReco_CClike_Neg = new TH2D("TrueVsReco_CClike_Neg","",nReco,redges,nTrue,tedges);
  TrueVsReco_CClike_NuMuCC = new TH2D("TrueVsReco_CClike_NuMuCC","",nReco,redges,nTrue,tedges);
  TrueVsReco_CClike_Pos_NuMuBarCC = new TH2D("TrueVsReco_CClike_Pos_NuMuBarCC","",nReco,redges,nTrue,tedges);
  TrueVsReco_CClike_Neg_NuMuCC = new TH2D("TrueVsReco_CClike_Neg_NuMuCC","",nReco,redges,nTrue,tedges);
  Reco_CClike = new TH1D("Reco_CClike","",nReco,redges);
  Reco_CClike_Pos = new TH1D("Reco_CClike_Pos","",nReco,redges);
  Reco_CClike_Neg = new TH1D("Reco_CClike_Neg","",nReco,redges);
  
  TrueVsReco_Fid_NueCC = new TH2D("TrueVsReco_Fid_NueCC","",nReco,redges,nTrue,tedges);
  TrueVsReco_Fid_NuTauCC = new TH2D("TrueVsReco_Fid_NuTauCC","",nReco,redges,nTrue,tedges);
  True_Fid_NueCC = new TH1D("True_Fid_NueCC","",nTrue,tedges);
  True_Fid_NuTauCC = new TH1D("True_Fid_NuTauCC","",nTrue,tedges);
  TrueVsReco_Fid_NueBarCC = new TH2D("TrueVsReco_Fid_NueBarCC","",nReco,redges,nTrue,tedges);
  TrueVsReco_Fid_NuTauBarCC = new TH2D("TrueVsReco_Fid_NuTauBarCC","",nReco,redges,nTrue,tedges);
  True_Fid_NueBarCC = new TH1D("True_Fid_NueBarCC","",nTrue,tedges);
  True_Fid_NuTauBarCC = new TH1D("True_Fid_NuTauBarCC","",nTrue,tedges);
  
  for(j=0;j<np;j++)
  {
    TrueVsRecoVsPID_NC[PIDlist[j]] = new TH3D(Form("TrueVsRecoVs%s_NC",PIDlist[j].c_str()),"",nPID,pedges,nReco,redges,nTrue,tedges);
    TrueVsRecoVsPID_NuMuToNuMu[PIDlist[j]] = new TH3D(Form("TrueVsRecoVs%s_NuMuToNuMu",PIDlist[j].c_str()),"",nPID,pedges,nReco,redges,nTrue,tedges);
    TrueVsRecoVsPID_BNueToNuMu[PIDlist[j]] = new TH3D(Form("TrueVsRecoVs%s_BNueToNuMu",PIDlist[j].c_str()),"",nPID,pedges,nReco,redges,nTrue,tedges);
    TrueVsRecoVsPID_BNueToBNue[PIDlist[j]] = new TH3D(Form("TrueVsRecoVs%s_BNueToBNue",PIDlist[j].c_str()),"",nPID,pedges,nReco,redges,nTrue,tedges);
    TrueVsRecoVsPID_NuMuToNue[PIDlist[j]] = new TH3D(Form("TrueVsRecoVs%s_NuMuToNue",PIDlist[j].c_str()),"",nPID,pedges,nReco,redges,nTrue,tedges);
    TrueVsRecoVsPID_NuMuToNuTau[PIDlist[j]] = new TH3D(Form("TrueVsRecoVs%s_NuMuToNuTau",PIDlist[j].c_str()),"",nPID,pedges,nReco,redges,nTrue,tedges);
    TrueVsRecoVsPID_BNueToNuTau[PIDlist[j]] = new TH3D(Form("TrueVsRecoVs%s_BNueToNuTau",PIDlist[j].c_str()),"",nPID,pedges,nReco,redges,nTrue,tedges);
    
    TrueVsRecoVsPID_NC_NoOscNue[PIDlist[j]] = new TH3D(Form("TrueVsRecoVs%s_NC_NoOscNue",PIDlist[j].c_str()),"",nPID,pedges,nReco,redges,nTrue,tedges);
    TrueVsRecoVsPID_NC_NoOscNuMu[PIDlist[j]] = new TH3D(Form("TrueVsRecoVs%s_NC_NoOscNuMu",PIDlist[j].c_str()),"",nPID,pedges,nReco,redges,nTrue,tedges);
    TrueVsRecoVsPID_NC_NoOscNueBar[PIDlist[j]] = new TH3D(Form("TrueVsRecoVs%s_NC_NoOscNueBar",PIDlist[j].c_str()),"",nPID,pedges,nReco,redges,nTrue,tedges);
    TrueVsRecoVsPID_NC_NoOscNuMuBar[PIDlist[j]] = new TH3D(Form("TrueVsRecoVs%s_NC_NoOscNuMuBar",PIDlist[j].c_str()),"",nPID,pedges,nReco,redges,nTrue,tedges);
  }
  
  TrueVsReco_Fid->Sumw2();
  TrueVsReco_Fid_NuMuCC->Sumw2();
  TrueVsReco_Fid_NuMuBarCC->Sumw2();
  TrueVsReco_CClike->Sumw2();
  TrueVsReco_CClike_Pos->Sumw2();
  TrueVsReco_CClike_Neg->Sumw2();
  TrueVsReco_CClike_NuMuCC->Sumw2();
  TrueVsReco_CClike_Pos_NuMuBarCC->Sumw2();
  TrueVsReco_CClike_Neg_NuMuCC->Sumw2();
  Reco_CClike->Sumw2();
  Reco_CClike_Pos->Sumw2();
  Reco_CClike_Neg->Sumw2();
  TrueVsReco_Fid_NueCC->Sumw2();
  TrueVsReco_Fid_NuTauCC->Sumw2();
  True_Fid_NueCC->Sumw2();
  True_Fid_NuTauCC->Sumw2();
  TrueVsReco_Fid_NueBarCC->Sumw2();
  TrueVsReco_Fid_NuTauBarCC->Sumw2();
  True_Fid_NueBarCC->Sumw2();
  True_Fid_NuTauBarCC->Sumw2();
  for(j=0;j<np;j++)
  {
    TrueVsRecoVsPID_NC[PIDlist[j]]->Sumw2();
    TrueVsRecoVsPID_NuMuToNuMu[PIDlist[j]]->Sumw2();
    TrueVsRecoVsPID_BNueToNuMu[PIDlist[j]]->Sumw2();
    TrueVsRecoVsPID_BNueToBNue[PIDlist[j]]->Sumw2();
    TrueVsRecoVsPID_NuMuToNue[PIDlist[j]]->Sumw2();
    TrueVsRecoVsPID_NuMuToNuTau[PIDlist[j]]->Sumw2();
    TrueVsRecoVsPID_BNueToNuTau[PIDlist[j]]->Sumw2();
    TrueVsRecoVsPID_NC_NoOscNue[PIDlist[j]]->Sumw2();
    TrueVsRecoVsPID_NC_NoOscNuMu[PIDlist[j]]->Sumw2();
    TrueVsRecoVsPID_NC_NoOscNueBar[PIDlist[j]]->Sumw2();
    TrueVsRecoVsPID_NC_NoOscNuMuBar[PIDlist[j]]->Sumw2();
  }
  
  delete [] redges;
  delete [] tedges;
  delete [] pedges;
  
  return;
}
