#include "NueAna/MultiBinAna/SysHists2D.h"

SysHists2D::SysHists2D(std::string name)
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

SysHists2D::~SysHists2D()
{
}
void SysHists2D::SetRecoBins(int n, double *edges)
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
void SysHists2D::SetTrueBins(int n, double *edges)
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
void SysHists2D::SetPIDBins(int n, double *edges)
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
void SysHists2D::AddPID(std::string pid)
{
  PIDlist.push_back(pid);
  
  return;
}
void SysHists2D::InitializeHists()
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
  
  ND_TrueVsReco = new TH2D("ND_TrueVsReco","",nReco,redges,nTrue,tedges);
  FD_TrueVsReco = new TH2D("FD_TrueVsReco","",nReco,redges,nTrue,tedges);
  ND_TrueVsReco_ErrorHist = new TH2D("ND_TrueVsReco_ErrorHist","",nReco,redges,nTrue,tedges);
  FD_TrueVsReco_ErrorHist = new TH2D("FD_TrueVsReco_ErrorHist","",nReco,redges,nTrue,tedges);
  
  ND_RecoVsPID.clear();
  FD_RecoVsPID.clear();
  ND_RecoVsPID_ErrorHist.clear();
  FD_RecoVsPID_ErrorHist.clear();
  
  TH2D *h;
  int np = PIDlist.size();
  for(i=0;i<np;i++)
  {
    h = new TH2D(Form("ND_RecoVs%s",PIDlist[i].c_str()),"",nPID,pedges,nReco,redges);
    ND_RecoVsPID[PIDlist[i]] = h;
    
    h = new TH2D(Form("FD_RecoVs%s",PIDlist[i].c_str()),"",nPID,pedges,nReco,redges);
    FD_RecoVsPID[PIDlist[i]] = h;
    
    h = new TH2D(Form("ND_RecoVs%s_ErrorHist",PIDlist[i].c_str()),"",nPID,pedges,nReco,redges);
    ND_RecoVsPID_ErrorHist[PIDlist[i]] = h;
    
    h = new TH2D(Form("FD_RecoVs%s_ErrorHist",PIDlist[i].c_str()),"",nPID,pedges,nReco,redges);
    FD_RecoVsPID_ErrorHist[PIDlist[i]] = h;
  }
  
  ND_TrueVsReco->Sumw2();
  FD_TrueVsReco->Sumw2();
  ND_TrueVsReco_ErrorHist->Sumw2();
  FD_TrueVsReco_ErrorHist->Sumw2();
  for(i=0;i<np;i++)
  {
    ND_RecoVsPID[PIDlist[i]]->Sumw2();
    FD_RecoVsPID[PIDlist[i]]->Sumw2();
    ND_RecoVsPID_ErrorHist[PIDlist[i]]->Sumw2();
    FD_RecoVsPID_ErrorHist[PIDlist[i]]->Sumw2();
  }
  
  delete [] redges;
  delete [] tedges;
  delete [] pedges;
  
  return;
}
