#define ErrorCalc_Joint_C

#include "NueAna/MultiBinAna/ErrorCalc_Joint.h"

ErrorCalc_Joint::ErrorCalc_Joint()
{
  Init = false;
  InitSys = false;
  InitDecompSys = false;
  UseGrid = false;
  
  return;
}
ErrorCalc_Joint::~ErrorCalc_Joint()
{
}

void ErrorCalc_Joint::AddExtrap(Extrapolate2D* /*E*/)
{
  cout << "For Joint Fit, use AddExtrapRHC or AddExtrapFHC instead.  Doing nothing..." << endl;
  return;
}

void ErrorCalc_Joint::AddExtrapFHC(Extrapolate2D *E)
{
  Extrap.push_back(E);
  //Make sure to indicate the Extrapolation type (0 for FHC):
  ExtrapType.push_back(0);
  return;
}
void ErrorCalc_Joint::AddExtrapRHC(Extrapolate2D *E)
{
  Extrap.push_back(E);
  //Make sure to indicate the Extrapolation type (1 for RHC):
  ExtrapType.push_back(1);
  return;
}


void ErrorCalc_Joint::Initialize()
{
  if (Init) return;
  
  int nbinsFHC = 0;
  int nbinsRHC = 0;

  cout.precision(5);
  cout<<fixed;
  
  unsigned int ie;
  for(ie=0;ie<Extrap.size();ie++)
  {
    Extrap[ie]->GetPrediction();
    if(Extrap[ie]->GetOscFlag())
    {
      cout<<"Error in ErrorCalc_Joint::Initialize(): You have called Extrapolate2D::OscillatePrediction() before initializing the systematics!  Quitting..."<<endl;
      return;
    }
    //Get number of bins for each type:
    if (ExtrapType[ie]==0) nbinsFHC = Extrap[ie]->GetNReco()*Extrap[ie]->GetNPID();
    if (ExtrapType[ie]==1) nbinsRHC = Extrap[ie]->GetNReco()*Extrap[ie]->GetNPID();
  }
  
  gROOT->cd();
  int nbins = nbinsFHC + nbinsRHC;//Used to be Extrap[0]->GetNReco()*Extrap[0]->GetNPID();
  CovMatrix = new TH2D("CovMatrix","",nbins,-0.5,nbins-0.5,nbins,-0.5,nbins-0.5);
  CovMatrix_Decomp = new TH2D("CovMatrix_Decomp","",nbins,-0.5,nbins-0.5,nbins,-0.5,nbins-0.5);
  

  //Double-check compatibility with joint fit later, but this should be fine for now
  string st;
  for(ie=0;ie<Extrap.size();ie++)
  {

    if (ExtrapType[ie]==0) nbins = nbinsFHC;
    if (ExtrapType[ie]==1) nbins = nbinsRHC;

      st = string(Background::AsString(Background::kNC));  
      GridPred[Background::kNC].push_back(new TH1D(Form("GridPred_%s_%i",st.c_str(),ie),"",nbins,-0.5,nbins-0.5));
      st = string(Background::AsString(Background::kNuMuCC));  
      GridPred[Background::kNuMuCC].push_back(new TH1D(Form("GridPred_%s_%i",st.c_str(),ie),"",nbins,-0.5,nbins-0.5));
      st = string(Background::AsString(Background::kBNueCC));  
      GridPred[Background::kBNueCC].push_back(new TH1D(Form("GridPred_%s_%i",st.c_str(),ie),"",nbins,-0.5,nbins-0.5));
      st = string(Background::AsString(Background::kNuTauCC));  
      GridPred[Background::kNuTauCC].push_back(new TH1D(Form("GridPred_%s_%i",st.c_str(),ie),"",nbins,-0.5,nbins-0.5));
      st = string(Background::AsString(Background::kNueCC));  
      GridPred[Background::kNueCC].push_back(new TH1D(Form("GridPred_%s_%i",st.c_str(),ie),"",nbins,-0.5,nbins-0.5));
   }

  Init = true;
  
  return;
}

void ErrorCalc_Joint::ReadSysFiles_FNExtrap_AllRuns(int n)
{

  int i;
  
  string plf,pln,mnf,mnn,std,systname,filetag;
  systname = FNExtrap_SystName.at(n);
  plf = FNExtrap_FarPlusTag.at(n);
  pln = FNExtrap_NearPlusTag.at(n);
  mnf = FNExtrap_FarMinusTag.at(n);
  mnn = FNExtrap_NearMinusTag.at(n);
  std = FNExtrap_StdTag.at(n);
  filetag=FNExtrap_FileTag.at(n);
  int flag = FNExtrap_Flag.at(n);
  
  unsigned int ie;

  vector<string> PreFiles;
  
  for(ie=0;ie<Extrap.size();ie++)
  {
    if (ExtrapType[ie]==0) PreFiles.push_back(SysFileDir+"/SysFile2D_FHC"+filetag+"_Pre_1.root");
    if (ExtrapType[ie]==1) PreFiles.push_back(SysFileDir+"/SysFile2D_RHC"+filetag+"_Pre_1.root");
  }

  unsigned int j;
  vector<Background::Background_t> bgs;
  bgs.push_back(Background::kNC);
  bgs.push_back(Background::kNuMuCC);
  bgs.push_back(Background::kBNueCC);
  
  string name;
  int np,nr;
  double *r,*p;
  TFile *fpre;
  TH2D *farstd=0,*farpl=0,*farmn=0;
  TH2D *nearstd=0,*nearpl=0,*nearmn=0;
  TH2D *totalfarstd[2]={0},*totalfarpl[2]={0},*totalfarmn[2]={0};
  TH2D *totalnearstd[2]={0},*totalnearpl[2]={0},*totalnearmn[2]={0};
  TH2D *fnpl[2]={0}, *fnmn[2]={0};
  double npot_pre,fpot_pre;
  TTree *treepre;
  double nPOTNear,nPOTFar;
  int nfiles[2]={0};
   
  //RBTNOTE: THIS SHOULD BE DONE SEPARATELY FOR RHC AND FHC:
  //Addendum: Why not just do it separately for each extrapolation?
  bool preExist = false;
  for (ie=0; ie<PreFiles.size(); ie++){
    if (!gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[ie].c_str()))) preExist = true;
  }

  //flag==1 means systematics related to numuCC extrapolation for signal/tau prediction - by definition 0 for these components
  if(flag==1 || (!preExist))//if all three Pre files do NOT exist
  {
    for(ie=0;ie<Extrap.size();ie++)
    {

      nr = Extrap[ie]->GetNReco();
      np = Extrap[ie]->GetNPID();
      r = new double[nr+1];
      p = new double[np+1];
      for(i=0;i<nr+1;i++)
	{
	  r[i] = Extrap[ie]->Pred_TotalBkgd->GetYaxis()->GetBinLowEdge(i+1);
	}
      for(i=0;i<np+1;i++)
	{
	  p[i] = Extrap[ie]->Pred_TotalBkgd->GetXaxis()->GetBinLowEdge(i+1);
	}

      FN_NC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FN_NC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      FN_NC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FN_NC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      
      FN_NuMuCC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FN_NuMuCC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      FN_NuMuCC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FN_NuMuCC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      
      FN_BNueCC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FN_BNueCC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      FN_BNueCC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FN_BNueCC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      
      ND_NC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_ND_NC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      ND_NC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_ND_NC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      
      ND_NuMuCC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_ND_NuMuCC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      ND_NuMuCC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_ND_NuMuCC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      
      ND_BNueCC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_ND_BNueCC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      ND_BNueCC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_ND_BNueCC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      
      FD_NC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FD_NC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      FD_NC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FD_NC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      
      FD_NuMuCC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FD_NuMuCC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      FD_NuMuCC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FD_NuMuCC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      
      FD_BNueCC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FD_BNueCC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      FD_BNueCC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FD_BNueCC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    }
    
    return;
  }

  int ietag=-1;
  
  for(j=0;j<bgs.size();j++)
  {

    //Set number of FHC and RHC files to 0
    nfiles[0]=0;
    nfiles[1]=0;

    for(ie=0;ie<Extrap.size();ie++)
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[ie].c_str())))//if file doesn't exist
      {
        continue;
      }

      //Get a tag of 0 for FHC, 1 for RHC
      ietag = ExtrapType[ie];

      nr = Extrap[ie]->GetNReco();
      np = Extrap[ie]->GetNPID();
      r = new double[nr+1];
      p = new double[np+1];
      for(i=0;i<nr+1;i++)
	{
	  r[i] = Extrap[ie]->Pred_TotalBkgd->GetYaxis()->GetBinLowEdge(i+1);
	}
      for(i=0;i<np+1;i++)
	{
	  p[i] = Extrap[ie]->Pred_TotalBkgd->GetXaxis()->GetBinLowEdge(i+1);
	}

      //Scaling:
      nPOTFar=Extrap[ie]->GetFarPOT();
      nPOTNear=Extrap[ie]->GetNearPOT();
      ///---

      fpre = new TFile(gSystem->ExpandPathName(PreFiles[ie].c_str()),"READ");
      
      treepre = (TTree*)fpre->Get("paramtree");
      treepre->SetBranchAddress("nearPOT",&npot_pre);
      treepre->SetBranchAddress("farPOT",&fpot_pre);
      treepre->GetEntry(0);
      
      name = string(Background::AsString(bgs[j])) + "_" + std + "_Presel/FD_RecoVs" + Extrap[ie]->GetPID();
      farstd = (TH2D*)fpre->Get(name.c_str());
      farstd->Scale(nPOTFar/fpot_pre);
      MBH.Rebin2DHist(farstd,np,p,nr,r);
      if(nfiles[ietag]==0) totalfarstd[ietag] = (TH2D*)farstd->Clone();
      else totalfarstd[ietag]->Add(farstd);
      
      name = string(Background::AsString(bgs[j])) + "_" + plf + "_Presel/FD_RecoVs" + Extrap[ie]->GetPID();
      if(plf==std)
      {
        farpl = (TH2D*)farstd->Clone(name.c_str());
      }
      else
      {
        farpl = (TH2D*)fpre->Get(name.c_str());
        farpl->Scale(nPOTFar/fpot_pre);
        MBH.Rebin2DHist(farpl,np,p,nr,r);
      }
      if(nfiles[ietag]==0) totalfarpl[ietag] = (TH2D*)farpl->Clone();
      else totalfarpl[ietag]->Add(farpl);
      
      name = string(Background::AsString(bgs[j])) + "_" + mnf + "_Presel/FD_RecoVs" + Extrap[ie]->GetPID();
      if(mnf==std)
      {
        farmn = (TH2D*)farstd->Clone(name.c_str());
      }
      else if(mnf==plf)
      {
        farmn = (TH2D*)farpl->Clone(name.c_str());
      }
      else
      {
        farmn = (TH2D*)fpre->Get(name.c_str());
        farmn->Scale(nPOTFar/fpot_pre);
        MBH.Rebin2DHist(farmn,np,p,nr,r);
      }
      if(nfiles[ietag]==0) totalfarmn[ietag] = (TH2D*)farmn->Clone();
      else totalfarmn[ietag]->Add(farmn);
      
      name = string(Background::AsString(bgs[j])) + "_" + std + "_Presel/ND_RecoVs" + Extrap[ie]->GetPID();
      nearstd = (TH2D*)fpre->Get(name.c_str());
      nearstd->Scale(nPOTNear/npot_pre);
      MBH.Rebin2DHist(nearstd,np,p,nr,r);
      if(nfiles[ietag]==0) totalnearstd[ietag] = (TH2D*)nearstd->Clone();
      else totalnearstd[ietag]->Add(nearstd);
      
      name = string(Background::AsString(bgs[j])) + "_" + pln + "_Presel/ND_RecoVs" + Extrap[ie]->GetPID();
      if(pln==std)
      {
        nearpl = (TH2D*)nearstd->Clone(name.c_str());
      }
      else
      {
        nearpl = (TH2D*)fpre->Get(name.c_str());
        nearpl->Scale(nPOTNear/npot_pre);
        MBH.Rebin2DHist(nearpl,np,p,nr,r);
      }
      if(nfiles[ietag]==0) totalnearpl[ietag] = (TH2D*)nearpl->Clone();
      else totalnearpl[ietag]->Add(nearpl);
      
      name = string(Background::AsString(bgs[j])) + "_" + mnn + "_Presel/ND_RecoVs" + Extrap[ie]->GetPID();
      if(mnn==std)
      {
        nearmn=(TH2D*)nearstd->Clone(name.c_str());
      }
      else if(mnn==pln)
      {
        nearmn=(TH2D*)nearpl->Clone(name.c_str());
      }
      else
      {
        nearmn = (TH2D*)fpre->Get(name.c_str());
        nearmn->Scale(nPOTNear/npot_pre);
        MBH.Rebin2DHist(nearmn,np,p,nr,r);
      }
      if(nfiles[ietag]==0) totalnearmn[ietag] = (TH2D*)nearmn->Clone();
      else totalnearmn[ietag]->Add(nearmn);
      
      nfiles[ietag]++;
    }
    
    for(ie=0;ie<Extrap.size();ie++)
    {

      //Get a tag of 0 for FHC, 1 for RHC
      ietag = ExtrapType[ie];

      if(bgs[j]==Background::kNC)
      {
        ND_NC_Plus1Sigma[systname].push_back((TH2D*)totalnearpl[ietag]->Clone(Form("%s_ND_NC_Plus1",systname.c_str())));
        ND_NC_Plus1Sigma[systname][ie]->Add(totalnearstd[ietag],-1.0);
        ND_NC_Plus1Sigma[systname][ie]->Divide(totalnearstd[ietag]);
        ND_NC_Minus1Sigma[systname].push_back((TH2D*)totalnearmn[ietag]->Clone(Form("%s_ND_NC_Minus1",systname.c_str())));
        ND_NC_Minus1Sigma[systname][ie]->Add(totalnearstd[ietag],-1.0);
        ND_NC_Minus1Sigma[systname][ie]->Divide(totalnearstd[ietag]);
        if(plf==mnf && pln==mnn)
        {
          ND_NC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
        
        FD_NC_Plus1Sigma[systname].push_back((TH2D*)totalfarpl[ietag]->Clone(Form("%s_FD_NC_Plus1",systname.c_str())));
        FD_NC_Plus1Sigma[systname][ie]->Add(totalfarstd[ietag],-1.0);
        FD_NC_Plus1Sigma[systname][ie]->Divide(totalfarstd[ietag]);
        FD_NC_Minus1Sigma[systname].push_back((TH2D*)totalfarmn[ietag]->Clone(Form("%s_FD_NC_Minus1",systname.c_str())));
        FD_NC_Minus1Sigma[systname][ie]->Add(totalfarstd[ietag],-1.0);
        FD_NC_Minus1Sigma[systname][ie]->Divide(totalfarstd[ietag]);
        if(plf==mnf && pln==mnn)
        {
          FD_NC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
      }
      else if(bgs[j]==Background::kNuMuCC)
      {
        ND_NuMuCC_Plus1Sigma[systname].push_back((TH2D*)totalnearpl[ietag]->Clone(Form("%s_ND_NuMuCC_Plus1",systname.c_str())));
        ND_NuMuCC_Plus1Sigma[systname][ie]->Add(totalnearstd[ietag],-1.0);
        ND_NuMuCC_Plus1Sigma[systname][ie]->Divide(totalnearstd[ietag]);
        ND_NuMuCC_Minus1Sigma[systname].push_back((TH2D*)totalnearmn[ietag]->Clone(Form("%s_ND_NuMuCC_Minus1",systname.c_str())));
        ND_NuMuCC_Minus1Sigma[systname][ie]->Add(totalnearstd[ietag],-1.0);
        ND_NuMuCC_Minus1Sigma[systname][ie]->Divide(totalnearstd[ietag]);
        if(plf==mnf && pln==mnn)
        {
          ND_NuMuCC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
        
        FD_NuMuCC_Plus1Sigma[systname].push_back((TH2D*)totalfarpl[ietag]->Clone(Form("%s_FD_NuMuCC_Plus1",systname.c_str())));
        FD_NuMuCC_Plus1Sigma[systname][ie]->Add(totalfarstd[ietag],-1.0);
        FD_NuMuCC_Plus1Sigma[systname][ie]->Divide(totalfarstd[ietag]);
        FD_NuMuCC_Minus1Sigma[systname].push_back((TH2D*)totalfarmn[ietag]->Clone(Form("%s_FD_NuMuCC_Minus1",systname.c_str())));
        FD_NuMuCC_Minus1Sigma[systname][ie]->Add(totalfarstd[ietag],-1.0);
        FD_NuMuCC_Minus1Sigma[systname][ie]->Divide(totalfarstd[ietag]);
        if(plf==mnf && pln==mnn)
        {
          FD_NuMuCC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
      }
      else if(bgs[j]==Background::kBNueCC)
      {
        ND_BNueCC_Plus1Sigma[systname].push_back((TH2D*)totalnearpl[ietag]->Clone(Form("%s_ND_BNueCC_Plus1",systname.c_str())));
        ND_BNueCC_Plus1Sigma[systname][ie]->Add(totalnearstd[ietag],-1.0);
        ND_BNueCC_Plus1Sigma[systname][ie]->Divide(totalnearstd[ietag]);
        ND_BNueCC_Minus1Sigma[systname].push_back((TH2D*)totalnearmn[ietag]->Clone(Form("%s_ND_BNueCC_Minus1",systname.c_str())));
        ND_BNueCC_Minus1Sigma[systname][ie]->Add(totalnearstd[ietag],-1.0);
        ND_BNueCC_Minus1Sigma[systname][ie]->Divide(totalnearstd[ietag]);
        if(plf==mnf && pln==mnn)
        {
          ND_BNueCC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
        
        FD_BNueCC_Plus1Sigma[systname].push_back((TH2D*)totalfarpl[ietag]->Clone(Form("%s_FD_BNueCC_Plus1",systname.c_str())));
        FD_BNueCC_Plus1Sigma[systname][ie]->Add(totalfarstd[ietag],-1.0);
        FD_BNueCC_Plus1Sigma[systname][ie]->Divide(totalfarstd[ietag]);
        FD_BNueCC_Minus1Sigma[systname].push_back((TH2D*)totalfarmn[ietag]->Clone(Form("%s_FD_BNueCC_Minus1",systname.c_str())));
        FD_BNueCC_Minus1Sigma[systname][ie]->Add(totalfarstd[ietag],-1.0);
        FD_BNueCC_Minus1Sigma[systname][ie]->Divide(totalfarstd[ietag]);
        if(plf==mnf && pln==mnn)
        {
          FD_BNueCC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
      }
    }
    
    //Make the total histograms for both FHC(itype=0) and RHC(itype=1)
    for (int itype=0; itype<2; itype++){
      if (nfiles[itype]>0){
	if (totalfarstd[itype]) totalfarstd[itype]->Divide(totalnearstd[itype]);
	if (totalfarpl[itype]) totalfarpl[itype]->Divide(totalnearpl[itype]);
	if (totalfarmn[itype]) totalfarmn[itype]->Divide(totalnearmn[itype]);
	
	fnpl[itype] = (TH2D*)totalfarpl[itype]->Clone("fnpl");
	fnpl[itype]->Add(totalfarstd[itype],-1.);
	fnpl[itype]->Divide(totalfarstd[itype]);
	
	fnmn[itype] = (TH2D*)totalfarmn[itype]->Clone("fnmn");
	fnmn[itype]->Add(totalfarstd[itype],-1.);
	fnmn[itype]->Divide(totalfarstd[itype]);
	
	if(plf==mnf && pln==mnn)
	  {
	    fnmn[itype]->Scale(-1.);
	  }
      }
    }

    for(ie=0;ie<Extrap.size();ie++)
    {

      //Get a tag of 0 for FHC, 1 for RHC
      ietag = ExtrapType[ie];

      if(bgs[j]==Background::kNC)
      {
        FN_NC_Plus1Sigma[systname].push_back((TH2D*)fnpl[ietag]->Clone(Form("%s_FN_NC_Plus1",systname.c_str())));
        FN_NC_Minus1Sigma[systname].push_back((TH2D*)fnmn[ietag]->Clone(Form("%s_FN_NC_Minus1",systname.c_str())));
      }
      if(bgs[j]==Background::kNuMuCC)
      {
        FN_NuMuCC_Plus1Sigma[systname].push_back((TH2D*)fnpl[ietag]->Clone(Form("%s_FN_NuMuCC_Plus1",systname.c_str())));
        FN_NuMuCC_Minus1Sigma[systname].push_back((TH2D*)fnmn[ietag]->Clone(Form("%s_FN_NuMuCC_Minus1",systname.c_str())));
      }
      if(bgs[j]==Background::kBNueCC)
      {
        FN_BNueCC_Plus1Sigma[systname].push_back((TH2D*)fnpl[ietag]->Clone(Form("%s_FN_BNueCC_Plus1",systname.c_str())));
        FN_BNueCC_Minus1Sigma[systname].push_back((TH2D*)fnmn[ietag]->Clone(Form("%s_FN_BNueCC_Minus1",systname.c_str())));
      }
    }
  }
  
  return;
}

void ErrorCalc_Joint::ReadSysFiles_Appearance_AllRuns(int n)
{

  int i;
  int ietag=-1;

  string plf,pln,mnf,mnn,std,systname,filetag;
  systname = FNExtrap_SystName.at(n);
  plf = FNExtrap_FarPlusTag.at(n);
  pln = FNExtrap_NearPlusTag.at(n);
  mnf = FNExtrap_FarMinusTag.at(n);
  mnn = FNExtrap_NearMinusTag.at(n);
  std = FNExtrap_StdTag.at(n);
  filetag=FNExtrap_FileTag.at(n);
  int flag = FNExtrap_Flag.at(n);
  
  vector<string> FidFiles;
  vector<string> CClikeFiles;
  vector<string> PreFiles;
  

  unsigned int ie;
  for(ie=0;ie<Extrap.size();ie++)
  {
  
    if(flag==1)//systs related to NuMuCC-like prediction
    {
    if (ExtrapType[ie]==0) FidFiles.push_back(SysFileDir+"/SysFile2D_FHC"+filetag+"_Fid_1.root");
    if (ExtrapType[ie]==1) FidFiles.push_back(SysFileDir+"/SysFile2D_RHC"+filetag+"_Fid_1.root");
    if (ExtrapType[ie]==0) CClikeFiles.push_back(SysFileDir+"/SysFile2D_FHC"+filetag+"_CCLike_1.root");
    if (ExtrapType[ie]==1) CClikeFiles.push_back(SysFileDir+"/SysFile2D_RHC"+filetag+"_CCLike_1.root");
    }
    else if(flag==2)//uncertainty in selected tau or signal FD MC
    {
    if (ExtrapType[ie]==0) PreFiles.push_back(SysFileDir+"/SysFile2D_FHC"+filetag+"_Pre_1.root");
    if (ExtrapType[ie]==1) PreFiles.push_back(SysFileDir+"/SysFile2D_RHC"+filetag+"_Pre_1.root");
    
    }
  }
  
  unsigned int j;
  vector<Background::Background_t> app;
  app.push_back(Background::kNueCC);
  app.push_back(Background::kNuTauCC);
  
  string name;
  int np,nr,nt;
  double *r,*p,*t;
  TFile *fpre,*ffid,*fcc;
  double npot_pre,fpot_pre;
  double npot_fid,fpot_fid;
  double npot_cc,fpot_cc;
  TTree *treepre,*treefid,*treecc;
  double nPOTNear,nPOTFar;
  TH2D *prestd=0,*prepl=0,*premn=0;
  TH2D *totalprestd[2]={0},*totalprepl[2]={0},*totalpremn[2]={0};
  double sum,temp;
  int ir,it,ip;
  TH1D *ccstd[2]={0},*ccpl[2]={0},*ccmn[2]={0};
  TH1D *data=0;
  TH2D *farfidccstd=0,*farcclikestd=0,*nearcclikestd=0;
  TH2D *farfidccpl=0,*farcclikepl=0,*nearcclikepl=0;
  TH2D *farfidccmn=0,*farcclikemn=0,*nearcclikemn=0;
  TH1D *totaldata[2]={0};
  TH2D *totalfarfidccstd[2]={0},*totalfarcclikestd[2]={0},*totalnearcclikestd[2]={0};
  TH2D *totalfarfidccpl[2]={0},*totalfarcclikepl[2]={0},*totalnearcclikepl[2]={0};
  TH2D *totalfarfidccmn[2]={0},*totalfarcclikemn[2]={0},*totalnearcclikemn[2]={0};
  TH2D *h2;
  int nfiles[2]={0};
  double sum_plus_tau,sum_minus_tau,sum_plus_nue,sum_minus_nue;
  string ccname = "CC";
    
  if(flag==0)
  {
    for(ie=0;ie<Extrap.size();ie++)
    {

      nr = Extrap[ie]->GetNReco();
      np = Extrap[ie]->GetNPID();
      nt = Extrap[ie]->GetNTrue();
      r = new double[nr+1];
      p = new double[np+1];
      t = new double[nt+1];
      for(i=0;i<nr+1;i++)
	{
	  r[i] = Extrap[ie]->Pred_TotalBkgd->GetYaxis()->GetBinLowEdge(i+1);
	}
      for(i=0;i<np+1;i++)
	{
	  p[i] = Extrap[ie]->Pred_TotalBkgd->GetXaxis()->GetBinLowEdge(i+1);
	}
      for(i=0;i<nt+1;i++)
	{
	  t[i] = Extrap[ie]->Pred_3D[Extrapolate2D::qNuMuToNuTau]->GetZaxis()->GetBinLowEdge(i+1);
	}


       //create empty histograms for overall signal/tau shifts
      NueCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NueCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NuTauCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NuTauCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      
       //creating empty numuCC shifts
      Pred_CC_Fid_Plus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Plus1",systname.c_str(),ie+1),"",nt,t));
      Pred_CC_Fid_Minus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Minus1",systname.c_str(),ie+1),"",nt,t));

    }
  }
  
  if(flag==1)
  {
    //create empty histograms for overall signal/tau shifts
    for(ie=0;ie<Extrap.size();ie++)
    {

      nr = Extrap[ie]->GetNReco();
      np = Extrap[ie]->GetNPID();
      r = new double[nr+1];
      p = new double[np+1];

      for(i=0;i<nr+1;i++)
	{
	  r[i] = Extrap[ie]->Pred_TotalBkgd->GetYaxis()->GetBinLowEdge(i+1);
	}
      for(i=0;i<np+1;i++)
	{
	  p[i] = Extrap[ie]->Pred_TotalBkgd->GetXaxis()->GetBinLowEdge(i+1);
	}
      NueCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NueCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NuTauCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NuTauCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    }
    

    //RBTNOTE: DO THIS SEPARATELY FOR RHC AND FHC!
    bool preExist = false;
    for (ie=0; ie<FidFiles.size(); ie++){
      if (!gSystem->AccessPathName(gSystem->ExpandPathName(FidFiles[ie].c_str()))) preExist = true;
    }
    for (ie=0; ie<CClikeFiles.size(); ie++){
      if (!gSystem->AccessPathName(gSystem->ExpandPathName(CClikeFiles[ie].c_str()))) preExist = true;
    }

    //If none of the files exists
    if (!preExist)
    {
      for(ie=0;ie<Extrap.size();ie++)
      {

	nt = Extrap[ie]->GetNTrue();
	t = new double[nt+1];
	for(i=0;i<nt+1;i++)
	  {
	    t[i] = Extrap[ie]->Pred_3D[Extrapolate2D::qNuMuToNuTau]->GetZaxis()->GetBinLowEdge(i+1);
	  }
	
        Pred_CC_Fid_Plus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Plus1",systname.c_str(),ie+1),"",nt,t));
        Pred_CC_Fid_Minus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Minus1",systname.c_str(),ie+1),"",nt,t));
      }
      
      return;
    }
    nfiles[0]=0;
    nfiles[1]=0;
    for(ie=0;ie<Extrap.size();ie++)
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(CClikeFiles[ie].c_str())) && gSystem->AccessPathName(gSystem->ExpandPathName(FidFiles[ie].c_str())))//if files don't exist
      { 
        continue;
      }
     
      ietag = ExtrapType[ie];

      //
      nPOTFar=Extrap[ie]->GetFarPOT();
      nPOTNear=Extrap[ie]->GetNearPOT();

      fcc = new TFile(gSystem->ExpandPathName(CClikeFiles[ie].c_str()),"READ");
      treecc = (TTree*)fcc->Get("paramtree");
      treecc->SetBranchAddress("nearPOT",&npot_cc);
      treecc->SetBranchAddress("farPOT",&fpot_cc);
      treecc->GetEntry(0);
      
      ffid = new TFile(gSystem->ExpandPathName(FidFiles[ie].c_str()),"READ");
      treefid = (TTree*)ffid->Get("paramtree");
      treefid->SetBranchAddress("nearPOT",&npot_fid);
      treefid->SetBranchAddress("farPOT",&fpot_fid);
      treefid->GetEntry(0);
      

      data = (TH1D*)Extrap[ie]->NDData_Reco_CClike->Clone(Form("data_%i",ie));
      //should scale data! don't need to?
      if(nfiles[ietag]==0) totaldata[ietag] = (TH1D*)data->Clone();
      else totaldata[ietag]->Add(data);
      

      name = "NuMuCC_" + std + "_Fid/FD_TrueVsReco";
      farfidccstd = (TH2D*)ffid->Get(name.c_str());
      farfidccstd->Scale(nPOTFar/fpot_fid);
      if(nfiles[ietag]==0) totalfarfidccstd[ietag] = (TH2D*)farfidccstd->Clone();
      else totalfarfidccstd[ietag]->Add(farfidccstd);
      
      name = "NuMuCC_" + std + "_" + ccname + "/FD_TrueVsReco";
      farcclikestd = (TH2D*)fcc->Get(name.c_str());
      name = "NC_" + std + "_" + ccname + "/FD_TrueVsReco";
      h2 = (TH2D*)fcc->Get(name.c_str());
      farcclikestd->Add(h2);
      name = "BNueCC_" + std + "_" + ccname + "/FD_TrueVsReco";
      h2 = (TH2D*)fcc->Get(name.c_str());
      farcclikestd->Add(h2);
      farcclikestd->Scale(nPOTFar/fpot_cc);
      if(nfiles[ietag]==0) totalfarcclikestd[ietag] = (TH2D*)farcclikestd->Clone();
      else totalfarcclikestd[ietag]->Add(farcclikestd);
      
      name = "NuMuCC_" + std + "_" + ccname + "/ND_TrueVsReco";
      nearcclikestd = (TH2D*)fcc->Get(name.c_str());
      name = "NC_" + std + "_" + ccname + "/ND_TrueVsReco";
      h2 = (TH2D*)fcc->Get(name.c_str());
      nearcclikestd->Add(h2);
      name = "BNueCC_" + std + "_" + ccname + "/ND_TrueVsReco";
      h2 = (TH2D*)fcc->Get(name.c_str());
      nearcclikestd->Add(h2);
      nearcclikestd->Scale(nPOTNear/npot_cc);
      if(nfiles[ietag]==0) totalnearcclikestd[ietag] = (TH2D*)nearcclikestd->Clone();
      else totalnearcclikestd[ietag]->Add(nearcclikestd);
      

      name = "NuMuCC_" + plf + "_Fid/FD_TrueVsReco";
      if(plf==std)
      {
        farfidccpl = (TH2D*)farfidccstd->Clone(name.c_str());
      }
      else
      {
        farfidccpl = (TH2D*)ffid->Get(name.c_str());
        farfidccpl->Scale(nPOTFar/fpot_fid);
      }
      if(nfiles[ietag]==0) totalfarfidccpl[ietag] = (TH2D*)farfidccpl->Clone();
      else totalfarfidccpl[ietag]->Add(farfidccpl);
      
      name = "NuMuCC_" + plf + "_" + ccname + "/FD_TrueVsReco";
      if(plf==std)
      {
        farcclikepl = (TH2D*)farcclikestd->Clone(name.c_str());
      }
      else
      {
        farcclikepl = (TH2D*)fcc->Get(name.c_str());
        name = "NC_" + plf + "_" + ccname + "/FD_TrueVsReco";
        h2 = (TH2D*)fcc->Get(name.c_str());
        farcclikepl->Add(h2);
        name = "BNueCC_" + plf + "_" + ccname + "/FD_TrueVsReco";
        h2 = (TH2D*)fcc->Get(name.c_str());
        farcclikepl->Add(h2);
        farcclikepl->Scale(nPOTFar/fpot_cc);
      }
      if(nfiles[ietag]==0) totalfarcclikepl[ietag] = (TH2D*)farcclikepl->Clone();
      else totalfarcclikepl[ietag]->Add(farcclikepl);
      
      name = "NuMuCC_" + pln + "_" + ccname + "/ND_TrueVsReco";
      if(pln==std)
      {
        nearcclikepl = (TH2D*)nearcclikestd->Clone(name.c_str());
      }
      else
      {
        nearcclikepl = (TH2D*)fcc->Get(name.c_str());
        name = "NC_" + pln + "_" + ccname + "/ND_TrueVsReco";
        h2 = (TH2D*)fcc->Get(name.c_str());
        nearcclikepl->Add(h2);
        name = "BNueCC_" + pln + "_" + ccname + "/ND_TrueVsReco";
        h2 = (TH2D*)fcc->Get(name.c_str());
        nearcclikepl->Add(h2);
        nearcclikepl->Scale(nPOTNear/npot_cc);
      }
      if(nfiles[ietag]==0) totalnearcclikepl[ietag] = (TH2D*)nearcclikepl->Clone();
      else totalnearcclikepl[ietag]->Add(nearcclikepl);
     
      name = "NuMuCC_" + mnf + "_Fid/FD_TrueVsReco";
      if(mnf==std)
      {
        farfidccmn = (TH2D*)farfidccstd->Clone(name.c_str());
      }
      else if(mnf==plf)
      {
        farfidccmn = (TH2D*)farfidccpl->Clone(name.c_str());
      }
      else
      {
        farfidccmn = (TH2D*)ffid->Get(name.c_str());
        farfidccmn->Scale(nPOTFar/fpot_fid);
      }
      if(nfiles[ietag]==0) totalfarfidccmn[ietag] = (TH2D*)farfidccmn->Clone();
      else totalfarfidccmn[ietag]->Add(farfidccmn);
      
      name = "NuMuCC_" + mnf + "_" + ccname + "/FD_TrueVsReco";
      if(mnf==std)
      {
        farcclikemn = (TH2D*)farcclikestd->Clone(name.c_str());
      }
      else if(mnf==plf)
      {
        farcclikemn = (TH2D*)farcclikepl->Clone(name.c_str());
      }
      else
      {
        farcclikemn = (TH2D*)fcc->Get(name.c_str());
        name = "NC_" + mnf + "_" + ccname + "/FD_TrueVsReco";
        h2 = (TH2D*)fcc->Get(name.c_str());
        farcclikemn->Add(h2);
        name = "BNueCC_" + mnf + "_" + ccname + "/FD_TrueVsReco";
        h2 = (TH2D*)fcc->Get(name.c_str());
        farcclikemn->Add(h2);
        farcclikemn->Scale(nPOTFar/fpot_cc);
      }
      if(nfiles[ietag]==0) totalfarcclikemn[ietag] = (TH2D*)farcclikemn->Clone();
      else totalfarcclikemn[ietag]->Add(farcclikemn);
      
      name = "NuMuCC_" + mnn + "_" + ccname + "/ND_TrueVsReco";
      if(mnn==std)
      {
        nearcclikemn = (TH2D*)nearcclikestd->Clone(name.c_str());
      }
      else if(mnn==pln)
      {
        nearcclikemn = (TH2D*)nearcclikepl->Clone(name.c_str());
      }
      else
      {
        nearcclikemn = (TH2D*)fcc->Get(name.c_str());
        name = "NC_" + mnn + "_" + ccname + "/ND_TrueVsReco";
        h2 = (TH2D*)fcc->Get(name.c_str());
        nearcclikemn->Add(h2);
        name = "BNueCC_" + mnn + "_" + ccname + "/ND_TrueVsReco";
        h2 = (TH2D*)fcc->Get(name.c_str());
        nearcclikemn->Add(h2);
        nearcclikemn->Scale(nPOTNear/npot_cc);
      }
      if(nfiles[ietag]==0) totalnearcclikemn[ietag] = (TH2D*)nearcclikemn->Clone();
      else totalnearcclikemn[ietag]->Add(nearcclikemn);
      
      nfiles[ietag]++;
    }

    //Make both types of file
    for (ietag=0; ietag<2; ietag++){
      if (nfiles[ietag]>0){
      for (ie=0; ie<Extrap.size(); ie++){
	if (ExtrapType[ie]==ietag) {
	  nt = Extrap[ie]->GetNTrue();
	  t = new double[nt+1];
	  for(i=0;i<nt+1;i++) t[i] = Extrap[ie]->Pred_3D[Extrapolate2D::qNuMuToNuTau]->GetZaxis()->GetBinLowEdge(i+1);
	  break;
	}
      }

      ccstd[ietag] = new TH1D(Form("ccstd%i",ietag),"",nt,t);
      for(it=0;it<nt;it++)
	{
	  sum=0;
	  for(ir=0;ir<totalfarcclikestd[ietag]->GetNbinsX();ir++)
	    {
	      if(totalnearcclikestd[ietag]->Integral(ir+1,ir+1,1,nt)>0)
		{
		  sum += (totaldata[ietag]->GetBinContent(ir+1)/totalnearcclikestd[ietag]->Integral(ir+1,ir+1,1,nt))*totalfarcclikestd[ietag]->GetBinContent(ir+1,it+1);
		}
	    }
	  temp=0;
	  if(totalfarcclikestd[ietag]->Integral(1,totalfarcclikestd[ietag]->GetNbinsX(),it+1,it+1)>0)
	    {
	      temp = totalfarfidccstd[ietag]->Integral(1,totalfarfidccstd[ietag]->GetNbinsX(),it+1,it+1)*(sum/totalfarcclikestd[ietag]->Integral(1,totalfarcclikestd[ietag]->GetNbinsX(),it+1,it+1));
	    }
	  ccstd[ietag]->SetBinContent(it+1,temp);
	}
      
      ccpl[ietag] = new TH1D(Form("ccpl%i",ietag),"",nt,t);
      for(it=0;it<nt;it++)
	{
	  sum=0;
	  for(ir=0;ir<totalfarcclikepl[ietag]->GetNbinsX();ir++)
	    {
	      if(totalnearcclikepl[ietag]->Integral(ir+1,ir+1,1,nt)>0)
		{
		  sum += (totaldata[ietag]->GetBinContent(ir+1)/totalnearcclikepl[ietag]->Integral(ir+1,ir+1,1,nt))*totalfarcclikepl[ietag]->GetBinContent(ir+1,it+1);
		}
	    }
	  temp=0;
	  if(totalfarcclikepl[ietag]->Integral(1,totalfarcclikepl[ietag]->GetNbinsX(),it+1,it+1)>0)
	    {
	      temp = totalfarfidccpl[ietag]->Integral(1,totalfarfidccpl[ietag]->GetNbinsX(),it+1,it+1)*(sum/totalfarcclikepl[ietag]->Integral(1,totalfarcclikepl[ietag]->GetNbinsX(),it+1,it+1));
	    }
	  ccpl[ietag]->SetBinContent(it+1,temp);
	}
      ccmn[ietag] = new TH1D(Form("ccmn%i",ietag),"",nt,t);
      for(it=0;it<nt;it++)
	{
	  sum=0;
	  for(ir=0;ir<totalfarcclikemn[ietag]->GetNbinsX();ir++)
	    {
	      if(totalnearcclikemn[ietag]->Integral(ir+1,ir+1,1,nt)>0)
		{
		  sum += (totaldata[ietag]->GetBinContent(ir+1)/totalnearcclikemn[ietag]->Integral(ir+1,ir+1,1,nt))*totalfarcclikemn[ietag]->GetBinContent(ir+1,it+1);
		}
	    }
	  temp=0;
	  if(totalfarcclikemn[ietag]->Integral(1,totalfarcclikemn[ietag]->GetNbinsX(),it+1,it+1)>0)
	    {
	      temp = totalfarfidccmn[ietag]->Integral(1,totalfarfidccmn[ietag]->GetNbinsX(),it+1,it+1)*(sum/totalfarcclikemn[ietag]->Integral(1,totalfarcclikemn[ietag]->GetNbinsX(),it+1,it+1));
	    }
	  ccmn[ietag]->SetBinContent(it+1,temp);
	}

      ccpl[ietag]->Add(ccstd[ietag],-1.);
      ccpl[ietag]->Divide(ccstd[ietag]);
      ccmn[ietag]->Add(ccstd[ietag],-1.);
      ccmn[ietag]->Divide(ccstd[ietag]);

      if(plf==mnf && pln==mnn)
	{
	  ccmn[ietag]->Scale(-1.);
	}
      }
    }

    for(ie=0;ie<Extrap.size();ie++)
      {
	ietag = ExtrapType[ie];

	Pred_CC_Fid_Plus1Sigma[systname].push_back((TH1D*)ccpl[ietag]->Clone(Form("%s_%i_Pred_CC_Fid_Plus1",systname.c_str(),ie+1)));
	Pred_CC_Fid_Minus1Sigma[systname].push_back((TH1D*)ccmn[ietag]->Clone(Form("%s_%i_Pred_CC_Fid_Minus1",systname.c_str(),ie+1)));
      }
    

      for(ie=0;ie<Extrap.size();ie++)
	{
	  nr = Extrap[ie]->GetNReco();
	  np = Extrap[ie]->GetNPID();
	  nt = Extrap[ie]->GetNTrue();

	  for(ip=0;ip<np;ip++)
	    {
	      for(ir=0;ir<nr;ir++)
		{
		  sum_plus_tau=0;
		  sum_minus_tau=0;
		  sum_plus_nue=0;
		  sum_minus_nue=0;
		  for(it=0;it<nt;it++)
		    {

		      sum_plus_tau += (Pred_CC_Fid_Plus1Sigma[systname][ie]->GetBinContent(it+1)*Extrap[ie]->Pred_3D[Extrapolate2D::qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1));
		      
		      sum_minus_tau += (Pred_CC_Fid_Minus1Sigma[systname][ie]->GetBinContent(it+1)*Extrap[ie]->Pred_3D[Extrapolate2D::qNuMuToNuTau]->GetBinContent(ip+1,ir+1,it+1));
		      
		      sum_plus_nue += (Pred_CC_Fid_Plus1Sigma[systname][ie]->GetBinContent(it+1)*Extrap[ie]->Pred_3D[Extrapolate2D::qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1));
		      
		      sum_minus_nue += (Pred_CC_Fid_Minus1Sigma[systname][ie]->GetBinContent(it+1)*Extrap[ie]->Pred_3D[Extrapolate2D::qNuMuToNue]->GetBinContent(ip+1,ir+1,it+1));
		    }

		  //make fractional errors again
		  sum_plus_tau = sum_plus_tau/Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip+1,ir+1);
		  sum_minus_tau = sum_minus_tau/Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip+1,ir+1);
		  sum_plus_nue = sum_plus_nue/Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1);
		  sum_minus_nue = sum_minus_nue/Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1);

		  NuTauCC_MC_Plus1Sigma[systname][ie]->SetBinContent(ip+1,ir+1,sum_plus_tau);
		  NuTauCC_MC_Minus1Sigma[systname][ie]->SetBinContent(ip+1,ir+1,sum_minus_tau);
		  NueCC_MC_Plus1Sigma[systname][ie]->SetBinContent(ip+1,ir+1,sum_plus_nue);
		  NueCC_MC_Minus1Sigma[systname][ie]->SetBinContent(ip+1,ir+1,sum_minus_nue);

		}
	    }
	}

  }
  
  if(flag==2)
    {

    //create empty histograms for numuCC shift
    for(ie=0;ie<Extrap.size();ie++)
    {
      nt = Extrap[ie]->GetNTrue();
      t = new double[nt+1];
      for(i=0;i<nt+1;i++)
	{
	  t[i] = Extrap[ie]->Pred_3D[Extrapolate2D::qNuMuToNuTau]->GetZaxis()->GetBinLowEdge(i+1);
	}
      Pred_CC_Fid_Plus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Plus1",systname.c_str(),ie+1),"",nt,t));
      Pred_CC_Fid_Minus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Minus1",systname.c_str(),ie+1),"",nt,t));

    }

    //RBTNOTE: DO THIS SEPARATELY FOR FHC AND RHC
    bool preExist = false;

    for (ie=0; ie<PreFiles.size(); ie++){
      if (!gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[ie].c_str()))) preExist = true;
    }

    if (!preExist)
    {
      for(ie=0;ie<Extrap.size();ie++)
      {

	nr = Extrap[ie]->GetNReco();
	np = Extrap[ie]->GetNPID();

        NueCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
        NueCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
        NuTauCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
        NuTauCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      }
      
      return;
    }

    for(j=0;j<app.size();j++)
    {

      nfiles[0]=0;
      nfiles[1]=0;

      for(ie=0;ie<Extrap.size();ie++)
      {
        if(gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[ie].c_str())))//if file doesn't exist
        {
          continue;
        }

	nr = Extrap[ie]->GetNReco();
	np = Extrap[ie]->GetNPID();
	nt = Extrap[ie]->GetNTrue();
	r = new double[nr+1];
	p = new double[np+1];
	for(i=0;i<nr+1;i++)
	  {
	    r[i] = Extrap[ie]->Pred_TotalBkgd->GetYaxis()->GetBinLowEdge(i+1);
	  }
	for(i=0;i<np+1;i++)
	  {
	    p[i] = Extrap[ie]->Pred_TotalBkgd->GetXaxis()->GetBinLowEdge(i+1);
	  }
        
	ietag=ExtrapType[ie];

	nPOTFar=Extrap[ie]->GetFarPOT();
	nPOTNear=Extrap[ie]->GetNearPOT();

        fpre = new TFile(gSystem->ExpandPathName(PreFiles[ie].c_str()),"READ");
        treepre = (TTree*)fpre->Get("paramtree");
        treepre->SetBranchAddress("nearPOT",&npot_pre);
        treepre->SetBranchAddress("farPOT",&fpot_pre);
        treepre->GetEntry(0);

        name = string(Background::AsString(app[j])) + "_" + std + "_Presel/FD_RecoVs" + Extrap[ie]->GetPID();

        prestd = (TH2D*)fpre->Get(name.c_str());
        prestd->Scale(nPOTFar/fpot_pre);
        MBH.Rebin2DHist(prestd,np,p,nr,r);
        if(nfiles[ietag]==0) totalprestd[ietag] = (TH2D*)prestd->Clone();
        else totalprestd[ietag]->Add(prestd);

        name = string(Background::AsString(app[j])) + "_" + plf + "_Presel/FD_RecoVs" + Extrap[ie]->GetPID();
        if(plf==std)
        {
          prepl = (TH2D*)prestd->Clone(name.c_str());
        }
        else
        {
          prepl = (TH2D*)fpre->Get(name.c_str());
          prepl->Scale(nPOTFar/fpot_pre);
          MBH.Rebin2DHist(prepl,np,p,nr,r);
        }
        if(nfiles[ietag]==0) totalprepl[ietag] = (TH2D*)prepl->Clone();
        else totalprepl[ietag]->Add(prepl);
        
        name = string(Background::AsString(app[j])) + "_" + mnf + "_Presel/FD_RecoVs" + Extrap[ie]->GetPID();
        if(mnf==std)
        {
          premn = (TH2D*)prestd->Clone(name.c_str());
        }
        else if(mnf==plf)
        {
          premn = (TH2D*)prepl->Clone(name.c_str());
        }
        else
        {
          premn = (TH2D*)fpre->Get(name.c_str());
          premn->Scale(nPOTFar/fpot_pre);
          MBH.Rebin2DHist(premn,np,p,nr,r);
        }
        if(nfiles[ietag]==0) totalpremn[ietag] = (TH2D*)premn->Clone();
        else totalpremn[ietag]->Add(premn);
        
        nfiles[ietag]++;
      }
      

      for (ietag=0; ietag<2; ietag++){
	if (nfiles[ietag]>0){
	totalprepl[ietag]->Add(totalprestd[ietag],-1.);
	totalprepl[ietag]->Divide(totalprestd[ietag]);
	totalpremn[ietag]->Add(totalprestd[ietag],-1.);
	totalpremn[ietag]->Divide(totalprestd[ietag]);
	
	if(plf==mnf && pln==mnn)
	  {
	    totalpremn[ietag]->Scale(-1.);
	  }
	}
      }

      for(ie=0;ie<Extrap.size();ie++)
      {
	ietag=ExtrapType[ie];

        if(app[j]==Background::kNueCC)
        {
          NueCC_MC_Plus1Sigma[systname].push_back((TH2D*)totalprepl[ietag]->Clone(Form("%s_%i_NueCC_MC_Plus1",systname.c_str(),ie+1)));
          NueCC_MC_Minus1Sigma[systname].push_back((TH2D*)totalpremn[ietag]->Clone(Form("%s_%i_NueCC_MC_Minus1",systname.c_str(),ie+1)));
        }
        else if(app[j]==Background::kNuTauCC)
        {
          NuTauCC_MC_Plus1Sigma[systname].push_back((TH2D*)totalprepl[ietag]->Clone(Form("%s_%i_NuTauCC_MC_Plus1",systname.c_str(),ie+1)));
          NuTauCC_MC_Minus1Sigma[systname].push_back((TH2D*)totalpremn[ietag]->Clone(Form("%s_%i_NuTauCC_MC_Minus1",systname.c_str(),ie+1)));
        }
      }
    }
  }
  
  return;
}
void ErrorCalc_Joint::ReadSpecialFiles(int n)
{


  string systname = SpecialSyst_SystName.at(n);
  TH2D *h;
  unsigned int ie;
  int i,nr,np,nt;
  double *r,*p,*t;
  
  //create empty histograms for all shifts
  for(ie=0;ie<Extrap.size();ie++)
  {

    nr = Extrap[ie]->GetNReco();
    np = Extrap[ie]->GetNPID();
    nt = Extrap[ie]->GetNTrue();
    r = new double[nr+1];
    p = new double[np+1];
    t = new double[nt+1];
    for(i=0;i<nr+1;i++)
      {
	r[i] = Extrap[ie]->Pred_TotalBkgd->GetYaxis()->GetBinLowEdge(i+1);
      }
    for(i=0;i<np+1;i++)
      {
	p[i] = Extrap[ie]->Pred_TotalBkgd->GetXaxis()->GetBinLowEdge(i+1);
      }
    for(i=0;i<nt+1;i++)
      {
	t[i] = Extrap[ie]->Pred_3D[Extrapolate2D::qNuMuToNuTau]->GetZaxis()->GetBinLowEdge(i+1);
      }

    FN_NC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FN_NC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
    FN_NC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FN_NC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    
    FN_NuMuCC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FN_NuMuCC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
    FN_NuMuCC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FN_NuMuCC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    
    FN_BNueCC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FN_BNueCC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
    FN_BNueCC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FN_BNueCC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    
    ND_NC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_ND_NC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
    ND_NC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_ND_NC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    
    ND_NuMuCC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_ND_NuMuCC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
    ND_NuMuCC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_ND_NuMuCC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    
    ND_BNueCC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_ND_BNueCC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
    ND_BNueCC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_ND_BNueCC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    
    FD_NC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FD_NC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
    FD_NC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FD_NC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    
    FD_NuMuCC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FD_NuMuCC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
    FD_NuMuCC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FD_NuMuCC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    
    FD_BNueCC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FD_BNueCC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
    FD_BNueCC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_FD_BNueCC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    
    NueCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
    NueCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    NuTauCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
    NuTauCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    
    Pred_CC_Fid_Plus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Plus1",systname.c_str(),ie+1),"",nt,t));
    Pred_CC_Fid_Minus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Minus1",systname.c_str(),ie+1),"",nt,t));
  }

  int flag = SpecialSyst_Flag.at(n);

  if(SpecialSyst_Value.at(n)>0)//if value is set, use that number in all bins
  {

    for(unsigned int ie=0;ie<Extrap.size();ie++)
      {
	
	h = (TH2D*)Extrap[ie]->Pred_TotalBkgd->Clone("htemp");
	h->Reset();
	for(int ip=0;ip<Extrap[ie]->GetNPID();ip++)
	  {
	    for(int ir=0;ir<Extrap[ie]->GetNReco();ir++)
	      {
		h->SetBinContent(ip+1,ir+1,SpecialSyst_Value.at(n));
	      }
	  }
	
      if(flag==0 || flag==1)//for NC,NuMuCC,BNueCC
      {
        FN_NC_Plus1Sigma[systname][ie]->Add(h);
        FN_NC_Minus1Sigma[systname][ie]->Add(h,-1);
        
        FN_NuMuCC_Plus1Sigma[systname][ie]->Add(h);
        FN_NuMuCC_Minus1Sigma[systname][ie]->Add(h,-1);
        
        FN_BNueCC_Plus1Sigma[systname][ie]->Add(h);
        FN_BNueCC_Minus1Sigma[systname][ie]->Add(h,-1);
        
        FD_BNueCC_Plus1Sigma[systname][ie]->Add(h);
        FD_BNueCC_Minus1Sigma[systname][ie]->Add(h,-1);
      }
      if(flag==0 || flag==2)
      {
        NuTauCC_MC_Plus1Sigma[systname][ie]->Add(h);
        NuTauCC_MC_Minus1Sigma[systname][ie]->Add(h,-1);
      }
      if(flag==0 || flag==3)
      {
        NueCC_MC_Plus1Sigma[systname][ie]->Add(h);
        NueCC_MC_Minus1Sigma[systname][ie]->Add(h,-1);
      }
      if(flag==4)//for NC and NuMuCC (fake HOO)
      {
        FN_NC_Plus1Sigma[systname][ie]->Add(h);
        FN_NC_Minus1Sigma[systname][ie]->Add(h,-1);
        
        FN_NuMuCC_Plus1Sigma[systname][ie]->Add(h);
        FN_NuMuCC_Minus1Sigma[systname][ie]->Add(h,-1);

      }


    }
    
    return;
  }
  
  //if value is 0, read from a file
  string filetag,histname;
  filetag=SpecialSyst_FileTag.at(n);
  histname = SpecialSyst_HistName.at(n);
  
  TFile *f;
  
  for(unsigned int ie=0;ie<Extrap.size();ie++)
  {
    if(ExtrapType[ie]==0&&gSystem->AccessPathName(gSystem->ExpandPathName(Form("%s/%s_FHC%i.root",SysFileDir.c_str(),filetag.c_str(),Extrap[ie]->GetRunPeriod()))))
    {
      continue;
    }
    if(ExtrapType[ie]==1&&gSystem->AccessPathName(gSystem->ExpandPathName(Form("%s/%s_RHC%i.root",SysFileDir.c_str(),filetag.c_str(),Extrap[ie]->GetRunPeriod()))))
    {
      continue;
    }
    
    if (ExtrapType[ie]==0) f = new TFile(gSystem->ExpandPathName(Form("%s/%s_FHC%i.root",SysFileDir.c_str(),filetag.c_str(),Extrap[ie]->GetRunPeriod())),"READ");
    if (ExtrapType[ie]==1) f = new TFile(gSystem->ExpandPathName(Form("%s/%s_RHC%i.root",SysFileDir.c_str(),filetag.c_str(),Extrap[ie]->GetRunPeriod())),"READ");
  
    h = (TH2D*)f->Get(histname.c_str());
    
    if(flag==0 || flag==1)//for NC,NuMuCC,BNueCC
    {
      FN_NC_Plus1Sigma[systname][ie]->Add(h);
      FN_NC_Minus1Sigma[systname][ie]->Add(h,-1);
        
      FN_NuMuCC_Plus1Sigma[systname][ie]->Add(h);
      FN_NuMuCC_Minus1Sigma[systname][ie]->Add(h,-1);
        
      FN_BNueCC_Plus1Sigma[systname][ie]->Add(h);
      FN_BNueCC_Minus1Sigma[systname][ie]->Add(h,-1);
        
      FD_BNueCC_Plus1Sigma[systname][ie]->Add(h);
      FD_BNueCC_Minus1Sigma[systname][ie]->Add(h,-1);
    }
    if(flag==0 || flag==2)
    {
      NuTauCC_MC_Plus1Sigma[systname][ie]->Add(h);
      NuTauCC_MC_Minus1Sigma[systname][ie]->Add(h,-1);
    }
    if(flag==0 || flag==3)
    {
      NueCC_MC_Plus1Sigma[systname][ie]->Add(h);
      NueCC_MC_Minus1Sigma[systname][ie]->Add(h,-1);
    }
    if(flag==4)//for NC and NuMuCC (fake HOO)
      {
	FN_NC_Plus1Sigma[systname][ie]->Add(h);
	FN_NC_Minus1Sigma[systname][ie]->Add(h,-1);

	FN_NuMuCC_Plus1Sigma[systname][ie]->Add(h);
	FN_NuMuCC_Minus1Sigma[systname][ie]->Add(h,-1);
      }

  }
  
  return;
}

void ErrorCalc_Joint::AddCovarianceMatrix(string systname,string file,string histname,int flag, int isRHC)
{
  if(flag>7)
  {
    cout<<"Error in AddCovarianceMatrix(): 'flag' should be <=7.  Not adding this covariance matrix."<<endl;
    return;
  }
  
  ExtraCovariance_SystName.push_back(systname);
  ExtraCovariance_File.push_back(file);
  ExtraCovariance_HistName.push_back(histname);
  ExtraCovariance_Flag.push_back(flag);
  ExtraCovariance_IsRHC.push_back(isRHC);
  
  return;
}


void ErrorCalc_Joint::ReadCovarianceFiles(int n)
{

  int matrixtype = ExtraCovariance_IsRHC.at(n);

  string systname = ExtraCovariance_SystName.at(n);
  string file = ExtraCovariance_File.at(n);
  string hist = ExtraCovariance_HistName.at(n);
//   int flag = ExtraCovariance_Flag.at(n);
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(file.c_str())))
  {
    cout<<"Can't read "<<file<<endl;
    return;
  }
  
  TFile *f = new TFile(file.c_str(),"READ");
  TH2D *h = (TH2D*)f->Get(hist.c_str());

  unsigned int ie=0;
  int nbins=0;

  //Bins for different parts of joint fit:
  for (ie=0; ie<Extrap.size(); ie++){
    if (ExtrapType[ie]==0&&matrixtype==0) nbins = Extrap[ie]->GetNPID()*Extrap[ie]->GetNReco();
    if (ExtrapType[ie]==1&&matrixtype==1) nbins = Extrap[ie]->GetNPID()*Extrap[ie]->GetNReco();
  }
  
  if(h->GetNbinsX() != nbins)
  {
    cout<<"Error in ReadCovarianceFiles(): # bins from matrix in file doesn't match extrapolation bins!"<<endl;
    return;
  }
  
  ExtraCovariance[systname] = h;
  ExtraCovarianceTag[systname] = matrixtype;
  
  return;
}

void ErrorCalc_Joint::CalculateSystErrorMatrixExtrap()
{
  //calculates the SYSTEMATIC part of the error matrix
  if(Extrap.size()==0)
  {
    cout<<"Error in ErrorCalc_Joint::CalculateSystErrorMatrix(): You must add an Extrapolate2D object.  Quitting..."<<endl;
    return;
  }
    
  Initialize();

  if(!Init) return;
  
  InitializeSys();

  //cout << "Calc matrix" << endl;

  string systname;
  
  unsigned int ie=0;

  int nPIDfhc=0;
  int nPIDrhc=0;
  int nRecofhc=0;
  int nRecorhc=0;

  for (ie=0; ie<Extrap.size(); ie++){
    if (ExtrapType[ie]==0) {
      nPIDfhc = Extrap[ie]->GetNPID();
      nRecofhc = Extrap[ie]->GetNReco();
    }
    if (ExtrapType[ie]==1) {
      nPIDrhc = Extrap[ie]->GetNPID();
      nRecorhc = Extrap[ie]->GetNReco();
    }
  }

  //int nPID = Extrap[0]->GetNPID();
  //int nbins = Extrap[0]->GetNReco()*Extrap[0]->GetNPID();
                                                                                                                   //int nt = Extrap[0]->GetNTrue();
  
  std::map<string, vector<TH2D*> >::iterator systiter = FN_NC_Plus1Sigma.begin();
  std::map<string, vector<TH2D*> >::iterator last  = FN_NC_Plus1Sigma.end();
  
  int i,j;
  int ip,ir,jp,jr;
  double element;
  double di,dj;
  bool addi, addj;

  CovMatrix->Reset();
  
  int nbinsFHC = nPIDfhc*nRecofhc;
  int nbinsRHC = nPIDrhc*nRecorhc;

  int nbins = nbinsFHC + nbinsRHC;

  for(i=0;i<nbins;i++)
  {

    if (i<nbinsFHC){
      ir = int(i/nPIDfhc);
      ip = i%nPIDfhc;
    }      
    if (i>=nbinsFHC){
      ir = int((i-nbinsFHC)/nPIDrhc);
      ip = (i-nbinsFHC)%nPIDrhc;
    }      


    for(j=0;j<nbins;j++)
      {

	if (j<nbinsFHC){
	  jr = int(j/nPIDfhc);
	  jp = j%nPIDfhc;
	}      
	if (j>=nbinsFHC){
	  jr = int((j-nbinsFHC)/nPIDrhc);
	  jp = (j-nbinsFHC)%nPIDrhc;
	}
	
	element = 0;
	systiter = FN_NC_Plus1Sigma.begin();
	while(systiter!=last)
	  {

	    systname = systiter->first;
	    
	    di=0;
	    dj=0;

	    for(ie=0;ie<Extrap.size();ie++)
	      {
		addi=false;
		addj=false;

		//Only add the element if bin is correct for that particular extrapolation
		if ((ExtrapType[ie]==0&&i<nbinsFHC)||(ExtrapType[ie]==1&&i>=nbinsFHC)) addi = true;
		if ((ExtrapType[ie]==0&&j<nbinsFHC)||(ExtrapType[ie]==1&&j>=nbinsFHC)) addj = true;


		if (addi) di += FN_NC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kNC]->GetBinContent(ip+1,ir+1) + FN_NuMuCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kNuMuCC]->GetBinContent(ip+1,ir+1) + NuTauCC_MC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip+1,ir+1) + NueCC_MC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1);
		
		if (addj) dj += FN_NC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*Extrap[ie]->Pred[Background::kNC]->GetBinContent(jp+1,jr+1) + FN_NuMuCC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*Extrap[ie]->Pred[Background::kNuMuCC]->GetBinContent(jp+1,jr+1) + NuTauCC_MC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(jp+1,jr+1) + NueCC_MC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(jp+1,jr+1);
		  
		if(Extrap[ie]->GetFNforBeamNue())
		  {
		    if (addi) di += FN_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
		    if (addj) dj += FN_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(jp+1,jr+1);
		  }
		else
		  {
		    if (addi) di += FD_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
		    if (addj) dj += FD_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(jp+1,jr+1);
		  }
		
	      } 
	    
	    element += (di*dj);
	    
	    systiter++;
	  }
	

	//Add an additional element (1% of stat error) to force matrix to be non-singular
	double totstat=0;

	for(ie=0;ie<Extrap.size();ie++)
	  {
	    addi=false;
	    addj=false;
	    
	    //Only add the element if bin is correct for that particular extrapolation
	    if ((ExtrapType[ie]==0&&i<nbinsFHC)||(ExtrapType[ie]==1&&i>=nbinsFHC)) addi = true;
	    if ((ExtrapType[ie]==0&&j<nbinsFHC)||(ExtrapType[ie]==1&&j>=nbinsFHC)) addj = true;
	    if ((i==j)&&addi&&addj) totstat += Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1) 
				      + Extrap[ie]->Pred_TotalBkgd->GetBinContent(ip+1,ir+1);
	  }
	
	element += (totstat*0.01);
	
	CovMatrix->SetBinContent(i+1,j+1,element);

      }


  }
  
  std::map<string,TH2D*>::iterator coviter = ExtraCovariance.begin();
  std::map<string,TH2D*>::iterator covlast  = ExtraCovariance.end();
  double ni,nj;
  int k=0;
  int flag;
  int ietype=-1;

  while(coviter!=covlast)
  {
    systname = coviter->first;
    flag = ExtraCovariance_Flag.at(k);

    ietype=ExtraCovarianceTag[systname];

    for(i=0;i<nbins;i++)
    {
      
      if (i<nbinsFHC){
	ir = int(i/nPIDfhc);
	ip = i%nPIDfhc;
      }      
      if (i>=nbinsFHC){
	ir = int((i-nbinsFHC)/nPIDrhc);
	ip = (i-nbinsFHC)%nPIDrhc;
      }      

      ni=0;
      for(ie=0;ie<Extrap.size();ie++)
      {

	//Only add the element if bin is correct for that particular extrapolation
	if ((ietype==0&&ExtrapType[ie]==0&&i<nbinsFHC)
	    ||(ietype==1&&ExtrapType[ie]==1&&i>=nbinsFHC)) {
	  
	  if(flag==0 || flag==1 || flag==4 || flag==5)
	    {
	      ni+=Extrap[ie]->Pred[Background::kNC]->GetBinContent(ip+1,ir+1);
	    }
	  if(flag==0 || flag==1 || flag==4 || flag==6)
	    {
	      ni+=Extrap[ie]->Pred[Background::kNuMuCC]->GetBinContent(ip+1,ir+1);
	    }
	  if(flag==0 || flag==1 || flag==7)
	    {
	      ni+=Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
	    }
	  if(flag==0 || flag==2)
	    {
	      ni+=Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip+1,ir+1);
	    }
	  if(flag==0 || flag==3)
	    {
	      ni+=Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1);
	    }
	}
      }
      
      for(j=0;j<nbins;j++)
      {
        
	if (j<nbinsFHC){
	  jr = int(j/nPIDfhc);
	  jp = j%nPIDfhc;
	}      
	if (j>=nbinsFHC){
	  jr = int((j-nbinsFHC)/nPIDrhc);
	  jp = (j-nbinsFHC)%nPIDrhc;
	}      

        nj=0;
	  
	for(ie=0;ie<Extrap.size();ie++)
	  {
	    
	    //Only add the element if bin is correct for that particular extrapolation
	    if ((ietype==0&&ExtrapType[ie]==0&&j<nbinsFHC)
		||(ietype==1&&ExtrapType[ie]==1&&j>=nbinsFHC)) {

	      if(flag==0 || flag==1 || flag==4 || flag==5)
		{
		  nj+=Extrap[ie]->Pred[Background::kNC]->GetBinContent(jp+1,jr+1);
		}
	      if(flag==0 || flag==1 || flag==4 || flag==6)
		{
		  nj+=Extrap[ie]->Pred[Background::kNuMuCC]->GetBinContent(jp+1,jr+1);
		}
	      if(flag==0 || flag==1 || flag==7)
		{
		  nj+=Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(jp+1,jr+1);
		}
	      if(flag==0 || flag==2)
		{
		  nj+=Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(jp+1,jr+1);
		}
	      if(flag==0 || flag==3)
		{
		  nj+=Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(jp+1,jr+1);
		}
	    }
	}

        element = CovMatrix->GetBinContent(i+1,j+1);
        element += (ni*nj*ExtraCovariance[systname]->GetBinContent(i+1,j+1));
        CovMatrix->SetBinContent(i+1,j+1,element);
      }
    }
    
    coviter++;
    k++;
  }
 
  return;
}

void ErrorCalc_Joint::CalculateSystErrorMatrixGrid()
{
  //calculates the SYSTEMATIC part of the error matrix
  if(Extrap.size()==0)
  {
    cout<<"Error in ErrorCalc::CalculateSystErrorMatrixGrid(): You must add an Extrapolate2D object.  Quitting..."<<endl;
    return;
  }

  Initialize();
  if(!Init) return;
  
  InitializeSys();

  string systname;
  
  int nPIDfhc=0;
  int nPIDrhc=0;
  int nRecofhc=0;
  int nRecorhc=0;
  unsigned int ie=0;

  for (ie=0; ie<Extrap.size(); ie++){

    if (ExtrapType[ie]==0) {
      nPIDfhc = Extrap[ie]->GetNPID();
      nRecofhc = Extrap[ie]->GetNReco();
    }
    if (ExtrapType[ie]==1) {
      nPIDrhc = Extrap[ie]->GetNPID();
      nRecorhc = Extrap[ie]->GetNReco();
    }
  }

  int nbinsFHC = nPIDfhc*nRecofhc;
  int nbinsRHC = nPIDrhc*nRecorhc;
  
  int nbins = nbinsFHC + nbinsRHC;

  std::map<string, vector<TH2D*> >::iterator systiter = FN_NC_Plus1Sigma.begin();
  std::map<string, vector<TH2D*> >::iterator last  = FN_NC_Plus1Sigma.end();
  

  int i,j;
  int ii,jj;
  int ip,ir,jp,jr;
  double element;
  double di,dj;
  bool addi, addj;
  
//  unsigned int ietmp=0;
//  unsigned int ieFHC=0;
//  unsigned int ieRHC=0;
  CovMatrix->Reset();
  
  for(ii=0;ii<nbins;ii++)
  {

    if (ii<nbinsFHC){
      i = ii;
      ir = int(ii/nPIDfhc);
      ip = ii%nPIDfhc;
    }      
    if (ii>=nbinsFHC){
      i = ii-nbinsFHC;
      ir = int((ii-nbinsFHC)/nPIDrhc);
      ip = (ii-nbinsFHC)%nPIDrhc;
    }      

    
    for(jj=0;jj<nbins;jj++)
    {
	if (jj<nbinsFHC){
	  j = jj;
	  jr = int(jj/nPIDfhc);
	  jp = jj%nPIDfhc;
	}      
	if (jj>=nbinsFHC){
	  j = jj-nbinsFHC;
	  jr = int((jj-nbinsFHC)/nPIDrhc);
	  jp = (jj-nbinsFHC)%nPIDrhc;
	}
      
      element = 0;
      systiter = FN_NC_Plus1Sigma.begin();
      while(systiter!=last)
      {
        systname = systiter->first;
        
        di=0;
        dj=0;
        for(ie=0;ie<Extrap.size();ie++)
        {

	  addi=false;
	  addj=false;
	  
	  //Only add the element if bin is correct for that particular extrapolation
	  if ((ExtrapType[ie]==0&&ii<nbinsFHC)||(ExtrapType[ie]==1&&ii>=nbinsFHC)) addi = true;
	  if ((ExtrapType[ie]==0&&jj<nbinsFHC)||(ExtrapType[ie]==1&&jj>=nbinsFHC)) addj = true;


          if (addi) di += FN_NC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*GridPred[Background::kNC][ie]->GetBinContent(i+1) + FN_NuMuCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*GridPred[Background::kNuMuCC][ie]->GetBinContent(i+1) + NuTauCC_MC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*GridPred[Background::kNuTauCC][ie]->GetBinContent(i+1) + NueCC_MC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*GridPred[Background::kNueCC][ie]->GetBinContent(i+1);
          
          if (addj) dj += FN_NC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*GridPred[Background::kNC][ie]->GetBinContent(j+1) + FN_NuMuCC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*GridPred[Background::kNuMuCC][ie]->GetBinContent(j+1) + NuTauCC_MC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*GridPred[Background::kNuTauCC][ie]->GetBinContent(j+1) + NueCC_MC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*GridPred[Background::kNueCC][ie]->GetBinContent(j+1);
          
          if(Extrap[ie]->GetFNforBeamNue())
          {
            if (addi) di += FN_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*GridPred[Background::kBNueCC][ie]->GetBinContent(i+1);
            
            if (addj) dj += FN_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*GridPred[Background::kBNueCC][ie]->GetBinContent(j+1);
            
          }
          else
          {
            if (addi) di += FD_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*GridPred[Background::kBNueCC][ie]->GetBinContent(i+1);
            
            if (addj) dj += FD_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*GridPred[Background::kBNueCC][ie]->GetBinContent(j+1);
            
          }
        }
        
        element += (di*dj);
        
        systiter++;
      }

      //Add an additional element (1% of stat error) to force matrix to be non-singular
	double totstat=0;

	for(ie=0;ie<Extrap.size();ie++)
	  {
	    addi=false;
	    addj=false;
	    
	    //Only add the element if bin is correct for that particular extrapolation
	    if ((ExtrapType[ie]==0&&ii<nbinsFHC)||(ExtrapType[ie]==1&&ii>=nbinsFHC)) addi = true;
	    if ((ExtrapType[ie]==0&&jj<nbinsFHC)||(ExtrapType[ie]==1&&jj>=nbinsFHC)) addj = true;
	    if ((ii==jj)&&addi&&addj) totstat += GridPred[Background::kNC][ie]->GetBinContent(i+1) + GridPred[Background::kNuMuCC][ie]->GetBinContent(i+1) + GridPred[Background::kBNueCC][ie]->GetBinContent(i+1) + GridPred[Background::kNuTauCC][ie]->GetBinContent(i+1) + GridPred[Background::kNueCC][ie]->GetBinContent(i+1);
				      
	  }
	
	element += (totstat*0.01);
	
	CovMatrix->SetBinContent(ii+1,jj+1,element);
    }
  }
  
  std::map<string,TH2D*>::iterator coviter = ExtraCovariance.begin();
  std::map<string,TH2D*>::iterator covlast  = ExtraCovariance.end();

  double ni,nj;
  int k=0;
  int flag;
  int ietype=-1;

  k=0;
  while(coviter!=covlast)
  {
    systname = coviter->first;
    flag = ExtraCovariance_Flag.at(k);

    ietype=ExtraCovarianceTag[systname];

    for(ii=0;ii<nbins;ii++)
    {
      
      if (ii<nbinsFHC)  i = ii;
      if (ii>=nbinsFHC) i = ii-nbinsFHC;
      
      ni=0;
      for(ie=0;ie<Extrap.size();ie++)
      {
	//Only add the element if bin is correct for that particular extrapolation
	if ((ietype==0&&ExtrapType[ie]==0&&ii<nbinsFHC)
	    ||(ietype==1&&ExtrapType[ie]==1&&ii>=nbinsFHC)) {
	  
	  if(flag==0 || flag==1 || flag==4 || flag==5)
	    {
	      ni+=GridPred[Background::kNC][ie]->GetBinContent(i+1);
	    }
	  if(flag==0 || flag==1 || flag==4 || flag==6)
	    {
	      ni+=GridPred[Background::kNuMuCC][ie]->GetBinContent(i+1);
	    }
	  if(flag==0 || flag==1 || flag==7)
	    {
	      ni+=GridPred[Background::kBNueCC][ie]->GetBinContent(i+1);
	    }
	  if(flag==0 || flag==2)
	    {
	      ni+=GridPred[Background::kNuTauCC][ie]->GetBinContent(i+1);
	    }
	  if(flag==0 || flag==3)
	    {
	      ni+=GridPred[Background::kNueCC][ie]->GetBinContent(i+1);
	    }
	}
      }
      
      for(jj=0;jj<nbins;jj++)
	{
	  
	if (jj<nbinsFHC)  j = jj;
	if (jj>=nbinsFHC) j = jj-nbinsFHC;
	
        nj=0;
        for(ie=0;ie<Extrap.size();ie++)
	  {
	    //Only add the element if bin is correct for that particular extrapolation
	    if ((ietype==0&&ExtrapType[ie]==0&&jj<nbinsFHC)
		||(ietype==1&&ExtrapType[ie]==1&&jj>=nbinsFHC)) {

	      if(flag==0 || flag==1 || flag==4 || flag==5)
		{
		  nj+=GridPred[Background::kNC][ie]->GetBinContent(j+1);
		}
	      if(flag==0 || flag==1 || flag==4 || flag==6)
		{
		  nj+=GridPred[Background::kNuMuCC][ie]->GetBinContent(j+1);
		}
	      if(flag==0 || flag==1 || flag==7)
		{
		  nj+=GridPred[Background::kBNueCC][ie]->GetBinContent(j+1);
		}
	      if(flag==0 || flag==2)
		{
		  nj+=GridPred[Background::kNuTauCC][ie]->GetBinContent(j+1);
		}
	      if(flag==0 || flag==3)
		{
		  nj+=GridPred[Background::kNueCC][ie]->GetBinContent(j+1);
		}
	    }
	  }

	//Added protection:
	if ((ietype==0&ii<nbinsFHC&&jj<nbinsFHC)
	    ||(ietype==1&&ii>=nbinsFHC&&jj>=nbinsFHC)){
	  element = CovMatrix->GetBinContent(ii+1,jj+1);
	  element += (ni*nj*ExtraCovariance[systname]->GetBinContent(i+1,j+1));
	  CovMatrix->SetBinContent(ii+1,jj+1,element);
	}


	}
    }
    
    coviter++;
    k++;
  }
  
  return;
}


 void ErrorCalc_Joint::CalculateHOOError()
{
  //For the joint fit, we're just going to initialize this to 0.
  
  Initialize();
  if(!Init) return;
  
  CovMatrix_Decomp->Reset();

  
  return;
}

void ErrorCalc_Joint::SetGridPred(int nbins, vector< vector<double> > nc, vector< vector<double> > cc, vector< vector<double> > bnue, vector< vector<double> > tau, vector< vector<double> > sig)
{
  if(Extrap.size()==0)
  {
    cout<<"Error in ErrorCalc::SetGridPred(): You must add an Extrapolate2D object.  Quitting..."<<endl;
    return;
  }

  unsigned int ie;
  int nbinsFHC=0;
  int nbinsRHC=0;

  for (ie=0;ie<Extrap.size();ie++){
    if (ExtrapType[ie]==0) nbinsFHC=GridPred[Background::kNC][ie]->GetNbinsX();
    if (ExtrapType[ie]==1) nbinsRHC=GridPred[Background::kNC][ie]->GetNbinsX();
  }

  if (nbins != (nbinsFHC + nbinsRHC)) cout << "Warning!  SetGridPred mismatch!" << endl;

  Initialize();
  if(!Init) return;
  
  int i;
  int ii;  

  unsigned int ieFHC=0;
  unsigned int ieRHC=0;
  unsigned int ietmp=0;

  for(ie=0;ie<Extrap.size();ie++)
  {
    GridPred[Background::kNC][ie]->Reset();
    GridPred[Background::kNuMuCC][ie]->Reset();
    GridPred[Background::kBNueCC][ie]->Reset();
    GridPred[Background::kNuTauCC][ie]->Reset();
    GridPred[Background::kNueCC][ie]->Reset();
   
    for(i=0;i<nbins;i++)
    {
      if (ExtrapType[ie]==0) { ii=i; ietmp = ieFHC; }
      if (ExtrapType[ie]==1) { ii=i-nbinsFHC; ietmp = ieRHC; }

      if ((ExtrapType[ie]==0&&i<nbinsFHC)||
	  (ExtrapType[ie]==1&&i>=nbinsFHC)){

	GridPred[Background::kNC][ie]->SetBinContent(ii+1,nc[i][ietmp]);
	GridPred[Background::kNuMuCC][ie]->SetBinContent(ii+1,cc[i][ietmp]);
	GridPred[Background::kBNueCC][ie]->SetBinContent(ii+1,bnue[i][ietmp]);
	GridPred[Background::kNuTauCC][ie]->SetBinContent(ii+1,tau[i][ietmp]);
	GridPred[Background::kNueCC][ie]->SetBinContent(ii+1,sig[i][ietmp]);
      }
    }

    if (ExtrapType[ie]==0) ieFHC++;
    if (ExtrapType[ie]==1) ieRHC++;

  }
  
  return;
}
