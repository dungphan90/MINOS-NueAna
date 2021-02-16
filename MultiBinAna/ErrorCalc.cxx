#define ErrorCalc_C

#include "NueAna/MultiBinAna/ErrorCalc.h"

ErrorCalc::ErrorCalc()
{
  Init = false;
  InitSys = false;
  InitDecompSys = false;
  UseGrid = false;

  reg = 0.5;
    
  return;
}
ErrorCalc::~ErrorCalc()
{
}
void ErrorCalc::AddExtrap(Extrapolate2D *E)
{
  Extrap.push_back(E);
  return;
}
void ErrorCalc::AddFNExtrapSyst(string systname,string farplustag, string nearplustag, string farminustag, string nearminustag,string stdtag,string filetag,int flag,bool useAllRuns)
{
  if(flag>2)
  {
    cout<<"Error in AddFNExtrapSyst(): 'flag' should be <=2.  Not adding this systematic."<<endl;
    return;
  }
  FNExtrap_SystName.push_back(systname);
  FNExtrap_FarPlusTag.push_back(farplustag);
  FNExtrap_NearPlusTag.push_back(nearplustag);
  FNExtrap_FarMinusTag.push_back(farminustag);
  FNExtrap_NearMinusTag.push_back(nearminustag);
  FNExtrap_StdTag.push_back(stdtag);
  FNExtrap_FileTag.push_back(filetag);
  FNExtrap_Flag.push_back(flag);
  FNExtrap_UseAllRuns.push_back(useAllRuns);
  
  return;
}
void ErrorCalc::AddSpecialFNSyst(string systname,double err,string histname,string filetag,int flag)
{
  if(flag>4)
  {
    cout<<"Error in AddSpecialFNSyst(): 'flag' should be <=4.  Not adding this systematic."<<endl;
    return;
  }
  SpecialSyst_SystName.push_back(systname);
  SpecialSyst_FileTag.push_back(filetag);
  SpecialSyst_HistName.push_back(histname);
  SpecialSyst_Value.push_back(err);
  SpecialSyst_Flag.push_back(flag);
  
  return;
}
void ErrorCalc::AddCovarianceMatrix(string systname,string file,string histname,int flag)
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
  
  return;
}
void ErrorCalc::Initialize()
{
  if(Init) return;
  
  cout.precision(5);
  cout<<fixed;
  
  unsigned int ie;
  for(ie=0;ie<Extrap.size();ie++)
  {
    Extrap[ie]->GetPrediction();
    if(Extrap[ie]->GetOscFlag())
    {
      cout<<"Error in ErrorCalc::Initialize(): You have called Extrapolate2D::OscillatePrediction() before initializing the systematics!  Quitting..."<<endl;
      return;
    }
  }
  
  gROOT->cd();
  int nbins = Extrap[0]->GetNReco()*Extrap[0]->GetNPID();
  CovMatrix = new TH2D("CovMatrix","",nbins,-0.5,nbins-0.5,nbins,-0.5,nbins-0.5);
  CovMatrix_Decomp = new TH2D("CovMatrix_Decomp","",nbins,-0.5,nbins-0.5,nbins,-0.5,nbins-0.5);
  
  string st;
  for(ie=0;ie<Extrap.size();ie++)
  {
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
void ErrorCalc::InitializeSys()
{
  if(InitSys) return;

  std::cout << " ErrorCalc::InitializeSys " << std::endl; 
  for(unsigned int i=0;i<FNExtrap_SystName.size();i++)
  {

    if(!FNExtrap_UseAllRuns[i])
    {
      std::cout << " ErrorCalc::InitializeSys ReadSysFiles_FNExtrap("<<i<<") " << std::endl; 
      ReadSysFiles_FNExtrap(i);
      std::cout << " ErrorCalc::InitializeSys ReadSysFiles_Appearance("<<i<<")" << std::endl; 
      ReadSysFiles_Appearance(i);
    }
    else
    {
      std::cout << " ErrorCalc::InitializeSys ReadSysFiles_FNExtrap_AllRuns("<<i<<")" << std::endl; 
      ReadSysFiles_FNExtrap_AllRuns(i);
      std::cout << " ErrorCalc::InitializeSys ReadSysFiles_Appearance_AllRuns("<<i<<")" << std::endl; 
      ReadSysFiles_Appearance_AllRuns(i);
    }
  }
  for(unsigned int i=0;i<SpecialSyst_SystName.size();i++)
  {
    std::cout << " ErrorCalc::InitializeSys ReadSpecialFiles("<<i<<")" << std::endl; 
    ReadSpecialFiles(i);
  }
  for(unsigned int i=0;i<ExtraCovariance_SystName.size();i++)
  {
    std::cout << " ErrorCalc::InitializeSys ReadCovarianceFiles("<<i<<")" << std::endl; 
    ReadCovarianceFiles(i);
  }
  
  InitSys = true;
  
  return;
}
void ErrorCalc::InitializeDecompSys()
{
  if(InitDecompSys) return;
  
  NDCovMatrix=0;
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(NDData_infile.c_str())))
  {
    cout<<"Error in InitializeDecompSys(): Failure to read ND Data file."<<endl;
    return;
  }
  
  if(Extrap.size()==0) 
  {
    cout<<"Error in InitializeDecompSys(): No Extrapolate2D objects have been added."<<endl;
    return;
  }
  
  unsigned int npid = Extrap[0]->GetNPID();//number pid bins
  unsigned int nenergy = Extrap[0]->GetNReco();//number of reco energy bins
  int ncomp = 3;//number of background components
  unsigned int nrun = Extrap.size();//number of runs
  int ntot=npid*nenergy*ncomp*nrun;
  
  TFile *f = new TFile(gSystem->ExpandPathName(NDData_infile.c_str()),"READ");
  string name = "CovMatrix_"+Extrap[0]->GetNDDataPID();
  
  if(f->Read(name.c_str())==0)
  {
    cout<<"Warning in InitializeDecompSys(): Can't read ND covariance matrix from file.  Setting to zero..."<<endl;
    NDCovMatrix = new TH2D("NDCovMatrix","",ntot,-0.5,ntot-0.5,ntot,-0.5,ntot-0.5);
    InitDecompSys=true;
    gROOT->cd();
    return;
  }
  
  NDCovMatrix = (TH2D*)f->Get(name.c_str());
  
  int n = NDCovMatrix->GetNbinsX();//dimension of error matrix = pid bins x energy bins x number of runs x number of components
  if(n!=ntot)
  {
    cout<<"Error in InitializeDecompSys(): The number of elements in the ND covariance matrix doesn't seem to match your setup, i.e. # PID bins x # energy bins x # runs x 3(NC,CC,Bnue)."<<endl;
    return;
  }
  
  int i,j;
  int ipid,ienergy,irun,icomp;
  
  Background::Background_t bg = Background::kNC;
  vector<double> N;//vector of near detector data with same binning as covariance matrix
  for(i=0;i<n;i++)
  {
    icomp = int(i/(npid*nenergy*nrun));
    irun = int((i - icomp*npid*nenergy*nrun)/(nenergy*npid));
    ipid = int((i - icomp*npid*nenergy*nrun - irun*nenergy*npid)/nenergy);
    ienergy = i%nenergy;
    
    if(icomp==0) bg = Background::kNC;
    else if(icomp==1) bg = Background::kNuMuCC;
    else if(icomp==2) bg = Background::kBNueCC;
    
    N.push_back(Extrap[irun]->NDData[bg]->GetBinContent(ipid+1,ienergy+1));
  }
  
  double temp;
  
  //make the ND covariance matrix a fractional covariance - this way we don't need to do (F/N)(Delta NData) to get the uncertainty, we do the equivalent thing: FPred*(Delta NData/NData) = (F/N)(NData)*(Delta NData/NData) = (F/N)(Delta NData)
  //this way it can be used with grid predictions for FC (where we have only the predictions, not F/N ratio)
  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      if (N[i]!=0&&N[j]!=0) {temp = NDCovMatrix->GetBinContent(i+1,j+1)/(N[i]*N[j]);}
      else temp = 0;
      NDCovMatrix->SetBinContent(i+1,j+1,temp);
    }
  }
  
  InitDecompSys=true;
  
  gROOT->cd();
  
  return;
}
TH2D* ErrorCalc::GetShifts(Background::Background_t bg,bool plus1sigma,string systname,int extrap)
{
  TH2D *hist = new TH2D("hist","hist",1,0,1,1,0,1);
  Initialize();
  if(!Init) return hist;
  
  InitializeSys();
  
  if(FN_NC_Plus1Sigma.count(systname)==0)
  {
    cout<<"Unrecognized systname."<<endl;
    return hist;
  }
  if(unsigned(extrap)>=Extrap.size())
  {
    cout<<"Unrecognized extrap number"<<endl;
    return hist;
  }
  
  if(bg==Background::kNC)
  {
    if(plus1sigma==true)
    {
      return FN_NC_Plus1Sigma[systname][extrap];
    }
    else
    {
      return FN_NC_Minus1Sigma[systname][extrap];
    }
  }
  else if(bg==Background::kNuMuCC)
  {
    if(plus1sigma==true)
    {
      return FN_NuMuCC_Plus1Sigma[systname][extrap];
    }
    else
    {
      return FN_NuMuCC_Minus1Sigma[systname][extrap];
    }
  }
  else if(bg==Background::kBNueCC)
  {
    if(plus1sigma==true)
    {
      return FN_BNueCC_Plus1Sigma[systname][extrap];
    }
    else
    {
      return FN_BNueCC_Minus1Sigma[systname][extrap];
    }
  }
  else if(bg==Background::kNuTauCC)
  {
    if(plus1sigma==true)
    {
      return NuTauCC_MC_Plus1Sigma[systname][extrap];
    }
    else
    {
      return NuTauCC_MC_Minus1Sigma[systname][extrap];
    }
  }
  else if(bg==Background::kNueCC)
  {
    if(plus1sigma==true)
    {
      return NueCC_MC_Plus1Sigma[systname][extrap];
    }
    else
    {
      return NueCC_MC_Minus1Sigma[systname][extrap];
    }
  }
  else
  {
    cout<<"Unknown background type."<<endl;
    return hist;
  }
}
void ErrorCalc::ReadSysFiles_FNExtrap(int n)
{
  int i;
  int k;  /// SG 

  string plf,pln,mnf,mnn,std,systname,filetag;
  systname = FNExtrap_SystName.at(n);
  plf = FNExtrap_FarPlusTag.at(n);
  pln = FNExtrap_NearPlusTag.at(n);
  mnf = FNExtrap_FarMinusTag.at(n);
  mnn = FNExtrap_NearMinusTag.at(n);
  std = FNExtrap_StdTag.at(n);
  filetag=FNExtrap_FileTag.at(n);
  int flag = FNExtrap_Flag.at(n);
  
  vector<string> PreFiles;
  
  unsigned int ie;
  for(ie=0;ie<Extrap.size();ie++)
  {
    PreFiles.push_back(SysFileDir+"/SysFile2D_"+filetag+"_Pre_"+Form("%i",Extrap[ie]->GetRunPeriod())+".root");
    if(gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[ie].c_str())))
    {
      if(flag!=1) cout<<"Warning: "<<PreFiles[ie]<<" doesn't exist."<<endl;
    }
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
  TH2D *farstd,*farpl,*farmn;
  TH2D *nearstd,*nearpl,*nearmn;
  TH2D *fnpl, *fnmn;
  double npot_pre,fpot_pre;
  TTree *treepre;
  double nPOTNear,nPOTFar;
  
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
    nPOTFar=Extrap[ie]->GetFarPOT();
    nPOTNear=Extrap[ie]->GetNearPOT();

    ///SG. Regularisation Histogram 
    TH2D *hreg = new TH2D("hreg", "", np,p, nr,r);
    for(i=1;i<=np;i++)   {      
      for(k=1;k<=nr;k++)    {	
	hreg->SetBinContent(i,k, reg);
      } // for np
    } // for nr
   
    //----------------------

    //flag==1 means systematics related to numuCC extrapolation for signal/tau prediction - by definition 0 for these components
    if(flag==1 || gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[ie].c_str())))//if Pre file does NOT exist
    {

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
   
      delete hreg;

      continue;
    }
    
    fpre = new TFile(gSystem->ExpandPathName(PreFiles[ie].c_str()),"READ");
    
    treepre = (TTree*)fpre->Get("paramtree");
    treepre->SetBranchAddress("nearPOT",&npot_pre);
    treepre->SetBranchAddress("farPOT",&fpot_pre);
    treepre->GetEntry(0);
    
    for(j=0;j<bgs.size();j++)
    {
      name = string(Background::AsString(bgs[j])) + "_" + std + "_Presel/FD_RecoVs" + Extrap[0]->GetPID();
      farstd = (TH2D*)fpre->Get(name.c_str());
      farstd->Scale(nPOTFar/fpot_pre);
      MBH.Rebin2DHist(farstd,np,p,nr,r);

      // Regularise
      farstd->Add(hreg);

      name = string(Background::AsString(bgs[j])) + "_" + plf + "_Presel/FD_RecoVs" + Extrap[0]->GetPID();
      if(plf==std)
      {
        farpl = (TH2D*)farstd->Clone(name.c_str());
      }
      else
      {
        farpl = (TH2D*)fpre->Get(name.c_str());
        farpl->Scale(nPOTFar/fpot_pre);
        MBH.Rebin2DHist(farpl,np,p,nr,r);
	farpl->Add(hreg);
      }
      

      name = string(Background::AsString(bgs[j])) + "_" + mnf + "_Presel/FD_RecoVs" + Extrap[0]->GetPID();
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
	farmn->Add(hreg);
      }
      
      name = string(Background::AsString(bgs[j])) + "_" + std + "_Presel/ND_RecoVs" + Extrap[0]->GetPID();
      nearstd = (TH2D*)fpre->Get(name.c_str());
      nearstd->Scale(nPOTNear/npot_pre);
      MBH.Rebin2DHist(nearstd,np,p,nr,r);
      
      // Regularise
      nearstd->Add(hreg);

      name = string(Background::AsString(bgs[j])) + "_" + pln + "_Presel/ND_RecoVs" + Extrap[0]->GetPID();
      if(pln==std)
      {
        nearpl = (TH2D*)nearstd->Clone(name.c_str());
      }
      else
      {
        nearpl = (TH2D*)fpre->Get(name.c_str());
        nearpl->Scale(nPOTNear/npot_pre);
        MBH.Rebin2DHist(nearpl,np,p,nr,r);
	nearpl->Add(hreg);
      }
      
      name = string(Background::AsString(bgs[j])) + "_" + mnn + "_Presel/ND_RecoVs" + Extrap[0]->GetPID();
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
	nearmn->Add(hreg);
      }
      
      if(bgs[j]==Background::kNC)
      {
        ND_NC_Plus1Sigma[systname].push_back((TH2D*)nearpl->Clone(Form("%s_%i_ND_NC_Plus1",systname.c_str(),ie+1)));
	//        ND_NC_Plus1Sigma[systname][ie]->Add(nearstd,-1.0);
        ND_NC_Plus1Sigma[systname][ie]->Add(nearmn,-1.0);  // SG 
	ND_NC_Plus1Sigma[systname][ie]->Scale(0.5);        // SG 
	ND_NC_Plus1Sigma[systname][ie]->Divide(nearstd);
        ND_NC_Minus1Sigma[systname].push_back((TH2D*)nearmn->Clone(Form("%s_%i_ND_NC_Minus1",systname.c_str(),ie+1)));
        ND_NC_Minus1Sigma[systname][ie]->Add(nearstd,-1.0);
	ND_NC_Minus1Sigma[systname][ie]->Divide(nearstd);
        if(plf==mnf && pln==mnn)
        {
          ND_NC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
        
        FD_NC_Plus1Sigma[systname].push_back((TH2D*)farpl->Clone(Form("%s_%i_FD_NC_Plus1",systname.c_str(),ie+1)));
	//        FD_NC_Plus1Sigma[systname][ie]->Add(farstd,-1.0);
        FD_NC_Plus1Sigma[systname][ie]->Add(farmn,-1.0);  // SG
        FD_NC_Plus1Sigma[systname][ie]->Scale(0.5);       // SG
	FD_NC_Plus1Sigma[systname][ie]->Divide(farstd);
        FD_NC_Minus1Sigma[systname].push_back((TH2D*)farmn->Clone(Form("%s_%i_FD_NC_Minus1",systname.c_str(),ie+1)));
        FD_NC_Minus1Sigma[systname][ie]->Add(farstd,-1.0);
	FD_NC_Minus1Sigma[systname][ie]->Divide(farstd);
        if(plf==mnf && pln==mnn)
        {
          FD_NC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
      }
      else if(bgs[j]==Background::kNuMuCC)
      {

        ND_NuMuCC_Plus1Sigma[systname].push_back((TH2D*)nearpl->Clone(Form("%s_%i_ND_NuMuCC_Plus1",systname.c_str(),ie+1)));
	//        ND_NuMuCC_Plus1Sigma[systname][ie]->Add(nearstd,-1.0);
        ND_NuMuCC_Plus1Sigma[systname][ie]->Add(nearmn,-1.0);    //SG
	ND_NuMuCC_Plus1Sigma[systname][ie]->Scale(0.5);          //SG
	ND_NuMuCC_Plus1Sigma[systname][ie]->Divide(nearstd);
        ND_NuMuCC_Minus1Sigma[systname].push_back((TH2D*)nearmn->Clone(Form("%s_%i_ND_NuMuCC_Minus1",systname.c_str(),ie+1)));
        ND_NuMuCC_Minus1Sigma[systname][ie]->Add(nearstd,-1.0);
	ND_NuMuCC_Minus1Sigma[systname][ie]->Divide(nearstd);
        if(plf==mnf && pln==mnn)
        {
          ND_NuMuCC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
        
        FD_NuMuCC_Plus1Sigma[systname].push_back((TH2D*)farpl->Clone(Form("%s_%i_FD_NuMuCC_Plus1",systname.c_str(),ie+1)));	
	//        FD_NuMuCC_Plus1Sigma[systname][ie]->Add(farstd,-1.0);	
        FD_NuMuCC_Plus1Sigma[systname][ie]->Add(farmn,-1.0);	// SG
        FD_NuMuCC_Plus1Sigma[systname][ie]->Scale(0.5);		// SG
	FD_NuMuCC_Plus1Sigma[systname][ie]->Divide(farstd);

        FD_NuMuCC_Minus1Sigma[systname].push_back((TH2D*)farmn->Clone(Form("%s_%i_FD_NuMuCC_Minus1",systname.c_str(),ie+1)));
	//        FD_NuMuCC_Minus1Sigma[systname][ie]->Add(farstd,-1.0);
        FD_NuMuCC_Minus1Sigma[systname][ie]->Add(farmn,-1.0);   // SG
        FD_NuMuCC_Minus1Sigma[systname][ie]->Scale(0.5);        // SG
	FD_NuMuCC_Minus1Sigma[systname][ie]->Divide(farstd);
        if(plf==mnf && pln==mnn)
        {
          FD_NuMuCC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
      }
      else if(bgs[j]==Background::kBNueCC)
      {
        ND_BNueCC_Plus1Sigma[systname].push_back((TH2D*)nearpl->Clone(Form("%s_%i_ND_BNueCC_Plus1",systname.c_str(),ie+1)));
	//        ND_BNueCC_Plus1Sigma[systname][ie]->Add(nearstd,-1.0);
        ND_BNueCC_Plus1Sigma[systname][ie]->Add(nearmn,-1.0); // SG
        ND_BNueCC_Plus1Sigma[systname][ie]->Scale(0.5);       // SG
	ND_BNueCC_Plus1Sigma[systname][ie]->Divide(nearstd);
        ND_BNueCC_Minus1Sigma[systname].push_back((TH2D*)nearmn->Clone(Form("%s_%i_ND_BNueCC_Minus1",systname.c_str(),ie+1)));
        ND_BNueCC_Minus1Sigma[systname][ie]->Add(nearstd,-1.0);
	ND_BNueCC_Minus1Sigma[systname][ie]->Divide(nearstd);
        if(plf==mnf && pln==mnn)
        {
          ND_BNueCC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
        
        FD_BNueCC_Plus1Sigma[systname].push_back((TH2D*)farpl->Clone(Form("%s_%i_FD_BNueCC_Plus1",systname.c_str(),ie+1)));
	//        FD_BNueCC_Plus1Sigma[systname][ie]->Add(farstd,-1.0);
        FD_BNueCC_Plus1Sigma[systname][ie]->Add(farmn,-1.0);   // SG
        FD_BNueCC_Plus1Sigma[systname][ie]->Scale(0.5);        // SG
	FD_BNueCC_Plus1Sigma[systname][ie]->Divide(farstd);
        FD_BNueCC_Minus1Sigma[systname].push_back((TH2D*)farmn->Clone(Form("%s_%i_FD_BNueCC_Minus1",systname.c_str(),ie+1)));
	FD_BNueCC_Minus1Sigma[systname][ie]->Add(farstd,-1.0);
	FD_BNueCC_Minus1Sigma[systname][ie]->Divide(farstd);
        if(plf==mnf && pln==mnn)
        {
          FD_BNueCC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
      }
      
      farstd->Divide(nearstd);
      farpl->Divide(nearpl);
      farmn->Divide(nearmn);
      
      fnpl = (TH2D*)farpl->Clone("fnpl");
      //      fnpl->Add(farstd,-1.);
      fnpl->Add(farmn,-1.);  // SG
      fnpl->Scale(0.5);      // SG
      fnpl->Divide(farstd);
      
      fnmn = (TH2D*)farmn->Clone("fnmn");
      fnmn->Add(farstd,-1.);
      fnmn->Divide(farstd);
      
      if(plf==mnf && pln==mnn)
      {
        fnmn->Scale(-1.);
      }
      
      if(bgs[j]==Background::kNC)
      {
        FN_NC_Plus1Sigma[systname].push_back((TH2D*)fnpl->Clone(Form("%s_%i_FN_NC_Plus1",systname.c_str(),ie+1)));
        FN_NC_Minus1Sigma[systname].push_back((TH2D*)fnmn->Clone(Form("%s_%i_FN_NC_Minus1",systname.c_str(),ie+1)));
      }
      if(bgs[j]==Background::kNuMuCC)
      {
        FN_NuMuCC_Plus1Sigma[systname].push_back((TH2D*)fnpl->Clone(Form("%s_%i_FN_NuMuCC_Plus1",systname.c_str(),ie+1)));
        FN_NuMuCC_Minus1Sigma[systname].push_back((TH2D*)fnmn->Clone(Form("%s_%i_FN_NuMuCC_Minus1",systname.c_str(),ie+1)));
      }
      if(bgs[j]==Background::kBNueCC)
      {
        FN_BNueCC_Plus1Sigma[systname].push_back((TH2D*)fnpl->Clone(Form("%s_%i_FN_BNueCC_Plus1",systname.c_str(),ie+1)));
        FN_BNueCC_Minus1Sigma[systname].push_back((TH2D*)fnmn->Clone(Form("%s_%i_FN_BNueCC_Minus1",systname.c_str(),ie+1)));
      }
    }
    delete hreg;
  }

  return;
}
void ErrorCalc::ReadSysFiles_FNExtrap_AllRuns(int n)
{

  std::cout << " ErrorCalc::ReadSysFiles_FNExtrap_AllRuns --> " << n << std::endl;
  int i;
  int k; // SG

  string plf,pln,mnf,mnn,std,systname,filetag;
  systname = FNExtrap_SystName.at(n);
  std::cout << " ErrorCalc::ReadSysFiles_FNExtrap_AllRuns  systname  " <<systname << std::endl;
  plf = FNExtrap_FarPlusTag.at(n);
  pln = FNExtrap_NearPlusTag.at(n);
  mnf = FNExtrap_FarMinusTag.at(n);
  mnn = FNExtrap_NearMinusTag.at(n);
  std = FNExtrap_StdTag.at(n);
  filetag=FNExtrap_FileTag.at(n);
  int flag = FNExtrap_Flag.at(n);
  
  vector<string> PreFiles;
  
  unsigned int ie;
  for(ie=0;ie<3;ie++)
  {
    PreFiles.push_back(SysFileDir+"/SysFile2D_"+filetag+"_Pre_"+Form("%i",ie+1)+".root");
    if(gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[ie].c_str())))
    {
      if(flag!=1) cout<<"Warning: "<<PreFiles[ie]<<" doesn't exist."<<endl;
    }
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
  TH2D *totalfarstd=0,*totalfarpl=0,*totalfarmn=0;
  TH2D *totalnearstd=0,*totalnearpl=0,*totalnearmn=0;
  TH2D *fnpl=0, *fnmn=0;
  double npot_pre,fpot_pre;
  TTree *treepre;
  double nPOTNear,nPOTFar;
  int nfiles=0;
  
  nr = Extrap[0]->GetNReco();
  np = Extrap[0]->GetNPID();
  r = new double[nr+1];
  p = new double[np+1];
  for(i=0;i<nr+1;i++)
  {
    r[i] = Extrap[0]->Pred_TotalBkgd->GetYaxis()->GetBinLowEdge(i+1);
  }
  for(i=0;i<np+1;i++)
  {
    p[i] = Extrap[0]->Pred_TotalBkgd->GetXaxis()->GetBinLowEdge(i+1);
  }
  
  ///SG. Regularisation Histogram 
  TH2D *hreg = new TH2D("hreg", "", np,p, nr,r);
  for(i=1;i<=np;i++)   {      
    for(k=1;k<=nr;k++)    {	      
      hreg->SetBinContent(i,k, reg);
    } // for np
  } // for nr


  std::cout << " ErrorCalc::ReadSysFiles_FNExtrap_AllRuns  flag==1 ?  " << std::endl;
  //flag==1 means systematics related to numuCC extrapolation for signal/tau prediction - by definition 0 for these components
  if(flag==1 || (gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[0].c_str())) &&  gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[1].c_str())) && gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[2].c_str()))))//if all three Pre files do NOT exist
  {
    for(ie=0;ie<Extrap.size();ie++)
    {

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

    delete hreg;
    return;
  }
  
  std::cout << " ErrorCalc::ReadSysFiles_FNExtrap_AllRuns  for bgs  " << std::endl;
  for(j=0;j<bgs.size();j++)
  {
    nfiles=0;
    for(ie=0;ie<3;ie++)
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[ie].c_str())))//if file doesn't exist
      {
        continue;
      }

      std::cout << gSystem->ExpandPathName(PreFiles[ie].c_str()) << std::endl;
      fpre = new TFile(gSystem->ExpandPathName(PreFiles[ie].c_str()),"READ");

      std::cout << " -- Tree --  paramtree" << std::endl;
      treepre = (TTree*)fpre->Get("paramtree");
      std::cout << " -- Tree --  nearPOT" << std::endl;
      treepre->SetBranchAddress("nearPOT",&npot_pre);
      std::cout << " -- Tree --  farPOT" << std::endl;
      treepre->SetBranchAddress("farPOT",&fpot_pre);
      std::cout << " -- Tree --  GetEntry(0)" << std::endl;
      treepre->GetEntry(0);
      
      
      std::cout << " -- Extrap["<<ie<<"]->GetFarPOT() " << std::endl;
      nPOTFar=Extrap[ie]->GetFarPOT();
      std::cout << "nPOTFar " << nPOTFar << std::endl;
      std::cout << " -- Extrap["<<ie<<"]->GetNearPOT() "  << std::endl;
      nPOTNear=Extrap[ie]->GetNearPOT();
      std::cout << "nPOTNear " << nPOTNear << std::endl;      
      name = string(Background::AsString(bgs[j])) + "_" + std + "_Presel/FD_RecoVs" + Extrap[0]->GetPID();

      std::cout << "hname " << name << std::endl;
      farstd = (TH2D*)fpre->Get(name.c_str());
      farstd->Scale(nPOTFar/fpot_pre);      
      std::cout << "Rebon2DHist " << name << "  np "<< np << "  nr " << nr << std::endl;
      MBH.Rebin2DHist(farstd,np,p,nr,r);

      if(nfiles==0) totalfarstd = (TH2D*)farstd->Clone();
      else totalfarstd->Add(farstd);

      name = string(Background::AsString(bgs[j])) + "_" + plf + "_Presel/FD_RecoVs" + Extrap[0]->GetPID();
      std::cout << "hname " << name << std::endl;
      if(plf==std)
      {
        farpl = (TH2D*)farstd->Clone(name.c_str());
      }
      else
      {
        farpl = (TH2D*)fpre->Get(name.c_str());
        farpl->Scale(nPOTFar/fpot_pre);
	std::cout << "Rebon2DHist " << name << "  np "<< np << "  nr " << nr << std::endl;
        MBH.Rebin2DHist(farpl,np,p,nr,r);
      }
      if(nfiles==0) totalfarpl = (TH2D*)farpl->Clone();
      else totalfarpl->Add(farpl);
      
      name = string(Background::AsString(bgs[j])) + "_" + mnf + "_Presel/FD_RecoVs" + Extrap[0]->GetPID();
      std::cout << "hname " << name << std::endl;
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
	std::cout << "Rebon2DHist " << name << "  np "<< np << "  nr " << nr << std::endl;
        MBH.Rebin2DHist(farmn,np,p,nr,r);
      }
      if(nfiles==0) totalfarmn = (TH2D*)farmn->Clone();
      else totalfarmn->Add(farmn);
      
      name = string(Background::AsString(bgs[j])) + "_" + std + "_Presel/ND_RecoVs" + Extrap[0]->GetPID();
      std::cout << "hname " << name << std::endl;
      nearstd = (TH2D*)fpre->Get(name.c_str());
      nearstd->Scale(nPOTNear/npot_pre);
      std::cout << "Rebon2DHist " << name << "  np "<< np << "  nr " << nr << std::endl;
      MBH.Rebin2DHist(nearstd,np,p,nr,r);
      if(nfiles==0) totalnearstd = (TH2D*)nearstd->Clone();
      else totalnearstd->Add(nearstd);


      name = string(Background::AsString(bgs[j])) + "_" + pln + "_Presel/ND_RecoVs" + Extrap[0]->GetPID();
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
      if(nfiles==0) totalnearpl = (TH2D*)nearpl->Clone();
      else totalnearpl->Add(nearpl);
      
      name = string(Background::AsString(bgs[j])) + "_" + mnn + "_Presel/ND_RecoVs" + Extrap[0]->GetPID();
      std::cout << "hname " << name << std::endl;
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
	std::cout << "Rebon2DHist " << name << "  np "<< np << "  nr " << nr << std::endl;
        MBH.Rebin2DHist(nearmn,np,p,nr,r);
      }
      if(nfiles==0) totalnearmn = (TH2D*)nearmn->Clone();
      else totalnearmn->Add(nearmn);
      
      nfiles++;
    }
  
    //Regularise
    totalfarstd->Add(hreg);      
    totalfarpl->Add(hreg);
    totalfarmn->Add(hreg);
    
    totalnearstd->Add(hreg);
    totalnearpl->Add(hreg);
    totalnearmn->Add(hreg);
    ///////////////////////
 

    std::cout << " ErrorCalc::ReadSysFiles_FNExtrap_AllRuns  for ExtarpSize  " << std::endl;  
    for(ie=0;ie<Extrap.size();ie++)
    {
      if(bgs[j]==Background::kNC)
      {
        ND_NC_Plus1Sigma[systname].push_back((TH2D*)totalnearpl->Clone(Form("%s_ND_NC_Plus1",systname.c_str())));
	//        ND_NC_Plus1Sigma[systname][ie]->Add(totalnearstd,-1.0);
        ND_NC_Plus1Sigma[systname][ie]->Add(totalnearmn,-1.0);   // SG 
        ND_NC_Plus1Sigma[systname][ie]->Scale(0.5);              // SG 
	ND_NC_Plus1Sigma[systname][ie]->Divide(totalnearstd);
        ND_NC_Minus1Sigma[systname].push_back((TH2D*)totalnearmn->Clone(Form("%s_ND_NC_Minus1",systname.c_str())));
        ND_NC_Minus1Sigma[systname][ie]->Add(totalnearstd,-1.0);
	ND_NC_Minus1Sigma[systname][ie]->Divide(totalnearstd);
        if(plf==mnf && pln==mnn)
        {
          ND_NC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
        
        FD_NC_Plus1Sigma[systname].push_back((TH2D*)totalfarpl->Clone(Form("%s_FD_NC_Plus1",systname.c_str())));
	//        FD_NC_Plus1Sigma[systname][ie]->Add(totalfarstd,-1.0);
	FD_NC_Plus1Sigma[systname][ie]->Add(totalfarmn,-1.0);   // SG
	FD_NC_Plus1Sigma[systname][ie]->Scale(0.5);             // SG
	FD_NC_Plus1Sigma[systname][ie]->Divide(totalfarstd);
        FD_NC_Minus1Sigma[systname].push_back((TH2D*)totalfarmn->Clone(Form("%s_FD_NC_Minus1",systname.c_str())));
        FD_NC_Minus1Sigma[systname][ie]->Add(totalfarstd,-1.0);
	FD_NC_Minus1Sigma[systname][ie]->Divide(totalfarstd);
        if(plf==mnf && pln==mnn)
        {
          FD_NC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
      }
      else if(bgs[j]==Background::kNuMuCC)
      {
        ND_NuMuCC_Plus1Sigma[systname].push_back((TH2D*)totalnearpl->Clone(Form("%s_ND_NuMuCC_Plus1",systname.c_str())));
	//        ND_NuMuCC_Plus1Sigma[systname][ie]->Add(totalnearstd,-1.0);
        ND_NuMuCC_Plus1Sigma[systname][ie]->Add(totalnearmn,-1.0);  // SG
        ND_NuMuCC_Plus1Sigma[systname][ie]->Scale(0.5);             // SG
	ND_NuMuCC_Plus1Sigma[systname][ie]->Divide(totalnearstd);
        ND_NuMuCC_Minus1Sigma[systname].push_back((TH2D*)totalnearmn->Clone(Form("%s_ND_NuMuCC_Minus1",systname.c_str())));
        ND_NuMuCC_Minus1Sigma[systname][ie]->Add(totalnearstd,-1.0);
	ND_NuMuCC_Minus1Sigma[systname][ie]->Divide(totalnearstd);
        if(plf==mnf && pln==mnn)
        {
          ND_NuMuCC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
        
        FD_NuMuCC_Plus1Sigma[systname].push_back((TH2D*)totalfarpl->Clone(Form("%s_FD_NuMuCC_Plus1",systname.c_str())));
	//        FD_NuMuCC_Plus1Sigma[systname][ie]->Add(totalfarstd,-1.0);
        FD_NuMuCC_Plus1Sigma[systname][ie]->Add(totalfarmn,-1.0);      // SG
        FD_NuMuCC_Plus1Sigma[systname][ie]->Scale(0.5);                // SG
	FD_NuMuCC_Plus1Sigma[systname][ie]->Divide(totalfarstd);
        FD_NuMuCC_Minus1Sigma[systname].push_back((TH2D*)totalfarmn->Clone(Form("%s_FD_NuMuCC_Minus1",systname.c_str())));
        FD_NuMuCC_Minus1Sigma[systname][ie]->Add(totalfarstd,-1.0);
	FD_NuMuCC_Minus1Sigma[systname][ie]->Divide(totalfarstd);
        if(plf==mnf && pln==mnn)
        {
          FD_NuMuCC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
      }
      else if(bgs[j]==Background::kBNueCC)
      {
        ND_BNueCC_Plus1Sigma[systname].push_back((TH2D*)totalnearpl->Clone(Form("%s_ND_BNueCC_Plus1",systname.c_str())));
	//        ND_BNueCC_Plus1Sigma[systname][ie]->Add(totalnearstd,-1.0);
        ND_BNueCC_Plus1Sigma[systname][ie]->Add(totalnearmn,-1.0);   // SG
        ND_BNueCC_Plus1Sigma[systname][ie]->Scale(0.5);              // SG
	ND_BNueCC_Plus1Sigma[systname][ie]->Divide(totalnearstd);
        ND_BNueCC_Minus1Sigma[systname].push_back((TH2D*)totalnearmn->Clone(Form("%s_ND_BNueCC_Minus1",systname.c_str())));
        ND_BNueCC_Minus1Sigma[systname][ie]->Add(totalnearstd,-1.0);
	ND_BNueCC_Minus1Sigma[systname][ie]->Divide(totalnearstd);
        if(plf==mnf && pln==mnn)
        {
          ND_BNueCC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
        
        FD_BNueCC_Plus1Sigma[systname].push_back((TH2D*)totalfarpl->Clone(Form("%s_FD_BNueCC_Plus1",systname.c_str())));
	//        FD_BNueCC_Plus1Sigma[systname][ie]->Add(totalfarstd,-1.0);
        FD_BNueCC_Plus1Sigma[systname][ie]->Add(totalfarmn,-1.0);    // SG
        FD_BNueCC_Plus1Sigma[systname][ie]->Scale(0.5);              // SG
	FD_BNueCC_Plus1Sigma[systname][ie]->Divide(totalfarstd);
        FD_BNueCC_Minus1Sigma[systname].push_back((TH2D*)totalfarmn->Clone(Form("%s_FD_BNueCC_Minus1",systname.c_str())));
        FD_BNueCC_Minus1Sigma[systname][ie]->Add(totalfarstd,-1.0);
	FD_BNueCC_Minus1Sigma[systname][ie]->Divide(totalfarstd);
        if(plf==mnf && pln==mnn)
        {
          FD_BNueCC_Minus1Sigma[systname][ie]->Scale(-1.);
        }
      }
    }
    

    totalfarstd->Divide(totalnearstd);
    totalfarpl->Divide(totalnearpl);
    totalfarmn->Divide(totalnearmn);

    
    fnpl = (TH2D*)totalfarpl->Clone("fnpl");
    //    fnpl->Add(totalfarstd,-1.);
    fnpl->Add(totalfarmn,-1.);  // SG
    fnpl->Scale(0.5);           // SG
    fnpl->Divide(totalfarstd);
    
    fnmn = (TH2D*)totalfarmn->Clone("fnmn");
    fnmn->Add(totalfarstd,-1.);
    fnmn->Divide(totalfarstd);
    
    if(plf==mnf && pln==mnn)
    {
      fnmn->Scale(-1.);
    }
   
    std::cout << " ErrorCalc::ReadSysFiles_FNExtrap_AllRuns  for ExtrapSize  " << std::endl;
    for(ie=0;ie<Extrap.size();ie++)
    {
      if(bgs[j]==Background::kNC)
      {
        FN_NC_Plus1Sigma[systname].push_back((TH2D*)fnpl->Clone(Form("%s_FN_NC_Plus1",systname.c_str())));
        FN_NC_Minus1Sigma[systname].push_back((TH2D*)fnmn->Clone(Form("%s_FN_NC_Minus1",systname.c_str())));
      }
      if(bgs[j]==Background::kNuMuCC)
      {
        FN_NuMuCC_Plus1Sigma[systname].push_back((TH2D*)fnpl->Clone(Form("%s_FN_NuMuCC_Plus1",systname.c_str())));
        FN_NuMuCC_Minus1Sigma[systname].push_back((TH2D*)fnmn->Clone(Form("%s_FN_NuMuCC_Minus1",systname.c_str())));
      }
      if(bgs[j]==Background::kBNueCC)
      {
        FN_BNueCC_Plus1Sigma[systname].push_back((TH2D*)fnpl->Clone(Form("%s_FN_BNueCC_Plus1",systname.c_str())));
        FN_BNueCC_Minus1Sigma[systname].push_back((TH2D*)fnmn->Clone(Form("%s_FN_BNueCC_Minus1",systname.c_str())));
      }
    }
  }
  
  delete hreg;  
  return;
}
void ErrorCalc::ReadSysFiles_Appearance(int n)
{
  int i;
  int k; // SG

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
      FidFiles.push_back(SysFileDir+"/SysFile2D_"+filetag+"_Fid_"+Form("%i",Extrap[ie]->GetRunPeriod())+".root");
      CClikeFiles.push_back(SysFileDir+"/SysFile2D_"+filetag+"_CCLike_"+Form("%i",Extrap[ie]->GetRunPeriod())+".root");
      if(gSystem->AccessPathName(gSystem->ExpandPathName(FidFiles[ie].c_str())))
      {
	cout<<"Warning: "<<FidFiles[ie]<<" doesn't exist."<<endl;
      }
      if(gSystem->AccessPathName(gSystem->ExpandPathName(CClikeFiles[ie].c_str())))
      {
	cout<<"Warning: "<<CClikeFiles[ie]<<" doesn't exist."<<endl;
      }
    }
    if(flag==2)//uncertainty in selected tau or signal FD MC
    {
      PreFiles.push_back(SysFileDir+"/SysFile2D_"+filetag+"_Pre_"+Form("%i",Extrap[ie]->GetRunPeriod())+".root");
      if(gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[ie].c_str())))
      {
	cout<<"Warning: "<<PreFiles[ie]<<" doesn't exist."<<endl;
      }
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
  TH2D *prestd,*prepl,*premn;
  double sum,temp;
  int ir,it,ip;
  TH1D *ccstd,*ccpl,*ccmn;
  TH1D *data;
  TH2D *farfidccstd,*farcclikestd,*nearcclikestd;
  TH2D *farfidccpl,*farcclikepl,*nearcclikepl;
  TH2D *farfidccmn,*farcclikemn,*nearcclikemn;
  TH2D *h2;
  double sum_plus_tau,sum_minus_tau,sum_plus_nue,sum_minus_nue;
  string ccname = "CC";
  
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
    nPOTFar=Extrap[ie]->GetFarPOT();
    nPOTNear=Extrap[ie]->GetNearPOT();

    ///SG. Regularisation Histogram 
    TH2D *hreg = new TH2D("hreg", "", np,p, nr,r);
    for(i=1;i<=np;i++)   {      
      for(k=1;k<=nr;k++)    {	
	hreg->SetBinContent(i,k, reg);
      } // for np
    } // for nr

    //------------------------------------
  
    if(flag==0)
    {
      //create empty histograms for overall signal/tau shifts
      NueCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NueCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NuTauCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NuTauCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      
      //creating empty numuCC shifts
      Pred_CC_Fid_Plus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Plus1",systname.c_str(),ie+1),"",nt,t));
      Pred_CC_Fid_Minus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Minus1",systname.c_str(),ie+1),"",nt,t));
    }
    if(flag==1)
    {
      //create empty histograms for overall signal/tau shifts
      NueCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NueCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NuTauCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NuTauCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      
      if(gSystem->AccessPathName(gSystem->ExpandPathName(CClikeFiles[ie].c_str())) || gSystem->AccessPathName(gSystem->ExpandPathName(FidFiles[ie].c_str())))//if CClike or Fid file does NOT exist
      {
        Pred_CC_Fid_Plus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Plus1",systname.c_str(),ie+1),"",nt,t));
        Pred_CC_Fid_Minus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Minus1",systname.c_str(),ie+1),"",nt,t));
        
	delete hreg;

        continue;
      }
     
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
      
      name = "NuMuCC_" + std + "_Fid/FD_TrueVsReco";
      farfidccstd = (TH2D*)ffid->Get(name.c_str());
      farfidccstd->Scale(nPOTFar/fpot_fid);
      
      name = "NuMuCC_" + std + "_" + ccname + "/FD_TrueVsReco";
      farcclikestd = (TH2D*)fcc->Get(name.c_str());
      name = "NC_" + std + "_" + ccname + "/FD_TrueVsReco";
      h2 = (TH2D*)fcc->Get(name.c_str());
      farcclikestd->Add(h2);
      name = "BNueCC_" + std + "_" + ccname + "/FD_TrueVsReco";
      h2 = (TH2D*)fcc->Get(name.c_str());
      farcclikestd->Add(h2);
      farcclikestd->Scale(nPOTFar/fpot_cc);
      
      name = "NuMuCC_" + std + "_" + ccname + "/ND_TrueVsReco";
      nearcclikestd = (TH2D*)fcc->Get(name.c_str());
      name = "NC_" + std + "_" + ccname + "/ND_TrueVsReco";
      h2 = (TH2D*)fcc->Get(name.c_str());
      nearcclikestd->Add(h2);
      name = "BNueCC_" + std + "_" + ccname + "/ND_TrueVsReco";
      h2 = (TH2D*)fcc->Get(name.c_str());
      nearcclikestd->Add(h2);
      nearcclikestd->Scale(nPOTNear/npot_cc);
      
      //Regularise
      if(farfidccstd->Integral()  >0)  farfidccstd->Add(hreg);      
      if(farcclikestd->Integral() >0)  farcclikestd->Add(hreg);
      if(nearcclikestd->Integral()>0)  nearcclikestd->Add(hreg);

      ccstd = new TH1D("ccstd","",nt,t);
      for(it=0;it<nt;it++)
      {
        sum=0;
        for(ir=0;ir<farcclikestd->GetNbinsX();ir++)
        {
          if(nearcclikestd->Integral(ir+1,ir+1,1,nt)>0)
          {
            sum += (data->GetBinContent(ir+1)/nearcclikestd->Integral(ir+1,ir+1,1,nt))*farcclikestd->GetBinContent(ir+1,it+1);
          }
        }
        temp=0;
        if(farcclikestd->Integral(1,farcclikestd->GetNbinsX(),it+1,it+1)>0)
        {
          temp = farfidccstd->Integral(1,farfidccstd->GetNbinsX(),it+1,it+1)*(sum/farcclikestd->Integral(1,farcclikestd->GetNbinsX(),it+1,it+1));
        }
        ccstd->SetBinContent(it+1,temp);
      }
      
      name = "NuMuCC_" + plf + "_Fid/FD_TrueVsReco";
      if(plf==std)
      {
        farfidccpl = (TH2D*)farfidccstd->Clone(name.c_str());
      }
      else
      {
        farfidccpl = (TH2D*)ffid->Get(name.c_str());
        farfidccpl->Scale(nPOTFar/fpot_fid);
	if(farfidccpl->Integral()>0)     farfidccpl->Add(hreg);
      }
      
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
	if(farcclikepl->Integral()>0)     farcclikepl->Add(hreg);
      }
      
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
	if(nearcclikepl->Integral()>0)     nearcclikepl->Add(hreg);
      }
      
      ccpl = new TH1D("ccpl","",nt,t);
      for(it=0;it<nt;it++)
      {
        sum=0;
        for(ir=0;ir<farcclikepl->GetNbinsX();ir++)
        {
          if(nearcclikepl->Integral(ir+1,ir+1,1,nt)>0)
          {
            sum += (data->GetBinContent(ir+1)/nearcclikepl->Integral(ir+1,ir+1,1,nt))*farcclikepl->GetBinContent(ir+1,it+1);
          }
        }
        temp=0;
        if(farcclikepl->Integral(1,farcclikepl->GetNbinsX(),it+1,it+1)>0)
        {
          temp = farfidccpl->Integral(1,farfidccpl->GetNbinsX(),it+1,it+1)*(sum/farcclikepl->Integral(1,farcclikepl->GetNbinsX(),it+1,it+1));
        }
        ccpl->SetBinContent(it+1,temp);
      }
      
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
	if(farfidccmn->Integral()>0)   farfidccmn->Add(hreg);
      }
      
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
	if(farcclikemn->Integral()>0)   farcclikemn->Add(hreg);
      }
      
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
	if(nearcclikemn->Integral()>0)   nearcclikemn->Add(hreg);
      }
      
      ccmn = new TH1D("ccmn","",nt,t);
      for(it=0;it<nt;it++)
      {
        sum=0;
        for(ir=0;ir<farcclikemn->GetNbinsX();ir++)
        {
          if(nearcclikemn->Integral(ir+1,ir+1,1,nt)>0)
          {
            sum += (data->GetBinContent(ir+1)/nearcclikemn->Integral(ir+1,ir+1,1,nt))*farcclikemn->GetBinContent(ir+1,it+1);
          }
        }
        temp=0;
        if(farcclikemn->Integral(1,farcclikemn->GetNbinsX(),it+1,it+1)>0)
        {
          temp = farfidccmn->Integral(1,farfidccmn->GetNbinsX(),it+1,it+1)*(sum/farcclikemn->Integral(1,farcclikemn->GetNbinsX(),it+1,it+1));
        }
        ccmn->SetBinContent(it+1,temp);
      }
      
      //      ccpl->Add(ccstd,-1.);
      ccpl->Add(ccmn,-1.);  // SG 
      ccpl->Scale(0.5);     // SG 
      ccpl->Divide(ccstd);
      ccmn->Add(ccstd,-1.);
      ccmn->Divide(ccstd);
      
      if(plf==mnf && pln==mnn)
      {
        ccmn->Scale(-1.);
      }
      
      Pred_CC_Fid_Plus1Sigma[systname].push_back((TH1D*)ccpl->Clone(Form("%s_%i_Pred_CC_Fid_Plus1",systname.c_str(),ie+1)));
      Pred_CC_Fid_Minus1Sigma[systname].push_back((TH1D*)ccmn->Clone(Form("%s_%i_Pred_CC_Fid_Minus1",systname.c_str(),ie+1)));
      
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
    
    if(flag==2)
    {
      //creating empty numuCC shifts
      Pred_CC_Fid_Plus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Plus1",systname.c_str(),ie+1),"",nt,t));
      Pred_CC_Fid_Minus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Minus1",systname.c_str(),ie+1),"",nt,t));
      
      if(gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[ie].c_str())))//if Pre file does NOT exist
      {
        NueCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
        NueCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
        NuTauCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
        NuTauCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
   
	delete hreg;

        continue;
      }
   
      fpre = new TFile(gSystem->ExpandPathName(PreFiles[ie].c_str()),"READ");
      treepre = (TTree*)fpre->Get("paramtree");
      treepre->SetBranchAddress("nearPOT",&npot_pre);
      treepre->SetBranchAddress("farPOT",&fpot_pre);
      treepre->GetEntry(0);
      
      for(j=0;j<app.size();j++)
      {
        name = string(Background::AsString(app[j])) + "_" + std + "_Presel/FD_RecoVs" + Extrap[0]->GetPID();
        prestd = (TH2D*)fpre->Get(name.c_str());
        prestd->Scale(nPOTFar/fpot_pre);
        MBH.Rebin2DHist(prestd,np,p,nr,r);
	
        name = string(Background::AsString(app[j])) + "_" + plf + "_Presel/FD_RecoVs" + Extrap[0]->GetPID();
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
      
        name = string(Background::AsString(app[j])) + "_" + mnf + "_Presel/FD_RecoVs" + Extrap[0]->GetPID();
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
        

	// Regularise ------
	prestd->Add(hreg);
	prepl->Add(hreg);
	premn->Add(hreg);
	//------------------
        
	//	prepl->Add(prestd,-1.);
	prepl->Add(premn,-1.);  // SG 
	prepl->Scale(0.5);      // SG
        prepl->Divide(prestd);
        premn->Add(prestd,-1.);
        premn->Divide(prestd);
        
        if(plf==mnf && pln==mnn)
        {
          premn->Scale(-1.);
        }
        
        if(app[j]==Background::kNueCC)
        {
          NueCC_MC_Plus1Sigma[systname].push_back((TH2D*)prepl->Clone(Form("%s_%i_NueCC_MC_Plus1",systname.c_str(),ie+1)));
          NueCC_MC_Minus1Sigma[systname].push_back((TH2D*)premn->Clone(Form("%s_%i_NueCC_MC_Minus1",systname.c_str(),ie+1)));
        }
        else if(app[j]==Background::kNuTauCC)
        {
          NuTauCC_MC_Plus1Sigma[systname].push_back((TH2D*)prepl->Clone(Form("%s_%i_NuTauCC_MC_Plus1",systname.c_str(),ie+1)));
          NuTauCC_MC_Minus1Sigma[systname].push_back((TH2D*)premn->Clone(Form("%s_%i_NuTauCC_MC_Minus1",systname.c_str(),ie+1)));
        }
      }
    }
    
    delete hreg;
  }
  
  return;
}
void ErrorCalc::ReadSysFiles_Appearance_AllRuns(int n)
{
  int i;
  int k; // SG

  string plf,pln,mnf,mnn,std,systname,filetag;
  systname = FNExtrap_SystName.at(n);

  std::cout <<  " ErrorCalc::ReadSysFiles_Appearance_AllRuns " <<  systname << std::endl;
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
      std::cout <<  " ErrorCalc::ReadSysFiles_Appearance_AllRuns  flag==1" << std::endl;
      FidFiles.push_back(SysFileDir+"/SysFile2D_"+filetag+"_Fid_"+Form("%i",Extrap[ie]->GetRunPeriod())+".root");
      CClikeFiles.push_back(SysFileDir+"/SysFile2D_"+filetag+"_CCLike_"+Form("%i",Extrap[ie]->GetRunPeriod())+".root");
      if(gSystem->AccessPathName(gSystem->ExpandPathName(FidFiles[ie].c_str())))
      {
        cout<<"Warning: "<<FidFiles[ie]<<" doesn't exist."<<endl;
      }
      if(gSystem->AccessPathName(gSystem->ExpandPathName(CClikeFiles[ie].c_str())))
      {
	cout<<"Warning: "<<CClikeFiles[ie]<<" doesn't exist."<<endl;
      }
    }
    else if(flag==2)//uncertainty in selected tau or signal FD MC
    {
      std::cout <<  " ErrorCalc::ReadSysFiles_Appearance_AllRuns  flag==2" << std::endl;
      PreFiles.push_back(SysFileDir+"/SysFile2D_"+filetag+"_Pre_"+Form("%i",Extrap[ie]->GetRunPeriod())+".root");
      if(gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[ie].c_str())))
      {
	cout<<"Warning: "<<PreFiles[ie]<<" doesn't exist."<<endl;
      }
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
  TH2D *totalprestd=0,*totalprepl=0,*totalpremn=0;
  double sum,temp;
  int ir,it,ip;
  TH1D *ccstd=0,*ccpl=0,*ccmn=0;
  TH1D *data=0;
  TH2D *farfidccstd=0,*farcclikestd=0,*nearcclikestd=0;
  TH2D *farfidccpl=0,*farcclikepl=0,*nearcclikepl=0;
  TH2D *farfidccmn=0,*farcclikemn=0,*nearcclikemn=0;
  TH1D *totaldata=0;
  TH2D *totalfarfidccstd=0,*totalfarcclikestd=0,*totalnearcclikestd=0;
  TH2D *totalfarfidccpl=0,*totalfarcclikepl=0,*totalnearcclikepl=0;
  TH2D *totalfarfidccmn=0,*totalfarcclikemn=0,*totalnearcclikemn=0;
  TH2D *h2;
  int nfiles;
  double sum_plus_tau,sum_minus_tau,sum_plus_nue,sum_minus_nue;
  string ccname = "CC";
  
  nr = Extrap[0]->GetNReco();
  np = Extrap[0]->GetNPID();
  nt = Extrap[0]->GetNTrue();
  std::cout <<  " ErrorCalc::ReadSysFiles_Appearance_AllRuns  nr " << nr << "  np " << np << "  nt " << nt  << std::endl;
  r = new double[nr+1];
  p = new double[np+1];
  t = new double[nt+1];
  for(i=0;i<nr+1;i++)
  {
    r[i] = Extrap[0]->Pred_TotalBkgd->GetYaxis()->GetBinLowEdge(i+1);
  }
  for(i=0;i<np+1;i++)
  {
    p[i] = Extrap[0]->Pred_TotalBkgd->GetXaxis()->GetBinLowEdge(i+1);
  }
  for(i=0;i<nt+1;i++)
  {
    t[i] = Extrap[0]->Pred_3D[Extrapolate2D::qNuMuToNuTau]->GetZaxis()->GetBinLowEdge(i+1);
  }

  ///SG. Regularisation Histogram 
  TH2D *hreg = new TH2D("hreg", "", np,p, nr,r);
  for(i=1;i<=np;i++)   {      
    for(k=1;k<=nr;k++)    {	
      hreg->SetBinContent(i,k, reg);
    } // for np
  } // for nr

  //----------------------

  if(flag==0)
  {
    std::cout <<  " ErrorCalc::ReadSysFiles_Appearance_AllRuns  flag==0" << std::endl;
    for(ie=0;ie<Extrap.size();ie++)
    {
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
    std::cout <<  " ErrorCalc::ReadSysFiles_Appearance_AllRuns  flag==1" << std::endl;
    //create empty histograms for overall signal/tau shifts
    for(ie=0;ie<Extrap.size();ie++)
    {
      std::cout <<  " ErrorCalc::ReadSysFiles_Appearance_AllRuns  flag==1 --NueCC NuTauCC  ie " << ie << std::endl;
      NueCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NueCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NuTauCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
      NuTauCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
    }
    
    if(gSystem->AccessPathName(gSystem->ExpandPathName(CClikeFiles[0].c_str())) && gSystem->AccessPathName(gSystem->ExpandPathName(CClikeFiles[1].c_str())) && gSystem->AccessPathName(gSystem->ExpandPathName(CClikeFiles[2].c_str())) && gSystem->AccessPathName(gSystem->ExpandPathName(FidFiles[0].c_str())) && gSystem->AccessPathName(gSystem->ExpandPathName(FidFiles[1].c_str())) && gSystem->AccessPathName(gSystem->ExpandPathName(FidFiles[2].c_str())))//if none of these files exist
    {
      for(ie=0;ie<Extrap.size();ie++)
      {
	std::cout <<  " ErrorCalc::ReadSysFiles_Appearance_AllRuns  flag==1 -- PredCC  ie " << ie << std::endl;
        Pred_CC_Fid_Plus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Plus1",systname.c_str(),ie+1),"",nt,t));
        Pred_CC_Fid_Minus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Minus1",systname.c_str(),ie+1),"",nt,t));
      }
     
      delete hreg;
      return;
    }
    
    nfiles=0;
    for(ie=0;ie<3;ie++)
    {
      if(gSystem->AccessPathName(gSystem->ExpandPathName(CClikeFiles[ie].c_str())) && gSystem->AccessPathName(gSystem->ExpandPathName(FidFiles[ie].c_str())))//if files don't exist
      { 
        continue;
      }
      
      std::cout << "  -----  " <<  ie  << " " << gSystem->ExpandPathName(CClikeFiles[ie].c_str()) << std::endl;

      fcc = new TFile(gSystem->ExpandPathName(CClikeFiles[ie].c_str()),"READ");
      treecc = (TTree*)fcc->Get("paramtree");
      treecc->SetBranchAddress("nearPOT",&npot_cc);
      treecc->SetBranchAddress("farPOT",&fpot_cc);
      treecc->GetEntry(0);

      std::cout << gSystem->ExpandPathName(FidFiles[ie].c_str()) << std::endl;

      ffid = new TFile(gSystem->ExpandPathName(FidFiles[ie].c_str()),"READ");

      std::cout << "paramtree" << std::endl;
      treefid = (TTree*)ffid->Get("paramtree");
      std::cout << "paramtree nearPOT" << std::endl;
      treefid->SetBranchAddress("nearPOT",&npot_fid);
      std::cout << "paramtree farPOT" << std::endl;
      treefid->SetBranchAddress("farPOT",&fpot_fid);
      std::cout << "paramtree GetEntry" << std::endl;
      treefid->GetEntry(0);
      
      std::cout << Form("data_%i",ie) << std::endl;
      std::cout << Extrap[ie]->NDData_Reco_CClike->GetName() << std::endl;
      data = (TH1D*)Extrap[ie]->NDData_Reco_CClike->Clone(Form("data_%i",ie));
      //should scale data!
      std::cout << " nfiles " << nfiles << std::endl;  
      if(nfiles==0) totaldata = (TH1D*)data->Clone();
      else totaldata->Add(data);
      
      std::cout << "Extrap[ie]->GetFarPOT()" << std::endl;
      nPOTFar=Extrap[ie]->GetFarPOT();
      std::cout << "Extrap[ie]->GetNearPOT()" << std::endl;
      nPOTNear=Extrap[ie]->GetNearPOT();
      
      std::cout << " POT Far Near " << nPOTFar << " " << nPOTNear << std::endl;
      name = "NuMuCC_" + std + "_Fid/FD_TrueVsReco";
      farfidccstd = (TH2D*)ffid->Get(name.c_str());
      farfidccstd->Scale(nPOTFar/fpot_fid);
      if(nfiles==0) totalfarfidccstd = (TH2D*)farfidccstd->Clone();
      else totalfarfidccstd->Add(farfidccstd);
      
      name = "NuMuCC_" + std + "_" + ccname + "/FD_TrueVsReco";
      farcclikestd = (TH2D*)fcc->Get(name.c_str());
      name = "NC_" + std + "_" + ccname + "/FD_TrueVsReco";
      h2 = (TH2D*)fcc->Get(name.c_str());
      farcclikestd->Add(h2);
      name = "BNueCC_" + std + "_" + ccname + "/FD_TrueVsReco";
      h2 = (TH2D*)fcc->Get(name.c_str());
      farcclikestd->Add(h2);
      farcclikestd->Scale(nPOTFar/fpot_cc);
      if(nfiles==0) totalfarcclikestd = (TH2D*)farcclikestd->Clone();
      else totalfarcclikestd->Add(farcclikestd);
      
      name = "NuMuCC_" + std + "_" + ccname + "/ND_TrueVsReco";
      nearcclikestd = (TH2D*)fcc->Get(name.c_str());
      name = "NC_" + std + "_" + ccname + "/ND_TrueVsReco";
      h2 = (TH2D*)fcc->Get(name.c_str());
      nearcclikestd->Add(h2);
      name = "BNueCC_" + std + "_" + ccname + "/ND_TrueVsReco";
      h2 = (TH2D*)fcc->Get(name.c_str());
      nearcclikestd->Add(h2);
      nearcclikestd->Scale(nPOTNear/npot_cc);
      if(nfiles==0) totalnearcclikestd = (TH2D*)nearcclikestd->Clone();
      else totalnearcclikestd->Add(nearcclikestd);
      
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
      if(nfiles==0) totalfarfidccpl = (TH2D*)farfidccpl->Clone();
      else totalfarfidccpl->Add(farfidccpl);
      
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
      if(nfiles==0) totalfarcclikepl = (TH2D*)farcclikepl->Clone();
      else totalfarcclikepl->Add(farcclikepl);
      
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
      if(nfiles==0) totalnearcclikepl = (TH2D*)nearcclikepl->Clone();
      else totalnearcclikepl->Add(nearcclikepl);
      
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
      if(nfiles==0) totalfarfidccmn = (TH2D*)farfidccmn->Clone();
      else totalfarfidccmn->Add(farfidccmn);
      
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
      if(nfiles==0) totalfarcclikemn = (TH2D*)farcclikemn->Clone();
      else totalfarcclikemn->Add(farcclikemn);
      
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
      if(nfiles==0) totalnearcclikemn = (TH2D*)nearcclikemn->Clone();
      else totalnearcclikemn->Add(nearcclikemn);
      
      nfiles++;
    }
    
    ccstd = new TH1D("ccstd","",nt,t);
    for(it=0;it<nt;it++)
    {
      sum=0;
      for(ir=0;ir<totalfarcclikestd->GetNbinsX();ir++)
      {
        if(totalnearcclikestd->Integral(ir+1,ir+1,1,nt)>0)
        {
          sum += (totaldata->GetBinContent(ir+1)/totalnearcclikestd->Integral(ir+1,ir+1,1,nt))*totalfarcclikestd->GetBinContent(ir+1,it+1);
        }
      }
      temp=0;
      if(totalfarcclikestd->Integral(1,totalfarcclikestd->GetNbinsX(),it+1,it+1)>0)
      {
        temp = totalfarfidccstd->Integral(1,totalfarfidccstd->GetNbinsX(),it+1,it+1)*(sum/totalfarcclikestd->Integral(1,totalfarcclikestd->GetNbinsX(),it+1,it+1));
      }
      ccstd->SetBinContent(it+1,temp);
    }
    
    ccpl = new TH1D("ccpl","",nt,t);
    for(it=0;it<nt;it++)
    {
      sum=0;
      for(ir=0;ir<totalfarcclikepl->GetNbinsX();ir++)
      {
        if(totalnearcclikepl->Integral(ir+1,ir+1,1,nt)>0)
        {
          sum += (totaldata->GetBinContent(ir+1)/totalnearcclikepl->Integral(ir+1,ir+1,1,nt))*totalfarcclikepl->GetBinContent(ir+1,it+1);
        }
      }
      temp=0;
      if(totalfarcclikepl->Integral(1,totalfarcclikepl->GetNbinsX(),it+1,it+1)>0)
      {
        temp = totalfarfidccpl->Integral(1,totalfarfidccpl->GetNbinsX(),it+1,it+1)*(sum/totalfarcclikepl->Integral(1,totalfarcclikepl->GetNbinsX(),it+1,it+1));
      }
      ccpl->SetBinContent(it+1,temp);
    }
    
    ccmn = new TH1D("ccmn","",nt,t);
    for(it=0;it<nt;it++)
    {
      sum=0;
      for(ir=0;ir<totalfarcclikemn->GetNbinsX();ir++)
      {
        if(totalnearcclikemn->Integral(ir+1,ir+1,1,nt)>0)
        {
          sum += (totaldata->GetBinContent(ir+1)/totalnearcclikemn->Integral(ir+1,ir+1,1,nt))*totalfarcclikemn->GetBinContent(ir+1,it+1);
        }
      }
      temp=0;
      if(totalfarcclikemn->Integral(1,totalfarcclikemn->GetNbinsX(),it+1,it+1)>0)
      {
        temp = totalfarfidccmn->Integral(1,totalfarfidccmn->GetNbinsX(),it+1,it+1)*(sum/totalfarcclikemn->Integral(1,totalfarcclikemn->GetNbinsX(),it+1,it+1));
      }
      ccmn->SetBinContent(it+1,temp);
    }

    /// ??????? 
    //    ccpl->Add(ccstd,-1.);
    ccpl->Add(ccmn,-1.);  // SG
    ccpl->Scale(0.5);     // SG 
    ccpl->Divide(ccstd);
    ccmn->Add(ccstd,-1.);
    ccmn->Divide(ccstd);
    /// ???????

    if(plf==mnf && pln==mnn)
    {
      ccmn->Scale(-1.);
    }
    
    for(ie=0;ie<Extrap.size();ie++)
    {
      Pred_CC_Fid_Plus1Sigma[systname].push_back((TH1D*)ccpl->Clone(Form("%s_%i_Pred_CC_Fid_Plus1",systname.c_str(),ie+1)));
      Pred_CC_Fid_Minus1Sigma[systname].push_back((TH1D*)ccmn->Clone(Form("%s_%i_Pred_CC_Fid_Minus1",systname.c_str(),ie+1)));
    }
    
    for(ie=0;ie<Extrap.size();ie++)
    {
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
    std::cout <<  " ErrorCalc::ReadSysFiles_Appearance_AllRuns  flag==2" << std::endl;
    //create empty histograms for numuCC shift
    for(ie=0;ie<Extrap.size();ie++)
    {
      Pred_CC_Fid_Plus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Plus1",systname.c_str(),ie+1),"",nt,t));
      Pred_CC_Fid_Minus1Sigma[systname].push_back(new TH1D(Form("%s_%i_Pred_CC_Fid_Minus1",systname.c_str(),ie+1),"",nt,t));
    }
    
    if(gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[0].c_str())) && gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[1].c_str())) && gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[2].c_str())))//if none of these files exist
    {
      for(ie=0;ie<Extrap.size();ie++)
      {
        NueCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
        NueCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NueCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
        NuTauCC_MC_Plus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Plus1",systname.c_str(),ie+1),"",np,p,nr,r));
        NuTauCC_MC_Minus1Sigma[systname].push_back(new TH2D(Form("%s_%i_NuTauCC_MC_Minus1",systname.c_str(),ie+1),"",np,p,nr,r));
      }
      
      delete hreg;

      return;
    }
    
    for(j=0;j<app.size();j++)
    {
      nfiles=0;
      for(ie=0;ie<3;ie++)
      {
        if(gSystem->AccessPathName(gSystem->ExpandPathName(PreFiles[ie].c_str())))//if file doesn't exist
        {
          continue;
        }
        
        fpre = new TFile(gSystem->ExpandPathName(PreFiles[ie].c_str()),"READ");
        treepre = (TTree*)fpre->Get("paramtree");
        treepre->SetBranchAddress("nearPOT",&npot_pre);
        treepre->SetBranchAddress("farPOT",&fpot_pre);
        treepre->GetEntry(0);

        nPOTFar=Extrap[ie]->GetFarPOT();
	nPOTNear=Extrap[ie]->GetNearPOT();
        
        name = string(Background::AsString(app[j])) + "_" + std + "_Presel/FD_RecoVs" + Extrap[0]->GetPID();
        prestd = (TH2D*)fpre->Get(name.c_str());
        prestd->Scale(nPOTFar/fpot_pre);
        MBH.Rebin2DHist(prestd,np,p,nr,r);
        if(nfiles==0) totalprestd = (TH2D*)prestd->Clone();
        else totalprestd->Add(prestd);
        
        name = string(Background::AsString(app[j])) + "_" + plf + "_Presel/FD_RecoVs" + Extrap[0]->GetPID();
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
        if(nfiles==0) totalprepl = (TH2D*)prepl->Clone();
        else totalprepl->Add(prepl);
        
        name = string(Background::AsString(app[j])) + "_" + mnf + "_Presel/FD_RecoVs" + Extrap[0]->GetPID();
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
        if(nfiles==0) totalpremn = (TH2D*)premn->Clone();
        else totalpremn->Add(premn);
        
        nfiles++;
      }
    
      //Regularise 
      totalprestd->Add(hreg);
      totalprepl->Add(hreg);
      totalpremn->Add(hreg);
      //-----------------------
    
      //      totalprepl->Add(totalprestd,-1.);
      totalprepl->Add(totalpremn,-1.);   // SG
      totalprepl->Scale(0.5);            // SG
      totalprepl->Divide(totalprestd);
      totalpremn->Add(totalprestd,-1.);
      totalpremn->Divide(totalprestd);
      
      if(plf==mnf && pln==mnn)
      {
        totalpremn->Scale(-1.);
      }
      
      for(ie=0;ie<Extrap.size();ie++)
      {
        if(app[j]==Background::kNueCC)
        {
          NueCC_MC_Plus1Sigma[systname].push_back((TH2D*)totalprepl->Clone(Form("%s_%i_NueCC_MC_Plus1",systname.c_str(),ie+1)));
          NueCC_MC_Minus1Sigma[systname].push_back((TH2D*)totalpremn->Clone(Form("%s_%i_NueCC_MC_Minus1",systname.c_str(),ie+1)));
        }
        else if(app[j]==Background::kNuTauCC)
        {
          NuTauCC_MC_Plus1Sigma[systname].push_back((TH2D*)totalprepl->Clone(Form("%s_%i_NuTauCC_MC_Plus1",systname.c_str(),ie+1)));
          NuTauCC_MC_Minus1Sigma[systname].push_back((TH2D*)totalpremn->Clone(Form("%s_%i_NuTauCC_MC_Minus1",systname.c_str(),ie+1)));
        }
      }
    }
  }
  delete hreg;
  
  return;
}
void ErrorCalc::ReadSpecialFiles(int n)
{
  string systname = SpecialSyst_SystName.at(n);
  TH2D *h;
  unsigned int ie;
  int i,nr,np,nt;
  double *r,*p,*t;
  nr = Extrap[0]->GetNReco();
  np = Extrap[0]->GetNPID();
  nt = Extrap[0]->GetNTrue();
  r = new double[nr+1];
  p = new double[np+1];
  t = new double[nt+1];
  for(i=0;i<nr+1;i++)
  {
    r[i] = Extrap[0]->Pred_TotalBkgd->GetYaxis()->GetBinLowEdge(i+1);
  }
  for(i=0;i<np+1;i++)
  {
    p[i] = Extrap[0]->Pred_TotalBkgd->GetXaxis()->GetBinLowEdge(i+1);
  }
  for(i=0;i<nt+1;i++)
  {
    t[i] = Extrap[0]->Pred_3D[Extrapolate2D::qNuMuToNuTau]->GetZaxis()->GetBinLowEdge(i+1);
  }
  
  //create empty histograms for all shifts
  for(ie=0;ie<Extrap.size();ie++)
  {
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
    h = (TH2D*)Extrap[0]->Pred_TotalBkgd->Clone("htemp");
    h->Reset();
    for(int ip=0;ip<Extrap[0]->GetNPID();ip++)
    {
      for(int ir=0;ir<Extrap[0]->GetNReco();ir++)
      {
        h->SetBinContent(ip+1,ir+1,SpecialSyst_Value.at(n));
      }
    }
    
    for(unsigned int ie=0;ie<Extrap.size();ie++)
    {
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
    if(gSystem->AccessPathName(gSystem->ExpandPathName(Form("%s/%s_%i.root",SysFileDir.c_str(),filetag.c_str(),Extrap[ie]->GetRunPeriod()))))
    {
      cout<<"Warning: "<<Form("%s/%s_%i.root",SysFileDir.c_str(),filetag.c_str(),Extrap[ie]->GetRunPeriod())<<" doesn't exist."<<endl;
      continue;
    }
    
    f = new TFile(gSystem->ExpandPathName(Form("%s/%s_%i.root",SysFileDir.c_str(),filetag.c_str(),Extrap[ie]->GetRunPeriod())),"READ");
  
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
void ErrorCalc::ReadCovarianceFiles(int n)
{
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
  
  int nbins = Extrap[0]->GetNPID()*Extrap[0]->GetNReco();
  
  cout << " ErrorCalc::ReadCovarianceFiles -  file - hist  " << file << " " << hist << endl; 
  cout << " ErrorCalc::ReadCovarianceFiles -  nbins " <<  nbins << " " << Extrap[0]->GetNPID() << " " << Extrap[0]->GetNReco() << endl;
  if(h->GetNbinsX() != nbins)
  {
    cout<<"Error in ReadCovarianceFiles(): # bins from matrix in file doesn't match extrapolation bins!"<<endl;
    cout << "h->GetNbinsX() -  nbins " << h->GetNbinsX() << " " <<  nbins << endl;
    return;
  }
  
  ExtraCovariance[systname] = h;
  
  return;
}
void ErrorCalc::CalculateFNExtrapError_SingleBin()
{
  Initialize();
  if(!Init) return;
  
  InitializeSys();
  
  string systname;
  
  unsigned int ie;
  TH2D *ch1,*ch2,*ch3;
  double bksum=0,pl=0,mn=0,bnue=0,bpl=0,bmn=0;
  double ncpl=0,ncmn=0,ccpl=0,ccmn=0;
  double nc=0,cc=0;
  double psum=0.,nsum=0.;
  
  std::map<string, vector<TH2D*> >::iterator systiter = FN_NC_Plus1Sigma.begin();
  std::map<string, vector<TH2D*> >::iterator last  = FN_NC_Plus1Sigma.end();
  
  cout<<"\\begin{table}[htp]"<<endl;
  cout<<"\\caption{}"<<endl;
  cout<<"\\centering"<<endl;
  cout<<"\\begin{tabular}{c|cc|cc|cc|cc}"<<endl;
  cout<<"Syst.&NC&&NuMuCC&&BNueCC&&Total&\\\\"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  
  while(systiter!=last)
  {
    systname = systiter->first;
    
    pl=0;
    mn=0;
    bpl=0;
    bmn=0;
    ncpl=0;
    ncmn=0;
    ccpl=0;
    ccmn=0;
    bksum=0;
    bnue=0;
    cc=0;
    nc=0;
    for(ie=0;ie<Extrap.size();ie++)
    {
      ch1 = (TH2D*)Extrap[ie]->Pred[Background::kNC]->Clone();
      ch1->Multiply(FN_NC_Plus1Sigma[systname][ie]);
      ncpl+=ch1->Integral();
      ch2 = (TH2D*)Extrap[ie]->Pred[Background::kNuMuCC]->Clone();
      ch2->Multiply(FN_NuMuCC_Plus1Sigma[systname][ie]);
      ccpl+=ch2->Integral();  cout <<" ccpl " << ch2->Integral() << endl;
      ch1->Add(ch2);
      if(Extrap[ie]->GetFNforBeamNue())
      {
	ch3 = (TH2D*)Extrap[ie]->Pred[Background::kBNueCC]->Clone();
	ch3->Multiply(FN_BNueCC_Plus1Sigma[systname][ie]);
	bpl+=ch3->Integral();
	ch1->Add(ch3);
      }
      pl += ch1->Integral();
      
      ch1 = (TH2D*)Extrap[ie]->Pred[Background::kNC]->Clone();
      ch1->Multiply(FN_NC_Minus1Sigma[systname][ie]);
      ncmn+=ch1->Integral();
      ch2 = (TH2D*)Extrap[ie]->Pred[Background::kNuMuCC]->Clone();
      ch2->Multiply(FN_NuMuCC_Minus1Sigma[systname][ie]);
      ccmn+=ch2->Integral();
      ch1->Add(ch2);
      if(Extrap[ie]->GetFNforBeamNue())
      {
	ch3 = (TH2D*)Extrap[ie]->Pred[Background::kBNueCC]->Clone();
	ch3->Multiply(FN_BNueCC_Minus1Sigma[systname][ie]);
	bmn+=ch3->Integral();
	ch1->Add(ch3);
      }
      mn += ch1->Integral();
      
      bksum+=(Extrap[ie]->Pred[Background::kNC]->Integral()+Extrap[ie]->Pred[Background::kNuMuCC]->Integral());
      cc+=Extrap[ie]->Pred[Background::kNuMuCC]->Integral();
      nc+=Extrap[ie]->Pred[Background::kNC]->Integral();
      if(Extrap[ie]->GetFNforBeamNue())
      {
	bksum+=Extrap[ie]->Pred[Background::kBNueCC]->Integral();
	bnue+=Extrap[ie]->Pred[Background::kBNueCC]->Integral();
      }
    }
    
    if(Extrap[0]->GetFNforBeamNue())
    {

      cout<<systname<<"&";
      cout<<100.*ncpl/nc<<"\\% &";
      if(TMath::Abs(ncpl+ncmn)<1e-10) cout<<"-&";
      else cout<<100.*ncmn/nc<<"\\% &";
      cout<< "==>> "<<100.*ccpl/cc<<"\\% &";
      if(TMath::Abs(ccpl+ccmn)<1e-10) cout<<"-&";
      else cout <<100.*ccmn/cc<<"\\% &";
      cout<<100.*bpl/bnue<<"\\% &";
      if(TMath::Abs(bpl+bmn)<1e-10) cout<<"-&";
      else cout<<100.*bmn/bnue<<"\\% &";
      cout<<"==>>>>"<<100.*pl/bksum<<"\\% &";
      if(TMath::Abs(pl+mn)<1e-10) cout<<"- \\\\"<<endl;
      else cout<<100.*mn/bksum<<"\\% \\\\"<<endl;
    }
    else
    {
      cout<<systname<<"&";
      cout<<100.*ncpl/nc<<"\\% &";
      if(TMath::Abs(ncpl+ncmn)<1e-10) cout<<"-&";
      else cout<<100.*ncmn/nc<<"\\% &";
      cout<<100.*ccpl/cc<<"\\% &";
      if(TMath::Abs(ccpl+ccmn)<1e-10) cout<<"-&";
      else cout<<100.*ccmn/cc<<"\\% &";
      cout<<100.*pl/bksum<<"\\% &";
      if(TMath::Abs(pl+mn)<1e-10) cout<<"- \\\\"<<endl;
      else cout<<100.*mn/bksum<<"\\% \\\\"<<endl;
    }
    
    if((pl>0 && mn>0) || (pl<0 && mn<0))
    {
      if(TMath::Abs(pl)>TMath::Abs(mn))
      {
        mn = -1.*pl;
      }
      else
      {
        pl = -1.*mn;
      }
    }
    
    if(pl>0) psum+=(pl*pl);
    else nsum+=(pl*pl);
    if(mn>0) psum+=(mn*mn);
    else nsum+=(mn*mn);
    
    systiter++;
  }
  
  psum = sqrt(psum);
  nsum = -1.*sqrt(nsum);
  
  cout<<"\\hline"<<endl;
  cout<<"Total Extrap&&&&&&&"<<100.*psum/bksum<<"\\% &"<<100.*nsum/bksum<<"\\% \\\\"<<endl;
  cout<<"\\end{tabular}"<<endl;
  cout<<"\\end{table}"<<endl;
  
  return;
}
void ErrorCalc::CalculateAppearanceExtrapError_SingleBin(Background::Background_t bg)
{
  if(bg!=Background::kNueCC && bg!=Background::kNuTauCC)
  {
    cout<<"ErrorCalc::CalculateAppearanceExtrapError() works only for NueCC or NuTauCC.  Quitting..."<<endl;
    return;
  }
  
  Initialize();
  if(!Init) return;
  
  InitializeSys();
  
//   TH3D *h;
  TH2D *g;
  double pl=0,mn=0,psum=0,nsum=0,sum=0;
  unsigned int ie;
  string systname;
  
  if(bg==Background::kNuTauCC) cout<<"Tau Table"<<endl;
  if(bg==Background::kNueCC) cout<<"Nue Table"<<endl;
  
  cout<<"\\begin{table}[htp]"<<endl;
  cout<<"\\caption{}"<<endl;
  cout<<"\\centering"<<endl;
  cout<<"\\begin{tabular}{c|cc}"<<endl;
  cout<<"Syst.&"<<string(Background::AsString(bg))<<"&\\\\"<<endl;
  cout<<"\\hline"<<endl;
  cout<<"\\hline"<<endl;
  
  psum=0;nsum=0;
  
  std::map<string, vector<TH2D*> >::iterator sysiter2;
  std::map<string, vector<TH2D*> >::iterator last2;
  if(bg==Background::kNueCC)
  {
    sysiter2 = NueCC_MC_Plus1Sigma.begin();
    last2  = NueCC_MC_Plus1Sigma.end();
  }
  else
  {
    sysiter2 = NuTauCC_MC_Plus1Sigma.begin();
    last2  = NuTauCC_MC_Plus1Sigma.end();
  }
  
  while(sysiter2!=last2)
  {
    systname = sysiter2->first;
  
    pl=0;
    mn=0;
    sum=0;
    for(ie=0;ie<Extrap.size();ie++)
    {
      if(bg==Background::kNuTauCC)
      {
        g = (TH2D*)Extrap[ie]->Pred[Background::kNuTauCC]->Clone();
        g->Multiply(NuTauCC_MC_Plus1Sigma[systname][ie]);
        pl += g->Integral();
      
        g = (TH2D*)Extrap[ie]->Pred[Background::kNuTauCC]->Clone();
        g->Multiply(NuTauCC_MC_Minus1Sigma[systname][ie]);
        mn += g->Integral();
      
        sum += Extrap[ie]->Pred[Background::kNuTauCC]->Integral();
      }
      if(bg==Background::kNueCC)
      {
        g = (TH2D*)Extrap[ie]->Pred[Background::kNueCC]->Clone();
        g->Multiply(NueCC_MC_Plus1Sigma[systname][ie]);
        pl += g->Integral();
      
        g = (TH2D*)Extrap[ie]->Pred[Background::kNueCC]->Clone();
        g->Multiply(NueCC_MC_Minus1Sigma[systname][ie]);
        mn += g->Integral();
      
        sum += Extrap[ie]->Pred[Background::kNueCC]->Integral();
      }
    }
    
    cout<<systname<<"&";
    cout<<100.*pl/sum<<"\\% &";
    if(TMath::Abs(pl+mn)<1e-10) cout<<"- \\\\"<<endl;
    else cout<<100.*mn/sum<<"\\% \\\\"<<endl;
    
    if((pl>0 && mn>0) || (pl<0 && mn<0))
    {
      if(TMath::Abs(pl)>TMath::Abs(mn))
      {
        mn = -1.*pl;
      }
      else
      {
        pl = -1.*mn;
      }
    }
    
    if(pl>0) psum+=(pl*pl);
    else nsum+=(pl*pl);
    if(mn>0) psum+=(mn*mn);
    else nsum+=(mn*mn);
    
    sysiter2++;
  }
  
  //no longer needed due to change in how Pred_CC_Fid_Plus/Minus1Sigma is used
//   std::map<string, vector<TH1D*> >::iterator sysiter3 = Pred_CC_Fid_Plus1Sigma.begin();
//   std::map<string, vector<TH1D*> >::iterator last3 = Pred_CC_Fid_Plus1Sigma.end();
//   
//   int nr,np,nt,it;
//   
//   while(sysiter3!=last3)
//   {
//     systname = sysiter3->first;
//   
//     pl=0;
//     mn=0;
//     sum=0;
//     for(ie=0;ie<Extrap.size();ie++)
//     {
//       nr = Extrap[ie]->GetNReco();
//       np = Extrap[ie]->GetNPID();
//       nt = Extrap[ie]->GetNTrue();
//       
//       if(bg==Background::kNuTauCC)
//       {
//         h = (TH3D*)Extrap[ie]->Pred_3D[Extrapolate2D::qNuMuToNuTau]->Clone();
//         h->Add(Extrap[ie]->Pred_3D[Extrapolate2D::qBNueToNuTau]);
//         sum += h->Integral();
//       }
//       else
//       {
//         h = (TH3D*)Extrap[ie]->Pred_3D[Extrapolate2D::qNuMuToNue]->Clone();
//         sum += h->Integral();
//       }
//       
//       for(it=0;it<nt;it++)
//       {
//         pl+=h->Integral(1,np,1,nr,it+1,it+1)*Pred_CC_Fid_Plus1Sigma[systname][ie]->GetBinContent(it+1);
//         
//         mn+=h->Integral(1,np,1,nr,it+1,it+1)*Pred_CC_Fid_Minus1Sigma[systname][ie]->GetBinContent(it+1);
//       }
//     }
//     
//     cout<<systname<<"&";
//     cout<<100.*pl/sum<<"\\% &";
//     if(TMath::Abs(pl+mn)<1e-10) cout<<"- \\\\"<<endl;
//     else cout<<100.*mn/sum<<"\\% \\\\"<<endl;
//     
//     if((pl>0 && mn>0) || (pl<0 && mn<0))
//     {
//       if(TMath::Abs(pl)>TMath::Abs(mn))
//       {
//         mn = -1.*pl;
//       }
//       else
//       {
//         pl = -1.*mn;
//       }
//     }
//     
//     if(pl>0) psum+=(pl*pl);
//     else nsum+=(pl*pl);
//     if(mn>0) psum+=(mn*mn);
//     else nsum+=(mn*mn);
//     
//     sysiter3++;
//   }
  
  psum = sqrt(psum);
  nsum = -1.*sqrt(nsum);
  
  cout<<"\\hline"<<endl;
  cout<<"Total Extrap&"<<100.*psum/sum<<"\\% &"<<100.*nsum/sum<<"\\% \\\\"<<endl;
  cout<<"\\end{tabular}"<<endl;
  cout<<"\\end{table}"<<endl;
  
  return;
}
void ErrorCalc::CalculateHOOError_SingleBin()
{
  //this function assumes the error matrix is normalized to the same POT that the denominator of the F/N ratio is (i.e.1e19)
  Initialize();
  if(!Init) return;
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(NDData_infile.c_str())))
  {


    cout << "CalculateHOOError_SingleBin " << endl;
    cout<<"Failure to read ND Data file."<<endl;
    return;
  }
  
  TFile *f = new TFile(gSystem->ExpandPathName(NDData_infile.c_str()),"READ");
  string name = "CovMatrix_"+Extrap[0]->GetNDDataPID();
  TH2D *VN = (TH2D*)f->Get(name.c_str());
  
  int ip,ir;
  unsigned int ie;
  int nr,np;
  int n = VN->GetNbinsX();
  vector<double> A;
  
  for(ie=0;ie<Extrap.size();ie++)
  {
    nr = Extrap[ie]->GetNReco();
    np = Extrap[ie]->GetNPID();
    for(ir=0;ir<nr;ir++)
    {
      for(ip=0;ip<np;ip++)
      {
        A.push_back(Extrap[ie]->FNRatio[Background::kNC]->GetBinContent(ip+1,ir+1));
      }
    }
  }
  
  for(ie=0;ie<Extrap.size();ie++)
  {
    nr = Extrap[ie]->GetNReco();
    np = Extrap[ie]->GetNPID();
    for(ir=0;ir<nr;ir++)
    {
      for(ip=0;ip<np;ip++)
      {
	A.push_back(Extrap[ie]->FNRatio[Background::kNuMuCC]->GetBinContent(ip+1,ir+1));
      }
    }
  }
  
  for(ie=0;ie<Extrap.size();ie++)
  {
    nr = Extrap[ie]->GetNReco();
    np = Extrap[ie]->GetNPID();
    for(ir=0;ir<nr;ir++)
    {
      for(ip=0;ip<np;ip++)
      {
	A.push_back(Extrap[ie]->FNRatio[Background::kBNueCC]->GetBinContent(ip+1,ir+1));
      }
    }
  }
  
  double VF = 0;
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<n;j++)
    {
      VF += (VN->GetBinContent(i+1,j+1)*A[i]*A[j]);
    }
  }
  
  double bksum=0;
  for(ie=0;ie<Extrap.size();ie++)
  {
    bksum+=(Extrap[ie]->Pred[Background::kNC]->Integral()+Extrap[ie]->Pred[Background::kNuMuCC]->Integral()+Extrap[ie]->Pred[Background::kBNueCC]->Integral());
  }
  
  cout<<"Total syst (NC+NuMuCC+BNueCC) due to HOO Input: "<<bksum<<" +- "<<100.*TMath::Sqrt(VF)/bksum<<"%"<<endl;
  
  return;
}
void ErrorCalc::CalculateMRCCError_SingleBin()
{
  //see DocDB- (Josh's syst note from the first analysis), section 6.2 for the explanation of this calculation.  Correlation between the CC and NC components is taken into account.
  
  Initialize();
  if(!Init) return;
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(NDData_infile.c_str())))
  {

    cout << "CalculateMRCCError_SingleBin " << endl;
    cout<<"Failure to read ND Data file."<<endl;
    return;
  }
  
  TFile *f = new TFile(gSystem->ExpandPathName(NDData_infile.c_str()),"READ");
  
  string p = Extrap[0]->GetNDDataPID();
  unsigned int ie;
  TH2D *temp = 0;
  
  TH2D *NDData_NuMuCC_Errors[3];
  TH2D *NDData_BNueCC_Errors[3];
  TH2D *NDData_Total[3];
  TH2D *NDData_NuMuCC[3];
  TH2D *NDData_BNueCC[3];
  for(ie=0;ie<3;ie++)
  {
    NDData_NuMuCC_Errors[ie] = (TH2D*)f->Get(Form("NDData_NuMuCC_Errors_%s_Run%i",p.c_str(),ie+1));
    NDData_BNueCC_Errors[ie] = (TH2D*)f->Get(Form("NDData_BNueCC_Errors_%s_Run%i",p.c_str(),ie+1));
    NDData_Total[ie] = (TH2D*)f->Get(Form("NDData_Total_%s_Run%i",p.c_str(),ie+1));
    NDData_NuMuCC[ie] = (TH2D*)f->Get(Form("NDData_NuMuCC_%s_Run%i",p.c_str(),ie+1));
    NDData_BNueCC[ie] = (TH2D*)f->Get(Form("NDData_BNueCC_%s_Run%i",p.c_str(),ie+1));
  }
  
  double numucc_syst=0;
  double bnuecc_syst=0;
  double total_stat=0;
  double numucc_stat=0;
  double bnuecc_stat=0;
  
  int ip,ir;
  int np,nr;
  double t;
  for(ie=0;ie<Extrap.size();ie++)
  {
    temp = (TH2D*)Extrap[ie]->FNRatio[Background::kNuMuCC]->Clone("temp");
    temp->Add(Extrap[ie]->FNRatio[Background::kNC],-1.0);
    temp->Multiply(NDData_NuMuCC_Errors[ie]);
    numucc_syst += temp->Integral();
    
    temp = (TH2D*)Extrap[ie]->FNRatio[Background::kBNueCC]->Clone("temp");
    temp->Add(Extrap[ie]->FNRatio[Background::kNC],-1.0);
    temp->Multiply(NDData_BNueCC_Errors[ie]);
    bnuecc_syst += temp->Integral();
    
    nr = Extrap[ie]->GetNReco();
    np = Extrap[ie]->GetNPID();
    for(ip=0;ip<np;ip++)
    {
      for(ir=0;ir<nr;ir++)
      {
        total_stat += (Extrap[ie]->FNRatio[Background::kNC]->GetBinContent(ip+1,ir+1)*Extrap[ie]->FNRatio[Background::kNC]->GetBinContent(ip+1,ir+1)*NDData_Total[ie]->GetBinError(ip+1,ir+1)*NDData_Total[ie]->GetBinError(ip+1,ir+1));
        
        t = Extrap[ie]->FNRatio[Background::kNuMuCC]->GetBinContent(ip+1,ir+1) - Extrap[ie]->FNRatio[Background::kNC]->GetBinContent(ip+1,ir+1);
        numucc_stat += (t*t*NDData_NuMuCC[ie]->GetBinError(ip+1,ir+1)*NDData_NuMuCC[ie]->GetBinError(ip+1,ir+1));
        
        t = Extrap[ie]->FNRatio[Background::kBNueCC]->GetBinContent(ip+1,ir+1) - Extrap[ie]->FNRatio[Background::kNC]->GetBinContent(ip+1,ir+1);
        bnuecc_stat += (t*t*NDData_BNueCC[ie]->GetBinError(ip+1,ir+1)*NDData_BNueCC[ie]->GetBinError(ip+1,ir+1));
      }
    }
  }
  
  double totalsyst = sqrt(numucc_syst*numucc_syst + bnuecc_syst*bnuecc_syst + total_stat + numucc_stat + bnuecc_stat);
  
  double bksum=0;
  for(ie=0;ie<Extrap.size();ie++)
  {
    bksum+=(Extrap[ie]->Pred[Background::kNC]->Integral()+Extrap[ie]->Pred[Background::kNuMuCC]->Integral()+Extrap[ie]->Pred[Background::kBNueCC]->Integral());
  }
  
  cout<<"Total syst (NC+NuMuCC+BNueCC) due to MRCC Input: "<<bksum<<" +- "<<100.*totalsyst/bksum<<"%"<<endl;
  
  return;
}
void ErrorCalc::FNRatioError(string ratioerrfilename)
{
  Initialize();
  if(!Init) return;
  
  InitializeSys();
  
  vector<TH2D*> FN_NC_Plus_Total;
  vector<TH2D*> FN_NC_Minus_Total;
  vector<TH2D*> FN_NuMuCC_Plus_Total;
  vector<TH2D*> FN_NuMuCC_Minus_Total;
  vector<TH2D*> FN_BNueCC_Plus_Total;
  vector<TH2D*> FN_BNueCC_Minus_Total;
  
  int i,ip,ir;
  unsigned int ie;
  
  int nr = Extrap[0]->GetNReco();
  int np = Extrap[0]->GetNPID();
  double *r = new double[nr+1];
  double *p = new double[np+1];
  for(i=0;i<nr+1;i++)
  {
    r[i] = Extrap[0]->Pred_TotalBkgd->GetYaxis()->GetBinLowEdge(i+1);
  }
  for(i=0;i<np+1;i++)
  {
    p[i] = Extrap[0]->Pred_TotalBkgd->GetXaxis()->GetBinLowEdge(i+1);
  }
  
  for(ie=0;ie<Extrap.size();ie++)
  {
    FN_NC_Plus_Total.push_back(new TH2D(Form("FN_NC_Plus_Total_Run%i",ie+1),"",np,p,nr,r));
    FN_NC_Minus_Total.push_back(new TH2D(Form("FN_NC_Minus_Total_Run%i",ie+1),"",np,p,nr,r));
    FN_NuMuCC_Plus_Total.push_back(new TH2D(Form("FN_NuMuCC_Plus_Total_Run%i",ie+1),"",np,p,nr,r));
    FN_NuMuCC_Minus_Total.push_back(new TH2D(Form("FN_NuMuCC_Minus_Total_Run%i",ie+1),"",np,p,nr,r));
    FN_BNueCC_Plus_Total.push_back(new TH2D(Form("FN_BNueCC_Plus_Total_Run%i",ie+1),"",np,p,nr,r));
    FN_BNueCC_Minus_Total.push_back(new TH2D(Form("FN_BNueCC_Minus_Total_Run%i",ie+1),"",np,p,nr,r));
  }
  
  string systname;
  double pl,mn,temp;
  
  std::map<string, vector<TH2D*> >::iterator systiter = FN_NC_Plus1Sigma.begin();
  std::map<string, vector<TH2D*> >::iterator last  = FN_NC_Plus1Sigma.end();
  
  for(ie=0;ie<Extrap.size();ie++)
  {
    for(ip=0;ip<np;ip++)
    {
      for(ir=0;ir<nr;ir++)
      {
	pl=0;
	mn=0;
	systiter = FN_NC_Plus1Sigma.begin();
	while(systiter!=last)
        {
          systname = systiter->first;
	  
	  temp = FN_NC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1);
	  if(temp>0) pl+=temp*temp;
	  else mn+=temp*temp;
	  
	  temp = FN_NC_Minus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1);
	  if(temp>0) pl+=temp*temp;
	  else mn+=temp*temp;
	  
	  systiter++;
	}
	FN_NC_Plus_Total[ie]->SetBinContent(ip+1,ir+1,sqrt(pl));
	FN_NC_Minus_Total[ie]->SetBinContent(ip+1,ir+1,-1.*sqrt(mn));
	
	pl=0;
	mn=0;
	systiter = FN_NC_Plus1Sigma.begin();
	while(systiter!=last)
        {
          systname = systiter->first;
	  
	  temp = FN_NuMuCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1);
	  if(temp>0) pl+=temp*temp;
	  else mn+=temp*temp;
	  
	  temp = FN_NuMuCC_Minus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1);
	  if(temp>0) pl+=temp*temp;
	  else mn+=temp*temp;
	  
	  systiter++;
	}
	FN_NuMuCC_Plus_Total[ie]->SetBinContent(ip+1,ir+1,sqrt(pl));
	FN_NuMuCC_Minus_Total[ie]->SetBinContent(ip+1,ir+1,-1.*sqrt(mn));
	
	pl=0;
	mn=0;
	systiter = FN_NC_Plus1Sigma.begin();
	while(systiter!=last)
        {
          systname = systiter->first;
	  
	  temp = FN_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1);
	  if(temp>0) pl+=temp*temp;
	  else mn+=temp*temp;
	  
	  temp = FN_BNueCC_Minus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1);
	  if(temp>0) pl+=temp*temp;
	  else mn+=temp*temp;
	  
	  systiter++;
	}
	FN_BNueCC_Plus_Total[ie]->SetBinContent(ip+1,ir+1,sqrt(pl));
	FN_BNueCC_Minus_Total[ie]->SetBinContent(ip+1,ir+1,-1.*sqrt(mn));
      }
    }
  }
  
  TFile *f = new TFile(gSystem->ExpandPathName(ratioerrfilename.c_str()),"RECREATE");
  for(ie=0;ie<Extrap.size();ie++)
  {
    FN_NC_Plus_Total[ie]->Write();
    FN_NC_Minus_Total[ie]->Write();
    FN_NuMuCC_Plus_Total[ie]->Write();
    FN_NuMuCC_Minus_Total[ie]->Write();
    FN_BNueCC_Plus_Total[ie]->Write();
    FN_BNueCC_Minus_Total[ie]->Write();
  }
  f->Close();
  
  return;
}
void ErrorCalc::NDSystematics(string nderrfilename)
{
  Initialize();
  if(!Init) return;
  
  InitializeSys();
  
  vector<TH2D*> ND_NC_Plus_Total;
  vector<TH2D*> ND_NC_Minus_Total;
  vector<TH2D*> ND_NuMuCC_Plus_Total;
  vector<TH2D*> ND_NuMuCC_Minus_Total;
  vector<TH2D*> ND_BNueCC_Plus_Total;
  vector<TH2D*> ND_BNueCC_Minus_Total;
  vector<TH2D*> ND_Total_Plus_Total;
  vector<TH2D*> ND_Total_Minus_Total;
  TH2D *ND_Total_Plus_Total_AllRuns;
  TH2D *ND_Total_Minus_Total_AllRuns;
  
  int i,ip,ir;
  unsigned int ie;
  
  int nr = Extrap[0]->GetNReco();
  int np = Extrap[0]->GetNPID();
  double *r = new double[nr+1];
  double *p = new double[np+1];
  for(i=0;i<nr+1;i++)
  {
    r[i] = Extrap[0]->Pred_TotalBkgd->GetYaxis()->GetBinLowEdge(i+1);
  }
  for(i=0;i<np+1;i++)
  {
    p[i] = Extrap[0]->Pred_TotalBkgd->GetXaxis()->GetBinLowEdge(i+1);
  }
  
  for(ie=0;ie<Extrap.size();ie++)
  {
    ND_NC_Plus_Total.push_back(new TH2D(Form("ND_NC_Plus_Total_Run%i",ie+1),"",np,p,nr,r));
    ND_NC_Minus_Total.push_back(new TH2D(Form("ND_NC_Minus_Total_Run%i",ie+1),"",np,p,nr,r));
    ND_NuMuCC_Plus_Total.push_back(new TH2D(Form("ND_NuMuCC_Plus_Total_Run%i",ie+1),"",np,p,nr,r));
    ND_NuMuCC_Minus_Total.push_back(new TH2D(Form("ND_NuMuCC_Minus_Total_Run%i",ie+1),"",np,p,nr,r));
    ND_BNueCC_Plus_Total.push_back(new TH2D(Form("ND_BNueCC_Plus_Total_Run%i",ie+1),"",np,p,nr,r));
    ND_BNueCC_Minus_Total.push_back(new TH2D(Form("ND_BNueCC_Minus_Total_Run%i",ie+1),"",np,p,nr,r));
    ND_Total_Plus_Total.push_back(new TH2D(Form("ND_Total_Plus_Total_Run%i",ie+1),"",np,p,nr,r));
    ND_Total_Minus_Total.push_back(new TH2D(Form("ND_Total_Minus_Total_Run%i",ie+1),"",np,p,nr,r));
  }
  ND_Total_Plus_Total_AllRuns = new TH2D("ND_Total_Plus_Total_AllRuns","",np,p,nr,r);
  ND_Total_Minus_Total_AllRuns = new TH2D("ND_Total_Minus_Total_AllRuns","",np,p,nr,r);
  
  string systname;
  double pl=0,mn=0,temp=0;
  double bksum=0,cc=0,nc=0,bnue=0;
  double ccpl=0,ncpl=0,bpl=0;
  double ccmn=0,ncmn=0,bmn=0;
  double psum=0,nsum=0;
  double ccpsum=0,ccnsum=0;
  double ncpsum=0,ncnsum=0;
  double bpsum=0,bnsum=0;
  TH2D *ch1,*ch2,*ch3;
  
  std::map<string, vector<TH2D*> >::iterator systiter = ND_NC_Plus1Sigma.begin();
  std::map<string, vector<TH2D*> >::iterator last  = ND_NC_Plus1Sigma.end();
  
  for(ie=0;ie<Extrap.size();ie++)
  {
    for(ip=0;ip<np;ip++)
    {
      for(ir=0;ir<nr;ir++)
      {
        pl=0;
        mn=0;
        systiter = ND_NC_Plus1Sigma.begin();
        while(systiter!=last)
        {
          systname = systiter->first;
          
          temp = ND_NC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1);
          if(temp>0) pl+=temp*temp;
          else mn+=temp*temp;
          
          temp = ND_NC_Minus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1);
          if(temp>0) pl+=temp*temp;
          else mn+=temp*temp;
          
          systiter++;
        }
        ND_NC_Plus_Total[ie]->SetBinContent(ip+1,ir+1,sqrt(pl));
        ND_NC_Minus_Total[ie]->SetBinContent(ip+1,ir+1,-1.*sqrt(mn));
        
        pl=0;
        mn=0;
        systiter = ND_NC_Plus1Sigma.begin();
        while(systiter!=last)
        {
          systname = systiter->first;
          
          temp = ND_NuMuCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1);
          if(temp>0) pl+=temp*temp;
          else mn+=temp*temp;
          
          temp = ND_NuMuCC_Minus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1);
          if(temp>0) pl+=temp*temp;
          else mn+=temp*temp;
          
          systiter++;
        }
        ND_NuMuCC_Plus_Total[ie]->SetBinContent(ip+1,ir+1,sqrt(pl));
        ND_NuMuCC_Minus_Total[ie]->SetBinContent(ip+1,ir+1,-1.*sqrt(mn));
        
        pl=0;
        mn=0;
        systiter = ND_NC_Plus1Sigma.begin();
        while(systiter!=last)
        {
          systname = systiter->first;
          
          temp = ND_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1);
          if(temp>0) pl+=temp*temp;
          else mn+=temp*temp;
          
          temp = ND_BNueCC_Minus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1);
          if(temp>0) pl+=temp*temp;
          else mn+=temp*temp;
          
          systiter++;
        }
        ND_BNueCC_Plus_Total[ie]->SetBinContent(ip+1,ir+1,sqrt(pl));
        ND_BNueCC_Minus_Total[ie]->SetBinContent(ip+1,ir+1,-1.*sqrt(mn));
        
        pl=0;
        mn=0;
        systiter = ND_NC_Plus1Sigma.begin();
        while(systiter!=last)
        {
          systname = systiter->first;
          
          temp = ND_NC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->NDMC[Background::kNC]->GetBinContent(ip+1,ir+1) + ND_NuMuCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->NDMC[Background::kNuMuCC]->GetBinContent(ip+1,ir+1) + ND_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->NDMC[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
          if(temp>0) pl+=temp*temp;
          else mn+=temp*temp;
          
          temp = ND_NC_Minus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->NDMC[Background::kNC]->GetBinContent(ip+1,ir+1) + ND_NuMuCC_Minus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->NDMC[Background::kNuMuCC]->GetBinContent(ip+1,ir+1) + ND_BNueCC_Minus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->NDMC[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
          if(temp>0) pl+=temp*temp;
          else mn+=temp*temp;
          
          systiter++;
        }
        temp = Extrap[ie]->NDMC[Background::kNC]->GetBinContent(ip+1,ir+1) + Extrap[ie]->NDMC[Background::kNuMuCC]->GetBinContent(ip+1,ir+1) + Extrap[ie]->NDMC[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
        
        ND_Total_Plus_Total[ie]->SetBinContent(ip+1,ir+1,sqrt(pl)/temp);
        ND_Total_Minus_Total[ie]->SetBinContent(ip+1,ir+1,-1.*sqrt(mn)/temp);
      }
    }
  }
  
  for(ip=0;ip<np;ip++)
  {
    for(ir=0;ir<nr;ir++)
    {
      pl=0;
      mn=0;
      systiter = ND_NC_Plus1Sigma.begin();
      while(systiter!=last)
      {
        systname = systiter->first;
        temp=0;
        for(ie=0;ie<Extrap.size();ie++)
        {
          temp += (ND_NC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->NDMC[Background::kNC]->GetBinContent(ip+1,ir+1) + ND_NuMuCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->NDMC[Background::kNuMuCC]->GetBinContent(ip+1,ir+1) + ND_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->NDMC[Background::kBNueCC]->GetBinContent(ip+1,ir+1));
        }
        if(temp>0) pl+=temp*temp;
        else mn+=temp*temp;
        
        temp=0;
        for(ie=0;ie<Extrap.size();ie++)
        { 
          temp += (ND_NC_Minus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->NDMC[Background::kNC]->GetBinContent(ip+1,ir+1) + ND_NuMuCC_Minus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->NDMC[Background::kNuMuCC]->GetBinContent(ip+1,ir+1) + ND_BNueCC_Minus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->NDMC[Background::kBNueCC]->GetBinContent(ip+1,ir+1));
        }
        if(temp>0) pl+=temp*temp;
        else mn+=temp*temp;
        
        systiter++;
      }
      
      temp=0;
      for(ie=0;ie<Extrap.size();ie++)
      {
        temp += (Extrap[ie]->NDMC[Background::kNC]->GetBinContent(ip+1,ir+1) + Extrap[ie]->NDMC[Background::kNuMuCC]->GetBinContent(ip+1,ir+1) + Extrap[ie]->NDMC[Background::kBNueCC]->GetBinContent(ip+1,ir+1));
      }
      ND_Total_Plus_Total_AllRuns->SetBinContent(ip+1,ir+1,sqrt(pl)/temp);
      ND_Total_Minus_Total_AllRuns->SetBinContent(ip+1,ir+1,-1.*sqrt(mn)/temp);
    }
  }
  
  TFile *f = new TFile(gSystem->ExpandPathName(nderrfilename.c_str()),"RECREATE");
  for(ie=0;ie<Extrap.size();ie++)
  {
    ND_NC_Plus_Total[ie]->Write();
    ND_NC_Minus_Total[ie]->Write();
    ND_NuMuCC_Plus_Total[ie]->Write();
    ND_NuMuCC_Minus_Total[ie]->Write();
    ND_BNueCC_Plus_Total[ie]->Write();
    ND_BNueCC_Minus_Total[ie]->Write();
    ND_Total_Plus_Total[ie]->Write();
    ND_Total_Minus_Total[ie]->Write();
  }
  ND_Total_Plus_Total_AllRuns->Write();
  ND_Total_Minus_Total_AllRuns->Write();
  
  f->Close();
  
  cout<<"ND Syst Table"<<endl;
  
  ofstream myfile("Tables/ndtable.tex");
  myfile.precision(2);
  myfile<<fixed;
  
  myfile<<"\\begin{table}[htp]"<<endl;
  myfile<<"\\caption{}"<<endl;
  myfile<<"\\centering"<<endl;
  myfile<<"\\begin{tabular}{c|cc|cc|cc|cc}"<<endl;
  myfile<<"Syst.&NC&&NuMuCC&&BNueCC&&Total&\\\\"<<endl;
  myfile<<"\\hline"<<endl;
  myfile<<"\\hline"<<endl;
  
  systiter = ND_NC_Plus1Sigma.begin();
  psum=0;
  nsum=0;
  ccpsum=0;
  ccnsum=0;
  ncpsum=0;
  ncnsum=0;
  bpsum=0;
  bnsum=0;
  while(systiter!=last)
  {
    systname = systiter->first;
    
    pl=0;
    mn=0;
    bpl=0;
    bmn=0;
    ncpl=0;
    ncmn=0;
    ccpl=0;
    ccmn=0;
    bksum=0;
    bnue=0;
    cc=0;
    nc=0;
    for(ie=0;ie<Extrap.size();ie++)
    {
      ch1 = (TH2D*)Extrap[ie]->NDMC[Background::kNC]->Clone();
      ch1->Multiply(ND_NC_Plus1Sigma[systname][ie]);
      ncpl+=ch1->Integral();
      ch2 = (TH2D*)Extrap[ie]->NDMC[Background::kNuMuCC]->Clone();
      ch2->Multiply(ND_NuMuCC_Plus1Sigma[systname][ie]);
      ccpl+=ch2->Integral();
      ch1->Add(ch2);
      ch3 = (TH2D*)Extrap[ie]->NDMC[Background::kBNueCC]->Clone();
      ch3->Multiply(ND_BNueCC_Plus1Sigma[systname][ie]);
      bpl+=ch3->Integral();
      ch1->Add(ch3);
      pl += ch1->Integral();
      
      ch1 = (TH2D*)Extrap[ie]->NDMC[Background::kNC]->Clone();
      ch1->Multiply(ND_NC_Minus1Sigma[systname][ie]);
      ncmn+=ch1->Integral();
      ch2 = (TH2D*)Extrap[ie]->NDMC[Background::kNuMuCC]->Clone();
      ch2->Multiply(ND_NuMuCC_Minus1Sigma[systname][ie]);
      ccmn+=ch2->Integral();
      ch1->Add(ch2);
      ch3 = (TH2D*)Extrap[ie]->NDMC[Background::kBNueCC]->Clone();
      ch3->Multiply(ND_BNueCC_Minus1Sigma[systname][ie]);
      bmn+=ch3->Integral();
      ch1->Add(ch3);
      mn += ch1->Integral();
      
      bksum+=(Extrap[ie]->NDMC[Background::kNuMuCC]->Integral() + Extrap[ie]->NDMC[Background::kNC]->Integral() + Extrap[ie]->NDMC[Background::kBNueCC]->Integral());
      cc+=Extrap[ie]->NDMC[Background::kNuMuCC]->Integral();
      nc+=Extrap[ie]->NDMC[Background::kNC]->Integral();
      bnue+=Extrap[ie]->NDMC[Background::kBNueCC]->Integral();
    }
    
    myfile<<systname<<"&";
    myfile<<100.*ncpl/nc<<"\\% &";
    if(TMath::Abs(ncpl+ncmn)<1e-20) myfile<<"-&";
    else myfile<<100.*ncmn/nc<<"\\% &";
    myfile<<100.*ccpl/cc<<"\\% &";
    if(TMath::Abs(ccpl+ccmn)<1e-20) myfile<<"-&";
    else myfile<<100.*ccmn/cc<<"\\% &";
    myfile<<100.*bpl/bnue<<"\\% &";
    if(TMath::Abs(bpl+bmn)<1e-20) myfile<<"-&";
    else myfile<<100.*bmn/bnue<<"\\% &";
    myfile<<100.*pl/bksum<<"\\% &";
    if(TMath::Abs(pl+mn)<1e-20) myfile<<"- \\\\"<<endl;
    else myfile<<100.*mn/bksum<<"\\% \\\\"<<endl;
    
    if((pl>0 && mn>0) || (pl<0 && mn<0))
    {
      if(TMath::Abs(pl)>TMath::Abs(mn))
      {
        mn = -1.*pl;
      }
      else
      {
        pl = -1.*mn;
      }
    }
    if((ccpl>0 && ccmn>0) || (ccpl<0 && ccmn<0))
    {
      if(TMath::Abs(ccpl)>TMath::Abs(ccmn))
      {
        ccmn = -1.*ccpl;
      }
      else
      {
        ccpl = -1.*ccmn;
      }
    }
    if((ncpl>0 && ncmn>0) || (ncpl<0 && ncmn<0))
    {
      if(TMath::Abs(ncpl)>TMath::Abs(ncmn))
      {
        ncmn = -1.*ncpl;
      }
      else
      {
        ncpl = -1.*ncmn;
      }
    }
    if((bpl>0 && bmn>0) || (bpl<0 && bmn<0))
    {
      if(TMath::Abs(bpl)>TMath::Abs(bmn))
      {
        bmn = -1.*bpl;
      }
      else
      {
        bpl = -1.*bmn;
      }
    }
    
    if(pl>0) psum+=(pl*pl);
    else nsum+=(pl*pl);
    if(mn>0) psum+=(mn*mn);
    else nsum+=(mn*mn);
    
    if(ccpl>0) ccpsum+=(ccpl*ccpl);
    else ccnsum+=(ccpl*ccpl);
    if(ccmn>0) ccpsum+=(ccmn*ccmn);
    else ccnsum+=(ccmn*ccmn);
    
    if(ncpl>0) ncpsum+=(ncpl*ncpl);
    else ncnsum+=(ncpl*ncpl);
    if(ncmn>0) ncpsum+=(ncmn*ncmn);
    else ncnsum+=(ncmn*ncmn);
    
    if(bpl>0) bpsum+=(bpl*bpl);
    else bnsum+=(bpl*bpl);
    if(bmn>0) bpsum+=(bmn*bmn);
    else bnsum+=(bmn*bmn);
    
    systiter++;
  }
  
  psum = sqrt(psum);
  nsum = -1.*sqrt(nsum);
  
  ccpsum = sqrt(ccpsum);
  ccnsum = -1.*sqrt(ccnsum);
  
  ncpsum = sqrt(ncpsum);
  ncnsum = -1.*sqrt(ncnsum);
  
  bpsum = sqrt(bpsum);
  bnsum = -1.*sqrt(bnsum);
  
  myfile<<"\\hline"<<endl;
  myfile<<"Total Extrap&"<<100.*ncpsum/nc<<"&"<<100.*ncnsum/nc<<"&"<<100.*ccpsum/cc<<"&"<<100.*ccnsum/cc<<"&"<<100.*bpsum/bnue<<"&"<<100.*bnsum/bnue<<"&"<<100.*psum/bksum<<"\\% &"<<100.*nsum/bksum<<"\\% \\\\"<<endl;
  myfile<<"\\end{tabular}"<<endl;
  myfile<<"\\end{table}"<<endl;
  myfile.close();
  
  return;
}
void ErrorCalc::FDSystematics()
{
  Initialize();
  if(!Init) return;
  
  InitializeSys();
  
  int i;
  unsigned int ie;
  
  int nr = Extrap[0]->GetNReco();
  int np = Extrap[0]->GetNPID();
  double *r = new double[nr+1];
  double *p = new double[np+1];
  for(i=0;i<nr+1;i++)
  {
    r[i] = Extrap[0]->Pred_TotalBkgd->GetYaxis()->GetBinLowEdge(i+1);
  }
  for(i=0;i<np+1;i++)
  {
    p[i] = Extrap[0]->Pred_TotalBkgd->GetXaxis()->GetBinLowEdge(i+1);
  }
  
  string systname;
  double pl=0,mn=0;
  double bksum=0,cc=0,nc=0,bnue=0;
  double ccpl=0,ncpl=0,bpl=0;
  double ccmn=0,ncmn=0,bmn=0;
  double psum=0,nsum=0;
  double ccpsum=0,ccnsum=0;
  double ncpsum=0,ncnsum=0;
  double bpsum=0,bnsum=0;
  TH2D *ch1,*ch2,*ch3;
  
  std::map<string, vector<TH2D*> >::iterator systiter = FD_NC_Plus1Sigma.begin();
  std::map<string, vector<TH2D*> >::iterator last  = FD_NC_Plus1Sigma.end();
  
  cout<<"FD Syst Table"<<endl;
  
  ofstream myfile("Tables/fdtable.tex");
  myfile.precision(2);
  myfile<<fixed;
  
  myfile<<"\\begin{table}[htp]"<<endl;
  myfile<<"\\caption{"<<Extrap[0]->GetPID()<<", Far Detector MC}"<<endl;
  myfile<<"\\centering"<<endl;
  myfile<<"\\begin{tabular}{c|cc|cc|cc|cc}"<<endl;
  myfile<<"Syst.&NC&&NuMuCC&&BNueCC&&Total&\\\\"<<endl;
  myfile<<"\\hline"<<endl;
  myfile<<"\\hline"<<endl;
  
  systiter = FD_NC_Plus1Sigma.begin();
  psum=0;
  nsum=0;
  ccpsum=0;
  ccnsum=0;
  ncpsum=0;
  ncnsum=0;
  bpsum=0;
  bnsum=0;
  while(systiter!=last)
  {
    systname = systiter->first;
    
    pl=0;
    mn=0;
    bpl=0;
    bmn=0;
    ncpl=0;
    ncmn=0;
    ccpl=0;
    ccmn=0;
    bksum=0;
    bnue=0;
    cc=0;
    nc=0;
    for(ie=0;ie<Extrap.size();ie++)
    {
      ch1 = (TH2D*)Extrap[ie]->FDMC[Background::kNC]->Clone();
      ch1->Multiply(FD_NC_Plus1Sigma[systname][ie]);
      ncpl+=ch1->Integral();
      ch2 = (TH2D*)Extrap[ie]->FDMC[Background::kNuMuCC]->Clone();
      ch2->Multiply(FD_NuMuCC_Plus1Sigma[systname][ie]);
      ccpl+=ch2->Integral();
      ch1->Add(ch2);
      ch3 = (TH2D*)Extrap[ie]->FDMC[Background::kBNueCC]->Clone();
      ch3->Multiply(FD_BNueCC_Plus1Sigma[systname][ie]);
      bpl+=ch3->Integral();
      ch1->Add(ch3);
      pl += ch1->Integral();
      
      ch1 = (TH2D*)Extrap[ie]->FDMC[Background::kNC]->Clone();
      ch1->Multiply(FD_NC_Minus1Sigma[systname][ie]);
      ncmn+=ch1->Integral();
      ch2 = (TH2D*)Extrap[ie]->FDMC[Background::kNuMuCC]->Clone();
      ch2->Multiply(FD_NuMuCC_Minus1Sigma[systname][ie]);
      ccmn+=ch2->Integral();
      ch1->Add(ch2);
      ch3 = (TH2D*)Extrap[ie]->FDMC[Background::kBNueCC]->Clone();
      ch3->Multiply(FD_BNueCC_Minus1Sigma[systname][ie]);
      bmn+=ch3->Integral();
      ch1->Add(ch3);
      mn += ch1->Integral();
      
      bksum+=(Extrap[ie]->FDMC[Background::kNuMuCC]->Integral() + Extrap[ie]->FDMC[Background::kNC]->Integral() + Extrap[ie]->FDMC[Background::kBNueCC]->Integral());
      cc+=Extrap[ie]->FDMC[Background::kNuMuCC]->Integral();
      nc+=Extrap[ie]->FDMC[Background::kNC]->Integral();
      bnue+=Extrap[ie]->FDMC[Background::kBNueCC]->Integral();
    }
    
    myfile<<systname<<"&";
    myfile<<100.*ncpl/nc<<"\\% &";
    if(TMath::Abs(ncpl+ncmn)<1e-20) myfile<<"-&";
    else myfile<<100.*ncmn/nc<<"\\% &";
    myfile<<100.*ccpl/cc<<"\\% &";
    if(TMath::Abs(ccpl+ccmn)<1e-20) myfile<<"-&";
    else myfile<<100.*ccmn/cc<<"\\% &";
    myfile<<100.*bpl/bnue<<"\\% &";
    if(TMath::Abs(bpl+bmn)<1e-20) myfile<<"-&";
    else myfile<<100.*bmn/bnue<<"\\% &";
    myfile<<100.*pl/bksum<<"\\% &";
    if(TMath::Abs(pl+mn)<1e-20) myfile<<"- \\\\"<<endl;
    else myfile<<100.*mn/bksum<<"\\% \\\\"<<endl;
    
    if((pl>0 && mn>0) || (pl<0 && mn<0))
    {
      if(TMath::Abs(pl)>TMath::Abs(mn))
      {
        mn = -1.*pl;
      }
      else
      {
        pl = -1.*mn;
      }
    }
    if((ccpl>0 && ccmn>0) || (ccpl<0 && ccmn<0))
    {
      if(TMath::Abs(ccpl)>TMath::Abs(ccmn))
      {
        ccmn = -1.*ccpl;
      }
      else
      {
        ccpl = -1.*ccmn;
      }
    }
    if((ncpl>0 && ncmn>0) || (ncpl<0 && ncmn<0))
    {
      if(TMath::Abs(ncpl)>TMath::Abs(ncmn))
      {
        ncmn = -1.*ncpl;
      }
      else
      {
        ncpl = -1.*ncmn;
      }
    }
    if((bpl>0 && bmn>0) || (bpl<0 && bmn<0))
    {
      if(TMath::Abs(bpl)>TMath::Abs(bmn))
      {
        bmn = -1.*bpl;
      }
      else
      {
        bpl = -1.*bmn;
      }
    }
    
    if(pl>0) psum+=(pl*pl);
    else nsum+=(pl*pl);
    if(mn>0) psum+=(mn*mn);
    else nsum+=(mn*mn);
    
    if(ccpl>0) ccpsum+=(ccpl*ccpl);
    else ccnsum+=(ccpl*ccpl);
    if(ccmn>0) ccpsum+=(ccmn*ccmn);
    else ccnsum+=(ccmn*ccmn);
    
    if(ncpl>0) ncpsum+=(ncpl*ncpl);
    else ncnsum+=(ncpl*ncpl);
    if(ncmn>0) ncpsum+=(ncmn*ncmn);
    else ncnsum+=(ncmn*ncmn);
    
    if(bpl>0) bpsum+=(bpl*bpl);
    else bnsum+=(bpl*bpl);
    if(bmn>0) bpsum+=(bmn*bmn);
    else bnsum+=(bmn*bmn);
    
    systiter++;
  }
  
  psum = sqrt(psum);
  nsum = -1.*sqrt(nsum);
  
  ccpsum = sqrt(ccpsum);
  ccnsum = -1.*sqrt(ccnsum);
  
  ncpsum = sqrt(ncpsum);
  ncnsum = -1.*sqrt(ncnsum);
  
  bpsum = sqrt(bpsum);
  bnsum = -1.*sqrt(bnsum);
  
  myfile<<"\\hline"<<endl;
  myfile<<"Total Extrap&"<<100.*ncpsum/nc<<"&"<<100.*ncnsum/nc<<"&"<<100.*ccpsum/cc<<"&"<<100.*ccnsum/cc<<"&"<<100.*bpsum/bnue<<"&"<<100.*bnsum/bnue<<"&"<<100.*psum/bksum<<"\\% &"<<100.*nsum/bksum<<"\\% \\\\"<<endl;
  myfile<<"\\end{tabular}"<<endl;
  myfile<<"\\end{table}"<<endl;
  myfile.close();
  
  return;
}
void ErrorCalc::CalculateSystErrorMatrix()
{
  if(UseGrid)
  {
    CalculateSystErrorMatrixGrid();
  }
  else
  {
    CalculateSystErrorMatrixExtrap();
  }
  return;
}
void ErrorCalc::CalculateSystErrorMatrixExtrap()
{
  //calculates the SYSTEMATIC part of the error matrix
  if(Extrap.size()==0)
  {
    cout<<"Error in ErrorCalc::CalculateSystErrorMatrix(): You must add an Extrapolate2D object.  Quitting..."<<endl;
    return;
  }
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(NDData_infile.c_str())))
  {

    cout << " CalculateSystErrorMatrixExtrap " << endl;
    cout<<"Failure to read ND Data file."<<endl;
    //SG 
    cout << " ----- > " << NDData_infile.c_str() << " " << gSystem->ExpandPathName(NDData_infile.c_str()) << endl;

    return;
  }
  
  Initialize();
  if(!Init) return;
  
  InitializeSys();
  
  string systname;
  
  int nPID = Extrap[0]->GetNPID();
  int nbins = Extrap[0]->GetNReco()*Extrap[0]->GetNPID();
//   int nt = Extrap[0]->GetNTrue();
  
  std::map<string, vector<TH2D*> >::iterator systiter = FN_NC_Plus1Sigma.begin();
  std::map<string, vector<TH2D*> >::iterator last  = FN_NC_Plus1Sigma.end();
  
  unsigned int ie;
  int i,j;
  int ip,ir,jp,jr;
  double element;
  double di,dj;
  
  CovMatrix->Reset();
  
  for(i=0;i<nbins;i++)
  {
    ir = int(i/nPID);
    ip = i%nPID;
      
    for(j=0;j<nbins;j++)
    {
      jr = int(j/nPID);
      jp = j%nPID;
      
      element = 0;
      systiter = FN_NC_Plus1Sigma.begin();
      while(systiter!=last)
      {
        systname = systiter->first;
        
        di=0;
        dj=0;
        for(ie=0;ie<Extrap.size();ie++)
        {
          di += FN_NC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kNC]->GetBinContent(ip+1,ir+1) + FN_NuMuCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kNuMuCC]->GetBinContent(ip+1,ir+1) + NuTauCC_MC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip+1,ir+1) + NueCC_MC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1);
          
          dj += FN_NC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*Extrap[ie]->Pred[Background::kNC]->GetBinContent(jp+1,jr+1) + FN_NuMuCC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*Extrap[ie]->Pred[Background::kNuMuCC]->GetBinContent(jp+1,jr+1) + NuTauCC_MC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(jp+1,jr+1) + NueCC_MC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(jp+1,jr+1);
          
          if(Extrap[0]->GetFNforBeamNue())
          {
          di += FN_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
          
          dj += FN_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(jp+1,jr+1);
            
          }
          else
          {
            di += FD_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
            
            dj += FD_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(jp+1,jr+1);
            
	     }
        }
        
        element += (di*dj);
        
        systiter++;
      }
      
      CovMatrix->SetBinContent(i+1,j+1,element);
    }
  }
  
  std::map<string,TH2D*>::iterator coviter = ExtraCovariance.begin();
  std::map<string,TH2D*>::iterator covlast  = ExtraCovariance.end();
  double ni,nj;
  //  int k=0;
  unsigned  int k=0;
  int flag;
  
  k=0;
  //  while(coviter!=covlast)
  while(k<ExtraCovariance.size())
  {
    //    systname = coviter->first;
    systname = ExtraCovariance_SystName.at(k);
    flag = ExtraCovariance_Flag.at(k);
    for(i=0;i<nbins;i++)
    {
      ir = int(i/nPID);
      ip = i%nPID;
      
      ni=0;
      for(ie=0;ie<Extrap.size();ie++)
      {
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
      
      for(j=0;j<nbins;j++)
      {
        jr = int(j/nPID);
        jp = j%nPID;
        
        nj=0;
        for(ie=0;ie<Extrap.size();ie++)
        {
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
void ErrorCalc::CalculateSystErrorMatrixGrid()
{
  //calculates the SYSTEMATIC part of the error matrix
  if(Extrap.size()==0)
  {
    cout<<"Error in ErrorCalc::CalculateSystErrorMatrixGrid(): You must add an Extrapolate2D object.  Quitting..."<<endl;
    return;
  }
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(NDData_infile.c_str())))
  {

    cout << " CalculateSystErrorMatrixGrid " << endl;
    cout<<"Failure to read ND Data file."<<endl;
    return;
  }
  
  Initialize();
  if(!Init) return;
  
  InitializeSys();
  
  string systname;
  
  int nPID = Extrap[0]->GetNPID();
  int nbins = Extrap[0]->GetNReco()*Extrap[0]->GetNPID();
  
  std::map<string, vector<TH2D*> >::iterator systiter = FN_NC_Plus1Sigma.begin();
  std::map<string, vector<TH2D*> >::iterator last  = FN_NC_Plus1Sigma.end();
  
  unsigned int ie;
  int i,j;
  int ip,ir,jp,jr;
  double element;
  double di,dj;
  
  CovMatrix->Reset();
  
  for(i=0;i<nbins;i++)
  {
    ir = int(i/nPID);
    ip = i%nPID;
    
    for(j=0;j<nbins;j++)
    {
      jr = int(j/nPID);
      jp = j%nPID;
      
      element = 0;
      systiter = FN_NC_Plus1Sigma.begin();
      while(systiter!=last)
      {
        systname = systiter->first;
        
        di=0;
        dj=0;
        for(ie=0;ie<Extrap.size();ie++)
        {
          di += FN_NC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*GridPred[Background::kNC][ie]->GetBinContent(i+1) + FN_NuMuCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*GridPred[Background::kNuMuCC][ie]->GetBinContent(i+1) + NuTauCC_MC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*GridPred[Background::kNuTauCC][ie]->GetBinContent(i+1) + NueCC_MC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*GridPred[Background::kNueCC][ie]->GetBinContent(i+1);
          
          dj += FN_NC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*GridPred[Background::kNC][ie]->GetBinContent(j+1) + FN_NuMuCC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*GridPred[Background::kNuMuCC][ie]->GetBinContent(j+1) + NuTauCC_MC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*GridPred[Background::kNuTauCC][ie]->GetBinContent(j+1) + NueCC_MC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*GridPred[Background::kNueCC][ie]->GetBinContent(j+1);
          
          if(Extrap[0]->GetFNforBeamNue())
          {
            di += FN_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*GridPred[Background::kBNueCC][ie]->GetBinContent(i+1);
            
            dj += FN_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*GridPred[Background::kBNueCC][ie]->GetBinContent(j+1);
            
          }
          else
          {
            di += FD_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*GridPred[Background::kBNueCC][ie]->GetBinContent(i+1);
            
            dj += FD_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(jp+1,jr+1)*GridPred[Background::kBNueCC][ie]->GetBinContent(j+1);
            
          }
        }
        
        element += (di*dj);
        
        systiter++;
      }
      
      CovMatrix->SetBinContent(i+1,j+1,element);
    }
  }
  
  std::map<string,TH2D*>::iterator coviter = ExtraCovariance.begin();
  std::map<string,TH2D*>::iterator covlast  = ExtraCovariance.end();
  double ni,nj;
  int k=0;
  int flag;
  
  k=0;
  while(coviter!=covlast)
  {
    systname = coviter->first;
    flag = ExtraCovariance_Flag.at(k);
    for(i=0;i<nbins;i++)
    {
      
      ni=0;
      for(ie=0;ie<Extrap.size();ie++)
      {
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
      
      for(j=0;j<nbins;j++)
      {
        nj=0;
        for(ie=0;ie<Extrap.size();ie++)
        {
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
void ErrorCalc::MakeFracError_Lists()
{
  //calculates Fractional Errors for the Standard Likelihood
  if(Extrap.size()==0)
  {
    cout<<"Error in ErrorCalc::CalculateSystErrorMatrix(): You must add an Extrapolate2D object.  Quitting..."<<endl;
    return;
  }
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(NDData_infile.c_str())))
  {
    cout << " MakeFracError_Lists " << endl;
    cout<<"Failure to read ND Data file."<<endl;
    return;
  }

  Initialize();
  if(!Init) return;
  
  InitializeSys();

  string systname;
  
  int nPID = Extrap[0]->GetNPID();
  int nbins = Extrap[0]->GetNReco()*Extrap[0]->GetNPID();
//   int nt = Extrap[0]->GetNTrue();
  
  std::map<string, vector<TH2D*> >::iterator systiter = FN_NC_Plus1Sigma.begin();
  std::map<string, vector<TH2D*> >::iterator last  = FN_NC_Plus1Sigma.end();
  
  unsigned int ie;
  int i;
  unsigned int j;
  int ip,ir;
  double element_bkg;
  double element_sig;
  double tot_bkg, tot_sig;
  double dnc, dnumu, dbnue;
  double dnue, dtau;
  double *bins = new double[nbins+1];;

  //binning for new histograms:

  for (int k=0; k<nbins+1; k++){
    bins[k] = Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetBinLowEdge(k+1);
  }

 //May or may not use this:
  //SysBkgd_Plus1Sigma.clear();
  //SysSig_Plus1Sigma.clear();

  i=0;
  j=0;

  systiter = FN_NC_Plus1Sigma.begin();
  while(systiter!=last)
    {
      systname = systiter->first;

      if (j>=SysBkgd_Plus1Sigma.size()){
      SysBkgd_Plus1Sigma.push_back(new TH1D(Form("bkg_sys_%s",systname.c_str()),"",nbins,bins));
      } else SysBkgd_Plus1Sigma[j]->Reset();

      if (j>=SysSig_Plus1Sigma.size()){
      SysSig_Plus1Sigma.push_back(new TH1D(Form("sig_sys_%s",systname.c_str()),"",nbins,bins));
      } else SysSig_Plus1Sigma[j]->Reset();

      if (SysBkgd_Plus1Sigma.size() != SysSig_Plus1Sigma.size()){
	cout << "ERROR: Got out of phase somehow in FracErr Calc." << endl;
	return;
      }

  
      for(i=0;i<nbins;i++)
	{
	  ir = int(i/nPID);
	  ip = i%nPID;
      
      
	  element_sig = 0;
	  element_bkg = 0;

	  dnc=0;        
	  dnumu=0;        
	  dbnue=0;        
	  dnue=0;
	  dtau=0;

	  tot_bkg=0;
	  tot_sig=0;

	  element_sig=0;
	  element_bkg=0;


	  for(ie=0;ie<Extrap.size();ie++)
	    {
	      dnc += FN_NC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kNC]->GetBinContent(ip+1,ir+1);
	      dnumu += FN_NuMuCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kNuMuCC]->GetBinContent(ip+1,ir+1);
	      dtau += NuTauCC_MC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip+1,ir+1);
	      dnue += NueCC_MC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1);
	      
          
	      if(Extrap[0]->GetFNforBeamNue())
	        {
	        dbnue += FN_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
	        
	       }
	      else
	      {
	         dbnue += FD_BNueCC_Plus1Sigma[systname][ie]->GetBinContent(ip+1,ir+1)*Extrap[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
	     	  
	      }

	      //Calculate total prediction in that in bin, too:	      
	      tot_bkg += Extrap[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1);
	      tot_sig += Extrap[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1);
	    }       
        
	  if (tot_sig>0){
	    element_sig = dnue/tot_sig;
	  } else element_sig = 0;

	  if (tot_bkg>0){
	    element_bkg = (dnc + dnumu + dtau + dbnue)/tot_bkg;
	  } else element_bkg = 0;

	  SysBkgd_Plus1Sigma[j]->SetBinContent(i+1,element_bkg);
	  SysSig_Plus1Sigma[j]->SetBinContent(i+1,element_sig);

	}

      j++;
      systiter++;
      

    }

  
  return;
}

void ErrorCalc::CalculateHOOError()
{
  //this function assumes the error matrix is normalized to the same POT that the denominator of the F/N ratio is (i.e.1e19)
  
  //note!!!! the binning convention used in the covariance matrix in the HOOHE files is:
  //bin = ienergy + ipid*Nenergy + irun*Nenergy*Npid + icomponent*Nenergy*Npid*Nrun
  //the FD prediction bin number convention is bin = ipid + ienergy*Npid
  
  Initialize();
  if(!Init) return;
  
  InitializeDecompSys();
  
  if(!InitDecompSys)
  {
    cout<<"Error in setting up CalculateHOOError().  Quitting..."<<endl;
    return;
  } 
  
  int n = NDCovMatrix->GetNbinsX();//dimension of error matrix = pid bins x energy bins x number of runs x number of components
  
  unsigned int npid = Extrap[0]->GetNPID();//number pid bins
  unsigned int nenergy = Extrap[0]->GetNReco();//number of reco energy bins
//   int ncomp = 3;//number of background components
  unsigned int nrun = Extrap.size();//number of runs
  
  unsigned int i;
  int j,k,m;
  int ipid,ienergy,irun,icomp;
  unsigned int ipred;
  Background::Background_t bg = Background::kNC;
  double temp;
  
  vector< vector<double> > A;//[npid*nenergy] by [npid*nenergy*nrun*ncomp] matrix of FD predictions to go between the covariance matrix binning to the FD prediction binning
  for(i=0;i<npid*nenergy;i++)
  {
    A.push_back( vector<double>() );
    for(j=0;j<n;j++)
    {
      icomp = int(j/(npid*nenergy*nrun));
      irun = int((j - icomp*npid*nenergy*nrun)/(nenergy*npid));
      ipid = int((j - icomp*npid*nenergy*nrun - irun*nenergy*npid)/nenergy);
      ienergy = j%nenergy;
      
      ipred = ipid + ienergy*npid;//FD prediction bin corresponding to HOOHE covariance matrix bin j
      
      if(icomp==0) bg = Background::kNC;
      else if(icomp==1) bg = Background::kNuMuCC;
      else if(icomp==2) bg = Background::kBNueCC;
      
      A[i].push_back(0);//start out with zeros
      if(ipred==i) //if this FD prediction bin in the matrix A (index i) corresponds to the FD prediction bin from the jth bin of the Decomp covariance matrix (given by ipred)
      {
        if(UseGrid) A[i][j] = GridPred[bg][irun]->GetBinContent(ipred+1);
        else A[i][j] = Extrap[irun]->Pred[bg]->GetBinContent(ipid+1,ienergy+1);
      }
    }
  }
  
//   cout<<"(NC+NuMuCC+BNueCC) background prediction:"<<endl;
//   vector<double> F;
//   for(i=0;i<npid*nenergy;i++)
//   {
//     F.push_back(0);
//     for(j=0;j<n;j++)
//     {
//       F[i] += A[i][j];
//     }
//     cout<<F[i]<<endl;
//   }
  
  CovMatrix_Decomp->Reset();
  for(unsigned int i=0;i<npid*nenergy;i++)
  {
    for(unsigned int j=0;j<npid*nenergy;j++)
    {
      temp=0;
      for(k=0;k<n;k++)
      {
        for(m=0;m<n;m++)
        {
          temp+=(NDCovMatrix->GetBinContent(k+1,m+1)*A[i][k]*A[j][m]);
        }
      }
      CovMatrix_Decomp->SetBinContent(i+1,j+1,temp);
    }
  }
  
  return;
}
void ErrorCalc::SetGridPred(int nbins,vector< vector<double> > nc, vector< vector<double> > cc, vector< vector<double> > bnue, vector< vector<double> > tau, vector< vector<double> > sig)
{
  if(Extrap.size()==0)
  {
    cout<<"Error in ErrorCalc::SetGridPred(): You must add an Extrapolate2D object.  Quitting..."<<endl;
    return;
  }
  
  Initialize();
  if(!Init) return;
  
  int i;
  unsigned int ie;
  
  if(nbins!=GridPred[Background::kNC][0]->GetNbinsX())
  {
    cout<<"Warning in SetGridPred(): nbins != GridPred[]->GetNbinsX()!!!"<<endl;
  }
  
  for(ie=0;ie<Extrap.size();ie++)
  {
    GridPred[Background::kNC][ie]->Reset();
    GridPred[Background::kNuMuCC][ie]->Reset();
    GridPred[Background::kBNueCC][ie]->Reset();
    GridPred[Background::kNuTauCC][ie]->Reset();
    GridPred[Background::kNueCC][ie]->Reset();
    
    for(i=0;i<nbins;i++)
    {
      GridPred[Background::kNC][ie]->SetBinContent(i+1,nc[i][ie]);
      GridPred[Background::kNuMuCC][ie]->SetBinContent(i+1,cc[i][ie]);
      GridPred[Background::kBNueCC][ie]->SetBinContent(i+1,bnue[i][ie]);
      GridPred[Background::kNuTauCC][ie]->SetBinContent(i+1,tau[i][ie]);
      GridPred[Background::kNueCC][ie]->SetBinContent(i+1,sig[i][ie]);
    }
  }
  
  return;
}
