#define NueFit2D_Joint_C

#include "NueAna/MultiBinAna/NueFit2D_Joint.h"


NueFit2D_Joint::NueFit2D_Joint()
{

  printMatrix=false;

  NObs=0;
  ErrCalc=0;
  ErrorMatrix=0;
  InvErrorMatrix=0;
  ExternalErrorMatrix=0;

  SetOutputFile();
  SetNExperiments();
  
  SetGridNorm();
  GridScale_Normal=1.;
  GridScale_Inverted=1.;
  GridScale_Normal_FHC=1.;
  GridScale_Inverted_FHC=1.;
  GridScale_Normal_RHC=1.;
  GridScale_Inverted_RHC=1.;
  
  nBins=0;
  nBinsFHC=0;
  nBinsRHC=0;
  
  IncludeOscParErrs = false;
  
  FormatChi2Hists();
  
  SetNDeltaSteps();
  SetNSinSq2Th13Steps();
  SetDeltaRange();
  SetSinSq2Th13Range();
  SetepsetauRange();
  SetNepsetauSteps();
  
  SetFitMethod();
  
  potScaleFHC=1.0;
  potScaleRHC=1.0;

  Theta12 = 0.59365;
  Theta12Unc = 0.041;
  Theta23 = 0.7854;
  Theta23Unc = 0.122;
  Theta13 = 0.1594;
  Theta13Unc = 0.011;
  DeltaMSq32 = 2.32e-3;
  DeltaMSq32Unc = 0.11e-3;
  DeltaMSq21 = 8.0e-5;
  DeltaMSq21Unc = 0.6e-5;
  delta = 0;
  
  //NSI parameters:
  Eps_ee = 0.0;
  Eps_emu = 0.0;
  Eps_etau = 0.0;
  Eps_etau_unc = 0.01;
  Eps_mumu = 0.0;
  Eps_mutau = 0.0; 
  Eps_tautau = 0.0; 
  Delta_emu = 0.0;
  Delta_etau = 0.0;
  Delta_etau_unc = 0.01;
  Delta_mutau = 0.0;

  return;
}
NueFit2D_Joint::~NueFit2D_Joint()
{
}

//RBNOTE: This is here purely for diagnostic purposes
void NueFit2D_Joint::CalculateDlnLMatrix(TH2D *SysMatrix=0, TH2D *HOOHEMatrix=0){
  //Calculate the covariance matrix for the "bin by bin" method

  //Make new inverted matrix if not done already
  if (InvErrorMatrix==0){
    if (SysMatrix){
      InvErrorMatrix = (TH2D*)SysMatrix->Clone("InvErrorMatrix");
    } else if (HOOHEMatrix) {
      InvErrorMatrix = (TH2D*)HOOHEMatrix->Clone("InvErrorMatrix");
    } else {
      cout << "Didn't give me any matrices to invert!" << endl;
      return;
    }
  }
  //Reset everything to zero:
  InvErrorMatrix->Reset();
  int totbins = nBins;
  int i=0;
  int j=0;
  int k=0;
  double scalei=1.0;
  double scalej=1.0;
  int nFHC=0;
  int nRHC=0;

  if (ExtrapFHC.size()>0) nFHC = ExtrapFHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ExtrapRHC.size()>0) nRHC = ExtrapRHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (printMatrix) cout << "PRE-INVERSION" << endl;
  //Add together error elements into a single array
  double *Varray = new double[totbins*totbins];
  for(i=0;i<totbins;i++)
  {

    if (i<nFHC) { scalei=potScaleFHC; }
    else { scalei=potScaleRHC; }

    for(j=0;j<totbins;j++)
    {

    if (j<nFHC) { scalej=potScaleFHC; }
    else { scalej=potScaleRHC; }

      k = i*totbins + j;
      Varray[k]=0;
      if (SysMatrix!=0) {
	Varray[k] += scalei*scalej*SysMatrix->GetBinContent(i+1,j+1);
      }

      if (HOOHEMatrix!=0) {
	Varray[k] += HOOHEMatrix->GetBinContent(i+1,j+1);
      }

      Varray[k] = Varray[k]*1.0e0;

    }

 }

  //Hand elements to a TMatrix
  TMatrixD *Vmatrix = new TMatrixD(totbins,totbins,Varray);

  if (printMatrix) Vmatrix->Print();
  //Vmatrix->Print();

  cout.precision(20);
  //Insert it (determ found for debugging purposes)!
  double determ;
  TMatrixD Vinvmatrix = Vmatrix->Invert(&determ);
  if (printMatrix) cout << "DETERM=" << determ << endl;
  cout.precision(5);

  if (printMatrix) Vinvmatrix.Print();

  //Extract the array of the inverted matrix:
  Double_t *Vinvarray = Vinvmatrix.GetMatrixArray();

  //Turn it into the inverted covariance matrix
  if (printMatrix) cout << "INVERTED" << endl;
  //make Vinv out of array
  for(i=0;i<totbins;i++)
  {
    for(j=0;j<totbins;j++)
    {
      InvErrorMatrix->SetBinContent(i+1,j+1,Vinvarray[i*totbins + j]*1.0e0);
    }
  }

  delete [] Varray;
  
  if (printMatrix) cout << "-------------" << endl;
  return;

}

void NueFit2D_Joint::AddExtrap(Extrapolate2D* /*E*/)
{

  cout << "Can't add generic extrapolation for a joint fit.  Use FHC/RHC options instead.  No extrapolation was added." << endl;
  return;

}

void NueFit2D_Joint::AddExtrapFHC(Extrapolate2D *E)
{
  ExtrapFHC.push_back(E);
  return;

}

void NueFit2D_Joint::AddExtrapRHC(Extrapolate2D *E)
{
  ExtrapRHC.push_back(E);
  return;

}

void NueFit2D_Joint::CombineFHCRHC()
{

  Bkgd->Reset();
  Sig->Reset();

  int nFHC=0;
  int nRHC=0;
 
  if (ExtrapFHC.size()>0) nFHC = ExtrapFHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ExtrapRHC.size()>0) nRHC = ExtrapRHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  unsigned int ie=0;
  TH1D* tmpsig = (TH1D*)Sig->Clone("tmpsig");
  tmpsig->Reset();
  TH1D* tmpbkg = (TH1D*)Bkgd->Clone("tmpbkg");
  tmpbkg->Reset();

  for(ie=0;ie<ExtrapFHC.size();ie++)
    {
      tmpbkg->Reset();
      tmpsig->Reset();
      for (int nb = 0; nb<nFHC; nb++){
	tmpbkg->SetBinContent(nb+1,potScaleFHC*(ExtrapFHC[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(nb+1)));
	tmpsig->SetBinContent(nb+1,potScaleFHC*(ExtrapFHC[ie]->Pred_Signal_VsBinNumber->GetBinContent(nb+1)));
      }
      Bkgd->Add(tmpbkg);
      Sig->Add(tmpsig);
    }

  for(ie=0;ie<ExtrapRHC.size();ie++)
    {
      tmpbkg->Reset();
      tmpsig->Reset();
      for (int nb = 0; nb<nRHC; nb++){
	tmpbkg->SetBinContent(nb+1+nFHC,potScaleRHC*(ExtrapRHC[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(nb+1)));
	tmpsig->SetBinContent(nb+1+nFHC,potScaleRHC*(ExtrapRHC[ie]->Pred_Signal_VsBinNumber->GetBinContent(nb+1)));
      }
      Bkgd->Add(tmpbkg);
      Sig->Add(tmpsig);

    }

}

void NueFit2D_Joint::ReadGridFiles()
{
  double fpfhc;
  double fprhc;
  unsigned int i;
  TTree *temp=0;
  unsigned int j;
  
  for(i=0;i<nBins;i++)
  {
    grid_bin_oscparerr.push_back(0);
  }
 

  TFile *fnorm;
  GridTree_Normal.clear();
  GridTree_2_Normal.clear();
  nPts_Normal = 0;
  if(!gSystem->AccessPathName(gSystem->ExpandPathName(GridFileName_Normal.c_str())))//if file exists
  {
    fnorm = new TFile(gSystem->ExpandPathName(GridFileName_Normal.c_str()),"READ");
    for(i=0;i<nBins;i++)
    {
      temp = (TTree*)fnorm->Get(Form("Bin_%i",i));
      temp->SetName(Form("Bin_%i_Normal",i));
      temp->SetBranchAddress("Background",&grid_background);
      temp->SetBranchAddress("Signal",&grid_signal);
      temp->SetBranchAddress("Delta",&grid_delta);
      temp->SetBranchAddress("Th13Axis",&grid_sinsq2th13);
      if(IncludeOscParErrs)
      {
        temp->SetBranchAddress("DNExp_DOscPars",&grid_bin_oscparerr[i]);
      }
      GridTree_Normal.push_back(temp);
      
      GridTree_2_Normal.push_back( vector<TTree*>() );

      if (i<nBinsFHC){
	for(j=0;j<ExtrapFHC.size();j++)
	  {
	    temp = (TTree*)fnorm->Get(Form("Bin_%i_Run_%i",i,j));
	    temp->SetName(Form("Bin_%i_Run_%i_Normal",i,j));
	    temp->SetBranchAddress("NC",&grid_nc);
	    temp->SetBranchAddress("NuMuCC",&grid_numucc);
	    temp->SetBranchAddress("BNueCC",&grid_bnuecc);
	    temp->SetBranchAddress("NuTauCC",&grid_nutaucc);
	    temp->SetBranchAddress("Signal",&grid_nue);
	    temp->SetBranchAddress("Delta",&grid_delta);
	    temp->SetBranchAddress("Th13Axis",&grid_sinsq2th13);
	    GridTree_2_Normal[i].push_back(temp);
	  }
      } else if (i>=nBinsFHC){
	for(j=0;j<ExtrapRHC.size();j++)
	  {
	    temp = (TTree*)fnorm->Get(Form("Bin_%i_Run_%i",i,j));
	    temp->SetName(Form("Bin_%i_Run_%i_Normal",i,j));
	    temp->SetBranchAddress("NC",&grid_nc);
	    temp->SetBranchAddress("NuMuCC",&grid_numucc);
	    temp->SetBranchAddress("BNueCC",&grid_bnuecc);
	    temp->SetBranchAddress("NuTauCC",&grid_nutaucc);
	    temp->SetBranchAddress("Signal",&grid_nue);
	    temp->SetBranchAddress("Delta",&grid_delta);
	    temp->SetBranchAddress("Th13Axis",&grid_sinsq2th13);
	    GridTree_2_Normal[i].push_back(temp);
	  }
      }
    }
    nPts_Normal = GridTree_Normal[0]->GetEntries();
    
    paramtree_Normal = (TTree*)fnorm->Get("paramtree");
    paramtree_Normal->SetName("paramtree_Normal");
    paramtree_Normal->SetBranchAddress("farPOTFHC",&fpfhc);
    paramtree_Normal->SetBranchAddress("farPOTRHC",&fprhc);
    paramtree_Normal->SetBranchAddress("Theta12",&grid_n_th12);
    paramtree_Normal->SetBranchAddress("Theta23",&grid_n_th23);
    paramtree_Normal->SetBranchAddress("DeltaMSq23",&grid_n_dm2_32);
    paramtree_Normal->SetBranchAddress("DeltaMSq12",&grid_n_dm2_21);
    paramtree_Normal->SetBranchAddress("Theta13",&grid_n_th13);
    paramtree_Normal->GetEntry(0);
    
    if(GridNorm>0)
    {
      GridScale_Normal_FHC = GridNorm/fpfhc;
      GridScale_Normal_RHC = GridNorm/fprhc;
    }
  }
  else
  {
    cout<<"Grid file (normal hierarchy) doesn't exist."<<endl;
    return;
  }

  TFile *finvt;
  GridTree_Inverted.clear();
  nPts_Inverted = 0;
  if(!gSystem->AccessPathName(gSystem->ExpandPathName(GridFileName_Inverted.c_str())))//if file exists
  {
    finvt = new TFile(gSystem->ExpandPathName(GridFileName_Inverted.c_str()),"READ");
    for(i=0;i<nBins;i++)
    {
      temp = (TTree*)finvt->Get(Form("Bin_%i",i));
      temp->SetName(Form("Bin_%i_Inverted",i));
      temp->SetBranchAddress("Background",&grid_background);
      temp->SetBranchAddress("Signal",&grid_signal);
      temp->SetBranchAddress("Delta",&grid_delta);
      temp->SetBranchAddress("Th13Axis",&grid_sinsq2th13);
      if(IncludeOscParErrs)
      {
        temp->SetBranchAddress("DNExp_DOscPars",&grid_bin_oscparerr[i]);
      }
      GridTree_Inverted.push_back(temp);
      
      GridTree_2_Inverted.push_back( vector<TTree*>() );

      if (i<nBinsFHC){      
	for(j=0;j<ExtrapFHC.size();j++)
	  {
	    temp = (TTree*)finvt->Get(Form("Bin_%i_Run_%i",i,j));
	    temp->SetName(Form("Bin_%i_Run_%i_Inverted",i,j));
	    temp->SetBranchAddress("NC",&grid_nc);
	    temp->SetBranchAddress("NuMuCC",&grid_numucc);
	    temp->SetBranchAddress("BNueCC",&grid_bnuecc);
	    temp->SetBranchAddress("NuTauCC",&grid_nutaucc);
	    temp->SetBranchAddress("Signal",&grid_nue);
	    temp->SetBranchAddress("Delta",&grid_delta);
	    temp->SetBranchAddress("Th13Axis",&grid_sinsq2th13);
	    GridTree_2_Inverted[i].push_back(temp);
	  }
      } else if (i>=nBinsFHC){
	for(j=0;j<ExtrapRHC.size();j++)
	  {
	    temp = (TTree*)finvt->Get(Form("Bin_%i_Run_%i",i,j));
	    temp->SetName(Form("Bin_%i_Run_%i_Inverted",i,j));
	    temp->SetBranchAddress("NC",&grid_nc);
	    temp->SetBranchAddress("NuMuCC",&grid_numucc);
	    temp->SetBranchAddress("BNueCC",&grid_bnuecc);
	    temp->SetBranchAddress("NuTauCC",&grid_nutaucc);
	    temp->SetBranchAddress("Signal",&grid_nue);
	    temp->SetBranchAddress("Delta",&grid_delta);
	    temp->SetBranchAddress("Th13Axis",&grid_sinsq2th13);
	    GridTree_2_Inverted[i].push_back(temp);
	  }
      }
    }
    nPts_Inverted = GridTree_Inverted[0]->GetEntries();
    
    paramtree_Inverted = (TTree*)finvt->Get("paramtree");
    paramtree_Inverted->SetName("paramtree_Inverted");
    paramtree_Inverted->SetBranchAddress("farPOTFHC",&fpfhc);
    paramtree_Inverted->SetBranchAddress("farPOTRHC",&fprhc);
    paramtree_Inverted->SetBranchAddress("Theta12",&grid_i_th12);
    paramtree_Inverted->SetBranchAddress("Theta23",&grid_i_th23);
    paramtree_Inverted->SetBranchAddress("DeltaMSq23",&grid_i_dm2_32);
    paramtree_Inverted->SetBranchAddress("DeltaMSq12",&grid_i_dm2_21);
    paramtree_Inverted->SetBranchAddress("Theta13",&grid_i_th13);
    paramtree_Inverted->GetEntry(0);
    

    if(GridNorm>0)
    {
      GridScale_Inverted_FHC = GridNorm/fpfhc;
      GridScale_Inverted_RHC = GridNorm/fprhc;
    }
  }
  else
  {
    cout<<"Grid file (inverted hierarchy) doesn't exist."<<endl;
    return;
  }
  
  return;
}


void NueFit2D_Joint::DefineBinDlnLMinuit(){
  //Function sets the maximum size of Minuit.

  //Clear things first:

  int npar = 0;

  if (nBins != 0){
    npar = nBins;
  } else if (NObs != 0){
      npar = NObs->GetNbinsX();
  } else if ((ExtrapFHC.size()+ExtrapRHC.size())!=0){
    npar = ExtrapFHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX()
     + ExtrapRHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();

  } else {
    cout << "ERROR: Add extrapolation or NObs before initializing Minuit size"
	 << endl;
    return;
  }

  //Make new minuit:
  minuit = new TMinuit(npar);
  minuit->SetPrintLevel(-1);

  double arglist[1];
  int ierflg=0;

  arglist[0] = 1.e-5;
  minuit->mnexcm("SET EPS",arglist,1,ierflg);


}

void NueFit2D_Joint::DefineStdDlnLMinuit(){
  //Function sets the maximum size of Minuit.

  //Clear things first:
  //minuit->mncler();

  int npar = 0;

  //Number of systematics inputted:
  if (FracErr_Bkgd_List.size()!=0){
    npar = FracErr_Bkgd_List.size();
  } 

  int nb = 0;
  //Number of bins (for HOOHE):
  if (nBins != 0){
    nb = nBins;
  } else if (NObs != 0){
    nb = NObs->GetNbinsX();
  } else if ((ExtrapFHC.size()+ExtrapRHC.size())!=0){
    nb = ExtrapFHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX()
     + ExtrapRHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  } else {
    cout << "ERROR: Add extrapolation or NObs before initializing Minuit size"
	 << endl;
    return;
  }

  npar += nb;

  //make new minuit
  minuit = new TMinuit(npar+nb);
  minuit->SetPrintLevel(-1);

}

void NueFit2D_Joint::RunDataGrid(string filenorm,string fileinvt)
{
  if(NObs==0)
  {
    cout<<"NObs not set.  Quitting..."<<endl;
    return;
  }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return;
  }
  for(unsigned int ie=0;ie<ExtrapFHC.size();ie++)
  {
    ExtrapFHC[ie]->GetPrediction();
  }
  for(unsigned int ie=0;ie<ExtrapRHC.size();ie++)
  {
    ExtrapRHC[ie]->GetPrediction();
  }
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
  {
    cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
    return;
  }
  if(ErrCalc==0)
  {
    cout<<"No ErrorCalc object set!  Quitting..."<<endl;
    return;
  }
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();
  
  nBins = NObs->GetNbinsX();
  nBinsFHC=0;
  nBinsRHC=0;  
  if (ExtrapFHC.size()>0) nBinsFHC = ExtrapFHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ExtrapRHC.size()>0) nBinsRHC = ExtrapRHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (nBins!=(nBinsFHC+nBinsRHC)) {
    cout << "ERROR: Bin mismatch between NObs and expected distribution.  Aborting..." << endl;
    return;
  }

  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  ReadGridFiles();
  
  if(nPts_Normal==0 || nPts_Inverted==0) return;
  
  TH1D *nexp_bkgd = new TH1D("nexp_bkgd","",nBins,-0.5,nBins-0.5);
  TH1D *nexp_signal = new TH1D("nexp_signal","",nBins,-0.5,nBins-0.5);
  TH1D *nexp = new TH1D("nexp","",nBins,-0.5,nBins-0.5);
  
  int i;
  unsigned int j,k;
  
  double chi2data,chi2min;
  
  vector< vector<double> > nc,numucc,bnuecc,nutaucc,sig;
  for(j=0;j<nBins;j++)
  {
    nc.push_back( vector<double>() );
    numucc.push_back( vector<double>() );
    bnuecc.push_back( vector<double>() );
    nutaucc.push_back( vector<double>() );
    sig.push_back( vector<double>() );

  if (j<nBinsFHC){
      for(k=0;k<ExtrapFHC.size();k++)
	{
	  nc[j].push_back(0);
	  numucc[j].push_back(0);
	  bnuecc[j].push_back(0);
	  nutaucc[j].push_back(0);
	  sig[j].push_back(0);
	}
    } else if (j>=nBinsFHC){
      for(k=0;k<ExtrapRHC.size();k++)
	{
	  nc[j].push_back(0);
	  numucc[j].push_back(0);
	  bnuecc[j].push_back(0);
	  nutaucc[j].push_back(0);
	  sig[j].push_back(0);
	}
    }
  }
  
  Bkgd = (TH1D*)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)NObs->Clone("Sig");
  Sig->Reset();
  
  //normal hierarchy
  
  ofstream myfile;
  myfile.open(gSystem->ExpandPathName(filenorm.c_str()));
  
  for(i=0;i<nPts_Normal;i++)
  {
    if(i%100==0) cout<<100.*i/nPts_Normal<<"% complete for normal hierarchy"<<endl;
    
    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();
    
    for(j=0;j<nBins;j++)
    {
      GridTree_Normal[j]->GetEntry(i);
      if (j<nBinsFHC){
	nexp_bkgd->SetBinContent(j+1,grid_background*GridScale_Normal_FHC);
	nexp_signal->SetBinContent(j+1,grid_signal*GridScale_Normal_FHC);
 	for(k=0;k<ExtrapFHC.size();k++)
	  {
	    GridTree_2_Normal[j][k]->GetEntry(i);
	    nc[j][k] = grid_nc*GridScale_Normal_FHC;
	    numucc[j][k] = grid_numucc*GridScale_Normal_FHC;
	    bnuecc[j][k] = grid_bnuecc*GridScale_Normal_FHC;
	    nutaucc[j][k] = grid_nutaucc*GridScale_Normal_FHC;
	    sig[j][k] = grid_nue*GridScale_Normal_FHC;
	  }
      } else if (j>=nBinsFHC){
	nexp_bkgd->SetBinContent(j+1,grid_background*GridScale_Normal_RHC);
	nexp_signal->SetBinContent(j+1,grid_signal*GridScale_Normal_RHC);
 	for(k=0;k<ExtrapRHC.size();k++)
	  {
	    GridTree_2_Normal[j][k]->GetEntry(i);
	    nc[j][k] = grid_nc*GridScale_Normal_RHC;
	    numucc[j][k] = grid_numucc*GridScale_Normal_RHC;
	    bnuecc[j][k] = grid_bnuecc*GridScale_Normal_RHC;
	    nutaucc[j][k] = grid_nutaucc*GridScale_Normal_RHC;
	    sig[j][k] = grid_nue*GridScale_Normal_RHC;
	  }
      }
    }
    nexp->Add(nexp_bkgd,nexp_signal,1,1);
    ErrCalc->SetGridPred(nBinsFHC+nBinsRHC,nc,numucc,bnuecc,nutaucc,sig);
    
    chi2min=GetMinLikelihood(grid_delta,true);
    
    ErrCalc->SetUseGrid(true);//use grid predictions set up above
    chi2data = 1e10;
    if(FitMethod==0)
    {
      chi2data = PoissonChi2(nexp) - chi2min;
    }
    else if(FitMethod==1)
    {
      chi2data = ScaledChi2(Bkgd,Sig) - chi2min;
    }
    else if(FitMethod==2)
    {
      chi2data = StandardChi2(nexp) - chi2min;
    }
    else if(FitMethod==3)
    {
      //Likelihood: "Standard" (N syst, N nuisance)
      //Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data = StandardLikelihood() - chi2min;
    }
    else if(FitMethod==4)
    {
      //Likelihood: Bin by Bin Calculation of Systematics
      //Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data = BinLikelihood() - chi2min;
    }
    else
    {
      cout<<"Error in RunMultiBinFC(): Unknown 'FitMethod'."<<endl;
    }
    
    myfile << grid_sinsq2th13 << " " << grid_delta << " " << nexp_signal->Integral() << " ";
    for(j=0;j<nBins;j++)
    {
      myfile << nexp_signal->GetBinContent(j+1) << " ";
    }
    myfile << chi2data << endl;
  }
  
  myfile.close();
  
  //inverted hierarchy
  
  myfile.open(gSystem->ExpandPathName(fileinvt.c_str()));
  
  for(i=0;i<nPts_Inverted;i++)
  {
    if(i%100==0) cout<<100.*i/nPts_Inverted<<"% complete for inverted hierarchy"<<endl;
    
    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();
    
    for(j=0;j<nBins;j++)
    {
      GridTree_Inverted[j]->GetEntry(i);
     if (j<nBinsFHC){
	nexp_bkgd->SetBinContent(j+1,grid_background*GridScale_Inverted_FHC);
	nexp_signal->SetBinContent(j+1,grid_signal*GridScale_Inverted_FHC);
	for(k=0;k<ExtrapFHC.size();k++)
	  {
	    GridTree_2_Inverted[j][k]->GetEntry(i);
	    nc[j][k] = grid_nc*GridScale_Inverted_FHC;
	    numucc[j][k] = grid_numucc*GridScale_Inverted_FHC;
	    bnuecc[j][k] = grid_bnuecc*GridScale_Inverted_FHC;
	    nutaucc[j][k] = grid_nutaucc*GridScale_Inverted_FHC;
	    sig[j][k] = grid_nue*GridScale_Inverted_FHC;
	  }
      } else if (j>=nBinsFHC){
	nexp_bkgd->SetBinContent(j+1,grid_background*GridScale_Inverted_RHC);
	nexp_signal->SetBinContent(j+1,grid_signal*GridScale_Inverted_RHC);
	for(k=0;k<ExtrapRHC.size();k++)
	  {
	    GridTree_2_Inverted[j][k]->GetEntry(i);
	    nc[j][k] = grid_nc*GridScale_Inverted_RHC;
	    numucc[j][k] = grid_numucc*GridScale_Inverted_RHC;
	    bnuecc[j][k] = grid_bnuecc*GridScale_Inverted_RHC;
	    nutaucc[j][k] = grid_nutaucc*GridScale_Inverted_RHC;
	    sig[j][k] = grid_nue*GridScale_Inverted_RHC;
	  }
      }
    }
    nexp->Add(nexp_bkgd,nexp_signal,1,1);
    ErrCalc->SetGridPred(nBins,nc,numucc,bnuecc,nutaucc,sig);
    
    chi2min=GetMinLikelihood(grid_delta,false);
    
    ErrCalc->SetUseGrid(true);//use grid predictions set up above
    chi2data = 1e10;
    if(FitMethod==0)
    {
      chi2data = PoissonChi2(nexp) - chi2min;
    }
    else if(FitMethod==1)
    {
      chi2data = ScaledChi2(Bkgd,Sig) - chi2min;
    }
    else if(FitMethod==2)
    {
      chi2data = StandardChi2(nexp) - chi2min;
    }
    else if(FitMethod==3)
    {
      //Likelihood: "Standard" (N syst, N nuisance)
      //Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data = StandardLikelihood() - chi2min;
    }
    else if(FitMethod==4)
    {
      //Likelihood: Bin by Bin Calculation of Systematics
      //Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data = BinLikelihood() - chi2min;
    }
    else
    {
      cout<<"Error in RunMultiBinFC(): Unknown 'FitMethod'."<<endl;
    }
    
    myfile << grid_sinsq2th13 << " " << grid_delta << " " << nexp_signal->Integral() << " ";
    for(j=0;j<nBins;j++)
    {
      myfile << nexp_signal->GetBinContent(j+1) << " ";
    }
    myfile << chi2data << endl;
  }
  
  myfile.close();
  
  return;
}

void NueFit2D_Joint::RunMultiBinPseudoExpts(bool Print)
{
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
  {
    cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
    return;
  }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return;
  }
  for(unsigned int ie=0;ie<ExtrapFHC.size();ie++)
  {
    ExtrapFHC[ie]->GetPrediction();
  }
  for(unsigned int ie=0;ie<ExtrapRHC.size();ie++)
  {
    ExtrapRHC[ie]->GetPrediction();
  }
  if(ErrCalc==0)
  {
    cout<<"Need to set ErrorCalc object!  Quitting..."<<endl;
    return;
  }
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();

  nBinsFHC=0;
  nBinsRHC=0;  
  if (ExtrapFHC.size()>0) nBinsFHC = ExtrapFHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ExtrapRHC.size()>0) nBinsRHC = ExtrapRHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  nBins = nBinsFHC + nBinsRHC;

  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  TH2D *Error4Expts = new TH2D("Error4Expts","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  ReadGridFiles();
  
  if(nPts_Normal==0 || nPts_Inverted==0) return;
  
  gRandom->SetSeed(0);
  
  int i,u;
  unsigned int j,k;
  TH1D *nexp_bkgd = new TH1D("nexp_bkgd","",nBins,-0.5,nBins-0.5);
  TH1D *nexp_signal = new TH1D("nexp_signal","",nBins,-0.5,nBins-0.5);
  TH1D *nexp = new TH1D("nexp","",nBins,-0.5,nBins-0.5);
  double delchi2,chi2min;
  TH1D *chi2hist = new TH1D("chi2hist","",110000,-10,100);
  double ele;
  int noff;
  
  vector< vector<double> > nc,numucc,bnuecc,nutaucc,sig;
  for(j=0;j<nBins;j++)
  {
    nc.push_back( vector<double>() );
    numucc.push_back( vector<double>() );
    bnuecc.push_back( vector<double>() );
    nutaucc.push_back( vector<double>() );
    sig.push_back( vector<double>() );
    if (j<nBinsFHC){
      for(k=0;k<ExtrapFHC.size();k++)
	{
	  nc[j].push_back(0);
	  numucc[j].push_back(0);
	  bnuecc[j].push_back(0);
	  nutaucc[j].push_back(0);
	  sig[j].push_back(0);
	}
    } else if (j>=nBinsFHC){
      for(k=0;k<ExtrapRHC.size();k++)
	{
	  nc[j].push_back(0);
	  numucc[j].push_back(0);
	  bnuecc[j].push_back(0);
	  nutaucc[j].push_back(0);
	  sig[j].push_back(0);
	}
    }
  }
  
  Bkgd = (TH1D*)nexp->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)nexp->Clone("Sig");
  Sig->Reset();
  
  ofstream myfile;
  string file,ofile;
  
  //normal hierarchy
  
  TFile *f = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  
  if(Print)
  {
    ofile = gSystem->ExpandPathName(outFileName.c_str());
    file = ofile.substr(0,ofile.length()-5) + "_Normal.dat";
    myfile.open(gSystem->ExpandPathName(file.c_str()));
  }
  
  for(i=0;i<nPts_Normal;i++)
  {
    cout<<"point "<<(i+1)<<"/"<<nPts_Normal<<" (normal hierarchy)"<<endl;
    
    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();
    
    for(j=0;j<nBins;j++)
    {
      GridTree_Normal[j]->GetEntry(i);
      if (j<nBinsFHC){
	nexp_bkgd->SetBinContent(j+1,grid_background*GridScale_Normal_FHC);
	nexp_signal->SetBinContent(j+1,grid_signal*GridScale_Normal_FHC);
 	for(k=0;k<ExtrapFHC.size();k++)
	  {
	    GridTree_2_Normal[j][k]->GetEntry(i);
	    nc[j][k] = grid_nc*GridScale_Normal_FHC;
	    numucc[j][k] = grid_numucc*GridScale_Normal_FHC;
	    bnuecc[j][k] = grid_bnuecc*GridScale_Normal_FHC;
	    nutaucc[j][k] = grid_nutaucc*GridScale_Normal_FHC;
	    sig[j][k] = grid_nue*GridScale_Normal_FHC;
	  }
      } else if (j>=nBinsFHC){
	nexp_bkgd->SetBinContent(j+1,grid_background*GridScale_Normal_RHC);
	nexp_signal->SetBinContent(j+1,grid_signal*GridScale_Normal_RHC);
 	for(k=0;k<ExtrapRHC.size();k++)
	  {
	    GridTree_2_Normal[j][k]->GetEntry(i);
	    nc[j][k] = grid_nc*GridScale_Normal_RHC;
	    numucc[j][k] = grid_numucc*GridScale_Normal_RHC;
	    bnuecc[j][k] = grid_bnuecc*GridScale_Normal_RHC;
	    nutaucc[j][k] = grid_nutaucc*GridScale_Normal_RHC;
	    sig[j][k] = grid_nue*GridScale_Normal_RHC;
	  }
      }
    }


    nexp->Add(nexp_bkgd,nexp_signal,1,1);
    ErrCalc->SetGridPred(nBinsFHC+nBinsRHC,nc,numucc,bnuecc,nutaucc,sig);
    
    Error4Expts->Reset();
    ErrCalc->SetUseGrid(true);
    ErrCalc->CalculateSystErrorMatrix();
    Error4Expts->Add(ErrCalc->CovMatrix);
    //ErrCalc->CalculateHOOError();
    //Error4Expts->Add(ErrCalc->CovMatrix_Decomp);

    if(IncludeOscParErrs)
    {
      noff=0;
      for(j=0;j<nBins;j++)
      {
	ele=Error4Expts->GetBinContent(j+1,j+1);
        ele+=(grid_bin_oscparerr[j]*grid_bin_oscparerr[j]*nexp->GetBinContent(j+1)*nexp->GetBinContent(j+1));
	Error4Expts->SetBinContent(j+1,j+1,ele);
        
	for(k=0;k<nBins;k++)
	{
	  if(k>j)
	  {
	    ele=Error4Expts->GetBinContent(j+1,k+1);
            ele+=(grid_bin_oscparerr[j]*grid_bin_oscparerr[k]*nexp->GetBinContent(j+1)*nexp->GetBinContent(k+1));
	    Error4Expts->SetBinContent(j+1,k+1,ele);
	    
	    ele=Error4Expts->GetBinContent(k+1,j+1);
            ele+=(grid_bin_oscparerr[j]*grid_bin_oscparerr[k]*nexp->GetBinContent(j+1)*nexp->GetBinContent(k+1));
	    Error4Expts->SetBinContent(k+1,j+1,ele);
            
	    noff++;
	  }
	}
      }
    }
    
    chi2hist->Reset();
    chi2hist->SetName(Form("Chi2Hist_Normal_%i",i));
    
    for(u=0;u<NumExpts;u++)
    {
      cout<<"expt "<<(u+1)<<"/"<<NumExpts<<endl;
      
      GenerateOneCorrelatedExp(nexp,Error4Expts);
      if(Print)
      {
        myfile << grid_sinsq2th13 << " " << grid_delta << " ";
        for(j=0;j<nBins;j++)
        {
          myfile << NObs->GetBinContent(j+1) << " ";
        }
        myfile << endl;
      }
      
      chi2min=GetMinLikelihood(grid_delta,true);
      
      ErrCalc->SetUseGrid(true);//will use the grid predictions set above
      delchi2 = 1e10;
      if(FitMethod==0)
      {
        delchi2 = PoissonChi2(nexp) - chi2min;
      }
      else if(FitMethod==1)
      {
        delchi2 = ScaledChi2(nexp_bkgd,nexp_signal) - chi2min;
      }
      else if(FitMethod==2)
      {
        delchi2 = StandardChi2(nexp) - chi2min;
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd);
        Sig->Reset();
        Sig->Add(nexp_signal);
        delchi2 = StandardLikelihood() - chi2min;
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd);
        Sig->Reset();
        Sig->Add(nexp_signal);

        delchi2 = BinLikelihood() - chi2min;
      }
      else
      {
        cout<<"Error in RunMultiBinPseudoExpts(): Unknown 'FitMethod'."<<endl;
      }
      chi2hist->Fill(delchi2);
      
    }
    f->cd();
    chi2hist->Write();
    f->Close();
    
    f = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"UPDATE");
  }
  
  if(Print) myfile.close();
  
  for(j=0;j<nBins;j++)
  {
    GridTree_Normal[j]->Write();
    
    if (j<nBinsFHC){
      for(k=0;k<ExtrapFHC.size();k++)
	{
	  GridTree_2_Normal[j][k]->Write();
	}
    } else if (j>=nBinsFHC){
      for(k=0;k<ExtrapRHC.size();k++)
	{
	  GridTree_2_Normal[j][k]->Write();
	}
    }
  }
  
  f->Close();
  
  //inverted hierarchy
  
  f = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"UPDATE");
  
  if(Print)
  {
    ofile = gSystem->ExpandPathName(outFileName.c_str());
    file = ofile.substr(0,ofile.length()-5) + "_Inverted.dat";
    myfile.open(gSystem->ExpandPathName(file.c_str()));
  }
  
  for(i=0;i<nPts_Inverted;i++)
  {
    cout<<"point "<<(i+1)<<"/"<<nPts_Inverted<<" (inverted hierarchy)"<<endl;
    
    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();
    
    for(j=0;j<nBins;j++)
    {
      GridTree_Inverted[j]->GetEntry(i);
      if (j<nBinsFHC){
	nexp_bkgd->SetBinContent(j+1,grid_background*GridScale_Inverted_FHC);
	nexp_signal->SetBinContent(j+1,grid_signal*GridScale_Inverted_FHC);
	for(k=0;k<ExtrapFHC.size();k++)
	  {
	    GridTree_2_Inverted[j][k]->GetEntry(i);
	    nc[j][k] = grid_nc*GridScale_Inverted_FHC;
	    numucc[j][k] = grid_numucc*GridScale_Inverted_FHC;
	    bnuecc[j][k] = grid_bnuecc*GridScale_Inverted_FHC;
	    nutaucc[j][k] = grid_nutaucc*GridScale_Inverted_FHC;
	    sig[j][k] = grid_nue*GridScale_Inverted_FHC;
	  }
      } else if (j>=nBinsFHC){
	nexp_bkgd->SetBinContent(j+1,grid_background*GridScale_Inverted_RHC);
	nexp_signal->SetBinContent(j+1,grid_signal*GridScale_Inverted_RHC);
	for(k=0;k<ExtrapRHC.size();k++)
	  {
	    GridTree_2_Inverted[j][k]->GetEntry(i);
	    nc[j][k] = grid_nc*GridScale_Inverted_RHC;
	    numucc[j][k] = grid_numucc*GridScale_Inverted_RHC;
	    bnuecc[j][k] = grid_bnuecc*GridScale_Inverted_RHC;
	    nutaucc[j][k] = grid_nutaucc*GridScale_Inverted_RHC;
	    sig[j][k] = grid_nue*GridScale_Inverted_RHC;
	  }
      }
    }

    nexp->Add(nexp_bkgd,nexp_signal,1,1);
    ErrCalc->SetGridPred(nBinsFHC+nBinsRHC,nc,numucc,bnuecc,nutaucc,sig);
    
    Error4Expts->Reset();
    ErrCalc->SetUseGrid(true);
    ErrCalc->CalculateSystErrorMatrix();
    Error4Expts->Add(ErrCalc->CovMatrix);
    //ErrCalc->CalculateHOOError();
    //Error4Expts->Add(ErrCalc->CovMatrix_Decomp);
    if(IncludeOscParErrs)
    {
      noff=0;
      for(j=0;j<nBins;j++)
      {
        ele=Error4Expts->GetBinContent(j+1,j+1);
        ele+=(grid_bin_oscparerr[j]*grid_bin_oscparerr[j]*nexp->GetBinContent(j+1)*nexp->GetBinContent(j+1));
        Error4Expts->SetBinContent(j+1,j+1,ele);
        
        for(k=0;k<nBins;k++)
        {
          if(k>j)
          {
            ele=Error4Expts->GetBinContent(j+1,k+1);
            ele+=(grid_bin_oscparerr[j]*grid_bin_oscparerr[k]*nexp->GetBinContent(j+1)*nexp->GetBinContent(k+1));
            Error4Expts->SetBinContent(j+1,k+1,ele);
            
            ele=Error4Expts->GetBinContent(k+1,j+1);
            ele+=(grid_bin_oscparerr[j]*grid_bin_oscparerr[k]*nexp->GetBinContent(j+1)*nexp->GetBinContent(k+1));
            Error4Expts->SetBinContent(k+1,j+1,ele);
            
            noff++;
          }
        }
      }
    }
    
    chi2hist->Reset();
    chi2hist->SetName(Form("Chi2Hist_Inverted_%i",i));
    
    for(u=0;u<NumExpts;u++)
    {
      cout<<"expt "<<(u+1)<<"/"<<NumExpts<<endl;
      
      GenerateOneCorrelatedExp(nexp,Error4Expts);
      if(Print)
      {
        myfile << grid_sinsq2th13 << " " << grid_delta << " ";
        for(j=0;j<nBins;j++)
        {
          myfile << NObs->GetBinContent(j+1) << " ";
        }
        myfile << endl;
      }
      
      chi2min=GetMinLikelihood(grid_delta,false);
      
      ErrCalc->SetUseGrid(true);//will use the grid predictions set above
      delchi2 = 1e10;
      if(FitMethod==0)
      {
        delchi2 = PoissonChi2(nexp) - chi2min;
      }
      else if(FitMethod==1)
      {
        delchi2 = ScaledChi2(nexp_bkgd,nexp_signal) - chi2min;
      }
      else if(FitMethod==2)
      {
        delchi2 = StandardChi2(nexp) - chi2min;
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd);
        Sig->Reset();
        Sig->Add(nexp_signal);
        delchi2 = StandardLikelihood() - chi2min;
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd);
        Sig->Reset();
        Sig->Add(nexp_signal);
        delchi2 = BinLikelihood() - chi2min;
      }
      else
      {
        cout<<"Error in RunMultiBinPseudoExpts(): Unknown 'FitMethod'."<<endl;
      }
      chi2hist->Fill(delchi2);
    }
    f->cd();
    chi2hist->Write();
    f->Close();
    f = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"UPDATE");
  }
  
  if(Print) myfile.close();
  
  for(j=0;j<nBins;j++)
  {
    GridTree_Inverted[j]->Write();
    if (j<nBinsFHC){
      for(k=0;k<ExtrapFHC.size();k++)
	{
	  GridTree_2_Inverted[j][k]->Write();
	}
    } else if (j>=nBinsFHC){
      for(k=0;k<ExtrapRHC.size();k++)
	{
	  GridTree_2_Inverted[j][k]->Write();
	}
    }
  }
  
  f->Close();
  
  return;
}

double NueFit2D_Joint::GetSensitivityAt(double delta, bool normalhier)
{

  if(NObs==0)
  {
    cout<<"NObs not set.  Quitting..."<<endl;
    return 0;
  }
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
  {
    cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
    return 0;
  }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return 0;
  }
  
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();
  
  Bkgd = (TH1D*)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)NObs->Clone("Sig");
  Sig->Reset();
  TH1D *NExp = (TH1D*)NObs->Clone("NExp");
  NExp->Reset();
  
  unsigned int ie;
  for(ie=0;ie<ExtrapFHC.size();ie++) ExtrapFHC[ie]->GetPrediction();
  for(ie=0;ie<ExtrapRHC.size();ie++) ExtrapRHC[ie]->GetPrediction();

  Int_t i;
  
  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  Double_t ss2th13 = 0;
  Double_t delchi2 = 0,delchi2prev = 0;
  Double_t contourlvl = 2.71;
  Double_t prev = 0;
  Double_t best = 0;
  Double_t Th13increment = 0;
  if(nSinSq2Th13Steps>0) Th13increment = (SinSq2Th13High - SinSq2Th13Low)/(nSinSq2Th13Steps);
  double sens=-1;
  double chi2[3],val[3];
  TGraph *g3;
  TF1* fit;
  double minchi2;
  double limit = 0;
  
  for(ie=0;ie<ExtrapRHC.size();ie++)
  {
    ExtrapRHC[ie]->SetDeltaCP(delta);
    if(!normalhier) ExtrapRHC[ie]->InvertMassHierarchy();
  }

  for(ie=0;ie<ExtrapFHC.size();ie++)
  {
    ExtrapFHC[ie]->SetDeltaCP(delta);
    if(!normalhier) ExtrapFHC[ie]->InvertMassHierarchy();
  }
  
  contourlvl = 2.71;
  
  best=-1.;
  minchi2=100;
  for(i=0;i<3;i++)
  {
    chi2[i]=-1.;
    val[i]=-1.;
  }
  for(i=0;i<nSinSq2Th13Steps+1;i++)
  {
    chi2[0] = chi2[1];
    chi2[1] = chi2[2];
    val[0] = val[1];
    val[1] = val[2];
    ss2th13 = i*Th13increment + SinSq2Th13Low;

    for(ie=0;ie<ExtrapFHC.size();ie++)
    {
      ExtrapFHC[ie]->SetSinSq2Th13(ss2th13);
      ExtrapFHC[ie]->OscillatePrediction();
    }

    for(ie=0;ie<ExtrapRHC.size();ie++)
    {
      ExtrapRHC[ie]->SetSinSq2Th13(ss2th13);
      ExtrapRHC[ie]->OscillatePrediction();
    }

    Bkgd->Reset();
    Sig->Reset();
    NExp->Reset();

    //Turn extrapolations into a prediction:
    CombineFHCRHC();
    NExp->Add(Bkgd);
    NExp->Add(Sig);
    
    chi2[2] = 1e10;
    if(FitMethod==0)
    {
      chi2[2] = PoissonChi2(NExp);
    }
    else if(FitMethod==1)
    {
      chi2[2] = ScaledChi2(Bkgd,Sig);
    }
    else if(FitMethod==2)
    {
      chi2[2] = StandardChi2(NExp);
    }
    else if(FitMethod==3)
    {
      //Likelihood: "Standard" (N syst, N nuisance)
      //Calculate the likelihood (x2 for chi)
      chi2[2] = StandardLikelihood();
    }
    else if(FitMethod==4)
    {
      //Likelihood: Bin by Bin Calculation of Systematics
      //Calculate the likelihood (x2 for chi)
      chi2[2] = BinLikelihood();
      //cout << ss2th13 << " -> " << chi2[2] << endl;
    }
    else
    {
      cout<<"Error in GetSensitivityAt(): Unknown 'FitMethod'."<<endl;
    }
    
    val[2] = ss2th13;
    
    if(i<2) continue;
    
    if(i<3 && chi2[2]>chi2[1] && chi2[1]>chi2[0])//first three points are increasing, first point is minimum.
    {
      best = val[0];
      minchi2 = chi2[0];
      cout<<"minimum at 1st point: "<<minchi2<<" at "<<best<<endl;
      break; //CHANGE BACK
    }
    
    if(chi2[2]>chi2[1] && chi2[0]>chi2[1])//found minimum
    {
      g3 = new TGraph(3, val, chi2);
      fit = new TF1("pol2", "pol2");
      g3->Fit(fit, "Q");//fit to second order polynominal
      if(fit->GetParameter(2) > 0)//if the x^2 term is nonzero
      {
        best = -fit->GetParameter(1)/(2*fit->GetParameter(2));//the location of the minimum is -p1/(2*p2)
        minchi2 = fit->GetParameter(0) + fit->GetParameter(1)*best + fit->GetParameter(2)*best*best;
        cout<<"minimum with fit: "<<minchi2<<" at "<<best<<endl;
      }
      else//if the x^2 term is zero, then just use the minimum you got by scanning
      {
        best = val[1];
        minchi2 = chi2[1];
        cout<<"minimum with scan: "<<minchi2<<" at "<<best<<endl;
      }
      break; //CHANGE BACK
    }

    //if (ss2th13>0.11) break; //DELETE ME
  }
  
  limit = -1;
  delchi2prev = 1000;
  prev = 0;
  for(i=0;i<nSinSq2Th13Steps+1;i++)
  {
    ss2th13 = i*Th13increment + SinSq2Th13Low;
    for(ie=0;ie<ExtrapFHC.size();ie++)
    {
      ExtrapFHC[ie]->SetSinSq2Th13(ss2th13);
      ExtrapFHC[ie]->OscillatePrediction();
    }
    for(ie=0;ie<ExtrapRHC.size();ie++)
    {
      ExtrapRHC[ie]->SetSinSq2Th13(ss2th13);
      ExtrapRHC[ie]->OscillatePrediction();
    }
    
    Bkgd->Reset();
    Sig->Reset();
    NExp->Reset();
    
    CombineFHCRHC();
    NExp->Add(Bkgd);
    NExp->Add(Sig);
    
    delchi2 = 1e10;
    if(FitMethod==0)
    {
      delchi2 = PoissonChi2(NExp) - minchi2;
    }
    else if(FitMethod==1)
    {
      delchi2 = ScaledChi2(Bkgd,Sig) - minchi2;
    }
    else if(FitMethod==2)
    {
      delchi2 = StandardChi2(NExp) - minchi2;
    }
    else if(FitMethod==3)
    {
      //Likelihood: "Standard" (N syst, N nuisance)
      //Calculate the likelihood (x2 for chi)
      delchi2 = StandardLikelihood() - minchi2;
    }
    else if(FitMethod==4)
    {
      //Likelihood: Bin by Bin Calculation of Systematics
      //Calculate the likelihood (x2 for chi)
      delchi2 = BinLikelihood() - minchi2;
      //cout << "ss2th13 = " << ss2th13 << " Chi2=" << delchi2 << endl;
    }
    else
    {
      cout<<"Error in GetSensitivityAt(): Unknown 'FitMethod'."<<endl;
    }

    //RBT: How to run the separate options

//     if (opt==0){
//       //Chi2 standard
//       chi2 = StandardChi2(NExp);
//     } else if (opt==1){
//       //Likelihood: Bin by Bin Calculation of Systematics
//       //Calculate the likelihood (x2 for chi)
//       chi2 = BinLikelihood();
//     } else if (opt==2){
//       //Likelihood: "Standard" (N syst, N nuisance)
//       //Calculate the likelihood (x2 for chi)
//       chi2 = StandardLikelihood();
//     }

    if(delchi2>contourlvl && delchi2prev<contourlvl && TMath::Abs(delchi2prev-delchi2)<0.1)
    {
      limit = ss2th13 + ((ss2th13-prev)/(delchi2 - delchi2prev))*(contourlvl - delchi2);
      sens = limit;
      break;
    }
    delchi2prev = delchi2;
    prev = ss2th13;
  }
  
  return sens;
}

void NueFit2D_Joint::RunDeltaChi2Contour(int cl)
{
  if(NObs==0)
  {
    cout<<"NObs not set.  Quitting..."<<endl;
    return;
  }
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
  {
    cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
    return;
  }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return;
  }
  
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();
  
  Bkgd = (TH1D*)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)NObs->Clone("Sig");
  Sig->Reset();
  NExp = (TH1D*)NObs->Clone("NExp");
  NExp->Reset();
  
  unsigned int ie;
  for(ie=0;ie<ExtrapFHC.size();ie++) ExtrapFHC[ie]->GetPrediction();
  for(ie=0;ie<ExtrapRHC.size();ie++) ExtrapRHC[ie]->GetPrediction();
  
  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  Int_t i,j;
  
  Double_t ss2th13 = 0;
  Double_t delta = 0;
  Double_t Deltaincrement = 0;
  if(nDeltaSteps>0) Deltaincrement = (DeltaHigh - DeltaLow)/(nDeltaSteps);
  Double_t Th13increment = 0;
  if(nSinSq2Th13Steps>0) Th13increment = (SinSq2Th13High - SinSq2Th13Low)/(nSinSq2Th13Steps);
  
  double chi2[3],val[3];
  TGraph *g3;
  TF1* fit;
  double best;
  double minchi2;
  
  double limit;
  double delchi2;
  double sprev,delchi2prev;
  double contourlvl = 0;
  if(cl==0)//90% CL
  {
    contourlvl = 2.71;
  }
  else if(cl==1)//68.3% cL
  {
    contourlvl = 1.0;
  }
  else
  {
    cout<<"Error in RunDeltaChi2Contour(): Input value should be 0 or 1 for 90% or 68.3%.  Quitting..."<<endl;
    return;
  }
  
  cout<<"Seeking ";
  if(cl==0) cout<<"90% ";
  else cout<<"68% ";
  cout<<" CL upper limit"<<endl;
  
  vector<double> deltapts;
  vector<double> bestfit_norm;
  vector<double> bestfit_invt;
  vector<double> limit_norm;
  vector<double> limit_invt;
  
  for(j=0;j<nDeltaSteps+1;j++)
  {
    delta = j*Deltaincrement*TMath::Pi() + DeltaLow;
    for(ie=0;ie<ExtrapFHC.size();ie++) ExtrapFHC[ie]->SetDeltaCP(delta);
    for(ie=0;ie<ExtrapRHC.size();ie++) ExtrapRHC[ie]->SetDeltaCP(delta);

    deltapts.push_back(delta/TMath::Pi());
    
    best=-1.;
    minchi2=100;
    for(i=0;i<3;i++)
    {
      chi2[i]=-1.;
      val[i]=-1.;
    }
    for(i=0;i<nSinSq2Th13Steps+1;i++)
    {
      chi2[0] = chi2[1];
      chi2[1] = chi2[2];
      val[0] = val[1];
      val[1] = val[2];
      ss2th13 = i*Th13increment + SinSq2Th13Low;

      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetSinSq2Th13(ss2th13);
        ExtrapFHC[ie]->OscillatePrediction();
      }
      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetSinSq2Th13(ss2th13);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      chi2[2] = 1e10;
      if(FitMethod==0)
      {
        chi2[2] = PoissonChi2(NExp);
      }
      else if(FitMethod==1)
      {
        chi2[2] = ScaledChi2(Bkgd,Sig);
      }
      else if(FitMethod==2)
      {
        chi2[2] = StandardChi2(NExp);
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        chi2[2] = StandardLikelihood();
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        chi2[2] = BinLikelihood();
      }
      else
      {
        cout<<"Error in RunDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
      
      val[2] = ss2th13;
      
      if(i<2) continue;
      
      if(i<3 && chi2[2]>chi2[1] && chi2[1]>chi2[0])//first three points are increasing, first point is minimum.
      {
	best = val[0];
	minchi2 = chi2[0];
	cout<<"minimum at 1st point: "<<minchi2<<" at "<<best<<endl;
	break;
      }
      
      if(chi2[2]>chi2[1] && chi2[0]>chi2[1])//found minimum
      {
        g3 = new TGraph(3, val, chi2);
        fit = new TF1("pol2", "pol2");
        g3->Fit(fit, "Q");//fit to second order polynominal
        if(fit->GetParameter(2) > 0)//if the x^2 term is nonzero
        {
          best = -fit->GetParameter(1)/(2*fit->GetParameter(2));//the location of the minimum is -p1/(2*p2)
          minchi2 = fit->GetParameter(0) + fit->GetParameter(1)*best + fit->GetParameter(2)*best*best;
          cout<<"minimum with fit: "<<minchi2<<" at "<<best<<endl;
        }
        else//if the x^2 term is zero, then just use the minimum you got by scanning
        {
          best = val[1];
          minchi2 = chi2[1];
          cout<<"minimum with scan: "<<minchi2<<" at "<<best<<endl;
        }
        break;
      }
    }
    bestfit_norm.push_back(best);
    
    limit = -1.;
    delchi2prev = 1000;
    sprev = 0;
    for(i=0;i<nSinSq2Th13Steps+1;i++)
    {
      ss2th13 = i*Th13increment + SinSq2Th13Low;

      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetSinSq2Th13(ss2th13);
        ExtrapFHC[ie]->OscillatePrediction();
      }
      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetSinSq2Th13(ss2th13);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      delchi2 = 1e10;
      if(FitMethod==0)
      {
        delchi2 = PoissonChi2(NExp) - minchi2;
      }
      else if(FitMethod==1)
      {
        delchi2 = ScaledChi2(Bkgd,Sig) - minchi2;
      }
      else if(FitMethod==2)
      {
        delchi2 = StandardChi2(NExp) - minchi2;
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        delchi2 = StandardLikelihood() - minchi2;
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        delchi2 = BinLikelihood() - minchi2;
      }
      else
      {
        cout<<"Error in RunDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
      
      if(i==1) continue;
      
      if(delchi2>contourlvl && delchi2prev<contourlvl)
      {
        limit = ss2th13 + ((ss2th13-sprev)/(delchi2 - delchi2prev))*(contourlvl - delchi2);
        cout<<delchi2prev<<", "<<sprev<<", "<<delchi2<<", "<<ss2th13<<", "<<limit<<endl;
        break;
      }
      delchi2prev = delchi2;
      sprev = ss2th13;
    }
    limit_norm.push_back(limit);
  }
  
  for(ie=0;ie<ExtrapFHC.size();ie++) ExtrapFHC[ie]->InvertMassHierarchy();
  for(ie=0;ie<ExtrapRHC.size();ie++) ExtrapRHC[ie]->InvertMassHierarchy();

  
  for(j=0;j<nDeltaSteps+1;j++)
  {
    delta = j*Deltaincrement*TMath::Pi() + DeltaLow;
    for(ie=0;ie<ExtrapFHC.size();ie++) ExtrapFHC[ie]->SetDeltaCP(delta);
    for(ie=0;ie<ExtrapRHC.size();ie++) ExtrapRHC[ie]->SetDeltaCP(delta);
    
    best=-1.;
    minchi2=100;
    for(i=0;i<3;i++)
    {
      chi2[i]=-1.;
      val[i]=-1.;
    }
    for(i=0;i<nSinSq2Th13Steps+1;i++)
    {
      chi2[0] = chi2[1];
      chi2[1] = chi2[2];
      val[0] = val[1];
      val[1] = val[2];
      ss2th13 = i*Th13increment + SinSq2Th13Low;

      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetSinSq2Th13(ss2th13);
        ExtrapFHC[ie]->OscillatePrediction();
      }

      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetSinSq2Th13(ss2th13);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      chi2[2] = 1e10;
      if(FitMethod==0)
      {
        chi2[2] = PoissonChi2(NExp);
      }
      else if(FitMethod==1)
      {
        chi2[2] = ScaledChi2(Bkgd,Sig);
      }
      else if(FitMethod==2)
      {
        chi2[2] = StandardChi2(NExp);
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        chi2[2] = StandardLikelihood();
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        chi2[2] = BinLikelihood();
      }
      else
      {
        cout<<"Error in RunDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
      
      val[2] = ss2th13;
      
      if(i<2) continue;
      
      if(i<3 && chi2[2]>chi2[1] && chi2[1]>chi2[0])//first three points are increasing, first point is minimum.
      {
	best = val[0];
	minchi2 = chi2[0];
	cout<<"minimum at 1st point: "<<minchi2<<" at "<<best<<endl;
	break;
      }
      
      if(chi2[2]>chi2[1] && chi2[0]>chi2[1])//found minimum
      {
        g3 = new TGraph(3, val, chi2);
        fit = new TF1("pol2", "pol2");
        g3->Fit(fit, "Q");//fit to second order polynominal
        if(fit->GetParameter(2) > 0)//if the x^2 term is nonzero
        {
          best = -fit->GetParameter(1)/(2*fit->GetParameter(2));//the location of the minimum is -p1/(2*p2)
          minchi2 = fit->GetParameter(0) + fit->GetParameter(1)*best + fit->GetParameter(2)*best*best;
          cout<<"minimum with fit: "<<minchi2<<" at "<<best<<endl;
        }
        else//if the x^2 term is zero, then just use the minimum you got by scanning
        {
          best = val[1];
          minchi2 = chi2[1];
          cout<<"minimum with scan: "<<minchi2<<" at "<<best<<endl;
        }
        break;
      }
    }
    bestfit_invt.push_back(best);
    
    limit = -1.;
    delchi2prev = 1000;
    sprev = 0;
    for(i=0;i<nSinSq2Th13Steps+1;i++)
    {
      ss2th13 = i*Th13increment + SinSq2Th13Low;

      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetSinSq2Th13(ss2th13);
        ExtrapFHC[ie]->OscillatePrediction();
      }
      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetSinSq2Th13(ss2th13);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      delchi2 = 1e10;
      if(FitMethod==0)
      {
        delchi2 = PoissonChi2(NExp) - minchi2;
      }
      else if(FitMethod==1)
      {
        delchi2 = ScaledChi2(Bkgd,Sig) - minchi2;
      }
      else if(FitMethod==2)
      {
        delchi2 = StandardChi2(NExp) - minchi2;
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        delchi2 = StandardLikelihood() - minchi2;
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        delchi2 = BinLikelihood() - minchi2;
      }
      else
      {
        cout<<"Error in RunDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
      
      if(i==1) continue;
      
      if(delchi2>contourlvl && delchi2prev<contourlvl)
      {
        limit = ss2th13 + ((ss2th13-sprev)/(delchi2 - delchi2prev))*(contourlvl - delchi2);
        cout<<delchi2prev<<", "<<sprev<<", "<<delchi2<<", "<<ss2th13<<", "<<limit<<endl;
        break;
      }
      delchi2prev = delchi2;
      sprev = ss2th13;
    }
    limit_invt.push_back(limit);
  }
  
  double *s_best_n = new double[nDeltaSteps+1];
  double *s_best_i = new double[nDeltaSteps+1];
  double *s_limit_n = new double[nDeltaSteps+1];
  double *s_limit_i = new double[nDeltaSteps+1];
  double *d = new double[nDeltaSteps+1];
  for(i=0;i<nDeltaSteps+1;i++)
  {
    d[i] = deltapts.at(i);
    s_best_n[i] = bestfit_norm.at(i);
    s_best_i[i] = bestfit_invt.at(i);
    s_limit_n[i] = limit_norm.at(i);
    s_limit_i[i] = limit_invt.at(i);
  }
  
  TGraph *gn = new TGraph(nDeltaSteps+1,s_best_n,d);
  gn->SetMarkerStyle(20);
  gn->SetTitle("");
  gn->GetYaxis()->SetTitle("#delta");
  gn->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gn->GetXaxis()->SetLimits(0,0.6);
  gn->SetLineWidth(4);
  gn->SetMaximum(2);
  gn->SetName("BestFit_Normal");
  
  TGraph *gi = new TGraph(nDeltaSteps+1,s_best_i,d);
  gi->SetMarkerStyle(20);
  gi->SetTitle("");
  gi->GetYaxis()->SetTitle("#delta");
  gi->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gi->GetXaxis()->SetLimits(0,0.6);
  gi->SetLineWidth(4);
  gi->SetLineStyle(2);
  gi->SetMaximum(2);
  gi->SetName("BestFit_Inverted");
  
  TGraph *gn_limit = new TGraph(nDeltaSteps+1,s_limit_n,d);
  gn_limit->SetMarkerStyle(20);
  gn_limit->SetTitle("");
  gn_limit->GetYaxis()->SetTitle("#delta");
  gn_limit->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gn_limit->GetXaxis()->SetLimits(0,0.6);
  gn_limit->SetLineWidth(4);
  gn_limit->SetLineColor(kBlue);
  gn_limit->SetMarkerColor(kBlue);
  gn_limit->SetMaximum(2);
  gn_limit->SetName("Limit_Normal");
  
  TGraph *gi_limit = new TGraph(nDeltaSteps+1,s_limit_i,d);
  gi_limit->SetMarkerStyle(20);
  gi_limit->SetTitle("");
  gi_limit->GetYaxis()->SetTitle("#delta");
  gi_limit->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  gi_limit->GetXaxis()->SetLimits(0,0.6);
  gi_limit->SetLineWidth(4);
  gi_limit->SetLineColor(kRed);
  gi_limit->SetMarkerColor(kRed);
  gi_limit->SetMaximum(2);
  gi_limit->SetName("Limit_Inverted");
  
  if(cl==0) cout<<"90% ";
  if(cl==1) cout<<"68% ";
  cout<<"confidence level limit = "<<limit_norm.at(0)<<", "<<limit_invt.at(0)<<endl;
  
  TFile *fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  gn->Write();
  gi->Write();
  gn_limit->Write();
  gi_limit->Write();
  fout->Write();
  fout->Close();
  
  delete [] s_best_n;
  delete [] s_best_i;
  delete [] s_limit_n;
  delete [] s_limit_i;
  delete [] d;
  
  return;
}

double NueFit2D_Joint::GetLikelihood(double t12,double t23,double t13,double dm2_32,double dm2_21,double delta)
{
  if(NObs==0)
  {
    cout<<"NObs not set.  Quitting..."<<endl;
    return 0;
  }
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
  {
    cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
    return 0;
  }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return 0;
  }
  
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();
  
  Bkgd = (TH1D*)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)NObs->Clone("Sig");
  Sig->Reset();
  TH1D *NExp = (TH1D*)NObs->Clone("NExp");
  NExp->Reset();
  
  unsigned int ie;
  for(ie=0;ie<ExtrapFHC.size();ie++)
  {
    ExtrapFHC[ie]->GetPrediction();
    ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,t12);
    ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,t23);
    ExtrapFHC[ie]->SetOscPar(OscPar::kTh13,t13);
    ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,dm2_32);
    ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,dm2_21);
    ExtrapFHC[ie]->SetOscPar(OscPar::kDelta,delta);
    ExtrapFHC[ie]->OscillatePrediction();
  }
  for(ie=0;ie<ExtrapRHC.size();ie++)
  {
    ExtrapRHC[ie]->GetPrediction();
    ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,t12);
    ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,t23);
    ExtrapRHC[ie]->SetOscPar(OscPar::kTh13,t13);
    ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,dm2_32);
    ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,dm2_21);
    ExtrapRHC[ie]->SetOscPar(OscPar::kDelta,delta);
    ExtrapRHC[ie]->OscillatePrediction();
  }
  
  Int_t i;
  
  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  Double_t ss2th13 = 0;
  Double_t delchi2 = 0;
  Double_t best = 0;
  Double_t Th13increment = 0;
  if(nSinSq2Th13Steps>0) Th13increment = (SinSq2Th13High - SinSq2Th13Low)/(nSinSq2Th13Steps);
  double chi2[3],val[3];
  TGraph *g3;
  TF1* fit;
  double minchi2;
  
  best=-1.;
  minchi2=100;
  for(i=0;i<3;i++)
  {
    chi2[i]=-1.;
    val[i]=-1.;
  }
  for(i=0;i<nSinSq2Th13Steps+1;i++)
  {
    chi2[0] = chi2[1];
    chi2[1] = chi2[2];
    val[0] = val[1];
    val[1] = val[2];
    ss2th13 = i*Th13increment + SinSq2Th13Low;

    for(ie=0;ie<ExtrapFHC.size();ie++)
    {
      ExtrapFHC[ie]->SetSinSq2Th13(ss2th13);
      ExtrapFHC[ie]->OscillatePrediction();
    }

    for(ie=0;ie<ExtrapRHC.size();ie++)
    {
      ExtrapRHC[ie]->SetSinSq2Th13(ss2th13);
      ExtrapRHC[ie]->OscillatePrediction();
    }

    Bkgd->Reset();
    Sig->Reset();
    NExp->Reset();

    CombineFHCRHC();
    NExp->Add(Bkgd);
    NExp->Add(Sig);
    
    chi2[2] = 1e10;
    if(FitMethod==0)
    {
      chi2[2] = PoissonChi2(NExp);
    }
    else if(FitMethod==1)
    {
      chi2[2] = ScaledChi2(Bkgd,Sig);
    }
    else if(FitMethod==2)
    {
      chi2[2] = StandardChi2(NExp);
    }
    else if(FitMethod==3)
    {
      //Likelihood: "Standard" (N syst, N nuisance)
      //Calculate the likelihood (x2 for chi)
      chi2[2] = StandardLikelihood();
    }
    else if(FitMethod==4)
    {
      //Likelihood: Bin by Bin Calculation of Systematics
      //Calculate the likelihood (x2 for chi)
      chi2[2] = BinLikelihood();
      //cout << "ss2th13 = " << ss2th13 << " Chi2=" << chi2[2] << endl;
    }
    else
    {
      cout<<"Error in GetLikelihood(): Unknown 'FitMethod'."<<endl;
    }
    
    val[2] = ss2th13;
    
    if(i<2) continue;
    
    if(i<3 && chi2[2]>chi2[1] && chi2[1]>chi2[0])//first three points are increasing, first point is minimum.
    {
      best = val[0];
      minchi2 = chi2[0];
      cout<<"minimum at 1st point: "<<minchi2<<" at "<<best<<endl;
      break;
    }
    
    if(chi2[2]>chi2[1] && chi2[0]>chi2[1])//found minimum
    {
      g3 = new TGraph(3, val, chi2);
      fit = new TF1("pol2", "pol2");
      g3->Fit(fit, "Q");//fit to second order polynominal
      if(fit->GetParameter(2) > 0)//if the x^2 term is nonzero
      {
        best = -fit->GetParameter(1)/(2*fit->GetParameter(2));//the location of the minimum is -p1/(2*p2)
        minchi2 = fit->GetParameter(0) + fit->GetParameter(1)*best + fit->GetParameter(2)*best*best;
        cout<<"minimum with fit: "<<minchi2<<" at "<<best<<endl;
      }
      else//if the x^2 term is zero, then just use the minimum you got by scanning
      {
        best = val[1];
        minchi2 = chi2[1];
        cout<<"minimum with scan: "<<minchi2<<" at "<<best<<endl;
      }
      break;
    }
  }
  
  for(ie=0;ie<ExtrapFHC.size();ie++)
  {
    ExtrapFHC[ie]->SetOscPar(OscPar::kTh13,t13);
    ExtrapFHC[ie]->OscillatePrediction();
  }
  for(ie=0;ie<ExtrapRHC.size();ie++)
  {
    ExtrapRHC[ie]->SetOscPar(OscPar::kTh13,t13);
    ExtrapRHC[ie]->OscillatePrediction();
  }

  Bkgd->Reset();
  Sig->Reset();
  NExp->Reset();

  CombineFHCRHC();
  NExp->Add(Bkgd);
  NExp->Add(Sig);
  
  delchi2 = 1e10;
  if(FitMethod==0)
  {
    delchi2 = PoissonChi2(NExp) - minchi2;
  }
  else if(FitMethod==1)
  {
    delchi2 = ScaledChi2(Bkgd,Sig) - minchi2;
  }
  else if(FitMethod==2)
  {
    delchi2 = StandardChi2(NExp) - minchi2;
  }
  else if(FitMethod==3)
  {
    //Likelihood: "Standard" (N syst, N nuisance)
    //Calculate the likelihood (x2 for chi)
    delchi2 = StandardLikelihood() - minchi2;
  }
  else if(FitMethod==4)
  {
    //Likelihood: Bin by Bin Calculation of Systematics
    //Calculate the likelihood (x2 for chi)
    delchi2 = BinLikelihood() - minchi2;
  }
  else
  {
    cout<<"Error in GetLikelihood(): Unknown 'FitMethod'."<<endl;
  }
  
  cout<<"delchi2 = "<<delchi2<<endl;
  
  return delchi2;
}

double NueFit2D_Joint::GetMinLikelihood(double delta,bool normalhier)
{
  if(NObs==0)
  {
    cout<<"NObs not set.  Quitting..."<<endl;
    return 0;
  }
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
  {
    cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
    return 0;
  }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return 0;
  }
 
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();
  if(FitMethod==3 || FitMethod==4)
  {
    if(ErrCalc!=0) ErrCalc->SetUseGrid(false);
  }
  
  Bkgd = (TH1D*)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)NObs->Clone("Sig");
  Sig->Reset();
  TH1D *NExp = (TH1D*)NObs->Clone("NExp");
  NExp->Reset();
  
  int nBinsFHC=0;
  int nBinsRHC=0;
  if (ExtrapFHC.size()>0) nBinsFHC = ExtrapFHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ExtrapRHC.size()>0) nBinsRHC = ExtrapRHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  nBins = nBinsFHC + nBinsRHC;  
  
  unsigned int ie;
  for(ie=0;ie<ExtrapFHC.size();ie++)
  {

    ExtrapFHC[ie]->GetPrediction();
    if(normalhier)
    {
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,grid_n_th12);
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,grid_n_th23);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,grid_n_dm2_21);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,grid_n_dm2_32);
    }
    else
    {
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,grid_i_th12);
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,grid_i_th23);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,grid_i_dm2_21);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,grid_i_dm2_32);
    }
    ExtrapFHC[ie]->SetOscPar(OscPar::kDelta,delta);
    ExtrapFHC[ie]->OscillatePrediction();
  }

  for(ie=0;ie<ExtrapRHC.size();ie++)
  {

    ExtrapRHC[ie]->GetPrediction();
    if(normalhier)
    {
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,grid_n_th12);
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,grid_n_th23);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,grid_n_dm2_21);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,grid_n_dm2_32);
    }
    else
    {
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,grid_i_th12);
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,grid_i_th23);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,grid_i_dm2_21);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,grid_i_dm2_32);
    }
    ExtrapRHC[ie]->SetOscPar(OscPar::kDelta,delta);
    ExtrapRHC[ie]->OscillatePrediction();
  }
  Int_t i;
  
  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  Double_t ss2th13 = 0;
  TF1* fit;
  double minchi2;
  
  const int ns=30;
  double max=0.6;
  double width = max/ns;
  double th13[ns],c2[ns];
  for(i=0;i<ns;i++)
  {

    ss2th13 = i*width;
    th13[i] = ss2th13;

    for(ie=0;ie<ExtrapFHC.size();ie++)
      {
	ExtrapFHC[ie]->SetSinSq2Th13(ss2th13);
	ExtrapFHC[ie]->OscillatePrediction();
      }
    for(ie=0;ie<ExtrapRHC.size();ie++)
      {
	ExtrapRHC[ie]->SetSinSq2Th13(ss2th13);
	ExtrapRHC[ie]->OscillatePrediction();
      }
    Bkgd->Reset();
    Sig->Reset();
    NExp->Reset();
 
  //Turn extrapolations into a prediction:
    CombineFHCRHC();
    NExp->Add(Bkgd);
    NExp->Add(Sig);
 
    c2[i] = 1e10;
    if(FitMethod==0)
    {
      c2[i] = PoissonChi2(NExp);
    }
    else if(FitMethod==1)
    {
      c2[i] = ScaledChi2(Bkgd,Sig);
    }
    else if(FitMethod==2)
    {
      c2[i] = StandardChi2(NExp);
    }
    else if(FitMethod==3)
    {
      //Likelihood: "Standard" (N syst, N nuisance)
      //Calculate the likelihood (x2 for chi)
      c2[i] = StandardLikelihood();
    }
    else if(FitMethod==4)
    {
      //Likelihood: Bin by Bin Calculation of Systematics
      //Calculate the likelihood (x2 for chi)
      c2[i] = BinLikelihood();
    }
    else
    {
      cout<<"Error in GetMinLikelihood(): Unknown 'FitMethod'."<<endl;
    }
  }
  
  TGraph *g = new TGraph(ns,th13,c2);
  fit = new TF1("pol6", "pol6");
  g->Fit(fit, "Q");//fit to second order polynominal
  minchi2 = fit->GetMinimum(0,max);
  
  return minchi2;
}
double NueFit2D_Joint::GetMinLikelihood_Delta(bool normalhier)
{
  if(NObs==0)
  {
    cout<<"NObs not set.  Quitting..."<<endl;
    return 0;
  }
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
  {
    cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
    return 0;
  }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return 0;
  }
  
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();
  if(FitMethod==3 || FitMethod==4)
  {
    if(ErrCalc!=0) ErrCalc->SetUseGrid(false);
  }
  
  Bkgd = (TH1D*)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)NObs->Clone("Sig");
  Sig->Reset();
  TH1D *NExp = (TH1D*)NObs->Clone("NExp");
  NExp->Reset();
  
  int nBinsFHC=0;
  int nBinsRHC=0;
  if (ExtrapFHC.size()>0) nBinsFHC = ExtrapFHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ExtrapRHC.size()>0) nBinsRHC = ExtrapRHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  nBins = nBinsFHC + nBinsRHC;  
  
  unsigned int ie;
  for(ie=0;ie<ExtrapFHC.size();ie++)
  {
    
    ExtrapFHC[ie]->GetPrediction();
    if(normalhier)
    {
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh13,grid_n_th13);
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,grid_n_th12);
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,grid_n_th23);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,grid_n_dm2_21);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,grid_n_dm2_32);
    }
    else
    {
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh13,grid_i_th13);
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,grid_i_th12);
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,grid_i_th23);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,grid_i_dm2_21);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,grid_i_dm2_32);
    }
    ExtrapFHC[ie]->OscillatePrediction();
  }
  
  for(ie=0;ie<ExtrapRHC.size();ie++)
  {
    
    ExtrapRHC[ie]->GetPrediction();
    if(normalhier)
    {
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh13,grid_n_th13);
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,grid_n_th12);
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,grid_n_th23);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,grid_n_dm2_21);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,grid_n_dm2_32);
    }
    else
    {
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh13,grid_i_th13);
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,grid_i_th12);
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,grid_i_th23);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,grid_i_dm2_21);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,grid_i_dm2_32);
    }
    ExtrapRHC[ie]->OscillatePrediction();
  }
  Int_t i;
  
  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  Double_t delta = 0;
  //TF1* fit;
  double minchi2,mindelta;
  Double_t increment = 0;
  if(nDeltaSteps>0) increment = 2*TMath::Pi()/(nDeltaSteps);
  const int nd=nDeltaSteps+1;
  double *d = new double[nDeltaSteps+1];
  double *c2 = new double[nDeltaSteps+1];
  for(i=0;i<nd;i++)
  {
    delta = i*increment;
    d[i] = delta;
    
    for(ie=0;ie<ExtrapFHC.size();ie++)
    {
      ExtrapFHC[ie]->SetDeltaCP(delta);
      ExtrapFHC[ie]->OscillatePrediction();
    }
    for(ie=0;ie<ExtrapRHC.size();ie++)
    {
      ExtrapRHC[ie]->SetDeltaCP(delta);
      ExtrapRHC[ie]->OscillatePrediction();
    }
    Bkgd->Reset();
    Sig->Reset();
    NExp->Reset();
    
    //Turn extrapolations into a prediction:
    CombineFHCRHC();
    NExp->Add(Bkgd);
    NExp->Add(Sig);
    
    c2[i] = 1e10;
    if(FitMethod==0)
    {
      c2[i] = PoissonChi2(NExp);
    }
    else if(FitMethod==1)
    {
      c2[i] = ScaledChi2(Bkgd,Sig);
    }
    else if(FitMethod==2)
    {
      c2[i] = StandardChi2(NExp);
    }
    else if(FitMethod==3)
    {
      //Likelihood: "Standard" (N syst, N nuisance)
      //Calculate the likelihood (x2 for chi)
      c2[i] = StandardLikelihood();
    }
    else if(FitMethod==4)
    {
      //Likelihood: Bin by Bin Calculation of Systematics
      //Calculate the likelihood (x2 for chi)
      c2[i] = BinLikelihood();
    }
    else
    {
      cout<<"Error in GetMinLikelihood_Delta(): Unknown 'FitMethod'."<<endl;
    }
  }
  
  //can we actually fit this to a nice function?
  //   TGraph *g = new TGraph(ns,th13,c2);
  //   fit = new TF1("pol6", "pol6");
  //   g->Fit(fit, "Q");//fit to second order polynominal
  //   minchi2 = fit->GetMinimum(0,max);
  
  minchi2=1000;
  mindelta=0;
  for(i=0;i<nd;i++)
  {
    if(c2[i]<minchi2)
    {
      minchi2 = c2[i];
      mindelta = d[i];
    }
  }
  
  cout<<"min delta = "<<mindelta<<" , min chi2 = "<<minchi2<<endl;
  
  delete [] d;
  delete [] c2;
  
  return minchi2;
}
void NueFit2D_Joint::RunMultiBinPseudoExpts_MHDeltaFit(bool Print)
{
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
  {
    cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
    return;
  }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return;
  }
  for(unsigned int ie=0;ie<ExtrapFHC.size();ie++)
  {
    ExtrapFHC[ie]->GetPrediction();
  }
  for(unsigned int ie=0;ie<ExtrapRHC.size();ie++)
  {
    ExtrapRHC[ie]->GetPrediction();
  }
  if(ErrCalc==0)
  {
    cout<<"Need to set ErrorCalc object!  Quitting..."<<endl;
    return;
  }
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();
  
  nBinsFHC=0;
  nBinsRHC=0;  
  if (ExtrapFHC.size()>0) nBinsFHC = ExtrapFHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ExtrapRHC.size()>0) nBinsRHC = ExtrapRHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  nBins = nBinsFHC + nBinsRHC;
  
  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  TH2D *Error4Expts = new TH2D("Error4Expts","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  ReadGridFiles();
  
  if(nPts_Normal==0 || nPts_Inverted==0) return;
  
  gRandom->SetSeed(0);
  
  int i,u;
  unsigned int j,k;
  TH1D *nexp_bkgd_n = new TH1D("nexp_bkgd_n","",nBins,-0.5,nBins-0.5);
  TH1D *nexp_signal_n = new TH1D("nexp_signal_n","",nBins,-0.5,nBins-0.5);
  TH1D *nexp_n = new TH1D("nexp_n","",nBins,-0.5,nBins-0.5);
  
  TH1D *nexp_bkgd_i = new TH1D("nexp_bkgd_i","",nBins,-0.5,nBins-0.5);
  TH1D *nexp_signal_i = new TH1D("nexp_signal_i","",nBins,-0.5,nBins-0.5);
  TH1D *nexp_i = new TH1D("nexp_i","",nBins,-0.5,nBins-0.5);
  
  double chi2_norm,chi2_invt,chi2min_norm,chi2min_invt;
  double delchi2_norm, delchi2_invt,delchi2_bestmh;
  
  TH1D *chi2hist_deltafit = new TH1D("chi2hist_deltafit","",110000,-10,100);
  TH1D *chi2hist_mhfit = new TH1D("chi2hist_mhfit","",110000,-10,100);
  double ele;
  int noff;
  
  vector< vector<double> > nc_n,numucc_n,bnuecc_n,nutaucc_n,sig_n;
  vector< vector<double> > nc_i,numucc_i,bnuecc_i,nutaucc_i,sig_i;
  for(j=0;j<nBins;j++)
  {
    nc_n.push_back( vector<double>() );
    numucc_n.push_back( vector<double>() );
    bnuecc_n.push_back( vector<double>() );
    nutaucc_n.push_back( vector<double>() );
    sig_n.push_back( vector<double>() );
    
    nc_i.push_back( vector<double>() );
    numucc_i.push_back( vector<double>() );
    bnuecc_i.push_back( vector<double>() );
    nutaucc_i.push_back( vector<double>() );
    sig_i.push_back( vector<double>() );
    
    if (j<nBinsFHC){
      for(k=0;k<ExtrapFHC.size();k++)
      {
        nc_n[j].push_back(0);
        numucc_n[j].push_back(0);
        bnuecc_n[j].push_back(0);
        nutaucc_n[j].push_back(0);
        sig_n[j].push_back(0);
        
        nc_i[j].push_back(0);
        numucc_i[j].push_back(0);
        bnuecc_i[j].push_back(0);
        nutaucc_i[j].push_back(0);
        sig_i[j].push_back(0);
      }
    } else if (j>=nBinsFHC){
      for(k=0;k<ExtrapRHC.size();k++)
      {
        nc_n[j].push_back(0);
        numucc_n[j].push_back(0);
        bnuecc_n[j].push_back(0);
        nutaucc_n[j].push_back(0);
        sig_n[j].push_back(0);
        
        nc_i[j].push_back(0);
        numucc_i[j].push_back(0);
        bnuecc_i[j].push_back(0);
        nutaucc_i[j].push_back(0);
        sig_i[j].push_back(0);
      }
    }
  }
  
  Bkgd = (TH1D*)nexp_n->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)nexp_n->Clone("Sig");
  Sig->Reset();
  
  ofstream myfile;
  string file,ofile;
  
  //normal hierarchy
  
  TFile *f = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  
  if(Print)
  {
    ofile = gSystem->ExpandPathName(outFileName.c_str());
    file = ofile.substr(0,ofile.length()-5) + "_Normal.dat";
    myfile.open(gSystem->ExpandPathName(file.c_str()));
  }
  
  for(i=0;i<nPts_Normal;i++)
  {
    cout<<"point "<<(i+1)<<"/"<<nPts_Normal<<" (normal hierarchy)"<<endl;
    
    nexp_bkgd_n->Reset();
    nexp_signal_n->Reset();
    nexp_n->Reset();
    
    nexp_bkgd_i->Reset();
    nexp_signal_i->Reset();
    nexp_i->Reset();
    
    for(j=0;j<nBins;j++)
    {
      GridTree_Normal[j]->GetEntry(i);
      if (j<nBinsFHC){
        nexp_bkgd_n->SetBinContent(j+1,grid_background*GridScale_Normal_FHC);
        nexp_signal_n->SetBinContent(j+1,grid_signal*GridScale_Normal_FHC);
        for(k=0;k<ExtrapFHC.size();k++)
        {
          GridTree_2_Normal[j][k]->GetEntry(i);
          nc_n[j][k] = grid_nc*GridScale_Normal_FHC;
          numucc_n[j][k] = grid_numucc*GridScale_Normal_FHC;
          bnuecc_n[j][k] = grid_bnuecc*GridScale_Normal_FHC;
          nutaucc_n[j][k] = grid_nutaucc*GridScale_Normal_FHC;
          sig_n[j][k] = grid_nue*GridScale_Normal_FHC;
        }
      } else if (j>=nBinsFHC){
        nexp_bkgd_n->SetBinContent(j+1,grid_background*GridScale_Normal_RHC);
        nexp_signal_n->SetBinContent(j+1,grid_signal*GridScale_Normal_RHC);
        for(k=0;k<ExtrapRHC.size();k++)
        {
          GridTree_2_Normal[j][k]->GetEntry(i);
          nc_n[j][k] = grid_nc*GridScale_Normal_RHC;
          numucc_n[j][k] = grid_numucc*GridScale_Normal_RHC;
          bnuecc_n[j][k] = grid_bnuecc*GridScale_Normal_RHC;
          nutaucc_n[j][k] = grid_nutaucc*GridScale_Normal_RHC;
          sig_n[j][k] = grid_nue*GridScale_Normal_RHC;
        }
      }
    }
    for(j=0;j<nBins;j++)
    {
      GridTree_Inverted[j]->GetEntry(i);
      if (j<nBinsFHC){
        nexp_bkgd_i->SetBinContent(j+1,grid_background*GridScale_Inverted_FHC);
        nexp_signal_i->SetBinContent(j+1,grid_signal*GridScale_Inverted_FHC);
        for(k=0;k<ExtrapFHC.size();k++)
        {
          GridTree_2_Inverted[j][k]->GetEntry(i);
          nc_i[j][k] = grid_nc*GridScale_Inverted_FHC;
          numucc_i[j][k] = grid_numucc*GridScale_Inverted_FHC;
          bnuecc_i[j][k] = grid_bnuecc*GridScale_Inverted_FHC;
          nutaucc_i[j][k] = grid_nutaucc*GridScale_Inverted_FHC;
          sig_i[j][k] = grid_nue*GridScale_Inverted_FHC;
        }
      } else if (j>=nBinsFHC){
        nexp_bkgd_i->SetBinContent(j+1,grid_background*GridScale_Inverted_RHC);
        nexp_signal_i->SetBinContent(j+1,grid_signal*GridScale_Inverted_RHC);
        for(k=0;k<ExtrapRHC.size();k++)
        {
          GridTree_2_Inverted[j][k]->GetEntry(i);
          nc_i[j][k] = grid_nc*GridScale_Inverted_RHC;
          numucc_i[j][k] = grid_numucc*GridScale_Inverted_RHC;
          bnuecc_i[j][k] = grid_bnuecc*GridScale_Inverted_RHC;
          nutaucc_i[j][k] = grid_nutaucc*GridScale_Inverted_RHC;
          sig_i[j][k] = grid_nue*GridScale_Inverted_RHC;
        }
      }
    }
    nexp_n->Add(nexp_bkgd_n,nexp_signal_n,1,1);
    nexp_i->Add(nexp_bkgd_i,nexp_signal_i,1,1);
    
    //want to generate pseudos for normal hierarchy
    ErrCalc->SetGridPred(nBins,nc_n,numucc_n,bnuecc_n,nutaucc_n,sig_n);
    
    Error4Expts->Reset();
    ErrCalc->SetUseGrid(true);
    ErrCalc->CalculateSystErrorMatrix();
    Error4Expts->Add(ErrCalc->CovMatrix);
    //ErrCalc->CalculateHOOError();
    //Error4Expts->Add(ErrCalc->CovMatrix_Decomp);
    
    if(IncludeOscParErrs)
    {
      noff=0;
      for(j=0;j<nBins;j++)
      {
        GridTree_Normal[j]->GetEntry(i);
        ele=Error4Expts->GetBinContent(j+1,j+1);
        ele+=(grid_bin_oscparerr[j]*grid_bin_oscparerr[j]*nexp_n->GetBinContent(j+1)*nexp_n->GetBinContent(j+1));
        Error4Expts->SetBinContent(j+1,j+1,ele);
        
        for(k=0;k<nBins;k++)
        {
          if(k>j)
          {
            ele=Error4Expts->GetBinContent(j+1,k+1);
            ele+=(grid_bin_oscparerr[j]*grid_bin_oscparerr[k]*nexp_n->GetBinContent(j+1)*nexp_n->GetBinContent(k+1));
            Error4Expts->SetBinContent(j+1,k+1,ele);
            
            ele=Error4Expts->GetBinContent(k+1,j+1);
            ele+=(grid_bin_oscparerr[j]*grid_bin_oscparerr[k]*nexp_n->GetBinContent(j+1)*nexp_n->GetBinContent(k+1));
            Error4Expts->SetBinContent(k+1,j+1,ele);
            
            noff++;
          }
        }
      }
    }
    
    chi2hist_deltafit->Reset();
    chi2hist_deltafit->SetName(Form("Chi2Hist_DeltaFit_Normal_%i",i));
    
    chi2hist_mhfit->Reset();
    chi2hist_mhfit->SetName(Form("Chi2Hist_MHFit_Normal_%i",i));
    
    for(u=0;u<NumExpts;u++)
    {
      cout<<"expt "<<(u+1)<<"/"<<NumExpts<<endl;
      
      GenerateOneCorrelatedExp(nexp_n,Error4Expts);
      if(Print)
      {
        myfile << grid_sinsq2th13 << " " << grid_delta << " ";
        for(j=0;j<nBins;j++)
        {
          myfile << NObs->GetBinContent(j+1) << " ";
        }
      }
      
      chi2min_norm=GetMinLikelihood_Delta(true);
      chi2min_invt=GetMinLikelihood_Delta(false);
      
      ErrCalc->SetGridPred(nBins,nc_n,numucc_n,bnuecc_n,nutaucc_n,sig_n);
      ErrCalc->SetUseGrid(true);
      chi2_norm = 1e10;
      if(FitMethod==0)
      {
        chi2_norm = PoissonChi2(nexp_n);
      }
      else if(FitMethod==1)
      {
        chi2_norm = ScaledChi2(nexp_bkgd_n,nexp_signal_n);
      }
      else if(FitMethod==2)
      {
        chi2_norm = StandardChi2(nexp_n);
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_n);
        Sig->Reset();
        Sig->Add(nexp_signal_n);
        chi2_norm = StandardLikelihood();
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_n);
        Sig->Reset();
        Sig->Add(nexp_signal_n);
        chi2_norm = BinLikelihood();
      }
      else
      {
        cout<<"Error in RunMultiBinPseudoExpts_MHDeltaFit(): Unknown 'FitMethod'."<<endl;
      }
      
      ErrCalc->SetGridPred(nBins,nc_i,numucc_i,bnuecc_i,nutaucc_i,sig_i);
      ErrCalc->SetUseGrid(true);
      chi2_invt = 1e10;
      if(FitMethod==0)
      {
        chi2_invt = PoissonChi2(nexp_i);
      }
      else if(FitMethod==1)
      {
        chi2_invt = ScaledChi2(nexp_bkgd_i,nexp_signal_i);
      }
      else if(FitMethod==2)
      {
        chi2_invt = StandardChi2(nexp_i);
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_i);
        Sig->Reset();
        Sig->Add(nexp_signal_i);
        chi2_invt = StandardLikelihood();
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_i);
        Sig->Reset();
        Sig->Add(nexp_signal_i);
        chi2_invt = BinLikelihood();
      }
      else
      {
        cout<<"Error in RunMultiBinPseudoExpts_MHDeltaFit(): Unknown 'FitMethod'."<<endl;
      }
      
      if(Print)
      {
        myfile << chi2_norm << " " << chi2_invt << " " << chi2min_norm << " " << chi2min_invt << endl;
      }
      
      delchi2_norm = chi2_norm - chi2min_norm;
      delchi2_invt = chi2_invt - chi2min_invt;
      delchi2_bestmh = delchi2_norm;
      if(delchi2_invt<delchi2_bestmh) delchi2_bestmh = delchi2_invt;
      
      chi2hist_deltafit->Fill(delchi2_norm);
      chi2hist_mhfit->Fill(delchi2_norm - delchi2_bestmh);
    }
    f->cd();
    chi2hist_deltafit->Write();
    chi2hist_mhfit->Write();
    f->Close();
    
    f = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"UPDATE");
  }
  
  if(Print) myfile.close();
  
  for(j=0;j<nBins;j++)
  {
    GridTree_Normal[j]->Write();
    
    if (j<nBinsFHC){
      for(k=0;k<ExtrapFHC.size();k++)
      {
        GridTree_2_Normal[j][k]->Write();
      }
    } else if (j>=nBinsFHC){
      for(k=0;k<ExtrapRHC.size();k++)
      {
        GridTree_2_Normal[j][k]->Write();
      }
    }
  }
  
  f->Close();
  
  //inverted hierarchy
  
  f = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"UPDATE");
  
  if(Print)
  {
    ofile = gSystem->ExpandPathName(outFileName.c_str());
    file = ofile.substr(0,ofile.length()-5) + "_Inverted.dat";
    myfile.open(gSystem->ExpandPathName(file.c_str()));
  }
  
  for(i=0;i<nPts_Inverted;i++)
  {
    cout<<"point "<<(i+1)<<"/"<<nPts_Inverted<<" (inverted hierarchy)"<<endl;
    
    nexp_bkgd_n->Reset();
    nexp_signal_n->Reset();
    nexp_n->Reset();
    
    nexp_bkgd_i->Reset();
    nexp_signal_i->Reset();
    nexp_i->Reset();
    
    for(j=0;j<nBins;j++)
    {
      GridTree_Normal[j]->GetEntry(i);
      if (j<nBinsFHC){
        nexp_bkgd_n->SetBinContent(j+1,grid_background*GridScale_Normal_FHC);
        nexp_signal_n->SetBinContent(j+1,grid_signal*GridScale_Normal_FHC);
        for(k=0;k<ExtrapFHC.size();k++)
        {
          GridTree_2_Normal[j][k]->GetEntry(i);
          nc_n[j][k] = grid_nc*GridScale_Normal_FHC;
          numucc_n[j][k] = grid_numucc*GridScale_Normal_FHC;
          bnuecc_n[j][k] = grid_bnuecc*GridScale_Normal_FHC;
          nutaucc_n[j][k] = grid_nutaucc*GridScale_Normal_FHC;
          sig_n[j][k] = grid_nue*GridScale_Normal_FHC;
        }
      } else if (j>=nBinsFHC){
        nexp_bkgd_n->SetBinContent(j+1,grid_background*GridScale_Normal_RHC);
        nexp_signal_n->SetBinContent(j+1,grid_signal*GridScale_Normal_RHC);
        for(k=0;k<ExtrapRHC.size();k++)
        {
          GridTree_2_Normal[j][k]->GetEntry(i);
          nc_n[j][k] = grid_nc*GridScale_Normal_RHC;
          numucc_n[j][k] = grid_numucc*GridScale_Normal_RHC;
          bnuecc_n[j][k] = grid_bnuecc*GridScale_Normal_RHC;
          nutaucc_n[j][k] = grid_nutaucc*GridScale_Normal_RHC;
          sig_n[j][k] = grid_nue*GridScale_Normal_RHC;
        }
      }
    }
    for(j=0;j<nBins;j++)
    {
      GridTree_Inverted[j]->GetEntry(i);
      if (j<nBinsFHC){
        nexp_bkgd_i->SetBinContent(j+1,grid_background*GridScale_Inverted_FHC);
        nexp_signal_i->SetBinContent(j+1,grid_signal*GridScale_Inverted_FHC);
        for(k=0;k<ExtrapFHC.size();k++)
        {
          GridTree_2_Inverted[j][k]->GetEntry(i);
          nc_i[j][k] = grid_nc*GridScale_Inverted_FHC;
          numucc_i[j][k] = grid_numucc*GridScale_Inverted_FHC;
          bnuecc_i[j][k] = grid_bnuecc*GridScale_Inverted_FHC;
          nutaucc_i[j][k] = grid_nutaucc*GridScale_Inverted_FHC;
          sig_i[j][k] = grid_nue*GridScale_Inverted_FHC;
        }
      } else if (j>=nBinsFHC){
        nexp_bkgd_i->SetBinContent(j+1,grid_background*GridScale_Inverted_RHC);
        nexp_signal_i->SetBinContent(j+1,grid_signal*GridScale_Inverted_RHC);
        for(k=0;k<ExtrapRHC.size();k++)
        {
          GridTree_2_Inverted[j][k]->GetEntry(i);
          nc_i[j][k] = grid_nc*GridScale_Inverted_RHC;
          numucc_i[j][k] = grid_numucc*GridScale_Inverted_RHC;
          bnuecc_i[j][k] = grid_bnuecc*GridScale_Inverted_RHC;
          nutaucc_i[j][k] = grid_nutaucc*GridScale_Inverted_RHC;
          sig_i[j][k] = grid_nue*GridScale_Inverted_RHC;
        }
      }
    }
    nexp_n->Add(nexp_bkgd_n,nexp_signal_n,1,1);
    nexp_i->Add(nexp_bkgd_i,nexp_signal_i,1,1);
    
    //want to generate pseudos for inverted hierarchy
    ErrCalc->SetGridPred(nBins,nc_i,numucc_i,bnuecc_i,nutaucc_i,sig_i);
    Error4Expts->Reset();
    ErrCalc->SetUseGrid(true);
    ErrCalc->CalculateSystErrorMatrix();
    Error4Expts->Add(ErrCalc->CovMatrix);
    //ErrCalc->CalculateHOOError();
    //Error4Expts->Add(ErrCalc->CovMatrix_Decomp);
    if(IncludeOscParErrs)
    {
      noff=0;
      for(j=0;j<nBins;j++)
      {
        ele=Error4Expts->GetBinContent(j+1,j+1);
        ele+=(grid_bin_oscparerr[j]*grid_bin_oscparerr[j]*nexp_i->GetBinContent(j+1)*nexp_i->GetBinContent(j+1));
        Error4Expts->SetBinContent(j+1,j+1,ele);
        
        for(k=0;k<nBins;k++)
        {
          if(k>j)
          {
            ele=Error4Expts->GetBinContent(j+1,k+1);
            ele+=(grid_bin_oscparerr[j]*grid_bin_oscparerr[k]*nexp_i->GetBinContent(j+1)*nexp_i->GetBinContent(k+1));
            Error4Expts->SetBinContent(j+1,k+1,ele);
            
            ele=Error4Expts->GetBinContent(k+1,j+1);
            ele+=(grid_bin_oscparerr[j]*grid_bin_oscparerr[k]*nexp_i->GetBinContent(j+1)*nexp_i->GetBinContent(k+1));
            Error4Expts->SetBinContent(k+1,j+1,ele);
            
            noff++;
          }
        }
      }
    }
    
    chi2hist_deltafit->Reset();
    chi2hist_deltafit->SetName(Form("Chi2Hist_DeltaFit_Inverted_%i",i));
    
    chi2hist_mhfit->Reset();
    chi2hist_mhfit->SetName(Form("Chi2Hist_MHFit_Inverted_%i",i));
    
    for(u=0;u<NumExpts;u++)
    {
      cout<<"expt "<<(u+1)<<"/"<<NumExpts<<endl;
      
      GenerateOneCorrelatedExp(nexp_i,Error4Expts);
      if(Print)
      {
        myfile << grid_sinsq2th13 << " " << grid_delta << " ";
        for(j=0;j<nBins;j++)
        {
          myfile << NObs->GetBinContent(j+1) << " ";
        }
      }
      
      chi2min_norm=GetMinLikelihood_Delta(true);
      chi2min_invt=GetMinLikelihood_Delta(false);
      
      ErrCalc->SetGridPred(nBins,nc_n,numucc_n,bnuecc_n,nutaucc_n,sig_n);
      ErrCalc->SetUseGrid(true);
      chi2_norm = 1e10;
      if(FitMethod==0)
      {
        chi2_norm = PoissonChi2(nexp_n);
      }
      else if(FitMethod==1)
      {
        chi2_norm = ScaledChi2(nexp_bkgd_n,nexp_signal_n);
      }
      else if(FitMethod==2)
      {
        chi2_norm = StandardChi2(nexp_n);
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_n);
        Sig->Reset();
        Sig->Add(nexp_signal_n);
        chi2_norm = StandardLikelihood();
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_n);
        Sig->Reset();
        Sig->Add(nexp_signal_n);
        chi2_norm = BinLikelihood();
      }
      else
      {
        cout<<"Error in RunMultiBinPseudoExpts_MHDeltaFit(): Unknown 'FitMethod'."<<endl;
      }
      
      ErrCalc->SetGridPred(nBins,nc_i,numucc_i,bnuecc_i,nutaucc_i,sig_i);
      ErrCalc->SetUseGrid(true);
      chi2_invt = 1e10;
      if(FitMethod==0)
      {
        chi2_invt = PoissonChi2(nexp_i);
      }
      else if(FitMethod==1)
      {
        chi2_invt = ScaledChi2(nexp_bkgd_i,nexp_signal_i);
      }
      else if(FitMethod==2)
      {
        chi2_invt = StandardChi2(nexp_i);
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_i);
        Sig->Reset();
        Sig->Add(nexp_signal_i);
        chi2_invt = StandardLikelihood();
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        Bkgd->Reset();
        Bkgd->Add(nexp_bkgd_i);
        Sig->Reset();
        Sig->Add(nexp_signal_i);
        chi2_invt = BinLikelihood();
      }
      else
      {
        cout<<"Error in RunMultiBinPseudoExpts_MHDeltaFit(): Unknown 'FitMethod'."<<endl;
      }
      
      if(Print)
      {
        myfile << chi2_norm << " " << chi2_invt << " " << chi2min_norm << " " << chi2min_invt << endl;
      }
      
      delchi2_norm = chi2_norm - chi2min_norm;
      delchi2_invt = chi2_invt - chi2min_invt;
      delchi2_bestmh = delchi2_invt;
      if(delchi2_norm<delchi2_bestmh) delchi2_bestmh = delchi2_norm;
      
      chi2hist_deltafit->Fill(delchi2_invt);
      chi2hist_mhfit->Fill(delchi2_invt - delchi2_bestmh);
    }
    f->cd();
    chi2hist_deltafit->Write();
    chi2hist_mhfit->Write();
    f->Close();
    f = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"UPDATE");
  }
  
  if(Print) myfile.close();
  
  for(j=0;j<nBins;j++)
  {
    GridTree_Inverted[j]->Write();
    if (j<nBinsFHC){
      for(k=0;k<ExtrapFHC.size();k++)
      {
        GridTree_2_Inverted[j][k]->Write();
      }
    } else if (j>=nBinsFHC){
      for(k=0;k<ExtrapRHC.size();k++)
      {
        GridTree_2_Inverted[j][k]->Write();
      }
    }
  }
  
  f->Close();
  
  return;
}
void NueFit2D_Joint::RunMultiBinFC_MHDeltaFit()
{
  if(NObs==0)
  {
    cout<<"NObs not set.  Quitting..."<<endl;
    return;
  }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return;
  }
  for(unsigned int ie=0;ie<ExtrapFHC.size();ie++)
  {
    ExtrapFHC[ie]->GetPrediction();
  }
  for(unsigned int ie=0;ie<ExtrapRHC.size();ie++)
  {
    ExtrapRHC[ie]->GetPrediction();
  }
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
  {
    cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
    return;
  }
  if(ErrCalc==0)
  {
    cout<<"No ErrorCalc object set!  Quitting..."<<endl;
    return;
  }
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();
  
  nBinsFHC=0;
  nBinsRHC=0;  
  if (ExtrapFHC.size()>0) nBinsFHC = ExtrapFHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ExtrapRHC.size()>0) nBinsRHC = ExtrapRHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  nBins = nBinsFHC + nBinsRHC;
  
  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  ReadGridFiles();
  SetupChi2Hists();
  
  if(nPts_Normal==0 || nPts_Inverted==0) return;
  
  TH1D *nexp_bkgd = new TH1D("nexp_bkgd","",nBins,-0.5,nBins-0.5);
  TH1D *nexp_signal = new TH1D("nexp_signal","",nBins,-0.5,nBins-0.5);
  TH1D *nexp = new TH1D("nexp","",nBins,-0.5,nBins-0.5);
  
  int i;
  unsigned int j,k;
  
  if(gSystem->AccessPathName(gSystem->ExpandPathName(PseudoExpFile.c_str())))
  {
    cout<<"Pseudo-experiment file doesn't exist."<<endl;
    return;
  }
  
  TFile *f = new TFile(gSystem->ExpandPathName(PseudoExpFile.c_str()),"READ");
  
  double *chi2data_NH = new double[nPts_Normal];
  double *chi2data_IH = new double[nPts_Inverted];
  double *chi2min_NH = new double[nPts_Normal];
  double *chi2min_IH = new double[nPts_Inverted];
  
  int chi2databin;
  TH1D *chi2hist = (TH1D*)f->Get("Chi2Hist_MHFit_Normal_0");
  double chi2binwidth = chi2hist->GetBinWidth(1);
  double chi2start = chi2hist->GetXaxis()->GetBinLowEdge(1);
  double frac;
  delete chi2hist;
  
  double *deltapts = new double[nPts_Normal];
  
  vector< vector<double> > nc,numucc,bnuecc,nutaucc,sig;
  for(j=0;j<nBins;j++)
  {
    nc.push_back( vector<double>() );
    numucc.push_back( vector<double>() );
    bnuecc.push_back( vector<double>() );
    nutaucc.push_back( vector<double>() );
    sig.push_back( vector<double>() );
    
    if (j<nBinsFHC){
      for(k=0;k<ExtrapFHC.size();k++)
      {
        nc[j].push_back(0);
        numucc[j].push_back(0);
        bnuecc[j].push_back(0);
        nutaucc[j].push_back(0);
        sig[j].push_back(0);
      }
    } else if (j>=nBinsFHC){
      for(k=0;k<ExtrapRHC.size();k++)
      {
        nc[j].push_back(0);
        numucc[j].push_back(0);
        bnuecc[j].push_back(0);
        nutaucc[j].push_back(0);
        sig[j].push_back(0);
      }
    }
  }
  
  Bkgd = (TH1D*)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)NObs->Clone("Sig");
  Sig->Reset();
  
  //normal hierarchy
  
  for(i=0;i<nPts_Normal;i++)
  {
    if(i%100==0) cout<<100.*i/nPts_Normal<<"% complete for normal hierarchy"<<endl;
    
    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();
    
    for(j=0;j<nBins;j++)
    {
      GridTree_Normal[j]->GetEntry(i);
      if (j<nBinsFHC){
        nexp_bkgd->SetBinContent(j+1,grid_background*GridScale_Normal_FHC);
        nexp_signal->SetBinContent(j+1,grid_signal*GridScale_Normal_FHC);
        for(k=0;k<ExtrapFHC.size();k++)
        {
          GridTree_2_Normal[j][k]->GetEntry(i);
          nc[j][k] = grid_nc*GridScale_Normal_FHC;
          numucc[j][k] = grid_numucc*GridScale_Normal_FHC;
          bnuecc[j][k] = grid_bnuecc*GridScale_Normal_FHC;
          nutaucc[j][k] = grid_nutaucc*GridScale_Normal_FHC;
          sig[j][k] = grid_nue*GridScale_Normal_FHC;
        }
      } else if (j>=nBinsFHC){
        nexp_bkgd->SetBinContent(j+1,grid_background*GridScale_Normal_RHC);
        nexp_signal->SetBinContent(j+1,grid_signal*GridScale_Normal_RHC);
        for(k=0;k<ExtrapRHC.size();k++)
        {
          GridTree_2_Normal[j][k]->GetEntry(i);
          nc[j][k] = grid_nc*GridScale_Normal_RHC;
          numucc[j][k] = grid_numucc*GridScale_Normal_RHC;
          bnuecc[j][k] = grid_bnuecc*GridScale_Normal_RHC;
          nutaucc[j][k] = grid_nutaucc*GridScale_Normal_RHC;
          sig[j][k] = grid_nue*GridScale_Normal_RHC;
        }
      }
    }
    nexp->Add(nexp_bkgd,nexp_signal,1,1);
    ErrCalc->SetGridPred(nBins,nc,numucc,bnuecc,nutaucc,sig);
    
    chi2min_NH[i]=GetMinLikelihood_Delta(true);
    
    ErrCalc->SetUseGrid(true);//use grid predictions set up above
    chi2data_NH[i] = 1e10;
    if(FitMethod==0)
    {
      chi2data_NH[i] = PoissonChi2(nexp);
    }
    else if(FitMethod==1)
    {
      chi2data_NH[i] = ScaledChi2(Bkgd,Sig);
    }
    else if(FitMethod==2)
    {
      chi2data_NH[i] = StandardChi2(nexp);
    }
    else if(FitMethod==3)
    {
      //Likelihood: "Standard" (N syst, N nuisance)
      //Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data_NH[i] = StandardLikelihood();
    }
    else if(FitMethod==4)
    {
      //Likelihood: Bin by Bin Calculation of Systematics
      //Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data_NH[i] = BinLikelihood();
    }
    else
    {
      cout<<"Error in RunMultiBinFC_MHDeltaFit(): Unknown 'FitMethod'."<<endl;
    }
    
    deltapts[i] = grid_delta/TMath::Pi();
  }
  
  //inverted hierarchy
  
  for(i=0;i<nPts_Inverted;i++)
  {
    if(i%100==0) cout<<100.*i/nPts_Inverted<<"% complete for inverted hierarchy"<<endl;
    
    nexp_bkgd->Reset();
    nexp_signal->Reset();
    nexp->Reset();
    
    for(j=0;j<nBins;j++)
    {
      GridTree_Inverted[j]->GetEntry(i);
      if (j<nBinsFHC){
        nexp_bkgd->SetBinContent(j+1,grid_background*GridScale_Inverted_FHC);
        nexp_signal->SetBinContent(j+1,grid_signal*GridScale_Inverted_FHC);
        for(k=0;k<ExtrapFHC.size();k++)
        {
          GridTree_2_Inverted[j][k]->GetEntry(i);
          nc[j][k] = grid_nc*GridScale_Inverted_FHC;
          numucc[j][k] = grid_numucc*GridScale_Inverted_FHC;
          bnuecc[j][k] = grid_bnuecc*GridScale_Inverted_FHC;
          nutaucc[j][k] = grid_nutaucc*GridScale_Inverted_FHC;
          sig[j][k] = grid_nue*GridScale_Inverted_FHC;
        }
      } else if (j>=nBinsFHC){
        nexp_bkgd->SetBinContent(j+1,grid_background*GridScale_Inverted_RHC);
        nexp_signal->SetBinContent(j+1,grid_signal*GridScale_Inverted_RHC);
        for(k=0;k<ExtrapRHC.size();k++)
        {
          GridTree_2_Inverted[j][k]->GetEntry(i);
          nc[j][k] = grid_nc*GridScale_Inverted_RHC;
          numucc[j][k] = grid_numucc*GridScale_Inverted_RHC;
          bnuecc[j][k] = grid_bnuecc*GridScale_Inverted_RHC;
          nutaucc[j][k] = grid_nutaucc*GridScale_Inverted_RHC;
          sig[j][k] = grid_nue*GridScale_Inverted_RHC;
        }
      }
    }
    nexp->Add(nexp_bkgd,nexp_signal,1,1);
    ErrCalc->SetGridPred(nBins,nc,numucc,bnuecc,nutaucc,sig);
    
    chi2min_IH[i]=GetMinLikelihood_Delta(false);
    
    ErrCalc->SetUseGrid(true);//use grid predictions set up above
    chi2data_IH[i] = 1e10;
    if(FitMethod==0)
    {
      chi2data_IH[i] = PoissonChi2(nexp);
    }
    else if(FitMethod==1)
    {
      chi2data_IH[i] = ScaledChi2(Bkgd,Sig);
    }
    else if(FitMethod==2)
    {
      chi2data_IH[i] = StandardChi2(nexp);
    }
    else if(FitMethod==3)
    {
      //Likelihood: "Standard" (N syst, N nuisance)
      //Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data_IH[i] = StandardLikelihood();
    }
    else if(FitMethod==4)
    {
      //Likelihood: Bin by Bin Calculation of Systematics
      //Calculate the likelihood (x2 for chi)
      Bkgd->Reset();
      Bkgd->Add(nexp_bkgd);
      Sig->Reset();
      Sig->Add(nexp_signal);
      chi2data_IH[i] = BinLikelihood();
    }
    else
    {
      cout<<"Error in RunMultiBinFC_MHDeltaFit(): Unknown 'FitMethod'."<<endl;
    }
  }
  
  if(nPts_Normal!=nPts_Inverted) cout<<"different number of points"<<endl;
  
  double *deltafit_nh = new double[nPts_Normal];
  double *deltafit_ih = new double[nPts_Inverted];
  double *mhfit_nh = new double[nPts_Normal];
  double *mhfit_ih = new double[nPts_Inverted];
  double delchi2,delchi2_norm,delchi2_invt,delchi2_bestmh;
  
  //delta fit
  for(i=0;i<nPts_Normal;i++)
  {
    //normal
    delchi2 = chi2data_NH[i] - chi2min_NH[i];
    chi2databin = int((delchi2-chi2start)/chi2binwidth)+1;
    
    TH1D *chi2hist = (TH1D*)f->Get(Form("Chi2Hist_DeltaFit_Normal_%i",i));
    
    if(chi2hist->Integral()<1)
    {
      cout<<"Warning, chi2hist is empty."<<endl;
      frac=0;
    }
    else
    {
      frac = chi2hist->Integral(1,chi2databin-1)/chi2hist->Integral();
    }
    if(chi2databin==1) frac=0;
    
    deltafit_nh[i] = frac;
    
    //inverted
    delchi2 = chi2data_IH[i] - chi2min_IH[i];
    chi2databin = int((delchi2-chi2start)/chi2binwidth)+1;
    
    chi2hist = (TH1D*)f->Get(Form("Chi2Hist_DeltaFit_Inverted_%i",i));
    
    if(chi2hist->Integral()<1)
    {
      cout<<"Warning, chi2hist is empty."<<endl;
      frac=0;
    }
    else
    {
      frac = chi2hist->Integral(1,chi2databin-1)/chi2hist->Integral();
    }
    if(chi2databin==1) frac=0;
    
    deltafit_ih[i] = frac;
    
    delete chi2hist;
  }
  
  //mh fit
  for(i=0;i<nPts_Normal;i++)
  {
    delchi2_norm = chi2data_NH[i] - chi2min_NH[i];
    delchi2_invt = chi2data_IH[i] - chi2min_IH[i];
    delchi2_bestmh = delchi2_invt;
    if(delchi2_norm<delchi2_bestmh) delchi2_bestmh = delchi2_norm;
    cout<<"NH: chi2min = "<<chi2min_NH[i]<<", IH: chi2min = "<<chi2min_IH[i]<<endl;
    
    delchi2 = delchi2_norm - delchi2_bestmh;
    chi2databin = int((delchi2-chi2start)/chi2binwidth)+1;
    
    TH1D *chi2hist = (TH1D*)f->Get(Form("Chi2Hist_MHFit_Normal_%i",i));
    
    if(chi2hist->Integral()<1)
    {
      cout<<"Warning, chi2hist is empty."<<endl;
      frac=0;
    }
    else
    {
      frac = chi2hist->Integral(1,chi2databin-1)/chi2hist->Integral();
    }
    if(chi2databin==1) frac=0;
    
    mhfit_nh[i] = frac;
    
    delchi2 = delchi2_invt - delchi2_bestmh;
    chi2databin = int((delchi2-chi2start)/chi2binwidth)+1;
    
    chi2hist = (TH1D*)f->Get(Form("Chi2Hist_MHFit_Inverted_%i",i));
    
    if(chi2hist->Integral()<1)
    {
      cout<<"Warning, chi2hist is empty."<<endl;
      frac=0;
    }
    else
    {
      frac = chi2hist->Integral(1,chi2databin-1)/chi2hist->Integral();
    }
    if(chi2databin==1) frac=0;
    
    mhfit_ih[i] = frac;
    
    delete chi2hist;
  }

  f->Close();
  
  TGraph *g_delta_NH = new TGraph(nPts_Normal,deltapts,deltafit_nh);
  g_delta_NH->SetName("delta_NH");
  
  TGraph *g_delta_IH = new TGraph(nPts_Normal,deltapts,deltafit_ih);
  g_delta_IH->SetName("delta_IH");
  
  TGraph *g_MH_NH = new TGraph(nPts_Normal,deltapts,mhfit_nh);
  g_MH_NH->SetName("mh_NH");
  
  TGraph *g_MH_IH = new TGraph(nPts_Normal,deltapts,mhfit_ih);
  g_MH_IH->SetName("mh_IH");
  
  TFile *fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  g_delta_NH->Write();
  g_delta_IH->Write();
  g_MH_NH->Write();
  g_MH_IH->Write();
  fout->Close();
  
  delete [] chi2data_NH;
  delete [] chi2data_IH;
  delete [] chi2min_NH;
  delete [] chi2min_IH;
  delete [] deltafit_nh;
  delete [] deltafit_ih;
  delete [] mhfit_nh;
  delete [] mhfit_ih;
  
  return;
}







////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//ADAM'S JUNK BEGINS HERE!
void NueFit2D_Joint::DoDeltaFit()
{
  //probably don't want to set central values in the function call
  //make separate functions that also set the uncertainties, ie
  //NueFit2D::SetTheta13(theta13 central value, +1 sigma error, -1sigam error)
    
  if(NObs==0)
    {
      cout<<"NObs not set.  Quitting..."<<endl;
      return;
    }
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
    {
      cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
      return;
    }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
    {
      cout<<"No Extrapolate2D input.  Quitting..."<<endl;
      return;
    }
  
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();
  
  Bkgd = (TH1D*)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)NObs->Clone("Sig");
  Sig->Reset();
  TH1D *NExp = (TH1D*)NObs->Clone("NExp");
  NExp->Reset();
  
  unsigned int ie;
  for(ie=0;ie<ExtrapFHC.size();ie++)
    {
      ExtrapFHC[ie]->GetPrediction();
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,Theta12);
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,Theta23);
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh13,Theta13);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq32);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq21);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDelta,delta);
      ExtrapFHC[ie]->OscillatePrediction();
    }
  for(ie=0;ie<ExtrapRHC.size();ie++)
    {
      ExtrapRHC[ie]->GetPrediction();
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,Theta12);
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,Theta23);
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh13,Theta13);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq32);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq21);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDelta,delta);
      ExtrapRHC[ie]->OscillatePrediction();
    }
  
  Int_t i;  
  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  //Double_t ss2th13 = 0;
  Double_t delchi2 = 0;
  Double_t best = 0;
  Double_t ttbest = 0;
  Double_t otbest = 0;
  Double_t dmsbest = 0;
  //
  double chi2[3],val[3]; 
  TGraph *g3;
  TF1* fit;
  double minchi2;
  double global_minchi2;
  
  //find global minimum  
  //find minimum in delta
  best=-1.;
  minchi2=100;
  for(i=0;i<3;i++)
    {
      chi2[i]=-1.;
      val[i]=-1.;
    }
  
  double inc = 100.0;
  double Deltaincrement = (1/inc)*TMath::Pi(); 
  double db;
  for(i=0; i<2*inc+1 ; i++)//i scans over all delta values so we can find the best fit delta
    {
      chi2[0] = chi2[1];
      chi2[1] = chi2[2];
      val[0] = val[1];
      val[1] = val[2];
      
      db = i*Deltaincrement + 0.5;
      cout << "It: " << i << " with delta: " << db << endl;
      for(ie=0;ie<ExtrapFHC.size();ie++)
	{
	  ExtrapFHC[ie]->SetOscPar(OscPar::kDelta,db);
	  ExtrapFHC[ie]->OscillatePrediction();
	}
      
      for(ie=0;ie<ExtrapRHC.size();ie++)
	{
	  ExtrapRHC[ie]->SetOscPar(OscPar::kDelta,db);
	  ExtrapRHC[ie]->OscillatePrediction();
	}
      
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      
      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      chi2[2] = 1e10;
      if(FitMethod==0)
	{
	  chi2[2] = PoissonChi2(NExp);
	}
      else if(FitMethod==1)
	{
	  chi2[2] = ScaledChi2(Bkgd,Sig);
	}
      else if(FitMethod==2)
	{
	  chi2[2] = StandardChi2(NExp);
	}
      else if(FitMethod==3)
	{
	  //Likelihood: "Standard" (N syst, N nuisance)
	  //Calculate the likelihood (x2 for chi)
	  chi2[2] = StandardLikelihood();
	}
      else if(FitMethod==4)
	{
	  //Likelihood: Bin by Bin Calculation of Systematics
	  //Calculate the likelihood (x2 for chi)
	  chi2[2] = BinLikelihood();
	  //cout << "ss2th13 = " << ss2th13 << " Chi2=" << chi2[2] << endl;
	}
      else
	{
	  cout<<"Error in GetLikelihood(): Unknown 'FitMethod'."<<endl;
	}
      
      val[2] = db;
      
      if(i<2) continue;
      
      if(i<3 && chi2[2]>chi2[1] && chi2[1]>chi2[0])//first three points are increasing, first point is minimum.
	{
	  best = val[0];
	  minchi2 = chi2[0];
	  cout<<"minimum at 1st point: "<<minchi2<<" at "<<best<<endl;
	  break;
	}
      
      if(chi2[2]>chi2[1] && chi2[0]>chi2[1])//found minimum
	{
	  g3 = new TGraph(3, val, chi2);
	  fit = new TF1("pol2", "pol2");
	  g3->Fit(fit, "Q");//fit to second order polynominal
	  if(fit->GetParameter(2) > 0)//if the x^2 term is nonzero
	    {
	      best = -fit->GetParameter(1)/(2*fit->GetParameter(2));//the location of the minimum is -p1/(2*p2)
	      minchi2 = fit->GetParameter(0) + fit->GetParameter(1)*best + fit->GetParameter(2)*best*best;
	      cout<<"minimum with fit: "<<minchi2<<" at "<<best<<endl;
	    }
	  else//if the x^2 term is zero, then just use the minimum you got by scanning
	    {
	      best = val[1];
	      minchi2 = chi2[1];
	      cout<<"minimum with scan: "<<minchi2<<" at "<<best<<endl;
	    }
	  break;
	}
    }
  
  double f;
  minchi2=1e10;
  inc = 2.0;//20.0;
  int d,a,r,m;
  Deltaincrement = 0.01; //*TMath::Pi();
  double Aincrement =2*Theta23Unc/inc;
  double Rincrement =2*Theta13Unc/inc;
  double Mincrement = 0.01e-3;
  double das, aas, ras, mas;
  double dbest; 
  double penpar;
  
  double ass, rss;
  //NEW
  das = 0.0;
  ass = 0.3;
  rss = 0.0167;
  mas = -2.70e-3;
  
  //Delta Increment
  for(d=0;d<=10;d++){
    das = best - (5)*Deltaincrement + d*Deltaincrement;
    
    //Delta increment for Andy's combined surface
    //das = d*(2*TMath::Pi())/32.0;
    
    //Theta 23 increment
    for(a=0;a<=inc;a++){
      aas = Theta23 - (inc/2)*Aincrement + a*Aincrement;
      
      //Theta23 increment used for Andy's combined surface
      //for(a=0;a<8;a++){
      //ass = 0.3 + a*0.06;
      //aas = TMath::ASin(TMath::Sqrt(ass)); 
      
      //Theta13 increment
      for(r=0;r<=inc;r++){
	ras = Theta13 - (inc/2)*Rincrement+ (r)*Rincrement;
	
	//Theta13 increment used for Andy's combined surface  
	//for(r=0;r<7;r++){
	//rss = 0.0167 + r*(0.0317-0.0167)/6;
	// ras = TMath::ASin(TMath::Sqrt(rss));
	
	//Mass increment
	for(m=0;m<=10;m++){
	  mas = DeltaMSq32 - (5)*Mincrement + (m)*Mincrement;
	  
          //Mass increment for Andy's combined surface
          //for(m=0;m<7;m++){
          //mas = -2.70e-3 + m*0.1e-3;
	  
	  cout << "Finding minchi2 for: Delta = " << das << " Theta23= " << aas << " Theta13 = " << ras << " with dms= " << mas << endl;
	  for(ie=0;ie<ExtrapFHC.size();ie++)
	    {
	      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,mas);
	      ExtrapFHC[ie]->SetOscPar(OscPar::kTh13,ras);
	      ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,aas); 
	      ExtrapFHC[ie]->SetOscPar(OscPar::kDelta,das);
	      //Th12 and DeltaM12 set for Andy's surface
	      //ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,0.587);
	      //ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,7.54e-5);
	      ExtrapFHC[ie]->OscillatePrediction();
	    }
	  
	  for(ie=0;ie<ExtrapRHC.size();ie++)
	    {
	      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,mas);
	      ExtrapRHC[ie]->SetOscPar(OscPar::kTh13,ras);
	      ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,aas);
	      ExtrapRHC[ie]->SetOscPar(OscPar::kDelta,das);
	      //Th12 and DeltaM12 set for Andy's surface
	      //ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,0.587);
	      //ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,7.54e-5);
	      ExtrapRHC[ie]->OscillatePrediction();
	    }
	  
	  Bkgd->Reset();
	  Sig->Reset();
	  NExp->Reset();
	  
	  CombineFHCRHC();
	  NExp->Add(Bkgd);
	  NExp->Add(Sig);
	  
	  penpar = ((TMath::Sin(2*aas))*(TMath::Sin(2*aas)) - ((TMath::Sin(2*Theta23))*(TMath::Sin(2*Theta23))))*((TMath::Sin(2*aas))*(TMath::Sin(2*aas)) - ((TMath::Sin(2*Theta23))*(TMath::Sin(2*Theta23))))/(0.036*0.036);
	  
	  f = BinLikelihood() + penpar + (ras - Theta13)*(ras-Theta13)/(Theta13Unc*Theta13Unc) + (mas - DeltaMSq32)*(mas-DeltaMSq32)/(DeltaMSq32Unc*DeltaMSq32Unc);
	  
	  //outfile for Andy's surface - to swap back, only use BinLikelihood above and remove penalty terms
	  //outfile << mas << " " << ass << " " << rss << " " << das << " " << f << endl;    
	  if(f<minchi2){
	    minchi2=f;
	    dbest = das;
	    ttbest = aas;
	    otbest = ras;
	    dmsbest = mas;
	  }
	}//m
      }//r
    }//a
  }//d
  
  // gm2 for -2LogL   
  global_minchi2 = 0; 
  
  cout << "Final Points - Theta13: " << otbest << " Theta23: " << ttbest << " with DMS: " << dmsbest << endl;  

  double plotinc = 20;
  int nn = int(2*plotinc+1); // number of values in vector deltaval. We'' use the same number to fill TH1* histogram.
  double deltaval[41];
  for(int i=0;i<nn;i++){
    deltaval[i]=i*(1/plotinc)*TMath::Pi();
  } 
  
  double deltatwologl[41];
  
  double dp;
  for(i=0; i<2*plotinc+1; i++)//looping over delta again, but now only the points we want to display on the delta plot
    {
      dp=deltaval[i];
      cout << "It: " << i << " with  delta: " << dp << endl;
      for(ie=0;ie<ExtrapFHC.size();ie++)
	{
	  ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,dmsbest);
	  ExtrapFHC[ie]->SetOscPar(OscPar::kTh13,otbest);
	  ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,ttbest);
	  ExtrapFHC[ie]->SetOscPar(OscPar::kDelta,dp);
	  ExtrapFHC[ie]->OscillatePrediction();
	}
      for(ie=0;ie<ExtrapRHC.size();ie++)
	{
	  ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,dmsbest);
	  ExtrapRHC[ie]->SetOscPar(OscPar::kTh13,otbest);
	  ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,ttbest);
	  ExtrapRHC[ie]->SetOscPar(OscPar::kDelta,dp);
	  ExtrapRHC[ie]->OscillatePrediction();
	}
      
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      
      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      delchi2 = 1e10;
      if(FitMethod==0)
	{
	  delchi2 = PoissonChi2(NExp) - global_minchi2;
	}
      else if(FitMethod==1)
	{
	  delchi2 = ScaledChi2(Bkgd,Sig) - global_minchi2;
	}
      else if(FitMethod==2)
	{
	  delchi2 = StandardChi2(NExp) - global_minchi2;
	}
      else if(FitMethod==3)
	{
	  //Likelihood: "Standard" (N syst, N nuisance)
	  //Calculate the likelihood (x2 for chi)
	  delchi2 = StandardLikelihood() - global_minchi2;
	}
      else if(FitMethod==4)
	{
	  //Likelihood: Bin by Bin Calculation of Systematics
	  //Calculate the likelihood (x2 for chi)
	  delchi2 = BinLikelihood() - global_minchi2;
	}
      else
	{
	  cout<<"Error in GetLikelihood(): Unknown 'FitMethod'."<<endl;
	}
      deltatwologl[i] = delchi2;
    }//end loop over delta values
  
  //make a TGraph
  TFile *w = new TFile("DeltaOutput.root","RECREATE");
  TGraph *g = new TGraph(nn,deltaval,deltatwologl);
  //g->Draw("al");
  g->SetName("deltasweep");
  g->Write();
  w->Close();
  //outfile for surface 
  //outfile.close();
  return;
}

///////////////////////////////Adam NSI LAND////////////////////////////////////////////////////////////////////////

void NueFit2D_Joint::DoNSIFitQuick()
{
  //probably don't want to set central values in the function call
  //make separate functions that also set the uncertainties, ie
  //NueFit2D::SetTheta13(theta13 central value, +1 sigma error, -1sigam error)
  if(NObs==0)
    {
      cout<<"NObs not set.  Quitting..."<<endl;
      return;
    }
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
    {
      cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
      return;
    }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
    {
      cout<<"No Extrapolate2D input.  Quitting..."<<endl;
      return;
    }
  
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();
  
  Bkgd = (TH1D*)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)NObs->Clone("Sig");
  Sig->Reset();
  TH1D *NExp = (TH1D*)NObs->Clone("NExp");
  NExp->Reset();
  
  unsigned int ie;
  for(ie=0;ie<ExtrapFHC.size();ie++)
    {
      ExtrapFHC[ie]->GetPrediction();
      /*
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,Theta12);
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,Theta23);
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh13,Theta13);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq32);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq21);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDelta,delta);

      // Add NSI parameters:
      ExtrapFHC[ie]->SetOscPar(OscPar::kEps_ee, Eps_ee);
      ExtrapFHC[ie]->SetOscPar(OscPar::kEps_emu, Eps_emu);
      ExtrapFHC[ie]->SetOscPar(OscPar::kEps_etau, Eps_etau);
      ExtrapFHC[ie]->SetOscPar(OscPar::kEps_mumu, Eps_mumu);
      ExtrapFHC[ie]->SetOscPar(OscPar::kEps_mutau, Eps_mutau); 
      ExtrapFHC[ie]->SetOscPar(OscPar::kEps_tautau, Eps_tautau); 
      ExtrapFHC[ie]->SetOscPar(OscPar::kDelta_emu, Delta_emu);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDelta_etau, Delta_etau);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDelta_mutau, Delta_mutau);
      //Oscillate prediction
      ExtrapFHC[ie]->OscillatePrediction();
      */
    }
  for(ie=0;ie<ExtrapRHC.size();ie++)
    {
      ExtrapRHC[ie]->GetPrediction();
      /*
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,Theta12);
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,Theta23);
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh13,Theta13);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq32);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq21);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDelta,delta);      
      // Add NSI parameters:
      ExtrapRHC[ie]->SetOscPar(OscPar::kEps_ee, Eps_ee);
      ExtrapRHC[ie]->SetOscPar(OscPar::kEps_emu, Eps_emu);
      ExtrapRHC[ie]->SetOscPar(OscPar::kEps_etau, Eps_etau);
      ExtrapRHC[ie]->SetOscPar(OscPar::kEps_mumu, Eps_mumu);
      ExtrapRHC[ie]->SetOscPar(OscPar::kEps_mutau, Eps_mutau); 
      ExtrapRHC[ie]->SetOscPar(OscPar::kEps_tautau, Eps_tautau); 
      ExtrapRHC[ie]->SetOscPar(OscPar::kDelta_emu, Delta_emu);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDelta_etau, Delta_etau);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDelta_mutau, Delta_mutau);
      // Oscillate Prediction
      ExtrapRHC[ie]->OscillatePrediction();
      */
    }
  
  Int_t i;  
  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  //Double_t ss2th13 = 0;
  Double_t delchi2 = 0;
  Double_t best = 0;
  //
  double chi2[601],val[601]; 
  TGraph *g3;
  TF1* fit;
  double minchi2;
  double global_minchi2;
  
  //find global minimum  
  //find minimum in delta
  best=-1.;
  minchi2=1.0e6;

  double inc = 600.0;
  for(i=0;i<=inc;i++)
    {
      chi2[i]=-1.;
      val[i]=-1.;
    }
  
  double epsetauincrement = (12.0/inc);
  double nn = 601; 
  double epsetaulog[601];
  double epsetauval[601];
  for(int thing=0;thing<=inc;thing++){
    epsetauval[thing]=thing*epsetauincrement-6.0;
  }

  double db;
  for(i=0; i<inc+1 ; i++)//i scans over all delta values so we can find the best fit delta
    {
      /*
      chi2[0] = chi2[1];
      chi2[1] = chi2[2];
      val[0] = val[1];
      val[1] = val[2];
      */
      db = i*epsetauincrement - 6.0;
  
      cout << "It: " << i << " with eps_etau: " << db << endl;
      for(ie=0;ie<ExtrapFHC.size();ie++)
	{
	  ExtrapFHC[ie]->SetOscPar(OscPar::kEps_etau,db);
	  ExtrapFHC[ie]->OscillatePrediction();
	}
      
      for(ie=0;ie<ExtrapRHC.size();ie++)
	{
	  ExtrapRHC[ie]->SetOscPar(OscPar::kEps_etau,db);
	  ExtrapRHC[ie]->OscillatePrediction();
	}
     
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      
      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      

      if(FitMethod==0)
	{
	  chi2[i] = PoissonChi2(NExp);
	}
      else if(FitMethod==1)
	{
	  chi2[i] = ScaledChi2(Bkgd,Sig);
	}
      else if(FitMethod==2)
	{
	  chi2[i] = StandardChi2(NExp);
	}
      else if(FitMethod==3)
	{
	  //Likelihood: "Standard" (N syst, N nuisance)
	  //Calculate the likelihood (x2 for chi)
	  chi2[i] = StandardLikelihood();
	}
      else if(FitMethod==4)
	{
	  //Likelihood: Bin by Bin Calculation of Systematics
	  //Calculate the likelihood (x2 for chi)
	  chi2[i] = BinLikelihood();
	  //cout << "ss2th13 = " << ss2th13 << " Chi2=" << chi2[2] << endl;
	}
      else
	{
	  cout<<"Error in GetLikelihood(): Unknown 'FitMethod'."<<endl;
	}
      
      val[i] = db;


      /*      
      if(i<2) continue;
      
      if(i<3 && chi2[2]>chi2[1] && chi2[1]>chi2[0])//first three points are increasing, first point is minimum.
	{
	  best = val[0];
	  minchi2 = chi2[0];
	  cout<<"minimum at 1st point: "<<minchi2<<" at "<<best<<endl;
	  break;
	}
      
      if(chi2[2]>chi2[1] && chi2[0]>chi2[1])//found minimum
	{
	  g3 = new TGraph(3, val, chi2);
	  fit = new TF1("pol2", "pol2");
	  g3->Fit(fit, "Q");//fit to second order polynominal
	  if(fit->GetParameter(2) > 0)//if the x^2 term is nonzero
	    {
	      best = -fit->GetParameter(1)/(2*fit->GetParameter(2));//the location of the minimum is -p1/(2*p2)
	      minchi2 = fit->GetParameter(0) + fit->GetParameter(1)*best + fit->GetParameter(2)*best*best;
	      cout<<"minimum with fit: "<<minchi2<<" at "<<best<<endl;
	    }
	  else//if the x^2 term is zero, then just use the minimum you got by scanning
	    {
	      best = val[1];
	      minchi2 = chi2[1];
	      cout<<"minimum with scan: "<<minchi2<<" at "<<best<<endl;
	    }
	  break;
	}
      */
    }
  int idx;
 
  for(int bst=0; bst<=inc; bst++){
    if(chi2[bst]<minchi2){
      cout << "Chi2 of " <<chi2[bst]<< " at " <<val[bst] <<endl;
      minchi2 = chi2[bst];
      cout << "Minchi2 is " <<minchi2<<endl;
      best = val[bst];
      idx = bst;
    }
  }

  if(idx>1 && idx<inc-1){
      double va[3]; va[0]=val[idx-2]; va[1]=val[idx-1]; va[2]=val[idx];
      double ch[3]; ch[0]=chi2[idx-2]; ch[1]=chi2[idx-1]; ch[2]=chi2[idx];
      g3 = new TGraph(3, va, ch);
      fit = new TF1("pol2", "pol2");
      g3->Fit(fit, "Q");//fit to second order polynominal
          if(fit->GetParameter(2) > 0)//if the x^2 term is nonzero
	    {
	      best = -fit->GetParameter(1)/(2*fit->GetParameter(2));//the location of the minimum is -p1/(2*p2)
	      minchi2 = fit->GetParameter(0) + fit->GetParameter(1)*best + fit->GetParameter(2)*best*best;
	      cout<<"minimum with fit: "<<minchi2<<" at "<<best<<endl;
	    }
	  else//if the x^2 term is zero, then just use the minimum you got by scanning
	    {
	      cout<<"Keeping scan best of: "<<minchi2<<" at "<<best<<endl;
	    }
	 
    }
 

  cout << "Fit 1 Complete" << endl;
  minchi2=1e10;
 global_minchi2=0;
 

  double dp;
  for(i=0; i<inc+1; i++){
    //looping over eps_Etau again, but now only the points we want to display on the eps plot
    
      dp=epsetauval[i];
      cout << "It: " << i << " with  eps_etau: " << dp << endl;
      for(ie=0;ie<ExtrapFHC.size();ie++)
	{
	  ExtrapFHC[ie]->SetOscPar(OscPar::kEps_etau,dp);
	  ExtrapFHC[ie]->OscillatePrediction();
	}
      for(ie=0;ie<ExtrapRHC.size();ie++)
	{
	  ExtrapRHC[ie]->SetOscPar(OscPar::kEps_etau,dp);
	  ExtrapRHC[ie]->OscillatePrediction();
	}
      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();
      
      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      delchi2 = 1e10;
      if(FitMethod==0)
	{
	  delchi2 = PoissonChi2(NExp) - global_minchi2;
	}
      else if(FitMethod==1)
	{
	  delchi2 = ScaledChi2(Bkgd,Sig) - global_minchi2;
	}
      else if(FitMethod==2)
	{
	  delchi2 = StandardChi2(NExp) - global_minchi2;
	}
      else if(FitMethod==3)
	{
	  //Likelihood: "Standard" (N syst, N nuisance)
	  //Calculate the likelihood (x2 for chi)
	  delchi2 = StandardLikelihood() - global_minchi2;
	}
      else if(FitMethod==4)
	{
	  //Likelihood: Bin by Bin Calculation of Systematics
	  //Calculate the likelihood (x2 for chi)
	  delchi2 = BinLikelihood() - global_minchi2;
	}
      else
	{
	  cout<<"Error in GetLikelihood(): Unknown 'FitMethod'."<<endl;
	}
      epsetaulog[i] = delchi2;
    }//end loop over delta values
  cout <<"Fit 2 Complete" << endl;
  //make a TGraph
  TFile *w = new TFile("NSIQuick.root","RECREATE");
  cout <<"Outfile Opened" << endl;
  TGraph *g = new TGraph(nn,epsetauval,epsetaulog);
  g->SetName("etausweep");
  g->GetXaxis()->SetTitle("#epsilon_{e#tau}");
  g->GetYaxis()->SetTitle("-2lnL");
  g->Write();
  w->Close();
   return;
}
///////////////////////////////////////////////////////////////////////////////////////
void NueFit2D_Joint::DoNSIFitForever()
{
  //probably don't want to set central values in the function call
  //make separate functions that also set the uncertainties, ie
  //NueFit2D::SetTheta13(theta13 central value, +1 sigma error, -1sigam error)
    
  if(NObs==0)
    {
      cout<<"NObs not set.  Quitting..."<<endl;
      return;
    }
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
    {
      cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
      return;
    }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
    {
      cout<<"No Extrapolate2D input.  Quitting..."<<endl;
      return;
    }
  
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();
  
  Bkgd = (TH1D*)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)NObs->Clone("Sig");
  Sig->Reset();
  TH1D *NExp = (TH1D*)NObs->Clone("NExp");
  NExp->Reset();
  
  unsigned int ie;
  for(ie=0;ie<ExtrapFHC.size();ie++)
    {
      ExtrapFHC[ie]->GetPrediction();
    }
  for(ie=0;ie<ExtrapRHC.size();ie++)
    {
      ExtrapRHC[ie]->GetPrediction();
    }
  
  
  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  Int_t i,j,k,q,s;
  i = 1;
  Int_t idx = 0;
  Double_t epsetau;
  Double_t delta = 0;
  Double_t Deltaincrement = 0;
  if(nDeltaSteps>0) Deltaincrement = (DeltaHigh - DeltaLow)/(nDeltaSteps);
  Double_t epsetauincrement = 0;
  if(nepsetauSteps>0) epsetauincrement = (epsetauHigh - epsetauLow)/(nepsetauSteps);
  int arraysize = (nDeltaSteps+1)*(nepsetauSteps+1);  
  cout << arraysize << endl;
  double x[18045],y[18045],chi[18045]; 
  //double x[53489],y[53489],chi[53489];  
  TF1* fit;
  double minchi2;
  double mc2;
  
  if(NHIH==1){
  for(ie=0;ie<ExtrapFHC.size();ie++) ExtrapFHC[ie]->InvertMassHierarchy();
  for(ie=0;ie<ExtrapRHC.size();ie++) ExtrapRHC[ie]->InvertMassHierarchy();
  }//DoIH function sets inverted hierarchy
  
 for(j=0;j<nDeltaSteps+1;j++)
  {
    delta = j*Deltaincrement*TMath::Pi() + DeltaLow*TMath::Pi();
    for(ie=0;ie<ExtrapFHC.size();ie++){
      ExtrapFHC[ie]->SetDelta_etau(delta);
    }
    for(ie=0;ie<ExtrapRHC.size();ie++){
      ExtrapRHC[ie]->SetDelta_etau(delta);
    }
    
    minchi2=100;
    mc2=1000;

    for(q=0;q<=nepsetauSteps;q++){
      if(idx>arraysize){
	cout << "The code is not doing what you think it is doing, check arraysize and idx." << endl;
       }
      epsetau = epsetauLow+q*epsetauincrement;
      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetEps_etau(epsetau);
        ExtrapFHC[ie]->OscillatePrediction();
      }
      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetEps_etau(epsetau);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      if(FitMethod==0)
      {
        x[idx] = epsetau;
        y[idx] = delta/TMath::Pi();
        chi[idx] = PoissonChi2(NExp);
      }
      else if(FitMethod==1)
      {
        x[idx] = epsetau;
        y[idx] = delta/TMath::Pi();
        chi[idx] = ScaledChi2(Bkgd,Sig);
      }
      else if(FitMethod==2)
      {
        x[idx] = epsetau;
        y[idx] = delta/TMath::Pi();
        chi[idx] = StandardChi2(NExp);
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        x[idx] = epsetau;
        y[idx] = delta/TMath::Pi();
        chi[idx] = StandardLikelihood();
      }
      else if(FitMethod==4)

      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        //cout << idx << " " << epsetau << " " << delta << endl;
        x[idx] = epsetau;
        y[idx] = delta/TMath::Pi();
        chi[idx] = BinLikelihood();
      }
      else
      {
        cout<<"Error in RunNSIDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
    
    idx++;    
  
    }//epsetauloop 
 
  }//deltaloop
 
 double dbest[1];
 double epsbest[1];
 for(k=0;k<arraysize;k++){
   if(chi[k]<mc2){
     mc2=chi[k];
     dbest[0]=y[k];
     epsbest[0]=x[k];
   }
 }
 cout << "Subtracting minchi value of: " << mc2 << " from coordinate ("<<epsbest[0]<<","<<dbest[0]<<"Pi)"<< endl;

 for(s=0;s<arraysize;s++){
   chi[s]=chi[s]-mc2;
   if(abs(chi[s])<0.00001){
     cout <<"I have a zero!"<<endl;
   }
   
 }

  //make a TGraph
  TFile *w = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  TGraph *bp = new TGraph(i,epsbest,dbest);
  TGraph2D *g = new TGraph2D(arraysize,x,y,chi);
  g->GetYaxis()->SetTitle("#delta_{e#tau}/#pi");
  g->GetXaxis()->SetTitle("#epsilon_{e#tau}");
  g->GetZaxis()->SetTitle("-2#DeltalnL");
  bp->SetName("bp");
  g->SetName("ets");
  g->Write();
  bp->Write();
  w->Close();
  return;
}

void NueFit2D_Joint::GetNSIFFX2()
{
  //probably don't want to set central values in the function call
  //make separate functions that also set the uncertainties, ie
  //NueFit2D::SetTheta13(theta13 central value, +1 sigma error, -1sigam error)
    
  if(NObs==0)
    {
      cout<<"NObs not set.  Quitting..."<<endl;
      return;
    }
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
    {
      cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
      return;
    }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
    {
      cout<<"No Extrapolate2D input.  Quitting..."<<endl;
      return;
    }
  
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();
  
  Bkgd = (TH1D*)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)NObs->Clone("Sig");
  Sig->Reset();
  TH1D *NExp = (TH1D*)NObs->Clone("NExp");
  NExp->Reset();
  
  unsigned int ie;
  for(ie=0;ie<ExtrapFHC.size();ie++)
    {
      ExtrapFHC[ie]->GetPrediction();
    }
  for(ie=0;ie<ExtrapRHC.size();ie++)
    {
      ExtrapRHC[ie]->GetPrediction();
    }
  
  
  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  Int_t i,j,k,q,s;
  i = 1;
  Int_t idx = 0;
  Double_t epsetau;
  Double_t delta = 0;
  Double_t Deltaincrement = 0;
  if(nDeltaSteps>0) Deltaincrement = (DeltaHigh - DeltaLow)/(nDeltaSteps);
  Double_t epsetauincrement = 0;
  if(nepsetauSteps>0) epsetauincrement = (epsetauHigh - epsetauLow)/(nepsetauSteps);
  int arraysize = 1;  
  cout << arraysize << endl;
  double x[1],y[1],chi[1]; 
   TF1* fit;
  double minchi2;
  double mc2;
  
  if(NHIH==1){
  for(ie=0;ie<ExtrapFHC.size();ie++) ExtrapFHC[ie]->InvertMassHierarchy();
  for(ie=0;ie<ExtrapRHC.size();ie++) ExtrapRHC[ie]->InvertMassHierarchy();
  }//DoIH function sets inverted hierarchy
  
 for(j=0;j<nDeltaSteps+1;j++)
  {
    delta = j*Deltaincrement*TMath::Pi() + DeltaLow*TMath::Pi();
    for(ie=0;ie<ExtrapFHC.size();ie++){
      ExtrapFHC[ie]->SetDelta_etau(delta);
    }
    for(ie=0;ie<ExtrapRHC.size();ie++){
      ExtrapRHC[ie]->SetDelta_etau(delta);
    }
    
    minchi2=100;
    mc2=1000;

    for(q=0;q<=nepsetauSteps;q++){
      if(idx>arraysize){
	cout << "The code is not doing what you think it is doing, check arraysize and idx." << endl;
       }
      epsetau = epsetauLow+q*epsetauincrement;
      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetEps_etau(epsetau);
        ExtrapFHC[ie]->OscillatePrediction();
      }
      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetEps_etau(epsetau);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      if(FitMethod==0)
      {
        x[idx] = epsetau;
        y[idx] = delta/TMath::Pi();
        chi[idx] = PoissonChi2(NExp);
      }
      else if(FitMethod==1)
      {
        x[idx] = epsetau;
        y[idx] = delta/TMath::Pi();
        chi[idx] = ScaledChi2(Bkgd,Sig);
      }
      else if(FitMethod==2)
      {
        x[idx] = epsetau;
        y[idx] = delta/TMath::Pi();
        chi[idx] = StandardChi2(NExp);
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        x[idx] = epsetau;
        y[idx] = delta/TMath::Pi();
        chi[idx] = StandardLikelihood();
      }
      else if(FitMethod==4)

      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        //cout << idx << " " << epsetau << " " << delta << endl;
        x[idx] = epsetau;
        y[idx] = delta/TMath::Pi();
        chi[idx] = BinLikelihood();
      }
      else
      {
        cout<<"Error in RunNSIDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
    
    idx++;    
  
    }//epsetauloop 
 
  }//deltaloop
 
 double dbest[1];
 double epsbest[1];
 for(k=0;k<arraysize;k++){
   if(chi[k]<mc2){
     mc2=chi[k];
     dbest[0]=y[k];
     epsbest[0]=x[k];
   }
 }
 cout << "Subtracting minchi value of: " << mc2 << " from coordinate ("<<epsbest[0]<<","<<dbest[0]<<"Pi)"<< endl;

 for(s=0;s<arraysize;s++){
   chi[s]=chi[s]-mc2;
   if(abs(chi[s])<0.00001){
     cout <<"I have a zero!"<<endl;
   }
   
 }

}


////////////////////// Adam's NSI Function Rewrite //////////////////
///////////////FUNCTION IS OBSOLETE - USE DoNSIFitForever////////////
void NueFit2D_Joint::RunNSIDeltaChi2Contour(int cl)
{
  if(NObs==0)
  {
    cout<<"NObs not set.  Quitting..."<<endl;
    return;
  }
  if(FitMethod==1 && (FracErr_Bkgd==0 || FracErr_Sig==0))
  {
    cout<<"FracErr_Bkgd and FracErr_Sig need to be set for ScaledChi2.  Quitting..."<<endl;
    return;
  }
  if((ExtrapFHC.size()+ExtrapRHC.size())==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return;
  }
  
  if(FitMethod==3) DefineStdDlnLMinuit();
  if(FitMethod==4) DefineBinDlnLMinuit();
  
  Bkgd = (TH1D*)NObs->Clone("Bkgd");
  Bkgd->Reset();
  Sig = (TH1D*)NObs->Clone("Sig");
  Sig->Reset();
  NExp = (TH1D*)NObs->Clone("NExp");
  NExp->Reset();
  
  unsigned int ie;
  for(ie=0;ie<ExtrapFHC.size();ie++) ExtrapFHC[ie]->GetPrediction();
  for(ie=0;ie<ExtrapRHC.size();ie++) ExtrapRHC[ie]->GetPrediction();
  
  if(ErrorMatrix==0) ErrorMatrix = new TH2D("ErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  if(InvErrorMatrix==0) InvErrorMatrix = new TH2D("InvErrorMatrix","",nBins,-0.5,nBins-0.5,nBins,-0.5,nBins-0.5);
  
  Int_t i,j,k,q,s;
  
  Double_t epsetau;
  Double_t delta = 0;
  Double_t Deltaincrement = 0;
  if(nDeltaSteps>0) Deltaincrement = (DeltaHigh - DeltaLow)/(nDeltaSteps);
  Double_t epsetauincrement = 0;
  if(nepsetauSteps>0) epsetauincrement = (epsetauHigh - epsetauLow)/(nepsetauSteps);
  
  double chi2[3],val[3];
  double ch[61],va[61];
  TGraph *g3;
  TF1* fit;
  double best;
  double vbt;
  double minchi2;
  double mc2;
  
  double limit;
  double delchi2;
  double sprev,delchi2prev;
  double contourlvl = 0;
  if(cl==0)//90% CL
  {
    contourlvl = 2.71;
  }
  else if(cl==1)//68.3% cL
  {
    contourlvl = 1.0;
  }
  else if(cl==2)//95% cL
  {
    contourlvl = 3.84;
  }
  else if(cl==3)//99% cL
  {
    contourlvl = 6.63;
  }
  else
  {
    cout<<"Error in RunNSIDeltaChi2Contour(): Input value should be 0,1,2,or 3 for 90%, 68.3%, 95%, or 99%.  Quitting..."<<endl;
    return;
  }
  
  cout<<"Seeking ";
  if(cl==0) cout<<"90% ";
  else if(cl==1)cout<<"68% ";
  else if(cl==2)cout<<"95% ";
  else if(cl==3)cout<<"99% ";
  else cout<< "??? ";
  cout<<" CL limits"<<endl;
  cout<<"PhaseBroom set to "<<PhaseBroom<<endl;
  
  vector<double> deltapts;
  vector<double> bestfit_norm;
  vector<double> bestfit_invt;
  vector<double> limit_norm;
  vector<double> limit_invt;
  vector<double> limit_normlow;
  vector<double> limit_invtlow;
  
  for(j=0;j<nDeltaSteps+1;j++)
  {
    delta = j*Deltaincrement*TMath::Pi() + DeltaLow;
    for(ie=0;ie<ExtrapFHC.size();ie++){
      ExtrapFHC[ie]->SetDeltaCP((1.0-PhaseBroom)*delta);
      ExtrapFHC[ie]->SetDelta_etau(PhaseBroom*delta);
    }
    for(ie=0;ie<ExtrapRHC.size();ie++){
      ExtrapRHC[ie]->SetDeltaCP((1.0-PhaseBroom)*delta);
      ExtrapRHC[ie]->SetDelta_etau(PhaseBroom*delta);
    }
    deltapts.push_back(delta/TMath::Pi());
    
    best=-1.;
    minchi2=100;
    mc2=1000;
    for(i=0;i<3;i++)
    {
      chi2[i]=-1.;
      val[i]=-1.;
    }
    for(k=0;k<=60;k++){
      ch[k]=-1.;
      va[k]=-1.;
    }
    for(q=0;q<=60;q++){
      epsetau = epsetauLow+q*0.1;
      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetEps_etau(epsetau);
        ExtrapFHC[ie]->OscillatePrediction();
      }
      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetEps_etau(epsetau);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      if(FitMethod==0)
      {
        ch[q] = PoissonChi2(NExp);
      }
      else if(FitMethod==1)
      {
        ch[q] = ScaledChi2(Bkgd,Sig);
      }
      else if(FitMethod==2)
      {
        ch[q] = StandardChi2(NExp);
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        ch[q] = StandardLikelihood();
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        ch[q] = BinLikelihood();
      }
      else
      {
        cout<<"Error in RunNSIDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
     
      va[q]=epsetau;
    }//coarse sweep

    for(s=0;s<61;s++){
      if(ch[s]<mc2){
	mc2=ch[s];
        vbt=va[s];
      }
    }//find junk
    cout << "Found CBF for Delta: "<<delta<< " at " <<vbt<< " with chi2 of "<<mc2<< endl;
    for(i=0;i<=nepsetauSteps/5;i++)
    {
      chi2[0] = chi2[1];
      chi2[1] = chi2[2];
      val[0] = val[1];
      val[1] = val[2];
      epsetau = vbt + (-(nepsetauSteps/10)+i)*epsetauincrement;
 
      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetEps_etau(epsetau);
        ExtrapFHC[ie]->OscillatePrediction();
      }
      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetEps_etau(epsetau);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      chi2[2] = 1e10;
      if(FitMethod==0)
      {
        chi2[2] = PoissonChi2(NExp);
      }
      else if(FitMethod==1)
      {
        chi2[2] = ScaledChi2(Bkgd,Sig);
      }
      else if(FitMethod==2)
      {
        chi2[2] = StandardChi2(NExp);
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        chi2[2] = StandardLikelihood();
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        chi2[2] = BinLikelihood();
      }
      else
      {
        cout<<"Error in RunNSIDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
      
      val[2] = epsetau;
      
      if(i<2) continue;
      
      if(i<3 && chi2[2]>chi2[1] && chi2[1]>chi2[0])//first three points are increasing, first point is minimum.
      {
	best = val[0];
	minchi2 = chi2[0];
	cout<<"minimum at 1st point: "<<minchi2<<" at "<<best<<endl;
	break;
      }
      
      if(chi2[2]>chi2[1] && chi2[0]>chi2[1])//found minimum
      {
        g3 = new TGraph(3, val, chi2);
        fit = new TF1("pol2", "pol2");
        g3->Fit(fit, "Q");//fit to second order polynominal
        if(fit->GetParameter(2) > 0)//if the x^2 term is nonzero
        {
          best = -fit->GetParameter(1)/(2*fit->GetParameter(2));//the location of the minimum is -p1/(2*p2)
          minchi2 = fit->GetParameter(0) + fit->GetParameter(1)*best + fit->GetParameter(2)*best*best;
          cout<<"minimum with fit: "<<minchi2<<" at "<<best<<endl;
        }
        else//if the x^2 term is zero, then just use the minimum you got by scanning
        {
          best = val[1];
          minchi2 = chi2[1];
          cout<<"minimum with scan: "<<minchi2<<" at "<<best<<endl;
        }
        break;
      }
    }
    bestfit_norm.push_back(best);
    
    limit = -1.;
    delchi2prev = 1000;
    sprev = best;
    for(i=0;i<(epsetauHigh-best)/epsetauincrement;i++)
    {
      //epsetau = i*epsetauincrement + epsetauLow;
      epsetau = i*epsetauincrement + best;
      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetEps_etau(epsetau);
        ExtrapFHC[ie]->OscillatePrediction();
      }
      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetEps_etau(epsetau);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      delchi2 = 1e10;
      if(FitMethod==0)
      {
        delchi2 = PoissonChi2(NExp) - minchi2;
      }
      else if(FitMethod==1)
      {
        delchi2 = ScaledChi2(Bkgd,Sig) - minchi2;
      }
      else if(FitMethod==2)
      {
        delchi2 = StandardChi2(NExp) - minchi2;
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        delchi2 = StandardLikelihood() - minchi2;
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        delchi2 = BinLikelihood() - minchi2;
      }
      else
      {
        cout<<"Error in RunNSIDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
      
      if(i==1) continue;
      
      if(delchi2>contourlvl && delchi2prev<contourlvl)
      {
        limit = epsetau + ((epsetau-sprev)/(delchi2 - delchi2prev))*(contourlvl - delchi2);
        cout<<delchi2prev<<", "<<sprev<<", "<<delchi2<<", "<<epsetau<<", "<<limit<<endl;
        break;
      }
      delchi2prev = delchi2;
      sprev = epsetau;
    }
    limit_norm.push_back(limit);

    limit = -1.;
    delchi2prev = 1000;
    sprev = best;
    for(i=0;i<(best-epsetauLow)/epsetauincrement;i++)
    {
      //epsetau = i*epsetauincrement + epsetauLow;
      epsetau = -i*epsetauincrement + best;
      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetEps_etau(epsetau);
        ExtrapFHC[ie]->OscillatePrediction();
      }
      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetEps_etau(epsetau);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      delchi2 = 1e10;
      if(FitMethod==0)
      {
        delchi2 = PoissonChi2(NExp) - minchi2;
      }
      else if(FitMethod==1)
      {
        delchi2 = ScaledChi2(Bkgd,Sig) - minchi2;
      }
      else if(FitMethod==2)
      {
        delchi2 = StandardChi2(NExp) - minchi2;
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        delchi2 = StandardLikelihood() - minchi2;
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        delchi2 = BinLikelihood() - minchi2;
      }
      else
      {
        cout<<"Error in RunNSIDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
      
      if(i==1) continue;
      
      if(delchi2>contourlvl && delchi2prev<contourlvl)
      {
        limit = epsetau + ((epsetau-sprev)/(delchi2 - delchi2prev))*(contourlvl - delchi2);
        cout<<delchi2prev<<", "<<sprev<<", "<<delchi2<<", "<<epsetau<<", "<<limit<<endl;
        break;
      }
      delchi2prev = delchi2;
      sprev = epsetau;
    }
    limit_normlow.push_back(limit);
  }
  
  for(ie=0;ie<ExtrapFHC.size();ie++) ExtrapFHC[ie]->InvertMassHierarchy();
  for(ie=0;ie<ExtrapRHC.size();ie++) ExtrapRHC[ie]->InvertMassHierarchy();

  
  for(j=0;j<nDeltaSteps+1;j++)
  {
    delta = j*Deltaincrement*TMath::Pi() + DeltaLow;
    for(ie=0;ie<ExtrapFHC.size();ie++){
      ExtrapFHC[ie]->SetDeltaCP((1.0-PhaseBroom)*delta);
      ExtrapFHC[ie]->SetDelta_etau(PhaseBroom*delta);
    }
    for(ie=0;ie<ExtrapRHC.size();ie++){
      ExtrapRHC[ie]->SetDeltaCP((1.0-PhaseBroom)*delta);
      ExtrapRHC[ie]->SetDelta_etau(PhaseBroom*delta);
    }

    best=-1.;
    minchi2=100;
    mc2=1000;
    for(i=0;i<3;i++)
    {
      chi2[i]=-1.;
      val[i]=-1.;
    }
    for(k=0;k<=60;k++){
      ch[k]=-1.;
      va[k]=-1.;
    }
    for(q=0;q<=60;q++){
      epsetau = epsetauLow+q*0.1;
      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetEps_etau(epsetau);
        ExtrapFHC[ie]->OscillatePrediction();
      }
      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetEps_etau(epsetau);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      if(FitMethod==0)
      {
        ch[q] = PoissonChi2(NExp);
      }
      else if(FitMethod==1)
      {
        ch[q] = ScaledChi2(Bkgd,Sig);
      }
      else if(FitMethod==2)
      {
        ch[q] = StandardChi2(NExp);
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        ch[q] = StandardLikelihood();
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        ch[q] = BinLikelihood();
      }
      else
      {
        cout<<"Error in RunNSIDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
     
      va[q]=epsetau;
    }//coarse sweep

    for(s=0;s<61;s++){
      if(ch[s]<mc2){
	mc2=ch[s];
        vbt=va[s];
      }
    }//find junk
    cout << "Found CBF for Delta: "<<delta<< " at " <<vbt<< " with chi2 of "<<mc2<< endl;
    for(i=0;i<=nepsetauSteps/5;i++)
    {
      chi2[0] = chi2[1];
      chi2[1] = chi2[2];
      val[0] = val[1];
      val[1] = val[2];
      epsetau = vbt + (-(nepsetauSteps/10)+i)*epsetauincrement;

      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetEps_etau(epsetau);
        ExtrapFHC[ie]->OscillatePrediction();
      }

      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetEps_etau(epsetau);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      chi2[2] = 1e10;
      if(FitMethod==0)
      {
        chi2[2] = PoissonChi2(NExp);
      }
      else if(FitMethod==1)
      {
        chi2[2] = ScaledChi2(Bkgd,Sig);
      }
      else if(FitMethod==2)
      {
        chi2[2] = StandardChi2(NExp);
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        chi2[2] = StandardLikelihood();
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        chi2[2] = BinLikelihood();
      }
      else
      {
        cout<<"Error in RunDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
      
      val[2] = epsetau;
      
      if(i<2) continue;
      
      if(i<3 && chi2[2]>chi2[1] && chi2[1]>chi2[0])//first three points are increasing, first point is minimum.
      {
	best = val[0];
	minchi2 = chi2[0];
	cout<<"minimum at 1st point: "<<minchi2<<" at "<<best<<endl;
	break;
      }
      
      if(chi2[2]>chi2[1] && chi2[0]>chi2[1])//found minimum
      {
        g3 = new TGraph(3, val, chi2);
        fit = new TF1("pol2", "pol2");
        g3->Fit(fit, "Q");//fit to second order polynominal
        if(fit->GetParameter(2) > 0)//if the x^2 term is nonzero
        {
          best = -fit->GetParameter(1)/(2*fit->GetParameter(2));//the location of the minimum is -p1/(2*p2)
          minchi2 = fit->GetParameter(0) + fit->GetParameter(1)*best + fit->GetParameter(2)*best*best;
          cout<<"minimum with fit: "<<minchi2<<" at "<<best<<endl;
        }
        else//if the x^2 term is zero, then just use the minimum you got by scanning
        {
          best = val[1];
          minchi2 = chi2[1];
          cout<<"minimum with scan: "<<minchi2<<" at "<<best<<endl;
        }
        break;
      }
    }
    bestfit_invt.push_back(best);
    
    limit = -1.;
    delchi2prev = 1000;
    sprev = best;
    for(i=0;i<(epsetauHigh-best)/epsetauincrement;i++)
    {
      //epsetau = i*epsetauincrement + epsetauLow;
      epsetau = i*epsetauincrement + best;
      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetEps_etau(epsetau);
        ExtrapFHC[ie]->OscillatePrediction();
      }
      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetEps_etau(epsetau);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      delchi2 = 1e10;
      if(FitMethod==0)
      {
        delchi2 = PoissonChi2(NExp) - minchi2;
      }
      else if(FitMethod==1)
      {
        delchi2 = ScaledChi2(Bkgd,Sig) - minchi2;
      }
      else if(FitMethod==2)
      {
        delchi2 = StandardChi2(NExp) - minchi2;
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        delchi2 = StandardLikelihood() - minchi2;
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        delchi2 = BinLikelihood() - minchi2;
      }
      else
      {
        cout<<"Error in RunNSIDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
      
      if(i==1) continue;
      
      if(delchi2>contourlvl && delchi2prev<contourlvl)
      {
        limit = epsetau + ((epsetau-sprev)/(delchi2 - delchi2prev))*(contourlvl - delchi2);
        cout<<delchi2prev<<", "<<sprev<<", "<<delchi2<<", "<<epsetau<<", "<<limit<<endl;
        break;
      }
      delchi2prev = delchi2;
      sprev = epsetau;
    }
    limit_invt.push_back(limit);


    limit = -1.;
    delchi2prev = 1000;
    sprev = best;
    for(i=0;i<(best-epsetauLow)/epsetauincrement;i++)
    {
      //epsetau = i*epsetauincrement + epsetauLow;
      epsetau = -i*epsetauincrement + best;
      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ExtrapFHC[ie]->SetEps_etau(epsetau);
        ExtrapFHC[ie]->OscillatePrediction();
      }
      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ExtrapRHC[ie]->SetEps_etau(epsetau);
        ExtrapRHC[ie]->OscillatePrediction();
      }

      Bkgd->Reset();
      Sig->Reset();
      NExp->Reset();

      CombineFHCRHC();
      NExp->Add(Bkgd);
      NExp->Add(Sig);
      
      delchi2 = 1e10;
      if(FitMethod==0)
      {
        delchi2 = PoissonChi2(NExp) - minchi2;
      }
      else if(FitMethod==1)
      {
        delchi2 = ScaledChi2(Bkgd,Sig) - minchi2;
      }
      else if(FitMethod==2)
      {
        delchi2 = StandardChi2(NExp) - minchi2;
      }
      else if(FitMethod==3)
      {
        //Likelihood: "Standard" (N syst, N nuisance)
        //Calculate the likelihood (x2 for chi)
        delchi2 = StandardLikelihood() - minchi2;
      }
      else if(FitMethod==4)
      {
        //Likelihood: Bin by Bin Calculation of Systematics
        //Calculate the likelihood (x2 for chi)
        delchi2 = BinLikelihood() - minchi2;
      }
      else
      {
        cout<<"Error in RunNSIDeltaChi2Contour(): Unknown 'FitMethod'."<<endl;
      }
      
      if(i==1) continue;
      
      if(delchi2>contourlvl && delchi2prev<contourlvl)
      {
        limit = epsetau + ((epsetau-sprev)/(delchi2 - delchi2prev))*(contourlvl - delchi2);
        cout<<delchi2prev<<", "<<sprev<<", "<<delchi2<<", "<<epsetau<<", "<<limit<<endl;
        break;
      }
      delchi2prev = delchi2;
      sprev = epsetau;
    }
    limit_invtlow.push_back(limit);
  }
  
  double *s_best_n = new double[nDeltaSteps+1];
  double *s_best_i = new double[nDeltaSteps+1];
  double *s_limit_n = new double[nDeltaSteps+1];
  double *s_limit_i = new double[nDeltaSteps+1];
  double *s_limit_nlow = new double[nDeltaSteps+1];
  double *s_limit_ilow = new double[nDeltaSteps+1];
  double *d = new double[nDeltaSteps+1];
  for(i=0;i<nDeltaSteps+1;i++)
  {
    d[i] = deltapts.at(i);
    s_best_n[i] = bestfit_norm.at(i);
    s_best_i[i] = bestfit_invt.at(i);
    s_limit_n[i] = limit_norm.at(i);
    s_limit_i[i] = limit_invt.at(i);
    
    s_limit_nlow[i] = limit_normlow.at(i);
    s_limit_ilow[i] = limit_invtlow.at(i);
  }
  
  TGraph *gn = new TGraph(nDeltaSteps+1,s_best_n,d);
  gn->SetMarkerStyle(20);
  gn->SetTitle("");
  gn->GetYaxis()->SetTitle("#delta");
  gn->GetXaxis()->SetTitle("#epsilon_{e#tau}");
  gn->GetXaxis()->SetLimits(-3.0,3.0);
  gn->SetLineWidth(4);
  gn->SetMaximum(2);
  gn->SetName("BestFit_Normal");
  
  TGraph *gi = new TGraph(nDeltaSteps+1,s_best_i,d);
  gi->SetMarkerStyle(20);
  gi->SetTitle("");
  gi->GetYaxis()->SetTitle("#delta");
  gi->GetXaxis()->SetTitle("#epsilon_{e#tau}");
  gi->GetXaxis()->SetLimits(-3.0,3.0);
  gi->SetLineWidth(4);
  gi->SetLineStyle(2);
  gi->SetMaximum(2);
  gi->SetName("BestFit_Inverted");
  
  TGraph *gn_limit = new TGraph(nDeltaSteps+1,s_limit_n,d);
  gn_limit->SetMarkerStyle(20);
  gn_limit->SetTitle("");
  gn_limit->GetYaxis()->SetTitle("#delta");
  gn_limit->GetXaxis()->SetTitle("#epsilon_{e#tau}");
  gn_limit->GetXaxis()->SetLimits(-3.0,3.0);
  gn_limit->SetLineWidth(4);
  gn_limit->SetLineColor(kBlue);
  gn_limit->SetMarkerColor(kBlue);
  gn_limit->SetMaximum(2);
  gn_limit->SetName("Limit_Normal");

 TGraph *gn_limitlow = new TGraph(nDeltaSteps+1,s_limit_nlow,d);
  gn_limitlow->SetMarkerStyle(20);
  gn_limitlow->SetTitle("");
  gn_limitlow->GetYaxis()->SetTitle("#delta");
  gn_limitlow->GetXaxis()->SetTitle("#epsilon_{e#tau}");
  gn_limitlow->GetXaxis()->SetLimits(-3.0,3.0);
  gn_limitlow->SetLineWidth(4);
  gn_limitlow->SetLineColor(kBlue);
  gn_limitlow->SetMarkerColor(kBlue);
  gn_limitlow->SetMaximum(2);
  gn_limitlow->SetName("Limit_NormalLow");
  
  TGraph *gi_limit = new TGraph(nDeltaSteps+1,s_limit_i,d);
  gi_limit->SetMarkerStyle(20);
  gi_limit->SetTitle("");
  gi_limit->GetYaxis()->SetTitle("#delta");
  gi_limit->GetXaxis()->SetTitle("#epsilon_{e#tau}");
  gi_limit->GetXaxis()->SetLimits(-3.0,3.0);
  gi_limit->SetLineWidth(4);
  gi_limit->SetLineColor(kRed);
  gi_limit->SetMarkerColor(kRed);
  gi_limit->SetMaximum(2);
  gi_limit->SetName("Limit_Inverted");
 
 TGraph *gi_limitlow = new TGraph(nDeltaSteps+1,s_limit_ilow,d);
  gi_limitlow->SetMarkerStyle(20);
  gi_limitlow->SetTitle("");
  gi_limitlow->GetYaxis()->SetTitle("#delta");
  gi_limitlow->GetXaxis()->SetTitle("#epsilon_{e#tau}");
  gi_limitlow->GetXaxis()->SetLimits(-3.0,3.0);
  gi_limitlow->SetLineWidth(4);
  gi_limitlow->SetLineColor(kRed);
  gi_limitlow->SetMarkerColor(kRed);
  gi_limitlow->SetMaximum(2);
  gi_limitlow->SetName("Limit_InvertedLow");
 
  if(cl==0) cout<<"90% ";
  if(cl==1) cout<<"68% ";
  if(cl==2) cout<<"95% ";
  if(cl==3) cout<<"99% ";
  cout<<"confidence level limit = "<<limit_norm.at(0)<<", "<<limit_invt.at(0)<<endl;
  
  TFile *fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  gn->Write();
  gi->Write();
  gn_limit->Write();
  gi_limit->Write();
  gn_limitlow->Write();
  gi_limitlow->Write();
  fout->Write();
  fout->Close();
  
  delete [] s_best_n;
  delete [] s_best_i;
  delete [] s_limit_n;
  delete [] s_limit_i;
  delete [] s_limit_nlow;
  delete [] s_limit_ilow;
  delete [] d;
  
  return;
}






//////////////////// NSI Fitting functions //////////////////////////



////////////////////// End of NSI fitting functions ///////////////////////
