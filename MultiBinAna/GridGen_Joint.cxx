#define GridGenJoint_C

#include "NueAna/MultiBinAna/GridGen_Joint.h"

GridGen_Joint::GridGen_Joint()
{
  SetOutputFile();
  
  SetNDeltaSteps();
  SetNSinSq2Th13Steps();
  SetDeltaRange();
  SetSinSq2Th13Range();
  SetNormalHierarchy();
  
  SetTheta12();
  SetTheta23();
  SetAbsValDeltaMSq23();
  SetDeltaMSq12();
  
  SetNExperiments();
  
  FreezeTheta23(false);
  
  return;
}
GridGen_Joint::~GridGen_Joint()
{
}

void GridGen_Joint::AddExtrapFHC(Extrapolate2D* E)
{
  ExtrapFHC.push_back(E);
  return;
}
void GridGen_Joint::AddExtrapRHC(Extrapolate2D* E)
{
  ExtrapRHC.push_back(E);
  return;
}


void GridGen_Joint::RunMultiBinOscParErrs(string s)
{
  if((ExtrapRHC.size() + ExtrapFHC.size())==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return;
  }
  
  unsigned int ie;
  for(ie=0;ie<ExtrapFHC.size();ie++)
  {
    ExtrapFHC[ie]->GetPrediction();
  }
  for(ie=0;ie<ExtrapRHC.size();ie++)
  {
    ExtrapRHC[ie]->GetPrediction();
  }
  
  for(ie=0;ie<ExtrapFHC.size();ie++)
  {
    ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,Theta12);
    ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,Theta23);
    ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
    ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);
    if(!NormalHier) ExtrapFHC[ie]->InvertMassHierarchy();
    ExtrapFHC[ie]->OscillatePrediction();
  }
    for(ie=0;ie<ExtrapRHC.size();ie++)
  {
    ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,Theta12);
    ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,Theta23);
    ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
    ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);
    if(!NormalHier) ExtrapRHC[ie]->InvertMassHierarchy();
    ExtrapRHC[ie]->OscillatePrediction();
  }
  int i,k;
  unsigned int j;

  int nbins=0;
  int nbinsFHC=0;
  int nbinsRHC=0;
  if (ExtrapFHC.size()>0) nbinsFHC = ExtrapFHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ExtrapRHC.size()>0) nbinsRHC = ExtrapRHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  nbins = nbinsRHC + nbinsFHC;

  int nPIDFHC=0;
  int nPIDRHC=0;
 if (ExtrapFHC.size()>0) nPIDFHC = ExtrapFHC[0]->GetNPID();
 if (ExtrapRHC.size()>0) nPIDRHC = ExtrapRHC[0]->GetNPID();
  double delta,t13axis;
  vector<double> sig, bkgd;
  vector< vector<double> > nc,numucc,bnue,tau,nue;
  vector<double> oscparerr;
  vector<double> oscparerr_offdiag;
  int noff;
  
  for(i=0;i<nbins;i++)
  {
    sig.push_back(0);
    bkgd.push_back(0);
    oscparerr.push_back(0);
    nc.push_back( vector<double>() );
    numucc.push_back( vector<double>() );
    bnue.push_back( vector<double>() );
    tau.push_back( vector<double>() );
    nue.push_back( vector<double>() );
    for(k=0;k<nbins;k++)
    {
      if(k>i)
      {
        oscparerr_offdiag.push_back(0);
      }
    }
    if (i<nbinsFHC){
      for(j=0;j<ExtrapFHC.size();j++)
	{
	  nc[i].push_back(0);
	  numucc[i].push_back(0);
	  bnue[i].push_back(0);
	  tau[i].push_back(0);
	  nue[i].push_back(0);
	} 
    } else if (i>=nbinsFHC){
      for(j=0;j<ExtrapRHC.size();j++)
	{
	  nc[i].push_back(0);
	  numucc[i].push_back(0);
	  bnue[i].push_back(0);
	  tau[i].push_back(0);
	  nue[i].push_back(0);
	}
    }
  }
  
  vector<TTree*> ftree;
  vector< vector<TTree*> > ftree2;
  TTree *ttmp;
  noff=0;
  for(i=0;i<nbins;i++)
  {
    ttmp = new TTree(Form("Bin_%i",i),Form("Bin_%i",i));
    ttmp->Branch("Delta",&delta,"Delta/D");
    ttmp->Branch("Th13Axis",&t13axis,"Th13Axis/D");
    ttmp->Branch("Signal",&sig[i],"Signal/D");
    ttmp->Branch("Background",&bkgd[i],"Background/D");
    ttmp->Branch("DNExp_DOscPars",&oscparerr[i],"DNExp_DOscPars/D");
    for(k=0;k<nbins;k++)
    {
      if(k>i)
      {
        ttmp->Branch(Form("Bin_%i_Bin_%i",i,k),&oscparerr_offdiag[noff],Form("Bin_%i_Bin_%i/D",i,k));
	noff++;
      }
    }
    ftree.push_back(ttmp);
    
    ttmp->Reset();
    
    ftree2.push_back( vector<TTree*>() );
    
   if (i<nbinsFHC){
 
     for(j=0;j<ExtrapFHC.size();j++)
       {
	 ttmp = new TTree(Form("Bin_%i_Run_%i",i,j),Form("Bin_%i_Run_%i",i,j));
	 ttmp->Branch("Delta",&delta,"Delta/D");
	 ttmp->Branch("Th13Axis",&t13axis,"Th13Axis/D");
	 ttmp->Branch("Signal",&nue[i][j],"Signal/D");
	 ttmp->Branch("NC",&nc[i][j],"NC/D");
	 ttmp->Branch("NuMuCC",&numucc[i][j],"NuMuCC/D");
	 ttmp->Branch("BNueCC",&bnue[i][j],"BNueCC/D");
	 ttmp->Branch("NuTauCC",&tau[i][j],"NuTauCC/D");
	 
	 ftree2[i].push_back(ttmp);
	 
	 ttmp->Reset();
       }
   } else if (i>=nbinsFHC){

     for(j=0;j<ExtrapRHC.size();j++)
       {
	 ttmp = new TTree(Form("Bin_%i_Run_%i",i,j),Form("Bin_%i_Run_%i",i,j));
	 ttmp->Branch("Delta",&delta,"Delta/D");
	 ttmp->Branch("Th13Axis",&t13axis,"Th13Axis/D");
	 ttmp->Branch("Signal",&nue[i][j],"Signal/D");
	 ttmp->Branch("NC",&nc[i][j],"NC/D");
	 ttmp->Branch("NuMuCC",&numucc[i][j],"NuMuCC/D");
	 ttmp->Branch("BNueCC",&bnue[i][j],"BNueCC/D");
	 ttmp->Branch("NuTauCC",&tau[i][j],"NuTauCC/D");
	 
	 ftree2[i].push_back(ttmp);
	 
	 ttmp->Reset();
       }
   }
  }
  
  double delta_increment = 0;
  if(nDeltaSteps>0) delta_increment = (DeltaHigh - DeltaLow)/(nDeltaSteps);
  double ssq2th13_increment = 0;
  if(nSinSq2Th13Steps>0) ssq2th13_increment = (SinSq2Th13High - SinSq2Th13Low)/(nSinSq2Th13Steps);
  
  int id,is,l,u;
  int ip,ir;
  double theta23,theta12,dm21,dm32,ssq2th13;
  vector<double> nexp,nobs;
  vector<TH1D*> delnexphist;
  vector<TH1D*> delnexphist_offdiag;
  for(i=0;i<nbins;i++)
  {
    delnexphist.push_back(new TH1D(Form("delnexphist_%i",i),"",400,-1,1));
    nexp.push_back(0);
    nobs.push_back(0);
    for(k=0;k<nbins;k++)
    {
      if(k>i)
      {
	delnexphist_offdiag.push_back(new TH1D(Form("delnexphist_%i_%i",i,k),"",400,-1,1));
      }
    }
  }
  
  gRandom->SetSeed(0);
  
  TFile *f = new TFile(gSystem->ExpandPathName(s.c_str()),"RECREATE");//save the osc par err distributions to this file
  
  l=0;
  for(id=0;id<nDeltaSteps+1;id++)
  {
    delta = (id*delta_increment + DeltaLow)*TMath::Pi();
    
    for(is=0;is<nSinSq2Th13Steps+1;is++)
    {
      t13axis = (is*ssq2th13_increment + SinSq2Th13Low);
      ssq2th13 = t13axis/(2.*TMath::Sin(Theta23)*TMath::Sin(Theta23));
      
      //get nominal prediction
      for(ie=0;ie<ExtrapFHC.size();ie++)
	{
	  ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,Theta12);
	  ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,Theta23);
	  ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
	  ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);//note that DeltaMSq23 is always positive
	  if(!NormalHier) ExtrapFHC[ie]->InvertMassHierarchy();
	  ExtrapFHC[ie]->SetDeltaCP(delta);
	  ExtrapFHC[ie]->SetSinSq2Th13(ssq2th13);
	  ExtrapFHC[ie]->OscillatePrediction();
	}
      
      for(ie=0;ie<ExtrapRHC.size();ie++)
	{
	  ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,Theta12);
	  ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,Theta23);
	  ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
	  ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);//note that DeltaMSq23 is always positive
	  if(!NormalHier) ExtrapRHC[ie]->InvertMassHierarchy();
	  ExtrapRHC[ie]->SetDeltaCP(delta);
	  ExtrapRHC[ie]->SetSinSq2Th13(ssq2th13);
	  ExtrapRHC[ie]->OscillatePrediction();
	}

      for(i=0;i<nbins;i++)
      {
        sig[i]=0;
        bkgd[i]=0;

	if (i<nbinsFHC){

	  ir = int(i/nPIDFHC);
	  ip = i%nPIDFHC;

	  for(ie=0;ie<ExtrapFHC.size();ie++)
	    {
	      bkgd[i] += (ExtrapFHC[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1));
	      sig[i] += (ExtrapFHC[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1));
	      
	      nc[i][ie] = ExtrapFHC[ie]->Pred[Background::kNC]->GetBinContent(ip+1,ir+1);
	      numucc[i][ie] = ExtrapFHC[ie]->Pred[Background::kNuMuCC]->GetBinContent(ip+1,ir+1);
	      bnue[i][ie] = ExtrapFHC[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
	      tau[i][ie] = ExtrapFHC[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip+1,ir+1);
	      nue[i][ie] = ExtrapFHC[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1);
	    }
	} else if (i>=nbinsFHC){

	  ir = int((i-nbinsFHC)/nPIDRHC);
	  ip = (i-nbinsFHC)%nPIDRHC;

	  for(ie=0;ie<ExtrapRHC.size();ie++)
	    {
	      bkgd[i] += (ExtrapRHC[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1-nbinsFHC));
	      sig[i] += (ExtrapRHC[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1-nbinsFHC));
	      
	      nc[i][ie] = ExtrapRHC[ie]->Pred[Background::kNC]->GetBinContent(ip+1,ir+1);
	      numucc[i][ie] = ExtrapRHC[ie]->Pred[Background::kNuMuCC]->GetBinContent(ip+1,ir+1);
	      bnue[i][ie] = ExtrapRHC[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
	      tau[i][ie] = ExtrapRHC[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip+1,ir+1);
	      nue[i][ie] = ExtrapRHC[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1);
	    }
	}

        nexp[i] = sig[i]+bkgd[i];
      }
      
      //do pseudo experiments
      noff=0;
      for(i=0;i<nbins;i++)
      {
        delnexphist[i]->Reset();
        delnexphist[i]->SetName(Form("DeltaNexp_%i_%i_Diag_%i",id,is,i));
	for(k=0;k<nbins;k++)
	{
	  if(k>i)
	  {
	    delnexphist_offdiag[noff]->Reset();
	    delnexphist_offdiag[noff]->SetName(Form("DeltaNexp_%i_%i_OffDiag_%i_%i",id,is,i,k));
	    noff++;
	  }
	}
      }

      //Do this NumExpts times
      for(u=0;u<NumExpts;u++)
      {
        theta12 = Theta12 + AsymGaus(dTheta12_dn,dTheta12_up);
        dm21 = DeltaMSq12 + AsymGaus(dDeltaMSq12_dn,dDeltaMSq12_up);
        dm32 = DeltaMSq23 + AsymGaus(dDeltaMSq23_dn,dDeltaMSq23_up);
        //theta23 = DrawTheta23(dTheta23_up);
        theta23 = Theta23 + AsymGaus(dTheta23_dn,dTheta23_up);
        
        ssq2th13 = t13axis/(2.*TMath::Sin(theta23)*TMath::Sin(theta23));

	for(ie=0;ie<ExtrapFHC.size();ie++)
        {
          ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,theta12);
          ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,dm21);
          ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,dm32);
          if(!NormalHier) ExtrapFHC[ie]->InvertMassHierarchy();
          ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,theta23);
          ExtrapFHC[ie]->SetSinSq2Th13(ssq2th13);
          ExtrapFHC[ie]->OscillatePrediction();
        }
	for(ie=0;ie<ExtrapRHC.size();ie++)
        {
          ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,theta12);
          ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,dm21);
          ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,dm32);
          if(!NormalHier) ExtrapRHC[ie]->InvertMassHierarchy();
          ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,theta23);
          ExtrapRHC[ie]->SetSinSq2Th13(ssq2th13);
          ExtrapRHC[ie]->OscillatePrediction();
        }
	
	noff=0;
	for(i=0;i<nbins;i++)
	{
	  nobs[i]=0;
	  
	  if (i<nbinsFHC){
	    for(ie=0;ie<ExtrapFHC.size();ie++)
	      {
		nobs[i] += (ExtrapFHC[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1));
		nobs[i] += (ExtrapFHC[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1));
	      }
	  } else if (i>=nbinsFHC){
	    for(ie=0;ie<ExtrapRHC.size();ie++)
	      {
		nobs[i] += (ExtrapRHC[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1-nbinsFHC));
		nobs[i] += (ExtrapRHC[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1-nbinsFHC));
	      }
	  }
          delnexphist[i]->Fill((nobs[i]-nexp[i])/(nexp[i]));
	}

	for(i=0;i<nbins;i++)
	{
	  for(k=0;k<nbins;k++)
	  {
	    if(k>i)
	    {
	      delnexphist_offdiag[noff]->Fill((nobs[i]-nexp[i])*(nobs[k]-nexp[k])/(nexp[i]*nexp[k]));
	      noff++;
	    }
	  }
	}
      }
      
      noff=0;
      for(i=0;i<nbins;i++)
      {
	oscparerr[i] = delnexphist[i]->GetRMS();
	delnexphist[i]->Write();
	for(k=0;k<nbins;k++)
	{
	  if(k>i)
	  {
	    oscparerr_offdiag[noff] = delnexphist_offdiag[noff]->GetRMS();
            if(delnexphist_offdiag[noff]->GetMean()<0) oscparerr_offdiag[noff] = -1.*oscparerr_offdiag[noff];
	    delnexphist_offdiag[noff]->Write();
	    noff++;
	  }
	}
      }
      f->Close();
      
      gROOT->cd("/");
      
      for(i=0;i<nbins;i++)
      {
	if (i<nbinsFHC){
	  for(ie=0;ie<ExtrapFHC.size();ie++)
	    {
	      ftree2[i][ie]->Fill();
	    }
	} else if (i>=nbinsFHC){
	  for(ie=0;ie<ExtrapRHC.size();ie++)
	    {
	      ftree2[i][ie]->Fill();
	    }
	}
        ftree[i]->Fill();
      }
      
      f = new TFile(gSystem->ExpandPathName(s.c_str()),"UPDATE");
      
      if(l%100==0) cout<<100.*l/((nDeltaSteps+1)*(nSinSq2Th13Steps+1))<<"% complete"<<endl;
      l++;
    }
  }
  
  f->Close();
  
  double nPOTNearFHC,nPOTFarFHC;
  double nPOTNearRHC,nPOTFarRHC;
  
  TTree *paramtree = new TTree("paramtree","paramtree");
  paramtree->Branch("nearPOTFHC",&nPOTNearFHC,"nearPOTFHC/D");
  paramtree->Branch("farPOTFHC",&nPOTFarFHC,"farPOTFHC/D");
  paramtree->Branch("nearPOTRHC",&nPOTNearRHC,"nearPOTRHC/D");
  paramtree->Branch("farPOTRHC",&nPOTFarRHC,"farPOTRHC/D");
  paramtree->Branch("Theta12",&Theta12,"Theta12/D");
  paramtree->Branch("Theta23",&Theta23,"Theta23/D");
  paramtree->Branch("DeltaMSq23",&dm32,"DeltaMSq23/D");
  paramtree->Branch("DeltaMSq12",&DeltaMSq12,"DeltaMSq12/D");
  
  dm32 = DeltaMSq23;
  if(!NormalHier) dm32 = -1.*DeltaMSq23;
  
  for(ie=0;ie<ExtrapFHC.size();ie++)
  {

    nPOTNearRHC=0;
    nPOTFarRHC=0;
    nPOTNearFHC = ExtrapFHC[ie]->GetNearPOT();
    nPOTFarFHC = ExtrapFHC[ie]->GetFarPOT();
    paramtree->Fill();
  }
  for(ie=0;ie<ExtrapRHC.size();ie++)
  {
    nPOTNearFHC=0;
    nPOTFarFHC=0;
    nPOTNearRHC=0;
    nPOTFarRHC=0;
    nPOTNearRHC = ExtrapRHC[ie]->GetNearPOT();
    nPOTFarRHC = ExtrapRHC[ie]->GetFarPOT();
    paramtree->Fill();
  }
  
  TFile *fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  for(i=0;i<nbins;i++)
  {
    ftree[i]->Write();
    if (i<nbinsFHC){
      for(ie=0;ie<ExtrapFHC.size();ie++)
	{
	  ftree2[i][ie]->Write();
	} 
    } else if (i>=nbinsFHC){
      for(ie=0;ie<ExtrapRHC.size();ie++)
	{
	  ftree2[i][ie]->Write();
	} 
    }
  }
  paramtree->Write();
  fout->Close();
  
  return;
}
void GridGen_Joint::RunMultiBin_VaryTheta13(string s)
{
  if((ExtrapRHC.size() + ExtrapFHC.size())==0)
  {
    cout<<"No Extrapolate2D input.  Quitting..."<<endl;
    return;
  }
  
  unsigned int ie;
  for(ie=0;ie<ExtrapFHC.size();ie++)
  {
    ExtrapFHC[ie]->GetPrediction();
  }
  for(ie=0;ie<ExtrapRHC.size();ie++)
  {
    ExtrapRHC[ie]->GetPrediction();
  }
  
  for(ie=0;ie<ExtrapFHC.size();ie++)
  {
    ExtrapFHC[ie]->SetOscPar(OscPar::kTh13,Theta13);
    ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,Theta12);
    ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,Theta23);
    ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
    ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);
    if(!NormalHier) ExtrapFHC[ie]->InvertMassHierarchy();
    ExtrapFHC[ie]->OscillatePrediction();
  }
  for(ie=0;ie<ExtrapRHC.size();ie++)
  {
    ExtrapRHC[ie]->SetOscPar(OscPar::kTh13,Theta13);
    ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,Theta12);
    ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,Theta23);
    ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
    ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);
    if(!NormalHier) ExtrapRHC[ie]->InvertMassHierarchy();
    ExtrapRHC[ie]->OscillatePrediction();
  }
  int i,k;
  unsigned int j;
  
  int nbins=0;
  int nbinsFHC=0;
  int nbinsRHC=0;
  if (ExtrapFHC.size()>0) nbinsFHC = ExtrapFHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  if (ExtrapRHC.size()>0) nbinsRHC = ExtrapRHC[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  nbins = nbinsRHC + nbinsFHC;
  
  int nPIDFHC=0;
  int nPIDRHC=0;
  if (ExtrapFHC.size()>0) nPIDFHC = ExtrapFHC[0]->GetNPID();
  if (ExtrapRHC.size()>0) nPIDRHC = ExtrapRHC[0]->GetNPID();
  double delta,t13axis;
  vector<double> sig, bkgd;
  vector< vector<double> > nc,numucc,bnue,tau,nue;
  vector<double> oscparerr;
  vector<double> oscparerr_offdiag;
  int noff;
  
  t13axis = TMath::Sin(2*Theta13)*TMath::Sin(2*Theta13);
  
  for(i=0;i<nbins;i++)
  {
    sig.push_back(0);
    bkgd.push_back(0);
    oscparerr.push_back(0);
    nc.push_back( vector<double>() );
    numucc.push_back( vector<double>() );
    bnue.push_back( vector<double>() );
    tau.push_back( vector<double>() );
    nue.push_back( vector<double>() );
    for(k=0;k<nbins;k++)
    {
      if(k>i)
      {
        oscparerr_offdiag.push_back(0);
      }
    }
    if (i<nbinsFHC){
      for(j=0;j<ExtrapFHC.size();j++)
      {
        nc[i].push_back(0);
        numucc[i].push_back(0);
        bnue[i].push_back(0);
        tau[i].push_back(0);
        nue[i].push_back(0);
      } 
    } else if (i>=nbinsFHC){
      for(j=0;j<ExtrapRHC.size();j++)
      {
        nc[i].push_back(0);
        numucc[i].push_back(0);
        bnue[i].push_back(0);
        tau[i].push_back(0);
        nue[i].push_back(0);
      }
    }
  }
  
  vector<TTree*> ftree;
  vector< vector<TTree*> > ftree2;
  TTree *ttmp;
  noff=0;
  for(i=0;i<nbins;i++)
  {
    ttmp = new TTree(Form("Bin_%i",i),Form("Bin_%i",i));
    ttmp->Branch("Delta",&delta,"Delta/D");
    ttmp->Branch("Th13Axis",&t13axis,"Th13Axis/D");
    ttmp->Branch("Signal",&sig[i],"Signal/D");
    ttmp->Branch("Background",&bkgd[i],"Background/D");
    ttmp->Branch("DNExp_DOscPars",&oscparerr[i],"DNExp_DOscPars/D");
    for(k=0;k<nbins;k++)
    {
      if(k>i)
      {
        ttmp->Branch(Form("Bin_%i_Bin_%i",i,k),&oscparerr_offdiag[noff],Form("Bin_%i_Bin_%i/D",i,k));
        noff++;
      }
    }
    ftree.push_back(ttmp);
    
    ttmp->Reset();
    
    ftree2.push_back( vector<TTree*>() );
    
    if (i<nbinsFHC){
      
      for(j=0;j<ExtrapFHC.size();j++)
      {
        ttmp = new TTree(Form("Bin_%i_Run_%i",i,j),Form("Bin_%i_Run_%i",i,j));
        ttmp->Branch("Delta",&delta,"Delta/D");
        ttmp->Branch("Th13Axis",&t13axis,"Th13Axis/D");
        ttmp->Branch("Signal",&nue[i][j],"Signal/D");
        ttmp->Branch("NC",&nc[i][j],"NC/D");
        ttmp->Branch("NuMuCC",&numucc[i][j],"NuMuCC/D");
        ttmp->Branch("BNueCC",&bnue[i][j],"BNueCC/D");
        ttmp->Branch("NuTauCC",&tau[i][j],"NuTauCC/D");
        
        ftree2[i].push_back(ttmp);
        
        ttmp->Reset();
      }
    } else if (i>=nbinsFHC){
      
      for(j=0;j<ExtrapRHC.size();j++)
      {
        ttmp = new TTree(Form("Bin_%i_Run_%i",i,j),Form("Bin_%i_Run_%i",i,j));
        ttmp->Branch("Delta",&delta,"Delta/D");
        ttmp->Branch("Th13Axis",&t13axis,"Th13Axis/D");
        ttmp->Branch("Signal",&nue[i][j],"Signal/D");
        ttmp->Branch("NC",&nc[i][j],"NC/D");
        ttmp->Branch("NuMuCC",&numucc[i][j],"NuMuCC/D");
        ttmp->Branch("BNueCC",&bnue[i][j],"BNueCC/D");
        ttmp->Branch("NuTauCC",&tau[i][j],"NuTauCC/D");
        
        ftree2[i].push_back(ttmp);
        
        ttmp->Reset();
      }
    }
  }
  
  double delta_increment = 0;
  if(nDeltaSteps>0) delta_increment = (DeltaHigh - DeltaLow)/(nDeltaSteps);
  
  int id,l,u;
  int ip,ir;
  double theta23,theta12,dm21,dm32,theta13;
  vector<double> nexp,nobs;
  vector<TH1D*> delnexphist;
  vector<TH1D*> delnexphist_offdiag;
  for(i=0;i<nbins;i++)
  {
    delnexphist.push_back(new TH1D(Form("delnexphist_%i",i),"",400,-1,1));
    nexp.push_back(0);
    nobs.push_back(0);
    for(k=0;k<nbins;k++)
    {
      if(k>i)
      {
        delnexphist_offdiag.push_back(new TH1D(Form("delnexphist_%i_%i",i,k),"",400,-1,1));
      }
    }
  }
  
  gRandom->SetSeed(0);
  
  TFile *f = new TFile(gSystem->ExpandPathName(s.c_str()),"RECREATE");//save the osc par err distributions to this file
  
  l=0;
  for(id=0;id<nDeltaSteps+1;id++)
  {
    delta = (id*delta_increment + DeltaLow)*TMath::Pi();
    
    
    //get nominal prediction
    for(ie=0;ie<ExtrapFHC.size();ie++)
    {
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh13,Theta13);
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,Theta12);
      ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,Theta23);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
      ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);//note that DeltaMSq23 is always positive
      if(!NormalHier) ExtrapFHC[ie]->InvertMassHierarchy();
      ExtrapFHC[ie]->SetDeltaCP(delta);
      ExtrapFHC[ie]->OscillatePrediction();
    }
    
    for(ie=0;ie<ExtrapRHC.size();ie++)
    {
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh13,Theta13);
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,Theta12);
      ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,Theta23);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,DeltaMSq12);
      ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,DeltaMSq23);//note that DeltaMSq23 is always positive
      if(!NormalHier) ExtrapRHC[ie]->InvertMassHierarchy();
      ExtrapRHC[ie]->SetDeltaCP(delta);
      ExtrapRHC[ie]->OscillatePrediction();
    }
    
    for(i=0;i<nbins;i++)
    {
      sig[i]=0;
      bkgd[i]=0;
      
      if (i<nbinsFHC){
        
        ir = int(i/nPIDFHC);
        ip = i%nPIDFHC;
        
        for(ie=0;ie<ExtrapFHC.size();ie++)
        {
          bkgd[i] += (ExtrapFHC[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1));
          sig[i] += (ExtrapFHC[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1));
          
          nc[i][ie] = ExtrapFHC[ie]->Pred[Background::kNC]->GetBinContent(ip+1,ir+1);
          numucc[i][ie] = ExtrapFHC[ie]->Pred[Background::kNuMuCC]->GetBinContent(ip+1,ir+1);
          bnue[i][ie] = ExtrapFHC[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
          tau[i][ie] = ExtrapFHC[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip+1,ir+1);
          nue[i][ie] = ExtrapFHC[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1);
        }
      } else if (i>=nbinsFHC){
        
        ir = int((i-nbinsFHC)/nPIDRHC);
        ip = (i-nbinsFHC)%nPIDRHC;
        
        for(ie=0;ie<ExtrapRHC.size();ie++)
        {
          bkgd[i] += (ExtrapRHC[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1-nbinsFHC));
          sig[i] += (ExtrapRHC[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1-nbinsFHC));
          
          nc[i][ie] = ExtrapRHC[ie]->Pred[Background::kNC]->GetBinContent(ip+1,ir+1);
          numucc[i][ie] = ExtrapRHC[ie]->Pred[Background::kNuMuCC]->GetBinContent(ip+1,ir+1);
          bnue[i][ie] = ExtrapRHC[ie]->Pred[Background::kBNueCC]->GetBinContent(ip+1,ir+1);
          tau[i][ie] = ExtrapRHC[ie]->Pred[Background::kNuTauCC]->GetBinContent(ip+1,ir+1);
          nue[i][ie] = ExtrapRHC[ie]->Pred[Background::kNueCC]->GetBinContent(ip+1,ir+1);
        }
      }
      
      nexp[i] = sig[i]+bkgd[i];
     }
      
      //do pseudo experiments
      noff=0;
      for(i=0;i<nbins;i++)
      {
        delnexphist[i]->Reset();
        delnexphist[i]->SetName(Form("DeltaNexp_%i_Diag_%i",id,i));
        for(k=0;k<nbins;k++)
        {
          if(k>i)
          {
            delnexphist_offdiag[noff]->Reset();
            delnexphist_offdiag[noff]->SetName(Form("DeltaNexp_%i_OffDiag_%i_%i",id,i,k));
            noff++;
          }
        }
      }
      
      //Do this NumExpts times
      for(u=0;u<NumExpts;u++)
      {
        theta13 = Theta13 + AsymGaus(dTheta13_dn,dTheta13_up);
        theta12 = Theta12 + AsymGaus(dTheta12_dn,dTheta12_up);
        dm21 = DeltaMSq12 + AsymGaus(dDeltaMSq12_dn,dDeltaMSq12_up);
        dm32 = DeltaMSq23 + AsymGaus(dDeltaMSq23_dn,dDeltaMSq23_up);
        theta23 = Theta23 + AsymGaus(dTheta23_dn,dTheta23_up);
        
        for(ie=0;ie<ExtrapFHC.size();ie++)
        {
          ExtrapFHC[ie]->SetOscPar(OscPar::kTh13,theta13);
          ExtrapFHC[ie]->SetOscPar(OscPar::kTh12,theta12);
          ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM12,dm21);
          ExtrapFHC[ie]->SetOscPar(OscPar::kDeltaM23,dm32);
          if(!NormalHier) ExtrapFHC[ie]->InvertMassHierarchy();
          ExtrapFHC[ie]->SetOscPar(OscPar::kTh23,theta23);
          ExtrapFHC[ie]->OscillatePrediction();
        }
        for(ie=0;ie<ExtrapRHC.size();ie++)
        {
          ExtrapRHC[ie]->SetOscPar(OscPar::kTh13,theta13);
          ExtrapRHC[ie]->SetOscPar(OscPar::kTh12,theta12);
          ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM12,dm21);
          ExtrapRHC[ie]->SetOscPar(OscPar::kDeltaM23,dm32);
          if(!NormalHier) ExtrapRHC[ie]->InvertMassHierarchy();
          ExtrapRHC[ie]->SetOscPar(OscPar::kTh23,theta23);
          ExtrapRHC[ie]->OscillatePrediction();
        }
        
        noff=0;
        for(i=0;i<nbins;i++)
        {
          nobs[i]=0;
          
          if (i<nbinsFHC){
            for(ie=0;ie<ExtrapFHC.size();ie++)
            {
              nobs[i] += (ExtrapFHC[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1));
              nobs[i] += (ExtrapFHC[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1));
            }
          } else if (i>=nbinsFHC){
            for(ie=0;ie<ExtrapRHC.size();ie++)
            {
              nobs[i] += (ExtrapRHC[ie]->Pred_TotalBkgd_VsBinNumber->GetBinContent(i+1-nbinsFHC));
              nobs[i] += (ExtrapRHC[ie]->Pred_Signal_VsBinNumber->GetBinContent(i+1-nbinsFHC));
            }
          }
          delnexphist[i]->Fill((nobs[i]-nexp[i])/(nexp[i]));
        }
        
        for(i=0;i<nbins;i++)
        {
          for(k=0;k<nbins;k++)
          {
            if(k>i)
            {
              delnexphist_offdiag[noff]->Fill((nobs[i]-nexp[i])*(nobs[k]-nexp[k])/(nexp[i]*nexp[k]));
              noff++;
            }
          }
        }
      }
      
      noff=0;
      for(i=0;i<nbins;i++)
      {
        oscparerr[i] = delnexphist[i]->GetRMS();
        delnexphist[i]->Write();
        for(k=0;k<nbins;k++)
        {
          if(k>i)
          {
            oscparerr_offdiag[noff] = delnexphist_offdiag[noff]->GetRMS();
            if(delnexphist_offdiag[noff]->GetMean()<0) oscparerr_offdiag[noff] = -1.*oscparerr_offdiag[noff];
            delnexphist_offdiag[noff]->Write();
            noff++;
          }
        }
      }
      f->Close();
      
      gROOT->cd("/");
      
      for(i=0;i<nbins;i++)
      {
        if (i<nbinsFHC){
          for(ie=0;ie<ExtrapFHC.size();ie++)
          {
            ftree2[i][ie]->Fill();
          }
        } else if (i>=nbinsFHC){
          for(ie=0;ie<ExtrapRHC.size();ie++)
          {
            ftree2[i][ie]->Fill();
          }
        }
        ftree[i]->Fill();
      }
      
      f = new TFile(gSystem->ExpandPathName(s.c_str()),"UPDATE");
      
      if(l%100==0) cout<<100.*l/(nDeltaSteps+1)<<"% complete"<<endl;
      l++;
    
  }
  
  f->Close();
  
  double nPOTNearFHC,nPOTFarFHC;
  double nPOTNearRHC,nPOTFarRHC;
  
  TTree *paramtree = new TTree("paramtree","paramtree");
  paramtree->Branch("nearPOTFHC",&nPOTNearFHC,"nearPOTFHC/D");
  paramtree->Branch("farPOTFHC",&nPOTFarFHC,"farPOTFHC/D");
  paramtree->Branch("nearPOTRHC",&nPOTNearRHC,"nearPOTRHC/D");
  paramtree->Branch("farPOTRHC",&nPOTFarRHC,"farPOTRHC/D");
  paramtree->Branch("Theta13",&Theta13,"Theta13/D");
  paramtree->Branch("Theta12",&Theta12,"Theta12/D");
  paramtree->Branch("Theta23",&Theta23,"Theta23/D");
  paramtree->Branch("DeltaMSq23",&dm32,"DeltaMSq23/D");
  paramtree->Branch("DeltaMSq12",&DeltaMSq12,"DeltaMSq12/D");
  
  dm32 = DeltaMSq23;
  if(!NormalHier) dm32 = -1.*DeltaMSq23;
  
  for(ie=0;ie<ExtrapFHC.size();ie++)
  {
    
    nPOTNearRHC=0;
    nPOTFarRHC=0;
    nPOTNearFHC = ExtrapFHC[ie]->GetNearPOT();
    nPOTFarFHC = ExtrapFHC[ie]->GetFarPOT();
    paramtree->Fill();
  }
  for(ie=0;ie<ExtrapRHC.size();ie++)
  {
    nPOTNearFHC=0;
    nPOTFarFHC=0;
    nPOTNearRHC=0;
    nPOTFarRHC=0;
    nPOTNearRHC = ExtrapRHC[ie]->GetNearPOT();
    nPOTFarRHC = ExtrapRHC[ie]->GetFarPOT();
    paramtree->Fill();
  }
  
  TFile *fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()),"RECREATE");
  for(i=0;i<nbins;i++)
  {
    ftree[i]->Write();
    if (i<nbinsFHC){
      for(ie=0;ie<ExtrapFHC.size();ie++)
      {
        ftree2[i][ie]->Write();
      } 
    } else if (i>=nbinsFHC){
      for(ie=0;ie<ExtrapRHC.size();ie++)
      {
        ftree2[i][ie]->Write();
      } 
    }
  }
  paramtree->Write();
  fout->Close();
  
  return;
}
