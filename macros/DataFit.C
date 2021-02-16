{
  
  TFile save("DataFit.root","RECREATE");
  
  Int_t nNuMuFiles;
  Int_t nNueFiles;
  Int_t nNuTauFiles;
  Int_t nNearFiles;
  Int_t nChallengeNearFiles;

  Double_t NuMuFilesPOT;
  Double_t NueFilesPOT;
  Double_t NuTauFilesPOT;
  Double_t NearFilesPOT;
  Double_t ChallengeNearFilesPOT;

  Double_t MDCChallengePOT = 7.4e20;  
  Double_t MDCNearToFar = 1;

  TFile *systematicFile;

  TFile *fmdcfile = NULL;
  TFile *fdatafile = NULL;

  if(true){   //cuts
    nNuMuFiles = 5;
    NuMuFilesPOT = 5.*6.5e20;
    nNueFiles = 5;
    NueFilesPOT = 5.*6.5e20;
    nNuTauFiles = 5;
    NuTauFilesPOT = 5.*6.5e20;
    nNearFiles = 76;
    NearFilesPOT = 76.*550*2.4e13;
    MDCNearToFar = 5.6093e-05;
    systematicFile = new TFile("sysHys/moneyplot_jbms.root","READ");
    nChallengeNearFiles = 246;
    ChallengeNearFilesPOT = 246.*550*2.4e13;
    fmdcfile = new TFile("SensitivityFile_MDC_cuts.root","READ");
    fdatafile = new TFile("SensitivityFile_Data_cuts.root","READ");
  }

  if(!fdatafile) {
    cout << "No data file. Quitting." << endl;
    return;
  }

  TH2F *systematicHist_allSys = (TH2F*) systematicFile->Get("hfdvnd");  
  TH2F *systematicHist_oscSys = (TH2F*) systematicFile->Get("hfdvnd");
  Double_t systematicHistNorm = 7.4e20; //Trish's normalisation

  Int_t nThetaPoints = 40;
  Double_t *theta13 = new Double_t[nThetaPoints];
  Int_t nDeltaPoints = 60;;
  Double_t *delta23 = new Double_t[nDeltaPoints];
  for(int i=0;i<nThetaPoints;i++){
    theta13[i] = 0.000 + 0.005*Double_t(i);
  }
  for(int i=0;i<nDeltaPoints;i++){
    delta23[i] = 0.001 + 0.0001*Double_t(i);
  }  

  TH2F *nhist = (TH2F*) fmdcfile.Get("sensitivity/nearHist");
  if(nhist->Integral()==0) nhist==NULL; 

  TH2F *fhist[40][60];
  for(int i=0;i<nThetaPoints;i++){
    for(int j=0;j<nDeltaPoints;j++){
      char name[256];
      if(fmdcfile){
	sprintf(name,"sensitivity/farHist%.0f_%.0f",
		1000.*theta13[i],10000.*delta23[j]);
	fhist[i][j] = (TH2F*) fmdcfile.Get(name);
      }
      else fhist[i][j] = NULL;
    }
  }

  TF1 *probgaus = new TF1("probgaus","gaus",-10,10);
  probgaus->SetParameters(1./TMath::Sqrt(2.*TMath::Pi()),0,1);

  TH2F *fardatahist = (TH2F*) fdatafile->Get("sensitivity/farHist0_25");
  TH1D *dataHist = fardatahist->ProjectionY();
  dataHist->Rebin(4);  
  TH2F *neardatahist = (TH2F*) fdatafile->Get("sensitivity/nearHist");
  if(neardatahist->Integral()==0) neardatahist = NULL;
  
  Double_t farnorm = 1;
  if(neardatahist && nhist) {
    farnorm = neardatahist->Integral()/nhist->Integral();
  }
 
  TH2F *Chi2Surf = new TH2F("Chi2Surf","Chi2Surf",
			    nThetaPoints,theta13[0]-0.0025,
			    theta13[nThetaPoints-1]+0.0025,
			    nDeltaPoints,delta23[0]-0.00005,
			    delta23[nDeltaPoints-1]+0.00005);
  
  Int_t bestdel = 0;
  Int_t bestthe = 0;
  Float_t bestChi2 = 9999;

  //analysis
  for(int i=0;i<nDeltaPoints;i++){
    for(int j=0;j<nThetaPoints;j++){
      TH2F *theFar = fhist[j][i];
      TH1D *mcHist = theFar->ProjectionY("mcHist",2,2);
      TH1D *bkgHist = theFar->ProjectionY("bkgHist");
      bkgHist->Add(mcHist,-1);
      mcHist->Rebin(4);
      bkgHist->Rebin(4);
      
      Float_t chi2 = 0;
      for(int k=1;k<=dataHist->GetNbinsX();k++){
	double data = dataHist->GetBinContent(k);
	double mc = mcHist->GetBinContent(k);
	double bkg = bkgHist->GetBinContent(k)*farnorm;
	if(data==0||mc+bkg==0) chi2 += 2.*(mc+bkg-data);
	else chi2 += 2.*(mc+bkg-data)+2.*data*TMath::Log(data/(mc+bkg));
      }
      Chi2Surf->Fill(theta13[j],delta23[i],chi2);
      if(chi2<bestChi2){
	bestdel = i;
	bestthe = j;
	bestChi2 = chi2;	
      }
      if(i==14) {
	if(j==0) {
	  cout << delta23[i] << " "  << theta13[j] << endl;
	  cout << "Data Events = " << dataHist->Integral() << endl;
	  cout << "Expected Background using simple ND Normalisation = " 
	       << bkgHist->Integral()*farnorm << endl;
	  Float_t mean = 0;
	  Float_t sigma = 0;
	  if(neardatahist){
	    Int_t bin = systematicHist_allSys->GetXaxis()
	      ->FindBin(neardatahist->Integral());
	    TH1D *tempHist = systematicHist_allSys
	      ->ProjectionY("tempHist",bin,bin);
	    mean = tempHist->GetMean();
	    sigma = tempHist->GetRMS();
	    cout << "Expected Mean from Sys Study = " << mean << endl;
	    cout << "Expected sigma from Sys Study = " << sigma << endl;
	  }
	  else {
	    mean = bkgHist->Integral()*farnorm;
	    sigma = 0.04 * mean; //assume 4% systematic
	  }
	  Float_t ts = ( (dataHist->Integral() - mean)/
			 sqrt(dataHist->Integral() + 
			      sigma*sigma));
	  cout << "Test Statistic = " << ts << endl;
	  cout << "Prob of this #events or greater given Null Hypothesis = "
	       << probgaus->Integral(ts,10.) << endl;
	}
      }
      delete mcHist;
      delete bkgHist;
    }
  }

  cout << bestChi2 << endl;
  cout << bestthe << " " << bestdel << endl;
  dataHist->SetName("Data");
  TH2F *theFar = fhist[bestthe][bestdel];
  TH1D *bestmc = theFar->ProjectionY("mcHist",2,2);
  bestmc->Rebin(4);
  bestmc->SetName("BestFit_Nue");
  TH1D *bestbkg = theFar->ProjectionY("bkgHist");
  bestbkg->Rebin(4);
  bestbkg->SetName("BestFit_Background");
  
  //end job
  save.cd();
  Chi2Surf->Write();
  dataHist->Write();
  bestmc->Write();
  bestbkg->Write();
  save.Close();

  systematicFile->Close();

  delete [] theta13;
  delete [] delta23;

}
