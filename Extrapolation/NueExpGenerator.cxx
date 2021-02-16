#include "NueExpGenerator.h"
#include "NueAna/NueStandard.h"
#include "TRandom3.h"
#include "NueGenConfig.h"
#include <cmath>

NueExpGenerator::NueExpGenerator(){ 
  outFile = "default.root";
  fSeed = -10;
  fNumberOfExp = 10000;
  fScale = 1.0;
  fMinTotal = 0.0;

  fOffSet = -1;
  fReadUntil = -1;
  fMethod = 0;
  fFitMethod = 0;
  fObservation = 35;
}

void NueExpGenerator::Run(string errinput, string data, string output)
{
  SetOutputFile(output);
  NueGenConfig* gc = new NueGenConfig(errinput);
  gc->AddInputFile(data);
  if(!gc->CheckConfig()) return;

  if(fOffSet >= 0) gc->SetOffset(fOffSet);
  if(fReadUntil > 0) gc->SetReadNumber(fReadUntil);

  double* temp = gc->GetErrors();
  for(int i = 0; i < 6; i++) fErrors[i] = temp[i];

  if(fSeed > -1) fRand.SetSeed((int) fSeed);

  int count = 0;

  while(gc->LoadNextNumberSet(fNumbers)){
     if(fScale != 1){
         ReScale(fNumbers);
         if(count == 0)   fMinTotal *= fScale;
     }
     double *temp = gc->GetOscPar();
     for(int i = 0; i < 4; i++) fOscPar[i] = temp[i];
     Position* pos = GenerateExperimentSet(fNumbers, fErrors);
     DetermineCuts(pos);
     WriteToFile(pos); 
     CleanPos(pos);
     delete pos;
     int temp2 = 1;
     if(fNumberOfExp < 1000*10000) temp2 = int(1000*10000.0/fNumberOfExp);
     if(count%(temp2) == 0) cout<<count<<" positions generated."<<endl;
     count++;
  }

  WriteToFile(0);
}

void NueExpGenerator::CleanPos(Position* pos)
{
   delete pos->nevent;
   delete pos->deltaLog;
   delete pos->fDirectory;
//   delete pos;
}

bool NueExpGenerator::DetermineCuts(Position* pos)
{
  TH1D* hist = pos->deltaLog;

  double total = hist->GetSum();
  double count = 0.0;
  bool found68 = false;
  bool found90 = false;
  bool found3s = false;

  for(int i = 0; i < hist->GetNbinsX(); i++){
    count += hist->GetBinContent(i);
    if(count/total > 0.68 && !found68){ found68 = true; pos->sixtyeight = hist->GetBinCenter(i);}
    if(count/total > 0.90 && !found90){ found90 = true; pos->ninety = hist->GetBinCenter(i);}
    if(count/total > 0.9973 && !found3s){ found3s = true; pos->threesigma = hist->GetBinCenter(i);}

    if(count/total > 0.9973) break;
  }

  return found90;
}

void NueExpGenerator::ReScale(double *num)
{
  for(int i = 0; i < 5; i++){
     num[i] = fScale*num[i];
  }
}

Position* NueExpGenerator::GenerateExperimentSet(double* num, double* err)
{
   float nominal = 0.0;
//   int nevent = 0;
   double totBG = 0.0;

   Position* pos = new Position;
      
   for(int i = 0; i < 5; i++){
      pos->fEventNo[i] = num[i];      
      nominal += num[i];
      if(i != 2) totBG += num[i];
   }
   fNumbers[5] = nominal;

   char ID[100];
   int one  = (int) rint(fOscPar[0]*1e3);
   int two  =  (int) rint(1e2*fOscPar[1]);
   int three = (int) rint(1e4*fOscPar[2]);
   if(fOscPar[3] > 0) 
    sprintf(ID, "Pos_%03d_%03d_%02d_N", one, two, three);
   else
    sprintf(ID, "Pos_%03d_%03d_%02d_I", one, two, three);

//   cout<<int(fOscPar[0]*1e3)<<"  "<<rint(1e2*fOscPar[1])<<"  "<<rint(1e4*fOscPar[2])<<"  "
//       <<fOscPar[0]*1e3<<"  "<<1e2*fOscPar[1]<<"  "<<1e4*fOscPar[2]<<endl;
//   cout<<"Name: "<<ID<<endl;
   pos->fDirectory = new TDirectory(ID, ID);
   pos->fDirectory->cd();

   TString n1 = "nevent";
   TString n2 = "DeltaLogLikely";

   pos->id = string(ID);
   int max = 199;
   if(nominal > 60) max = int(3*nominal);
   pos->nevent = new TH1D(n1, n1, max+1, -1, max);
   pos->deltaLog = new TH1D(n2, n2, 50000, 0, 1000);

   TString n3[6];
   n3[0] = n1+"nc";  n3[1] = n1+"cc";  n3[2] = n1+"nue";
   n3[3] = n1+"tau"; n3[4] = n1+"bnue"; n3[5] = n1 + "tot";

//   for(int i = 0; i < 6; i++){
//     pos->sysevents[i] = new TH1D(n3[i], n3[i], max+1, -1, max);
//     statevents[i] = new TH1D(n1, n1, max+1, -1, max);
//   }


   gDirectory->cd("/");


   //fMethod 0, 1, 2
   // 0  shift the observation but not expectation
   // 1  shift the observation and expectation together
   // 2  Poissson(obs) and Gauss(exp) 

//   cout<<err[5]<<"  "<<totBG<<"  "<<err[2]<<"   "<<num[2]<<endl; 
   for(int k = 0; k < fNumberOfExp; k++)
   {
     double shift = fRand.Gaus(0, err[5]);
     if(shift < -1){       cout<<"Wow "<<shift<<" "<<err[5]<<endl; shift = -1; }

     double shiftedBg = totBG*(1+shift); 
     
     double bgSeed = shiftedBg;
     if(fMethod >= 2) bgSeed = totBG;

     int bgObs =  fRand.Poisson(bgSeed);

     shift = fRand.Gaus(0, err[ClassType::nue]);
     if(shift < -1) {     cout<<"Wow "<<shift<<" "<<err[ClassType::nue]<<endl; shift = -1; }
 
     double shiftedSig = num[ClassType::nue]*(1+shift);
     
     double sigSeed = shiftedSig;
     if(fMethod >= 2) sigSeed = num[ClassType::nue];
  
     int sigObs = fRand.Poisson(sigSeed);

     int nObs = sigObs + bgObs;
     
     double nexp = nominal;
     double bexp, sexp;
     bexp = sexp = 0.0;
     if(fMethod == 1 || fMethod >= 2){
          nexp = shiftedBg + shiftedSig;
          bexp = shiftedBg;
          sexp = shiftedSig;
     }

//     cout<<fMethod<<"  "<<nexp<<"  "<<nominal<<endl;

     pos->nevent->Fill(nObs);
     if(k%5000000 == 0 && k > 0) cout<<k<<"/"<<fNumberOfExp<<endl;

     double dll = 0;
 
     // Next we implement the special condition for n == nObs
     if(nObs == fObservation){
         bexp = totBG;
         sexp = num[2];
         nexp = bexp + sexp;
     }   

     if(fFitMethod == 0) dll = CalculateDeltaLog(nexp, nObs, bexp);
     else{     
        double errBg = err[5]*err[5]*bexp*bexp;;
        double errK = err[2]*err[2]*sexp*sexp/(num[2]*num[2]);
       dll = CalculateChi2withBestFit(bexp, num[2], sexp/num[2], errBg,  errK, nObs);
     }
//     CalculateDeltaLog(nexp, nObs, shiftedBg);

     pos->deltaLog->Fill(dll);
   }

//   cout<<"I'm going home"<<endl;
   return pos;
}

double NueExpGenerator::CalculateDeltaLog(double nexp, double nobs, double b)
{   
  //I assume that if the observation has a value greater than the minimum number
  // of events, than it can be predicted perfectly
  double bound = fMinTotal;
  if(b > -0.5) bound = b; 

  double chi2 = 0;
  
  if(nobs < 1e-2){
     chi2 = 2*nexp;  //protection against nobs = 0
     if(nobs < bound)
     chi2 -= 2*(bound);
  }
  else{
     chi2 = 2*(nexp - nobs + nobs*TMath::Log(nobs/nexp));
  
     if(nobs < bound)
       chi2 -= 2*(bound - nobs + nobs*TMath::Log(nobs/bound));
   }

//   if(chi2 < -1e-10) cout<<"really: "<< 2*(nexp - nobs + nobs*TMath::Log(nobs/nexp))<<"  "<<nexp<<" "<<fMinTotal<<"  "<<chi2<<endl;

  if(nexp+1e-4 < fMinTotal){
//    cout<<"The expected number is less then the minimum number of events - "
//        <<" something is messed up, please check: min: "<<fMinTotal<<" nexp: "<<nexp<<endl;
  }
  
  return chi2;
}

double NueExpGenerator::CalculateChi2withBestFit(double b, double s, double k, double errBg, double errK, int n)
{
    //First calculate the best fit to the particular point mu given these values of b,s,k
    // These are the same formulae used in the analytic feldman cousins

    double sol_A = 1 + errK/errBg*s*s;
    double sol_B = errBg + s*k - 2*errK*s*s*b/errBg - b + s*s*errK; 
    double sol_C = -n*errBg + s*k*errBg - s*s*errK*b - s*k*b + s*s*b*b*errK/errBg;

    double rad = sol_B*sol_B-4*sol_A*sol_C;
    if(rad < 0) std::cout<<"negative radical - not sure about this...."<<std::endl;

   double betaHat = (-sol_B + TMath::Sqrt(rad))/(2*sol_A);
   double kHat = k - s*errK/errBg*(b-betaHat);

   double nexp = betaHat + s*kHat;
   double chisq = 2*(nexp - n + n*TMath::Log(n/nexp)) 
                             + (betaHat-b)*(betaHat-b)/errBg + (kHat - k)*(kHat-k)/errK;

   double bBest = 0.0;
     // Working out  chi^2 BF 
     double chiBF = 0.0;
    
     if(n > b){
          chiBF = 0.0;
     }else{
       double  bBest = 0.5*(b - errBg + TMath::Sqrt( (b-errBg)*(b - errBg) + 4*n*errBg));
        // then we have no signal at best fit point

        chiBF = 2*(bBest - n + n*TMath::Log(n/bBest)) + (b-bBest)*(b-bBest)/errBg;
    }

   if(chiBF > chisq-1e-10 && chiBF > 0) 
       cout<<n<<"  "<<chisq<<"  "<<chiBF<<"  "<<bBest<<"  "<<betaHat<<"  "<<nexp<<"  "<<kHat<<"  "<<k<<"  "<<errBg<<"  "<<errK<<"  "<<s<<endl;
           
    return (chisq - chiBF);
}


void NueExpGenerator::WriteToFile(Position *pos)
{
  static TFile* file = 0;
  static TTree* tree = 0;
  static char name[256];
  static int count = 0;
  static double sixty;
  static double ninety;
  static double threesig;
                                                                                
  if(file == 0){
     file = new TFile(outFile.c_str(),"RECREATE");
     tree = new TTree("PositionTree","PositionTree");
     tree->Branch("NC",&fNumbers[0],"NC/D");
     tree->Branch("NuMuCC",&fNumbers[1],"NuMuCC/D");
     tree->Branch("NueCC",&fNumbers[2], "NueCC/D");
     tree->Branch("NuTauCC",&fNumbers[3],"NuTauCC/D");
     tree->Branch("BNueCC",&fNumbers[4],"BNueCC/D");
     tree->Branch("Total", &fNumbers[5], "Total/D");

//     tree->Branch("NCErr",&fErrors[0],"NCErr/D");
//     tree->Branch("NuMuCCErr",&fErrors[1],"NuMuCCErr/D");
//     tree->Branch("NueCCErr",&fErrors[2], "NueCCErr/D");
//     tree->Branch("NuTauCCErr",&fErrors[3],"NuTauCCErr/D");
//     tree->Branch("BNueCCErr",&fErrors[4],"BNueCCErr/D");

     tree->Branch("Th13",&fOscPar[0],"Th13/D");
     tree->Branch("DeltaCP",&fOscPar[1],"DeltaCP/D");
     tree->Branch("DeltaM32",&fOscPar[2], "DeltaM32/D");
     tree->Branch("Hierarchy",&fOscPar[3],"Hierarchy/D");
     tree->Branch("ParamID", name, "ParamID/C");

     tree->Branch("SixtyEightCut", &sixty, "SixtyEightCut/D");
     tree->Branch("NinetyCut", &ninety, "NinetyCut/D");
     tree->Branch("ThreeSigma",&threesig, "ThreeSigma/D");
//     tree->SetAutoDelete(kFALSE);

  }
  file->cd();
 
  if(pos != 0){   
     string fnh_name = pos->id;
     sprintf(name, "%s", fnh_name.c_str());
     TDirectory *filedir = file->mkdir(fnh_name.c_str());
     filedir->cd();
     pos->nevent->Write();
     pos->deltaLog->Write();     
//     for(int i = 0; i < 6; i++) pos->sysevents[i]->Write();
     file->cd();
     sixty = pos->sixtyeight;
     ninety = pos->ninety;
     threesig = pos->threesigma;
//     cout<<sixty<<"  "<<ninety<<endl;
     tree->Fill();
     count++;
  }
  else{ 
    cout<<"Writing the trees"<<endl;
    tree->Write();

    TTree* tree2 = 0;
    tree2 = new TTree("GenerationTree","GenerationTree");
    tree2->Branch("MinTotal",&fMinTotal,"MinTotal/D");
    tree2->Branch("NumExperiment",&fNumberOfExp,"NumExperiment/I");
    tree2->Branch("NCErr",&fErrors[0],"NCErr/D");
    tree2->Branch("NuMuCCErr",&fErrors[1],"NuMuCCErr/D");
    tree2->Branch("NueCCErr",&fErrors[2], "NueCCErr/D");
    tree2->Branch("NuTauCCErr",&fErrors[3],"NuTauCCErr/D");
    tree2->Branch("BNueCCErr",&fErrors[4],"BNueCCErr/D");
    tree2->Branch("TotalBGErr",&fErrors[5],"TotalBGErr/D");
    tree2->Fill();

    tree2->Write(); 
    file->Close();
  }
}   
                                                                                    

