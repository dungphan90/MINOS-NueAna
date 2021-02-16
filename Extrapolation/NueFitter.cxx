
// THis will be the actual engine that handles a Full Extrapolation
#include "NueAna/Extrapolation/NueFitter.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <iostream>

using namespace std;

NueFitter::NueFitter()
{
   fOutFile = "DefaultOut.root";
   fMinTotalEvents = -10;
   fFitMethod = 0;
}

bool NueFitter::LoadFiles()
{
   for(unsigned int i = 0; i < fFileList.size(); i++)
   {
      TFile* file = new TFile(fFileList[i].c_str(), "READ");
      TTree* Pos;
      file->GetObject("PositionTree", Pos);

      double dm32, delta, th13, massH;
      char name[100];
 
      double nevent, onesigma, ninety, threesigma, signal;
   
      Pos->SetBranchAddress("Th13", &th13);
      Pos->SetBranchAddress("DeltaM32", &dm32);
      Pos->SetBranchAddress("DeltaCP",&delta);
      Pos->SetBranchAddress("Hierarchy",&massH);
      Pos->SetBranchAddress("ParamID", &name);
      Pos->SetBranchAddress("Total", &nevent);
      Pos->SetBranchAddress("SixtyEightCut", &onesigma);
      Pos->SetBranchAddress("NinetyCut", &ninety);
      Pos->SetBranchAddress("ThreeSigma", &threesigma);
      Pos->SetBranchAddress("NueCC", &signal);
 
      Point temp;
       
      for(int j = 0; j < Pos->GetEntries(); j++){ 
        // For each position in this particular file
        Pos->GetEntry(j);
        temp.cutVals.clear();
        temp.nExpected = nevent;
        temp.nSignal = signal;
        double thVal = TMath::Sin(2*th13);
        thVal *= thVal;        

        temp.th13 = thVal;

        for(unsigned int k = 0; k < fContourLevels.size(); k++){
           Chi2Cut tempChi; 
           tempChi.percent = fContourLevels[k]; tempChi.cutval = 0;
           if(fContourLevels[k] == 68){
             tempChi.percent = 68; tempChi.cutval = onesigma; 
           }
           if(fContourLevels[k] == 90){
             tempChi.percent = 90; tempChi.cutval = ninety;
           }
           if(fContourLevels[k] == 99.73){
             tempChi.percent = 99.73; tempChi.cutval = threesigma;
           }

           if(tempChi.cutval == 0){
            cout<<"Why am I here... "<<fContourLevels[k]<<endl;
            //Clever function that takes the deltaLog histogram
            // from the file and figures out the cut values 
            // for now its the Generator code, but not really 
            //efficient if multiple cuts here
            TH1D* hist = 0;
            
            TString hName = string(name) + "/DeltaLogLikely";
            file->GetObject(hName, hist);

            double total = hist->GetSum();
            double count = 0.0;
            for(int m = 0; m < hist->GetNbinsX(); m++){
              count += hist->GetBinContent(m);
              if(count/total > fContourLevels[k]/100.0){
                 tempChi.cutval = hist->GetBinLowEdge(m);
                 m = hist->GetNbinsX();
              }
            }
            delete hist;
           }

           temp.cutVals.push_back(tempChi);
        } //end of loop over contours 

        //store it and move on

        vector<Chi2Cut>(temp.cutVals).swap(temp.cutVals);

        if(massH > 0) 
           fDeltaM32MapNormal[dm32][delta].push_back(temp);
        else
           fDeltaM32MapInvert[dm32][delta].push_back(temp);

        if(temp.th13 < 0 || temp.th13> 1)
          cout<<"Pushing "<<dm32<<" "<<delta<<"  "<<temp.th13<<"  "<<massH<<endl;
//        cout<<temp.cutVals.capacity()<<"  "<<fDeltaM32MapNormal[dm32][delta].capacity()<<"  "<<endl;
        if(j%10000 == 0) cout<<"Loaded "<<j<<" points..."<<endl;
     }
     
     //Check this file against the others
     TTree* fileInfo;
     file->GetObject("GenerationTree", fileInfo);
                                                         
     const int NERR = 6;
     double error[NERR];
     double eventMin;
                                                                                
     fileInfo->SetBranchAddress("MinTotal", &eventMin);
     fileInfo->SetBranchAddress("NCErr", &error[0]);
     fileInfo->SetBranchAddress("NuMuCCErr",&error[1]);
     fileInfo->SetBranchAddress("NueCCErr",&error[2]);
     fileInfo->SetBranchAddress("NuTauCCErr", &error[3]); 
     fileInfo->SetBranchAddress("BNueCCErr", &error[4]);
     fileInfo->SetBranchAddress("TotalBGErr", &error[5]);

     fileInfo->GetEntry(0);
     if(fMinTotalEvents < 0){
       fMinTotalEvents = eventMin;
       for(int m = 0; m < NERR; m++) fErrors[m] = error[m];
     }else{
       bool good = true;
//       good = good && (fMinTotalEvents == eventMin);
       if(fMinTotalEvents > eventMin) fMinTotalEvents = eventMin;
       cout<<fMinTotalEvents<<"  "<<eventMin<<endl;
       for(int m = 0; m < NERR; m++) {
          if(fErrors[m] > 0) cout<<fErrors[m]<<"  "<<error[m]<<endl;
          good = good && (fErrors[m] == error[m]);
       }

       if(!good)
         cout<<"Mismatch between files - check the GenerationTrees!"<<endl;
     }
  
     delete fileInfo;
     delete Pos;
     delete file;
   }
   std::cout<<"Finished Loading all files"<<std::endl;

   return true;
}

bool NueFitter::InitializeFittingHistograms()
{
   double dmDV = fDeltaM32MapNormal.begin()->first;

   std::map<double, std::vector<Point> >::iterator delBeg = fDeltaM32MapNormal[dmDV].begin();
   std::map<double, std::vector<Point> >::iterator delEnd = fDeltaM32MapNormal[dmDV].end();

   vector<double> deltaPos;
   while(delBeg != delEnd){
     deltaPos.push_back(delBeg->first);
     delBeg++;
   }
   
   vector<double> sin22th13Pos;
   std::vector<Point>::iterator thBeg = fDeltaM32MapNormal[dmDV][0].begin();
   std::vector<Point>::iterator thEnd = fDeltaM32MapNormal[dmDV][0].end();
                                                                                
   while(thBeg != thEnd){
     sin22th13Pos.push_back(thBeg->th13);
     thBeg++;
   }

   //Loaded in all the positions - now convert to array 
   // with offsets
   sort(deltaPos.begin(), deltaPos.end());
   sort(sin22th13Pos.begin(), sin22th13Pos.end());

   int nDelta = deltaPos.size();
   int nTh13 = sin22th13Pos.size();

   if(nDelta <= 1){
      cout<<"Only zero/one delta value??"<<std::endl;
      return false;
   }
   if(nTh13 <= 1){
      cout<<"Only zero/one theta value??"<<std::endl;
      return false;
   }

   double* deltaArray = new double[nDelta+1];
   double* sinthArray = new double[nTh13+1];   

   double offset = 0.0;
   for(int i = 1; i < nDelta; i++){
     double prev = deltaPos[i-1];
     double curr = deltaPos[i];
     offset = (curr - prev)/2.0;
     if(i == 1) deltaArray[0] = prev - offset;
     deltaArray[i] = curr - offset;
   }
   deltaArray[nDelta] = deltaArray[nDelta-1] + 2*offset;
 
   offset = 0.0;
   for(int i = 1; i < nTh13; i++){
     double prev = sin22th13Pos[i-1];
     double curr = sin22th13Pos[i];
     cout<<prev<<"  "<<curr<<endl;
     offset = (curr - prev)/2.0;
     if(i == 1) sinthArray[0] = prev - offset;
     sinthArray[i] = curr - offset;
   }
   sinthArray[nTh13] = sinthArray[nTh13-1] + 2*offset;

   TString n1 = "Chi2HistN";
   TString n2 = "Chi2HistI";

   fFitPointChi2N = new TH2D(n1, n1, nTh13, sinthArray, nDelta, deltaArray);
   fFitPointChi2I = new TH2D(n2, n2, nTh13, sinthArray, nDelta, deltaArray);

   fFitPointChi2N->SetDirectory(0);
   fFitPointChi2I->SetDirectory(0);   
   
   for(unsigned int i = 0; i < fContourLevels.size(); i++)
   {
      char num[3];
      sprintf(num, "%d", i);
      TString one = "Chi2Hist" + string(num) + "N";
      TString two = "Chi2Hist" + string(num) + "I";

      TH2D* temp = new TH2D(one, one, nTh13, sinthArray, nDelta, deltaArray);
      fContourHistsN.push_back(temp);
      temp  = new TH2D(two, two, nTh13, sinthArray, nDelta, deltaArray);
      fContourHistsI.push_back(temp);

      fContourHistsN[i]->SetDirectory(0);
      fContourHistsI[i]->SetDirectory(0);
   }

   cout<<"Finished Initializing Histograms"<<endl;
   return true;
}

double NueFitter::CalculateDeltaLog(double nexp, double nobs)
{
  //I assume that if the observation has a value greater than the minimum number  // of events, than it can be predicted perfectly
                                                                                
  double chi2 = 2*(nexp - nobs + nobs*TMath::Log(nobs/nexp));
                                                                                
  if(nobs < fMinTotalEvents)
    chi2 -= 2*(fMinTotalEvents - nobs + nobs*TMath::Log(nobs/fMinTotalEvents));
                                                                                
  if(nexp+1e-4 < fMinTotalEvents){
    cout<<"The expected number is less then the minimum number of events - "
        <<" something is messed up, please check: min: "<<fMinTotalEvents<<" nexp: "<<nexp<<endl;
  }
     
  return chi2;
}

double NueFitter::CalculateDeltaLog(double nexp, double nobs, double nsig)
{ 
    //First calculate the best fit to the particular point mu given these values of b,s,k
    // These are the same formulae used in the analytic feldman cousins

    double b = nexp - nsig;
    double s = nsig;
    double k = 1.0;
    double errK = 0.0774118*0.0774118;
    double errBg = b*b*0.073739*0.073739;
    double n = nobs;

    double sol_A = 1 + errK/errBg*s*s;
    double sol_B = errBg + s*k - 2*errK*s*s*b/errBg - b + s*s*errK;
    double sol_C = -n*errBg + s*k*errBg - s*s*errK*b - s*k*b + s*s*b*b*errK/errBg;

    double rad = sol_B*sol_B-4*sol_A*sol_C;
    if(rad < 0) std::cout<<"negative radical - not sure about this...."<<std::endl;

   double betaHat = (-sol_B + TMath::Sqrt(rad))/(2*sol_A);
   double kHat = k - s*errK/errBg*(b-betaHat);

   nexp = betaHat + s*kHat;
   double chisq = 2*(nexp - n + n*TMath::Log(n/nexp))
                             + (betaHat-b)*(betaHat-b)/errBg + (kHat - k)*(kHat-k)/errK;

   double chiBF = 0.0;

  if(n > b){
          chiBF = 0.0;
     }else{
       double  bBest = 0.5*(b - errBg + TMath::Sqrt( (b-errBg)*(b - errBg) + 4*n*errBg));
        // then we have no signal at best fit point

        chiBF = 2*(bBest - n + n*TMath::Log(n/bBest)) + (b-bBest)*(b-bBest)/errBg;
    }

//   if(chiBF > chisq+1e-10)
//       cout<<n<<"  "<<chisq<<"  "<<chiBF<<"  "<<bBest<<"  "<<betaHat<<"  "<<nexp<<"  "<<kHat<<"  "<<k<<"  "<<errBg<<"  "<<errK<<"  "<<s<<endl;

    return (chisq - chiBF);
}

bool NueFitter::BuildContourMaps()
{
  //This fills the maps for the various chi2 levels

  double dmDV = fDeltaM32MapNormal.begin()->first;

  std::map<double, std::vector<Point> >::iterator delBeg = fDeltaM32MapNormal[dmDV].begin();
  std::map<double, std::vector<Point> >::iterator delEnd = fDeltaM32MapNormal[dmDV].end();

  while(delBeg != delEnd){     
    double delta = delBeg->first;
    for(unsigned int i = 0; i < delBeg->second.size(); i++){
       // at this point
       double ss2th = delBeg->second[i].th13;
  
       int nChi2 = delBeg->second[i].cutVals.size();
       for(int j = 0; j < nChi2; j++){
         double val = delBeg->second[i].cutVals[j].cutval;
         fContourHistsN[j]->Fill(ss2th, delta, val);
       }
    }    
    delBeg++;
  }

  dmDV = fDeltaM32MapInvert.begin()->first;
  delBeg = fDeltaM32MapInvert[dmDV].begin();
  delEnd = fDeltaM32MapInvert[dmDV].end();
                                                                                
  while(delBeg != delEnd){
    for(unsigned int i = 0; i < delBeg->second.size(); i++){
       // at this point
       double delta = delBeg->first;
       double ss2th = delBeg->second[i].th13;
                                                                                
       int nChi2 = delBeg->second[i].cutVals.size();
       for(int j = 0; j < nChi2; j++){
         double val = delBeg->second[i].cutVals[j].cutval;
         fContourHistsI[j]->Fill(ss2th, delta, val);
       }
    }
    delBeg++;
  }

  cout<<"Completed contour maps"<<endl;
  return true;
}
 
bool NueFitter::PerformFit(double input, string outName)
{
   PrepareFit();
   return FitInput(input, outName);
}

bool NueFitter::PrepareFit()
{
   LoadFiles();
   if(!InitializeFittingHistograms()) return false;
   BuildContourMaps();

   return true;
}

bool NueFitter::FitInput(double input, string outName)
{
   if(outName != "dummy.root") SetOutputFile(outName);

   BuildChi2Map(input);
   WriteToFile();

   return true;
}

bool NueFitter::BuildChi2Map(double val)
{
  double dmDV = fDeltaM32MapNormal.begin()->first;
                                                                                
  std::map<double, std::vector<Point> >::iterator delBeg = fDeltaM32MapNormal[dmDV].begin();
  std::map<double, std::vector<Point> >::iterator delEnd = fDeltaM32MapNormal[dmDV].end();
  fFitPointChi2N->Reset();
  fFitPointChi2I->Reset();
                                                                                
  while(delBeg != delEnd){
    for(unsigned int i = 0; i < delBeg->second.size(); i++){
       // at this point
       double delta = delBeg->first;
       double ss2th = delBeg->second[i].th13;
       double nexp  = delBeg->second[i].nExpected;
       double nsig =  delBeg->second[i].nSignal;

       double chi2 = 0;
       if(fFitMethod == 0) chi2 = CalculateDeltaLog(nexp, val);
       if(fFitMethod == 1) chi2 = CalculateDeltaLog(nexp, val, nsig);

       fFitPointChi2N->Fill(ss2th, delta, chi2);
    }
    delBeg++;
  }

  delBeg = fDeltaM32MapInvert[dmDV].begin();
  delEnd = fDeltaM32MapInvert[dmDV].end();
                                                                                
  while(delBeg != delEnd){
    for(unsigned int i = 0; i < delBeg->second.size(); i++){
       // at this point
       double delta = delBeg->first;
       double ss2th = delBeg->second[i].th13;
       double nexp  = delBeg->second[i].nExpected;
       double nsig =  delBeg->second[i].nSignal;

       double chi2 = 0;
       if(fFitMethod == 0) chi2 = CalculateDeltaLog(nexp, val);
       if(fFitMethod == 1) chi2 = CalculateDeltaLog(nexp, val, nsig);
                                                                                
       fFitPointChi2I->Fill(ss2th, delta, chi2);
    }
    delBeg++;
  }

  cout<<"Completed chi2 map"<<endl;
  return true;
}

/*
bool NueFitter::BuildChi2Map(double val,  
    std::map<double, std::map<double, std::vector<Point> > > & map, 
    TH2D* hist){
*/

         
void NueFitter::WriteToFile()
{
  TFile* file = 0;  

  file = new TFile(fOutFile.c_str(),"RECREATE");
  file->cd();

  fFitPointChi2N->Write(); 
  fFitPointChi2I->Write();
                                                                                
  for(unsigned int j = 0; j < fContourLevels.size(); j++){
    fContourHistsN[j]->Write();
    fContourHistsI[j]->Write();
  }
  file->Close();
}

