#include <string>
#include "NueAna/NueAnaTools/NueConvention.h"
#include "NueSensitivity.h"
#include "TMath.h"
#include "TChain.h"
#include "NueAna/NuePOT.h"
#include <iostream>
#include <fstream>

NueSensitivity::NueSensitivity(){
  t13Step = 0.001;
  dStep = 0.02;
  dm23Step = 0.025;
  fBaseLine = 735;
  fObservedIsSet = false;
}

void NueSensitivity::Run(std::string input, std::string output, double outPOT)
{
   //Loading in the data
   nsc = new NueSenseConfig(input);
   if(!nsc->CheckConfig()) return;
   std::cout<<"Configuration file read and confirmed"<<std::endl;

   fMeasuredPOT = outPOT;
   // if just numbers roll forward, if taking form a chain pull out
   //  just enough information to do the oscillations, if hist grab it

   fMethod = nsc->GetDataMethod();
   fNumConfig = nsc->GetNumberConfig();

   Initialize();    //Load in any values and prepare listings
   std::cout<<"Initialization complete"<<std::endl;

   if(fMethod != 4)  RunStandardApproach();
   else              RunFromGrid();
 //Now we have all the data points
   // now we write it out to file
   WriteToFile(output);
}

void NueSensitivity::RunStandardApproach()
{
   NSCDataParam DPst13 = nsc->GetSinS2Th13();
   NSCDataParam DPdelta = nsc->GetDelta();
   NSCDataParam DPdm23 = nsc->GetDeltaMS23();

   fZeroValBg = new float[fNumConfig];
   fZeroValSig   = new float[fNumConfig];

   //int total = 
   SetupChi2Histograms(DPst13, DPdelta, DPdm23);

   SetOscParamBase(0.0024, 0.0, 0, 1);

   //Determine the number of background at Ue3 = 0
   Oscillate(0.0024, 0.0, 0, 1);
   for(int i = 0; i < fNumConfig; i++){
      fZeroValBg[i] = fZeroValSig[i] = 0.0;
      GetPoint(i, fZeroValBg[i], fZeroValSig[i]);
      std::cout<<"For configuration "<<i<<" starting with "
          <<fZeroValBg[i]<<", "<<fZeroValSig[i]<<" events."<<std::endl;
   } 

   //And now the big run over all the numbers
   bool d2,d3;  //some bools to keep the loops right
   float bg, sig;
   double pi = 3.1415926;
   int count = 0;

   for(double sins2t13 = DPst13.start; sins2t13 <= DPst13.end; sins2t13 += t13Step)
   {
      double Th13 = TMath::ASin(TMath::Sqrt(sins2t13))/2.;
      fOscCalc.SetOscParam(OscPar::kTh13, Th13);

      d2 = DPdelta.isfixed;
      for(double delta = DPdelta.start; delta <= DPdelta.end || d2 ; delta += dStep)
      {
         d3 = DPdm23.isfixed;

         fOscCalc.SetOscParam(OscPar::kDelta, delta*pi);
         for(double dm23 = DPdm23.start; dm23 <= DPdm23.end || d3; dm23+= dm23Step)
         {
            int sign = 1;
            fOscCalc.SetOscParam(OscPar::kDeltaM23, dm23*1e-3);

            Oscillate(dm23*1e-3, Th13, delta*pi, sign);   //Change the values
            for(int i = 0; i < fNumConfig; i++){
	      GetPoint(i, bg, sig);
              float chi2 = CalculateFitChi2(i, bg, sig);
              chi2n[i]->Fill(sins2t13, delta, dm23, chi2);
	    }

            if(fMethod != NueSenseConfig::kAllNumbers){
              sign = -1;
              fOscCalc.SetOscParam(OscPar::kDeltaM23, sign*dm23*1e-3);
              Oscillate(dm23*1e-3, Th13, delta*pi, sign);   //Change the values

              for(int i = 0; i < fNumConfig; i++){
                GetPoint(i, bg, sig);
                float chi2 = CalculateFitChi2(i, bg, sig);
                chi2i[i]->Fill(sins2t13, delta, dm23, chi2);
              }
            }
            if(DPdm23.isfixed){ d3 = false;  dm23 = DPdm23.end + 1; }
            count++;
            if(count%100000 == 0) std::cout<<"On iteration "<<count<<std::endl;
         }//end of dm23
         if(DPdelta.isfixed){ d2 = false;  delta = DPdelta.end + 1; }
      }//end of delta
   }//end of sins2theta13
}

void NueSensitivity::RunFromGrid()
{ 
   SetupGridRun();
   std::cout<<"Setup is Complete"<<std::endl;

   double dmDV = fDeltaM32MapNormal.begin()->first;

   std::map<double, std::vector<Point> >::iterator delBeg = fDeltaM32MapNormal[dmDV].begin();
   std::map<double, std::vector<Point> >::iterator delEnd = fDeltaM32MapNormal[dmDV].end();

   while(delBeg != delEnd){
     double delta = delBeg->first;
     for(unsigned int i = 0; i < delBeg->second.size(); i++){
        // at this point
        double ss2th = delBeg->second[i].th13;

        double signal = delBeg->second[i].nsignal;
        double bg = delBeg->second[i].nbg;

        for(int j = 0; j < fNumConfig; j++){
          double val = CalculateFitChi2(j, bg, signal);
          chi2n[j]->Fill(ss2th, delta, dmDV, val);
        }
    }
    delBeg++;
  }
  dmDV = fDeltaM32MapInvert.begin()->first;
  delBeg = fDeltaM32MapInvert[dmDV].begin();
  delEnd = fDeltaM32MapInvert[dmDV].end();
  
  while(delBeg != delEnd){
    double delta = delBeg->first;
    for(unsigned int i = 0; i < delBeg->second.size(); i++){
       // at this point
       double ss2th = delBeg->second[i].th13;
       double signal = delBeg->second[i].nsignal;
       double bg = delBeg->second[i].nbg;

       for(int j = 0; j < fNumConfig; j++){
          double val = CalculateFitChi2(j, bg, signal);
          chi2i[j]->Fill(ss2th, delta, dmDV, val);
       }
    }
    delBeg++;
  }
}

void NueSensitivity::SetupGridRun()
{ 
  std::vector<std::string> datafiles = nsc->GetDataFiles();
  float inPOT = nsc->GetPOT();
  float scale = fMeasuredPOT/inPOT;
 
  double dmVal = 0.0; 
  for(unsigned int i = 0; i < datafiles.size(); i++){
     double th13, delta, mass, ni, bg[5], sig, dummy;
     std::ifstream ins(datafiles[i].c_str());

     ins>>th13>>delta>>mass>>ni>>bg[0]>>bg[1]>>bg[2]>>bg[3]>>sig>>dummy;
     if(i == 0) dmVal = mass;
     while(!ins.eof()){
       Point temp;
       double temp2 = TMath::Sin(2*th13);

       temp.th13 = temp2*temp2;
       temp.nsignal = sig*scale;
       temp.nbg = (dummy-sig)*scale;

       if(ni > 0)  fDeltaM32MapNormal[mass][delta].push_back(temp);
       else        fDeltaM32MapInvert[mass][delta].push_back(temp);

       ins>>th13>>delta>>mass>>ni>>bg[0]>>bg[1]>>bg[2]>>bg[3]>>sig>>dummy;
     }
     ins.clear();
   }

   std::map<double, std::vector<Point> >::iterator delBeg = fDeltaM32MapNormal[dmVal].begin();
   std::map<double, std::vector<Point> >::iterator delEnd = fDeltaM32MapNormal[dmVal].end();

   std::vector<double> deltaM2Pos;
   std::map<double, std::map<double, std::vector<Point> > >::iterator dmBeg = fDeltaM32MapNormal.begin();
   std::map<double, std::map<double, std::vector<Point> > >::iterator dmEnd = fDeltaM32MapNormal.end();

   while(dmBeg != dmEnd){
     deltaM2Pos.push_back(dmBeg->first);
     dmBeg++;
   }

   std::vector<double> deltaPos;
   while(delBeg != delEnd){
     deltaPos.push_back(delBeg->first);
     delBeg++;
   }

   std::vector<double> sin22th13Pos;
   std::vector<Point>::iterator thBeg = fDeltaM32MapNormal[dmVal][0].begin();
   std::vector<Point>::iterator thEnd = fDeltaM32MapNormal[dmVal][0].end();

   while(thBeg != thEnd){
     sin22th13Pos.push_back(thBeg->th13);
     thBeg++;
   }

   sort(deltaM2Pos.begin(), deltaM2Pos.end());
   sort(deltaPos.begin(), deltaPos.end());
   sort(sin22th13Pos.begin(), sin22th13Pos.end());

   int nDM2 = deltaM2Pos.size();
   int nDelta = deltaPos.size();
   int nTh13 = sin22th13Pos.size();

   if(nDelta <= 1){
      std::cout<<"Only zero/one delta value??"<<std::endl;
      return;
   }
   if(nTh13 <= 1){
      std::cout<<"Only zero/one theta value??"<<std::endl;
      return;
   }

   double* deltaM2Array = new double[nDM2+1];
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
     offset = (curr - prev)/2.0;
     if(i == 1) sinthArray[0] = prev - offset;
     sinthArray[i] = curr - offset;
   }
   sinthArray[nTh13] = sinthArray[nTh13-1] + 2*offset;

   offset = 0.00002;
   deltaM2Array[0] = deltaM2Pos[0];
   for(int i = 1; i < nDM2; i++){
     double prev = deltaM2Pos[i-1];
     double curr = deltaM2Pos[i];
     offset = (curr - prev)/2.0;
     if(i == 1) deltaM2Array[0] = prev - offset;
     deltaM2Array[i] = curr - offset;
   }
   deltaM2Array[nDM2] = deltaM2Array[0] + 2*offset;

   TString n1 = "Chi2HistN";
   TString n2 = "Chi2HistI";

   chi2n = new TH3D*[fNumConfig];
   chi2i = new TH3D*[fNumConfig];

   for(int i = 0; i < fNumConfig; i++){
     char temp[20];
     sprintf(temp, "chi2normal%d", i);
     chi2n[i] = new TH3D(temp, temp, nTh13, sinthArray, nDelta, deltaArray, nDM2, deltaM2Array);
     sprintf(temp, "chi2inverted%d", i);
     chi2i[i] = new TH3D(temp, temp, nTh13, sinthArray, nDelta, deltaArray, nDM2, deltaM2Array);

     chi2n[i]->SetDirectory(0);
     chi2i[i]->SetDirectory(0);
   }
   std::cout<<"Chi2 Histograms created ("<<nTh13<<", "<<nDelta<<", "<<nDM2<<"): Running Grid Method"<<std::endl;

   //So now I have built all the histograms as appropriate 
   // - now I just loop over the data

}

void NueSensitivity::GetPoint(int i, float &bg, float &sig)
{
   NSCErrorParam err = nsc->GetErrorConfig(i);  

   float bgt =  currentVal[ClassType::NC]  * (1.0 + err.scale[ClassType::NC])
             + currentVal[ClassType::numu]* (1.0 + err.scale[ClassType::numu])
             + currentVal[ClassType::bnue]* (1.0 + err.scale[ClassType::bnue])
             + currentVal[ClassType::nutau]*(1.0 + err.scale[ClassType::nutau]);

   bg = bgt;
   sig = currentVal[ClassType::nue]*(1.0 + err.scale[ClassType::nue]);
}

void NueSensitivity::Initialize()
{  
   float inPOT = nsc->GetPOT();
   float scale = fMeasuredPOT/inPOT;

   //Initilize the background
   if(fMethod == NueSenseConfig::kAllNumbers || fMethod == NueSenseConfig::kNueHist){
     currentVal[ClassType::NC] = nsc->GetNumber(ClassType::NC)*scale;
     currentVal[ClassType::numu] = nsc->GetNumber(ClassType::numu)*scale;
     currentVal[ClassType::bnue] = nsc->GetNumber(ClassType::bnue)*scale;
     currentVal[ClassType::nutau] = nsc->GetNumber(ClassType::nutau)*scale;
   }

   if(fMethod == NueSenseConfig::kAllNumbers){
     float nuenumber = nsc->GetNumber(ClassType::nue);
     if(nsc->ShouldDeOsc()){
        //Assumed originally Oscillated to:
        //  pow(TMath::Sin(theta23),2)*4.*UE32*(1-UE32)*pow(oscterm,2);
        float UE32 = nsc->GetOldUe3Square();
        float deoscweight = UE32*(1-UE32);
        signuenum = nuenumber/deoscweight*scale;
     }else{
        signuenum = nuenumber*scale;
     }
   }

   if(fMethod == NueSenseConfig::kNueHist){
      TFile inf(nsc->GetNueHistFile().c_str());
      TH1D *h = dynamic_cast<TH1D*>(inf.Get(nsc->GetNueHistName().c_str()));
      if(h == NULL){
        TH1D *h2 = dynamic_cast<TH1D*>(inf.Get(nsc->GetNueHistName().c_str()));
        if(h2 == NULL)
          std::cout<<"Failure to load signal histogram"<<std::endl;
        //This isn't really setup yet, but ideally I would convert it over to the right type
      }

      TH1D* hist = dynamic_cast<TH1D*>(h->Clone());
      hist->SetDirectory(0);

      if(nsc->ShouldDeOsc()){
        hist->Reset();

        float UE32 = nsc->GetOldUe3Square();
        float dms23 = nsc->GetOldDeltaMSquare();                                                
        for(int i = 0; i < h->GetNbinsX()+1; i++){
          float E = h->GetBinCenter(i);
          float num = h->GetBinContent(i);

          double invKmToeV = 1.97e-10; //convert 1/km to eV
          double Delta13Old = dms23*1e-3*735/(4.*E*1e9*invKmToeV);        
          double unweight = UE32*(1-UE32)*4*(1.0/2.0)
                             *TMath::Power(TMath::Sin(Delta13Old),2);

          hist->Fill(E, num/unweight*scale);
        }
     }else{
       hist->Scale(scale);
     }

     signuehist = dynamic_cast<TH1D*>(hist->Clone());
     signuehist->SetDirectory(0);
     reweight = dynamic_cast<TH1D*>(signuehist->Clone());
     reweight->SetDirectory(0);

     fBinCenter = new double[reweight->GetNbinsX()+1];
     for(int i = 0; i < reweight->GetNbinsX()+1; i++){
         fBinCenter[i] = reweight->GetBinCenter(i);
     }
     
   }

   if(fMethod == NueSenseConfig::kAnaNueFiles)  
      LoadEventsFromFile();


   fDeltaMS12 = nsc->GetDeltaMS12();
   float temp = nsc->GetSinS2Th12();
   fTh12 = TMath::ASin(TMath::Sqrt(temp))/2.;
   temp = nsc->GetSinS2Th23();
   fTh23 = TMath::ASin(TMath::Sqrt(temp))/2.;
   
   fDensity = nsc->GetDensity();

}

void NueSensitivity::LoadEventsFromFile()
{
  for(int i = 0; i < nsc->GetNumFiles(); i++){
   
     TChain *selected = new TChain("ana_nue");
     TChain *pots = new TChain("pottree");
    
     selected->Add(nsc->GetFile(i).c_str());
     pots->Add(nsc->GetFile(i).c_str());
 
     int nonOscFlavor, nuFlavor, type;
     float trueE;
     double skzpweight;        

     selected->SetMakeClass(1);
     selected->SetBranchStatus("*", 0);
     selected->SetBranchStatus("mctrue*", 1);
     selected->SetBranchStatus("fluxweight.totbeamweight", 1);
     selected->SetBranchAddress("mctrue.nonOscNuFlavor",&nonOscFlavor);
     selected->SetBranchAddress("mctrue.nuFlavor",&nuFlavor);
     selected->SetBranchAddress("mctrue.fNueClass",&type);
     selected->SetBranchAddress("mctrue.nuEnergy", &trueE);
     selected->SetBranchAddress("fluxweight.totbeamweight",&skzpweight);
        
     NuePOT *np;                                                             
     pots->SetBranchAddress("NuePOT", &np);
     double filePOT, scale;
     filePOT = 0.0;                                                                     

     for(int j=0; j < pots->GetEntries(); j++) {
         pots->GetEntry(j);
         filePOT += np->pot;
     }
     scale = filePOT/fMeasuredPOT;
        
     for(int j = 0; j < selected->GetEntries(); j++){
        selected->GetEntry(j);
        mininfo newInfo;
        newInfo.energy = trueE;
        newInfo.NuFlavor = nuFlavor;
        newInfo.NuFlavorBeforeOsc = nonOscFlavor;
        newInfo.nuClass = type;
        if(skzpweight < 0) skzpweight = 1.0;
        newInfo.weight = skzpweight*scale;
        mcInfo.push_back(newInfo);
     }  //End of this file

     delete selected;  
     delete pots; 
  }  //Done loading all files and POT Normalization is taken care of
}

float NueSensitivity::Oscillate(double dm23, double Ue32, double delta, int sign)
{
  if(fMethod == NueSenseConfig::kAllNumbers)
     return OscillateNumber(Ue32);

  if(fMethod == NueSenseConfig::kNueHist)
     return OscillateHist(dm23, Ue32, delta, sign);
   
  if(fMethod == NueSenseConfig::kAnaNueFiles)
     return OscillateFile(dm23, Ue32, delta, sign);

  return -1;
}

float NueSensitivity::OscillateNumber(double Ue32)
{
//  double Ue3 = TMath::Sin(t13);

  //Assumed originally Oscillated to:
  //  pow(TMath::Sin(theta23),2)*4.*UE32*(1-UE32)*pow(oscterm,2);
  currentVal[ClassType::nue] = signuenum*Ue32*(1-Ue32);
                                                                                           
  return currentVal[ClassType::nue];
}

float NueSensitivity::OscillateHist(double /*dm23*/, double /*th13*/, double /*delta*/, int /*sign*/)
{
  reweight->Reset("ICE");

  for(int i = 0; i < signuehist->GetNbinsX()+1; i++){
     float E = fBinCenter[i];
     if(E > 10) break;
     float num = signuehist->GetBinContent(i);
     float weight = 1.0;
     if(num > 0)    weight =OscillateMatter(12,14,E); //,dm23,th13,delta,sign);
                                                                                
     reweight->Fill(E, num*weight);
  }
 
  currentVal[ClassType::nue] = reweight->GetSum();
//  std::cout<<dm23<<"  "<<Ue32<<"  "<<delta<<"  "<<sign<<"  "<<reweight->GetSum()<<std::endl;
  return reweight->GetSum();
}                                                                                

float NueSensitivity::OscillateFile(double /*dm23*/, double /*Ue32*/, double /*delta*/, int /*sign*/)
{
  float total[5];

  for(int i = 0; i < 5; i++) total[i] = 0;

  for(unsigned int i = 0; i < mcInfo.size(); i++){
    float E = mcInfo[i].energy;
    int inu = mcInfo[i].NuFlavor;
    int inunoosc = mcInfo[i].NuFlavorBeforeOsc;
    int nuClass = mcInfo[i].nuClass;
    double weight = mcInfo[i].weight;

    double oscweight = 
      OscillateMatter(inu,inunoosc,E);

    total[nuClass] += weight*oscweight;
  }

  for(int i = 0; i < 5; i++){
    currentVal[i] = total[i];
  }

  return currentVal[ClassType::nue];
}

float NueSensitivity::CalculateChi2(int i, float bg, float sig)
{
  NSCErrorParam err = nsc->GetErrorConfig(i);

  float nObs = sig+bg;
  float nZero = fZeroValBg[i] + fZeroValSig[i];

  double bgerr = err.bg_systematic*fZeroValBg[i];
  double sigerr = err.sig_systematic*fZeroValSig[i];

  float syserrsq = bgerr*bgerr + sigerr*sigerr;
 
  float errScale = nObs/(nObs + syserrsq);
                                                                                
  float chi2 = 2*(nZero - nObs + nObs*TMath::Log(nObs/nZero)) * errScale;
  return chi2;
}

float NueSensitivity::CalculateFitChi2(int i, float bg, float sig)
{
  NSCErrorParam err = nsc->GetErrorConfig(i);
                                                                                
  float nZero = sig+bg;
  float nObs = fObserved;
  if(!fObservedIsSet) nObs = fZeroValBg[i] + fZeroValSig[i];

  double bgerr = err.bg_systematic*bg;
  double sigerr = err.sig_systematic*sig;

  float syserrsq = bgerr*bgerr + sigerr*sigerr;                                               
  float errScale = nZero/(nZero + syserrsq);
                                                                                
  float chi2 = 2*(nZero - nObs + nObs*TMath::Log(nObs/nZero)) * errScale;
  return chi2;
}


int NueSensitivity::SetupChi2Histograms(NSCDataParam DPst13, NSCDataParam DPdelta, NSCDataParam DPdm23)
{
   chi2n = new TH3D*[fNumConfig];
   chi2i = new TH3D*[fNumConfig];
                                                                                           
   int t13Bins = int((DPst13.end - DPst13.start)/t13Step) + 1;
   int delBins = int((DPdelta.end - DPdelta.start)/dStep) + 1;
   int dm23Bins = int((DPdm23.end - DPdm23.start)/dm23Step) + 1;
   
   double xstart, xend, ystart, yend, zstart, zend;
   xstart = xend = ystart = yend = zstart = zend = 0;
                                                                                        
   if(DPst13.isfixed){ // i have no idea why you are running the code
   }else{
      xstart = DPst13.start - t13Step/2; xend = xstart + t13Bins*t13Step;
   }
   if(DPdelta.isfixed){
      ystart = DPdelta.start - 3*dStep/2;
      delBins = 5;
   }else{
      ystart = DPdelta.start - dStep/2;
   }
   yend = ystart + delBins*dStep;

   if(DPdm23.isfixed){
      zstart = DPdm23.start - 3*dm23Step/2;
      dm23Bins = 5;
   }else{
      zstart = DPdm23.start - dm23Step/2;
   }
   zend = zstart + dm23Bins*dm23Step;

   for(int i = 0; i < fNumConfig; i++){
      char temp[20];
      sprintf(temp, "chi2normal%d", i);
      chi2n[i] = new TH3D(temp, temp, t13Bins, xstart, xend, delBins,
                             ystart, yend, dm23Bins, zstart, zend);
      sprintf(temp, "chi2inverted%d", i);
      chi2i[i] = new TH3D(temp, temp, t13Bins, xstart, xend, delBins,
                             ystart, yend, dm23Bins, zstart, zend);
   }
   std::cout<<"Chi2 Histograms created ("<<t13Bins<<", "<<delBins<<", "<<dm23Bins<<"): " 
            <<t13Bins*delBins*dm23Bins<<" iterations to perform"<<std::endl;
   return t13Bins*delBins*dm23Bins;
}

void NueSensitivity::WriteToFile(std::string file)
{
   TFile out(file.c_str(), "RECREATE");
   out.cd();

   for(int i = 0; i < fNumConfig; i++){
      chi2n[i]->Write();
      chi2i[i]->Write();
   }
   TTree *info, *err;

   ProduceTree(info, err);
   info->Write();
   err->Write();
   out.Close();
}

double NueSensitivity::OscillateMatter(int nuFlavor, int nonOscNuFlavor,
                                    float Energy, float dm2, float th13,
                                    float delta, int hierarchy)
{
 Double_t x[1] = {};
  x[0] = Energy;
  Double_t dm2_12 = fDeltaMS12*1e-5; //best fit SNO
  Double_t dm2_23 = dm2;
                                                                                
                                                                                
  Double_t par[9] = {0};
  par[0] = fBaseLine;
  par[1] = fTh23;
  par[2] = fTh12;
  par[3] = th13;
  par[4] = hierarchy*dm2_23;
  par[5] = dm2_12;
  par[6] = fDensity; //standard rock density
  par[7] = delta;
  par[8] = 1;
  if(nonOscNuFlavor < 0) par[8] = -1;

  //std::cout<<"About to call "<<dm2<<"  "<<ss13<<"  "<<delta<<std::endl;
  fOscCalc.SetOscParam(par);
  return fOscCalc.Oscillate(nuFlavor, nonOscNuFlavor, Energy);
}

double NueSensitivity::OscillateMatter(int nuFlavor, int nonOscNuFlavor,
                                    float Energy)
{

  if(nonOscNuFlavor > 0) fOscCalc.SetOscParam(OscPar::kNuAntiNu, 1);
  if(nonOscNuFlavor < 0) fOscCalc.SetOscParam(OscPar::kNuAntiNu, -1);

  return fOscCalc.Oscillate(nuFlavor, nonOscNuFlavor, Energy);
}

void NueSensitivity::SetOscParamBase( float dm2, float ss13,
                                    float delta, int hierarchy){

  Double_t dm2_12 = fDeltaMS12*1e-5; //best fit SNO
  Double_t dm2_23 = dm2;
                                                                                
  Double_t par[9] = {0};
  par[OscPar::kL] = fBaseLine;
  par[OscPar::kTh23] = fTh23;
  par[OscPar::kTh12] = fTh12;
  par[OscPar::kTh13] = ss13; // TMath::ASin(TMath::Sqrt(ss2th13))/2.;
  par[OscPar::kDeltaM23] = hierarchy*dm2_23;
  par[OscPar::kDeltaM12] = dm2_12;
  par[OscPar::kDensity] = fDensity; //standard rock density
  par[OscPar::kDelta] = delta;
  par[OscPar::kNuAntiNu] = 1;
                                                                                
//  std::cout<<"About to call "<<dm2<<"  "<<ss13<<"  "<<delta<<std::endl;
  fOscCalc.SetOscParam(par);
}

void NueSensitivity::ProduceTree(TTree * &infotree, TTree * &errtree)
{
  infotree = new TTree("OscPar","OscPar");

  double Th23, Th12, dm2_12;
  dm2_12 = nsc->GetDeltaMS12()*1e-5;
  Th12 = nsc->GetSinS2Th12();
  Th23 = nsc->GetSinS2Th23();

  infotree->Branch("L",&fBaseLine, "L/F");
  infotree->Branch("Sin2_2Theta23",&Th23, "Sin2_2Theta23/D");
  infotree->Branch("Sin2_2Theta12",&Th12, "Sin2_2Theta12/D");
  infotree->Branch("DeltaMS12",&dm2_12, "DeltaMS12/D");
  infotree->Branch("Density",&fDensity, "Density/D");

  NSCDataParam DPst13 = nsc->GetSinS2Th13();
  NSCDataParam DPdelta = nsc->GetDelta();
  NSCDataParam DPdm23 = nsc->GetDeltaMS23();

  infotree->Branch("DeltaMS23_start",&DPst13.start, "DeltaMS23_start/F");
  infotree->Branch("DeltaMS23_end",&DPst13.end, "DeltaMS23_end/F");
  infotree->Branch("Delta_start",&DPdelta.start, "Delta_start/F");
  infotree->Branch("Delta_end",&DPdelta.end, "Delta_end/F");
  infotree->Branch("Sin2_2Theta13_start",&DPdm23.start, "Sin2_2Theta13_start/F");
  infotree->Branch("Sin2_2Theta13_end",&DPdm23.end, "Sin2_2Theta13_end/F");

  infotree->Fill();

  errtree = new TTree("ErrPar","ErrPar");
  
  double bgsys, sigsys, nc, numu, nue, nutau, bnue;
  errtree->Branch("BG_systematic",&bgsys,"BG_systematic/D");
  errtree->Branch("SIG_systematic",&sigsys, "SIG_systematic/D");
  errtree->Branch("nc_scale",&nc, "nc_scale/D");
  errtree->Branch("numu_scale",&numu, "numu_scale/D");
  errtree->Branch("nue_scale",&nue, "nue_scale/D");
  errtree->Branch("nutau_scale",&nutau, "nutau_scale/D");
  errtree->Branch("bnue_scale",&bnue, "bnue_scale/D");

  for(int i = 0; i < nsc->GetNumberConfig(); i++){
      NSCErrorParam err = nsc->GetErrorConfig(i);
      bgsys = err.bg_systematic;
      sigsys = err.sig_systematic;
      nc = err.scale[ClassType::NC];
      numu = err.scale[ClassType::numu];
      bnue =err.scale[ClassType::bnue];
      nutau = err.scale[ClassType::nutau];
      nue = err.scale[ClassType::nue];
      
      errtree->Fill();
  }
}




