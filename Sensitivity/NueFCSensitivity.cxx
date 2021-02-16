#include <string>
#include "NueAna/NueAnaTools/NueConvention.h"
#include "NueFCSensitivity.h"
#include "TMath.h"
#include "TChain.h"
#include "NueAna/NuePOT.h"
#include <iostream>
#include <fstream>
#include "TMinuit.h"

double Poisson(double mean, int n, bool numOnly)
{
   double numerator = 1;
   if(mean > 100 || n > 100){
      // use roots protected function -
      if(!numOnly) return TMath::Poisson(n, mean);

      // else     pull a trick to prevent explosions
      double logN = n*TMath::Log(mean) - mean;
      numerator = TMath::Exp(logN);
   } 
   else{
    numerator = TMath::Power(mean,n)*TMath::Exp(-mean);   
   }
   double denom = 1.0;
   if(!numOnly) denom = TMath::Factorial(n);
   return numerator/denom;
} 
 
double Gaussian(double x, double mean, double sig)
{ 
  double front = 1/TMath::Sqrt(2*3.1415926*sig);
  double exp = (x-mean)*(x-mean)/(2*sig);
  return front*TMath::Exp(-exp);
} 

void WrapperFunction(Int_t & /*npar*/, Double_t * /*gin*/, Double_t &f, Double_t *par, Int_t /*iflag*/)
{
   NueFCSensitivity * fcs = dynamic_cast<NueFCSensitivity *>(gMinuit->GetObjectFit());

//   std::cout<<"In the wrapper"<<par[0]<<"  "<<par[1]<<std::endl;
   f = fcs->MinimizationFunction(par);
}

ClassImp(NueFCSensitivity)

NueFCSensitivity::NueFCSensitivity(){
  t13Step = 0.001;
  dStep = 0.02;
  dm23Step = 0.025;
  fBaseLine = 735;
  fObservedIsSet = false;
  fDeltaStart = -1;
  fDeltaEnd = 3;
  fTheta23 = -1e4;
  for(int i = 0; i < 6; i++)  fFixedFitPar[i] = 0.0;
}


NueFCSensitivity::~NueFCSensitivity(){
}

double NueFCSensitivity::MinimizationFunction(double* par)
{
   int n = (int) fFixedFitPar[0];
   double s = fFixedFitPar[1];
   double b = fFixedFitPar[2];
   double sigBg = fFixedFitPar[3];
   double k = fFixedFitPar[4];
   double sigK = fFixedFitPar[5];

//   std::cout<<n<<"  "<<b<<"  "<<s<<std::endl;

   double bFit = par[0];
   double kFit = par[1];

   double tot = s*kFit +bFit;
   if(tot < 0) std::cout<<s<<" "<<kFit<<"  "<<bFit<<" serious error "<<std::endl;

   double lnF = -n*TMath::Log(s*kFit + bFit) + (s*kFit + bFit) + 0.5*1/sigBg*(bFit-b)*(bFit-b);
          lnF +=  0.5*1/sigK*(kFit-k)*(kFit-k);

//   std::cout<<"res: "<<kFit<<"  "<<bFit<<"  "<<lnF<<std::endl;
   fResult = lnF;

   return  lnF;
}

void NueFCSensitivity::TestCode()
{

   double results[3];
   fObserved = 35;
   FindBestFit(35, 1, 26.5, 4, 1, 0.07, results);

   std::cout<<results[0]<<"  "<<results[1]<<"  "<<results[2]<<std::endl;
  
   TMinuit* min = gMinuit;
   min->DefineParameter(0,"bfit",35,0.01, 0,60);
   min->DefineParameter(1,"kfit",1,0.01, 0,22);

   FindBestFit(35, 1, 26.5, 4, 1, 0.07, results);

   std::cout<<results[0]<<"  "<<results[1]<<"  "<<results[2]<<std::endl;

   min->DefineParameter(0,"bfit",35,3, 0,60);
   min->DefineParameter(1,"kfit",1,3, 0,22);
 
   FindBestFit(35, 1, 26.5, 4, 1, 0.07, results);
 
   std::cout<<results[0]<<"  "<<results[1]<<"  "<<results[2]<<std::endl;



}


void NueFCSensitivity::SetDeltaRange(double start, double end)
{
   if(start > -1e-3) fDeltaStart = start;
   if(end > -1e-3)   fDeltaEnd = end;
}

void NueFCSensitivity::Run(std::string input, std::string output, double outPOT)
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

void NueFCSensitivity::RunStandardApproach()
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

void NueFCSensitivity::RunFromGrid()
{
   SetupGridRun();
   std::cout<<"Setup is Complete"<<std::endl;  
 
   double dmDV = fDeltaM32MapNormal.begin()->first;

   std::map<double, std::vector<Point> >::iterator delBeg = fDeltaM32MapNormal[dmDV].begin();
   std::map<double, std::vector<Point> >::iterator delEnd = fDeltaM32MapNormal[dmDV].end();

   int ndCount = 0;
   double th23scale = 1.0;

   if(fTheta23 > -100) th23scale = 2*TMath::Sin(fTheta23)*TMath::Sin(fTheta23);

   while(delBeg != delEnd){
     double delta = delBeg->first;
     if(delta - fDeltaStart > -1e-4 && delta - fDeltaEnd < -1e-4){ 
       if(ndCount%10 == 0) std::cout<<"Running: "<<delta<<"  "<<fDeltaStart<<"  "<<fDeltaEnd<<std::endl;

       for(unsigned int i = 0; i < delBeg->second.size(); i++){
          // at this point
          double ss2th = delBeg->second[i].th13 * th23scale;
           
 
          double signal = delBeg->second[i].nsignal; 
          double bg = delBeg->second[i].nbg;

//          if(ss2th < 0.426) continue;
//          std::cout<<ss2th<<"  "<<signal<<"  "<<std::endl;
        // For each 
          for(int j = 0; j < fNumConfig; j++){
            double omega = EvaluateOmega(signal, bg, j);
//          double val = CalculateFitChi2(j, bg, signal);
            chi2n[j]->Fill(ss2th, delta, dmDV, omega);
            nsigN[j]->Fill(ss2th, delta, dmDV, signal);
            nbgN[j]->Fill(ss2th, delta, dmDV, bg);
          }
       }
       ndCount++;
    }
    delBeg++;
  }

   std::cout<<"Finished Normal"<<std::endl;

  dmDV = fDeltaM32MapInvert.begin()->first;
  delBeg = fDeltaM32MapInvert[dmDV].begin();
  delEnd = fDeltaM32MapInvert[dmDV].end();

  ndCount = 0;

  while(delBeg != delEnd){
    double delta = delBeg->first;
    if(delta - fDeltaStart > -1e-4 && delta - fDeltaEnd < -1e-4){
      for(unsigned int i = 0; i < delBeg->second.size(); i++){
         // at this point
         double ss2th = delBeg->second[i].th13 * th23scale;
         double signal = delBeg->second[i].nsignal;
         double bg = delBeg->second[i].nbg;

         for(int j = 0; j < fNumConfig; j++){
            double omega =  EvaluateOmega(signal, bg, j);
            // CalculateFitChi2(j, bg, signal);
            chi2i[j]->Fill(ss2th, delta, dmDV, omega);
            nsigI[j]->Fill(ss2th, delta, dmDV, signal);
            nbgI[j]->Fill(ss2th, delta, dmDV, bg);
         }  
      }
      ndCount++;
      if(ndCount%10 == 0) std::cout<<"Finished inv delta = "<<delBeg->first<<std::endl;
    }
    delBeg++;
    
  }
}

double NueFCSensitivity::EvaluateOmega(double signal, double background, int i)
{
   int n0 = fObserved;
   double b0 = background;
   double s0 = signal;
   double k0 = 1.0;

   NSCErrorParam err = nsc->GetErrorConfig(i);

   double errK = 1.0;
   if(s0 > 0) errK = (err.sig_systematic)*(err.sig_systematic);
   if(errK < 1e-5){ errK = 1.0;  std::cout<<"can't do perfect signal"<<std::endl; }
   double errBg = err.bg_systematic*b0*err.bg_systematic*b0;
 
   double rank0 = CalculateRank(n0, s0, b0, errBg, k0, errK);

   double results[3];
   FindBestFit(n0, s0, b0, errBg, k0, errK, results);

   double betaHat = results[0];
   double kHat = results[1];

   errBg = (betaHat*err.bg_systematic)*(betaHat*err.bg_systematic);
   errK = (err.sig_systematic)*(err.sig_systematic)*kHat*kHat;

//   std::cout<<n0<<"  "<<s0<<"  "<<b0<<"  "<<k0<<"  "<<betaHat<<" "<<kHat<<"  "<<errBg<<" "<<errK<<std::endl;
//   std::cout<<s0*kHat + betaHat<<"  "<<n0<<std::endl;

   double omega = 0;
   double omegaBar = Poisson(s0*kHat + betaHat, n0);

   int nBase = int(s0*kHat + betaHat);

   int nHigh = nBase;
   int nLow  = nBase - 1;

   double delta = 1.0;
   bool filled = false;

   const int NUMB = 30;
   const int NUMK = 30;

   static double bVal[NUMB];
   static double kVal[NUMK];
   static double errBVar[NUMB];
   static double errKVar[NUMK];

   static bool first = true;
   static double scale = 1.0;

   double bStart = 0.4*b0;   double bStop = 2.2*b0;
   double kStart = 0.4;   double kStop = 1.6;

   if(first){
      for(int i = 0; i < NUMB; i++) { 
          bVal[i] = bStart + i*(bStop - bStart)/float(NUMB); 
          errBVar[i] = err.bg_systematic*err.bg_systematic*bVal[i]*bVal[i];
      }
  
      for(int i = 0; i < NUMK; i++) { 
          kVal[i] = kStart + i*(kStop - kStart)/float(NUMK);
          errKVar[i] = (err.sig_systematic)*(err.sig_systematic)*kVal[i]*kVal[i];
       }
      first = false; 
      scale = (bStop-bStart)*(kStop-kStart)/(NUMK*NUMB);
   }
   double bkProb[NUMB][NUMK];

   /*   so lets be clear about this, the values for these gaussians are the same
    but its the rank that changes for any given value of n (the values themselve differ with mu)
    so the first time through i calculate the contribution for each point in the space
    then when looping through just look it up in a giant array             */

   // better yet - I can save cycles by building the arrays separately on the first pass 

   double bGauss[NUMB];
   double kGauss[NUMK];

   for(int i = 0; i < NUMB; i++) { 
      bGauss[i] = Gaussian(betaHat, bVal[i], errBg);
   }
   for(int i = 0; i < NUMB; i++) { 
      kGauss[i] = Gaussian(kHat, kVal[i], errK);
   }    

   bool risingH = false, risingL = false;
   bool doneH = false;  // some calc savers
   bool doneL = false;  // some calc savers

   double ThreshHold = 1e-5;
   double pfThresh = ThreshHold*1e-2;

   double slip = 0;   // Error on the numerical integral

   while(delta > ThreshHold || (1 - (omega + omegaBar) > 2*ThreshHold) )
   {
      if(nHigh == n0) nHigh++;
      if(nLow == n0) nLow--;

      double dOmegaH = 0.0, dOmegaBarH = 0.0;
      double dOmegaL = 0.0, dOmegaBarL = 0.0;

      double PrefactorH = Poisson(s0*kHat+betaHat, nHigh);
      double PrefactorL = Poisson(s0*kHat+betaHat, nLow);

      if(PrefactorH > pfThresh) risingH = true;
      if(PrefactorL > pfThresh) risingL = true;
      if(PrefactorH < pfThresh && risingH) doneH = true;
      if(PrefactorL < pfThresh && risingL) doneL = true;
      if(doneH && doneL){
          // just a sanity check that the code doesn't think its done too early, this will cause an inf loop
          if(1 - (omega + omegaBar) < 1.5*slip) break;
          std::cout<<"Thats unexpected:  "<<1 - (omega + omegaBar)<<"  "<<slip<<std::endl;
      }

      if(PrefactorH > pfThresh || PrefactorL > pfThresh){  // No need to loop if contribution too small
        for(int i = 0; i < NUMB; i++) //  0 to infinity
        {
          double b = bVal[i];
          double eb = errBVar[i];
          for(int j = 0 ; j < NUMK; j++) //  -inf to inf
          {
              double k = kVal[j];
              if(!filled) bkProb[i][j] = bGauss[i]*kGauss[j];
 
              double val = bkProb[i][j];   
              double eK = errKVar[j];
              if(val < ThreshHold*1e-4) continue; // no point in using these contributions 
              
              double rank = 1.0;

              if(PrefactorH > pfThresh){
                rank = CalculateRank(nHigh, s0, b, eb, k, eK);
                if(rank > rank0) dOmegaH += val*scale;
                else          dOmegaBarH += val*scale; 
              }
 
              if(PrefactorL > pfThresh){
                if(nLow >= 0){
                  rank = CalculateRank(nLow, s0, b, eb, k, eK);
                  if(rank > rank0) dOmegaL += val*scale;
                  else          dOmegaBarL += val*scale; 
                }
              }
          }
        }
        if(!filled) filled = true;          
        if(PrefactorH > pfThresh && 1 - (dOmegaH+ dOmegaBarH) > slip) slip = 1 - (dOmegaH+ dOmegaBarH);
        if(PrefactorL > pfThresh && 1 - (dOmegaL+ dOmegaBarL) > slip) slip = 1 - (dOmegaL+ dOmegaBarL);
      }
      delta = TMath::Max(PrefactorH*dOmegaH + PrefactorL*dOmegaL, PrefactorH*dOmegaBarH + PrefactorL*dOmegaBarL);
      omega += PrefactorH*dOmegaH + PrefactorL*dOmegaL;
      omegaBar += PrefactorH*dOmegaBarH + PrefactorL*dOmegaBarL;
      nHigh++; nLow--;
   }

   if(slip > ThreshHold){
     std::cout<<"Integral error past threshold:  "<<omega<<"  "<<omegaBar<<"  "<<omega+omegaBar<<"  "<<slip<<std::endl;
     std::cout<<n0<<"  "<<s0<<"  "<<b0<<"  "<<k0<<"  "<<betaHat<<" "<<kHat<<"  "<<errBg<<" "<<errK<<std::endl;
     std::cout<<s0*kHat + betaHat<<"  "<<n0<<std::endl;
   }
   return omega;
}


double NueFCSensitivity::CalculateRank(int n, double s, double b, double errBg, double k, double errK)
{
   double results[3]; 
   //Calculate the numerator
   FindBestFit(n, s, b, errBg, k, errK, results);

//   double numerator = results[2];
   
   double sBest, kBest, bBest;

   // Now calculate the denominator 
   if(n >= b) { // then we have n = sk+b, Beta = b; k = k0// don't have to do anything  
        sBest = (n-b)/k; kBest = k; bBest = b;
   }
   else{   //if(n < b) {
        bBest = 0.5*(b - errBg + TMath::Sqrt( (b-errBg)*(b - errBg) + 4*n*errBg)); 
        // then we have no signal at best fit point
        sBest = 0;  kBest = k;
    }

   double betaHat = results[0];
   double kHat = results[1];
   double logN = n*TMath::Log(s*kHat+betaHat) - (s*kHat+betaHat) - (betaHat-b)*(betaHat-b)/(2*errBg)
                      - (kHat-k)*(kHat-k)/(2*errK);
   double logD = n*TMath::Log(sBest*kBest+bBest) - (sBest*kBest+bBest) - (bBest-b)*(bBest-b)/(2*errBg)
                      - (kBest-k)*(kBest-k)/(2*errK);

    double rank = 1.0;
    if(logD != 0) { rank = TMath::Exp(logN-logD); }
    else std::cout<<"odd"<<std::endl;

    return rank;
}

bool NueFCSensitivity::FindBestFit(int n,  double s, double b, double errBg, double k, double errK, double* res)
{
//  At the moment this function is solvable in ~closed form (we ignore the dependance of errBg, errK on khat, betaHat
//  this lets us avoid having a numerical minimization, but the second half of this function is the code that allows for that

   double sol_A = 1 + errK/errBg*s*s;
   double sol_B = errBg + s*k - 2*errK*s*s*b/errBg - b + s*s*errK;
   double sol_C = -n*errBg + s*k*errBg - s*s*errK*b - s*k*b + s*s*b*b*errK/errBg;

   double rad = sol_B*sol_B-4*sol_A*sol_C;
   if(rad < 0) std::cout<<"negative radical - not sure about this...."<<std::endl;

   double betaHat = (-sol_B + TMath::Sqrt(rad))/(2*sol_A);

   //solving for kHat given betaHat
   double kHat = k - s*errK/errBg*(b-betaHat);

   res[0] = betaHat;
   res[1] = kHat;

/*
   fFixedFitPar[0] = n;
   fFixedFitPar[1] = s;
   fFixedFitPar[2] = b;
   fFixedFitPar[3] = errBg;
   fFixedFitPar[4] = k;
   fFixedFitPar[5] = errK;
*/
//   res[2] = fResult = Poisson(s*kHat+betaHat, n, true)*Gaussian(betaHat, b, errBg)*Gaussian(kHat, k, errK);
/*
   static bool first = true;
   static TMinuit *min = 0;

   if(first){    
      min = new TMinuit(2);
      gMinuit = min;
      min->SetFCN(WrapperFunction);
      min->SetObjectFit(this);
      min->DefineParameter(0,"bfit",fObserved,1, 0,60);
      min->DefineParameter(1,"kfit",1,0.99, 0,22);

      const double ERRDEF=1.;
      min->SetErrorDef(ERRDEF);
      min->SetPrintLevel(-1);
      first = false;
     std::cout<<"Built"<<std::endl;
   }
   Double_t arglist[10];
   Int_t ierflg = 0;

   arglist[0] = 10000;       //max calls
   arglist[1] = 0.01;         //tolerance

   min->mnexcm("SIMPLEX",arglist,2,ierflg);

   double errs[2];
   for(int i=0;i<2;i++){
     min->GetParameter(i,res[i],errs[i]);
//     std::cout<<res[i]<<"  "<<errs[i]<<std::endl;
   }
*/
//   res[2] = fResult;

   return true;
}

void NueFCSensitivity::SetupGridRun()
{
  std::vector<std::string> datafiles = nsc->GetDataFiles();
  float inPOT = nsc->GetPOT();
  float scale = fMeasuredPOT/inPOT;

  NSCErrorParam err = nsc->GetErrorConfig(0);
  bool first = true;

  double hold = 2.32e-3;

  for(unsigned int i = 0; i < datafiles.size(); i++){
     std::cout<<"Openning file: "<<datafiles[i]<<std::endl;
     double th13, delta, mass, ni, bg[5], sig, dummy;
     std::ifstream ins(datafiles[i].c_str());

     ins>>th13>>delta>>mass>>ni>>bg[0]>>bg[1]>>bg[2]>>bg[3]>>sig>>dummy;
     hold = mass;
     while(!ins.eof()){
       Point temp;
       double temp2 = TMath::Sin(2*th13);

//       bg[0] *=  1.0 + err.scale[ClassType::NC];
//       bg[1] *=  1.0 + err.scale[ClassType::numu];
//       bg[2] *=  1.0 + err.scale[ClassType::bnue];
//       bg[3] *=  1.0 + err.scale[ClassType::nutau];
//       sig *=  1.0 + err.scale[ClassType::nue];
       dummy = bg[0]+bg[1]+bg[2]+bg[3]+sig;
       temp.th13 = temp2*temp2;
       temp.nsignal = sig*scale;
       temp.nbg = (dummy-sig)*scale;
 
       if(i == 0 && first){ std::cout<< temp.nbg<<std::endl; first = false; }

       if(ni > 0)  fDeltaM32MapNormal[mass][delta].push_back(temp);
       else        fDeltaM32MapInvert[mass][delta].push_back(temp);

       ins>>th13>>delta>>mass>>ni>>bg[0]>>bg[1]>>bg[2]>>bg[3]>>sig>>dummy;
     }
     ins.clear();
   }
   

   std::map<double, std::vector<Point> >::iterator delBeg = fDeltaM32MapNormal[hold].begin();
   std::map<double, std::vector<Point> >::iterator delEnd = fDeltaM32MapNormal[hold].end();

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
   std::vector<Point>::iterator thBeg = fDeltaM32MapNormal[hold][0].begin();
   std::vector<Point>::iterator thEnd = fDeltaM32MapNormal[hold][0].end();

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
   double th23scale = 1.0;
   if(fTheta23 > -100) th23scale = 2*TMath::Sin(fTheta23)*TMath::Sin(fTheta23);

   for(int i = 1; i < nTh13; i++){
     double prev = sin22th13Pos[i-1]*th23scale;
     double curr = sin22th13Pos[i]*th23scale;
     offset = (curr - prev)/2.0;
     if(i == 1) sinthArray[0] = prev - offset;
     sinthArray[i] = curr - offset;
   }
   sinthArray[nTh13] = sinthArray[nTh13-1] + 2*offset;

   offset = 0.00002;
   deltaM2Array[0] = deltaM2Pos[0] - 0.00001;
   deltaM2Array[1] = deltaM2Pos[0] + 0.00001; 
/*   for(int i = 1; i < nDM2; i++){
     double prev = deltaM2Pos[i-1];
     double curr = deltaM2Pos[i];
     offset = (curr - prev)/2.0;
     if(i == 1) deltaM2Array[0] = prev - offset;
     deltaM2Array[i] = curr - offset;
   }*/ 
//   deltaM2Array[nDM2] = deltaM2Array[0] + 2*offset;

   TString n1 = "Chi2HistN";
   TString n2 = "Chi2HistI";

   chi2n = new TH3D*[fNumConfig];
   chi2i = new TH3D*[fNumConfig];
   nsigN = new TH3D*[fNumConfig];
   nsigI = new TH3D*[fNumConfig];
   nbgN = new TH3D*[fNumConfig];
   nbgI = new TH3D*[fNumConfig];  

   for(int i = 0; i < fNumConfig; i++){
     char temp[20];
     sprintf(temp, "chi2normal%d", i);
     chi2n[i] = new TH3D(temp, temp, nTh13, sinthArray, nDelta, deltaArray, nDM2, deltaM2Array);
     sprintf(temp, "chi2inverted%d", i);
     chi2i[i] = new TH3D(temp, temp, nTh13, sinthArray, nDelta, deltaArray, nDM2, deltaM2Array);
 
     chi2n[i]->SetDirectory(0);
     chi2i[i]->SetDirectory(0);
     sprintf(temp, "signormal%d", i);
     nsigN[i] = new TH3D(temp, temp, nTh13, sinthArray, nDelta, deltaArray, nDM2, deltaM2Array);

     sprintf(temp, "siginverted%d", i);
     nsigI[i] = new TH3D(temp, temp, nTh13, sinthArray, nDelta, deltaArray, nDM2, deltaM2Array);

     sprintf(temp, "bgnormal%d", i);
     nbgN[i] = new TH3D(temp, temp, nTh13, sinthArray, nDelta, deltaArray, nDM2, deltaM2Array);

     sprintf(temp, "bginverted%d", i);
     nbgI[i] = new TH3D(temp, temp, nTh13, sinthArray, nDelta, deltaArray, nDM2, deltaM2Array);

     nsigN[i]->SetDirectory(0);
     nsigI[i]->SetDirectory(0);
     nbgN[i]->SetDirectory(0);
     nbgI[i]->SetDirectory(0);
   }
   std::cout<<"Chi2 Histograms created ("<<nTh13<<", "<<nDelta<<", "<<nDM2<<"): Running Grid Method"<<std::endl;

   //So now I have built all the histograms as appropriate 
   // - now I just loop over the data

}


void NueFCSensitivity::GetPoint(int i, float &bg, float &sig)
{
   NSCErrorParam err = nsc->GetErrorConfig(i);  

   float bgt =  currentVal[ClassType::NC]  * (1.0 + err.scale[ClassType::NC])
             + currentVal[ClassType::numu]* (1.0 + err.scale[ClassType::numu])
             + currentVal[ClassType::bnue]* (1.0 + err.scale[ClassType::bnue])
             + currentVal[ClassType::nutau]*(1.0 + err.scale[ClassType::nutau]);

   bg = bgt;
   sig = currentVal[ClassType::nue]*(1.0 + err.scale[ClassType::nue]);
}

void NueFCSensitivity::Initialize()
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

void NueFCSensitivity::LoadEventsFromFile()
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

float NueFCSensitivity::Oscillate(double dm23, double Ue32, double delta, int sign)
{
  if(fMethod == NueSenseConfig::kAllNumbers)
     return OscillateNumber(Ue32);

  if(fMethod == NueSenseConfig::kNueHist)
     return OscillateHist(dm23, Ue32, delta, sign);
   
  if(fMethod == NueSenseConfig::kAnaNueFiles)
     return OscillateFile(dm23, Ue32, delta, sign);

  return -1;
}

float NueFCSensitivity::OscillateNumber(double Ue32)
{
//  double Ue3 = TMath::Sin(t13);

  //Assumed originally Oscillated to:
  //  pow(TMath::Sin(theta23),2)*4.*UE32*(1-UE32)*pow(oscterm,2);
  currentVal[ClassType::nue] = signuenum*Ue32*(1-Ue32);
                                                                                           
  return currentVal[ClassType::nue];
}

float NueFCSensitivity::OscillateHist(double /*dm23*/, double /*th13*/, double /*delta*/, int /*sign*/)
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

float NueFCSensitivity::OscillateFile(double dm23, double Ue32, double delta, int sign)
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
      OscillateMatter(inu,inunoosc,E,dm23,Ue32,delta,sign);

    total[nuClass] += weight*oscweight;
  }

  for(int i = 0; i < 5; i++){
    currentVal[i] = total[i];
  }

  return currentVal[ClassType::nue];
}

float NueFCSensitivity::CalculateChi2(int i, float bg, float sig)
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

float NueFCSensitivity::CalculateFitChi2(int i, float bg, float sig)
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


int NueFCSensitivity::SetupChi2Histograms(NSCDataParam DPst13, NSCDataParam DPdelta, NSCDataParam DPdm23)
{
   chi2n = new TH3D*[fNumConfig];
   chi2i = new TH3D*[fNumConfig];

   nsigN = new TH3D*[fNumConfig];
   nsigI = new TH3D*[fNumConfig];
   nbgN = new TH3D*[fNumConfig];
   nbgI = new TH3D*[fNumConfig];
                                                                                           
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

     chi2n[i]->SetDirectory(0);
     chi2i[i]->SetDirectory(0);

     sprintf(temp, "signormal%d", i);
     nsigN[i] = new TH3D(temp, temp, t13Bins, xstart, xend, delBins,
                             ystart, yend, dm23Bins, zstart, zend);

     sprintf(temp, "siginverted%d", i);
     nsigI[i] = new TH3D(temp, temp, t13Bins, xstart, xend, delBins,
                             ystart, yend, dm23Bins, zstart, zend);

     sprintf(temp, "bgnormal%d", i);
     nbgN[i] = new TH3D(temp, temp, t13Bins, xstart, xend, delBins,
                             ystart, yend, dm23Bins, zstart, zend);
 
     sprintf(temp, "bginverted%d", i); 
     nbgI[i] = new TH3D(temp, temp, t13Bins, xstart, xend, delBins,
                             ystart, yend, dm23Bins, zstart, zend);
   }
   std::cout<<"Chi2 Histograms created ("<<t13Bins<<", "<<delBins<<", "<<dm23Bins<<"): " 
            <<t13Bins*delBins*dm23Bins<<" iterations to perform"<<std::endl;
   return t13Bins*delBins*dm23Bins;
}

void NueFCSensitivity::WriteToFile(std::string file)
{
   TFile out(file.c_str(), "RECREATE");
   out.cd();

   for(int i = 0; i < fNumConfig; i++){
      chi2n[i]->Write();
      chi2i[i]->Write();
      nsigN[i]->Write();
      nbgN[i]->Write();
      nsigI[i]->Write();
      nbgI[i]->Write();
   }
   TTree *info, *err;

   ProduceTree(info, err);
   info->Write();
   err->Write();
   out.Close();
}

double NueFCSensitivity::OscillateMatter(int nuFlavor, int nonOscNuFlavor,
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

double NueFCSensitivity::OscillateMatter(int nuFlavor, int nonOscNuFlavor,
                                    float Energy)
{

  if(nonOscNuFlavor > 0) fOscCalc.SetOscParam(OscPar::kNuAntiNu, 1);
  if(nonOscNuFlavor < 0) fOscCalc.SetOscParam(OscPar::kNuAntiNu, -1);

  return fOscCalc.Oscillate(nuFlavor, nonOscNuFlavor, Energy);
}

void NueFCSensitivity::SetOscParamBase( float dm2, float ss13,
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

void NueFCSensitivity::ProduceTree(TTree * &infotree, TTree * &errtree)
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




