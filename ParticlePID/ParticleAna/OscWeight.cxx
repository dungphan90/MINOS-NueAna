#include "OscWeight.h"
#include "TMath.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std; 


void OscWeight::GetOscParam(double *par)
{
   for(int i = 0; i < int(OscPar::kNumParameters); i++){
     par[i] = fOscPar[i];
   }
}


OscWeight::OscWeight()
{
  for(int i = 0; i < 9; i++) fOscPar[i] = 0;  
  double par[9] = {0};
  SetOscParam(par, true);
}

double OscWeight::Oscillate(int nuFlavor, int nonOscNuFlavor, double Energy)
{
   double L = fOscPar[OscPar::kL];
   double dm2 = fOscPar[OscPar::kDeltaM23];
//   double theta23 = fOscPar[OscPar::kTh23];
                                                                                         
   double oscterm = TMath::Sin(1.269*dm2*L/Energy);
   oscterm = oscterm*oscterm;
                                                                                         
   double sin2th23 = fSin2Th[OscPar::kTh23];
   double sinsq2th23 = sin2th23*sin2th23;
   double sinth23 = fSinTh[OscPar::kTh23];
   double costh23 = fCosTh[OscPar::kTh23];
   double UE32 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
                                                                                         
   double pmt=(1-UE32)*(1-UE32)*sinsq2th23*oscterm;
   double pme=sinth23*sinth23*4.*UE32*(1-UE32)*oscterm;
   double pmm=1.-pmt-pme;
                                                                                         
   double pet=4*(1-UE32)*UE32*oscterm*costh23*costh23;
   double pem=sinth23*sinth23*4.*UE32*(1-UE32)*oscterm;
   double pee=1.-pet-pem;
                                                                                         
                                                                                         
   if(abs(nonOscNuFlavor)==14){
      if(abs(nuFlavor)==12)       { return pme; }
      else if(abs(nuFlavor)==14)  { return pmm; }
      else if(abs(nuFlavor)==16)  { return pmt; }
   }
   else if(abs(nonOscNuFlavor)==12){
      if(abs(nuFlavor)==12)       { return pee; }
      else if(abs(nuFlavor)==14)  { return pem; }
      else if(abs(nuFlavor)==16)  { return pet; }
   }
   else{
     std::cout<<"I don't know what to do with "<<nonOscNuFlavor
         <<" "<<nuFlavor<<" "<<pee<<std::endl;
   }
   return 0.;
}
   
double OscWeight::Oscillate(int nuFlavor, int nonOscNuFlavor, float Energy,
                              float L, float dm2, float theta23, float UE32)
{
   double par[9] = {0};
   for(int i = 0; i < 9; i++) {par[i] = fOscPar[i];}
   par[OscPar::kL] = L;
   par[OscPar::kDeltaM23] = dm2;
   par[OscPar::kTh23] = theta23;

   SetOscParam(par);

   double oscterm = TMath::Sin(1.269*dm2*L/Energy);
   oscterm = oscterm*oscterm;

   double sin2th23 = fSin2Th[OscPar::kTh23];
   double sinsq2th23 = sin2th23*sin2th23;
   double sinth23 = fSinTh[OscPar::kTh23];
   double costh23 = fCosTh[OscPar::kTh23];

   double pmt=(1-UE32)*(1-UE32)*sinsq2th23*oscterm;
   double pme=sinth23*sinth23*4.*UE32*(1-UE32)*oscterm;
   double pmm=1.-pmt-pme;

   double pet=4*(1-UE32)*UE32*oscterm*costh23*costh23;
   double pem=sinth23*sinth23*4.*UE32*(1-UE32)*oscterm;
   double pee=1.-pet-pem;


   if(abs(nonOscNuFlavor)==14){
      if(abs(nuFlavor)==12)       { return pme; }
      else if(abs(nuFlavor)==14)  { return pmm; }
      else if(abs(nuFlavor)==16)  { return pmt; }
   }
   else if(abs(nonOscNuFlavor)==12){
      if(abs(nuFlavor)==12)       { return pee; }
      else if(abs(nuFlavor)==14)  { return pem; }
      else if(abs(nuFlavor)==16)  { return pet; }
   }
   else{
     std::cout<<"I don't know what to do with "<<nonOscNuFlavor
         <<" "<<nuFlavor<<" "<<pee<<std::endl;
   }
   return 0.;
}



double OscWeight::OscillateMatter(int nuFlavor, int nonOscNuFlavor,
                                    double Energy, double * par)
{
  SetOscParam(par);
  return OscillateMatter(nuFlavor, nonOscNuFlavor, Energy); 
}


double OscWeight::OscillateMatter(int nuFlavor, int nonOscNuFlavor,
                                    double Energy)
{
  double x[1] = {};
  x[0] = Energy;

  double L = fOscPar[OscPar::kL];
  double cos_th23 = fCosTh[OscPar::kTh23];
  double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
  double dm2 = fOscPar[OscPar::kDeltaM23];

  if(abs(nonOscNuFlavor)==14) {
    if(abs(nuFlavor)==12)      return OscWeight::ElecAppear(x);
    else if(abs(nuFlavor)==14) return OscWeight::MuSurvive(x);
    else if(abs(nuFlavor)==16) return OscWeight::MuToTau(x);
  }
  if(abs(nonOscNuFlavor)==12) {
    if(abs(nuFlavor)==14) return OscWeight::ElecAppear(x);

    double sin_dm = sin(1.27*dm2*L/x[0]);
    double eToTau = sinsq_2th13*cos_th23*cos_th23*sin_dm*sin_dm;

    if(abs(nuFlavor)==12) return 1 - OscWeight::ElecAppear(x) - eToTau;
    else if(abs(nuFlavor)==16) return (eToTau);
  }

  return 0;
}

void OscWeight::SetOscParam(OscPar::OscPar_t pos, double val)
{
  //cout<<"setting  "<<pos<<"  "<<val<<endl;
  if(fabs(fOscPar[pos] - val) > 1e-9){
    fOscPar[pos] = val;
    if(pos < 3){
       fSinTh[pos] = sin(val);
       fCosTh[pos] = cos(val);
       fSin2Th[pos] = sin(2*val);
    }
    if(pos == OscPar::kDensity){
      double ne = OscPar::z_a*OscPar::A_av*val; //electron density #/cm^{3}
      fElecDensity = ne*OscPar::invCmToeV*OscPar::invCmToGeV*OscPar::invCmToGeV;
      //electron density with units Gev^{2} eV
      //Gev^{2} to cancel with GeV^{-2} in Gf
                                                                                
       fA = OscPar::root2*OscPar::Gf*fElecDensity; //eV
       fSinAL = sin(fA*fOscPar[OscPar::kL]/(OscPar::invKmToeV*2.));
    }
  }
}

void OscWeight::SetOscParam(double *par, bool force)
{
   for(int i = 0; i < int(OscPar::kNumParameters); i++){
     if(fOscPar[i] != par[i] || force){
        fOscPar[i] = par[i];
        if(i < 3){
            fSinTh[i] = sin(par[i]);
            fCosTh[i] = cos(par[i]);
            fSin2Th[i] = sin(2*par[i]);
        }
        if(i == OscPar::kDensity){
          double ne = OscPar::z_a*OscPar::A_av*par[i]; //electron density #/cm^{3}
          fElecDensity = ne*OscPar::invCmToeV*OscPar::invCmToGeV*OscPar::invCmToGeV;
          //electron density with units Gev^{2} eV
          //Gev^{2} to cancel with GeV^{-2} in Gf
          
          fA = OscPar::root2*OscPar::Gf*fElecDensity; //eV
          fSinAL = sin(fA*par[OscPar::kL]/(OscPar::invKmToeV*2.));
        }
     }
   }
}

double OscWeight::ElecAppear(double *x)
{

  //x[0] = E

  double E = x[0]; //energy
  double L = fOscPar[OscPar::kL];  //baseline
  double plusminus = int(fOscPar[OscPar::kNuAntiNu]);

  double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
  double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];                   

  double cos_th23 = fCosTh[OscPar::kTh23];
  double sin_th23 = fSinTh[OscPar::kTh23];
  double cos_th13 = fCosTh[OscPar::kTh13];

  double sin_2th23 = fSin2Th[OscPar::kTh23];
  double sin_2th13 = fSin2Th[OscPar::kTh13];
  double sin_2th12 = fSin2Th[OscPar::kTh12];

  double d_cp = fOscPar[OscPar::kDelta];

  double dmsq_23 = fOscPar[OscPar::kDeltaM23];
  double dmsq_12 = fOscPar[OscPar::kDeltaM12];
  double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}

  //double Delta23 = dmsq_23/(2.*E*1e9); //eV
  double Delta12 = dmsq_12/(2.*E*1e9);
  double Delta13 = dmsq_13/(2.*E*1e9);

  double B = TMath::Abs(fA - plusminus*Delta13); //eV
  double J = cos_th13*sin_2th12*sin_2th13*sin_2th23;

  double sin_BL = sin(B*L/(OscPar::invKmToeV*2.));
  double sin_AL = fSinAL;


  double p1 = sin_th23*sin_th23*sinsq_2th13*(Delta13/B * Delta13/B)
    *sin_BL*sin_BL;

  double p2 = 0;
  double p3 = 0;
  if(fA !=0){
    p2 = cos_th23*cos_th23*sinsq_2th12*(Delta12/fA * Delta12/fA)
      *sin_AL*sin_AL;

    p3 = J*Delta12*Delta13*sin_AL*sin_BL
      *cos(plusminus*d_cp + Delta13*L/(OscPar::invKmToeV*2.))/(fA*B);
  }

  //std::cout<<dmsq_13<<"  "<<dmsq_12<<"  "<<dmsq_23<<"  "<<Delta13<<"  "
    //       <<fA<<"  "<<B<<std::endl;


  if(p1+p2+p3>1) return 1;
  return p1+p2+p3;
}

double OscWeight::MuToTau(double *x)
{

  //x[0] = E

  double E = x[0]; //energy
  double L = fOscPar[OscPar::kL];  //baseline

  double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];

  double cos_th13 = fCosTh[OscPar::kTh13];

  double dmsq_23 = fOscPar[OscPar::kDeltaM23];

  double Delta23 = dmsq_23/(2.*E*1e9); //eV
  double sin_delta = sin(Delta23*L/(OscPar::invKmToeV*2.));  
 
  double p1 = sinsq_2th23*cos_th13*cos_th13*cos_th13*cos_th13
    *sin_delta*sin_delta;  //numu->nutau


  return p1;

}


double OscWeight::MuSurvive(double *x)
{
  double p1 = 1. - OscWeight::MuToTau(x) - OscWeight::ElecAppear(x);
  if(p1 < 0) p1 = 0;
  return p1;
}

