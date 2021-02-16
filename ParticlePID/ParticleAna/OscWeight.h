#ifndef OSCWEIGHT_H
#define OSCWEIGHT_H

/*
namespace OscPar
{
   typedef enum EOscPar{
     kTh12 = 0,
     kTh23 = 1,
     kTh13 = 2,
     kDeltaM23 = 3,
     kDeltaM12 = 4,
     kDelta = 5,
     kDensity = 6,
     kL = 7,
     kNuAntiNu = 8,
     kUnknown = 9
   } OscPar_t;


#if !defined(__CINT__) || defined(__MAKECINT__)
   // Numbers pulled from 2006 PDG pg 97
   static const double z_a = 0.5; //average Z/A
   static const double A_av = 6.0221415e23; //avogadro's number
   static const double invCmToeV = 1.97326968e-5; //convert 1/cm to eV
   static const double invCmToGeV = 1.97326968e-14; //convert 1/cm to GeV
   static const double invKmToeV = 1.97326968e-10; //convert 1/km to eV
   static const double Gf = 1.166371e-5; //fermi constant (GeV^{-2})
   static const double root2 = 1.41421356;
#endif

}
*/
 // the above is in OscProb/OscCalc.h ... use that so theres only one defintion
#include "OscProb/OscCalc.h"


class OscWeight
{
  public:

   OscWeight();

   void SetOscParam(double * par, bool force = false);  
   void SetOscParam(OscPar::OscPar_t pos, double val); 

	void GetOscParam(double *par);

   double Oscillate(int nuFlavor, int nonOscNuFlavor, float Energy,
                      float L, float dm2, float theta23, float UE32);

   double Oscillate(int nuFlavor, int nonOscNuFlavor, double Energy);

   double OscillateMatter(int nuFlavor, int nonOscNuFlavor, double Energy, double * par);
   double OscillateMatter(int nuFlavor, int nonOscNuFlavor, double Energy);

   double ElecAppear(double *x);
   double MuToTau(double *x);
   double MuSurvive(double *x);

  private:

   double fOscPar[10];

   double fSinTh[3];
   double fCosTh[3];
   double fSin2Th[3];
   double fElecDensity;
   double fA;
   double fSinAL;

};


#endif


