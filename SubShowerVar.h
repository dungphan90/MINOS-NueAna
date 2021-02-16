///////////////////////////////////////////////////////////////////////////
// 
// SubShowerVar: variables from Chris and Hai's SubShower/Shower package
//
//
///////////////////////////////////////////////////////////////////////////
#ifndef SUBSHOWERVAR_H
#define SUBSHOWERVAR_H

#include "TObject.h"

class SubShowerVar : public TObject
{

 public:
  SubShowerVar();
  virtual ~SubShowerVar();

  //virtual void Draw(Option_t *option);
  //virtual void Print(Option_t *option) const;
  void Zero();
  void Reset();

  //SubShowerVar variables
  Int_t ncluster;
  Int_t nclusterU;
  Int_t nclusterV;
  Int_t nPhysClusterU;
  Int_t nPhysClusterV;
  Int_t nstp0U;
  Int_t nstp0V;
  Double_t E2to1U;
  Double_t E2to1V;
  Double_t PHAvgIDU;
  Double_t PHAvgIDV;
  Double_t PHAvgProbEMU;
  Double_t PHAvgProbEMV;
  Double_t PHFracRMSU;
  Double_t PHFracRMSV;
  Double_t PHAvgDevU;
  Double_t PHAvgDevV;
  Double_t pid;

  Int_t nclusterU_he;    //number of U clusters with p>=0.5 GeV
  Int_t nclusterV_he;    //number of V clusters with p>=0.5 GeV

  //Numbers  
  Int_t nEMU;
  Int_t nEMV;
  Int_t nHadU;
  Int_t nHadV;
  Int_t nTrkU;
  Int_t nTrkV;
  Int_t nXTkU;
  Int_t nXTkV;
  Int_t nHalU;
  Int_t nHalV;

  //Numbers  
  Int_t nEMU_he;
  Int_t nEMV_he;
  Int_t nHadU_he;
  Int_t nHadV_he;
  Int_t nTrkU_he;
  Int_t nTrkV_he;
  Int_t nXTkU_he;
  Int_t nXTkV_he;
  Int_t nHalU_he;
  Int_t nHalV_he;

  //total energy:
  Double_t eEMU;
  Double_t eEMV;
  Double_t eHadU;
  Double_t eHadV;
  Double_t eTrkU;
  Double_t eTrkV;
  Double_t eXTkU;
  Double_t eXTkV;
  Double_t eHalU;
  Double_t eHalV;

  //effective number:
  Double_t wEMU;
  Double_t wEMV;
  Double_t wHadU;
  Double_t wHadV;
  Double_t wTrkU;
  Double_t wTrkV;
  Double_t wXTkU;
  Double_t wXTkV;
  Double_t wHalU;
  Double_t wHalV;

 private:

  ClassDef(SubShowerVar,7)
};

#endif// SUBSHOWERVAR_H
