#ifndef NUESYSTEMATIC_H
#define NUESYSTEMATIC_H
#include "NueAna/Extrapolation/Systematic.h"
#include "NueAna/NueAnaTools/Selection.h"
#include "NueAna/NueRecord.h"
#include "MCReweight/NeugenWeightCalculator.h"
#include "MCReweight/SKZPWeightCalculator.h"
#include "MCReweight/Zfluk.h"
#include "MCReweight/Zbeam.h"
#include "TTree.h"
#include "NueAna/Extrapolation/Background.h"
#include "OscProb/OscCalc.h"
 
class NueSystematic
{
  
 public:

  NueSystematic(std::string);
  ~NueSystematic();
  void Clear() {fSysList.clear();}
  
  void AddSystematic(Systematic::Systematic_t sys,Double_t val) {fSysList[sys] = val;}
  Double_t UpdateRecord(NueRecord*,Selection::Selection_t);
  Double_t UpdateRecord(NueRecord*,Selection::Selection_t,Background::Background_t bg);
  Double_t GetAppearanceWeight(NueRecord *rec, Background::Background_t bg);

  void MakeBranches(TTree*);
  Char_t* GetName() {return fName;}

  void SetSKZPParams(std::string cfg);
  std::string GetSKZPParams();
  void SetOscParams(double * par);
  void SetOscParams(Double_t,Double_t,Double_t,
		    Double_t,Double_t,Double_t deltaCP=0,Int_t massH=1);
  void GetOscParams(Double_t &,Double_t &,Double_t &,
		    Double_t &,Double_t &,Double_t &,Int_t &);


  Double_t GetSysValue(Systematic::Systematic_t sys);
  void SetSysValue(Systematic::Systematic_t sys, Double_t val);

  Double_t DoOscCalc(NueRecord *,Double_t);
  Double_t DoNeugenCalc(NueRecord *,Systematic::Systematic_t,Double_t);

 private:
  void Init();

  Double_t DoSKZPCalc(NueRecord *,Double_t);
//  Double_t DoNeugenCalc(NueRecord *,Systematic::Systematic_t,Double_t);
  Double_t DoCalibShift(NueRecord *,Systematic::Systematic_t,Double_t);
  Double_t DoShwDevCalc(NueRecord *,Double_t,Selection::Selection_t);
  Double_t DoTauProd(NueRecord*,Double_t);
  Double_t DoPIDSkew(NueRecord*,Double_t,Selection::Selection_t);
  Double_t DoNormCalc(NueRecord*,Double_t);
  Double_t DoNCScale(NueRecord *record,Double_t val, Selection::Selection_t sel);
  Double_t DoCCShwEnergyScale(NueRecord *record,Double_t val, Selection::Selection_t sel);


  Char_t fName[256];

  //Zbeam fBeamWeighter;
  //Zfluk fHadProdWeighter;  
  //NeugenWeightCalculator *fNWC;
  static SKZPWeightCalculator* skzp;  

  std::map<Systematic::Systematic_t,Double_t> fSysList;
  
//  Double_t  fBeamWeightPars[6];
//  Double_t  fHadProdPars[7];
  std::string skzpcfg;

  Double_t  fTheta23;
  Double_t  fTheta13;
  Double_t  fTheta12;
  Double_t  fDeltaMSq23;
  Double_t  fDeltaMSq12;
  Double_t  fDeltaCP;
  Int_t     fMassHierarchy;
  Double_t  fTempDouble;

  OscCalc fOscCalc;

};
#endif //NUESYSTEMATIC_H
