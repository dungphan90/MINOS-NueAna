////////////////////////////////////////////////////////////////////////
// $Id: NueStandard.h,v 1.58 2016/05/31 21:27:23 wingmc Exp $
//
// The NueStandard namespace is intended to streamline cuts 
//  and other nue related values
//  Big ones include:
//       Data Quality Cuts
//       Fiducial Volume Cuts
//       Preselection Cuts
//
// Author: Josh Boehm 
// Created: October 2, 2007
////////////////////////////////////////////////////////////////////////
#ifndef  NUESTANDARD
#define NUESTANDARD
                                                                                
// This file contains a lot of the standards used in the NueAnalysis
#include "TObject.h"
#include <iostream>
#include "NueRecord.h"
#include "NueAna/NueAnaTools/Selection.h"
#include "OscProb/OscCalc.h"
                                                                                
namespace NueStandard
{
  static const int NRUNPERIODS=2;
  static const double pot_fd[NRUNPERIODS]={0.398,0.602};
  static const double pot_nd[NRUNPERIODS]={0.383,0.617};
  static const double pot_ndmrcc[NRUNPERIODS]={0.393161,0.606839};
  static const double pot_fdmrcc[NRUNPERIODS]={0.398,0.602};
  static const double pot_ndhornoff[NRUNPERIODS]={0.509,0.491};//doc-3429

  static const Double_t kNormalizedNearPOT = 1.0e7; //in units of e12
  static const Double_t kNormalizedFarPOT_Run1 = 1.20963e8;  //1.20933e8; //in units of e12
  static const Double_t kNormalizedFarPOT_Run2 = 1.93577e8; //1.9348e8; //in units of e12
//Real Run3
  static const Double_t kNormalizedFarPOT_Run3 = 3.88063e8; //3.87245e8; //in units of e12
//Run 3Prime
//  static const Double_t kNormalizedFarPOT_Run3 = 5.01668e8; //in units of e12

  static const Double_t kNormalizedFarPOT_Run4 = 0.0810616e8; //0.088e8; //in units of e12
  static const Double_t kNormalizedFarPOT_Run4RHC = 1.6997e8; //in units of e12
  static const Double_t kNormalizedFarPOT_Run5 = 0.449622e8;  //0.450e8; //in units of e12
  static const Double_t kNormalizedFarPOT_Run6 = 0.61379e8;   //0.611e8; //in units of e12
  static const Double_t kNormalizedFarPOT_Run7RHC = 1.23052e8; //in units of e12
  static const Double_t kNormalizedFarPOT_Run8 = 0.122894e8; // in units of e12
  static const Double_t kNormalizedFarPOT_Run9RHC = 0.407265e8; //in units of e12
  static const Double_t kNormalizedFarPOT_Run10prelim = 1.51499e8 + 0.611644e8; // in units of e12
  static const Double_t kNormalizedFarPOT_Run10 = 2.34614e8; // in units of e12

  static const Double_t kNormalizedFarPOT = kNormalizedFarPOT_Run1 + kNormalizedFarPOT_Run2 + kNormalizedFarPOT_Run3; //in units of e12

  static const Double_t kNormalizedFarPOT_PerI = kNormalizedFarPOT_Run1 + kNormalizedFarPOT_Run2 + kNormalizedFarPOT_Run3 + kNormalizedFarPOT_Run4 + kNormalizedFarPOT_Run5 + kNormalizedFarPOT_Run6; // in units of e12

  static const Double_t kNormalizedFarPOT_RHC = kNormalizedFarPOT_Run4RHC + kNormalizedFarPOT_Run7RHC + kNormalizedFarPOT_Run9RHC; //in units of e12

  static const Double_t kNormalizedFarPOT_PerII = kNormalizedFarPOT_RHC; //in units of e12

  static const Double_t kNormalizedFarPOT_PerIII = kNormalizedFarPOT_Run8 + kNormalizedFarPOT_Run10; //in units of e12

//Begin MINOS+ consts here

  static const Double_t kNormalizedFarPOT_Run11 = 2.96697e8;//in units of e12
  static const Double_t kNormalizedFarPOT_Run12 = -5743.8; //Change me into something in e12 units
  static const Double_t kNormalizedFarPOT_Run13 = -5743.8; //Change me into something in e12 units



  static OscCalc fOscGen;

  static const Double_t predwts[5][5] = {{1.00232,1.00311,0.981931,0.981639,1.01361},
					 {0.946369,1.03716,0.96031,0.80668,0.666657},
					 {0.927385,1.00622,0.937925,0.794856,0.614556},
					 {1.20908,1.0726,1.10532,0.938343,1.12485},
					 {0.972292,0.972605,0.972527,0.972485,0.972984}};

  static const Double_t predwts_nuecc[5][5] = {{1.02611,1.02567,1.0064,1.00467,1.03611},
					       {0.98242,0.982655,0.962493,0.96204,0.992857},
					       {0.950664,0.952551,0.9306,0.931326,0.962668},
					       {0.94193,0.945058,0.922595,0.923665,0.954923},
					       {0.972311,0.972081,0.952034,0.951971,0.983061}};
  
  static const Double_t predwts_nc[5][5] = {{0.863239,0.988428,0.931521,0.824146,0.708488},
					    {0.921802,1.10491,1.03858,0.827883,0.636916},
					    {0.956744,1.13551,0.963699,0.798763,0.622532},
					    {0.983896,1.04438,0.925654,0.764158,0.761439},
					    {1.02799,0.951312,0.939436,0.833891,0.899809}};
  
  static const Double_t predwts_numucc[5][5] = {{0.905015,1.06799,0.996018,0.899026,0.777587},
						{0.91316,1.0588,0.974574,0.832987,0.591156},
						{0.90616,0.933719,0.865338,0.66765,0.530119},
						{0.930668,0.938328,0.7928,0.712665,0.584523},
						{0.995441,0.801244,0.845764,0.695091,0.798452}};
  
  static const Double_t predwts_bnuecc[5][5] = {{1.07792,1.02552,1.12685,1.04827,1.15549},
						{1.10753,0.998143,1.34637,0.821027,1.20424},
						{0.821103,0.950048,1.20129,0.757376,1.08648},
						{1.03969,1.18587,0.887482,1.06532,1.07671},
						{1.26007,1.00518,0.972437,1.10064,1.05491}};

  static const Double_t predwts_nutaucc[5][5] = {{0.989861,0.992486,0.990626,0.986435,0.989023},
						 {0.986954,0.986549,0.98681,0.988837,0.988621},
						 {0.981899,0.982569,0.983938,0.982339,0.982898},
						 {0.975608,0.973232,0.977278,0.978017,0.977017},
						 {0.982638,0.984615,0.98071,0.98373,0.987868}};
  
  bool PassesDataQuality(NueRecord *nr);
  bool PassesNearDataQuality(int gbm, float cc, int st);
  bool PassesFarDataQuality(float li, int rc, int dpfddq,
                                 float tmin, double spillt);
  bool PassesFarDataTiming(NueRecord *nr);
    
  bool IsGoodRun(NueRecord *nr);
  bool IsGoodNearRun(NueRecord *nr);
  bool IsGoodFarRun(NueRecord *nr);

  bool PassesPOTStandards(NueRecord *nr);
  bool PassesCosmicCut(NueRecord *nr);
  bool PassesCosmicCutFunction(NueRecord *nr);
  void FillCosmicCut(NueRecord *nr);
  bool IsInFid(NueRecord *nr);

    bool PassesTrackPlaneCut(NueRecord *nr);
    bool PassesTrackPlaneCut(int trkplane);
    bool PassesTrackLikePlaneCut(NueRecord *nr);
    bool PassesTrackLikePlaneCut(int tlp);
    bool PassesLowEnergyCut(NueRecord *nr);
    bool PassesLowEnergyCut(float energy);
    bool PassesHighEnergyCut(NueRecord *nr);
    bool PassesHighEnergyCut(float energy);
    bool IsLargestEvent(NueRecord *nr);

    bool PassesMinPlaneCut(NueRecord *nr);
    bool PassesMinPlaneCut(int planes);
    bool PassesShowerCut(NueRecord *nr);
    bool PassesShowerCut(int nshw);

    bool PassesPreSelection(NueRecord *nr);
    bool PassesPreSelection(int tp, int tlp, float energy);
    bool PassesNonHEPreSelection(NueRecord *nr);
    bool PassesNonHEPreSelection(int tp, int tlp, float energy);
    bool PassesPreSelectionTrackCuts(NueRecord *nr);
    bool PassesPreSelectionTrackCuts(int trkplane, int tlp);
    bool PassesPreSelectionBasicCuts(NueRecord *nr);

    bool PassesSysPreSelectionNoHE(NueRecord *nr);
    bool PassesSysPreSelectionNoHE(int tp, int tlp, float energy);
    bool PassesSysPreSelection(NueRecord *nr);
    bool PassesSysPreSelection(int tp, int tlp, float energy);

    bool PassesMRCCFiducial(NueRecord *nr);
    bool PassesMRCCPreSelection(NueRecord *nr);
    bool PassesMRCCPreSelection(float bestComp, int fitp, double pid);
    
    bool PassesMREFiducial(NueRecord *nr);
    bool PassesMREPreSelection(NueRecord *nr);
    bool PassesMREPreSelection(float bestComp, int fitp, double pid);
 
    bool PassesPIDSelection(NueRecord*nr, Selection::Selection_t sel); 
    bool PassesSelection(NueRecord*nr, Selection::Selection_t sel);
    bool PassesSelection(NueRecord *nr, Selection::Selection_t sel, Selection::Selection_t sel2);

    bool PassesNCCleaningCuts(NueRecord* nr);
    
    bool PassesParticlePIDCut(NueRecord* nr);
    bool PassesParticlePIDPreselectionCut(NueRecord *nr);

    double GetPIDValue(NueRecord *nr, Selection::Selection_t sel);

    bool PassesCCSelection(NueRecord *nr);

    double GetRPWBeamWeight(NueRecord *nr, bool ismrcc=false);
    double GetRPWBeamWeight(std::vector<double> weights, std::vector<double> pots);

    void SetDefaultOscParam();
    void SetDefaultOscParamNoNue();    
    void SetOscParamBestFitANN();
    void SetOscParamBestFitkNN();

    void SetOscNoMatter();
    void GetOscParam(double *par);
    void SetOscParam(double *par);
    void SetOscParam(OscPar::OscPar_t pos, double val);
    double GetOscWeight(int nuFlavor, int nonOsc, double E);
    double GetNSIOscWeight(int nuFlavor,int nonOsc,double E);
    double GetLSNDOscWeight(int nuFlavor,int nonOsc,double E);
    void FillDefaultOscParam(double* par);

    Bool_t IsRun1(NueRecord* nr);
    Bool_t IsRun2(NueRecord* nr);
    Bool_t IsRun3(NueRecord* nr);
    Bool_t IsRun3Prime(NueRecord* nr);
    Bool_t IsRun3NotPrime(NueRecord* nr);
    Bool_t IsRun4(NueRecord* nr);
    Bool_t IsRun4RHC(NueRecord* nr);
    Bool_t IsRun5(NueRecord* nr);
    Bool_t IsRun6(NueRecord* nr);
    Bool_t IsRun7RHC(NueRecord* nr);
    Bool_t IsRun8(NueRecord* nr);
    Bool_t IsRun9RHC(NueRecord* nr);
    Bool_t IsRun10(NueRecord* nr);
    Bool_t IsRun11(NueRecord* nr);
    Bool_t IsRun12(NueRecord* nr);
    Bool_t IsRun13(NueRecord* nr);
    Bool_t IsSpecialRun(NueRecord* nr);

    Double_t GetIntensityBeamWeight(NueRecord *nr);
    Double_t GetSKZPBeamWeight(NueRecord *nr);
    Double_t GetMCWeights(NueRecord *nr);
    Double_t GetNDDataWeights(NueRecord *nr);
    Double_t GetPredWeights_DO_NOT_USE(NueRecord *nr); //For 5 LEM bin x 1 E bin extrapolation
    Double_t GetPredWeights(NueRecord *nr);

    void ModifyANNPID(NueRecord *nr);
    Double_t Calc4thAnaANN(NueRecord *nr, Selection::Selection_t sel);

    void SetE50PID(NueRecord *nr);

    //NSI SetNSI
    void SetNSI(bool nsiflag=true);

    //LSND SetLSND
    void SetLSND(bool lsndflag=true);
}

#endif
