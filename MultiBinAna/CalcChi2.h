#pragma once

#include "TROOT.h"
#include "TMatrixD.h"

class TFile;
class TH1D;
class TH2D;
class TString;

//---------------------------------------------------------------------------------
struct params
{
  double Dm232;
  double Dm221;
  double th23;
  double th12;
  double th13;
  double deltaCP;
  double Dm241;
  double th24;
  double th34;
  double th14;
  double delta14;
  double delta24;
};
//---------------------------------------------------------------------------------
void PrintParStatus(params my_pars);
//---------------------------------------------------------------------------------
class CalcChi2
{
 public:
  CalcChi2();
  void SetDataFit(bool doData) {doDataFit = doData;}
  double Chi2(double Dm232, double Dm221, double th23,
	      double th12, double th13, double deltaCP,
	      double Dm241, double th24, double th34,
	      double th14, double delta14, double delta24);
  
 private:
  void Zombie(TFile* f);
  void LoadInputHistograms(TString fileName);
  void GenerateUnOscillatedSpectra();
  void GenerateOscillatedSpectra(params my_pars);
  TH1D* GetTwoDetSpectrum(TH1D* hND, TH1D* hFD);
  TMatrixD* ScaleCovarianceMatrix(TH1D* pred, TMatrixD* mtx);
  Double_t PenaltyTermDm232(Double_t dm232);
  Double_t PenaltyTermNuisance(Double_t par, Double_t mean, Double_t sigma);
  Double_t ChiSqFunction(TH1D* rPred, TH1D* rData, TMatrixD* CoVarInvert);
  Double_t ComparePredWithData(TH1D* predCC, 
			       TH1D* dataCC,
			       TMatrixD* CoVarCC_inverted, 
			       Double_t Dm2);
  TH1D* CreateTotalSpectrum(params my_pars,
			    TH2D* TrueNC,
			    TH2D* NuMu,
			    TH2D* BeamNue,
			    TH2D* AppNue,
			    TH2D* AppNuTau,
			    double baseline);
  double FourFlavourNuMuToNuSProbability
    (const double energy, 
     double dm232, const double theta23, double dm221, 
     double dm243, const double theta12, 
     const double theta13, const double theta14,
     const double theta24, const double theta34,
     double delta1, double delta2, double delta3,
     const double baseline);
  double FourFlavourDisappearanceWeight
    (const double energy, 
     double dm232, const double theta23, double dm221, 
     double dm243, const double theta12, 
     const double theta13, const double theta14,
     const double theta24, const double theta34,
     double delta1, double delta2, double delta3,
     const double baseline);
  double FourFlavourNuESurvivalProbability
    (const double energy, 
     double dm232, const double theta23, double dm221, 
     double dm243, const double theta12, 
     const double theta13, const double theta14,
     const double theta24, const double theta34,
     double delta1, double delta2, double delta3,
     const double baseline);
  double FourFlavourNuMuToNuEProbability
    (const double energy, 
     double dm232, const double theta23, double dm221, 
     double dm243, const double theta12, 
     const double theta13, const double theta14,
     const double theta24, const double theta34,
     double delta1, double delta2, double delta3,
     const double baseline);
  double FourFlavourNuMuToNuTauProbability
    (const double energy, 
     double dm232, const double theta23, double dm221, 
     double dm243, const double theta12, 
     const double theta13, const double theta14,
     const double theta24, const double theta34,
     double delta1, double delta2, double delta3,
     const double baseline);
  TH1D* CreateSpectrumComponent(params my_pars, TString OscType, TH2D* oscDummy, Double_t baseline);
  
  double Penalty_nuisance;
  double Penalty_dm232;

  double totalChi2_CC;
  double totalChi2;

  TFile *InFile;

  TH1D  *dataCC;
  TH1D  *fakeDataCC;

  bool got_beamOptics;
  bool got_inputHistos;
  bool got_covMx;
  bool got_unOscHistos;
  bool doDataFit;
  
  // Near MC RecoVtrue MINOS
  TH2D *NDCC_TrueNC_minos, *NDCC_NuMu_minos, *NDCC_BeamNue_minos, *NDCC_AppNue_minos, *NDCC_AppNuTau_minos;
  // Far MC RecoVtrue MINOS
  TH2D *FDCC_TrueNC_minos, *FDCC_NuMu_minos, *FDCC_BeamNue_minos, *FDCC_AppNue_minos, *FDCC_AppNuTau_minos;
  // MC (Unoscillated) MINOS
  TH1D *FDUnOscCC_MC_minos, *NDUnOscCC_MC_minos;
  // MC (Oscillated) MINOS
  TH1D *FDOscCC_MC_minos, *NDOscCC_MC_minos;
  // Data MINOS
  TH1D *FD_dataCC_minos, *ND_dataCC_minos;
  // Near MC RecoVtrue MINOS+
  TH2D *NDCC_TrueNC_minosPlus, *NDCC_NuMu_minosPlus, *NDCC_BeamNue_minosPlus, *NDCC_AppNue_minosPlus, *NDCC_AppNuTau_minosPlus;
  // Far MC RecoVtrue MINOS+
  TH2D *FDCC_TrueNC_minosPlus, *FDCC_NuMu_minosPlus, *FDCC_BeamNue_minosPlus, *FDCC_AppNue_minosPlus, *FDCC_AppNuTau_minosPlus;
  // MC (Unoscillated) MINOS+
  TH1D *FDUnOscCC_MC_minosPlus, *NDUnOscCC_MC_minosPlus;
  // MC (Oscillated) MINOS+
  TH1D *FDOscCC_MC_minosPlus, *NDOscCC_MC_minosPlus;
  // Data MINOS+
  TH1D *FD_dataCC_minosPlus, *ND_dataCC_minosPlus;
  // MC (Unoscillated) Joint MINOS/MINOS+
  TH1D *FDUnOscCC_MC, *NDUnOscCC_MC;
  // MC (Oscillated) Joint MINOS/MINOS+
  TH1D *FDOscCC_MC, *NDOscCC_MC;
  // MC (Oscillated & Systematicall fit) Joint MINOS/MINOS+
  TH1D *FDSystOscCC_MC, *NDSystOscCC_MC;
  // Data Joint MINOS/MINOS+
  TH1D *FD_dataCC, *ND_dataCC;
  //Covariance Matrices -- relative variance
  TMatrixD* CoVarCC_relative;
  TMatrixD* CoVarCC_inverted;
};

