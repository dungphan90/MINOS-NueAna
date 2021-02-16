#ifndef GridGenSterile_h
#define GridGenSterile_h

#include "NueAna/Extrapolation/Background.h"
#include "NueAna/MultiBinAna/Extrapolate2D.h"
#include "OscProb/OscCalc.h"
#include "TFile.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TTree.h"
#include <vector>
#include <string>
#include <time.h>

class GridGenSterile {
public:
  GridGenSterile();
  virtual ~GridGenSterile();

  void SetOutputFile(string s = "GridOut.root") { outFileName = s; }
  void SetNExperiments(int n = 2000) { NumExpts = n; };
  void AddExtrap(Extrapolate2D *E);

  void SetNormalHierarchy(Bool_t normal = true) { NormalHier = normal; }
  void SetTheta12(double val = 0.58725, double errup = 0.01409, double errdn = 0.01409);
  void SetTheta13(double val = 0.14819, double errup = 0.00240, double errdn = 0.00240);
  void SetTheta23(double val = 0.82143, double errup = 0.02306, double errdn = 0.02807);
  void SetAbsValDeltaMSq23(double val = 0.00255, double errup = 0.00004, double errdn = 0.00004);
  void SetDeltaMSq12(double val = 7.53E-5, double errup = 0.18E-5, double errdn = 0.18E-5);

  void SetSinSqTh14Grid(double ssqth14) {SinSqTh14Grid = ssqth14;}
  void SetSinSqTh24Grid(double ssqth24) {SinSqTh24Grid = ssqth24;}
  void SetDMSQ41Grid(double dmsq41) {DMSQ41Grid = dmsq41;}

  virtual void RunMultiBinOscParErrsSterileFit(std::string s = "OscParErrDistributions.root");

protected:
  double AsymGaus(double sigma_minus, double sigma_plus);

private:
  std::string outFileName;

  std::vector<Extrapolate2D *> Extrap;

  // These PDG 3-flavor values will remain fixed
  Double_t Theta12;
  Double_t Theta13;
  Double_t Theta23;
  Double_t DeltaMSq23;
  Double_t DeltaMSq12;
  Double_t dTheta13_up, dTheta13_dn;
  Double_t dTheta12_up, dTheta12_dn;
  Double_t dTheta23_up, dTheta23_dn;
  Double_t dDeltaMSq23_up, dDeltaMSq23_dn;
  Double_t dDeltaMSq12_up, dDeltaMSq12_dn;

  double SinSqTh14Grid;
  double SinSqTh24Grid;
  double DMSQ41Grid;

  Bool_t NormalHier;

  Int_t NumExpts;
};

#endif