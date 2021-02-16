#ifndef NueSterileFit_h
#define NueSterileFit_h

#include "NueAna/MultiBinAna/ErrorCalc.h"
#include "NueAna/MultiBinAna/Extrapolate2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TMinuit.h"
#include "TObject.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TSystem.h"

#include "TRandom3.h"
#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "TDecompQRH.h"

#ifndef ROOT_TVectorD
#include "TVectorD.h"
#endif

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <time.h>

using namespace std;
using namespace TMath;


class MultiVarGauss {
public:
  MultiVarGauss(TH1 *hMean, TH2 *hCov);
  MultiVarGauss(TH1 *hMean, TMatrixD *hCov);
  ~MultiVarGauss();

  TH1 *Generate();

private:
  void DecomposeCov();

  TRandom3 *_generator;
  Int_t _npar;
  TVectorD _bestValues;
  TH1 *_hVaried;

  TMatrixDSym _cov;
  Double_t _det;
  TMatrixD _UT;
};

class NueSterileFit : public TObject {
public:
  NueSterileFit();
  virtual ~NueSterileFit();

  virtual void Run2DSterileSlice();
  virtual void Run3FlavorDeltaCPFit();

  virtual void SetIncludeOscParErrs() {IncludeOscParErrs = true;}
  virtual void ReadGridFiles();
  virtual void SetGridFiles(std::string snorm, std::string sinvt);
  void SetGridNorm(double n = 1.) {GridNorm = n;}

  void SetNExperiments(unsigned int n = 1000) {NumExpts = n;}
  virtual void GenerateOneCorrelatedExp(TH1D *nexp, std::vector<TMatrixD>& Matrices);
  virtual void RunMultiBinPseudoExpts(bool Print);
  std::vector<TMatrixD> GenerateDiagonalizedCovarianceMatrix(TH1D *nexp, TH2D *err);
  void PrintMatrix(TMatrixD& matrix);
  double BinLikeComparisonSterileNeutrinoFitFeldmanCousins(std::vector<double> npar);
  double DoBinMinParamSterileNeutrinoFitFeldmanCousins();
  double FindBestFitForPseudoExperiment();

  void AddError(ErrorCalc *Err); // Set the errorcalc object.
  virtual void AddExtrap(Extrapolate2D *E);
  void SetNObs(TH1D *n); // x axis is 2D bin number, content is number of events in that bin

  void SetOutputFile(std::string s = "SensOut.root");

  void SetNSinSqTh14Steps(int n = 100) { nSinSqTh14Steps = n; };
  void SetNSinSqTh24Steps(int n = 100) { nSinSqTh24Steps = n; };
  void SetSinSqTh14Range(double l = 0.0, double h = 1.0);
  void SetSinSqTh24Range(double l = 0.0, double h = 1.0);

  void SetDeltaCP13(double dcp_low, double dcp_high, unsigned int nSteps) {
    DeltaCP13Low  = dcp_low;
    DeltaCP13High = dcp_high;
    nDeltaCP13Steps = nSteps;
  }

  void SetFitMethod(int m = 0) {
    // 0 = poisson, 1 = scaled chi2, 2 = standard chi2, 3 = standard (nuisance
    // param for each syst) likelihood, 4 = bin by bin likelihood

    FitMethod = m;
  }

  void DefTheta13(double t13, double t13unc) {
    Theta13 = t13;
    Theta13Unc = t13unc;
    return;
  }

  void DefTheta23(double t23, double t23unc) {
    Theta23 = t23;
    Theta23Unc = t23unc;
    return;
  }

  void DefTheta12(double t12, double t12unc) {
    Theta12 = t12;
    Theta12Unc = t12unc;
    return;
  }

  void DefDMS32(double dms32, double dms32unc) {
    DeltaMSq32 = dms32;
    DeltaMSq32Unc = dms32unc;
    return;
  }

  void DefDMS21(double dms21, double dms21unc) {
    DeltaMSq21 = dms21;
    DeltaMSq21Unc = dms21unc;
    return;
  }

  void DefDeltaCP(double deltaCP) {
    Delta = deltaCP;
    return;
  }

  void DefDMS41(double dms41) {
    DeltaMSq41 = dms41;
    return;
  }

  double StdLikeComparison(std::vector<double> npar);
  double BinLikeComparison(std::vector<double> npar);

private:
  Double_t Theta12;
  Double_t Theta12Unc;
  Double_t Theta23;
  Double_t Theta23Unc;
  Double_t Theta13;
  Double_t Theta13Unc;
  Double_t DeltaMSq32;
  Double_t DeltaMSq32Unc;
  Double_t DeltaMSq21;
  Double_t DeltaMSq21Unc;
  Double_t Delta;

  double DeltaMSq41;

  ErrorCalc *ErrCalc;
  TH2D *ErrorMatrix;
  TH2D *InvErrorMatrix;
  TH2D *ExternalErrorMatrix;
  std::vector<TH1D*> FracErr_Bkgd_List;
  std::vector<TH1D*> FracErr_Sig_List;
  TH1D *FracErr_Bkgd;
  TH1D *FracErr_Sig;

  std::vector<Extrapolate2D*> Extrap;

  TH1D *NObs;
  TH1D *Bkgd;
  TH1D *Sig;
  TH1D *NExp;
  unsigned int nBins;

  TMinuit *minuit;
  int FitMethod;

  int nDeltaCP13Steps;
  double DeltaCP13Low;
  double DeltaCP13High;

  int nSinSqTh14Steps;
  int nSinSqTh24Steps;
  double SinSqTh14Low, SinSqTh14High;
  double SinSqTh24Low, SinSqTh24High;

  virtual void DefineStdDlnLMinuit();
  virtual void DefineBinDlnLMinuit();

  double PoissonChi2(TH1D *nexp);
  double ScaledChi2(TH1D *nexp_bkgd, TH1D *nexp_signal);

  double StandardChi2(TH1D *nexp);
  void CalculateErrorMatrixInv(TH1D *nexp);

  double StandardLikelihood();
  double DoStdMinParam();

  double BinLikelihood();
  double DoBinMinParam();
  virtual void CalculateDlnLMatrix(TH2D *SystMatrix, TH2D *HOOHEMatrix);

  virtual bool VariateParameter(std::vector<OscPar::EOscPar> ParameterEnumToVariates, std::vector<std::pair<double, double> > ParameterValueToVariates, std::vector<unsigned int> numberOfGridPoint, int idx);

  std::string outFileName;

  bool IncludeOscParErrs;
  std::string GridFileName_Normal;
  std::string GridFileName_Inverted;
  std::vector<TTree*> GridTree_Normal;
  std::vector<TTree*> GridTree_Inverted;
  std::vector<std::vector<TTree*> > GridTree_2_Normal;
  std::vector<std::vector<TTree*> > GridTree_2_Inverted;
  TTree *paramtree_Normal;
  TTree *paramtree_Inverted;
  int nPts_Normal;
  int nPts_Inverted;
  double GridNorm;
  double GridScale_Normal;
  double GridScale_Inverted;

  unsigned int NumExpts;

  double grid_background;
  double grid_signal;
  double grid_sinsqth14;
  double grid_sinsqth24;
  double grid_dmsq41;
  double grid_delta13;
  double grid_delta14;
  double grid_delta24;
  double grid_th34;
  std::vector<double> grid_bin_oscparerr;
  double grid_nc;
  double grid_numucc;
  double grid_bnuecc;
  double grid_nutaucc;
  double grid_nue;
  double grid_normal_th12;
  double grid_normal_th13;
  double grid_normal_th23;
  double grid_normal_dm2_32;
  double grid_normal_dm2_21;
  double grid_inverted_th12;
  double grid_inverted_th13;
  double grid_inverted_th23;
  double grid_inverted_dm2_32;
  double grid_inverted_dm2_21;

  bool FixTheta14Theta24SterileNeutrinoFitFeldmanCousins;

  ClassDef(NueSterileFit, 1)
};

#endif
