#ifndef ErrorCalc_h
#define ErrorCalc_h

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include "NueAna/Extrapolation/Background.h"
#include "NueAna/MultiBinAna/Extrapolate2D.h"
#include "NueAna/MultiBinAna/MultiBinAnaHelper.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMatrixD.h"
#include "TSystem.h"
#include "TTree.h"

using namespace std;

class ErrorCalc {
 public:
  ErrorCalc();
  virtual ~ErrorCalc();

  virtual void AddExtrap(Extrapolate2D *E);
  void AddFNExtrapSyst(string systname, string farplustag, string nearplustag,
                       string farminustag, string nearminustag,
                       string nominaltag, string filetag, int flag = 0,
                       bool useAllRuns = true);
  void SetSysFileDir(string dir) { SysFileDir = dir; };
  void AddSpecialFNSyst(string systname, double err = 0,
                        string histname = "NULL", string filetag = "NULL",
                        int flag = 0);
  void AddCovarianceMatrix(string systname, string file, string histname,
                           int flag = 0);
  void SetNDDataFile(string s) { NDData_infile = s; };

  void CalculateFNExtrapError_SingleBin();
  void CalculateAppearanceExtrapError_SingleBin(Background::Background_t bg);

  void CalculateSystErrorMatrix();

  void CalculateHOOError_SingleBin();
  void CalculateMRCCError_SingleBin();

  virtual void CalculateHOOError();

  void FNRatioError(string ratioerrfilename = "FNRatioErrs.root");
  void NDSystematics(string nderrfilename = "NDSysts.root");
  void FDSystematics();

  TH2D *CovMatrix;
  TH2D *CovMatrix_Decomp;

  void MakeFracError_Lists();
  std::vector<TH1D *> SysBkgd_Plus1Sigma;
  std::vector<TH1D *> SysSig_Plus1Sigma;

  virtual void SetGridPred(int nbins, vector<vector<double> > nc,
                           vector<vector<double> > cc,
                           vector<vector<double> > bnue,
                           vector<vector<double> > tau,
                           vector<vector<double> > sig);
  virtual void SetUseGrid(bool grid = false) { UseGrid = grid; };

  TH2D *GetShifts(Background::Background_t bg, bool plus1sigma, string systname,
                  int extrap);

 protected:
  virtual void Initialize();
  void InitializeSys();
  void InitializeDecompSys();

  void ReadSysFiles_FNExtrap(int n);    // i.e. LisaFiles
  void ReadSysFiles_Appearance(int n);  // i.e. LisaFiles
  virtual void ReadSysFiles_FNExtrap_AllRuns(int n);  // i.e. LisaFiles, combining runs 1,2,&3
  virtual void ReadSysFiles_Appearance_AllRuns(int n);  // i.e. LisaFiles, combining runs 1,2,&3

  virtual void ReadSpecialFiles(int n);
  virtual void ReadCovarianceFiles(int n);

  virtual void
  CalculateSystErrorMatrixGrid();  // function to use for FC when using grid
                                   // files instead of oscillating on the fly
  virtual void CalculateSystErrorMatrixExtrap();  // usual function

  vector<Extrapolate2D *> Extrap;

  MultiBinAnaHelper MBH;

  Bool_t Init;
  Bool_t InitSys;
  Bool_t InitDecompSys;

  string SysFileDir;

  string NDData_infile;

  TH2D *NDCovMatrix;  // covariance matrix from decomposition results

  vector<string> FNExtrap_FarPlusTag;
  vector<string> FNExtrap_NearPlusTag;
  vector<string> FNExtrap_FarMinusTag;
  vector<string> FNExtrap_NearMinusTag;
  vector<string> FNExtrap_StdTag;
  vector<string> FNExtrap_SystName;
  vector<string> FNExtrap_FileTag;
  vector<int> FNExtrap_Flag;  // 0 for ONLY nc,numucc,bnue; 1 for numu
                              // extrapolation contribution to tau/signal, 2 for
                              // overall MC change to tau/signal
  vector<bool> FNExtrap_UseAllRuns;

  vector<string> SpecialSyst_SystName;
  vector<string> SpecialSyst_FileTag;
  vector<string> SpecialSyst_HistName;
  vector<double> SpecialSyst_Value;
  vector<int> SpecialSyst_Flag;  // 0 for all components, 1 for NC+NuMuCC+BNueCC
                                 // only, 2 for NuTauCC only, 3 for NueCC only, 4
                                 // for NuMuCC+NC (fake HOO error)

  vector<string> ExtraCovariance_SystName;
  vector<string> ExtraCovariance_File;
  vector<string> ExtraCovariance_HistName;
  vector<int>
      ExtraCovariance_Flag;  // 0 for all components, 1 for NC+NuMuCC+BNueCC
                             // only, 2 for NuTauCC only, 3 for NueCC only, 4 for
                             // NuMuCC+NC (fake HOO error), 5 for NC only, 6 for
                             // NuMuCC only, 7 for BNueCC only

  std::map<string, vector<TH2D *> > FN_NC_Plus1Sigma;
  std::map<string, vector<TH2D *> > FN_NC_Minus1Sigma;
  std::map<string, vector<TH2D *> > FN_NuMuCC_Plus1Sigma;
  std::map<string, vector<TH2D *> > FN_NuMuCC_Minus1Sigma;
  std::map<string, vector<TH2D *> >
      FN_BNueCC_Plus1Sigma;  // if using f/n method for beam nues
  std::map<string, vector<TH2D *> > FN_BNueCC_Minus1Sigma;

  std::map<string, vector<TH2D *> > FD_NC_Plus1Sigma;
  std::map<string, vector<TH2D *> > ND_NC_Plus1Sigma;
  std::map<string, vector<TH2D *> > FD_NuMuCC_Plus1Sigma;
  std::map<string, vector<TH2D *> > ND_NuMuCC_Plus1Sigma;
  std::map<string, vector<TH2D *> > FD_BNueCC_Plus1Sigma;
  std::map<string, vector<TH2D *> > ND_BNueCC_Plus1Sigma;

  std::map<string, vector<TH2D *> > FD_NC_Minus1Sigma;
  std::map<string, vector<TH2D *> > ND_NC_Minus1Sigma;
  std::map<string, vector<TH2D *> > FD_NuMuCC_Minus1Sigma;
  std::map<string, vector<TH2D *> > ND_NuMuCC_Minus1Sigma;
  std::map<string, vector<TH2D *> > FD_BNueCC_Minus1Sigma;
  std::map<string, vector<TH2D *> > ND_BNueCC_Minus1Sigma;

  std::map<string, vector<TH1D *> > Pred_CC_Fid_Plus1Sigma;
  std::map<string, vector<TH1D *> > Pred_CC_Fid_Minus1Sigma;

  std::map<string, vector<TH2D *> > NuTauCC_MC_Plus1Sigma;
  std::map<string, vector<TH2D *> > NuTauCC_MC_Minus1Sigma;
  std::map<string, vector<TH2D *> > NueCC_MC_Plus1Sigma;
  std::map<string, vector<TH2D *> > NueCC_MC_Minus1Sigma;

  std::map<string, TH2D *> ExtraCovariance;

  std::map<Background::Background_t, vector<TH1D *> > GridPred;

  bool UseGrid;

  //////
  Double_t reg;

};

#endif
