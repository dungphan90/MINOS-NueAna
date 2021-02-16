#ifndef Extrapolate2D_h
#define Extrapolate2D_h

#include "TROOT.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TMath.h"
#include "NueAna/MultiBinAna/MultiBinAnaHelper.h"
#include <vector>
#include "OscProb/OscCalc.h"
#include "OscProb/NuOscProbCalc.h"
#include "NueAna/Extrapolation/Background.h"
#include "TH3D.h"
#include <map>
#include "TFile.h"
#include <iostream>
#include "TTree.h"
#include "TSystem.h"

using namespace std;

class Extrapolate2D{
  public:
    Extrapolate2D();
    virtual ~Extrapolate2D();

    typedef enum OscType {
      qNC = 0,
      qNuMuToNuMu = 1,
      qBNueToNuMu = 2,
      qBNueToBNue = 3,
      qNuMuToNue = 4,
      qNuMuToNuTau = 5,
      qBNueToNuTau = 6,
      qNuMuNC = 7,
      qBNueNC = 8
    } OscType_t;

    void GetPrediction();//Gets the prediction at the nominal values (from the input file) of the oscillation parameters.  To get the prediction with different parameter values, call this, then Oscillate Prediction().
    virtual void OscillatePrediction();//this only changes the 2D prediction histograms (not the 3D ones!)

    void SetRecoBins(Int_t n,Double_t *e);
    void SetPIDBins(Int_t n,Double_t *e);
    void SetTrueBins(Int_t n,Double_t *e);
    void SetNearPOT(double pot=1e19) { nPOTNear = pot;};
    void SetFarPOT(double pot=7e20)  { nPOTFar = pot;};
    void SetOutputFile(string filename);//will write the predictions with the nominal values (from the input file) of the oscillation parameters
    void SetPID(string pid="ANN11");
    void SetNDDataPID(string pid="ANN11");//in case the ND Data input files have a different pid identifier than the FN input files
    void SetFNforBeamNue() { FNforBeamNue = true; };
    void SetRunPeriod(int run=1) { RunPeriod = run; };
    void SetReadExtrapFromFile(string s);
    void SetPrintResult() { PrintResult = true; };

    void SetOscPar(OscPar::OscPar_t par, double val);
    void SetSinSq2Th13(double val);
    void SetSinSq2Th14(double val);
    void SetSinSqTh14(double val);
    void SetSinSqTh24(double val);
    void SetDeltaCP(double val);
    void InvertMassHierarchy();

    void SetNDDataFile(string s);
    void SetFNFile(string s);
    void SetMREFile(string s);
    void SetXSecFile(string s);

    int GetNPID() { return nPID; };
    int GetNReco() { return nReco; };
    int GetNTrue() { return nTrue; };
    double GetNearPOT() { return nPOTNear; };
    double GetFarPOT() { return nPOTFar; };
    string GetPID() { return PID; };
    string GetNDDataPID() { return PID_NDData; };
    bool GetFNforBeamNue() { return FNforBeamNue; };
    int GetRunPeriod() { return RunPeriod; };
    TH2D* GetPredHist(Background::Background_t bg) { return Pred[bg]; };
    bool GetOscFlag() { return Oscillated; };

    std::map<OscType_t, TH3D*> Pred_3D;
    std::map<OscType_t, TH3D*> Pred_Nu_3D;
    std::map<OscType_t, TH3D*> Pred_NuBar_3D;
    std::map<Background::Background_t, TH2D*> Pred;
    TH2D *Pred_NueCCSignal;
    TH2D *Pred_NueBarCCSignal;
    TH2D *Pred_TotalBkgd;
    TH1D *Pred_CC_Fid;
    TH1D *Pred_NuMuCC_Fid;
    TH1D *Pred_NuMuBarCC_Fid;

    TH1D *Pred_TotalBkgd_VsBinNumber;//for input to NueFit2D
    TH1D *Pred_Signal_VsBinNumber;//for input to NueFit2D

    std::map<Background::Background_t, TH2D *> FNRatio;

    //ND data input hists
    std::map<Background::Background_t,TH2D*> NDData;
    TH2D *NDData_Total;
    TH1D *NDData_Reco_CClike;
    TH1D *NDData_Reco_CClike_Pos;
    TH1D *NDData_Reco_CClike_Neg;

    std::map<Background::Background_t,TH2D*> NDMC;
    std::map<Background::Background_t,TH2D*> FDMC;
    //NSI
    void SetEps_ee(double val);
    void SetEps_emu(double val);
    void SetEps_etau(double val);
    void SetEps_mumu(double val);
    void SetEps_mutau(double val);
    void SetEps_tautau(double val);

    void SetDelta_emu(double val);
    void SetDelta_etau(double val);
    void SetDelta_mutau(double val);

    void SetTheta14(double val);
    void SetTheta24(double val);
    void SetTheta34(double val);
    void SetDm41(double val);
    void SetDelta14(double val);
    void SetDelta24(double val);

   void SetOscMethod(int m=0){ OscMethod = m; }; // 0 = StdNueOscCalc, 1 = osc.Oscillate, 2 = osc.OscillateNSI, 3 = osc.OscillateLSND


  protected:

    virtual void ReadNDDataFile();
    virtual void ReadFNFile();
    void ReadMREFile();
    void ReadXSecFile();
    virtual void ReadFiles();
    void ReadExtrap();

    void InitializeOscCalc();
    virtual void SetupPredHists();
    void GetNoOscCCPrediction();
    virtual void FarNearPred(Background::Background_t bg);
    virtual void AppearancePred(Background::Background_t bg);
    virtual void BeamNuePred();
    void SetNominalOscProb();

    void Set1DPredHists();//takes 2D Pred histograms and makes 1D (content vs 2D bin number) histograms.  useful for NueFit2D and GridGen.

    virtual void RebinInputHists();

    void WriteToFile();

    void OscillatePrediction_SepNuNuBar();
    void OscillatePrediction_NoSep();

    string PID;
    string PID_NDData;

    string outFileName;

    int nReco;
    int nTrue;
    int nPID;
    int nRecoCC;

    vector<double> RecoEdges;
    vector<double> PIDEdges;
    vector<double> TrueEdges;
    vector<double> RecoCCEdges;

    Bool_t RebinE;
    Bool_t RebinP;

    Double_t nPOTNear,nPOTFar;

    vector<Background::Background_t> FDComponents;

    Bool_t FNforBeamNue;

    Int_t RunPeriod;

    MultiBinAnaHelper MBH;

  public:
    OscCalc osc;

  protected:
    //input files
    string NDData_infile;
    string FN_infile;
    string MRE_infile;
    string XSec_infile;

    bool ReadMREFromFile;

    //input hists, FD
    map<OscType_t,TH3D*> FD_TrueVsRecoVsPID;
    map<OscType_t,TH3D*> FD_Nu_TrueVsRecoVsPID;
    map<OscType_t,TH3D*> FD_NuBar_TrueVsRecoVsPID;

    TH2D *FD_TrueVsReco_Fid;
    TH2D *FD_TrueVsReco_Fid_NuMuCC;
    TH2D *FD_TrueVsReco_Fid_NuMuBarCC;

    TH2D *FD_TrueVsReco_CClike;
    TH2D *FD_TrueVsReco_CClike_Pos;
    TH2D *FD_TrueVsReco_CClike_Neg;

    TH2D *FD_TrueVsReco_CClike_NuMuCC;
    TH2D *FD_TrueVsReco_CClike_Pos_NuMuBarCC;
    TH2D *FD_TrueVsReco_CClike_Neg_NuMuCC;

    TH2D *FD_TrueVsReco_Fid_NueCC;
    TH2D *FD_TrueVsReco_Fid_NuTauCC;
    TH1D *FD_True_Fid_NueCC;
    TH1D *FD_True_Fid_NuTauCC;

    TH2D *FD_TrueVsReco_Fid_NueBarCC;
    TH2D *FD_TrueVsReco_Fid_NuTauBarCC;
    TH1D *FD_True_Fid_NueBarCC;
    TH1D *FD_True_Fid_NuTauBarCC;

    //input hists, ND
    map<OscType_t,TH3D*> ND_TrueVsRecoVsPID;

    TH2D *ND_TrueVsReco_Fid;
    TH2D *ND_TrueVsReco_Fid_NuMuCC;
    TH2D *ND_TrueVsReco_Fid_NuMuBarCC;

    TH2D *ND_TrueVsReco_CClike;
    TH2D *ND_TrueVsReco_CClike_Pos;
    TH2D *ND_TrueVsReco_CClike_Neg;

    TH2D *ND_TrueVsReco_CClike_NuMuCC;
    TH2D *ND_TrueVsReco_CClike_Pos_NuMuBarCC;
    TH2D *ND_TrueVsReco_CClike_Neg_NuMuCC;

    //derived hists for extrapolation
    std::map<Background::Background_t,TH2D*> ND_DataOverMC_RecoVsPID;

    TH1D *ND_Reco_CClike;
    TH1D *ND_Reco_CClike_Pos;
    TH1D *ND_Reco_CClike_Neg;

    TH1D *FD_Reco_CClike;
    TH1D *FD_Reco_CClike_Pos;
    TH1D *FD_Reco_CClike_Neg;

    TH2D *FD_Reco2True_CClike;
    TH2D *FD_Reco2True_CClike_Pos;
    TH2D *FD_Reco2True_CClike_Neg;

    TH1D *FD_Purity_CC;//purity of (numu+numubar) in CClike sample
    TH1D *FD_Eff_CC;
    TH1D *FD_Purity_NuMuCC_Neg;//purity of numu in CClike-neg sample
    TH1D *FD_Eff_NuMuCC_Neg;
    TH1D *FD_Purity_NuMuBarCC_Pos;//purity of numubar in CClike-pos sample
    TH1D *FD_Eff_NuMuBarCC_Pos;

    std::map<Background::Background_t,TH2D*> FD_True2Reco_Fid;
    std::map<Background::Background_t,TH2D*> FD_Eff;
    std::map<Background::Background_t,TH2D*> FD_True2Reco_NuBar_Fid;
    std::map<Background::Background_t,TH2D*> FD_Eff_NuBar;


    //MRE efficiency ratio
    TH2D *MREEffRatio;

    //cross sections
    TH1F* fNuMuCCXSec[5];
    TH1F* fNuTauCCXSec[5];
    TH1F* fNueCCXSec[5];
    TH1F* fNuMuBarCCXSec[5];
    TH1F* fNuTauBarCCXSec[5];
    TH1F* fNueBarCCXSec[5];
    std::map<Background::Background_t,TH1D*> XSecWeight;
    std::map<Background::Background_t,TH1D*> XSecWeight_NuBar;



    //these values store the osc. params used to create the 3D prediction histograms...if you write the output to a file, these values are written to a tree along with the 3D histograms.
    Double_t Theta13;
    Double_t Theta12;
    Double_t Theta23;
    Double_t DeltaMSq23;
    Double_t DeltaMSq12;
    Double_t DeltaCP;
    //
    Double_t Eps_ee;
    Double_t Eps_emu;
    Double_t Eps_etau;
    Double_t Eps_mumu;
    Double_t Eps_mutau;
    Double_t Eps_tautau;
    Double_t Delta_emu;
    Double_t Delta_etau;
    Double_t Delta_mutau;
    //
    Double_t Theta14;
    Double_t Theta24;
    Double_t Theta34;
    Double_t Dm41;
    Double_t Delta14;
    Double_t Delta24;

    int OscMethod;

    std::map<OscType_t,TH1D*> NominalOscProb;//for neutrinos
    std::map<OscType_t,TH1D*> NominalOscProb_NuBar;//for antineutrinos

    Bool_t ReadExtrapFromFile;
    string inFileName;

    bool PrintResult;

    bool WriteOutput;

    bool Init;

    bool ReadError;

    TFile *ExtrapFile;

    bool Oscillated;//flag to indicated prediction hists have been oscillated from their default values

    bool UseSeparateNuNuBar;//if separate nu/nubar distributions are included in the extrapolation input file, can oscillate nu and nubars separately; if not, do it the old way

    TH3D* FD_TrueVsRecoVsPID_NC_NueFrac;
    TH3D* FD_TrueVsRecoVsPID_NC_NueBarFrac;
    TH3D* FD_TrueVsRecoVsPID_NC_NuMuFrac;
    TH3D* FD_TrueVsRecoVsPID_NC_NuMuBarFrac;

};


#endif
