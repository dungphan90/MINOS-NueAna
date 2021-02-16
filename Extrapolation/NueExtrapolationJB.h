#ifndef NUEEXTRAPOLATIONJB_H
#define NUEEXTRAPOLATIONJB_H

// THis will be the actual engine that handles a Full Extrapolation
#include <vector>
#include "NueData.h"
#include "NueAna/Extrapolation/NueExtrapHelper.h"
#include "NueAna/Extrapolation/FNHistsM2.h"

class NueExtrapolationJB: public NueExtrapHelper
{
   public:

    NueExtrapolationJB();
    virtual ~NueExtrapolationJB() {};

    void Reset() {};
    void Initialize();

    void AddChain(TChain *ch, double POT, bool isNue = true);
    bool LoadFiles();

    void AddNueSystematic(NueSystematic*);
    void SetSeparation(string );

    void PrepareExtrapHistograms(Selection::Selection_t);
    TH1* GetPrediction(Background::Background_t bg);

    void RunSystematicStudy(Selection::Selection_t sel);
    void RunSystematicStudy(vector<Selection::Selection_t> &sel);
    void RunExtrapStudy(Selection::Selection_t sel);


    void SetOutputFile(string name);
    void WriteToFile(Selection::Selection_t sel);
 
    void UseMCAsData(bool in);
    void UseMCAsCCData(bool in);

    void SetNueRecoBins(int, double, double);
    void SetNueRecoBins(int, double*);
    void SetTrueBins(int, double, double);
    void SetTrueBins(int, double*);
    void SetCCRecoBins(int, double, double);
    void SetCCRecoBins(int, double*);

    void SetDataFileName(string name);
    void SetExtrapMethod(int method) { fExtrapMethod = method;};

  private:
    void ResetExtrapHistograms();
    void InitializeExtrapHistograms();
    void InitializePredictionHistograms();
    void InitializeNeugen();

    void MakePrediction();

    void MakeDataHistograms(Selection::Selection_t sel);
    void LoadDataHistograms(string fname, Selection::Selection_t sel);

    void BuildAppTrueHistExact(Background::Background_t bg, TH1D* trueHist);
    void BuildAppTrueHistFast(Background::Background_t bg, TH1D* trueHist);
    TH1* CreateOscHist(TH2* hist, int nonOsc, int nuFlavor);
                                                                                                                


    TH1 *GetCCRatio(string hist, Background::Background_t bg);
    TH1 *GetFNRatio(Background::Background_t bg);
 


    vector<NueData*> fData;
    std::map<TChain*, int> fChainDataMap;  //Only interesting for setting
                                           //    up loading

    NueSystematic* fCurrentSys;
    vector<NueSystematic*> fSystematics;   //Used for Sys Study or current

    std::map<Background::Background_t, FNHistsM2*> fHistMap;  //Helper Hists
    std::map<Background::Background_t, FNHistsM2*> fDataMap;  //Data Hists
    std::map<Background::Background_t, TH1D*> fPredMap;

    std::string inDataFileName;
    std::string outFileName;
    double fTargetPOT;
    bool fUseMCAsData;
    bool fUseMCAsCCData;

    Int_t     fNZBins;
    Double_t *fZBins;

    double fFarCCPOT;
    double fNearCCPOT;

    double fNearDataPOT;
    double fNearCCDataPOT;

    TH1F* fNuMuXSec;
    TH1F* fNuTauXSec;
    TH1F* fNueXSec;
    TH1F* fNuMuBarXSec;
    TH1F* fNuTauBarXSec;
    TH1F* fNueBarXSec;
  
    TH1F* fNuMuCCXSec[5];
    TH1F* fNuTauCCXSec[5];
    TH1F* fNueCCXSec[5];

    TH1F* fNuMuBarCCXSec[5];
    TH1F* fNuTauBarCCXSec[5];
    TH1F* fNueBarCCXSec[5];

    double fXSecPos[2000];

    TH1D* fTauEff;
    TH1D* fMREEff;
    std::string fSeparation;
    int fExtrapMethod; // 1 - Accurate, 2 - Fast, 
                       // 3 - Crazy fast, only use if only oscillation changes
                       //         are being made (extrap hists only made once)
};

inline   void NueExtrapolationJB::UseMCAsData(bool in) {fUseMCAsData = in; }
inline   void NueExtrapolationJB::UseMCAsCCData(bool in) {fUseMCAsCCData = in; }
inline   void NueExtrapolationJB::SetSeparation(string in) {fSeparation = in;}
inline   void NueExtrapolationJB::SetDataFileName(string name) {inDataFileName = name;}
#endif
