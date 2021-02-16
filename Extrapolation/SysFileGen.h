#ifndef SYSFILEGEN_H
#define SYSFILEGEN_H

// THis will be the actual engine that handles a Full Extrapolation
#include <vector>
#include "NueData.h"
#include "NueAna/Extrapolation/NueExtrapHelper.h"
#include "NueAna/Extrapolation/SysHists.h"

class SysFileGen: public NueExtrapHelper
{
   public:

    SysFileGen();
    virtual ~SysFileGen() {};

    void Reset() {};
    void Initialize();

    void AddSystematic(string sys) { fSystematics.push_back(sys);};
    void PrepareExtrapHistograms(Selection::Selection_t);
    TH1* GetPrediction(Background::Background_t bg);

    void RunSystematicStudy(Selection::Selection_t sel);

    void SetOutputFile(string name);
    void WriteToFile(Selection::Selection_t sel);

    void SetNearPOT(double pot) { fNearPOT = pot;};
    void SetFarPOT(double pot)  { fFarPOT = pot;};
 
    void SetNueRecoBins(int, double, double);
    void SetNueRecoBins(int, double*);
    void SetTrueBins(int, double, double);
    void SetTrueBins(int, double*);
    void InitializeSysHists(SysHists* one);

    void SetDataFileName(string name);
    void SetOscParams(Double_t theta12,Double_t theta23,Double_t theta13,
                                 Double_t deltam12,Double_t deltam23,Double_t deltaCP,
                                 Int_t massH);

  private:
    void ResetHistograms();
    void InitializeHistograms();
    void PrepareHistograms(Selection::Selection_t sel);

    void FillSysHistograms();

    vector<string> fSystematics;   //Used for Sys Study or current
    string fCurrentSys;
  
    std::map<Background::Background_t, SysHists*> fHistMap;  //Helper Hists
    
    std::string outFileName;
    double fTargetPOT;

    Double_t fTheta12;
    Double_t fTheta23;
    Double_t fTheta13;
    Double_t fDeltaMSq23;
    Double_t fDeltaMSq12;
    Double_t fDeltaCP;
    Double_t fMassHierarchy;

};

#endif
