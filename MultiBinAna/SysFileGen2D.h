#ifndef SYSFILEGEN2D_H
#define SYSFILEGEN2D_H

#include <map>
#include <vector>
#include <string>
#include "NueAna/MultiBinAna/SysHists2D.h"
#include "NueAna/NueStandard.h"
#include "NueAna/NueRecord.h"
#include "NueAna/Extrapolation/Background.h"
#include "NueAna/NueAnaTools/Selection.h"
#include "OscProb/OscCalc.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TDirectory.h"
#include "TCollection.h"
#include "TH2D.h"
#include "TObject.h"

using namespace std;

class SysFileGen2D{
   public:

    SysFileGen2D();
    virtual ~SysFileGen2D() {};
    
    void Initialize();
    void AddSystematic(std::string sys) { fSysNames.push_back(sys);};//one should be 'Nominal', another could be e.g. 'Xtalk'
    void SetOutputFile(std::string name="DefaultOut.root");
    void SetNearPOT(double pot=1e19) { fNearPOT = pot;};
    void SetFarPOT(double pot=3.14e20)  { fFarPOT = pot;};
    void AddPID(std::string pid);
    void SetCutLevel(Selection::Selection_t s=Selection::kPre) { CutLevel = s; };
    
    double GetNueRecoEnergy(NueRecord *nr);
    double GetCCRecoEnergy(NueRecord *nr);
    std::string GetCurrentSysName() { return fCurrentSysName; };
    
    void ResetHistograms();//initializes hists if it hasn't been done or resets them
    void WriteToFile();
    virtual void FillHistograms();//empty
    void RunSystematics();//for each element in fSysNames: initialize/reset all histograms, call FillHistograms, and write the hists to the file
    
    std::map<Background::Background_t, SysHists2D*> fSysHistsMap;

  private:
    void InitializeHistograms();
    void FinalizeHistograms();//add PID overflow bin to last PID bin
    void SetDefaultOsc();
    
    vector<std::string> fSysNames;
    std::string fCurrentSysName;
    Selection::Selection_t CutLevel;
    std::string outFileName;
    double fNearPOT;
    double fFarPOT;
    vector<std::string> PIDlist;
    
    TFile *fout;
    TTree *paramtree;
    
    double Theta12;
    double Theta13;
    double Theta23;
    double DeltaMSq23;
    double DeltaMSq12;
    double DeltaCP;
    
};

#endif
