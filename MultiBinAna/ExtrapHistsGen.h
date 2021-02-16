#ifndef EXTRAPHISTSGEN_H
#define EXTRAPHISTSGEN_H

#include <map>
#include <vector>
#include <string>
#include "ExtrapHists2D.h"
#include "NueAna/NueStandard.h"
#include "NueAna/NueRecord.h"
#include "NueAna/NuePOT.h"
#include "NueAna/Extrapolation/Background.h"
#include "NueAna/NueAnaTools/Selection.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "OscProb/OscCalc.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TDirectory.h"
#include "TCollection.h"
#include "TH2D.h"
#include "TObject.h"
#include <fstream>
#include "TChain.h"
#include "TSystem.h"

using namespace std;

class ExtrapHistsGen{
   public:

    ExtrapHistsGen();
    virtual ~ExtrapHistsGen() {};
    
    void SetOutputFile(string name="DefaultOut.root");
    void AddPID(string pidname,Selection::Selection_t pidsel);
    void SetFileList(string list="NULL");
    void FillHistograms();
    
    void SetRecoBins(int n, double *edges);
    void SetPIDBins(int n, double *edges);
    void SetTrueBins(int n, double *edges);
    
    void SetRHC() {isRHC = true;};
    
  private:
    void Initialize();
    void InitializeHistograms();
    void FinalizeHistograms();//add PID overflow bin to last PID bin
    void SetDataType(string type="NULL");
    void WriteToFile();
    
    double GetNueRecoEnergy(NueRecord *nr);
    double GetCCRecoEnergy(NueRecord *nr);
    
    ExtrapHists2D *ExHists;
    ExtrapHists2D *ExHists_Nu;
    ExtrapHists2D *ExHists_NuBar;
    
    string outFileName;
    string FileList;
    string DataType;
    
    map<string,Selection::Selection_t> PIDlist;
    
    TFile *fout;
    TTree *paramtree;
    
    double Theta12;
    double Theta13;
    double Theta23;
    double DeltaMSq23;
    double DeltaMSq12;
    double DeltaCP;
    
    int nReco;
    int nTrue;
    int nPID;
    
    vector<double> RecoEdges;
    vector<double> PIDEdges;
    vector<double> TrueEdges;
    
    bool isRHC;//need to set this in order to call RHCNueEnergyCorrection()
    
};

#endif
