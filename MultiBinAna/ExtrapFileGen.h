#ifndef EXTRAPFILEGEN_H
#define EXTRAPFILEGEN_H

#include <map>
#include <vector>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TDirectory.h"
#include "TCollection.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TObject.h"
#include <fstream>
#include "TSystem.h"
#include "TROOT.h"
#include "NueAna/MultiBinAna/ExtrapHists2D.h"
#include <iostream>

using namespace std;

class ExtrapFileGen{
   public:

    ExtrapFileGen();
    virtual ~ExtrapFileGen() {};
    
    void SetOutputFile(string name="DefaultOut.root");
    void SetNearPOT(double pot=1e19) { fNearPOT = pot;};
    void SetFarPOT(double pot=7e20)  { fFarPOT = pot;};
    void AddPID(string pidname);
    void SetFileLists(string nddata="NULL",string ndmc="NULL",string fdmcbeam="NULL",string fdmcnue="NULL",string fdmctau="NULL");
    void Run();
    void SetFakeData() {fakedata=true;};//use MC instead of data for ND CC-like selection

  private:
    void Initialize();
    void InitializeHistograms();
    void WriteToFile();
    
    ExtrapHists2D *ExHists_NDMC;
    ExtrapHists2D *ExHists_NDData;
    ExtrapHists2D *ExHists_FDMC;
    ExtrapHists2D *ExHists_Nu_FDMC;
    ExtrapHists2D *ExHists_NuBar_FDMC;
    
    string outFileName;
    string FileList_NDData;
    string FileList_NDMC;
    string FileList_FDMC_Beam;
    string FileList_FDMC_Nue;
    string FileList_FDMC_Tau;
    
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
    
    int nReco;
    int nTrue;
    int nPID;
    
    vector<double> RecoEdges;
    vector<double> PIDEdges;
    vector<double> TrueEdges;
    
    bool fakedata;
    
    
};

#endif
