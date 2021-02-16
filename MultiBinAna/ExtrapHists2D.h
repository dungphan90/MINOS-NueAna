#ifndef EXTRAPHISTS2D_H
#define EXTRAPHISTS2D_H

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TDirectory.h"
#include <map>
#include <vector>
#include <string>

using namespace std;

class ExtrapHists2D
{

 public:
   
   ExtrapHists2D(std::string);
   ~ExtrapHists2D();
   
   TDirectory *fDirectory;
   
   TH2D *TrueVsReco_Fid;
   TH2D *TrueVsReco_Fid_NuMuCC;
   TH2D *TrueVsReco_Fid_NuMuBarCC;
   TH2D *TrueVsReco_CClike;
   TH2D *TrueVsReco_CClike_Pos;
   TH2D *TrueVsReco_CClike_Neg;
   TH2D *TrueVsReco_CClike_NuMuCC;
   TH2D *TrueVsReco_CClike_Pos_NuMuBarCC;
   TH2D *TrueVsReco_CClike_Neg_NuMuCC;
   TH1D *Reco_CClike;
   TH1D *Reco_CClike_Pos;
   TH1D *Reco_CClike_Neg;
   
   TH2D *TrueVsReco_Fid_NueCC;
   TH2D *TrueVsReco_Fid_NuTauCC;
   TH1D *True_Fid_NueCC;
   TH1D *True_Fid_NuTauCC;
   TH2D *TrueVsReco_Fid_NueBarCC;
   TH2D *TrueVsReco_Fid_NuTauBarCC;
   TH1D *True_Fid_NueBarCC;
   TH1D *True_Fid_NuTauBarCC;
   
   map<string,TH3D*> TrueVsRecoVsPID_NC;
   map<string,TH3D*> TrueVsRecoVsPID_NuMuToNuMu;
   map<string,TH3D*> TrueVsRecoVsPID_BNueToNuMu;
   map<string,TH3D*> TrueVsRecoVsPID_BNueToBNue;
   map<string,TH3D*> TrueVsRecoVsPID_NuMuToNue;
   map<string,TH3D*> TrueVsRecoVsPID_NuMuToNuTau;
   map<string,TH3D*> TrueVsRecoVsPID_BNueToNuTau;
   
   map<string,TH3D*> TrueVsRecoVsPID_NC_NoOscNue;
   map<string,TH3D*> TrueVsRecoVsPID_NC_NoOscNuMu;
   map<string,TH3D*> TrueVsRecoVsPID_NC_NoOscNueBar;
   map<string,TH3D*> TrueVsRecoVsPID_NC_NoOscNuMuBar;
   
   void SetRecoBins(int n, double *edges);
   void SetPIDBins(int n, double *edges);
   void SetTrueBins(int n, double *edges);
   
   void SetPOT(double p);
   
   void AddPID(std::string pid);
   void InitializeHists();
   
   int GetRecoBins() { return nReco; };
   int GetPIDBins() { return nPID; };
   int GetTrueBins() { return nTrue; };
   double GetPOT() { return POT; };
   
  private:
    
    int nReco;
    int nTrue;
    int nPID;
    
    vector<double> RecoEdges;
    vector<double> PIDEdges;
    vector<double> TrueEdges;
    
    vector<std::string> PIDlist;
    
    double POT;

};
#endif //SYSHISTS2D_H
