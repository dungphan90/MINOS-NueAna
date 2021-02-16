#ifndef SYSHISTS2D_H
#define SYSHISTS2D_H

#include "TH1D.h"
#include "TH2D.h"
#include "TDirectory.h"
#include <map>
#include <vector>
#include <string>

using namespace std;

class SysHists2D
{

 public:
   
   SysHists2D(std::string);
   ~SysHists2D();
   
   TDirectory *fDirectory;
   
   TH2D* ND_TrueVsReco;//x = reco energy, y = true energy
   TH2D* FD_TrueVsReco;
   
   std::map<std::string, TH2D*> ND_RecoVsPID;//x = pid, y = reco energy
   std::map<std::string, TH2D*> FD_RecoVsPID;
   
   TH2D *ND_TrueVsReco_ErrorHist;
   TH2D* FD_TrueVsReco_ErrorHist;
   
   std::map<std::string, TH2D*> ND_RecoVsPID_ErrorHist;
   std::map<std::string, TH2D*> FD_RecoVsPID_ErrorHist;
   
   void SetRecoBins(int n, double *edges);
   void SetPIDBins(int n, double *edges);
   void SetTrueBins(int n, double *edges);
   
   void AddPID(std::string pid);
   void InitializeHists();
   
   int GetRecoBins() { return nReco; };
   int GetPIDBins() { return nPID; };
   int GetTrueBins() { return nTrue; };
   
  private:
    
    int nReco;
    int nTrue;
    int nPID;
    
    vector<double> RecoEdges;
    vector<double> PIDEdges;
    vector<double> TrueEdges;
    
    vector<std::string> PIDlist;

};
#endif //SYSHISTS2D_H
