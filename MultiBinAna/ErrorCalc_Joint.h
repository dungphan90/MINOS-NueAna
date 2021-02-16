#ifndef ErrorCalc_Joint_h
#define ErrorCalc_Joint_h

#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "NueAna/MultiBinAna/MultiBinAnaHelper.h"
#include <map>
#include <string>
#include <vector>
#include "NueAna/Extrapolation/Background.h"
#include "TFile.h"
#include "TH3D.h"
#include "TTree.h"
#include "NueAna/MultiBinAna/Extrapolate2D.h"
#include "NueAna/MultiBinAna/ErrorCalc.h"
#include "TSystem.h"
#include <fstream>
#include "TMatrixD.h"

using namespace std;

class ErrorCalc_Joint: public ErrorCalc {
  
  public:
    ErrorCalc_Joint();
    virtual ~ErrorCalc_Joint();
    
    void AddExtrap(Extrapolate2D* E);
    void AddExtrapFHC(Extrapolate2D* E);
    void AddExtrapRHC(Extrapolate2D* E);

    void AddCovarianceMatrix(string systname,string file,string histname,int flag, int isRHC);

    void CalculateHOOError();
    void SetGridPred(int nbins, vector< vector<double> > nc, vector< vector<double> > cc, vector< vector<double> > bnue, vector< vector<double> > tau, vector< vector<double> > sig);

  private:
    void Initialize();
    
    void ReadSysFiles_FNExtrap_AllRuns(int n);//i.e. LisaFiles, combining runs 1,2,&3
    void ReadSysFiles_Appearance_AllRuns(int n);//i.e. LisaFiles, combining runs 1,2,&3
    
    void ReadSpecialFiles(int n);
    void ReadCovarianceFiles(int n);
    
    void CalculateSystErrorMatrixExtrap();//usual function
    void CalculateSystErrorMatrixGrid();
    vector<Extrapolate2D*> Extrap;
    vector<int> ExtrapType;
    
    std::map<string,TH2D*> ExtraCovariance;
    std::map<string,int> ExtraCovarianceTag;
    vector<int> ExtraCovariance_IsRHC;
    

};


#endif
