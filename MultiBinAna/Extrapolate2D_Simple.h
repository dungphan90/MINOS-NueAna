#ifndef Extrapolate2D_Simple_h
#define Extrapolate2D_Simple_h

#include "TROOT.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TMath.h"
#include "NueAna/MultiBinAna/MultiBinAnaHelper.h"
#include <vector>
#include "OscProb/OscCalc.h"
#include "NueAna/Extrapolation/Background.h"
#include "TH3D.h"
#include <map>
#include "TFile.h"
#include <iostream>
#include "TTree.h"
#include "NueAna/MultiBinAna/Extrapolate2D.h"

using namespace std;

class Extrapolate2D_Simple: public Extrapolate2D
{
  public:
    Extrapolate2D_Simple();
    virtual ~Extrapolate2D_Simple();
    
    void SetUseInputWithNo3DHists() { UseInputWithNo3D = true; };
    
    void OscillatePrediction();
    
  private:
    
    void ReadNDDataFile();
    void ReadFNFile();
    void ReadFNFileNo3D();
    void ReadFiles();
    
    void SetupPredHists();
    void FarNearPred(Background::Background_t bg);
    void AppearancePred(Background::Background_t bg);
    void BeamNuePred();
    
    void RebinInputHists();
    
    void OscillatePrediction_SepNuNuBar();
    void OscillatePrediction_NoSep();
    
    bool UseInputWithNo3D;
         
    //input hists, FD
    map<Background::Background_t,TH2D*> FD_RecoVsPID;
    
    //input hists, ND
    map<Background::Background_t,TH2D*> ND_RecoVsPID;
    
    //ND data input hists
    TH2D *NDMC_Total;
    TH2D *NDMC_BNueCC;
    TH2D *NDMC_NuMuCC;
    TH2D *NDMC_NC;
    TH2D *NDData_Ratio;
    
};


#endif
