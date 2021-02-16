#ifndef NueFit2D_Joint_h
#define NueFit2D_Joint_h

#include "TROOT.h"
#include "NueAna/MultiBinAna/NueFit2D.h"
#include "NueAna/MultiBinAna/Extrapolate2D.h"
#include "NueAna/MultiBinAna/ErrorCalc.h"
#include "NueAna/MultiBinAna/ErrorCalc_Joint.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include <vector>
#include "TMatrixD.h"
#include <string>
#include "TMinuit.h"
#include "TObject.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include <iostream>
#include "TSystem.h"
#include "TRandom.h"
#include "TF1.h"
#include "TMatrixDEigen.h"

using namespace std;

class NueFit2D_Joint: public NueFit2D {
  public:
    NueFit2D_Joint();
    virtual ~NueFit2D_Joint();

    ///////////Moved to NueFit2D.h
    /*
    //Delta Fit Additions
    Double_t Theta12;  
    Double_t Theta12Unc;
    Double_t Theta23;
    Double_t Theta23Unc; 
    Double_t Theta13;
    Double_t Theta13Unc;
    Double_t DeltaMSq32;
    Double_t DeltaMSq32Unc;
    Double_t DeltaMSq21;
    Double_t DeltaMSq21Unc;
    Double_t delta;
    //NSI
    double Eps_ee;
    double Eps_emu;
    double Eps_etau;
    double Eps_etau_unc;
    double Eps_mumu;
    double Eps_mutau; 
    double Eps_tautau; 
    double Delta_emu;
    double Delta_etau;
    double Delta_etau_unc;
    double Delta_mutau;


    void DefTheta13(double t13, double t13unc){
      Theta13 = t13;
      Theta13Unc = t13unc;
      return;
    }

    void DefTheta23(double t23, double t23unc){
      Theta23 = t23;
      Theta23Unc = t23unc;
      return;
    }

    void DefTheta12(double t12, double t12unc){
      Theta12 = t12;
      Theta12Unc = t12unc;
      return;
    }

    void DefDMS32(double dms32, double dms32unc){
      DeltaMSq32 = dms32;
      DeltaMSq32Unc = dms32unc;
      return;
    }

    void DefDMS21(double dms21, double dms21unc){
      DeltaMSq21 = dms21;
      DeltaMSq21Unc = dms21unc;
      return;
    }

    void DefDeltaCP(double deltaCP){
      delta = deltaCP;
      return;
    }

    // NSI


    void DefEps_ee(double eps_ee){
      Eps_ee = eps_ee;
      return;
    }

    void DefEps_emu(double eps_emu){
      Eps_emu = eps_emu;
      return;
    }

    void DefEps_etau(double eps_etau, double eps_etau_unc){
      Eps_etau = eps_etau;
      Eps_etau_unc = eps_etau_unc;
      return;
    }

    void DefDelta_etau(double delta_etau, double delta_etau_unc){
      Delta_etau = delta_etau;
      Delta_etau_unc = delta_etau_unc;
      return;
    }   
    */
    ///////////////////////////////////////////////////////////////////   

    //for this analysis only
    void ScalePredictionFHC(double scale=1.0) { potScaleFHC=scale; };
    void ScalePredictionRHC(double scale=1.0) { potScaleRHC=scale; };

    //Multi-bin FC functions:
    void RunMultiBinPseudoExpts(bool Print=true);
    void ReadGridFiles();
    void RunDataGrid(string filenorm,string fileinvt);
   
    //Various fitting functions:
    double GetSensitivityAt(double delta=0, bool normalhier=true);//returns sin^2 2th13 90% CL upper limit for given values of delta and the hierarchy 
    void RunDeltaChi2Contour(int cl=0);//0 = 90%, 1 = 68.4%
    double GetLikelihood(double t12=0.6,double t23=0.785398,double t13=0,double dm2_32=2.32e-3,double dm2_21=7.59e-5,double delta=0);//gets the likelihood value for NObs and prediction at these osc par values
    void DoDeltaFit();

    //Extrapolations added for joint fit:
    void AddExtrap(Extrapolate2D* E);
    void AddExtrapRHC(Extrapolate2D* E);
    void AddExtrapFHC(Extrapolate2D* E);
    void CombineFHCRHC();
    
    void RunMultiBinPseudoExpts_MHDeltaFit(bool Print=true);
    void RunMultiBinFC_MHDeltaFit();

    //Add NSI Fitting subroutine
    void DoNSIFitQuick();
    void DoNSIFitForever();
    void GetNSIFFX2();
    void RunNSIDeltaChi2Contour(int cl=0);//0 = 90%, 1 = 68.4%   
    void DoIH(int doihinput=0){
      NHIH = doihinput;
    }
    void SetPhaseBroom(double pb=0.0){ 
      if(pb>=0.0 && pb <=1.0){ 
       PhaseBroom = pb; //0 cycles delta_CP, 1 cylces delta_etau, inbetween cycles by fraction
      }
      else{
       PhaseBroom = 0.0;
      }
    }

  private:

    unsigned int nBinsFHC;
    unsigned int nBinsRHC;

    double PhaseBroom;
    int NHIH;

    double GridScale_Normal_FHC;
    double GridScale_Inverted_FHC;
    double GridScale_Normal_RHC;
    double GridScale_Inverted_RHC;
 
    bool printMatrix;
    
    //remove this later
    void CalculateDlnLMatrix(TH2D *SystMatrix, TH2D *HOOHEMatrix);

    double potScaleFHC;
    double potScaleRHC;

    void DefineStdDlnLMinuit();
    void DefineBinDlnLMinuit();

    double GetMinLikelihood(double delta=0,bool normalhier=true);//used only within FC code - gets the best fit likelihood value for NObs; other parameters fixed at values in grid files and given delta
    
    double GetMinLikelihood_Delta(bool normalhier=true);//used only within FC code - gets the best fit likelihood value for NObs; other parameters fixed at values in grid files (including theta13).  Fit is along delta.
    
    vector<Extrapolate2D*> ExtrapRHC;
    vector<Extrapolate2D*> ExtrapFHC;
    ClassDef(NueFit2D_Joint,1)
};


#endif
