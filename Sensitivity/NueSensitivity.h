#ifndef NUESENSITIVITY_H
#define NUESENSITIVITY_H
                                                                                
#include <string>
#include "NueAna/NueAnaTools/NueConvention.h"
#include "NueSenseConfig.h"
#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TTree.h"
#include "OscProb/OscCalc.h"
#include <vector>
#include <map>

struct mininfo{
  float energy;
  int NuFlavor;
  int NuFlavorBeforeOsc;
  int nuClass;
  double weight;
};

struct Point{
  double th13;
  double nsignal;
  double nbg;
};

class NueSensitivity
{
  public:
     NueSensitivity();
     ~NueSensitivity();
                                                                                
     void Run(std::string input, std::string output, double outPOT = 4.0);

     void WriteToFile(std::string file);
     float Oscillate(double dm23, double t13, double delta, int sign);

     void RunFromGrid();
     void RunStandardApproach();

     void SetObserved(int in) {fObserved = in; fObservedIsSet = true;}

  private:
     void Initialize();
     void SetupGridRun();
     void LoadEventsFromFile();
     void GetPoint(int i, float &bg, float &sig);
     float OscillateNumber(double Ue32);
     float OscillateHist(double dm23, double Ue32, double delta, int sign);
     float OscillateFile(double dm23, double Ue32, double delta, int sign);
     float CalculateChi2(int i, float bg, float sig);
     float CalculateFitChi2(int i, float bg, float sig);

     int SetupChi2Histograms(NSCDataParam DPst13, NSCDataParam DPdelta, NSCDataParam DPdm23);
     void ProduceTree(TTree * &infotree, TTree * &errtree);

     double OscillateMatter(int nuFlavor, int nonOscNuFlavor,
                                    float Energy);

     double OscillateMatter(int nuFlavor, int nonOscNuFlavor,
                                    float Energy, float dm2, float th13,
                                    float delta, int hierarchy);

     void SetOscParamBase( float dm2, float ss13, float delta, int hierarchy);
 
     float signuenum;    //the unoscillated number of signal nue
     TH1D* signuehist;   //the unoscillated energy spectrum of signal
     TH1D* reweight;     //a clone of signuehist for use in reweighting

     double* fBinCenter;
     
     std::vector<mininfo> mcInfo;
     double fMeasuredPOT;    

     float currentVal[5];  //value of the params at a particular point
                           //-> no scaling

     float* fZeroValBg;    // Number of background events at Ue3 = 0 
     float* fZeroValSig;   // Number of signal events at the Ue3 = Zero
     //Now technically the above two could be evaluatied for each deltam^2_23 
     //  but since the only contribution to this order is from the solar term
     //   we don't have to worry about it, just have one for each error set

     float fBaseLine;

     int fMethod;
     int fNumConfig;
     NueSenseConfig* nsc;

     TH3D** chi2n;  //histograms for normal heirarchy
     TH3D** chi2i;  //histograms for inverted

     float t13Step;
     float dStep;
     float dm23Step;

     //Additional OscPar
     double fDeltaMS12;
     double fTh12;
     double fTh23;
     double fDensity;

     std::map<double, std::map<double, std::vector<Point> > > fDeltaM32MapNormal;
     std::map<double, std::map<double, std::vector<Point> > > fDeltaM32MapInvert;

     OscCalc fOscCalc;

     int fObserved;
     bool fObservedIsSet;
};

#endif
