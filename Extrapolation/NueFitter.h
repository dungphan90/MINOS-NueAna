#ifndef NUEFITTER_H
#define NUEFITTER_H

// This will be the actual engine that handles fitting
//
//  Assumes that it will read in the files produced during the massive simulation stage
//  allow it to figure out what contours are desired and all that

#include <vector>
#include <map>
#include "TH1D.h"
#include "TH2D.h"

struct Chi2Cut
{
   double percent;
   double cutval;
};

struct Point{
  double th13;
  double nExpected;
  double nSignal;
  std::vector<Chi2Cut> cutVals;
};

class NueFitter
{
   public:
    NueFitter();
    ~NueFitter();

    void Reset();
    void Initialize();

    void AddFile(std::string file) { fFileList.push_back(file); };
    void AddContour(double cont) {fContourLevels.push_back(cont);};
    void SetOutputFile(std::string out) {fOutFile = out;};
                                                                                
    bool PerformFit(double input, std::string outName = "dummy.root");
    bool PrepareFit();
    bool FitInput(double input, std::string outName);    

    bool LoadFiles();
    bool BuildChi2Map(double val);
    bool BuildContourMaps();
    void WriteToFile();
    void SetFitMethod(int meth) {fFitMethod = meth;};
 
  private:
    bool InitializeFittingHistograms();
    double CalculateDeltaLog(double nexp, double nobs);
    double CalculateDeltaLog(double nexp, double nobs, double nsig);

    std::vector<std::string> fFileList;
    std::vector<double> fContourLevels;

    std::map<double, std::map<double, std::vector<Point> > > fDeltaM32MapNormal;
    std::map<double, std::map<double, std::vector<Point> > > fDeltaM32MapInvert;

//    bool PerformFit(double input, vector<Point> &row);

    TH2D* fFitPointChi2N;    
    TH2D* fFitPointChi2I;

    std::vector<TH2D*> fContourHistsN;
    std::vector<TH2D*> fContourHistsI;

    double fMinTotalEvents;
    double fErrors[10];

    std::string fOutFile;

    int fFitMethod;   // 0 - standard chi2, 1 - fit chi2
};

#endif
