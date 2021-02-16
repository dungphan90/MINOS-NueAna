#ifndef GridGenJoint_h
#define GridGenJoint_h

#include "TROOT.h"
#include "NueAna/MultiBinAna/Extrapolate2D.h"
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "OscProb/OscCalc.h"
#include "NueAna/Extrapolation/Background.h"
#include "TRandom.h"
#include "NueAna/MultiBinAna/GridGen.h"

using namespace std;

class GridGen_Joint : public GridGen {
  public:
    GridGen_Joint();
    virtual ~GridGen_Joint();
    
    void AddExtrapFHC(Extrapolate2D* E);
    void AddExtrapRHC(Extrapolate2D* E);
    void RunMultiBinOscParErrs(string s = "OscParErrDistributions.root");
    void RunMultiBin_VaryTheta13(string s = "OscParErrDistributions.root");
    
  private:
 
    vector<Extrapolate2D*> ExtrapRHC;
    vector<Extrapolate2D*> ExtrapFHC;
    ClassDef(GridGen_Joint,1)    
};

#endif
