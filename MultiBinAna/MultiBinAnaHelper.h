#ifndef MultiBinAnaHelper_h
#define MultiBinAnaHelper_h

#include "TH2D.h"
#include <iostream>
#include "TCanvas.h"
#include "TMath.h"
#include "TH1D.h"
#include <string>
#include "TH3D.h"

using namespace std;

class MultiBinAnaHelper {
  public :
    MultiBinAnaHelper();
    virtual ~MultiBinAnaHelper();
    void Rebin2DHist(TH2D *h,Int_t nx,Double_t *x,Int_t ny,Double_t *y);
    void Rebin3DHist(TH3D*& h,Int_t nx,Double_t *x,Int_t ny,Double_t *y,Int_t nz,Double_t *z);
    void CopyYtoZ(TH2D* h2, TH3D*& h3);
};


#endif
