/// $Id: AngClusterFitAna.h,v 1.7 2006/09/19 15:29:48 boehm Exp $ 
///
/// class AngClusterFitAna
///
/// NueAna package
///
/// Purpose: Fit a EM shower profile to the 3D hits contained in 
///          the primary cluster determined by AngClusterAna. 
///
/// Author: Alex Sousa <asousa@minos.phy.tufts.edu>
///
/// Created on: Fri May 06 2005

#ifndef ANGCLUSTERFITANA_H
#define ANGCLUSTERFITANA_H

#include <deque>
#include <vector>
#include "TObject.h"
#include "TMatrixD.h"
#include "TVector3.h"
#include "TF1.h"
#include "TMinuit.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/AngClusterFit.h"
#include "Conventions/PlaneView.h"

class NtpSRRecord;
class TH1F;
class TH3F;
class TPolyMarker3D;
class TPad;

class AngClusterFitAna : public NueAnaBase
{

public:

    // Typedefs an enumerations:
    typedef std::deque<Int_t> DeqInt_t;
    typedef DeqInt_t::iterator IterDeqInt_t;       
    typedef DeqInt_t::size_type SizeDeqInt_t;
  
    typedef std::deque<Float_t> DeqFloat_t;
    typedef DeqFloat_t::iterator IterDeqFloat_t;       

    typedef std::vector<Int_t> VecInt_t;
    typedef VecInt_t::iterator IterVecInt_t;    

    typedef std::vector<Float_t> VecFloat_t;
    typedef VecFloat_t::iterator IterVecFloat_t;    

    typedef std::deque< std::deque <Int_t> > DeqDeqInt_t;
    typedef DeqDeqInt_t::iterator IterDeqDeqInt_t;
    typedef DeqDeqInt_t::size_type SizeDeqDeqInt_t;
 
    typedef std::deque <TMatrixD> DeqTMatrixD;
    typedef DeqTMatrixD::iterator IterDeqTMatrixD; 
    
    typedef std::deque <TPolyMarker3D *> DeqTPoly; 
    typedef DeqTPoly::iterator IterDeqTPoly; 
    

    // Con/de-structors:
    AngClusterFitAna(AngClusterFit &acf);
    ~AngClusterFitAna();

    TH1F *fLongHitHisto;  //Longitudinal profile histogram.
    TH1F *fLongHitHistoRad;  //Longitudinal profile histogram.
    TH1F *fTransHitHisto; //Transverse profile histogram.
    TH1F *fTransHitHistoRad; //Transverse profile histogram.

    TH3F *fScatterAxis; //Create 3D plot axis.

    TF1 *fLongFit;  //Longitudinal Fitting function
    TF1 *fTransFit; //Transverse Fitting function

    void Set3DHit(DeqFloat_t &x, DeqFloat_t &y, DeqFloat_t &z, DeqFloat_t &e);

    void SetAngCluster(Int_t &primShow
                       , DeqDeqInt_t &clusterMap
                       , TVector3 &primDir);
    
    void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);
    
    void FitShower(Float_t totalHitEnergy);

    void TransVarCalc(TH1F *transHisto);

    void Draw(TPad *pad);
  
    void Reset();

private:
     Double_t GetMaximum(TF1* efit, Double_t xmin=0, Double_t xmax=0);
    AngClusterFit &fAngClusterFit;

    // copy constructor, assignment:
    // AngClusterFitAna(const AngClusterFitAna& rhs); // copy constructor
    // AngClusterFitAna& operator=(const AngClusterFitAna& rhs); // assignment

    //To be set to the referenced HitCalcAna Hit member vectors.
    DeqFloat_t fX; 
    DeqFloat_t fY; 
    DeqFloat_t fZ; 
    DeqFloat_t fE; 

    //To be set to referenced AngClusterAna cluster quantities.
    Int_t fPrimShow;
    DeqDeqInt_t fClusterMap;
    TVector3 fPrimDir;
};                              // end of class AngClusterFitAna

#endif  // ANGCLUSTERFITANA_H
