/// $Id: AngClusterFit.h,v 1.7 2005/06/03 17:14:27 asousa Exp $ 
///
/// class AngClusterFit
///
/// NueAna package
///
/// Purpose: Return 3D Hit angular clustering shower fitting variables. 
///
/// Author: Alex Sousa <asousa@minos.phy.tufts.edu>
///
/// Created on: Fri May 06 2005


#ifndef ANGCLUSTERFIT_H
#define ANGCLUSTERFIT_H

#include "TObject.h"

class TH1F;
class TH2F;

class AngClusterFit : public TObject
{

public:

    // Con/de-structors:
    AngClusterFit();
    ~AngClusterFit();

    virtual void Draw(Option_t *option);
    virtual void Print(Option_t *option) const; 

    void Reset();

    //Cluster shower fit variables
    Float_t fACluFitParA;
    Float_t fACluFitParB;
    Float_t fACluFitParLongE0;
    Float_t fACluFitShwMax;
    Float_t fACluFitE0EnergyRatio;
    
    Float_t fACluFitParL1;
    Float_t fACluFitParL2;
    Float_t fACluFitParC12;
    Float_t fACluFitParTransE0;
 
    Float_t fACluFitLongChiSq;
    Int_t   fACluFitLongConv;
    Float_t fACluFitLongNDF;


    Float_t fACluFitTransChiSq;
    Int_t   fACluFitTransConv;
    Float_t fACluFitTransNDF;

    Float_t fACluFitAsymPeak;
    Float_t fACluFitAsymVert;
    Float_t fACluFitMolRadPeak;
    Float_t fACluFitMolRadVert;
    Float_t fACluFitMean;
    Float_t fACluFitRMS;
    Float_t fACluFitSkew;
    Float_t fACluFitKurt;

private:

    // copy constructor, assignment:
    // AngClusterFit(const AngClusterFit& rhs); // copy constructor
    // AngClusterFit& operator=(const AngClusterFit& rhs); // assignment

    ClassDef(AngClusterFit,3);

};                              // end of class AngClusterFit

#endif  // ANGCLUSTERFIT_H

