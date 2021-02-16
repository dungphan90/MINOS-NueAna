/// $Id: HitCalc.h,v 1.3 2005/05/19 18:18:50 vahle Exp $
///
/// class HitCalc
///
/// NueAna package
///
/// Purpose: Calculate 3D Hits, Hit angular clustering and 
///          Tufts Analysis variables.
///
/// Author: Alex Sousa <asousa@minos.phy.tufts.edu>
///
/// Created on: Tue Mar 29 2005

#ifndef HITCALC_H
#define HITCALC_H

#include "TObject.h"

class TH1F;
class TH2F;

class HitCalc : public TObject
{
    
public:
    // Con/de-structors: 
    HitCalc();
    //    HitCalc(const HitCalc& rhs); // copy constructor
    virtual ~HitCalc();
    
    virtual void Draw(Option_t *option);
    virtual void Print(Option_t *option) const;
    void Reset();
   
    //HitCalc variables
    Float_t fHitTotalEnergy;
    Float_t fHitTransEnergy;
    Float_t fHitLongEnergy;
    Float_t fHitTransCMEnergy;
    Float_t fHitTransEnergyRatio;
    Float_t fHitLongEnergyRatio;
    Float_t fHitTransLongEnergyRatio;
    Float_t fHitTransCMEnergyRatio;
    Float_t fHitFarMomBalance;
    Float_t fHitPeakMomBalance;
    Float_t fHitFarAngle;
    Float_t fHitPeakAngle;


private:
    
    // copy constructor, assignment:
    HitCalc& operator=(const HitCalc& rhs); // assignment
    
    ClassDef(HitCalc,1);
        
};                              // end of class HitCalc

#endif  // HITCALC_H
