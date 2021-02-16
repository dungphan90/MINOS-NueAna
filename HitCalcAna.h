/// $Id: HitCalcAna.h,v 1.10 2005/05/23 18:44:07 asousa Exp $
///
/// class HitCalcAna
///
/// NueAna package
///
/// Purpose: Calculate 3D Hit variables passed by HitCalc object. 
///
/// Author: Alex Sousa <asousa@minos.phy.tufts.edu>
///
/// Created on: Mon Apr 11 2005

#ifndef HITCALCANA_H
#define HITCALCANA_H

#include <deque>
#include <vector>
#include "TObject.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/HitCalc.h"
#include "Conventions/PlaneView.h"

class NtpSRRecord;
class NtpSREvent;

class HitCalcAna : public NueAnaBase
{
public:

    // Typedefs
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
 
 
    // Con/de-structors:
    HitCalcAna(HitCalc &hc);
    virtual ~HitCalcAna();
    
    // Member variables
    DeqFloat_t fX;
    DeqFloat_t fY;
    DeqFloat_t fZ;
    DeqFloat_t fE;
    IterDeqFloat_t fHitIndex;

    void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);

    //Calculate 3D Hits for each event
    Int_t ComputeHits(RecRecordImp<RecCandHeader> *srobj, NtpSREvent *event);

    void Get3DHit(DeqFloat_t &x, DeqFloat_t &y, DeqFloat_t &z, DeqFloat_t &e);

    void Reset();
    
private:
    
    HitCalc &fHitCalc;
    
    // copy constructor, assignment:
    // HitCalcAna(const HitCalcAna& rhs); // copy constructor
    // HitCalcAna& operator=(const HitCalcAna& rhs); // assignment
    

    void VarCalc();        // Calculate Hit variables 
    

};                         // end of class HitCalcAna

#endif  // HITCALCANA_H
