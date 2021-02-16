/// $Id: MdaDiscrimAna.h,v 1.6 2006/10/30 19:19:24 asousa Exp $ 
///
/// class MdaDiscrimAna
///
/// NueAna package
///
/// Purpose: Calculate PID from SAS discriminant functions and perform 
///          posterior event classification. 
///
/// Author: Alex Sousa <a.sousa@physics.ox.ac.uk>
///
/// Created on:  Fri Feb 17 15:50:47 2006

#ifndef MDADISCRIMANA_H
#define MDADISCRIMANA_H

#include <deque>
#include <vector>
#include "TObject.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/MdaDiscrim.h"

class NueRecord;
//class MdaDiscrim;

class MdaDiscrimAna : public TObject
{

public:

    // Typedefs and enumerations:
    typedef std::deque< std::deque <Int_t> > DeqDeqInt_t;
    typedef DeqDeqInt_t::iterator IterDeqDeqInt_t;
    typedef DeqDeqInt_t::size_type SizeDeqDeqInt_t;

    typedef std::deque<Double_t> DeqDouble_t;
    typedef DeqDouble_t::iterator IterDeqDouble_t;    
 
    typedef std::deque <TMatrixD> DeqTMatrixD;
    typedef DeqTMatrixD::iterator IterDeqTMatrixD;

    typedef std::deque <TVectorD> DeqTVectorD;
    typedef DeqTVectorD::iterator IterDeqTVectorD;
  
    // Con/de-structors:
    MdaDiscrimAna(NueRecord &nr, MdaDiscrim &md);
    ~MdaDiscrimAna();

    void Analyze();

    void Reset();
    
    Bool_t NeedsSpecialAttention(TString name
                                 , NueRecord *nr
                                 , Double_t &value);

    void SetMdaParams(Double_t threshC, std::string sasFN){
        threshCut=threshC;
        sasFileName=sasFN;
    }

    //Member Variables
    DeqTVectorD varPVector;

    DeqTVectorD discrimDiff;
   
    DeqDouble_t probClass;
    

    // Static member declaration.
    static bool isFilled;    
    bool isFillDone() const {return isFilled; }
        
private:

    NueRecord &nueRec;
    MdaDiscrim &fMdaDiscrim;

    // copy constructor, assignment:
    MdaDiscrimAna(const MdaDiscrimAna& rhs); // copy constructor
    MdaDiscrimAna& operator=(const MdaDiscrimAna& rhs); // assignment
  
    //Member functions
    Bool_t FillCalibArrays();

    void FillPVector();

    // Data members:
   
    void PIDCalc();
       
    void MdaClassify();    

    //Member Variables
    Double_t threshCut;
    std::string sasFileName;
    
    static Int_t nClass;

    static std::deque <TMatrixD> quadCoeff;
    static std::deque <TVectorD> linCoeff;
    static std::deque <TVectorD> constCoeff;
    static std::deque <TVectorD> meanVec;   

    static std::deque < std::deque <std::string> > typeClassVec;
    static std::deque < std::string> typeCoeffVec;
    static std::deque < std::string> varNameVec;


    ClassDef(MdaDiscrimAna,1)
     
};                              // end of class MdaDiscrimAna

#endif  // MDADISCRIMANA_H
