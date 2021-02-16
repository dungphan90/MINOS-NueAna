/// $Id: MdaDiscrim.h,v 1.2 2006/03/01 14:37:11 asousa Exp $ 
///
/// class MdaDiscrim
///
/// NueAna package
///
/// Purpose: Obtain PID from discriminant functions obtained using SAS. 
///
/// Author: Alex Sousa <a.sousa@physics.ox.ac.uk>
///
/// Created on: Fri Feb 17 12:17:14 2006

#ifndef MDADISCRIM_H
#define MDADISCRIM_H

#include "TObject.h"

class MdaDiscrim : public TObject
{

public:

    // Typedefs an enumerations:
    
    // Con/de-structors:
    MdaDiscrim();
    virtual ~MdaDiscrim();

    virtual void Draw(Option_t *option);
    virtual void Print(Option_t *option) const; 
    
    void Reset();

    //MdaDiscrim Variable
    Double_t fMdaPIDnue;
    Double_t fMdaPIDnc;
    Double_t fMdaPIDnumu;
    Double_t fMdaPIDnutau;
    
    Int_t   fMdaClass;



private:

    MdaDiscrim& operator=(const MdaDiscrim& rhs); // assignment

    ClassDef(MdaDiscrim,1);

    
};                              // end of class MdaDiscrim

#endif  // MDADISCRIM_H
