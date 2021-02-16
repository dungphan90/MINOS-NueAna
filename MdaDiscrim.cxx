/// $Id: MdaDiscrim.cxx,v 1.2 2006/03/01 14:37:11 asousa Exp $
///
/// class MdaDiscrim
///
/// NueAna package
///
/// Purpose:  Obtain PID from discriminant functions obtained using SAS. 
///
/// Author: Alex Sousa <a.sousa@physics.ox.ac.uk>
///
/// Created on: Fri Feb 17 12:17:14 2006

#include <iostream>
#include "NueAna/MdaDiscrim.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"


ClassImp(MdaDiscrim)

MdaDiscrim::MdaDiscrim():
    fMdaPIDnue(ANtpDefVal::kDouble),
    fMdaPIDnc(ANtpDefVal::kDouble),
    fMdaPIDnumu(ANtpDefVal::kDouble),
    fMdaPIDnutau(ANtpDefVal::kDouble),
    fMdaClass(ANtpDefVal::kInt)
{}

MdaDiscrim::~MdaDiscrim()
{}

void MdaDiscrim::Draw(Option_t */*option*/)
{

/// To be filled later.

}

void MdaDiscrim::Print(Option_t */*option*/) const
{
    std::cout<<"There's a lot of stuff to print here!"<<std::endl;
}

void MdaDiscrim::Reset()
{ 

    fMdaPIDnue=ANtpDefVal::kDouble;
    fMdaPIDnc=ANtpDefVal::kDouble;
    fMdaPIDnumu=ANtpDefVal::kDouble;
    fMdaPIDnutau=ANtpDefVal::kDouble;
    fMdaClass=ANtpDefVal::kInt;
   
}
