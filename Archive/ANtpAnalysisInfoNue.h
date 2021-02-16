/**
 *
 * $Id: ANtpAnalysisInfoNue.h,v 1.1 2006/09/18 21:41:22 boehm Exp $
 *
 * \class ANtpAnalysisInfoNue
 *
 * \package NueAna
 *
 * \brief 
 *
 * Contact: M. Sanchez
 *
 * Created on: Mon Apr 18 17:35:12 2005
 *
 */

#ifndef ANTPANALYSISINFONUE_H
#define ANTPANALYSISINFONUE_H

#include "AnalysisNtuples/ANtpAnalysisInfo.h"
#include "AnalysisNtuples/ANtpRecoInfo.h"

class NtpSREvent;

class ANtpAnalysisInfoNue : public ANtpRecoInfo
{

public:

    ANtpAnalysisInfoNue();
    ~ANtpAnalysisInfoNue();


    void Reset();

    
    Float_t dpCCPID;
    Double_t nsCCPID;

private:

    ClassDef(ANtpAnalysisInfoNue, 3) //ANtpAnalysisInfoNue

};                              // end of class ANtpAnalysisInfoNue

#endif  // ANTPANALYSISINFONUE_H
