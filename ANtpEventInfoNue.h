/**
 *
 * $Id: ANtpEventInfoNue.h,v 1.15 2008/07/17 21:35:10 boehm Exp $
 *
 * \class ANtpEventInfoNue
 *
 * \package NueAna
 *
 * \brief A class to hold information about an event in the Nue analysis ntuple
 *
 * Contact: M. Sanchez 4/2005
 *
 * Created on: Fri Apr  8 17:00:51 2005
 *
 */

#ifndef ANTPEVENTINFONUE_H
#define ANTPEVENTINFONUE_H

#include "AnalysisNtuples/ANtpEventInfo.h"

class NtpSREvent;


class ANtpEventInfoNue : public ANtpEventInfo
{

public:

    ANtpEventInfoNue();
    ~ANtpEventInfoNue();

    void Reset();

    Double_t timeLength;
    Float_t phMip;
    Float_t phMeu;
    Float_t phNueGeV;
    Int_t triggerPass;

    Int_t hotch; //0=good, 1=bad

    Int_t triggerSource;
    Float_t triggerTime;
    Int_t spillType;
    Float_t coilStatus;
    Float_t coilCurrent;
    Float_t coilCurrentSM2;
    Int_t coilQuality;    // CoilTools IsOK
    Int_t coilDirection;   
    Float_t liTime;
    Bool_t dmxStatus;
    Float_t eventSignalFull;
    Float_t eventSignalPartial;
    Float_t eventTimeMax;
    Float_t eventTimeMin;
    Int_t rcBoundary;
    Int_t largestEvent;
    Int_t daveFDDataQuality;
    Int_t passCosmicCut;

private:

    ClassDef(ANtpEventInfoNue, 10) //ANtpEventInfoNue
   
};                              // end of class ANtpEventInfoNue

#endif  // ANTPEVENTINFONUE_H
