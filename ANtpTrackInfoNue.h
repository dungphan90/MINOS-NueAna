/**
 *
 * $Id: ANtpTrackInfoNue.h,v 1.8 2007/05/22 19:12:54 boehm Exp $
 *
 * \class ANtpTrackInfoNue
 *
 * \package NueAna
 *
 * \brief 
 *
 * Contact: M. Sanchez 4/2005
 *
 * Created on: Thu Apr 14 17:37:59 2005
 *
 */

#ifndef ANTPTRACKINFONUE_H
#define ANTPTRACKINFONUE_H


#include "AnalysisNtuples/ANtpTrackInfo.h"

class NtpSREvent;
class NtpSRTrack;

class ANtpTrackInfoNue : public ANtpTrackInfo
{

public:

    ANtpTrackInfoNue();
    ~ANtpTrackInfoNue();

    void Reset();


    Int_t trklikePlanes;
    Float_t trklikeRatio;
    Float_t pulseHeightRatio;
    Float_t phMip;
    Float_t phMeu;
    Float_t phNueGeV;
    Float_t phCCGeV;
    Float_t trackSignalFull;
    Float_t trackSignalPartial;
    Int_t deltaUVVtx;
    Int_t muonEnergyMethod;
private:

    ClassDef(ANtpTrackInfoNue, 5) //ANtpTrackInfoNue

    
};                              // end of class ANtpTrackInfoNue

#endif  // ANTPTRACKINFONUE_H
