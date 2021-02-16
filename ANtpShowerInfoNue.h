/**
 *
 * $Id: ANtpShowerInfoNue.h,v 1.10 2008/02/24 23:37:04 boehm Exp $
 *
 * \class ANtpShowerInfoNue
 *
 * \package NueAna
 *
 * \brief 
 *
 * Contact: M. Sanchez
 *
 * Created on: Thu Apr 14 18:36:34 2005
 *
 */

#ifndef ANTPSHOWERINFONUE_H
#define ANTPSHOWERINFONUE_H


#include "AnalysisNtuples/ANtpShowerInfo.h"

class NtpSREvent;
class NtpSRShower;


class ANtpShowerInfoNue : public ANtpShowerInfo
{

public:

    ANtpShowerInfoNue();
    ~ANtpShowerInfoNue();

    void Reset();

    Float_t stripRatio;
    Float_t planeRatio;
    Float_t pulseHeightRatio;
    Float_t phMip;
    Float_t phMeu;
    Float_t phNueGeV;
    Float_t phCCGeV;
    Float_t phNCGeV;

    //sphericity measures:
    Float_t EValUZ0;
    Float_t EValUZ1;
    Float_t EVecUZ0[2];
    Float_t EVecUZ1[2];
    Float_t EValVZ0;
    Float_t EValVZ1;
    Float_t EVecVZ0[2];
    Float_t EVecVZ1[2];

                                                                                                                                                     
    Int_t longestTrackGapU;
    Int_t longestTrackGapV;
    Int_t gapPlanesU;
    Int_t gapPlanesV;
    Int_t longestTrackGap;
    Int_t gapPlanes;


private:

    ClassDef(ANtpShowerInfoNue, 7) //ANtpShowerInfoNue
    
};                              // end of class ANtpShowerInfoNue

#endif  // ANTPSHOWERINFONUE_H
