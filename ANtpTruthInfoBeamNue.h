/**
 *
 * $Id: ANtpTruthInfoBeamNue.h,v 1.7 2008/07/17 21:35:10 boehm Exp $
 *
 * \class ANtpTruthkInfoBeamNue
 *
 * \package NueAna
 *
 * \brief 
 *
 * Contact: J. Boehm 4/2005
 *
 * Created on: Thu Apr 14 17:37:59 2005
 *
 */

#ifndef ANTPTRUTHINFOBEAMNUE_H
#define ANTPTRUTHINFOBEAMNUE_H


#include "AnalysisNtuples/ANtpTruthInfoBeam.h"

class NtpSREvent;
class NtpSRTrack;

class ANtpTruthInfoBeamNue : public ANtpTruthInfoBeam
{

public:

    ANtpTruthInfoBeamNue();
    ~ANtpTruthInfoBeamNue();

    void Reset();

    Float_t DirCosNeu;
    Float_t DirCosZ_pan;
    Int_t istruckq;   //PDG id of struck quark
    Int_t iflags;     // generator flags   
    Float_t sigmadiff; // differential cross section
    Int_t itg;        //PDG id of target nucleon, p=2212, n=2112

    Int_t fNueClass;
    Float_t fOscProb;
    Float_t fOscProbMatterNormal;
    Float_t fOscProbMatterInverted;
    Float_t fNueWeight;
 
    Float_t Baseline;
    Float_t Ue3Squared;
    Float_t DeltamSquared23;
    Float_t Theta23;
    Float_t Theta12;
    Float_t Theta13;
    Float_t DeltamSquared12;
    Float_t Density;
    Float_t Delta;
   
private:

    ClassDef(ANtpTruthInfoBeamNue, 7) //ANtpTrackInfoNue

    
};                              // end of class ANtpTrackInfoNue

#endif  // ANTPTRACKINFONUE_H
