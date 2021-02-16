/**
 *
 * $Id: AnalysisInfoNue.h,v 1.5 2012/02/14 22:29:31 whitehd Exp $
 *
 * \class AnalysisInfoNue
 *
 * \package NueAna
 *
 * \brief 
 *
 * Contact: J Boehm
 *
 *
 */

#ifndef ANALYSISINFONUE_H
#define ANALYSISINFONUE_H
#include "TObject.h"
                                                                                
class NtpSREvent;

class AnalysisInfoNue : public TObject
{
  public:

    AnalysisInfoNue();
    ~AnalysisInfoNue();


    void Reset();

  Int_t inFiducialVolume;     //flag for whether the neutrino interaction is
                              //in the fiducial volume
  Int_t isFullyContained;     //flag for whether the event is fully contained.
  Float_t dpCCPID;
  Double_t nsCCPID;
  Double_t roCCPID;
  Double_t abCCPID;
  
  
  //andy's variables
  float abcostheta;
  float abeventenergy;
  float abtrackcharge;
  float abtrackenergy;
  float abtracklikeplanes;
  float abtrackphfrac;
  float abtrackphmean;
  float abtrackqpsigmaqp;
  float aby;

  //rustem's variables
  float roNScientPlanes;
  float roMeanTrkSig;
  float roMeanRatio;
  float roTrkMeanVsWindow;
  float roNuMuBar;
  float roRelAng;

  bool isRHC;

  private:

    ClassDef(AnalysisInfoNue, 3) //AnalysisInfoNue

};                              // end of class AnalysisInfoNue

#endif  // ANALYSISINFONUE_H
