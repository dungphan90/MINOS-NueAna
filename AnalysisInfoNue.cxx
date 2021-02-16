/**
 *
 * $Id: AnalysisInfoNue.cxx,v 1.6 2012/02/14 22:29:31 whitehd Exp $
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


#include "NueAna/AnalysisInfoNue.h"
#include "MessageService/MsgService.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

#include <cassert>
#include <algorithm>

ClassImp(AnalysisInfoNue)

CVSID("$Id: AnalysisInfoNue.cxx,v 1.6 2012/02/14 22:29:31 whitehd Exp $");

AnalysisInfoNue::AnalysisInfoNue():
    inFiducialVolume(ANtpDefVal::kInt),
    isFullyContained(ANtpDefVal::kInt),
    dpCCPID(ANtpDefVal::kFloat),
    nsCCPID(ANtpDefVal::kDouble),
    roCCPID(ANtpDefVal::kDouble),
    abCCPID(ANtpDefVal::kDouble),
    abcostheta(ANtpDefVal::kFloat),
    abeventenergy(ANtpDefVal::kFloat),
    abtrackcharge(ANtpDefVal::kFloat),
    abtrackenergy(ANtpDefVal::kFloat),
    abtracklikeplanes(ANtpDefVal::kFloat),
    abtrackphfrac(ANtpDefVal::kFloat),
    abtrackphmean(ANtpDefVal::kFloat),
    abtrackqpsigmaqp(ANtpDefVal::kFloat),
    aby(ANtpDefVal::kFloat),
    roNScientPlanes(ANtpDefVal::kFloat),
    roMeanTrkSig(ANtpDefVal::kFloat),
    roMeanRatio(ANtpDefVal::kFloat),
    roTrkMeanVsWindow(ANtpDefVal::kFloat),
    roNuMuBar(ANtpDefVal::kFloat),
    roRelAng(ANtpDefVal::kFloat),
    isRHC(false)

{

       MSG("AnalysisInfoNue", Msg::kDebug) 
           << "AnalysisInfoNue::Constructor" << endl;

}

//----------------------------------------------------------------------
AnalysisInfoNue::~AnalysisInfoNue()
{

     MSG("AnalysisInfoNue", Msg::kDebug) 
        << "AnalysisInfoNue::Destructor" << endl;

}
//----------------------------------------------------------------------
void AnalysisInfoNue::Reset()
{
  inFiducialVolume = 0;
  isFullyContained = 0;
  dpCCPID =ANtpDefVal::kFloat;
  nsCCPID =ANtpDefVal::kDouble;
  roCCPID =ANtpDefVal::kDouble;
  abCCPID =ANtpDefVal::kDouble;

  abcostheta=ANtpDefVal::kFloat;
  abeventenergy=ANtpDefVal::kFloat;
  abtrackcharge=ANtpDefVal::kFloat;
  abtrackenergy=ANtpDefVal::kFloat;
  abtracklikeplanes=ANtpDefVal::kFloat;
  abtrackphfrac=ANtpDefVal::kFloat;
  abtrackphmean=ANtpDefVal::kFloat;
  abtrackqpsigmaqp=ANtpDefVal::kFloat;
  aby=ANtpDefVal::kFloat;
  roNScientPlanes =ANtpDefVal::kFloat;
  roMeanTrkSig = ANtpDefVal::kFloat; 
  roMeanRatio = ANtpDefVal::kFloat;
  roTrkMeanVsWindow = ANtpDefVal::kFloat;
  roNuMuBar = ANtpDefVal::kFloat; 
  roRelAng = ANtpDefVal::kFloat;
  isRHC = false;

  return;
}
