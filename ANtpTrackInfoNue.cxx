
/**
 *
 * $Id: ANtpTrackInfoNue.cxx,v 1.8 2007/05/22 17:05:47 boehm Exp $
 *
 * \class ANtpTrackInfoNue
 *
 * \package NueAna
 *
 * \brief 
 *
 * Contact: M. Sanchez 4/2005
 *
 * Created on: Thu Apr 14 18:00:03 2005
 *
 */


#include "NueAna/ANtpTrackInfoNue.h"
#include "MessageService/MsgService.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

#include <cassert>
#include <algorithm>

ClassImp(ANtpTrackInfoNue)

CVSID("$Id: ANtpTrackInfoNue.cxx,v 1.8 2007/05/22 17:05:47 boehm Exp $");

//-------------------------------------------------------------------
ANtpTrackInfoNue::ANtpTrackInfoNue() :
    trklikePlanes(ANtpDefVal::kInt),
    trklikeRatio(ANtpDefVal::kFloat),
    pulseHeightRatio(ANtpDefVal::kFloat),
    phMip(ANtpDefVal::kFloat),
    phMeu(ANtpDefVal::kFloat),
    phNueGeV(ANtpDefVal::kFloat),
    phCCGeV(ANtpDefVal::kFloat),
    trackSignalFull(ANtpDefVal::kFloat),
    trackSignalPartial(ANtpDefVal::kFloat),
    deltaUVVtx(ANtpDefVal::kInt),
    muonEnergyMethod(ANtpDefVal::kInt)
{

    MSG("ANtpTrackInfoNue", Msg::kDebug) << "ANtpTrackInfoNue::Constructor" << endl;

}

//----------------------------------------------------------------------
ANtpTrackInfoNue::~ANtpTrackInfoNue()
{

        MSG("ANtpTrackInfoNue", Msg::kDebug) << "ANtpTrackInfoNue::Destructor" << endl;

}
//----------------------------------------------------------------------
void ANtpTrackInfoNue::Reset()
{
        ANtpTrackInfo::Reset();
        trklikePlanes=ANtpDefVal::kInt;
        trklikeRatio=ANtpDefVal::kFloat;
        pulseHeightRatio=ANtpDefVal::kFloat;

        phMip = ANtpDefVal::kFloat;
        phMeu = ANtpDefVal::kFloat;

        phNueGeV = ANtpDefVal::kFloat;
        phCCGeV = ANtpDefVal::kFloat;

        trackSignalFull = ANtpDefVal::kFloat;
        trackSignalPartial = ANtpDefVal::kFloat;
        deltaUVVtx = ANtpDefVal::kInt;
        muonEnergyMethod = ANtpDefVal::kInt;

        return;
}
