/**
 *
 * $Id: ANtpEventInfoNue.cxx,v 1.13 2008/07/17 21:35:10 boehm Exp $
 *
 * \class ANtpEventInfoNue
 *
 * \package NueAna
 *
 * \brief 
 *
 * Contact: M. Sanchez 4/2005
 *
 * Created on: Fri Apr  8 17:26:32 2005
 *
 */


#include "NueAna/ANtpEventInfoNue.h"
#include "MessageService/MsgService.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

#include <cassert>
#include <algorithm>


ClassImp(ANtpEventInfoNue)

CVSID("$Id: ANtpEventInfoNue.cxx,v 1.13 2008/07/17 21:35:10 boehm Exp $");


//-------------------------------------------------------------------
ANtpEventInfoNue::ANtpEventInfoNue() :
    timeLength(ANtpDefVal::kDouble),
    phMip(ANtpDefVal::kFloat),
    phMeu(ANtpDefVal::kFloat),
    phNueGeV(ANtpDefVal::kFloat),
    triggerPass(ANtpDefVal::kInt),
    hotch(0),
    triggerSource(ANtpDefVal::kInt),
    triggerTime(ANtpDefVal::kFloat),
    spillType(ANtpDefVal::kInt),
    coilStatus(ANtpDefVal::kFloat),
    coilCurrent(ANtpDefVal::kFloat),
    coilCurrentSM2(ANtpDefVal::kFloat),
    coilQuality(ANtpDefVal::kInt),
    coilDirection(ANtpDefVal::kInt),
    liTime(ANtpDefVal::kFloat),
    dmxStatus(ANtpDefVal::kBool),
    eventSignalFull(ANtpDefVal::kFloat),
    eventSignalPartial(ANtpDefVal::kFloat),
    eventTimeMax(ANtpDefVal::kFloat),
    eventTimeMin(ANtpDefVal::kFloat),
    rcBoundary(ANtpDefVal::kInt),
    largestEvent(ANtpDefVal::kInt),
    daveFDDataQuality(ANtpDefVal::kInt)
    
{

        MSG("ANtpEventInfoNue", Msg::kDebug) << "ANtpEventInfoNue::Constructor" << endl;


}
//----------------------------------------------------------------------
ANtpEventInfoNue::~ANtpEventInfoNue()
{

        MSG("ANtpEventInfoNue", Msg::kDebug) << "ANtpEventInfoNue::Destructor" << endl;

}
//----------------------------------------------------------------------
void ANtpEventInfoNue::Reset()
{
    ANtpEventInfo::Reset();
    timeLength=ANtpDefVal::kDouble;
    phMip = ANtpDefVal::kFloat;
    phMeu = ANtpDefVal::kFloat;
    triggerPass = ANtpDefVal::kInt;
    phNueGeV = ANtpDefVal::kFloat;
    hotch = 0;
    triggerSource = ANtpDefVal::kInt;
    triggerTime = ANtpDefVal::kFloat;
    spillType = ANtpDefVal::kInt;
    coilStatus = ANtpDefVal::kFloat;
    coilCurrent = ANtpDefVal::kFloat;
    coilCurrentSM2 = ANtpDefVal::kFloat;
    coilQuality = ANtpDefVal::kInt;
    coilDirection = ANtpDefVal::kInt;
    liTime = ANtpDefVal::kFloat;
    dmxStatus = ANtpDefVal::kBool;
    eventSignalFull = ANtpDefVal::kFloat;
    eventSignalPartial=ANtpDefVal::kFloat;
    eventTimeMax =ANtpDefVal::kFloat;
    eventTimeMin =ANtpDefVal::kFloat;
    rcBoundary = ANtpDefVal::kInt;
    largestEvent = ANtpDefVal::kInt;
    daveFDDataQuality = ANtpDefVal::kInt;

        return;
}
