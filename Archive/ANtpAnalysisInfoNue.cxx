/**
 *
 * $Id: ANtpAnalysisInfoNue.cxx,v 1.1 2006/09/18 21:41:22 boehm Exp $
 *
 * \class ANtpAnalysisInfoNue
 *
 * \package NueAna
 *
 * \brief 
 *
 * Contact: M. Sanchez
 *
 * Created on: Mon Apr 18 17:41:57 2005
 *
 */


#include "NueAna/ANtpAnalysisInfoNue.h"
#include "MessageService/MsgService.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

#include <cassert>
#include <algorithm>

ClassImp(ANtpAnalysisInfoNue)

CVSID("$Id:");

ANtpAnalysisInfoNue::ANtpAnalysisInfoNue():
    dpCCPID(ANtpDefVal::kFloat),
    nsCCPID(ANtpDefVal::kFloat)
{

        MSG("ANtpAnalysisInfoNue", Msg::kDebug) << "ANtpAnalysisInfoNue::Constructor" << endl;

}

//----------------------------------------------------------------------
ANtpAnalysisInfoNue::~ANtpAnalysisInfoNue()
{

        MSG("ANtpAnalysisInfoNue", Msg::kDebug) << "ANtpAnalysisInfoNue::Destructor" << endl;

}
//----------------------------------------------------------------------
void ANtpAnalysisInfoNue::Reset()
{
  //ANtpAnalysisInfo::Reset();
	ANtpRecoInfo::Reset();
      dpCCPID =ANtpDefVal::kFloat;
      nsCCPID =ANtpDefVal::kFloat;


        return;
}
