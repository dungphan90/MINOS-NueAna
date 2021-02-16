/**
 *
 * $Id: Ann.cxx,v 1.6 2009/09/13 22:14:25 jjling Exp $
 *
 * \class Ann
 *
 * \package NueAna
 *
 * \brief Hold variables related to the AnnAna package
 **/
#include "NueAna/Ann.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "MessageService/MsgService.h"

CVSID("$Id: Ann.cxx,v 1.6 2009/09/13 22:14:25 jjling Exp $");

ClassImp(Ann)

Ann::Ann():
  pid(ANtpDefaultValue::kDouble),
  pid_30inp(ANtpDefaultValue::kDouble),
  pid_6inp(ANtpDefaultValue::kDouble),
  pid_11inp(ANtpDefaultValue::kDouble),
  pid_11inp_daikon04(ANtpDefaultValue::kDouble),
  pid_14inp_daikon04(ANtpDefaultValue::kDouble)
{}


Ann::~Ann()
{}

void Ann::Reset()
{
   pid = ANtpDefaultValue::kDouble;
   pid_30inp = ANtpDefaultValue::kDouble;
   pid_6inp = ANtpDefaultValue::kDouble;
   pid_11inp = ANtpDefaultValue::kDouble;
   pid_11inp_daikon04 = ANtpDefaultValue::kDouble;
   pid_14inp_daikon04 = ANtpDefaultValue::kDouble;
}

