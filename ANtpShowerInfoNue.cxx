
/**
 *
 * $Id: ANtpShowerInfoNue.cxx,v 1.9 2008/02/24 23:37:04 boehm Exp $
 *
 * \class ANtpShowerInfoNue
 *
 * \package NueAna
 *
 * \brief 
 *
 * Contact: M. Sanchez 4/2005
 *
 * Created on: Thu Apr 14 18:39:51 2005
 *
 */


#include "NueAna/ANtpShowerInfoNue.h"
#include "MessageService/MsgService.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

#include <cassert>
#include <algorithm>

ClassImp(ANtpShowerInfoNue)

CVSID("$Id: ANtpShowerInfoNue.cxx,v 1.9 2008/02/24 23:37:04 boehm Exp $");

//-------------------------------------------------------------------
ANtpShowerInfoNue::ANtpShowerInfoNue() :
  stripRatio(ANtpDefVal::kFloat),
  planeRatio(ANtpDefVal::kFloat),
  pulseHeightRatio(ANtpDefVal::kFloat),
  phMip(ANtpDefVal::kFloat),
  phMeu(ANtpDefVal::kFloat),
  phNueGeV(ANtpDefVal::kFloat),
  phCCGeV(ANtpDefVal::kFloat),
  phNCGeV(ANtpDefVal::kFloat),
  EValUZ0(ANtpDefVal::kFloat),
  EValUZ1(ANtpDefVal::kFloat),
  EValVZ0(ANtpDefVal::kFloat),
  EValVZ1(ANtpDefVal::kFloat),
  longestTrackGapU(ANtpDefVal::kInt),
  longestTrackGapV(ANtpDefVal::kInt),
  gapPlanesU(ANtpDefVal::kInt),
  gapPlanesV(ANtpDefVal::kInt),
  longestTrackGap(ANtpDefVal::kInt),
  gapPlanes(ANtpDefVal::kInt)
{
    MSG("ANtpShowerInfoNue", Msg::kDebug) << "ANtpShowerInfoNue::Constructor" << endl;
    EVecUZ0[0] = EVecUZ1[1] = EVecVZ0[0] = EVecVZ1[1] = ANtpDefVal::kFloat;
}
//----------------------------------------------------------------------
ANtpShowerInfoNue::~ANtpShowerInfoNue()
{

        MSG("ANtpShowerInfoNue", Msg::kDebug) << "ANtpShowerInfoNue::Destructor" << endl;

}
//----------------------------------------------------------------------
void ANtpShowerInfoNue::Reset()
{
        ANtpShowerInfo::Reset();
        stripRatio=ANtpDefVal::kFloat;
        planeRatio=ANtpDefVal::kFloat;
        pulseHeightRatio=ANtpDefVal::kFloat;
        phMeu =ANtpDefVal::kFloat;
        phMip = ANtpDefVal::kFloat;
        phNueGeV = ANtpDefVal::kFloat;
        phCCGeV = ANtpDefVal::kFloat;        
        phNCGeV = ANtpDefVal::kFloat;
	EValUZ0 = ANtpDefVal::kFloat;
	EValUZ1 = ANtpDefVal::kFloat;
	EValVZ0 = ANtpDefVal::kFloat;
	EValVZ1 = ANtpDefVal::kFloat;
	EVecUZ0[0] = EVecUZ1[1] = ANtpDefVal::kFloat;
	EVecVZ0[0] = EVecVZ1[1] = ANtpDefVal::kFloat;	
    longestTrackGapU = ANtpDefVal::kInt;
    longestTrackGapV = ANtpDefVal::kInt;
    gapPlanesU = ANtpDefVal::kInt;
    gapPlanesV = ANtpDefVal::kInt;
    longestTrackGap = ANtpDefVal::kInt;
    gapPlanes = ANtpDefVal::kInt;

        return;
}

