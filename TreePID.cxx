/**
 * \class TreePID
 *
 * \ingroup NueAna
 *
 * \brief Calculate PID from Tree or Cut methods and perform classification.
 *
 * Author: Mayly Sanchez (msanchez@physics.harvard.edu)
 *
 * \version $Revision: 1.4 $
 *
 * \date $Date: 2007/05/25 20:58:47 $
 *
 * Created on: Thu Mar 16 19:26:32 2006
 *
 * $Id: TreePID.cxx,v 1.4 2007/05/25 20:58:47 msanchez Exp $
 *
 */

#include "NueAna/TreePID.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

ClassImp(TreePID)

TreePID::TreePID():
    fPassBCuts(ANtpDefVal::kInt),
    fCutPID(ANtpDefVal::kInt),
    fCutPID1(ANtpDefVal::kInt),
    fCutPID2(ANtpDefVal::kInt),
    fCutPID3(ANtpDefVal::kInt),
    fCutClass(ANtpDefVal::kInt),    
    fCutClass1(ANtpDefVal::kInt),
    fCutClass2(ANtpDefVal::kInt),
    fCutClass3(ANtpDefVal::kInt)
{}

TreePID::~TreePID() 
{}

void TreePID::Draw(Option_t */*option*/) 
{}

void TreePID::Print(Option_t */*option*/) const
{}

void TreePID::Reset() 
{
    fPassBCuts=ANtpDefVal::kInt;
    fCutPID=ANtpDefVal::kInt;
    fCutPID1=ANtpDefVal::kInt;
    fCutPID2=ANtpDefVal::kInt;
    fCutPID3=ANtpDefVal::kInt;
    fCutClass=ANtpDefVal::kInt;
    fCutClass1=ANtpDefVal::kInt;
    fCutClass2=ANtpDefVal::kInt;
    fCutClass3=ANtpDefVal::kInt;
}
