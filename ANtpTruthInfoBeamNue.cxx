
/**
 *
 * $Id: ANtpTruthInfoBeamNue.cxx,v 1.6 2008/08/27 15:55:34 boehm Exp $
 *
 * \class ANtpTruthInfoBeamNue
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


#include "NueAna/ANtpTruthInfoBeamNue.h"
#include "MessageService/MsgService.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

#include <cassert>
#include <algorithm>

ClassImp(ANtpTruthInfoBeamNue)

CVSID("$Id: ANtpTruthInfoBeamNue.cxx,v 1.6 2008/08/27 15:55:34 boehm Exp $");

//-------------------------------------------------------------------
ANtpTruthInfoBeamNue::ANtpTruthInfoBeamNue() :
    DirCosNeu(ANtpDefVal::kFloat),
    DirCosZ_pan(ANtpDefVal::kFloat),
    fNueClass(ANtpDefVal::kInt),
    fOscProb(ANtpDefVal::kFloat),
    fNueWeight(ANtpDefVal::kFloat),
    Baseline(ANtpDefVal::kFloat),
    Ue3Squared(ANtpDefVal::kFloat),
    DeltamSquared23(ANtpDefVal::kFloat),
    Theta23(ANtpDefVal::kFloat)			
{

    Reset();
    MSG("ANtpTruthInfoBeamNue", Msg::kDebug) << "ANtpTruthInfoBeamNue::Constructor" << endl;

}

//----------------------------------------------------------------------
ANtpTruthInfoBeamNue::~ANtpTruthInfoBeamNue()
{

        MSG("ANtpTruthInfoBeamNue", Msg::kDebug) << "ANtpTruthInfoBeamNue::Destructor" << endl;

}
//----------------------------------------------------------------------
void ANtpTruthInfoBeamNue::Reset()
{
        ANtpTruthInfoBeam::Reset();

        DirCosNeu = ANtpDefVal::kFloat;
        DirCosZ_pan = ANtpDefVal::kFloat;

        fNueClass = ANtpDefVal::kInt;
        fOscProb = ANtpDefVal::kFloat;
        fNueWeight = ANtpDefVal::kFloat;
                                                                                
        Baseline = ANtpDefVal::kFloat;
	Ue3Squared = ANtpDefVal::kFloat;
	DeltamSquared23 = ANtpDefVal::kFloat;
	Theta23 = ANtpDefVal::kFloat;

        istruckq = ANtpDefVal::kInt;
        iflags = ANtpDefVal::kInt;
        sigmadiff = ANtpDefVal::kFloat;
  	itg = ANtpDefVal::kInt;

        fOscProbMatterNormal = ANtpDefVal::kFloat;
        fOscProbMatterInverted = ANtpDefVal::kFloat;
        fNueWeight = ANtpDefVal::kFloat;
                                                                                
        Theta12 = ANtpDefVal::kFloat; 
        Theta13 = ANtpDefVal::kFloat;
        DeltamSquared12 = ANtpDefVal::kFloat;
        Density = ANtpDefVal::kFloat;
        Delta = ANtpDefVal::kFloat;

        return;
}
