/**
 * \class NueBeamMonModule
 *
 * \ingroup NueAna
 *
 * \brief Fill in the beam monitor information in the NueAna ntuples
 *
 *
 * \author (last to touch it) $Author: rhatcher $
 *
 * \version $Revision: 1.5 $
 *
 * \date $Date: 2008/11/19 18:22:51 $
 *
 * Contact: mdier@bnl.gov
 *
 * Created on: Tue Feb 28 14:48:02 2006
 *
 * $Id: NueBeamMonModule.h,v 1.5 2008/11/19 18:22:51 rhatcher Exp $
 *
 */

#ifndef NUEBEAMMONMODULE_H
#define NUEBEAMMONMODULE_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif

#include "BeamDataUtil/BMSpillAna.h"
#include "Conventions/BeamType.h"

class NuePOT;

class NueBeamMonModule : public JobCModule
{

public:
    
    // Con/de-structors:
    NueBeamMonModule();
    ~NueBeamMonModule();

    const Registry& DefaultConfig() const;
    void Config(const Registry& r);

    // Analysis and Reconstruction methods
    JobCResult Reco(MomNavigator* mom);
    void EndJob();

    /// next two were declared but never implemented .. dangerous!
    // float CarrotMCPOT(int bt,Detector::Detector_t det);
    // float DaikonMCPOT(int bt,Detector::Detector_t det);

private:
    BMSpillAna fBMSpillAna;
//    double fTotPot;

//    std::string kOutputFile;
    std::string kPOTTreeName;
    std::string beamstring;
    bool firstevent;
    BeamType::BeamType_t beamtype;
    int zbeamtype;
    int lastrun;
    NuePOT *nuepot;
};                              // end of class NueBeamMonModule


#endif  // NUEBEAMMONMODULE_H
