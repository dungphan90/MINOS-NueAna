/**
 * \class TreePIDAna
 *
 * \ingroup NueAna
 *
 * \brief Calculate PID from Tree or Cut methods and perform classification.
 *
 * Author: Mayly Sanchez (msanchez@physics.harvard.edu)
 *
 * \version $Revision: 1.11 $
 *
 * \date $Date: 2008/11/19 18:22:51 $
 *
 * Created on: Thu Mar 16 20:10:10 2006
 *
 * $Id: TreePIDAna.cxx,v 1.11 2008/11/19 18:22:51 rhatcher Exp $
 *
 */

#include "NueAna/TreePIDAna.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "NueAna/NueRecord.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "MessageService/MsgService.h"

CVSID("$Id: TreePIDAna.cxx,v 1.11 2008/11/19 18:22:51 rhatcher Exp $");


ClassImp(TreePIDAna)

TreePIDAna::TreePIDAna(NueRecord &nr, TreePID &tp):
    nueRec(nr),
    fTreePID(tp)
{}

TreePIDAna::~TreePIDAna()
{}

void TreePIDAna::Analyze()
{
    NueRecord *nr=&nueRec;

    fTreePID.fPassBCuts=1;
    fTreePID.fCutPID=1;
    fTreePID.fCutClass=1;
    fTreePID.fCutPID1=1;
    fTreePID.fCutClass1=1;
    fTreePID.fCutPID2=1;
    fTreePID.fCutClass2=1;
    fTreePID.fCutPID3=1;
    fTreePID.fCutClass3=1;

      // Cut for no reco event
      if (nr->GetHeader().GetEventNo()<0) fTreePID.fPassBCuts=0;
      // Cut for Fiducial volume
      if (nr->anainfo.inFiducialVolume != 1) fTreePID.fPassBCuts=0;
      // Cut for Full Containment, note difference in Near/Far
      if (nr->GetHeader().GetVldContext().GetDetector()== Detector::kFar
          && nr->anainfo.isFullyContained != 1) fTreePID.fPassBCuts=0;
      if (nr->GetHeader().GetVldContext().GetDetector()== Detector::kNear
          && nr->anainfo.isFullyContained != 1 
          && nr->anainfo.isFullyContained != -2) fTreePID.fPassBCuts=0;


      Int_t kHiPlaneTrackCut=25;
      Int_t kPhProngCut=5000; //sigcor
      Int_t kMeuEnergyCut=150; //Meu
      Int_t kHiTrackLikeCut=18;
      Int_t kLoPlaneEventCut=-1;
      Int_t kHiEnergyCut=-1;//Meu
      Int_t kLoEnergyCut=-1;//Meu
      Int_t kHiEnergyShowerCut=-1; //gev
      Int_t kLoEnergyShowerCut=-1; //gev
      Int_t kLoPhNStripCut=-1;
      Int_t kLoPhNPlaneCut=-1;

      // Cut for max number of planes in Track
      if (nr->srtrack.planes>=kHiPlaneTrackCut&&kHiPlaneTrackCut>=0
          ) fTreePID.fPassBCuts=0;
      // Cut on max total event energy in Meu
      if(nr->srevent.phMip>kMeuEnergyCut&&kMeuEnergyCut>=0) fTreePID.fPassBCuts=0;
      // Cut on min Total pulse height per prong (sigcor)
      if(TMath::Max(nr->srtrack.pulseHeight,nr->srshower.pulseHeight)
	 <kPhProngCut&&kPhProngCut>=0) fTreePID.fPassBCuts=0;
      // Cut on max number of trklike planes in Track
      if (nr->srtrack.trklikePlanes>=kHiTrackLikeCut&&kHiTrackLikeCut>=0
         ) fTreePID.fPassBCuts=0;
      // Cut on min number of planes in Event
      if (nr->srevent.planes<kLoPlaneEventCut&&kLoPlaneEventCut>=0
          ) fTreePID.fPassBCuts=0;
      // Cut on max total event energy in Meu
      if (nr->srevent.phMip>=kHiEnergyCut
	  &&kHiEnergyCut>=0) fTreePID.fPassBCuts=0;
      // Cut on min total event energy in Meu
      if (nr->srevent.phMip<kLoEnergyCut
	  &&kLoEnergyCut>=0) fTreePID.fPassBCuts=0;
      // Cut on max shower energy (in gev) 
      if (nr->srevent.phNueGeV>=kHiEnergyShowerCut
	  &&kHiEnergyShowerCut>=0) fTreePID.fPassBCuts=0;
      // Cut on min shower energy (in gev) 
      if (nr->srevent.phNueGeV<kLoEnergyShowerCut
	  &&kLoEnergyShowerCut>=0) fTreePID.fPassBCuts=0;
      // Cut on min number of strips with above a threshold ph
      if (nr->shwfit.hiPhStripCount<kLoPhNStripCut
	  &&kLoPhNStripCut>=0) fTreePID.fPassBCuts=0;
      // Cut on min number of planes with above a threshold ph
      if (nr->shwfit.hiPhPlaneCount<kLoPhNPlaneCut
	  &&kLoPhNPlaneCut>=0) fTreePID.fPassBCuts=0;

      
      if((nr->srevent.planes<1.5||nr->srevent.planes>15.5) ||
         (nr->hitcalc.fHitLongEnergy<42.8||nr->hitcalc.fHitLongEnergy>315.)||
         (nr->shwfit.uv_molrad_vert<0.70||nr->shwfit.uv_molrad_vert>5.13)||
         (nr->shwfit.par_b<0.26||nr->shwfit.par_b>2.4)) fTreePID.fCutPID=0;

      if((nr->shwfit.uv_rms<1.28||nr->shwfit.uv_rms>3.28) ||
         (nr->hitcalc.fHitLongEnergy<35.0||nr->hitcalc.fHitLongEnergy>313.5)||
         (nr->shwfit.par_b<0.26||nr->shwfit.par_b>1.80)||
         (nr->shwfit.uv_kurt<21.67)||
         (nr->shwfit.shwmaxplane>44.5)||
         (nr->angclusterfit.fACluFitKurt<124.37)||
         (nr->srshower.planes>14.5)) fTreePID.fCutPID1=0;


     if((nr->srevent.planes<6.5||nr->srevent.planes>14.5) ||
         (nr->hitcalc.fHitLongEnergy<54.7||nr->hitcalc.fHitLongEnergy>298.2)||
         (nr->shwfit.uv_molrad_vert<1.71||nr->shwfit.uv_molrad_vert>6.06)||
         (nr->shwfit.par_b<0.32||nr->shwfit.par_b>1.31)) fTreePID.fCutPID2=0;


      if(!fTreePID.fPassBCuts||!fTreePID.fCutPID) fTreePID.fCutClass=0;
      if(!fTreePID.fPassBCuts||!fTreePID.fCutPID1) fTreePID.fCutClass1=0;
      if(!fTreePID.fPassBCuts||!fTreePID.fCutPID2) fTreePID.fCutClass2=0;

}
  
void TreePIDAna::Reset()
{}
