#ifndef NueMiniPID_H
#define NueMiniPID_H

#include <ctime>
#include <vector>
#include <map>
#include "Conventions/Detector.h"
#include "Conventions/BeamType.h"
#include "Conventions/ReleaseType.h"

#include "NueAna/NueAnaTools/Selection.h"

#include "TObject.h"

using namespace std;                                   

class NueRecord;

class NueMiniPID: public TObject
{
   public:
     NueMiniPID();
     NueMiniPID(BeamType::BeamType_t beam,
             Detector::Detector_t det,
             ReleaseType::Release_t rel);
     virtual ~NueMiniPID();

     BeamType::BeamType_t fBeam;
     Detector::Detector_t fDet;
     ReleaseType::Release_t fRelease;
     Selection::Selection_t fSelectionLevel;

	 //time stamp
     time_t timestamp;


     double fPOT;

     double trkRecoCCEnergy;
     double shwRecoCCEnergy;
     double evtRecoNueEnergy;
     double evtRecoMEUEnergy;
     double skzpWeight;
     double MCWeight;    
 
     double weight;

     //MC Variables

     double shiEpi0;
     double shiEmEnergy;
     double nuEnergy;
     double nueOscProb;
     double ParentType;     
     double nuDCosX;
     double nuDCosY;
     double nuDCosZ;
     double hadronicY;
     double w2;
     double q2;
     double bjorkenX;
     double targetPX;
     double targetPY;
     double targetPZ;
     double targetEnergy;
     double atomicNumber;
     double atomicWeight;
     double neugenStdXsec;

     // PID Variables
     double ann14;
     double ann2pe_daikon04;
     double ann2pe;
     double ann30;
     double ann6;
     double ssPID;
     double mcnnPID;
     double abCCPID;
     double roCCPID;

     double mri_abCCPID;
     double mri_roCCPID;


     double longest_s;
     double event_length;
     double event_energy;
     double particle_energy;

     //pids
     double pidA;
     double pidB;
     double pidC;
     double pidD;
     double pidE;
     double pidF;

     //mrcc
     double mrcc_s;
     double pars[14];

     double nueVtxX;
     double nueVtxY;
     double nueVtxZ;
     double vtxU;
     double vtxV;
     double vtxZ;


	 //nue mrcc vars
	 double mri_qp;
	 double mri_orig_cc_pid;
	 double mri_SigmaQP;
	 
	 float mri_best_complete;
     int mri_fitp;
    
    
    
     // Preselection variables
     int nshower;
     int contPlanes;
     int cosmicCut;
     int largestEvent;

     int trkPlanes;
     int trkEndPlane;
     int trkBegPlane;
     int trkLikePlanes;



     int hadronicFinalState;
     int initialState;
     int resonanceCode;
     int nuFlavor;
     int nonOscNuFlavor;
     int nueClass;
     int interactionType;
     int mcnnMatch;
     int cutPID;


     int ntrack;
     int trkPass;
     int endPlaneU;
     int endPlaneV;
     int deltaUVVtx;

     //MRCC variables
     int mri_trkPass;
     int gapPlanes;
     
     
     
     //PID variables
     //preselection
     int ntot;


     //nuecut summary
     //event header variables
     int run;
     int subrun;
     int event;
     int snarl;


     //PID variables
     //preselection
     bool infid;
     bool contained;
     
     //NueStandard cuts...
     bool passes_NueStandard_PassesDataQuality; 
	 bool passes_NueStandard_IsInFid;     
     bool passes_NueStandard_PassesPOTStandards;
     bool passes_NueStandard_PassesCosmicCut;
     bool passes_NueStandard_PassesNonHEPreSelection;
     bool passes_NueStandard_PassesPreSelection;    
   
  	 bool passes_NueStandard_PassesMinPlaneCut;
     bool passes_NueStandard_PassesShowerCut;
     bool passes_NueStandard_PassesTrackPlaneCut;
     bool passes_NueStandard_PassesTrackLikePlaneCut;
     bool passes_NueStandard_PassesLowEnergyCut;
     bool passes_NueStandard_PassesHighEnergyCut;

     bool passes_NueStandard_PassesMRCCFiducial;
     bool passes_NueStandard_PassesMRCCPreSelection;  



	 //NN variable bounds check    
     bool pass_var_check;
     bool pass_nvar_check;




private:
    ClassDef(NueMiniPID, 6)
};

#endif //NueMiniPID_H
