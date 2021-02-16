#ifndef NUEMINI_H
#define NUEMINI_H                                                               

#include <vector>
#include <map>
#include "Conventions/Detector.h"
#include "Conventions/BeamType.h"
#include "Conventions/ReleaseType.h"

#include "NueAna/NueAnaTools/Selection.h"

#include "TObject.h"

using namespace std;                                   

class NueRecord;

class NueMini: public TObject
{
   public:
     NueMini();
     NueMini(BeamType::BeamType_t beam,
             Detector::Detector_t det,
             ReleaseType::Release_t rel);
     virtual ~NueMini();


	 ///////////////////
	 ////!!!!!!
	 ////
	 ////  WARNING!... the items in this list are structured to optimize the filesize (and thus speed)
	 ////
	 ////  when adding new variables, group them in order of bytesize for a given type...
	 ////     so doubles, ints, bools  .....
	 ////
	 ////  for more info, see
	 ////
	 ////	http://root.cern.ch/root/roottalk/roottalk01/3900.html
	 ////
	 ////   --Steve Cavanaugh 
	 ////
	 ///////////////////
	 

	 /////MINOS Types...

     BeamType::BeamType_t fBeam;
     Detector::Detector_t fDet;
     ReleaseType::Release_t fRelease;
     Selection::Selection_t fSelectionLevel;


	 //maybe not in use?
     double fPOT;

	 //user placeholders for later use
	 double weight;

	 //reco energies
     double trkRecoCCEnergy;
     double shwRecoCCEnergy;
     double evtRecoNueEnergy;
     double evtRecoMEUEnergy;


     //MC Variables
     //MC base variables
     double shiEpi0;
     double shiEmEnergy;
     double nuEnergy;
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
     
     //MC calculated variables
     double nueOscProb;
     double skzpWeight;
    

     // PID Variables
     double annpid_11inp;
     double annpid_11inp_daikon04;
     double annpid_14inp_daikon04;
     double annpid;   
     double ssPID;
     double mcnnPID;
     double abCCPID;
     double roCCPID;



     //MRCC variables
     double mri_orig_abCCPID;
     double mri_orig_roCCPID;
     double mri_best_complete_phw;
     double mri_best_purity_phw;
     double mri_shwe;
     double mri_pmux;
     double mri_pmuy;
     double mri_pmuz;
     double mri_vtxx;
     double mri_vtxy;
     double mri_vtxz;





     //ANN 11 Variables:
     double shwfit_par_a;
     double shwfit_par_b;
     double shwfit_uv_molrad_peak_9s_2pe_dw;
     double shwfit_uv_rms_9s_2pe_dw;
     double mstvars_e4w_PLUS_mstvars_o4w;
     double fracvars_fract_road;
     double fracvars_fract_2_planes;
     double fracvars_fract_4_planes;
     double fracvars_fract_6_planes;
     double fracvars_fract_8_counters;
     double shwfit_LongE;

     //ANN 14 new added Variables:
     double fracvars_shw_max;
     double fracvars_shw_slp;
     double fracvars_dis2stp;
     double fracvars_fract_1_plane;
     double fracvars_fract_5_planes;
     double fracvars_fract_6_counters;
     double fracvars_fract_20_counters; 
     double event_length;

     //LEM variables:
     double mcnnv_fracCC;
     double mcnnv_ymean;
     double mcnnv_meanfracQmatched; 
     double mcnnv_var2;
     double mcnnv_fracCCy;
     double mcnnv_qtot;

     //SSPID variables:
     double ssvar1;
     double ssvar2;
     double ssvar3;
     double ssvar4;

     //event vertex:
     double vtxX;
     double vtxY;
     double vtxZ;

     //time stamp
     time_t timestamp;


     //MC Variables
     int nuFlavor;
     int nonOscNuFlavor;
     int nueClass;
     int hadronicFinalState;
     int interactionType;
     int initialState;
     int resonanceCode;
     
     
     // Preselection variables
     int nshower;
     int contPlanes;
     int cosmicCut;
     int largestEvent;

     int trkPlanes;
     int trkEndPlane;
     int trkBegPlane;
     int trkLikePlanes;

     // PID Variables
     int ntrack;
     int trkPass;
     int endPlaneU;
     int endPlaneV;
     int deltaUVVtx;
     int mcnnMatch;

     //MRCC variables
     int mri_trkPass;
     int gapPlanes;
     
     //event header variables
     int run;
     int subrun;
     int event;
     int snarl;
     
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












private:
    ClassDef(NueMini, 9)
};

#endif //NUEMINI_H
