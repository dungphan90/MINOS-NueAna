#ifndef NUEDATA_H
#define NUEDATA_H                                                               

#include <vector>
#include <map>
#include "Conventions/Detector.h"
#include "Conventions/BeamType.h"
#include "Conventions/ReleaseType.h"
#include "NueAna/NueMini.h"
#include "NueAna/NueMiniAna.h"

#include "NueAna/NueAnaTools/Selection.h"

using namespace std;                                   

class NueRecord;
class NueHeader;

class NueData
{
   public:
     NueData();
     NueData(BeamType::BeamType_t beam,   
             Detector::Detector_t det,
             ReleaseType::Release_t rel,
             bool isNueData=true);

     ~NueData();

     void AddEvent(NueMini *nm);
     void AddEvent(NueRecord *nr);
     void FillRecord(NueRecord *nr, int i);
     void SetupNueHeader(NueHeader &nh);

     double GetNeugenStdXsec(int i);
     void   SetNeugenStdXsec(double xsec, int i);


     BeamType::BeamType_t    GetBeamType()    const;
     Detector::Detector_t    GetDetector()    const; 
     ReleaseType::Release_t  GetRelease()     const;
     Selection::Selection_t  GetSelectionLevel() const;
     bool                    IsNueData()      const;
     bool                    IsData()         const;
     bool                    IsMR()           const;
     void                    SetIsMR(bool);

     void     SetPOT(double);
     void     AddPOT(double);
     double   GetPOT();

     int      NumberOfEntries() const;
     void     Clear();

   private:
     NueMini *fNM;
     NueMiniAna fNMA;

     BeamType::BeamType_t fBeam;
     Detector::Detector_t fDet;
     ReleaseType::Release_t fRelease;
     Selection::Selection_t fSelectionLevel;
     bool kIsNueData;   // True for Nue Selection, False -> CC selection
     bool kIsMR;        // takes the place of FoundMR

     double fPOT;

     vector<double> trkRecoCCEnergy;
     vector<double> shwRecoCCEnergy;
     vector<double> evtRecoNueEnergy;
     vector<double> evtRecoMEUEnergy;
                                                                                
     vector<double> skzpWeight;
     vector<double> weight;

     // these can go once I fix the Preselection Cuts                 
     vector<int> trkPlanes;
     vector<int> trkEndPlane;
     vector<int> trkBegPlane;
     vector<int> trkLikePlanes;
     vector<int> nshower;
     vector<int> contPlanes;


/*
     vector<double> vtxX;
     vector<double> vtxY;
     vector<double> vtxZ;
     vector<double> triggerTime;
     vector<double> spillType;
     vector<double> coilCurrent;
     vector<double> liTime;
     vector<double> eventTimeMin;
     vector<double> rcBoundary;
     vector<int> daveFDDataQuality;     
     vector<double> goodBeamMon;
     vector<double> bmonTimeDiff;
*/

     //MC Truth Values

     vector<double> shiEpi0;
     vector<double> shiEmEnergy;
     vector<double> nuEnergy;
     vector<int> nuFlavor;
     vector<int> nonOscNuFlavor;

     vector<double> nueOscProb;
     vector<int> nueClass;
//     vector<double> emShowerFraction;

//  If I get rid of the Neugen X-Sec Reweighting these can go
     vector<double> ParentType;
     vector<int> interactionType;
     vector<double> nuDCosX;
     vector<double> nuDCosY;
     vector<double> nuDCosZ;
//     vector<double> leptonMom;
     vector<double> hadronicY;
     vector<int> hadronicFinalState;
     vector<double> w2;
     vector<double> q2;
     vector<double> bjorkenX;
     vector<double> targetPX;
     vector<double> targetPY;
     vector<double> targetPZ;
     vector<double> targetEnergy;
     vector<double> atomicNumber;
     vector<double> atomicWeight;
     vector<int> initialState;
     vector<int> resonanceCode;
// End of X-Sec section

     //PID distributions

     vector<double> ann2pe_daikon04;
     vector<double> ann2pe;
     vector<double> ann30;
     vector<double> ann6;
     vector<double> ssPID;
     vector<double> mcnnPID;
     vector<int> mcnnMatch;
     vector<int> cutPID;

     //CC pid information 
     vector<double> abCCPID;
     vector<double> roCCPID;
     vector<int> ntrack;
     vector<int> trkPass;
     vector<int> endPlaneU;
     vector<int> endPlaneV;
     vector<int> deltaUVVtx;

     vector<double> mri_roCCPID;
     vector<double> mri_abCCPID;
     vector<int> mri_trkPass;
     vector<int> gapPlanes;
                                                                                
     vector<double> neugenStdXsec;
     vector<int> cosmicCut;
     vector<int> largestEvent;
};

inline BeamType::BeamType_t    NueData::GetBeamType()  const {return fBeam;}
inline Detector::Detector_t    NueData::GetDetector()  const {return fDet;}
inline ReleaseType::Release_t  NueData::GetRelease()   const {return fRelease;}
inline Selection::Selection_t  NueData::GetSelectionLevel() const {return fSelectionLevel;}

inline bool    NueData::IsNueData()  const {return kIsNueData;}
inline bool    NueData::IsData()     const {return ReleaseType::IsData(fRelease);}
inline bool    NueData::IsMR()       const {return kIsMR; }
inline void    NueData::SetIsMR(bool in)   { kIsMR = in;}
inline void    NueData::SetPOT(double in) {fPOT  = in;} 
inline void    NueData::AddPOT(double in) {fPOT += in;}
inline double  NueData::GetPOT()          {return fPOT;}
inline int   NueData::NumberOfEntries() const {return evtRecoNueEnergy.size();}


#endif //NUEDATA_H
