#ifndef SNTPHELPERS_H
#define SNTPHELPERS_H
#include <vector>
#include "Record/RecCandHeader.h"
#include "Record/RecRecordImp.h"

class TObject;
class NtpSRRecord;
class NtpMCRecord;
class NtpTHRecord;
class NtpSREvent;
class NtpSRTrack;
class NtpSRStrip;
class NtpSRShower;
class NtpSRCluster;
class NtpSRSlice;
class NtpMCTruth;
class NtpMCStdHep;
class NtpStRecord;
class NtpSRShieldSummary;
class NtpSRShieldStrip;
class NtpMRRecord;
class NtpMREvent;
class NtpMRTruth;

namespace SntpHelpers
{

  NtpSREvent *GetEvent(int evnt, RecRecordImp<RecCandHeader> *rri);
  NtpSRTrack *GetTrack(int trkn, RecRecordImp<RecCandHeader> *rri);
  NtpSRShower *GetShower(int shwn, RecRecordImp<RecCandHeader> *rri);
  NtpSRCluster *GetCluster(int clun, RecRecordImp<RecCandHeader> *rri);
  NtpSRStrip *GetStrip(int stpn, RecRecordImp<RecCandHeader> *rri);
  NtpSRShieldSummary *shieldSummary(RecRecordImp<RecCandHeader> *rri);
  NtpSRShieldStrip *GetShieldStrip(int stpn, RecRecordImp<RecCandHeader> *rri);  
  NtpSREvent *GetEvent(int evnt, NtpSRRecord *sr);
  NtpSRTrack *GetTrack(int trkn, NtpSRRecord *sr);
  NtpSRShower *GetShower(int shwn, NtpSRRecord *sr);
  NtpSRCluster *GetCluster(int clun, NtpSRRecord *sR);
  NtpSRSlice *GetSlice(int slcn, NtpSRRecord *sr);
  NtpSRStrip *GetStrip(int stpn, NtpSRRecord *sr);
  NtpSRShieldSummary *shieldSummary(NtpSRRecord *sr);
  NtpSRShieldStrip *GetShieldStrip(int stpn, NtpSRRecord *sr); 
  NtpSRShower *GetPrimaryShower(int evnt, RecRecordImp<RecCandHeader> *rri);

//  NtpSRShower *GetPrimaryShower(NtpSREvent* evt, RecRecordImp<RecCandHeader> *rri);
  
  NtpSREvent *GetEvent(int evnt, NtpStRecord *st);
  NtpSRTrack *GetTrack(int trkn, NtpStRecord *st);
  NtpSRTrack *GetPrimaryTrack(int evtn,NtpStRecord *st);
  NtpSRShower *GetShower(int shwn, NtpStRecord *st);
  NtpSRCluster *GetCluster(int clun, NtpStRecord *st);
  NtpSRSlice *GetSlice(int slcn, NtpStRecord *st);
  NtpSRStrip *GetStrip(int stpn, NtpStRecord *st);
  NtpSRShieldSummary *shieldSummary(NtpStRecord *st);
  NtpSRShieldStrip *GetShieldStrip(int stpn, NtpStRecord *st);
    
  int GetStripIndex(int stpn, TObject *obj);
  int GetTrackIndex(int trkn, NtpSREvent *event);
  int GetShowerIndex(int shwn, NtpSREvent *event);
  int GetClusterIndex(int clun,NtpSRShower *shower);
  
  NtpMCTruth *GetMCTruth(int evnt, NtpMCRecord *mc);
  NtpMCTruth *GetMCTruth(int evnt, NtpStRecord *st);
  std::vector<NtpMCStdHep *> GetStdHepArray(int evnt, NtpMCRecord *mc);
  std::vector<NtpMCStdHep *> GetStdHepArray(int evnt, NtpStRecord *st);
  
  int GetEvent2MCIndex(int evnt, NtpTHRecord *th);
  int GetEvent2MCIndex(int evnt, NtpStRecord *st);
 
  Int_t InPartialRegion(UShort_t plane, UShort_t strip); 

  NtpMREvent *GetMREvent(int ind,NtpMRRecord *mr);
  NtpMRTruth *GetMRTruth(int ind,NtpMRRecord *mr);
  void FillEventEnergy(float* ph0, float* ph1, int evtn, NtpStRecord* st, const
int SIZE);
}

#endif// SNTPHELPERS_H
