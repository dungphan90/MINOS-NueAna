#ifndef MUONREMOVALINFOANA_H
#define MUONREMOVALINFOANA_H

#include "NueAna/NueAnaBase.h"
#include "NueAna/MuonRemovalInfo.h"
#include "MuonRemoval/NtpMRRecord.h"
#include "Mad/MadDpID.h"
#include "Mad/MadNsID.h"
#include "Mad/MadAbID.h"
#include "PhysicsNtuple/Store/Interface.h"


class NtpStRecord;

class MuonRemovalInfoAna : public NueAnaBase
{

 public:

  MuonRemovalInfoAna(MuonRemovalInfo &sv);
  virtual ~MuonRemovalInfoAna();

  void Analyze(int evtn, RecRecordImp<RecCandHeader> *mrobj);
  void Analyze(int evtn, RecRecordImp<RecCandHeader> *mrobj,
	       RecRecordImp<RecCandHeader> *oldstobj);
  void Analyze(int evtn, NtpMRRecord *mrobj, NtpStRecord *oldstobj);

  float LoadROPID(int evtn);

 private:
  MuonRemovalInfo &fMuonRemovalInfo;
  
   static MadDpID dpid;
   static MadNsID nsid;
   static MadAbID abid;
                                                                                
    static bool readabidfile;

};

#endif// MUONREMOVALINFOANA_H
