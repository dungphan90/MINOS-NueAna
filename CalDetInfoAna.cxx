#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "NueAna/CalDetInfoAna.h"
#include "CalDetDST/UberRecord.h"
#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

CalDetInfoAna::CalDetInfoAna(CalDetInfo &sv):
  fCalDetInfo(sv)
{
}

CalDetInfoAna::~CalDetInfoAna()
{}

void CalDetInfoAna::Analyze(RecRecordImp<UberRecHeader> *uberrecord)
{
  fCalDetInfo.Reset();
  
  if(uberrecord==0) return;  
  UberRecord *ur=0;
  if(((ur=dynamic_cast<UberRecord *>(uberrecord))==0)) return;
    
  fCalDetInfo.Zero();
  
  fCalDetInfo.beamp = ur->GetHeader().GetBeamMomentum();
  fCalDetInfo.inct = ur->cpid.inct;
  fCalDetInfo.pid = ur->cpid.pid;
  fCalDetInfo.olchi2 = ur->cpid.olchi2;
  fCalDetInfo.p0stripmaxmip = ur->p0stripmaxmip;
  
}

			  
