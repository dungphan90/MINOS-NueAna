#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "Conventions/SimFlag.h"
#include "Conventions/Munits.h"
#include "MessageService/MsgService.h"
#include "BeamData/ana/Summary/BeamSummary.h"
#include "NueAna/BeamMonAna.h"
#include "NueAna/BeamMon.h"
#include <SpillTiming/SpillTimeFinder.h>

#include <BeamDataUtil/BeamMonSpill.h>
#include <BeamDataUtil/BDSpillAccessor.h>

CVSID("$ID:");


BeamMonAna::BeamMonAna(BeamMon &bmon):
  fBmon(bmon),
  bs(0)
{}

BeamMonAna::~BeamMonAna()
{}

void BeamMonAna::SetBeamSummary(BeamSummary *b)
{
    if (b) {
	bs=b;
	return;
    }
}

//void BeamMonAna::Analyze(int /*evtn*/, NtpSRRecord *srobj, 
//			 NtpMCRecord */*mc*/, NtpTHRecord */*th*/)
void BeamMonAna::Analyze(int /*evtn*/, RecRecordImp<RecCandHeader> *srobj)
{
  if(srobj==0){
    return;
  }
  if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
     ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
    return;
  }

  
  VldContext vc=srobj->GetHeader().GetVldContext();
  if(vc.GetSimFlag()!=SimFlag::kData){
    return;
  }


  VldTimeStamp vts = vc.GetTimeStamp();
  Int_t tsec = vts.GetSec();
  Int_t tnsec = vts.GetNanoSec();

  VldTimeStamp stnd = SpillTimeFinder::Instance().GetTimeOfNearestSpill(vc);

  if (bs) {
      bs->MatchSpillWithTime(tsec,tnsec);
      
      fBmon.bI=bs->beamIntensity;
      fBmon.hbw=bs->hBeamWidth;
      fBmon.vbw=bs->vBeamWidth;
      fBmon.hpos1=bs->hPosatTargetPM;
      fBmon.vpos1=bs->vPosatTargetPM;
      fBmon.hpos2=bs->hPosatTargetBPM;
      fBmon.vpos2=bs->vPosatTargetBPM;
      fBmon.htan=bs->tanHoriz;
      fBmon.vtan=bs->tanVert;
      fBmon.hornI=bs->hornPeakCurrent;
      fBmon.nuTarZ=bs->nuTarZ;
      fBmon.time=bs->timeStampD;
      fBmon.dt_bmst = ((double)(vts)) - fBmon.time;
      fBmon.dt_stnd = (double)(vts-stnd);
      return;
  }

  const BeamMonSpill* spill = BDSpillAccessor::Get().LoadSpill(vts);

      
  if (!spill) {
    MSG("BeamMonAna",Msg::kError)
	<< "No BeamMonSpill for " << vts << endl;
    return;
  }

  double tor = spill->fTortgt;
  if (tor == 0) tor = spill->fTrtgtd;
  if (tor == 0) tor = spill->fTr101d;
  if (tor == 0) tor = spill->fTor101;
  fBmon.bI = tor;

  fBmon.hbw = spill->fProfWidX;
  fBmon.vbw = spill->fProfWidY;

  fBmon.hpos1 = spill->fTargProfX;
  fBmon.vpos1 = spill->fTargProfY;

  double xbpm,ybpm,xrms,yrms;
  spill->BpmAtTarget(xbpm,ybpm,xrms,yrms);
  fBmon.hpos2 = xbpm;
  fBmon.vpos2 = ybpm;

  fBmon.htan = 0;
  fBmon.vtan = 0;

  fBmon.hornI = spill->fHornCur;

  fBmon.nuTarZ = 0;
  int beam_type = spill->GetStatusBits().beam_type;
  if (beam_type==4)
      fBmon.nuTarZ = 1.50 * Munits::meter;
  if (beam_type==5)
      fBmon.nuTarZ = 2.50 * Munits::meter;

  VldTimeStamp st = spill->SpillTime();
  fBmon.time = (double)st;
  fBmon.dt_stnd = (double)(vts-stnd);
  fBmon.dt_bmst = (double)(vts-st);
  
  MSG("BeamMonAna",Msg::kDebug)
      << vts << " - " << st << " = " << (double)(vts-st)
      << " : " << tor << " (" << fBmon.hpos1 << " x " << fBmon.vpos1 << ") "
      << "(" << xbpm << " x " << ybpm << ") " << fBmon.hornI << " " << beam_type << endl;

}

