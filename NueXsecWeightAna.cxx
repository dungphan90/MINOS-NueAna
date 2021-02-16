#include "StandardNtuple/NtpStRecord.h"
#include "MessageService/MsgService.h"
#include "Registry/Registry.h"
#include "MCReweight/MCReweight.h"
#include "MCReweight/MCEventInfo.h"
#include "MCReweight/ReweightHelpers.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/NueXsecWeightAna.h"

CVSID("$Id: NueXsecWeightAna.cxx,v 1.6 2007/03/12 14:05:34 cbs Exp $");


NueXsecWeightAna::NueXsecWeightAna(NueXsecWeight &nuexs):
  fNueXsecWeight(nuexs)
{}

NueXsecWeightAna::~NueXsecWeightAna()
{}

void NueXsecWeightAna::Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj)
{
  NtpStRecord *st=dynamic_cast<NtpStRecord *>(srobj);
  if(st==0){
    MSG("NueXsecWeightAna",Msg::kError)<<"Trying to do mc reweighting on an event"
					<<" that comes from a non NtpStRecord"<<endl
					<<"That's not a good idea"
				       <<"I haven't even tried to implement it"<<endl;
    return;
  }
  if(st->evthdr.nevent == 0 ) return;

  double gw = 1.0;

  if(ReleaseType::IsCarrot(release)){
    MCEventInfo ei;

    int thn = SntpHelpers::GetEvent2MCIndex(evtn,st);
    if(thn<0) {
      fNueXsecWeight.xsecweight=1.0;
      return;
    }
    
    ReweightHelpers::MCEventInfoFilla(&ei,st,thn);
    NuParent *np=0;
    if(ei.iresonance!=1005){
      gw = mcr->ComputeWeight(&ei,np,rwtreg);
    }
  }
  
  fNueXsecWeight.xsecweight=gw;
}
