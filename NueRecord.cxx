#include "MessageService/MsgService.h"
#include "NueAna/NueRecord.h"




ClassImp(NueRecord)
CVSID("$Id: NueRecord.cxx,v 1.34 2009/06/23 18:08:31 scavan Exp $");

NueRecord::NueRecord():
    RecRecordImp<NueHeader>(),
    shwfit(),
//   reco(),
    hitcalc(),
    angcluster(),
    angclusterfit(),
//   vertfind(),
//   vtx(),
//   emvars(),
    mstvars(),
    fracvars(),
    subshowervars(),
    highhitvars(),
    shieldrejvars(),
    ann(),
    anainfo(),
    srevent(),
    srshower(),
    srtrack(),
    mctrue(),
    bmon(),
    mdadiscrim(),
    treepid(),
    fluxinfo(),
    fluxweights(),
    xsecweights(),
    shi(),
    mri(),
    cdi(),
    timingvars(),
    mcnnv(),
    dtree(),
    eventq(),
    precord(),
    precordMRCC()
{
//   SetClearable(true);
   MSG("NueRecord",Msg::kDebug)<<"In NueRecord()"<<endl;
}

NueRecord::NueRecord(const NueHeader& head):
    RecRecordImp<NueHeader>(head),
    shwfit(),
//   reco(),
    hitcalc(),
    angcluster(),
    angclusterfit(),
//   vertfind(),
//   vtx(),
//   emvars(),
    mstvars(),
    fracvars(),
    subshowervars(),
    highhitvars(),
    shieldrejvars(),
    ann(),
    anainfo(),
    srevent(),
    srshower(),
    srtrack(),
    mctrue(),
    bmon(),
    mdadiscrim(),
    treepid(),
    fluxinfo(),
    fluxweights(),
    xsecweights(),
    shi(),
    mri(),
    cdi(),
    timingvars(),
    mcnnv(),
    dtree(),
    eventq(),
    precord(),
    precordMRCC()
{
//   SetClearable(true);
   MSG("NueRecord",Msg::kDebug)<<"In NueRecord(const NueHeader &)"<<endl;
}

NueRecord::NueRecord(const NueRecord &nr):
    RecRecordImp<NueHeader>(nr.GetHeader()),
    shwfit(nr.shwfit),
    hitcalc(nr.hitcalc),
    angcluster(nr.angcluster),
    angclusterfit(nr.angclusterfit),
    mstvars(nr.mstvars),
    fracvars(nr.fracvars),
    subshowervars(nr.subshowervars),
    highhitvars(nr.highhitvars),
    shieldrejvars(nr.shieldrejvars),
    ann(nr.ann),
    anainfo(nr.anainfo),
    srevent(nr.srevent),
    srshower(nr.srshower),
    srtrack(nr.srtrack),
    mctrue(nr.mctrue),
    bmon(nr.bmon),
    mdadiscrim(nr.mdadiscrim),    
    treepid(nr.treepid),
    fluxinfo(nr.fluxinfo),
    fluxweights(nr.fluxweights),
    xsecweights(nr.xsecweights),
    shi(nr.shi),
    mri(nr.mri),
    cdi(nr.cdi),
    timingvars(nr.timingvars),
    mcnnv(nr.mcnnv),
    dtree(nr.dtree),
    eventq(nr.eventq),
    precord(nr.precord),
    precordMRCC(nr.precordMRCC)
{
   MSG("NueRecord",Msg::kDebug)<<"In Copy constructor NueRecord"<<endl;
}

void NueRecord::Clear(Option_t* /* option */)
{
   Reset();
   fluxweights.Clear();
   shwfit.Clear();
}

void NueRecord::Reset()
{
    //Need some sort of header reset
    shwfit.Reset();
    hitcalc.Reset();
    angcluster.Reset();
    angclusterfit.Reset();
    mstvars.Reset();
    fracvars.Reset();
    subshowervars.Reset();
    highhitvars.Reset();
    ann.Reset();
    anainfo.Reset();
    srevent.Reset();
    srshower.Reset();
    srtrack.Reset();
    mctrue.Reset();
    bmon.Reset();
    mdadiscrim.Reset();
    treepid.Reset();
    fluxweights.Reset();
    xsecweights.Reset();
    shi.Reset();
    mri.Reset();
    cdi.Reset();
    timingvars.Reset();
    mcnnv.Reset();
    dtree.Reset();
    eventq.Reset();
    precord.Reset();
    precordMRCC.Reset();
}

NueRecord::~NueRecord()
{
   MSG("NueRecord",Msg::kDebug)<<"In ~NueRecord"<<endl;
}

