#include "MessageService/MsgService.h"
#include "MCNtuple/NtpMCTruth.h"
#include "StandardNtuple/NtpStRecord.h"
#include "MCReweight/MuParentHelper.h"
#include "MCReweight/NuParent.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/MCFluxInfoAna.h"

CVSID("$Id: MCFluxInfoAna.cxx,v 1.5 2007/03/01 16:38:49 rhatcher Exp $");

MCFluxInfoAna::MCFluxInfoAna(NtpMCFluxInfo &fluxinfo):
  fMCFluxInfo(fluxinfo),
  fluxfileneedsdebuggering(true),
  mupar(0)
{}

MCFluxInfoAna::~MCFluxInfoAna()
{
  MSG("MCFluxInfoAna",Msg::kDebug)<<"Destructing MCFluxInfoAna"<<endl;
  
  MSG("MCFluxInfoAna",Msg::kDebug)<<"Leaving ~MCFluxInfoAna"<<endl;
}

void MCFluxInfoAna::Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj)
{
  //only NtpStRecords have the flux info
  NtpStRecord *st = dynamic_cast<NtpStRecord *>(srobj);
  if(!st){
    MSG("MCFluxInfoAna",Msg::kError)<<"Trying to get flux info from a non NtpStRecord"<<endl
			<<"That's not a good idea"<<endl;
    return;
  }
  int thn = SntpHelpers::GetEvent2MCIndex(evtn,st);
  if(thn<0){
    return;
  }
  NtpMCTruth *mcrec = SntpHelpers::GetMCTruth(thn,st);
  MSG("MCFluxInfoAna",Msg::kDebug)<<"in MCFluxInfoAna::Analyze mcrec: "<<mcrec<<endl;
  if(mcrec==0){
      MSG("MCFluxInfoAna",Msg::kError)<<"No mcrec! cannot fill flux info: "<<mcrec<<endl;
      return;
  }
  Analyze(mcrec);
}

void MCFluxInfoAna::Analyze(NtpMCTruth *mcrec)
{

  if(mcrec==0){
    return;
  }

  MSG("MCFluxInfoAna",Msg::kDebug)<<"in analyze(NtpMCTruth *)"<<endl;
  MSG("MCFluxInfoAna",Msg::kDebug)<<"flux info before copying: "
  				  <<" tpx: "<<mcrec->flux.tpx<<" "
    				  <<" tpy: "<<mcrec->flux.tpy<<" "
    				  <<" tpz: "<<mcrec->flux.tpz<<endl;
  
  fMCFluxInfo=mcrec->flux;
  MSG("MCFluxInfoAna",Msg::kDebug)<<"copying flux info: "
  				  <<" tpx: "<<fMCFluxInfo.tpx<<" "
    				  <<" tpy: "<<fMCFluxInfo.tpy<<" "
    				  <<" tpz: "<<fMCFluxInfo.tpz<<endl;

  if(fluxfileneedsdebuggering){
    ResetMuParentInfo();
  }

}

void  MCFluxInfoAna::ResetMuParentInfo()
{
  if(!mupar){
    MSG("MCFluxInfoAna",Msg::kError)<<"If you want to reset the muparent info, you must provide a mu parent helper"<<endl;
    return;
  }
  NuParent par;
  if(fMCFluxInfo.tptype==-13 || fMCFluxInfo.tptype==13) {
    mupar->GetMuParent(fMCFluxInfo.fluxrun,fMCFluxInfo.fluxevtno,
		      fMCFluxInfo.tpx,fMCFluxInfo.tpy,fMCFluxInfo.tpz,par);
  

    fMCFluxInfo.tvx = par.GetX();
    fMCFluxInfo.tvy = par.GetY();
    fMCFluxInfo.tvz = par.GetZ();
    fMCFluxInfo.tpx = par.GetPx();
    fMCFluxInfo.tpy = par.GetPy();
    fMCFluxInfo.tpz = par.GetPz();
    fMCFluxInfo.tptype = par.GetPID();
    
    fMCFluxInfo.tgen = par.GetGen();
  }
}
