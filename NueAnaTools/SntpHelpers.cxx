#include "TClonesArray.h"
#include "MCNtuple/NtpMCRecord.h"
#include "MCNtuple/NtpMCTruth.h"
#include "MCNtuple/NtpMCStdHep.h"
#include "TruthHelperNtuple/NtpTHRecord.h"
#include "TruthHelperNtuple/NtpTHEvent.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "CandNtupleSR/NtpSRCluster.h"
#include "CandNtupleSR/NtpSRSlice.h"
#include "CandNtupleSR/NtpSRShieldStrip.h"
#include "CandNtupleSR/NtpSRShieldSummary.h"
#include "StandardNtuple/NtpStRecord.h"
#include "Record/RecRecordImp.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "MuonRemoval/NtpMRRecord.h"
#include "MuonRemoval/NtpMREvent.h"
#include "MuonRemoval/NtpMRTruth.h"
#include "MessageService/MsgService.h"

#include <iostream>

CVSID("$Id: SntpHelpers.cxx,v 1.9 2008/09/07 18:21:49 boehm Exp $");

using namespace std;

NtpSREvent *SntpHelpers::GetEvent(int evnt, RecRecordImp<RecCandHeader> *rri)
{
  NtpStRecord *st=dynamic_cast<NtpStRecord *>(rri);

  if(st!=0){
    return GetEvent(evnt,st);
  }

  NtpSRRecord *sr=dynamic_cast<NtpSRRecord *>(rri);
  if(sr!=0){
    return GetEvent(evnt,sr);
  }

  return 0;
}

NtpSRTrack *SntpHelpers::GetTrack(int trkn, RecRecordImp<RecCandHeader> *rri)
{
  NtpStRecord *st=dynamic_cast<NtpStRecord *>(rri);
  if(st!=0){
    return GetTrack(trkn,st);
  }

  NtpSRRecord *sr=dynamic_cast<NtpSRRecord *>(rri);
  if(sr!=0){
    return GetTrack(trkn,sr);
  }
  return 0;
}

NtpSRShower *SntpHelpers::GetShower(int shwn, RecRecordImp<RecCandHeader> *rri)
{
  NtpStRecord *st=dynamic_cast<NtpStRecord *>(rri);
  if(st!=0){
    return GetShower(shwn,st);
  }

  NtpSRRecord *sr=dynamic_cast<NtpSRRecord *>(rri);
  if(sr!=0){
    return GetShower(shwn,sr);
  }
  return 0;
}

NtpSRShower *SntpHelpers::GetPrimaryShower(int evnt, RecRecordImp<RecCandHeader> *rri)
{
   NtpSRShower *shower=0;
   if(rri==0){
      return shower;
   }

   NtpSREvent *event = SntpHelpers::GetEvent(evnt,rri);
/*   return SntpHelpers::GetPrimaryShower(event, rri);
}

NtpSRShower *GetPrimaryShower(NtpSREvent* evt, RecRecordImp<RecCandHeader> *rri)
{
   NtpSRShower *shower=0;
   if(rri==0){
      return shower;
   }
                                                                                
   NtpSREvent *event = evt;
*/
   if (event==0){
       return shower;
   }
                                                                                
   Int_t nshws = event->nshower;
   if (nshws){//at least one shower
     if(event->primshw < 0){
//        cout<<"Event has no primray shower but "<<nshws<<" showers."<<endl;
         return shower;
     }
     int index = event->shw[event->primshw];
//     cout<<"Taking shower at "<<event->primshw<<" points to  "<<index<<endl;
     if(index < 0) return shower;
     shower = SntpHelpers::GetShower(index,rri);
   }
   return shower;
}

NtpSRCluster *SntpHelpers::GetCluster(int clun, RecRecordImp<RecCandHeader> *rri)
{
  NtpStRecord *st=dynamic_cast<NtpStRecord *>(rri);
  if(st!=0){
    return GetCluster(clun,st);
  }

  NtpSRRecord *sr=dynamic_cast<NtpSRRecord *>(rri);
  if(sr!=0){
    return GetCluster(clun,sr);
  }
  return 0;
}

NtpSRStrip *SntpHelpers::GetStrip(int stpn, RecRecordImp<RecCandHeader> *rri)
{
  NtpStRecord *st=dynamic_cast<NtpStRecord *>(rri);
  if(st!=0){
    return GetStrip(stpn,st);
  }

  NtpSRRecord *sr=dynamic_cast<NtpSRRecord *>(rri);
  if(sr!=0){
    return GetStrip(stpn,sr);
  }
    return 0;
}

NtpSRShieldSummary *SntpHelpers::shieldSummary(RecRecordImp<RecCandHeader> *rri)
{
  NtpStRecord *st=dynamic_cast<NtpStRecord *>(rri);
  if(st!=0){
    return shieldSummary(st);
  }

  NtpSRRecord *sr=dynamic_cast<NtpSRRecord *>(rri);
  if(sr!=0){
    return shieldSummary(sr);
  }
    return 0;
}

NtpSRShieldStrip *SntpHelpers::GetShieldStrip(int stpn, RecRecordImp<RecCandHeader> *rri)
{
  NtpStRecord *st=dynamic_cast<NtpStRecord *>(rri);
  if(st!=0){
    return GetShieldStrip(stpn,st);
  }

  NtpSRRecord *sr=dynamic_cast<NtpSRRecord *>(rri);
  if(sr!=0){
    return GetShieldStrip(stpn,sr);
  }
    return 0;
}

NtpSREvent *SntpHelpers::GetEvent(int evnt,NtpSRRecord *sr)
{
   NtpSREvent *event=0;
   if(sr==0){
      return event;
   }
   if(evnt>=sr->evthdr.nevent){
      return event;
   }

   event = dynamic_cast<NtpSREvent *>((*sr->evt)[evnt]);

   return event;
}

NtpSRTrack *SntpHelpers::GetTrack(int trkn,NtpSRRecord *sr)
{
   NtpSRTrack *track=0;
   if(sr==0){
      return track;
   }
   if(trkn>=sr->evthdr.ntrack){
      return track;
   }

   track = dynamic_cast<NtpSRTrack *>((*sr->trk)[trkn]);

   return track;
}

NtpSRShower *SntpHelpers::GetShower(int shwn,NtpSRRecord *sr)
{
   NtpSRShower *shower=0;
   if(sr==0){
      return shower;
   }
   if(shwn>=sr->evthdr.nshower){
      return shower;
   }

   shower = dynamic_cast<NtpSRShower *>((*sr->shw)[shwn]);

   return shower;
}

NtpSRCluster *SntpHelpers::GetCluster(int clun,NtpSRRecord *sr)
{
   NtpSRCluster *cluster=0;
   if(sr==0){
      return cluster;
   }
   if(clun>=sr->evthdr.ncluster){
      return cluster;
   }

   cluster = dynamic_cast<NtpSRCluster *>((*sr->clu)[clun]);

   return cluster;
}

NtpSRSlice *SntpHelpers::GetSlice(int slcn,NtpSRRecord *sr)
{
   NtpSRSlice *slice=0;
   if(sr==0){
      return slice;
   }
   if(slcn>=sr->evthdr.nslice){
      return slice;
   }

   slice = dynamic_cast<NtpSRSlice *>((*sr->slc)[slcn]);

   return slice;
}

NtpSRStrip *SntpHelpers::GetStrip(int stpn,NtpSRRecord *sr)
{
   NtpSRStrip *strip=0;
   if(sr==0){
      return strip;
   }
   if(stpn>=(int)(sr->evthdr.nstrip)){
      return strip;
   }

   strip = dynamic_cast<NtpSRStrip *>((*sr->stp)[stpn]);

   return strip;
}

NtpSRShieldSummary *SntpHelpers::shieldSummary(NtpSRRecord *sr)
{
   NtpSRShieldSummary *shsm=0;
   if(sr==0){
      return shsm;
   }

   shsm = dynamic_cast<NtpSRShieldSummary *>(&sr->vetohdr);

   return shsm;
}

NtpSRShieldStrip *SntpHelpers::GetShieldStrip(int stpn,NtpSRRecord *sr)
{
   NtpSRShieldStrip *strip=0;
   if(sr==0){
      return strip;
   }
   int nshdigt = sr->vetohdr.ndigit[0]+sr->vetohdr.ndigit[1]+sr->vetohdr.ndigit[2];
   if(sr->vetohdr.ishit==0||stpn>=nshdigt){
      return strip;
   }
   strip = dynamic_cast<NtpSRShieldStrip *>((*sr->vetostp)[stpn]);

   return strip;
}

NtpSREvent *SntpHelpers::GetEvent(int evnt,NtpStRecord *st)
{
   NtpSREvent *event=0;
   if(st==0){
      return event;
   }
   if(evnt>=st->evthdr.nevent){
      return event;
   }

   event = dynamic_cast<NtpSREvent *>((*st->evt)[evnt]);

   return event;
}

NtpSRTrack *SntpHelpers::GetTrack(int trkn,NtpStRecord *st)
{
   NtpSRTrack *track=0;
   if(st==0){
      return track;
   }
   if(trkn>=st->evthdr.ntrack){
      return track;
   }


   track = dynamic_cast<NtpSRTrack *>((*st->trk)[trkn]);

   return track;
}

NtpSRTrack *SntpHelpers::GetPrimaryTrack(int evtn,NtpStRecord *st)
{
   NtpSRTrack *ntpTrack=0;
   NtpSRTrack *primaryTrack = 0;
   if(st==0){ return ntpTrack; }
   NtpSREvent* event = SntpHelpers::GetEvent(evtn,st);
  
   int length = 0;
   int tmpLength = 0;

   for(int i = 0; i < event->ntrack; i++){
     int index = event->trk[i];
     ntpTrack = SntpHelpers::GetTrack(index, st);
     if(!primaryTrack) primaryTrack = ntpTrack;
                                                                                
     tmpLength = TMath::Abs(ntpTrack->plane.end - ntpTrack->plane.beg);
                                                                                
     if(tmpLength > length){
       length = tmpLength;
       primaryTrack = ntpTrack;
     }
   }
                             
   return primaryTrack;
}



NtpSRShower *SntpHelpers::GetShower(int shwn,NtpStRecord *st)
{
   NtpSRShower *shower=0;
   if(st==0){
      return shower;
   }
   if(shwn>=st->evthdr.nshower){
      return shower;
   }


   shower = dynamic_cast<NtpSRShower *>((*st->shw)[shwn]);

   return shower;
}

NtpSRCluster *SntpHelpers::GetCluster(int clun,NtpStRecord *st)
{
   NtpSRCluster *cluster=0;
   if(st==0){
      return cluster;
   }
   if(clun>=st->evthdr.ncluster){
      return cluster;
   }

   cluster = dynamic_cast<NtpSRCluster *>((*st->clu)[clun]);

   return cluster;
}

NtpSRSlice *SntpHelpers::GetSlice(int slcn,NtpStRecord *st)
{
   NtpSRSlice *slice=0;
   if(st==0){
      return slice;
   }
   if(slcn>=st->evthdr.nslice){
      return slice;
   }

   slice = dynamic_cast<NtpSRSlice *>((*st->slc)[slcn]);

   return slice;
}

NtpSRStrip *SntpHelpers::GetStrip(int stpn,NtpStRecord *st)
{
   NtpSRStrip *strip=0;
   if(st==0){
      return strip;
   }
   if(stpn>=(int)(st->evthdr.nstrip) || stpn < 0){
      return strip;
   }

   strip = dynamic_cast<NtpSRStrip *>((*st->stp)[stpn]);

   return strip;
}

NtpSRShieldSummary *SntpHelpers::shieldSummary(NtpStRecord *st)
{
   NtpSRShieldSummary *shsm=0;
   if(st==0){
      return shsm;
   }

   shsm = dynamic_cast<NtpSRShieldSummary *>(&st->vetohdr);

   return shsm;
}

NtpSRShieldStrip *SntpHelpers::GetShieldStrip(int stpn,NtpStRecord *st)
{
   NtpSRShieldStrip *strip=0;
   if(st==0){
      return strip;
   }
   int nshdigt = st->vetohdr.ndigit[0]+st->vetohdr.ndigit[1]+st->vetohdr.ndigit[2];
   if(st->vetohdr.ishit==0||stpn>=nshdigt){
      return strip;
   }

   strip = dynamic_cast<NtpSRShieldStrip *>((*st->vetostp)[stpn]);

   return strip;
}

int SntpHelpers::GetStripIndex(int stpn, TObject *obj)
{
   //this may be a stupid way of doing this, but it's the only way I can think of
   //we may want to do this for other reco objects in the future
  if(obj==0){
    return -1;
  }

   const char* cn = obj->ClassName();

   if(strcmp(cn,"NtpSREvent")==0){
      NtpSREvent *event = dynamic_cast<NtpSREvent *>(obj);
      if(stpn>=event->nstrip){
	 return -1;
      }
      else{
	 return event->stp[stpn];
      }
   }
   else if(strcmp(cn,"NtpSRTrack")==0){
      NtpSRTrack *track = dynamic_cast<NtpSRTrack *>(obj);
      if(stpn>=track->nstrip){
	 return -1;
      }
      else{
	 return track->stp[stpn];
      }
   }
   else if(strcmp(cn,"NtpSRShower")==0){
      NtpSRShower *shower = dynamic_cast<NtpSRShower *>(obj);
      if(stpn>=shower->nstrip){
	 return -1;
      }
      else{
	 return shower->stp[stpn];
      }

   }
   else if(strcmp(cn,"NtpSRCluster")==0){
      NtpSRCluster *cluster = dynamic_cast<NtpSRCluster *>(obj);
      if(stpn>=cluster->nstrip){
	 return -1;
      }
      else{
	 return cluster->stp[stpn];
      }

   }else  if(strcmp(cn,"NtpSRSlice")==0){
      NtpSRSlice *slice = dynamic_cast<NtpSRSlice *>(obj);
      if(stpn>=slice->nstrip){
         return -1;
      }
      else{
         return slice->stp[stpn];
      }
   }
   else{
      return -1;
   }

   return -1;
}

int SntpHelpers::GetTrackIndex(int trkn, NtpSREvent *event)
{
   if(trkn>=event->ntrack){
	 return -1;
      }
   else{
      return event->trk[trkn];
   }
   return -1;
}

int SntpHelpers::GetShowerIndex(int shwn, NtpSREvent *event)
{
   if(shwn>=event->nshower){
      return -1;
   }
   else{
      return event->shw[shwn];
   }
   return -1;
}

int SntpHelpers::GetClusterIndex(int clun, NtpSRShower *shower)
{
   if(clun>=shower->ncluster){
      return -1;
   }
   else{
      return shower->clu[clun];
   }
   return -1;
}

int SntpHelpers::GetEvent2MCIndex(int evnt, NtpTHRecord *th)
{
   int thsize=th->thevt->GetEntries();
   if(evnt>=thsize){
      return -1;
   }
   NtpTHEvent *the = dynamic_cast<NtpTHEvent *>((*th->thevt)[evnt]);
   return the->neumc;
   
}

int SntpHelpers::GetEvent2MCIndex(int evnt, NtpStRecord *st)
{
   int thsize=st->thevt->GetEntries();
   if(evnt>=thsize){
      return -1;
   }
   NtpTHEvent *the = dynamic_cast<NtpTHEvent *>((*st->thevt)[evnt]);
   return the->neumc;
   
}

NtpMCTruth *SntpHelpers::GetMCTruth(int index, NtpMCRecord *mc)
{
   NtpMCTruth *mct=0;
   
   if(index>=mc->mc->GetEntries()){
      return mct;
   }
   if(index<0){
     return mct;
   }
   
   mct=dynamic_cast<NtpMCTruth *>((*mc->mc)[index]);
   return mct;

}

NtpMCTruth *SntpHelpers::GetMCTruth(int index, NtpStRecord *st)
{
   NtpMCTruth *mct=0;
   
   if(index>=st->mc->GetEntries()){
      return mct;
   }
   if(index<0){
     return mct;
   }
   
   mct=dynamic_cast<NtpMCTruth *>((*st->mc)[index]);
   return mct;

}

std::vector<NtpMCStdHep *> SntpHelpers::GetStdHepArray(int index, NtpMCRecord *mc)
{
   std::vector<NtpMCStdHep *> v;
   for(int i=0;i<mc->stdhep->GetEntries();i++){
      NtpMCStdHep *stdhep = dynamic_cast<NtpMCStdHep *>((*mc->stdhep)[i]);
      if(stdhep->mc==index){
	 v.push_back(stdhep);
      }
   }
   return v;
}

std::vector<NtpMCStdHep *> SntpHelpers::GetStdHepArray(int index, NtpStRecord *st)
{
   std::vector<NtpMCStdHep *> v;
   for(int i=0;i<st->stdhep->GetEntries();i++){
      NtpMCStdHep *stdhep = dynamic_cast<NtpMCStdHep *>((*st->stdhep)[i]);
      if(stdhep->mc==index){
	 v.push_back(stdhep);
      }
   }
   return v;
}

Int_t SntpHelpers::InPartialRegion(UShort_t plane, UShort_t strip){
  // figure out if this (plane,strip) corresponds to something in
  // the partially instrumented region
  //
  // this region is defined as:
  // v planes: (strip<=4 || strip>=67)
  // partial u: (strip==0 || strip=63)
  // full u: (strip<=26 || strip>=88) 
  //
  // if so, return 1
  // if not, return -1
  // if error, return 0


  // make a lookup ptype to hold the type of each plane
  // 1 = v partial   2 = u partial
  // 3 = v full   4 = u full
  // 0 = uninstrumented
  static bool first=true;
  static UShort_t ptype[282];
  if(first){
    ptype[0]=0;
    for(int i=1; i<=281; i++){
      if(i%2==0) ptype[i]=1; // a v plane
      else ptype[i]=2; // a u plane
      if((i-1)%5 == 0) ptype[i]+=2; // fully instrumented
      else if(i>120) ptype[i]=0; // not instrumented
    }
    first=false;
  }
  if(plane>281){
    //    std::cerr<<"InPartialRegion passed plane = "<<plane<<std::endl;
    return 0;
  }
  UShort_t pt = ptype[plane];
  
  Int_t result;
  switch(pt){
  case 1:
  case 3:
    if(strip<=4 || strip>=67) result=1;
    else result = -1;
    break;
  case 2:
    if(strip==0 || strip == 63) result=1;
    else result = -1;
    break;
  case 4:
    if(strip<=26 || strip>=88) result=1;
    else result = -1;
    break;
  case 0:
  default:
    result=0;
    break;
  }
  return result;

}

NtpMREvent *SntpHelpers::GetMREvent(int ind,NtpMRRecord *mr)
{
   NtpMREvent *mrevt=0;
   if(mr==0){
      return mrevt;
   }
   if(ind>=(int)(mr->mrhdr.nmrevt)){
      return mrevt;
   }
   mrevt = dynamic_cast<NtpMREvent *>((*mr->mrevt)[ind]);
   return mrevt;
}

NtpMRTruth *SntpHelpers::GetMRTruth(int ind,NtpMRRecord *mr)
{
   NtpMRTruth *mrtru=0;
   if(mr==0) return mrtru;
   if(mr->GetHeader().GetVldContext().GetSimFlag()!=4) return mrtru;
   if(ind>=(int)(mr->mrhdr.nmrevt)) return mrtru;
   mrtru = dynamic_cast<NtpMRTruth *>((*mr->mrtru)[ind]);
   return mrtru;
}

void SntpHelpers::FillEventEnergy(float* ph0, float* ph1, int evtn, NtpStRecord* st, const int SIZE)
{
   NtpSREvent *event = 0;
   event = SntpHelpers::GetEvent(evtn,st);
   if(event == 0) return;
                                                                                
   for(int j=0;j<event->ntrack;j++){
     int tindex = SntpHelpers::GetTrackIndex(j,event);
     NtpSRTrack *track = SntpHelpers::GetTrack(tindex,st);
     for(int k = 0; k < track->nstrip; k++){
        int index = SntpHelpers::GetStripIndex(k, track);
        if(index > SIZE){
          MSG("SntpHelper",Msg::kFatal)<<" evt mip array insufficiently large: "
              <<"index "<<index<<"/"<<SIZE<<" requested"<<std::endl;
        }
                                                                                
        float mips1 = track->stpph1mip[k];
        float mips0 = track->stpph0mip[k];
                                                                                
        if(ph0[index] < 0)   ph0[index] = mips0;
        if(ph1[index] < 0)   ph1[index] = mips1;
       }
    }
    for(int j=0;j<event->nshower;j++){
       int sindex = SntpHelpers::GetShowerIndex(j,event);
       NtpSRShower *shw = SntpHelpers::GetShower(sindex,st);
       for(int k = 0; k < shw->nstrip; k++){
          int index = SntpHelpers::GetStripIndex(k, shw);
          if(index > SIZE){
          MSG("SntpHelper",Msg::kFatal)<<" evt mip array insufficiently large: "
              <<"index "<<index<<"/"<<SIZE<<" requested"<<std::endl;
          }
                                                                                
          float mips1 = shw->stpph1mip[k];
          float mips0 = shw->stpph0mip[k];
                                                                                
          if(ph0[index] < 0)   ph0[index] = mips0;
          if(ph1[index] < 0)   ph1[index] = mips1;
       }
    }
    return;
}

