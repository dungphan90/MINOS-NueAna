#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "StandardNtuple/NtpStRecord.h"
#include "NueAna/EventFilter.h"
#include "NueAna/SntpHelpers.h"
#include "NueAna/NueRecord.h"

Bool_t EventFilter::PassesHiTrackCut(NtpSREvent *event, NtpStRecord *st, int cut){


    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}

    int longestTrack=0;
    for(int j=0;j<event->ntrack;j++){
        int tindex = SntpHelpers::GetTrackIndex(j,event);
        NtpSRTrack *track = SntpHelpers::GetTrack(tindex,st);
        if(longestTrack<track->plane.n){
	    longestTrack = track->plane.n;
        }
    }
    
    if(longestTrack<cut){passesCut=true;}
    return passesCut;
}

Bool_t EventFilter::PassesHiTrackCut(NtpSREvent *event, NtpSRRecord *st, int cut){


    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}

    int longestTrack=0;
    for(int j=0;j<event->ntrack;j++){
        int tindex = SntpHelpers::GetTrackIndex(j,event);
        NtpSRTrack *track = SntpHelpers::GetTrack(tindex,st);
        if(longestTrack<track->plane.n){
	    longestTrack = track->plane.n;
        }
    }
    
    if(longestTrack<cut){passesCut=true;}
    return passesCut;
}
Bool_t EventFilter::PassesTrackLikeCut(NtpSREvent *event, NtpStRecord *st, int cut){


    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}

    int longestTrack=0;
    for(int j=0;j<event->ntrack;j++){
        int tindex = SntpHelpers::GetTrackIndex(j,event);
        NtpSRTrack *track = SntpHelpers::GetTrack(tindex,st);
        if(longestTrack<track->plane.ntrklike){
	    longestTrack = track->plane.ntrklike;
        }
    }
    
    if(longestTrack<cut){passesCut=true;}
    return passesCut;
}

Bool_t EventFilter::PassesTrackLikeCut(NtpSREvent *event, NtpSRRecord *st, int cut){


    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}

    int longestTrack=0;
    for(int j=0;j<event->ntrack;j++){
        int tindex = SntpHelpers::GetTrackIndex(j,event);
        NtpSRTrack *track = SntpHelpers::GetTrack(tindex,st);
        if(longestTrack<track->plane.ntrklike){
	    longestTrack = track->plane.ntrklike;
        }
    }
    
    if(longestTrack<cut){passesCut=true;}
    return passesCut;
}
Bool_t EventFilter::PassesHiShowerCut(NtpSREvent *event, NtpStRecord *st, float cut){

    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    

    float highestShwEng=0;
    for(int j=0;j<event->nshower;j++){
        int index = SntpHelpers::GetShowerIndex(j,event);
        NtpSRShower *shower = SntpHelpers::GetShower(index,st);
        if(highestShwEng<shower->ph.gev){
            highestShwEng=shower->ph.gev;
        }
    }
    
    if(highestShwEng<cut){passesCut=true;}
    return passesCut;
}
Bool_t EventFilter::PassesHiShowerCut(NtpSREvent *event, NtpSRRecord *st, float cut){

    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    

    float highestShwEng=0;
    for(int j=0;j<event->nshower;j++){
        int index = SntpHelpers::GetShowerIndex(j,event);
        NtpSRShower *shower = SntpHelpers::GetShower(index,st);
        if(highestShwEng<shower->ph.gev){
            highestShwEng=shower->ph.gev;
        }
    }
    
    if(highestShwEng<cut){passesCut=true;}
    return passesCut;

}

Bool_t EventFilter::PassesLoShowerCut(NtpSREvent *event, NtpStRecord *st, float cut){

    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    
    int nshws = event->nshower;  //no. of showers
    if(nshws==0) {passesCut=false; return passesCut;}

    float highestShwEng=0;
    for(int j=0;j<event->nshower;j++){
        int index = SntpHelpers::GetShowerIndex(j,event);
        NtpSRShower *shower = SntpHelpers::GetShower(index,st);
        if(highestShwEng<shower->ph.gev){
            highestShwEng=shower->ph.gev;
        }
    }
    

    if(nshws==1&&highestShwEng>cut){passesCut=true;}
    else if(nshws>1&&highestShwEng>cut/2){passesCut=true;}
    return passesCut;
}
Bool_t EventFilter::PassesLoShowerCut(NtpSREvent *event, NtpSRRecord *st, float cut){

    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    
    int nshws = event->nshower;  //no. of showers
    if(nshws==0) {passesCut=false; return passesCut;}

    float highestShwEng=0;
    for(int j=0;j<event->nshower;j++){
        int index = SntpHelpers::GetShowerIndex(j,event);
        NtpSRShower *shower = SntpHelpers::GetShower(index,st);
        if(highestShwEng<shower->ph.gev){
            highestShwEng=shower->ph.gev;
        }
    }
    
    if(nshws==1&&highestShwEng>cut){passesCut=true;}
    else if(nshws>1&&highestShwEng>cut/2){passesCut=true;}
    return passesCut;

}




Bool_t EventFilter::PassesHiEnergyCut(NtpSREvent *event, float cut){
    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    if(event->ph.sigcor<cut) passesCut=true;
    return passesCut;
}


Bool_t EventFilter::PassesLoEnergyCut(NtpSREvent *event,  float cut){
    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    if(event->ph.sigcor>cut) passesCut=true;
    return passesCut;
}

Bool_t EventFilter::PassesLoEventCut(NtpSREvent *event,  int cut){
    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    if(event->plane.n>cut) passesCut=true;
    return passesCut;
}


Bool_t EventFilter::PassesAllCuts(NtpSREvent *event, NtpStRecord *st, int trkcut,int trklike, int loevtcut, float hiecut, float loecut, float hishwcut, float loshwcut){


    Bool_t passesCut=false;

    if(PassesHiTrackCut(event,st,trkcut)&&
       PassesTrackLikeCut(event,st,trklike)&&
       PassesHiEnergyCut(event,hiecut) && 
       PassesLoEnergyCut(event,loecut) && 
       PassesHiShowerCut(event,st,hishwcut) && 
       PassesLoShowerCut(event,st,loshwcut) && 
       PassesLoEventCut(event,loevtcut) ) passesCut=true;

    return passesCut;
}



Bool_t EventFilter::PassesAllCuts(NtpSREvent *event, NtpSRRecord *st, int trkcut, int trklike, int loevtcut, float hiecut, float loecut, float hishwcut, float loshwcut){

    Bool_t passesCut=false;

    if(PassesHiTrackCut(event,st,trkcut)&&
       PassesTrackLikeCut(event,st,trklike)&&
       PassesHiEnergyCut(event,hiecut) && 
       PassesLoEnergyCut(event,loecut) && 
       PassesHiShowerCut(event,st,hishwcut) && 
       PassesLoShowerCut(event,st,loshwcut) && 
       PassesLoEventCut(event,loevtcut) ) passesCut=true;

    return passesCut;
}






Bool_t EventFilter::PassesHiTrackCut(NueRecord *nr, int cut){


    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    if(nr->srtrack.planes<cut){passesCut=true;}
    return passesCut;
}

Bool_t EventFilter::PassesTrackLikeCut(NueRecord *nr, int cut){


    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    if(nr->srtrack.trklikePlanes<cut){passesCut=true;}
    return passesCut;
}

Bool_t EventFilter::PassesHiShowerCut(NueRecord *nr, float cut){

    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    if(nr->srshower.phNueGeV<cut){passesCut=true;}
    return passesCut;

}

Bool_t EventFilter::PassesLoShowerCut(NueRecord *nr, float cut){

    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    
    int nshws = nr->srevent.showers;  //no. of showers
    if(nshws==0) {passesCut=false; return passesCut;}

    if(nshws==1 && nr->srshower.phNueGeV > cut){passesCut=true;}
    else if(nshws>1 && nr->srshower.phNueGeV > cut/2.0){passesCut=true;}
    return passesCut;
}

Bool_t EventFilter::PassesHiEnergyCut(NueRecord *nr, float cut){
    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    if(nr->srevent.phMeu<cut) passesCut=true;
    return passesCut;
}


Bool_t EventFilter::PassesLoEnergyCut(NueRecord *nr,  float cut){
    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    if(nr->srevent.phMeu>cut) passesCut=true;
    return passesCut;
}

Bool_t EventFilter::PassesLoEventCut(NueRecord *nr,  int cut){
    Bool_t passesCut=false;
    if(cut<0) {passesCut=true; return passesCut;}
    if(nr->srevent.planes>cut) passesCut=true;
    return passesCut;
}


Bool_t EventFilter::PassesAllCuts(NueRecord *nr, int trkcut,int trklike, int loevtcut, float hiecut, float loecut, float hishwcut, float loshwcut){


    Bool_t passesCut=false;

    if(PassesHiTrackCut(nr,trkcut)&&
       PassesTrackLikeCut(nr,trklike)&&
       PassesHiEnergyCut(nr,hiecut) && 
       PassesLoEnergyCut(nr,loecut) && 
       PassesHiShowerCut(nr,hishwcut) && 
       PassesLoShowerCut(nr,loshwcut) && 
       PassesLoEventCut(nr,loevtcut) ) passesCut=true;

    return passesCut;
}

