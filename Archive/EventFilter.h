#ifndef EVENTFILTER_H
#define EVENTFILTER_H

class NtpSREvent;
class NtpStRecord;
class NtpSRRecord;
class NueRecord;

namespace EventFilter
{

    Bool_t PassesHiTrackCut(NtpSREvent *event, NtpStRecord *st, int cut);
    Bool_t PassesTrackLikeCut(NtpSREvent *event, NtpStRecord *st, int cut);
    Bool_t PassesHiShowerCut(NtpSREvent *event, NtpStRecord *st, float cut);
    Bool_t PassesLoShowerCut(NtpSREvent *event, NtpStRecord *st, float cut);
    Bool_t PassesAllCuts(NtpSREvent *event, NtpStRecord *st, int trkcut,int trklike, int loevtcut, float hiecut, float loecut, float hishwcut, float loshwcut);

    Bool_t PassesHiTrackCut(NtpSREvent *event, NtpSRRecord *st, int cut);
    Bool_t PassesTrackLikeCut(NtpSREvent *event, NtpSRRecord *st, int cut);
    Bool_t PassesHiShowerCut(NtpSREvent *event, NtpSRRecord *st, float cut);
    Bool_t PassesLoShowerCut(NtpSREvent *event, NtpSRRecord *st, float cut);
    Bool_t PassesAllCuts(NtpSREvent *event, NtpSRRecord *st, int trkcut, int trklike, int loevtcut, float hiecut, float loecut, float hishwcut, float loshwcut);


    Bool_t PassesHiEnergyCut(NtpSREvent *event, float cut);
    Bool_t PassesLoEnergyCut(NtpSREvent *event, float cut);
    Bool_t PassesLoEventCut(NtpSREvent *event, int cut);

    Bool_t PassesHiTrackCut(NueRecord *nr, int cut);
    Bool_t PassesTrackLikeCut(NueRecord *nr, int cut);
    Bool_t PassesHiShowerCut(NueRecord *nr, float cut);
    Bool_t PassesLoShowerCut(NueRecord *nr, float cut);
    Bool_t PassesAllCuts(NueRecord *nr, int trkcut,int trklike, int loevtcut, float hiecut, float loecut, float hishwcut, float loshwcut);

    Bool_t PassesHiEnergyCut(NueRecord *nr, float cut);
    Bool_t PassesLoEnergyCut(NueRecord *nr, float cut);
    Bool_t PassesLoEventCut(NueRecord *nr, int cut);
}
#endif// EVENTFILTER_H


