#include "NueAna/ParticlePID/ParticleFinder/EventQuality.h"

ClassImp(EventQuality)

EventQuality::EventQuality()
{
	Reset();
}

EventQuality::~EventQuality()
{}


void EventQuality::Reset()
{

	single_view_long_muon=0;
	single_view_primary_shower=0;
	
	foundlongmuon=0;
	foundprimaryshower=0;
	
}


bool EventQuality::IsQuality()
{
	bool isGood=true;

	isGood = isGood && single_view_long_muon==0;
	isGood = isGood && single_view_primary_shower==0;

	return isGood;
}

