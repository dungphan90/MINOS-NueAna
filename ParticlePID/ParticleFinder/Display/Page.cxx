#include "NueAna/ParticlePID/ParticleFinder/Display/Page.h"


Page::Page():myCanvas(0)
{}


Page::~Page()
{}


void Page::BuildDisplay(TCanvas *c)
{
	myCanvas=c;
}



void Page::DrawEvent(ParticleObjectHolder * /*poh*/, NtpStRecord * /*ntp*/,int /*ntpEvt*/)
{
	
	
}

