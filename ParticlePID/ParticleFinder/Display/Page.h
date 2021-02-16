#ifndef PAGE_H
#define PAGE_H

#include "TCanvas.h"
#include "StandardNtuple/NtpStRecord.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"


class Page
{
	public:
		Page();
		virtual ~Page();

		virtual void BuildDisplay(TCanvas *c);
		virtual void DrawEvent(ParticleObjectHolder * poh=0,  NtpStRecord *ntp=0,int ntpEvt=-1);
	
	protected:
	
		TCanvas * myCanvas;
	

};

#endif
