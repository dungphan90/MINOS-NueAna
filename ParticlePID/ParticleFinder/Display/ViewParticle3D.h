#ifndef VIEWPARTICLE3D_H
#define VIEWPARTICLE3D_H


#include "NueAna/ParticlePID/ParticleFinder/Display/Page.h"

#include "TCanvas.h"
#include "StandardNtuple/NtpStRecord.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleAna/PRecord.h"
#include "TH3F.h"
#include <string>
#include "TPad.h"

class ViewParticle3D : public Page
{
	public:
		ViewParticle3D();
		~ViewParticle3D();
		
		virtual void BuildDisplay(TCanvas *c);
		virtual void DrawEvent(ParticleObjectHolder * poh=0, NtpStRecord *ntp=0, int ntpEvt=-1);
		void DrawParticles(ParticleObjectHolder * poh);
		void DrawLegend();
		void DrawParticlesDetail(ParticleObjectHolder * poh, PRecord * precord);
		void DrawTextLine(std::string text, int ln);
		
	private:
		TH3F * t3;
		TPad *leg;
		TPad *main;
		TPad *detail;
		
		double meupergev;

};

#endif

