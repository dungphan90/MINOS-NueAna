#ifndef FITVIEW_H
#define FITVIEW_H

#include "NueAna/ParticlePID/ParticleFinder/Display/Page.h"
#include "TH2.h"
#include "TPad.h"

#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleFinder/ChainHelper.h"


class FitView : public Page
{
	
	public:
		FitView();
		~FitView();

		virtual void BuildDisplay(TCanvas *c);
		virtual void DrawEvent(ParticleObjectHolder * poh=0,  NtpStRecord *ntp=0,int ntpEvt=-1);
	
		void DrawAngles();
		
	private:
		TPad * padU;
		TPad * padV;
		TPad * padMC;
		TPad * padInfo;
	
		TH1D * viewTheta;
		TH1D * viewPhi;
		TH2D * view3D;
		ParticleObjectHolder * mypoh;



};

#endif

