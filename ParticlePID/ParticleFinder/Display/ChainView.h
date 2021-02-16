#ifndef CHAINVIEW_H
#define CHAINVIEW_H

#include "NueAna/ParticlePID/ParticleFinder/Display/Page.h"
#include "TH2.h"
#include "TPad.h"

#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleFinder/ChainHelper.h"

#include "TLatex.h"

class ChainView : public Page
{
	
	public:
		ChainView();
		~ChainView();

		virtual void BuildDisplay(TCanvas *c);
		virtual void DrawEvent(ParticleObjectHolder * poh=0,  NtpStRecord *ntp=0,int ntpEvt=-1);
	
		void DrawChains(int view);
		void DrawTrueVertex(int view);
		void DrawVertex(int view);
		void DrawClusters(int view, TH2 * h);
		void DrawInteractionDiagram(NtpStRecord *str,int index);
		void DrawBasicInfo(NtpStRecord *str,int index);

	private:
		TPad * padU;
		TPad * padV;
		TPad * padMC;
		TPad * padInfo;
	
		TH2D * viewU;
		TH2D * viewV;
		
		ParticleObjectHolder * mypoh;



		
		TLatex * info6;		
		TLatex * info7;
		TLatex * info8;
		TLatex * info9;    
  		TLatex * info10;
  		TLatex * info11;
  		TLatex * info12;
  		TLatex * info13;

};

#endif

