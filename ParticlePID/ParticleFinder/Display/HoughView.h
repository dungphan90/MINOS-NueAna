#ifndef HOUGHVIEW_H
#define HOUGHVIEW_H

#include "NueAna/ParticlePID/ParticleFinder/Display/Page.h"
#include "TH2.h"
#include "TPad.h"

#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleFinder/ChainHelper.h"
#include "NueAna/ParticlePID/ParticleFinder/PrimaryShowerFinder.h"
#include "TLatex.h"

class HoughView : public Page
{
	
	public:
		HoughView();
		~HoughView();

		virtual void BuildDisplay(TCanvas *c);
		virtual void DrawEvent(ParticleObjectHolder * poh=0,  NtpStRecord *ntp=0,int ntpEvt=-1);
	
		void DrawChains(int view);
		void DrawTrueVertex(int view);
		void DrawVertex(int view);
		void DrawClusters(int view, TH2 * h);


		void DrawHough(int view, TH2 * h);

	private:
		TPad * padU;
		TPad * padV;

	
		TH2D * viewU;
		TH2D * viewV;

		TH2D * houghviewU;
		TH2D * houghviewV;

		ParticleObjectHolder * mypoh;


		void ClearHit(TH2D * his, double below_val, int curx, int cury);
		void SaveHitGroup(TH2D * his, TH2D * saver,  double save_val,double with_val, int curx, int cury);
		void GetPeakAreaAverage(double &x, double &y,double &val, int & cnt, int curx, int cury, TH2D * hist);
};

#endif

