#include "NueAna/ParticlePID/ParticleFinder/Display/FitView.h"
#include "TPolyLine.h"
#include "TMarker.h"
#include "TArrow.h"
#include "TMath.h"

#include "Conventions/SimFlag.h"
#include <vector>
#include "MCNtuple/NtpMCStdHep.h"
#include "MCNtuple/NtpMCRecord.h"
#include "MCNtuple/NtpMCTruth.h"

#include "TruthHelperNtuple/NtpTHEvent.h"

using namespace std;

FitView::FitView():Page()
{
	viewTheta=0;
	viewPhi=0;
	padU=0;
	padV=0;

}


FitView::~FitView()
{}


void FitView::BuildDisplay(TCanvas *c)
{
	myCanvas=c;
	padV=new TPad("chainview_padv","chainview_padu",0,0,0.7,0.5);
	padU=new TPad("chainview_padu","chainview_padv",0,0.5,0.7,1.0);
	padMC=new TPad("chainview_padmc","chainview_padmc",0.71,0.0,1.0,0.8);
	padInfo=new TPad("chainview_padinfo","chainview_padinfo",0.71,0.8,1.0,1.0);
	myCanvas->cd();
	padU->Draw();
	padV->Draw();
	padMC->Draw();
	padInfo->Draw();
	
}


void FitView::DrawEvent(ParticleObjectHolder * poh,  NtpStRecord *ntp,int /*ntpEvt*/)
{
	padU->Clear();
	padV->Clear();
	padMC->Clear();
	padInfo->Clear();
	
	mypoh=poh;
	
	int isMC=0;
	if(poh->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC)isMC=1;
	
	double minz= mypoh->strips.minz;
	double minu= mypoh->strips.minu;
	double minv= mypoh->strips.minv;
	double maxz= mypoh->strips.maxz;
	double maxu= mypoh->strips.maxu;
	double maxv= mypoh->strips.maxv;
	

	
	if(isMC)
	{
		if(mypoh->mctrue.vtx_z>=0)minz=minz<mypoh->mctrue.vtx_z ? minz : mypoh->mctrue.vtx_z;
		maxz=maxz>mypoh->mctrue.vtx_z ? maxz : mypoh->mctrue.vtx_z;
		minu=minu<mypoh->mctrue.vtx_u ? minu : mypoh->mctrue.vtx_u;
		maxu=maxu>mypoh->mctrue.vtx_u ? maxu : mypoh->mctrue.vtx_u;
		minv=minv<mypoh->mctrue.vtx_v ? minv : mypoh->mctrue.vtx_v;
		maxv=maxv>mypoh->mctrue.vtx_v ? maxv : mypoh->mctrue.vtx_v;
	}
	



	minz -= 0.1;
	minu -= 0.1;
	minv -= 0.1;
	maxz += 0.1;
	maxu += 0.1;
	maxv += 0.1;
	
	//int nu=(int)((maxu-minu)/0.0412);
	//int nv=(int)((maxv-minv)/0.0412);
	//int nz=(int)((maxz-minz)/0.0354);


	if(viewTheta)delete viewTheta;
	if(viewPhi)delete viewPhi;

	viewTheta = new TH1D("hu","Theta",100,0,2*3.141592);
	viewPhi = new TH1D("hv","Phi", 100,0,3.141592);
	view3D = new TH2D("view3D","view3D",100,0,2*3.141592,100,0,3.141592);
//	viewTheta->SetStats(0);
//	viewPhi->SetStats(0);
	viewTheta->SetDirectory(0);
	viewPhi->SetDirectory(0);

	if(mypoh)
	{	
		DrawAngles();

	}
	
	if(ntp)
	{

	}
	
	padU->cd();
	viewTheta->Draw();
	padV->cd();
	viewPhi->Draw();
	padMC->cd();
	view3D->Draw("colz");
	
	
	padU->Update();
	padV->Update();
	padMC->Update();
	padInfo->Update();
	myCanvas->Update();
	
}

void FitView::DrawAngles()
{
	//find max particle index...
	double maxe=0;
	int maxi=-1;
	std::vector<Particle3D>  particles3d = mypoh->particles3d;	
	for(unsigned int i=0;i<particles3d.size();i++)
	{
		Particle3D * p3 = (Particle3D *)&particles3d[i];
		if(p3 <=0)continue;
		if(p3->sum_e>maxe)
		{
			maxe=p3->sum_e;
			maxi=i;
		}
	}

	if(maxi<0)return;

	Particle3D * p = &particles3d[maxi];
	
	for(unsigned int i=0;i<p->e.size();i++)
	{
		double u = p->u[i]-mypoh->event.vtx_u;
		double v = p->v[i]-mypoh->event.vtx_v;
		double z = p->z[i]-mypoh->event.vtx_z;
	
		double r = sqrt(u*u+v*v+z*z);
		double theta = atan(u/v);
		double phi = acos(z/r);

		phi=phi<0?phi+2*3.141592:phi;
		theta=theta<0?theta+2*3.141592:theta;
		
		viewPhi->Fill(phi,p->e[i]);
		viewTheta->Fill(theta,p->e[i]);
		view3D->Fill(theta,phi,p->e[i]);
	
	}
	
	



} 
