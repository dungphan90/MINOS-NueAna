#include "NueAna/ParticlePID/ParticleFinder/Display/ViewParticle3D.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle3D.h"

#include "TClonesArray.h"
#include "TMarker3DBox.h"
#include "TPolyLine3D.h"
#include "TLine.h"
#include "TLatex.h"
#include <map>


#include "NueAna/ParticlePID/ParticleAna/ParticlesAna.h"


ViewParticle3D::ViewParticle3D() :Page()
{
	myCanvas=0;

	meupergev = 25.;//rough value... need not be precise here
	
}


ViewParticle3D::~ViewParticle3D()
{}

		
void ViewParticle3D::BuildDisplay(TCanvas *c)
{
	myCanvas=c;
	t3 = new TH3F("th3","3D Particles",100,0,100,100,-10,10,100,-10,10);
	c->cd();
	
	main=new TPad("main","main",0,0.05,0.7,1);
	leg = new TPad("leg","leg",0,0,0.7,0.05);
	detail = new TPad("detail","detail",0.7,0,1,1);
	main->Draw();
	leg->Draw();
	detail->Draw();
	
	main->cd();
	t3->Draw();
}


void ViewParticle3D::DrawEvent(ParticleObjectHolder * poh, NtpStRecord * /*ntp*/,int /*ntpEvt*/)
{

	if(!poh)return;
 
 	PRecord * precord = new PRecord();
  	ParticlesAna pa;
  	pa.ana(poh,&precord->particles);



	if(!myCanvas)return;

	main->cd();
	t3->Draw();
	DrawParticles(poh);

	leg->cd();
	DrawLegend();
	
	detail->Clear();
	detail->cd();
	DrawParticlesDetail(poh,precord);
	myCanvas->Update();
}

void ViewParticle3D::DrawParticlesDetail(ParticleObjectHolder * poh, PRecord * precord)
{
	if(!poh)return;

	//TClonesArray * particles3d = poh->particles3d1;
	std::vector<Particle3D>  particles3d = poh->particles3d;
		
		
	int ln=1;
	std::string text="List of particles:";
	DrawTextLine(text,ln++);


	std::map<double,int> emap;
	for(unsigned int i=0;i<particles3d.size();i++)
	{
		Particle3D * p3 = (Particle3D *)&particles3d[i];
		if(p3 <=0)continue;
		emap.insert(std::pair<double,int>(p3->sum_e,i));
	}

	char t[200];	
	std::map<double,int>::reverse_iterator it;
	for(it=emap.rbegin();it!=emap.rend();it++)
	{
		Particle3D * p3 = (Particle3D *)&particles3d[it->second];

		if(p3->particletype==Particle3D::other)sprintf(t,"Other");
		if(p3->particletype==Particle3D::muon)sprintf(t,"Muon");
		if(p3->particletype==Particle3D::electron)sprintf(t,"Electron");		
		if(p3->particletype==Particle3D::proton)sprintf(t,"Proton");	
		if(p3->particletype==Particle3D::neutron)sprintf(t,"Neutron");	

		sprintf(t,"%s %2.1f GeV ",t,(p3->sum_e/meupergev));
	
		sprintf(t,"%s (%2.1f,%2.1f,%2.1f)",t,p3->start_u,p3->start_v,p3->start_z);
		sprintf(t,"%s (%2.1f,%2.1f,%2.1f)",t,p3->end_u,p3->end_v,p3->end_z);
	

		DrawTextLine(std::string(t),ln++);
	}
	
	DrawTextLine("",ln++);
	DrawTextLine("",ln++);
	DrawTextLine("Details",ln++);
	
	if(precord)
	{
		sprintf(t,"Total Visible Energy: %2.1f GeV",precord->particles.totvise/meupergev);
		DrawTextLine(std::string(t),ln++);	
	}
		
}

void ViewParticle3D::DrawTextLine(std::string text, int ln) 
{
    TLatex *tex;
    tex = new TLatex(0,1-ln*0.03,text.c_str());
    tex->SetTextFont(132);
    tex->SetTextSize(0.05);
    tex->Draw();
}	


void ViewParticle3D::DrawParticles(ParticleObjectHolder * poh)
{
	if(!poh)return;

//	TClonesArray * particles3d = poh->particles3d1;
	std::vector<Particle3D>  particles3d = poh->particles3d;

	double u0=10000;
	double v0=10000;
	double z0=10000;	
	double u1=-10000;
	double v1=-10000;
	double z1=-10000;


	for(unsigned int i=0;i<particles3d.size();i++)
	{
		Particle3D * p3 = (Particle3D *)&particles3d[i];
		if(p3 <=0)continue;
		
		if(p3->entries <1)continue;
		
		TPolyLine3D * line = new TPolyLine3D(p3->entries);
		int color=2;
		if(p3->particletype)
		{
			if(p3->particletype==Particle3D::muon)
			{
				color=3;
			}
			if(p3->particletype==Particle3D::electron)
			{
				color=4;
			}
			if(p3->particletype==Particle3D::proton)
			{
				color=5;
			}
			if(p3->particletype==Particle3D::neutron)
			{
				color=6;
			}
		}
		
		line->SetLineColor(color);	
		
		for(int j=0;j<p3->entries;j++)
		{
				line->SetPoint(j,p3->z[j],p3->u[j],p3->v[j]);


				double s = p3->e[j]/100*0.1;
				TMarker3DBox *m = new TMarker3DBox(p3->z[j],p3->u[j],p3->v[j],s,s,s,0,0);
				

				
				m->SetFillColor(color);
				m->SetLineColor(color);
				
				if(p3->shared[j]>0)m->SetLineColor(5);
				
				m->Draw("same");
			
				u0 = u0 < p3->u[j] ? u0 :  p3->u[j];
				v0 = v0 < p3->v[j] ? v0 :  p3->v[j];
				z0 = z0 < p3->z[j] ? z0 :  p3->z[j];		

			
				u1 = u1 > p3->u[j] ? u1 :  p3->u[j];
				v1 = v1 > p3->v[j] ? v1 :  p3->v[j];
				z1 = z1 > p3->z[j] ? z1 :  p3->z[j];
		
		}
	
		line->Draw();
	}


	double du = (u1-u0)*0.1;
	u1+=du;
	u0-=du;

	double dv = (v1-v0)*0.1;
	v1+=dv;
	v0-=dv;

	double dz = (z1-z0)*0.1;
	z1+=dz;
	z0-=dz;



	TMarker3DBox *m = new TMarker3DBox(poh->event.vtx_z, poh->event.vtx_u, poh->event.vtx_v,0.005,0.005,0.005,0,0);
	m->SetFillColor(1);
	m->SetLineColor(1);
	m->Draw("same");




//	view3d->SetRange(u0,v0,z0,u1,v1,z1);
//	view3d->ShowAxis();
//	view3d->ToggleRulers();
	
	t3->GetYaxis()->SetRangeUser(u0,u1);
	t3->GetZaxis()->SetRangeUser(v0,v1);
	t3->GetXaxis()->SetRangeUser(z0,z1);

	t3->SetStats(0);





}


void ViewParticle3D::DrawLegend()
{

leg->cd();
  //copied Chris' code
    TLine *line;
    line = new TLine(0.01,0.5,0.04,0.5);
    line->SetLineColor(4);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.11,0.5,0.14,0.5);
    line->SetLineColor(3);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.21,0.5,0.24,0.5);
    line->SetLineColor(5);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.31,0.5,0.34,0.5);
    line->SetLineColor(6);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.41,0.5,0.44,0.5);
    line->SetLineColor(2);
    line->SetLineWidth(2);
    line->Draw();
  


    TLatex *tex;
    tex = new TLatex(0.06,0.3,"e");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();
    tex = new TLatex(0.16,0.3,"#mu");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();
    tex = new TLatex(0.26,0.3,"p");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();
    tex = new TLatex(0.36,0.3,"n");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();
    tex = new TLatex(0.46,0.3,"other");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();


}

