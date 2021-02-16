#include "NueAna/ParticlePID/ParticleAna/PRecordAna.h"
#include "NueAna/ParticlePID/ParticleAna/EventAna.h"
#include "NueAna/ParticlePID/ParticleAna/MCTrueAna.h"
#include "NueAna/ParticlePID/ParticleAna/ParticlesAna.h"
#include "NueAna/ParticlePID/ParticleAna/TruthCompareAna.h"
#include "NueAna/ParticlePID/ParticleAna/MRCCAna.h"



PRecordAna::PRecordAna() 
{


}

PRecordAna::~PRecordAna()
{




}

void PRecordAna::ana(ParticleObjectHolder * h, PRecord * r,ParticleBeamMon  *  /*bmon*/)
{

	EventAna eventana;
	eventana.ana(h,&(r->event));


//	MCTrueAna mctrueana;
//	mctrueana.ana(h,&(r->mctrue),bmon);

	r->mctrue=h->mctrue;

	ParticlesAna particles;
	particles.ana(h,&(r->particles));
	

	TruthCompareAna a;
	a.ana(h,&r->truthcompare);
	
	MRCCAna ma;

        if((&(h->mrcc)!=0) && (&(h->mrcc->removedmuon)!=0))
        {
		
		ma.ana(h,&r->mrccinfo);
	}



}


