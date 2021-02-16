#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "Reporter/ParticleReportObject.h"

#include "Record/RecArrayAllocator.h"

ClassImp(ParticleObjectHolder)

ParticleObjectHolder::ParticleObjectHolder(): 
  RecRecordImp<RecCandHeader>()
{

this->Init();


}


ParticleObjectHolder::ParticleObjectHolder(const RecCandHeader& hdr) : 
  RecRecordImp<RecCandHeader>(hdr),  particles(0)
{
this->Init();

}

MRCCInfo * ParticleObjectHolder::GetMRCCInfo()
{
	if(!mrcc)mrcc=new MRCCInfo();
	return mrcc;
}


void ParticleObjectHolder::Init()
{
	SetClearable(true);

	particles = new TObjArray(0);
	particles->SetOwner(1);
	particletruth = new TObjArray( 0);
	particletruth->SetOwner(1);
	particlematch = new TObjArray(0);
	particlematch->SetOwner(1);
//	particles3d1 = new TObjArray(0);
//	particles3d1->SetOwner(1);
	mrcc=0;
}



ParticleObjectHolder::~ParticleObjectHolder()
{
	particles->Delete();
	delete particles;
	
	particletruth->Delete();
	delete particletruth;
	
	particlematch->Delete();
	delete particlematch;
	
//	particles3d1->Delete();
//	delete particles3d1;
	
			
	particles3d.clear();
	cluster_map.clear();	
	if(mrcc)delete mrcc;mrcc=0;
}

void ParticleObjectHolder::Clear(Option_t* /* option */) {

if(particles){particles->Clear("C");}
if(particletruth){particletruth->Clear("C");}
if(particlematch){particlematch->Clear("C");}
//if(particles3d1){particles3d1->Clear("C");}
if(mrcc)delete mrcc;mrcc=0;
}



void ParticleObjectHolder::AddParticle3D(Particle3D  p)
{

//	Particle3D *n = (Particle3D*)particles3d1->New(particles3d1->GetEntries());


//	Particle3D *n = new Particle3D(p);
//	particles3d1->Add(n);
		
	//(*n) = p;
	
	particles3d.push_back(p);

}


void ParticleObjectHolder::AddParticle(Particle * /*p*/){
//	particles->AddLast(&p);


}




/*
void ParticleObjectHolder::AddParticle(int pid, int bs, int es, int bp, int ep, double sume, double fite, double begt, double endt, double begz, double endz)
{

 ParticleObject p;

	p.particle_id=pid;
	p.begstrip=bs;
	p.endstrip=es;
	p.begplane=bp;
	p.endplane=ep;
	p.sumenergy=sume;
	p.fitenergy=fite;

//	particles.push_back(p);
*/

//	particles->AddLast(p);

	//ParticleObject *p = new ParticleObject();
/*
	ParticleObject *p = (ParticleObject*)particles->New(a);

	p->sumenergy=sume;
	p->endt=endt;
	p->begt=begt;
	p->endz=endz;
	p->begz=begz;
	p->begstrip=bs;
	p->endstrip=es;
	p->begplane=bp;
	p->endplane=ep;

printf("particle added e %f at p %d s %d t %f z %f\n",sume, bp, ep, begt, begz);

//	particles->AddLast(p);

//	particles[a]=p;

	a++;

}*/
