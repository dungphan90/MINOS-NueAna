#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedHit.h"

using namespace Managed;

ClassImp(ManagedHit)

int ManagedHit::idcounter=0;

ManagedHit::ManagedHit(int view,int plane, int strip, double z,double t,double e)
{
	id=idcounter++;

	this->view=view;
	this->z=z;
	this->t=t;
	this->e_original=e;
	this->plane=plane;
	this->strip=strip;
	e_remaining=e_original;
}

ManagedHit::~ManagedHit()
{}

void ManagedHit::AdvanceID()
{
	id=idcounter++;
}

void ManagedHit::ResetIDCounter()
{
	idcounter=0;
}

