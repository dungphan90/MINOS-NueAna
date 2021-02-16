#include "NueAna/ParticlePID/ParticleAna/MRCC.h"

ClassImp(MRCC)

MRCC::MRCC()
{
	Init();

}

MRCC::~MRCC()
{}


void MRCC::Clear(Option_t* /* option */) {
  // Purpose: Clear memory allocated to arrays so that record can
  // be reused.  

	Init();

}
 
 
 
void MRCC::Init() {
  // 
  // Purpose: Initialize ntuple TClonesArrays
  //
	hasMRCC=0;

	start_u=0;
	end_u=0;
	start_v=0;
	end_v=0;
	start_z=0;
	end_z=0;
		
	sum_e=0;

	muonfrac=0;
	entries=0;
	
	theta=0; 	
	phi=0; 	
				
	calibrated_energy=0;
	avg_rms_t=0;
	particle_s=0;
	stage=0;

}



