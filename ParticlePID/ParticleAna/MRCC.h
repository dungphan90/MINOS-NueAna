#ifndef MRCC_H
#define MRCC_H

#include "TObject.h"

class MRCC : public TObject
{

	public:
	
		MRCC();
		virtual ~MRCC();
		virtual void Clear(Option_t* option = "");


		double start_u;
		double end_u;
		double start_v;
		double end_v;
		double start_z;
		double end_z;

		double sum_e;

		double muonfrac;

		
		//description of angle of particle from vertex
		//should represent something like fit of first 5 hits, etc
		double theta; //spherical coords, angle in xy plane
		double phi; //spherical coords, angle off of z
				
		double calibrated_energy;
		double avg_rms_t;
		double particle_s;
		
		
		
		int hasMRCC;
		int entries;		
		int stage;

	private:
		void Init();

		
		ClassDef(MRCC,2)

};

#endif

