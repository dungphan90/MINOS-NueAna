#ifndef PRECORD_H
#define PRECORD_H

//#include "TNamed.h"

#ifndef RECRECORDIMP_H
#include "Record/RecRecordImp.h" // base class
#endif
#ifndef RECCANDHEADER_H
#include "Record/RecCandHeader.h" 
#endif

#include "NueAna/ParticlePID/ParticleAna/Event.h"
#include "NueAna/ParticlePID/ParticleFinder/MCTrue.h"
#include "NueAna/ParticlePID/ParticleAna/Particles.h"
#include "NueAna/ParticlePID/ParticleAna/TruthCompare.h"
#include "NueAna/ParticlePID/ParticleAna/MRCC.h"

class PRecord : public RecRecordImp<RecCandHeader>
{
	public:

		PRecord();
		PRecord(const RecCandHeader& header);
		virtual ~PRecord();
			
		void Reset();

		Particles particles;  
		Event event; 
		MCTrue mctrue; //!	
		TruthCompare truthcompare;  //!

		MRCC  mrccinfo;

	private:



	ClassDef(PRecord,3)

};

#endif

