#ifndef PARTICLEOBJECTHOLDER_H
#define PARTICLEOBJECTHOLDER_H


#ifndef RECRECORDIMP_H
#include "Record/RecRecordImp.h" // base class
#endif
#ifndef RECCANDHEADER_H
#include "Record/RecCandHeader.h" 
#endif

//#include "TNamed.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleObject.h"
//#include <vector>
#include "NueAna/ParticlePID/ParticleFinder/Particle.h"
#include "TObjArray.h"

#include "NueAna/ParticlePID/ParticleFinder/ChainHelper.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleEvent.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle3D.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedCluster.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ClusterManager.h"
#include "NueAna/ParticlePID/ParticleFinder/StripHolder.h"
#include "NueAna/ParticlePID/ParticleFinder/MCTrue.h"
#include "NueAna/ParticlePID/ParticleFinder/MRCCInfo.h"
#include "NueAna/ParticlePID/ParticleFinder/EventQuality.h"

#include "NueAna/ParticlePID/ParticleFinder/PrimaryShowerFinder.h"

#include <map>

//#include "ParticleTruthObject.h"

class ParticleObjectHolder : public RecRecordImp<RecCandHeader> //: public TNamed
{
	public:

		ParticleObjectHolder();
		ParticleObjectHolder(const RecCandHeader& header);
		virtual ~ParticleObjectHolder();

        //void AddParticle(int pid, int bs, int es, int bp, int ep, double sume, double fite, double begt, double endt, double begz, double endz);
		void AddParticle(Particle *p);
		
		
		void AddParticle3D(Particle3D  p);
		
//		std::vector<ParticleObject> particles;



  		virtual void Clear(Option_t* option = "");
		TObjArray * particles;
		TObjArray * particletruth;
		TObjArray * particlematch;
	//	TObjArray * particles3d1;


		MCTrue mctrue; 
		
		ParticleEvent  event;
		
		ChainHelper chu;//!
		ChainHelper chv;//!
		
		Managed::ClusterManager clusterer;//!
		
		Managed::ClusterSaver cluster_saver;//!
		
		EventQuality eventquality;
		
		std::vector<Particle3D>particles3d;


	    std::map<double, std::map<double, std::pair<double, int> > > cluster_map;	//!


		//Managed::ManagedCluster cluster;
		
		StripHolder strips; 
		
		PrimaryShowerFinder psf; //!
		
		MRCCInfo *mrcc;
		
		MRCCInfo * GetMRCCInfo();		
	private:
		void Init();


	ClassDef(ParticleObjectHolder,2)

};

#endif

