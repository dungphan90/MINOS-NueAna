#ifndef LongMuonFinder_H
#define LongMuonFinder_H

#include "NueAna/ParticlePID/ParticleFinder/Managed/ClusterManager.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle3D.h"
#include "NueAna/ParticlePID/ParticleFinder/Chain.h"
#include <vector>
#include "Conventions/Detector.h"

class LongMuonFinder
{
	public:
	
		LongMuonFinder(Detector::Detector_t d);
		~LongMuonFinder();
		
		//find a long muon particle, returning 0 if none found
		int FindLongMuon(Managed::ClusterManager *cm);

		
		Particle3D* GetParticle3D(){if(foundparticle)return foundparticle3d;return 0;};
		Chain * GetChain(int view){if(foundparticle){if(view==2)return chain_u;if(view==3)return chain_v;}return 0;};
	
	
		int FoundSingleViewLongMuon(){return single_view_long_muon;};
	
	private:
		//we want to check the found chain to make sure that it is muon-like	
		int CheckChainQuality(Chain *c, int view, int partialcount);

		//find a muon chain in a given view
		Chain * FindMuonChain(Managed::ClusterManager *cl, int view);
	
		void Reset();	
				
		//attempt to make a particle3d from the found chains
		void MakeParticle3D();
		
		void AbsorbMuonClusters(Chain *c,int view,double past_z);
		void RemoveNonMuonEnergy(Chain *c,int view,double past_z);

		//take a chain in each view and extrapolate 
		//  how many hits in each view are in the partially instrumented region
		std::pair<int,int> CountInPartiallyInstrumentedRegion(Chain *viewU, Chain *viewV);

		static Particle3D * foundparticle3d;
		int foundparticle;
		static Chain * chain_u;
		static Chain * chain_v;
		
		void DumpParticle3D();
			
		Managed::ClusterManager *cluster_manager;
		
			
		int single_view_long_muon;
		double FindIsolationZ();
		
		//find the appropriate cluster in the chain, and then merge the other clusters into it!
		//adjust the stored information about the cluster in the chain as well!
		void MergeChainClusters(Chain * ch, std::vector<int> *clusters);

		int CheckChainOverlap(Chain * chain_u, Chain * chain_v, double isolation_z);
		void ClearFrontVertex(Chain * chain_u, Chain * chain_v);
		Detector::Detector_t detector;
	
		int IsPartiallyInstrumented(double t, double z, int view);
	
};

#endif

