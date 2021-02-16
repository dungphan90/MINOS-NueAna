#ifndef FINDER_H
#define FINDER_H


#include <vector>
#include "NueAna/ParticlePID/ParticleFinder/ParticleObject.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle3D.h"
#include "NueAna/ParticlePID/ParticleFinder/ChainHelper.h"

//#include "NueAna/ParticlePID/ParticleFinder/ClusterHelper.h"
#include <map>
#include <string>

#include "TH2F.h"
#include "TCanvas.h"

#include "ShareHolder.h"

#include "NueAna/ParticlePID/ParticleFinder/Managed/ClusterManager.h"

class Finder
{
	public:

		Finder();	
		~Finder();

		void Reset();
		
		void AddStrip(int myplane, int mystrip, double myenergy, double st, double sz, int view);

		int GetStrips(){return plane.size();};

		void RecordLostHits(std::vector<Particle3D*> p3d);


		void SetPOH(ParticleObjectHolder *p){mypoh = p;};

		void Process(ParticleObjectHolder &p);

		void FindIsolatedHits();		


		void FindMuons(ChainHelper *ch);
		
		void FindNeutrons(std::vector<Particle3D*> & pout);

		void MergeChains(ChainHelper *ch);

		void SetTrueVertex(double myu, double myv, double myz){true_vtx_u = myu;true_vtx_v = myv;true_vtx_z = myz;};

		void FindVertex(ChainHelper * cu, ChainHelper * cv);

		void MakeChains(Managed::ClusterManager *cl, ChainHelper * ch, int view);
		
		void RemoveNonVertexPointingChains(ChainHelper *ch, int view);

		void Weave(ChainHelper * chu, ChainHelper * chv);
	
		void ClearXTalk();


		std::vector<Particle3D*> ProcessParticle3D(Particle3D * p3d);

		std::vector<Particle3D*> ProcessParticle3D1(Particle3D * p3d, int check_unique_muon=1);
		std::pair<Particle3D*,Particle3D*> StripMuon1(Particle3D * p3d);
		void RemoveSharedHits(std::vector<Particle3D*> pv);



		std::vector<Particle3D*> SetShared(std::vector<Particle3D*> p3v);

		void ShareHit(int view, int chain, int chainhit, std::vector<Particle3D*> shared);

		std::vector<std::pair<int,int> > matchViews(std::vector<foundpath> pu, std::vector<foundpath> pv);
		std::vector<Particle3D> shareEnergy(std::vector<Particle3D> p3v);
	
		std::pair<Particle3D*,Particle3D*> StripMuon(Particle3D * p3d);
		
		void FinalizeParticles3D(std::vector<Particle3D*>pout);
		
		double vtx_u;
		double vtx_v;
		double vtx_z;
		
		Managed::ClusterManager clustermanager_u;
		Managed::ClusterManager clustermanager_v;
		
		void SetClusterManagerU(Managed::ClusterManager &m){clustermanager_u=m;};
		void SetClusterManagerV(Managed::ClusterManager &m){clustermanager_v=m;};
	
		void SetMRCC(int i){DoMRCC=i;};

		void SetMEUperGeV(double d){meupergev=d;};
		

	
	private:

		ParticleObjectHolder *  mypoh;


		std::vector<int> view;
		std::vector<int> plane;
		std::vector<int> strip;
		std::vector<double> t;
		std::vector<double> z;

		std::vector<double> energy;



		std::vector<Particle> particles;

		std::multimap<double,int> sorter_map;

		std::map<int, std::map<int, int> > loc_map;  //plane, strip, index in vectors
		
        /*  RWH 2014-02-17 these seem never to be used anywhere
		double maxz;
		double minz;
		double maxu;
		double minu;
		double maxv;
		double minv;
        */

		double true_vtx_u;
		double true_vtx_v;
		double true_vtx_z;
		

		
		std::vector<foundpath> GetPaths(ChainHelper*ch);
		void DumpPaths(std::vector<foundpath> a);
		
		Particle3D Make3DParticle(std::vector<int>upath, std::vector<int>vpath, ChainHelper * chu, ChainHelper * chv,int multu, int multv);
		
		void DumpParticle3D(std::vector<Particle3D*> p3d);
		
		double meupergev;
		
		ShareHolder shareholder;
		int DoMRCC;
		
		void SetStatus(Chain *c, int status);
		
};

#endif

