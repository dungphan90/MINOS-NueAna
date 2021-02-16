#ifndef ClusterSaver_H
#define ClusterSaver_H
#include "NueAna/ParticlePID/ParticleFinder/Managed/ClusterManager.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/HitManager.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedHit.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedCluster.h"


#include <vector>
#include <map>

#include "TObject.h"

namespace Managed
{

class ClusterSaver :public TObject
{
	public:
		ClusterSaver();
		~ClusterSaver();
		


		void Reset();
			
		
		Managed::ManagedCluster * GetCluster(int cid);
		
		void FillClusterMap(std::map<double, std::map<double, std::pair<double, int> > > * cluster_map	);

		std::map<double, std::map<double, int>  > cluster_map;	//loc_z, sum_e_t_, cluster_id
		std::map<double, std::map<double, int>  > cluster_map_u;	//loc_z, sum_e_t_, cluster_id
		std::map<double, std::map<double, int>  > cluster_map_v;	//loc_z, sum_e_t_, cluster_id

		std::vector<Managed::ManagedCluster> clusters;
		
		int SaveCluster(Managed::ManagedCluster *cluster);

		std::map<double, std::map<double, int>  > * GetClusterMap(int view=0);

		
		double maxz;
		double minz;
		double maxt;
		double mint;
		double minu;
		double maxu;
		double minv;
		double maxv;		
		int nClusters;
		
		void DumpClusters();
		
		int save_id; //holds the next cluster save id;

		//get a map of plane, strip, energy for all strips in use
		//split strips in clusters are merged here...
		std::map<std::pair<int,int>, double> GetStripEnergy();
		void recomputeBounds();
	
	private:

		int needMapRebuild;
		void RebuildClusterMaps();
		std::vector<int>clusters_to_delete;
		
			
	ClassDef(ClusterSaver,1);

};

}

#endif

