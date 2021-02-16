#ifndef ClusterManager_H
#define ClusterManager_H
#include "NueAna/ParticlePID/ParticleFinder/Managed/HitManager.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedHit.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedCluster.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ClusterSaver.h"

#include <vector>
#include <map>

#include "TObject.h"

namespace Managed
{

class ClusterSaver;

class ClusterManager : public TObject
{
	public:
		ClusterManager();
		~ClusterManager();
		

		void MakeClusters(double min_cluster_e=0.4, double min_strip_e=0.4, double max_t_skip=0.025);

		void Reset();
		void AddStrip(int myplane, int mystrip, double myenergy, double st, double sz, int view);
		
		void LoadHits();
		
		std::vector<int> GetPlanes();
		std::vector<int> GetStrips(int plane);
		
		Managed::ManagedCluster * GetCluster(int cid);
		Managed::ManagedCluster * GetSavedCluster(int cid);
		
		void SplitEnergy(Managed::ManagedCluster * cluster,double energy, int matchenergy);		
		int MergeClusters(int cid1, int cid2);
		
		void FillClusterMap(std::map<double, std::map<double, std::pair<double, int> > > * cluster_map	);

		std::vector<int> FindClustersInZ(double z,int view);


		std::map<double, std::map<double, int>  > cluster_map;	//loc_z, sum_e_t_, cluster_id
		std::map<double, std::map<double, int>  > cluster_map_u;	//loc_z, sum_e_t_, cluster_id
		std::map<double, std::map<double, int>  > cluster_map_v;	//loc_z, sum_e_t_, cluster_id

		std::map<double, std::map<double, int>  > * GetClusterMap(int view=0);

		std::vector<Managed::ManagedCluster> clusters;
		
		//get a vector of indicies in the clusters vector corresponding to a given view
		std::vector<int>  GetViewIndex(int view);
		
		double maxz;
		double minz;
		double maxu;
		double minu;
		double maxv;
		double minv;
				
		void SetHitManager(Managed::HitManager *h){hitmanager=h;};
		HitManager * GetHitManager(){return hitmanager;};
		
		void SetClusterSaver(Managed::ClusterSaver *cs){clustersaver=cs;};
		Managed::ClusterSaver * GetClusterSaver(){return clustersaver;};
		
		void DumpClusters();
		
		int SaveCluster(int cluster_id, double energy_to_save=0, int status=0);
		
		void ClearInUse(){inuse.clear();};
		
		void MarkUsed(int id){inuse.push_back(id);};
		
	private:
		int clusters_are_made;
		std::map<int, std::map<int,int> > loc_map;
		std::vector<ManagedHit> hits;
		HitManager *hitmanager;
		ClusterSaver *clustersaver;
		
		void EraseCluster(int cid);
		
		void AdjustCluster(Managed::ManagedCluster *cluster);
		
		ClassDef(ClusterManager,1);
		
		int needMapRebuild;
		void RebuildClusterMaps();
		
		std::vector<int>inuse;
		
		std::vector<int>clusters_to_delete;
};

}

#endif

