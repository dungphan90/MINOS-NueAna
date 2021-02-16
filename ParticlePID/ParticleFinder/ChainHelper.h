#ifndef CHAINHELPER_H
#define CHAINHELPER_H

#include <vector>
#include "NueAna/ParticlePID/ParticleFinder/Chain.h"
#include <map>
#include "Conventions/Detector.h"

#include "TObject.h"

struct foundpath
{

	double start_t;
	double start_z;
	double end_t;
	double end_z;
	double energy;
	double muonlike;
	int entries;
	std::vector<int> path;

	foundpath():start_t(0.0),start_z(0.0),end_t(0.0),end_z(0.0),energy(0.0),muonlike(0.0),entries(0){};

};



class ChainHelper :public TObject
{
	public:
	
		ChainHelper();
		~ChainHelper();
		
		
		std::vector<Chain > working;
		std::vector<Chain > finished;
		
		std::vector<Chain > particles;
		
		int parents;
		int muonlikechains;
		int emlikechains;
		int totalchains;
		
		void insert(double it, double iz, double ie, int my_cluster_id);
		
		Chain * split(Chain  * d);
		
		Chain * SplitAt(Chain * c, double splitpointz);
		
		Chain  AttachAt(Chain *mother, Chain * daughter, double z); //returns other daughter... put that in the right place so its not lost!
		
		void AdjustLevel(Chain * c, int level);
		void AttachAsChild(Chain * parent, Chain * child);

		void Reset();  //only resets vectors, does not reset helper parameters....
		
		
		void print_finished();
		void finish();
		
		void SetDetector(Detector::Detector_t d){detector=d;};

		void matchChains();  //connect chains that point to eachother that were not previously matched due to large gaps....
		
		
		std::vector<int> FindMaxPath();
		std::pair< std::vector<int>, double> FindMaxPath(Chain * c);
		
		
		int NewChain();
		
		void ChangeHitId(int oldid, int newid);
		
		std::vector<int> GetAllChildren(Chain *c);
			
		void add_to_plane(double it, double iz, double ie, int cluster_id);
		void process_plane();
		
		std::vector<double> pending_t;
		std::vector<double> pending_z;
		std::vector<double> pending_e;
		std::vector<int> pending_cluster_id;
		
		Chain* GetChain(int id);
		void DeleteChain(int id);
		
		double max_z_gap;
		double max_slope;
		
		int found_max_path;

		std::vector<int> maxpath;

		void AddFinishedChain(Chain c);
		
		
		std::vector<foundpath> FindPaths(Chain*c);
		
		
		double vtx_z;
		double vtx_t;
		
	private:


		std::map<int, int> idhelper;

		Detector::Detector_t detector;
		ClassDef(ChainHelper,1);

};


struct LTId
{
	bool operator()(Chain a, Chain b) 
	{
		return a.myId < b.myId;
	}
};




#endif

