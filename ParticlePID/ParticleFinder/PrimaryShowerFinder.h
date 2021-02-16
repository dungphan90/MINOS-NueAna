#ifndef PrimaryShowerFinder_H
#define PrimaryShowerFinder_H

#include "NueAna/ParticlePID/ParticleFinder/Managed/ClusterManager.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle3D.h"
#include "NueAna/ParticlePID/ParticleFinder/Chain.h"
#include <vector>

#include "TH2D.h"

#include "TObject.h"
#include "NueAna/ParticlePID/ParticleFinder/Chain.h"
#include "NueAna/ParticlePID/ParticleFinder/HoughLine.h"

#include "NueAna/ParticlePID/ParticleFinder/ChainHelper.h"

class PrimaryShowerFinder : public TObject
{
	public:
	
		PrimaryShowerFinder();
		~PrimaryShowerFinder();
		
		//find a long Shower particle, returning 0 if none found
		int FindPrimaryShower(Managed::ClusterManager *cm);

		
		Particle3D* GetParticle3D(){if(foundparticle)return foundparticle3d;return 0;};
		Chain * GetChain(int view){if(foundparticle){if(view==2)return chain_u;if(view==3)return chain_v;}return 0;};
	
	
		int FoundSingleViewPrimaryShower(){return single_view_long_shower;};

		static void LoadCloseHits(HoughLine * hl, Managed::ClusterManager *cm, int view, double max_tdist,int limit_to_current_size=0);
	
		TH2D * GetHoughMap(int view){if(view==2)return houghmapU; if(view==3)return houghmapV; return 0;};
		std::vector<HoughLine> * GetHoughLines(int view){if(view==2)return &houghlinesU;if(view==3)return &houghlinesV; return 0;};
	
		//find the z position that has the first intersection of hough lines
		//returns 1 if intersection found
		int FindFirstIntersection(int view, double &t, double &z);


		//merge similar houghlines which roughly have the same hits
		void Bundle(int view);
	
		static TH2D * intU;
		static TH2D * intV;	

		static ChainHelper *chu;
		static ChainHelper *chv;
		
		int ran;
		
		
		void SetVertex(double u, double v, double z){vtx_u=u;vtx_v=v;vtx_z=z;foundvertex=1; foundvertexU=1;foundvertexV=1;};
		
	
	private:
		//we want to check the found chain to make sure that it is Shower-like	
		int CheckChainQuality(Chain *c);

		//find a Shower chain in a given view
		Chain * FindShowerChain(Managed::ClusterManager *cl, int view);
	
		void Reset(int reset_vertex=1);	
				
		//attempt to make a particle3d from the found chains
		void MakeParticle3D();
		
		void MakeHoughMap(int view);
		void SaveHitGroup(TH2D * his, TH2D * saver, double save_val, double with_val, int curx, int cury);
		void GetPeakAreaAverage(double &x, double &y,double &val, int & cnt, int curx, int cury, TH2D * hist, std::vector<double> &peakvalue, std::vector<int> &peakstillgood, std::vector<int> &peakbinx,std::vector<int>  &peakbiny);
	
	//	void GetPeakAreaAverage(double &xv, double &yv, double &val, std::vector<int> &peakbinx, std::vector<int> &peakbiny, std::vector<double> &peakx, std::vector<double> &peaky, std::vector<double> &peakvalue, std::vector<int> &peakstillgood);


		void ExpandShowerChain(Chain * chain,int view, double start_z, double end_z);
		Chain * ExtractMuon(Chain *c);



		static Particle3D * foundparticle3d;
		int foundparticle;
		static Chain * chain_u;
		static Chain * chain_v;
		
		void DumpParticle3D();

		void DumpHoughLines(int view);

		void MakeChains(int view);
			
		Managed::ClusterManager *cluster_manager;
		Chain *  MakeShowerChain(int view);
			
		int single_view_long_shower;
		
		//find the appropriate cluster in the chain, and then merge the other clusters into it!
		//adjust the stored information about the cluster in the chain as well!
		void MergeChainClusters(Chain * ch, std::vector<int> *clusters);

		int CheckChainOverlap(Chain * chain_u, Chain * chain_v);
		void ClearFrontVertex(Chain * chain_u, Chain * chain_v);

	
		static TH2D * houghmapU;
		static TH2D * houghmapV;
		
		
		HoughLine *SplitHoughLine(Chain *c,HoughLine *hl, Managed::ClusterManager *mycm=0);
		
		std::vector<HoughLine> houghlinesU;
		std::vector<HoughLine> houghlinesV;
	
		std::vector<std::pair<int,int> >houghlineMatch;
		
		void FindHoughLineMatches();
	
		void CleanHoughLines(int view,std::vector<HoughLine>*hhh);
			
		//vertex to use when finding chains
		//either from the long muon finder or from the first chain found	
		double vtx_u;
		double vtx_v;
		double vtx_z;
		int foundvertex;
		int foundvertexU;
		int foundvertexV;


		ClassDef(PrimaryShowerFinder,1);
		
		



};


#endif

