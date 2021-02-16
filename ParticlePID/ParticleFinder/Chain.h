#ifndef CHAIN_H
#define CHAIN_H

#include "TObject.h"
#include <vector>

class Chain : public TObject
{
	public:
	
		Chain();
		~Chain();
		
		std::vector<double>t;
		std::vector<double>z;
		std::vector<double>e;
		std::vector<int>cluster_id;
		
		std::vector<int> children;
		
		double start_t;
		double end_t;
		double start_z;
		double end_z;
	
		static int lastid;
		
		int particletype;
	
		double sum_e;
		
		int available;
		
		double weighted_t;
		
		int entries;

		int good_slope();
		
		void Recalc();
		
		static void ResetCounter();
		void updateMuonFrac(double it, double iz, double ie);
		void updateNumPeaks(double it, double iz, double ie);
		void updateEMLike(double it, double iz, double ie);

		double avg_slope;
		double avg_offset;
		double last_slope; //last 2
		
		double front_slope; //first 4 slope
		double back_slope; //last 4 slope
		double front_offset;
		double back_offset;
		
		void add_to_back(double it, double iz, double ie, int my_cluster_id);  //allways assuming adding from front to back
		void insert(double it, double iz, double ie, int my_cluster_id); //add to the proper location


		int parentChain;
		int myId;
		int level;
		
		double muonfrac;
		double interior_muonfrac;
		double emlike;
		int num_peaks;
		int strict_decreasing;
		int strict_increasing;
		
		
		double muon_threshold_max;
		double muon_threshold_min;
		int   muon_min_hits; 
		
		void Reverse();  //cause the chain to reverse itself
		void ClearHits();		
		
		void PrintChain();


		// t=a+bz
		double a;
		double b;		

		double interpolate(double z);//returns best guess at t
		static double interpolate(double z0, double t0, double z1, double t1, double z2, double t2,  double z);//interpolate after fitting a parobola

	private:
		
		double sum_z;
		double sum_t;
		double sum_z_z;
		double sum_z_t;
		
		int muonlike;
		
		double lastpeake;
		int lastpeakprob;

		ClassDef(Chain,1)

		

};

#endif

