#ifndef MAKEPLOTS_H
#define MAKEPLOTS_H

#include "TChain.h"
#include "TH1F.h"
#include "TPad.h"

#include <vector>
#include <string>

#define nCut 4

#define nHistos 22

class MakePlots
{
	public:
		MakePlots();
		~MakePlots(){};
		void FillHistos(int set);
		void SetupHistos();

		//count pots
		std::vector<double> pots;
		
		
		
		
		
		void Run();
		
		void SavePlots(int norm=0);
		void SaveRatioPlots();
		
		void SetRatioPair(int set1, int set2, const char* name){
			ratiopair.push_back(std::pair<int,int>(set1,set2));
			rationame.push_back(name);
		};
		
		void AddSet(const char* filename, const char* typesname, int ismc, int ismrcc=0);		
		
		void SetOutDir(const char*o);
		void SetOutFile(const char*o);
		
		
		void SetPots(int chain, double p){if(chain>-1 && chain<nchains)pots[chain]=p;};
		
	private:
		void SortOutStats(TPad *pad,float fracPadx=0.3,float fracPady = 0.5, 
		  int numPerCol=4,float startx=0.98,float starty=0.98);
		void SetBranches(TChain*c);
	
		void FillHistos(int type,int cut);


		std::vector<std::vector<std::vector<TH1F *> > >histos;

		std::vector<std::pair<int,int> >ratiopair;
		std::vector<std::string>rationame;

		char names[nHistos][200];
		std::vector<int> ismc;
		std::vector<int> ismrcc;
		
		char cutnames[nCut][200];
		
		std::vector<TChain *> chains;
		std::vector<std::string> typesname;
		int nchains;
		
		char outdir[200];
		char outfile[200];
		
		Double_t         particles_primary_long_e;		
	Double_t         particles_longest_s_particle_s;
	Double_t         particles_elec_vise;
	Double_t         particles_primary_phi;
	Double_t         particles_prim_par_b;
	Double_t         particles_mol_rad_r;
	Double_t         particles_longest_z;
	Double_t         particles_emfrac;
	Int_t            particles_ntot;
	Double_t         particles_total_long_e;
	Double_t         particles_weighted_phi;		

	Double_t		 particles_largest_particle_e;
	Double_t	     particles_totvise;
		
	Double_t		 particles_prim_vise;
	Double_t		 event_visenergy;
	Double_t		 event_max_z;
	Double_t		 event_min_z;
	Int_t			 event_nstrips;
	Int_t			 event_inFiducial;
	Int_t			 event_contained;

	Double_t		 event_pidA;
	Double_t		 event_pidB;
	Double_t		 event_pidC;
	Double_t		 event_pidD;

	Int_t			 mctrue_type;
	Double_t		 mctrue_oscprob;
	Double_t	  	 mctrue_totbeamweight;
	Int_t			 mctrue_iresonance;

	Double_t		 trainweight;
	Int_t 			 type;
	
	Double_t		 largest_frac;
	Double_t		 reco_frac;

	Int_t			 isnc;
	Int_t			 isnue;
	Int_t			 iscc;
	Int_t			 istau;
	Int_t 			 isbeamve;
	Int_t			 isndis;
	Int_t 			 isdis;	

	Double_t		 particles_rms_r;
	Double_t		 particles_prim_par_a;
	Double_t		 particles_prim_par_e0;
	Double_t		 particles_prim_par_chisq;
	Double_t		 particles_largest_particle_peakdiff;
	Double_t		 particles_largest_particle_cmp_chisq;
	Int_t		 particles_largest_particlecmp_ndf;

	Double_t		particles_longest_particle_type;
	Double_t		  mrccinfo_particle_s;
	Double_t		mrccinfo_sum_e;
};

#endif
