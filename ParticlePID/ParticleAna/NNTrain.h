#ifndef NNTrain_H
#define NNTrain_H
#include "TMultiLayerPerceptron.h"
#include "OscProb/OscCalc.h"
#include "NueAna/NueRecord.h"
#include "TF1.h"
#include <string>
#include <vector>
using std::string;

#include "TROOT.h"


class OscCalc;

class NNTrain
{
	public:
		NNTrain();
		~NNTrain();
		
		void SetInputFile(int type, string file){trainfile[type].clear();trainfile[type].push_back(file);};
		void AddInputFile(int type, string file){trainfile[type].push_back(file);};
		void Train(int steps=20, int update=5, int method=1,string form="5",double ncemcut=0.3,double veemcut=0.2);
		
		void SetTrainFile(string file,int needToWeight=0){inputFile=file;this->needToWeight=needToWeight;};
		
		void SetTrainContinue(int v){traincontinue=v;};

		void MakeTestTree();
		
		void SetDelta(double v){delta=v;};
		void SetEpsilon(double v){epsilon=v;};
		void SetEta(double v){eta=v;};
		void SetEtaDecay(double v){etadecay=v;};
		void SetTau(double v){tau=v;};
		void ResetTrainParams();
		TMultiLayerPerceptron * mlp;

		TTree * MakeTrainTree(int makeTestTree=0);

		void FillTreePid(string file,string outfile="nnout.root", string MLPfile="MLP.root"); 



	private:
		int traincontinue;
		std::vector<string> trainfile[3];
		string inputFile;
		double delta;
		double epsilon;
		double eta;
		double etadecay;
		double tau;

		void SetBranches(TTree*t);

		double ncemcut;
		double veemcut;

//vars

	
	Double_t particles_prim_cmp_chisq;
	Int_t particles_prim_cmp_ndf;
	Double_t particles_prim_peakdiff;	

	Double_t particles_largest_particle_cmp_chisq;
	Int_t particles_largest_particle_cmp_ndf;
	Double_t particles_largest_particle_peakdiff;	

	Double_t particles_longest_s_particle_cmp_chisq;
	Int_t particles_longest_s_particle_cmp_ndf;
	Double_t particles_longest_s_particle_peakdiff;	
				
	Double_t         particles_primary_long_e;		
	Double_t         particles_longest_s_particle_s;
	Double_t         particles_elec_vise;
	Double_t         particles_primary_phi;
	Double_t         particles_mol_rad_r;
	Double_t         particles_longest_z;
	Double_t         particles_emfrac;
	Int_t            particles_ntot;
	Double_t         particles_total_long_e;
	Double_t         particles_weighted_phi;		
	Float_t			 mctrue_nuenergy;
	
	Int_t			 mctrue_inu;
	Int_t			 mctrue_inunoosc;
	
	Double_t		 particles_frac_particle_2;

	Double_t		 particles_largest_particle_e;
	Double_t	     particles_totvise;
	
	Double_t		 particles_longest_s_particle_e;
	Double_t		 particles_rms_r;
		
	Double_t		 particles_prim_vise;
	Double_t		 event_visenergy;
	Double_t		 event_max_z;
	Double_t		 event_min_z;
	Int_t		 	 particles_longest_particle_type;
	
	Int_t			 event_nstrips;
	Int_t			 event_inFiducial;

	Int_t			 mctrue_type;
	Float_t		 mctrue_oscprob;
	Double_t	  	 mctrue_totbeamweight;
	Int_t			 mctrue_iresonance;
	Double_t 		 mctrue_visible_energy;
	
	Double_t		 trainweight;
	Int_t 			 type;
	
	Double_t		 largest_frac;
	Double_t		 reco_frac;
	Double_t                 ntot_lsps;

	Double_t			 isnc;
	Double_t			 isnue;
	Double_t			 iscc;
	Double_t			 istau;
	Double_t 			 isbeamve;
	Double_t			 isndis;
	Double_t 			 isdis;
	Double_t			 particles_nelec;

	Double_t 			length_z;

	Double_t particles_longest_s_particle_par_b;
	Double_t particles_longest_s_particle_par_e0;
	Double_t particles_longest_s_particle_par_a;
	Double_t particles_longest_s_particle_par_chisq;
	Double_t longest_ae0;

	Double_t particles_prim_par_b;
	Double_t particles_prim_par_e0;
	Double_t particles_prim_par_a;
	Double_t particles_prim_par_chisq;
	Double_t prim_ae0;
	
	Double_t truthcompare_emenergy;
	
	Double_t largest_cmp_chisqndf;
	Float_t	 trueEMFrac;

	Double_t emfrac;

	Double_t weight;

  	Float_t nuerec_shwfit_par_a;
  	Float_t nuerec_shwfit_par_b;
	Float_t nuerec_shwfit_uv_molrad_peak_9s_2pe_dw;
  	Float_t nuerec_shwfit_uv_rms_9s_2pe_dw;
  	Float_t nuerec_mstvars_e4w;
  	Float_t nuerec_mstvars_o4w;
	Float_t mstvar_combine;

    Float_t nuerec_fracvars_fract_2_planes;
    Float_t nuerec_fracvars_fract_4_planes;
    Float_t nuerec_fracvars_fract_6_planes;
    Float_t nuerec_fracvars_fract_8_counters;
    Float_t nuerec_fracvars_fract_road;
    Float_t nuerec_shwfit_LongE;
	Float_t mctrue_emShowerFraction;
	
	Double_t tweight;

		OscCalc fOscCalc;
		float fBaseLine;
		double fDeltaMS12;
		double fTh12;
		double fTh23;
		double fDensity;
		void SetOscParamBase( float dm2, float ss13, float delta, int hierarchy);

	Double_t pars[30];
	NueRecord *nr;
	int needToWeight;
	
	TF1 * f2;
	float OscillationProb(TF1* f2, int ntype, float NuE, float sinth23=0.5, float sin2th13=0.15);

	double osc(double nuEnergy, int interactionType, int nonOscFlavor, int oscFlavor);
};


#endif



