#ifndef PARTICLE3D_H
#define PARTICLE3D_H
#include "TObject.h"
#include "ParticleType.h"
#include <vector>





class Particle3D :public TObject
{




	public:
	
	
		typedef enum EParticle3DType
		{
			other=0,
			electron=11,
			muon=13,
			proton=2212,
			photon=22,
			neutron=2112
		}Particle3DType;

	
		Particle3D();
		virtual ~Particle3D();
		
		void ResetHits();// only reset the hits and related info... do not change particle type or other externally calculated values
		
		std::vector<double>u;
		std::vector<double>v;		
		std::vector<double>z;
		std::vector<double>e;
		std::vector<int>chain;		
		std::vector<int>chainhit;	
		std::vector<int>view;
		std::vector<double>rms_t;
		
		std::vector<int>shared;  //mark if its the same chain/chainhit as used in another particle3d
		std::vector<int>lockshared; //mark if its shared, but the energy in this particle3d is fixed (ie: already an extracted muon....)
		std::vector<ParticleType> types;	
				
		double start_u;
		double end_u;
		double start_v;
		double end_v;
		double start_z;
		double end_z;
		double sum_e;
		double muonfrac;

		//description of angle of particle from vertex
		//should represent something like fit of first 5 hits, etc
		double theta; //spherical coords, angle in xy plane
		double phi; //spherical coords, angle off of z			
		
		double emfit_a;
		double emfit_b;
		double emfit_e0;
		double emfit_a_err;
		double emfit_b_err;
		double emfit_e0_err;
		double emfit_prob;
		double emfit_chisq;
		double emfit_ndf;
		
		
		double pred_e_a;
		double pred_g_a;
		double pred_b;
		double pred_e0;
		double pred_e_chisq;
		double pred_e_ndf;
		double pred_g_chisq;
		double pred_g_ndf;
	
		double pre_over;
		double pre_under;
		double post_over;
		double post_under;	
		
				
		double calibrated_energy;
		
		
		double avg_rms_t;

	
		double pp_chisq;
		double pp_p;

		Particle3DType particletype;
		int entries;
		int numshared;	
		int pp_ndf;
		int pp_igood;
		
		double cmp_chisq;
		int cmp_ndf;
		double peakdiff;				
				
		void SetEnergy(int ichain,int ichainhit,int inview, double energy);
		double GetNextEnergy(int chain,int chainhit,int inview);
		double GetPreviousEnergy(int chain,int chainhit,int inview);
		
		
		void SetShared(int ichain, int ichainhit);
		void UnsetShared(int ichain, int ichainhit);
		
		void finalize();
		
		int ShareLocked(int ichain, int ichainhit);
		
		void add_to_back(double iu, double iv, double iz, double ie, int chain, int chainhit, int view, double irms_t);  //allways assuming adding from front to back

		int hasShared(){return numshared;};
		
		void Clean();
		

		
	


		

		
	private:
		
	

		
		double muon_threshold_max;
		double muon_threshold_min;
		int muon_min_hits; 
		int muonlike;	
		
		
			ClassDef(Particle3D,2)
				

};

#endif

