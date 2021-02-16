#ifndef PARTICLES_H
#define PARTICLES_H

#include "TObject.h"

class Particles : public TObject
{

	public:
	
		Particles();
		virtual ~Particles();

 		virtual void Clear(Option_t* option = "");


		
		double emfrac;
		double totcale;
		double totvise;
		
		double elec_vise;
		double muon_vise;
		double prot_vise;
		double neut_vise;
		double other_vise;


		double  longest_s_particle_s;
		double  longest_s_particle_e;
		double  longest_s_particle_cal_e;
		double  longest_s_particle_type;
		double  longest_s_particle_avg_rms;
		double  longest_s_particle_avg_e;
		
		double  longest_s_particle_z;
		double  longest_s_particle_par_b;
		double  longest_s_particle_par_a;
		double  longest_s_particle_par_e0;
		double  longest_s_particle_par_a_err;
		double  longest_s_particle_par_b_err;
		double  longest_s_particle_par_e0_err;
		double  longest_s_particle_par_prob;
		double  longest_s_particle_par_chisq;
		double  longest_s_particle_par_ndf;		
		double  longest_s_long_e;


		double largest_particle_s;
		double largest_particle_e;
		double largest_particle_cal_e;
		double largest_particle_type;
		double largest_particle_avg_rms;
		double largest_particle_avg_e;
		
		double largest_particle_z;
		double largest_particle_par_b;
		double largest_particle_par_a;
		double largest_particle_par_e0;
		double largest_particle_par_a_err;
		double largest_particle_par_b_err;
		double largest_particle_par_e0_err;
		double largest_particle_par_prob;
		double largest_particle_par_chisq;
		double largest_particle_par_ndf;		
		
		double rough_primary_theta_z;
		double primary_phi;
		double primary_theta;
		
		double weighted_phi;
		double weighted_theta;
		
		
		double prim_par_a;
		double prim_par_b;
		double prim_par_e0;
		double prim_par_a_err;
		double prim_par_b_err;
		double prim_par_e0_err;
		double prim_par_prob;
		double prim_par_chisq;
		double prim_par_ndf;
		double prim_vise;
	
		double prim_pred_e_a;
		double prim_pred_g_a;
		double prim_pred_b;
		double prim_pred_e0;
		double prim_pred_e_chisq;
		double prim_pred_e_ndf;
		double prim_pred_g_chisq;
		double prim_pred_g_ndf;

		double prim_cmp_chisq;
		int prim_cmp_ndf;
		double prim_peakdiff;	

		double largest_particle_cmp_chisq;
		int largest_particle_cmp_ndf;
		double largest_particle_peakdiff;	

		double longest_s_particle_cmp_chisq;
		int longest_s_particle_cmp_ndf;
		double longest_s_particle_peakdiff;	

			
		double prim_pred_pre_over;		
		double prim_pred_pre_under;		
		double prim_pred_post_over;		
		double prim_pred_post_under;	
	
		
		double primary_long_e;
		double total_long_e;
		double total_long_e_frac;
		
		double pointing_phi;
		double pointing_theta;
		double pointing_r;
		double pointing_u;
		double pointing_v;
		double pointing_z;
		
		
		double mol_rad_r;
		double rms_r;
		
		double frac_particle_2;
		
		double elec_muon_asym;
		double elec_muon_asym_cale_weight;
		double elec_muon_asym_vise_weight;
		
		double elec_other_asym;
		double elec_other_asym_cale_weight;
		double elec_other_asym_vise_weight;
		
		double longest_z;
		double longest_particle_type;
		double longest_particle_vise;
		double longest_particle_avg_rms;
		
		double length_mean;
		double length_rms;
		double length_weighted_mean;
		double length_weighted_rms;
		
		double maxe_phi;
		double maxe_theta;
		double maxe_phi_rms;
		double maxe_theta_rms;
		
	
		double prim_pp_chisq;//!
		double prim_pp_p;//!
		
		
		int prim_pp_ndf;//!
		int prim_pp_igood;//!

		int ntot;
		int nelec;
		int nmuon;
		int nprot;
		int nother;
		int nneut;
		
		int nshort;
		int nmed;
		int nlong;
		
		
		
		

  
	private:
		void Init();


        ClassDef(Particles,2)
				
		
		
};

#endif

