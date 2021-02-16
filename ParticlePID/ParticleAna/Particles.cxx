#include "NueAna/ParticlePID/ParticleAna/Particles.h"

ClassImp(Particles)


Particles::Particles() 
{
	Init();
}

Particles::~Particles()
{}






void Particles::Clear(Option_t* /* option */) {
  // Purpose: Clear memory allocated to arrays so that record can
  // be reused.  

	Init();

}
 
 
 
void Particles::Init() {
  
	emfrac=0;
	totcale=0;
	totvise=0;
	
	elec_vise=0;
	muon_vise=0;
	prot_vise=0;
	neut_vise=0;
	other_vise=0;


	 longest_s_particle_s=0;
	 longest_s_particle_e=0;
	 longest_s_particle_cal_e=0;
	 longest_s_particle_type=0;
	 longest_s_particle_avg_rms=0;
	 longest_s_particle_avg_e=0;
	
	 longest_s_particle_z=0;
	 longest_s_particle_par_b=0;
	 longest_s_particle_par_a=0;
	 longest_s_particle_par_e0=0;
	 longest_s_particle_par_a_err=0;
	 longest_s_particle_par_b_err=0;
	 longest_s_particle_par_e0_err=0;
	 longest_s_particle_par_prob=0;
	 longest_s_particle_par_chisq=0;
	 longest_s_particle_par_ndf=0;	
	 longest_s_long_e=0;


	largest_particle_s=0;
	largest_particle_e=0;
	largest_particle_cal_e=0;
	largest_particle_type=0;
	largest_particle_avg_rms=0;
	largest_particle_avg_e=0;
	
	largest_particle_z=0;
	largest_particle_par_b=0;
	largest_particle_par_a=0;
	largest_particle_par_e0=0;
	largest_particle_par_a_err=0;
	largest_particle_par_b_err=0;
	largest_particle_par_e0_err=0;
	largest_particle_par_prob=0;
	largest_particle_par_chisq=0;
	largest_particle_par_ndf=0;	
	
	rough_primary_theta_z=0;
	primary_phi=0;
	primary_theta=0;
	
	weighted_phi=0;
	weighted_theta=0;
	
	
	prim_par_a=0;
	prim_par_b=0;
	prim_par_e0=0;
	prim_par_a_err=0;
	prim_par_b_err=0;
	prim_par_e0_err=0;
	prim_par_prob=0;
	prim_par_chisq=0;
	prim_par_ndf=0;
	prim_vise=0;
	
	prim_pred_e_a=0;
	prim_pred_g_a=0;
	prim_pred_b=0;
	prim_pred_e0=0;
	prim_pred_e_chisq=0;
	prim_pred_e_ndf=0;
	prim_pred_g_chisq=0;
	prim_pred_g_ndf=0;
		
		
	prim_pred_pre_over=0;	
	prim_pred_pre_under=0;	
	prim_pred_post_over=0;	
	prim_pred_post_under=0;	
	
	
	primary_long_e=0;
	total_long_e=0;
	total_long_e_frac=0;
	
	pointing_phi=0;
	pointing_theta=0;
	pointing_r=0;
	pointing_u=0;
	pointing_v=0;
	pointing_z=0;
	
	
	mol_rad_r=0;
	rms_r=0;
	
	frac_particle_2=0;
	
	elec_muon_asym=0;
	elec_muon_asym_cale_weight=0;
	elec_muon_asym_vise_weight=0;
	
	elec_other_asym=0;
	elec_other_asym_cale_weight=0;
	elec_other_asym_vise_weight=0;
	
	longest_z=0;
	longest_particle_type=0;
	longest_particle_vise=0;
	longest_particle_avg_rms=0;
	
	length_mean=0;
	length_rms=0;
	length_weighted_mean=0;
	length_weighted_rms=0;
	
	maxe_phi=0;
	maxe_theta=0;
	maxe_phi_rms=0;
	maxe_theta_rms=0;
	
	
	prim_pp_chisq=0;
	prim_pp_p=0;
	
	
	prim_pp_ndf=0;
	prim_pp_igood=0;
	
	ntot=0;
	nelec=0;
	nmuon=0;
	nprot=0;
	nother=0;
	nneut=0;
	
	nshort=0;
	nmed=0;
	nlong=0;


	prim_cmp_chisq=0;
	prim_cmp_ndf=0;
	prim_peakdiff=0;	

	largest_particle_cmp_chisq=0;
	largest_particle_cmp_ndf=0;
	largest_particle_peakdiff=0;	

	longest_s_particle_cmp_chisq=0;
	longest_s_particle_cmp_ndf=0;
	longest_s_particle_peakdiff=0;	


  
}




