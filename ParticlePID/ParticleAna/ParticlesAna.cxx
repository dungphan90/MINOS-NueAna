#include "NueAna/ParticlePID/ParticleAna/ParticlesAna.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle3D.h"
#include "TMath.h"


TH1F * ParticlesAna::rhist=0;
TH1F * ParticlesAna::viewTheta=0;
TH1F * ParticlesAna::viewPhi=0;


ParticlesAna::ParticlesAna()
{
	if(!rhist)rhist=new TH1F("rhist","rhist",1000,0,8);
	
	if(!viewTheta)viewTheta = new TH1F("viewTheta","Theta",100,0,2*3.141592);
	if(!viewPhi)viewPhi = new TH1F("viewPhi","Phi", 100,0,3.141592);
}

ParticlesAna::~ParticlesAna()
{
//	if(rhist){
//		delete rhist;
//		rhist=0;
//	}
}


void ParticlesAna::ana(ParticleObjectHolder * poh, Particles * e)
{

	rhist->Reset();
	viewTheta->Reset();
	viewPhi->Reset();


//	printf("------------------\nParticlesAna\n");

	

	int ntot = (int)poh->particles3d.size();//1->GetEntries();
//	int ntot = poh->particles3d.size();
	if(ntot==0)return;

	e->ntot=ntot;
	
	double largest_e=0;
	int largest_e_idx=-1;
	int largest_electron_idx=-1;
	double largest_elec_e=0;
	int second_largest_e_idx=-1;
	double second_largest_e=0;
	
	double sum_extheta=0;
	double sum_exphi=0;
	double sum_e=0;
	
	double pointing_u=0;
	double pointing_v=0;
	double pointing_z=0;
	
	double emvise=0;
	double emcale=0;
	
	double largest_s=0;
	
	//printf("\n\ntot %d\n",ntot);
	
	
	for(unsigned int i=0;i<poh->particles3d.size();i++)
	{
		
		Particle3D * p = (Particle3D *)&(poh->particles3d[i]);//1->At(i);
		if(!p)continue;
		
		if(p->entries <1)continue;
	
		//printf("particle %d with e %f\n",i,p->sum_e);	
		
		if(p->particletype==Particle3D::neutron)//neutrons mess it up... they are only 2d
		{
			e->nneut++;
			e->totcale+=p->calibrated_energy;
			e->totvise+=p->sum_e;
			e->neut_vise+=sum_e;
			

			//but do allow neutrons to be second largest energy particle
	                if(p->sum_e < largest_e && p->sum_e > second_largest_e)
        	        {
               	        	second_largest_e_idx=i;
               	        	second_largest_e=p->sum_e;
                		//printf("neut\n");
			}

		
			continue;
		}else if(p->particletype==Particle3D::electron)
		{
			e->nelec++;
			e->elec_vise+=p->sum_e;
			
			if(p->sum_e>largest_elec_e)
			{
				largest_elec_e=p->sum_e;
				largest_electron_idx=i;
			
			}
			
		}else if(p->particletype==Particle3D::muon) 
		{
			e->nmuon++;
			e->muon_vise+=p->sum_e;
		}else if(p->particletype==Particle3D::proton)
		{
			e->nprot++;
			e->prot_vise+=p->sum_e;
		}
		else if(p->particletype==Particle3D::other)
		{
			e->other_vise+=p->sum_e;
			e->nother++;
			
			
			//for now... consider it an electron
			if(p->sum_e>largest_elec_e)
			{
				largest_elec_e=p->sum_e;
				largest_electron_idx=i;
			
			}
			
		}else
		{
			std::cout<<"particle of unknown type! "<<p->particletype<<"\n";
			continue;
		}
		
		double plen = p->end_z-p->start_z;
		plen=plen<0?-plen:plen;
		
		
		e->length_mean+=plen;
		e->length_weighted_mean+=plen*p->sum_e;
		
		if( e->longest_z < plen)
		{
			e->longest_z = plen;
			e->longest_particle_type=p->particletype;
			e->longest_particle_vise=p->sum_e;
			e->longest_particle_avg_rms=p->avg_rms_t; 
		}
		

				
		double particle_s=0;		
		for(int j=0;j<p->entries;j++)
		{
			double u=p->u[j]-poh->event.vtx_u;
			double v=p->v[j]-poh->event.vtx_v;
			
		//	printf(" ??? %f %f %f\n",u,v,p->e[j]);
			
			if(p->e[j]>0 && p->e[j]<1000000)
			rhist->Fill(TMath::Sqrt(u*u+v*v), p->e[j]);		
			
					
			e->rms_r += u*u+v*v;
			
			if(j>0)
			{
				double du = p->u[j]-p->u[j-1];
				double dv = p->v[j]-p->v[j-1];
				double dz = p->z[j]-p->z[j-1];
				
				double dd=sqrt(du*du+dv*dv+dz*dz);
				if(dd<100)particle_s+=dd; //in case we have a single bogus value, we will skip it
			}
		
		}
		
		
		
		if(p->particletype == Particle3D::electron || p->particletype==Particle3D::muon)
		{
			e->elec_muon_asym += (p->particletype==Particle3D::electron?+1:-1);
			e->elec_muon_asym_cale_weight+= (p->particletype==Particle3D::electron?+1:-1) * p->calibrated_energy;
			e->elec_muon_asym_vise_weight+= (p->particletype==Particle3D::electron?+1:-1) * p->sum_e;
			emvise+=p->sum_e;
			emcale+=p->calibrated_energy;
			
		}
		
		e->elec_other_asym= (p->particletype==Particle3D::electron?+1:-1);
		e->elec_other_asym_cale_weight= (p->particletype==Particle3D::electron?+1:-1) * p->calibrated_energy;
		e->elec_other_asym_vise_weight= (p->particletype==Particle3D::electron?+1:-1) * p->sum_e;

		
		
		
		
		
		
		e->totcale+=p->calibrated_energy;
		e->totvise+=p->sum_e;
		
	//	printf("part %d sume %f\n", i,p->sum_e);
	//
	
		double du = p->end_u-p->start_u;
		double dv = p->end_v-p->start_v;
		double dz = p->end_z-p->start_z;
		double r = TMath::Sqrt(du*du + dv*dv+dz*dz);
	
		if(r)
		{	
			pointing_u+=du/r*p->sum_e;
			pointing_v+=dv/r*p->sum_e;
			pointing_z+=dz/r*p->sum_e;
		}

		double longe = p->sum_e * TMath::Sin(p->phi);
	
		e->total_long_e+=longe;
		
		if(p->calibrated_energy==0 && p->sum_e>0)p->calibrated_energy=p->sum_e*60;
			
		if(p->sum_e < largest_e && p->sum_e > second_largest_e)
		{
			second_largest_e_idx=i;
			second_largest_e=p->sum_e;			
		}	

	
		if(p->sum_e > largest_e)
		{

			second_largest_e_idx=largest_e_idx;
			second_largest_e=largest_e;
			largest_e_idx=i;
			largest_e = p->sum_e;
			e->largest_particle_avg_rms=p->avg_rms_t;
			e->primary_long_e=longe;
	
			e->largest_particle_s = particle_s;
			e->largest_particle_z=p->end_z-p->start_z;
			e->largest_particle_par_b = p->emfit_b;
			e->largest_particle_par_a = p->emfit_a;
			e->largest_particle_par_e0 = p->emfit_e0;
			e->largest_particle_par_a_err = p->emfit_a_err;
			e->largest_particle_par_b_err = p->emfit_b_err;
			e->largest_particle_par_e0_err = p->emfit_e0_err;
			e->largest_particle_par_prob = p->emfit_prob;
			e->largest_particle_par_chisq = p->emfit_chisq;
			e->largest_particle_par_ndf = p->emfit_ndf;

			e->largest_particle_cmp_chisq = p->cmp_chisq;
			e->largest_particle_cmp_ndf = p->cmp_ndf;
			e->largest_particle_peakdiff = p->peakdiff;


		}
	

		//printf("largest %d %f 2nd %d %f\n",largest_e_idx,largest_e,second_largest_e_idx,second_largest_e);

	
		
		if(particle_s > largest_s)
		{
		
			largest_s=particle_s;
			
			
			e->longest_s_particle_s=particle_s;		

			
		
			
			e->longest_s_particle_e=p->sum_e;
			e->longest_s_particle_cal_e=p->calibrated_energy;
			e->longest_s_particle_type=p->particletype;
			e->longest_s_particle_avg_e=p->sum_e/particle_s;
		


			e->longest_s_particle_avg_rms=p->avg_rms_t;
			e->longest_s_long_e=longe;

			e->longest_s_particle_z=p->end_z-p->start_z;
			e->longest_s_particle_par_b = p->emfit_b;
			e->longest_s_particle_par_a = p->emfit_a;
			e->longest_s_particle_par_e0 = p->emfit_e0;
			e->longest_s_particle_par_a_err = p->emfit_a_err;
			e->longest_s_particle_par_b_err = p->emfit_b_err;
			e->longest_s_particle_par_e0_err = p->emfit_e0_err;
			e->longest_s_particle_par_prob = p->emfit_prob;
			e->longest_s_particle_par_chisq = p->emfit_chisq;
			e->longest_s_particle_par_ndf = p->emfit_ndf;		


			e->longest_s_particle_cmp_chisq = p->cmp_chisq;
			e->longest_s_particle_cmp_ndf = p->cmp_ndf;
			e->longest_s_particle_peakdiff = p->peakdiff;


	

		
		}


		
		
		if(p->entries < 10)e->nshort++;
		if(p->entries >= 10 && p->entries < 20) e->nmed++;
		if(p->entries >=20)e->nlong++;
	
	
		sum_extheta+=p->sum_e * p->theta;
		sum_exphi+=p->sum_e * p->phi;
		sum_e+=p->sum_e;
	}



	e->pointing_u=pointing_u;
	e->pointing_v=pointing_v;
	e->pointing_z=pointing_z;

	if(emcale>0)e->elec_muon_asym_cale_weight/=emcale;else e->elec_muon_asym_cale_weight=0;
	if(emvise>0)e->elec_muon_asym_vise_weight/=emvise;else e->elec_muon_asym_vise_weight=0;
	
	if(e->totcale>0)e->elec_other_asym_cale_weight/=e->totcale; else e->elec_other_asym_cale_weight=0;
	if(e->totvise>0)e->elec_other_asym_vise_weight/=e->totvise; else e->elec_other_asym_vise_weight=0;


	if(e->totvise>0)e->length_weighted_mean/=e->totvise;else e->length_weighted_mean=0;



	int tcount=0;
	for(int i=0;i<ntot;i++)
	{
		
		Particle3D * p = (Particle3D *)&(poh->particles3d[i]);//1->At(i);
		if(p <=0)continue;
		
		if(p->entries <1)continue;
		if(p->particletype==Particle3D::neutron)continue;
		double plen = p->end_z-p->start_z;
							
		e->length_weighted_rms+=(plen-e->length_weighted_mean)*(plen-e->length_weighted_mean);
		e->length_rms+=plen*plen;
		tcount++;
	}				
					
	if(tcount>0)e->length_weighted_rms=sqrt(e->length_weighted_rms/tcount);else e->length_weighted_rms=0;
	if(tcount>0)e->length_rms=sqrt(e->length_rms/tcount);else e->length_rms=0;



	e->pointing_r=TMath::Sqrt(pointing_u*pointing_u*pointing_v*pointing_v*pointing_z*pointing_z);
	
	if(e->pointing_r>0)
	{
		e->pointing_phi=TMath::ACos(pointing_z/e->pointing_r);
		e->pointing_theta=TMath::ATan2(pointing_u,pointing_v);
	}


	if(e->totcale>0)
	e->total_long_e_frac = e->total_long_e/e->totcale;

	if(sum_e>0)
	{
		e->weighted_phi=sum_exphi/sum_e;
		e->weighted_theta=sum_extheta/sum_e;
	}

	if(largest_e_idx>-1)
	{
		
		Particle3D * p = (Particle3D *)&(poh->particles3d[largest_e_idx]);
		if(p)
		{
			
		e->largest_particle_e=p->sum_e;
		e->largest_particle_cal_e=p->calibrated_energy;
		e->largest_particle_type=p->particletype;
		
		
		e->rough_primary_theta_z = p->start_z-p->end_z ? TMath::ATan(TMath::Sqrt((p->end_u-p->start_u)*(p->end_u-p->start_u)+(p->end_v-p->start_v)*(p->end_v-p->start_v))/TMath::Abs(p->start_z-p->end_z)) : 1.5707;//if its vertical, don't compute it
      	e->primary_phi=p->phi;
      	e->primary_theta=p->theta;
		
		
		e->largest_particle_avg_e= p->end_z-p->start_z ? p->sum_e / (p->end_z-p->start_z) : p->sum_e;
		

		if(second_largest_e_idx>-1)
		{
			Particle3D * p2 = (Particle3D *)&(poh->particles3d[second_largest_e_idx]);
			if(p2)
			{
				if(e->largest_particle_e>0)
				e->frac_particle_2 =   p2->sum_e / e->largest_particle_e;
		//		printf("fp2 %f from %f / %f\n",e->frac_particle_2,p2->sum_e ,e->largest_particle_e);
			}
		}
		}





		
		for(unsigned int i=0;i<p->e.size();i++)
		{
			double u = p->u[i]-poh->event.vtx_u;
			double v = p->v[i]-poh->event.vtx_v;
			double z = p->z[i]-poh->event.vtx_z;
			if(!u || !v || !z)continue;
	
			double r = sqrt(u*u+v*v+z*z);
			double theta = atan(u/v);
			double phi = acos(z/r);
			
			phi=phi<0?phi+2*3.141592:phi;
			theta=theta<0?theta+2*3.141592:theta;
		
			viewPhi->Fill(phi,p->e[i]);
			viewTheta->Fill(theta,p->e[i]);
		}

		e->maxe_phi=viewPhi->GetMean();
		e->maxe_theta=viewTheta->GetMean();
		e->maxe_phi_rms=viewPhi->GetRMS();
		e->maxe_theta_rms=viewTheta->GetRMS();


	}
	
	
	
	if(largest_electron_idx>-1)
	{
		Particle3D * p = (Particle3D *)&(poh->particles3d[largest_e_idx]);
		if(p)
		{
			e->emfrac = p->sum_e / e->totvise;
			e->prim_par_a = p->emfit_a;
			e->prim_par_b = p->emfit_b;
			e->prim_par_e0 = p->emfit_e0;
			e->prim_par_a_err = p->emfit_a_err;
			e->prim_par_b_err = p->emfit_b_err;
			e->prim_par_e0_err = p->emfit_e0_err;
			e->prim_par_prob = p->emfit_prob;
			e->prim_par_chisq=p->emfit_chisq;
			e->prim_par_ndf=p->emfit_ndf;
			e->prim_vise = p->sum_e;
			
			
			e->prim_pred_e_a=p->pred_e_a;
			e->prim_pred_g_a=p->pred_g_a;
			e->prim_pred_b=p->pred_b;
			e->prim_pred_e0=p->pred_e0;
			e->prim_pred_e_chisq=p->pred_e_chisq;
			e->prim_pred_e_ndf=p->pred_e_ndf;
			e->prim_pred_g_chisq=p->pred_g_chisq;
			e->prim_pred_g_ndf=p->pred_g_ndf;
		
			e->prim_pred_pre_over=p->pre_over;
			e->prim_pred_pre_under=p->pre_under;
			e->prim_pred_post_over=p->post_over;	
			e->prim_pred_post_under=p->post_under;		
			
			e->prim_pp_chisq=p->pp_chisq;
			e->prim_pp_ndf=p->pp_ndf;
			e->prim_pp_igood=p->pp_igood;
			e->prim_pp_p=p->pp_p;	

			e->prim_cmp_chisq = p->cmp_chisq;
			e->prim_cmp_ndf = p->cmp_ndf;
			e->prim_peakdiff = p->peakdiff;

	//	printf("%f %d %f\n",p->cmp_chisq,p->cmp_ndf, p->peakdiff);

		}
	}	
		
	

	double rhist_e = rhist->Integral();
	double rhist_partial_e=0;
	
//	printf("rhiste %f \n",rhist_e);
	
	if(rhist_e>0)
	for(int i=0;i<rhist->GetNbinsX();i++)
	{
			if(rhist_partial_e/rhist_e < 0.9)
			{
				rhist_partial_e+=rhist->GetBinContent(i);
				
	//			printf ("-- %f %f -- %f\n",rhist_partial_e,rhist->GetBinContent(i),(double)i / (double) rhist->GetNbinsX());
				
				e->mol_rad_r=(double)i / (double) rhist->GetNbinsX(); ///right now its normalized!
			}

	}
	

	e->rms_r = TMath::Sqrt(e->rms_r);
	
}


