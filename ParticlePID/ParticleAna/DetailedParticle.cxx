#include "NueAna/ParticlePID/ParticleAna/DetailedParticle.h"
#include "TH3F.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle3D.h"

#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"

#include "MCNtuple/NtpMCStdHep.h"

#include "NueAna/ParticlePID/ParticleAna/TruthCompareAna.h"
#include "MessageService/MsgService.h"

DetailedParticle::DetailedParticle(const char*fname)
{
	input=new TChain("PO");
	inputPA=new TChain("PA");
	input->Add(fname);
	inputPA->Add(fname);
	poh=new ParticleObjectHolder();
	prec=new PRecord();
	input->SetBranchAddress("ParticleObjectHolder",&poh);
	inputPA->SetBranchAddress("PRecord",&prec);
}

DetailedParticle::~DetailedParticle()
{}


int DetailedParticle::PassesPreselec()
{
	int pass=1;
	
//	pass = pass && prec->event.inFiducial==1;


	pass = pass && ( (prec->event.vtx_z>0.5 && prec->event.vtx_z<14.5) || (prec->event.vtx_z>16.5 && prec->event.vtx_z<29.5) );

	double r = sqrt(prec->event.vtx_u*prec->event.vtx_u+prec->event.vtx_v*prec->event.vtx_v);
	
	pass = pass && (r>0.25 && r<3.9);
	
	pass = pass && prec->event.nstrips>10;
	pass = pass && prec->event.nstrips<50;
	pass = pass && prec->particles.longest_s_particle_s<0.9;
	pass = pass && prec->particles.emfrac>0.4;
	pass = pass && prec->particles.elec_vise/25.>1;
	pass = pass && prec->particles.primary_long_e<40;
	pass = pass && prec->particles.largest_particle_e/prec->particles.totvise>0.65;
/*	pass = pass && prec->event.visenergy/25.-prec->particles.totvise/25>0;
	pass = pass && prec->event.visenergy/25.-prec->particles.totvise/25<0.5;
*/	
	
	/*
pass = pass && prec->particles.ntot>0;
pass = pass && prec->particles.longest_s_particle_s<1.1;
pass = pass && prec->particles.elec_vise<200;
pass = pass && prec->particles.nelec>0;
pass = pass && prec->particles.largest_particle_type==11;
pass = pass && prec->event.nstrips>10;
pass = pass && prec->event.nstrips<55;
pass = pass && prec->particles.largest_particle_avg_rms>0.004;
*/

/*
	pass = pass && prec->event.nstrips>10 && prec->event.nstrips<50;	
	pass = pass && prec->particles.ntot>0;
	pass = pass &&  prec->event.vtx_z-prec->event.min_z<0.35;
	pass = pass &&  prec->event.nstrips>15 && prec->event.nstrips<45;
	pass = pass &&  prec->particles.rms_r<1.1;
	pass = pass &&  prec->particles.prim_par_b>0.1&&prec->particles.prim_par_b<1.2;
	pass = pass &&  prec->particles.ntot/prec->event.visenergy<0.14;
	pass = pass &&  prec->particles.mol_rad_r<0.0425;
	pass = pass &&  prec->particles.largest_particle_avg_rms<0.017;
	pass = pass &&  prec->particles.largest_particle_avg_rms>0.005;
	pass = pass &&  prec->particles.ntot<8 ;
	pass = pass &&  prec->particles.longest_z<0.8;
	pass = pass &&  prec->particles.primary_long_e<60;
	
*/	

	
/*	
	pass = pass && prec->particles.rough_primary_theta_z<0.7;
	pass = pass && prec->particles.primary_long_e<15;
	
	pass = pass && sqrt((prec->event.max_u-prec->event.min_u)*(prec->event.max_u-prec->event.min_u)
			+(prec->event.max_v-prec->event.min_v)*(prec->event.max_v-prec->event.min_v))>0.3;
			
	pass = pass && prec->particles.largest_particle_avg_rms>0.0045;
	pass = pass && prec->particles.largest_particle_e/prec->event.visenergy>0.85;

	pass = pass && prec->particles.total_long_e<25;
	pass = pass && prec->particles.mol_rad_r>0.001;
	pass = pass && prec->particles.mol_rad_r<0.035;

	pass = pass && prec->particles.nelec>0;
	pass = pass && prec->particles.ntot<4;
	pass = pass && prec->particles.prim_par_b>0.1;
	pass = pass && prec->particles.prim_par_b<1;

	pass = pass && prec->particles.longest_z<1;
	pass = pass && prec->particles.longest_z>0.1;
	pass = pass && prec->particles.longest_particle_type<14;
*/
	return pass;
}

void DetailedParticle::Run()
{

	
	TH2F * hist_sumdist[5];
	TH2F * hist_weighted_sumdist[5]; 
	TH1F * hist_max_dist[5]; 
	TH1F * hist_max_dist_per_gev[5]; 	


	TH1F * hist_weighted_area[5]; 	


	TH1F * hist_first_peak[5];
	TH1F * hist_second_peak[5];
	TH1F * hist_twopeak_diff[5];
	
	TH1F * hist_back_peak_dist[5];
	
	TH1F * hist_view_e_asym[5];		


	TH2F * hist_prim_angle_vs_prim_eres[5];		
	TH2F * hist_prim_angle_vs_prim_true[5];		
	
	TH2F * hist_sum_curv_vs_part_e[5];
	TH2F * hist_pu_pz[5];
	TH2F * hist_pv_pz[5];
	TH2F * hist_pt_pz[5];
	TH2F * hist_pu_pv[5];
	
	for(int i=0;i<5;i++)
	{
		char tmp[200];
		
		sprintf(tmp,"sumdist_%d",i);
		hist_sumdist[i] = new TH2F(tmp,tmp,100,0, 50,20,0,20);

		sprintf(tmp,"weighted_sumdist_%d",i);		
		hist_weighted_sumdist[i] = new TH2F(tmp,tmp,100,0, 50,20,0,20);

		sprintf(tmp,"max_dist_%d",i);			
		hist_max_dist[i] = new TH1F(tmp,tmp,100,0,20);
		
		sprintf(tmp,"max_dist_per_gev_%d",i);		
		hist_max_dist_per_gev[i] = new TH1F(tmp,tmp,100,0,1);
	
		sprintf(tmp,"weighted_area_%d",i);		
		hist_weighted_area[i] = new TH1F(tmp,tmp,100,0,20);

		sprintf(tmp,"first_peak_%d",i);		
		hist_first_peak[i] = new TH1F(tmp,tmp,100,0,2);

		sprintf(tmp,"second_peak_%d",i);		
		hist_second_peak[i] = new TH1F(tmp,tmp,100,0,2);
		
		sprintf(tmp,"twopeak_diff_peak_%d",i);		
		hist_twopeak_diff[i] = new TH1F(tmp,tmp,100,0,2);
				
		sprintf(tmp,"back_peak_dist_%d",i);		
		hist_back_peak_dist[i] = new TH1F(tmp,tmp,100,0,2);	
		
		sprintf(tmp,"view_e_asym_%d",i);		
		hist_view_e_asym[i] = new TH1F(tmp,tmp,100,0,1);	

		sprintf(tmp,"prim_angle_vs_prim_eres_%d",i);
		hist_prim_angle_vs_prim_eres[i] = new TH2F(tmp,tmp,100,0, 3,100,-1,1);		
		sprintf(tmp,"prim_angle_vs_prim_true_%d",i);
		hist_prim_angle_vs_prim_true[i] = new TH2F(tmp,tmp,100,0, 3,100,0,10);	


		sprintf(tmp,"sum_curv_vs_part_e_%d",i);
		hist_sum_curv_vs_part_e[i]=new TH2F(tmp,tmp,1000,0,50,100,0,10);
		sprintf(tmp,"pu_pz_%d",i);
		hist_pu_pz[i]=new TH2F(tmp,tmp,100,-10,10,100,0,10);
		sprintf(tmp,"pv_pz_%d",i);
		hist_pv_pz[i]=new TH2F(tmp,tmp,100,-10,10,100,0,10);
		sprintf(tmp,"pt_pz_%d",i);
		hist_pt_pz[i]=new TH2F(tmp,tmp,100,0,10,100,0,10);

		sprintf(tmp,"pu_pv_%d",i);
		hist_pu_pv[i]=new TH2F(tmp,tmp,100,-10,10,100,-10,10);
				
	}
	




	int ent =input->GetEntries();
	int entpa =inputPA->GetEntries();
	printf("%d %d total entries\n",ent,entpa);
	
	
	double typecnt[5];
	for(int i=0;i<5;i++)typecnt[i]=0;
	
	
	int saved=0;
	int ji=0;
	input->GetEntry(0);
	for(int ent_idx=0;ent_idx<entpa;ent_idx++)
	{
		//input->GetEntry(i);
		/*if(i==0)*/inputPA->GetEntry(ent_idx);
		
		while(poh->GetHeader().GetRun() != prec->GetHeader().GetRun() || poh->GetHeader().GetSnarl()!=prec->GetHeader().GetSnarl()  || poh->GetHeader().GetEvent()!=prec->GetHeader().GetEvent())
			input->GetEntry(ji++);
			
/*
		int j=0;
		while(poh->GetHeader().GetRun() != prec->GetHeader().GetRun())
		{
			j++;
			inputPA->GetEntry(i+j);
	//					printf("sync a %d %d\n",poh->GetHeader().GetSnarl(),prec->GetHeader().GetSnarl());
		}
		

		while(poh->GetHeader().GetSnarl()<prec->GetHeader().GetSnarl() && poh->GetHeader().GetRun()==prec->GetHeader().GetRun())
		{
			i++;
			input->GetEntry(i);
	//				printf("sync b %d %d\n",poh->GetHeader().GetSnarl(),prec->GetHeader().GetSnarl());
		}
		
		j=0;
		while(poh->GetHeader().GetSnarl()>prec->GetHeader().GetSnarl() && poh->GetHeader().GetRun()==prec->GetHeader().GetRun())
		{
			j++;
			inputPA->GetEntry(i+j);
	//				printf("sync c %d %d\n",poh->GetHeader().GetSnarl(),prec->GetHeader().GetSnarl());
		}
		
		if(poh->GetHeader().GetRun()!=prec->GetHeader().GetRun())continue;
		*/
		saved++;
		
/*
		while(poh->GetHeader().GetSnarl()>prec->GetHeader().GetSnarl() && poh->GetHeader().GetRun()==prec->GetHeader().GetRun())
		{
			printf("sync %d %d\n",poh->GetHeader().GetSnarl(),prec->GetHeader().GetSnarl());
			j++;
			if(i+j>inputPA->GetEntries()-1)break;
			
			inputPA->GetEntry(i+j);
		}
		
		
		if(i+j>inputPA->GetEntries()-1)break;
		
		if(poh->GetHeader().GetSnarl()<prec->GetHeader().GetSnarl() || poh->GetHeader().GetRun()!=prec->GetHeader().GetRun())
		{
			printf("sync %d %d\n",poh->GetHeader().GetSnarl(),prec->GetHeader().GetSnarl());

			continue;
		}
*/		
		
		if(!PassesPreselec())continue;
	
	
	
		double weight = poh->mctrue.totbeamweight*poh->mctrue.oscprob*3.25/6.5/5;
		
		double paweight = prec->mctrue.totbeamweight*prec->mctrue.oscprob*3.25/6.5/5;
		
		if(fabs(weight-paweight)>0.001)printf("weight mismatch!\n");
		
		int type=poh->mctrue.type;
	
		if(type<0 || type>5)continue;
		
		if(type==0)weight=weight/3.;
	
	
	
	//	MsgService * ms = MsgService::Instance();
	//	ms->GetStream("TruthCompareAna")->SetLogLevel(Msg::kDebug);
	//	TruthCompareAna a;
	//	a.ana(poh,&prec->truthcompare);
	
	/*
		////calculate expected particles.......
		printf("\nstdhep....run %d snarl %d event %d\n",prec->GetHeader().GetRun() , prec->GetHeader().GetSnarl() , prec->GetHeader().GetEvent());

		double true_elec_p4[4];//largest true election
		true_elec_p4[3]=0;
		int true_elec_idx=-1;
		
		for(unsigned int i=0;i<prec->mctrue.stdhep.size();i++)
		{
			NtpMCStdHep *mc = &(prec->mctrue.stdhep[i]);
			
			if(mc->IstHEP!=1)continue;
			
			printf("found %d with e %f at %d  (%f %f %f)  p3 (%f %f %f)\n",mc->IdHEP,sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass),mc->IstHEP,(mc->vtx[0]+mc->vtx[1])/sqrt(2.0),(-mc->vtx[0]+mc->vtx[1])/sqrt(2.0),mc->vtx[2],(mc->p4[0]+mc->p4[1])/sqrt(2.0),(-mc->p4[0]+mc->p4[1])/sqrt(2.0),mc->p4[2]);
			
			//need lepton energy from cc
			//need sum of pi0 energy
			//number of pi0
			
			//number of pip,pim
			
			if(fabs(mc->IdHEP)==11 && true_elec_p4[3]<mc->p4[3])
			{
				true_elec_p4[0]=mc->p4[0];
				true_elec_p4[1]=mc->p4[1];
				true_elec_p4[2]=mc->p4[2];
				true_elec_p4[3]=mc->p4[3];
				true_elec_idx=i;
			}
		
		}
		printf("\n");
		
		
		///////////////////////
	*/
	
	
		//calculate view momentums of found particles relative to reco vertex
		
		double pu=0;
		double pv=0;
		double pz_u=0;
		double pz_v=0;
		double pz=0;
		double pt=0;			
		for(int j=0;j<(int)poh->particles3d.size();j++)
		{

			Particle3D * p = &poh->particles3d[j];
			if(!p)continue;
			
			for(int i=0;i<p->entries;i++)
			{
				if(p->view[i]==2)
				{
					double dt=p->u[i]-poh->event.vtx_u;
					double dz=p->z[i]-poh->event.vtx_z;
					double r = sqrt(dt*dt+dz*dz);
					if(r>0)pu+= p->e[i]/25.*dt/r;
					if(r>0)pz_u+= p->e[i]/25.*dz/r;
				}

				if(p->view[i]==2)
				{
					double dt=p->v[i]-poh->event.vtx_v;
					double dz=p->z[i]-poh->event.vtx_z;
					double r = sqrt(dt*dt+dz*dz);
					if(r>0)pv+= p->e[i]/25.*dt/r;
					if(r>0)pz_v+= p->e[i]/25.*dz/r;
				}			
			
				double du=p->u[i]-poh->event.vtx_u;
				double dv=p->v[i]-poh->event.vtx_v;
				double dz=p->z[i]-poh->event.vtx_z;
				double r = sqrt(du*du+dv*dv+dz*dz);
				double dt=sqrt(du*du+dv*dv);
				if(r>0)pt+=p->e[i]/25.*dt/r;
				if(r>0)pz+=p->e[i]/25.*dz/r;
			
			}
			
		}
		
		hist_pu_pz[type]->Fill(pu,pz_u,weight);
		hist_pv_pz[type]->Fill(pv,pz_v,weight);
		hist_pt_pz[type]->Fill(pt,pz,weight);
		hist_pu_pv[type]->Fill(pu,pv,weight);
		
		
		int longest_id=-1;
		double longZ=0;
		
		int largest_id=-1;
		double largeE=0;
				
		for(int j=0;j<(int)poh->particles3d.size();j++)
		{

			Particle3D * p = &poh->particles3d[j];
			if(!p)continue;
			

			
			if(p->sum_e>largeE)
			{
				largest_id=j;
				largeE=p->sum_e;
			}
			
			if(p->end_z-p->start_z>longZ)
			{
				longest_id=j;
				longZ=p->end_z-p->start_z;
			}
			

			//////////
			double sum_z2=0;
			double sum_z=0;
			double sum_u=0;
			double sum_zu=0;			
			double sum_v=0;
			double sum_zv=0;
			int n=0;	
						
			for(int k=0;k<p->entries;k++)
			{
				sum_z2+=p->z[k]*p->z[k];
				sum_z+=p->z[k];
				sum_zu+=p->z[k]*p->u[k];
				sum_u+=p->u[k];
				sum_zv+=p->z[k]*p->v[k];
				sum_v+=p->v[k];	
				n++;			
			}
			
			//double off_u=  (sum_u*sum_z2 - sum_z*sum_zu)/(n*sum_z2-sum_z*sum_z);
			double slope_u=  (n*sum_zu - sum_z*sum_u)/(n*sum_z2-sum_z*sum_z);
			
			
			//double off_v=  (sum_v*sum_z2 - sum_z*sum_zv)/(n*sum_z2-sum_z*sum_z);
			double slope_v=  (n*sum_zv - sum_z*sum_v)/(n*sum_z2-sum_z*sum_z);
			

	

			double dz=p->end_z-p->start_z;
			double du=dz*slope_u;
			double dv=dz*slope_v;			
			
			double r2 = du*du+dv*dv+dz*dz;
			r2=sqrt(r2);
			//double pu= p->sum_e/25. * du/r2;
			//double pv= p->sum_e/25. * dv/r2;
			//double pz= p->sum_e/25. * dz/r2;
			
			//printf("particle %d with e %f  (%f %f %f) p3 (%f %f %f)\n",p->particletype,p->sum_e/25.,p->start_u,p->start_v,p->start_z,pu,pv,pz);
	
	

		}
	
	
	
		
		if(largest_id>-1)
		{
			Particle3D * p = &poh->particles3d[largest_id];
		

			//calculate particle curv
			double curv=0;
			for(int i=1;i<p->entries;i++)
			{
				double du = p->u[i]-p->u[i-1];
				double dv = p->v[i]-p->v[i-1];
				double dz = p->z[i]-p->z[i-1];
				double ds = sqrt(du*du+dv*dv+dz*dz);
				if(dz>0)curv+=ds/dz;
			}
			hist_sum_curv_vs_part_e[type]->Fill(curv,p->sum_e/25.,weight);
	

			double sum_z2=0;
			double sum_z=0;
			double sum_u=0;
			double sum_zu=0;			
			double sum_v=0;
			double sum_zv=0;
			int n=0;	
						
			for(int i=0;i<p->entries;i++)
			{
				sum_z2+=p->z[i]*p->z[i];
				sum_z+=p->z[i];
				sum_zu+=p->z[i]*p->u[i];
				sum_u+=p->u[i];
				sum_zv+=p->z[i]*p->v[i];
				sum_v+=p->v[i];	
				n++;			
			}
			
			double off_u=  (sum_u*sum_z2 - sum_z*sum_zu)/(n*sum_z2-sum_z*sum_z);
			double slope_u=  (n*sum_zu - sum_z*sum_u)/(n*sum_z2-sum_z*sum_z);
			
			
			double off_v=  (sum_v*sum_z2 - sum_z*sum_zv)/(n*sum_z2-sum_z*sum_z);
			double slope_v=  (n*sum_zv - sum_z*sum_v)/(n*sum_z2-sum_z*sum_z);
			
/*
			////////////////////
			double dz=p->end_z-p->start_z;
			double du=dz*slope_u;
			double dv=dz*slope_v;			
			
			double r2 = du*du+dv*dv+dz*dz;
			r2=sqrt(r2);
			double pu= p->sum_e/25. * du/r2;
			double pv= p->sum_e/25. * dv/r2;
			double pz= p->sum_e/25. * dz/r2;


			if(p->particletype==Particle3D::electron)
			{
				printf("largest elec %d with e %f  (%f %f %f) p3 (%f %f %f)\n",p->particletype,p->sum_e/25.,p->start_u,p->start_v,p->start_z,pu,pv,pz);
				if(true_elec_idx>-1)
				{
					NtpMCStdHep *mc = &(prec->mctrue.stdhep[true_elec_idx]);
			
					printf("compare to true elec %d with e %f at %d  (%f %f %f)  p3 (%f %f %f)\n",mc->IdHEP,sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass),mc->IstHEP,(mc->vtx[0]+mc->vtx[1])/sqrt(2.0),(-mc->vtx[0]+mc->vtx[1])/sqrt(2.0),mc->vtx[2],(mc->p4[0]+mc->p4[1])/sqrt(2.0),(-mc->p4[0]+mc->p4[1])/sqrt(2.0),mc->p4[2]);		
					//calculate the angle!!!!
					
					//angle between p vectors....
					
					double true_pu=(mc->p4[0]+mc->p4[1])/sqrt(2.0);
					double true_pv=(-mc->p4[0]+mc->p4[1])/sqrt(2.0);
					double true_pz=mc->p4[2];
					
					double reco_r = sqrt(pu*pu+pv*pv+pz*pz);	
					double true_r = sqrt(true_pu*true_pu+true_pv*true_pv+true_pz*true_pz);
					
					double angle = acos( (pu*true_pu+pv*true_pv+pz*true_pz) / (reco_r*true_r) );
					printf("angle between: %f\n",angle);
					
					hist_prim_angle_vs_prim_eres[type]->Fill(angle,(sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass)-p->sum_e/25.)/(sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass)),weight);		
					hist_prim_angle_vs_prim_true[type]->Fill(angle,sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass));	
					
					
				}
			}
			
			/////////////////////
	*/		
			//calculate rms
			
			double sumdist=0;
			double weighted_sumdist=0;
			double sume=0;
			
			for(int i=0;i<p->entries;i++)
			{
				double pu = slope_u*p->z[i]+off_u;	
				double pv = slope_v*p->z[i]+off_v;
				
				sumdist += sqrt((pu*p->u[i])*(pu*p->u[i])+(pv*p->v[i])*(pv*p->v[i]));
				
				//printf("pt %f %f true %f %f\n",pu,pv,p->u[i],p->v[i]);
				
				weighted_sumdist += p->e[i]*sqrt((pu*p->u[i])*(pu*p->u[i])+(pv*p->v[i])*(pv*p->v[i]));
				sume+=p->e[i];	
			}	
			weighted_sumdist/=sume;
			sumdist/=n;		
			hist_sumdist[type]->Fill(sumdist,n,weight);
			hist_weighted_sumdist[type]->Fill(weighted_sumdist,n,weight);	
		
		
			double max_e=0;
			double max_e_dist=0;
			for(int i=0;i<p->entries;i++)
			{
				if(max_e<p->e[i])
				{
					max_e=p->e[i];
					max_e_dist=sqrt((p->u[i]-p->u[0])*(p->u[i]-p->u[0]) + (p->v[i]-p->v[0])*(p->v[i]-p->v[0]) + (p->z[i]-p->z[0])*(p->z[i]-p->z[0])   );
				
				
				}
			}
		
		
			hist_max_dist[type]->Fill(max_e_dist,weight);
			hist_max_dist_per_gev[type]->Fill(max_e_dist/sume*25.0,weight);
	
	

	
	
	
			///hist weighted area
			//////////
			
			//draw a line from the first to the last point and calculation the energy weighted r^2 from that line
			
			int ent=p->entries-1;
			double slopeu = (p->u[ent]-p->u[0])/(p->z[ent]-p->z[0]);
			double slopev = (p->v[ent]-p->v[0])/(p->z[ent]-p->z[0]);			
			
			double offu = p->u[0]-p->z[0]*slopeu;
			double offv = p->v[0]-p->z[0]*slopev;			
		
			double sum_dr=0;
			
			for(int i=0;i<p->entries;i++)
			{
				double eu=p->z[i]*slopeu+offu;
				double ev=p->z[i]*slopev+offv;
				
				double dr = sqrt((eu-p->u[i])*(eu-p->u[i]) + (ev-p->v[i])*(ev-p->v[i]));

				sum_dr+=dr*p->e[i];
			}
			
			sum_dr/=sume;
			hist_weighted_area[type]->Fill(sum_dr,weight);
			
			/////////
		
			int fp=-1;
			double fpe=0;
			int sp=-1;
			//second peak distance
			for(int i=1;i<p->entries;i++)
			{		
			
				if(i==1 &&p->e[0]>p->e[1])fp=0;
				
				if(i<p->entries-1)
				if(p->e[i-1]<p->e[i] && p->e[i+1]<p->e[i])//peak
				{
					if(fp<0){
						fp=i;
						fpe=p->e[i];
					}
					else if(p->e[i] *0.5 > fpe)
					{
						sp=i;
						break;
					}
				
				}
			
			
				if(i==p->entries-1 && p->e[i-1]<p->e[i])
				{
					if(fp>-1)sp=i;
					else fp=i;
				}
			
			}
			
			if(fp>-1)hist_first_peak[type]->Fill(p->z[fp]-p->z[0],weight);
				else hist_first_peak[type]->Fill(0.,weight);
			if(sp>-1)hist_second_peak[type]->Fill(p->z[sp]-p->z[0],weight);
				else hist_second_peak[type]->Fill(0.,weight);
			if(fp>-1 && sp>-1)hist_twopeak_diff[type]->Fill(p->z[sp]-p->z[fp],weight);	
				else hist_twopeak_diff[type]->Fill(0.,weight);



		/////////////sizemax diff
			double maxe=0;
			for(int i=0;i<p->entries;i++)
			{		
				if(maxe<p->e[i])maxe=p->e[i];
			}
			
			int pkidx=-1;
			for(int i=p->entries-1;i>-1;i--)
			{		
				if(maxe*0.7<p->e[i])
				{
					pkidx=i;
					break;
				}
			}			
		
	//		if(pkidx>-1 && p->z[pkidx]-p->z[0]>0.3)continue;
		
			if(pkidx>-1)
				hist_back_peak_dist[type]->Fill(p->z[pkidx]-p->z[0],weight);
			else
				hist_back_peak_dist[type]->Fill(0.,weight);
		/////////////


			//view e asym
			
			double eu=0;
			double ev=0;
			for(int i=0;i<p->entries;i++)
			{	
				if(p->view[i]==2)eu+=p->e[i];
				if(p->view[i]==3)ev+=p->e[i];
			}
			
			hist_view_e_asym[type]->Fill(fabs(eu-ev)/(eu+ev),weight);


		}



			typecnt[type]+=weight;



		if(ent_idx%10000==0)printf("Entry %d\n",ent_idx);
	}


	printf("used %d entries\n",saved);



printf("values\n");
for(int i=0;i<5;i++)
	printf("%d   %2.2f\n",i,typecnt[i]);
	
printf("FOM %2.2f\n",typecnt[1]/sqrt(typecnt[0]+typecnt[2]+typecnt[3]+typecnt[4]));


	TFile *f = new TFile("detailplots.root","RECREATE");
	f->cd();
	for(int i=0;i<5;i++)
	{
		hist_sumdist[i]->Write();
		hist_weighted_sumdist[i]->Write();
	
		hist_max_dist[i]->Write();
		hist_max_dist_per_gev[i]->Write(); 	
		
		hist_weighted_area[i]->Write();

		hist_first_peak[i]->Write();
		hist_second_peak[i]->Write();
		hist_twopeak_diff[i]->Write();
		hist_back_peak_dist[i]->Write();
		hist_view_e_asym[i]->Write();

		hist_prim_angle_vs_prim_eres[i]->Write();		
		hist_prim_angle_vs_prim_true[i]->Write();	


		hist_sum_curv_vs_part_e[i]->Write();	
		hist_pu_pz[i]->Write();	
		hist_pv_pz[i]->Write();	
		hist_pt_pz[i]->Write();	
		hist_pu_pv[i]->Write();

	}
	f->Close();


}



/*

TH3F * d_e_x_elec = new TH3F("d_e_x_elec","d_e_x_elec",100,0,1000,100,0,1000,100,0,10);
	TH3F * d_e_x_muon = new TH3F("d_e_x_muon","d_e_x_muon",100,0,1000,100,0,1000,100,0,10);
	TH3F * d_e_x_prot = new TH3F("d_e_x_prot","d_e_x_prot",100,0,1000,100,0,1000,100,0,10);
	TH3F * d_e_x_neut = new TH3F("d_e_x_neut","d_e_x_neut",100,0,1000,100,0,1000,100,0,10);		
	TH3F * d_e_x_other = new TH3F("d_e_x_other","d_e_x_other",100,0,1000,100,0,1000,100,0,10);
	
	TH2F * par_b_elec = new TH2F("par_b_elec","par_b_elec",100,0,10,100,0,100);
	TH2F * par_b_other = new TH2F("par_b_other","par_b_other",100,0,10,100,0,100);
	TH2F * par_b_else = new TH2F("par_b_else","par_b_else",100,0,10,100,0,100);
	
	
	TH2F * par_b_elec_rms = new TH2F("par_b_elec_rms","par_b_elec_rms",100,0,10,1000,0,0.02);
	TH2F * par_b_other_rms = new TH2F("par_b_other_rms","par_b_other_rms",100,0,10,1000,0,0.02);

	TH2F * tmax_vs_e = new TH2F("tmax_vs_e","tmax_vs_e",100,0,10,1000,0,10);
	
	TH2F * fitvscale = new TH2F("fitvscale","fitvscale",100,0,10000,100,0,10000);
	
	TH1F * fitvscaleres = new TH1F("fitvscaleres","fitvscaleres",100,-0.2,0.2);
	
	//MINE
	TH2F * neut_tail_z = new TH2F("neut_tail_z", "neut_tail_z", 50, -1.5,0.5,100,0,200);
	TH3F * neut_tail_uv = new TH3F("neut_tail_uv", "neut_tail_uv",50, -3,2,50,-1,1,100,0,200);
	int event_n_elec, event_n_neut;
	double n_tail_z, n_tail_u, n_tail_v, nvtxZ, evtxZ, nvtxU, evtxU, nvtxV, evtxV;
	//end MINE
	
	int ent =input->GetEntries();
	printf("%d total entries\n",ent);
	for(int i=0;i<ent;i++)
	{
		input->GetEntry(i);
		//printf("snarl %d\n",poh->GetHeader().GetSnarl());
		
		
		int type = (poh->mctrue.iaction>0)*(abs(poh->mctrue.inu)/2-5);
		
		if(truetype>-1)if(truetype!=type)continue;
		
//		if(poh->mctrue.iaction!=0)continue;
//		if(poh->mctrue.iresonance!=1001)continue;
//		if(poh->mctrue.inu!=12)continue;
	
//		if(poh->particles3d1->GetEntries()>1)continue;
		
		Particle3D * maxParticle=0;

		//MINE
		event_n_elec = 0;
		event_n_neut = 0;
		//end MINE
		
				
		for(int j=0;j<poh->particles3d.size();j++)
		{

			Particle3D * p = &poh->particles3d[j];
			if(!p)continue;
			
			if(maxParticle==0 || maxParticle->sum_e<p->sum_e)maxParticle=p;
			
			double dex=p->sum_e / (p->end_z-p->start_z);
			//double dex=p->sum_e / (double)(p->entries);
			
			switch(p->particletype)
			{
				case Particle3D::neutron:
					d_e_x_neut->Fill(dex,p->sum_e,(p->end_z-p->start_z));
					
					//MINE
					event_n_neut++;
					nvtxZ = p->start_z;
					nvtxU = p->start_u;
					nvtxV = p->start_v;
					//end MINE					
					break;
				case Particle3D::proton:
					d_e_x_prot->Fill(dex,p->sum_e,(p->end_z-p->start_z));
					break;
				case Particle3D::muon:
					d_e_x_muon->Fill(dex,p->sum_e,(p->end_z-p->start_z));
					break;
				case Particle3D::electron:
					d_e_x_elec->Fill(dex,p->sum_e,(p->end_z-p->start_z));
					if(p->emfit_b>0)
					{
						par_b_elec->Fill(p->emfit_b,p->sum_e);	
						par_b_elec_rms->Fill(p->emfit_b,p->avg_rms_t);
					}
					
					//MINE
					event_n_elec++;
					evtxZ = p->end_z;
					evtxU = p->end_u;
					evtxV = p->end_v;
					//end MINE
					break;
				case Particle3D::other:
					d_e_x_other->Fill(dex,p->sum_e,(p->end_z-p->start_z));
					if(p->emfit_b>0)
					{
						par_b_other->Fill(p->emfit_b,p->sum_e);
						par_b_other_rms->Fill(p->emfit_b,p->avg_rms_t);	
					}
					break;
					
				//else
					if(p->emfit_b>0)par_b_else->Fill(p->emfit_b,p->sum_e);	
					
					
			
			}
			
			if(p->emfit_e0>0)
			{
				fitvscale->Fill(p->emfit_e0,p->calibrated_energy);
				fitvscaleres->Fill((p->emfit_e0-p->calibrated_energy)/p->calibrated_energy);
			}
			
			//MINE
			if(event_n_neut ==1 && event_n_elec ==1){
					n_tail_z = nvtxZ - evtxZ;
					n_tail_u = nvtxU - evtxZ;
					n_tail_v = nvtxV - evtxV;
					neut_tail_z->Fill(n_tail_z,p->sum_e);
					neut_tail_uv->Fill(n_tail_u, n_tail_v,p->sum_e);
			}
			//end MINE
		}

		if(maxParticle!=0)	if(maxParticle->emfit_b>0 && maxParticle->sum_e>0)tmax_vs_e->Fill((maxParticle->emfit_a-1)/maxParticle->emfit_b,TMath::Log(maxParticle->sum_e));//,poh->event.oscprob);


		if(i%10000==0)printf("Entry %d\n",i);
	}

	TFile *f = new TFile("detailplots.root","RECREATE");
	f->cd();
	d_e_x_elec->Write();
	d_e_x_muon->Write();
	d_e_x_prot->Write();
	d_e_x_neut->Write();
	d_e_x_other->Write();
	par_b_other->Write();
	par_b_elec->Write();
	par_b_else->Write();
	par_b_other_rms->Write();
	par_b_elec_rms->Write();
	tmax_vs_e->Write();
	//MINE
	neut_tail_z->Write();
	neut_tail_uv->Write();
	fitvscale->Write();
	//end MINE
	f->Close();

*/



