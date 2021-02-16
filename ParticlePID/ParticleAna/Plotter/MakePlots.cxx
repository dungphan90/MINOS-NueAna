#include "NueAna/ParticlePID/ParticleAna/Plotter/MakePlots.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TLegend.h"
#include "TList.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TLine.h"
#include "TFile.h"

MakePlots::MakePlots()
{
	nchains=0;
	sprintf(outdir,"plots");
	sprintf(outfile,"plots.root");
}

void MakePlots::AddSet(const char* filename, const char* tn, int mc, int mrcc)
{
	nchains=chains.size();
	chains.push_back(new TChain("PA"));
	int added = chains[nchains]->Add(filename);
	if(!added)
	{
		printf("error adding chain %s '%s' mc %d\n",filename,tn, mc);
		delete chains[nchains];
		return;
	}
	printf("adding chain %s '%s' mc %d\n",filename,tn, mc);
	char tmpchar[200];
	sprintf(tmpchar,"%s",tn);
	typesname.push_back(tmpchar);
	ismc.push_back(mc);
	ismrcc.push_back(mrcc);
	
	TChain * POTchain=new TChain("pottree");
	POTchain->Add(filename);
	int ent=POTchain->GetEntries();
	pots.push_back(0);
	for(int i=0;i<ent;i++)
	{
		POTchain->GetEntry(i);
		pots[nchains]+=POTchain->GetLeaf("pot")->GetValue();
	}
	printf("%s has %f POTs \n",tn,pots[nchains]);


	nchains++;
	std::vector<std::vector<TH1F *> > a;
	histos.push_back(a);
}

void MakePlots::Run()
{

//	AddSet("/minos/data/dogwood0/near/data/PO*.root","near data",0);
//	AddSet("/minos/data/dogwood0/near/mc/PO*.root","near mc",1);
	

/*	TChain * near_mc=new TChain("PA");
	TChain * near_data=new TChain("PA");
	
	TChain * near_mc_POT=new TChain("pottree");
	TChain * near_data_POT=new TChain("pottree");
	
	/ *near_data->Add("/minos/data/dogwood0/near/data/PO*.root");
	near_data_POT->Add("/minos/data/dogwood0/near/data/PO*.root");

	near_mc->Add("/minos/data/dogwood0/near/mc/PO*.root");
	near_mc_POT->Add("/minos/data/dogwood0/near/mc/PO*.root");
	*/
	
/*	near_data->Add("/data/minos/ManagedPIDCVS/PO-N00011449_0000.spill.sntp.dogwood0.0-Trim-NN.root");
	near_data_POT->Add("/data/minos/ManagedPIDCVS/PO-N00011449_0000.spill.sntp.dogwood0.0-Trim-NN.root");

	near_mc->Add("/data/minos/ManagedPIDCVS/PO-n13037040_0000_L010185N_D04.sntp.dogwood0-Trim-NN.root");
	near_mc_POT->Add("/data/minos/ManagedPIDCVS/PO-n13037040_0000_L010185N_D04.sntp.dogwood0-Trim-NN.root");
	
	
	
	for(int i=0;i<nType;i++)pots[i]=0;
	
	int ent;
	
	ent=near_data_POT->GetEntries();
	for(int i=0;i<ent;i++)
	{
		near_data_POT->GetEntry(i);
		pots[0]+=near_data_POT->GetLeaf("pot")->GetValue();
	}
	printf("%f near data POTs\n",pots[0]);
	ismc[0]=0;

	ent=near_mc_POT->GetEntries();
	for(int i=0;i<ent;i++)
	{
		near_mc_POT->GetEntry(i);
		pots[1]+=near_mc_POT->GetLeaf("pot")->GetValue();
	}
	printf("%f near mc POTs\n",pots[1]);
	ismc[1]=1;
*/	

	if(!nchains)
	{
		printf("You must specify some input sets...\n");
		return;
	}
	
	printf("outdir: %s\n",outdir);
	
	
	//define cut names
	sprintf(cutnames[0],"No Cuts");
	sprintf(cutnames[1],"Fiducial");
	sprintf(cutnames[2],"Preselection");
	sprintf(cutnames[3],"PID");


	SetupHistos();
	
	for(int i=0;i<nchains;i++)
		FillHistos(i);


}

void MakePlots::SetOutDir(const char*o)
{

	sprintf(outdir,"%s",o);
}


void MakePlots::SetOutFile(const char*o)
{

	sprintf(outfile,"%s",o);
}


void MakePlots::SetupHistos()
{
	char tmp[200];
	
	for(int i=0;i<(int)histos.size();i++)	
	for(int j=0;j<nCut;j++)	
	{
		std::vector<TH1F*> a;
		histos[i].push_back(a);
		
		sprintf(tmp,"particles_longest_z_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,10,0,2));
		sprintf(names[0],"particles_longest_z");	

		sprintf(tmp,"particles_longest_s_particle_s_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,2));
		sprintf(names[1],"particles_longest_s_particle_s");			

		sprintf(tmp,"particles_primary_long_e_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,10,0,2));
		sprintf(names[2],"particles_primary_long_e");
		
		sprintf(tmp,"particles_elec_vise_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,10));
		sprintf(names[3],"particles_elec_vise");				

		sprintf(tmp,"particles_primary_phi_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,1.6));
		sprintf(names[4],"particles_primary_phi");	
			
		sprintf(tmp,"particles_prim_par_b_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,1.5));		
		sprintf(names[5],"particles_prim_par_b");	
			
		sprintf(tmp,"particles_mol_rad_r_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,0.2));
		sprintf(names[6],"particles_mol_rad_r");	
			
		sprintf(tmp,"particles_emfrac_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,1.1));	
		sprintf(names[7],"particles_emfrac");	
		
		sprintf(tmp,"particles_ntot_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,20));
		sprintf(names[8],"particles_ntot");	
			
		sprintf(tmp,"particles_total_long_e_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,200));	
		sprintf(names[9],"particles_total_long_e");	
		
		sprintf(tmp,"particles_weighted_phi_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,1.6));	
		sprintf(names[10],"particles_weighted_phi");	
		
		sprintf(tmp,"largest_frac_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,50,0,2));	
		sprintf(names[11],"largest_frac");	
		
		sprintf(tmp,"reco_frac_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,50,0,2));	
		sprintf(names[12],"reco_frac");	

		sprintf(tmp,"particles_totvise_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,20));	
		sprintf(names[13],"particles_totvise");	
		
		sprintf(tmp,"event_visenergy_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,20));	
		sprintf(names[14],"event_visenergy");			


		sprintf(tmp,"particles_rms_r_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,4));	
		sprintf(names[15],"particles_rms_r");	
		
		sprintf(tmp,"particles_prim_par_a_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,10));	
		sprintf(names[16],"particles_prim_par_a");	


		sprintf(tmp,"particles_prim_par_e0_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,10));	
		sprintf(names[17],"particles_prim_par_e0");	

		sprintf(tmp,"particles_prim_par_chisq_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0.1,500));	
		sprintf(names[18],"particles_prim_par_chisq");	

		sprintf(tmp,"particles_largest_particle_peakdiff_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,-5,5));	
		sprintf(names[19],"particles_largest_particle_peakdiff");									

		sprintf(tmp,"largest_cmp_chisqndf_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,2));	
		sprintf(names[20],"largest_cmp_chisqndf");	
		
		sprintf(tmp,"event_nstrips_%d_%d",i,j);
		histos[i][j].push_back(new TH1F(tmp,tmp,20,0,80));	
		sprintf(names[21],"event_nstrips");	

	}

	for(int j=0;j<nCut;j++)	
		for(int i=0;i<(int)histos.size();i++)
			for(int k=0;k<nHistos;k++)
				histos[i][j][k]->Sumw2();
}

void MakePlots::SetBranches(TChain*c)
{
	c->SetMakeClass(1);
	c->SetBranchStatus("*",0);
	c->SetBranchStatus("particles.primary_long_e",1);
	c->SetBranchStatus("particles.longest_s_particle_s",1);
	c->SetBranchStatus("particles.elec_vise",1);
	c->SetBranchStatus("particles.primary_phi",1);
	c->SetBranchStatus("particles.prim_par_b",1);
	c->SetBranchStatus("particles.mol_rad_r",1);
	c->SetBranchStatus("particles.longest_z",1);
	c->SetBranchStatus("particles.emfrac",1);
	c->SetBranchStatus("particles.ntot",1);
	c->SetBranchStatus("particles.total_long_e",1);
	c->SetBranchStatus("particles.weighted_phi",1);
	c->SetBranchStatus("particles.largest_particle_e",1);
	c->SetBranchStatus("particles.totvise",1);
	c->SetBranchStatus("particles.prim_vise",1);
	c->SetBranchStatus("event.visenergy",1);
	c->SetBranchStatus("event.nstrips",1);
	c->SetBranchStatus("event.inFiducial",1);
	c->SetBranchStatus("event.max_z",1);
	c->SetBranchStatus("event.min_z",1);
	c->SetBranchStatus("mctrue.type",1);
	c->SetBranchStatus("mctrue.oscprob",1);
	c->SetBranchStatus("mctrue.totbeamweight",1);
	c->SetBranchStatus("mctrue.iresonance",1);
	c->SetBranchStatus("event.contained",1);			
	c->SetBranchStatus("event.pidA",1);		
	c->SetBranchStatus("event.pidB",1);		
	c->SetBranchStatus("event.pidC",1);		
	c->SetBranchStatus("event.pidD",1);
	c->SetBranchStatus("mrccinfo.particle_s",1);
	c->SetBranchStatus("mrccinfo.sum_e",1);

	c->SetBranchStatus("particles.rms_r",1);
	c->SetBranchStatus("particles.prim_par_a",1);
	c->SetBranchStatus("particles.prim_par_e0",1);
	c->SetBranchStatus("particles.prim_par_chisq",1);
	c->SetBranchStatus("particles.largest_particle_peakdiff",1);
	c->SetBranchStatus("particles.largest_particle_cmp_chisq",1);
	c->SetBranchStatus("particles.largest_particle_cmp_ndf",1);	
	c->SetBranchStatus("particles.longest_particle_type",1);

	
	c->SetBranchAddress("particles.primary_long_e",&particles_primary_long_e);
	c->SetBranchAddress("particles.longest_s_particle_s",&particles_longest_s_particle_s);
	c->SetBranchAddress("particles.elec_vise",&particles_elec_vise);
	c->SetBranchAddress("particles.primary_phi",&particles_primary_phi);
	c->SetBranchAddress("particles.prim_par_b",&particles_prim_par_b);
	c->SetBranchAddress("particles.mol_rad_r",&particles_mol_rad_r);
	c->SetBranchAddress("particles.longest_z",&particles_longest_z);
	c->SetBranchAddress("particles.emfrac",&particles_emfrac);
	c->SetBranchAddress("particles.ntot",&particles_ntot);
	c->SetBranchAddress("particles.total_long_e",&particles_total_long_e);
	c->SetBranchAddress("particles.weighted_phi",&particles_weighted_phi);
	c->SetBranchAddress("particles.largest_particle_e",&particles_largest_particle_e);
	c->SetBranchAddress("particles.totvise",&particles_totvise);
	c->SetBranchAddress("particles.prim_vise",&particles_prim_vise);
	c->SetBranchAddress("event.visenergy",&event_visenergy);
	c->SetBranchAddress("event.nstrips",&event_nstrips);
	c->SetBranchAddress("event.inFiducial",&event_inFiducial);
	c->SetBranchAddress("event.max_z",&event_max_z);
	c->SetBranchAddress("event.min_z",&event_min_z);
	c->SetBranchAddress("mctrue.type",&mctrue_type);
	c->SetBranchAddress("mctrue.oscprob",&mctrue_oscprob);
	c->SetBranchAddress("mctrue.totbeamweight",&mctrue_totbeamweight);
	c->SetBranchAddress("mctrue.iresonance",&mctrue_iresonance);
	c->SetBranchAddress("event.contained",&event_contained);
	c->SetBranchAddress("event.pidA",&event_pidA);
	c->SetBranchAddress("event.pidB",&event_pidB);
	c->SetBranchAddress("event.pidC",&event_pidC);
	c->SetBranchAddress("event.pidD",&event_pidD);
	c->SetBranchAddress("mrccinfo.particle_s",&mrccinfo_particle_s);
	c->SetBranchAddress("mrccinfo.sum_e",&mrccinfo_sum_e);
	
	c->SetBranchAddress("particles.rms_r",&particles_rms_r);
	c->SetBranchAddress("particles.prim_par_a",&particles_prim_par_a);
	c->SetBranchAddress("particles.prim_par_e0",&particles_prim_par_e0);
	c->SetBranchAddress("particles.prim_par_chisq",&particles_prim_par_chisq);
	c->SetBranchAddress("particles.largest_particle_peakdiff",&particles_largest_particle_peakdiff);
	c->SetBranchAddress("particles.largest_particle_cmp_chisq",&particles_largest_particle_cmp_chisq);
	c->SetBranchAddress("particles.largest_particle_cmp_ndf",&particles_largest_particlecmp_ndf);	
	c->SetBranchAddress("particles.longest_particle_type",&particles_longest_particle_type);
	

}

void MakePlots::SaveRatioPlots()
{
	if(ratiopair.size()<1)return;



	const int nHist=ratiopair.size();
	std::vector<TH1F*>  histos[nCut][nHistos];
	
	
	for(int j=0;j<nCut;j++)
	for(int k=0;k<nHistos;k++)
	for(int i=0;i<nHist;i++)
	{
		histos[j][k].push_back((TH1F*)0);
		
		printf("making %s with %s and %s %d %d\n",rationame[i].c_str(),this->histos[ratiopair[i].first][j][k]->GetName(),this->histos[ratiopair[i].second][j][k]->GetName(),j,k);
		histos[j][k][i]=(TH1F*)this->histos[ratiopair[i].first][j][k]->Clone(rationame[i].c_str());
		histos[j][k][i]->SetDirectory(0);
		histos[j][k][i]->Divide((TH1F*)this->histos[ratiopair[i].second][j][k]);
		printf("made %s max %f\n",histos[j][k][i]->GetName(),histos[j][k][i]->GetMaximum());
	}

	TFile *of = TFile::Open(outfile,"RECREATE");
	of->cd();
	for(int k=0;k<nHistos;k++)
	for(int j=0;j<nCut;j++)
	for(int i=0;i<nHist;i++)
		histos[j][k][i]->Write();
	
	of->Close();


	TCanvas c1;
	TPad p1("p1","",0,0.4,1,1);
	TPad p2("p2","",0,0,1,0.4);
	c1.Draw();
	p1.Draw();
	p2.Draw();
	for(int k=0;k<nHistos;k++)
	for(int j=0;j<nCut;j++)
	{
			p1.cd();

			double maxb=0;
			for(int i=0;i<nHist;i++)
			{
				//char tmpchar[200];
				double maxt=histos[j][k][i]->GetMaximum();
				if(maxt>maxb)maxb=maxt;
			}
			
			for(int i=0;i<nHist;i++)
			{
		//		histos[j][k][i]->GetYaxis()->SetRangeUser(0.5,1.5);//0,maxb*1.2);
			
				if(i==0)
				{
					histos[j][k][i]->Draw("E");
					continue;
				}
				histos[j][k][i]->SetLineColor(1+i);
				histos[j][k][i]->Draw("SAMESE");

			}
			
	
			p2.cd();
			
			double extt=0;
			std::vector<TH1F*> tmp;
			
			if(nHist>1)
			{
			for(int i=1;i<nHist;i++)
			{
				tmp.push_back((TH1F*)0);
				tmp[i-1]= (TH1F*)histos[j][k][0]->Clone("tmp");
				tmp[i-1]->Add(histos[j][k][i],-1);	
				tmp[i-1]->Divide(histos[j][k][i]);
				tmp[i-1]->SetLineColor(i+1);
				
				int minbin=tmp[i-1]->GetMinimumBin();
				int maxbin=tmp[i-1]->GetMaximumBin();
				
				double min=tmp[i-1]->GetBinContent(minbin);
				double max=tmp[i-1]->GetBinContent(maxbin);
			
				double ext = TMath::Abs(min);
				ext = TMath::Abs(max) < ext?ext:TMath::Abs(max);
				if(ext>extt)extt=ext;
			}
			extt=extt>0.2?0.2:extt;
			extt*=1.25;
			
			for(int i=0;i<nHist-1;i++)
			{
				tmp[i]->GetYaxis()->SetRangeUser(-extt,extt);
			
				tmp[i]->SetStats(0);
				tmp[i]->SetTitle("");
				if(i==0)
				{
					tmp[i]->Draw();
					continue;
				}
				tmp[i]->Draw("same");
			}
			
			TLine l(tmp[0]->GetXaxis()->GetXmin(),0,tmp[0]->GetXaxis()->GetXmax(),0);
			l.Draw("same");
			}
			
			SortOutStats(&c1);
			
			char tmpchar[200];
			printf("%d  %s/%s_%d.gif\n",k,outdir,names[k],j);
			sprintf(tmpchar,"%s/%s_%d.gif",outdir,names[k],j);
			c1.SaveAs(tmpchar);
			for(int i=0;i<(int)tmp.size()-1;i++){delete tmp[i];tmp[i]=0;}
	}
}



void MakePlots::SavePlots(int norm)
{
	TFile *of = TFile::Open(outfile,"RECREATE");
	of->cd();


	
	
	for(int k=0;k<nHistos;k++)
	for(int j=0;j<nCut;j++)
	for(int i=0;i<(int)histos.size();i++)
	{
		if(norm && histos[i][j][k]->Integral())histos[i][j][k]->Scale(1./histos[i][j][k]->Integral());
		histos[i][j][k]->Write();
	}
	of->Close();


	TCanvas c1;
	TPad p1("p1","",0,0.4,1,1);
	TPad p2("p2","",0,0,1,0.4);
	c1.Draw();
	p1.Draw();
	p2.Draw();
	for(int k=0;k<nHistos;k++)
	for(int j=0;j<nCut;j++)
	{
			p1.cd();

			double maxb=0;
			for(int i=0;i<(int)histos.size();i++)
			{
				char tmpchar[200];
				sprintf(tmpchar,"%s  %s",names[k],cutnames[j]);
				histos[i][j][k]->SetTitle(tmpchar);
				double maxt=histos[i][j][k]->GetMaximum();
				histos[i][j][k]->SetXTitle(names[k]);
				if(maxt>maxb)maxb=maxt;
			}
			
			
			
			for(int i=0;i<(int)histos.size();i++)
			{
				histos[i][j][k]->GetYaxis()->SetRangeUser(0,maxb*1.2);
				histos[i][j][k]->SetName(typesname[i].c_str());
				if(i==0)
				{
					histos[i][j][k]->Draw("E");
					continue;
				}
				histos[i][j][k]->SetLineColor(1+i);
				histos[i][j][k]->Draw("SAMESE");

			}
			
			
			
			

			

	
			p2.cd();
			
			double extt=0;
			std::vector<TH1F *> tmp;
			for(int i=1;i<(int)histos.size();i++)
			{
				tmp.push_back((TH1F*)histos[0][j][k]->Clone("tmp"));
				tmp[i-1]->Add(histos[i][j][k],-1);	
				tmp[i-1]->Divide(histos[i][j][k]);
				tmp[i-1]->SetLineColor(i+1);
				
				int minbin=tmp[i-1]->GetMinimumBin();
				int maxbin=tmp[i-1]->GetMaximumBin();
				
				double min=tmp[i-1]->GetBinContent(minbin);
				double max=tmp[i-1]->GetBinContent(maxbin);
			
				double ext = TMath::Abs(min);
				ext = TMath::Abs(max) < ext?ext:TMath::Abs(max);
				if(ext>extt)extt=ext;
			}
			extt=extt>0.2?0.2:extt;
			extt*=1.25;
			
			for(int i=0;i<(int)histos.size()-1;i++)
			{
				tmp[i]->GetYaxis()->SetRangeUser(-extt,extt);
			
				tmp[i]->SetStats(0);
				tmp[i]->SetTitle("");
				if(i==0)
				{
					tmp[i]->Draw();
					continue;
				}
				tmp[i]->Draw("same");
			}
			
			TLine l(tmp[0]->GetXaxis()->GetXmin(),0,tmp[0]->GetXaxis()->GetXmax(),0);
			l.Draw("same");
			
			
			SortOutStats(&c1);
			
			char tmpchar[200];
			printf("%d  %s/%s_%d.gif\n",k,outdir,names[k],j);
			sprintf(tmpchar,"%s/%s_%d.gif",outdir,names[k],j);
			c1.SaveAs(tmpchar);
			for(int i=0;i<(int)histos.size()-1;i++)delete tmp[i];
	}
	
	

	
}



void MakePlots::FillHistos(int set)
{
	TChain *c=chains[set];

	printf("make plots set %d\n",set);

	SetBranches(c);

	if(!pots[set])return;
	int ent=c->GetEntries();
	printf("%d entries\n",ent);
	for(int i=0;i<ent;i++)
	{
	//	if(i>1000)break;//testing
	
		if(i%10000==0)printf("entry %d\n",i);
		c->GetEntry(i);
		
		//is this an mrcc file?  if so do the mrcc cut
		if(ismrcc[set])
		{
			if(mrccinfo_particle_s<3)continue;//no long muon track
			event_visenergy-= (mrccinfo_sum_e*0.03799819+0.4803261); //adjust!
		}
		
		
		
		FillHistos(set,0);
		
		
		if(event_inFiducial!=1 || event_contained!=1)continue;
		FillHistos(set,1);
		
		//presel
		if(particles_ntot<1)continue;
		if(particles_longest_s_particle_s>1.2 ||particles_longest_s_particle_s<0.2)continue;
		//broken in mrcc files if(event_max_z-event_min_z>0.9 || event_max_z-event_min_z < 0.2)continue;
		//need to use particle total energy rather than strip energy for mrcc...
		if(particles_totvise<0.2 || particles_totvise>8)continue;

		//if(event_nstrips<10 || event_nstrips>40)continue;

		//don't want a long muon..
		if(particles_longest_particle_type>12)continue;		


		if(particles_largest_particle_cmp_chisq>30)continue;



		FillHistos(set,2);

		if(event_pidB<0.5)continue;
	
	//	if(event_pidA<0.5)continue;
		FillHistos(set,3);

	}

	
	for(int j=0;j<nCut;j++)
	for(int k=0;k<nHistos;k++)
	{
		if(pots[set])histos[set][j][k]->Scale(3.25e8/pots[set]);
	}

}


void MakePlots::FillHistos(int type,int cut)
{
	if(!pots[type])return;
	
	double weight = 1;
	if(ismc[type]==1)weight*=mctrue_totbeamweight*mctrue_oscprob;
	
	histos[type][cut][0]->Fill(particles_longest_z,weight);
	histos[type][cut][1]->Fill(particles_longest_s_particle_s,weight);
	histos[type][cut][2]->Fill(particles_primary_long_e,weight);
	histos[type][cut][3]->Fill(particles_elec_vise/25.,weight);
	histos[type][cut][4]->Fill(particles_primary_phi,weight);
	if(particles_prim_par_e0)histos[type][cut][5]->Fill(particles_prim_par_b,weight);
	histos[type][cut][6]->Fill(particles_mol_rad_r,weight);
	histos[type][cut][7]->Fill(particles_emfrac,weight);
	histos[type][cut][8]->Fill(particles_ntot,weight);
	histos[type][cut][9]->Fill(particles_total_long_e,weight);
	histos[type][cut][10]->Fill(particles_weighted_phi,weight);
	histos[type][cut][11]->Fill((particles_largest_particle_e*0.03379042+0.5866994)/particles_totvise,weight);
	histos[type][cut][12]->Fill(particles_prim_vise/25./event_visenergy,weight);
	histos[type][cut][13]->Fill(particles_totvise,weight);
	histos[type][cut][14]->Fill(event_visenergy,weight);

	histos[type][cut][15]->Fill(particles_rms_r,weight);
	if(particles_prim_par_e0)histos[type][cut][16]->Fill(particles_prim_par_a,weight);
	if(particles_prim_par_e0)histos[type][cut][17]->Fill(particles_prim_par_e0/25./60.,weight);
	histos[type][cut][18]->Fill(particles_prim_par_chisq,weight);
	histos[type][cut][19]->Fill(particles_largest_particle_peakdiff/25.,weight);
	histos[type][cut][20]->Fill(particles_largest_particle_cmp_chisq/particles_largest_particlecmp_ndf,weight);
	histos[type][cut][21]->Fill(event_nstrips,weight);	
	
}



void MakePlots::SortOutStats(TPad *pad,float fracPadx,float fracPady,int numPerCol,
		  float startx,float starty){

  pad->Update();
  pad->Modified();

  TList *theList = pad->GetListOfPrimitives();
  TIterator *iter = theList->MakeIterator();
  TObject *ob;
  int numStats = 0;
  while((ob = iter->Next())){
    if(ob->InheritsFrom("TPad")) {
      TPad *p = (TPad*) ob;
      SortOutStats(p,fracPadx,fracPady,numPerCol,startx,starty);
    }
    else if(ob->InheritsFrom("TH1")||ob->InheritsFrom("TH2")
	    ||ob->InheritsFrom("TGraph")||ob->InheritsFrom("TProfile")){
      if(ob->FindObject("stats")) numStats+=1;
    }
    else if(ob->InheritsFrom("TLegend")){
      numStats+=1;
    }
  }
  
  if(numStats!=0){
    
    float xwidth = fracPadx/float(int((numStats-1)/numPerCol)+1);
    float ywidth = 0.2;
    if(numStats>=numPerCol) ywidth = fracPady/numPerCol;
    else ywidth = fracPady/float(numStats%numPerCol);
    float xgap = 0.01;
    float ygap = 0.01;
    float x2 = startx;
    float y2 = starty;
    float x1 = x2-xwidth;
    float y1 = y2-ywidth;
    int cnt = 0;
    
    iter->Reset();
    while((ob = iter->Next())){
      int col = 1;
      if(ob->InheritsFrom("TH1")){
	TH1 *h1 = (TH1*) ob;
	col = h1->GetLineColor();
      }
      else if(ob->InheritsFrom("TProfile")){
	TProfile *p1 = (TProfile*) ob;
	col = p1->GetLineColor();
      }
      else if(ob->InheritsFrom("TGraph")){
	TGraph *g1 = (TGraph*) ob;
	col = g1->GetLineColor();
      }

      TPaveStats *stat = 0;
      if((stat = (TPaveStats*) ob->FindObject("stats"))){
	stat->SetTextColor(col);
	stat->SetX1NDC(x1-int(cnt/numPerCol)*(xwidth+xgap));
	stat->SetX2NDC(x2-int(cnt/numPerCol)*(xwidth+xgap));
	stat->SetY1NDC(y1-float(cnt%numPerCol)*(ywidth+ygap));
	stat->SetY2NDC(y2-float(cnt%numPerCol)*(ywidth+ygap));
	cnt+=1;
      }
      else if(ob->InheritsFrom("TLegend")){
	TLegend *leg = (TLegend*) ob;
	leg->SetX1NDC(x1-int(cnt/numPerCol)*(xwidth+xgap));
	leg->SetX2NDC(x2-int(cnt/numPerCol)*(xwidth+xgap));
	leg->SetY1NDC(y1-float(cnt%numPerCol)*(ywidth+ygap));
	leg->SetY2NDC(y2-float(cnt%numPerCol)*(ywidth+ygap));
	cnt+=1;
      }
    }
  }
}


