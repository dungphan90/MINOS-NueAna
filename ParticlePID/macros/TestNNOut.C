#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"


void TestNNOut(string fname)
{
	TChain *c = new TChain("TrainTree");
	int added = c->Add(fname.c_str());
	if(!added)
	{
		printf("unable to add files %s\n",fname.c_str());
		exit(1);
	}
	
	int fNueClass;
	float oscprob;
	double totbeamweight;
	double pid;
	double weight;

	int resonance;
	double recoe;
	float truee;
	
	c->SetMakeClass(1);
	c->SetBranchAddress("mctrue_fNueClass",&fNueClass);
	c->SetBranchAddress("mctrue_oscprob",&oscprob);
	c->SetBranchAddress("mctrue_totbeamweight",&totbeamweight);
	c->SetBranchAddress("pid",&pid);
	c->SetBranchAddress("weight",&weight);
	c->SetBranchAddress("resonanceCode",&resonance);
	c->SetBranchAddress("mctrue_nuEnergy",&truee);
	c->SetBranchAddress("visenergy",&recoe);
	
	c->SetBranchStatus("*",0);
	c->SetBranchStatus("mctrue_fNueClass",1);
	c->SetBranchStatus("mctrue_oscprob",1);
	c->SetBranchStatus("mctrue_totbeamweight",1);
	c->SetBranchStatus("pid",1);
	c->SetBranchStatus("weight",1);
	c->SetBranchStatus("resonanceCode",1);
	c->SetBranchStatus("mctrue_nuEnergy",1);
	c->SetBranchStatus("visenergy",1);
	
	int ent = c->GetEntries();
	

	
	TFile * fout = TFile::Open("hists.root","RECREATE");
	TH1D * pidtype[5];
	TH3D * ehist[5];
	TH3D * rprh[5];
	TH3D * rpth[5];
	for(int i=0;i<5;i++)
	{
		char tmp[200];
		sprintf(tmp,"pid_type_%d",i);
		pidtype[i]=new TH1D(tmp,tmp,140,-0.2,1.2);
		sprintf(tmp,"e_res_true_reco_%d",i);
		ehist[i]=new TH3D(tmp,tmp,5,1,6,1000,0,20,1000,0,20);	
		sprintf(tmp,"e_res_pid_reco_%d",i);
		rprh[i]=new TH3D(tmp,tmp,5,1,6,1000,-0.2,1.2,1000,0,20);
		sprintf(tmp,"e_res_pid_true_%d",i);
		rpth[i]=new TH3D(tmp,tmp,5,1,6,1000,-0.2,1.2,1000,0,20);
	}
	
	TH1D * hpid = new TH1D("pid","pid",140,-0.2,1.2);

	
	printf("%d entries\n",ent);
	
	for(int i=0;i<ent;i++)
	{
		c->GetEntry(i);
		
		double w = weight*oscprob*totbeamweight;
		if(fNueClass<0||fNueClass>4)continue;
		pidtype[fNueClass]->Fill(pid,w);
		hpid->Fill(pid,w);

		//calibrate
		recoe=recoe*0.0387301+ 0.489987;
             
            
        


		resonance=resonance-1000;	
		ehist[fNueClass]->Fill(resonance,truee,recoe,w);
		rprh[fNueClass]->Fill(resonance,pid,recoe,w);
		rpth[fNueClass]->Fill(resonance,pid,truee,w);
	
		
	}
	

	
	TH1D * pidbg = (TH1D*)pidtype[0]->Clone("pidbg");
	pidbg->Add(pidtype[1]);
	pidbg->Add(pidtype[3]);
	pidbg->Add(pidtype[4]);

	TH1D * pidsig = (TH1D*)pidtype[2]->Clone("pidsig");	

	TH1D * fom = new TH1D("fom","fom",140,-0.2,1.2);
	TH1D * superfom = new TH1D("superfom","superfom",140,-0.2,1.2);	
	
	for(int i=1;i<141;i++)
	{
		double s=pidsig->Integral(i,141);
		double b=pidbg->Integral(i,141);
		
		double f = b ? s/sqrt(b):0;
		fom->SetBinContent(i,f);
		
		double sf = b ? s/sqrt(b+b*b*0.01) : 0;
		superfom->SetBinContent(i,sf);
	
	}

	for(int i=0;i<5;i++)
	{
		pidtype[i]->Write();
		ehist[i]->Write();
		rprh[i]->Write();
		rpth[i]->Write();
	}


	fom->Write();
	superfom->Write();
	pidbg->Write();
	pidsig->Write();
	hpid->Write();
	
	
	printf("MaxFom %f at %f\n",fom->GetMaximum(),fom->GetBinLowEdge(fom->GetMaximumBin()));


		
	printf("MaxSuperFom %f at %f\n",superfom->GetMaximum(),superfom->GetBinLowEdge(superfom->GetMaximumBin()));

	int pidcut=superfom->GetMaximumBin();

	printf("at maxfom\n");
	printf("Signal %f\n",pidtype[2]->Integral(pidcut,141));
        printf("BG %f\n",pidbg->Integral(pidcut,141));
	printf("NC %f\n",pidtype[0]->Integral(pidcut,141));
        printf("CC %f\n",pidtype[1]->Integral(pidcut,141));
        printf("Tau %f\n",pidtype[3]->Integral(pidcut,141));
        printf("Beam %f\n",pidtype[4]->Integral(pidcut,141));


	fout->Close();

	
}
