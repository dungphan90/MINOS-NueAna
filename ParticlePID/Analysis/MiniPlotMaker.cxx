#include "NueAna/ParticlePID/Analysis/MiniPlotMaker.h"

#include <fstream>
#include <math.h>
#include "OscProb/OscCalc.h"
#include "NueAna/NueStandard.h"
#include "TTree.h"

 #include <exception>

//from greg's nue how to
//sept 1 2009
// in units of 1
double MiniPlotMaker::potratios[3]={1.20378e20,1.92724e20,3.78478e20}; 
double MiniPlotMaker::potratios_mrcc[3]={4.207220e+19,1.416530e+20,2.565563e+20};

int MiniPlotMaker::AddInputFiles(string s,int group)
{
	infileCount[group] += infiles[group]->Add(s.c_str()); 
	infilePOTTreeCount[group] += infilesPOTTree[group]->Add(s.c_str());
	infileForPOTCountingCount[group] += infilesForPOTCounting[group]->Add(s.c_str());
	 
	return infileCount[group];
}


int MiniPlotMaker::GetRunPeriod(time_t ts)
{
		
	if(ts<1145036420)return 0;
	if(ts<1190028111)return 1;
		
	return 2;
}

void MiniPlotMaker::CountPots(TChain *c, TChain *d, double rp[3])
{
	//set up chain
	double tpot=0;
	c->SetMakeClass(1);
	c->SetBranchStatus("*",0);
	c->SetBranchStatus("pot",1);
	c->SetBranchAddress("pot",&tpot);


	int nsnarls=0;
        c->SetBranchStatus("nsnarls",1);
	c->SetBranchAddress("nsnarls",&nsnarls);
	

	d->SetBranchStatus("*",1);
	d->SetBranchAddress("NueMiniPID",&nm);
	d->SetBranchStatus("*",0);
	d->SetBranchStatus("timestamp",1);
	

	rp[0]=0;
	rp[1]=0;
	rp[2]=0;
	
	int ent=c->GetEntries();
	
	Long64_t *offsets=d->GetTreeOffset();
	
	int snarls=0;
	
	if(!isMC || detector==0)
	{
		for(int i=0;i<ent;i++)
		{
			c->GetEntry(i);
			//check the timestamp of the first entry from this tree
			d->GetEntry(offsets[i]);
			
	
			snarls+=nsnarls;
			rp[GetRunPeriod(nm->timestamp)]+=tpot;
		}
		printf("estimating %d snarls...\n",snarls);
		return;
	}
	
	//if we are here... we are doing near mc
	//near mc has entries from all runs in a single file!
	
	printf("counting pots for near mc\n");
	
	

	d->SetBranchStatus("fDet",1);
	d->SetBranchStatus("fBeam",1);
	d->SetBranchStatus("fRelease",1);
	d->SetBranchStatus("snarl",1);
	d->SetBranchStatus("fPOT",1);
	d->SetBranchStatus("passes_NueStandard_PassesDataQuality",1);

	ent = d->GetEntries();
	
	int lastsnarl=-1;
	for(int i=0;i<ent;i++)
	{
		d->GetEntry(i);
		snarls++;
		//if(nm->snarl==lastsnarl)continue;
		
		if(!nm->passes_NueStandard_PassesDataQuality)continue;
		if(!nm->fPOT)continue;
		//a near mc will have a fPOT of -9999 when we should be counting it...
		
		if(isMC)
		{
			rp[GetRunPeriod(nm->timestamp)]+=
				MCInfo::GetMCPoT((Detector::Detector_t)nm->fDet,
					(BeamType::BeamType_t)nm->fBeam,
					(ReleaseType::Release_t)nm->fRelease);
		}
				
		lastsnarl=nm->snarl;
		
	}
	
	d->SetMakeClass(0);
	d->ResetBranchAddresses();
	d->SetBranchStatus("*",0);

        printf("estimating %d snarls...\n",snarls);

	
}




MiniPlotMaker::~MiniPlotMaker()
{
}

void MiniPlotMaker::MakePlots()
{

	//check the input counts to make sure we don't have bogus files...
	int fcheck = infileCount[0] == infilePOTTreeCount[0];
	fcheck = fcheck && infileCount[0] == infilePOTTreeCount[0];
	fcheck = fcheck && infileCount[0] == infilePOTTreeCount[0];

	if(!fcheck)
	{
		printf("NueMini/POTTree count mismatch!\n");
		return;
	}

	PrintConfigs();

	outfile=0;
	
	if(outputFile!="" && !overwrite)outfile=new TFile(outputFile.c_str(),"UPDATE");	
	if(outputFile!="" && overwrite)outfile=new TFile(outputFile.c_str(),"RECREATE");
	if(!outfile)
	{
		printf("error opening output file %s\n",outputFile.c_str());
		return;
	}

	outfile->SetCompressionLevel(9);

	TDirectory *workingDir = MakeDirectory(outfile);
	workingDir->cd();

	//make histos.....
	MakeHistos();



	//compute pots
	double weight[3][3];
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)weight[i][j]=0;


	fstream file_list("filelist",fstream::out | fstream::app);

	for(int i=0;i<3;i++)
	{
		CountPots(infilesPOTTree[i],
			infilesForPOTCounting[i],	
			weight[i]);
		
		weight[i][0]*=1e12;
		weight[i][1]*=1e12;
		weight[i][2]*=1e12;
		
		//save a record for the tree!
		for(int j=0;j<3;j++)
			inputPOT[i][j]=weight[i][j];	
		
		
		char tmp[500];
		printf("group %d has %e pots %e %e %e\n",i,(weight[i][0]+weight[i][1]+weight[i][2]),	
			weight[i][0],weight[i][1],weight[i][2]);

		sprintf(tmp,"group %d has %e pots %e %e %e\n",i,(weight[i][0]+weight[i][1]+weight[i][2]),	
			weight[i][0],weight[i][1],weight[i][2]);
			
		file_list<<tmp;	
			
		if(isMC)
		{	
		
			//in the near, we are normalizing each run to the same...
			if(detector==1)
			{
				if(weight[i][0])weight[i][0]=wantPOTs/weight[i][0];
				if(weight[i][1])weight[i][1]=wantPOTs/weight[i][1];
				if(weight[i][2])weight[i][2]=wantPOTs/weight[i][2];						
			}else{
		
				if(weight[i][0] && weight[i][1] && weight[i][2])
				{
					if(!isMRCC){
						weight[i][0]=1./weight[i][0]*potratios[0];
						weight[i][1]=1./weight[i][1]*potratios[1];
						weight[i][2]=1./weight[i][2]*potratios[2];
					}else{
                                                weight[i][0]=1./weight[i][0]*potratios_mrcc[0];
                                                weight[i][1]=1./weight[i][1]*potratios_mrcc[1];
                                                weight[i][2]=1./weight[i][2]*potratios_mrcc[2];
					}
				}
			
				//a testing scenario when we only have a single run of data
				if(weight[i][0] && !weight[i][1] && !weight[i][2])
				{
					weight[i][0]=1./weight[i][0];
				}
			}		
			
		}else{
			//in data, we already have events with proper weight ratios
			//but we normalize the near
			if(detector==1)
			{
				if(weight[i][0])weight[i][0]=wantPOTs/weight[i][0];
				if(weight[i][1])weight[i][1]=wantPOTs/weight[i][1];
				if(weight[i][2])weight[i][2]=wantPOTs/weight[i][2];
			}

			//far data... don't normalize at all!
			if(detector==0)
			{
				printf("FAR DATA!!!! not normalizing!\n");
				weight[i][0]=1.;
				weight[i][1]=1.;
				weight[i][2]=1.;
			}

		}
	}
	file_list.close();


	//loop over each set of chains

	for(int i=0;i<3;i++)
	{
		if(!weight[i])continue; //chain is empty...
		
		printf("processing chain %d\n",i);
		SetUpChain(infiles[i]);



		
		if(recoType==0)
		{
			if(!isMRCC)
				ProcessChainNormal(infiles[i],weight[i]);
			else
				ProcessChainNormalMRCC(infiles[i],weight[i]);
			
		}
		if(recoType==1)ProcessChainParticlePID(infiles[i],weight[i]);		
		
		
		if(detector==0 && recoType==1 && isMC==1 && isMRCC==0)
			ProcessChain2Pairs(infiles[i],weight[i]);
	
	}


	//save histos
	SaveHistos(workingDir);

	outfile->Write();
	outfile->Close();
}

void MiniPlotMaker::SaveHistos(TDirectory *td)
{
	
	TDirectory *base = td;

	for(int rp=0;rp<4;rp++)
	{

		td=DirectoryHelpers::GetDirectory(base,GetRunName(rp),2);

	if( detector==0 && recoType==1 && isMC==1 && isMRCC==0)
	{
		string x = GetDirectoryString();
		string xx="/";
		vector<string> aa;
		DirectoryHelpers::Tokenize(x,aa,"/");
		for(int i=0;i<(int)aa.size()-2;i++)
			xx+=aa[i]+"/";
		xx+="Combined/"+GetRunName(rp);
		TDirectory *sp=DirectoryHelpers::GetDirectory(outfile,xx,2);
		
		sp->cd();
		char tmp[200];
		for(int i=0;i<6;i++)
		{
			sprintf(tmp,"pid_ParticlePID_vs_ANN11%s",GetNueClassSuffix(i).c_str());
			hist_pids_PID_ann[rp][i]->SetDirectory(sp);
			hist_pids_PID_ann[rp][i]->SetName(tmp);
			hist_pids_PID_ann[rp][i]->SetTitle(tmp);		
		
		}
		
	}
		
		
		

	td->cd();
	
	int nCuts=GetNCuts();
	if(isMRCC)nCuts*=4;
	for(int i=0;i<nCuts;i++)
	{
		string cutname=GetCutName(i%GetNCuts());
		TDirectory *sdir = td->GetDirectory(cutname.c_str());
		if(!sdir)td->mkdir(cutname.c_str());
		sdir = td->GetDirectory(cutname.c_str());	
		sdir->cd();

		if(isMRCC)
		{
			cutname=GetMRCCQPCut(i/GetNCuts());
			TDirectory *sdirt = sdir->GetDirectory(cutname.c_str());
			if(!sdirt)sdir->mkdir(cutname.c_str());
			sdir = sdir->GetDirectory(cutname.c_str());	
			sdir->cd();			
		}
		
		
		char tmp[200];
		int start = isMC ? 0 : 5;
		for(int j=start;j<6;j++)
		{
			if(detector==1 && j==2)continue;//no signal in near
			if(detector==1 && j==3)continue;//no tau in near

			sprintf(tmp,"recoE%s",GetNueClassSuffix(j).c_str());
			hist_recoE[rp][i][j]->SetDirectory(sdir);
			hist_recoE[rp][i][j]->SetName(tmp);
			hist_recoE[rp][i][j]->SetTitle(tmp);
			
		//	if(i==0&&j==0&&rp==3)hist_recoE[rp][i][j]->Dump();

			if(isMC)
			{
				sprintf(tmp,"resCode%s",GetNueClassSuffix(j).c_str());
				hist_resCode[rp][i][j]->SetDirectory(sdir);
				hist_resCode[rp][i][j]->SetName(tmp);
				hist_resCode[rp][i][j]->SetTitle(tmp);		

				sprintf(tmp,"nuEnergy%s",GetNueClassSuffix(j).c_str());
				hist_nuEnergy[rp][i][j]->SetDirectory(sdir);
				hist_nuEnergy[rp][i][j]->SetName(tmp);
				hist_nuEnergy[rp][i][j]->SetTitle(tmp);	
			}
			
			
				sprintf(tmp,"ann11%s",GetNueClassSuffix(j).c_str());
				hist_ann11[rp][i][j]->SetDirectory(sdir);
				hist_ann11[rp][i][j]->SetName(tmp);
				hist_ann11[rp][i][j]->SetTitle(tmp);	
			
				sprintf(tmp,"ann11_firebird%s",GetNueClassSuffix(j).c_str());
				hist_ann11_firebird[rp][i][j]->SetDirectory(sdir);
				hist_ann11_firebird[rp][i][j]->SetName(tmp);
				hist_ann11_firebird[rp][i][j]->SetTitle(tmp);							

                                sprintf(tmp,"ann14%s",GetNueClassSuffix(j).c_str());
                                hist_ann14[rp][i][j]->SetDirectory(sdir);
                                hist_ann14[rp][i][j]->SetName(tmp);
                                hist_ann14[rp][i][j]->SetTitle(tmp);    

			
			
				sprintf(tmp,"ntot%s",GetNueClassSuffix(j).c_str());
				hist_ntot[rp][i][j]->SetDirectory(sdir);
				hist_ntot[rp][i][j]->SetName(tmp);
				hist_ntot[rp][i][j]->SetTitle(tmp);					

				sprintf(tmp,"pidA%s",GetNueClassSuffix(j).c_str());
				hist_pidA[rp][i][j]->SetDirectory(sdir);
				hist_pidA[rp][i][j]->SetName(tmp);
				hist_pidA[rp][i][j]->SetTitle(tmp);					

				sprintf(tmp,"pidB%s",GetNueClassSuffix(j).c_str());
				hist_pidB[rp][i][j]->SetDirectory(sdir);
				hist_pidB[rp][i][j]->SetName(tmp);
				hist_pidB[rp][i][j]->SetTitle(tmp);					

				sprintf(tmp,"pidC%s",GetNueClassSuffix(j).c_str());
				hist_pidC[rp][i][j]->SetDirectory(sdir);
				hist_pidC[rp][i][j]->SetName(tmp);
				hist_pidC[rp][i][j]->SetTitle(tmp);					

				sprintf(tmp,"pidD%s",GetNueClassSuffix(j).c_str());
				hist_pidD[rp][i][j]->SetDirectory(sdir);
				hist_pidD[rp][i][j]->SetName(tmp);
				hist_pidD[rp][i][j]->SetTitle(tmp);					

				sprintf(tmp,"pidE%s",GetNueClassSuffix(j).c_str());
				hist_pidE[rp][i][j]->SetDirectory(sdir);
				hist_pidE[rp][i][j]->SetName(tmp);
				hist_pidE[rp][i][j]->SetTitle(tmp);					

				sprintf(tmp,"pidF%s",GetNueClassSuffix(j).c_str());
				hist_pidF[rp][i][j]->SetDirectory(sdir);
				hist_pidF[rp][i][j]->SetName(tmp);
				hist_pidF[rp][i][j]->SetTitle(tmp);					


				sprintf(tmp,"pidA_recoE%s",GetNueClassSuffix(j).c_str());
				hist_pidA_recoE[rp][i][j]->SetDirectory(sdir);
				hist_pidA_recoE[rp][i][j]->SetName(tmp);
				hist_pidA_recoE[rp][i][j]->SetTitle(tmp);		
							
				sprintf(tmp,"pidB_recoE%s",GetNueClassSuffix(j).c_str());
				hist_pidB_recoE[rp][i][j]->SetDirectory(sdir);
				hist_pidB_recoE[rp][i][j]->SetName(tmp);
				hist_pidB_recoE[rp][i][j]->SetTitle(tmp);	
								
				sprintf(tmp,"pidC_recoE%s",GetNueClassSuffix(j).c_str());
				hist_pidC_recoE[rp][i][j]->SetDirectory(sdir);
				hist_pidC_recoE[rp][i][j]->SetName(tmp);
				hist_pidC_recoE[rp][i][j]->SetTitle(tmp);			
						
				sprintf(tmp,"pidD_recoE%s",GetNueClassSuffix(j).c_str());
				hist_pidD_recoE[rp][i][j]->SetDirectory(sdir);
				hist_pidD_recoE[rp][i][j]->SetName(tmp);
				hist_pidD_recoE[rp][i][j]->SetTitle(tmp);			
						
				sprintf(tmp,"pidE_recoE%s",GetNueClassSuffix(j).c_str());
				hist_pidE_recoE[rp][i][j]->SetDirectory(sdir);
				hist_pidE_recoE[rp][i][j]->SetName(tmp);
				hist_pidE_recoE[rp][i][j]->SetTitle(tmp);				
					
				sprintf(tmp,"pidF_recoE%s",GetNueClassSuffix(j).c_str());
				hist_pidF_recoE[rp][i][j]->SetDirectory(sdir);
				hist_pidF_recoE[rp][i][j]->SetName(tmp);
				hist_pidF_recoE[rp][i][j]->SetTitle(tmp);					
				
				//save file size by only saving these where needed	
				if(recoType==1 ||
					(recoType==0 && (
						i%GetNCuts()==0 || //fiducial				
						i%GetNCuts()==1 || //preselection
						i%GetNCuts()==4 //pidF
					))
				){
					//we don't need these for systematic tests (saves space)
					if(!systematic)
					for(int k=0;k<14;k++)
					{
						sprintf(tmp,"pars_%d%s",k,GetNueClassSuffix(j).c_str());
						hist_pars[k][rp][i][j]->SetDirectory(sdir);
						hist_pars[k][rp][i][j]->SetName(tmp);
						hist_pars[k][rp][i][j]->SetTitle(tmp);		
					}
				}
		
		
				////////
				//lisa histos
//only need these for fid and presel (save a lot of space!)

if(i%GetNCuts()<2)
{				
				if(isMC)MoveOverFlow(hist_lisa_trueE_recoE[rp][i][j]);
				MoveOverFlow(hist_lisa_recoE_pidF[rp][i][j]);
				MoveOverFlow(hist_lisa_recoE_ann11[rp][i][j]);
				MoveOverFlow(hist_lisa_recoE_ann11_firebird[rp][i][j]);
				MoveOverFlow(hist_lisa_recoE_ann14[rp][i][j]);
				
				if(isMC)
				{
					sprintf(tmp,"lisa_trueE_recoE%s",GetNueClassSuffix(j).c_str());
					hist_lisa_trueE_recoE[rp][i][j]->SetDirectory(sdir);
					hist_lisa_trueE_recoE[rp][i][j]->SetName(tmp);
					hist_lisa_trueE_recoE[rp][i][j]->SetTitle(tmp);
				}
				
				
				sprintf(tmp,"lisa_recoE_pidF%s",GetNueClassSuffix(j).c_str());
				hist_lisa_recoE_pidF[rp][i][j]->SetDirectory(sdir);
				hist_lisa_recoE_pidF[rp][i][j]->SetName(tmp);
				hist_lisa_recoE_pidF[rp][i][j]->SetTitle(tmp);		
				
				sprintf(tmp,"lisa_recoE_ann11%s",GetNueClassSuffix(j).c_str());	
				hist_lisa_recoE_ann11[rp][i][j]->SetDirectory(sdir);
				hist_lisa_recoE_ann11[rp][i][j]->SetName(tmp);
				hist_lisa_recoE_ann11[rp][i][j]->SetTitle(tmp);		
				
				sprintf(tmp,"lisa_recoE_ann11_firebird%s",GetNueClassSuffix(j).c_str());	
				hist_lisa_recoE_ann11_firebird[rp][i][j]->SetDirectory(sdir);
				hist_lisa_recoE_ann11_firebird[rp][i][j]->SetName(tmp);
				hist_lisa_recoE_ann11_firebird[rp][i][j]->SetTitle(tmp);		
				
				sprintf(tmp,"lisa_recoE_ann14%s",GetNueClassSuffix(j).c_str());	
				hist_lisa_recoE_ann14[rp][i][j]->SetDirectory(sdir);
				hist_lisa_recoE_ann14[rp][i][j]->SetName(tmp);
				hist_lisa_recoE_ann14[rp][i][j]->SetTitle(tmp);		
}
				//lisa histos	
				////////	
				
				
		}


		//lisas tree for sample information
		//...not everything is being filled now...
		char selection[256];
		sprintf(selection,"%s",GetCutName(i%GetNCuts()).c_str());

		double par[100]={0};
		NueStandard::GetOscParam(par);
		double Theta12 = par[OscPar::kTh12];
		double Theta13 = par[OscPar::kTh13];
		double Theta23 = par[OscPar::kTh23];
		double DeltaMSq23 = par[OscPar::kDeltaM23];
		double DeltaMSq12 = par[OscPar::kDeltaM12];
		double DeltaCP = par[OscPar::kDelta];
		
		
		double fFarPOT = 0;
		double fNearPOT = 0;
		if(detector==0)fFarPOT = wantPOTs;
		if(detector==1)fNearPOT = wantPOTs;
		

	    TTree * paramtree = new TTree("paramtree","paramtree");
	    paramtree->Branch("Selection",selection,"Selection/C");
 	    paramtree->Branch("nearPOT",&fNearPOT,"nearPOT/D");
	    paramtree->Branch("farPOT",&fFarPOT,"farPOT/D");
	    paramtree->Branch("Theta12",&Theta12,"Theta12/D");
	    paramtree->Branch("Theta13",&Theta13,"Theta13/D");
	    paramtree->Branch("Theta23",&Theta23,"Theta23/D");
 		paramtree->Branch("DeltaMSq23",&DeltaMSq23,"DeltaMSq23/D");
	    paramtree->Branch("DeltaMSq12",&DeltaMSq12,"DeltaMSq12/D");
	    paramtree->Branch("DeltaCP",&DeltaCP,"DeltaCP/D");
	    paramtree->Branch("inputPOT0[3]",inputPOT[0]);
            paramtree->Branch("inputPOT1[3]",inputPOT[1]);
            paramtree->Branch("inputPOT2[3]",inputPOT[2]);


		paramtree->SetDirectory(sdir);
		paramtree->Fill();
		//lisa's tree
		/////////

	
	}
	}
}

//move the overflow bins into the normal bin region....
//this is for lisa's histogram fitting
void MiniPlotMaker::MoveOverFlow(TH2D *h)
{
	int nx = h->GetNbinsX();
	int ny = h->GetNbinsY();
	
	for(int i=1;i<nx+1;i++)
	{
		int obin = h->GetBin(i,ny+1);
		int tbin = h->GetBin(i,ny);
		h->AddBinContent(tbin,h->GetBinContent(obin));
		h->AddBinContent(obin,-h->GetBinContent(obin));
	
	}


	for(int i=1;i<ny+1;i++)
	{
		int obin = h->GetBin(nx+1,i);
		int tbin = h->GetBin(nx,i);
		h->AddBinContent(tbin,h->GetBinContent(obin));
		h->AddBinContent(obin,-h->GetBinContent(obin));
	}

	//get the top diagonal!
	int obin = h->GetBin(nx+1,ny+1);
	int tbin = h->GetBin(nx,ny);
	h->AddBinContent(tbin,h->GetBinContent(obin));
	h->AddBinContent(obin,-h->GetBinContent(obin));
	
}

string MiniPlotMaker::GetRunName(int j)
{
	switch(j)
	{
		case 0:return "Run1";
		case 1:return "Run2";
		case 2:return "Run3";
		case 3:return "All";
	};

	return "?";
}

//get a suffix for the histogram name
string MiniPlotMaker::GetNueClassSuffix(int j)
{
	switch(j)
	{
		case 0:
			return "_nc";
		case 1:
			return "_cc";
		case 2:
			return "_sig";
		case 3:
			return "_tau";
		case 4:
			return "_beam";
		case 5:
			return "";
	}

	return "?";
}

int MiniPlotMaker::GetNCuts()
{
	switch (recoType)
	{
		case 0:
			return 6;
		case 1:
			return 8;
	
	};
	return 0;
}

string MiniPlotMaker::GetCutName(int cut)
{
	switch (recoType)
	{
		case 0:
			switch (cut)
			{
				case 0:
					return "fiducial";
				case 1:
					return "presel";
				case 2:
					return "ann11";
				case 3:
					return "ann11_firebird";
				case 4:
					return "pidF";
				case 5:
					return "ann14";
			
			};
		case 1:
			switch (cut)
			{
				case 0:
					return "fiducial";
				case 1:
					return "presel";
				case 2:
					return "pidA";
				case 3:
					return "pidB";			
				case 4:
					return "pidC";
				case 5:
					return "pidD";		
				case 6:
					return "pidE";
				case 7:
					return "pidF";										
			};
		
	};

	return "?";
}


void MiniPlotMaker::SetUpChain(TChain * c)
{

	c->SetBranchStatus("*",1);

	c->SetBranchAddress("NueMiniPID",&nm);

	c->SetBranchStatus("*",0);


 	c->SetBranchStatus("fBeam",1);
 	c->SetBranchStatus("fDet",1);
 	c->SetBranchStatus("fRelease",1);
 	
 	
	c->SetBranchStatus("run",1);
 	c->SetBranchStatus("snarl",1);
 	c->SetBranchStatus("event",1);
 

    c->SetBranchStatus("nueClass",1);
    c->SetBranchStatus("nueOscProb",1);
    c->SetBranchStatus("weight",1);
    c->SetBranchStatus("skzpWeight",1);
    
    c->SetBranchStatus("MCWeight",1);
   	
	c->SetBranchStatus("event_energy",1);
	c->SetBranchStatus("resonanceCode",1);	
	c->SetBranchStatus("nuEnergy",1);	
	
	
	c->SetBranchStatus("timestamp",1);


	if(isMRCC)
	{
		c->SetBranchStatus("mri_qp",1);
		c->SetBranchStatus("mri_orig_cc_pid",1);
		c->SetBranchStatus("mri_SigmaQP",1);
	 
		c->SetBranchStatus("mri_best_complete",1);
		c->SetBranchStatus("mri_fitp",1);
	
	}

	if(recoType==0)
	{
		//preselects
		c->SetBranchStatus("passes_NueStandard_PassesDataQuality",1);          
	    c->SetBranchStatus("passes_NueStandard_IsInFid",1);         
    	c->SetBranchStatus("passes_NueStandard_PassesNonHEPreSelection",1);    
    	c->SetBranchStatus("passes_NueStandard_PassesPreSelection",1);

   		c->SetBranchStatus("evtRecoNueEnergy",1);

       
        c->SetBranchStatus("nshower",1);
        c->SetBranchStatus("contPlanes",1);
        c->SetBranchStatus("cosmicCut",1);
        c->SetBranchStatus("largestEvent",1);
        c->SetBranchStatus("trkPlanes",1);
        c->SetBranchStatus("trkLikePlanes",1);

		//pids
 		c->SetBranchStatus("ann2pe_daikon04",1);
    	c->SetBranchStatus("ann2pe",1);	
		c->SetBranchStatus("pidF",1);

	c->SetBranchStatus("ann14",1);

		if(isMRCC)
		{
			c->SetBranchStatus("passes_NueStandard_PassesMRCCFiducial",1);
			c->SetBranchStatus(
				"passes_NueStandard_PassesMRCCPreSelection",1);  
		}

	}
	
	if(1|| recoType==1 )
	{
		//preselects
		
		c->SetBranchStatus("infid",1);
        c->SetBranchStatus("contained",1);
		c->SetBranchStatus("pars[14]",1);
		
		c->SetBranchStatus("passes_NueStandard_PassesDataQuality",1);          
	    c->SetBranchStatus("passes_NueStandard_IsInFid",1);         
    	c->SetBranchStatus("passes_NueStandard_PassesNonHEPreSelection",1);    
    	c->SetBranchStatus("passes_NueStandard_PassesPreSelection",1);

	    c->SetBranchStatus("longest_s",1);
    	c->SetBranchStatus("event_length",1);
    	c->SetBranchStatus("ntot",1);
   		c->SetBranchStatus("event_energy",1);

		//pids
		c->SetBranchStatus("pidA",1);
		c->SetBranchStatus("pidB",1);
		c->SetBranchStatus("pidC",1);
		c->SetBranchStatus("pidD",1);
		c->SetBranchStatus("pidE",1);
		c->SetBranchStatus("pidF",1);
		
		
 	  	c->SetBranchStatus("pass_var_check",1);
		
		if(isMRCC)
		{
			c->SetBranchStatus("mrcc_s",1);
		
		}
	}

	if(systematic)
	{

		c->SetBranchStatus("nuEnergy",1);
                c->SetBranchStatus("targetEnergy",1);  
                c->SetBranchStatus("targetPX",1);
                c->SetBranchStatus("targetPY",1);
                c->SetBranchStatus("targetPZ",1);
                c->SetBranchStatus("hadronicY",1);
                c->SetBranchStatus("bjorkenX",1);
                c->SetBranchStatus("q2",1);
                c->SetBranchStatus("w2",1);
                c->SetBranchStatus("interactionType",1);
                c->SetBranchStatus("nuFlavor",1);
                c->SetBranchStatus("resonanceCode",1);
                c->SetBranchStatus("initialState",1);
                c->SetBranchStatus("atomicNumber",1);
                c->SetBranchStatus("atomicWeight",1);
                c->SetBranchStatus("hadronicFinalState",1);
	}

}



void MiniPlotMaker::ProcessChain2Pairs(TChain * c, double chainweight[])
{
	printf("(in ProcessChain2Pairs)\n");
	//do the loop
 	c->SetBranchStatus("ann2pe_daikon04",1);
    c->SetBranchStatus("ann2pe",1);		
	
	//int ent=c->GetEntries();
	//if(maxEntries && maxEntries < ent)ent=maxEntries;
	
	//int p100 = ent/10;
	
	
	
	//printf("running on %d entries...\n",ent);
	
	for(int i=0;c->GetEntry(i);i++)
	{
	///	if(i%p100==0)printf("%.2f%%\n",(double)i/(double)ent*100.);
		if(i%100000==0)printf("%d\n",i);
	//	c->GetEntry(i);
		ApplySystematic();


		int myrp=GetRunPeriod(nm->timestamp);
		
		//dq cut
		if(!nm->passes_NueStandard_PassesDataQuality)continue;
			
		
		//fiducial cut
		if(!nm->passes_NueStandard_IsInFid)continue;

		if(isMRCC)
		{
			if(!nm->passes_NueStandard_PassesMRCCFiducial) continue;   
		}

				
		//preselection cut for NUE
		if(!nm->passes_NueStandard_PassesPreSelection)continue;
		
		if(isMRCC)
		{
			if(!nm->passes_NueStandard_PassesMRCCPreSelection)continue;
		}
		/////////////
		


		//preselection cut for ParticlePID
				int pre = 1;
		
		pre = pre && nm->event_length>0.1 && nm->event_length<1.2;
		pre = pre && nm->longest_s>0.1 && nm->longest_s<1.2;
		pre = pre && nm->ntot>0;
		pre = pre && nm->event_energy>0.5 && nm->event_energy<8;	
		
		if(!pre)continue;
		//make sure inputs are good!
		if(!nm->pass_var_check)continue;
		/////////////




	if(isMC && (nm->nueClass<0 || nm->nueClass>4))
	{
//		printf("unknown nue class\n");
		return;
	}
	
	if(!isMC)nm->nueClass=0;
	
	double totweight = chainweight[myrp];
	
	//MCWeight has already been "fixed" (set to 1 for data) by ApplySystematic()
	totweight*=nm->MCWeight;
	
	
	if(isMC)
	{
		if(detector==0)totweight*=nm->nueOscProb;
	}



	for(int rpc=0;rpc<2;rpc++)
	{
		int rp=rpc==0 ? 3 : myrp;

		hist_pids_PID_ann[rp][nm->nueClass]->Fill(nm->pidC,nm->ann2pe_daikon04,totweight);
		hist_pids_PID_ann[rp][5]->Fill(nm->pidC,nm->ann2pe_daikon04,totweight);
	}
	
	}
}

void MiniPlotMaker::ProcessChainNormal(TChain * c, double chainweight[])
{


	
	//do the loop
	
//	int ent=c->GetEntries();
//	if(maxEntries && maxEntries < ent)ent=maxEntries;
	
//	int p100 = ent/10;
	
	
//	printf("running on %d entries...\n",ent);
	
	for(int i=0;c->GetEntry(i);i++)
	{
//		if(i%p100==0)printf("%.2f%%\n",(double)i/(double)ent*100.);
		if(i%100000==0)printf("%d\n",i);
	
//		c->GetEntry(i);
		ApplySystematic();


		int rp=GetRunPeriod(nm->timestamp);
		//dq cut
		
		
		if(!nm->passes_NueStandard_PassesDataQuality)continue;
	
	
		//fiducial cut
		if(!nm->passes_NueStandard_IsInFid)continue;

		if(isMRCC)
		{
			if(!nm->passes_NueStandard_PassesMRCCFiducial) continue;   
		}

		
		FillHistos(0, chainweight[rp], rp);
		
		//preselection cut
		if(!nm->passes_NueStandard_PassesPreSelection)continue;
		
		if(isMRCC)
		{
			if(!nm->passes_NueStandard_PassesMRCCPreSelection)continue;
		}
		
		
		FillHistos(1, chainweight[rp], rp);
		
		//pid cuts..
		
		if(nm->ann2pe_daikon04>0.7)
			FillHistos(2, chainweight[rp], rp);
		
		if(nm->ann2pe>0.7)
			FillHistos(3, chainweight[rp], rp);

                if(nm->ann14>0.7)
                        FillHistos(5, chainweight[rp], rp);


		if(nm->pidF>0.7)
		{

			int passfid=1;
					//fiducial cut
			if(!nm->infid 
				|| !nm->contained 
				|| !nm->passes_NueStandard_IsInFid
			)passfid=0;
		
			//preselection cut
			int pre = 1;
		
			pre = pre && nm->event_length>0.1 && nm->event_length<1.2;
			pre = pre && nm->longest_s>0.1 && nm->longest_s<1.2;
			pre = pre && nm->ntot>0;
			pre = pre && nm->event_energy>0.5 && nm->event_energy<8;	
			
			bool passPP = true;
			passPP = passPP && passfid;
			passPP = passPP && pre;
		
			if(passPP)FillHistos(4, chainweight[rp], rp);

		}	
	
	}


}



void MiniPlotMaker::ProcessChainNormalMRCC(TChain * c, double chainweight[])
{


	
	//do the loop
	
//	int ent=c->GetEntries();
//	if(maxEntries && maxEntries < ent)ent=maxEntries;
	
//	int p100 = ent/10;
	
	
//	printf("running on %d entries...\n",ent);
	
	for(int i=0;c->GetEntry(i);i++)
	{
//		if(i%p100==0)printf("%.2f%%\n",(double)i/(double)ent*100.);
		if(i%100000==0)printf("%d\n",i);
	
//		c->GetEntry(i);
		ApplySystematic();


		int rp=GetRunPeriod(nm->timestamp);
		//dq cut


	for(int qpcut=0;qpcut<4;qpcut++)
	{		

		
		if(!nm->passes_NueStandard_PassesDataQuality)continue;
	
	
		//fiducial cut
		if(!nm->passes_NueStandard_IsInFid)continue;


		if(isMRCC)
		{
			if(!nm->passes_NueStandard_PassesMRCCFiducial) continue; 
			
			
			//do the qp cut;
			bool passesQP=false;
			switch(qpcut)
			{
				case 0:  //no cut
					passesQP=true;break;
				case 1: // numu         i.e. qp<=0 and ccpid>-0.1
					if(nm->mri_qp<=0)passesQP=true;break; 
				case 2: // antinumu01   i.e. qp>0 and ccpid>-0.1 
					if(nm->mri_qp>0)passesQP=true;break; 
				case 3: // antinumu04   i.e. qp>0 and ccpid>0.4 and SigmaQP/QP<0.3 
					if(nm->mri_qp>0 && nm->mri_orig_cc_pid>0.4 && 
						(nm->mri_SigmaQP/nm->mri_qp)<0.3)passesQP=true;break;
			}
			 
				 
			if(!passesQP)continue;
		}

	
		int ncuts = GetNCuts();
		
		FillHistos((qpcut*ncuts)+0, chainweight[rp], rp);
		
		//preselection cut
		if(!nm->passes_NueStandard_PassesPreSelection)continue;
		
		if(isMRCC)
		{
			if(!nm->passes_NueStandard_PassesMRCCPreSelection)continue;
		}
		
		
		FillHistos((qpcut*ncuts)+1, chainweight[rp], rp);
		
		//pid cuts..
		
		if(nm->ann2pe_daikon04>0.7)
			FillHistos((qpcut*ncuts)+2, chainweight[rp], rp);
		
		if(nm->ann2pe>0.7)
			FillHistos((qpcut*ncuts)+3, chainweight[rp], rp);

		if(nm->ann14>0.7)
			FillHistos((qpcut*ncuts)+5, chainweight[rp], rp);


		if(nm->pidF>0.7)
		{

			int passfid=1;
					//fiducial cut
			if(!nm->infid 
				|| !nm->contained 
				|| !nm->passes_NueStandard_IsInFid
			)passfid=0;
		
			//preselection cut
			int pre = 1;
		
			pre = pre && nm->event_length>0.1 && nm->event_length<1.2;
			pre = pre && nm->longest_s>0.1 && nm->longest_s<1.2;
			pre = pre && nm->ntot>0;
			pre = pre && nm->event_energy>0.5 && nm->event_energy<8;	
			
			bool passPP = true;
			passPP = passPP && passfid;
			passPP = passPP && pre;
		
			if(passPP)FillHistos((qpcut*ncuts)+4, chainweight[rp], rp);

		}	
	}
	}


}


void MiniPlotMaker::ProcessChainParticlePID(TChain * c, double chainweight[])
{


	
	//do the loop
	
//	int ent=c->GetEntries();
//	if(maxEntries && maxEntries < ent)ent=maxEntries;
	
//	int p100 = ent/10;
	
	
//	printf("running on %d entries...\n",ent);
	
	for(int i=0;c->GetEntry(i);i++)
	{
//		if(i%p100==0)printf("%.2f%%\n",(double)i/(double)ent*100.);
		if(i%100000==0)printf("%d\n",i);
	
//		c->GetEntry(i);
		ApplySystematic();
		
		int rp=GetRunPeriod(nm->timestamp);

		
		
		//mrcc cut?
		if(isMRCC)
		{
			if(nm->mrcc_s<1.2)continue;   //require removed long muon to be longer than X
		
		}
		
		
		
		//dq cut
		if(!nm->passes_NueStandard_PassesDataQuality)continue;
			
		
		//fiducial cut
		if(!nm->infid 
			|| !nm->contained 
			|| !nm->passes_NueStandard_IsInFid
		)continue;
		
		FillHistos(0, chainweight[rp], rp);
		
		
		//preselection cut
		int pre = 1;
		
		pre = pre && nm->event_length>0.1 && nm->event_length<1.2;
		pre = pre && nm->longest_s>0.1 && nm->longest_s<1.2;
		pre = pre && nm->ntot>0;
		pre = pre && nm->event_energy>0.5 && nm->event_energy<8;	
		
		if(!pre)continue;
		
		
		FillHistos(1, chainweight[rp], rp);
		
		//pid cuts..
		
		//make sure inputs are good!
	//	if(!nm->pass_var_check)continue;
		
		if(nm->pidA>0.7)
			FillHistos(2, chainweight[rp], rp);
		
		if(nm->pidB>0.7)
			FillHistos(3, chainweight[rp], rp);

		if(nm->pidC>0.7)
			FillHistos(4, chainweight[rp], rp);

		if(nm->pidD>0.7)
			FillHistos(5, chainweight[rp], rp);
			
		if(nm->pidE>0.7)
			FillHistos(6, chainweight[rp], rp);
			
		if(nm->pidF>0.7)
		{
			FillHistos(7, chainweight[rp], rp);		
			if(nm->run==13037001)printf("%d %d\n",nm->snarl,nm->event);
		}
	
	}


}



void MiniPlotMaker::MakeHistos()
{

	///////
	//lisa histos params
	
	vector<double>RecoEdges;
	vector<double>TrueEdges;
	vector<double>PIDEdges;
	
	int nReco = 100;
	for(int i=0;i<nReco+1;i++)
	{
		RecoEdges.push_back(i*1.);
	}
	int nTrue = 1200;
	for(int i=0;i<nTrue+1;i++)
	{
		TrueEdges.push_back(i*0.1);
	}
	int nPID = 20;
	for(int i=0;i<nPID+1;i++)
	{
		PIDEdges.push_back(i*0.05);
	}
	
	
	double *redges = new double[nReco+1];
	double *tedges = new double[nTrue+1];
	double *pedges = new double[nPID+1];

	for(int i=0;i<nReco+1;i++)
	{
		redges[i] = RecoEdges[i];
	}
	for(int i=0;i<nTrue+1;i++)
	{
		tedges[i] = TrueEdges[i];
	}
	for(int i=0;i<nPID+1;i++)
	{
		pedges[i] = PIDEdges[i];
	}
	//lisa histos params
	///////


	int cuts=GetNCuts();
	
	if(isMRCC)cuts*=4;

	int start= isMC ? 0 :5;

	char tmp[500];
	
	for(int rp=0;rp<4;rp++)
	{

	for(int j=start;j<6;j++)
	{

			sprintf(tmp,"pids_PID_ann_%d_%d",rp,j);
			hist_pids_PID_ann[rp][j] = new TH2D(tmp,tmp,1000,-0.5,1.5,1000,-0.5,1.5);
//			hist_pids_PID_ann[rp][j]->Sumw2();
		
	for(int i=0;i<cuts;i++)
	{

		
		sprintf(tmp,"recoE_%d_%d_%d",rp,i,j);
		hist_recoE[rp][i][j]=new TH1D(tmp,tmp,20,0,10);
//		hist_recoE[rp][i][j]->Sumw2();


		///////
		//lisa histos
		if(isMC)
		{
			sprintf(tmp,"lisa_trueE_recoE_%d_%d_%d",rp,i,j);
			hist_lisa_trueE_recoE[rp][i][j] = new TH2D(tmp,tmp,nReco,redges,nTrue,tedges);
//			hist_lisa_trueE_recoE[rp][i][j]->Sumw2();
		}
		

		sprintf(tmp,"lisa_recoE_pidF_%d_%d_%d",rp,i,j);
		hist_lisa_recoE_pidF[rp][i][j] = new TH2D(tmp,tmp,nPID,pedges,nReco,redges);
//		hist_lisa_recoE_pidF[rp][i][j]->Sumw2();
		
               

		sprintf(tmp,"lisa_recoE_ann11_%d_%d_%d",rp,i,j);
		hist_lisa_recoE_ann11[rp][i][j] = new TH2D(tmp,tmp,nPID,pedges,nReco,redges);
//		hist_lisa_recoE_ann11[rp][i][j]->Sumw2();
		
              


		sprintf(tmp,"lisa_recoE_ann11_firebird_%d_%d_%d",rp,i,j);
		hist_lisa_recoE_ann11_firebird[rp][i][j] = new TH2D(tmp,tmp,nPID,pedges,nReco,redges);
//		hist_lisa_recoE_ann11_firebird[rp][i][j]->Sumw2();
		
             


		sprintf(tmp,"lisa_recoE_ann14_%d_%d_%d",rp,i,j);
		hist_lisa_recoE_ann14[rp][i][j] = new TH2D(tmp,tmp,nPID,pedges,nReco,redges);
//		hist_lisa_recoE_ann14[rp][i][j]->Sumw2();

            

		//lisa histos
		///////


		if(isMC)
		{
			sprintf(tmp,"resCode_%d_%d_%d",rp,i,j);
			hist_resCode[rp][i][j]=new TH1D(tmp,tmp,5,1001,1006);
//			hist_resCode[rp][i][j]->Sumw2();
		
			sprintf(tmp,"nuEnergy_%d_%d_%d",rp,i,j);
			hist_nuEnergy[rp][i][j]=new TH1D(tmp,tmp,20,0,10);
//			hist_nuEnergy[rp][i][j]->Sumw2();
		}
		
		
			sprintf(tmp,"ann11_%d_%d_%d",rp,i,j);
			hist_ann11[rp][i][j]=new TH1D(tmp,tmp,150,-0.5,1.5);
//			hist_ann11[rp][i][j]->Sumw2();
		
			sprintf(tmp,"ann11_firebird_%d_%d_%d",rp,i,j);
			hist_ann11_firebird[rp][i][j]=new TH1D(tmp,tmp,150,-0.5,1.5);
//			hist_ann11_firebird[rp][i][j]->Sumw2();				
	
                        sprintf(tmp,"ann14_%d_%d_%d",rp,i,j);
                        hist_ann14[rp][i][j]=new TH1D(tmp,tmp,150,-0.5,1.5);
  //                      hist_ann14[rp][i][j]->Sumw2();

	
		
			sprintf(tmp,"ntot_%d_%d_%d",rp,i,j);
			hist_ntot[rp][i][j]=new TH1D(tmp,tmp,100,0,100);
//			hist_ntot[rp][i][j]->Sumw2();			

			sprintf(tmp,"pidA_%d_%d_%d",rp,i,j);
			hist_pidA[rp][i][j]=new TH1D(tmp,tmp,150,-0.5,1.5);
//			hist_pidA[rp][i][j]->Sumw2();	

			sprintf(tmp,"pidB_%d_%d_%d",rp,i,j);
			hist_pidB[rp][i][j]=new TH1D(tmp,tmp,150,-0.5,1.5);
//			hist_pidB[rp][i][j]->Sumw2();	
			
			sprintf(tmp,"pidC_%d_%d_%d",rp,i,j);
			hist_pidC[rp][i][j]=new TH1D(tmp,tmp,150,-0.5,1.5);
//			hist_pidC[rp][i][j]->Sumw2();	
			
			sprintf(tmp,"pidD_%d_%d_%d",rp,i,j);
			hist_pidD[rp][i][j]=new TH1D(tmp,tmp,150,-0.5,1.5);
//			hist_pidD[rp][i][j]->Sumw2();	
			
			sprintf(tmp,"pidE_%d_%d_%d",rp,i,j);
			hist_pidE[rp][i][j]=new TH1D(tmp,tmp,150,-0.5,1.5);
//			hist_pidE[rp][i][j]->Sumw2();	
			
			sprintf(tmp,"pidF_%d_%d_%d",rp,i,j);
			hist_pidF[rp][i][j]=new TH1D(tmp,tmp,150,-0.5,1.5);
//			hist_pidF[rp][i][j]->Sumw2();				


			sprintf(tmp,"pidA_recoE_%d_%d_%d",rp,i,j);
			hist_pidA_recoE[rp][i][j]=new TH2D(tmp,tmp,150,-0.5,1.5,20,0,10);
//			hist_pidA_recoE[rp][i][j]->Sumw2();	
			
			sprintf(tmp,"pidB_recoE_%d_%d_%d",rp,i,j);
			hist_pidB_recoE[rp][i][j]=new TH2D(tmp,tmp,150,-0.5,1.5,20,0,10);
//			hist_pidB_recoE[rp][i][j]->Sumw2();	
			
			sprintf(tmp,"pidC_recoE_%d_%d_%d",rp,i,j);
			hist_pidC_recoE[rp][i][j]=new TH2D(tmp,tmp,150,-0.5,1.5,20,0,10);
//			hist_pidC_recoE[rp][i][j]->Sumw2();	
			
			sprintf(tmp,"pidD_recoE_%d_%d_%d",rp,i,j);
			hist_pidD_recoE[rp][i][j]=new TH2D(tmp,tmp,150,-0.5,1.5,20,0,10);
//			hist_pidD_recoE[rp][i][j]->Sumw2();	
			
			sprintf(tmp,"pidE_recoE_%d_%d_%d",rp,i,j);
			hist_pidE_recoE[rp][i][j]=new TH2D(tmp,tmp,150,-0.5,1.5,20,0,10);
//			hist_pidE_recoE[rp][i][j]->Sumw2();	
			
			sprintf(tmp,"pidF_recoE_%d_%d_%d",rp,i,j);
			hist_pidF_recoE[rp][i][j]=new TH2D(tmp,tmp,150,-0.5,1.5,20,0,10);
//			hist_pidF_recoE[rp][i][j]->Sumw2();	



			for(int k=0;k<14;k++)
			{
				int nbin=1;
				double min=0;
				double max=1;
				switch(k)
				{
					case 0:
						nbin=300;
						min=0;
						max=3;
						break;
					case 1:
						nbin=200;
						min=0;
						max=0.1;
						break;
					case 2:
						nbin=100;
						min=0;
						max=1;
						break;
					case 3:
						nbin=50;
						min=0;
						max=50;
						break;
					case 4:
						nbin=100;
						min=0;
						max=1.2;
						break;
					case 5:
						nbin=100;
						min=0;
						max=1;
						break;
					case 6:
						nbin=100;
						min=0;
						max=1.5;
						break;
					case 7:
						nbin=100;
						min=0;
						max=8000;
						break;
					case 8:
						nbin=100;
						min=0;
						max=500;
						break;
					case 9:
						nbin=100;
						min=-20;
						max=20;
						break;
					case 10:
						nbin=100;
						min=0;
						max=1.5;
						break;
					case 11:
						nbin=100;
						min=0;
						max=6;
						break;
					case 12:
						nbin=100;
						min=0;
						max=50;
						break;
					case 13:
						nbin=100;
						min=0;
						max=0.005;
						break;
				}
				sprintf(tmp,"pars_%d_%d_%d_%d",rp,k,i,j);
				hist_pars[k][rp][i][j]=new TH1D(tmp,tmp,nbin,min,max);
//				hist_pars[k][rp][i][j]->Sumw2();	
			}

		
		
		
	}
	}
	}
	



	delete [] redges;
	delete [] tedges;
	delete [] pedges;

	
}


void MiniPlotMaker::FillHistos(int cut, double weight, int myrp)
{
	if(isMC && (nm->nueClass<0 || nm->nueClass>4))
	{
		printf("unknown nue class\n");
		return;
	}
	
	//do all and also the correct run
	for(int rpc =0;rpc<2;rpc++)
	{	
	
	int rp=rpc ==0 ? 3 : myrp;


	if(!isMC)nm->nueClass=0;
	
	double totweight = weight;
	double osctotweight = totweight;
	if(isMC)
	{
		totweight*=nm->MCWeight;
		osctotweight*=nm->MCWeight;
		if(detector==0)osctotweight*=nm->nueOscProb;
	}


	
	//if we are near, we want things normalized to 1e19
	if(detector==1)
	{
		if(rp<3) //specific run period
			;// already normalized 
		else //for all
		{
			totweight/=3.;	
			osctotweight/=3.;
		}
	}
	
	

	if(isMC)
	{
	
		hist_resCode[rp][cut][5]->Fill(nm->resonanceCode,osctotweight);
		if(isMC)hist_resCode[rp][cut][nm->nueClass]->Fill(nm->resonanceCode,osctotweight);

		hist_nuEnergy[rp][cut][5]->Fill(nm->nuEnergy,osctotweight);
		if(isMC)hist_nuEnergy[rp][cut][nm->nueClass]->Fill(nm->nuEnergy,osctotweight);
		
		///////
		//lisa histos
		//don't oscillate!
		hist_lisa_trueE_recoE[rp][cut][5]->Fill(nm->evtRecoNueEnergy, nm->nuEnergy,totweight);
		if(isMC)hist_lisa_trueE_recoE[rp][cut][nm->nueClass]->Fill(nm->evtRecoNueEnergy, nm->nuEnergy,totweight);
		//lisa histos
		///////
	}
	

	///////
	//lisa histos	

	double recoE=0;
        if(recoType==0 )
        	recoE=nm->evtRecoNueEnergy;
        if(recoType==1 ) // || cut==4 )
                recoE=nm->event_energy;


	hist_lisa_recoE_pidF[rp][cut][5]->Fill(nm->pidF, recoE,osctotweight);
	if(isMC)hist_lisa_recoE_pidF[rp][cut][nm->nueClass]->Fill(nm->pidF, recoE,osctotweight);
		
	hist_lisa_recoE_ann11[rp][cut][5]->Fill(nm->ann2pe_daikon04, recoE,osctotweight);
	if(isMC)hist_lisa_recoE_ann11[rp][cut][nm->nueClass]->Fill(nm->ann2pe_daikon04, recoE,osctotweight);
		
	hist_lisa_recoE_ann11_firebird[rp][cut][5]->Fill(nm->ann2pe, recoE,osctotweight);
	if(isMC)hist_lisa_recoE_ann11_firebird[rp][cut][nm->nueClass]->Fill(nm->ann2pe, recoE,osctotweight);
		
	hist_lisa_recoE_ann14[rp][cut][5]->Fill(nm->ann14, recoE,osctotweight);
	if(isMC)hist_lisa_recoE_ann14[rp][cut][nm->nueClass]->Fill(nm->ann14, recoE,osctotweight);
	//lisa histos
	///////

	if(recoType==0)
	{
	
		hist_ann11[rp][cut][5]->Fill(nm->ann2pe_daikon04,osctotweight);
		if(isMC)hist_ann11[rp][cut][nm->nueClass]->Fill(nm->ann2pe_daikon04,osctotweight);
	
		hist_ann11_firebird[rp][cut][5]->Fill(nm->ann2pe,osctotweight);
		if(isMC)hist_ann11_firebird[rp][cut][nm->nueClass]->Fill(nm->ann2pe,osctotweight);	

		hist_ann14[rp][cut][5]->Fill(nm->ann14,osctotweight);
		if(isMC)hist_ann14[rp][cut][nm->nueClass]->Fill(nm->ann14,osctotweight);


	}

	
	if(recoType==0 )//&& cut<4)
	{
		hist_recoE[rp][cut][5]->Fill(nm->evtRecoNueEnergy,osctotweight);
		if(isMC)hist_recoE[rp][cut][nm->nueClass]->Fill(nm->evtRecoNueEnergy,osctotweight);
	}
	
	
	if(recoType==1 ) // || cut==4 )
	{
		hist_recoE[rp][cut][5]->Fill(nm->event_energy,osctotweight);
		if(isMC)hist_recoE[rp][cut][nm->nueClass]->Fill(nm->event_energy,osctotweight);	
	}
	
	if( 1 || recoType==1 || (recoType==0 && isMC==1 && detector==1))  
	//for debugging with anna's results... want a detail of variables in near mc standard hornON normal
	{			
		hist_ntot[rp][cut][5]->Fill(nm->ntot,osctotweight);
		if(isMC)hist_ntot[rp][cut][nm->nueClass]->Fill(nm->ntot,osctotweight);
		
		
		hist_pidA[rp][cut][5]->Fill(nm->pidA,osctotweight);
		if(isMC)hist_pidA[rp][cut][nm->nueClass]->Fill(nm->pidA,osctotweight);	
			
		hist_pidB[rp][cut][5]->Fill(nm->pidB,osctotweight);
		if(isMC)hist_pidB[rp][cut][nm->nueClass]->Fill(nm->pidB,osctotweight);		
		
		hist_pidC[rp][cut][5]->Fill(nm->pidC,osctotweight);
		if(isMC)hist_pidC[rp][cut][nm->nueClass]->Fill(nm->pidC,osctotweight);		
		
		hist_pidD[rp][cut][5]->Fill(nm->pidD,osctotweight);
		if(isMC)hist_pidD[rp][cut][nm->nueClass]->Fill(nm->pidD,osctotweight);		
		
		hist_pidE[rp][cut][5]->Fill(nm->pidE,osctotweight);
		if(isMC)hist_pidE[rp][cut][nm->nueClass]->Fill(nm->pidE,osctotweight);		
		
		hist_pidF[rp][cut][5]->Fill(nm->pidF,osctotweight);
		if(isMC)hist_pidF[rp][cut][nm->nueClass]->Fill(nm->pidF,osctotweight);
		
		
		hist_pidA_recoE[rp][cut][5]->Fill(nm->pidA,nm->event_energy,osctotweight);
		if(isMC)hist_pidA_recoE[rp][cut][nm->nueClass]->Fill(nm->pidA,nm->event_energy,osctotweight);
		
		hist_pidB_recoE[rp][cut][5]->Fill(nm->pidB,nm->event_energy,osctotweight);
		if(isMC)hist_pidB_recoE[rp][cut][nm->nueClass]->Fill(nm->pidB,nm->event_energy,osctotweight);
		
		hist_pidC_recoE[rp][cut][5]->Fill(nm->pidC,nm->event_energy,osctotweight);
		if(isMC)hist_pidC_recoE[rp][cut][nm->nueClass]->Fill(nm->pidC,nm->event_energy,osctotweight);
		
		hist_pidD_recoE[rp][cut][5]->Fill(nm->pidD,nm->event_energy,osctotweight);
		if(isMC)hist_pidD_recoE[rp][cut][nm->nueClass]->Fill(nm->pidD,nm->event_energy,osctotweight);
		
		hist_pidE_recoE[rp][cut][5]->Fill(nm->pidE,nm->event_energy,osctotweight);
		if(isMC)hist_pidE_recoE[rp][cut][nm->nueClass]->Fill(nm->pidE,nm->event_energy,osctotweight);
		
		hist_pidF_recoE[rp][cut][5]->Fill(nm->pidF,nm->event_energy,osctotweight);
		if(isMC)hist_pidF_recoE[rp][cut][nm->nueClass]->Fill(nm->pidF,nm->event_energy,osctotweight);
		
				

		for(int i=0;i<14;i++)
		{
			double pv=nm->pars[i];
			
			//we don't want to draw all values...
			switch(i)
			{
				case 2:
					if(pv==0)continue;
					break;
				case 6:
					if(pv==0)continue;
					if(pv>1.4499)continue;
					if(pv==0.5)continue;
					break;
				case 7:
					if(pv==0)continue;
					break;
				case 8:
					if(pv==0)continue;
					break;
				case 9:
					if(pv==0)continue;
					break;
				case 10:
					if(pv==0)continue;
					break;
				case 11:
					if(pv<1.5001)continue;
					if(pv>4.9999)continue;
					if(pv==3)continue;
					if(pv==0)continue;
					break;
				case 13:
					if(pv==0)continue;
					break;
			
			}
		
			hist_pars[i][rp][cut][5]->Fill(pv,osctotweight);
			if(isMC)hist_pars[i][rp][cut][nm->nueClass]->Fill(pv,osctotweight);				
		}

	}


}

}



string MiniPlotMaker::GetDirectoryString()
{
	string stage="stage0";
	
	
	stage+="/"+GetSystematic();
	stage+="/"+GetDetectorString();
	stage+="/"+GetMCString();	
	stage+="/"+GetMRCCString();

	if(detector==1)//need to show the hornoff/on 
	{
		stage+="/"+GetHornString();
	}
	
	
	
	stage+="/"+GetRecoTypeString();
	return stage;

}


TDirectory * MiniPlotMaker::MakeDirectory(TFile *f)
{


	string stage=GetDirectoryString();

	return DirectoryHelpers::GetDirectory(f,stage,2);
}


vector<string> MiniPlotMaker::GetFileList(TChain *c)
{
	vector<string> a;
	
	TObjArray *fileElements=c->GetListOfFiles();
    TIter next(fileElements);
    TChainElement *chEl=0;
    while (( chEl=(TChainElement*)next() )) {
         a.push_back(chEl->GetTitle());
    }

	return a;
}

string MiniPlotMaker::GetDirectory()
{
	string reco = GetRecoTypeString();
	string mc = GetMCString();
	string mrcc = GetMRCCString();
	string det = GetDetectorString();

	return "/"+det+"/"+mc+"/"+mrcc+"/"+reco;
}

string MiniPlotMaker::GetRecoTypeString()
{
	switch (recoType)
	{
		case 0:
			return "normal";
		case 1:
			return "ParticlePID";
	};
	return "?";
}


string MiniPlotMaker::GetMCString()
{
	switch(isMC)
	{
		case 0:
			return "data";
		case 1:
			return "MC";
	};
	return "?";
}


string MiniPlotMaker::GetMRCCString()
{
	switch (isMRCC)
	{
		case 0:
			return "standard";
		case 1:
			return "MRCC";
	};
	return "?";
}


string MiniPlotMaker::GetDetectorString()
{
	switch (detector)
	{
		case 0:
			return "far";
		case 1:
			return "near";
	};
	return "?";
}

string MiniPlotMaker::GetHornString()
{
	switch (horn)
	{
		case 0:
			return "horn_on";
		case 1:
			return "horn_off";
	};
	return "?";
}


		
string MiniPlotMaker::GetMRCCQPCut(int cut)
{
	switch(cut)
	{
		case 0:return "mrcc_all";
		case 1:return "mrcc_numu";
		case 2:return "mrcc_antinumu01";
		case 3:return "mrcc_antinumu04";
	}
	return "?";
}
	


void MiniPlotMaker::SetSystematic(int t)
{
	systematic=t;
	
	if(!t)return;/// normal
	
	//otherwise...set up the necessary helpers
	
	

	if(!skzp)skzp=new SKZPWeightCalculator("DetXs",true);
	bool needToAdd=(!mcr) || (!nwc);
	if(!mcr)mcr=&MCReweight::Instance();
	if(!nwc)nwc=new NeugenWeightCalculator();
	if(needToAdd)mcr->AddWeightCalculator(nwc);
	if(!rwtconfig)rwtconfig = new Registry();


}

MiniPlotMaker::MiniPlotMaker()
{

	NueStandard::SetDefaultOscParam();	

	skzp=0;
	mcr=0;
	nwc=0;
	rwtconfig=0;

	lastReweightType=-1;

	//we will need to manually place the histos....
	TH1::AddDirectory(0);
	TH1::SetDefaultSumw2(true);

	printf("MiniPlotMaker....\n");

	for(int i=0;i<3;i++)
	{
		infiles[i]=new TChain("NueMiniPID");
		infilesPOTTree[i]=new TChain("pottree");
		infilesForPOTCounting[i]=new TChain("NueMiniPID");
	}
	
	SetWantPOTs();
	
	recoType=0;
	isMC=0;
	isMRCC=0;
	detector=0;
	horn=0;
	systematic=0;
	
	overwrite=1;
	
	outputFile="";
	
	maxEntries=0;
	
	for(int i=0;i<3;i++)
	{
		infileCount[i]=0;
		infilePOTTreeCount[i]=0;
		infileForPOTCountingCount[i]=0;
	}
	
	nm = new NueMiniPID();


	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
		inputPOT[i][j]=0;

}

void MiniPlotMaker::PrintConfigs()
{
	printf("configuration:\n");
	printf("normalizing to %e POTs\n",wantPOTs);
	printf("output file %s\n",outputFile.c_str());
	printf("working directory %s\n",GetDirectory().c_str());
	
	string setup = GetRecoTypeString() + ' ' + GetMCString() + ' ' + 
		GetMRCCString() + ' ' + GetDetectorString() + ' ' + GetHornString() + ' ' + GetSystematic();

	printf("%s\n",setup.c_str());
	
	fstream file_list("filelist",ios::out);
	
	char tmp[500];
	sprintf(tmp,"normalizing to %e POTs\n",wantPOTs);
	file_list<<tmp;
	sprintf(tmp,"output file %s\n",outputFile.c_str());
	file_list<<tmp;
	sprintf(tmp,"working directory %s\n",GetDirectory().c_str());
	file_list<<tmp;
	sprintf(tmp,"%s\n",setup.c_str());
	file_list<<tmp;	
	
	
	file_list<<"input files:\n";
	for(int i=0;i<3;i++)
	{
		if(!infileCount[i])continue;
		
		file_list<< "-------------------\nGroup "<<i<<"\n";
		vector<string> filelist=GetFileList(infiles[i]);
		file_list<<"  with "<<filelist.size()<<" files\n";
		for(int i=0;i<(int)filelist.size();i++)
			file_list<<filelist[i].c_str()<<"\n";
	
		file_list<<"\n\n";
	}
	
	file_list.close();
}


string MiniPlotMaker::GetSystematic()
{
	switch (systematic)
	{
		case 0:
			return "full";
		case 1:
			return "shower_e_plus_11percent";
		case 2:
			return "shower_e_minus_11percent";
		case 3:
			return "flux_error_weight_increased";
		case 4:
			return "flux_error_weight_decreased";
		case 5:
			return "cross_section_qema_plus_15percent";
		case 6:
			return "cross_section_qema_minus_15percent";
		case 7:
			return "cross_section_resma_plus_15percent";
		case 8:
			return "cross_section_resma_minus_15percent";
		case 9:
			return "cross_section_kno_plus_50percent";
		case 10:
			return "cross_section_kno_minus_50percent";
	};
	return "?";


}

void MiniPlotMaker::ApplySystematic()
{
		
	//calibrate ParticlePID energy!
		
	if(detector==0)
	{
		float offset =   0.489987;
       	float slope  =  0.0387301;
        nm->event_energy = nm->event_energy*slope+offset;
	}else if(detector==1)
	{
		float offset = 0.4803261;
		float slope = 0.03799819;
		nm->event_energy = nm->event_energy*slope+offset;
	}

	switch (systematic)
	{
		case 0: //none
			return;
			
		case 1: //shower_e_plus_11percent
				nm->event_energy*=1.11;
				nm->evtRecoNueEnergy*=1.11;
				nm->passes_NueStandard_PassesPreSelection=PassesNuePreselection();
			return;
			
		case 2: //shower_e_minus_11percent
				nm->event_energy*=0.89;
				nm->evtRecoNueEnergy*=0.89;
                                nm->passes_NueStandard_PassesPreSelection=PassesNuePreselection();

			return;
			
		case 3: //flux_error_weight_increased
				nm->MCWeight=(isMC?nm->MCWeight:1)*skzp->GetFluxError(1,2,nm->nuFlavor,nm->nuEnergy);	
			return;
			
		case 4: //flux_error_weight_decreased
				nm->MCWeight=(isMC?nm->MCWeight:1)*(2-skzp->GetFluxError(1,2,nm->nuFlavor,nm->nuEnergy));
			return;
			
		case 5: //cross_section_qema_plus_15percent
				nm->MCWeight=(isMC?nm->MCWeight:1)*GetReweight(0);
			return;
			
		case 6: //cross_section_qema_minus_15percent
				nm->MCWeight=(isMC?nm->MCWeight:1)*GetReweight(1);
			return;
			
		case 7: //cross_section_resma_plus_15percent
				nm->MCWeight=(isMC?nm->MCWeight:1)*GetReweight(2);
			return;
			
		case 8: //cross_section_resma_minus_15percent
				nm->MCWeight=(isMC?nm->MCWeight:1)*GetReweight(3);
			return;
			
		case 9: //cross_section_kno_plus_50percent
				nm->MCWeight=(isMC?nm->MCWeight:1)*GetReweight(4);
			return;
			
		case 10: //cross_section_kno_minus_50percent
				nm->MCWeight=(isMC?nm->MCWeight:1)*GetReweight(5);
			return;
	
	};
}


double MiniPlotMaker::GetReweight(int type)
{

	bool reset=false;
	if(lastReweightType!=type)
	{
		rwtconfig->UnLockValues();
		rwtconfig->UnLockKeys();
		rwtconfig->Clear();
		rwtconfig->Set("neugen:use_scale_factors",1);
	
		switch(type)
		{
			case 0:
				rwtconfig->Set("neugen:ma_qe",1.15);break;
			case 1:
	 			    rwtconfig->Set("neugen:ma_qe",0.85);break;
			case 2:
		      	    rwtconfig->Set("neugen:ma_res",1.15);break;
			case 3:
				rwtconfig->Set("neugen:ma_res",0.85);break;
			case 4:
				rwtconfig->Set("neugen:scale_kno_all",1.5);break;
			case 5:
				rwtconfig->Set("neugen:scale_kno_all",0.5);break;
		
		};
		rwtconfig->LockValues();
		rwtconfig->LockKeys();
		lastReweightType=type;
		reset=true;
	}
	
	
    MCEventInfo ei;

    ////////////////fill truth info for xsec reweighting/////////////
    ei.UseStoredXSec(false);
    ei.nuE=nm->nuEnergy;
    ei.nuPx=nm->nuDCosX*sqrt(nm->nuEnergy*nm->nuEnergy);
    ei.nuPy=nm->nuDCosY*sqrt(nm->nuEnergy*nm->nuEnergy);
    ei.nuPz=nm->nuDCosZ*sqrt(nm->nuEnergy*nm->nuEnergy);
    ei.tarE=nm->targetEnergy;
    ei.tarPx=nm->targetPX;
    ei.tarPy=nm->targetPY;
    ei.tarPz=nm->targetPZ;
    ei.y=nm->hadronicY;
    ei.x=nm->bjorkenX;
    ei.q2=nm->q2;
    ei.w2=nm->w2;
    ei.iaction=nm->interactionType;
    ei.inu=nm->nuFlavor;
    ei.iresonance=nm->resonanceCode;
    ei.initial_state=nm->initialState;
    ei.nucleus=ReweightHelpers::FindNucleusNumber((int)nm->atomicNumber,
                                                    (int)nm->atomicWeight);
    ei.had_fs=nm->hadronicFinalState;

  if(ei.inu<-9998 ||
     (ei.iresonance==1003 && TMath::Abs(ei.had_fs)<200) ) return 1;


    NuParent *nuparent = 0;
    double xsw=1;

    if(reset)xsw = mcr->ComputeWeight(&ei,nuparent,rwtconfig);
   else xsw = mcr->ComputeWeight(&ei,nuparent);


	return xsw;
}


bool MiniPlotMaker::PassesNuePreselection()
{
	bool pass=1;
	pass = pass && (nm->evtRecoNueEnergy>1 && nm->evtRecoNueEnergy<8);
	pass = pass && (nm->nshower>0);
	pass = pass && (nm->contPlanes>4);
	pass = pass && (nm->cosmicCut==1);
	pass = pass && (detector!=0 || nm->largestEvent==1);
	pass = pass && (nm->trkPlanes<25);
	pass = pass && (nm->trkLikePlanes<16);
	return pass;
}


