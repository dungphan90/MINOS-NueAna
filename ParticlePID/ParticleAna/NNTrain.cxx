////////////////////////
//
//  this is for training and testing NN pids...
//
//
/*
		train like this:
		NNTrain n
		n.SetInputFile("/minos/data/dogwood0/official/far_MC_dogwood1_withXtalkFilterOn_Griffin_PreProduction_v00/testset/AnaNue*")
		n.SetTrainFile("nout.root")
		n.Train(20,5,6,"14:7");

		test like this:
		NNTrain n
		n.SetTrainFile("trainfile.root") 
		n.SetInputFile(0,"/minos/data/dogwood0/official/far_MC_dogwood1_withXtalkFilterOn_Griffin_PreProduction_v00/testset2/AnaNue-f210*")
		n.SetInputFile(1,"/minos/data/dogwood0/official/far_MC_dogwood1_withXtalkFilterOn_Griffin_PreProduction_v00/testset2/AnaNue-f213*")
		n.SetInputFile(2,"/minos/data/dogwood0/official/far_MC_dogwood1_withXtalkFilterOn_Griffin_PreProduction_v00/testset2/AnaNue-f214*")
		n.MakeTestTree() 
*/
//
////////////////////////



#include "NueAna/ParticlePID/ParticleAna/NNTrain.h"
#include <TTree.h>
#include <TMath.h>

#include "TMLPAnalyzer.h"
#include <TLegend.h>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TRandom.h"

#include "NueAna/ParticlePID/ParticleAna/PRecord.h"
#include "NueAna/NueStandard.h"




//comment out dodraw line for grid machines which hate graphics...
//#define dodraw 1

NNTrain::NNTrain()
{
	for(int i=0;i<3;i++)
		trainfile[i].clear();
	ResetTrainParams();
	mlp=0;
	needToWeight=0;
		
	//osc
	fBaseLine = 735;
  	fDeltaMS12= 8.0;
	fTh12= 0.816;
	fTh23=3.141592/4.0;
	fDensity=2.75;
	SetOscParamBase(0.00243, TMath::ASin(TMath::Sqrt(0.15))/2.0, 0, 1);	
	
	nr=new NueRecord();


	traincontinue=0;

  f2=new TF1("f2","sin(1.267*[0]*[1]/x)*sin(1.267*[0]*[1]/x)",0.,120.);
  float delm2=0.0024;
  float Ltd=735.;  
  f2->SetParameters(delm2,Ltd);

}


NNTrain::~NNTrain()
{}


void NNTrain::ResetTrainParams()
{
	delta=0;
	epsilon=0;
	eta=0;
	etadecay=0;
	tau=0;
}

//make a test tree will all but
// rejected events (due to presel, etc)
void NNTrain::MakeTestTree()
{
	
	MakeTrainTree(1);

}

//attempts to open file to find a sim tree
//then attempt to fill tree with pid from MLP
void NNTrain::FillTreePid(string file,string outfile, string MLPfile)
{

	TFile * mlpfile = TFile::Open(MLPfile.c_str());
	TMultiLayerPerceptron * mlp = 0;
	if(mlpfile)
	{
		mlp=(TMultiLayerPerceptron*)mlpfile->Get("MLP")->Clone();
		mlpfile->Close();
	}


	if(!mlp)
	{
		printf("unable to load MLP\n");
		exit(1);
	}


	TFile * fin = TFile::Open(file.c_str(), "UPDATE");
	if(!fin)
	{
		printf("no input file specified\n");
		exit(1);
	}

	TTree * tree = (TTree*)fin->Get("TrainTree");
	if(!tree)
	{
		printf("no traintree tree found in input file!\n");
		exit(1);
	}

	TTree * descriptor = (TTree*)fin->Get("descriptor");
	if(!descriptor)
	{
		printf("no descriptor tree found in input file!\n");
		exit(1);
	}


	double pid;
	double pars[100];

	TFile * fout = TFile::Open(outfile.c_str(),"RECREATE");
	fout->cd();
	TTree * treeout = (TTree*)tree->CloneTree(0);
	TTree * descriptorout = (TTree*)descriptor->CloneTree();

	if(treeout->GetBranch("pid"))
		treeout->SetBranchAddress("pid",&pid);
	else
		treeout->Branch("pid",&pid);


	tree->SetBranchAddress("pars",&pars);
	treeout->SetBranchAddress("pars",&pars);

	int ent = tree->GetEntries();
	printf("Filling pid for %d entries\n",ent);
	for(int i=0;i<ent;i++)
	{
		tree->GetEntry(i);
		pid=mlp->Evaluate(0,pars);
		treeout->Fill();
	}
	treeout->Write("TrainTree",TObject::kOverwrite);
	descriptorout->Write();
	fin->Close();
	fout->Close();

}

TTree * NNTrain::MakeTrainTree(int makeTestTree)  //0 to bypass random osc selection
{
	int files=0;
	for(int i=0;i<3;i++)
		files+=trainfile[i].size();	
	if(files<1)
	{
		printf("Specify a input file!\n");
		return 0;
	}

	double filePOTs[3];
	for(int i=0;i<3;i++)
	{
		filePOTs[i]=0;
		TChain *pt = new TChain("pottree");
		for(unsigned int j=0;j<trainfile[i].size();j++)
			pt->Add(trainfile[i][j].c_str());
		
		int ent = pt->GetEntries();
		for(int j=0;j<ent;j++)
		{
			pt->GetEntry(j);
			filePOTs[i]+= pt->GetLeaf("pot")->GetValue();	
		}
	
		printf("file type %d has %f POTs\n",i,filePOTs[i]);
	}



	TChain *c = new TChain("ana_nue");

	c->SetBranchAddress("NueRecord",&nr);

	//there is a problem with turning off all branches together "*"... maybe 
	//because of the header?.... so turn them all off manually!
	c->SetBranchStatus("*",1);
	
	c->SetBranchStatus("shwfit*",0);
	c->SetBranchStatus("hitcalc*",0);
	c->SetBranchStatus("angcluster*",0);
	c->SetBranchStatus("angclusterfit*",0);
	c->SetBranchStatus("mstvars*",0);
	c->SetBranchStatus("fracvars*",0);
	c->SetBranchStatus("subshowervars*",0);
	c->SetBranchStatus("highhitvars*",0);
	c->SetBranchStatus("shieldrejvars*",0);
	c->SetBranchStatus("ann*",0);
	c->SetBranchStatus("anainfo*",0);
	c->SetBranchStatus("srevent*",0);
	c->SetBranchStatus("srshower*",0);
	c->SetBranchStatus("srtrack*",0);
	c->SetBranchStatus("mctrue*",0);
 	//c->SetBranchStatus("bmon*",0);
	c->SetBranchStatus("mdadiscri*",0);;
	c->SetBranchStatus("treepid*",0);
	c->SetBranchStatus("fluxinfo*",0);
	c->SetBranchStatus("fluxweights*",0);
	c->SetBranchStatus("xsecweights*",0);
	c->SetBranchStatus("shi*",0);
	c->SetBranchStatus("mri*",0);
	c->SetBranchStatus("cdi*",0);
	c->SetBranchStatus("timingvars*",0);
	c->SetBranchStatus("mcnnv*",0);
	c->SetBranchStatus("dtree*",0);
	//c->SetBranchStatus("eventq*",0);
	c->SetBranchStatus("precord*",0);
	c->SetBranchStatus("precordMRCC*",0);
	
	
	c->SetBranchStatus("precord.particles*",1);
	c->SetBranchStatus("precord.event*",1);
	
	c->SetBranchStatus("fracvars.fract_2_planes",1);
	c->SetBranchStatus("fracvars.fract_4_planes",1);	
	c->SetBranchStatus("fracvars.fract_6_planes",1);	
	c->SetBranchStatus("fracvars.fract_8_counters",1);	
	c->SetBranchStatus("fracvars.fract_road",1);

	c->SetBranchStatus("shwfit.par_a",1);
	c->SetBranchStatus("shwfit.par_b",1);
	c->SetBranchStatus("shwfit.LongE",1);
	c->SetBranchStatus("shwfit.uv_molrad_peak_9s_2pe_dw",1);
	c->SetBranchStatus("shwfit.uv_rms_9s_2pe_dw",1);
	
	c->SetBranchStatus("mstvars.e4w",1);
	c->SetBranchStatus("mstvars.o4w",1);	
	
	c->SetBranchStatus("fluxweights.totbeamweight",1);
	c->SetBranchStatus("fluxweights.RPtotbeamweight",1);
	
	c->SetBranchStatus("mctrue.fNueClass",1);
	c->SetBranchStatus("mctrue.fOscProb",1);
	c->SetBranchStatus("mctrue.resonanceCode",1);
	c->SetBranchStatus("mctrue.nuEnergy",1);
	c->SetBranchStatus("mctrue.emShowerFraction",1);
	c->SetBranchStatus("mctrue.trueVisibleE",1);
	c->SetBranchStatus("mctrue.nuFlavor",1);
	c->SetBranchStatus("mctrue.nonOscNuFlavor",1);
	c->SetBranchStatus("mctrue.interactionType",1);
	
	c->SetBranchStatus("shwfit.contPlaneCount050",1);
	c->SetBranchStatus("shwfit.UVSlope",1);
	c->SetBranchStatus("srtrack.endX",1);
	c->SetBranchStatus("srtrack.endY",1);
	c->SetBranchStatus("srtrack.endZ",1);
	c->SetBranchStatus("srtrack.vtxX",1);
	c->SetBranchStatus("srtrack.vtxY",1);
	c->SetBranchStatus("srtrack.vtxZ",1);
	
	c->SetBranchStatus("srtrack.endPlane",1);
	c->SetBranchStatus("srtrack.begPlane",1);
	c->SetBranchStatus("srtrack.trklikePlanes",1);
	c->SetBranchStatus("srevent.tracks",1);
	c->SetBranchStatus("srevent.largestEvent",1);
	c->SetBranchStatus("srevent.showers",1);
	c->SetBranchStatus("srevent.phNueGeV",1);
	c->SetBranchStatus("srevent.vtxX",1);
	c->SetBranchStatus("srevent.vtxY",1);
	c->SetBranchStatus("srevent.vtxZ",1);
	
	
	
	
	
	Double_t pid=0;
	TMultiLayerPerceptron *mlp=0;
	if(makeTestTree==1)
	{
		TFile *fmlp = TFile::Open("MLP.root");
		if(fmlp)mlp = (TMultiLayerPerceptron*) fmlp->Get("MLP")->Clone();
		if(!mlp)
		{
			printf("NO MLP FOUND when making test tree\n");
		//	return 0;
		
		}
		if(fmlp)fmlp->Close();
		
	}
	
	
	
	TFile *fout=TFile::Open(inputFile.c_str(),"RECREATE");
	fout->cd();
	TTree *sim = new TTree("TrainTree","TrainTree");
	
	
	sim->Branch("mctrue_oscprob",&mctrue_oscprob);
	sim->Branch("mctrue_totbeamweight",&mctrue_totbeamweight);
	sim->Branch("type",&type);
	sim->Branch("mctrue_fNueClass",&mctrue_type);
	sim->Branch("weight",&weight);
	sim->Branch("pars",&pars,"pars[30]/D");
	
	sim->Branch("mctrue_nuEnergy",&mctrue_nuenergy);	
	sim->Branch("visenergy",&event_visenergy);
	sim->Branch("resonanceCode",&mctrue_iresonance);
	sim->Branch("tweight",&tweight);
	sim->Branch("emfrac",&emfrac);

	Int_t preselect;
	sim->Branch("strict_preselection",&preselect);	
	
	if(mlp && makeTestTree==1)
	{
		sim->Branch("pid",&pid);	
	}
	
	
	for(int i=0;i<20;i++)pars[i]=0;
	
	for(unsigned int k=0;k<3;k++)
	for(unsigned int i=0;i<trainfile[k].size();i++)
	{ 
		int added =c->Add(trainfile[k][i].c_str());
		if(!added)
		{
			printf("input file %s bad\n",trainfile[k][i].c_str());
			return 0;
		}else{
			printf("added %d files from %s\n",added,trainfile[k][i].c_str());
		}
	}


	c->GetEntry(0);
	int det = nr->GetHeader().GetVldContext().GetDetector();
	int mc = nr->GetHeader().GetVldContext().GetSimFlag();
	printf("detector %d mc %d\n",det,mc);



	NueStandard::SetDefaultOscParam();
	NueStandard::SetOscNoMatter(); 


	//dump the osc params...
	double oscparams[20];
	NueStandard::GetOscParam(oscparams);
	if(det==Detector::kFar)
	for(int oi = 0; oi < int(OscPar::kNumParameters); oi++){
    		printf("oscpar %d: %f\n",oi,oscparams[oi]);
    	}



	int npars=0;//so we can automate the adjustment for training...

	int ent = c->GetEntries();
	printf("%d total entries\n",ent);

	//determine max osc*totbeam for each type if needed
	double type_weight[5];
	for(int i=0;i<5;i++)type_weight[i]=1;
	if(!makeTestTree && det==Detector::kFar && mc==SimFlag::kMC)
	{
		for(int i=0;i<5;i++)type_weight[i]=0;
		for(int i=0;i<ent;i++)
		{
			c->GetEntry(i);
			if(nr->mctrue.fNueClass<0||nr->mctrue.fNueClass>4)continue;


			double w = nr->mctrue.fOscProb*nr->fluxweights.totbeamweight;
			if(w>type_weight[nr->mctrue.fNueClass])type_weight[nr->mctrue.fNueClass]=w;
		}
		for(int i=0;i<5;i++)printf("type_weight %d %f\n",i,type_weight[i]);
		
	}


	int failedpval=0;
	int failedpp=0;
	int failedspp=0;
	int faileddq=0;





	int require_nue_presel=0;
	int require_pid_presel=0;
	int require_strict_pid_presel=1;
	int checkNueVars=0;
	int checkParticleVars=1;

	for(int i=0;i<ent;i++)
	{
		c->GetEntry(i);
		if(i%10000==0)printf("%2.2f%%...\n",100.*i/ent);

	        det = nr->GetHeader().GetVldContext().GetDetector();
	        mc = nr->GetHeader().GetVldContext().GetSimFlag();

		if(mc!=SimFlag::kMC)
			nr->mctrue.fNueClass=0;	

		mctrue_type=nr->mctrue.fNueClass;
	       
		
		if(mctrue_type<0||mctrue_type>4)continue;



		mctrue_nuenergy = nr->mctrue.nuEnergy;
		mctrue_iresonance = nr->mctrue.resonanceCode;
		event_visenergy = nr->precord.event.visenergy;
		emfrac = nr->mctrue.emShowerFraction;
		
		///////////////////////
		//// calculate file and beam weights

				
		double pots=0;
		if(det==Detector::kNear)pots=filePOTs[0];
		else
		{


		if(nr->mctrue.interactionType==0)
		{
			pots=filePOTs[0]+filePOTs[1]+filePOTs[2];
		}else if (abs(nr->mctrue.nonOscNuFlavor)==14&&
					abs(nr->mctrue.nuFlavor)==14)
		{
			pots=filePOTs[0];
		}else if (abs(nr->mctrue.nonOscNuFlavor)==12&&
    				abs(nr->mctrue.nuFlavor)==12)
    		{
    			pots=filePOTs[0];
    		}else if (abs(nr->mctrue.nonOscNuFlavor)==12&&
    				abs(nr->mctrue.nuFlavor)==14)
    		{
    			pots=filePOTs[2];
   	 	}else if (abs(nr->mctrue.nonOscNuFlavor)==14&&
    				abs(nr->mctrue.nuFlavor)==12)
    		{
    			pots=filePOTs[2];
   	 	}else if (abs(nr->mctrue.nonOscNuFlavor)==14&&
    				abs(nr->mctrue.nuFlavor)==16)
	    	{
	    		pots=filePOTs[1];
  	  	}else if (abs(nr->mctrue.nonOscNuFlavor)==12&&
    				abs(nr->mctrue.nuFlavor)==16)
   	 	{
    			pots=filePOTs[1];
    		}
		}


		
		if(!pots)continue;


		
		mctrue_totbeamweight= (mc==SimFlag::kData) ? 1 :
			NueStandard::GetRPWBeamWeight(nr);		

		if(det == Detector::kNear && mc == SimFlag::kData)
			if(!NueStandard::PassesDataQuality(nr)){faileddq++;continue;}	

	
		//// end file and beam weights
		///////////////////////
		
		

		///////////////////////
		//// reoscillation!	
				
		mctrue_oscprob=1;

		if(det == Detector::kFar && mc == SimFlag::kMC)
		{
			mctrue_oscprob = (nr->mctrue.fNueClass==0) ? 1 :
		//	osc( nr->mctrue.nuEnergy, nr->mctrue.fNueClass,nr->mctrue.nonOscNuFlavor,nr->mctrue.nuFlavor);
			//fOscCalc.Oscillate(nr->mctrue.nuFlavor, nr->mctrue.nonOscNuFlavor, nr->mctrue.nuEnergy);

			NueStandard::GetOscWeight(nr->mctrue.nuFlavor,nr->mctrue.nonOscNuFlavor,nr->mctrue.nuEnergy); 
			if(mctrue_oscprob<-1000)continue;
		}


		//// end reoscillation!	
		///////////////////////




		///////////////////////		
		//// trainweight
		
				
		tweight = 0;
	
		if(det == Detector::kFar)
		{		
		if(nr->mctrue.interactionType==0)
		{
			tweight=filePOTs[0]/(filePOTs[0]+filePOTs[2]);
		}else if (abs(nr->mctrue.nonOscNuFlavor)==14&&
				abs(nr->mctrue.nuFlavor)==14)
		{
			tweight=mctrue_oscprob;
		}else if (abs(nr->mctrue.nonOscNuFlavor)==12&&
 					abs(nr->mctrue.nuFlavor)==12)
	 	{
   		}else if (abs(nr->mctrue.nonOscNuFlavor)==12&&
   				abs(nr->mctrue.nuFlavor)==14)
   		{
   			tweight=mctrue_oscprob;
	 	}else if (abs(nr->mctrue.nonOscNuFlavor)==14&&
   				abs(nr->mctrue.nuFlavor)==12)
   		{
   			tweight=mctrue_oscprob/0.075;
	  	}else if (abs(nr->mctrue.nonOscNuFlavor)==14&&
 				abs(nr->mctrue.nuFlavor)==16)
   		{
   	 	}else if (abs(nr->mctrue.nonOscNuFlavor)==12&&
   				abs(nr->mctrue.nuFlavor)==16)
   		{
   		}
		}
			
			
		//make the sample realistic!
		if(!makeTestTree && det==Detector::kFar && mc == SimFlag::kMC)
		{	
		//	double nueweight=1.3;//20;
		//	if(mctrue_type==2)tweight= tweight*20.;
			
		
			if( tweight<gRandom->Uniform()/**type_weight[mctrue_type]*/)continue;
		}


		//// end trainweight
		///////////////////////		




		weight=3.25*(1e8)/pots;





		PRecord * pr = &nr->precord;


		trueEMFrac = nr->mctrue.emShowerFraction;
		//trueEMFrac = truthcompare_emenergy;//10000*truthcompare_emenergy/mctrue_visible_energy;//
		//truthcompare_emenergy/event_visenergy*25.;
		if(trueEMFrac>1)trueEMFrac=1;
		




		///////
		//compute precord variables
		
		int pass_precord_preselect=1;
		preselect=1;
		
		length_z = pr->event.max_z - pr->event.min_z;
		if(pr->particles.totvise==0)pass_precord_preselect=0;
		largest_frac=pr->particles.totvise ? pr->particles.largest_particle_e/pr->particles.totvise : 1;
	//	reco_frac=pr->particles.prim_vise/pr->event.visenergy;
		ntot_lsps=pr->particles.ntot*pr->particles.longest_s_particle_s;		
	//	longest_ae0=pr->particles.longest_s_particle_par_a/pr->particles.longest_s_particle_par_e0;

		prim_ae0=pr->particles.prim_par_e0 ? pr->particles.prim_par_a/pr->particles.prim_par_e0 : 0;

		//
		//////
		
		if(pr->event.inFiducial!=1){pass_precord_preselect=0;preselect=0;}
		if(pr->event.contained!=1){pass_precord_preselect=0;preselect=0;} 
		if(pr->particles.ntot<1){pass_precord_preselect=0;preselect=0;}

		
		//weak preselection for training
		if(pr->particles.longest_s_particle_s>2)pass_precord_preselect=0;
		if(pr->particles.longest_s_particle_s<0.0 )pass_precord_preselect=0;
		
		if(pr->event.max_z-pr->event.min_z>2)pass_precord_preselect=0;
		if(pr->event.max_z-pr->event.min_z<0.0)pass_precord_preselect=0;
		
		if(pr->event.visenergy/25.<0.0 )pass_precord_preselect=0;
		if(pr->event.visenergy/25.>10.)pass_precord_preselect=0;


		if(!pass_precord_preselect&&require_pid_presel){failedpp++;continue;}

		//strict preselection
		
		if(pr->particles.longest_s_particle_s>1.2)preselect=0;
		if(pr->particles.longest_s_particle_s<0.1)preselect=0;
		
		if(pr->event.max_z-pr->event.min_z>1.2)preselect=0;
		if(pr->event.max_z-pr->event.min_z<0.1)preselect=0;
		
		if(pr->event.visenergy/25.<0.5)preselect=0;
		if(pr->event.visenergy/25.>8)preselect=0;

		if(!preselect&&require_strict_pid_presel){failedspp++;continue;}

//		if(!pass_precord_preselect)...;

//		if(particles_largest_particle_cmp_chisq>30)...;

	//	if(event_nstrips<10 || event_nstrips>40)...;

		//don't want a long muon..
//		if(particles_longest_particle_type>12)...;		
		
		


		int pass = 1;

		//for nr use only!
		
		pass = pass && NueStandard::IsInFid(nr);
		pass = pass && NueStandard::PassesPreSelection(nr);
		
		if(!pass && require_nue_presel)continue;




		
 
		if (nr->shwfit.par_a<-1000) pass=0;//{printf("nspa\n");pass=0;}
		if (nr->shwfit.par_b<-1000)pass=0;//{printf("nspb\n");pass=0;}
		if (nr->shwfit.uv_molrad_peak_9s_2pe_dw<-1000)pass=0;//{printf("nsum\n");pass=0;}
		if (nr->shwfit.uv_rms_9s_2pe_dw<-1000) pass=0;//{printf("nsur\n");pass=0;}
		if (nr->mstvars.e4w<-1000) pass=0;//{printf("nme\n");pass=0;}
		if (nr->mstvars.e4w>500)pass=0;// {printf("nme\n");pass=0;}
		if (nr->mstvars.o4w<-1000)pass=0;//{printf("nmo\n");pass=0;}
		if (nr->mstvars.o4w>500) pass=0;//{printf("nmo\n");pass=0;}
		if (nr->fracvars.fract_2_planes<-1000)pass=0;// {printf("ff2\n");pass=0;}
		if (nr->fracvars.fract_4_planes<-1000)pass=0;// {printf("ff4\n");pass=0;}
		if (nr->fracvars.fract_6_planes<-1000)pass=0;// {printf("ff6\n");pass=0;}
		if (nr->fracvars.fract_8_counters<-1000) pass=0;//{printf("ff8\n");pass=0;}
		if (nr->fracvars.fract_road<-1000) pass=0;//{printf("ffr\n");pass=0;}
		if (nr->shwfit.LongE<-1000)pass=0;// {printf("nle\n");pass=0;}
		if (nr->shwfit.LongE>1000)pass=0;//{printf("nle\n");pass=0;}	
		if(!pass && checkNueVars)continue;//trainweight=0;
		
		
		mstvar_combine = nr->mstvars.e4w+nr->mstvars.o4w;
	
	
	
	
		
		//////////////////////
		////fix up the variables
		pass=1;
		if(pr->particles.longest_s_particle_s<0 || pr->particles.longest_s_particle_s>6)pass=0;//{printf("plsp\n");pass=0;}
	//"@particles_elec_vise, "
		//ok particles_primary_phi
		//ok particles_mol_rad_r, "
		if(pr->particles.longest_z<0 || pr->particles.longest_z>6)pass=0;//{printf("plz\n");pass=0;}

	//ok"@particles_emfrac, "
		if(pr->particles.ntot<0 || pr->particles.ntot>50)pass=0;//{printf("pn\n");pass=0;}

	//"@particles_total_long_e, "
		//ok "@particles_weighted_phi, "
		//ok "@largest_frac, "
	//"@reco_frac, "
		if(pr->particles.rms_r<0 || pr->particles.rms_r>100)pass=0;
	//"@particles_prim_par_a, "
		//ok "@particles_prim_par_b, "
		if(pr->particles.prim_par_e0<0 || pr->particles.prim_par_e0>40e3)pass=0;//{printf("ppe0\n");pass=0;}
		if(pr->particles.prim_par_chisq<0 || pr->particles.prim_par_chisq>1000)pass=0;//{printf("ppc\n");pass=0;}
		if(pr->particles.largest_particle_peakdiff<-200 || pr->particles.largest_particle_peakdiff>200)pass=0;//{printf("plpp\n");pass=0;}
		
		

 		//ok "@largest_cmp_chisqndf, "
 		//ok "@particles_frac_particle_2, "
// 	"@particles_primary_long_e, "
//	"@particles_longest_s_particle_e, "
//	"@ntot_lsps, "
//	"@particles_longest_s_particle_par_a, "	
//	"@particles_longest_s_particle_par_b, "	
//	"@particles_longest_s_particle_par_e0, "	
//	"@particles_longest_s_particle_par_chisq, "
		//"@longest_ae0, "	
		//"@prim_ae0 "

//	"@particles_prim_cmp_chisq, "
//	"@particles_prim_cmp_ndf, "
//	"@particles_prim_peakdiff, "
//	"@particles_largest_particle_cmp_chisq, "
//	"@particles_largest_particle_cmp_ndf, "
//	"@particles_longest_s_particle_cmp_chisq, "
//	"@particles_longest_s_particle_cmp_ndf, "
//	"@particles_longest_s_particle_peakdiff "
	if(length_z<0 || length_z>6)pass=0;//{printf("lz\n");pass=0;}

	
	
	//ok "@shwfit_par_a, "
	if(nr->shwfit.par_b<0 || nr->shwfit.par_b>6)pass=0;//{printf("fspb\n");pass=0;}
	//ok "@shwfit_uv_molrad_peak_9s_2pe_dw, "
	//ok "@shwfit_uv_rms_9s_2pe_dw, "
	//ok "@mstvar_combine, "
	if(mstvar_combine<-1000 || nr->shwfit.par_b>1000)pass=0;//{printf("mvc\n");pass=0;}
	//ok "@fracvars_fract_2_planes, "
	//ok "@fracvars_fract_4_planes, "
	//ok "@fracvars_fract_6_planes, "
	//ok "@fracvars_fract_8_counters, "
	//ok "@fracvars_fract_road, "
	if(nr->shwfit.LongE<0 || nr->shwfit.LongE>1200)pass=0;//{printf("sle\n");pass=0;}

		if(!pass && checkParticleVars){failedpval++;continue;}



		
		
		isdis=0;
		if(nr->mctrue.resonanceCode==1003)isdis=1;
		isndis=isdis?0:1;
		
		
		type=1;
		mctrue_type=nr->mctrue.fNueClass;
		if(mctrue_type!=2)type=0;
		
		isnue=mctrue_type==2;
		isnc=mctrue_type==0;
		iscc=mctrue_type==1;
		istau=mctrue_type==3;
		isbeamve=mctrue_type==4;
		



		largest_cmp_chisqndf = pr->particles.largest_particle_cmp_ndf ? pr->particles.largest_particle_cmp_chisq / pr->particles.largest_particle_cmp_ndf : 0;

		if(largest_cmp_chisqndf>3)largest_cmp_chisqndf=3;

		




		for(int i=0;i<11;i++)
			pars[i]=0;
		

		
		int z=0;	


/*		pars[z++]=nr->shwfit.par_a;
		pars[z++]=nr->shwfit.par_b;
		pars[z++]=nr->shwfit.uv_molrad_peak_9s_2pe_dw;
		pars[z++]=nr->shwfit.uv_rms_9s_2pe_dw;
		pars[z++]=mstvar_combine;
		pars[z++]=nr->fracvars.fract_2_planes;
		pars[z++]=nr->fracvars.fract_4_planes;
		pars[z++]=nr->fracvars.fract_6_planes;
		pars[z++]=nr->fracvars.fract_8_counters;
		pars[z++]=nr->fracvars.fract_road;
		pars[z++]=nr->shwfit.LongE;
*/


		pars[z++]=pr->particles.longest_s_particle_s;
		//pars[z++]=pr->particles.primary_phi;
		pars[z++]=pr->particles.mol_rad_r;
		//pars[z++]=pr->particles.longest_z;
		pars[z++]=pr->particles.emfrac;
		pars[z++]=pr->particles.ntot;
		pars[z++]=pr->particles.weighted_phi;
		pars[z++]=largest_frac;

		//pars[z++]=pr->particles.rms_r;

		pars[z++]=pr->particles.prim_par_b;
		pars[z++]=pr->particles.prim_par_e0;
		pars[z++]=pr->particles.prim_par_chisq;
		//pars[z++]=pr->particles.largest_particle_peakdiff;
		pars[z++]=largest_cmp_chisqndf;

		//pars[z++]=pr->particles.frac_particle_2;
		//pars[z++]=length_z;
                pars[z++]=pr->particles.prim_par_a;
		//pars[z++]=pr->event.nstrips;
		pars[z++]=pr->event.nclusters;
		pars[z++]=prim_ae0;
		//pars[z++]=ntot_lsps;


		//pars[z++]=nr->shwfit.uv_molrad_peak_9s_2pe_dw;
		//pars[z++]=nr->shwfit.uv_rms_9s_2pe_dw;







		npars=z;

		
		if(makeTestTree==1)
		{
			
			if(mlp)pid=mlp->Evaluate(0,pars);
		}
		
		sim->Fill();
	}


	printf("failed %d on dq\n",faileddq);
	printf("failed %d on pvals\n",failedpval);
	printf("failed %d on particle pre\n",failedpp);
	printf("failed %d on strict particle pre\n",failedspp);

	sim->Write("TrainTree", TObject::kOverwrite);
	TTree * descriptor = new TTree("descriptor","descriptor");
	descriptor->Branch("npars",&npars);
	descriptor->Fill();
	descriptor->Write("descriptor", TObject::kOverwrite);


	fout->Close();
	return sim;
}



void NNTrain::Train(int steps, int update, int method, string form,double ncemcut,double veemcut)
{
	double pots=6.5*1e8*99;

	//attempt to load previously made tree!
	this->ncemcut=ncemcut;
	this->veemcut=veemcut;

	TTree *simold=0;
	TTree *descriptorin=0;
	TFile * fin = TFile::Open(inputFile.c_str());
	if(fin)simold=(TTree*)fin->Get("TrainTree");
	if(fin)descriptorin=(TTree*)fin->Get("descriptor");
	
	
	
	if(!fin)
	{
		simold = MakeTrainTree();
		TFile * fin = TFile::Open(inputFile.c_str());
		if(fin)simold=(TTree*)fin->Get("TrainTree");
	}



	if(!simold)
	{
		printf("sim tree missing\n");
		exit(1);
	}

	if(!descriptorin)
	{
		printf("descriptor tree missing\n");
		exit(1);
	}

    TFile * f = new TFile("nntemp1.root","RECREATE");
    f->cd();
    
    
        
    TTree * sim=0;    
    TTree * descriptor = descriptorin->CloneTree();
    
    if(needToWeight)
    {
    	printf("need to reweight this tree...\n");
		sim=simold->CloneTree(0);
        SetBranches(sim);
        SetBranches(simold);
        
        int ent=simold->GetEntries();
        printf("looking over %d entries...\n",ent);
        int sav=0;
        int stype[5];
        double sweight[5];
        for(int i=0;i<5;i++){stype[i]=0;sweight[i]=0;}
        
       	//determine max osc*totbeam for each type if needed
		double type_weight[5];
		for(int i=0;i<5;i++)type_weight[i]=0;

		for(int i=0;i<ent;i++)
		{
			simold->GetEntry(i);
			
			if(mctrue_type<0||mctrue_type>4)continue;
			
			
			double w = mctrue_oscprob*mctrue_totbeamweight;
			if(w>type_weight[mctrue_type])type_weight[mctrue_type]=w;
		}
		for(int i=0;i<5;i++)printf("type_weight %d %f\n",i,type_weight[i]);

        
        for(int i=0;i<ent;i++)
        {
        	simold->GetEntry(i);
        	

			
        //	double nueweight=1.3;
        //	if(mctrue_type==2)tweight= tweight*20.;
		
			if(  tweight<gRandom->Uniform()/**type_weight[mctrue_type]*/)continue;
        	
        	sim->Fill();
        	sav++;
        	stype[mctrue_type]++;
        	sweight[mctrue_type]+=tweight;

        }
        printf("using %d entries...\n",sav);
        for(int i=0;i<5;i++)
        {
        	printf("type: %d  cnt: %d  sumweight %f\n",i,stype[i],sweight[i]);
        }
	}else{
		sim=simold->CloneTree(-1,"fast");
        SetBranches(sim);
	}


  	char structure[1000];

	int npars = (int)descriptor->GetLeaf("npars")->GetValue();

	sprintf(structure,"%s","");
	for(int i=0;i<npars-1;i++)sprintf(structure,"%s@pars[%d], ",structure,i);
	sprintf(structure,"%s@pars[%d] ",structure,npars-1);
 	sprintf(structure,"%s:%s:type",structure,form.c_str()); //isnc,isnue,iscc,istau,isbeamve,isndis,isdis
 
 if(mlp)delete mlp;mlp=0;
 mlp = new TMultiLayerPerceptron(structure,sim,"(Entry$+1)%2","(Entry$)%2");
 
 
 if(method==1)mlp->SetLearningMethod(TMultiLayerPerceptron::kSteepestDescent);
  if(method==2) mlp->SetLearningMethod(TMultiLayerPerceptron::kStochastic);
    if(method==3) mlp->SetLearningMethod(TMultiLayerPerceptron::kBatch);

 if(method==4)mlp->SetLearningMethod(TMultiLayerPerceptron::kRibierePolak);
 if(method==5)mlp->SetLearningMethod(TMultiLayerPerceptron::kFletcherReeves);
   if(method==6)mlp->SetLearningMethod(TMultiLayerPerceptron::kBFGS);
   
   
   
 //  if(delta)mlp->SetDelta(delta);
   //if(epsilon)mlp->SetEpsilon(epsilon);
 //  if(eta)mlp->SetEta(eta);
 //  if(etadecay)mlp->SetEtaDecay(etadecay);
 //  if(tau)mlp->SetTau(tau);
   
   
 //  mlp->SetEventWeight("trueEMFrac");
 //mlp->SetEventWeight("emfrac");
   
	//mlp->SetEventWeight("trainweight");
   

   	char tmpa[200];
   	#ifdef dodraw
   		sprintf(tmpa,"+,text,graph,update=%d",update);
    #else
   		sprintf(tmpa,"+,text,update=%d",update);    
    #endif
    
    
	if(!traincontinue)
	{
		mlp->Randomize();
		mlp->DumpWeights("ini.out");
	}
	mlp->LoadWeights("ini.out");

    
    mlp->Train(steps, tmpa);

	mlp->DumpWeights("mlp.out");
    
    
    printf("training done\n");
    
   //mlp->Export("test","python");
   // Use TMLPAnalyzer to see what it looks for

	#ifdef dodraw
   TCanvas* mlpa_canvas = new TCanvas("mlpa_canvas","Network analysis");
   mlpa_canvas->Divide(2,2);
   TMLPAnalyzer ana(mlp);
   // Initialisation
   ana.GatherInformations();
   // output to the console
   ana.CheckNetwork();
   mlpa_canvas->cd(1);
   // shows how each variable influences the network
   ana.DrawDInputs();
   mlpa_canvas->cd(2);
   // shows the network structure
   mlp->Draw();
   mlpa_canvas->cd(3);
   // draws the resulting network
   ana.DrawNetwork(0,"type==1","type==0");
   mlpa_canvas->cd(4);
   #endif
 



	double wantpots= 3.25;  //in e20
	double potscale = (wantpots*1e8)/pots;

printf("b\n");

	int npidbins=250;

  TH1F *bg = new TH1F("bgh", "NN output", npidbins, -.5, 1.5);
   TH1F *sign = new TH1F("sigh", "NN output", npidbins, -.5, 1.5);

  TH2F *bg_recoE = new TH2F("bg_recoE", "bg_recoE", npidbins, -.5, 1.5, 20,0,10);
   TH2F *sig_recoE = new TH2F("sig_recoE", "sig_recoE", npidbins, -.5, 1.5, 20,0,10);
   
   
   bg->SetDirectory(0);
   sign->SetDirectory(0);
   //Double_t params[50];
  
  	Double_t 		 pid=0;
  	Double_t 		 oscweight=0;
 /* 	Double_t 		 pid0=0;
  	Double_t 		 pid1=0;
  	Double_t 		 pid2=0;
  	Double_t 		 pid3=0;
  	Double_t 		 pid4=0;
	Double_t 		 ndis=0;
  	Double_t 		 dis=0;	
   */
   printf("storing %d entries\n",(int)sim->GetEntries());
   TBranch * pidbranch = sim->Branch("pid",&pid);
   TBranch * oscweightbranch = sim->Branch("oscweight",&oscweight);   
  /* TBranch * pid0branch = sim->Branch("pid0",&pid0);
   TBranch * pid1branch = sim->Branch("pid1",&pid1);
   TBranch * pid2branch = sim->Branch("pid2",&pid2);    
   TBranch * pid3branch = sim->Branch("pid3",&pid3);
   TBranch * pid4branch = sim->Branch("pid4",&pid4);
	TBranch * ndisbranch = sim->Branch("ndis",&ndis);
   TBranch * disbranch = sim->Branch("dis",&dis); 
   */
   
   for (int i = 0; i < sim->GetEntries(); i++) {
   
       sim->GetEntry(i);


      pid=mlp->Evaluate(0, pars);


      if(type==0)bg->Fill(pid,mctrue_oscprob * mctrue_totbeamweight * weight);
      if(type==1)sign->Fill(pid,mctrue_oscprob * mctrue_totbeamweight * weight);
      
     if(type==0) bg_recoE->Fill(pid,event_visenergy/25.,mctrue_oscprob * mctrue_totbeamweight * weight);
     if(type==1) sig_recoE->Fill(pid,event_visenergy/25.,mctrue_oscprob * mctrue_totbeamweight * weight);  
         
     oscweight=mctrue_oscprob * mctrue_totbeamweight * potscale;    
	oscweightbranch->Fill();
	
	
	
	pidbranch->Fill();
	 }
 

       
      
   
   sim->Write("", TObject::kOverwrite);
   

   bg->SetLineColor(kBlue);
   bg->SetFillStyle(3008);   bg->SetFillColor(kBlue);
   sign->SetLineColor(kRed);
   sign->SetFillStyle(3003); sign->SetFillColor(kRed);
   bg->SetStats(0);
   sign->SetStats(0);
  #ifdef dodraw
   bg->Draw();
   sign->Draw("same");
   TLegend *legend = new TLegend(.75, .80, .95, .95);
   legend->AddEntry(bg, "Background");
   legend->AddEntry(sign, "Signal");
   legend->Draw();
   mlpa_canvas->cd(0);
#endif




	sign->Write();
	bg->Write();
	

	char tmp[200];
	sprintf(tmp,"fom @ %2.2f pots",wantpots);	
	
	TH1F * hz = new TH1F("fom",tmp,npidbins,-0.5,1.5);
	for(int i=0;i<npidbins;i++)
	{
		 double d = bg->Integral(i,npidbins);
		if(d>0)
		{
			 double n = sign->Integral(i,npidbins);
			hz->SetBinContent(i,n/sqrt(d));
		
		}
	
	}
	hz->Write();


	double bg_sys_frac=0.1;
	double sig_sys_frac=0.1;

	sprintf(tmp,"multibin superfom @ %2.2f pots",wantpots);	
	TH1F * hz1 = new TH1F("fom_mb_super",tmp,npidbins,-0.5,1.5);
	for(int i=0;i<npidbins;i++)
	{
		double d = bg->Integral(i,npidbins);
		if(d<=0)continue;
		{
			 //double n = sign->Integral(i,npidbins);
			 
		double sum_sigbg=0;
		double sum_sig=0;
		double sum_bg=0;
		double sum_sigsig=0;
		
		TH1D *sE = sig_recoE->ProjectionY("sE",i,npidbins);
		TH1D *bE = bg_recoE->ProjectionY("bE",i,npidbins);
		
		for(int k=1;k<sE->GetNbinsX()+1;k++)
		{
			double sig=sE->GetBinContent(k);
			double nbg=bE->GetBinContent(k);
			sum_sig+=sig;
			sum_bg+=nbg;
			sum_sigbg+=nbg*nbg*bg_sys_frac*bg_sys_frac;
			sum_sigsig+=sig*sig*sig_sys_frac*sig_sys_frac;
		}
		
		double orig_sum_fom1=0;
		if(sum_bg)orig_sum_fom1 = sum_sig/sqrt(sum_sig+sum_bg+sum_sigbg+sum_sigsig);
			 
			 
			 
			hz1->SetBinContent(i,orig_sum_fom1);
		
		}
	
	}
	hz1->Write();
	


	sprintf(tmp,"fom @ %2.2f pots",wantpots);	
	
	TH1F * hz2 = new TH1F("superfom",tmp,npidbins,-0.5,1.5);
	for(int i=0;i<npidbins;i++)
	{
		 double d = bg->Integral(i,npidbins);
		if(d>0)
		{
			 double n = sign->Integral(i,npidbins);
			hz2->SetBinContent(i,!d ? 0 : n/sqrt(n+d+d*d*bg_sys_frac*bg_sys_frac+n*n*sig_sys_frac*sig_sys_frac));
		
		}
	
	}
	hz2->Write();


	sprintf(tmp,"fom @ %2.2f pots",wantpots);	
	
	TH1F *hz2a = new TH1F("simplesuperfom",tmp,npidbins,-0.5,1.5);
	for(int i=0;i<npidbins;i++)
	{
		 double d = bg->Integral(i,npidbins);
		if(d>0)
		{
			 double n = sign->Integral(i,npidbins);
			hz2a->SetBinContent(i,!d ? 0 : n/sqrt(d+d*d*bg_sys_frac*bg_sys_frac));
		
		}
	
	}
	hz2a->Write();



	
	
	int mb = hz->GetMaximumBin();
	
	printf("Max fom at %2.2f pots is %f sig %f back %f cut above %f\n",wantpots,hz->GetMaximum(),sign->Integral(mb,npidbins),bg->Integral(mb,npidbins),hz->GetBinLowEdge(mb));

	int mb1 = hz1->GetMaximumBin();
	
	printf("Max multibin super fom at %2.2f pots is %f sig %f back %f cut above %f\n",wantpots,hz1->GetMaximum(),sign->Integral(mb1,npidbins),bg->Integral(mb1,npidbins),hz1->GetBinLowEdge(mb1));

	int mb2 = hz2->GetMaximumBin();
	
	printf("Max super fom at %2.2f pots is %f sig %f back %f cut above %f\n",wantpots,hz2->GetMaximum(),sign->Integral(mb2,npidbins),bg->Integral(mb2,npidbins),hz2->GetBinLowEdge(mb2));

	printf("Max standard super fom at %2.2f pots is %f sig %f back %f cut above %f\n",wantpots,hz2a->GetMaximum(),sign->Integral(mb2,npidbins),bg->Integral(mb2,npidbins),hz2a->GetBinLowEdge(mb2));


	TFile * fmlp = TFile::Open("MLP.root","RECREATE");
	fmlp->cd();
	mlp->Write("MLP");
	

	#ifdef dodraw
		mlpa_canvas->SaveAs("trainImage.eps");
	#endif	

		
	fmlp->Close();
if(fin)fin->Close();
f->Close();



}

void NNTrain::SetBranches(TTree *sim)
{
	sim->SetMakeClass(1);
	sim->SetBranchStatus("*",1);
	sim->SetBranchAddress("weight",&weight);	
	sim->SetBranchAddress("type",&type);
	sim->SetBranchAddress("pars",&pars);
	sim->SetBranchAddress("mctrue_oscprob",&mctrue_oscprob);
	sim->SetBranchAddress("mctrue_fNueClass",&mctrue_type);
	sim->SetBranchAddress("mctrue_totbeamweight",&mctrue_totbeamweight);
	sim->SetBranchAddress("mctrue_nuEnergy",&mctrue_nuenergy);
	sim->SetBranchAddress("visenergy",&event_visenergy);
	sim->SetBranchAddress("resonanceCode",&mctrue_iresonance);
	sim->SetBranchAddress("tweight",&tweight); 

}


void NNTrain::SetOscParamBase( float dm2, float ss13,
                                    float delta, int hierarchy){

  Double_t dm2_12 = fDeltaMS12*1e-5; //best fit SNO
  Double_t dm2_23 = dm2;
                                                                                
  Double_t par[9] = {0};
  par[OscPar::kL] = fBaseLine;
  par[OscPar::kTh23] = fTh23;
  par[OscPar::kTh12] = fTh12;
  par[OscPar::kTh13] = ss13; // TMath::ASin(TMath::Sqrt(ss2th13))/2.;
  par[OscPar::kDeltaM23] = hierarchy*dm2_23;
  par[OscPar::kDeltaM12] = dm2_12;
  par[OscPar::kDensity] = fDensity; //standard rock density
  par[OscPar::kDelta] = delta;
  par[OscPar::kNuAntiNu] = 1;
                                                                                
//  std::cout<<"About to call "<<dm2<<"  "<<ss13<<"  "<<delta<<std::endl;
  fOscCalc.SetOscParam(par);
}




double NNTrain::osc(double nuEnergy, int interactionType, int nonOscFlavor, int oscFlavor)
{



	if(nuEnergy<0)nuEnergy=-nuEnergy;
	if(oscFlavor<0)oscFlavor=-oscFlavor;
	if(nonOscFlavor<0)nonOscFlavor=-nonOscFlavor;

    if (interactionType==0){//NC
      return OscillationProb(f2,0,nuEnergy);
    }
    if (nonOscFlavor==14&&oscFlavor==14){//numu->numu
      return OscillationProb(f2,1,nuEnergy);
    }
    if (nonOscFlavor==12&&oscFlavor==12){//nue->nue
      return OscillationProb(f2,4,nuEnergy);

    }
    if (nonOscFlavor==12&&oscFlavor==14){//nue->numu
      return OscillationProb(f2,6,nuEnergy);
    }
    if (nonOscFlavor==14&&oscFlavor==12){//numu->nue
      return OscillationProb(f2,2,nuEnergy);///0.075;;//for training?
    }
    if (nonOscFlavor==14&&oscFlavor==16){//numu->nutau
      return OscillationProb(f2,3,nuEnergy);

    }
    if (nonOscFlavor==12&&oscFlavor==16){//nue->nutau
      return OscillationProb(f2,5,nuEnergy);

    }
    
  //  printf("? %d %d %f\n",nonOscFlavor,oscFlavor,nuEnergy);
  return 0;  
}    


float NNTrain::OscillationProb(TF1* f2, int ntype, float NuE, float sinth23, float sin2th13) {

  // sinth23 : sine square theta 23
  // sin2th13 : sine square 2 theta 13
  //
  float OscProb = 0 ;
  float NumuToNutau ;
  float NumuToNue ;
  float NueSurvival ;
  float NumuSurvival ;
  float NueToNutau  ;
  float NueToNumu;
  
  if(NuE<0)NuE=-NuE;

  if (ntype==0)
  {
  	OscProb = 1 ;
	return OscProb;
  }




	if (ntype==4)
	{
		  NueSurvival = 1.- sin2th13*f2->Eval(NuE) ;
		OscProb = NueSurvival ;
		return OscProb;
	}
	
	if (ntype==5)
	{
		  NueToNutau = (1.-sinth23)*sin2th13*f2->Eval(NuE) ;
		OscProb = NueToNutau ;
		return OscProb;
	}	
	  
	
  NumuToNue = sinth23*sin2th13*f2->Eval(NuE) ;
  NueToNumu = NumuToNue;
  if (ntype==6){ OscProb = NueToNumu; return OscProb;}
  if (ntype==2){ OscProb = NumuToNue; return OscProb;}
  
  
  NumuToNutau = 4.*sinth23*(1.-sinth23)*pow(1-sin2th13/4,2) ;
  NumuToNutau *= f2->Eval(NuE) ;

if (ntype==3){ OscProb = NumuToNutau; return OscProb;}



	if (ntype==1 )
	{
	  NumuSurvival = 1. - NumuToNutau - NumuToNue ;
		OscProb = NumuSurvival ;  
  	
  		return OscProb;
  	}


  //NueSurvival = 1.;

  //NumuSurvival = 1.;




  
/*
 // if (ntype==0) OscProb = 1 ;
 // else if (ntype==1) OscProb = NumuSurvival ;
 // else if (ntype==2) OscProb = NumuToNue ;
 // else if (ntype==3) OscProb = NumuToNutau ;
//  else if (ntype==4) OscProb = NueSurvival ;
//else if (ntype==5) OscProb = NueToNutau ;
 // else if (ntype==6) OscProb = NueToNumu;
*/
//  cout<<"Period = "<<f2->Eval(fabs(NuE))<<endl ;
//  cout<<"Oscillation probability = "<<OscProb<<endl ;

  return OscProb ;
}


