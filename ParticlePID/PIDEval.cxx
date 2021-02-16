#include "NueAna/ParticlePID/PIDEval.h"

#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "Conventions/Detector.h"

#include "TTree.h"



JOBMODULE(PIDEval, "PIDEval",
          "does the pid eval for ParticleFinder");
CVSID("$Id: PIDEval.cxx,v 1.9 2011/02/17 21:15:15 whitehd Exp $");


TMultiLayerPerceptron * PIDEval::mlp1=0;
TMultiLayerPerceptron * PIDEval::mlp2=0;
TMultiLayerPerceptron * PIDEval::mlp3=0;
TMultiLayerPerceptron * PIDEval::mlp4=0;
TMultiLayerPerceptron * PIDEval::mlp5=0;
TMultiLayerPerceptron * PIDEval::mlp6=0;

int PIDEval::first=1;

PIDEval::PIDEval()
{
	if(mlp1==0 && first) //only try once...
    {
		LoadMLP();
		first=0;
    }
    
	infile="";
	outfile="";
	cnt=0;
	total=1;
	
	
	//osc

/*	fBaseLine = 735;
  	fDeltaMS12= 8.0;
	fTh12= 0.816;
	fTh23=3.141592/4.0;
	fDensity=2.75;
	SetOscParamBase(0.00243, TMath::ASin(TMath::Sqrt(0.15))/2.0, 0, 1);	
  */		   
  
}

PIDEval::~PIDEval()
{}

void PIDEval::LoadMLP()
{
  TFile *mlpfile0 = TFile::Open("/minos/app/nue/Releases/Griffin/Data/ParticlePID/finalMLPs/MLP_13inp_20:9.root");
  if (mlpfile0==0)
    {
      printf("Error in PIDEval::LoadMLP() 1-- Bad file name\n");
    }
  if(mlpfile0)mlp1 = (TMultiLayerPerceptron*) mlpfile0->Get("MLP");
  if (mlp1==0)
    {
      printf("Error in PIDEval::LoadMLP() 1-- File does not contain mlp\n");
    }else{
		printf("loaded ParticlePID from %s\n",mlpfile0->GetName());
  		mlp1=(TMultiLayerPerceptron*)mlp1->Clone();
  		mlpfile0->Close();
  	}

	TFile *mlpfile1 = TFile::Open("/minos/app/nue/Releases/Griffin/Data/ParticlePID/finalMLPs/MLP_13inp_25:10.root");
  if (mlpfile1==0)
    {
      printf("Error in PIDEval::LoadMLP() 2-- Bad file name\n");
      
    }
  if(mlpfile1)mlp2 = (TMultiLayerPerceptron*) mlpfile1->Get("MLP");
  if (mlp2==0)
    {
      printf("Error in PIDEval::LoadMLP() 2-- File does not contain mlp\n");
    }else{
                printf("loaded ParticlePID from %s\n",mlpfile1->GetName());
  		mlp2=(TMultiLayerPerceptron*)mlp2->Clone();
  		mlpfile1->Close();
  	}
  
  
	TFile *mlpfile2 = TFile::Open("/minos/app/nue/Releases/Griffin/Data/ParticlePID/finalMLPs/MLP_14inp_16:9.root");
  if (mlpfile2==0)
    {
      printf("Error in PIDEval::LoadMLP() 3-- Bad file name\n");
      
    }
  if(mlpfile2)mlp3 = (TMultiLayerPerceptron*) mlpfile2->Get("MLP");
  if (mlp3==0)
    {
      printf("Error in PIDEval::LoadMLP() 3-- File does not contain mlp\n");
    }else{
                printf("loaded ParticlePID from %s\n",mlpfile2->GetName());
  		mlp3=(TMultiLayerPerceptron*)mlp3->Clone();
  		mlpfile2->Close();
  	}  
  
  
  
  
	TFile *mlpfile3 = TFile::Open("/minos/app/nue/Releases/Griffin/Data/ParticlePID/finalMLPs/MLP_14inp_18:11.root");

  if (mlpfile3==0)
    {
      printf("Error in PIDEval::LoadMLP() 4-- Bad file name\n");
      
    }
  if(mlpfile3)mlp4 = (TMultiLayerPerceptron*) mlpfile3->Get("MLP");
  if (mlp4==0)
    {
      printf("Error in PIDEval::LoadMLP() 4-- File does not contain mlp\n");
    }else{
                printf("loaded ParticlePID from %s\n",mlpfile3->GetName());
  		mlp4=(TMultiLayerPerceptron*)mlp4->Clone();
  		mlpfile3->Close();
  	}    
  
 
  
	TFile *mlpfile4 = TFile::Open("/minos/app/nue/Releases/Griffin/Data/ParticlePID/finalMLPs/MLP_14inp_18:6.root");

  if (mlpfile4==0)
    {
      printf("Error in PIDEval::LoadMLP() 5-- Bad file name\n");
      
    }
  if(mlpfile4)mlp5 = (TMultiLayerPerceptron*) mlpfile4->Get("MLP");
  if (mlp5==0)
    {
      printf("Error in PIDEval::LoadMLP() 5-- File does not contain mlp\n");
    }else{
                printf("loaded ParticlePID from %s\n",mlpfile4->GetName());
  		mlp5=(TMultiLayerPerceptron*)mlp5->Clone();
  		mlpfile4->Close();
  	}  
  	
  	
 
  
	TFile *mlpfile5 = TFile::Open("/minos/app/nue/Releases/Griffin/Data/ParticlePID/finalMLPs/MLP_14inp_23:6.root");

  if (mlpfile5==0)
    {
      printf("Error in PIDEval::LoadMLP() 6-- Bad file name\n");
      
    }
  if(mlpfile5)mlp6 = (TMultiLayerPerceptron*) mlpfile5->Get("MLP");
  if (mlp6==0)
    {
      printf("Error in PIDEval::LoadMLP() 6-- File does not contain mlp\n");
    }else{
                printf("loaded ParticlePID from %s\n",mlpfile5->GetName());
  		mlp6=(TMultiLayerPerceptron*)mlp6->Clone();
  		mlpfile5->Close();
  	}    	  


    if(! (mlp1 || mlp2 || mlp3 || mlp4 || mlp5 || mlp6) )
    {           
        printf("ParticlePID files are missing! directory must have moved\nexiting...\n");
        exit(1);
    }


  
  
}

void PIDEval::ana(Particles * particles, Event * event, NueRecord* /*nr*/)
{

	if(cnt%10000==0)printf("%d ... \n",cnt);//%2.2f%%\n",cnt,100.*((double)cnt)/((double)total));
	cnt++;



	event->pidA=-1;
        event->pidB=-1;
        event->pidC=-1;
        event->pidD=-1;
        event->pidE=-1;	
        event->pidF=-1;



	//double reco_frac=particles->prim_vise/event->visenergy;
	//double length_z = event->max_z - event->min_z;	
    
    double largest_frac=particles->totvise ? 
    	particles->largest_particle_e/particles->totvise : 1;
    double prim_ae0=particles->prim_par_e0 ? 
    	particles->prim_par_a/particles->prim_par_e0 : 0;
 	double largest_cmp_chisqndf = particles->largest_particle_cmp_ndf ?
 		particles->largest_particle_cmp_chisq / 
 			particles->largest_particle_cmp_ndf : 0;
 	
 	

	//check ranges
	if(particles->longest_s_particle_s<0 || 
		particles->longest_s_particle_s>6)return;
	//removing restriction on longest_z since it is not in the nn
	//if(particles->longest_z<0 || 
	//	particles->longest_z>6)return;
	if(particles->ntot<0 || 
		particles->ntot>50)return;
	//removing restriction on rms_r since it is not in the nn
	//if(particles->rms_r<0 || 
	//	particles->rms_r>100)return;
	if(particles->prim_par_e0<0 || 
		particles->prim_par_e0>40e3)return;
	//changing max on prim_par_chisq from 1000 to 10000 for flexibility
	if(particles->prim_par_chisq<0 || 
		particles->prim_par_chisq>10000)return;
	if(particles->largest_particle_peakdiff<-200 || 
		particles->largest_particle_peakdiff>200)return;		
	///////

    
    

	double pars[30];//30? Pick numer based on appropriate neural network.
	
	int z=0;
  	
    if(mlp1 || mlp2)
    {

		z=0;
		

		pars[z++]=particles->longest_s_particle_s;
		pars[z++]=particles->mol_rad_r;
		pars[z++]=particles->emfrac;
		pars[z++]=particles->ntot;
		pars[z++]=particles->weighted_phi;
		pars[z++]=largest_frac;
		pars[z++]=particles->prim_par_b;
		pars[z++]=particles->prim_par_e0;
		pars[z++]=particles->prim_par_chisq;
		pars[z++]=largest_cmp_chisqndf;
       		pars[z++]=particles->prim_par_a;
		pars[z++]=event->nclusters;
		pars[z++]=prim_ae0;

 	 	for(int i=z;i<30;i++)pars[i]=0;
  	
 
		event->pidA = mlp1 ? mlp1->Evaluate(0, pars) : -1;

		event->pidB = mlp2 ? mlp2->Evaluate(0, pars) : -1;
 
 	}
 	
 	


  	
    if(mlp3 || mlp4 || mlp5 || mlp6)
    {

		z=0;
		

		pars[z++]=particles->longest_s_particle_s;
		pars[z++]=particles->mol_rad_r;
		pars[z++]=particles->emfrac;
		pars[z++]=particles->ntot;
		pars[z++]=particles->weighted_phi;
		pars[z++]=largest_frac;
		pars[z++]=particles->prim_par_b;
		pars[z++]=particles->prim_par_e0;
		pars[z++]=particles->prim_par_chisq;
		pars[z++]=particles->largest_particle_peakdiff;
		pars[z++]=largest_cmp_chisqndf;
        	pars[z++]=particles->prim_par_a;
		pars[z++]=event->nclusters;
		pars[z++]=prim_ae0;

 	 	for(int i=z;i<30;i++)pars[i]=0;
  	
 
		event->pidC = mlp3 ? mlp3->Evaluate(0, pars) : -1;

		event->pidD = mlp4 ? mlp4->Evaluate(0, pars) : -1;

		event->pidE = mlp5 ? mlp5->Evaluate(0, pars) : -1;

		event->pidF = mlp6 ? mlp6->Evaluate(0, pars) : -1;
		 
 	}
 	
 	 
    
    

}




//......................................................................

void PIDEval::BeginJob()
{
  ///
  /// (Document me!)
  ///

/*	TFile *f = TFile::Open(infile.c_str());
	if(!f)return;
	TTree * pt = (TTree*)f->Get("PA");
	if(pt)
	{
		total=pt->GetEntries();
		printf("%d total PRecords\n",total);
	}
	f->Close();
*/

}

//......................................................................

void PIDEval::EndJob()
{
  ///
  /// (Document me!)
  ///


	//copy pot tree?
/*	if(infile!="" && outfile!="")
	{
		//copy the pot tree

   		TDirectory *savedir = gDirectory;
        
   		TFile *fpf = dynamic_cast<TFile *>(gROOT->GetListOfFiles()->FindObject(outfile.c_str()));
   		if(fpf){

			TFile *f = TFile::Open(infile.c_str());
			TTree * pt = (TTree*)f->Get("pottree");
			fpf->cd();
			if(pt)
				pt->CloneTree()->Write("pottree");
			else
				printf("NO POTTREE FOUND...!\n");
			savedir->cd();
		}
	
	}
*/


  if(kPOTTreeName=="")return; //don't want to copy pottree...
  

   printf("PIDEval copying POTtree from input file\n");  
   TFile *fpf = dynamic_cast<TFile *>(gROOT->GetListOfFiles()->FindObject(kPOTTreeName.c_str()));
   if(fpf){
        printf("got outfile\n");
   }


   TFile *fpi = dynamic_cast<TFile *>(gROOT->GetListOfFiles()->FindObject(kPOTTreeNameIn.c_str()));
   if(fpi){
        printf("got infile\n");
        TTree *pt = (TTree*)fpi->Get("pottree");
        if(pt)
        {
                printf("got pottree\n");
                fpf->cd();
                pt->CloneTree()->Write();
        }
   }




}

//......................................................................

JobCResult PIDEval::Reco(MomNavigator* mom)
{
  ///
  /// (Document me!)
  ///


	std::vector<TObject* > hft =( mom->GetFragmentList("PRecord"));

	for(unsigned int s =0;s<hft.size();s++)
	{
		PRecord *n=0;
		n=dynamic_cast<PRecord *>(hft[s]);
		ana(&(n->particles),&(n->event),0);
		
		//recalc_oscprob(n);
		//adjRecoE(n);
	}

	if(hft.size()<1) //no PRecords.. so look for nuerecords...
	{
		std::vector<TObject* > nrv =( mom->GetFragmentList("NueRecord"));
		NueRecord * nr =0;
		for(unsigned int s=0;s<nrv.size();s++)
		{
			nr = dynamic_cast<NueRecord *>(nrv[s]);		
			PRecord *n=&nr->precord;
			if(n)ana(&(n->particles),&(n->event),nr);
			n=&nr->precordMRCC;
			if(n)ana(&(n->particles),&(n->event),nr);
			
	//		recalc_oscprob(nr);
		}
	}

	return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

//......................................................................

const Registry& PIDEval::DefaultConfig() const
{
  ///
  /// Supply the default configuration for the module
  ///

  static Registry r; // Default configuration for module

  // Set name of config
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());
  r.Set("InFileName",""); //don't set when used with a real reco chain
  r.Set("OutFileName","");
  // Set values in configuration

  r.Set("POTTreeName","");
  r.Set("POTTreeNameIn","");






  r.UnLockValues();
  r.LockValues();


  return r;
}

//......................................................................

void PIDEval::Config(const Registry& r)
{
  ///
  /// Configure the module given the Registry r
  ///
  
 
    const char* tmps;
    if(r.Get("InFileName",tmps)){infile=tmps;}
    if(r.Get("OutFileName",tmps)){outfile=tmps;}  


    if(r.Get("POTTreeName",tmps))kPOTTreeName=tmps;
    if(r.Get("POTTreeNameIn",tmps))kPOTTreeNameIn=tmps;

}

//......................................................................

void PIDEval::Reset()
{
  ///
  /// (Document me!)
  ///
}

////////////////////////////////////////////////////////////////////////

/*
void PIDEval::adjRecoE(PRecord * p)
{
	if(p->GetHeader().GetVldContext().GetDetector()==Detector::kNear)
	{
	
		p->event.visenergy=p->event.visenergy*0.03799819+0.4803261;
		p->particles.totvise=p->particles.totvise*0.03958938+0.4161953;
	
	}else if(p->GetHeader().GetVldContext().GetDetector()==Detector::kFar)
	{
		p->event.visenergy=p->event.visenergy*0.0387301+0.489987;
		p->particles.totvise=p->particles.totvise*0.03379042+0.5866994;
	
	}
	

}
*/

/*
void PIDEval::recalc_oscprob(NueRecord * p)
{
	if(p->GetHeader().GetVldContext().GetDetector()==Detector::kNear)
	{
//		p->mctrue.fOscProb=1.;
		return;
	}


	
	double newoscprob =  fOscCalc.Oscillate(p->mctrue.nuFlavor, p->mctrue.nonOscNuFlavor, p->mctrue.nuEnergy);
	
//	if(p->mctrue.fNueClass==0)newoscprob/=3.; //adjust for equal file types
	
//	printf("%f %f\n",p->mctrue.oscprob,newoscprob);
	p->mctrue.fOscProb = newoscprob;
	
	/ *
	p->mctrue.osc_L=fBaseLine;
	p->mctrue.osc_dm2=fDeltaMS12;
	p->mctrue.osc_sinth23=fTh12;
	p->mctrue.osc_sin2th13=fTh23;
* /
}
*/
/*
void PIDEval::SetOscParamBase( float dm2, float ss13,
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
*/
