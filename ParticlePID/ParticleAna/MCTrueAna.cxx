#include "NueAna/ParticlePID/ParticleAna/MCTrueAna.h"
#include "TMath.h"

TF1 * MCTrueAna::osceq=0;
SKZPWeightCalculator * MCTrueAna::skzpCalc=0;


MCTrueAna::MCTrueAna()
{
	osc_dm2=0.0024;
 	osc_L=735.; 
	osc_sinth23=sin(3.1415926535/4.0);
	osc_sin2th13=0.15;


	if(!osceq)
 	{
 		osceq = new TF1("f2","sin(1.267*[0]*[1]/x)*sin(1.267*[0]*[1]/x)",0.,120.);

 		osceq->SetParameters(osc_dm2,osc_L);
	}
	
	if(!skzpCalc)skzpCalc=new SKZPWeightCalculator("DetXs",true);
}

MCTrueAna::~MCTrueAna()
{
}


void MCTrueAna::ana(ParticleObjectHolder * poh, MCTrue * e, ParticleBeamMon * bmon)
{

	////fluxweights/////
	int ismc=bmon->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC;
	if(!ismc)return; //can't do anything with data!

	VldContext evt_vldc = poh->GetHeader().GetVldContext();
	int det=evt_vldc.GetDetector(); 

	e->vtx_u=poh->mctrue.vtx_u;
	e->vtx_v=poh->mctrue.vtx_v;
	e->vtx_z=poh->mctrue.vtx_z;
			
	e->inu=poh->mctrue.inu;
	e->iresonance=poh->mctrue.iresonance;
	e->iaction=poh->mctrue.iaction;
	e->inunoosc=poh->mctrue.inunoosc;
		
	e->nuenergy=poh->mctrue.nuenergy;
//	e->visenergy=poh->mctrue.visenergy;
			
	


	
	/////oscillation/////////
	//only in far!

	if(det==Detector::kFar)
	{
	int otype=0;
	if(e->iaction==1)
	{
		if(abs(e->inu)==12)
		{
			if(abs(e->inunoosc)==12) otype=4;
			else otype=2;
		}	
		if(abs(e->inu)==14)
		{
			if(abs(e->inunoosc)==12) otype=6;
			else otype=1;
		}	
		if(abs(e->inu)==16)
		{
			if(abs(e->inunoosc)==12) otype=5;
			if(abs(e->inunoosc)==14) otype=3;
		}	
	}

	e->oscprob=OscillationProb(osceq, otype, e->nuenergy, osc_sinth23, osc_sin2th13); 
	e->osc_L=osc_L;
	e->osc_dm2=osc_dm2;
	e->osc_sinth23=osc_sinth23;
	e->osc_sin2th13=osc_sin2th13;
	}
	///////
	
	
	e->type = e->iaction * (abs(e->inu)/2-5);
	if(e->type==1 && (e->inu == e ->inunoosc)) e->type=4; //bnue
	


	
	

   	BeamType::BeamType_t beam=bmon->beamtype;
    NtpMCFluxInfo fi= poh->mctrue.flux;
  	double pt = sqrt(fi.tpx*fi.tpx+fi.tpy*fi.tpy);
	double pz = 1.*fi.tpz;
	int tptype = fi.tptype;
 	int zbeam = BeamType::ToZarko(beam);
	e->totbeamweight = skzpCalc->GetBeamWeight(det,zbeam,tptype,pt,pz,e->nuenergy,e->inu);

	//////

	//e->trainweight = e->oscprob * e->totbeamweight ;



}



	    
		


//Oscillation code from Tingjun  7/15/08
double MCTrueAna::OscillationProb(TF1* f2, int ntype, double NuE, double sinth23, double sin2th13) {

 // sinth23 : sine square theta 23
 // sin2th13 : sine square 2 theta 13
 //
 double OscProb = 0 ;
 double NumuToNutau ;
 double NumuToNue ;
 double NueSurvival ;
 double NumuSurvival ;
 double NueToNutau  ;
 double NueToNumu;

 NumuToNutau = 4.*sinth23*(1.-sinth23)*pow(1-sin2th13/4,2) ;
 NumuToNutau *= osceq->Eval(TMath::Abs(NuE)) ;

 NumuToNue = sinth23*sin2th13*f2->Eval(TMath::Abs(NuE)) ;

 NueSurvival = 1.- sin2th13*f2->Eval(TMath::Abs(NuE)) ;
 //NueSurvival = 1.;

 NumuSurvival = 1. - NumuToNutau - NumuToNue ;
 //NumuSurvival = 1.;

 NueToNutau = (1.-sinth23)*sin2th13*f2->Eval(TMath::Abs(NuE)) ;

 NueToNumu = NumuToNue;

 if (ntype==0) OscProb = 1 ;
 else if (ntype==1) OscProb = NumuSurvival ;
 else if (ntype==2) OscProb = NumuToNue ;
 else if (ntype==3) OscProb = NumuToNutau ;
 else if (ntype==4) OscProb = NueSurvival ;
 else if (ntype==5) OscProb = NueToNutau ;
 else if (ntype==6) OscProb = NueToNumu;
//  cout<<"Period = "<<f2->Eval(TMath::Abs(NuE))<<endl ;
//  cout<<"Oscillation probability = "<<OscProb<<endl ;

 return OscProb ;
} 

