
//uncomment next line to draw the fits for debuggins!
#define EMDEBUG 1

#include "MessageService/MsgService.h"

#include "TMinuit.h"
#include "TMath.h"

#include "ShwFit.h"

#ifdef EMDEBUG	
#include "TCanvas.h"
#include "TLatex.h"
#endif


#include <math.h>

#include <iostream>
#include <fstream>


#include "Plex/PlexPlaneId.h"
#include "UgliGeometry/UgliGeomHandle.h"
#include "UgliGeometry/UgliScintPlnHandle.h"
#include "Plex/PlexStripEndId.h"
#include "UgliGeometry/UgliStripHandle.h"

#include "UgliGeometry/UgliSteelPlnHandle.h"



CVSID("$Id: ShwFit.cxx,v 1.5 2009/07/02 18:31:32 scavan Exp $");



TH1F * ShwFit::lenepl=0;
TF1 * ShwFit::efit=0;
TH1F * ShwFit::lenepl3d=0;

std::map<double,int> * ShwFit::planemap=0;


/*
static Double_t shwfunc(Double_t *x, Double_t *par)
{

  // function to be fitted for em showers par[0]=a par[1]=E0
  Double_t Rsteel = 1.445; // for steel
//    Double_t R = 1.46676;
  Double_t Rscint=0.0241; //for scint only

	Double_t R = Rscint+Rsteel;


 Double_t xx=R*x[0];

if(xx==0)return 0;//to fix a bug..... we don't expect anything here anyways..

//	Double_t xx=R*(x[0]+par[3]);

 // cout<<"In shwfunc "<<xx<<" par[0] "<<par[0]<<" par[1] "<<par[1]<<" par[2] "<<par[2]<<endl;

  Double_t lnf = TMath::Log(Rscint*par[2]*par[1])+(par[0]-1)*TMath::Log(xx*par[1])-
     par[1]*xx-TMath::LnGamma(par[0]);
     
 //     Double_t lnf = TMath::Log(par[2]*par[1])+(par[0]-1)*TMath::Log(xx*par[1])-
  //   par[1]*xx-TMath::LnGamma(par[0]);


//	xx+=Rscint;
//	 Double_t lnfh = TMath::Log(par[2]*par[1])+(par[0]-1)*TMath::Log(xx*par[1])-
  //   par[1]*xx-TMath::LnGamma(par[0]);


  Double_t f = exp(lnf);
	
	//Double_t fh=exp(lnfh);
	
	//return (fh+f)/2*Rscint;


//  Double_t f = 1.46676*par[2]*TMath::Power(par[1],par[0])*
//               TMath::Power(xx,par[0]-1)*
//               TMath::Exp(-par[1]*(xx))/TMath::Gamma(par[0]);
  return f;
}
*/

ShwFit::ShwFit()
{
     	debug=0;
       	if(MsgService::Instance()->IsActive("ShwFit",Msg::kDebug))debug=1;


	Reset();
   //const Int_t LHISTBINS=100;
  // if(lenepl)delete lenepl;
   //	 lenepl = new TH1F("lenepl","longitudinal energy by plane",LHISTBINS,0.0,LHISTBINS);
	//lenepl->SetDirectory(0);
	if(!planemap)BuildPlaneMap();



	
	
	if(debug)MSG("ShwFit",Msg::kError)<<"ParticleFinder ShwFit is in debug mode... there will be a purposeful memory leak\n";
	
	
}

void ShwFit::BuildPlaneMap()
{
	//try to open plane location file ... for faster running!
	
	ifstream infile("planemap.txt");
	if(infile.is_open())
	{
		MSG("ShwFit",Msg::kDebug)<<"loading shwfit planemap from file!....  ";
	
		planemap=new std::map<double,int>;
		string line;
		int cnt=0;
		while(!infile.eof())
		{
			getline(infile,line);
			int planei=0;
			double planez=0;
			sscanf(line.c_str(),"%lf %d",&planez, &planei);
			planemap->insert(make_pair(planez,planei));	
			cnt++;
		}
		
		MSG("ShwFit",Msg::kDebug)<<"loaded "<<cnt<<" planes....";
	
		infile.close();
	}else{


		VldContext vldc(Detector::kFar , SimFlag::kData , VldTimeStamp());
		UgliGeomHandle geo(vldc);       

		planemap=new std::map<double,int>;

		//iterate over the planes, and insert the upper z and plane number into the map
		std::vector<UgliSteelPlnHandle>  plns = geo.GetSteelPlnHandleVector();
 
 		ofstream outfile("planemap.txt");
 		
 		
 		for(unsigned int i=0;i<plns.size();i++)
 		{
 			planemap->insert(make_pair(plns[i].GetZ0(),i));
 	
 			if(outfile.is_open())
 			{
 				outfile << (double)plns[i].GetZ0() << " "<<i<<"\n";
 			}
 		}
 		
 		outfile.close();
 	
 	}
 	


}


int ShwFit::GetPlaneFromZ(double z)
{
	int plane=0;
	
	std::map<double,int>::iterator i = planemap->lower_bound(z);
	
	if(i!=planemap->end()) plane=i->second;

	
	return plane;
}


ShwFit::~ShwFit()
{

	if(!debug)
	{
		if(lenepl){delete lenepl;lenepl=0;}
		if(lenepl3d){delete lenepl3d;lenepl3d=0;}
		if(efit){delete efit;efit=0;}
	}

}



void ShwFit::Reset()
{

	if(!debug)
	{
		if(lenepl){delete lenepl;lenepl=0;}
		if(lenepl3d){delete lenepl3d;lenepl3d=0;}
		if(efit){delete efit;efit=0;}
	}

	cmp_chisq=0;
	cmp_ndf = 0;
	peakdiff=0;
		
	pp_chisq=0;
	pp_ndf=0;
	pp_igood=0;
	pp_p=0;
	
	lenepl=0;
	par_a=0;
	par_b=0;
	par_e0=0;
	par_a_err=0;
	par_b_err=0;
	par_e0_err=0;	
	chisq=0;
	ndf=0;
	prob=0;
	shwmax=0;
	shwmaxplane=0;
	conv=0;
	
	
	pred_e_a=0;
	pred_g_a=0;
	pred_b=0;
	pred_e0=0;
	pred_e_chisq=0;
	pred_e_ndf=0;
	pred_g_chisq=0;
	pred_g_ndf=0;	
	

	detector=0;
	
	v_ph.clear();
	v_u.clear();
	v_v.clear();
	v_z.clear();
	
	zshift=0;
	
	pre_over=0;
	pre_under=0;
	post_over=0;
	post_under=0;
	
}

void ShwFit::Insert(double ph, double z)
{
	lenepl->Fill(z+0.5,ph);
}

Double_t ShwFit::GetMaximumX(TF1* efit, Double_t xmin, Double_t xmax)
{                           
   Double_t fXmin = 0.001;
   Double_t fXmax = 100;
   Double_t fNpx = 100;                                                    

   Double_t x,y;
   if (xmin >= xmax) {xmin = fXmin; xmax = fXmax;}
   Double_t dx = (xmax-xmin)/fNpx;
   Double_t xxmax = xmin;
   Double_t yymax = efit->Eval(xmin+dx);
   for (Int_t i=0;i<fNpx;i++) {
      x = xmin + (i+0.5)*dx;
      y = efit->Eval(x);
      if (y > yymax) {xxmax = x; yymax = y;}
   }
   if (dx < 1.e-9*(fXmax-fXmin)) return TMath::Min(xmax,xxmax);
   else return GetMaximumX(efit, TMath::Max(xmin,xxmax-dx), TMath::Min(xmax,xxmax+dx));
}


//function for fitting to 3d info.... histogram is in radlengths!
static Double_t shwfunc3d(Double_t *x, Double_t *par)
{

	// function to be fitted for em showers par[0]=a par[1]=E0
	//Double_t Rsteel = 1.445; // for steel
	//    Double_t R = 1.46676;
	Double_t Rscint=0.0241; //for scint only

	//Double_t R = Rscint+Rsteel;
	// Double_t xx=R*x[0];
	Double_t xx=x[0];//+par[3]; //uncomment for offset fitting

	//MSG("ShwFit",Msg::kWarning)<<"at "<<xx<<endl;

	if(xx<=0 || xx >1e9 )
	{
		//cout<<"badx "<<xx<<"\n";

		//	MSG("ShwFit",Msg::kError)<<"removing hit at "<<xx<<endl;
		return 0;//to fix a bug..... we don't expect anything here anyways..
	}


	//	Double_t xx=R*(x[0]+par[3]);

 	// cout<<"In shwfunc "<<xx<<" par[0] "<<par[0]<<" par[1] "<<par[1]<<" par[2] "<<par[2]<<endl;

	if(TMath::Abs(par[0])<1e-6)return 0;

	Double_t lnf = TMath::Log(Rscint*par[2]*par[1])+(par[0]-1)*TMath::Log(xx*par[1])-
	par[1]*xx-TMath::LnGamma(par[0]);

	
     
	//     Double_t lnf = TMath::Log(par[2]*par[1])+(par[0]-1)*TMath::Log(xx*par[1])-
	//   par[1]*xx-TMath::LnGamma(par[0]);


	//	xx+=Rscint;
	//	 Double_t lnfh = TMath::Log(par[2]*par[1])+(par[0]-1)*TMath::Log(xx*par[1])-
	//   par[1]*xx-TMath::LnGamma(par[0]);

	//cout <<"eval "<<lnf<<" ";
	Double_t f = exp(lnf);
	

	//    cout <<f<<"\n";
 
	//if the number is close to 0... there are FPEs on minos machines, 	probably due to lost precision 
	if(f<1e-10)return 0;



	//Double_t fh=exp(lnfh);
	
	//return (fh+f)/2*Rscint;


	//  Double_t f = 1.46676*par[2]*TMath::Power(par[1],par[0])*
	//               TMath::Power(xx,par[0]-1)*
	//               TMath::Exp(-par[1]*(xx))/TMath::Gamma(par[0]);
	return f;
}


void ShwFit::SetDetector(int mydet)
{
	detector=mydet;
}

void ShwFit::Insert3d(double ph, double u, double v, double z)
{
	v_ph.push_back(ph);
	v_u.push_back(u);
	v_v.push_back(v);
	v_z.push_back(z);
}



void ShwFit::Fit3d(int psf)
{

	double Rsteel = 1.445; // for steel
  	double Rscint=0.0241; //for scint only
	double R = Rscint+Rsteel;
	
	double offset=Rsteel / 2; //offset from front of scint plane to count as 0 radiation lengths...

/*
	//build a vector measuring distance in radiation lengths from first hit
	VldContext vldc(Detector::kFar , SimFlag::kData , VldTimeStamp());
	UgliGeomHandle geo(vldc);       
*/
	std::vector<int> planes;
	
	for(unsigned int i=0;i<v_z.size();i++)
	{
/*		PlexPlaneId p;
		p = geo.GetPlaneIdFromZ(v_z[i]); 	
		planes.push_back(p.GetPlane());
*/
		planes.push_back(GetPlaneFromZ(v_z[i]));	
	}
	

	std::vector<double> rdlen;	
	double rddist=offset;	
	rdlen.push_back(rddist);
		
	int sm1=0;
	int sm2=0;
	
	if(v_z[0]<15.5)sm1=1;
	if(v_z[0]>=15.5)sm2=1;
	
	for(unsigned int i=1;i<v_z.size();i++)
	{
		double dz = v_z[i]-v_z[i-1];
		double du = v_u[i]-v_u[i-1];
		double dv = v_v[i]-v_v[i-1];
	
		double tdist = sqrt(dz*dz +du*du +dv*dv);

		if(planes[i]==planes[i-1]) //in the same scint plane
		{
			rddist+=tdist * Rscint;
		}else{
			int dp = planes[i]-planes[i-1];
			double dt = sqrt(du*du+dv*dv);
			rddist += sqrt( R*R*dp*dp +  dt*dt);
		}
		rdlen.push_back(rddist);
		
		if(v_z[i]<15.5)sm1=1;
		if(v_z[i]>=15.5)sm2=1;	
	}
	
	
	if(sm1 && sm2)
	{
		MSG("ShwFit",Msg::kError)<<"Attempting to fit between super modules! exiting\n";
		return;	
	}
	
	
	//for(unsigned int i=0;i<v_z.size();i++)
		//std::cout << "z " << v_z[i] << " plane " <<planes[i] <<" rdlen "<<rdlen[i] <<endl;	
	

	//decide what size histo to make, and make it	
	int nbins = v_z.size()*2;	
	double rlmax = rdlen[rdlen.size()-1]*2;

	//don't delete the histo if we want to be viewing it later
	//this will leak... so just don't run this module in debug mode
	//except when debugging!
	if(!debug) if (lenepl3d){delete lenepl3d;lenepl3d=0;}
	
	if(debug && psf)
		lenepl3d = new TH1F("lenepl3dPSF","longitudinal energy by plane",nbins,0.0,rlmax);
	else
		lenepl3d = new TH1F("lenepl3d","longitudinal energy by plane",nbins,0.0,rlmax);
	//hide the histo unless we are debugging the primary shower
	if(!debug)	lenepl3d->SetDirectory(0);


	//try a smear fill...		 
	for(unsigned int i=0;i<rdlen.size();i++) 
	{
		double ph = v_ph[i];
		ph/=R;
		//cout<<"adding to hist "<<rdlen[i]<<" e "<<ph<<"\n";	
		if(i==0)
		{
			lenepl3d->Fill(rdlen[i],ph);//*0.8);
			//lenepl3d->Fill(rdlen[i]+lenepl3d->GetBinWidth(1),ph*0.2);
			//lenepl3d->Fill(rdlen[i]+lenepl3d->GetBinWidth(1)*2,ph*0.1);
		}else{
			//lenepl3d->Fill(rdlen[i]+lenepl3d->GetBinWidth(1),ph*0.2);
			lenepl3d->Fill(rdlen[i],ph);//*0.60);
			//lenepl3d->Fill(rdlen[i]-lenepl3d->GetBinWidth(1),ph*0.2);
			
			//lenepl3d->Fill(rdlen[i]+lenepl3d->GetBinWidth(1),ph*0.05);
			//lenepl3d->Fill(rdlen[i]-lenepl3d->GetBinWidth(1),ph*0.05);
		}
	}
	
	///////////////////////
	//calculate expected values here....
		///see if we need a z shift!



 	
 	
	double chisq_elec=-1;
	double chisq_gamma=-1;

	//what should a be for this energy?
	//a = 1 +b(lnE/Ec -0.5)  (or+0.5 for gamma)  b=0.5
	//Ec is about 36.64 MeV assuming steel has z 26 and scint is about z=7

	//calculate total e
	double etot=0;
	for(int i=0;i<(int)v_ph.size();i++)
		etot+=v_ph[i];
	
	//make it MeV with roughly 25 mip/gev

	etot=etot/23*1000;// *2; 

	if(etot<1e-5)return;

	pred_b=0.55;
	pred_e0=etot;


	double a=1+0.5*(log(etot/36.64)-0.5);
	TF1* predElec = 0;
	if(psf)predElec = new TF1("predElecPSF",shwfunc3d,0.1,rlmax,3);
	else predElec = new TF1("predElec",shwfunc3d,0.1,rlmax,3);
   	predElec->SetParNames("a","b","e0");
   	predElec->SetParameter(0,a);
   	predElec->SetParameter(1,0.5); 	
   	predElec->SetParameter(2,etot);
 	predElec->SetLineColor(2);

	pred_e_a=a;
	
	MSG("ShwFit",Msg::kDebug)<<"elec a b e0 "<<a<<" "<<0.5<<" "<<etot<<"\n";


	double maxX = GetMaximumX(predElec,0,rlmax);
	
	int maxbin = lenepl3d->GetMaximumBin();
	double maxH = lenepl3d->GetBinCenter(maxbin);
	
	MSG("ShwFit",Msg::kDebug)<<"pred max "<<maxX<<" hist max "<<maxH<<"\n";
	
/*	
////////////shifting	
	//do we need to reset the offset of the histogram?
	if(fabs(maxH-maxX)>0.05)
	{
		zshift = maxX-maxH;

		//lenepl3d->Reset();
		if(lenepl3d)delete lenepl3d;lenepl3d=0;
		lenepl3d = new TH1F("lenepl3d","longitudinal energy by plane",nbins,0.0,rlmax);
		lenepl3d->SetDirectory(0);
		
		//try a smear fill...		 
		for(unsigned int i=0;i<rdlen.size();i++) 
		{
			double ph = v_ph[i];
			ph/=R;
		
			if(i==0)
			{
				lenepl3d->Fill(rdlen[i]+zshift,ph);// *0.8);
				//lenepl3d->Fill(rdlen[i]+zshift+lenepl3d->GetBinWidth(1),ph*0.2);
				//lenepl3d->Fill(rdlen[i]+zshift+lenepl3d->GetBinWidth(1)*2,ph*0.1);
			}else{
			//	lenepl3d->Fill(rdlen[i]+zshift+lenepl3d->GetBinWidth(1),ph*0.2);
				lenepl3d->Fill(rdlen[i]+zshift,ph);// *0.60);
			//	lenepl3d->Fill(rdlen[i]+zshift-lenepl3d->GetBinWidth(1),ph*0.2);
			
				//lenepl3d->Fill(rdlen[i]+zshift+lenepl3d->GetBinWidth(1),ph*0.05);
				//lenepl3d->Fill(rdlen[i]+zshift-lenepl3d->GetBinWidth(1),ph*0.05);
			}
		}	
		
		//store zshift in terms of actual distance, not radiation length
		zshift/=R; //to get number of planes
		zshift*=0.0354;//average plane width
	
		int maxbin = lenepl3d->GetMaximumBin();
		double maxH = lenepl3d->GetBinCenter(maxbin);
		printf("new hist max %f\n",maxH);
	}
/////////////
*/


	a=1+0.5*(log(etot/36.64)+0.5); //gamma
	TF1* predGamma = 0;
	if(psf) predGamma = new TF1("predGammaPSF",shwfunc3d,0.1,rlmax,3);
	else predGamma = new TF1("predGamma",shwfunc3d,0.1,rlmax,3);
   	predGamma->SetParNames("a","b","e0");
   	predGamma->SetParameter(0,a);
   	predGamma->SetParameter(1,0.5); 	
   	predGamma->SetParameter(2,etot);
 	predGamma->SetLineColor(3);

	pred_g_a=a;
	MSG("ShwFit",Msg::kDebug)<<"gamma a b e0 "<<a<<" "<<0.5<<" "<<etot<<"\n";
  
	//calculate chi sq's ....

  
	chisq_elec=0;
  	chisq_gamma=0;
  	int dc_e=0;
  	int dc_g=0;
  	for(int i=1;i<lenepl3d->GetNbinsX()+1;i++)
  	{
  		double val = lenepl3d->GetBinContent(i);
  		double x = lenepl3d->GetBinCenter(i);
 		MSG("ShwFit",Msg::kDebug)<<"v "<<val<<" x "<<x<<" eval "<<predElec->Eval(x)<<" "<<predGamma->Eval(x)<<"\n";
  
  		//	double predval = predElec->Integral( lenepl3d->GetBinLowEdge(i),lenepl3d->GetBinLowEdge(i+1),(double *)0,1e-3);
  
		double bincenter = lenepl3d->GetBinCenter(i);
		double predval = predElec->Eval(bincenter);
  	
  	//	if(predval>1e-4 && val>0)
  		if(predval)chisq_elec+=(val-predval)*(val-predval)/(predval);
  		//else dc_e++;
  	
  		//predval = predGamma->Integral( lenepl3d->GetBinLowEdge(i),lenepl3d->GetBinLowEdge(i+1),(double *)0,1e-3);
		predval = predGamma->Eval(bincenter);
  	
  		//if(predval>1e-4 && val>0)
  		if(predval)
  			chisq_gamma+=(val-predval)*(val-predval)/(predval);
  		//else dc_g++;
  
	}

	pred_e_chisq=chisq_elec;
	pred_e_ndf=lenepl3d->GetNbinsX()-dc_e;
	pred_g_chisq=chisq_gamma;
	pred_g_ndf=lenepl3d->GetNbinsX()-dc_g;  
  
	chisq_elec/=lenepl3d->GetNbinsX()-dc_e;
	chisq_gamma/=lenepl3d->GetNbinsX()-dc_g;
  
 	post_over=0;
 	post_under=0;
 	pre_over=0;
 	pre_under=0;
	int m = lenepl3d->GetMaximumBin();
 



	//////////////////////////

	//do the fit

	int hmin=0;
   	int npare=3;
   
   	double pulseheight = lenepl3d->Integral();
	if(!debug)	if(efit){delete efit;efit=0;}

	efit = new TF1("efit",shwfunc3d,hmin+0.001,rlmax,npare);
 	efit->SetParNames("a","b","e0","zshift");
 	
 	//efit->SetParLimits(0,hmin+0.001,rlmax);
 	//efit->SetParLimits(1,0.001,20000); 	
 	//efit->SetParLimits(2,0+0.001,1000000);
  
        efit->SetParameters(3,0.5,300); 
 	efit->SetParLimits(0,1.5,5);
   	efit->SetParLimits(1,0.01,1.5); 	
   	efit->SetParLimits(2,0,1000000);
  
       

 	//efit->SetParameters(3,-5,5);
   
  	//efit->SetParameters(3.,0.5,pulseheight *60*25,0); //gues e as sum of ph in projection....
   
	TCanvas * c=0; 
	
	if(debug){
  		MSG("ShwFit",Msg::kDebug)<<(psf!=0)<<" "<<(c!=0)<<"\n";
		if(psf){
  			c = new TCanvas("psfcan","psfcan");
  			c->cd();
			lenepl3d->Draw();
		}
		MSG("ShwFit",Msg::kDebug)<<(psf!=0)<<" "<<(c!=0)<<"\n";
	}

	if(!debug){		
	   	lenepl3d->Fit(efit,"NWLLQ"); ///fast version //N do not draw
	}else{
	    if(psf)
 		   	lenepl3d->Fit(efit,"WLLV+"); ///detail version
		else
			lenepl3d->Fit(efit,"NWLLV"); ///detail version
	}	


	TH1F * predhist= 0;

	if(psf){
		predhist = new TH1F("predhistPSF","longitudinal energy by plane",nbins,0.0,rlmax);
	//	predhist->SetDirectory(0);
		predhist->SetLineColor(2);
	}else{
		predhist = new TH1F("predhist","longitudinal energy by plane",nbins,0.0,rlmax);
		predhist->SetDirectory(0);
		predhist->SetLineColor(2);
	}
	
	
  	for(int i=1;i<lenepl3d->GetNbinsX()+1;i++)
 	{

		double x = lenepl3d->GetBinCenter(i);
   		//	double pred = predElec->Integral(lenepl3d->GetBinLowEdge(i),lenepl3d->GetBinLowEdge(i+1),(double *)0,1e-3);
		double pred = predElec->Eval(lenepl3d->GetBinCenter(i));
   
		double val = lenepl3d->GetBinContent(i);
   	
		if(predhist)predhist->Fill(x,pred);
  	
		if(i>m)
		{
			if(pred<val)
				post_over+=val-pred;
			else
				post_under+=pred-val;
		}else{
			if(pred<val)
				pre_over+=val-pred;
			else
				pre_under+=pred-val;
		}
	}

	if(debug&&psf&&predhist)predhist->Draw("same");

	//draw the predicted fit here...
	if(debug){
		MSG("ShwFit",Msg::kDebug)<<(psf!=0)<<" "<<(c!=0)<<"\n";

		if(psf && c){
			c->cd();
			predElec->Draw("same");
			predGamma->Draw("same");
			MSG("ShwFit",Msg::kDebug)<<"ASDFASDFASDF\n";
		}
	}

	pp_chisq=0;
	pp_ndf=0;
	pp_igood=0;
	pp_p=0;
	if(psf){
			//require nonempty histograms to prevent endless loop!
	
		if(predhist->Integral()>0.01 && lenepl3d->Integral()>0.01)
		{
			MSG("ShwFit",Msg::kDebug)<<"int pred "<<(predhist->Integral()>0)<<" shw "<< (lenepl3d->Integral()>0)<<"\n";	
			//pp_p=lenepl3d->Chi2TestX(predhist,pp_chisq,pp_ndf,pp_igood,"WW",(double*)0);
		}
	}
	
	//lets directly compare the predicted histogram with the original one
	cmp_chisq=0;
	cmp_ndf = 0;
	peakdiff=0;
	if(predhist)
	{
		cmp_ndf=predhist->GetNbinsX();
		for(int k=1;k<cmp_ndf+1;k++)
		{
			double oc = lenepl3d->GetBinContent(k);
			double pc = predhist->GetBinContent(k);
			if(pc)cmp_chisq+=(oc-pc)*(oc-pc)/pc;
			else if(oc)cmp_ndf--;
		}
		if(cmp_chisq>0)cmp_chisq=sqrt(cmp_chisq);
		
		int maxbin = predhist->GetMaximumBin();
		if(maxbin>1)
		{
			peakdiff+= lenepl3d->GetBinContent(maxbin-1)-predhist->GetBinContent(maxbin-1);
			peakdiff+= lenepl3d->GetBinContent(maxbin)-predhist->GetBinContent(maxbin);
			peakdiff+= lenepl3d->GetBinContent(maxbin+1)-predhist->GetBinContent(maxbin+1);					
		}else if(maxbin==1)
		{
			peakdiff+= lenepl3d->GetBinContent(maxbin)-predhist->GetBinContent(maxbin);
			peakdiff+= lenepl3d->GetBinContent(maxbin+2)-predhist->GetBinContent(maxbin+2);
			peakdiff+= lenepl3d->GetBinContent(maxbin+1)-predhist->GetBinContent(maxbin+1);					
		}
	
	//	printf("-- %f %d %f\n",cmp_chisq,cmp_ndf, peakdiff);


	}
	



	MSG("ShwFit",Msg::kDebug)<<" STATUS: "<<gMinuit->fCstatu
			<<" "<<efit->GetParameter(0)
			<<" "<<efit->GetNDF()<<" "<<pulseheight<<endl;

	string fitstatus = (string)(gMinuit->fCstatu);
	MSG("ShwFit",Msg::kDebug)<<" STATUS: "<<fitstatus
			<<", "<<efit->GetParameter(0)
			<<", "<<efit->GetNDF()<<" "<<pulseheight<<endl;

	if(fitstatus=="CONVERGED "&&efit->GetParameter(0)<29.9&&
		efit->GetNDF()>0&&pulseheight>0){

		MSG("ShwFit",Msg::kDebug)<<" filling vars"<<endl;

		par_a=efit->GetParameter(0);
		par_b=efit->GetParameter(1);
		par_e0=efit->GetParameter(2);
  
		par_a_err=efit->GetParError(0);
		par_b_err=efit->GetParError(1);
		par_e0_err=efit->GetParError(2);
      
		chisq=efit->GetChisquare();
     
		ndf = efit->GetNDF();
    

		//ndf should be equal to the number of bins... because we are giving 0 bins weight 1
		ndf = nbins;

		if(debug && psf){
			c->cd();
			TLatex l;
			l.SetTextAlign(12);
			l.SetNDC(1);
	
			char txt[100];
			sprintf(txt,"\\chi^2=%f",efit->GetChisquare());
			l.DrawLatex(0.5,0.8,txt);    
   
  	 		sprintf(txt,"NDF=%f",ndf);
			l.DrawLatex(0.5,0.75,txt);    
  	 
 	 		sprintf(txt,"\\chi^2/NDF=%f",efit->GetChisquare()/ndf);
			l.DrawLatex(0.5,0.7,txt);    

 	 		sprintf(txt,"par a %f \\pm %f",par_a,par_a_err);
			l.DrawLatex(0.5,0.65,txt);    
		
 	 		sprintf(txt,"par b %f \\pm %f",par_b,par_b_err);
			l.DrawLatex(0.5,0.6,txt);    

 	 		sprintf(txt,"e0 %f \\pm %f",par_e0,par_e0_err);
			l.DrawLatex(0.5,0.55,txt);    
		
			sprintf(txt,"elec \\chi^2/NDF = %f ",chisq_elec);
			l.DrawLatex(0.5,0.50,txt);    

			sprintf(txt,"gamma \\chi^2/NDF = %f ",chisq_gamma);
			l.DrawLatex(0.5,0.45,txt);  

			sprintf(txt,"pre over %f under %f",pre_over,pre_under);
			l.DrawLatex(0.5,0.40,txt);  
			sprintf(txt,"post over %f under %f",post_over,post_under);
			l.DrawLatex(0.5,0.35,txt);  
		
			sprintf(txt,"%f --- %f %d %f", pp_p ,pp_chisq,pp_ndf,pp_chisq/pp_ndf);
			l.DrawLatex(0.5,0.3,txt); 	

			if(cmp_ndf)
			{
				sprintf(txt,"cmp c2 %f ndf %d c2ndf %f ", cmp_chisq,cmp_ndf,cmp_chisq/(double)cmp_ndf);
				l.DrawLatex(0.5,0.25,txt); 				
			}

			sprintf(txt,"pd %f ",peakdiff);
			l.DrawLatex(0.5,0.2,txt); 	
			
		}


     
		shwmax=1.*GetMaximumX(efit);
		prob=efit->GetProb();
	
		shwmaxplane=(int)(lenepl3d->
			GetBinContent(lenepl3d->FindBin(shwmax)));
		conv=1;
	}


	if(!psf || !debug)
	{
		delete predElec;predElec=0;
		delete predGamma;predGamma=0;
		delete predhist;predhist=0;
	}	

	if(!psf)
	{
		delete efit;efit=0;

	}
}



