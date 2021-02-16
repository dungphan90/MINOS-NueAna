/// $Id: AngClusterFitAna.cxx,v 1.21 2007/03/01 16:38:49 rhatcher Exp $ 
///
/// class AngClusterFitAna
///
/// NueAna package
///
/// Purpose: Fit a EM shower profile to the 3D hits contained in 
///          the primary cluster determined by AngClusterAna.. 
///
/// Author: Alex Sousa <asousa@minos.phy.tufts.edu>
///
/// Created on: Fri May 06 2005

#include <iostream>
#include <algorithm>
#include <TCanvas.h>
#include <TRotation.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "MessageService/MsgService.h"
#include "NueAna/AngClusterFitAna.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

CVSID("$Id: AngClusterFitAna.cxx,v 1.21 2007/03/01 16:38:49 rhatcher Exp $");

using std::cout;
using std::endl;

ClassImp(AngClusterFitAna)

const Float_t kDetectorPitch=0.05949; //Munits::meters

const Float_t kStripWidth=0.041; //Munits::meters

const Float_t kRadLengthToPlane=1.54; //X0->planes conversion factor

const Float_t kMolRadDet=0.0391; //Characteristic Moliere Radius 
                                 //of the Minos Detectors

const Float_t kMeuToGeV=1.0;  //Conversion factor from MEU to GeV

AngClusterFitAna::AngClusterFitAna(AngClusterFit &acf) :
    fLongHitHisto(0),
    fLongHitHistoRad(0),
    fTransHitHisto(0),    
    fTransHitHistoRad(0),
    fScatterAxis(0),
    fLongFit(0),
    fTransFit(0),
    fAngClusterFit(acf)
{
        
}

AngClusterFitAna::~AngClusterFitAna()
{
    //Prevent memory leaks
    if(fLongHitHisto!=0){
        fLongHitHisto->Delete();
        fLongHitHisto=0;
    }
    
    if(fLongHitHistoRad !=0){
        fLongHitHistoRad->Delete();
    }
    if(fTransHitHisto !=0){
        fTransHitHisto->Delete();
        fTransHitHisto=0;
    }
    
    if(fTransHitHistoRad !=0){
        fTransHitHistoRad->Delete();
        fTransHitHistoRad=0;
    }

    if(fScatterAxis !=0){
        fScatterAxis->Delete();
        fScatterAxis=0;
    }

    if(fLongFit !=0){
        fLongFit->Delete();
        fLongFit=0;
    }

    if(fTransFit !=0){
        fTransFit->Delete();
        fTransFit=0;
    }

}

void AngClusterFitAna::Set3DHit(DeqFloat_t &x
                                , DeqFloat_t &y
                                , DeqFloat_t &z
                                , DeqFloat_t &e){
    
    if(x.size()!=0 && y.size()!=0 && z.size()!=0 && e.size()!=0){

        fX=x;
        fY=y;
        fZ=z;
        fE=e;  
    
    }

}

void AngClusterFitAna::SetAngCluster(Int_t &primShow
                                     , DeqDeqInt_t &clusterMap
                                     , TVector3 &primDir)
{

    if(clusterMap.size()!=0){

        fClusterMap=clusterMap;
        fPrimShow  =primShow;
        fPrimDir   =primDir;

    } else fPrimShow=ANtpDefVal::kInt; //Will trigger stopping of 
                                       //further computations.
    
}


static Double_t shwLongFunc(Double_t *x, Double_t *par)
{

    // Longitudinal profile to be fitted to em showers 
    //xx=Distance from the vertex in units of X0 
    //par[0]=a, par[1]=b, par[2]=E0
   
    Float_t xx=x[0];
    
    Double_t lnf = TMath::Log(kRadLengthToPlane*par[2]*par[1])+(par[0]-1)
        *TMath::Log(xx*par[1])-par[1]*xx-TMath::LnGamma(par[0]);
    
    Double_t f = TMath::Exp(lnf);

    return f;
}

static Double_t shwTransFunc(Double_t *x, Double_t *par)
{
    //Transverse fit profile
    //x[0]=Radial distance from shower axis in units of Moliere radius
    //par[0]=Lambda1 attenuation length for high energy core
    //par[1]=Lambda2 attenuation length for shower halo
    //par[2]=Relative weight between lambda1 and lambda2
    //par[3]=E0 Energy normalization

    Float_t xx=x[0];
    Double_t f=par[3]*(TMath::Exp(-TMath::Sqrt(xx/par[0]))
                       +par[2]*TMath::Exp(-xx/par[1]));    
    return f;
}

void AngClusterFitAna::Analyze(int evtn
                            ,RecRecordImp<RecCandHeader> *srobj)
{

    if (ANtpDefVal::IsDefault(fPrimShow)){
        
        MSG("AngClusterFitAna",Msg::kInfo)<< "No primary cluster " 
                                          << " found in event."
                                          << endl << endl;
        fAngClusterFit.Reset();
        return;
    }

    if(fZ.size() < 3 || fZ.size() > 1000) return;  

    if(fPrimDir.X()==0.0 && fPrimDir.Y()==0.0 && fPrimDir.Z()==0) return;
    
    MSG("AngClusterFitAna",Msg::kDebug)<<"In AngClusterFitAna::Analyze"<<endl;
    MSG("AngClusterFitAna",Msg::kDebug)<<"On Snarl "
                                      <<srobj->GetHeader().GetSnarl()
                                      <<" event "<<evtn<<endl;
    
    if(srobj==0){
        return;
    }
    if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
       ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
        return;
    }

    Float_t cosZAngle=fPrimDir.CosTheta(); //Angle between 
                                           //Primary Shower dir 
                                           //and the Z axis.
    
    const Int_t kLHistBin=40;
    const Int_t kTHistBin=40;
    
    const Float_t kTHistXLim=static_cast<Float_t>(kTHistBin)*kStripWidth;

    char lhh[100];
    char thh[100];
    char lhhr[100];
    char thhr[100];

    sprintf(lhh ,"LongHitHisto_%d_%d"    ,srobj->GetHeader().GetSnarl(),evtn);
    sprintf(lhhr,"LongHitHistoRad_%d_%d" ,srobj->GetHeader().GetSnarl(),evtn);
    sprintf(thh ,"TransHitHisto_%d_%d"   ,srobj->GetHeader().GetSnarl(),evtn);
    sprintf(thhr,"TransHitHistoRad_%d_%d",srobj->GetHeader().GetSnarl(),evtn);
   

    fLongHitHisto  = new TH1F(lhh,"longitudinal energy hit profile"
                              ,kLHistBin,0.0,kLHistBin);
    
    fLongHitHistoRad  = new TH1F(lhhr,"longitudinal energy hit profile(X0)"
                                 ,kLHistBin
                                 ,0.0
                                 ,kLHistBin*kRadLengthToPlane/cosZAngle);
   
    fTransHitHisto = new TH1F(thh,"trasverse energy hit profile",
			      kTHistBin,-kTHistXLim,kTHistXLim);

    fTransHitHistoRad = new TH1F(thhr,"trasverse energy hit profile(Rm)",
                                 kTHistBin,0,kTHistBin);


    TVector3 unitDir; //Unit vector parallel to PrimShow dir
    unitDir=fPrimDir.Unit();

    TRotation r,a;

    r.SetToIdentity();
    a.SetToIdentity();
    
    if (TMath::Abs(fPrimDir.Z()) > 0.0){  //Do not rotate Z axis if shower dir
                                          //goes along it.

        r.SetZAxis(unitDir); //Rotates the Z axis to be parallel to the 
                             //primDir vector
       
        a=r.Inverse(); //Inverting the rotation means when one
                       //applies the transformation to the hit 
                       //vector, one obtains the primDir coordinates
    }
    
    TVector3 hitTemp;

    DeqDeqInt_t cluTemp; 
    cluTemp.push_back(fClusterMap.at(fPrimShow)); //Get PrimShow cluster

    Float_t rPrime; //Hit Radius around the shower direction. 

    Float_t totalHitEnergy=0.0;

    //Loop over primary shower hits and fill Long and Trans histograms
    for (IterDeqDeqInt_t cluIter = cluTemp.begin()
             ;cluIter != cluTemp.end()
             ; ++cluIter ){
        for (IterDeqInt_t hitIter =  cluIter->begin()
                 ;hitIter < cluIter->end()
                 ; ++hitIter){
            
            hitTemp.SetX(fX.at(*hitIter));
            hitTemp.SetY(fY.at(*hitIter));
            hitTemp.SetZ(fZ.at(*hitIter));
     
            totalHitEnergy+=(fE.at(*hitIter))*kMeuToGeV;

            Float_t nPlane=hitTemp.Z()/(kDetectorPitch);
            
            Float_t nPlaneDir=nPlane/(cosZAngle);

            Float_t radDist=nPlaneDir*kRadLengthToPlane; //Longitudinal hit 
                                                         //distance from 
                                                         //vertex in X0 units.

            //Transform hit coordinate system to Zprime
            hitTemp.Transform(a);

            rPrime=TMath::Sqrt(hitTemp.X()*hitTemp.X()
                               +hitTemp.Y()*hitTemp.Y());

            Float_t rPrimeRad=rPrime/(kMolRadDet); //Transverse hit distance 
                                                   //from Zprime in Moliere 
                                                   //radius units

            fLongHitHisto    ->Fill(nPlane,fE.at(*hitIter)*kMeuToGeV);
            
            fLongHitHistoRad ->Fill(radDist,fE.at(*hitIter)*kMeuToGeV);

            fTransHitHisto   ->Fill(hitTemp.X(),fE.at(*hitIter)*kMeuToGeV);

            fTransHitHistoRad->Fill(rPrimeRad,fE.at(*hitIter)*kMeuToGeV);

        }//for (hitIter)
    }//for(cluIter)

    //Do the Longitudinal and Transverse fits
    FitShower(totalHitEnergy);

    //fill transverse variables
    TransVarCalc(fTransHitHisto);

    MSG("AngClusterFitAna",Msg::kDebug)<<"Leaving AngClusterFitAna::Analyze"
                                       <<endl;
}

void AngClusterFitAna::FitShower(Float_t totalHitEnergy)
{
	
    //Configure fit functions
    Int_t lMin=0;
    Int_t lMax=100;
    Int_t tMin=0;
    Int_t tMax=100;
    Int_t nParLong=3;
    Int_t nParTrans=4;
    char lfn[100];
    char tfn[100];
    sprintf(lfn,"lfit_%d_%d",lMin,lMax);
    sprintf(tfn,"tfit_%d_%d",tMin,tMax);
    
    fLongFit  = new TF1(lfn,shwLongFunc,lMin+0.001,lMax,nParLong);
    fTransFit = new TF1(tfn,shwTransFunc,tMin+0.001,tMax,nParTrans);
   
    fLongFit->SetParNames("a","b","E0");
    fLongFit->SetParLimits(0,lMin+0.001,lMax);
    fLongFit->SetParLimits(1,0.001,20000);
    fLongFit->SetParLimits(2,0.001,20000);

    fTransFit->SetParNames("l1,l2,c12,E0");
    fTransFit->SetParLimits(0,tMin+0.001,tMax);
    fTransFit->SetParLimits(1,tMin+0.001,tMax);
    fTransFit->SetParLimits(2,0.001,50);
    fTransFit->SetParLimits(3,0.001,20000);

    //Do the Longitudinal fit
    fLongFit->SetParameters(3.,0.5,totalHitEnergy);
    fLongHitHistoRad->Fit(fLongFit,"RLQ0+");

    string fitStatus = (string)(gMinuit->fCstatu);
    MSG("AngClusterFitAna",Msg::kDebug)<<"LongFit status:  "<<fitStatus<<endl;
	
    if(fitStatus=="CONVERGED " 
       && fLongFit->GetParameter(0)<29.9
       && fLongFit->GetNDF()>0
       && totalHitEnergy > 0.0){ 
        
        fAngClusterFit.fACluFitParA     =fLongFit->GetParameter(0);
        fAngClusterFit.fACluFitParB     =fLongFit->GetParameter(1);
        fAngClusterFit.fACluFitParLongE0=fLongFit->GetParameter(2);
        fAngClusterFit.fACluFitShwMax   =GetMaximum(fLongFit);
        fAngClusterFit.fACluFitLongChiSq=fLongFit->GetChisquare();
        fAngClusterFit.fACluFitLongConv =1;
        fAngClusterFit.fACluFitLongNDF
             =fAngClusterFit.fACluFitLongChiSq/fLongFit->GetNDF();
        fAngClusterFit.fACluFitE0EnergyRatio
            =fAngClusterFit.fACluFitParLongE0/totalHitEnergy;
    }//if(fitStatus)

    //Do the transverse fit
    fTransFit->SetParameters(1.0,1.0,1.0,totalHitEnergy);
    fTransHitHistoRad->Fit(fTransFit,"RLQ0+");

    fitStatus = (string)(gMinuit->fCstatu);
    MSG("AngClusterFitAna",Msg::kDebug)<<"TransFit status:  "<<fitStatus<<endl;
	

    if(fitStatus=="CONVERGED " 
       && fTransFit->GetNDF()>0 ){ 
        fAngClusterFit.fACluFitParL1     =fTransFit->GetParameter(0);
        fAngClusterFit.fACluFitParL2     =fTransFit->GetParameter(1);
        fAngClusterFit.fACluFitParC12    =fTransFit->GetParameter(2);
        fAngClusterFit.fACluFitParTransE0=fTransFit->GetParameter(3);
        fAngClusterFit.fACluFitTransChiSq=fTransFit->GetChisquare();
        fAngClusterFit.fACluFitTransConv =1;
        fAngClusterFit.fACluFitTransNDF
             =fAngClusterFit.fACluFitTransChiSq/fTransFit->GetNDF();
    }

    MSG("AngClusterFitAna",Msg::kDebug)<< "total energy " 
                                       << totalHitEnergy 
                                       << " from long histo " 
                                       << fLongHitHistoRad->Integral() 
                                       << " from trans histo " 
                                       << fTransHitHistoRad->Integral() 
                                       << endl;
}

void AngClusterFitAna::TransVarCalc(TH1F *transHisto)
{
    Int_t kTransHistoBin = static_cast <Int_t> (transHisto->GetNbinsX());
    
    Int_t binMax   =transHisto->GetMaximumBin();
    Int_t binVert  =Int_t(((kTransHistoBin-1)/2)+1);

    Float_t peakBin=transHisto->Integral(binMax,binMax);
    Float_t vertBin=transHisto->Integral(binVert,binVert);
    Float_t total  =transHisto->Integral(1,kTransHistoBin);
  
    Float_t asymPeak= ANtpDefVal::kFloat;
    Float_t asymVert= ANtpDefVal::kFloat;

    MSG("AngClusterFitAna",Msg::kDebug)<<"total:  "
                                       <<total<<" binMax: "<<binMax
                                       <<" binVert "<<binVert
	                               <<"  peakBin: "<<peakBin
                                       <<" kTransHistBin "<<kTransHistoBin
                                       <<endl;

    MSG("AngClusterFitAna",Msg::kDebug)<<"HistoI1:  "
                                       <<transHisto->Integral(binMax+1
                                                              ,kTransHistoBin)
                                       <<endl;
    
    MSG("AngClusterFitAna",Msg::kDebug)<<"HistoI1:  "
                                       <<transHisto->Integral(1, binMax-1)
                                       <<endl;
    
    MSG("AngClusterFitAna",Msg::kDebug)<<"HistoI1:  "
                                       <<transHisto->Integral(binVert+1
                                                              ,kTransHistoBin)
                                       <<endl;
    
    MSG("AngClusterFitAna",Msg::kDebug)<<"HistoI1:  "
                                       <<transHisto->Integral(1, binVert-1)
                                       <<endl;
   


    
    if(TMath::Abs(total-peakBin) > 1e-5){
        MSG("AngClusterFitAna",Msg::kDebug)<<"Evaluating asymPeak:  "
                                           <<TMath::Abs(transHisto
                                                        ->Integral(binMax+1
                                                        ,kTransHistoBin)
                                                        -transHisto
                                                        ->Integral(1,binMax-1))
                                           <<" divided by "<<(total-peakBin)
                                           <<endl;
		
        asymPeak=(TMath::Abs(transHisto->Integral(binMax+1,kTransHistoBin)
                             -transHisto->Integral(1,binMax-1)))
                             /(total-peakBin);
    }

    
    if(TMath::Abs(total-vertBin) > 1e-5){
       MSG("AngClusterFitAna",Msg::kDebug)<<"Evaluating asymVert:  "<<
               TMath::Abs(transHisto->Integral(binVert+1,kTransHistoBin)
                                            -transHisto->Integral(1,binVert-1))
                  <<" divided by "<<(total-peakBin)<<endl;
		    
        asymVert=(TMath::Abs(transHisto->Integral(binVert+1,kTransHistoBin)
                             -transHisto->Integral(1,binVert-1)))
                             /(total-vertBin);
    }


    MSG("AngClusterFitAna",Msg::kDebug)<<"asymPeak:  "<<asymPeak
		                       <<" asymVert: "<<asymVert<<endl;
    
    Float_t ratio,molRadPeak=0,molRadVert=0;
    
    if(total){
        for(Int_t i=0; i<=binVert-1;i++){
            ratio=transHisto->Integral(binMax-i>0?binMax-i
                                       :1,binMax+i<kTransHistoBin
                                       ?binMax+i:kTransHistoBin)/total;
            if(ratio>0.90){molRadPeak=i+1; break;}
        }
        for(Int_t i=0; i<=binVert-1;i++){
            ratio=transHisto->Integral(binVert-i,binVert+i)/total;
            if(ratio>0.90){molRadVert=i+1; break;}
        }
        
    }
    
    Float_t mean=transHisto->GetMean();
    Float_t rms =transHisto->GetRMS();
    Float_t skew=0,kurt=0;
    Int_t nCount=0;
                                                   
    MSG("AngClusterFitAna",Msg::kDebug)<<"mean:  "<<mean<<" rms: "<<rms<<endl;
	
    
    for(Int_t i=1;i<=kTransHistoBin;i++){
        
        if(transHisto->GetBinContent(i)){
            skew=skew+(TMath::Power((transHisto->GetBinCenter(i)-mean),3)
                       *transHisto->GetBinContent(i));
            kurt=kurt+(TMath::Power((transHisto->GetBinCenter(i)-mean),4)
                       *transHisto->GetBinContent(i));
            nCount++;
        }
        
    }
    
    if(rms>1e-8 && nCount>1){
        skew=skew /(static_cast<Float_t>(nCount-1)*TMath::Power((rms),3.));
        kurt=(kurt/(static_cast<Float_t>(nCount-1)*TMath::Power((rms),4.)))-3.;
        if(TMath::Abs(skew) > 1.e8 ) skew = ANtpDefVal::kFloat;
        if(TMath::Abs(kurt) > 1.e8 ) kurt = ANtpDefVal::kFloat;
    }
    else{skew= ANtpDefVal::kFloat; kurt= ANtpDefVal::kFloat;}
    
    //Fill variables
    fAngClusterFit.fACluFitAsymPeak=asymPeak;
    fAngClusterFit.fACluFitAsymVert=asymVert;
    fAngClusterFit.fACluFitMolRadPeak=molRadPeak;
    fAngClusterFit.fACluFitMolRadVert=molRadVert;
    fAngClusterFit.fACluFitMean=mean;
    fAngClusterFit.fACluFitRMS=rms;
    fAngClusterFit.fACluFitSkew=skew;
    fAngClusterFit.fACluFitKurt=kurt;
   
    MSG("AngClusterFitAna",Msg::kDebug)
        <<"Leaving AngClusterFitAna::TransCalc"
        <<endl;
}

void AngClusterFitAna::Draw(TPad *pad)
{

    if (ANtpDefVal::IsDefault(fPrimShow)){
        
        MSG("AngClusterFitAna",Msg::kInfo)<< "No primary cluster " 
                                          << " found in event."
                                          << endl << endl;
        fAngClusterFit.Reset();
        return;
    }

    if(fZ.size() < 3 || fZ.size() > 1000) return;  
    if(fPrimDir.X()==0.0 && fPrimDir.Y()==0.0 && fPrimDir.Z()==0) return;

    DeqTPoly scatterCluster;
    TPolyMarker3D *scatterTemp;

    DeqDeqInt_t scatTemp; 

    Int_t nClu=0;
    Int_t nHit=0;

    //Loop over all cluster hits and fill polymarkers
    for (IterDeqDeqInt_t cluIter = fClusterMap.begin()
             ;cluIter != fClusterMap.end()
             ; ++cluIter ){
        scatTemp.push_back(fClusterMap.at(nClu));
        scatterTemp = new TPolyMarker3D(scatTemp.size(),1);
        
        for (IterDeqInt_t hitIter =  cluIter->begin()
                 ;hitIter < cluIter->end()
                 ; ++hitIter){
            
            scatterTemp->SetPoint(nHit,fZ.at(*hitIter)
                                  ,fX.at(*hitIter)
                                  ,fY.at(*hitIter));
            nHit++;
        }//for (hitIter)
        scatterCluster.push_back(scatterTemp);
        nHit=0;
        nClu++;
    }//for(cluIter)

    // Find min and max values;
    Float_t minX=static_cast <Float_t> (*min_element(fX.begin(),fX.end()));
    Float_t maxX=static_cast <Float_t> (*max_element(fX.begin(),fX.end()));
    Float_t minY=static_cast <Float_t> (*min_element(fY.begin(),fY.end()));
    Float_t maxY=static_cast <Float_t> (*max_element(fY.begin(),fY.end()));
    Float_t minZ=static_cast <Float_t> (*min_element(fZ.begin(),fZ.end()));
    Float_t maxZ=static_cast <Float_t> (*max_element(fZ.begin(),fZ.end()));
    
    if(fScatterAxis) delete fScatterAxis;
    fScatterAxis = new TH3F("cluHisto","Scatter plot of all hits"
                           , 10,minZ-0.1,maxZ+0.2
                           , 10,minX-0.2,maxX+0.2
                           , 10,minY-0.2,maxY+0.2);
    
    fScatterAxis->SetXTitle("Z(m)");
    fScatterAxis->SetYTitle("X(m)");
    fScatterAxis->SetZTitle("Y(m)");
    
    TPolyMarker3D *evtVertex;
    evtVertex = new TPolyMarker3D(1,8);
    evtVertex->SetMarkerColor(1);
    evtVertex->SetPoint(0,0.,0.,0.);
    
    TPolyLine3D *lineCluster;
    lineCluster= new TPolyLine3D(2);
    
    lineCluster->SetPoint(0,0.,0.,0.);
    lineCluster->SetPoint(1,fPrimDir.Z(),fPrimDir.X(),fPrimDir.Y());
    
    pad->cd(1);
    
    fScatterAxis->Draw();

    Int_t polyColor=1;

    for (IterDeqTPoly polyIter = scatterCluster.begin()
             ;polyIter < scatterCluster.end()
             ; ++polyIter){
        (*polyIter)->SetMarkerStyle(7);
        (*polyIter)->SetMarkerColor(polyColor);
        (*polyIter)->Draw("sames");
        polyColor++;   
    }//for(polyIter)

    evtVertex->Draw("sames");
    lineCluster->SetLineWidth(2);
    lineCluster->Draw("sames");

    pad->cd(3);
    fLongHitHistoRad->Draw();
    fLongFit->Draw("sames");
    pad->cd(2);
    
    fTransHitHistoRad->Draw();
    fTransFit->Draw("sames");
    pad->cd(4);
    fTransHitHisto->Draw();
    pad->cd();
    pad->Modified();
    
    //Clearing ROOT objects
    evtVertex->Delete();
    scatterTemp->Delete();
    lineCluster->Delete();
    
}

void AngClusterFitAna::Reset()
{


}

Double_t AngClusterFitAna::GetMaximum(TF1* efit, Double_t xmin, Double_t xmax)
{
   Double_t fXmin = 0.001;
   Double_t fXmax = 100;
   Double_t fNpx = 40;
                                                                                
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
   else return GetMaximum(efit, TMath::Max(xmin,xxmax-dx), TMath::Min(xmax,xxmax+dx));
}

