/// $Id: HitCalcAna.cxx,v 1.27 2008/11/19 18:22:51 rhatcher Exp $
///
/// class HitCalcAna
///
/// NueAna package
///
/// Purpose: Calculate HitCalc variables passed by HitCalc object. 
///
/// Author: Alex Sousa <asousa@minos.phy.tufts.edu>
///
/// Created on: Mon Apr 11 2005

#include <iostream>
#include <algorithm>
#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "MessageService/MsgService.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/HitCalcAna.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "VertexFinder/NtpVtxFinder/NtpVtxFinder.h"

CVSID("$Id: HitCalcAna.cxx,v 1.27 2008/11/19 18:22:51 rhatcher Exp $");

using std::cout;
using std::endl;

ClassImp(HitCalcAna)


//......................................................................


HitCalcAna::HitCalcAna(HitCalc &hc):
    fHitCalc(hc)
{
    
  
}

HitCalcAna::~HitCalcAna()
{

    
}

void HitCalcAna::Analyze(int evtn
                         ,RecRecordImp<RecCandHeader> *srobj)
{
    
    if(srobj==0){
        return;
    }
    if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
       ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
        return;
    }
    
    MSG("HitCalcAna",Msg::kInfo)<<"In HitCalcAna::Analyze"<<endl;
    MSG("HitCalcAna",Msg::kInfo)<<"On Snarl "<<srobj->GetHeader().GetSnarl()
                                 <<" event "<<evtn<<endl;
   
    NtpSREvent *event = SntpHelpers::GetEvent(evtn,srobj);
    if(!event){
        MSG("HitCalcAna",Msg::kError)<<"Couldn't get event "<<evtn
                                     <<" from Snarl "
                                     <<srobj->GetHeader().GetSnarl()<<endl;
        return;
    }
    
    Int_t hitTotal;
    
    //Match strips into 3D Hits.
    hitTotal=ComputeHits(srobj,event);

    if (hitTotal > 3){
        //Calculate simple hit variables and fill branch.
        VarCalc();

    }
    
}

void HitCalcAna::Get3DHit(DeqFloat_t &x
                           , DeqFloat_t &y
                           , DeqFloat_t &z
                           , DeqFloat_t &e)
{

    if(fX.size()==0 || fY.size()==0 || fZ.size()==0 || fE.size()==0){
        MSG("HitCalcAna",Msg::kWarning)<< "3D Hits not calculated for event,"
                                       << "stopping any further"
                                       << " Hit related processing." 
                                       << endl;            
        return;
    }

    //Reference the 3D Hit vectors.
    x=fX;
    y=fY;
    z=fZ;
    e=fE;

}


Int_t HitCalcAna::ComputeHits(RecRecordImp<RecCandHeader> *srobj
                              ,NtpSREvent *event)
{
 
    if(srobj==0){
        return 1;
    }
    if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
       ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
        return 1;
    }
    
    // Clear 3D Hit Coordinate vectors
    fX.clear();
    fY.clear();
    fZ.clear();
    fE.clear();

    //Determine detector type from Validity Context
    Detector::Detector_t fDetectorType = 
        srobj->GetHeader().GetVldContext().GetDetector();
    
    //Scale sigcor to pe.
    //    const Float_t kNFScale = 1.0/100.0; // 92.3/65.2; 

    Float_t totalVisEnerg;
    
    // Auxiliary hit variables
    Float_t  hitX;
    Float_t  hitY;
    Float_t  hitZ;
    Float_t  hitE;
      
    //Strip vectors
    DeqFloat_t stripUV;
    DeqFloat_t stripZ;
    DeqFloat_t stripE;

    DeqFloat_t stripUVTemp;
    DeqFloat_t stripZTemp;
    DeqFloat_t stripETemp;

    Int_t nPlaneSM = 0; // account for U and V even-odd switch 
                        // from FD SM1 to SM2 

    Int_t nMax = 0;     //Total number of Hits found
    Int_t nHit = 0;     //Number of hits counter
    
    Float_t vertexX = event->vtx.x; //Munits::meters
    Float_t vertexY = event->vtx.y; //Munits::meters
    Float_t vertexZ = event->vtx.z; //Munits::meters
   
    if(ReleaseType::IsCedar(release)){
      NtpSRTrack *primaryTrack = 0;
      NtpSRTrack *ntpTrack = 0;
      Int_t index = 0;

      for(Int_t i = 0; i < event->ntrack; ++i){
        index = event->trk[i];
        ntpTrack = SntpHelpers::GetTrack(index, srobj);
        if(!primaryTrack) primaryTrack = ntpTrack;

        if(ntpTrack->plane.n > primaryTrack->plane.n)
          primaryTrack = ntpTrack;
      }

      NtpStRecord* st = dynamic_cast<NtpStRecord *>(srobj);
      NtpVtxFinder vtxf;  
      vtxf.SetTargetEvent(st, event, primaryTrack);
      if(vtxf.FindVertex() > 0){
         vertexX = vtxf.VtxX();
         vertexY = vtxf.VtxY();
         vertexZ = vtxf.VtxZ();
      }
    }
 
    Int_t event_end   = event -> plane.end;
    Int_t event_start = event -> plane.beg;

    if(fDetectorType == Detector::kNear){
        if (event_end >= 121) {
            MSG("HitCalcAna",Msg::kInfo)<< "ND event ends outside"
                                        << " calorimeter region in plane " 
                                        << event_end 
                                        <<". Truncated hit matching." 
                                        << endl;            
            event_end = 121;
        }
        if (event_start >= 121) {
            MSG("HitCalcAna",Msg::kError)<< "ND event starts outside"
                                         << " calorimeter region in " 
                                         << event_start <<" Stopping." 
                                         << endl;
            return 1;
        }
    }
    
    totalVisEnerg = 0.0;
    nHit=0;
    
    Float_t stpEDiff = 0.0;

    //Loop over all planes in the event
    for(Int_t nPlane = event_start; 
        nPlane < event_end; nPlane++){
    
        //Reset the strip vectors for a new plane
        stripUVTemp.clear();
        stripZTemp.clear();
        stripETemp.clear();

        stripUV.clear();
        stripZ.clear();
        stripE.clear();	  

        // Loop over all strips in the event
        for(Int_t iStp = 0; iStp < event->nstrip; iStp++){ 
            Int_t index = SntpHelpers::GetStripIndex(iStp,event);
            NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
        
            if(!strip){
                MSG("HitCalcAna",Msg::kDebug)<< "Couldn't get strip "<< index
                                             << " in snarl "
                                             << srobj->GetHeader().GetSnarl()
                                             << " next"<<endl;
                continue;
            }
 	    if(!evtstp0mip){
	       MSG("HitCalcAna",Msg::kError)<<"No mip strip information"<<endl;
	       continue;
	    }


            float charge = evtstp0mip[index] + evtstp1mip[index];;
            
            totalVisEnerg += charge;
            
            //Eliminate crosstalk strips 
            if(TMath::Abs(strip->z-vertexZ) < 4. && (charge)>200.0/sigcormeu){
                
                if(strip->plane == nPlane){
                    stripUVTemp.push_back(strip->tpos);
                    stripZTemp.push_back(strip ->z);
                    stripETemp.push_back(charge);
                }
	
                if(strip->plane == (nPlane+1)){
                    stripUV.push_back(strip->tpos);
                    stripZ.push_back(strip ->z);
                    stripE.push_back(charge);
                }                
            }
        } // iStp


        // Loop over all strips and match hits with strips in 
        // consecutive planes
        for(UInt_t nStrip = 0; nStrip < stripZ.size(); nStrip++) {
            for(UInt_t nStripTemp = 0; nStripTemp < stripZTemp.size()
                    ; nStripTemp++){
                
                MSG("HitCalcAna",Msg::kDebug)<< "StripUV(" << nStrip 
                                             << ") = " 
                                             << stripUV.at(nStrip) << "  "
                                             << "StripUVTemp(" << nStripTemp 
                                             << ") = " 
                                             << stripUVTemp.at(nStripTemp) 
                                             << endl;
                MSG("HitCalcAna",Msg::kDebug)<< "StripZ(" << nStrip 
                                             << ") = " 
                                             << stripZ.at(nStrip) << "  "
                                             << "StripZTemp(" << nStripTemp 
                                             << ") = " 
                                             << stripZTemp.at(nStripTemp) 
                                             << endl;

                MSG("HitCalcAna",Msg::kDebug)<< "StripE(" << nStrip 
                                             << ") = " 
                                             << stripE.at(nStrip) << "  "
                                             << "StripETemp(" << nStripTemp 
                                             << ") = " 
                                             << stripETemp.at(nStripTemp) 
                                             << endl;                     

                //Only match hits if there are two strips in 
                //consecutive planes without unphysical values. 
                if(stripUV.at(nStrip) != stripUVTemp.at(nStripTemp) 
                   && (stripUV.at(nStrip) > -99. 
                       && stripUVTemp.at(nStripTemp) > -99.)){
                    
                    if(fDetectorType == Detector::kFar){ //far detector
                        stpEDiff=10; //sigcormeu (1GeV~25.76 sigcormeu)
                        
                        if(nPlane >= 0 && nPlane <= 249) 
                            nPlaneSM = nPlane+1;  //u - even and v - odd
                        if(nPlane >= 250 && nPlane <= 485) 
                            nPlaneSM = nPlane;  //u - odd and v - even
                    }
                    
                    if(fDetectorType == Detector::kNear||
                       fDetectorType==Detector::kCalDet){
                        nPlaneSM = nPlane; // near detector
                        stpEDiff=10; //sigcormeu
                    }
                    // Only match hits if adjacent strips have reasonably 
                    // similar energy deposition.
                    
//		     cout << "StpEDiff              " << stpEDiff << endl;
                    
                    
                    if(TMath::Abs(stripE.at(nStrip)
                                  -stripETemp.at(nStripTemp)) < stpEDiff ){
                        
                        MSG("HitCalcAna",Msg::kDebug)<< "Found Hit " 
                                                     << nHit << endl;

                        hitX = pow((-1.0),nPlaneSM)*(sqrt(2.0)/2.0)
                            *(stripUV.at(nStrip)-stripUVTemp.at(nStripTemp));
                        hitY = (sqrt(2.0)/2.0)
                            *(stripUV.at(nStrip)+stripUVTemp.at(nStripTemp));
                        hitZ = (stripZ.at(nStrip)
                                +stripZTemp.at(nStripTemp))/2.0;
                        hitE = (stripE.at(nStrip)
                                +stripETemp.at(nStripTemp))/2.0;
                     
                        fX.push_back(hitX-vertexX);
                        fY.push_back(hitY-vertexY);
                        fZ.push_back(hitZ-vertexZ);
                        fE.push_back(hitE);
                        
                        MSG("HitCalcAna",Msg::kDebug)<< "Vertex Pos(X,Y,Z)= (" 
                                                     << vertexX << ","
                                                     << vertexY << ","
                                                     << vertexZ << ")"
                                                     << endl << endl;
                        
                        MSG("HitCalcAna",Msg::kDebug)<< "fX[" 
                                                     << nHit <<"]= "
                                                     << fX[nHit]
                                                     << endl << endl;
                        
                        MSG("HitCalcAna",Msg::kDebug)<< "fY[" 
                                                     << nHit <<"]= "
                                                     << fY[nHit]
                                                 << endl << endl;
                        
                        MSG("HitCalcAna",Msg::kDebug)<< "fZ[" 
                                                     << nHit <<"]= "
                                                     << fZ[nHit]
                                                     << endl << endl;
                        
                        MSG("HitCalcAna",Msg::kDebug)<< "fE[" 
                                                     << nHit <<"]= "
                                                     << fE[nHit]
                                                     << endl << endl;
                        
                        nHit++;
                        
                    } else continue;
                }
            } // nStripTemp
        } // nStrip
    } // nPlane 
    
    nMax = nHit;

    MSG("HitCalcAna",Msg::kInfo)<< "Found " 
                                << nMax 
                                << " 3D Hits." 
                                << endl;

    return nMax;
}


void HitCalcAna::Reset()
{
 

}

void HitCalcAna::VarCalc()
{
    
    //Variable declaration
    
    Float_t totalHitEnerg;
    Float_t arg;
    Float_t argTheta;
    Float_t eTrans;
    Float_t eTransX;
    Float_t eTransY;
    Float_t sumETransX;
    Float_t sumETransY;
    Float_t eLong;
    Float_t sTrans;
    
    Float_t phi;
    Float_t phiMax;
    Float_t phiEMax;
    Float_t thetaLMax;
    Float_t thetaEMax;
    Float_t s;
    Float_t sE;
    
    Float_t length;
    Float_t lengthMax;
    Float_t eMax;
    
    Float_t sumOppSide;
    Float_t sumSameSide;
    Float_t ratioFar;
    
    Float_t sumOppSideE;
    Float_t sumSameSideE;
    Float_t ratioEnerg;
    
    Int_t nFar;
    Int_t nEMax;
    
    Float_t zMax;
    Float_t zL;
    
    VecFloat_t theta;
    theta.assign(fZ.size(),0.0);
    
    //////////////////////////////////
    
    sumETransX=0.0;
    sumETransY=0.0;
    eLong=0.0;
    sTrans=0.0;
    length=0.0;
    lengthMax=0.0;
    eMax=0.0;
    phiMax=0.0;
    phiEMax=0.0;
    thetaLMax=0.0;
    thetaEMax=0.0;
    nFar=0;
    nEMax=0;
    zMax=0;
    totalHitEnerg=0.0;
    
    Float_t epsilon=0.00001;
    
    //Preliminary Calc of Total Hit Energy and Hit with highest z
    for (UInt_t n = 0; n < fZ.size(); n++){
        if(fZ.at(n) >= 0.0){      //Get rid of Hits behind the vertex
            totalHitEnerg += fE.at(n);
            if(fZ.at(n) > zMax) zMax = fZ.at(n);
        }
    }
    
    
    for (UInt_t n = 0; n < fZ.size(); n++){ //BIG Loop over all 3D Hits
        if(fZ.at(n) >= 0.0){      //Get rid of Hits behind the vertex
            if(fZ.at(n) == 0){ fZ.at(n) += epsilon; } //Protect against nans
            
            arg      = fX.at(n)*fX.at(n)+fY.at(n)*fY.at(n);            
            argTheta = TMath::Sqrt(arg)/fZ.at(n);
            theta.at(n) = TMath::ATan(argTheta);
            phi      = TMath::ACos(fX.at(n)/(argTheta*fZ.at(n) + epsilon));

            MSG("HitCalcAna",Msg::kDebug)<< "Intermediate VarCalc Monitor " 
                                         << endl
                                         << "arg = " << arg << endl
                                         << "argTheta = " << argTheta 
                                         << endl
                                         << "fE[" << n << "] = " << fE.at(n) 
                                         << endl
                                         << "fX = " << fX.at(n) << endl
                                         << "fY = " << fY.at(n) << endl
                                         << "fZ = " << fZ.at(n) << endl
                                         << "theta =  " << theta.at(n) 
                                         << endl << endl;
            
            eTransX = fE.at(n)*fX.at(n)*TMath::Cos(theta.at(n))/fZ.at(n);
            eTransY = fE.at(n)*fY.at(n)*TMath::Cos(theta.at(n))/fZ.at(n);
            sumETransX += eTransX;
            sumETransY += eTransY;
            
            sTrans += (eTransX*eTransX+eTransY*eTransY); //Get square of the 
                                                         //transverse energy

            eLong  += fE.at(n)*TMath::Cos(theta.at(n)); //Get Longitudinal 
                                                        //Energy
            
            MSG("HitCalcAna",Msg::kDebug)<< "Intermediate VarCalc Monitor " 
                                         << endl
                                         << "sumETransX[" << n << "] = " 
                                         << sumETransX << endl
                                         << "sumETransY[" << n << "] = " 
                                         << sumETransY << endl
                                         << "eLong = " << eLong 
                                         << endl << endl;
 
            //Calculate Farthest Hit
            length = TMath::Sqrt(arg+fZ.at(n)*fZ.at(n));
            if(length > lengthMax){
                if ( fZ.at(n) > 0.8*zMax 
                     && fZ.at(n) <= zMax){ //Try to avoid 
                                           //far noise hits 
                                           //causing very high 
                                           //theta values.
                    lengthMax = length;
                    nFar = n;
                    zL = fZ.at(n);
                    phiMax = phi;
                    thetaLMax = theta.at(n);
                }
            }
            
            //Calculate Highest Energy Hit:
            if(fE.at(n) > eMax){
                eMax = fE.at(n);
                nEMax = n;
                phiEMax = phi;
                thetaEMax = theta.at(n);
            }
        }
    }
    
    eTrans = TMath::Sqrt(sumETransX*sumETransX
                         +sumETransY*sumETransY); //Transverse 
                                                  //energy
    
    sumSameSide  = 0.0;
    sumOppSide   = 0.0;
    sumSameSideE = 0.0;
    sumOppSideE  = 0.0;
    
    for (UInt_t n=0; n < fZ.size(); n++){
        if(fZ.at(n) >= 0.0){
            //Calculate Ratio Far variable
            s = fX.at(n)*TMath::Cos(phiMax)+fY.at(n)*TMath::Sin(phiMax);
            
            if(s >= 0.0) sumSameSide += fE.at(n)*TMath::Sin(theta.at(n));
            else         sumOppSide  += fE.at(n)*TMath::Sin(theta.at(n));
            
            //Calculate Ratio variable for showers.
            sE = fX.at(n)*TMath::Cos(phiEMax)+fY.at(n)*TMath::Sin(phiEMax);
            
            if(sE >= 0.0) sumSameSideE += fE.at(n)*TMath::Sin(theta.at(n));
            else          sumOppSideE  += fE.at(n)*TMath::Sin(theta.at(n));

        }
    }
    
    if (sumSameSide > 0.0){
        ratioFar = sumOppSide/sumSameSide;
    }else ratioFar = ANtpDefVal::kFloat;
    
    if (sumSameSideE > 0.0){
        ratioEnerg = sumOppSideE/sumSameSideE;
    }else ratioEnerg = ANtpDefVal::kFloat;
    
    if (ratioFar   > 50.0) {ratioFar   = ANtpDefVal::kFloat;}
    if (ratioEnerg > 50.0) {ratioEnerg = ANtpDefVal::kFloat;}
    
    MSG("HitCalcAna",Msg::kInfo)<< "Calculated Hit variables "  << endl
                                << "eTotal = " << totalHitEnerg 
                                << " eTrans = " << eTrans 
                                << " eLong = " << eLong 
                                << " sTrans = " << sTrans 
                                << " farMomBalance = " << ratioFar
                                << " peakMomBalance = " << ratioEnerg 
                                << " farAngle = " << thetaLMax 
                                << " peakAngle = " << thetaEMax 
                                << " zMax = " << zMax << endl;
    
    Float_t transLongRatio=ANtpDefVal::kFloat;
    if(eLong > 0.0){transLongRatio=eTrans/eLong;}
    
    //Fill the variables
    if (totalHitEnerg > 0.0){
        fHitCalc.fHitTotalEnergy          = totalHitEnerg;    
        fHitCalc.fHitTransEnergy          = eTrans;
        fHitCalc.fHitLongEnergy           = eLong;
        fHitCalc.fHitTransCMEnergy        = sTrans;    
        fHitCalc.fHitTransEnergyRatio     = eTrans/totalHitEnerg;    
        fHitCalc.fHitLongEnergyRatio      = eLong/totalHitEnerg;    
        fHitCalc.fHitTransLongEnergyRatio = transLongRatio;
        fHitCalc.fHitTransCMEnergyRatio   = sTrans/totalHitEnerg;    
        fHitCalc.fHitFarMomBalance        = ratioFar;
        fHitCalc.fHitPeakMomBalance       = ratioEnerg;
        fHitCalc.fHitFarAngle             = thetaLMax;
        fHitCalc.fHitPeakAngle            = thetaEMax;
    }
    else {
        MSG("HitCalcAna",Msg::kInfo)<< "No hits found in Event" << endl;
        fHitCalc.fHitTotalEnergy          = ANtpDefVal::kFloat;    
        fHitCalc.fHitTransEnergy          = ANtpDefVal::kFloat;
        fHitCalc.fHitLongEnergy           = ANtpDefVal::kFloat;
        fHitCalc.fHitTransCMEnergy        = ANtpDefVal::kFloat;    
        fHitCalc.fHitTransEnergyRatio     = ANtpDefVal::kFloat;    
        fHitCalc.fHitLongEnergyRatio      = ANtpDefVal::kFloat;    
        fHitCalc.fHitTransLongEnergyRatio = ANtpDefVal::kFloat;
        fHitCalc.fHitTransCMEnergyRatio   = ANtpDefVal::kFloat;    
        fHitCalc.fHitFarMomBalance        = ANtpDefVal::kFloat;
        fHitCalc.fHitPeakMomBalance       = ANtpDefVal::kFloat;
        fHitCalc.fHitFarAngle             = ANtpDefVal::kFloat;
        fHitCalc.fHitPeakAngle            = ANtpDefVal::kFloat;
        
    }
}
