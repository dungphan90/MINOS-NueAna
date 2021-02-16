/// $Id: AngClusterAna.cxx,v 1.24 2007/06/06 15:35:30 boehm Exp $
///
/// class AngClusterAna
///
/// NueAna package
///
/// Purpose: Cluster 3D Hits using Cos(x) vs Cos(y) representations.
///
/// Author: Alex Sousa <asousa@minos.phy.tufts.edu>
///
/// Created on: Thu Apr 28 2005

#include <iostream>
#include <algorithm>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "MessageService/MsgService.h"
#include "NueAna/AngClusterAna.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "VertexFinder/NtpVtxFinder/NtpVtxFinder.h"



CVSID("$Id: AngClusterAna.cxx,v 1.24 2007/06/06 15:35:30 boehm Exp $");

using std::cout;
using std::endl;

ClassImp(AngClusterAna)

const Float_t kDetectorPitch=0.05949; //Munits::meter
const Float_t kStripWidth=0.041; //Munits::meters

const Float_t kMeuToGeV=1.0;  //Conversion factor from MEU to GeV

//......................................................................

AngClusterAna::AngClusterAna(AngCluster &ac):
    fAngCluster(ac)
{


}

AngClusterAna::~AngClusterAna()
{

}

void AngClusterAna::Analyze(int evtn
                            ,RecRecordImp<RecCandHeader> *srobj)
{

    MSG("AngClusterAna",Msg::kInfo)<<"In AngClusterAna::Analyze"<<endl;
    MSG("AngClusterAna",Msg::kDebug)<<"On Snarl "<<srobj->GetHeader().GetSnarl()
                                 <<" event "<<evtn<<endl;

    if(srobj==0){
        return;
    }
    if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
       ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
        return;
    }


    if(fZ.size() >=3 && fZ.size() <= 1000){

        //Find Clusters
        FindCluster();

        //Calculate cluster variables and fill branch.
        ClusterVarCalc();
    }

    WeightedEnergy(evtn, srobj);
}

void AngClusterAna::Set3DHit(DeqFloat_t &x
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

void AngClusterAna::GetAngCluster(Int_t &primShow
                                  , DeqDeqInt_t &clusterMap
                                  , TVector3 &primDir)
{
    
    if(fClusterMap.size()==0){
        MSG("AngClusterAna",Msg::kWarning)<< "Primary cluster not found for " 
                                       << "event, stopping any further"
                                       << " clustering related processing"
                                       << endl;            
        return;
    }

    if(fPrimDir.X()==0.0 && fPrimDir.Y()==0.0 && fPrimDir.Z()==0) return;

    //Reference relevant quantities.
    primShow  =fPrimShow;
    clusterMap=fClusterMap;
    primDir   =fPrimDir;
}



void AngClusterAna::FindCentroidBlob(TMatrixD & grid
                                  ,const Int_t nBinX
                                  ,const Int_t nBinY
                                  ,Int_t indX
                                  ,Int_t indY
                                  ,TMatrixD & blobCentro)
{
    //Recursively find centroid agglomerations in the centroid histogram
    // No approximations are used, all centroid aggregates are found.
    
    //Make sure we did not go over the edges.    
    if (indX < 0 || indY < 0 || indX >= nBinX || indY >= nBinY)
        return;
    
    //Finish recursion if the next grid point is empty.
    if(TMath::Nint(grid(indX,indY)) != 1)
        return;
    
    //else
    grid(indX,indY) = 0.0;    //Mark as used
    blobCentro(indX,indY)=1.0;
  
    MSG("AngClusterAna",Msg::kDebug)<< "Adding point(" << indX
                                    <<"," << indY << ") to blob " 
                                    << endl;
  
    FindCentroidBlob(grid, nBinX, nBinY, indX - 1, indY, blobCentro);
    FindCentroidBlob(grid, nBinX, nBinY, indX + 1, indY, blobCentro);
    FindCentroidBlob(grid, nBinX, nBinY, indX, indY - 1, blobCentro);
    FindCentroidBlob(grid, nBinX, nBinY, indX, indY + 1, blobCentro);
    
}

void AngClusterAna::FindCluster()
{

   
    //Array of TVector3 hit objects containing
    //normal and XZY and ZYX CS rotations  
    std::vector<TVector3> hitVecArrayX(fZ.size());
    std::vector<TVector3> hitVecArrayY(fZ.size());
    
    //Vectors to be used in clustering
    VecFloat_t cosX;
    VecFloat_t cosY;
    cosX.assign(fZ.size(),0.0);
    cosY.assign(fZ.size(),0.0);
    
    fTotalEventHitEnergy=0.0;

    for (UInt_t n=0; n < hitVecArrayX.size(); n++){ //Loop over all 3D Hits
        
        (hitVecArrayX.at(n)).SetZ(fX.at(n));
        (hitVecArrayX.at(n)).SetY(fY.at(n));
        (hitVecArrayX.at(n)).SetX(fZ.at(n));

        (hitVecArrayY.at(n)).SetX(fX.at(n));
        (hitVecArrayY.at(n)).SetZ(fY.at(n));
        (hitVecArrayY.at(n)).SetY(fZ.at(n));

        cosX.at(n)=(hitVecArrayX.at(n)).CosTheta();
        cosY.at(n)=(hitVecArrayY.at(n)).CosTheta();

        //Calculate total 3D hit event energy  
        fTotalEventHitEnergy+=fE.at(n)*kMeuToGeV;
        
    }

    //============================================================
    //Start of cluster finding algorithm (nearest-neighbor method)
    //Do it in CosTheta X-Y space for now
    
    TMatrixD distHit(fZ.size(),fZ.size());   //Holds all distances 
                                             //between hits
    TMatrixD usedDist(fZ.size(),fZ.size());  //Keep track of already 
                                             //used distances
    
    VecFloat_t usedHit; //keep track of already used hits
    usedHit.assign(fZ.size(),0.0);
    
    Int_t nMin;  //index of smallest inter-hit distance
    Int_t kMin;
    Float_t minDist=999999.;

    //Calculate all inter-hit distances ommitting reciprocals
    for (UInt_t n=0; n < fZ.size(); n++){
        for (UInt_t k=0; k < fZ.size(); k++){
            distHit(n,k) =0.0;
            usedDist(n,k)=0.0;
            if (k > n){

                distHit(n,k)= TMath::Sqrt((cosX.at(k)-cosX.at(n))
                                          *(cosX.at(k)-cosX.at(n))
                                          +(cosY.at(k)-cosY.at(n))
                                          *(cosY.at(k)-cosY.at(n)));
               
                if (distHit(n,k) < minDist){
                    minDist=distHit(n,k);
                    nMin=n;
                    kMin=k;
                }
                
            } //if
        }//for(k)
    }//for(n)
    
    MSG("AngClusterAna",Msg::kDebug)<< "minDist[" 
                                    << nMin << "][" << kMin << "] = " 
                                    << distHit(nMin,kMin) << endl;


    Float_t radius=0.5; //limit radius for nearest neighbor 

    DeqInt_t aggHitSize; //Keeps track of the size of each aggregate 

    DeqDeqInt_t aggHitArray; //Container filled with all the 
                             //hit aggregate sets.

    UInt_t aggHitSizeMax=0;

    Int_t ahaIndex=0;
    Int_t ahaIndexMax=0;

    for (UInt_t n=0; n < fZ.size(); n++){
        for (UInt_t k=(n+1); k < fZ.size(); k++){
            if (k > n){
                //Take the (n,k) pair as an aggregate seed 
                if (TMath::Nint(usedDist(n,k)) == 0
                    && distHit(n,k) < radius
                    && distHit(n,k) > 0){
                    
                    usedDist(n,k)=1.0;
                    
                    DeqInt_t aggHit; //aggregate hit array
                                       
                    aggHit.push_back(n); //Store first hit from seed
                    aggHit.push_back(k); //Store second hit
                    
                    //Look for neighbor hits above the upper hit of 
                    //the seed pair
                    for (UInt_t kTemp=(n+1); kTemp < fZ.size(); kTemp++){
                        if (TMath::Nint(usedDist(n,kTemp)) == 0){
                            usedDist(n,kTemp)=1.0;
                            if (kTemp != k 
                                && distHit(n,kTemp) > 0.0 
                                && distHit(n,kTemp) < radius){
                                
                                aggHit.push_back(kTemp);
                                usedDist(n,kTemp)=1.0;

                            }//if(kTemp !=k)
                                                       
                        }else continue; //if (usedDist(n,kTemp)) 
                  
                    }//for(kTemp)
                    
                    //Look for neighbor hits below the lower hit of 
                    //the seed pair
                    for (UInt_t nTemp=0; nTemp < k; nTemp++){
                        if (TMath::Nint(usedDist(nTemp,k)) == 0){
                            usedDist(nTemp,k)=1.0;
                            if (nTemp != n 
                                && distHit(nTemp,k) > 0.0 
                                && distHit(nTemp,k) < radius){
                                
                                aggHit.push_back(nTemp);
                                usedDist(nTemp,k)=1.0;

                            }//if(nTemp !=n)

                        }else continue; //if (usedDist(nTemp,k)) 
                  
                    }//for(nTemp)

                    aggHitSize.push_back(aggHit.size()); //Record aggreg. size
                    aggHitArray.push_back(aggHit);       //Store hit aggregate

                    if (aggHit.size() > aggHitSizeMax){
                        aggHitSizeMax=aggHit.size();
                        ahaIndexMax=ahaIndex;
                    }
                    ahaIndex++;
                }//if(distHit < radius)
                
            }//if(k > n)

        }//for(k)
        
    }//for(n)
    
    
    MSG("AngClusterAna",Msg::kInfo)<< "Found "
                                   << aggHitArray.size() 
                                   << " aggregates."
                                   << endl;
    
    MSG("AngClusterAna",Msg::kInfo)<< "Largest aggregate contains "
                                   << aggHitSizeMax 
                                   << " Hits and has Index "
                                   << ahaIndexMax
                                   << endl;
 
    //Cluster determination:
    //Eliminate all overlapping aggregates in favor of the largest one.
    //Method: Calculate and histogram aggregate centroids. Run a 
    //recursion blob finder function on the histogram and determine 
    //high centroid density region bounds. From all the centroids
    //consistent with a region boundary, choose the one 
    //corresponding to the aggregate with most hits and call 
    //the aggregate a cluster.

    const Int_t kBinCentroX=10;
    const Int_t kBinCentroY=10;

    const Float_t kMinCentroX =-1.0;
    const Float_t kMaxCentroX = 1.0;
    const Float_t kMinCentroY =-1.0;
    const Float_t kMaxCentroY = 1.0;

    Float_t centSumX=0.0;
    Float_t centSumY=0.0;
    Float_t centValX=0.0;
    Float_t centValY=0.0;
    
    // Write out centroid histogram for debugging purposes 
    Bool_t kDebug=0;
    
    TH2F *centroHisto;
    centroHisto = new TH2F("centroHisto","CosX vs CosY centroid histogram"
                           ,kBinCentroX,kMinCentroX,kMaxCentroX
                           ,kBinCentroY,kMinCentroY,kMaxCentroY);

    centroidX.clear();
    centroidY.clear();

    for (IterDeqDeqInt_t ahaIter = aggHitArray.begin()
             ;ahaIter != aggHitArray.end()
             ; ++ahaIter ){
        centSumX=0.0;
        centSumY=0.0;
        for (IterDeqInt_t ahIter = ahaIter->begin()
                 ;ahIter < ahaIter->end()
                 ; ++ahIter){
            
            centSumX+=cosX.at((*ahIter));
            centSumY+=cosY.at((*ahIter)); 
        }//for(ahIter)
        
        centValX=centSumX/(static_cast <Float_t> (ahaIter->size()));
        centValY=centSumY/(static_cast <Float_t> (ahaIter->size()));
        
        centroidX.push_back(centValX);
        centroidY.push_back(centValY);
        
        centroHisto->Fill(centValX,centValY);
        
    }//for(ahaIter)

        
    TMatrixD gridHisto(kBinCentroX,kBinCentroY); //Centroid Grid
                                                 //for blob finding 

    const Float_t kMinCentroid=1.0; //Minimum number of centroids 
                                    //to be found in a  bin that 
                                    // will be part of a blob;
    
    if (kMinCentroid < 1.0){
        MSG("AngClusterAna",Msg::kFatal)<< "Attempting to cluster "
                                        << "empty space. Check "
                                        << "value of kMinCentroid"
                                        << endl;
    }
    
    Float_t centroStepX
        =(kMaxCentroX-kMinCentroX)/static_cast <Float_t> (kBinCentroX);
    Float_t centroStepY
        =(kMaxCentroY-kMinCentroY)/static_cast <Float_t> (kBinCentroY);

    //Map each centroid to bin position.
    VecInt_t centroBinX;
    centroBinX.assign(aggHitSize.size(),0);
    VecInt_t centroBinY;
    centroBinY.assign(aggHitSize.size(),0);

    //Populate blob grid.
    //Reminder: ROOT bin indexing starts at 1
    for (Int_t cBinX=0; cBinX < kBinCentroX; cBinX++){
      for (Int_t cBinY=0; cBinY < kBinCentroY; cBinY++){

          gridHisto(cBinX,cBinY)=0.0;
          
          //Find bin coordinates for each centroid
          for(UInt_t ah=0; ah < aggHitArray.size(); ah++){

              if (centroidX.at(ah) > (centroHisto->GetXaxis()
                                      ->GetBinLowEdge(cBinX+1)) 
                  && centroidX.at(ah) < centroStepX
                                        +(centroHisto
                                          ->GetXaxis()
                                          ->GetBinLowEdge(cBinX+1)) 
                  && centroidY.at(ah) > (centroHisto->GetYaxis()
                                         ->GetBinLowEdge(cBinY+1)) 
                  && centroidY.at(ah) < centroStepY
                                        +(centroHisto
                                          ->GetYaxis()
                                          ->GetBinLowEdge(cBinY+1))){
              
                  centroBinX.at(ah)=cBinX;
                  centroBinY.at(ah)=cBinY;

              }//if(centroid)
                  
          }//for(ah)


          //Ignore bins with low centroid density
          Int_t cBin=centroHisto->GetBin(cBinX+1,cBinY+1);
          if (centroHisto->GetBinContent(cBin) >= kMinCentroid ){
          
              gridHisto(cBinX,cBinY)=1.0;

          }//if(content)
      }//for(cBinY)
    }//for(cBinX)

    TMatrixD blobCentro(kBinCentroX,kBinCentroY); //Contains all centroids
                                                  //on a given blob

    DeqTMatrixD blobArray; //Contains all the found blobs 

    //Find blobs recursively in the centroid histogram
    for (Int_t nX=0; nX < kBinCentroX; nX++){
        for (Int_t nY=0; nY < kBinCentroY; nY++){
            blobCentro(nX,nY)=0.0;
            if(TMath::Nint(gridHisto(nX,nY))==1.0){
                
                FindCentroidBlob(gridHisto
                                 ,kBinCentroX
                                 ,kBinCentroY
                                 ,nX
                                 ,nY
                                 ,blobCentro);

                blobArray.push_back(blobCentro);
                if(kDebug) blobCentro.Print();
                blobCentro.Zero();
            }//if(gridHisto)
        }//for(nX)
    }//for(nY)

    
    Int_t sizeMax=-1;   //keep track of max agg size within a blob
    Int_t cluSelect=-1; //index of max agg size in a blob

    VecInt_t cluAggMap; //Correspondance between final cluster and
                        //largest Hit aggregate for that cluster
    cluAggMap.assign(aggHitArray.size(),0);

    Int_t cluSizeMax=0;   //keep track of maximum size over all blobs
    Int_t cluPrimShow=-1; //index of largest aggregate over all blobs
    Int_t aggIndexMax=-1;

    fClusterSize.clear(); //Prepare for the new event.
    fClusterMap.clear();

    Int_t blobIndex=0;

    //Find the largest aggregate for each blob and call that a cluster
    for (IterDeqTMatrixD baIter = blobArray.begin()
             ;baIter != blobArray.end()
             ; ++baIter ){

        sizeMax=-1;
        cluSelect=-1;

        for (Int_t i = 0; i < kBinCentroX; i++){
            for (Int_t j = 0; j < kBinCentroY; j++){
                if(TMath::Nint((*baIter)(i,j))==1){

                    for(UInt_t ah=0; ah < aggHitArray.size(); ah++){

                        if (centroBinX.at(ah)==i && centroBinY.at(ah)==j){

                            if(aggHitSize.at(ah) > sizeMax){
                                sizeMax=aggHitSize.at(ah);
                                cluSelect=ah;
                            }//if(aggHitSize)

                        }//if(centroBin)

                    }//for(ah)

                }//if(*baIter(i,j)==1)
            }//for(i)
        }//for(i)

        if(cluSelect >=0){
            cluAggMap.at(blobIndex)=cluSelect;

            //Fill cluster size member deque;
            fClusterSize.push_back(aggHitSize.at(cluSelect));

            //Fill cluster map member container
            fClusterMap.push_back(aggHitArray.at(cluSelect));

            //Determination of largest cluster index ("primary shower")
            if (aggHitSize.at(cluSelect) > cluSizeMax){
                cluSizeMax=aggHitSize.at(cluSelect);
                cluPrimShow=blobIndex;
                aggIndexMax=cluSelect;
            }

            MSG("AngClusterAna",Msg::kDebug)<< "Cluster "
                                       << blobIndex
                                       << " largest agg is "
                                       << cluSelect
                                       << " with "
                                       << sizeMax
                                       << " Hits"
                                       << endl;

            blobIndex++;

        } else{
            MSG("AngClusterAna",Msg::kWarning)<< "Unable to find "
                                           << " largest cluster."
                                           << " Quitting."
                                           << endl << endl;
            break;
        }

    }//for(baIter)



    if(cluPrimShow >= 0){
        fPrimShow=cluPrimShow;
        MSG("AngClusterAna",Msg::kInfo)<< "Primary Cluster "
                                       << "found in blob "
                                       << cluPrimShow << "/"
                                       << blobArray.size()
                                       << " contains "
                                       << aggHitSize.at
                                       (cluAggMap.at(cluPrimShow))
                                       << " hits and was found in aggregate "
                                       << aggIndexMax
                                       << endl;
    }else {

        MSG("AngClusterAna",Msg::kInfo)<< "No primary cluster "
                                       << " found in event."
                                       << endl << endl;
        if(centroHisto)
            delete centroHisto;

        fPrimShow=ANtpDefVal::kInt;

        return;
    }


    //Find all high energy hits close to the vertex (low angular resolution)
    DeqInt_t highEnergyHits;

    for (UInt_t n=0; n < fZ.size(); n++){
        if(fZ.at(n)/kDetectorPitch > 0. && fZ.at(n)/kDetectorPitch < 4.){
            if(TMath::Sqrt(fX.at(n)*fX.at(n)+fY.at(n)*fY.at(n))
               /kStripWidth < 4.){

                if(fE.at(n) > 5.0){
                    highEnergyHits.push_back(n);
		    //Add energy of close hits to total event energy
		    fTotalEventHitEnergy+=fE.at(n)*kMeuToGeV;
                    MSG("AngClusterAna",Msg::kDebug)<< "(X,Y,Z,E) HighHits = "
                                                   << "("
                                                   << fX.at(n) << ","
                                                   << fY.at(n) << ","
                                                   << fZ.at(n) << ","
                                                   << fE.at(n) << ")"
                                                   << endl;
                }//for(fE.at(n))
            }//for(fX.at(n))
        } //for(fZ.at(n))
    }// for(n)

    //Update the primary shower cluster with the high energy hits above
    DeqDeqInt_t clusterMapTemp;
    DeqInt_t    hitTemp;
    Int_t       cluIndexTemp=0;

    for (IterDeqDeqInt_t cluIter = fClusterMap.begin()
             ;cluIter != fClusterMap.end()
             ; ++cluIter ){

        hitTemp=*cluIter;

        if (cluIndexTemp==fPrimShow){

            for (UInt_t n=0; n < highEnergyHits.size(); n++){
                hitTemp.push_back(highEnergyHits.at(n));
            }//for(n)

        }//if(fPrimShow)

        clusterMapTemp.push_back(hitTemp);
        cluIndexTemp++;

    }//for(cluIter)

    //Fill member variables
    fClusterMap.clear();
    fClusterMap=clusterMapTemp;

    fClusterSize.at(fPrimShow)+=highEnergyHits.size();

    fNCluster=fClusterMap.size();

    Int_t cluIndex=0;

    fClusterID.clear();

    fClusterID.assign(fZ.size(),-1); //hits outside any cluster get -1.

    //Build the clusterID vector (returns the cluster one hit belongs to.)
    for (IterDeqDeqInt_t cluIter = fClusterMap.begin()
             ;cluIter != fClusterMap.end()
             ; ++cluIter ){
        for (IterDeqInt_t hitIter = cluIter->begin()
                 ;hitIter < cluIter->end()
                 ; ++hitIter){

            fClusterID.at((*hitIter))=cluIndex;

        }//for (hitIter)

        cluIndex++;

    }//for(cluIter)



    //Delete centroid histogram
    if(centroHisto)
        delete centroHisto;

}

void AngClusterAna::ClusterVarCalc()
{
    if (ANtpDefVal::IsDefault(fPrimShow)){

        MSG("AngClusterAna",Msg::kInfo)<< "No primary cluster "
                                       << " found in event."
                                       << endl << endl;
        fAngCluster.Reset();
        return;
    }

    //Get primary cluster centroid coordinates
    Float_t cosXCent=centroidX.at(fPrimShow);
    Float_t cosYCent=centroidY.at(fPrimShow);

    Float_t temp = 1.0-(cosXCent*cosXCent)-(cosYCent*cosYCent);
    if(temp < 0) temp = 0;
    Float_t cosZCent =TMath::Sqrt(temp);

    //Define histograms
    Int_t zBin =21;
    Int_t rBin =21;
    Float_t zMin =-0.05;
    Float_t rMin =-0.05;
    Float_t zMax =1.05;
    Float_t rMax =1.05;

    TH2F *hRZPrime;
    hRZPrime = new TH2F("RZPrime","r vs zprime"
                        , zBin, zMin, zMax
                        , rBin, rMin, rMax);

    TH2F *hRZ;
    hRZ = new TH2F("RZ","r vs z"
                   ,zBin, zMin, zMax
                   ,rBin, rMin, rMax);

    TH1F *hZPrime;
    hZPrime = new TH1F("ZPrime","ELongPrime"
                     ,zBin, zMin, zMax);
    
    TH1F *hRPrime;
    hRPrime = new TH1F("RPrime","ETransPrime"
                       ,rBin, rMin, rMax);
    
    TH1F *hZAxis;
    hZAxis = new TH1F("ZAxis","ELong"
                    ,zBin, zMin, zMax);
    
    TH1F *hRAxis;
    hRAxis = new TH1F("RAxis","ETrans"
                      ,rBin, rMin, rMax);

    
    //Iterate over all hits in primary cluster and 
    //calculate XYZ position centroid
    Float_t zPrime;
    Float_t distVertexSq;

    Float_t rPrime;
    Float_t rAxis;

    Float_t centroX=0;
    Float_t centroY=0;
    Float_t centroZ=0;

    //Needed to calculate cluster centroid weighted by distance to vertex.
    Float_t centroWeightX=0;
    Float_t centroWeightY=0;
    Float_t centroWeightZ=0;

    Float_t  totalHitEnergy=0.0;

    Float_t epsilon=0.0000001;

    DeqDeqInt_t cluTemp; 
    cluTemp.push_back(fClusterMap.at(fPrimShow)); //Select the 
                                                  //primary cluster
  
    for (IterDeqDeqInt_t cluIter = cluTemp.begin()
             ;cluIter != cluTemp.end()
             ; ++cluIter ){
        for (IterDeqInt_t hitIter =  cluIter->begin()
                 ;hitIter < cluIter->end()
                 ; ++hitIter){
            
            zPrime=fX.at(*hitIter)*cosXCent
                +fY.at(*hitIter)*cosYCent
                +fZ.at(*hitIter)*cosZCent;
            
            distVertexSq=fX.at(*hitIter)*fX.at(*hitIter)
                +fY.at(*hitIter)*fY.at(*hitIter)
                +fZ.at(*hitIter)*fZ.at(*hitIter);
            
            rPrime=TMath::Sqrt(
                TMath::Abs(distVertexSq-(zPrime*zPrime)));
            
            rAxis=TMath::Sqrt(
                TMath::Abs(distVertexSq-fZ.at(*hitIter)*fZ.at(*hitIter)));

            hRZPrime->Fill(zPrime,rPrime,fE.at(*hitIter));
            hRZ     ->Fill(fZ.at(*hitIter),rAxis,fE.at(*hitIter));
         
            hZPrime ->Fill(zPrime,fE.at(*hitIter));
            hRPrime ->Fill(rPrime,fE.at(*hitIter));

            hZAxis  ->Fill(fZ.at(*hitIter),fE.at(*hitIter));
            hRAxis  ->Fill(rAxis,fE.at(*hitIter));

            centroX+=fX.at(*hitIter);
            centroY+=fY.at(*hitIter);
            centroZ+=fZ.at(*hitIter);

            //Total Energy of primary cluster
            totalHitEnergy+=(fE.at(*hitIter))*kMeuToGeV;

            MSG("AngClusterAna",Msg::kDebug)<< "(X,Y,Z) = " << "(" 
                                            << fX.at(*hitIter) << ","
                                            << fY.at(*hitIter) << ","
                                            << fZ.at(*hitIter) << ")"
                                            << endl;
            
            MSG("AngClusterAna",Msg::kDebug)<< "centroTemp =" << "("
                                            << centroX << ","
                                            << centroY << ","
                                            << centroZ << ")"
                                            << endl;
           
            centroWeightX+=fX.at(*hitIter)*fX.at(*hitIter);
            centroWeightY+=fY.at(*hitIter)*fY.at(*hitIter);
            centroWeightZ+=fZ.at(*hitIter)*fZ.at(*hitIter);
            

        }//for (hitIter)
    }//for(cluIter)

    MSG("AngClusterAna",Msg::kDebug)<< "number of hits = " 
                                    << fClusterSize.at(fPrimShow)
                                    << endl;
        
    if (centroX == 0.0) centroX+=epsilon;
    if (centroY == 0.0) centroY+=epsilon;
    if (centroZ == 0.0) centroZ+=epsilon;

    centroWeightX/=centroX;
    centroWeightY/=centroY;
    centroWeightZ/=centroZ;

    centroX/=(fClusterSize.at(fPrimShow));
    centroY/=(fClusterSize.at(fPrimShow));
    centroZ/=(fClusterSize.at(fPrimShow));
  
    //Fill branch variables
    fAngCluster.fACluRmsShwAxis=hRPrime->GetRMS();
    fAngCluster.fACluRmsZAxis  =hRAxis->GetRMS();
    //     fAngCluster.fACluShwDirX   =centroWeightX;
    //     fAngCluster.fACluShwDirY   =centroWeightY;
    //     fAngCluster.fACluShwDirZ   =centroWeightZ;

    fAngCluster.fACluPrimEnergy     =totalHitEnergy;
    if (fTotalEventHitEnergy > 0.0)
        fAngCluster.fACluPrimEnergyRatio=totalHitEnergy/fTotalEventHitEnergy;

    fAngCluster.fACluShwDirX   =centroX;
    fAngCluster.fACluShwDirY   =centroY;
    fAngCluster.fACluShwDirZ   =centroZ;

    //Fill primary direction vector.
    fPrimDir.SetX(centroX);
    fPrimDir.SetY(centroY);
    fPrimDir.SetZ(centroZ);

    //Print variables
    MSG("AngClusterAna",Msg::kInfo)<< "fACluRmsShwAxis = "
                                   << fAngCluster.fACluRmsShwAxis
                                   << " fACluRmsZAxis = "
                                   << fAngCluster.fACluRmsZAxis
                                   << " fACluShwDirX = "
                                   << fAngCluster.fACluShwDirX
                                   << " fACluShwDirY = "
                                   << fAngCluster.fACluShwDirY
                                   << " fACluShwDirZ = "
                                   << fAngCluster.fACluShwDirZ
                                   << endl;

    //delete histograms
    if(hRZPrime)
        delete hRZPrime;

    if(hRZ)
        delete hRZ;

     if(hZPrime)
         delete hZPrime;

     if(hRPrime)
         delete hRPrime;

     if(hZAxis)
         delete hZAxis;

     if(hRAxis)
         delete hRAxis;

}

void AngClusterAna::WeightedEnergy(int evtn, RecRecordImp<RecCandHeader> *srobj){
  if(srobj==0){
    return;
  }
  if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
     ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
    return;
  }

  MSG("ShwfitAna",Msg::kDebug)<<"In ShwfitAna::Analyze"<<endl;
  MSG("ShwfitAna",Msg::kDebug)<<"On Snarl "<<srobj->GetHeader().GetSnarl()
			      <<" event "<<evtn<<endl;

//  Reset(srobj->GetHeader().GetSnarl(),evtn);

  NtpSREvent *event = SntpHelpers::GetEvent(evtn,srobj);

  if(!event){
      MSG("ShwfitAna",Msg::kError)<<"Couldn't get event "<<evtn
				   <<" from Snarl "<<srobj->GetHeader().GetSnarl()<<endl;
      return;
   }

  Float_t cent_x = fAngCluster.fACluShwDirX;
  Float_t cent_y = fAngCluster.fACluShwDirY;
  Float_t cent_z = fAngCluster.fACluShwDirZ;
  
  Float_t cent_u = TMath::Sqrt(2)*(cent_x-cent_y);
  Float_t cent_v = TMath::Sqrt(2)*(cent_x+cent_y);

  Int_t MaxStripU = 0;
  Int_t MaxStripV = 0;
  Float_t MaxStripU_ph = 0;
  Float_t MaxStripV_ph = 0;

  Float_t Weighted_ph0 = 0.0;
  Float_t Weighted_ph1 = 0.0;
  Float_t Weighted_ph2 = 0.0;
  Float_t Weighted_ph3 = 0.0;

  Int_t vtxPlane = event->vtx.plane;
  Float_t vtxU = event->vtx.u;
  Float_t vtxV = event->vtx.v;

  if(ReleaseType::IsCedar(release)){
    NtpStRecord* st = dynamic_cast<NtpStRecord *>(srobj);
    NtpVtxFinder vtxf(evtn, st);
    if(vtxf.FindVertex() > 0){
       vtxPlane = vtxf.VtxPlane();
       vtxU = vtxf.VtxU();
       vtxV = vtxf.VtxV();
    }
  }

  for(int i=0;i<event->nstrip;i++){
    Int_t index = SntpHelpers::GetStripIndex(i,event);
    NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
    if(!strip){
      continue;
    }
    if(!evtstp0mip){
       MSG("AngCluster",Msg::kError)<<"No mip strip information"<<endl;
       continue;
    }


    Float_t Strip_ph_meu = evtstp0mip[index] + evtstp1mip[index];
    //(strip->ph0.sigcor+strip->ph1.sigcor)/sigcormeu;

    if(strip->plane>=event->vtx.plane && strip->plane<=event->vtx.plane+(cent_z/0.0354) ){
      if(strip->planeview==PlaneView::kU){
        if(strip->tpos-event->vtx.u+cent_u*(strip->plane-event->vtx.plane)/(cent_z/0.0354)<.0205){
          if( Strip_ph_meu>MaxStripU_ph ){MaxStripU_ph=Strip_ph_meu; MaxStripU = i;}
        }
      }else if(strip->planeview==PlaneView::kV){
        if(strip->tpos-event->vtx.v+cent_v*(strip->plane-event->vtx.plane)/(cent_z/0.0354)<.0205){
          if( Strip_ph_meu>MaxStripV_ph ){MaxStripV_ph=Strip_ph_meu; MaxStripV = i;}
        }
      }
    }
  }

  Int_t Main_index_u = SntpHelpers::GetStripIndex(MaxStripU,event);
  NtpSRStrip *Main_strip_u = SntpHelpers::GetStrip(Main_index_u,srobj);
  Int_t Main_index_v = SntpHelpers::GetStripIndex(MaxStripV,event);
  NtpSRStrip *Main_strip_v = SntpHelpers::GetStrip(Main_index_v,srobj);

  for(int i=0;i<event->nstrip;i++){
    Int_t index = SntpHelpers::GetStripIndex(i,event);
    NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
    if(!strip){
      continue;
    }
    if(!evtstp0mip){
       MSG("MSTCalcAna",Msg::kError)<<"No mip strip information"<<endl;
       continue;
    }

    float charge = evtstp0mip[index] + evtstp1mip[index];;
    //  = (strip->ph0.sigcor+strip->ph1.sigcor)/sigcormeu


    Float_t Strip_ph_meu0 = charge;
    Float_t Strip_ph_meu1 = charge;
    Float_t Strip_ph_meu2 = charge;
    Float_t Strip_ph_meu3 = charge;

    if(strip->planeview==PlaneView::kU){
      if(Main_strip_u->plane-strip->plane>0){
        Strip_ph_meu0*=TMath::Sqrt(2*TMath::Abs(Main_strip_u->plane-strip->plane)*0.0354) + 40*TMath::Sqrt(TMath::Abs(Main_strip_u->tpos-strip->tpos));
        Strip_ph_meu1*=TMath::Sqrt(2*TMath::Abs(Main_strip_u->plane-strip->plane)*0.0354) + 20*TMath::Abs(Main_strip_u->tpos-strip->tpos);
        Strip_ph_meu2*=TMath::Sqrt( (5*TMath::Abs(Main_strip_u->plane-strip->plane)*0.0354) + 40*TMath::Abs(Main_strip_u->tpos-strip->tpos) );
        Strip_ph_meu3*=TMath::Sqrt( (2*TMath::Abs(Main_strip_u->plane-strip->plane)*0.0354) + 40*TMath::Abs(Main_strip_u->tpos-strip->tpos) );

      }else{
        Strip_ph_meu0*=TMath::Sqrt(TMath::Abs(Main_strip_u->plane-strip->plane)*0.0354) + 40*TMath::Sqrt(TMath::Abs(Main_strip_u->tpos-strip->tpos));
        Strip_ph_meu1*=TMath::Sqrt(TMath::Abs(Main_strip_u->plane-strip->plane)*0.0354) + 20*TMath::Abs(Main_strip_u->tpos-strip->tpos);
        Strip_ph_meu2*=TMath::Sqrt( (TMath::Abs(Main_strip_u->plane-strip->plane)*0.0354) + 40*TMath::Abs(Main_strip_u->tpos-strip->tpos) );
        Strip_ph_meu3*=TMath::Sqrt( (TMath::Abs(Main_strip_u->plane-strip->plane)*0.0354) + 40*TMath::Abs(Main_strip_u->tpos-strip->tpos) );

      }

    }else if(strip->planeview==PlaneView::kV){
      if(Main_strip_v->plane-strip->plane>0){
        Strip_ph_meu0*=TMath::Sqrt(2*TMath::Abs(Main_strip_v->plane-strip->plane)*0.0354) + 40*TMath::Sqrt(TMath::Abs(Main_strip_v->tpos-strip->tpos));
        Strip_ph_meu1*=TMath::Sqrt(2*TMath::Abs(Main_strip_v->plane-strip->plane)*0.0354) + 20*TMath::Abs(Main_strip_v->tpos-strip->tpos);
        Strip_ph_meu2*=TMath::Sqrt( (5*TMath::Abs(Main_strip_v->plane-strip->plane)*0.0354) + 40*TMath::Abs(Main_strip_v->tpos-strip->tpos) );
        Strip_ph_meu3*=TMath::Sqrt( (2*TMath::Abs(Main_strip_v->plane-strip->plane)*0.0354) + 40*TMath::Abs(Main_strip_v->tpos-strip->tpos) );

      }else{
        Strip_ph_meu0*=TMath::Sqrt(TMath::Abs(Main_strip_v->plane-strip->plane)*0.0354) + 40*TMath::Sqrt(TMath::Abs(Main_strip_v->tpos-strip->tpos));
        Strip_ph_meu1*=TMath::Sqrt(TMath::Abs(Main_strip_v->plane-strip->plane)*0.0354) + 20*TMath::Abs(Main_strip_v->tpos-strip->tpos);
        Strip_ph_meu2*=TMath::Sqrt( (TMath::Abs(Main_strip_v->plane-strip->plane)*0.0354) + 40*TMath::Abs(Main_strip_v->tpos-strip->tpos) );
        Strip_ph_meu3*=TMath::Sqrt( (TMath::Abs(Main_strip_v->plane-strip->plane)*0.0354) + 40*TMath::Abs(Main_strip_v->tpos-strip->tpos) );
      }

    }

    Weighted_ph0+=Strip_ph_meu0;
    Weighted_ph1+=Strip_ph_meu1;
    Weighted_ph2+=Strip_ph_meu2;
    Weighted_ph3+=Strip_ph_meu3;
  }

  if(MaxStripU==0 && MaxStripV==0){Weighted_ph0=-100;Weighted_ph1=-100;Weighted_ph2=-100;Weighted_ph3=-100;}
  fAngCluster.weightedPH0 = Weighted_ph0;
  fAngCluster.weightedPH1 = Weighted_ph1;
  fAngCluster.weightedPH2 = Weighted_ph2;
  fAngCluster.weightedPH3 = Weighted_ph3;

}


void AngClusterAna::Reset()
{


}
