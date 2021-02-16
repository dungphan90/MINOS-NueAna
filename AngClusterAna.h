/// $Id: AngClusterAna.h,v 1.8 2007/05/18 21:01:09 danche Exp $
///
/// class AngCluster
///
/// NueAna package
///
/// Purpose: Cluster 3D Hits using Cos(x)vsCos(y) representations.
///
/// Author: Alex Sousa <asousa@minos.phy.tufts.edu>
///
/// Created on: Thu Apr 28 2005


#ifndef ANGCLUSTERANA_H
#define ANGCLUSTERANA_H

#include <deque>
#include <vector>
#include "TObject.h"
#include "TMatrixD.h"
#include "TVector3.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/AngCluster.h"
#include "Conventions/PlaneView.h"

class NtpSRRecord;
class NtpSREvent;

class AngClusterAna : public NueAnaBase
{

public:

    // Typedefs an enumerations:
    typedef std::deque<Int_t> DeqInt_t;
    typedef DeqInt_t::iterator IterDeqInt_t;
    typedef DeqInt_t::size_type SizeDeqInt_t;

    typedef std::deque<Float_t> DeqFloat_t;
    typedef DeqFloat_t::iterator IterDeqFloat_t;

    typedef std::vector<Int_t> VecInt_t;
    typedef VecInt_t::iterator IterVecInt_t;

    typedef std::vector<Float_t> VecFloat_t;
    typedef VecFloat_t::iterator IterVecFloat_t;

    typedef std::deque< std::deque <Int_t> > DeqDeqInt_t;
    typedef DeqDeqInt_t::iterator IterDeqDeqInt_t;
    typedef DeqDeqInt_t::size_type SizeDeqDeqInt_t;

    typedef std::deque <TMatrixD> DeqTMatrixD;
    typedef DeqTMatrixD::iterator IterDeqTMatrixD;


    // Con/de-structors:
    AngClusterAna(AngCluster &ac);
    ~AngClusterAna();

    void Set3DHit(DeqFloat_t &x, DeqFloat_t &y, DeqFloat_t &z, DeqFloat_t &e);

    void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);

    void WeightedEnergy(int evtn, RecRecordImp<RecCandHeader> *srobj);

    void Reset();

    void GetAngCluster(Int_t &primShow
                       , DeqDeqInt_t &clusterMap
                       , TVector3 &primDir);

    //Member Variables
    VecInt_t fClusterID; //Vector that assigns each hit to a cluster
    DeqInt_t fClusterSize; //Contains number of hits in each cluster

    Int_t fNCluster; //Total number of clusters
    Int_t fPrimShow; //Largest cluster index

    DeqDeqInt_t fClusterMap; // Contains list of hits for each cluster

    TVector3 fPrimDir;  //Primary Shower Direction

    Float_t fTotalEventHitEnergy; //Total event hit energy

private:

    AngCluster &fAngCluster;

    // copy constructor, assignment:
    AngClusterAna(const AngClusterAna& rhs); // copy constructor
    AngClusterAna& operator=(const AngClusterAna& rhs); // assignment

    DeqFloat_t centroidX;
    DeqFloat_t centroidY;

    //To be set to the referenced HitCalcAna Hit member vectors.
    DeqFloat_t fX; 
    DeqFloat_t fY; 
    DeqFloat_t fZ; 
    DeqFloat_t fE; 

    void FindCluster();    // Find Clusters in hit views
    
    void ClusterVarCalc(); // Calculate cluster 
                                           //related variables
    
    void FindCentroidBlob(TMatrixD & grid
                          ,const Int_t nBinX
                          ,const Int_t nBinY
                          ,Int_t indX
                          ,Int_t indY
                          ,TMatrixD & blobCentro); //Find accumulations 
                                                   //in Centroid Histo. 
    
};                              // end of class AngClusterAna

#endif  // ANGCLUSTERANA_H
