///////////////////////////////////////////////////////////////////////////
// 
// MuonRemovalInfo
//
///////////////////////////////////////////////////////////////////////////
#ifndef MUONREMOVALINFO_H
#define MUONREMOVALINFO_H

#include "TObject.h"

class MuonRemovalInfo : public TObject
{
  
 public:
  MuonRemovalInfo();
  virtual ~MuonRemovalInfo();
  
  void Reset();
  void Zero();

  //retained digits/strips:
  Int_t   ndigit;
  Int_t   nstrip;
  //new event matching info
  Int_t   orig_event; //original event match
  Float_t best_purity;
  Float_t best_complete;
  Float_t elec_complete;
  Float_t best_purity_phw;
  Float_t best_complete_phw;
  Float_t elec_complete_phw;
  //original event info:
  Float_t orig_cc_pid;  //dpCCPID
//  Float_t orig_dpCCPID;
  Double_t orig_nsCCPID;
  Double_t orig_roCCPID;
  Double_t orig_abCCPID;

  Int_t   nrmstp;
  //removed track info:
  Float_t vtxx;
  Float_t vtxy;
  Float_t vtxz;
  Float_t vtxp;
  Int_t   npln;
  Float_t prng;
  Float_t pcrv;
  Float_t pvdx;
  Float_t pvdy;
  Float_t pvdz;
  Int_t   fitp;
  Int_t   endc;
  Int_t   pass;
  Float_t pmux;
  Float_t pmuy;
  Float_t pmuz;
  Float_t shwe;
  Int_t   mxpl;
  Float_t qp;
  Float_t SigmaQP;
  Float_t vtxdistance;  // track vetext distance to the edge
  Float_t endx;  //track end x
  Float_t endy;  //track end y
  Float_t endz;  //track end z
  Float_t enddistance;  // track end distance to the edge
  Int_t endp;  //track end plane
  Float_t zenith;  // Zenith angle  
  Float_t azimuth; // azi


  Float_t mrmpmux;
  Float_t mrmpmuy;
  Float_t mrmpmuz;
  Float_t mrmQ2;
  Float_t mrmEshw;

  //original shower info:
  UShort_t origShwPlanes;
  UShort_t origShwBegPlane;
  UShort_t origShwEndPlane;
  Int_t origShwStrips;
  Int_t origShwVtxPlane;
  Float_t origShwVtxX;
  Float_t origShwVtxY;
  Float_t origShwVtxZ;
  
  //original event purity/completeness:
  Float_t orig_evt_purity;
  Float_t orig_evt_complete;
  //truth digit info:
  Int_t   nMuonDig;
  Int_t   nMuonDigRetained;
  Int_t   nShwDig;
  Int_t   nShwDigRetained;
  Int_t   nShwDigAtVtx;
  Int_t   nShwDigRetainedAtVtx;
  Float_t nShwPE;
  Float_t nShwPERetained;
  Float_t nShwPEAtVtx;
  Float_t nShwPERetainedAtVtx;
  Int_t   nRetained;
  Int_t   nRetainedMuon;
  Int_t   nRetainedShw;
  Int_t   nRetainedBoth;
  Float_t nPERetained;
  Float_t nPERetainedMuon;
  Float_t nPERetainedShw;
  Float_t nPERetainedBoth;
  Int_t   nRejected;
  Int_t   nRejectedMuon;
  Int_t   nRejectedShw;
  Int_t   nRejectedBoth;
  Int_t   nRejShw;
  Int_t   nRejShwMaxTrk;
  Int_t   nRejShwFakeTrk;
  Int_t   nRejShwMix;

  ClassDef(MuonRemovalInfo,8)
};

#endif// MUONREMOVALINFO_H
