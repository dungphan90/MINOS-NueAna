#include "NueAna/MuonRemovalInfo.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

ClassImp(MuonRemovalInfo)

MuonRemovalInfo::MuonRemovalInfo():
  ndigit(ANtpDefVal::kInt),nstrip(ANtpDefVal::kInt),
  orig_event(ANtpDefVal::kInt),best_purity(ANtpDefVal::kDouble),
  best_complete(ANtpDefVal::kDouble),
  elec_complete(ANtpDefVal::kDouble),
  best_purity_phw(ANtpDefVal::kDouble),
  best_complete_phw(ANtpDefVal::kDouble),
  elec_complete_phw(ANtpDefVal::kDouble),
  orig_cc_pid(ANtpDefVal::kDouble),
  orig_nsCCPID(ANtpDefVal::kDouble),
  orig_roCCPID(ANtpDefVal::kDouble),
  orig_abCCPID(ANtpDefVal::kDouble),
  nrmstp(ANtpDefVal::kInt), 
  vtxx(ANtpDefVal::kDouble),
  vtxy(ANtpDefVal::kDouble),vtxz(ANtpDefVal::kDouble),
  vtxp(ANtpDefVal::kDouble),npln(ANtpDefVal::kInt),
  prng(ANtpDefVal::kDouble),pcrv(ANtpDefVal::kDouble),
  pvdx(ANtpDefVal::kDouble),pvdy(ANtpDefVal::kDouble),
  pvdz(ANtpDefVal::kDouble),fitp(ANtpDefVal::kInt),
  endc(ANtpDefVal::kInt),pass(ANtpDefVal::kInt),
  pmux(ANtpDefVal::kDouble),pmuy(ANtpDefVal::kDouble),
  pmuz(ANtpDefVal::kDouble),shwe(ANtpDefVal::kDouble),mxpl(ANtpDefVal::kInt), 
  qp(ANtpDefVal::kDouble),
  SigmaQP(ANtpDefVal::kDouble),
  vtxdistance(ANtpDefVal::kDouble), 
  endx(ANtpDefVal::kDouble), 
  endy(ANtpDefVal::kDouble),
  endz(ANtpDefVal::kDouble),
  enddistance(ANtpDefVal::kDouble), 
  endp(ANtpDefVal::kInt),
  zenith(ANtpDefVal::kDouble),
  azimuth(ANtpDefVal::kDouble),
  mrmpmux(ANtpDefVal::kDouble),mrmpmuy(ANtpDefVal::kDouble),
  mrmpmuz(ANtpDefVal::kDouble),mrmQ2(ANtpDefVal::kDouble),
  mrmEshw(ANtpDefVal::kInt), 
  origShwPlanes(ANtpDefVal::kInt),
  origShwBegPlane(ANtpDefVal::kInt),
  origShwEndPlane(ANtpDefVal::kInt),
  origShwStrips(ANtpDefVal::kInt),
  origShwVtxPlane(ANtpDefVal::kInt),
  origShwVtxX(ANtpDefVal::kFloat),
  origShwVtxY(ANtpDefVal::kFloat),
  origShwVtxZ(ANtpDefVal::kFloat),
  orig_evt_purity(ANtpDefVal::kDouble),orig_evt_complete(ANtpDefVal::kDouble),
  nMuonDig(ANtpDefVal::kInt),
  nMuonDigRetained(ANtpDefVal::kInt),nShwDig(ANtpDefVal::kInt),
  nShwDigRetained(ANtpDefVal::kInt),nShwDigAtVtx(ANtpDefVal::kInt),
  nShwDigRetainedAtVtx(ANtpDefVal::kInt),nShwPE(ANtpDefVal::kDouble),
  nShwPERetained(ANtpDefVal::kDouble),nShwPEAtVtx(ANtpDefVal::kDouble),
  nShwPERetainedAtVtx(ANtpDefVal::kDouble),nRetained(ANtpDefVal::kInt),
  nRetainedMuon(ANtpDefVal::kInt),nRetainedShw(ANtpDefVal::kInt),
  nRetainedBoth(ANtpDefVal::kInt),nPERetained(ANtpDefVal::kDouble),
  nPERetainedMuon(ANtpDefVal::kDouble),nPERetainedShw(ANtpDefVal::kDouble),
  nPERetainedBoth(ANtpDefVal::kDouble),nRejected(ANtpDefVal::kInt),
  nRejectedMuon(ANtpDefVal::kInt),nRejectedShw(ANtpDefVal::kInt),
  nRejectedBoth(ANtpDefVal::kInt),nRejShw(ANtpDefVal::kInt),
  nRejShwMaxTrk(ANtpDefVal::kInt),nRejShwFakeTrk(ANtpDefVal::kInt),
  nRejShwMix(ANtpDefVal::kInt)
{}

MuonRemovalInfo::~MuonRemovalInfo(){}

void MuonRemovalInfo::Reset()
{
  ndigit = ANtpDefVal::kInt;
  nstrip = ANtpDefVal::kInt;
  orig_event = ANtpDefVal::kInt;
  best_purity = ANtpDefVal::kDouble;
  best_complete = ANtpDefVal::kDouble;
  elec_complete = ANtpDefVal::kDouble;
  best_purity_phw = ANtpDefVal::kDouble;
  best_complete_phw = ANtpDefVal::kDouble;
  elec_complete_phw = ANtpDefVal::kDouble;
  orig_cc_pid = ANtpDefVal::kDouble;
  orig_nsCCPID = ANtpDefVal::kDouble;
  orig_roCCPID = ANtpDefVal::kDouble;
  orig_abCCPID = ANtpDefVal::kDouble; 
  nrmstp = ANtpDefVal::kInt;
  vtxx = ANtpDefVal::kDouble;
  vtxy = ANtpDefVal::kDouble;
  vtxz = ANtpDefVal::kDouble;
  vtxp = ANtpDefVal::kDouble;
  npln = ANtpDefVal::kInt;
  prng = ANtpDefVal::kDouble;
  pcrv = ANtpDefVal::kDouble;
  pvdx = ANtpDefVal::kDouble;
  pvdy = ANtpDefVal::kDouble;
  pvdz = ANtpDefVal::kDouble;
  fitp = ANtpDefVal::kInt;
  endc = ANtpDefVal::kInt;
  pass = ANtpDefVal::kInt;
  pmux = ANtpDefVal::kDouble;
  pmuy = ANtpDefVal::kDouble;
  pmuz = ANtpDefVal::kDouble;
  shwe = ANtpDefVal::kDouble;
  mxpl = ANtpDefVal::kInt;
  qp = ANtpDefVal::kDouble;
  SigmaQP = ANtpDefVal::kDouble;
  vtxdistance = ANtpDefVal::kDouble; 
  endx = ANtpDefVal::kDouble;
  endy = ANtpDefVal::kDouble;
  endz = ANtpDefVal::kDouble;
  enddistance = ANtpDefVal::kDouble; 
  endp = ANtpDefVal::kInt;
  zenith = ANtpDefVal::kDouble;
  azimuth = ANtpDefVal::kDouble;
  mrmpmux = ANtpDefVal::kDouble;
  mrmpmuy = ANtpDefVal::kDouble;
  mrmpmuz = ANtpDefVal::kDouble;
  mrmQ2 = ANtpDefVal::kDouble;
  mrmEshw = ANtpDefVal::kDouble;
  origShwPlanes = ANtpDefVal::kInt;
  origShwBegPlane = ANtpDefVal::kInt;
  origShwEndPlane = ANtpDefVal::kInt;
  origShwStrips = ANtpDefVal::kInt;
  origShwVtxPlane = ANtpDefVal::kInt;
  origShwVtxX = ANtpDefVal::kFloat;
  origShwVtxY = ANtpDefVal::kFloat;
  origShwVtxZ = ANtpDefVal::kFloat;
  orig_evt_purity = ANtpDefVal::kDouble;
  orig_evt_complete = ANtpDefVal::kDouble;
  nMuonDig = ANtpDefVal::kInt;
  nMuonDigRetained = ANtpDefVal::kInt;
  nShwDig = ANtpDefVal::kInt;
  nShwDigRetained = ANtpDefVal::kInt;
  nShwDigAtVtx = ANtpDefVal::kInt;
  nShwDigRetainedAtVtx = ANtpDefVal::kInt;
  nShwPE = ANtpDefVal::kDouble;
  nShwPERetained = ANtpDefVal::kDouble;
  nShwPEAtVtx = ANtpDefVal::kDouble;
  nShwPERetainedAtVtx = ANtpDefVal::kDouble;
  nRetained = ANtpDefVal::kInt;
  nRetainedMuon = ANtpDefVal::kInt;
  nRetainedShw = ANtpDefVal::kInt;
  nRetainedBoth = ANtpDefVal::kInt;
  nPERetained = ANtpDefVal::kDouble;
  nPERetainedMuon = ANtpDefVal::kDouble;
  nPERetainedShw = ANtpDefVal::kDouble;
  nPERetainedBoth = ANtpDefVal::kDouble;
  nRejected = ANtpDefVal::kInt;
  nRejectedMuon = ANtpDefVal::kInt;
  nRejectedShw = ANtpDefVal::kInt;
  nRejectedBoth = ANtpDefVal::kInt;
  nRejShw = ANtpDefVal::kInt;
  nRejShwMaxTrk = ANtpDefVal::kInt;
  nRejShwFakeTrk = ANtpDefVal::kInt;
  nRejShwMix = ANtpDefVal::kInt;
}

void MuonRemovalInfo::Zero()
{
  ndigit = 0;
  nstrip = 0;
  orig_event = -1;
  best_purity = 0;
  best_complete = 0;
  elec_complete = 0;
  best_purity_phw = 0;
  best_complete_phw = 0;
  elec_complete_phw = 0;
  orig_cc_pid = 0;
  orig_nsCCPID = 0;
  orig_roCCPID = 0;
  orig_abCCPID = 0;  
  qp = 0;
  SigmaQP = 0;
  vtxdistance = 0; 
  endx = 0; 
  endy = 0;
  endz = 0;
  enddistance = 0; 
  endp = 0;
  zenith = 0;
  azimuth = 0;
  origShwPlanes = 0;
  origShwBegPlane = 0;
  origShwEndPlane = 0;
  origShwStrips = 0;
  origShwVtxPlane = 0;
  origShwVtxX = 0;
  origShwVtxY = 0;
  origShwVtxZ = 0;
  nrmstp = 0;
  vtxx = 0;
  vtxy = 0;
  vtxz = 0;
  vtxp = 0;
  npln = 0;
  prng = 0;
  pcrv = 0;
  pvdx = 0;
  pvdy = 0;
  pvdz = 0;
  fitp = 0;
  endc = 0;
  pass = 0;
  pmux = 0;
  pmuy = 0;
  pmuz = 0;
  shwe = 0;
  mxpl = 0;
  mrmpmux = 0;
  mrmpmuy = 0;
  mrmpmuz = 0;
  mrmQ2 = 0;
  mrmEshw = 0;
  orig_evt_purity = 0;
  orig_evt_complete = 0;
  nMuonDig = 0;
  nMuonDigRetained = 0;
  nShwDig = 0;
  nShwDigRetained = 0;
  nShwDigAtVtx = 0;
  nShwDigRetainedAtVtx = 0;
  nShwPE = 0;
  nShwPERetained = 0;
  nShwPEAtVtx = 0;
  nShwPERetainedAtVtx = 0;
  nRetained = 0;
  nRetainedMuon = 0;
  nRetainedShw = 0;
  nRetainedBoth = 0;
  nPERetained = 0;
  nPERetainedMuon = 0;
  nPERetainedShw = 0;
  nPERetainedBoth = 0;
  nRejected = 0;
  nRejectedMuon = 0;
  nRejectedShw = 0;
  nRejectedBoth = 0;
  nRejShw = 0;
  nRejShwMaxTrk = 0;
  nRejShwFakeTrk = 0;
  nRejShwMix = 0;
}
