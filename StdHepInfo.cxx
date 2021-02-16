#include "StdHepInfo.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

ClassImp(StdHepInfo)

StdHepInfo::StdHepInfo():
  nmult(ANtpDefVal::kInt),  nhad(ANtpDefVal::kInt),
  nfs(ANtpDefVal::kInt),nlep(ANtpDefVal::kInt),
  npip(ANtpDefVal::kInt),npim(ANtpDefVal::kInt),
  npi0(ANtpDefVal::kInt),nprot(ANtpDefVal::kInt),
  nneut(ANtpDefVal::kInt),nkaon(ANtpDefVal::kInt),
  ngeant(ANtpDefVal::kInt),
  nfs_he(ANtpDefVal::kInt),nlep_he(ANtpDefVal::kInt),
  npip_he(ANtpDefVal::kInt),npim_he(ANtpDefVal::kInt),
  npi0_he(ANtpDefVal::kInt),nprot_he(ANtpDefVal::kInt),
  nneut_he(ANtpDefVal::kInt),nkaon_he(ANtpDefVal::kInt),
  ngeant_he(ANtpDefVal::kInt),
  etot(ANtpDefVal::kDouble),elep(ANtpDefVal::kDouble),
  epip(ANtpDefVal::kDouble),epim(ANtpDefVal::kDouble),
  epi0(ANtpDefVal::kDouble),eprot(ANtpDefVal::kDouble),
  eneut(ANtpDefVal::kDouble),ekaon(ANtpDefVal::kDouble),
  egeant(ANtpDefVal::kDouble),
  pztot(ANtpDefVal::kDouble),pzlep(ANtpDefVal::kDouble),
  pzpip(ANtpDefVal::kDouble),pzpim(ANtpDefVal::kDouble),
  pzpi0(ANtpDefVal::kDouble),pzprot(ANtpDefVal::kDouble),
  pzneut(ANtpDefVal::kDouble),pzkaon(ANtpDefVal::kDouble),
  pzgeant(ANtpDefVal::kDouble),
  pttot(ANtpDefVal::kDouble),ptlep(ANtpDefVal::kDouble),
  ptpip(ANtpDefVal::kDouble),ptpim(ANtpDefVal::kDouble),
  ptpi0(ANtpDefVal::kDouble),ptprot(ANtpDefVal::kDouble),
  ptneut(ANtpDefVal::kDouble),ptkaon(ANtpDefVal::kDouble),
  ptgeant(ANtpDefVal::kDouble),
  wfs(ANtpDefVal::kDouble),wlep(ANtpDefVal::kDouble),
  wpip(ANtpDefVal::kDouble),wpim(ANtpDefVal::kDouble),
  wpi0(ANtpDefVal::kDouble),wprot(ANtpDefVal::kDouble),
  wneut(ANtpDefVal::kDouble),wkaon(ANtpDefVal::kDouble),
  wgeant(ANtpDefVal::kDouble),//finalID(ANtpDefVal::kDouble),
  lepnmult(ANtpDefVal::kInt),
  lepnfs(ANtpDefVal::kInt),lepnlep(ANtpDefVal::kInt),
  lepnpip(ANtpDefVal::kInt),lepnpim(ANtpDefVal::kInt),
  lepnpi0(ANtpDefVal::kInt),lepnprot(ANtpDefVal::kInt),
  lepnneut(ANtpDefVal::kInt),lepnkaon(ANtpDefVal::kInt),
  lepngeant(ANtpDefVal::kInt),
  lepnfs_he(ANtpDefVal::kInt),lepnlep_he(ANtpDefVal::kInt),
  lepnpip_he(ANtpDefVal::kInt),lepnpim_he(ANtpDefVal::kInt),
  lepnpi0_he(ANtpDefVal::kInt),lepnprot_he(ANtpDefVal::kInt),
  lepnneut_he(ANtpDefVal::kInt),lepnkaon_he(ANtpDefVal::kInt),
  lepngeant_he(ANtpDefVal::kInt),
  lepetot(ANtpDefVal::kDouble),lepelep(ANtpDefVal::kDouble),
  lepepip(ANtpDefVal::kDouble),lepepim(ANtpDefVal::kDouble),
  lepepi0(ANtpDefVal::kDouble),lepeprot(ANtpDefVal::kDouble),
  lepeneut(ANtpDefVal::kDouble),lepekaon(ANtpDefVal::kDouble),
  lepegeant(ANtpDefVal::kDouble),
  leppztot(ANtpDefVal::kDouble),leppzlep(ANtpDefVal::kDouble),
  leppzpip(ANtpDefVal::kDouble),leppzpim(ANtpDefVal::kDouble),
  leppzpi0(ANtpDefVal::kDouble),leppzprot(ANtpDefVal::kDouble),
  leppzneut(ANtpDefVal::kDouble),leppzkaon(ANtpDefVal::kDouble),
  leppzgeant(ANtpDefVal::kDouble),
  leppttot(ANtpDefVal::kDouble),lepptlep(ANtpDefVal::kDouble),
  lepptpip(ANtpDefVal::kDouble),lepptpim(ANtpDefVal::kDouble),
  lepptpi0(ANtpDefVal::kDouble),lepptprot(ANtpDefVal::kDouble),
  lepptneut(ANtpDefVal::kDouble),lepptkaon(ANtpDefVal::kDouble),
  lepptgeant(ANtpDefVal::kDouble),
  lepwfs(ANtpDefVal::kDouble),lepwlep(ANtpDefVal::kDouble),
  lepwpip(ANtpDefVal::kDouble),lepwpim(ANtpDefVal::kDouble),
  lepwpi0(ANtpDefVal::kDouble),lepwprot(ANtpDefVal::kDouble),
  lepwneut(ANtpDefVal::kDouble),lepwkaon(ANtpDefVal::kDouble),
  lepwgeant(ANtpDefVal::kDouble),lep2leptype(ANtpDefVal::kDouble),
  emfrac(ANtpDefVal::kDouble),emcount(ANtpDefVal::kInt),
  emenergy(ANtpDefVal::kDouble),epi0_total(ANtpDefVal::kDouble),
  epi0_intranuke(ANtpDefVal::kDouble),epi0_neugen(ANtpDefVal::kDouble),
  epi0_decay(ANtpDefVal::kDouble),epi0_abs(ANtpDefVal::kDouble),
  npi0_neugen(ANtpDefVal::kInt),baryonpt(ANtpDefVal::kDouble),
  baryonxf(ANtpDefVal::kDouble),totpt(ANtpDefVal::kDouble),
  e_total(ANtpDefVal::kDouble),
  em_total(ANtpDefVal::kDouble),
  npi0_jbhg(ANtpDefVal::kInt),
  nch_jbhg(ANtpDefVal::kInt),
  ntot_jbhg(ANtpDefVal::kInt),

  nfs_nu(ANtpDefVal::kInt), 
  nfs_nu_he(ANtpDefVal::kInt),
  etot_nu(ANtpDefVal::kDouble),
  pztot_nu(ANtpDefVal::kDouble),  
  pttot_nu(ANtpDefVal::kDouble),
  wfs_nu(ANtpDefVal::kDouble),
  lepnfs_nu(ANtpDefVal::kInt),
  lepnfs_nu_he(ANtpDefVal::kInt),
  lepetot_nu(ANtpDefVal::kDouble),
  leppztot_nu(ANtpDefVal::kDouble),  
  leppttot_nu(ANtpDefVal::kDouble),
  lepwfs_nu(ANtpDefVal::kDouble),
  n_nu_elec(ANtpDefVal::kInt),
  n_nu_elec_he(ANtpDefVal::kInt),
  e_nu_elec(ANtpDefVal::kDouble),
  pz_nu_elec(ANtpDefVal::kDouble),  
  pt_nu_elec(ANtpDefVal::kDouble),
  w_nu_elec(ANtpDefVal::kDouble),
  lepn_nu_elec(ANtpDefVal::kInt),
  lepn_nu_elec_he(ANtpDefVal::kInt),
  lepe_nu_elec(ANtpDefVal::kDouble),
  leppz_nu_elec(ANtpDefVal::kDouble),  
  leppt_nu_elec(ANtpDefVal::kDouble),
  lepw_nu_elec(ANtpDefVal::kDouble),
  n_nu_muon(ANtpDefVal::kInt),
  n_nu_muon_he(ANtpDefVal::kInt),
  e_nu_muon(ANtpDefVal::kDouble),
  pz_nu_muon(ANtpDefVal::kDouble),  
  pt_nu_muon(ANtpDefVal::kDouble),
  w_nu_muon(ANtpDefVal::kDouble),
  lepn_nu_muon(ANtpDefVal::kInt),
  lepn_nu_muon_he(ANtpDefVal::kInt),
  lepe_nu_muon(ANtpDefVal::kDouble),
  leppz_nu_muon(ANtpDefVal::kDouble),  
  leppt_nu_muon(ANtpDefVal::kDouble),
  lepw_nu_muon(ANtpDefVal::kDouble),
  n_nu_tau(ANtpDefVal::kInt),
  n_nu_tau_he(ANtpDefVal::kInt),
  e_nu_tau(ANtpDefVal::kDouble),
  pz_nu_tau(ANtpDefVal::kDouble),  
  pt_nu_tau(ANtpDefVal::kDouble),
  w_nu_tau(ANtpDefVal::kDouble),
  lepn_nu_tau(ANtpDefVal::kInt),
  lepn_nu_tau_he(ANtpDefVal::kInt),
  lepe_nu_tau(ANtpDefVal::kDouble),
  leppz_nu_tau(ANtpDefVal::kDouble),  
  leppt_nu_tau(ANtpDefVal::kDouble),
  lepw_nu_tau(ANtpDefVal::kDouble),
  nlep_elec(ANtpDefVal::kInt),
  nlep_elec_he(ANtpDefVal::kInt),
  elep_elec(ANtpDefVal::kDouble),
  pzlep_elec(ANtpDefVal::kDouble),  
  ptlep_elec(ANtpDefVal::kDouble),
  wlep_elec(ANtpDefVal::kDouble),
  lepnlep_elec(ANtpDefVal::kInt),
  lepnlep_elec_he(ANtpDefVal::kInt),
  lepelep_elec(ANtpDefVal::kDouble),
  leppzlep_elec(ANtpDefVal::kDouble),  
  lepptlep_elec(ANtpDefVal::kDouble),
  lepwlep_elec(ANtpDefVal::kDouble),
  nlep_muon(ANtpDefVal::kInt),
  nlep_muon_he(ANtpDefVal::kInt),
  elep_muon(ANtpDefVal::kDouble),
  pzlep_muon(ANtpDefVal::kDouble),  
  ptlep_muon(ANtpDefVal::kDouble),
  wlep_muon(ANtpDefVal::kDouble),
  lepnlep_muon(ANtpDefVal::kInt),
  lepnlep_muon_he(ANtpDefVal::kInt),
  lepelep_muon(ANtpDefVal::kDouble),
  leppzlep_muon(ANtpDefVal::kDouble),  
  lepptlep_muon(ANtpDefVal::kDouble),
  lepwlep_muon(ANtpDefVal::kDouble),
  nlep_tau(ANtpDefVal::kInt),
  nlep_tau_he(ANtpDefVal::kInt),
  elep_tau(ANtpDefVal::kDouble),
  pzlep_tau(ANtpDefVal::kDouble),  
  ptlep_tau(ANtpDefVal::kDouble),
  wlep_tau(ANtpDefVal::kDouble),
  lepnlep_tau(ANtpDefVal::kInt),
  lepnlep_tau_he(ANtpDefVal::kInt),
  lepelep_tau(ANtpDefVal::kDouble),
  leppzlep_tau(ANtpDefVal::kDouble),  
  lepptlep_tau(ANtpDefVal::kDouble),
  lepwlep_tau(ANtpDefVal::kDouble),
  nkaonplus(ANtpDefVal::kInt),
  nkaonplus_he(ANtpDefVal::kInt),
  ekaonplus(ANtpDefVal::kDouble),
  pzkaonplus(ANtpDefVal::kDouble),  
  ptkaonplus(ANtpDefVal::kDouble),
  wkaonplus(ANtpDefVal::kDouble),
  lepnkaonplus(ANtpDefVal::kInt),
  lepnkaonplus_he(ANtpDefVal::kInt),
  lepekaonplus(ANtpDefVal::kDouble),
  leppzkaonplus(ANtpDefVal::kDouble),  
  lepptkaonplus(ANtpDefVal::kDouble),
  lepwkaonplus(ANtpDefVal::kDouble),
  nkaonminus(ANtpDefVal::kInt),
  nkaonminus_he(ANtpDefVal::kInt),
  ekaonminus(ANtpDefVal::kDouble),
  pzkaonminus(ANtpDefVal::kDouble),  
  ptkaonminus(ANtpDefVal::kDouble),
  wkaonminus(ANtpDefVal::kDouble),
  lepnkaonminus(ANtpDefVal::kInt),
  lepnkaonminus_he(ANtpDefVal::kInt),
  lepekaonminus(ANtpDefVal::kDouble),
  leppzkaonminus(ANtpDefVal::kDouble),  
  lepptkaonminus(ANtpDefVal::kDouble),
  lepwkaonminus(ANtpDefVal::kDouble),
  nkaon0L(ANtpDefVal::kInt),
  nkaon0L_he(ANtpDefVal::kInt),
  ekaon0L(ANtpDefVal::kDouble),
  pzkaon0L(ANtpDefVal::kDouble),  
  ptkaon0L(ANtpDefVal::kDouble),
  wkaon0L(ANtpDefVal::kDouble),
  lepnkaon0L(ANtpDefVal::kInt),
  lepnkaon0L_he(ANtpDefVal::kInt),
  lepekaon0L(ANtpDefVal::kDouble),
  leppzkaon0L(ANtpDefVal::kDouble),  
  lepptkaon0L(ANtpDefVal::kDouble),
  lepwkaon0L(ANtpDefVal::kDouble),
  nkaon0S(ANtpDefVal::kInt),
  nkaon0S_he(ANtpDefVal::kInt),
  ekaon0S(ANtpDefVal::kDouble),
  pzkaon0S(ANtpDefVal::kDouble),  
  ptkaon0S(ANtpDefVal::kDouble),
  wkaon0S(ANtpDefVal::kDouble),
  lepnkaon0S(ANtpDefVal::kInt),
  lepnkaon0S_he(ANtpDefVal::kInt),
  lepekaon0S(ANtpDefVal::kDouble),
  leppzkaon0S(ANtpDefVal::kDouble),  
  lepptkaon0S(ANtpDefVal::kDouble),
  lepwkaon0S(ANtpDefVal::kDouble),
  nkaon0(ANtpDefVal::kInt),
  nkaon0_he(ANtpDefVal::kInt),
  ekaon0(ANtpDefVal::kDouble),
  pzkaon0(ANtpDefVal::kDouble),  
  ptkaon0(ANtpDefVal::kDouble),
  wkaon0(ANtpDefVal::kDouble),
  lepnkaon0(ANtpDefVal::kInt),
  lepnkaon0_he(ANtpDefVal::kInt),
  lepekaon0(ANtpDefVal::kDouble),
  leppzkaon0(ANtpDefVal::kDouble),  
  lepptkaon0(ANtpDefVal::kDouble),
  lepwkaon0(ANtpDefVal::kDouble),
  nutauevent(ANtpDefVal::kInt),
  nutauchannel(ANtpDefVal::kInt)
{}

StdHepInfo::~StdHepInfo(){}

void StdHepInfo::Reset()
{
  nmult   = ANtpDefVal::kInt;
  nhad    = ANtpDefVal::kInt;
  nfs_he  = ANtpDefVal::kInt;
  nfs     = ANtpDefVal::kInt;  
  nlep    = ANtpDefVal::kInt;
  npip    = ANtpDefVal::kInt;
  npim    = ANtpDefVal::kInt;
  npi0    = ANtpDefVal::kInt;
  nprot   = ANtpDefVal::kInt;
  nneut   = ANtpDefVal::kInt;
  nkaon   = ANtpDefVal::kInt;
  ngeant  = ANtpDefVal::kInt;
  etot    = ANtpDefVal::kDouble;
  elep    = ANtpDefVal::kDouble;
  epip    = ANtpDefVal::kDouble;
  epim    = ANtpDefVal::kDouble;
  epi0    = ANtpDefVal::kDouble;
  eprot   = ANtpDefVal::kDouble;
  eneut   = ANtpDefVal::kDouble;
  ekaon   = ANtpDefVal::kDouble;
  egeant  = ANtpDefVal::kDouble;
  pztot   = ANtpDefVal::kDouble;
  pzlep   = ANtpDefVal::kDouble;
  pzpip   = ANtpDefVal::kDouble;
  pzpim   = ANtpDefVal::kDouble;
  pzpi0   = ANtpDefVal::kDouble;
  pzprot  = ANtpDefVal::kDouble;
  pzneut  = ANtpDefVal::kDouble;
  pzkaon  = ANtpDefVal::kDouble;
  pzgeant = ANtpDefVal::kDouble;
  pttot   = ANtpDefVal::kDouble;
  ptlep   = ANtpDefVal::kDouble;
  ptpip   = ANtpDefVal::kDouble;
  ptpim   = ANtpDefVal::kDouble;
  ptpi0   = ANtpDefVal::kDouble;
  ptprot  = ANtpDefVal::kDouble;
  ptneut  = ANtpDefVal::kDouble;
  ptkaon  = ANtpDefVal::kDouble;
  ptgeant = ANtpDefVal::kDouble;
  wfs     = ANtpDefVal::kDouble;
  wlep    = ANtpDefVal::kDouble;
  wpip    = ANtpDefVal::kDouble;
  wpim    = ANtpDefVal::kDouble;
  wpi0    = ANtpDefVal::kDouble;
  wprot   = ANtpDefVal::kDouble;
  wneut   = ANtpDefVal::kDouble;
  wkaon   = ANtpDefVal::kDouble;
  wgeant  = ANtpDefVal::kDouble;
  //////// finalID = ANtpDefVal::kDouble;
  lepnmult   = ANtpDefVal::kInt;
  lepnfs_he  = ANtpDefVal::kInt;
  lepnfs     = ANtpDefVal::kInt;  
  lepnlep    = ANtpDefVal::kInt;
  lepnpip    = ANtpDefVal::kInt;
  lepnpim    = ANtpDefVal::kInt;
  lepnpi0    = ANtpDefVal::kInt;
  lepnprot   = ANtpDefVal::kInt;
  lepnneut   = ANtpDefVal::kInt;
  lepnkaon   = ANtpDefVal::kInt;
  lepngeant  = ANtpDefVal::kInt;
  lepetot    = ANtpDefVal::kDouble;
  lepelep    = ANtpDefVal::kDouble;
  lepepip    = ANtpDefVal::kDouble;
  lepepim    = ANtpDefVal::kDouble;
  lepepi0    = ANtpDefVal::kDouble;
  lepeprot   = ANtpDefVal::kDouble;
  lepeneut   = ANtpDefVal::kDouble;
  lepekaon   = ANtpDefVal::kDouble;
  lepegeant  = ANtpDefVal::kDouble;
  leppztot   = ANtpDefVal::kDouble;
  leppzlep   = ANtpDefVal::kDouble;
  leppzpip   = ANtpDefVal::kDouble;
  leppzpim   = ANtpDefVal::kDouble;
  leppzpi0   = ANtpDefVal::kDouble;
  leppzprot  = ANtpDefVal::kDouble;
  leppzneut  = ANtpDefVal::kDouble;
  leppzkaon  = ANtpDefVal::kDouble;
  leppzgeant = ANtpDefVal::kDouble;
  leppttot   = ANtpDefVal::kDouble;
  lepptlep   = ANtpDefVal::kDouble;
  lepptpip   = ANtpDefVal::kDouble;
  lepptpim   = ANtpDefVal::kDouble;
  lepptpi0   = ANtpDefVal::kDouble;
  lepptprot  = ANtpDefVal::kDouble;
  lepptneut  = ANtpDefVal::kDouble;
  lepptkaon  = ANtpDefVal::kDouble;
  lepptgeant = ANtpDefVal::kDouble;
  lepwfs     = ANtpDefVal::kDouble;
  lepwlep    = ANtpDefVal::kDouble;
  lepwpip    = ANtpDefVal::kDouble;
  lepwpim    = ANtpDefVal::kDouble;
  lepwpi0    = ANtpDefVal::kDouble;
  lepwprot   = ANtpDefVal::kDouble;
  lepwneut   = ANtpDefVal::kDouble;
  lepwkaon   = ANtpDefVal::kDouble;
  lepwgeant  = ANtpDefVal::kDouble;
  lep2leptype= ANtpDefVal::kDouble;
  emfrac     = ANtpDefVal::kDouble;
  emenergy   = ANtpDefVal::kDouble;
  emcount    = ANtpDefVal::kInt;
  epi0_total = ANtpDefVal::kDouble;
  epi0_intranuke = ANtpDefVal::kDouble;
  epi0_neugen= ANtpDefVal::kDouble;
  epi0_decay = ANtpDefVal::kDouble;
  epi0_abs   = ANtpDefVal::kDouble;
  npi0_neugen = ANtpDefVal::kInt;
  baryonpt   = ANtpDefVal::kDouble;
  baryonxf   = ANtpDefVal::kDouble;
  totpt   = ANtpDefVal::kDouble;

  e_total = ANtpDefVal::kDouble;
  em_total = ANtpDefVal::kDouble;
  npi0_jbhg = ANtpDefVal::kInt;
  nch_jbhg = ANtpDefVal::kInt;
  ntot_jbhg = ANtpDefVal::kInt;

  //added Anna vars:
  nfs_nu = ANtpDefVal::kInt; 
  nfs_nu_he = ANtpDefVal::kInt;
  etot_nu = ANtpDefVal::kDouble;
  pztot_nu = ANtpDefVal::kDouble;  
  pttot_nu = ANtpDefVal::kDouble;
  wfs_nu = ANtpDefVal::kDouble;
  lepnfs_nu = ANtpDefVal::kInt;
  lepnfs_nu_he = ANtpDefVal::kInt;
  lepetot_nu = ANtpDefVal::kDouble;
  leppztot_nu = ANtpDefVal::kDouble;  
  leppttot_nu = ANtpDefVal::kDouble;
  lepwfs_nu = ANtpDefVal::kDouble;
  n_nu_elec = ANtpDefVal::kInt;
  n_nu_elec_he = ANtpDefVal::kInt;
  e_nu_elec = ANtpDefVal::kDouble;
  pz_nu_elec = ANtpDefVal::kDouble;  
  pt_nu_elec = ANtpDefVal::kDouble;
  w_nu_elec = ANtpDefVal::kDouble;
  lepn_nu_elec = ANtpDefVal::kInt;
  lepn_nu_elec_he = ANtpDefVal::kInt;
  lepe_nu_elec = ANtpDefVal::kDouble;
  leppz_nu_elec = ANtpDefVal::kDouble;  
  leppt_nu_elec = ANtpDefVal::kDouble;
  lepw_nu_elec = ANtpDefVal::kDouble;
  n_nu_muon = ANtpDefVal::kInt;
  n_nu_muon_he = ANtpDefVal::kInt;
  e_nu_muon = ANtpDefVal::kDouble;
  pz_nu_muon = ANtpDefVal::kDouble;  
  pt_nu_muon = ANtpDefVal::kDouble;
  w_nu_muon = ANtpDefVal::kDouble;
  lepn_nu_muon = ANtpDefVal::kInt;
  lepn_nu_muon_he = ANtpDefVal::kInt;
  lepe_nu_muon = ANtpDefVal::kDouble;
  leppz_nu_muon = ANtpDefVal::kDouble;  
  leppt_nu_muon = ANtpDefVal::kDouble;
  lepw_nu_muon = ANtpDefVal::kDouble;
  n_nu_tau = ANtpDefVal::kInt;
  n_nu_tau_he = ANtpDefVal::kInt;
  e_nu_tau = ANtpDefVal::kDouble;
  pz_nu_tau = ANtpDefVal::kDouble;  
  pt_nu_tau = ANtpDefVal::kDouble;
  w_nu_tau = ANtpDefVal::kDouble;
  lepn_nu_tau = ANtpDefVal::kInt;
  lepn_nu_tau_he = ANtpDefVal::kInt;
  lepe_nu_tau = ANtpDefVal::kDouble;
  leppz_nu_tau = ANtpDefVal::kDouble;  
  leppt_nu_tau = ANtpDefVal::kDouble;
  lepw_nu_tau = ANtpDefVal::kDouble;
  nlep_elec = ANtpDefVal::kInt;
  nlep_elec_he = ANtpDefVal::kInt;
  elep_elec = ANtpDefVal::kDouble;
  pzlep_elec = ANtpDefVal::kDouble;  
  ptlep_elec = ANtpDefVal::kDouble;
  wlep_elec = ANtpDefVal::kDouble;
  lepnlep_elec = ANtpDefVal::kInt;
  lepnlep_elec_he = ANtpDefVal::kInt;
  lepelep_elec = ANtpDefVal::kDouble;
  leppzlep_elec = ANtpDefVal::kDouble;  
  lepptlep_elec = ANtpDefVal::kDouble;
  lepwlep_elec = ANtpDefVal::kDouble;
  nlep_muon = ANtpDefVal::kInt;
  nlep_muon_he = ANtpDefVal::kInt;
  elep_muon = ANtpDefVal::kDouble;
  pzlep_muon = ANtpDefVal::kDouble;  
  ptlep_muon = ANtpDefVal::kDouble;
  wlep_muon = ANtpDefVal::kDouble;
  lepnlep_muon = ANtpDefVal::kInt;
  lepnlep_muon_he = ANtpDefVal::kInt;
  lepelep_muon = ANtpDefVal::kDouble;
  leppzlep_muon = ANtpDefVal::kDouble;  
  lepptlep_muon = ANtpDefVal::kDouble;
  lepwlep_muon = ANtpDefVal::kDouble;
  nlep_tau = ANtpDefVal::kInt;
  nlep_tau_he = ANtpDefVal::kInt;
  elep_tau = ANtpDefVal::kDouble;
  pzlep_tau = ANtpDefVal::kDouble;  
  ptlep_tau = ANtpDefVal::kDouble;
  wlep_tau = ANtpDefVal::kDouble;
  lepnlep_tau = ANtpDefVal::kInt;
  lepnlep_tau_he = ANtpDefVal::kInt;
  lepelep_tau = ANtpDefVal::kDouble;
  leppzlep_tau = ANtpDefVal::kDouble;  
  lepptlep_tau = ANtpDefVal::kDouble;
  lepwlep_tau = ANtpDefVal::kDouble;
  nkaonplus = ANtpDefVal::kInt;
  nkaonplus_he = ANtpDefVal::kInt;
  ekaonplus = ANtpDefVal::kDouble;
  pzkaonplus = ANtpDefVal::kDouble;  
  ptkaonplus = ANtpDefVal::kDouble;
  wkaonplus = ANtpDefVal::kDouble;
  lepnkaonplus = ANtpDefVal::kInt;
  lepnkaonplus_he = ANtpDefVal::kInt;
  lepekaonplus = ANtpDefVal::kDouble;
  leppzkaonplus = ANtpDefVal::kDouble;  
  lepptkaonplus = ANtpDefVal::kDouble;
  lepwkaonplus = ANtpDefVal::kDouble;
  nkaonminus = ANtpDefVal::kInt;
  nkaonminus_he = ANtpDefVal::kInt;
  ekaonminus = ANtpDefVal::kDouble;
  pzkaonminus = ANtpDefVal::kDouble;  
  ptkaonminus = ANtpDefVal::kDouble;
  wkaonminus = ANtpDefVal::kDouble;
  lepnkaonminus = ANtpDefVal::kInt;
  lepnkaonminus_he = ANtpDefVal::kInt;
  lepekaonminus = ANtpDefVal::kDouble;
  leppzkaonminus = ANtpDefVal::kDouble;  
  lepptkaonminus = ANtpDefVal::kDouble;
  lepwkaonminus = ANtpDefVal::kDouble;
  nkaon0L = ANtpDefVal::kInt;
  nkaon0L_he = ANtpDefVal::kInt;
  ekaon0L = ANtpDefVal::kDouble;
  pzkaon0L = ANtpDefVal::kDouble;  
  ptkaon0L = ANtpDefVal::kDouble;
  wkaon0L = ANtpDefVal::kDouble;
  lepnkaon0L = ANtpDefVal::kInt;
  lepnkaon0L_he = ANtpDefVal::kInt;
  lepekaon0L = ANtpDefVal::kDouble;
  leppzkaon0L = ANtpDefVal::kDouble;  
  lepptkaon0L = ANtpDefVal::kDouble;
  lepwkaon0L = ANtpDefVal::kDouble;
  nkaon0S = ANtpDefVal::kInt;
  nkaon0S_he = ANtpDefVal::kInt;
  ekaon0S = ANtpDefVal::kDouble;
  pzkaon0S = ANtpDefVal::kDouble;  
  ptkaon0S = ANtpDefVal::kDouble;
  wkaon0S = ANtpDefVal::kDouble;
  lepnkaon0S = ANtpDefVal::kInt;
  lepnkaon0S_he = ANtpDefVal::kInt;
  lepekaon0S = ANtpDefVal::kDouble;
  leppzkaon0S = ANtpDefVal::kDouble;  
  lepptkaon0S = ANtpDefVal::kDouble;
  lepwkaon0S = ANtpDefVal::kDouble;
  nkaon0 = ANtpDefVal::kInt;
  nkaon0_he = ANtpDefVal::kInt;
  ekaon0 = ANtpDefVal::kDouble;
  pzkaon0 = ANtpDefVal::kDouble;  
  ptkaon0 = ANtpDefVal::kDouble;
  wkaon0 = ANtpDefVal::kDouble;
  lepnkaon0 = ANtpDefVal::kInt;
  lepnkaon0_he = ANtpDefVal::kInt;
  lepekaon0 = ANtpDefVal::kDouble;
  leppzkaon0 = ANtpDefVal::kDouble;  
  lepptkaon0 = ANtpDefVal::kDouble;
  lepwkaon0 = ANtpDefVal::kDouble;
  nutauevent = ANtpDefVal::kInt;
  nutauchannel = ANtpDefVal::kInt;

}

void StdHepInfo::Zero()
{
  nmult      = 0;
  nhad       = 0;
  nfs        = 0;  
  nlep       = 0;
  npip       = 0;
  npim       = 0;
  npi0       = 0;
  nprot      = 0;
  nneut      = 0;
  nkaon      = 0;
  ngeant     = 0;
  nfs_he     = 0;  
  nlep_he    = 0;
  npip_he    = 0;
  npim_he    = 0;
  npi0_he    = 0;
  nprot_he   = 0;
  nneut_he   = 0;
  nkaon_he   = 0;
  ngeant_he  = 0;
  etot       = 0.0;
  elep       = 0.0;
  epip       = 0.0;
  epim       = 0.0;
  epi0       = 0.0;
  eprot      = 0.0;
  eneut      = 0.0;
  ekaon      = 0.0;
  egeant     = 0.0;
  pztot      = 0.0;
  pzlep      = 0.0;
  pzpip      = 0.0;
  pzpim      = 0.0;
  pzpi0      = 0.0;
  pzprot     = 0.0;
  pzneut     = 0.0;
  pzkaon     = 0.0;
  pzgeant    = 0.0;
  pttot      = 0.0;
  ptlep      = 0.0;
  ptpip      = 0.0;
  ptpim      = 0.0;
  ptpi0      = 0.0;
  ptprot     = 0.0;
  ptneut     = 0.0;
  ptkaon     = 0.0;
  ptgeant    = 0.0;
  wfs        = 0.0;
  wlep       = 0.0;
  wpip       = 0.0;
  wpim       = 0.0;
  wpi0       = 0.0;
  wprot      = 0.0;
  wneut      = 0.0;
  wkaon      = 0.0;
  wgeant     = 0.0;
  ////// finalID    = 0.0;
  lepnmult      = 0;
  lepnfs        = 0;  
  lepnlep       = 0;
  lepnpip       = 0;
  lepnpim       = 0;
  lepnpi0       = 0;
  lepnprot      = 0;
  lepnneut      = 0;
  lepnkaon      = 0;
  lepngeant     = 0;
  lepnfs_he     = 0;  
  lepnlep_he    = 0;
  lepnpip_he    = 0;
  lepnpim_he    = 0;
  lepnpi0_he    = 0;
  lepnprot_he   = 0;
  lepnneut_he   = 0;
  lepnkaon_he   = 0;
  lepngeant_he  = 0;
  lepetot       = 0.0;
  lepelep       = 0.0;
  lepepip       = 0.0;
  lepepim       = 0.0;
  lepepi0       = 0.0;
  lepeprot      = 0.0;
  lepeneut      = 0.0;
  lepekaon      = 0.0;
  lepegeant     = 0.0;
  leppztot      = 0.0;
  leppzlep      = 0.0;
  leppzpip      = 0.0;
  leppzpim      = 0.0;
  leppzpi0      = 0.0;
  leppzprot     = 0.0;
  leppzneut     = 0.0;
  leppzkaon     = 0.0;
  leppzgeant    = 0.0;
  leppttot      = 0.0;
  lepptlep      = 0.0;
  lepptpip      = 0.0;
  lepptpim      = 0.0;
  lepptpi0      = 0.0;
  lepptprot     = 0.0;
  lepptneut     = 0.0;
  lepptkaon     = 0.0;
  lepptgeant    = 0.0;
  lepwfs        = 0.0;
  lepwlep       = 0.0;
  lepwpip       = 0.0;
  lepwpim       = 0.0;
  lepwpi0       = 0.0;
  lepwprot      = 0.0;
  lepwneut      = 0.0;
  lepwkaon      = 0.0;
  lepwgeant     = 0.0;
  lep2leptype   = 0.0;
  emenergy      = 0.0;
  emfrac        = 0.0;
  emcount       = 0;
  epi0_total    = 0.0;
  epi0_intranuke= 0.0;
  epi0_neugen   = 0.0;
  epi0_decay    = 0.0;
  epi0_abs      = 0.0;
  npi0_neugen   = 0;
  baryonpt      = 0.0;
  baryonxf      = 0.0;
  totpt      = 0.0;
  e_total = 0;
  em_total = 0;
  npi0_jbhg  = 0;
  nch_jbhg  = 0;
  ntot_jbhg  = 0;

  //added Anna vars:
  nfs_nu = 0; 
  nfs_nu_he = 0;
  etot_nu = 0.0;
  pztot_nu = 0.0;  
  pttot_nu = 0.0;
  wfs_nu = 0.0;
  lepnfs_nu = 0;
  lepnfs_nu_he = 0;
  lepetot_nu = 0.0;
  leppztot_nu = 0.0;  
  leppttot_nu = 0.0;
  lepwfs_nu = 0.0;
  n_nu_elec = 0;
  n_nu_elec_he = 0;
  e_nu_elec = 0.0;
  pz_nu_elec = 0.0;  
  pt_nu_elec = 0.0;
  w_nu_elec = 0.0;
  lepn_nu_elec = 0;
  lepn_nu_elec_he = 0;
  lepe_nu_elec = 0.0;
  leppz_nu_elec = 0.0;  
  leppt_nu_elec = 0.0;
  lepw_nu_elec = 0.0;
  n_nu_muon = 0;
  n_nu_muon_he = 0;
  e_nu_muon = 0.0;
  pz_nu_muon = 0.0;  
  pt_nu_muon = 0.0;
  w_nu_muon = 0.0;
  lepn_nu_muon = 0;
  lepn_nu_muon_he = 0;
  lepe_nu_muon = 0.0;
  leppz_nu_muon = 0.0;  
  leppt_nu_muon = 0.0;
  lepw_nu_muon = 0.0;
  n_nu_tau = 0;
  n_nu_tau_he = 0;
  e_nu_tau = 0.0;
  pz_nu_tau = 0.0;  
  pt_nu_tau = 0.0;
  w_nu_tau = 0.0;
  lepn_nu_tau = 0;
  lepn_nu_tau_he = 0;
  lepe_nu_tau = 0.0;
  leppz_nu_tau = 0.0;  
  leppt_nu_tau = 0.0;
  lepw_nu_tau = 0.0;
  nlep_elec = 0;
  nlep_elec_he = 0;
  elep_elec = 0.0;
  pzlep_elec = 0.0;  
  ptlep_elec = 0.0;
  wlep_elec = 0.0;
  lepnlep_elec = 0;
  lepnlep_elec_he = 0;
  lepelep_elec = 0.0;
  leppzlep_elec = 0.0;  
  lepptlep_elec = 0.0;
  lepwlep_elec = 0.0;
  nlep_muon = 0;
  nlep_muon_he = 0;
  elep_muon = 0.0;
  pzlep_muon = 0.0;  
  ptlep_muon = 0.0;
  wlep_muon = 0.0;
  lepnlep_muon = 0;
  lepnlep_muon_he = 0;
  lepelep_muon = 0.0;
  leppzlep_muon = 0.0;  
  lepptlep_muon = 0.0;
  lepwlep_muon = 0.0;
  nlep_tau = 0;
  nlep_tau_he = 0;
  elep_tau = 0.0;
  pzlep_tau = 0.0;  
  ptlep_tau = 0.0;
  wlep_tau = 0.0;
  lepnlep_tau = 0;
  lepnlep_tau_he = 0;
  lepelep_tau = 0.0;
  leppzlep_tau = 0.0;  
  lepptlep_tau = 0.0;
  lepwlep_tau = 0.0;
  nkaonplus = 0;
  nkaonplus_he = 0;
  ekaonplus = 0.0;
  pzkaonplus = 0.0;  
  ptkaonplus = 0.0;
  wkaonplus = 0.0;
  lepnkaonplus = 0;
  lepnkaonplus_he = 0;
  lepekaonplus = 0.0;
  leppzkaonplus = 0.0;  
  lepptkaonplus = 0.0;
  lepwkaonplus = 0.0;
  nkaonminus = 0;
  nkaonminus_he = 0;
  ekaonminus = 0.0;
  pzkaonminus = 0.0;  
  ptkaonminus = 0.0;
  wkaonminus = 0.0;
  lepnkaonminus = 0;
  lepnkaonminus_he = 0;
  lepekaonminus = 0.0;
  leppzkaonminus = 0.0;  
  lepptkaonminus = 0.0;
  lepwkaonminus = 0.0;
  nkaon0L = 0;
  nkaon0L_he = 0;
  ekaon0L = 0.0;
  pzkaon0L = 0.0;  
  ptkaon0L = 0.0;
  wkaon0L = 0.0;
  lepnkaon0L = 0;
  lepnkaon0L_he = 0;
  lepekaon0L = 0.0;
  leppzkaon0L = 0.0;  
  lepptkaon0L = 0.0;
  lepwkaon0L = 0.0;
  nkaon0S = 0;
  nkaon0S_he = 0;
  ekaon0S = 0.0;
  pzkaon0S = 0.0;  
  ptkaon0S = 0.0;
  wkaon0S = 0.0;
  lepnkaon0S = 0;
  lepnkaon0S_he = 0;
  lepekaon0S = 0.0;
  leppzkaon0S = 0.0;  
  lepptkaon0S = 0.0;
  lepwkaon0S = 0.0;
  nkaon0 = 0;
  nkaon0_he = 0;
  ekaon0 = 0.0;
  pzkaon0 = 0.0;  
  ptkaon0 = 0.0;
  wkaon0 = 0.0;
  lepnkaon0 = 0;
  lepnkaon0_he = 0;
  lepekaon0 = 0.0;
  leppzkaon0 = 0.0;  
  lepptkaon0 = 0.0;
  lepwkaon0 = 0.0;
  nutauevent = 0;
  nutauchannel = 0;



}
