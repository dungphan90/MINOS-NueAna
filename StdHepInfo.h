///////////////////////////////////////////////////////////////////////////
// 
// StdHepInfo
//
///////////////////////////////////////////////////////////////////////////
#ifndef STDHEPINFO_H
#define STDHEPINFO_H

#include "TObject.h"

class StdHepInfo : public TObject
{

 public:
  StdHepInfo();
  virtual ~StdHepInfo();

  void Reset();
  void Zero();

  //StdHepInfo variables
  //general:
  Int_t nmult;     //neugen multiplicity
  Int_t nhad;      //should be the same as nmult
  
  //numbers:
  Int_t nfs;       //number of final state particles
  Int_t nlep;      //number of leptons
  Int_t npip;      //number of pi+
  Int_t npim;      //number of pi-
  Int_t npi0;      //number of pi0
  Int_t nprot;     //number of protons
  Int_t nneut;     //number of neutrons
  Int_t nkaon;     //number of kaons (+/-/0)
  Int_t ngeant;    //number of geantinos

  Int_t nfs_he;    //number of final state particles with p>=0.5 GeV
  Int_t nlep_he;   //number of leptons with p>=0.5 GeV
  Int_t npip_he;   //number of pi+ with p>=0.5 GeV
  Int_t npim_he;   //number of pi- with p>=0.5 GeV
  Int_t npi0_he;   //number of pi0 with p>=0.5 GeV
  Int_t nprot_he;  //number of protons with p>=0.5 GeV
  Int_t nneut_he;  //number of neutrons with p>=0.5 GeV
  Int_t nkaon_he;  //number of kaons (+/-/0) with p>=0.5 GeV
  Int_t ngeant_he; //number of geantinos with p>=0.5 GeV

  //total energy:
  Double_t etot;      //total energy of final state particles
  Double_t elep;      //total energy of leptons
  Double_t epip;      //total energy of pi+
  Double_t epim;      //total energy of pi-
  Double_t epi0;      //total energy of pi0
  Double_t eprot;     //total energy of protons
  Double_t eneut;     //total energy of neutrons
  Double_t ekaon;     //total energy of kaons (+/-/0)
  Double_t egeant;    //total energy of geantinos

  //momenta:
  Double_t pztot;      //longitudinal momentum of final state particles
  Double_t pzlep;      //longitudinal momentum of leptons
  Double_t pzpip;      //longitudinal momentum of pi+
  Double_t pzpim;      //longitudinal momentum of pi-
  Double_t pzpi0;      //longitudinal momentum of pi0
  Double_t pzprot;     //longitudinal momentum of protons
  Double_t pzneut;     //longitudinal momentum of neutrons
  Double_t pzkaon;     //longitudinal momentum of kaons (+/-/0)
  Double_t pzgeant;    //longitudinal momentum of geantinos

  Double_t pttot;      //transverse momentum of final state particles
  Double_t ptlep;      //transverse momentum of leptons
  Double_t ptpip;      //transverse momentum of pi+
  Double_t ptpim;      //transverse momentum of pi-
  Double_t ptpi0;      //transverse momentum of pi0
  Double_t ptprot;     //transverse momentum of protons
  Double_t ptneut;     //transverse momentum of neutrons
  Double_t ptkaon;     //transverse momentum of kaons (+/-/0)
  Double_t ptgeant;    //transverse momentum of geantinos

  //effective number:
  Double_t wfs;       //effective number of final states
  Double_t wlep;      //effective number of leptons
  Double_t wpip;      //effective number of pi+
  Double_t wpim;      //effective number of pi-
  Double_t wpi0;      //effective number of pi0
  Double_t wprot;     //effective number of protons
  Double_t wneut;     //effective number of neutrons
  Double_t wkaon;     //effective number of kaons (+/-/0)
  Double_t wgeant;    //effective number of geantinos

  //added by steve cavanaugh

  ///  Double_t finalID;   //final ID of particle which is daughter of inital neutrino

  //StdHepInfo variables with leptonic parent... not all variables are in use at this time!
  //general:
  Int_t lepnmult;     //neugen multiplicity
  
  //numbers:
  Int_t lepnfs;       //number of final state particles
  Int_t lepnlep;      //number of leptons
  Int_t lepnpip;      //number of pi+
  Int_t lepnpim;      //number of pi-
  Int_t lepnpi0;      //number of pi0
  Int_t lepnprot;     //number of protons
  Int_t lepnneut;     //number of neutrons
  Int_t lepnkaon;     //number of kaons (+/-/0)
  Int_t lepngeant;    //number of geantinos

  Int_t lepnfs_he;    //number of final state particles with p>=0.5 GeV
  Int_t lepnlep_he;   //number of leptons with p>=0.5 GeV
  Int_t lepnpip_he;   //number of pi+ with p>=0.5 GeV
  Int_t lepnpim_he;   //number of pi- with p>=0.5 GeV
  Int_t lepnpi0_he;   //number of pi0 with p>=0.5 GeV
  Int_t lepnprot_he;  //number of protons with p>=0.5 GeV
  Int_t lepnneut_he;  //number of neutrons with p>=0.5 GeV
  Int_t lepnkaon_he;  //number of kaons (+/-/0) with p>=0.5 GeV
  Int_t lepngeant_he; //number of geantinos with p>=0.5 GeV

  //total energy:
  Double_t lepetot;      //total energy of final state particles
  Double_t lepelep;      //total energy of leptons
  Double_t lepepip;      //total energy of pi+
  Double_t lepepim;      //total energy of pi-
  Double_t lepepi0;      //total energy of pi0
  Double_t lepeprot;     //total energy of protons
  Double_t lepeneut;     //total energy of neutrons
  Double_t lepekaon;     //total energy of kaons (+/-/0)
  Double_t lepegeant;    //total energy of geantinos

  //momenta:
  Double_t leppztot;      //longitudinal momentum of final state particles
  Double_t leppzlep;      //longitudinal momentum of leptons
  Double_t leppzpip;      //longitudinal momentum of pi+
  Double_t leppzpim;      //longitudinal momentum of pi-
  Double_t leppzpi0;      //longitudinal momentum of pi0
  Double_t leppzprot;     //longitudinal momentum of protons
  Double_t leppzneut;     //longitudinal momentum of neutrons
  Double_t leppzkaon;     //longitudinal momentum of kaons (+/-/0)
  Double_t leppzgeant;    //longitudinal momentum of geantinos

  Double_t leppttot;      //transverse momentum of final state particles
  Double_t lepptlep;      //transverse momentum of leptons
  Double_t lepptpip;      //transverse momentum of pi+
  Double_t lepptpim;      //transverse momentum of pi-
  Double_t lepptpi0;      //transverse momentum of pi0
  Double_t lepptprot;     //transverse momentum of protons
  Double_t lepptneut;     //transverse momentum of neutrons
  Double_t lepptkaon;     //transverse momentum of kaons (+/-/0)
  Double_t lepptgeant;    //transverse momentum of geantinos

  //effective number:
  Double_t lepwfs;       //effective number of final states
  Double_t lepwlep;      //effective number of leptons
  Double_t lepwpip;      //effective number of pi+
  Double_t lepwpim;      //effective number of pi-
  Double_t lepwpi0;      //effective number of pi0
  Double_t lepwprot;     //effective number of protons
  Double_t lepwneut;     //effective number of neutrons
  Double_t lepwkaon;     //effective number of kaons (+/-/0)
  Double_t lepwgeant;    //effective number of geantinos


  Double_t lep2leptype;  //particle id of finalstate lepton e, u, t (which then decays) if applicable

  Double_t emfrac;       //fraction of em energy to total energy
  Int_t emcount;         //number of em particles in em energy
  Double_t emenergy;     //total em energy



    /////////end sc edit

  //calculate pi0 energy  T.Y. 09/06
  Double_t epi0_total;       //total pi0's in the final state
  Double_t epi0_intranuke;   //pi0's generated in the intranuclear interaction
  Double_t epi0_neugen;      //pi0's generated in neugen hadronization model
  Double_t epi0_decay;       //pi0's from hadron decay(Ks...)
  Double_t epi0_abs;         //pi0's absorbed in the intranuclear interaction
  Int_t    npi0_neugen;             //number of pi0's generated in neugen hadronization model
  Double_t baryonpt;         //pt of the baryon in the hadronic system
  Double_t baryonxf;         //xf of the baryon in the hadronic system
  Double_t totpt;            //overall pt


  Double_t e_total;    // Defined by Hugh/Josh for reweighting
  Double_t em_total;   // Defined by Hugh/Josh for reweighting
  Int_t    npi0_jbhg;  // Defined by Hugh/Josh for reweighting
  Int_t    nch_jbhg;   // Defined by Hugh/Josh for reweighting
  Int_t    ntot_jbhg;  // Defined by Hugh/Josh for reweighting


  //Now include Anna Vars for nutau analysis
  Int_t nfs_nu;       //number of neutrinos in final states
  Int_t nfs_nu_he;  //numbers of neutrinos above 0.5GeV
  Double_t etot_nu;  //total energy taken away by neutrinos in final states
  Double_t pztot_nu;  //total final state neutrino pz
  Double_t pttot_nu;  //total final state neutrino pt
  Double_t wfs_nu;   //effective final state neutrinos 
  Int_t lepnfs_nu;  //same as above, but counts only neutrinos coming from original neutrino parent
   Int_t lepnfs_nu_he;
   Double_t lepetot_nu;
   Double_t leppztot_nu;
   Double_t leppttot_nu;
   Double_t lepwfs_nu;

   Int_t n_nu_elec; //same as first block, but only counts electron neutrinos
   Int_t n_nu_elec_he;
   Double_t e_nu_elec;
   Double_t pz_nu_elec;
   Double_t pt_nu_elec;
   Double_t w_nu_elec;
   Int_t lepn_nu_elec;
   Int_t lepn_nu_elec_he;
   Double_t lepe_nu_elec;
   Double_t leppz_nu_elec;
   Double_t leppt_nu_elec;
   Double_t lepw_nu_elec;

   Int_t n_nu_muon; //same as first block, but only counts muon neutrinos
   Int_t n_nu_muon_he;
   Double_t e_nu_muon;
   Double_t pz_nu_muon;
   Double_t pt_nu_muon;
   Double_t w_nu_muon;
   Int_t lepn_nu_muon;
   Int_t lepn_nu_muon_he;
   Double_t lepe_nu_muon;
   Double_t leppz_nu_muon;
   Double_t leppt_nu_muon;
   Double_t lepw_nu_muon;

   Int_t n_nu_tau; //same as first block, but only counts tau neutrinos
   Int_t n_nu_tau_he;
   Double_t e_nu_tau;
   Double_t pz_nu_tau;
   Double_t pt_nu_tau;
   Double_t w_nu_tau;
   Int_t lepn_nu_tau;
   Int_t lepn_nu_tau_he;
   Double_t lepe_nu_tau;
   Double_t leppz_nu_tau;
   Double_t leppt_nu_tau;
   Double_t lepw_nu_tau;

  //now this is the block that counts all final states EXCEPT final state neutrinos
   Int_t nlep_elec; // this block counts electrons
   Int_t nlep_elec_he;
   Double_t elep_elec;
   Double_t pzlep_elec;
   Double_t ptlep_elec;
   Double_t wlep_elec;
   Int_t lepnlep_elec; //analogously to before, this counts electrons coming from original neutrino parent
   Int_t lepnlep_elec_he;
   Double_t lepelep_elec;
   Double_t leppzlep_elec;
   Double_t lepptlep_elec;
   Double_t lepwlep_elec;

   Int_t nlep_muon;  // this block counts muons
   Int_t nlep_muon_he;
   Double_t elep_muon;
   Double_t pzlep_muon;
   Double_t ptlep_muon;
   Double_t wlep_muon;
   Int_t lepnlep_muon;
   Int_t lepnlep_muon_he;
   Double_t lepelep_muon;
   Double_t leppzlep_muon;
   Double_t lepptlep_muon;
   Double_t lepwlep_muon;

   Int_t nlep_tau;   // this block counts taus in final states - this should always be 0 probably
   Int_t nlep_tau_he;
   Double_t elep_tau;
   Double_t pzlep_tau;
   Double_t ptlep_tau;
   Double_t wlep_tau;
   Int_t lepnlep_tau;
   Int_t lepnlep_tau_he;
   Double_t lepelep_tau;
   Double_t leppzlep_tau;
   Double_t lepptlep_tau;
   Double_t lepwlep_tau;

   Int_t nkaonplus; // this block counts K+
   Int_t nkaonplus_he;
   Double_t ekaonplus;
   Double_t pzkaonplus;
   Double_t ptkaonplus;
   Double_t wkaonplus;
   Int_t lepnkaonplus;
   Int_t lepnkaonplus_he;
   Double_t lepekaonplus;
   Double_t leppzkaonplus;
   Double_t lepptkaonplus;
   Double_t lepwkaonplus;

   Int_t nkaonminus;    // this block counts K-
   Int_t nkaonminus_he;
   Double_t ekaonminus;
   Double_t pzkaonminus;
   Double_t ptkaonminus;
   Double_t wkaonminus;
   Int_t lepnkaonminus;
   Int_t lepnkaonminus_he;
   Double_t lepekaonminus;
   Double_t leppzkaonminus;
   Double_t lepptkaonminus;
   Double_t lepwkaonminus;

   Int_t nkaon0L;    // this block counts K0L
   Int_t nkaon0L_he;
   Double_t ekaon0L;
   Double_t pzkaon0L;
   Double_t ptkaon0L;
   Double_t wkaon0L;
   Int_t lepnkaon0L;
   Int_t lepnkaon0L_he;
   Double_t lepekaon0L;
   Double_t leppzkaon0L;
   Double_t lepptkaon0L;
   Double_t lepwkaon0L;

   Int_t nkaon0S;   // this block counts K0S
   Int_t nkaon0S_he;
   Double_t ekaon0S;
   Double_t pzkaon0S;
   Double_t ptkaon0S;
   Double_t wkaon0S;
   Int_t lepnkaon0S;
   Int_t lepnkaon0S_he;
   Double_t lepekaon0S;
   Double_t leppzkaon0S;
   Double_t lepptkaon0S;
   Double_t lepwkaon0S;

   Int_t nkaon0;   // this block counts K0
   Int_t nkaon0_he;
   Double_t ekaon0;
   Double_t pzkaon0;
   Double_t ptkaon0;
   Double_t wkaon0;
   Int_t lepnkaon0;
   Int_t lepnkaon0_he;
   Double_t lepekaon0;
   Double_t leppzkaon0;
   Double_t lepptkaon0;
   Double_t lepwkaon0;

   Int_t nutauevent;   // This variable just says whether it is a nutau event: 0- not nutau, 1 = nutau CC, -1 = antinutau CC, 2 = nutau NC
   Int_t nutauchannel;  // gives the decay channel: 1=nutau->tau->muon, 2=nutau->tau->electron, 3=nutau->hadronic with pi0s; 4=nutau->hadronic without pi0s; 5=other;

  ClassDef(StdHepInfo,10)
};

#endif// STDHEPINFO_H
