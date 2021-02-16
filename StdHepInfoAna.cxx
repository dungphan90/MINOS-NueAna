#include "CandNtupleSR/NtpSREvent.h"
#include "StandardNtuple/NtpStRecord.h"
#include "MCNtuple/NtpMCTruth.h"
#include "MCNtuple/NtpMCStdHep.h"
#include "TruthHelperNtuple/NtpTHEvent.h"
#include "MessageService/MsgService.h"
#include "StdHepInfoAna.h"
#include "AnalysisNtuples/Module/ANtpRecoNtpManipulator.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

CVSID("$Id: StdHepInfoAna.cxx,v 1.17 2013/05/15 12:45:08 annah1 Exp $");

void CalcDaughters(int start, int end, double &E, double &EmE, int &npi0, int &nch, int &ntot, TClonesArray& stdhepArray);
void HughCode(NtpMCTruth* mctruth, TClonesArray& stdhepArray, StdHepInfo &shi);


StdHepInfoAna::StdHepInfoAna(StdHepInfo &shi):
  fStdHepInfo(shi)
{}

StdHepInfoAna::~StdHepInfoAna()
{}

bool isLeptonic(int parent, TClonesArray& stdhepArray);
bool checkParent(int parent, TClonesArray& stdhepArray);
void countpi0(int index, TClonesArray& stdhepArray, StdHepInfo &fStdHepInfo);
void calbaryonpt(int index, TClonesArray& stdhepArray, StdHepInfo &fStdHepInfo);

void StdHepInfoAna::Analyze(int event, RecRecordImp<RecCandHeader> *srobj)
{
  if(srobj==0){
    return;
  }
  NtpStRecord *st = dynamic_cast<NtpStRecord *>(srobj);
  Analyze(event,st);
}

void StdHepInfoAna::Analyze(int evtn, NtpStRecord *srobj)
{
  fStdHepInfo.Reset();
  if(srobj==0){
    return;
  }
  NtpMCTruth *mctruth = 0;
  NtpMCStdHep *ntpStdHep = 0;
  NtpTHEvent *thevent = 0;
  ANtpRecoNtpManipulator ntpManipulator(srobj);
  thevent = ntpManipulator.GetMCManipulator()->GetNtpTHEvent(evtn);
  if(thevent){
    mctruth = ntpManipulator.GetMCManipulator()->GetNtpMCTruth(thevent->neumc);
    if(mctruth){
      fStdHepInfo.Zero();
      TClonesArray& stdhepArray = *(srobj->stdhep);
      Int_t nStdHep = stdhepArray.GetEntries();
      HughCode(mctruth, stdhepArray, fStdHepInfo);
      Double_t maxmom = 0;
      Double_t maxke = 0;
      Double_t lepmaxmom = 0;
      Double_t lepmaxke = 0;
      Double_t maxmom_nu = 0;  //edit by Anna
      Double_t maxke_nu = 0;  //edit by Anna
      Double_t lepmaxmom_nu = 0;  //edit by Anna
      Double_t lepmaxke_nu = 0;  //edit by Anna
      Double_t initialNuEnergy=0;
      //count pi0's
      countpi0(mctruth->index,stdhepArray,fStdHepInfo);
      //calculate baryon pt
      calbaryonpt(mctruth->index,stdhepArray,fStdHepInfo);
      if (mctruth->w2>0) fStdHepInfo.baryonxf *= 2./sqrt(mctruth->w2);

      for(int i=mctruth->stdhep[0];i<=mctruth->stdhep[1]&&i<nStdHep;i++){
	ntpStdHep = dynamic_cast<NtpMCStdHep *>(stdhepArray[i]);

	//get initial neutrino energy:
	if(ntpStdHep->mc==mctruth->index &&ntpStdHep->IstHEP==0&& (abs(ntpStdHep->IdHEP)==12 || abs(ntpStdHep->IdHEP)==14 || abs(ntpStdHep->IdHEP)==16)){
	  initialNuEnergy=TMath::Abs(ntpStdHep->p4[3]);

	  	  if(ntpStdHep->IdHEP==16&&mctruth->iaction==1) fStdHepInfo.nutauevent=1;  //nutau CC event
	          if(ntpStdHep->IdHEP==-16&&mctruth->iaction==1) fStdHepInfo.nutauevent=-1; // anti-nutau CC event
	  	  if(abs(ntpStdHep->IdHEP)==16&&mctruth->iaction==0) fStdHepInfo.nutauevent=2; //nutau NC event
		  // fStdHepInfo.nutauevent=2; //nutau NC event
	}


	//check that stdhep entry corresponds to this mc event
	//count final state neutrinos only here - edit by Anna
	if(ntpStdHep->mc==mctruth->index && 
          (abs(ntpStdHep->IdHEP)==12 || abs(ntpStdHep->IdHEP)==14 || abs(ntpStdHep->IdHEP)==16) &&
           ntpStdHep->IstHEP==1){       //get final state neutrinos:  

	    Double_t energy = TMath::Abs(ntpStdHep->p4[3]);
	    Double_t pz = TMath::Abs(ntpStdHep->p4[2]);
	    Double_t pt = ( ntpStdHep->p4[0]*ntpStdHep->p4[0] +
			    ntpStdHep->p4[1]*ntpStdHep->p4[1] );
	    Double_t mom = TMath::Sqrt(pz*pz + pt);
	    Double_t ke = energy - ntpStdHep->mass;

	    if(mom>maxmom_nu) maxmom_nu = mom;
	    if(ke>maxke_nu) maxke_nu = ke;
	    if(ke>0.5) fStdHepInfo.nfs_nu_he += 1;
	    
	    fStdHepInfo.nfs_nu += 1;
	    fStdHepInfo.etot_nu += energy;
	    fStdHepInfo.pztot_nu += pz;
	    fStdHepInfo.pttot_nu += pt;
	    fStdHepInfo.wfs_nu += ke;

	    if (isLeptonic(ntpStdHep->index, stdhepArray)){
	      fStdHepInfo.lepnfs_nu += 1;
	      fStdHepInfo.lepetot_nu += energy;
	      fStdHepInfo.leppztot_nu += pz;
	      fStdHepInfo.leppttot_nu += pt;
	      fStdHepInfo.lepwfs_nu += ke;
	      
	      if(mom>lepmaxmom_nu) lepmaxmom_nu = mom;
	      if(ke>lepmaxke_nu) lepmaxke_nu = ke;
	      if(ke>0.5) fStdHepInfo.lepnfs_nu_he += 1;

	      if(ntpStdHep->IdHEP==12) {
		fStdHepInfo.lepn_nu_elec += 1;
		if(ke>0.5) fStdHepInfo.lepn_nu_elec_he += 1;
		fStdHepInfo.lepe_nu_elec += energy;
		fStdHepInfo.leppz_nu_elec += pz;
		fStdHepInfo.leppt_nu_elec += pt;
		fStdHepInfo.lepw_nu_elec += ke;
	      }
	      else if(ntpStdHep->IdHEP==14) {
		fStdHepInfo.lepn_nu_muon += 1;
		if(ke>0.5) fStdHepInfo.lepn_nu_muon_he += 1;
		fStdHepInfo.lepe_nu_muon += energy;
		fStdHepInfo.leppz_nu_muon += pz;
		fStdHepInfo.leppt_nu_muon += pt;
		fStdHepInfo.lepw_nu_muon += ke;
	      }
	      else if(ntpStdHep->IdHEP==16) {
		fStdHepInfo.lepn_nu_tau += 1;
		if(ke>0.5) fStdHepInfo.lepn_nu_tau_he += 1;
		fStdHepInfo.lepe_nu_tau += energy;
		fStdHepInfo.leppz_nu_tau += pz;
		fStdHepInfo.leppt_nu_tau += pt;
		fStdHepInfo.lepw_nu_tau += ke;
	      }
	    }  //end of isleptonic for final state neutrinos

	    
	      if(ntpStdHep->IdHEP==12) {
		fStdHepInfo.n_nu_elec += 1;
		if(ke>0.5) fStdHepInfo.n_nu_elec_he += 1;
		fStdHepInfo.e_nu_elec += energy;
		fStdHepInfo.pz_nu_elec += pz;
		fStdHepInfo.pt_nu_elec += pt;
		fStdHepInfo.w_nu_elec += ke;
	      }
	      else if(ntpStdHep->IdHEP==14) {
		fStdHepInfo.n_nu_muon += 1;
		if(ke>0.5) fStdHepInfo.n_nu_muon_he += 1;
		fStdHepInfo.e_nu_muon += energy;
		fStdHepInfo.pz_nu_muon += pz;
		fStdHepInfo.pt_nu_muon += pt;
		fStdHepInfo.w_nu_muon += ke;
	      }
	      else if(ntpStdHep->IdHEP==16) {
		fStdHepInfo.n_nu_tau += 1;
		if(ke>0.5) fStdHepInfo.n_nu_tau_he += 1;
		fStdHepInfo.e_nu_tau += energy;
		fStdHepInfo.pz_nu_tau += pz;
		fStdHepInfo.pt_nu_tau += pt;
		fStdHepInfo.w_nu_tau += ke;
	      }
	
	} //edit Anna edit for final state neutrinos
      
      
	//check that stdhep entry corresponds to this mc event
	//also don't include final state neutrinos in counting
	if(ntpStdHep->mc==mctruth->index && abs(ntpStdHep->IdHEP)!=12 &&
	   abs(ntpStdHep->IdHEP)!=14 && abs(ntpStdHep->IdHEP)!=16){
	  //look for neugen intermediate particle to get multiplicity
	  if(ntpStdHep->IstHEP==3 &&
	     (mctruth->iresonance!=1002 || 
	      TMath::Abs(ntpStdHep->IdHEP)!=15)){
	    Int_t hfs = ntpStdHep->IdHEP%1000;
	    fStdHepInfo.nmult = Int_t((hfs-128)/4)-20;
	  }
	

	  // look for pi0, e, and photons to make the emcount, emfrac, emenergy variables
	  //check final decays, but dont double count energies (dont count a particle if its parents energy is also to be counted
	  //i am using a recursive method of checking - i do not know how far the mc decay tree extends in IstHEP==205 
	  //it goes at least two decays deep in IstHEP==205, possible more
	  else if(ntpStdHep->IstHEP==205){
	    
	    if(abs(ntpStdHep->IdHEP)==11||abs(ntpStdHep->IdHEP)==22||abs(ntpStdHep->IdHEP)==111){
	      Double_t energy = TMath::Abs(ntpStdHep->p4[3]);
	      if(!checkParent(ntpStdHep->parent[0],stdhepArray)){
		if(ntpStdHep->parent[0]!=ntpStdHep->parent[1]){
		  if(!checkParent(ntpStdHep->parent[1],stdhepArray)){
		    fStdHepInfo.emenergy+=energy;
		    fStdHepInfo.emcount++;
		  }
		}else{
		  fStdHepInfo.emenergy+=energy;
		  fStdHepInfo.emcount++;
		}
	      }		
	    }
	  }





	  else if((ntpStdHep->IstHEP==1 ) && //get final state particles:  
		  ntpStdHep->IdHEP<1000000000){ //excluding nuclei
	    Double_t energy = TMath::Abs(ntpStdHep->p4[3]);
	    Double_t pz = TMath::Abs(ntpStdHep->p4[2]);
	    Double_t pt = ( ntpStdHep->p4[0]*ntpStdHep->p4[0] +
			    ntpStdHep->p4[1]*ntpStdHep->p4[1] );
	    Double_t mom = TMath::Sqrt(pz*pz + pt);
	    Double_t ke = energy - ntpStdHep->mass;

	    if(mom>maxmom) maxmom = mom;
	    if(ke>maxke) maxke = ke;
	    if(ke>0.5) fStdHepInfo.nfs_he += 1;
	    
	    fStdHepInfo.nfs += 1;
	    fStdHepInfo.etot += energy;
	    fStdHepInfo.pztot += pz;
	    fStdHepInfo.pttot += pt;
	    fStdHepInfo.wfs += ke;

	    
            // by steve cavanaugh
	   
	    //electrons, protons, and pi0 contribute to em frac
	    //decay particle of these which also fit 11,22,111 do not contribute as per code dealing with IstHEP===205 above
	    if(abs(ntpStdHep->IdHEP)==11||abs(ntpStdHep->IdHEP)==22||abs(ntpStdHep->IdHEP)==111){
	      fStdHepInfo.emenergy+=energy;
	      fStdHepInfo.emcount++;
	    }




	    //see if the current particle is a decay product of the original neutrino
	    if (isLeptonic(ntpStdHep->index, stdhepArray)){
	    
	      fStdHepInfo.lepnfs += 1;
	      fStdHepInfo.lepetot += energy;
	      fStdHepInfo.leppztot += pz;
	      fStdHepInfo.leppttot += pt;
	      fStdHepInfo.lepwfs += ke;
	      
	      if(mom>lepmaxmom) lepmaxmom = mom;
	      if(ke>lepmaxke) lepmaxke = ke;
	      if(ke>0.5) fStdHepInfo.lepnfs_he += 1;
	     

	      //if not a tau, check if electron or muon, if so set flags and energies
	      if(abs(ntpStdHep->IdHEP)==11||abs(ntpStdHep->IdHEP)==13||abs(ntpStdHep->IdHEP)==15){ // should i be checking tau also? probably only appears with IstHEP=2
		fStdHepInfo.lep2leptype =ntpStdHep->IdHEP;
	      }
	    
			 

	      //otherwise do standard particle contribution ...

	      if(abs(ntpStdHep->IdHEP)==11 || abs(ntpStdHep->IdHEP)==13 || 
		 abs(ntpStdHep->IdHEP)==15) {     
//ive probably already accounted for tau energyies from the daughter particles.... so i probably do not want them here!!
		fStdHepInfo.lepnlep += 1;
		if(ke>0.5) fStdHepInfo.lepnlep_he += 1;
		fStdHepInfo.lepelep += energy;
		fStdHepInfo.leppzlep += pz;
		fStdHepInfo.lepptlep += pt;
		fStdHepInfo.lepwlep += ke;
		//edit by Anna to find out which leptons are present                                                                                                      
		if(abs(ntpStdHep->IdHEP)==11){
		  fStdHepInfo.lepnlep_elec += 1;
		  if(ke>0.5) fStdHepInfo.lepnlep_elec_he += 1;
		  fStdHepInfo.lepelep_elec += energy;
		  fStdHepInfo.leppzlep_elec += pz;
		  fStdHepInfo.lepptlep_elec += pt;
		  fStdHepInfo.lepwlep_elec += ke;
		}
		if(abs(ntpStdHep->IdHEP)==13){
		  fStdHepInfo.lepnlep_muon += 1;
		  if(ke>0.5) fStdHepInfo.lepnlep_muon_he += 1;
		  fStdHepInfo.lepelep_muon += energy;
		  fStdHepInfo.leppzlep_muon += pz;
		  fStdHepInfo.lepptlep_muon += pt;
		  fStdHepInfo.lepwlep_muon += ke;
		}
		if(abs(ntpStdHep->IdHEP)==15){
		  fStdHepInfo.lepnlep_tau += 1;
		  if(ke>0.5) fStdHepInfo.lepnlep_tau_he += 1;
		  fStdHepInfo.lepelep_tau += energy;
		  fStdHepInfo.leppzlep_tau += pz;
		  fStdHepInfo.lepptlep_tau += pt;
		  fStdHepInfo.lepwlep_tau += ke;
		}
		//end of edit by Anna to find out which leptons are present   
	      }
	      else if(ntpStdHep->IdHEP==28) {
		fStdHepInfo.lepngeant += 1;
		if(ke>0.5) fStdHepInfo.lepngeant_he += 1;
		fStdHepInfo.lepegeant += energy;
		fStdHepInfo.leppzgeant += pz;
		fStdHepInfo.lepptgeant += pt;
		fStdHepInfo.lepwgeant += ke;
	      }
	      else if(ntpStdHep->IdHEP==2212) {
		fStdHepInfo.lepnprot += 1;
		if(ke>0.5) fStdHepInfo.lepnprot_he += 1;
		fStdHepInfo.lepeprot += energy;
		fStdHepInfo.leppzprot += pz;
		fStdHepInfo.lepptprot += pt;
		fStdHepInfo.lepwprot += ke;
	      }
	      else if(ntpStdHep->IdHEP==2112) {
		fStdHepInfo.lepnneut += 1;
		if(ke>0.5) fStdHepInfo.lepnneut_he += 1;
		fStdHepInfo.lepeneut += energy;
		fStdHepInfo.leppzneut += pz;
		fStdHepInfo.lepptneut += pt;
		fStdHepInfo.lepwneut += ke;
	      }
	      else if(ntpStdHep->IdHEP==211) {
		fStdHepInfo.lepnpip += 1;
		if(ke>0.5) fStdHepInfo.lepnpip_he += 1;
		fStdHepInfo.lepepip += energy;
		fStdHepInfo.leppzpip += pz;
		fStdHepInfo.lepptpip += pt;
		fStdHepInfo.lepwpip += ke;
	      }
	      else if(ntpStdHep->IdHEP==-211) {
		fStdHepInfo.lepnpim += 1;
		if(ke>0.5) fStdHepInfo.lepnpim_he += 1;
		fStdHepInfo.lepepim += energy;
		fStdHepInfo.leppzpim += pz;
		fStdHepInfo.lepptpim += pt;
		fStdHepInfo.lepwpim += ke;
	      }
	      else if(ntpStdHep->IdHEP==111) {
		fStdHepInfo.lepnpi0 += 1;
		if(ke>0.5) fStdHepInfo.lepnpi0_he += 1;
		fStdHepInfo.lepepi0 += energy;
		fStdHepInfo.leppzpi0 += pz;
		fStdHepInfo.lepptpi0 += pt;
		fStdHepInfo.lepwpi0 += ke;
	      }
	      else if(TMath::Abs(ntpStdHep->IdHEP)==321 || //K+/-
		      ntpStdHep->IdHEP==130 ||     //K0_L
		      ntpStdHep->IdHEP==310 ||     //K0_S
		      ntpStdHep->IdHEP==311 ) {    //K0
		fStdHepInfo.lepnkaon += 1;
		if(ke>0.5) fStdHepInfo.lepnkaon_he += 1;
		fStdHepInfo.lepekaon += energy;
		fStdHepInfo.leppzkaon += pz;
		fStdHepInfo.lepptkaon += pt;
		fStdHepInfo.lepwkaon += ke;

		//edit by Anna within the kaon loop to get the kaon species details                                                                                     
                if(ntpStdHep->IdHEP==321){      //K+                                                                                                                    
		  fStdHepInfo.lepnkaonplus += 1;
		  if(ke>0.5) fStdHepInfo.lepnkaonplus_he += 1;
		  fStdHepInfo.lepekaonplus += energy;
		  fStdHepInfo.leppzkaonplus += pz;
		  fStdHepInfo.lepptkaonplus += pt;
		  fStdHepInfo.lepwkaonplus += ke;
		}
                if(ntpStdHep->IdHEP==-321){  //K-                                                                                                                       
		  fStdHepInfo.lepnkaonminus += 1;
		  if(ke>0.5) fStdHepInfo.lepnkaonminus_he += 1;
		  fStdHepInfo.lepekaonminus += energy;
		  fStdHepInfo.leppzkaonminus += pz;
		  fStdHepInfo.lepptkaonminus += pt;
		  fStdHepInfo.lepwkaonminus += ke;
		}
                if(ntpStdHep->IdHEP==130){     //K0_L                                                                                                                   
		  fStdHepInfo.lepnkaon0L += 1;
		  if(ke>0.5) fStdHepInfo.lepnkaon0L_he += 1;
		  fStdHepInfo.lepekaon0L += energy;
		  fStdHepInfo.leppzkaon0L += pz;
		  fStdHepInfo.lepptkaon0L += pt;
		  fStdHepInfo.lepwkaon0L += ke;
		}
                if(ntpStdHep->IdHEP==310){    //K0_S                                                                                                                    
		  fStdHepInfo.lepnkaon0S += 1;
		  if(ke>0.5) fStdHepInfo.lepnkaon0S_he += 1;
		  fStdHepInfo.lepekaon0S += energy;
		  fStdHepInfo.leppzkaon0S += pz;
		  fStdHepInfo.lepptkaon0S += pt;
		  fStdHepInfo.lepwkaon0S += ke;
		}
                if(ntpStdHep->IdHEP==311){   //K0                                                                                                                       
		  fStdHepInfo.lepnkaon0 += 1;
		  if(ke>0.5) fStdHepInfo.lepnkaon0_he += 1;
		  fStdHepInfo.lepekaon0 += energy;
		  fStdHepInfo.leppzkaon0 += pz;
		  fStdHepInfo.lepptkaon0 += pt;
		  fStdHepInfo.lepwkaon0 += ke;
		}
		//end of edit by Anna to count individual kaon species                 

	      }


	    }
            //end edit by steve cavanaugh
	    

	    if(abs(ntpStdHep->IdHEP)==11 || abs(ntpStdHep->IdHEP)==13 || 
	       abs(ntpStdHep->IdHEP)==15) {
	      fStdHepInfo.nlep += 1;
	      if(ke>0.5) fStdHepInfo.nlep_he += 1;
	      fStdHepInfo.elep += energy;
	      fStdHepInfo.pzlep += pz;
	      fStdHepInfo.ptlep += pt;
	      fStdHepInfo.wlep += ke;
              //edit by Anna to find out which leptons are present                                                                                                      
              if(abs(ntpStdHep->IdHEP)==11){
		fStdHepInfo.nlep_elec += 1;
		if(ke>0.5) fStdHepInfo.nlep_elec_he += 1;
		fStdHepInfo.elep_elec += energy;
		fStdHepInfo.pzlep_elec += pz;
		fStdHepInfo.ptlep_elec += pt;
		fStdHepInfo.wlep_elec += ke;
              }
              if(abs(ntpStdHep->IdHEP)==13){
		fStdHepInfo.nlep_muon += 1;
		if(ke>0.5) fStdHepInfo.nlep_muon_he += 1;
		fStdHepInfo.elep_muon += energy;
		fStdHepInfo.pzlep_muon += pz;
		fStdHepInfo.ptlep_muon += pt;
		fStdHepInfo.wlep_muon += ke;
              }
              if(abs(ntpStdHep->IdHEP)==15){
		fStdHepInfo.nlep_tau += 1;
		if(ke>0.5) fStdHepInfo.nlep_tau_he += 1;
		fStdHepInfo.elep_tau += energy;
		fStdHepInfo.pzlep_tau += pz;
		fStdHepInfo.ptlep_tau += pt;
		fStdHepInfo.wlep_tau += ke;
              }
              //end of edit by Anna to find out which leptons are present                                                                                               
	    }
	    else if(ntpStdHep->IdHEP==28) {
	      fStdHepInfo.ngeant += 1;
	      if(ke>0.5) fStdHepInfo.ngeant_he += 1;
	      fStdHepInfo.egeant += energy;
	      fStdHepInfo.pzgeant += pz;
	      fStdHepInfo.ptgeant += pt;
	      fStdHepInfo.wgeant += ke;
	    }
	    else if(ntpStdHep->IdHEP==2212) {
	      fStdHepInfo.nprot += 1;
	      if(ke>0.5) fStdHepInfo.nprot_he += 1;
	      fStdHepInfo.eprot += energy;
	      fStdHepInfo.pzprot += pz;
	      fStdHepInfo.ptprot += pt;
	      fStdHepInfo.wprot += ke;
	    }
	    else if(ntpStdHep->IdHEP==2112) {
	      fStdHepInfo.nneut += 1;
	      if(ke>0.5) fStdHepInfo.nneut_he += 1;
	      fStdHepInfo.eneut += energy;
	      fStdHepInfo.pzneut += pz;
	      fStdHepInfo.ptneut += pt;
	      fStdHepInfo.wneut += ke;
	    }
	    else if(ntpStdHep->IdHEP==211) {
	      fStdHepInfo.npip += 1;
	      if(ke>0.5) fStdHepInfo.npip_he += 1;
	      fStdHepInfo.epip += energy;
	      fStdHepInfo.pzpip += pz;
	      fStdHepInfo.ptpip += pt;
	      fStdHepInfo.wpip += ke;
	    }
	    else if(ntpStdHep->IdHEP==-211) {
	      fStdHepInfo.npim += 1;
	      if(ke>0.5) fStdHepInfo.npim_he += 1;
	      fStdHepInfo.epim += energy;
	      fStdHepInfo.pzpim += pz;
	      fStdHepInfo.ptpim += pt;
	      fStdHepInfo.wpim += ke;
	    }
	    else if(ntpStdHep->IdHEP==111) {
	      fStdHepInfo.npi0 += 1;
	      if(ke>0.5) fStdHepInfo.npi0_he += 1;
	      fStdHepInfo.epi0 += energy;
	      fStdHepInfo.pzpi0 += pz;
	      fStdHepInfo.ptpi0 += pt;
	      fStdHepInfo.wpi0 += ke;
	    }
	    else if(TMath::Abs(ntpStdHep->IdHEP)==321 || //K+/-
		    ntpStdHep->IdHEP==130 ||     //K0_L
		    ntpStdHep->IdHEP==310 ||     //K0_S
		    ntpStdHep->IdHEP==311 ) {    //K0
	      fStdHepInfo.nkaon += 1;
	      if(ke>0.5) fStdHepInfo.nkaon_he += 1;
	      fStdHepInfo.ekaon += energy;
	      fStdHepInfo.pzkaon += pz;
	      fStdHepInfo.ptkaon += pt;
	      fStdHepInfo.wkaon += ke;

	      //edit by Anna within the kaon loop to get the kaon species details                                                                                     
	      if(ntpStdHep->IdHEP==321){      //K+                                                                                                                    
		fStdHepInfo.nkaonplus += 1;
		if(ke>0.5) fStdHepInfo.nkaonplus_he += 1;
		fStdHepInfo.ekaonplus += energy;
		fStdHepInfo.pzkaonplus += pz;
		fStdHepInfo.ptkaonplus += pt;
		fStdHepInfo.wkaonplus += ke;
	      }
	      if(ntpStdHep->IdHEP==-321){  //K-                                                                                                                       
		fStdHepInfo.nkaonminus += 1;
		if(ke>0.5) fStdHepInfo.nkaonminus_he += 1;
		fStdHepInfo.ekaonminus += energy;
		fStdHepInfo.pzkaonminus += pz;
		fStdHepInfo.ptkaonminus += pt;
		fStdHepInfo.wkaonminus += ke;
	      }
	      if(ntpStdHep->IdHEP==130){     //K0_L                                                                                                                   
		fStdHepInfo.nkaon0L += 1;
		if(ke>0.5) fStdHepInfo.nkaon0L_he += 1;
		fStdHepInfo.ekaon0L += energy;
		fStdHepInfo.pzkaon0L += pz;
		fStdHepInfo.ptkaon0L += pt;
		fStdHepInfo.wkaon0L += ke;
	      }
	      if(ntpStdHep->IdHEP==310){    //K0_S                                                                                                                    
		fStdHepInfo.nkaon0S += 1;
		if(ke>0.5) fStdHepInfo.nkaon0S_he += 1;
		fStdHepInfo.ekaon0S += energy;
		fStdHepInfo.pzkaon0S += pz;
		fStdHepInfo.ptkaon0S += pt;
		fStdHepInfo.wkaon0S += ke;
	      }
	      if(ntpStdHep->IdHEP==311){   //K0                                                                                                                       
		fStdHepInfo.nkaon0 += 1;
		if(ke>0.5) fStdHepInfo.nkaon0_he += 1;
		fStdHepInfo.ekaon0 += energy;
		fStdHepInfo.pzkaon0 += pz;
		fStdHepInfo.ptkaon0 += pt;
		fStdHepInfo.wkaon0 += ke;
	      }
              //end of edit by Anna to count individual kaon species                                                                                                    
	    }
	  }
	}
      }

      if(maxke>0){
	fStdHepInfo.wfs   /= maxke;
	fStdHepInfo.wgeant /= maxke;
	fStdHepInfo.wlep   /= maxke;
	fStdHepInfo.wprot  /= maxke;
	fStdHepInfo.wneut  /= maxke;
	fStdHepInfo.wpip   /= maxke;
	fStdHepInfo.wpim   /= maxke;
	fStdHepInfo.wpi0   /= maxke;
	fStdHepInfo.wkaon  /= maxke;
	fStdHepInfo.wlep_elec   /= maxke;   //Adding Anna Vars
	fStdHepInfo.wlep_muon   /= maxke;   //Adding Anna Vars
	fStdHepInfo.wlep_tau   /= maxke;   //Adding Anna Vars
	fStdHepInfo.wkaonplus  /= maxke;   //Adding Anna Vars
	fStdHepInfo.wkaonminus  /= maxke;   //Adding Anna Vars
	fStdHepInfo.wkaon0L  /= maxke;   //Adding Anna Vars
	fStdHepInfo.wkaon0S  /= maxke;   //Adding Anna Vars
	fStdHepInfo.wkaon0  /= maxke;   //Adding Anna Vars
      }

      if(lepmaxke>0){
	fStdHepInfo.lepwfs   /= lepmaxke;
	fStdHepInfo.lepwgeant /= lepmaxke;
	fStdHepInfo.lepwlep   /= lepmaxke;
	fStdHepInfo.lepwprot  /= lepmaxke;
	fStdHepInfo.lepwneut  /= lepmaxke;
	fStdHepInfo.lepwpip   /= lepmaxke;
	fStdHepInfo.lepwpim   /= lepmaxke;
	fStdHepInfo.lepwpi0   /= lepmaxke;
	fStdHepInfo.lepwkaon  /= lepmaxke;
	fStdHepInfo.lepwlep_elec   /= lepmaxke;   //Adding Anna Vars
	fStdHepInfo.lepwlep_muon   /= lepmaxke;   //Adding Anna Vars
	fStdHepInfo.lepwlep_tau   /= lepmaxke;   //Adding Anna Vars
	fStdHepInfo.lepwkaonplus  /= lepmaxke;   //Adding Anna Vars
	fStdHepInfo.lepwkaonminus  /= lepmaxke;   //Adding Anna Vars
	fStdHepInfo.lepwkaon0L  /= lepmaxke;   //Adding Anna Vars
	fStdHepInfo.lepwkaon0S  /= lepmaxke;   //Adding Anna Vars
	fStdHepInfo.lepwkaon0  /= lepmaxke;   //Adding Anna Vars
      }

      //Adding in Anna Vars
      if(maxke_nu>0){
	fStdHepInfo.wfs_nu   /= maxke_nu;
	fStdHepInfo.w_nu_elec   /= maxke_nu;
	fStdHepInfo.w_nu_muon   /= maxke_nu;
	fStdHepInfo.w_nu_tau   /= maxke_nu;
      }
      if(lepmaxke>0){
	fStdHepInfo.lepwfs_nu   /= lepmaxke_nu;
	fStdHepInfo.lepw_nu_elec   /= lepmaxke_nu;
	fStdHepInfo.lepw_nu_muon   /= lepmaxke_nu;
	fStdHepInfo.lepw_nu_tau   /= lepmaxke_nu;
      }
      //end Anna Vars    


      if(initialNuEnergy>0){
	fStdHepInfo.emfrac=fStdHepInfo.emenergy /initialNuEnergy;
	
      }else{
	fStdHepInfo.emfrac= ANtpDefVal::kDouble;
      }



      if(fStdHepInfo.pttot>0)
	fStdHepInfo.pttot = TMath::Sqrt(fStdHepInfo.pttot);
      if(fStdHepInfo.ptgeant>0) 
	fStdHepInfo.ptgeant = TMath::Sqrt(fStdHepInfo.ptgeant);
      if(fStdHepInfo.ptlep>0)
	fStdHepInfo.ptlep = TMath::Sqrt(fStdHepInfo.ptlep);
      if(fStdHepInfo.ptprot>0) 
	fStdHepInfo.ptprot = TMath::Sqrt(fStdHepInfo.ptprot);
      if(fStdHepInfo.ptneut>0) 
	fStdHepInfo.ptneut = TMath::Sqrt(fStdHepInfo.ptneut);
      if(fStdHepInfo.ptpip>0) 
	fStdHepInfo.ptpip = TMath::Sqrt(fStdHepInfo.ptpip);
      if(fStdHepInfo.ptpim>0) 
	fStdHepInfo.ptpim = TMath::Sqrt(fStdHepInfo.ptpim);
      if(fStdHepInfo.ptpi0>0) 
	fStdHepInfo.ptpi0 = TMath::Sqrt(fStdHepInfo.ptpi0);
      if(fStdHepInfo.ptkaon>0) 
	fStdHepInfo.ptkaon = TMath::Sqrt(fStdHepInfo.ptkaon);
      //adding Anna vars
      if(fStdHepInfo.pttot_nu>0)
	fStdHepInfo.pttot_nu = TMath::Sqrt(fStdHepInfo.pttot_nu);
      if(fStdHepInfo.pt_nu_elec>0)
	fStdHepInfo.pt_nu_elec = TMath::Sqrt(fStdHepInfo.pt_nu_elec);
      if(fStdHepInfo.pt_nu_muon>0)
	fStdHepInfo.pt_nu_muon = TMath::Sqrt(fStdHepInfo.pt_nu_muon);
      if(fStdHepInfo.pt_nu_tau>0)
	fStdHepInfo.pt_nu_tau = TMath::Sqrt(fStdHepInfo.pt_nu_tau);
      if(fStdHepInfo.ptlep_elec>0)
	fStdHepInfo.ptlep_elec = TMath::Sqrt(fStdHepInfo.ptlep_elec);
      if(fStdHepInfo.ptlep_muon>0)
	fStdHepInfo.ptlep_muon = TMath::Sqrt(fStdHepInfo.ptlep_muon);
      if(fStdHepInfo.ptlep_tau>0)
	fStdHepInfo.ptlep_tau = TMath::Sqrt(fStdHepInfo.ptlep_tau);
      if(fStdHepInfo.ptkaonplus>0) 
	fStdHepInfo.ptkaonplus = TMath::Sqrt(fStdHepInfo.ptkaonplus);
      if(fStdHepInfo.ptkaonminus>0) 
	fStdHepInfo.ptkaonminus = TMath::Sqrt(fStdHepInfo.ptkaonminus);
      if(fStdHepInfo.ptkaon0L>0) 
	fStdHepInfo.ptkaon0L = TMath::Sqrt(fStdHepInfo.ptkaon0L);
      if(fStdHepInfo.ptkaon0S>0) 
	fStdHepInfo.ptkaon0S = TMath::Sqrt(fStdHepInfo.ptkaon0S);
      if(fStdHepInfo.ptkaon0>0) 
	fStdHepInfo.ptkaon0 = TMath::Sqrt(fStdHepInfo.ptkaon0);
      //end Anna vars

      if(fStdHepInfo.leppttot>0)
	fStdHepInfo.leppttot = TMath::Sqrt(fStdHepInfo.leppttot);
      if(fStdHepInfo.lepptgeant>0) 
	fStdHepInfo.lepptgeant = TMath::Sqrt(fStdHepInfo.lepptgeant);
      if(fStdHepInfo.lepptlep>0)
	fStdHepInfo.lepptlep = TMath::Sqrt(fStdHepInfo.lepptlep);
      if(fStdHepInfo.lepptprot>0) 
	fStdHepInfo.lepptprot = TMath::Sqrt(fStdHepInfo.lepptprot);
      if(fStdHepInfo.lepptneut>0) 
	fStdHepInfo.lepptneut = TMath::Sqrt(fStdHepInfo.lepptneut);
      if(fStdHepInfo.lepptpip>0) 
	fStdHepInfo.lepptpip = TMath::Sqrt(fStdHepInfo.lepptpip);
      if(fStdHepInfo.lepptpim>0) 
	fStdHepInfo.lepptpim = TMath::Sqrt(fStdHepInfo.lepptpim);
      if(fStdHepInfo.lepptpi0>0) 
	fStdHepInfo.lepptpi0 = TMath::Sqrt(fStdHepInfo.lepptpi0);
      if(fStdHepInfo.lepptkaon>0) 
	fStdHepInfo.lepptkaon = TMath::Sqrt(fStdHepInfo.lepptkaon);
      //adding Anna vars
      if(fStdHepInfo.leppttot_nu>0)
	fStdHepInfo.leppttot_nu = TMath::Sqrt(fStdHepInfo.leppttot_nu);
      if(fStdHepInfo.leppt_nu_elec>0)
	fStdHepInfo.leppt_nu_elec = TMath::Sqrt(fStdHepInfo.leppt_nu_elec);
      if(fStdHepInfo.leppt_nu_muon>0)
	fStdHepInfo.leppt_nu_muon = TMath::Sqrt(fStdHepInfo.leppt_nu_muon);
      if(fStdHepInfo.leppt_nu_tau>0)
	fStdHepInfo.leppt_nu_tau = TMath::Sqrt(fStdHepInfo.leppt_nu_tau);
      if(fStdHepInfo.lepptlep_elec>0)
	fStdHepInfo.lepptlep_elec = TMath::Sqrt(fStdHepInfo.lepptlep_elec);
      if(fStdHepInfo.lepptlep_muon>0)
	fStdHepInfo.lepptlep_muon = TMath::Sqrt(fStdHepInfo.lepptlep_muon);
      if(fStdHepInfo.lepptlep_tau>0)
	fStdHepInfo.lepptlep_tau = TMath::Sqrt(fStdHepInfo.lepptlep_tau);
      if(fStdHepInfo.lepptkaonplus>0) 
	fStdHepInfo.lepptkaonplus = TMath::Sqrt(fStdHepInfo.lepptkaonplus);
      if(fStdHepInfo.lepptkaonminus>0) 
	fStdHepInfo.lepptkaonminus = TMath::Sqrt(fStdHepInfo.lepptkaonminus);
      if(fStdHepInfo.lepptkaon0L>0) 
	fStdHepInfo.lepptkaon0L = TMath::Sqrt(fStdHepInfo.lepptkaon0L);
      if(fStdHepInfo.lepptkaon0S>0) 
	fStdHepInfo.lepptkaon0S = TMath::Sqrt(fStdHepInfo.lepptkaon0S);
      if(fStdHepInfo.lepptkaon0>0) 
	fStdHepInfo.lepptkaon0 = TMath::Sqrt(fStdHepInfo.lepptkaon0);
      //end Anna vars

      //Now include Anna vars to determine the tau decay channel for this event:
      if(fStdHepInfo.nutauevent==1 || fStdHepInfo.nutauevent==-1){
	if(fStdHepInfo.lepnlep_elec==0 && fStdHepInfo.lepnlep_muon>0 ) fStdHepInfo.nutauchannel=1;
	else if(fStdHepInfo.lepnlep_elec>0 && fStdHepInfo.lepnlep_muon==0 ) fStdHepInfo.nutauchannel=2;
	else if(fStdHepInfo.lepnlep_elec==0 && fStdHepInfo.lepnlep_muon==0 && fStdHepInfo.lepnlep_elec==0 && fStdHepInfo.lepnpi0>0 ) fStdHepInfo.nutauchannel=3; 
	else if(fStdHepInfo.lepnlep_elec==0 && fStdHepInfo.lepnlep_muon==0 && fStdHepInfo.lepnlep_elec==0 && fStdHepInfo.lepnpi0==0 ) fStdHepInfo.nutauchannel=4; 
	else { 
        fStdHepInfo.nutauchannel=5; 
	}
       }
      //end Anna vars

    }
  }                 	
}




//this class will determine if a given particle has a parent which is from the original neutrino, pass in parent of current particle
bool isLeptonic(int parent, TClonesArray& stdhepArray)  
{

  NtpMCStdHep *ntpStdHepTemp = 0;  
  ntpStdHepTemp = dynamic_cast<NtpMCStdHep *>(stdhepArray[parent]);
  int parentID = ntpStdHepTemp->IdHEP;

  if ((abs(parentID)==12||abs(parentID)==14||abs(parentID)==16)&&ntpStdHepTemp->IstHEP==0)return true;  //related to the original neutrino!
  if (ntpStdHepTemp->parent[0]==-1||ntpStdHepTemp->parent[1]==-1)return false;  //at end of chain, not parent so exit

  if (isLeptonic(ntpStdHepTemp->parent[0],stdhepArray)) return true; 
  if (ntpStdHepTemp->parent[0]!=ntpStdHepTemp->parent[1])  //check each parent branch if they are different
    {
      if (isLeptonic(ntpStdHepTemp->parent[1],stdhepArray)) return true;
    }

  return false;

}


//this class will determine if a given particle has a parent with final state 1 or 205 of type electron, pi0, or photon
//returns false if no parent is wanted type
bool checkParent(int parent, TClonesArray& stdhepArray)  
{

  NtpMCStdHep *ntpStdHepTemp = 0;  
  ntpStdHepTemp = dynamic_cast<NtpMCStdHep *>(stdhepArray[parent]);
  int parentID = ntpStdHepTemp->IdHEP;

  if(!(ntpStdHepTemp->IstHEP==1||ntpStdHepTemp->IstHEP==205))return false;

  if ((abs(parentID)==22||abs(parentID)==111||abs(parentID)==11))return true;  //
  
  if (checkParent(ntpStdHepTemp->parent[0],stdhepArray)) return true; 
  if (ntpStdHepTemp->parent[0]!=ntpStdHepTemp->parent[1])  //check each parent branch if they are different
    {
      if (checkParent(ntpStdHepTemp->parent[1],stdhepArray)) return true;
    }

  return false;

}



void countpi0(int index, TClonesArray& stdhepArray, StdHepInfo &fStdHepInfo)
{
  //cout<<"counting pi0"<<endl;
  Int_t nStdHep = stdhepArray.GetEntries();
  NtpMCStdHep *ntpStdHep = 0;
  for(int i=0; i<nStdHep; i++){
    ntpStdHep = dynamic_cast<NtpMCStdHep *>(stdhepArray[i]);
    //check that stdhep entry corresponds to this mc event
    if(ntpStdHep->mc!=index) continue;
    
    if(ntpStdHep->IdHEP==111){//a pi0 is spotted
      //check if this pi0 comes from a resonance decay or DIS
      bool resdecay = true;
      for (int ip = ntpStdHep->parent[0]; ip<=ntpStdHep->parent[1]; ip++){
	if (ip==-1) {
	  resdecay = false;
	  break;
	}
	NtpMCStdHep *parent = dynamic_cast<NtpMCStdHep *>(stdhepArray[ip]);
	resdecay = resdecay && (parent->IstHEP==3);
	if (!resdecay) break;
      }
      if (resdecay) fStdHepInfo.epi0_neugen += TMath::Abs(ntpStdHep->p4[3]);

      Double_t etmp = 0; // count epi0 in the child list
      for (int ic = ntpStdHep->child[0]; ic<=ntpStdHep->child[1]; ic++){
	if (ic==-1) break;
	NtpMCStdHep *child = dynamic_cast<NtpMCStdHep *>(stdhepArray[ic]);
	if (child->IdHEP==111) etmp += TMath::Abs(child->p4[3]);
      }
      if (ntpStdHep->IstHEP!=205 && //not a pi0 decay
	  ntpStdHep->IstHEP!=1205 &&
	  ntpStdHep->child[0]!=-1 &&
	  ntpStdHep->child[1]!=-1){
	fStdHepInfo.epi0_abs += TMath::Abs(ntpStdHep->p4[3]) - etmp;
      }

      bool pi0pro = true;
      for (int ip = ntpStdHep->parent[0]; ip<=ntpStdHep->parent[1]; ip++){
	if (ip==-1) {
	  pi0pro = false;
	  break;
	}
 	NtpMCStdHep *parent = dynamic_cast<NtpMCStdHep *>(stdhepArray[ip]);     
	pi0pro = pi0pro&&(parent->IstHEP == 14 && parent->IdHEP != 111);
	if (!pi0pro) break;
      }
      if (pi0pro) {
	fStdHepInfo.epi0_intranuke += TMath::Abs(ntpStdHep->p4[3]);
	fStdHepInfo.epi0_total += TMath::Abs(ntpStdHep->p4[3]);
      }
      else if (ntpStdHep->IstHEP == 1){
	fStdHepInfo.epi0_total += TMath::Abs(ntpStdHep->p4[3]);
      }

      bool decay = true;
      for (int ip = ntpStdHep->parent[0]; ip<=ntpStdHep->parent[1]; ip++){
	if (ip==-1) {
	  decay = false;
	  break;
	}
	NtpMCStdHep *parent = dynamic_cast<NtpMCStdHep *>(stdhepArray[ip]);     
	decay = decay&&(parent->IstHEP == 1 &&
			(ntpStdHep->IstHEP == 205||ntpStdHep->IstHEP == 1205) 
			&& parent->IdHEP != 111);
	if (!decay) break;
      }
      if (decay) {
	fStdHepInfo.epi0_intranuke += TMath::Abs(ntpStdHep->p4[3]);
	fStdHepInfo.epi0_total += TMath::Abs(ntpStdHep->p4[3]);
      }
    }//find a pi0
    else if(ntpStdHep->IdHEP==22){//a gamma is spotted
      bool decay = true;
      for (int ip = ntpStdHep->parent[0]; ip<=ntpStdHep->parent[1]; ip++){
	if (ip==-1) {
	  decay = false;
	  break;
	}
	NtpMCStdHep *parent = dynamic_cast<NtpMCStdHep *>(stdhepArray[ip]);     
	decay = decay&&(parent->IstHEP == 1 && 
			(ntpStdHep->IstHEP == 205||ntpStdHep->IstHEP == 1205)
			&&parent->IdHEP != 111);
	if (!decay) break;
      }
      if (decay) {
	fStdHepInfo.epi0_intranuke += TMath::Abs(ntpStdHep->p4[3]);      
	fStdHepInfo.epi0_total += TMath::Abs(ntpStdHep->p4[3]);
      }
    }
  }//loop through the particle list
//  cout<<"pi0 summary "<<endl;
//  cout<<"epi0_total "<<fStdHepInfo.epi0_total<<endl;
//  cout<<"epi0_neugen "<<fStdHepInfo.epi0_neugen<<endl;
//  cout<<"epi0_intranuke "<<fStdHepInfo.epi0_intranuke<<endl;
//  cout<<"epi0_abs "<<fStdHepInfo.epi0_abs<<endl;

}

void calbaryonpt(int index, TClonesArray& stdhepArray, StdHepInfo &fStdHepInfo)
{

  vector<int> pid;
  vector<double> px;
  vector<double> py;
  vector<double> pz;
  vector<double> eng;
  vector<double> mass;

  double hspx = 0.;
  double hspy = 0.;
  double hspz = 0.;
  double hse  = 0.;
  double beta = 0.;
  double gamma = 0.;

  //cout<<"counting pi0"<<endl;
  Int_t nStdHep = stdhepArray.GetEntries();
  NtpMCStdHep *ntpStdHep = 0;
  for(int i=0; i<nStdHep; i++){
    ntpStdHep = dynamic_cast<NtpMCStdHep *>(stdhepArray[i]);
    //check that stdhep entry corresponds to this mc event
    if(ntpStdHep->mc!=index) continue;
    
    if(ntpStdHep->IstHEP==3){//resonance decay
      if (TMath::Abs(ntpStdHep->p4[3])>hse){
	hspx = ntpStdHep->p4[0];
	hspy = ntpStdHep->p4[1];
	hspz = ntpStdHep->p4[2];
	hse  = TMath::Abs(ntpStdHep->p4[3]);
      }
    }
    else {
      //check if this particle comes from a resonance decay or DIS
      bool resdecay = true;
      for (int ip = ntpStdHep->parent[0]; ip<=ntpStdHep->parent[1]; ip++){
	if (ip==-1) {
	  resdecay = false;
	  break;
	}
	NtpMCStdHep *parent = dynamic_cast<NtpMCStdHep *>(stdhepArray[ip]);
	resdecay = resdecay && (parent->IstHEP==3);
	if (!resdecay) break;
      }
      if (resdecay&&ntpStdHep->IdHEP) {
	pid.push_back(ntpStdHep->IdHEP);
	px.push_back(ntpStdHep->p4[0]);
	py.push_back(ntpStdHep->p4[1]);
	pz.push_back(ntpStdHep->p4[2]);
	eng.push_back(ntpStdHep->p4[3]);
	mass.push_back(ntpStdHep->mass);
      }
    }
  }

  double maxbaryone = 0;
  double Phs = sqrt(pow(hspx,2)+pow(hspy,2)+pow(hspz,2));

  if (Phs>=0.&&hse>Phs){
    beta = Phs/hse;
    gamma = hse/sqrt(pow(hse,2)-pow(Phs,2));
  }


  for (unsigned ipar = 0; ipar<pid.size(); ipar++){//loop through hadrons
    double ptot = sqrt(pow(px[ipar],2)+pow(py[ipar],2)+pow(pz[ipar],2));
    double Pz = 0;
    if (Phs>0) Pz = (px[ipar]*hspx+py[ipar]*hspy+pz[ipar]*hspz)/Phs;
    double pt = pow(ptot,2)-pow(Pz,2);
    Pz = gamma*(Pz-beta*eng[ipar]);
    if (pt>0.00000001) pt = sqrt(pt);
    else pt = 0;
    if (TMath::Abs(pid[ipar])==2212||TMath::Abs(pid[ipar])==2112){
      if (eng[ipar]>maxbaryone){
	maxbaryone = eng[ipar];
	fStdHepInfo.baryonpt = pt;
	fStdHepInfo.baryonxf = Pz;
      }
    }
    if (pid[ipar]==111) fStdHepInfo.npi0_neugen++;
    fStdHepInfo.totpt += pt;
  }
  fStdHepInfo.nhad = int(pid.size());
}

void HughCode(NtpMCTruth* mctruth, TClonesArray& stdhepArray, StdHepInfo& shi)
{
   //find the mc information for this event
   NtpMCStdHep* ntpStdHep = 0;
                                                                                                             
   bool foundFirst = false;
   int index = 0;
                                                                                                             
   double Etot = 0.0;
   double EMtot = 0.0;
   int npi0 = 0;
   int nch = 0;
   int ntot = 0;
                                                                                                             
   for(int i=mctruth->stdhep[0];i<=mctruth->stdhep[1]&&!foundFirst;i++){
      ntpStdHep = dynamic_cast<NtpMCStdHep *>(stdhepArray[i]);
                                                                                                             
      if(ntpStdHep->mc!=mctruth->index)  cout<<"Big Problem"<<endl;
                                                                                                             
      if(!foundFirst && ntpStdHep->IstHEP == 3) foundFirst = true;
      else continue;
      index = i;
   }
                                                                                                             
   if(foundFirst){
     int start = ntpStdHep->child[0];
     int end  =  ntpStdHep->child[1];
                                                                                                             
     CalcDaughters(start, end, Etot, EMtot, npi0, nch, ntot, stdhepArray);
                                                                                                             
     Etot = Etot - 0.935;
   }

   shi.e_total = Etot;
   shi.em_total = EMtot;
   shi.npi0_jbhg = npi0;
   shi.nch_jbhg = nch;
   shi.ntot_jbhg = ntot;
}
                                                                                                             
void CalcDaughters(int start, int end, double &E, double &EmE, int &npi0, int &nch, int &ntot, TClonesArray&
stdhepArray)
{
   NtpMCStdHep* ntpStdHep = 0;
                                                                                                             
   for(int j = start; j <= end; j++){
     ntpStdHep = dynamic_cast<NtpMCStdHep *>(stdhepArray[j]);
     int id = ntpStdHep->IdHEP;
     int ist = ntpStdHep->IstHEP;
     if(ist == 3){
       int st = ntpStdHep->child[0];
       int en =  ntpStdHep->child[1];
       CalcDaughters(st, en, E, EmE, npi0, nch, ntot, stdhepArray);
     }
                                                                                                             
     float Etemp = ntpStdHep->p4[3];
                                                                                                             
     E += Etemp;
     if(id == 111 || id == 22){
       EmE += Etemp;
       npi0++;
     }
     if(id == 211 || id == -211)   nch++;
     ntot++;
  }
}
                                                                                                             

