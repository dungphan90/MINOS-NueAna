////////////////////////////////////////////////////////////////////////
// $Id: NueReweightFast.C,v 1.2 2005/06/09 13:19:29 vahle Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle
////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "NueAna/NueRecord.h"
#include "NueAna/Reweight/NueRW.h"
#include "NueAna/NueRWHelpers.h"
#include "Conventions/DetectorType.h"
#include "MCReweight/NeugenWeightCalculator.h"
#include "MCReweight/MCReweight.h"
#include "MCReweight/MCEventInfo.h"
#include "MCReweight/ReweightHelpers.h"
#include "Registry/Registry.h"

using namespace std;

vector<vector<float> >ReadRandom(const char* fname);
void SetRandRow(unsigned int row, const vector<vector<float> > &rvec, NueRW *rw);
void FillMCEventInfo(const ANtpTruthInfoBeam &mc, MCEventInfo &ei);
const unsigned int NCOLS=25;


int NueReweightFast(const char* rname, TChain *nue, const char* outname)
{
//Run on multiple nd trees at once.  Nothing gets written out if there
//are no entries in the filtered tree, then it's hard to
//figure out theproper normalization!!!!!

   const int NFILES=nue->GetNtrees();

  //Set branch address to a NueRecord
  //this assumes we're reading in a FilteredPID TREE
  //it will loop over every event in the tree
  //it makes no PID cuts
  //you must be using a ana_nue_filt tree!
  NueRecord *nr = new NueRecord();
  nue->SetBranchAddress("NueRecord",&nr);
  const float NENTRIES=nue->GetEntries();
  if(NENTRIES<1){
    cout<<"ERROR No entries in chain"<<endl;
    return -1;
  }
  //make a mcr
  MCReweight mcr(MCReweight::Instance());
  // make a WeightCalculator and add it to the MCReweight singleton
  NeugenWeightCalculator *n=new NeugenWeightCalculator();
  mcr.AddWeightCalculator(n);
  
  //read in the list of parameter values to loop over
  vector<vector<float> >rvec = ReadRandom(rname);
  const unsigned int NROWS = rvec.size();
  
  cout<<"Events*Loops="<<NENTRIES*NROWS<<endl;

  //get first entry and fill rw's with initial info
  nue->GetEntry(0);  
  vector<NueRW *>rw(NROWS);
  int run = nr->GetHeader().GetRun();
  int srun =  nr->GetHeader().GetSubRun();
  DetectorType::Detector_t det=nr->GetHeader().GetVldContext().GetDetector();
  int iflv = (int)((run%100000-run%1000)/10000);
  NueRW::FileType_t ftype;
  if(iflv==0){
    ftype = NueRW::kBEAM;
  }
  else if(iflv==1){
    ftype = NueRW::kNUE;
  }
  else if(iflv==3){
    ftype = NueRW::kTAU;
  }
  else{
    ftype = NueRW::kUnknown;
  }
  for(unsigned int i=0;i<rw.size();i++){
    rw[i]=new NueRW();
    rw[i]->randrow=i;
    rw[i]->fRun = run;
    rw[i]->fSubRun = srun;
    rw[i]->fDet = det;
    rw[i]->fFileType=ftype;
    rw[i]->nfiles=NFILES;
    SetRandRow(i,rvec,rw[i]);
  }


  //loop over events
  int cntr=0;
  for(int z=0;z<NENTRIES;z++){
    if(z%1000==0){
      cout<<"On Entry "<<z<<endl;
    }
    nue->GetEntry(z);
    //make event regsitry
    //    Registry evreg;
    //    ReweightHelpers::EventRegistryFilla(&(nr->mctrue),evreg);
    //make an MCEventInfo to hold mctruth for this event
    MCEventInfo event;
    event.UseStoredXSec(true);
    FillMCEventInfo(nr->mctrue,event);
    //make an empty(for now) nuparent
    NuParent *np=0;

    //loop over reweight rows
    for(unsigned int i=0;i<rw.size();i++){
      cntr++;
      if(cntr%500000==0){
	cout<<"On loop "<<cntr<<endl;
      }
    //    for(unsigned int i=0;i<50;i++){
      if(nr->GetHeader().GetEventNo()==0){
	rw[i]->nsnarls++;
      }
      rw[i]->nevents++;   
      rw[i]->nacc++;      
      
      //fill the reweight config
      Registry rwtconfig;
      rwtconfig.UnLockValues();
      rwtconfig.UnLockKeys();
      rwtconfig.Clear();
      rwtconfig.Set("neugen:use_scale_factors",1);
      rwtconfig.Set("neugen:ma_qe",rw[i]->qel_ma);
      rwtconfig.Set("neugen:ma_res",rw[i]->res_ma);
      rwtconfig.Set("neugen:kno_r112",rw[i]->kno_r112);
      rwtconfig.Set("neugen:kno_r122",rw[i]->kno_r122);
      rwtconfig.Set("neugen:kno_r132",rw[i]->kno_r132);
      rwtconfig.Set("neugen:kno_r142",rw[i]->kno_r142);
      rwtconfig.Set("neugen:kno_r113",rw[i]->kno_r113);
      rwtconfig.Set("neugen:kno_r123",rw[i]->kno_r123);
      rwtconfig.Set("neugen:kno_r133",rw[i]->kno_r133);
      rwtconfig.Set("neugen:kno_r143",rw[i]->kno_r143);
      rwtconfig.LockValues();
      rwtconfig.LockKeys();
      
      //do reweighting--
      // now all that's left is to compute the weight as follows:
      float reweight = mcr.ComputeWeight(&event,np,&rwtconfig);
      
      //oscillate
      float oscprob = 1.;
      if(nr->GetHeader().GetVldContext().GetDetector()==DetectorType::kFar){
	float theta23 = TMath::ASin(sqrt(rw[i]->ss2th))/2.;
	oscprob = NueRWHelpers::Oscillate(&(nr->mctrue),735.,
					  rw[i]->dm2,theta23,rw[i]->UE32);
      }

      //      cout<<"oscprob is "<<oscprob<<endl;
      int ebin=rw[i]->FindEBin(nr->srevent.energyGeV);
      float scale=reweight*oscprob;

      int iaction = nr->mctrue.interactionType;
      int inu = nr->mctrue.nuFlavor;
      int inunoosc = nr->mctrue.nonOscNuFlavor;
      
      //fill the elements in rw
      if(iaction==0){
	//neutral current;
	rw[i]->nnc+=scale;
	rw[i]->nbg+=scale;
	rw[i]->nncE[ebin]+=scale;
	rw[i]->nbgE[ebin]+=scale;
      }
      else{
	if(abs(inu)==12){
	  if(abs(inunoosc)==12){
	     rw[i]->nnueb+=scale;
	     rw[i]->nbg+=scale;
	     rw[i]->nnuebE[ebin]+=scale;
	     rw[i]->nbgE[ebin]+=scale;
	  }
	  else if(abs(inunoosc)==14){
	    rw[i]->nsig+=scale;
	    rw[i]->nsigE[ebin]+=scale;
	  }
	}
	else if(abs(inu)==14){
	  rw[i]->nnumu+=scale;
	  rw[i]->nbg+=scale;
	  rw[i]->nnumuE[ebin]+=scale;
	  rw[i]->nbgE[ebin]+=scale;
	}
	else if(abs(inu)==16){
	  rw[i]->nnutau+=scale;
	  rw[i]->nbg+=scale;
	  rw[i]->nnutauE[ebin]+=scale;
	  rw[i]->nbgE[ebin]+=scale;
	}
      }
    }
  }

  //write the rw's to a tree.
  TFile *f = new TFile(outname,"RECREATE");
  TTree *tree = new TTree("rwtree","rwtree");
  NueRW *nrw = new NueRW();
  tree->Branch("NueRW","NueRW",&nrw,64000,99);
  delete nrw;
  nrw=0;
  
  for(unsigned int i=0;i<rw.size();i++){
    tree->SetBranchAddress("NueRW",&rw[i]);
    tree->Fill();
  }

  f->Write();
  f->Close();

  return 0;
}



//......................................................................
vector<vector<float> >ReadRandom(const char* fname)
{
  vector< vector<float> > rval;
  ifstream in(fname);
  if(!in){
    cout<<"ERROR Could not open random number file in priviate area"<<endl
	<<"Tried: "<<fname<<endl;
    return rval;
  }

  float r;
  in>>r;
  while(!in.eof()){
    vector<float> rands;
    rands.push_back(r);
    for(unsigned int i=1;i<NCOLS;i++){
      in>>r;
      rands.push_back(r);
    }
    if(rands.size()!=NCOLS){
      cout<<"ERROR row  only has "<<rands.size()<<" entries"<<endl;
    }
    rval.push_back(rands);
    in>>r;
  }
  
  return rval;
}

void SetRandRow(unsigned int row, const vector<vector<float> > &rvec, NueRW *rw)
{
  if(row>rvec.size()){
    cout<<"ERROR Can not fill row "<<row<<" rvec only has "
	<<rvec.size()<<" rows"<<endl;
    return;
  }

   rw->qel_ma=fabs(rvec[row][0]);
   rw->res_ma=fabs(rvec[row][1]);
   rw->coh_ma=fabs(rvec[row][2]);
   rw->qel_fa0=-1.*fabs(rvec[row][3]);
   rw->qel_eta=fabs(rvec[row][4]);
   rw->res_omega=fabs(rvec[row][5]);
   rw->res_z=fabs(rvec[row][6]);
   rw->coh_r0=fabs(rvec[row][7]);
   rw->coh_rei=fabs(rvec[row][8]);
   rw->kno_a1=fabs(rvec[row][9]);
   rw->kno_a2=fabs(rvec[row][10]);
   rw->kno_a3=fabs(rvec[row][11]);
   rw->kno_a4=fabs(rvec[row][12]);
   rw->kno_b=fabs(rvec[row][13]);

   //use same scale for all kno_r* variables
   rw->kno_r112=fabs(rvec[row][14]);
   rw->kno_r122=fabs(rvec[row][14]);
   rw->kno_r132=fabs(rvec[row][14]);
   rw->kno_r142=fabs(rvec[row][14]);   
   rw->kno_r113=fabs(rvec[row][14]);
   rw->kno_r123=fabs(rvec[row][14]);
   rw->kno_r133=fabs(rvec[row][14]);
   rw->kno_r143=fabs(rvec[row][14]);
   rw->dm2=fabs(rvec[row][22]);
   if(rvec[row][23]>1){
     rw->ss2th=1.-(rvec[row][23]-1.);
   }
   else{
     rw->ss2th=fabs(rvec[row][23]);
   }
   rw->UE32=fabs(rvec[row][24]);
   return;
}

//......................................................................

void FillMCEventInfo(const ANtpTruthInfoBeam &mc, MCEventInfo &ei)
{


ei.nuE=mc.nuEnergy;
ei.nuPx=mc.nuDCosX*sqrt(mc.nuEnergy*mc.nuEnergy);
ei.nuPy=mc.nuDCosY*sqrt(mc.nuEnergy*mc.nuEnergy);
ei.nuPz=mc.nuDCosZ*sqrt(mc.nuEnergy*mc.nuEnergy);

ei.tarE=mc.targetEnergy;
ei.tarPx=mc.targetPX;
ei.tarPy=mc.targetPY;
ei.tarPz=mc.targetPZ;

ei.y=mc.hadronicY;
ei.x=mc.bjorkenX;
ei.q2=mc.q2;
ei.w2=mc.w2;

ei.iaction=mc.interactionType;
ei.inu=mc.nuFlavor;
ei.iresonance=mc.resonanceCode;
ei.initial_state=mc.initialState;

   Int_t nucleus=1;
   switch ((int)(mc.atomicNumber)) {
   case 1:
      nucleus=0;   // free nucleon
      break;
   case 6:
      nucleus=274; // carbon
      break;
   case 8:
      nucleus=284; // oxygen
      break;
   case 26:
      nucleus=372; // mcon
      break;
   default:
      nucleus=1;   // unknown
      break;
   }
ei.nucleus=nucleus;
ei.had_fs=mc.hadronicFinalState;

}
