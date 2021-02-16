#include <iostream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "Conventions/DetectorType.h"
#include "NueAna/Reweight/NueRW.h" 

using namespace std;

void AddRWFiles(const char* fnames);
void MakeChallengeFile(const char* ndf, const char *fdf);

void AddRWFiles(const char* fnames)
{
   
//   const int NROWS=2675;
//   const int NROWS=5103;
   const int NROWS=4913;
   const float NEARPOTPERFILE=2.4e-7*550;
   const float FARPOTPERFILE=6.5;

   int nnfiles[NROWS];
   int nfbeam[NROWS];
   int nfnue[NROWS];
   int nfnutau[NROWS];

   NueRW *nrw[NROWS];
   NueRW *fbrw[NROWS];
   NueRW *fnrw[NROWS];
   NueRW *ftrw[NROWS];
   NueRW *totfrw[NROWS];
   for(int i=0;i<NROWS;i++){
     nrw[i]=new NueRW();
     nrw[i]->fDet=DetectorType::kNear;
     nrw[i]->fFileType=NueRW::kBEAM;
     nrw[i]->randrow=i;

     fbrw[i]=new NueRW();
     fbrw[i]->fDet=DetectorType::kFar;
     fbrw[i]->fFileType=NueRW::kBEAM;
     fbrw[i]->randrow=i;

     fnrw[i]=new NueRW();
     fnrw[i]->fDet=DetectorType::kFar;
     fnrw[i]->fFileType=NueRW::kNUE;
     fnrw[i]->randrow=i;

     ftrw[i]=new NueRW();
     ftrw[i]->fDet=DetectorType::kFar;
     ftrw[i]->fFileType=NueRW::kTAU;
     ftrw[i]->randrow=i;

     totfrw[i]=new NueRW();
     totfrw[i]->fDet=DetectorType::kFar;
     totfrw[i]->fFileType=NueRW::kAGG;
     totfrw[i]->randrow=i;

     nnfiles[i]=0;	
     nfbeam[i]=0;
     nfnue[i]=0;
     nfnutau[i]=0;
   }

   TChain *rwchain = new TChain("rwtree","rwtree");
   rwchain->Add(fnames);
   
   NueRW *looprw = new NueRW();
   rwchain->SetBranchAddress("NueRW",&looprw);
   //loop over chain entries
   for(int z=0;z<rwchain->GetEntries();z++){
     if(z%10000==0){
       cout<<"On entry "<<z<<endl;
     }
     rwchain->GetEntry(z);
     if(looprw->randrow>=NROWS||looprw->randrow<0){
       //cout<<"ERROR: randrow="<<fdl->randrow<<endl;
       continue;
     }
     if(looprw->fDet==0){
       continue;
     }
     if(looprw->fDet==DetectorType::kNear){
       *(nrw[looprw->randrow])=(*looprw)+(*(nrw[looprw->randrow]));
       nnfiles[looprw->randrow]+=looprw->nfiles;
     }
     else if(looprw->fDet==DetectorType::kFar){       
       if(looprw->fFileType==NueRW::kBEAM){
	 *(fbrw[looprw->randrow])=(*looprw)+(*(fbrw[looprw->randrow]));
	 nfbeam[looprw->randrow]+=looprw->nfiles;
       }
       else if(looprw->fFileType==NueRW::kNUE){
	 *(fnrw[looprw->randrow])=(*looprw)+(*(fnrw[looprw->randrow]));
	 nfnue[looprw->randrow]+=looprw->nfiles;
       }
       else if(looprw->fFileType==NueRW::kTAU){
	 *(ftrw[looprw->randrow])=(*looprw)+(*(ftrw[looprw->randrow]));
	 nfnutau[looprw->randrow]+=looprw->nfiles;
       }
     }
     else{
       cout<<"UNKNOWN FILE TYPE!!!!!!!!!!!!!"<<endl;
     }
   }

   //   char newfname[100];
   //   sprintf(newfname,"rwtree-%s.root",DET);
   //   TFile *f=new TFile(newfname,"RECREATE");
   TFile *f=new TFile("rwtree.root","RECREATE");
   TTree *nttree = new TTree("rwtree_ND_norm","rwtree_ND_norm");
   TTree *fttree = new TTree("rwtree_FD_norm","rwtree_FD_norm");
   TTree *ntree = new TTree("rwtree_ND","rwtree_ND");
   TTree *fbtree = new TTree("rwtree_FDBeam","rwtree_FDBeam");
   TTree *fntree = new TTree("rwtree_FDNue","rwtree_FDNue");
   TTree *ftautree = new TTree("rwtree_FDTau","rwtree_FDTau");
   NueRW *temp = new NueRW();
   nttree->Branch("NueRW","NueRW",&temp,64000,99);
   fttree->Branch("NueRW","NueRW",&temp,64000,99);
   ntree->Branch("NueRW","NueRW",&temp,64000,99);
   fbtree->Branch("NueRW","NueRW",&temp,64000,99);
   fntree->Branch("NueRW","NueRW",&temp,64000,99);
   ftautree->Branch("NueRW","NueRW",&temp,64000,99);

   for(int i=0;i<NROWS;i++){
     ntree->SetBranchAddress("NueRW",&(nrw[i]));
     fbtree->SetBranchAddress("NueRW",&(fbrw[i]));
     fntree->SetBranchAddress("NueRW",&(fnrw[i]));
     ftautree->SetBranchAddress("NueRW",&(ftrw[i]));

     nrw[i]->nfiles=nnfiles[i];
     ntree->Fill();
     fbrw[i]->nfiles=nfbeam[i];
     fbtree->Fill();
     fnrw[i]->nfiles=nfnue[i];
     fntree->Fill();
     ftrw[i]->nfiles=nfnutau[i];
     ftautree->Fill();

     //renormalize total trees to exposure
     //near's easy, just divide by POT*files
     if(nnfiles[i]!=0){
       *(nrw[i])=*(nrw[i])/(NEARPOTPERFILE*nnfiles[i]);
     }
     else{
       cout<<"Row "<<i<<" has 0 near files, not normalizing"<<endl;
     }
     nttree->SetBranchAddress("NueRW",&(nrw[i]));
     nttree->Fill();
     
     //far's harder because there could be different number of each file type

     if(nfbeam[i]==0||nfnue[i]==0||nfnutau[i]==0){
       cout<<"Missing some FD files, will not normalize FD"
	   <<" beam: "<<nfbeam[i]
	   <<" nue "<<nfnue[i]
	   <<" tau "<<nfnutau[i]<<endl;
     }
     else{
       *(fbrw[i])=*(fbrw[i])/(1.*FARPOTPERFILE*nfbeam[i]);
       *(fnrw[i])=*(fnrw[i])/(1.*FARPOTPERFILE*nfnue[i]);
       *(ftrw[i])=*(ftrw[i])/(1.*FARPOTPERFILE*nfnutau[i]);
       *(totfrw[i])=*(fbrw[i])+*(fnrw[i])+*(ftrw[i]);
       //but now the mf nuebeams are off since we get them from both beam and nue files
       //as is the total background since we added in too much nuebeams
       //first subtract the wrong nuebeam bg
       totfrw[i]->nbg=(totfrw[i]->nbg-totfrw[i]->nnueb);
       //figure out the right nuebeam bg
       totfrw[i]->nnueb=(fbrw[i]->nnueb+fnrw[i]->nnueb)/2.;
       //add back the right nuebeam bg
       totfrw[i]->nbg=(totfrw[i]->nbg+totfrw[i]->nnueb);
       
       //and now fix the by energy entries
       for(int j=0;j<totfrw[i]->EBINS;j++){
	 //first subtract the wrong nuebeam bg
	 totfrw[i]->nbgE[j]=(totfrw[i]->nbgE[j]-totfrw[i]->nnuebE[j]);
	 //figure out the right nuebeam bg
	 totfrw[i]->nnuebE[j]=(fbrw[i]->nnuebE[j]+fnrw[i]->nnuebE[j])/2.;
	 //add back the right nuebeam bg
	 totfrw[i]->nbgE[j]=(totfrw[i]->nbgE[j]+totfrw[i]->nnuebE[j]);
       }
     }

     fttree->SetBranchAddress("NueRW",&(totfrw[i]));
     fttree->Fill();
   }
   f->Write();
   
}

void MakeChallengeFile(const char* ndf, const char *fdf)
{
  TChain *nnue = new TChain("ana_nue_filt","ana_nue_filt");
  nnue->Add(ndf);

  TChain *fnue = new TChain("ana_nue_filt","ana_nue_filt");
  fnue->Add(fdf);

  TFile f("challengespec.root","RECREATE");
  f.cd();

  TH1F ndspec("ndspec","ndspec",20,0,20);
  nnue->Draw("NueRecord.srevent.energyGeV>>ndspec","","goff");
  TH1F fdh("fdspec","fdspec",20,0,20);
  fnue->Draw("NueRecord.srevent.energyGeV>>fdspec","","goff");
  
  f.Write();
}

