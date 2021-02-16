////////////////////////////////////////////////////////////////////////
// $Id: CompareMST.cxx,v 1.2 2008/11/19 18:22:51 rhatcher Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TFile.h"
#include "Conventions/Detector.h"
#include "NueAna/CompareMST.h"
#include "NueAna/NueRecord.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "MessageService/MsgService.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "HistMan/HistMan.h"

const int CALend=110;
const int SM1end=230;
const int SM2end=484;

JOBMODULE(CompareMST, "CompareMST",
          "Makes hists to compare MST variables between detectors and truth classes");
CVSID("$Id:");
//......................................................................

CompareMST::CompareMST():
  counter(0)
{}

//......................................................................

CompareMST::~CompareMST()
{}

//......................................................................
void CompareMST::BeginJob()
{

  if(counter%1000==0){
    cout<<"On entry "<<counter<<endl;
  }
  counter++;
 
 static HistMan *hm = new HistMan("mstcomp");
  for(int d=0;d<2;d++){
    char det[2];
    if(d==0){
	  sprintf(det,"f");
    }
    else if(d==1){
      sprintf(det,"n");
    }
    else{
      sprintf(det,"u");
    }
    for(int b=0;b<4;b++){
      char bgc[6];
      if(b==0){
	sprintf(bgc,"nue");
      }
      else if(b==1){
	sprintf(bgc,"numu");
      }
      else if(b==2){
	sprintf(bgc,"nutau");
      }
      else if(b==3){
	sprintf(bgc,"nc");
      }
      else{
	sprintf(bgc,"u");
      }
      for(int r=0;r<4;r++){
	char id[20];
	sprintf(id,"%s_%s_%d",det,bgc,r+1001);
	char nemsg[50];
	char newtot[50];
	char nentot[50];
	char new1[50];
	char nen1[50];
	char new1nn1[50];
	char nesmtot[50];
	char nesm1[50];
	char neal[50];
	char nep[50];
	char new4[50];
	char nen4[50];
	char new4n4[50];
	char nesm4w4[50];
	char newall[50];
	char nemall[50];
	char neb1[50];
	char neb25[50];
	char nebranch[50];

	sprintf(nemsg,"enmsg_%s",id);
	sprintf(newtot,"ewtot_%s",id);
	sprintf(nentot,"entot_%s",id);
	sprintf(new1,"ew1_%s",id);
	sprintf(nen1,"en1_%s",id);
	sprintf(new1nn1,"ew1nn1_%s",id);
	sprintf(nesmtot,"esmtot_%s",id);
	sprintf(nesm1,"esm1_%s",id);
	sprintf(neal,"eal_%s",id);
	sprintf(nep,"ep_%s",id);
	sprintf(new4,"ew4_%s",id);
	sprintf(nen4,"en4_%s",id);
	sprintf(new4n4,"ew4n4_%s",id);
	sprintf(nesm4w4,"esm4w4_%s",id);
	sprintf(newall,"ewall_%s",id);
	sprintf(nemall,"emall_%s",id);
	sprintf(neb1,"eb1_%s",id);
	sprintf(neb25,"eb25_%s",id);
	sprintf(nebranch,"enbranch_%s",id);

	/*
	TH1F *hemsg=hm->Book<TH1F>(nemsg,nemsg,11,-0.5,10.5);
	TH1F *hewtot=hm->Book<TH1F>(newtot,newtot,100,0,2000);
	TH1F *hentot=hm->Book<TH1F>(nentot,nentot,100,0,300);
	TH1F *hew1=hm->Book<TH1F>(new1,new1,100,0,2000);
	TH1F *hen1=hm->Book<TH1F>(nen1,nen1,100,0,300);
	TH1F *hew1nn1=hm->Book<TH1F>(new1nn1,new1nn1,100,0,50);
	TH1F *hesmtot=hm->Book<TH1F>(nesmtot,nesmtot,100,0,5000);
	TH1F *hesm1=hm->Book<TH1F>(nesm1,nesm1,100,0,5000);
	TH1F *heal=hm->Book<TH1F>(neal,neal,51,-0.01,1.01);
	TH1F *hep=hm->Book<TH1F>(nep,nep,50,0,100);
	TH1F *hew4=hm->Book<TH1F>(new4,new4,50,0,1000);
	TH1F *hen4=hm->Book<TH1F>(nen4,nen4,50,0,100);
	TH1F *hew4n4=hm->Book<TH1F>(new4n4,new4n4,100,0,50);
	TH1F *hesm4w4=hm->Book<TH1F>(nesm4w4,nesm4w4,100,0,5);
	TH1F *hewall=hm->Book<TH1F>(newall,newall,100,0,100);
	TH1F *hemall=hm->Book<TH1F>(nemall,nemall,100,0,100);
	TH1F *heb1=hm->Book<TH1F>(neb1,neb1,10,0,1);
	TH1F *heb25=hm->Book<TH1F>(neb25,neb25,10,0,1);
	TH1F *hebranch=hm->Book<TH1F>(nebranch,nebranch,21,-0.5,20.5);
	*/
	hm->Book<TH1F>(nemsg,nemsg,11,-0.5,10.5);
	hm->Book<TH1F>(newtot,newtot,100,0,2000);
	hm->Book<TH1F>(nentot,nentot,100,0,300);
	hm->Book<TH1F>(new1,new1,100,0,2000);
	hm->Book<TH1F>(nen1,nen1,100,0,300);
	hm->Book<TH1F>(new1nn1,new1nn1,100,0,50);
	hm->Book<TH1F>(nesmtot,nesmtot,100,0,5000);
	hm->Book<TH1F>(nesm1,nesm1,100,0,5000);
	hm->Book<TH1F>(neal,neal,51,-0.01,1.01);
	hm->Book<TH1F>(nep,nep,50,0,100);
	hm->Book<TH1F>(new4,new4,50,0,1000);
	hm->Book<TH1F>(nen4,nen4,50,0,100);
	hm->Book<TH1F>(new4n4,new4n4,100,0,50);
	hm->Book<TH1F>(nesm4w4,nesm4w4,100,0,5);
	hm->Book<TH1F>(newall,newall,100,0,100);
	hm->Book<TH1F>(nemall,nemall,100,0,100);
	hm->Book<TH1F>(neb1,neb1,10,0,1);
	hm->Book<TH1F>(neb25,neb25,10,0,1);
	hm->Book<TH1F>(nebranch,nebranch,21,-0.5,20.5);


	char nomsg[50];
	char nowtot[50];
	char nontot[50];
	char now1[50];
	char non1[50];
	char now1nn1[50];
	char nosmtot[50];
	char nosm1[50];
	char noal[50];
	char nop[50];
	char now4[50];
	char non4[50];
	char now4n4[50];
	char nosm4w4[50];
	char nowall[50];
	char nomall[50];
	char nob1[50];
	char nob25[50];
	char nobranch[50];

	sprintf(nomsg,"onmsg_%s",id);
	sprintf(nowtot,"owtot_%s",id);
	sprintf(nontot,"ontot_%s",id);
	sprintf(now1,"ow1_%s",id);
	sprintf(non1,"on1_%s",id);
	sprintf(now1nn1,"ow1nn1_%s",id);
	sprintf(nosmtot,"osmtot_%s",id);
	sprintf(nosm1,"osm1_%s",id);
	sprintf(noal,"oal_%s",id);
	sprintf(nop,"op_%s",id);
	sprintf(now4,"ow4_%s",id);
	sprintf(non4,"on4_%s",id);
	sprintf(now4n4,"ow4n4_%s",id);
	sprintf(nosm4w4,"osm4w4_%s",id);
	sprintf(nowall,"owall_%s",id);
	sprintf(nomall,"omall_%s",id);
	sprintf(nob1,"ob1_%s",id);
	sprintf(nob25,"ob25_%s",id);
	sprintf(nobranch,"onbranch_%s",id);
	/*
	TH1F *homsg=hm->Book<TH1F>(nomsg,nomsg,11,-0.5,10.5);
	TH1F *howtot=hm->Book<TH1F>(nowtot,nowtot,100,0,2000);
	TH1F *hontot=hm->Book<TH1F>(nontot,nontot,100,0,300);
	TH1F *how1=hm->Book<TH1F>(now1,now1,100,0,2000);
	TH1F *hon1=hm->Book<TH1F>(non1,non1,100,0,300);
	TH1F *how1nn1=hm->Book<TH1F>(now1nn1,now1nn1,100,0,50);
	TH1F *hosmtot=hm->Book<TH1F>(nosmtot,nosmtot,100,0,5000);
	TH1F *hosm1=hm->Book<TH1F>(nosm1,nosm1,100,0,5000);
	TH1F *hoal=hm->Book<TH1F>(noal,noal,51,-0.01,1.01);
	TH1F *hop=hm->Book<TH1F>(nop,nop,50,0,100);
	TH1F *how4=hm->Book<TH1F>(now4,now4,50,0,1000);
	TH1F *hon4=hm->Book<TH1F>(non4,non4,50,0,100);
	TH1F *how4n4=hm->Book<TH1F>(now4n4,now4n4,100,0,50);
	TH1F *hosm4w4=hm->Book<TH1F>(nosm4w4,nosm4w4,100,0,5);
	TH1F *howall=hm->Book<TH1F>(nowall,nowall,100,0,100);
	TH1F *homall=hm->Book<TH1F>(nomall,nomall,100,0,100);
	TH1F *hob1=hm->Book<TH1F>(nob1,nob1,10,0,1);
	TH1F *hob25=hm->Book<TH1F>(nob25,nob25,10,0,1);
	TH1F *hobranch=hm->Book<TH1F>(nobranch,nobranch,21,-0.5,20.5);
	*/	
	hm->Book<TH1F>(nomsg,nomsg,11,-0.5,10.5);
	hm->Book<TH1F>(nowtot,nowtot,100,0,2000);
	hm->Book<TH1F>(nontot,nontot,100,0,300);
	hm->Book<TH1F>(now1,now1,100,0,2000);
	hm->Book<TH1F>(non1,non1,100,0,300);
	hm->Book<TH1F>(now1nn1,now1nn1,100,0,50);
	hm->Book<TH1F>(nosmtot,nosmtot,100,0,5000);
	hm->Book<TH1F>(nosm1,nosm1,100,0,5000);
	hm->Book<TH1F>(noal,noal,51,-0.01,1.01);
	hm->Book<TH1F>(nop,nop,50,0,100);
	hm->Book<TH1F>(now4,now4,50,0,1000);
	hm->Book<TH1F>(non4,non4,50,0,100);
	hm->Book<TH1F>(now4n4,now4n4,100,0,50);
	hm->Book<TH1F>(nosm4w4,nosm4w4,100,0,5);
	hm->Book<TH1F>(nowall,nowall,100,0,100);
	hm->Book<TH1F>(nomall,nomall,100,0,100);
	hm->Book<TH1F>(nob1,nob1,10,0,1);
	hm->Book<TH1F>(nob25,nob25,10,0,1);
	hm->Book<TH1F>(nobranch,nobranch,21,-0.5,20.5);

      }
    }
  }
}


JobCResult CompareMST::Ana(const MomNavigator* mom)
{
   //get all NueRecords from mom 
   //may have more than one per go since mom reads in a snarl's worth of data
   //so, this is a little more complicated than just asking for a NueRecord
   TObject *obj=0;
   TIter objiter = mom->FragmentIter();
   while((obj=objiter.Next())){
      NueRecord *nr = static_cast<NueRecord *>(obj);
      if(nr){
	 MSG("CompareMST",Msg::kDebug)<<"Found a NueRecord in MOM"<<endl;
      }
      else{
	 MSG("CompareMST",Msg::kError)<<"Didn't find a NueRecord in MOM"<<endl;
	 continue;
      }
      MSG("CompareMST",Msg::kDebug)<<"Found a NueRecord "<<nr<<endl;
      
      //need at least 1 sucessfully reconstructed event
      if(nr->GetHeader().GetEventNo()<0){
	continue;
      }
      //only look at events for which mst actually gets calculated
      if(nr->GetHeader().GetTrackLength()>25){
	continue;
      }
      //only look at contained events (vertex contained, end doesn't go past 
      //calorimeter or sm gap
      if(nr->GetHeader().GetVldContext().GetDetector()==Detector::kFar){
	if(sqrt(pow(nr->srevent.vtxX,2)+pow(nr->srevent.vtxY,2))>3.87){
	  continue;
	}
	if(nr->srevent.vtxZ<1||(nr->srevent.vtxZ>14&&nr->srevent.vtxZ<17)||
	   nr->srevent.vtxZ>20){
	  continue;
	}
	if((nr->srevent.begPlane<SM1end&&nr->srevent.endPlane>SM1end-2)||
	   (nr->srevent.begPlane>SM1end&&nr->srevent.endPlane>SM2end-2)){
	  continue;
	}	  
      }
      else if(nr->GetHeader().GetVldContext().GetDetector()==Detector::kNear){
	if(sqrt(pow(nr->srevent.vtxX-1.5,2)+pow(nr->srevent.vtxY,2))>1.0){
	  continue;
	}
	if(nr->srevent.vtxZ<4||nr->srevent.vtxZ>6.5){
	  continue;
	}
	if(nr->srevent.endPlane>CALend-2){
	  continue;
	}	  
      }
      //only look at events where mst output makes sense, nn1>3
      if(nr->mstvars.enn1<=3||nr->mstvars.onn1<=3){
	continue;
      }
      
      //make id string
      Detector::Detector_t d = nr->GetHeader().GetVldContext().GetDetector();
      char evdet[2];
      if(d==Detector::kFar){
	sprintf(evdet,"f");
      }
      else if(d==Detector::kNear){
	sprintf(evdet,"n");
      }
      else{
	sprintf(evdet,"u");
      }
      
      char evbgc[6];
      if(nr->mctrue.interactionType==1){
	if(abs(nr->mctrue.nuFlavor)==12){
	  sprintf(evbgc,"nue");
	}
	else if(abs(nr->mctrue.nuFlavor)==14){
	  sprintf(evbgc,"numu");
	}
	else if(abs(nr->mctrue.nuFlavor)==16){
	  sprintf(evbgc,"nutau");
	}
	else{
	  sprintf(evbgc,"u");
	}
      }
      else{
	sprintf(evbgc,"nc");
      }

      char evid[20];
      sprintf(evid,"%s_%s_%d",evdet,evbgc,nr->mctrue.resonanceCode);

      //make histname strings
      char nemsg[50];
      char newtot[50];
      char nentot[50];
      char new1[50];
      char nen1[50];
      char new1nn1[50];
      char nesmtot[50];
      char nesm1[50];
      char neal[50];
      char nep[50];
      char new4[50];
      char nen4[50];
      char new4n4[50];
      char nesm4w4[50];
      char newall[50];
      char nemall[50];
      char neb1[50];
      char neb25[50];
      char nebranch[50];
      
      sprintf(nemsg,"enmsg_%s",evid);
      sprintf(newtot,"ewtot_%s",evid);
      sprintf(nentot,"entot_%s",evid);
      sprintf(new1,"ew1_%s",evid);
      sprintf(nen1,"en1_%s",evid);
      sprintf(new1nn1,"ew1nn1_%s",evid);
      sprintf(nesmtot,"esmtot_%s",evid);
      sprintf(nesm1,"esm1_%s",evid);
      sprintf(neal,"eal_%s",evid);
      sprintf(nep,"ep_%s",evid);
      sprintf(new4,"ew4_%s",evid);
      sprintf(nen4,"en4_%s",evid);
      sprintf(new4n4,"ew4n4_%s",evid);
      sprintf(nesm4w4,"esm4w4_%s",evid);
      sprintf(newall,"ewall_%s",evid);
      sprintf(nemall,"emall_%s",evid);
      sprintf(neb1,"eb1_%s",evid);
      sprintf(neb25,"eb25_%s",evid);
      sprintf(nebranch,"enbranch_%s",evid);

      char nomsg[50];
      char nowtot[50];
      char nontot[50];
      char now1[50];
      char non1[50];
      char now1nn1[50];
      char nosmtot[50];
      char nosm1[50];
      char noal[50];
      char nop[50];
      char now4[50];
      char non4[50];
      char now4n4[50];
      char nosm4w4[50];
      char nowall[50];
      char nomall[50];
      char nob1[50];
      char nob25[50];
      char nobranch[50];

      sprintf(nomsg,"onmsg_%s",evid);
      sprintf(nowtot,"owtot_%s",evid);
      sprintf(nontot,"ontot_%s",evid);
      sprintf(now1,"ow1_%s",evid);
      sprintf(non1,"on1_%s",evid);
      sprintf(now1nn1,"ow1nn1_%s",evid);
      sprintf(nosmtot,"osmtot_%s",evid);
      sprintf(nosm1,"osm1_%s",evid);
      sprintf(noal,"oal_%s",evid);
      sprintf(nop,"op_%s",evid);
      sprintf(now4,"ow4_%s",evid);
      sprintf(non4,"on4_%s",evid);
      sprintf(now4n4,"ow4n4_%s",evid);
      sprintf(nosm4w4,"osm4w4_%s",evid);
      sprintf(nowall,"owall_%s",evid);
      sprintf(nomall,"omall_%s",evid);
      sprintf(nob1,"ob1_%s",evid);
      sprintf(nob25,"ob25_%s",evid);
      sprintf(nobranch,"onbranch_%s",evid);

      //fill histograms      
      static HistMan *hm = new HistMan("mstcomp");
      //      cout<<"Printing hm "<<endl;
      //      hm->BaseFolder().ls();
      //      cout<<"******************************"<<endl;

      hm->Fill1d(nemsg,nr->mstvars.enmsg);
      hm->Fill1d(newtot,nr->mstvars.ewtot);
      hm->Fill1d(nentot,nr->mstvars.enntot);
      hm->Fill1d(new1,nr->mstvars.ew1);
      hm->Fill1d(nen1,nr->mstvars.enn1);
      if(nr->mstvars.enn1!=0){
	hm->Fill1d(new1nn1,nr->mstvars.ew1/nr->mstvars.enn1);
      }
      hm->Fill1d(nesmtot,nr->mstvars.esmtot);
      hm->Fill1d(nesm1,nr->mstvars.esm1);
      hm->Fill1d(neal,nr->mstvars.ealpha);
      hm->Fill1d(nep,nr->mstvars.eeprob);
      hm->Fill1d(new4,nr->mstvars.e4w);
      hm->Fill1d(nen4,nr->mstvars.e4nn);
      if(nr->mstvars.e4nn!=0){
	hm->Fill1d(new4n4,nr->mstvars.e4w/nr->mstvars.e4nn);
      }
      if(nr->mstvars.e4w!=0){
	hm->Fill1d(nesm4w4,nr->mstvars.e4sm/nr->mstvars.e4w);
      }
      for(int i=0;i<nr->mstvars.enn1;i++){
	hm->Fill1d(newall,nr->mstvars.eallw1[i]);
	hm->Fill1d(nemall,nr->mstvars.eallm1[i]);
      }
      hm->Fill1d(neb1,nr->mstvars.eb1);
      hm->Fill1d(neb25,nr->mstvars.eb25);
      hm->Fill1d(nebranch,nr->mstvars.enbranch);

      hm->Fill1d(nomsg,nr->mstvars.onmsg);
      hm->Fill1d(nowtot,nr->mstvars.owtot);
      hm->Fill1d(nontot,nr->mstvars.onntot);
      hm->Fill1d(now1,nr->mstvars.ow1);
      hm->Fill1d(non1,nr->mstvars.onn1);
      if(nr->mstvars.onn1!=0){
	hm->Fill1d(now1nn1,nr->mstvars.ow1/nr->mstvars.onn1);
      }
      hm->Fill1d(nosmtot,nr->mstvars.osmtot);
      hm->Fill1d(nosm1,nr->mstvars.osm1);
      hm->Fill1d(noal,nr->mstvars.oalpha);
      hm->Fill1d(nop,nr->mstvars.oeprob);
      hm->Fill1d(now4,nr->mstvars.o4w);
      hm->Fill1d(non4,nr->mstvars.o4nn);
      if(nr->mstvars.o4nn!=0){
	hm->Fill1d(now4n4,nr->mstvars.o4w/nr->mstvars.o4nn);
      }
      if(nr->mstvars.o4w!=0){
	hm->Fill1d(nosm4w4,nr->mstvars.o4sm/nr->mstvars.o4w);
      }
      for(int i=0;i<nr->mstvars.onn1;i++){
	hm->Fill1d(nowall,nr->mstvars.oallw1[i]);
	hm->Fill1d(nomall,nr->mstvars.oallm1[i]);
      }
      hm->Fill1d(nob1,nr->mstvars.ob1);
      hm->Fill1d(nob25,nr->mstvars.ob25);
      hm->Fill1d(nobranch,nr->mstvars.onbranch);



   }
   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

////////////////////////////////////////////////////////////////////////
void CompareMST::EndJob()
{}
