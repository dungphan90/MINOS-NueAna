#include "TList.h"
#include "TH2D.h"
#include "THStack.h"
#include "TVector3.h"
#include <deque>
#include "TCanvas.h"
#include "TNtuple.h"
#include "TPad.h"
#include "TEllipse.h"
#include "TText.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TButton.h"
#include "TLine.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TPolyLine.h"
#include "TDirectory.h"
#include "TList.h"

#include "NueAna/ParticlePID/ParticleFinder/Display/ParticleDisplay.h"
#include "NueAna/Display/SelectPad.h"
#include "NueAna/NueRecord.h"
//#include "NueAna/NuePID.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"

#include "CandSubShowerSR/ClusterType.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "CandNtupleSR/NtpSRCluster.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "CandNtupleSR/NtpSRShieldStrip.h"
#include "CandNtupleSR/NtpSRShowerPulseHeight.h"
#include "MCNtuple/NtpMCRecord.h"
#include "MCNtuple/NtpMCTruth.h"
#include "MCNtuple/NtpMCStdHep.h"
#include "TruthHelperNtuple/NtpTHRecord.h"
#include "TruthHelperNtuple/NtpTHEvent.h"
#include "StandardNtuple/NtpStRecord.h"
#include "Conventions/Detector.h"
#include "Conventions/SimFlag.h"

#include "Midad/Base/Mint.h"
#include "Midad/Base/PageDisplay.h"
#include "Midad/Base/CanvasSignals.h"
#include "Midad/Base/SteelOutline.h"
#include "Midad/Base/TimeHist.h"
#include "Midad/Gui/GuiButton.h"
#include "Midad/Gui/GuiBox.h"
#include "Midad/Gui/GuiTextEntry.h"
#include "Midad/Gui/GuiTextView.h"
#include "Midad/Gui/GuiTab.h"
#include "Midad/Gui/GuiLabel.h"
#include "Midad/Gui/GuiMainWindow.h"
#include <sigc++/sigc++.h>
#include <sigc++/class_slot.h>

#include "JobControl/JobC.h"
#include "MinosObjectMap/MomNavigator.h"

#include "DataUtil/GetDetectorType.h"
#include "DataUtil/CDL2STL.h"
#include "DataUtil/GetRunSnarlEvent.h"
#include "DataUtil/PlaneOutline.h"

#include "Calibrator/CalMIPCalibration.h"
#include "Record/RecRecordImp.h"
#include "Record/RecCandHeader.h"

#include "Plex/PlexStripEndId.h"

#include "RecoBase/PropagationVelocity.h"

#include "UgliGeometry/UgliStripHandle.h"
#include "UgliGeometry/UgliGeomHandle.h"

#include "CandFitShowerEM/FitterEM.h"


#include "NueAna/NueAnaTools/NueConvention.h"
#include "NueAna/ParticlePID/ParticleAna/ParticlesAna.h"

#include <fstream>
#include <cmath>

using namespace DataUtil;
using namespace SigC;

// Boiler plate for using the Message Service.
#include "MessageService/MsgService.h"
CVSID("$Id: ParticleDisplay.cxx,v 1.3 2009/06/23 22:42:25 scavan Exp $");

// Boiler plate needed for us to be a Job Module.
#include "JobControl/JobCModuleRegistry.h"
JOBMODULE(ParticleDisplay,"ParticleDisplay","Example of adding a Midad display in a Job Module\n");

#include "DatabaseInterface/DbiResultPtr.tpl"

const string evtpcode[] = {"N/A", "mu","e","NC","mu/NC?","e/NC?","???"};
const string topocode[] = {"N/A", "QE","RES","DIS","???"};


ParticleDisplay::ParticleDisplay()
  : clickbutton(0)
  , fSlice(-1)
  , fNumSlices(-1)
  , fEvent(-1)
  , fEvent_old(-1)
  , fNumEvents(-1)
  , fSnarl(-1)
  , RunNo_old(-1)
  , SubRunNo_old(-1)
  , imctruth(0)
  , ifixmcinfo(0)
  , kDPlaneCut(-1)
  , kLoPhNStripCut(-1)
  , kLoPhNPlaneCut(-1)
  , kPhStripCut(-1)
  , kPhPlaneCut(-1)
  , kCPhPlaneCut(-1.)
  , kPIDCut(-1)
  , kScanMode(0)
  , kTestMode(0)
  , kHideRunSnarl(0)
  , kDrawClu(0)
  , kIntReco(0)
  , foundmeu(false)
  , SIGCORRMEU(1.)
  , fHBox(0)
  , fVBox1(0)
  , fVBox2(0)
  , fVBox3(0)
  , fHBox1(0)
  , fHBox2(0)
  , fHBox3(0)
  , fHBox4(0)
  , fHBox5(0)
  , ievtp(0)
  , itopo(0)



  , selecevtp(0)
  , selectopo(0)
  , fRel(ReleaseType::kUnknown)
{	 
 
  for (int i = 0; i<7; i++){ iEvtp[i] = 0;}
  for (int i = 0; i<5; i++){ iTopo[i] = 0;}

  showStpTrueMu = false;
  showStpTrueShw = false;
  showStpScaled = false;
  showStpRecoTrk = false;
  showStpAll = true;


  showOrigEvent = true;
  showNewEvent = true;
}

ParticleDisplay::~ParticleDisplay()
{

}

void ParticleDisplay::BeginRun()
{
  cout<<"In ParticleDisplay beginrun()"<<endl;
  if (gMint ) this->BuildDisplay(); //must check to see if its already built!
}

JobCResult ParticleDisplay::Ana(const MomNavigator *mom)
{
  if (gMint) {                // paranoia
    if (&gMint->GetJobC().Mom != mom) {
      MSG("ParticleDisplay",Msg::kError) << "Module's mom and JobC's mom differ: "
				     << (void*)&gMint->GetJobC().Mom << " != "
				     << (void*) mom << endl;
    }
  }
  //get data from mom
  foundPOH=false;
  foundST=false;
 
  st = 0;
  poh = 0;



  VldContext vc;  
  
  hft = mom->GetFragmentList("ParticleObjectHolder");
  
  std::vector<TObject* > precord  = mom->GetFragmentList("PRecord");
  
  cout <<" there are "<<precord.size()<< " ana records\n";
  
  cout <<"particles in each: ";
  for(unsigned int i=0;i<precord.size();i++)
  	cout << ((PRecord*)precord[i])->particles.ntot;
  cout <<"\n";
  
  if(hft.size()>0)poh=(ParticleObjectHolder *)hft[0];
  
  st =(NtpStRecord*)mom->GetFragment("NtpStRecord","Primary");
  
  
  if(poh){
    foundPOH=true;
    vc=st->GetHeader().GetVldContext();

    string relName = st->GetTitle();
    string reco = relName.substr(0,relName.find_first_of("("));
    if(reco == "CEDAR") fRel = ReleaseType::kCedar;
    else fRel = ReleaseType::kBirch;
    if(vc.GetSimFlag() == SimFlag::kMC){
 //    NtpMCGenInfo* genInfo = &(st->mchdr.geninfo);
   //   if(strstr(genInfo->codename.c_str(), "daikon") == 0)
//	fRel += ReleaseType::kCarrot;
//      if(strstr(genInfo->codename.c_str(), "daikon") != 0)
//	fRel += ReleaseType::kDaikon;
    }
    if(vc.GetSimFlag() != SimFlag::kMC) fRel += ReleaseType::kData;
  }
  
  
    if(!foundPOH) {
    MSG("ParticleDisplay",Msg::kError)<<"Got Nothing to display"<<endl;
      return JobCResult::kFailed;
   }


   
  fEvent = 0;
  this->UpdateDisplay(1,1);
  
  return JobCResult::kAOK;
}

void ParticleDisplay::BuildDisplay()
{
  MSG("ParticleDisplay",Msg::kDebug)<<"In BuildDisplay"<<endl;


  //////////////////////////////////////////////
  //do some resets...
  //
  int run=0, subrun = 0, snarl = 0, evt = 0;
  DataUtil::GetRunSnarlEvent(&(gMint->GetJobC().Mom),run,snarl,evt);
  RecRecordImp<RecCandHeader> *rr = 
    dynamic_cast<RecRecordImp<RecCandHeader>*>
    ((gMint->GetJobC().Mom).GetFragment("RecRecordImp<RecCandHeader>"));
  if ( rr ) {
    subrun   = rr->GetHeader().GetSubRun();
    fDetectorType = rr->GetHeader().GetVldContext().GetDetector();
    fSimFlag = rr->GetHeader().GetVldContext().GetSimFlag();
  }
  RunNo = run;
  SubRunNo = subrun;

  fSnarl = -1; //reset snarl no at the beginning of each run
  if (RunNo_old==-1&&SubRunNo_old==-1){
    RunNo_old = RunNo;
    SubRunNo_old = SubRunNo;
  }
  //////////////////////////////////////////////////////////////////

  //generate main display

  const int width = 1000, height = 680;

  PageDisplay* pd = gMint->GetDisplay();
  //get Button Box on the left
  GuiBox *fBBox = pd->GetButtonBox();
  // No pre-existing display, so make one to our size
  if (!pd) {
    MSG("ParticleDisplay",Msg::kDebug)<<"No display, making one"<<endl;
    pd = gMint->SpawnDisplay(width,height);
  }

  //////////////////////////////////////////////////////////////////////

  //add buttons in the Button Box
  //"Next Event" "Prev Event"
  fNextEventbut = pd->AddButton("Next Event ");
  fNextEventbut->clicked.connect(slot_class(*this,&ParticleDisplay::NextEvent));
  fNextEventbut = pd->AddButton("Prev Event ");
  fNextEventbut->clicked.connect(slot_class(*this,&ParticleDisplay::PrevEvent));

  //Event No
  fEventNo = manage(new GuiLabel(*fBBox," "));
  fEventNo->SetLayoutHints(kLHintsExpandX);
  fBBox->Add(*fEventNo);
  fEventEntry = pd->AddEntry(" ");
  fEventEntry->activated.connect(slot_class(*this,&ParticleDisplay::GotoEvent));
  
  //"NextSelEvt" "PreSelEvt"
  fNextSelEvt = pd->AddButton("NextSelEvt");
  fNextSelEvt->clicked.connect(slot_class(*this,&ParticleDisplay::NextSelEvt));
  fPrevSelEvt = pd->AddButton("PrevSelEvt");
  fPrevSelEvt->clicked.connect(slot_class(*this,&ParticleDisplay::PrevSelEvt));

  if(!kTestMode){
    fMCTruth = pd->AddButton("MC Truth");
    fMCTruth->clicked.connect(slot_class(*this,&ParticleDisplay::showmctruth));

    fFixMCInfo = pd->AddButton("Fix MC Info");
    fFixMCInfo->clicked.connect(slot_class(*this,&ParticleDisplay::fixmcinfo));
  }


  ////////////////////////////////////////////////////

  //construct frame for voting and I/O facilities

  //main box fHBox
  if(!fHBox){
    fHBox = manage(new GuiBox(*pd,kHorizontalFrame));
    fHBox->SetHeight(140);
    pd->SetWidth(pd->GetWidth());
    fHBox->SetWidth(pd->GetWidth());
    fHBox->SetLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX);
    pd->Add(*fHBox);
  } 
  


  //reco info
  fHBox2 = manage(new GuiBox(*fHBox,kHorizontalFrame));
  fHBox->Add(*fHBox2);

  fRecoInfo = manage(new GuiTextView(*fHBox2,1,1));
  fHBox2->Add(*fRecoInfo);
  fRecoInfo->AddLine("Reco Information");



  ///////////////////////////////////////////////////////////////    

  //add user-defined canvas
  
  
    
  CanvasSignals* cs0 = 0;
  CanvasSignals* cs1 = 0;
  CanvasSignals* cs2 = 0;
  //CanvasSignals* cs3 = 0;
  cs0 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));
  cs1 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));
  cs2 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));  
  
/*  hitview.BuildDisplay(&cs0->GetCanvas());
      gPad->Modified();
    cs0->GetCanvas().Update();
*/
/*  fitview.BuildDisplay(&cs0->GetCanvas());
      gPad->Modified();
    cs0->GetCanvas().Update();
*/

  chainview.BuildDisplay(&cs1->GetCanvas());
      gPad->Modified();
    cs1->GetCanvas().Update();
    
  view3d.BuildDisplay(&cs2->GetCanvas());
        gPad->Modified();
    cs2->GetCanvas().Update();


    
  houghview.BuildDisplay(&cs0->GetCanvas());
        gPad->Modified();
    cs0->GetCanvas().Update();
  
    fEvent = 0;    
}

/*
static void clear_hist(TH2D* hist)
{
    hist->Reset();
    hist->GetXaxis()->UnZoom();
    hist->GetYaxis()->UnZoom();
    hist->SetMaximum(-1111);
    hist->SetMinimum(-1111);
}
*/
void ParticleDisplay::UpdateDisplay(bool /*foundnr*/, bool /*foundpid*/)
{
  //MSG("ParticleDisplay",Msg::kDebug)
  cout<<"In UpdateDisplay"<<endl;
  
  
   
  int run=0, subrun = 0, snarl = 0, evt = 0;
  DataUtil::GetRunSnarlEvent(&(gMint->GetJobC().Mom),run,snarl,evt);
  RecRecordImp<RecCandHeader> *rr = 
    dynamic_cast<RecRecordImp<RecCandHeader>*>
    ((gMint->GetJobC().Mom).GetFragment("RecRecordImp<RecCandHeader>"));
  if ( rr ) {
    subrun   = rr->GetHeader().GetSubRun();
    fDetectorType = rr->GetHeader().GetVldContext().GetDetector();
    fSimFlag = rr->GetHeader().GetVldContext().GetSimFlag();
  }
  RunNo = run;
  SubRunNo = subrun;







  if (imctruth||ifixmcinfo){
    imctruth = 0;
    fMCTruth->SetDown(false);
   
  }
  if (foundPOH) {
    fNumEvents = hft.size();//st->fHeader.nevent;
  }
 
  if (clickbutton == -2 && fNumEvents){
    fEvent = fNumEvents-1;
    clickbutton = -1;
  }
  else if (clickbutton == -2){
    gMint->Prev();
    return;
  }

  if (!GetEvent(fEvent)){
    if (foundPOH){
      MSG("ParticleDisplay",Msg::kError)<<"Couldn't get event "<<fEvent
					 <<" from Snarl "<<st->GetHeader().GetSnarl()<<endl;
      if (clickbutton == -1) gMint->Prev();
      if (clickbutton ==  1) gMint->Next();
      return;
    }

  }
 
  
  
  if(st)
  cout <<  st->GetHeader().GetSnarl()<<" "<<st->GetHeader().GetEvent()<<endl;
  
 
     fEventEntry->SetText(Form("%d",fEvent));
    fEventNo->SetText(Form("Event: %d(%d)",fEvent,fNumEvents));
 
 
 	double particleE=0;
	std::vector<Particle3D>  particles3d = poh->particles3d;
	for(unsigned int ik=0;ik<particles3d.size();ik++)
	{
		Particle3D * p3 = (Particle3D *)&particles3d[ik];
		if(p3 <=0)continue;
		particleE+=p3->sum_e;
	}




 
	cout <<"#################################\n";
	cout <<"###### RUN " << RunNo <<" SNARL "<<snarl << " EVENT "<<fEvent<<"\n";
	cout <<"###### totE " << poh->event.visenergy << " particleE "<< particleE;
	if(particleE>poh->event.visenergy) cout <<" !!!energy mismatch!";
	cout <<"\n";
	cout <<"###### Vertex U,V,Z "<< poh->event.vtx_u<<", "<<poh->event.vtx_v<<", "<<poh->event.vtx_z<<"\n";
	cout <<"#################################\n";
	



  poh=0;

  if(fEvent>-1 && (int)hft.size()>fEvent)poh=(ParticleObjectHolder*)hft[fEvent];

 // hitview.DrawEvent(poh,st,fEvent);
 // fitview.DrawEvent(poh,st,fEvent);  
  chainview.DrawEvent(poh,st,fEvent);
  view3d.DrawEvent(poh,st,fEvent);
   houghview.DrawEvent(poh,st,fEvent);
  
  return;

 
}

Int_t ParticleDisplay::GetEvent(Int_t /*evt*/){

  // Get Ugli for later.  Bail immediately if fail

/*
  if (foundST) {
    event = SntpHelpers::GetEvent(evt,st);
  }
  else if(foundSR){
    event = SntpHelpers::GetEvent(evt,sr);
  }
*/


  return 1;
}
  
    


void ParticleDisplay::NextEvent(){
  clickbutton = 1;

  if(fNumEvents>0){
    fEvent++;
    if(fEvent>=fNumEvents){
      fEvent = 0;
      gMint->Next();
      return;
    }
    this->UpdateDisplay(1, 1);
  }
  else {
    gMint->Next();
    return;
  }
}


void ParticleDisplay::PrevEvent(){

  if(fNumEvents>0){
    fEvent--;
    if(fEvent<0){
      clickbutton = -2;
      gMint->Prev();
      return;
    }
    this->UpdateDisplay(1, 1);
  }
  else {
    gMint->Prev();
    return;
  }
}

void ParticleDisplay::GotoEvent()
{
  if (!gMint) return;

  string eventno = fEventEntry->GetText();
  if (atoi(eventno.c_str())<fNumEvents && atoi(eventno.c_str())>=0){
    fEvent = atoi(eventno.c_str());
    this->UpdateDisplay(1, 1);
  }
  else {
    return;
  }
}

void ParticleDisplay::showmctruth()
{
  
//  if (imctruth) {
//    fMCTruth->SetDown(true);
//    return;
//  }
  if (fSimFlag!=SimFlag::kMC) return;
  if (ifixmcinfo) return;
  if (!imctruth){//not pressed yet
    imctruth = 1;
    fMCTruth->SetDown(true);

  }
  else {
    imctruth = 0;
    fMCTruth->SetDown(false);

  }
}
void ParticleDisplay::fixmcinfo()
{
  if (fSimFlag!=SimFlag::kMC) return;
  if (!ifixmcinfo){
    fFixMCInfo->SetDown(true);
    ifixmcinfo = 1;
    if (!imctruth){

    }
    else {
      fMCTruth->SetDown(false);
      imctruth = 0;
    }
  }
  else {
    fFixMCInfo->SetDown(false);
    ifixmcinfo = 0;

  }
}
  
const Registry& ParticleDisplay::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
   MSG("ParticleDisplay",Msg::kDebug)<<"In ParticleDisplay::DefaultConfig"<<endl;

  static Registry r = fCut.DefaultConfig();

  // Set name of config
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());

  // Set values in configuration
  r.UnLockValues();
  r.Set("DPlaneCut",-1);
  r.Set("LoPhNStripCut",-1);
  r.Set("LoPhNPlaneCut",-1);
  r.Set("PhStripCut",-1);
  r.Set("PhPlaneCut",-1);
  r.Set("ContPhPlaneCut",-1);

  r.Set("PIDCut",-1);
  r.Set("ScanMode",0);
  r.Set("TestMode",0);
  r.Set("HideRunSnarl",0);
  r.Set("DrawClu",0);
  r.Set("IntReco",0);
  r.LockValues();

  return r;
}

//......................................................................

void ParticleDisplay::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
  MSG("ParticleDisplay",Msg::kDebug)<<"In ParticleDisplay::Config"<<endl;

  fCut.Config(r);

  int imps;
  if(r.Get("DPlaneCut",imps)) { kDPlaneCut=imps;}
  if(r.Get("LoPhNStripCut",imps)) { kLoPhNStripCut=imps;}
  if(r.Get("LoPhNPlaneCut",imps)) { kLoPhNPlaneCut=imps;}

  if(r.Get("ScanMode",imps)) {kScanMode=imps;}
  if(r.Get("TestMode",imps)) {kTestMode=imps;}
  if(r.Get("HideRunSnarl",imps)) {kHideRunSnarl=imps;}
  if(r.Get("DrawClu",imps)) {kDrawClu=imps;}
  if(r.Get("IntReco",imps)) {kIntReco=imps;}
  double fmps;
  if(r.Get("PhStripCut",fmps)) { kPhStripCut=fmps;}
  if(r.Get("PhPlaneCut",fmps)) { kPhPlaneCut=fmps;}
  if(r.Get("ContPhPlaneCut",fmps)) { kCPhPlaneCut=fmps;}

  if(r.Get("PIDCut",fmps)) { kPIDCut=fmps;}
}

void ParticleDisplay::NextSelEvt(){
  if (snarlno.size()){
    /*
      cout<<*runnoitr<<" "<<*subrunnoitr<<" "
      <<fSnarl<<" "<<fEvent
      <<" "<<*snarlnoitr<<" "<<*eventnoitr<<endl;

    if(runnoitr<runno.end()&&runnoitr>=runno.begin()&&
       subrunnoitr<subrunno.end()&&subrunnoitr>=subrunno.begin()&&
       snarlnoitr<snarlno.end()&&snarlnoitr>=snarlno.begin()
       &&(fSnarl < *snarlnoitr || fSnarl == *snarlnoitr && 
       fEvent < *eventnoitr)){
    if(fSnarl == *snarlnoitr && fEvent < *eventnoitr){
      if ((!selecevtp||iEvtp[*typeitr])&&(!selectopo||iTopo[*topoitr])){
	fEvent = *eventnoitr;
	ievtp = *typeitr;
	itopo = *topoitr;
	//gMint->GoTo(*runnoitr,*snarlnoitr);
	if (*typeitr>=0&&*typeitr<7&&*topoitr>=0&&*topoitr<5){
	  string tmp = evtpcode[*typeitr] + " | " + topocode[*topoitr];
	  fTYPETOPO->SetText(tmp.c_str());
	}
	fComment->SetText((*commitr).c_str());
      }
    }
    else 
    */
    if (snarlnoitr != snarlno.end() && snarlnoitr != snarlno.end() - 1){    
      runnoitr++;
      subrunnoitr++;
      snarlnoitr++;
      eventnoitr++;
      typeitr++;
      topoitr++;
      commitr++;
      if ((!selecevtp||iEvtp[*typeitr])&&(!selectopo||iTopo[*topoitr])){
	fEvent = *eventnoitr;
	ievtp = *typeitr;
	itopo = *topoitr;
	gMint->GoTo(*runnoitr,*snarlnoitr);
	if (*typeitr>=0&&*typeitr<7&&*topoitr>=0&&*topoitr<5){
	  string tmp = evtpcode[*typeitr] + " | " + topocode[*topoitr];
	  fTYPETOPO->SetText(tmp.c_str());
	}
	fComment->SetText((*commitr).c_str());
      }
      else {this->NextSelEvt();}
    }
  }
}

void ParticleDisplay::PrevSelEvt(){
  if (snarlno.size()){
    /*
    cout<<*runnoitr<<" "<<*subrunnoitr<<" "
	<<fSnarl<<" "<<fEvent<<" "
	<<*snarlnoitr<<" "<<*eventnoitr<<endl;
    if(runnoitr<runno.end()&&runnoitr>=runno.begin()&&
       subrunnoitr<subrunno.end()&&subrunnoitr>=subrunno.begin()&&
       snarlnoitr<snarlno.end()&&snarlnoitr>=snarlno.begin()&&
       (fSnarl > *snarlnoitr || fSnarl == *snarlnoitr && 
	fEvent > *eventnoitr)){
      if ((!selecevtp||iEvtp[*typeitr])&&(!selectopo||iTopo[*topoitr])){
	fEvent = *eventnoitr;
	ievtp = *typeitr;
	itopo = *topoitr;
	gMint->GoTo(*runnoitr,*snarlnoitr);
	if (*typeitr>=0&&*typeitr<7&&*topoitr>=0&&*topoitr<5){
	  string tmp = evtpcode[*typeitr] + " | " + topocode[*topoitr];
	  fTYPETOPO->SetText(tmp.c_str());
	}
	fComment->SetText((*commitr).c_str());
      }
    }
    else */
    if (snarlnoitr != snarlno.begin() && snarlnoitr != snarlno.begin() -1){
      runnoitr--;
      subrunnoitr--;
      snarlnoitr--;
      eventnoitr--;
      typeitr--;
      topoitr--;
      commitr--;
      if ((!selecevtp||iEvtp[*typeitr])&&(!selectopo||iTopo[*topoitr])){
	fEvent = *eventnoitr;
	gMint->GoTo(*runnoitr,*snarlnoitr);
	if (*typeitr>=0&&*typeitr<7&&*topoitr>=0&&*topoitr<5){
	  string tmp = evtpcode[*typeitr] + " | " + topocode[*topoitr];
	  fTYPETOPO->SetText(tmp.c_str());
	}
	fComment->SetText((*commitr).c_str());
      }
      else {this->PrevSelEvt();}
    }
  }
}


void ParticleDisplay::GhostGraph(TGraph *gr)
{
  if(gr == 0) return;

  gr->SetMarkerSize(0);
  gr->SetMarkerStyle(8);
  gr->SetMarkerColor(kWhite);
}

void ParticleDisplay::ColorGraph(TGraph *gr, int color, 
                                   int style, Double_t size)
{
  if(gr == 0) return;

  gr->SetMarkerSize(size);
  gr->SetMarkerStyle(style);
  gr->SetMarkerColor(color);
}
