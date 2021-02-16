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

#include "NueAna/Display/NueDisplayModule.h"
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

#include "DataUtil/GetDetector.h"
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
#include "MuonRemoval/NtpMRRecord.h"
#include "MuonRemoval/NtpMREvent.h"
#include "MuonRemoval/NtpMRTruth.h"

#include "NueAna/NueAnaTools/NueConvention.h"

#include <fstream>
#include <cmath>

using namespace DataUtil;
using namespace SigC;

// Boiler plate for using the Message Service.
#include "MessageService/MsgService.h"
CVSID("$Id: NueDisplayModule.cxx,v 1.115 2010/07/26 22:31:32 mho Exp $");

// Boiler plate needed for us to be a Job Module.
#include "JobControl/JobCModuleRegistry.h"
JOBMODULE(NueDisplayModule,"NueDisplayModule","Example of adding a Midad display in a Job Module\n");

#include "DatabaseInterface/DbiResultPtr.tpl"

const string evtpcode[] = {"N/A", "mu","e","NC","mu/NC?","e/NC?","???"};
const string topocode[] = {"N/A", "QE","RES","DIS","???"};


NueDisplayModule::NueDisplayModule()
  : fCanvas0(0)
  , fUVview(0)
  , fButtonPad(0)
  , fHistPad(0)
  , fInfo0(0)
  , fStdHepCan(0)
  , fUZview(0)
  , fVZview(0)
  , fTrkUZ(0)
  , fTrkVZ(0)
  , fShwUZ(0)
  , fShwVZ(0)
  , fSlcUZ(0)
  , fSlcVZ(0)
  , fSteelOutline(0)
  , pu1_outline(0)
  , fu1_outline(0)
  , pv1_outline(0)
  , fv1_outline(0)
  , pu2_outline(0)
  , fu2_outline(0)
  , pv2_outline(0)
  , fv2_outline(0)
  , stdhepinfo()
  , shia(stdhepinfo)
  , fCanvas1(0)
  , fHistcolz(0)
  , fHistlego(0)
  , fUZcolz(0)
  , fVZcolz(0)
  , fUZlego(0)
  , fVZlego(0)
  , ifixc(0)
  , ifixl(0)
  , evthighest_z(30)
  , evtlowest_z(0)
  , evthighest_t0(4.05)
  , evtlowest_t0(-4.05)
  , evthighest_t1(4.05)
  , evtlowest_t1(-4.05)
  , fFracVar_plots(0)
  , fFracVar_info(0)
  , highest_z(30)
  , lowest_z(0)
  , highest_t0(4.05)
  , lowest_t0(-4.05)
  , highest_t1(4.05)
  , lowest_t1(-4.05)
  , fracvars()
  , fva(fracvars)
  , fShwfit_plots(0)
  , fShwfit_plots_sub(0)
  , fShwfit_info(0)
  , shwfit()
  , sfa(shwfit)
  , hitcalc()
  , hca(hitcalc)
  , fAngClusterFitAna_plots(0)
  , angcluster()
  , aca(angcluster)
  , angclusterfit()
  , acfa(angclusterfit)
  , fCanvas5(0)
  , ssGraphU(0)
  , ssGraphV(0)
  , cluLegU(0)
  , cluLegV(0)
  , fCanvas6(0)
  , fCanvas7(0)
  , fSelectPad1(0)
  , fSelectPad2(0)
  , fCanvas8(0)
  , TimeHst(0)
  , TimeHstTrk(0)
  , TimeHstShw(0)
  , TimeHstUV(0)
  , TimeHstTrkU(0)
  , TimeHstTrkV(0)
  , TimeHstShwU(0)
  , TimeHstShwV(0)
  , TimeHst2(0)
  , TimeHstTrk2(0)
  , TimeHstShw2(0)
  , TimeHst2UV(0)
  , TimeHstTrk2U(0)
  , TimeHstTrk2V(0)
  , TimeHstShw2U(0)
  , TimeHstShw2V(0)
  , gr_dtds(0)
  , fCanvas9(0)
  , vzEventOverlay(0)
  , uzEventOverlay(0)
  , uzSliceOverlay(0)
  , vzSliceOverlay(0)
  , leg(0)
  , fPlusMinusOne(0)
  , fPlusMinusTwo(0)
  , fFullSnarl(0)
  , kShowRange(-1)
  , fMRuShowAll(0)
  , fMRuShowMR(0)
  , fMRuShowOld(0)                                                       
  , fMRdShowAll(0)
  , fMRdShowTrueMu(0)
  , fMRdShowTrueShw(0)
  , fMRdShowScaled(0)
  , fMRdShowReco(0)
  , clickbutton(0)
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
  , SIGMAPMEU(1.)
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
  , iFileW(0)
  , hitlog(0)
  , icomm(0)
  , passfid(0)
  , passtrk(0)
  , passtrklike(0)
  , passshw(0)
  , preselec(0)
  , iIO(0)
  , selecevtp(0)
  , selectopo(0)
  , fRel(ReleaseType::kUnknown)
{	 
  ftrkshw = new TNtuple("trkshw","trkshw","x:y:type");
  info1 = new TLatex();
  info2 = new TLatex();
  info3 = new TLatex();
  info4 = new TLatex();
  info41 = new TLatex();
  info5 = new TLatex();
  info6 = new TLatex();
  info7 = new TLatex();
  info8 = new TLatex();
  info9 = new TLatex();
  info10 = new TLatex();
  info11 = new TLatex();
  info12 = new TLatex();
  info13 = new TLatex();
  mrInfo1 = new TLatex();
  mrInfo2 = new TLatex();
  mrInfo3 = new TLatex();
  mrInfo4 = new TLatex();
  mrInfo5 = new TLatex();
  mrInfo6 = new TLatex();
  mrInfo7 = new TLatex();
  mrInfo8 = new TLatex();
  mrInfo9 = new TLatex();


  mcvtx_u = new TMarker();
  mcvtx_u->SetMarkerStyle(2);
  mcvtx_v = new TMarker();
  mcvtx_v->SetMarkerStyle(2);
  srvtx_u = new TMarker();
  srvtx_u->SetMarkerStyle(5);
  srvtx_v = new TMarker();
  srvtx_v->SetMarkerStyle(5);
  srvtx_xy = new TMarker();
  srvtx_xy->SetMarkerStyle(5);
  fva.SetDisplay(1);

  for (int i = 0; i<7; i++){ iEvtp[i] = 0;}
  for (int i = 0; i<5; i++){ iTopo[i] = 0;}

  tfit_dt_ds_pos = new TF1("tfit_dt_ds_pos","(1./299792458.0)*x+[0]",0,30) ;
  tfit_dt_ds_neg = new TF1("tfit_dt_ds_neg","(-1./299792458.0)*x+[0]",0,30) ;

                                                                    
  showStpTrueMu = false;
  showStpTrueShw = false;
  showStpScaled = false;
  showStpRecoTrk = false;
  showStpAll = true;

  showMREvent = true;
  showOrigEvent = true;
  showNewEvent = true;
  
  //gStyle->SetTitleH(0.04);
  gStyle->SetPadTopMargin(0.1);
}

NueDisplayModule::~NueDisplayModule()
{
  if (!iIO){ //record decisions
    if (iFileW){
      if (!kScanMode){//nornal mode
	if (!icomm && hitlog) outfile<<endl;
	outfile<<endl;
	outfile.close();
      }
      else {//'fast scan' mode
	if (fSnarl>0&&!hitlog){
	  outfile<<RunNo_old<<" "<<SubRunNo_old<<" "<<fSnarl<<" "<<fEvent_old<<" "<<ievtp<<" "<<itopo<<"  "<<fComment->GetText()<<endl;
	  outfile<<endl;
	  outfile.close();
	}
	else if (fSnarl>0&&hitlog){
	  if(!icomm) outfile<<endl;
	  outfile<<endl;
	  outfile.close();
	}
      }
    }
  }
}

void NueDisplayModule::BeginRun()
{
  cout<<"In NueDisplayModule beginrun()"<<endl;
  if (gMint && !fCanvas0 && !fCanvas1) this->BuildDisplay();
}

JobCResult NueDisplayModule::Ana(const MomNavigator *mom)
{
  if (gMint) {                // paranoia
    if (&gMint->GetJobC().Mom != mom) {
      MSG("NueDisplayModule",Msg::kError) << "Module's mom and JobC's mom differ: "
				     << (void*)&gMint->GetJobC().Mom << " != "
				     << (void*) mom << endl;
    }
  }
  
  if (!fCanvas0||!fCanvas1) {
    // No canvas - we can't do anything, but otherwise nothing
    // to say that there is a problem, so just return "ok".
    MSG("NueDisplayModule",Msg::kWarning) << "I have no canvas!\n";
    return JobCResult::kAOK;
  }
  
  //get data from mom
  foundST=false;
  foundSR=false;
  foundMC=false;
  foundTH=false;
  foundMR=false;
  foundSTOld=false;
  st = 0;
  sr = 0;
  mc = 0;
  th = 0;
  mr = 0;
  stOld = 0;

  VldContext vc;  
  st=dynamic_cast<NtpStRecord *>(mom->GetFragment("NtpStRecord", "Primary"));
  if(!st)
  {
    MSG("NueDisplayModule",Msg::kError) << "Can't find primary NtpStRecord.... looking for any!\n";
    st=dynamic_cast<NtpStRecord *>(mom->GetFragment("NtpStRecord"));
  }


  if(st){
    foundST=true;
    vc=st->GetHeader().GetVldContext();

    //MHO - This doesn't work.
    /*    string relName = st->GetTitle();
    string reco = relName.substr(0,relName.find_first_of("("));
    if(reco == "CEDAR") fRel = ReleaseType::kCedar;
    else if(reco == "DOGWOOD") fRel = ReleaseType::kDogwood;
    else fRel = ReleaseType::kBirch;*/
    ReleaseType::Release_t relName = st->GetJobHistory().GetProdReleaseType();
    if (ReleaseType::IsBirch(relName)) fRel = ReleaseType::kBirch;
    if (ReleaseType::IsCedar(relName)) fRel = ReleaseType::kCedar;
    if (ReleaseType::IsDogwood(relName)) fRel = ReleaseType::kDogwood;

    if(vc.GetSimFlag() == SimFlag::kMC){
      NtpMCGenInfo* genInfo = &(st->mchdr.geninfo);
      if(strstr(genInfo->codename.c_str(), "daikon") == 0) {
	fRel += ReleaseType::kCarrot;
      }
      if(strstr(genInfo->codename.c_str(), "daikon") != 0) {
	fRel += ReleaseType::kDaikon;
      }
    }
    if(vc.GetSimFlag() != SimFlag::kMC) fRel += ReleaseType::kData;

    //check for MR:
    mr=dynamic_cast<NtpMRRecord *>(mom->GetFragment("NtpMRRecord"));
    if(mr){
      if(ReleaseType::IsBirch(fRel)) {
	stOld=dynamic_cast<NtpStRecord *>(mom->GetFragment("NtpStRecord",
							   "NtpStRecordOld"));
      }
      else if(ReleaseType::IsCedar(fRel)||ReleaseType::IsDogwood(fRel)){
	stOld=st;
	st=dynamic_cast<NtpStRecord *>(mom->GetFragment("NtpStRecord",
							"MuonRemoved"));
      }
      if(stOld) foundMR = true;
    }
  }
  else {
    fRel = ReleaseType::kUnknown;
    sr = dynamic_cast<NtpSRRecord *>(mom->GetFragment("NtpSRRecord"));
    if (sr) {
      foundSR = true;
      vc=sr->GetHeader().GetVldContext();        
    }
    mc = static_cast<NtpMCRecord *>(mom->GetFragment("NtpMCRecord"));
    if (mc) {
      foundMC = true;
      th = static_cast<NtpTHRecord *>(mom->GetFragment("NtpTHRecord"));
      if (th) foundTH = true;
    }
  }

  if(!foundSR&&!foundST) {
    MSG("NueDisplayModule",Msg::kError)<<"Got Nothing to display"<<endl;
      return JobCResult::kFailed;
   }
    

  // Load DB CalMIPCalibration and set it for packages that need it


   if(!foundmeu) 
   {
    //Setup and give context
    DbiResultPtr<CalMIPCalibration> fResPtr;
    fResPtr.NewQuery(vc,10);
    
    //Get Scale
    if(fResPtr.GetNumRows()>0)
    {
     const CalMIPCalibration* mipcal = fResPtr.GetRow(0);
     SIGMAPMEU = mipcal->GetScale();
     foundmeu=true;
    }
   }

   sfa.SetParams(SIGMAPMEU);
   sfa.SetCutParams(kDPlaneCut,kPhStripCut,kPhPlaneCut, kCPhPlaneCut);
   fva.SetParams(SIGMAPMEU);
   hca.SetParams(SIGMAPMEU);
   //aneia.SetParams(SIGMAPMEU);

   TObject *obj=0;
   TIter objiter = mom->FragmentIter();
   vector<NueRecord *> vr;
   //vector<NuePID *> vpid;
   while((obj=objiter.Next())){
    const char *cn=obj->ClassName();
    if(strcmp(cn,"NueRecord")==0){
      NueRecord *nr = dynamic_cast<NueRecord *>(obj);
      MSG("NueDisplayModule",Msg::kDebug)<<"Found a NueRecord in MOM"
					 <<" Snarl "<<nr->GetHeader().GetSnarl()
					 <<" Event "<<nr->GetHeader().GetEventNo()<<endl;
      vr.push_back(nr);
    }
    else if(strcmp(cn,"NuePID")==0){
      //NuePID *npid  = dynamic_cast<NuePID *>(obj);
      MSG("NueDisplayModule",Msg::kDebug)<<"Found a NuePID in MOM"
	//     <<" Snarl "<<npid->GetHeader().GetSnarl()
	//     <<" Event "<<npid->GetHeader().GetEventNo()<< "PID: " << npid->IsNue i
					 << endl;
      //vpid.push_back(npid);
    }
    else{
      continue;
    }
  }
  
  
  //so, mom will match up snarls for us,
  //but, we have to match up events for ourselves.
  MSG("NueDisplayModule",Msg::kDebug)<<"Starting to match records"<<endl;
  MSG("NueDisplayModule",Msg::kDebug)<<"found "<<vr.size()<<" NueR's"<<endl;
  foundvrmatch=false;
  for(unsigned int i=0;i<vr.size();i++){
    int event = vr[i]->GetHeader().GetEventNo();
    if(event==fEvent){
      nr=vr[i];
      foundvrmatch=true;
      break;
    }
  }	
  if(foundvrmatch) MSG("NueDisplayModule",Msg::kDebug)<<"Found vr"<<endl;
  foundpidmatch=false;
/*  for(unsigned int i=0;i<vpid.size();i++){
    int event = vpid[i]->GetHeader().GetEventNo();
    if(event==fEvent){
      pid=vpid[i];
      foundpidmatch=true;
      break;
    }
  }	
*/
  if(foundvrmatch) MSG("NueDisplayModule",Msg::kDebug)<<"Found pid"<<endl; 
  
  MSG("NueDisplayModule",Msg::kDebug)<<"Matches?"<<endl;
  if(foundpidmatch){
    MSG("NueDisplayModule",Msg::kDebug)<<"found pid"<<endl;
  }
  if(foundvrmatch){
    MSG("NueDisplayModule",Msg::kDebug)<<"found record"<<endl;
  }
  
  if (!iIO) fEvent = 0; //in the Write mode
  this->UpdateDisplay(foundvrmatch, foundpidmatch);
  
  return JobCResult::kAOK;
}

void NueDisplayModule::BuildDisplay()
{
  MSG("NueDisplayModule",Msg::kDebug)<<"In BuildDisplay"<<endl;

  //get detector/simflag/run information

  // Get Ugli for later.  Bail immediately if fail
  UgliGeomHandle ugh = gMint->GetUgliGeomHandle();
  if (! ugh.IsValid()) {
    MSG("NueDisplayModule",Msg::kWarning) << "Got invalid Ugli\n";
    //return;
  }

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
  
  //CanvasSignals* cs_reco = 0;
  CanvasSignals* cs0 = 0;
  CanvasSignals* cs1 = 0;
  CanvasSignals* cs2 = 0;
  CanvasSignals* cs3 = 0;
  CanvasSignals* cs4 = 0;
  CanvasSignals* cs5 = 0;
  CanvasSignals* cs6 = 0;
  CanvasSignals* cs7 = 0;
  CanvasSignals* cs8 = 0;
  CanvasSignals* cs9 = 0;
  CanvasSignals* mrC1 = 0;

  PageDisplay* pd = gMint->GetDisplay();
  //get Button Box on the left
  GuiBox *fBBox = pd->GetButtonBox();
  // No pre-existing display, so make one to our size
  if (!pd) {
    MSG("NueDisplayModule",Msg::kDebug)<<"No display, making one"<<endl;
    pd = gMint->SpawnDisplay(width,height);
  }

  //////////////////////////////////////////////////////////////////////

  //add buttons in the Button Box
  //"Next Event" "Prev Event"
  fNextEventbut = pd->AddButton("Next Event ");
  fNextEventbut->clicked.connect(slot_class(*this,&NueDisplayModule::NextEvent));
  fNextEventbut = pd->AddButton("Prev Event ");
  fNextEventbut->clicked.connect(slot_class(*this,&NueDisplayModule::PrevEvent));

  //Event No
  fEventNo = manage(new GuiLabel(*fBBox," "));
  fEventNo->SetLayoutHints(kLHintsExpandX);
  fBBox->Add(*fEventNo);
  fEventEntry = pd->AddEntry(" ");
  fEventEntry->activated.connect(slot_class(*this,&NueDisplayModule::GotoEvent));
  
  //preselection
  fCuts = pd->AddButton("Cuts: OFF");
  fCuts->clicked.connect(slot_class(*this,&NueDisplayModule::SetCuts));

  //I/O
  fIO = pd->AddButton("Write");
  fIO->clicked.connect(slot_class(*this,&NueDisplayModule::SetMode));

  //"NextSelEvt" "PreSelEvt"
  fNextSelEvt = pd->AddButton("NextSelEvt");
  fNextSelEvt->clicked.connect(slot_class(*this,&NueDisplayModule::NextSelEvt));
  fPrevSelEvt = pd->AddButton("PrevSelEvt");
  fPrevSelEvt->clicked.connect(slot_class(*this,&NueDisplayModule::PrevSelEvt));

  if(!kTestMode){
    fMCTruth = pd->AddButton("MC Truth");
    fMCTruth->clicked.connect(slot_class(*this,&NueDisplayModule::showmctruth));

    fFixMCInfo = pd->AddButton("Fix MC Info");
    fFixMCInfo->clicked.connect(slot_class(*this,&NueDisplayModule::fixmcinfo));
  }

  if(kIntReco){
    fIntRecoCalc2 = pd->AddButton("IR: Re-Draw ");
    fIntRecoCalc2->clicked.connect(slot_class(*this,
					     &NueDisplayModule::SetUpStripButtons));
    fIntRecoCalc1 = pd->AddButton("IR: Calculate");
    fIntRecoCalc1->clicked.connect(slot_class(*this,
					     &NueDisplayModule::IntRecoCalc));
  }

  fPrint = pd->AddButton("Print");
  fPrint->clicked.connect(slot_class(*this,&NueDisplayModule::PrintPads));

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
  
  //voting stuff
  fVBox1 = manage(new GuiBox(*fHBox,kVerticalFrame));
  fVBox1->SetLayoutHints(kLHintsExpandY);
  fHBox->Add(*fVBox1);

  //event type
  evtp1 =  manage(new GuiTextButton(*fVBox1,"mu"));
  evtp1->SetLayoutHints(kLHintsExpandX);
  fVBox1->Add(*evtp1);
  evtp1->clicked.connect(slot_class(*this,&NueDisplayModule::vote10));

  evtp2 =  manage(new GuiTextButton(*fVBox1,"e"));
  evtp2->SetLayoutHints(kLHintsExpandX);
  fVBox1->Add(*evtp2);
  evtp2->clicked.connect(slot_class(*this,&NueDisplayModule::vote20));

  evtp3 =  manage(new GuiTextButton(*fVBox1,"NC"));
  evtp3->SetLayoutHints(kLHintsExpandX);
  fVBox1->Add(*evtp3);
  evtp3->clicked.connect(slot_class(*this,&NueDisplayModule::vote30));

  evtp4 =  manage(new GuiTextButton(*fVBox1,"mu/NC?"));
  evtp4->SetLayoutHints(kLHintsExpandX);
  fVBox1->Add(*evtp4);
  evtp4->clicked.connect(slot_class(*this,&NueDisplayModule::vote40));

  evtp5 =  manage(new GuiTextButton(*fVBox1,"e/NC?"));
  evtp5->SetLayoutHints(kLHintsExpandX);
  fVBox1->Add(*evtp5);
  evtp5->clicked.connect(slot_class(*this,&NueDisplayModule::vote50));

  evtp6 =  manage(new GuiTextButton(*fVBox1,"???"));
  evtp6->SetLayoutHints(kLHintsExpandX);
  fVBox1->Add(*evtp6);
  evtp6->clicked.connect(slot_class(*this,&NueDisplayModule::vote60));

  fVBox2 = manage(new GuiBox(*fHBox,kVerticalFrame,2));
  fVBox2->SetLayoutHints(kLHintsExpandY);
  fHBox->Add(*fVBox2);

  //topology
  topo1 =  manage(new GuiTextButton(*fVBox2,"QE"));
  topo1->SetLayoutHints(kLHintsExpandX);
  fVBox2->Add(*topo1);
  topo1->clicked.connect(slot_class(*this,&NueDisplayModule::vote01));

  topo2 =  manage(new GuiTextButton(*fVBox2,"RES"));
  topo2->SetLayoutHints(kLHintsExpandX);
  fVBox2->Add(*topo2);
  topo2->clicked.connect(slot_class(*this,&NueDisplayModule::vote02));

  topo3 =  manage(new GuiTextButton(*fVBox2,"DIS"));
  topo3->SetLayoutHints(kLHintsExpandX);
  fVBox2->Add(*topo3);
  topo3->clicked.connect(slot_class(*this,&NueDisplayModule::vote03));

  topo4 =  manage(new GuiTextButton(*fVBox2,"???"));
  topo4->SetLayoutHints(kLHintsExpandX);
  fVBox2->Add(*topo4);
  topo4->clicked.connect(slot_class(*this,&NueDisplayModule::vote04));

  fEncoded = manage(new GuiLabel(*fVBox2," "));
  fEncoded->SetLayoutHints(kLHintsExpandX);
  fVBox2->Add(*fEncoded);
  //updateEncoded();

  vote = manage(new GuiTextButton(*fVBox2,"Log Details"));
  vote->SetLayoutHints(kLHintsExpandX);
  fVBox2->Add(*vote);
  vote->clicked.connect(slot_class(*this,&NueDisplayModule::logvote));

  fVBox3 = manage(new GuiBox(*fHBox,kVerticalFrame));
  fHBox->Add(*fVBox3);

  //make comments?
  fHBox3 = manage(new GuiBox(*fVBox3,kHorizontalFrame));
  fHBox3->SetLayoutHints(kLHintsExpandX);
  fVBox3->Add(*fHBox3);
  fComLab = manage(new GuiLabel(*fHBox3," Comments? "));
  fComLab->SetLayoutHints(kLHintsLeft);
  fHBox3->Add(*fComLab);
  fComment = manage(new GuiTextEntry(*fHBox3,""));
  fComment->SetLayoutHints(kLHintsExpandX);
  fHBox3->Add(*fComment);
  fComment->activated.connect(slot_class(*this,&NueDisplayModule::makecomment));

  //open file to read
  fHBox4 = manage(new GuiBox(*fVBox3,kHorizontalFrame));
  fHBox4->SetLayoutHints(kLHintsExpandX);
  fVBox3->Add(*fHBox4);
  fFileLab = manage(new GuiLabel(*fHBox4," File to Open? "));
  fFileLab->SetLayoutHints(kLHintsLeft);
  fHBox4->Add(*fFileLab);
  fFileEnt = manage(new GuiTextEntry(*fHBox4,""));
  fFileEnt->SetLayoutHints(kLHintsExpandX);
  fHBox4->Add(*fFileEnt);
  fFileEnt->activated.connect(slot_class(*this,&NueDisplayModule::OpenFileRead));

  //event type
  fHBox5 = manage(new GuiBox(*fVBox3,kHorizontalFrame));
  fHBox5->SetLayoutHints(kLHintsExpandX);
  fVBox3->Add(*fHBox5);
  fTypeLab = manage(new GuiLabel(*fHBox5," TYPE         "));
  fTypeLab->SetLayoutHints(kLHintsLeft);
  fHBox5->Add(*fTypeLab);
  Evtp1 = manage(new GuiCheckButton(*fHBox5,"mu"));
  Evtp1->SetLayoutHints(kLHintsLeft);
  fHBox5->Add(*Evtp1);
  Evtp1->clicked.connect(slot_class(*this,&NueDisplayModule::selec10));
  Evtp2 = manage(new GuiCheckButton(*fHBox5,"e"));
  Evtp2->SetLayoutHints(kLHintsLeft);
  fHBox5->Add(*Evtp2);
  Evtp2->clicked.connect(slot_class(*this,&NueDisplayModule::selec20));
  Evtp3 = manage(new GuiCheckButton(*fHBox5,"NC"));
  Evtp3->SetLayoutHints(kLHintsLeft);
  fHBox5->Add(*Evtp3);
  Evtp3->clicked.connect(slot_class(*this,&NueDisplayModule::selec30));
  Evtp4 = manage(new GuiCheckButton(*fHBox5,"mu/NC?"));
  Evtp4->SetLayoutHints(kLHintsLeft);
  fHBox5->Add(*Evtp4);
  Evtp4->clicked.connect(slot_class(*this,&NueDisplayModule::selec40));
  Evtp5 = manage(new GuiCheckButton(*fHBox5,"e/NC?"));
  Evtp5->SetLayoutHints(kLHintsLeft);
  fHBox5->Add(*Evtp5);
  Evtp5->clicked.connect(slot_class(*this,&NueDisplayModule::selec50));
  Evtp6 = manage(new GuiCheckButton(*fHBox5,"???"));
  Evtp6->SetLayoutHints(kLHintsLeft);
  fHBox5->Add(*Evtp6);
  Evtp6->clicked.connect(slot_class(*this,&NueDisplayModule::selec60));
  Evtp7 = manage(new GuiCheckButton(*fHBox5,"N/A"));
  Evtp7->SetLayoutHints(kLHintsLeft);
  fHBox5->Add(*Evtp7);
  Evtp7->clicked.connect(slot_class(*this,&NueDisplayModule::selec70));

  //topology
  fHBox6 = manage(new GuiBox(*fVBox3,kHorizontalFrame));
  fHBox6->SetLayoutHints(kLHintsExpandX);
  fVBox3->Add(*fHBox6);
  fTopoLab = manage(new GuiLabel(*fHBox6," TOPOLOGY "));
  fTopoLab->SetLayoutHints(kLHintsLeft);
  fHBox6->Add(*fTopoLab);
  Topo1 = manage(new GuiCheckButton(*fHBox6,"QE"));
  Topo1->SetLayoutHints(kLHintsLeft);
  fHBox6->Add(*Topo1);
  Topo1->clicked.connect(slot_class(*this,&NueDisplayModule::selec01));
  Topo2 = manage(new GuiCheckButton(*fHBox6,"RES"));
  Topo2->SetLayoutHints(kLHintsLeft);
  fHBox6->Add(*Topo2);
  Topo2->clicked.connect(slot_class(*this,&NueDisplayModule::selec02));
  Topo3 = manage(new GuiCheckButton(*fHBox6,"DIS"));
  Topo3->SetLayoutHints(kLHintsLeft);
  fHBox6->Add(*Topo3);
  Topo3->clicked.connect(slot_class(*this,&NueDisplayModule::selec03));
  Topo4 = manage(new GuiCheckButton(*fHBox6,"???"));
  Topo4->SetLayoutHints(kLHintsLeft);
  fHBox6->Add(*Topo4);
  Topo4->clicked.connect(slot_class(*this,&NueDisplayModule::selec04));
  Topo5 = manage(new GuiCheckButton(*fHBox6,"N/A"));
  Topo5->SetLayoutHints(kLHintsLeft);
  fHBox6->Add(*Topo5);
  Topo5->clicked.connect(slot_class(*this,&NueDisplayModule::selec05));

  fHBox7 = manage(new GuiBox(*fVBox3,kHorizontalFrame));
  fHBox7->SetLayoutHints(kLHintsExpandX);
  fVBox3->Add(*fHBox7);
  fFileR = manage(new GuiLabel(*fHBox7,""));
  fHBox7->Add(*fFileR);
  fTYPETOPO = manage(new GuiLabel(*fHBox7,""));
  fHBox7->Add(*fTYPETOPO);
  
  fHBox8 = manage(new GuiBox(*fVBox3,kHorizontalFrame));
  fHBox8->SetLayoutHints(kLHintsExpandX);
  fVBox3->Add(*fHBox8);
  fFileW = manage(new GuiLabel(*fHBox8,""));
  fHBox8->Add(*fFileW);

  fHBox9 = manage(new GuiBox(*fVBox3,kHorizontalFrame));
  fHBox9->SetLayoutHints(kLHintsExpandX);
  fVBox3->Add(*fHBox9);

  fNextEventbut1 = manage(new GuiTextButton(*fHBox9,"Next Event"));
  //fNextEventbut1->SetLayoutHints(kLHintsExpandX);
  fHBox9->Add(*fNextEventbut1);
  fNextEventbut1->clicked.connect(slot_class(*this,&NueDisplayModule::NextEvent));

  fPrevEventbut1 = manage(new GuiTextButton(*fHBox9,"Prev Event"));
  //fPrevEventbut1->SetLayoutHints(kLHintsExpandX);
  fHBox9->Add(*fPrevEventbut1);
  fPrevEventbut1->clicked.connect(slot_class(*this,&NueDisplayModule::PrevEvent));


  //reco info
  fHBox2 = manage(new GuiBox(*fHBox,kHorizontalFrame));
  fHBox->Add(*fHBox2);

  fRecoInfo = manage(new GuiTextView(*fHBox2,1,1));
  fHBox2->Add(*fRecoInfo);
  fRecoInfo->AddLine("Reco Information");



  ///////////////////////////////////////////////////////////////    

  //add user-defined canvas

  cs0 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));
  cs1 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));
  cs2 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));
  cs3 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));
  cs4 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));
  if(kDrawClu) cs5 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));
  if(kIntReco) {
    cs7 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));
    cs6 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));
  }
  cs8 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));
  cs9 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));
  mrC1 = dynamic_cast<CanvasSignals*>(pd->AddPage("UserCanvas"));

  if (!cs0||!cs1||!cs2||!cs3||!cs4||!cs8|!cs9|!mrC1) {
    MSG("NueDisplayModule",Msg::kWarning) << "Failed to get NueCanvas's CanvasSignals\n";
    return;
  }

    //fCanvas_reco = &cs_reco->GetCanvas();
    fCanvas0 = &cs0->GetCanvas();
    fCanvas1 = &cs1->GetCanvas();
    fCanvas2 = &cs2->GetCanvas();
    fCanvas3 = &cs3->GetCanvas();
    fCanvas4 = &cs4->GetCanvas();
    if(kDrawClu) {
      fCanvas5 = &cs5->GetCanvas();    
      fCanvas5->SetHighLightColor(0);
    }
    if(kIntReco) {
      fCanvas6 = &cs6->GetCanvas();
      fCanvas6->SetHighLightColor(0);
      fCanvas7 = &cs7->GetCanvas();
      fCanvas7->SetHighLightColor(0);
    }
    fCanvas8 = &cs8->GetCanvas();
    fCanvas9 = &cs9->GetCanvas();
    fMRCanvas1 = &mrC1->GetCanvas();
    
    //please set highlightcolor of all the canvases to be 0
    //otherwise buttons on the second page will be screwed
    fCanvas0->SetHighLightColor(0);
    fCanvas1->SetHighLightColor(0);
    fCanvas2->SetHighLightColor(0);
    fCanvas3->SetHighLightColor(0);
    fCanvas4->SetHighLightColor(0);
    fCanvas8->SetHighLightColor(0);
    fCanvas9->SetHighLightColor(0);
    fMRCanvas1->SetHighLightColor(0);

    //main canvas
    //fCanvas0->Divide(2,1);
    //fCanvas0->SetName("fCanvas0");
    //fCanvas1->Divide(2,1);
    fCanvas2->Divide(2,1);
    fCanvas3->Divide(2,2);
    fCanvas4->Divide(2,2);    
    fCanvas8->Divide(2,2);
    fCanvas8->SetName("fCanvas8");
    fCanvas9->Divide(2,2);
//    fMRCanvas1->Divide(2,2);

    if(kDrawClu) {
      fCanvas5->cd();
      fCanvas5->SetName("fCanvas5");
      TPad *temp_pad = new TPad("c5_1","c5_1",0.0,0.5,0.42,1);
      temp_pad->SetNumber(1); temp_pad->Draw();
      temp_pad = new TPad("c5_2","c5_2",0.42,0.5,0.84,1);
      temp_pad->SetNumber(2); temp_pad->Draw();
      temp_pad = new TPad("c5_3","c5_3",0.84,0.5,1.0,1);
      temp_pad->SetNumber(3); temp_pad->Draw();
      temp_pad = new TPad("c5_4","c5_4",0.0,0.0,0.42,0.5);
      temp_pad->SetNumber(4); temp_pad->Draw();
      temp_pad = new TPad("c5_5","c5_5",0.42,0.0,0.84,0.5);
      temp_pad->SetNumber(5); temp_pad->Draw();
      temp_pad = new TPad("c5_6","c5_6",0.84,0.0,1.0,0.5);
      temp_pad->SetNumber(6); temp_pad->Draw();
    }

    if(!kTestMode){
      fCanvas0->cd();
      fStdHepCan = new TPad("fStdHepCan","fStdHepCan",0.5,0.,1.,0.5);
      fStdHepCan->Draw();
      fStdHepCan->Range(0,0,1,1.1);
      fStdHepCan->cd();
      TPaveText *infoTex1 = new TPaveText(0.05,1.,0.2,1.07);
      infoTex1->AddText("Initial State");
      infoTex1->SetBorderSize(1);
      TPaveText *infoTex2 = new TPaveText(0.3,1.,0.45,1.07);
      infoTex2->AddText("Intermediate");
      infoTex2->SetBorderSize(1);
      TPaveText *infoTex3 = new TPaveText(0.55,1.,0.7,1.07);
      infoTex3->AddText("Final State");
      infoTex3->SetBorderSize(1);
      TPaveText *infoTex4 = new TPaveText(0.8,1.,0.95,1.07);
      infoTex4->AddText("Later Decays");
      infoTex4->SetBorderSize(1);
      infoTex1->SetTextSize(0.05);
      infoTex1->SetTextFont(12);
      infoTex1->SetTextColor(1);
      infoTex2->SetTextSize(0.05);
      infoTex2->SetTextFont(12);
      infoTex2->SetTextColor(1);
      infoTex3->SetTextSize(0.05);
      infoTex3->SetTextFont(12);
      infoTex3->SetTextColor(1);
      infoTex4->SetTextSize(0.05);
      infoTex4->SetTextFont(12);
      infoTex4->SetTextColor(1);
      infoTex1->SetName("infoTex1");
      infoTex2->SetName("infoTex2");
      infoTex3->SetName("infoTex3");
      infoTex4->SetName("infoTex4");
      infoTex1->Draw();
      infoTex2->Draw();
      infoTex3->Draw();
      infoTex4->Draw();
    }

    fCanvas0->cd();
    fButtonPad = new TPad("fButtonPad","Buttons",0,0,0.5,0.05);
    fButtonPad->Draw();
    fButtonPad->cd();
    //the following button can zoom in/out the histograms on the first Canvas. TButton has no access to the class members so I have to pass the parameters through the histogram option. 
    fZoom = new TButton("Zoom","float a1,a2,a3,a4,a5,a6;int i;sscanf(UZview->GetOption(),\"%f,%f,%f,%f,%f,%f,%d\",&a1,&a2,&a3,&a4,&a5,&a6,&i);if(i==0){UZview->GetXaxis()->SetRangeUser(a2,a1);UZview->GetYaxis()->SetRangeUser(a4,a3);VZview->GetXaxis()->SetRangeUser(a2,a1);VZview->GetYaxis()->SetRangeUser(a6,a5);char setoptions[100];sprintf(setoptions,\"%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d\",a1,a2,a3,a4,a5,a6,1);UZview->SetOption(setoptions);} else if(i==1){double max=UZview->GetMaximum();UZview->GetXaxis()->UnZoom();UZview->GetYaxis()->UnZoom();VZview->GetXaxis()->UnZoom();VZview->GetYaxis()->UnZoom();UZview->SetMaximum(max);VZview->SetMaximum(max);char setoptions[100];sprintf(setoptions,\"%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d\",a1,a2,a3,a4,a5,a6,0);UZview->SetOption(setoptions);}fHistPad->cd(1);gPad->Modified();fHistPad->cd(2);gPad->Modified();",0,0,0.1,1);
    fZoom->SetTextSize(0.5);
    fZoom->Draw();

    fCanvas0->cd();
    fHistPad = new TPad("fHistPad","UZ and VZ views",0,0.05,0.5,1);
    //fHistPad->SetTopMargin(1);
    fHistPad->Draw();
    
    fHistPad->cd();
    fHistPad->Divide(1,2);
    
    int trkcolor = 1;
    int shwcolor = 6;
    // U vs. Z view
    fHistPad->cd(1);
    fUZview = new TH2D("UZview","Strip vs. Plane, U view",
		       500,0,30,210,-4.05,4.05);
    fUZview->SetStats(false);
    fUZview->SetXTitle("z(m)");
    fUZview->SetYTitle("u(m)");
    fUZview->GetYaxis()->SetTitleOffset(0.6);
    fUZview->Draw("COLZ");
    srvtx_u->Draw();
    fTrkUZ = new TH2D("TrkUZ","track strips, U view",
		       500,0,30,210,-4.05,4.05);
    fTrkUZ->SetStats(false);
    fTrkUZ->SetLineColor(trkcolor);
    fShwUZ = new TH2D("ShwUZ","shower strips, U view",
		       500,0,30,210,-4.05,4.05);
    fShwUZ->SetStats(false);
    fShwUZ->SetLineColor(shwcolor);
    fShwUZ->SetLineWidth(2);
    fSlcUZ = new TH2D("SlcUZ","Strip vs. Plane, U view",
		       500,0,30,210,-4.05,4.05);
    fSlcUZ->SetStats(false);

    // V vs. Z view
    fHistPad->cd(2);
    fVZview = new TH2D("VZview","Strip vs. Plane, V view",
		       500,0,30,210,-4.05,4.05);
    fVZview->SetStats(false);
    fVZview->SetXTitle("z(m)");
    fVZview->SetYTitle("v(m)");
    fVZview->GetYaxis()->SetTitleOffset(0.6);
    fVZview->Draw("COLZ");
    srvtx_v->Draw();
    fTrkVZ = new TH2D("TrkVZ","track strips, V view",
		       500,0,30,210,-4.05,4.05);
    fTrkVZ->SetStats(false);
    fTrkVZ->SetLineColor(trkcolor);
    fShwVZ = new TH2D("ShwVZ","shower strips, V view",
		       500,0,30,210,-4.05,4.05);
    fShwVZ->SetLineColor(shwcolor);
    fShwVZ->SetStats(false);
    fShwVZ->SetLineWidth(2);
    fSlcVZ = new TH2D("SlcVZ","Strip vs. Plane, V view",
		       500,0,30,210,-4.05,4.05);
    fSlcVZ->SetStats(false);

    // U vs. V view
    fCanvas0->cd();
    fUVview = new TPad("UVview","U vs V view",0.5,0.5,0.725,1);
    fUVview->SetFillStyle(4000);
    
    float t[4];
    ugh.GetTransverseExtent(PlaneView::kU,t[0],t[1]);
    ugh.GetTransverseExtent(PlaneView::kV,t[2],t[3]);
    for (int ind=0; ind<4;++ind) t[ind] = TMath::Abs(t[ind]);
    int maxind = TMath::LocMax(4,t);
    float tsize = TMath::Sqrt(2.0) * t[maxind];
    
    fUVview->Range(-tsize,-tsize,tsize,tsize);
    
    fUVview->Draw();
    
    fInfo0 = new TPad("Info0","Basic information",0.725,0.5,1,1);
    fInfo0->Draw();
    //fInfo0->cd();
    //fInfo0->Divide(1,2);
    fUVview->cd();
    fSteelOutline = new SteelOutline(DataUtil::GetDetector(gMint->GetJobC().Mom));
    if (!fSteelOutline) MSG("NueDisplayModule",Msg::kWarning)<<"Can't get the stell outline."<<endl;
    if (fDetectorType == Detector::kFar){//draw veto shield,copied Jim's EVD code
      //UgliGeomHandle ugh = gMint->GetUgliGeomHandle();
      if (ugh.IsValid()){
	vector<UgliScintPlnHandle> usphv = ugh.GetScintPlnHandleVector();
	vector<UgliScintPlnHandle>::reverse_iterator rit, rdone = usphv.rend();
	VS = new TList;
	VS->SetOwner();
	float cx,cy; 
	StripEnd::StripEnd_t two_sided[] = { StripEnd::kEast,
					     StripEnd::kWest };
	StripEnd::StripEnd_t *side_list = two_sided;
   
	for (rit = usphv.rbegin(); rit != rdone; ++rit) {
	  if (! rit->GetPlexPlaneId().IsVetoShield()) break;
	  int nstrips = rit->NumberOfStrips();
	  for (int istrip=0; istrip < nstrips; istrip++) {
	    PlexStripEndId oneend(rit->GetPlexPlaneId(),istrip,side_list[0]);
	    UgliStripHandle ush = ugh.GetStripHandle(oneend);
	    TVector3 stripxyz0(ush.GlobalPos(-ush.GetHalfLength()));
	    TVector3 stripxyz1(ush.GlobalPos(ush.GetHalfLength()));
	    cx = 0.5*(stripxyz0[0]+stripxyz1[0]);
	    cy = 0.5*(stripxyz0[1]+stripxyz1[1]);     
	    Double_t dy[5] = {-ush.GetHalfThickness(),
			      -ush.GetHalfThickness(),
			      +ush.GetHalfThickness(),
			      +ush.GetHalfThickness(),
			      -ush.GetHalfThickness()};
	    Double_t dx[5] = {-ush.GetHalfWidth(),
			      +ush.GetHalfWidth(),
			      +ush.GetHalfWidth(),
			      -ush.GetHalfWidth(),
			      -ush.GetHalfWidth()};
	    if(rit->GetPlaneView()==PlaneView::kVSWallOnEdge){
	      for (int id=0;id<6;id++){
		double tmp=dy[id];
		dy[id]=dx[id];
		dx[id]=tmp;
	      }
	    }
	    Double_t rot=0.9;
	    if(rit->GetPlaneView()==PlaneView::kVSTopEastSlant){
	      double newdx;
	      double newdy;
	      for (int id=0;id<6;id++){
		newdx=dx[id]*cos(-rot)+dy[id]*sin(-rot);
		newdy=-dx[id]*sin(-rot)+dy[id]*cos(-rot);
		dx[id]=newdx;
		dy[id]=newdy;
	      }
	    }
	    if(rit->GetPlaneView()==PlaneView::kVSTopWestSlant){
	      double newdx;
	      double newdy;
	      for (int id=0;id<6;id++){
		newdx=dx[id]*cos(rot)+dy[id]*sin(rot);
		newdy=-dx[id]*sin(rot)+dy[id]*cos(rot);
		dx[id]=newdx;
		dy[id]=newdy;
	      }
	      
	    }
	    Double_t xbox[5]={cx+dx[0],cx+dx[1],cx+dx[2],cx+dx[3],cx+dx[4]};
	    Double_t ybox[5]={cy+dy[0],cy+dy[1],cy+dy[2],cy+dy[3],cy+dy[4]};
	    TPolyLine * VSstrip = new TPolyLine(5,xbox,ybox);
	    VSstrip->SetLineColor(13);
	    VSstrip->SetLineWidth(1);
	    VS->Add(VSstrip);
	  }
	}
      }
    }
    Double_t l_coil = 0.2;
    Double_t x_coil[5] = {0,l_coil,0,-1*l_coil,0};
    Double_t y_coil[5] = {l_coil,0,-1*l_coil,0,l_coil};
    coil = new TPolyLine(5,x_coil,y_coil);
    coil->SetLineWidth(2);

    if (fDetectorType == Detector::kNear && !pu1_outline){
      //from Mad....
      PlaneOutline po;
      Color_t colu=38;
      Color_t colv=46;
      po.GetOutline(PlaneView::kV, PlaneCoverage::kNearPartial,
		    pv1_outline, pv2_outline);
      po.GetOutline(PlaneView::kV, PlaneCoverage::kNearFull,
		    fv1_outline, fv2_outline);
      po.GetOutline(PlaneView::kU, PlaneCoverage::kNearPartial,
		    pu1_outline, pu2_outline);
      po.GetOutline(PlaneView::kU, PlaneCoverage::kNearFull,
		    fu1_outline, fu2_outline);
      pv1_outline->SetLineColor(colv);
      pu1_outline->SetLineColor(colu);

      fv1_outline->SetLineColor(colv);
      fu1_outline->SetLineColor(colu);
      fv2_outline->SetFillColor(16);
      fu2_outline->SetFillColor(16);
      fv2_outline->SetFillStyle(4020);
      fu2_outline->SetFillStyle(4020);
      
      pv1_outline->SetLineWidth(1);
      pu1_outline->SetLineWidth(1);
      fv1_outline->SetLineWidth(1);
      fu1_outline->SetLineWidth(1);
      fv2_outline->SetLineWidth(1);
      fu2_outline->SetLineWidth(1);

    }

    //magnified images canvas
    fCanvas1->cd();
    fBC = new TPad("fBC","Button Panel for colz display",0,0.05,0.05,1);
    fBC->Draw();
    fBC->cd();

    //add buttons to get rid of low ph hitss 
    //originally 0,1,3,5 pe = 0,60,180,300 sigcor (far) = 0,0.1062,0.3186,0.5410 MEU (far)
    //originally 0,1,3,5 pe = 0,100,300,500 sigcor (near) = 
    cbut[0] = new TButton("PH>0","UZcolz->SetMinimum(0);VZcolz->SetMinimum(0);fHistcolz->cd(1);gPad->Modified();fHistcolz->cd(2);gPad->Modified();for(int i=0;i<4;i++){TButton *cbut=(TButton*)fBC->FindObject(Form(\"cbut%d\",i));if(i==0) cbut->SetFillColor(4); else cbut->SetFillColor(5);cbut->Paint();}",0,0.31,1,.36);
    cbut[0]->SetName("cbut0");
    cbut[1] = new TButton("PH>0.1","UZcolz->SetMinimum(0.10);VZcolz->SetMinimum(0.1);fHistcolz->cd(1);gPad->Modified();fHistcolz->cd(2);gPad->Modified();for(int i=0;i<4;i++){TButton *cbut=(TButton*)fBC->FindObject(Form(\"cbut%d\",i));if(i==1) cbut->SetFillColor(4); else cbut->SetFillColor(5);cbut->Paint();}",0,0.24,1,0.29);
    cbut[1]->SetName("cbut1");
    cbut[2] = new TButton("PH>0.3","UZcolz->SetMinimum(0.30);VZcolz->SetMinimum(0.3);fHistcolz->cd(1);gPad->Modified();fHistcolz->cd(2);gPad->Modified();for(int i=0;i<4;i++){TButton *cbut=(TButton*)fBC->FindObject(Form(\"cbut%d\",i));if(i==2) cbut->SetFillColor(4); else cbut->SetFillColor(5);cbut->Paint();}",0,0.17,1,0.22);
    cbut[2]->SetName("cbut2");
    cbut[3] = new TButton("PH>0.5","UZcolz->SetMinimum(0.50);VZcolz->SetMinimum(0.5);fHistcolz->cd(1);gPad->Modified();fHistcolz->cd(2);gPad->Modified();for(int i=0;i<4;i++){TButton *cbut=(TButton*)fBC->FindObject(Form(\"cbut%d\",i));if(i==3) cbut->SetFillColor(4); else cbut->SetFillColor(5);cbut->Paint();}",0,0.10,1,0.15);
    cbut[3]->SetName("cbut3");
    //cbut[4] = new TButton("Fix","fixc();",0,0.6,1,0.65);
    //cbut[4]->SetName("cbut4");
    for (int icbut = 0; icbut<4; icbut++){
      if (icbut == 0){
	cbut[icbut]->SetFillColor(4);
      }
      else cbut[icbut]->SetFillColor(5);
      cbut[icbut]->SetTextSize(0.5);
      cbut[icbut]->Draw();
    }

    fCanvas1->cd();
    fHistcolz = new TPad("fHistcolz","UZ and VZ views",0.05,0.05,0.525,1);
    fHistcolz->Draw();
    fHistcolz->cd();
    fHistcolz->Divide(1,2);
    fHistcolz->cd(1);
    fUZcolz = new TH2D("UZcolz","Strip vs. Plane, U view",1,0,1,1,0,1);
    //fUZcolz->SetTitle("Strip vs. Plane, U view");
    //fUZcolz->SetBins(1,0,1,1,0,1);
    fUZcolz->SetStats(false);
    fUZcolz->SetXTitle("z(m)");
    fUZcolz->SetYTitle("u(m)");
    fUZcolz->GetYaxis()->SetTitleOffset(0.6);
    fUZcolz->Draw("colz");
    srvtx_u->Draw();
    //mcvtx_u->Draw();

    fHistcolz->cd(2);
    fVZcolz = new TH2D("VZcolz","Strip vs. Plane, V view",1,0,1,1,0,1);
    fVZcolz->SetStats(false);
    fVZcolz->SetXTitle("z(m)");
    fVZcolz->SetYTitle("v(m)");
    fVZcolz->GetYaxis()->SetTitleOffset(0.6);
    //fVZcolz->SetTitle("Strip vs. Plane, V view");
    //fVZcolz->SetBins(1,0,1,1,0,1);
    fVZcolz->Draw("colz");
    srvtx_v->Draw();
    //mcvtx_v->Draw();
    
    fCanvas1->cd();
    fHistlego = new TPad("fHistlego","UZ and VZ views",0.525,0.05,1,1);
    fHistlego->Draw();
    fHistlego->cd();
    fHistlego->Divide(1,2);
    fHistlego->cd(1);
    fUZlego = new TH2D("UZlego","Strip vs. Plane, U view",1,0,1,1,0,1);
    fUZlego->SetStats(false);
    //fUZlego->SetName("fUZlego");
    //fUZlego->SetTitle("Strip vs. Plane, U view");
    //fUZlego->SetBins(1,0,1,1,0,1);
    fUZlego->Draw("lego2");
    fHistlego->cd(2);
    fVZlego = new TH2D("VZlego","Strip vs. Plane, V view",1,0,1,1,0,1);
    fVZlego->SetStats(false);
    //fVZlego->SetName("fVZlego");
    //fVZlego->SetTitle("Strip vs. Plane, V view");
    //fVZlego->SetBins(1,0,1,1,0,1);
    fVZlego->Draw("lego2");

    fCanvas1->cd();
    fBL = new TPad("fBL","Button Panel for lego display",0.525,0,1,0.05);
    fBL->Draw();
    fBL->cd();

    //add buttons to rotate lego plots
    lbut[0] = new TButton("Pan view","fHistlego->cd(1);gPad->SetPhi(30);gPad->SetTheta(30);gPad->Modified();fHistlego->cd(2);gPad->SetPhi(30);gPad->SetTheta(30);gPad->Modified();for(int i=0;i<4;i++){TButton *lbut=(TButton*)fBL->FindObject(Form(\"lbut%d\",i));if(i==0) lbut->SetFillColor(4); else lbut->SetFillColor(5);lbut->Paint();}",0.1,0,0.22,1);
    lbut[0]->SetName("lbut0");
    lbut[1] = new TButton("Side view","fHistlego->cd(1);gPad->SetPhi(15);gPad->SetTheta(30);gPad->Modified();fHistlego->cd(2);gPad->SetPhi(15);gPad->SetTheta(30);gPad->Modified();for(int i=0;i<4;i++){TButton *lbut=(TButton*)fBL->FindObject(Form(\"lbut%d\",i));if(i==1) lbut->SetFillColor(4); else lbut->SetFillColor(5);lbut->Paint();}",0.23,0,0.35,1);
    lbut[1]->SetName("lbut1");
    lbut[2] = new TButton("Front view","fHistlego->cd(1);gPad->SetPhi(75);gPad->SetTheta(60);gPad->Modified();fHistlego->cd(2);gPad->SetPhi(75);gPad->SetTheta(60);gPad->Modified();for(int i=0;i<4;i++){TButton *lbut=(TButton*)fBL->FindObject(Form(\"lbut%d\",i));if(i==2) lbut->SetFillColor(4); else lbut->SetFillColor(5);lbut->Paint();}",0.36,0,0.48,1);
    lbut[2]->SetName("lbut2");
    lbut[3] = new TButton("Top view","fHistlego->cd(1);gPad->SetPhi(30);gPad->SetTheta(75);gPad->Modified();fHistlego->cd(2);gPad->SetPhi(30);gPad->SetTheta(75);gPad->Modified();for(int i=0;i<4;i++){TButton *lbut=(TButton*)fBL->FindObject(Form(\"lbut%d\",i));if(i==3) lbut->SetFillColor(4); else lbut->SetFillColor(5);lbut->Paint();}",0.49,0,0.61,1);
    lbut[3]->SetName("lbut3");
    //lbut[4] = new TButton("Fix","",0.64,0,0.76,1);
    for (int ilbut = 0; ilbut<4; ilbut++){
      if (ilbut == 0){
	lbut[ilbut]->SetFillColor(4);
      }
      else lbut[ilbut]->SetFillColor(5);
      lbut[ilbut]->SetTextSize(0.5);
      lbut[ilbut]->Draw();
    }

    fCanvas1->cd();
    fStdLeg = new TPad("fStdLeg","legend",0,0,0.525,0.05);
    fStdLeg->Draw();
    fStdLeg->cd();
    
    //copied Chris' code
    TLine *line;
    line = new TLine(0.01,0.5,0.04,0.5);
    line->SetLineColor(3);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.11,0.5,0.14,0.5);
    line->SetLineColor(4);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.21,0.5,0.24,0.5);
    line->SetLineColor(2);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.31,0.5,0.34,0.5);
    line->SetLineColor(28);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.41,0.5,0.44,0.5);
    line->SetLineColor(6);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.51,0.5,0.54,0.5);
    line->SetLineColor(7);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.61,0.5,0.64,0.5);
    line->SetLineColor(31);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.71,0.5,0.74,0.5);
    line->SetLineColor(9);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.81,0.5,0.84,0.5);
    line->SetLineColor(5);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.91,0.5,0.94,0.5);
    line->SetLineColor(1);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw();


    TLatex *tex;
    tex = new TLatex(0.06,0.3,"e");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();
    tex = new TLatex(0.16,0.3,"#mu");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();
    tex = new TLatex(0.26,0.3,"p");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();
    tex = new TLatex(0.36,0.3,"n");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();
    tex = new TLatex(0.46,0.3,"#pi^{+/-}");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();
    tex = new TLatex(0.56,0.3,"#pi^{0}");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();
    tex = new TLatex(0.66,0.3,"K^{+/-/0}");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();
    tex = new TLatex(0.76,0.3,"#gamma");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();
    tex = new TLatex(0.86,0.3,"tau");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();
    tex = new TLatex(0.96,0.3,"#nu");
    tex->SetTextFont(132);
    tex->SetTextSize(0.7);
    tex->Draw();

    gPad->Modified();
    fCanvas1->Update();

    fCanvas2->cd(1);
    fFracVar_plots = new TPad("fFracVar_plots","FracVar",0,0,1,1);
    fFracVar_plots->Draw();
    fFracVar_plots->Divide(1,2);
    fCanvas2->cd(2);
    fReco_plots = new TPad("fReco_plots","Reco",0,0,1,1);
    fReco_plots->Draw();
    fReco_plots->Divide(1,2);
    fReco_plots->cd(1);
    fSlcUZ->SetXTitle("z(m)");
    fSlcUZ->SetYTitle("u(m)");
    fSlcUZ->GetYaxis()->SetTitleOffset(0.6);
    fSlcUZ->Draw("colz");
    fShwUZ->SetMaximum(1);
    fShwUZ->SetMinimum(0);
    fShwUZ->Draw("box same");
    fTrkUZ->SetMaximum(1);
    fTrkUZ->SetMinimum(0);
    fTrkUZ->Draw("box same");

    fReco_plots->cd(2);
    fSlcVZ->SetXTitle("z(m)");
    fSlcVZ->SetYTitle("v(m)");
    fSlcVZ->GetYaxis()->SetTitleOffset(0.6);
    fSlcVZ->Draw("colz");
    fShwVZ->SetMaximum(1);
    fShwVZ->SetMinimum(0);
    fShwVZ->Draw("box same");
    fTrkVZ->SetMaximum(1);
    fTrkVZ->SetMinimum(0);
    fTrkVZ->Draw("box same");

    fShwfit_plots = fCanvas3;

    fAngClusterFitAna_plots = fCanvas4;

    if(kIntReco) {
      fCanvas6->cd();      
      fCanvas6->SetName("fCanvas6");
      fCanvas6->Divide(2,2);
      fCanvas6->SetFillColor(42);
      fCanvas6->cd(1);
      gPad->Range(0,0,1,1);
      fCanvas6->cd(2);
      gPad->Range(0,0,1,1);
      fCanvas6->cd(3);
      gPad->Range(0,0,1,1);
      fCanvas6->cd(4);
      gPad->Range(0,0,1,1);

      fCanvas6->cd(2);
      fIntRecoHistU = new TH1F("IntRecoHistU","Interactively Reconstructed Shower - Transverse View",1,0,1);
      fIntRecoHistU->SetXTitle("Transverse distance from event vertex (m)");
      fIntRecoHistU->SetYTitle("Summed PH (PE)");
      fIntRecoHistU->SetLineColor(4);
      fIntRecoHistU->SetStats(0);
      fIntRecoHistV = new TH1F("IntRecoHistV","Interactively Reconstructed Shower - Transverse View",1,0,1);
      fIntRecoHistV->SetXTitle("Transverse distance from event vertex (m)");
      fIntRecoHistV->SetYTitle("Summed PH (PE)");
      fIntRecoHistV->SetLineColor(6);
      fIntRecoHistV->SetStats(0);

      fPredictedHistU = new TH1F("PredictedHistU","Interactively Reconstructed Shower - Transverse View",1,0,1);
      fPredictedHistU->SetXTitle("Transverse distance from event vertex (m)");
      fPredictedHistU->SetYTitle("Summed PH (PE)");
      fPredictedHistU->SetLineColor(4);
      fPredictedHistU->SetStats(0); 
      fPredictedHistV = new TH1F("PredictedHistV","Interactively Reconstructed Shower - Transverse View",1,0,1);
      fPredictedHistV->SetXTitle("Transverse distance from event vertex (m)");
      fPredictedHistV->SetYTitle("Summed PH (PE)");
      fPredictedHistV->SetLineColor(6);
      fPredictedHistV->SetStats(0);

      fPredictedHistU->Draw();
      fPredictedHistV->Draw("sames");
      fIntRecoHistU->Draw("e1sames");
      fIntRecoHistV->Draw("e1sames");

      fCanvas6->cd(4);
      fIntRecoHistZ = new TH1F("IntRecoHistZ","Interactively Reconstructed Shower - Longitudinal View",1,0,1);
      fIntRecoHistZ->SetXTitle("Longitudinal distance from event vertex (m)");
      fIntRecoHistZ->SetYTitle("Summed PH (PE)");
      fIntRecoHistZ->SetStats(0);

      fPredictedHistZ = new TH1F("PredictedHistZ","Interactively Reconstructed Shower - Longitudinal View",1,0,1);
      fPredictedHistZ->SetXTitle("Longitudinal distance from event vertex (m)");
      fPredictedHistZ->SetYTitle("Summed PH (PE)");
      fPredictedHistZ->SetLineColor(2);
      fPredictedHistZ->SetStats(0);
      fPredictedHistZ->Draw();
      fIntRecoHistZ->Draw("e1sames");

      fUZPred = new TH2F("PredictedHistUZ","EM Prediction Based on Interactive Reconstruction - UZ View",1,0,1,1,0,1);
      fUZPred->SetXTitle("Plane");
      fUZPred->SetYTitle("U Strip");
      fUZPred->SetZTitle("Normalised PH");
      fUZPred->SetStats(0);
      fVZPred = new TH2F("PredictedHistVZ","EM Prediction Based on Interactive Reconstruction - VZ View",1,0,1,1,0,1);
      fVZPred->SetXTitle("Plane");
      fVZPred->SetYTitle("V Strip");
      fVZPred->SetZTitle("Normalised PH");      
      fVZPred->SetStats(0);

      fCanvas6->cd(1);
      fUZPred->Draw("COLZ");
      fCanvas6->cd(3);
      fVZPred->Draw("COLZ");

      fCanvas7->SetName("fCanvas7");
      fCanvas7->cd();
      fCanvas7->SetFillColor(42);
      fSelectPad1 = new SelectPad("fSelPad1","fSelPad1",
				  0.01,0.505,0.99,0.99);
      fSelectPad1->Draw();      
      fSelectPad2 = new SelectPad("fSelPad2","fSelPad2",
				  0.01,0.01,0.99,0.495);
      fSelectPad2->Draw();
      fSelectPad2->cd();
      fIntRecoDoSim = new TButton("Do Sim","TButton *but = (TButton*) fSelPad2->FindObject(\"IntRecoDoSimButton\");if(but->GetFillColor()==2) but->SetFillColor(3); else but->SetFillColor(2);but->Paint();",0.01,0.04,0.075,0.15);
      fIntRecoDoSim->SetName("IntRecoDoSimButton");
      fIntRecoDoSim->SetFillColor(2);
      fIntRecoDoSim->Draw();
    }

    fCanvas8->cd(1);
    TimeHst = new THStack("TimeHist",
			  "Digit Times (ns) Red - Trk, Blue - Shw");
    TimeHstTrk = new TH1F();
    TimeHstTrk->SetName("TimeHstTrk");
    TimeHstShw = new TH1F();
    TimeHstShw->SetName("TimeHstShw");

    TimeHst->Add(TimeHstTrk);
    TimeHst->Add(TimeHstShw);
    TimeHst->Draw("nostack");

    fCanvas8->cd(2);
    TimeHstUV = new THStack("TimeHistUV",
			    "UV Digit Times (ns) Red - Trk, Blue - Shw");
    TimeHstTrkU = new TH1F();
    TimeHstTrkU->SetName("TimeHstTrkU");
    TimeHstTrkV = new TH1F();
    TimeHstTrkV->SetName("TimeHstTrkV");
    TimeHstShwU = new TH1F();
    TimeHstShwU->SetName("TimeHstShwU");
    TimeHstShwV = new TH1F();
    TimeHstShwV->SetName("TimeHstShwV");

    TimeHstUV->Add(TimeHstTrkU);
    TimeHstUV->Add(TimeHstTrkV);
    TimeHstUV->Add(TimeHstShwU);
    TimeHstUV->Add(TimeHstShwV);
    TimeHstUV->Draw("nostack");
    
    fCanvas8->cd(3);
    TimeHst2 = new THStack("TimeHist2",
			   "Digit Times (ns) after prop. cor. Red - Trk, Blue - Shw");
    TimeHstTrk2 = new TH1F();
    TimeHstTrk2->SetName("TimeHstTrk2");
    TimeHstShw2 = new TH1F();
    TimeHstShw2->SetName("TimeHstShw2");
    
    TimeHst2->Add(TimeHstTrk2);
    TimeHst2->Add(TimeHstShw2);
    TimeHst2->Draw("nostack");

    TLatex *tmptextrk = new TLatex(0.6,0.85,"tmp");
    tmptextrk->SetName("tmptextrk");
    tmptextrk->SetTextColor(4);
    tmptextrk->Draw();

    fCanvas8->cd(4);
    TimeHst2UV = new THStack("TimeHist2UV",
			     "UV Digit Times (ns) after prop. cor. Red - Trk, Blue - Shw");
    TimeHstTrk2U = new TH1F();
    TimeHstTrk2U->SetName("TimeHstTrk2U");
    TimeHstTrk2V = new TH1F();
    TimeHstTrk2V->SetName("TimeHstTrk2V");
    TimeHstShw2U = new TH1F();
    TimeHstShw2U->SetName("TimeHstShw2U");
    TimeHstShw2V = new TH1F();
    TimeHstShw2V->SetName("TimeHstShw2V");

    TimeHst2UV->Add(TimeHstTrk2U);
    TimeHst2UV->Add(TimeHstTrk2V);
    TimeHst2UV->Add(TimeHstShw2U);
    TimeHst2UV->Add(TimeHstShw2V);
    TimeHst2UV->Draw("nostack");

  //////////////////////////////////////////////////////////////////

    fCanvas9->cd(2);

    for(int i =0; i < 6; i++) place[i] = new TBox(0,0,0,0);
    place[0]->SetFillColor(kMagenta);
    place[1]->SetFillColor(kRed);
    place[2]->SetFillColor(kBlack);
    place[3]->SetFillColor(kBlue);
    place[4]->SetFillColor(kCyan);
    place[5]->SetFillColor(kGreen);

    leg = new TLegend(0.8,0.49,0.99,0.99);
    leg->SetFillColor(0);
    leg->AddEntry(place[0], "Earlier", "f");
    leg->AddEntry(place[1], "Previous", "f");
    leg->AddEntry(place[2], "Current", "f");
    leg->AddEntry(place[3], "Next", "f");
    leg->AddEntry(place[4], "Later", "f");
    leg->Draw();

    fCanvas9->cd(3);

    char cmd1[100];  char cmd2[100];  char cmd3[100];
    sprintf(cmd1, "gROOT->ProcessLineSync(\"((NueDisplayModule*)%p)->AdjustOverlayRange(1);\");", (void*)this);
    sprintf(cmd2, "gROOT->ProcessLineSync(\"((NueDisplayModule*)%p)->AdjustOverlayRange(2);\");", (void*)this); 
    sprintf(cmd3, "gROOT->ProcessLineSync(\"((NueDisplayModule*)%p)->AdjustOverlayRange(-1);\");", (void*)this);
    
    fCanvas9->cd();
    fPlusMinusOne = new TButton("+-One",cmd1,0,0.01,0.04,0.05);
    fPlusMinusTwo = new TButton("+-Two",cmd2,0.,0.06,0.04,0.10);
    fFullSnarl = new TButton("All",cmd3,0,0.11,0.04,0.15);

    fPlusMinusOne->Draw();
    fPlusMinusTwo->Draw();
    fFullSnarl->Draw();
    ///////////////////////////////////////////////////////////

    //set up the buttons in MR pad1
    fMRCanvas1->cd();
    
    mrInfoPanel = new TPad("Info0","Basic information",0.0,0.0,0.149,1);
    mrInfoPanel->Draw();
    mrInfoPanel->cd();
   
    char mrcmd1[100];
    sprintf(mrcmd1, "gROOT->ProcessLineSync(\"((NueDisplayModule*)%p)->MRShowEvent(0);\");", (void*)this);
    fMRuShowAll = new TButton("New",mrcmd1,0,0.85,0.33,0.89);

    sprintf(mrcmd1, "gROOT->ProcessLineSync(\"((NueDisplayModule*)%p)->MRShowEvent(1);\");", (void*)this);
    fMRuShowMR  = new TButton("MR",mrcmd1,0.,0.90,0.33,0.94);

    sprintf(mrcmd1, "gROOT->ProcessLineSync(\"((NueDisplayModule*)%p)->MRShowEvent(2);\");", (void*)this);
    fMRuShowOld = new TButton("Orig",mrcmd1,0,0.95,0.33,0.99);

    fMRuShowAll->Draw();
    fMRuShowMR->Draw();
    fMRuShowOld->Draw();

    sprintf(mrcmd1, "gROOT->ProcessLineSync(\"((NueDisplayModule*)%p)->MRShowMRStrips(0);\");", (void*)this);
    fMRdShowAll = new TButton("All",mrcmd1, 0 ,0.21,0.33,0.25);
      
    sprintf(mrcmd1, "gROOT->ProcessLineSync(\"((NueDisplayModule*)%p)->MRShowMRStrips(1);\");", (void*)this);
    fMRdShowTrueMu = new TButton("TruMu",mrcmd1,0,0.16,0.33,0.20);

    sprintf(mrcmd1, "gROOT->ProcessLineSync(\"((NueDisplayModule*)%p)->MRShowMRStrips(2);\");", (void*)this);
    fMRdShowTrueShw = new TButton("TruShw",mrcmd1,0.0,0.11,0.33,0.15);

    sprintf(mrcmd1, "gROOT->ProcessLineSync(\"((NueDisplayModule*)%p)->MRShowMRStrips(3);\");", (void*)this);
    fMRdShowScaled = new TButton("Scaled",mrcmd1,0,0.06,0.33,0.10);

    sprintf(mrcmd1, "gROOT->ProcessLineSync(\"((NueDisplayModule*)%p)->MRShowMRStrips(4);\");", (void*)this);
    fMRdShowReco = new TButton("Reco",mrcmd1,0.0,0.01,0.33,0.05);

    fMRdShowAll->Draw();
    fMRdShowTrueMu->Draw();
    fMRdShowTrueShw->Draw();
    fMRdShowScaled->Draw();
    fMRdShowReco->Draw();

    mrInfoPanel->Draw();
    
    /* reusing the boxes from above
    place[0]->SetFillColor(kMagenta);
    place[1]->SetFillColor(kRed);
    place[2]->SetFillColor(kBlack);
    place[3]->SetFillColor(kBlue);
    place[4]->SetFillColor(kCyan);
    */

//    fMRCanvas1->cd(2);
    fMRUpperLeg = new TLegend(0.35,0.85,0.99,0.99);
    fMRUpperLeg->SetFillColor(0);
    fMRUpperLeg->AddEntry(place[3], "Original", "f");
    fMRUpperLeg->AddEntry(place[1], "MR Strips", "f");
    fMRUpperLeg->AddEntry(place[5], "New Event", "f");
    fMRUpperLeg->SetTextSize(0.1);
    fMRUpperLeg->SetBorderSize(1);
    fMRUpperLeg->Draw();

//    fMRCanvas1->cd(4);
    fMRLowerLeg = new TLegend(0.35,0.01,0.99,0.2);
    fMRLowerLeg->SetFillColor(0);
    fMRLowerLeg->AddEntry(place[1], "True Mu", "f");
    fMRLowerLeg->AddEntry(place[5], "True Shw", "f");
    fMRLowerLeg->AddEntry(place[3], "Scaled", "f");
    fMRLowerLeg->AddEntry(place[4], "Reco", "f");
    fMRLowerLeg->SetTextSize(0.1);
    fMRLowerLeg->SetBorderSize(1);
    fMRLowerLeg->Draw();
   
   ///////////////////////////////////////////////////

    fMRCanvas1->cd();
    mrGraphPad = new TPad("Info0","mrgraphs",0.15,0.0,1,1);
    mrGraphPad->Draw();
    mrGraphPad->Divide(2,2);
                                                                                
    mrGraphPad->Draw();


    fEvent = 0;    
}

static void clear_hist(TH2D* hist)
{
    hist->Reset();
    hist->GetXaxis()->UnZoom();
    hist->GetYaxis()->UnZoom();
    hist->SetMaximum(-1111);
    hist->SetMinimum(-1111);
}

void NueDisplayModule::UpdateDisplay(bool /*foundnr*/, bool /*foundpid*/)
{
  MSG("NueDisplayModule",Msg::kDebug)<<"In UpdateDisplay"<<endl;
  
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

  clear_hist(fUZview);
  clear_hist(fVZview);
  clear_hist(fTrkUZ);
  clear_hist(fTrkVZ);
  clear_hist(fShwUZ);
  clear_hist(fShwVZ);
  clear_hist(fSlcUZ);
  clear_hist(fSlcVZ);
  clear_hist(fUZcolz);
  clear_hist(fVZcolz);
  clear_hist(fUZlego);
  clear_hist(fVZlego);
  fUVview->Clear();
  fInfo0->Clear();
  ftrkshw->Reset();

  if (TimeHst){
    TimeHst->RecursiveRemove(TimeHstTrk);
    TimeHst->RecursiveRemove(TimeHstShw);
    delete TimeHst;
    TimeHst = NULL;
  }
  if(TimeHstUV){
    TimeHstUV->RecursiveRemove(TimeHstTrkU);
    TimeHstUV->RecursiveRemove(TimeHstTrkV);
    TimeHstUV->RecursiveRemove(TimeHstShwU);
    TimeHstUV->RecursiveRemove(TimeHstShwV);
    delete TimeHstUV;
    TimeHstUV = NULL;
  }

  TimeHstTrk->Reset();
  TimeHstTrk->SetLineColor(2);
  TimeHstTrkU->Reset();
  TimeHstTrkU->SetLineColor(2);
  TimeHstTrkU->SetFillColor(2);
  TimeHstTrkU->SetFillStyle(3004);
  TimeHstTrkV->Reset();
  TimeHstTrkV->SetLineColor(46);
  TimeHstTrkV->SetFillColor(46);
  TimeHstTrkV->SetFillStyle(3005);
  TimeHstShw->Reset();
  TimeHstShw->SetLineColor(4);
  TimeHstShwU->Reset();
  TimeHstShwU->SetLineColor(4);
  TimeHstShwU->SetFillColor(4);
  TimeHstShwU->SetFillStyle(3004);
  TimeHstShwV->Reset();
  TimeHstShwV->SetLineColor(38);
  TimeHstShwV->SetFillColor(38);
  TimeHstShwV->SetFillStyle(3005);

  if (TimeHst2){
    TimeHst2->RecursiveRemove(TimeHstTrk2);
    TimeHst2->RecursiveRemove(TimeHstShw2);
    delete TimeHst2;
    TimeHst2 = NULL;
  }
  if (TimeHst2UV){
    TimeHst2UV->RecursiveRemove(TimeHstTrk2U);
    TimeHst2UV->RecursiveRemove(TimeHstTrk2V);
    TimeHst2UV->RecursiveRemove(TimeHstShw2U);
    TimeHst2UV->RecursiveRemove(TimeHstShw2V);
    delete TimeHst2UV;
    TimeHst2UV = NULL;
  }

  TimeHstTrk2->Reset();
  TimeHstTrk2->SetLineColor(2);
  TimeHstTrk2U->Reset();
  TimeHstTrk2U->SetLineColor(2);
  TimeHstTrk2U->SetFillColor(2);
  TimeHstTrk2U->SetFillStyle(3004);
  TimeHstTrk2V->Reset();
  TimeHstTrk2V->SetLineColor(46);
  TimeHstTrk2V->SetFillColor(46);
  TimeHstTrk2V->SetFillStyle(3005);
  TimeHstShw2->Reset();
  TimeHstShw2->SetLineColor(4);
  TimeHstShw2U->Reset();
  TimeHstShw2U->SetLineColor(4);
  TimeHstShw2U->SetFillColor(4);
  TimeHstShw2U->SetFillStyle(3004);
  TimeHstShw2V->Reset();
  TimeHstShw2V->SetLineColor(38);
  TimeHstShw2V->SetFillColor(38);
  TimeHstShw2V->SetFillStyle(3005);

  if (gr_dtds){
    delete gr_dtds;
    gr_dtds = NULL;
  }

  if(kDrawClu){
    delete cluLegU; cluLegU = NULL;
    delete cluLegV; cluLegV = NULL;
    delete ssGraphU; ssGraphU = NULL;
    delete ssGraphV; ssGraphV = NULL;
  }

  if(kIntReco){
    for(UInt_t i=0;i<fStripButtonU.size();i++){
      delete fStripButtonU[i];
    }
    fStripButtonU.clear();
    fStpIndexMapU.clear();
    for(UInt_t i=0;i<fStripButtonV.size();i++){
      delete fStripButtonV[i];
    }
    fStripButtonV.clear();
    fStpIndexMapV.clear();    
    fIntRecoHistZ->Reset();
    fIntRecoHistU->Reset();
    fIntRecoHistV->Reset();
    fPredictedHistZ->Reset();
    fPredictedHistV->Reset();
    fPredictedHistU->Reset();
    fUZPred->Reset();
    fVZPred->Reset();
    for(int i=1;i<=4;i++){
      fCanvas6->cd(i);
      gPad->Modified();
    }
    fCanvas6->Update();
    fSelectPad1->cd(); gPad->Modified();
    fSelectPad2->cd(); gPad->Modified();
    fCanvas7->Update();
  }

  fUZcolz->SetMinimum(0);
  fVZcolz->SetMinimum(0);

  for (int i = 0; i<4; i++){
    if (!i) cbut[i]->SetFillColor(4);
    else cbut[i]->SetFillColor(5);
    cbut[i]->Paint();
  }
    
  for (int i = 0; i<4; i++){
    if (!i) lbut[i]->SetFillColor(4);
    else lbut[i]->SetFillColor(5);
    lbut[i]->Paint();
  }

  fHistlego->cd(1);
  gPad->SetPhi(30);
  gPad->SetTheta(30);
  gPad->Modified();
  fHistlego->cd(2);
  gPad->SetPhi(30);
  gPad->SetTheta(30);
  gPad->Modified();
    
  if (imctruth||ifixmcinfo){
    imctruth = 0;
    fMCTruth->SetDown(false);
    this->delmc();
  }
  if (foundST) {
    fNumEvents = st->evthdr.nevent;
  }
  else if (foundSR) {
    fNumEvents = sr->evthdr.nevent;
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
    if (foundST){
      MSG("NueDisplayModule",Msg::kError)<<"Couldn't get event "<<fEvent
					 <<" from Snarl "<<st->GetHeader().GetSnarl()<<endl;
      if (clickbutton == -1) gMint->Prev();
      if (clickbutton ==  1) gMint->Next();
      return;
    }
    else if(foundSR){
      MSG("NueDisplayModule",Msg::kError)<<"Couldn't get event "<<fEvent
					 <<" from Snarl "<<sr->GetHeader().GetSnarl()<<endl;
      if (clickbutton == -1) gMint->Prev();
      if (clickbutton ==  1) gMint->Next();
      return;
    } 
  }
  
  //Determine if event passes the general cuts  
  bool passcuts = PassCuts();
  // Determine if event passes pid cuts
  bool passpidcuts = true;
  if (foundpidmatch){
//    Int_t test=pid->IsNue; 
//    MSG("NueDisplayModule",Msg::kDebug)<<" PID: "<< test << " kPIDCut: "<<kPIDCut<<endl;
//    if (test<kPIDCut) passpidcuts = false;
  }
  // Determine if event passes shwfit cut (aka hiPh cut) 
  bool passShwfitCut = true;

  if(kLoPhNStripCut>0||kLoPhNPlaneCut>0){
    if (foundST){
      sfa.Analyze(fEvent,st);
    }
    else if (foundSR){
      sfa.Analyze(fEvent,sr);
    }    

    passShwfitCut = sfa.PassCuts(kLoPhNStripCut,kLoPhNPlaneCut);
    if(passShwfitCut) cout << " Passes Shwfit Count Cut " << endl; 
  }

  if (((passcuts&&passShwfitCut&&preselec==1)||!preselec)&&passpidcuts){ //passcuts
    //  if (((passcuts&&preselec==1)||!preselec)){

    //reseting voting stuff
    if (!iIO) {
      if (kScanMode&&clickbutton&&fSnarl>0){
	if (!iFileW){
	  this->OpenFileWrite();
	}
	if (!hitlog){//if 'Log Details' was not clicked
	  outfile<<RunNo_old<<" "<<SubRunNo_old<<" "<<fSnarl<<" "<<fEvent_old<<" "<<ievtp<<" "<<itopo<<"  "<<fComment->GetText()<<endl;
	  fComment->Clear();
	  if (RunNo_old!=RunNo || SubRunNo_old!=SubRunNo){
	    RunNo_old = RunNo;
	    SubRunNo_old = SubRunNo;
	  }
	}
      }
      ievtp = 0;
      itopo = 0;
    }
    updateEncoded();
    if (!icomm && hitlog) outfile<<endl;
    hitlog = 0;
    icomm = 0;
    fEvent_old = fEvent;
    Analyze(fEvent);

    GetBasicInfo();

    if(kDrawClu){
      fCanvas5->cd(1);
      fUZcolz->Draw("COLZ");
      fCanvas5->cd(4);
      fVZcolz->Draw("COLZ");
      fCanvas5->cd(2);
      if(ssGraphU) ssGraphU->Draw("AP");
      fCanvas5->cd(3);
      if(cluLegU) cluLegU->Draw();
      fCanvas5->cd(5);
      if(ssGraphV) ssGraphV->Draw("AP");
      fCanvas5->cd(6);
      if(cluLegV) cluLegV->Draw();    
      fCanvas5->Update();
      fCanvas5->Modified();
    }

    fHistPad->cd(1);
    gPad->Modified();
    fHistPad->cd(2);
    gPad->Modified();
    fUVview->cd();
    fSteelOutline->Draw();
    if (fDetectorType == Detector::kFar) {
      VS->Draw();
      coil->Draw();
    }
    if (fDetectorType == Detector::kNear){
      fu1_outline->Draw();
      fu2_outline->Draw();
      fv1_outline->Draw();
      fv2_outline->Draw();
      pu1_outline->Draw();
      //pu2_outline->Draw();
      pv1_outline->Draw();
      //pv2_outline->Draw();
    }
    ftrkshw->SetMarkerStyle(20);
    ftrkshw->SetMarkerColor(2);
    ftrkshw->SetMarkerSize(0.4);
    ftrkshw->Draw("y:x","type>0","same");
    ftrkshw->SetMarkerColor(3);
    ftrkshw->Draw("y:x","type<0","same");
    srvtx_xy->Draw("same");
    gPad->Modified();

    fInfo0->cd();
    info1->Draw();
    info2->Draw();
    info3->Draw();
    info4->Draw();
    info41->Draw();
    info5->Draw();
    //info6->Draw();
    gPad->Modified();
    fCanvas0->Update();

    fHistcolz->cd(1);
    gPad->Modified();
    fHistcolz->cd(2);
    gPad->Modified();
    fHistlego->cd(1);
    gPad->Modified();
    fHistlego->cd(2);
    gPad->Modified();
    fCanvas1->Update();

    fEventEntry->SetText(Form("%d",fEvent));
    fEventNo->SetText(Form("Event: %d(%d)",fEvent,fNumEvents));

    
    if (ifixmcinfo) {this->plotmc();}

    fReco_plots->cd(1);
    fSlcUZ->GetXaxis()->SetRangeUser(lowest_z,highest_z);
    fSlcUZ->GetYaxis()->SetRangeUser(lowest_t0,highest_t0);
    gPad->Modified();
    fReco_plots->cd(2);
    fSlcVZ->GetXaxis()->SetRangeUser(lowest_z,highest_z);
    fSlcVZ->GetYaxis()->SetRangeUser(lowest_t1,highest_t1);
    gPad->Modified();
    fCanvas2->Update();

    fCanvas8->cd(1);
    TimeHst->Draw("nostack");    
    gPad->Modified();
    fCanvas8->cd(2);
    TimeHstUV->Draw("nostack");
    gPad->Modified();

    fCanvas8->cd(3);
    Double_t phFrac = TimeHstTrk2->Integral();
    if(phFrac>0) phFrac = TimeHstTrk2->GetMaximum()/phFrac;
    cout << "Track Timing Peak Fraction = " << phFrac << endl;
    phFrac = TimeHstShw2->Integral();
    if(phFrac>0) phFrac = TimeHstShw2->GetMaximum()/phFrac;
    cout << "Shower Timing Peak Fraction = " << phFrac << endl;
    TimeHst2->Draw("nostack");
    gPad->Modified();

    fCanvas8->cd(4);
    TimeHst2UV->Draw("nostack");
    gPad->Modified();

    if(false){
      fCanvas8->cd(2);
      if (gr_dtds){
	gPad->Clear();
	gr_dtds->SetMarkerStyle(20);
	gr_dtds->SetMarkerSize(0.5);
	gr_dtds->Draw("ap");
	tfit_dt_ds_pos->SetLineWidth(2);
	tfit_dt_ds_pos->SetLineColor(2);
	tfit_dt_ds_pos->Draw("same");
	tfit_dt_ds_neg->SetLineWidth(2);
	tfit_dt_ds_neg->SetLineColor(4);
	tfit_dt_ds_neg->Draw("same");
	gPad->Modified();
      }
      else{
	gPad->Clear();
      }
    }
    fCanvas8->Update();
    
    MSG("NueDisplayModule",Msg::kDebug)<<"Alldone with update"<<endl;
  }
  else{ 
    
    MSG("NueDisplayModule",Msg::kInfo)<<"Event does not pass cuts" << endl;

    info1->Clear();
    char text1[100];
    sprintf(text1,"Event did not pass cuts");
    info1->SetText(0.01,0.95,text1);
    fInfo0->cd();
    info1->Draw();
    gPad->Modified();
    fCanvas0->Update();

    if (clickbutton == -1) PrevEvent();
    if (clickbutton ==  1) NextEvent();
    
    return;
  }
  clickbutton = 0;

  fCanvas9->cd();
  UpdateOverlayGraphs(fEvent);

  fCanvas9->cd(1);
  if(uzSliceOverlay){
    uzSliceOverlay->Draw("AP");
    uzSliceOverlay->GetXaxis()->SetTitle("z position (m)");
    uzSliceOverlay->GetYaxis()->SetTitle("u Position (m)");
  }
  gPad->Modified();
  fCanvas9->cd(2);
  if(vzSliceOverlay){
    vzSliceOverlay->Draw("AP");
    vzSliceOverlay->GetXaxis()->SetTitle("z position (m)");
    vzSliceOverlay->GetYaxis()->SetTitle("v Position (m)");
  }
  leg->Draw();

  gPad->Modified();  fCanvas9->cd(3);
  if(uzEventOverlay){
    uzEventOverlay->Draw("AP");
    uzEventOverlay->GetXaxis()->SetTitle("z position (m)");
    uzEventOverlay->GetYaxis()->SetTitle("u Position (m)");
  }
  gPad->Modified(); fCanvas9->cd(4);
  if(vzEventOverlay){
    vzEventOverlay->Draw("AP");
    vzEventOverlay->GetXaxis()->SetTitle("z position (m)");
    vzEventOverlay->GetYaxis()->SetTitle("v Position (m)");
  }
  gPad->Modified();

  fCanvas9->Update();
  fCanvas9->Modified();

  if(foundMR){
    //fMRCanvas1->cd();
   int response = UpdateMRGraphs(fEvent);

   DrawUpperMRCanvas();
   if(response > 0){
      DrawLowerMRCanvas();
   }else
     {
        mrInfoPanel->cd();
        mrInfo1->SetText(0.01,0.80,"There was no MR Match");
        mrInfo1->SetTextSize(0.06);
        mrInfo1->SetTextColor(kBlack);
        mrInfo1->Draw();
     }

   /*
     fMRCanvas1->cd();
     fMRuShowAll->Draw();
     fMRuShowMR->Draw();
     fMRuShowOld->Draw();
   */                                                                                                               
    fMRCanvas1->Update();
    fMRCanvas1->Modified();

  }
  else
  {
     mrInfoPanel->cd();
     mrInfo1->SetText(0.01,0.80,"This Canvas is for MR Record Display");
     mrInfo1->SetText(0.01,0.74," No such records found");

     mrInfo1->SetTextSize(0.06);
     mrInfo1->SetTextColor(kBlack);
     mrInfo1->Draw();
  }

}

Int_t NueDisplayModule::GetEvent(Int_t evt){

  // Get Ugli for later.  Bail immediately if fail
  UgliGeomHandle ugh = gMint->GetUgliGeomHandle();
  if (! ugh.IsValid()) {
    MSG("NueDisplayModule",Msg::kWarning) << "Got invalid Ugli\n";
  }

  if (foundST) {
    event = SntpHelpers::GetEvent(evt,st);
  }
  else if(foundSR){
    event = SntpHelpers::GetEvent(evt,sr);
  }

  if(!event) return 0;
 
  
  for(int k=0; k < PHSIZE; k++) { ph0[k] = ph1[k] = -1;}
  SntpHelpers::FillEventEnergy(ph0, ph1, evt, st, PHSIZE);

  shia.SetEventEnergyArray(ph0, ph1);
  fva.SetEventEnergyArray(ph0, ph1);
  sfa.SetEventEnergyArray(ph0, ph1);
  hca.SetEventEnergyArray(ph0, ph1);
  aca.SetEventEnergyArray(ph0, ph1);
  acfa.SetEventEnergyArray(ph0, ph1);
  

  fUZcolz->SetBins(60,event->vtx.z-0.3,event->vtx.z+1.,
		     30,event->vtx.u-0.618,event->vtx.u+0.618);
  fVZcolz->SetBins(60,event->vtx.z-0.3,event->vtx.z+1.,
		   30,event->vtx.v-0.618,event->vtx.v+0.618);
  
  fUZlego->SetBins(100,event->vtx.z-0.3,event->vtx.z+1.,
		   30,event->vtx.u-0.618,event->vtx.u+0.618);
  fVZlego->SetBins(100,event->vtx.z-0.3,event->vtx.z+1.,
		   30,event->vtx.v-0.618,event->vtx.v+0.618);
  
  
  if (foundST){//ST
    float evthighest_plane = 0;
    float evtlowest_plane = 500;
    float evthighest_strip0 = 0;
    float evtlowest_strip0 = 192;
    float evthighest_strip1 = 0;
    float evtlowest_strip1 = 192;
    
    evthighest_z = 0;
    evtlowest_z = 30.;
    evthighest_t0 = -4.0;
    evtlowest_t0 = 4.0;
    evthighest_t1 = -4.0;
    evtlowest_t1 = 4.0;    

    //loop over strips to fill histograms
    MSG("NueDisplayModule",Msg::kDebug)<<"Looping over strips "<<event->nstrip<<endl;
    for(int i=0;i<event->nstrip;i++){
      Int_t index = SntpHelpers::GetStripIndex(i,event);
      NtpSRStrip *strip = SntpHelpers::GetStrip(index,st);
      if(!strip){
	MSG("NueDisplayModule",Msg::kError)<<"Couldn't get strip "<<index<<" from event 0"
					   <<" in snarl "<<st->GetHeader().GetSnarl()
					   <<" something has gone horribly wrong, I'll just go"
					   <<" on to the next event"<<endl;
	return 0;
      }
      else{
//	MSG("NueDisplayModule",Msg::kDebug)<<"got strip"<<endl;
      }


	if(strip->ph0.pe+strip->ph1.pe<1e-5)continue;//empty strip

      int tempo_pln = strip->plane;
      int tempo_stp = strip->strip;
      float tempo_tpos = strip->tpos;
      if(tempo_pln<evtlowest_plane) {
	evtlowest_plane=tempo_pln;
	evtlowest_z=strip->z;
      }
      if(tempo_pln>evthighest_plane) {
	evthighest_plane=tempo_pln;
	evthighest_z=strip->z;
      }
      
      if(strip->planeview==PlaneView::kU){
//	MSG("NueDisplayModule",Msg::kDebug)<<"filling u"<<endl;
	fUZview->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					    strip->ph1.sigcor)/SIGMAPMEU);
	if (strip->z>=event->vtx.z-0.3&&strip->z<=event->vtx.z+1.&&
	    strip->tpos>=event->vtx.u-0.618&&strip->tpos<=event->vtx.u+0.618){
	  fUZcolz->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					      strip->ph1.sigcor)/SIGMAPMEU);
	  fUZlego->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					      strip->ph1.sigcor)/SIGMAPMEU);
	  fStpIndexMapU[fUZcolz->FindBin(strip->z,strip->tpos)] = strip->index;
	}
	if(tempo_tpos<evtlowest_t0) {
	  evtlowest_strip0=tempo_stp;
	  evtlowest_t0=tempo_tpos;
	}
	if(tempo_tpos>evthighest_t0) {
	  evthighest_strip0=tempo_stp;
	  evthighest_t0=tempo_tpos;
	}
      }
      else if(strip->planeview==PlaneView::kV){
//	MSG("NueDisplayModule",Msg::kDebug)<<"filling v"<<endl;
	fVZview->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					    strip->ph1.sigcor)/SIGMAPMEU);
	if (strip->z>=event->vtx.z-0.3&&strip->z<=event->vtx.z+1.&&
	    strip->tpos>=event->vtx.v-0.618&&strip->tpos<=event->vtx.v+0.618){
	  fVZcolz->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					      strip->ph1.sigcor)/SIGMAPMEU);
	  fVZlego->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					      strip->ph1.sigcor)/SIGMAPMEU);
	  fStpIndexMapV[fVZcolz->FindBin(strip->z,strip->tpos)] = strip->index;
	}
	if(tempo_tpos<evtlowest_t1) {
	  evtlowest_strip1=tempo_stp;
	  evtlowest_t1=tempo_tpos;
	}
	if(tempo_tpos>evthighest_t1) {
	  evthighest_strip1=tempo_stp;
	  evthighest_t1=tempo_tpos;
	}
      }
    }
     
    if(evtlowest_plane-5>=0) {
      evtlowest_plane-=5;
      evtlowest_z-=5.*0.06;
    }
    else {
      evtlowest_plane=0;
      evtlowest_z=0.;
    }

    if (fDetectorType == Detector::kNear){
      if (evtlowest_t0 - 3.*0.041>=-4.0){
	evtlowest_strip0-=3;
	evtlowest_t0-=3.*0.041;
      }
      else {
	evtlowest_strip0=0;
	evtlowest_t0=-4.0;
      }
    }
    else {
      if (evtlowest_strip0-3>=0){
	evtlowest_strip0-=3;
	evtlowest_t0-=3.*0.041;
      }
      else {
	evtlowest_strip0=0;
	evtlowest_t0=-4.0;
      }
    }

    if (fDetectorType == Detector::kNear){
      if (evtlowest_t1 - 3.*0.041>=-4.0){
	evtlowest_strip1-=3;
	evtlowest_t1-=3.*0.041;
      }
      else {
	evtlowest_strip1=0;
	evtlowest_t1=-4.0;
      }
    }
    else {
      if (evtlowest_strip1-3>=0){
	evtlowest_strip1-=3;
	evtlowest_t1-=3.*0.041;
      }
      else {
	evtlowest_strip1=0;
	evtlowest_t1=-4.0;
      }
    }

    
    if(evthighest_plane+5<=485) {
      evthighest_plane+=5;
      evthighest_z+=5.*0.06;
    }
    else {
      evthighest_plane=485;
      evthighest_z=30.;
    }
    
    if(evthighest_strip0+3<=191) {
      evthighest_strip0+=3;
      evthighest_t0+=3.*0.041;
    }
    else {
      evthighest_strip0=191;
      evthighest_t0=4.0;
    }
      
    if(evthighest_strip1+3<=191) {
      evthighest_strip1+=3;
      evthighest_t1+=3.*0.041;
    }
    else {
      evthighest_strip1=191;
      evthighest_t1=4.0;
    }   

    //I implant the 6 parameters in the option string so they can be retrived interactively later. This is a risky trick. I hope no one will actually make use of that "options".TY
    char setoptions[100];
    sprintf(setoptions,"%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d",evthighest_z,evtlowest_z,evthighest_t0,evtlowest_t0,evthighest_t1,evtlowest_t1,0);
    fUZview->SetOption(setoptions);

    Double_t histmax = TMath::Max(fUZview->GetMaximum(),fVZview->GetMaximum());
    cout<<"histmax "<<fUZview->GetMaximum()<<" "<<fVZview->GetMaximum()<<" "<<histmax<<endl;
    fUZview->SetMaximum(1.05*histmax);
    fVZview->SetMaximum(1.05*histmax);
    Double_t colzmax = fUZcolz->GetMaximum()>fVZcolz->GetMaximum()?fUZcolz->GetMaximum():fVZcolz->GetMaximum();
    fUZcolz->SetMaximum(1.05*colzmax);
    fVZcolz->SetMaximum(1.05*colzmax);
    Double_t legomax = fUZlego->GetMaximum()>fVZlego->GetMaximum()?fUZlego->GetMaximum():fVZlego->GetMaximum();
    fUZlego->SetMaximum(1.05*legomax);
    fVZlego->SetMaximum(1.05*legomax);
    
    //record track hits
    int ntrks = event->ntrack;
    if (ntrks){
      for (int i = 0; i<ntrks; i++){
	int index = SntpHelpers::GetTrackIndex(i,event);
	NtpSRTrack *track = SntpHelpers::GetTrack(index,st);
	for (int istp = 0; istp<track->nstrip; istp++){
	  float x = track->stpx[istp];
	  float y = track->stpy[istp];
	  ftrkshw->Fill(x,y,1);
	}
      }
    }

    //record veto shield hits
    if (st->vetohdr.ishit){//projection is within boundaries of shield plane
      TClonesArray* shieldstparray = st->vetostp; 
      for (int i = 0; i<shieldstparray->GetEntries(); i++){
	NtpSRShieldStrip* shieldstp = dynamic_cast<NtpSRShieldStrip *>(shieldstparray->At(i));
	ftrkshw->Fill(shieldstp->x,shieldstp->y,-1);
      }
    }

    //timing histogram
    //find range
    double tmin=0, tmax=0;
    bool first = true;
    NtpSRSlice *slice = SntpHelpers::GetSlice(event->slc,st);

    //shamelessly stolen from Mad
    if (slice){//slice
      float highest_plane = 0;
      float lowest_plane = 500;
      float highest_strip0 = 0;
      float lowest_strip0 = 192;
      float highest_strip1 = 0;
      float lowest_strip1 = 192;
    
      highest_z = 0;
      lowest_z = 30.;
      highest_t0 = -4.0;
      lowest_t0 = 4.0;
      highest_t1 = -4.0;
      lowest_t1 = 4.0;    

      for (int i = 0; i<slice->nstrip; i++){//loop over strips
	NtpSRStrip *strip = SntpHelpers::GetStrip(slice->stp[i],st);	
	double t = strip->time0;
	if (t>-100){
	  if (first){
	    tmin = tmax = t;
	    first = false;
	  }
	  if (t < tmin) tmin = t;
	  if (t > tmax) tmax = t;
	}
	t = strip->time1;
	if (t>-100){
	  if (first){
	    tmin = tmax = t;
	    first = false;
	  }
	  if (t < tmin) tmin = t;
	  if (t > tmax) tmax = t;
	}
	int tempo_pln = strip->plane;
	int tempo_stp = strip->strip;
	float tempo_tpos = strip->tpos;
	if(tempo_pln<lowest_plane) {
	  lowest_plane=tempo_pln;
	  lowest_z=strip->z;
	}
	if(tempo_pln>highest_plane) {
	  highest_plane=tempo_pln;
	  highest_z=strip->z;
	}

	if (strip->planeview==PlaneView::kU){
	  fSlcUZ->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					     strip->ph1.sigcor)/SIGMAPMEU);
	  if(tempo_tpos<lowest_t0) {
	    lowest_strip0=tempo_stp;
	    lowest_t0=tempo_tpos;
	  }
	  if(tempo_tpos>highest_t0) {
	    highest_strip0=tempo_stp;
	    highest_t0=tempo_tpos;
	  }
	}
	else if (strip->planeview==PlaneView::kV){
	  fSlcVZ->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					     strip->ph1.sigcor)/SIGMAPMEU);
	  if(tempo_tpos<lowest_t1) {
	    lowest_strip1=tempo_stp;
	    lowest_t1=tempo_tpos;
	  }
	  if(tempo_tpos>highest_t1) {
	    highest_strip1=tempo_stp;
	    highest_t1=tempo_tpos;
	  }
	}
      }
      if(lowest_plane-10>=0) {
	lowest_plane-=10;
	lowest_z-=10.*0.06;
      }
      else {
	lowest_plane=0;
	lowest_z=0.;
      }
  
      if(lowest_strip0-5>=0) {
	lowest_strip0-=5;
	lowest_t0-=5.*0.041;
      }
      else {
	lowest_strip0=0;
	lowest_t0=-4.0;
      }
      
      if(lowest_strip1-5>=0) {
	lowest_strip1-=5;
	lowest_t1-=5.*0.041;
      }
      else {
	lowest_strip1=0;
	lowest_t1=-4.0;
      }
      
      if(highest_plane+10<=485) {
	highest_plane+=10;
	highest_z+=10.*0.06;
      }
      else {
	highest_plane=485;
	highest_z=30.;
      }
      
      if(highest_strip0+5<=191) {
	highest_strip0+=5;
	highest_t0+=5.*0.041;
      }
      else {
	highest_strip0=191;
	highest_t0=4.0;
      }
      
      if(highest_strip1+5<=191) {
	highest_strip1+=5;
	highest_t1+=5.*0.041;
      }
      else {
	highest_strip1=191;
	highest_t1=4.0;
      }      
    }
    // give some buffer at either end...
    tmin -= 50e-9;
    tmax += 20e-9;

    double eps = 1.0e-8;
    if (tmin == tmax) { tmin -= eps; tmax += eps; }
    
    // for far detector, suppress display of pre-trigger time interval
    if (fDetectorType == Detector::kFar){
      if ((tmax - tmin)*1e9>1500) tmin = tmax - 1500e-9;
    }

    TimeHstTrk->SetBins(100,tmin/1e-9,tmax/1e-9);
    TimeHstTrkU->SetBins(100,tmin/1e-9,tmax/1e-9);
    TimeHstTrkV->SetBins(100,tmin/1e-9,tmax/1e-9);
    TimeHstShw->SetBins(100,tmin/1e-9,tmax/1e-9);
    TimeHstShwU->SetBins(100,tmin/1e-9,tmax/1e-9);
    TimeHstShwV->SetBins(100,tmin/1e-9,tmax/1e-9);

    TimeHstTrk2->SetBins(50,0,100);
    TimeHstTrk2U->SetBins(50,0,100);
    TimeHstTrk2V->SetBins(50,0,100);
    TimeHstShw2->SetBins(50,0,100);
    TimeHstShw2U->SetBins(50,0,100);
    TimeHstShw2V->SetBins(50,0,100);

    //record track hits
    //int ntrks = event->ntrack;
    if (ntrks){//if(ntrks)
      int trkidx = -1;
      int trkplanes = -1;
      for (int i = 0; i<ntrks; i++){
	int index = SntpHelpers::GetTrackIndex(i,event);
	NtpSRTrack *track = SntpHelpers::GetTrack(index,st);
	if (track->plane.n>trkplanes){
	  trkplanes = track->plane.n;
	  trkidx = index;
	}
	for (int istp = 0; istp<track->nstrip; istp++){
	  //ftrkshw->Fill(track->stpx[istp],track->stpy[istp],1);
	  NtpSRStrip *strip = SntpHelpers::GetStrip(track->stp[istp],st);
	  TimeHstTrk->Fill(track->stpt0[istp]/1e-9,strip->ph0.pe);
	  TimeHstTrk->Fill(track->stpt1[istp]/1e-9,strip->ph1.pe);
	  TimeHstTrk2->Fill((track->stpt0[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph0.pe);
	  TimeHstTrk2->Fill((track->stpt1[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph1.pe);
	  //record track strip information
	  if(strip->planeview==PlaneView::kU){
	    fTrkUZ->Fill(strip->z,strip->tpos,100000);//paranoia

	    TimeHstTrkU->Fill(track->stpt0[istp]/1e-9,strip->ph0.pe);
	    TimeHstTrkU->Fill(track->stpt1[istp]/1e-9,strip->ph1.pe);
	    TimeHstTrk2U->Fill((track->stpt0[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph0.pe);
	    TimeHstTrk2U->Fill((track->stpt1[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph1.pe);
	  }
	  else if(strip->planeview==PlaneView::kV){
	    fTrkVZ->Fill(strip->z,strip->tpos,100000);

	    TimeHstTrkV->Fill(track->stpt0[istp]/1e-9,strip->ph0.pe);
	    TimeHstTrkV->Fill(track->stpt1[istp]/1e-9,strip->ph1.pe);
	    TimeHstTrk2V->Fill((track->stpt0[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph0.pe);
	    TimeHstTrk2V->Fill((track->stpt1[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph1.pe);
	  }
	}
      }
      vector<double> spathLength;
      vector<double> st0;
      if (trkidx != -1 && fDetectorType == Detector::kNear){
	NtpSRTrack *track = SntpHelpers::GetTrack(trkidx,st);
	for (int istp = 0; istp<track->nstrip; istp++){
	  NtpSRStrip *strip = SntpHelpers::GetStrip(track->stp[istp],st);
	  spathLength.push_back(track->ds-track->stpds[istp]);
	  PlexStripEndId seid(Detector::kNear,strip->plane,strip->strip,StripEnd::kWest);
	  UgliStripHandle stripHandle = ugh.GetStripHandle(seid);
	  float halfLength = stripHandle.GetHalfLength();
	  const TVector3 ghitxyz(track->stpx[istp],track->stpy[istp],track->stpz[istp]);
	  TVector3 lhitxyz = stripHandle.GlobalToLocal(ghitxyz);
	  float fiberDist = (halfLength - lhitxyz.x() + stripHandle.ClearFiber(StripEnd::kWest) + stripHandle.WlsPigtail(StripEnd::kWest));
	  //using strip time, I don't understand track time TJ
	  //st0.push_back(track->stpt1[istp]-fiberDist/PropagationVelocity::Velocity());
	  st0.push_back(strip->time1-fiberDist/PropagationVelocity::Velocity());
	  //cout<<strip->plane<<" "<<strip->strip<<" "<<track->ds-track->stpds[istp]<<Form(" %.9f %.9f %.9f",track->stpt1[istp],strip->time1,strip->time1-fiberDist/PropagationVelocity::Velocity())<<endl;
	}
	gr_dtds = new TGraph(st0.size(),&spathLength[0],&st0[0]);
	tfit_dt_ds_pos -> SetParameter(0,0.0) ;
	gr_dtds->Fit("tfit_dt_ds_pos","QR") ;
	tfit_dt_ds_neg -> SetParameter(0,0.0) ;
	gr_dtds->Fit("tfit_dt_ds_neg","QR") ;

	float trms1 = 0;
	float trms2 = 0;

	for (unsigned i = 0; i<st0.size(); i++){
	  double x,y;
	  gr_dtds->GetPoint(i,x,y);
	  trms1 += pow(y-tfit_dt_ds_pos->Eval(x),2);
	  trms2 += pow(y-tfit_dt_ds_neg->Eval(x),2);
	}
	trms1/=st0.size();
	trms2/=st0.size();
	trms1=sqrt(trms1)*1e9;
	trms2=sqrt(trms2)*1e9;
	char grtitle[100];
	sprintf(grtitle,"dt vs ds, rms+:%.1fns rms-:%.1fns (rms+)-(rms-):%.1fns",trms1,trms2,trms1-trms2);
	gr_dtds->SetTitle(grtitle);
	gr_dtds->SetName("gr_dtds");
      }
    }
    
    //record shower hits
    int nshws = event->nshower;
    if (nshws){
      for (int i = 0; i<nshws; i++){
	int index = SntpHelpers::GetShowerIndex(i,event);
	NtpSRShower *shower = SntpHelpers::GetShower(index,st);
	for (int istp = 0; istp<shower->nstrip; istp++){
	  NtpSRStrip *strip = SntpHelpers::GetStrip(shower->stp[istp],st);
	  if(ReleaseType::IsCedar(fRel)||ReleaseType::IsDogwood(fRel)){
	    TimeHstShw->Fill(shower->stpt0[istp]/1e-9,strip->ph0.pe);
	    TimeHstShw->Fill(shower->stpt1[istp]/1e-9,strip->ph1.pe);
	    TimeHstShw2->Fill((shower->stpt0[istp]-tmin)/1e-9,strip->ph0.pe);
	    TimeHstShw2->Fill((shower->stpt1[istp]-tmin)/1e-9,strip->ph1.pe);
	  }
	  else {
	    TimeHstShw->Fill(strip->time0/1e-9,strip->ph0.pe);
	    TimeHstShw->Fill(strip->time1/1e-9,strip->ph1.pe);
	    TimeHstShw2->Fill((strip->time0-tmin)/1e-9,strip->ph0.pe);
	    TimeHstShw2->Fill((strip->time1-tmin)/1e-9,strip->ph1.pe);
	  }
	  //record shower strip information
	  if(strip->planeview==PlaneView::kU){
	    fShwUZ->Fill(strip->z,strip->tpos,100000);//paranoia
	    if(ReleaseType::IsCedar(fRel)||ReleaseType::IsDogwood(fRel)){
	      TimeHstShwU->Fill(shower->stpt0[istp]/1e-9,strip->ph0.pe);
	      TimeHstShwU->Fill(shower->stpt1[istp]/1e-9,strip->ph1.pe);
	      TimeHstShw2U->Fill((shower->stpt0[istp]-tmin)/1e-9,strip->ph0.pe);
	      TimeHstShw2U->Fill((shower->stpt1[istp]-tmin)/1e-9,strip->ph1.pe);
	    }
	    else {
	      TimeHstShwU->Fill(strip->time0/1e-9,strip->ph0.pe);
	      TimeHstShwU->Fill(strip->time1/1e-9,strip->ph1.pe);
	      TimeHstShw2U->Fill((strip->time0-tmin)/1e-9,strip->ph0.pe);
	      TimeHstShw2U->Fill((strip->time1-tmin)/1e-9,strip->ph1.pe);
	    }
	  }
	  else if(strip->planeview==PlaneView::kV){
	    fShwVZ->Fill(strip->z,strip->tpos,100000);
	    if(ReleaseType::IsCedar(fRel)||ReleaseType::IsDogwood(fRel)){
	      TimeHstShwV->Fill(shower->stpt0[istp]/1e-9,strip->ph0.pe);
	      TimeHstShwV->Fill(shower->stpt1[istp]/1e-9,strip->ph1.pe);
	      TimeHstShw2V->Fill((shower->stpt0[istp]-tmin)/1e-9,strip->ph0.pe);
	      TimeHstShw2V->Fill((shower->stpt1[istp]-tmin)/1e-9,strip->ph1.pe);
	    }
	    else {
	      TimeHstShwV->Fill(strip->time0/1e-9,strip->ph0.pe);
	      TimeHstShwV->Fill(strip->time1/1e-9,strip->ph1.pe);
	      TimeHstShw2V->Fill((strip->time0-tmin)/1e-9,strip->ph0.pe);
	      TimeHstShw2V->Fill((strip->time1-tmin)/1e-9,strip->ph1.pe);
	    }
	  }
	}
      }
    }
    TimeHst = new THStack("TimeHist",
			  "Digit Times (ns) Red - Trk, Blue - Shw");
    TimeHst->Add(TimeHstTrk);
    TimeHst->Add(TimeHstShw);

    TimeHstUV = new THStack("TimeHistUV",
			    "UV Digit Times (ns) Red - Trk, Blue - Shw");
    TimeHstUV->Add(TimeHstTrkU);
    TimeHstUV->Add(TimeHstTrkV);
    TimeHstUV->Add(TimeHstShwU);
    TimeHstUV->Add(TimeHstShwV);

    TimeHst2 = new THStack("TimeHist2",
			   "Digit Times (ns) after prop. cor. Red - Trk, Blue - Shw");
    TimeHst2->Add(TimeHstTrk2);
    TimeHst2->Add(TimeHstShw2);

    TimeHst2UV = new THStack("TimeHist2UV",
			   "UV Digit Times (ns) after prop. cor. Red - Trk, Blue - Shw");
    TimeHst2UV->Add(TimeHstTrk2U);
    TimeHst2UV->Add(TimeHstTrk2V);
    TimeHst2UV->Add(TimeHstShw2U);
    TimeHst2UV->Add(TimeHstShw2V);
  }
  else if (foundSR){//SR
    float evthighest_plane = 0;
    float evtlowest_plane = 500;
    float evthighest_strip0 = 0;
    float evtlowest_strip0 = 192;
    float evthighest_strip1 = 0;
    float evtlowest_strip1 = 192;
    
    evthighest_z = 0;
    evtlowest_z = 30.;
    evthighest_t0 = -4.0;
    evtlowest_t0 = 4.0;
    evthighest_t1 = -4.0;
    evtlowest_t1 = 4.0;    
    
    //loop over strips to fill histograms
    MSG("NueDisplayModule",Msg::kDebug)<<"Looping over strips "<<event->nstrip<<endl;
    for(int i=0;i<event->nstrip;i++){
      Int_t index = SntpHelpers::GetStripIndex(i,event);
      NtpSRStrip *strip = SntpHelpers::GetStrip(index,sr);
      if(!strip){
	MSG("NueDisplayModule",Msg::kError)<<"Couldn't get strip "<<index<<" from event 0"
					   <<" in snarl "<<sr->GetHeader().GetSnarl()
					   <<" something has gone horribly wrong, I'll just go"
					   <<" on to the next event"<<endl;
	return 0;
      }
      else{
	MSG("NueDisplayModule",Msg::kDebug)<<"got strip"<<endl;
      }
      int tempo_pln = strip->plane;
      int tempo_stp = strip->strip;
      float tempo_tpos = strip->tpos;
      if(tempo_pln<evtlowest_plane) {
	evtlowest_plane=tempo_pln;
	evtlowest_z=strip->z;
      }
      if(tempo_pln>evthighest_plane) {
	evthighest_plane=tempo_pln;
	evthighest_z=strip->z;
      }
      
      if(strip->planeview==PlaneView::kU){
	MSG("NueDisplayModule",Msg::kDebug)<<"filling u"<<endl;
	fUZview->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					    strip->ph1.sigcor)/SIGMAPMEU);
	if (strip->z>=event->vtx.z-0.3&&strip->z<=event->vtx.z+1.&&
	    strip->tpos>=event->vtx.u-0.618&&strip->tpos<=event->vtx.u+0.618){
	  fUZcolz->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					      strip->ph1.sigcor)/SIGMAPMEU);
	  fUZlego->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					      strip->ph1.sigcor)/SIGMAPMEU);
	  fStpIndexMapU[fUZcolz->FindBin(strip->z,strip->tpos)] = strip->index;
	}
	if(tempo_tpos<evtlowest_t0) {
	  evtlowest_strip0=tempo_stp;
	  evtlowest_t0=tempo_tpos;
	}
	if(tempo_tpos>evthighest_t0) {
	  evthighest_strip0=tempo_stp;
	  evthighest_t0=tempo_tpos;
	}
      }
      else if(strip->planeview==PlaneView::kV){
	MSG("NueDisplayModule",Msg::kDebug)<<"filling v"<<endl;
	fVZview->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					    strip->ph1.sigcor)/SIGMAPMEU);
	if (strip->z>=event->vtx.z-0.3&&strip->z<=event->vtx.z+1.&&
	    strip->tpos>=event->vtx.v-0.618&&strip->tpos<=event->vtx.v+0.618){
	  fVZcolz->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					      strip->ph1.sigcor)/SIGMAPMEU);
	  fVZlego->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					      strip->ph1.sigcor)/SIGMAPMEU);
	  fStpIndexMapV[fVZcolz->FindBin(strip->z,strip->tpos)] = strip->index;
	}
	if(tempo_tpos<evtlowest_t1) {
	  evtlowest_strip1=tempo_stp;
	  evtlowest_t1=tempo_tpos;
	}
	if(tempo_tpos>evthighest_t1) {
	  evthighest_strip1=tempo_stp;
	  evthighest_t1=tempo_tpos;
	}
      }
    }
    if(evtlowest_plane-5>=0) {
      evtlowest_plane-=5;
      evtlowest_z-=5.*0.06;
    }
    else {
      evtlowest_plane=0;
      evtlowest_z=0.;
    }

    if (fDetectorType == Detector::kNear){
      if (evtlowest_t0 - 3.*0.041>=-4.0){
	evtlowest_strip0-=3;
	evtlowest_t0-=3.*0.041;
      }
      else {
	evtlowest_strip0=0;
	evtlowest_t0=-4.0;
      }
    }
    else {
      if (evtlowest_strip0-3>=0){
	evtlowest_strip0-=3;
	evtlowest_t0-=3.*0.041;
      }
      else {
	evtlowest_strip0=0;
	evtlowest_t0=-4.0;
      }
    }

    if (fDetectorType == Detector::kNear){
      if (evtlowest_t1 - 3.*0.041>=-4.0){
	evtlowest_strip1-=3;
	evtlowest_t1-=3.*0.041;
      }
      else {
	evtlowest_strip1=0;
	evtlowest_t1=-4.0;
      }
    }
    else {
      if (evtlowest_strip1-3>=0){
	evtlowest_strip1-=3;
	evtlowest_t1-=3.*0.041;
      }
      else {
	evtlowest_strip1=0;
	evtlowest_t1=-4.0;
      }
    }    
    
    if(evthighest_plane+5<=485) {
      evthighest_plane+=5;
      evthighest_z+=5.*0.06;
    }
    else {
      evthighest_plane=485;
      evthighest_z=30.;
    }
    
    if(evthighest_strip0+3<=191) {
      evthighest_strip0+=3;
      evthighest_t0+=3.*0.041;
    }
    else {
      evthighest_strip0=191;
      evthighest_t0=4.0;
    }
      
    if(evthighest_strip1+3<=191) {
      evthighest_strip1+=3;
      evthighest_t1+=3.*0.041;
    }
    else {
      evthighest_strip1=191;
      evthighest_t1=4.0;
    }       

    //I implant the 6 parameters in the option string so they can be retrived interactively later. This is a risky trick. I hope no one will actually make use of that "options".
    char setoptions[100];
    sprintf(setoptions,"%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d",evthighest_z,evtlowest_z,evthighest_t0,evtlowest_t0,evthighest_t1,evtlowest_t1,0);
    fUZview->SetOption(setoptions);

    //set histogram ranges
    Double_t histmax = TMath::Max(fUZview->GetMaximum(),fVZview->GetMaximum());
    fUZview->SetMaximum(1.05*histmax);
    fVZview->SetMaximum(1.05*histmax);
    Double_t colzmax = fUZcolz->GetMaximum()>fVZcolz->GetMaximum()?fUZcolz->GetMaximum():fVZcolz->GetMaximum();
    fUZcolz->SetMaximum(1.05*colzmax);
    fVZcolz->SetMaximum(1.05*colzmax);
    Double_t legomax = fUZlego->GetMaximum()>fVZlego->GetMaximum()?fUZlego->GetMaximum():fVZlego->GetMaximum();
    fUZlego->SetMaximum(1.05*legomax);
    fVZlego->SetMaximum(1.05*legomax);
    
    int ntrks = event->ntrack;
    if (ntrks){
      for (int i = 0; i<ntrks; i++){
	int index = SntpHelpers::GetTrackIndex(i,event);
	NtpSRTrack *track = SntpHelpers::GetTrack(index,sr);
	for (int istp = 0; istp<track->nstrip; istp++){
	  float x = track->stpx[istp];
	  float y = track->stpy[istp];
	  ftrkshw->Fill(x,y,1);
	}
      }
    }

    //record veto shield hits
    if (sr->vetohdr.ishit){//projection is within boundaries of shield plane 
      TClonesArray* shieldstparray = sr->vetostp; 
      for (int i = 0; i<shieldstparray->GetEntries(); i++){
	NtpSRShieldStrip* shieldstp = dynamic_cast<NtpSRShieldStrip *>(shieldstparray->At(i));
	ftrkshw->Fill(shieldstp->x,shieldstp->y,-1);
      }
    }
    //timing histogram
    //find range
    double tmin=0, tmax=0;
    bool first = true;
    NtpSRSlice *slice = SntpHelpers::GetSlice(event->slc,sr);
    if (slice){
      float highest_plane = 0;
      float lowest_plane = 500;
      float highest_strip0 = 0;
      float lowest_strip0 = 192;
      float highest_strip1 = 0;
      float lowest_strip1 = 192;
    
      highest_z = 0;
      lowest_z = 30.;
      highest_t0 = -4.0;
      lowest_t0 = 4.0;
      highest_t1 = -4.0;
      lowest_t1 = 4.0;    

      for (int i = 0; i<slice->nstrip; i++){
	NtpSRStrip *strip = SntpHelpers::GetStrip(slice->stp[i],sr);	
	double t = strip->time0;
	if (t>-100){
	  if (first){
	    tmin = tmax = t;
	    first = false;
	  }
	  if (t < tmin) tmin = t;
	  if (t > tmax) tmax = t;
	}
	t = strip->time1;
	if (t>-100){
	  if (first){
	    tmin = tmax = t;
	    first = false;
	  }
	  if (t < tmin) tmin = t;
	  if (t > tmax) tmax = t;
	}
	int tempo_pln = strip->plane;
	int tempo_stp = strip->strip;
	float tempo_tpos = strip->tpos;
	if(tempo_pln<lowest_plane) {
	  lowest_plane=tempo_pln;
	  lowest_z=strip->z;
	}
	if(tempo_pln>highest_plane) {
	  highest_plane=tempo_pln;
	  highest_z=strip->z;
	}

	if (strip->planeview==PlaneView::kU){
	  fSlcUZ->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					     strip->ph1.sigcor)/SIGMAPMEU);
	  if(tempo_tpos<lowest_t0) {
	    lowest_strip0=tempo_stp;
	    lowest_t0=tempo_tpos;
	  }
	  if(tempo_tpos>highest_t0) {
	    highest_strip0=tempo_stp;
	    highest_t0=tempo_tpos;
	  }
	}
	else if (strip->planeview==PlaneView::kV){
	  fSlcVZ->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + 
					     strip->ph1.sigcor)/SIGMAPMEU);
	  if(tempo_tpos<lowest_t1) {
	    lowest_strip1=tempo_stp;
	    lowest_t1=tempo_tpos;
	  }
	  if(tempo_tpos>highest_t1) {
	    highest_strip1=tempo_stp;
	    highest_t1=tempo_tpos;
	  }
	}
      }
      if(lowest_plane-10>=0) {
	lowest_plane-=10;
	lowest_z-=10.*0.06;
      }
      else {
	lowest_plane=0;
	lowest_z=0.;
      }
  
      if(lowest_strip0-5>=0) {
	lowest_strip0-=5;
	lowest_t0-=5.*0.041;
      }
      else {
	lowest_strip0=0;
	lowest_t0=-4.0;
      }
      
      if(lowest_strip1-5>=0) {
	lowest_strip1-=5;
	lowest_t1-=5.*0.041;
      }
      else {
	lowest_strip1=0;
	lowest_t1=-4.0;
      }
      
      if(highest_plane+10<=485) {
	highest_plane+=10;
	highest_z+=10.*0.06;
      }
      else {
	highest_plane=485;
	highest_z=30.;
      }
      
      if(highest_strip0+5<=191) {
	highest_strip0+=5;
	highest_t0+=5.*0.041;
      }
      else {
	highest_strip0=191;
	highest_t0=4.0;
      }
      
      if(highest_strip1+5<=191) {
	highest_strip1+=5;
	highest_t1+=5.*0.041;
      }
      else {
	highest_strip1=191;
	highest_t1=4.0;
      }      
    }
    // give some buffer at either end...
    tmin -= 50e-9;
    tmax += 20e-9;

    double eps = 1.0e-8;
    if (tmin == tmax) { tmin -= eps; tmax += eps; }
    
    // for far detector, suppress display of pre-trigger time interval
    if (fDetectorType == Detector::kFar){
      if ((tmax - tmin)*1e9>1500) tmin = tmax - 1500e-9;
    }

    TimeHstTrk->SetBins(100,tmin/1e-9,tmax/1e-9);
    TimeHstTrkU->SetBins(100,tmin/1e-9,tmax/1e-9);
    TimeHstTrkV->SetBins(100,tmin/1e-9,tmax/1e-9);
    TimeHstShw->SetBins(100,tmin/1e-9,tmax/1e-9);
    TimeHstShwU->SetBins(100,tmin/1e-9,tmax/1e-9);
    TimeHstShwV->SetBins(100,tmin/1e-9,tmax/1e-9);

    TimeHstTrk2->SetBins(50,0,100);
    TimeHstTrk2U->SetBins(50,0,100);
    TimeHstTrk2V->SetBins(50,0,100);
    TimeHstShw2->SetBins(50,0,100);
    TimeHstShw2U->SetBins(50,0,100);
    TimeHstShw2V->SetBins(50,0,100);

//    TH1F *h1 = (TH1F*)TimeHst->GetHists()->FindObject("TimeHstTrk");
//    TH1F *h2 = (TH1F*)TimeHst->GetHists()->FindObject("TimeHstShw");
//    h1->SetBins(100,tmin/1e-9,tmax/1e-9);
//    h2->SetBins(100,tmin/1e-9,tmax/1e-9);

    //record track hits
    //int ntrks = event->ntrack;
    if (ntrks){
      for (int i = 0; i<ntrks; i++){
	int index = SntpHelpers::GetTrackIndex(i,event);
	NtpSRTrack *track = SntpHelpers::GetTrack(index,sr);
	for (int istp = 0; istp<track->nstrip; istp++){
	  //ftrkshw->Fill(track->stpx[istp],track->stpy[istp],1);
	  NtpSRStrip *strip = SntpHelpers::GetStrip(track->stp[istp],sr);
	  TimeHstTrk->Fill(track->stpt0[istp]/1e-9,strip->ph0.pe);
	  TimeHstTrk->Fill(track->stpt1[istp]/1e-9,strip->ph1.pe);
	  TimeHstTrk2->Fill((track->stpt0[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph0.pe);
	  TimeHstTrk2->Fill((track->stpt1[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph1.pe);
	  //record track strip information
	  if(strip->planeview==PlaneView::kU){
	    fTrkUZ->Fill(strip->z,strip->tpos,100000);//paranoia

	    TimeHstTrkU->Fill(track->stpt0[istp]/1e-9,strip->ph0.pe);
	    TimeHstTrkU->Fill(track->stpt1[istp]/1e-9,strip->ph1.pe);
	    TimeHstTrk2U->Fill((track->stpt0[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph0.pe);
	    TimeHstTrk2U->Fill((track->stpt1[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph1.pe);
	  }
	  else if(strip->planeview==PlaneView::kV){
	    fTrkVZ->Fill(strip->z,strip->tpos,100000);

	    TimeHstTrkV->Fill(track->stpt0[istp]/1e-9,strip->ph0.pe);
	    TimeHstTrkV->Fill(track->stpt1[istp]/1e-9,strip->ph1.pe);
	    TimeHstTrk2V->Fill((track->stpt0[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph0.pe);
	    TimeHstTrk2V->Fill((track->stpt1[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph1.pe);	    
	  } 
	}
      }
    }

    //record shower hits
    int nshws = event->nshower;
    if (nshws){
      for (int i = 0; i<nshws; i++){
	int index = SntpHelpers::GetShowerIndex(i,event);
	NtpSRShower *shower = SntpHelpers::GetShower(index,sr);
	for (int istp = 0; istp<shower->nstrip; istp++){
	  NtpSRStrip *strip = SntpHelpers::GetStrip(shower->stp[istp],sr);
	  TimeHstShw->Fill(strip->time0/1e-9,strip->ph0.pe);
	  TimeHstShw->Fill(strip->time1/1e-9,strip->ph1.pe);
	  TimeHstShw2->Fill((strip->time0-tmin)/1e-9,strip->ph0.pe);
	  TimeHstShw2->Fill((strip->time1-tmin)/1e-9,strip->ph1.pe);
	  //record shower strip information
	  if(strip->planeview==PlaneView::kU){
	    fShwUZ->Fill(strip->z,strip->tpos,100000);//paranoia

	    TimeHstShwU->Fill(strip->time0/1e-9,strip->ph0.pe);
	    TimeHstShwU->Fill(strip->time1/1e-9,strip->ph1.pe);
	    TimeHstShw2U->Fill((strip->time0-tmin)/1e-9,strip->ph0.pe);
	    TimeHstShw2U->Fill((strip->time1-tmin)/1e-9,strip->ph1.pe);
	  }
	  else if(strip->planeview==PlaneView::kV){
	    fShwVZ->Fill(strip->z,strip->tpos,100000);

	    TimeHstShwV->Fill(strip->time0/1e-9,strip->ph0.pe);
	    TimeHstShwV->Fill(strip->time1/1e-9,strip->ph1.pe);
	    TimeHstShw2V->Fill((strip->time0-tmin)/1e-9,strip->ph0.pe);
	    TimeHstShw2V->Fill((strip->time1-tmin)/1e-9,strip->ph1.pe);
	  } 
	}
      }
    }
    TimeHst = new THStack("TimeHist","Digit Times (ns)");
    TimeHst->Add(TimeHstTrk);
    TimeHst->Add(TimeHstShw);

    TimeHstUV = new THStack("TimeHistUV","UV Digit Times (ns)");
    TimeHstUV->Add(TimeHstTrkU);
    TimeHstUV->Add(TimeHstTrkV);
    TimeHstUV->Add(TimeHstShwU);
    TimeHstUV->Add(TimeHstShwV);

    TimeHst2 = new THStack("TimeHist2","Digit Times (ns)");
    TimeHst2->Add(TimeHstTrk2);
    TimeHst2->Add(TimeHstShw2);

    TimeHst2UV = new THStack("TimeHist2UV","UV Digit Times (ns)");
    TimeHst2UV->Add(TimeHstTrk2U);
    TimeHst2UV->Add(TimeHstTrk2V);
    TimeHst2UV->Add(TimeHstShw2U);
    TimeHst2UV->Add(TimeHstShw2V);
  }  
  
  srvtx_u->SetX(event->vtx.z);
  srvtx_u->SetY(event->vtx.u);
  srvtx_v->SetX(event->vtx.z);
  srvtx_v->SetY(event->vtx.v);
  srvtx_xy->SetX(event->vtx.x);
  srvtx_xy->SetY(event->vtx.y);
  
  if(kDrawClu) FillClusterGraphs();  
  if(kIntReco) SetUpStripButtons();
  
  return 1;
  }
  
  void NueDisplayModule::FillClusterGraphs()
  {
    /*
      RecRecordImp<RecCandHeader> *rr = 
      dynamic_cast<RecRecordImp<RecCandHeader>*>
      ((gMint->GetJobC().Mom).GetFragment("RecRecordImp<RecCandHeader>"));
    */

  ssGraphU = new TMultiGraph();
  ssGraphU->SetName("usubshowers");
  ssGraphU->SetTitle("Transverse Position vs Z - U View");
  cluLegU = new TLegend(0.05,0.05,0.95,0.95,"  Key: ID (P_{EM})");
  cluLegU->SetBorderSize(0);
  cluLegU->SetTextSize(0.08);
  
  ssGraphV = new TMultiGraph();
  ssGraphV->SetName("vsubshowers");
  ssGraphV->SetTitle("Transverse Position vs Z - V View");
  cluLegV = new TLegend(0.05,0.05,0.95,0.95,"  Key: ID (P_{EM})");
  cluLegV->SetBorderSize(0);
  cluLegV->SetTextSize(0.08);

  int nUclus = 0;
  int col0 = 1;
  int nVclus = 0;
  int col1 = 1;
  
  for(int i=0;i<event->nshower;i++){  //loop over showers
    int shower_index = SntpHelpers::GetShowerIndex(i,event);
    NtpSRShower *shower = SntpHelpers::GetShower(shower_index,st);  
    if(shower==0) continue;

    int numclustp=0;
    int numclustp0=0;
    int numclustp1=0;
    
    for(int j=0;j<shower->ncluster;j++){ //loop over clusters
      int cluster_index = SntpHelpers::GetClusterIndex(j,shower);
      NtpSRCluster *cluster = SntpHelpers::GetCluster(cluster_index,st);
      numclustp+=cluster->nstrip;
      if(cluster->planeview==2) nUclus+=1;
      else if(cluster->planeview==3) nVclus+=1;
    }

    float *clu_tpos = new float[numclustp];
    float *clu_tpos0 = new float[numclustp];
    float *clu_z0 = new float[numclustp];
    float *clu_tpos1 = new float[numclustp];
    float *clu_z1 = new float[numclustp];
    double *clu_z = new double[numclustp];
    int count = 0;

    for(int j=0;j<shower->ncluster;j++){ //loop over clusters
      int cluster_index = SntpHelpers::GetClusterIndex(j,shower);
      NtpSRCluster *cluster = SntpHelpers::GetCluster(cluster_index,st);
      
      count=0;
      numclustp0=0;
      numclustp1=0;
      
      for(int k=0;k<cluster->nstrip;k++){
	int strip_index = SntpHelpers::GetStripIndex(k,cluster);
	NtpSRStrip *strip = SntpHelpers::GetStrip(strip_index,st);
	clu_tpos[count] = strip->tpos;
	clu_z[count]    = strip->z;
	if(strip->planeview==2){
	  clu_tpos0[numclustp0] = clu_tpos[count];
	  clu_z0[numclustp0] = clu_z[count];
	  numclustp0+=1;
	}
	else {
	  clu_tpos1[numclustp1] = clu_tpos[count];
	  clu_z1[numclustp1] = clu_z[count];
	  numclustp1+=1;
	}
	count++;
      }
    
      if(numclustp0>0){
	TGraph *temp = new TGraph(numclustp0,clu_z0,clu_tpos0);
	if(col0==10) col0+=1;
	temp->SetMarkerColor(col0);
	temp->SetMarkerSize(0.6);
	temp->SetMarkerStyle(21);
	if(cluster->id==2 ||
	   cluster->id==4) temp->SetMarkerStyle(25);
	ssGraphU->Add(temp);
	col0+=1;
	char ssnom[256];
	if(cluster->id==0) {
	  sprintf(ssnom,"EM (%.2f)",cluster->probem);
	}
	else if(cluster->id==1){
	  sprintf(ssnom,"HAD (%.2f)",cluster->probem);
	}
	else {
	  sprintf(ssnom,"%s",
		  ClusterType::AsString(ClusterType::
					EClusterType(cluster->id)));
	}
	cluLegU->AddEntry(temp,ssnom,"p");
      }
      if(numclustp1>0){
	if(col1==10) col1+=1;
	TGraph *temp = new TGraph(numclustp1,clu_z1,clu_tpos1);
	temp->SetMarkerColor(col1);
	temp->SetMarkerSize(0.6);
	temp->SetMarkerStyle(21);
	if(cluster->id==2 ||
	   cluster->id==4) temp->SetMarkerStyle(25);
	ssGraphV->Add(temp);
	col1+=1;
	char ssnom[256];
	if(cluster->id==0){
	  sprintf(ssnom,"EM (%.2f)",cluster->probem);
	}
	else if(cluster->id==1){
	  sprintf(ssnom,"HAD (%.2f)",cluster->probem);
	}
	else {
	  sprintf(ssnom,"%s",
		  ClusterType::AsString(ClusterType::
				      EClusterType(cluster->id)));
	}
	cluLegV->AddEntry(temp,ssnom,"p");
      }
    }
                                                                                
    delete [] clu_tpos;
    delete [] clu_z;
    delete [] clu_tpos0;
    delete [] clu_z0;
    delete [] clu_tpos1;
    delete [] clu_z1;
  }
  if(nUclus==0) {
    delete ssGraphU;
    ssGraphU=NULL;
    delete cluLegU;
    cluLegU = NULL;
  }
  if(nVclus==0) {
    delete ssGraphV;
    ssGraphV=NULL;
    delete cluLegV;
    cluLegV = NULL;
  }

}



void NueDisplayModule::DrawInteractionDiagram(Int_t index){
  //modified by steve cavanaugh aug 9, 2006 to take into account double particles in mc with IstHEP=2 
  //with modifications for readability
  
  fStdHepCan->cd();
  fStdHepCan->Range(0,0,1,1.1);
  
  vector<NtpMCStdHep *> hep;
  if (foundST) hep = SntpHelpers::GetStdHepArray(index, st);
  else hep = SntpHelpers::GetStdHepArray(index, mc);



  Int_t nStdHep = int(hep.size());
  Int_t *indicesToUse = new Int_t[nStdHep];
  Int_t *parent = new Int_t[nStdHep];
  Int_t *parent1 = new Int_t[nStdHep];
  Int_t incomingNeutrino = -1;
  Int_t cnt = 0;



  for(int i=0;i<nStdHep;i++){    
    if(hep[i]->mc==index) {
      indicesToUse[cnt] = i;



      //parent[i] = hep[i]->parent[0];

      //in the case where we have more than 1 event per snarl we might 
      //have indices which refer to the actual index, not the one for the particular event.... 
      //so make sure we match this up properly
      for(int j=0;j<nStdHep;j++){
        if ((Int_t)hep[j]->index==hep[i]->parent[0]){
          parent[i]=j;
          break;
        }
      }

      
      for(int j=0;j<nStdHep;j++){
        if ((Int_t)hep[j]->index==hep[i]->parent[1]){
          parent1[i]=j;
          break;
        }
      }
	



      if(hep[i]->IstHEP==0){
	if(abs(hep[i]->IdHEP)==12||abs(hep[i]->IdHEP)==14
	   ||abs(hep[i]->IdHEP)==16) {
	  incomingNeutrino=i;


	}
	parent[i] = -1; //don't draw arrows to initial state particles
	parent1[i] = -1;
      }
      cnt++;
    }
    else{
      parent[i] = -1;
      parent1[i] = -1;
    }
  }


  




  //make arrows and markers
  TArrow *arrow[1000];
  TMarker *marker[1000];
  for(int i=0;i<nStdHep;i++) {
    arrow[i] = new TArrow(0,0,0,0,.01,"|>"); 
    arrow[i]->SetLineWidth(1);
    arrow[i]->SetLineStyle(3);
    arrow[i]->SetFillColor(39);
    arrow[i]->SetLineColor(39);
    marker[i] = new TMarker(0,0,24);
    marker[i]->SetMarkerSize(.2);
    marker[i]->SetMarkerColor(38);
  }

  //now loop through valid stdhep entries and fill variables  
  Float_t Available[5] = {0.9,0.7,0.7,0.7};  
  
  //cout<<"nStdHep: "<<nStdHep<<"\n";

  for(int i=0;i<cnt;i++){
    
    int toUse = indicesToUse[i];
    if(hep[toUse]->IstHEP==999){
      parent[i]=-1;
      continue;
    }



    //    cout<<toUse<<"  "<<hep[toUse]->IstHEP<<"   "<<nStdHep<<endl;
    //cout<<parent[toUse]<<"   "<<parent1[toUse]<<endl;



    
    
    //sc - avoid any double taus' charms' etc, but set the proper parent/child 
    //(remove the particle between true parent and child)

    if((hep[toUse]->child[0]==hep[toUse]->child[1]) &&
       (hep[toUse]->IstHEP==2 || hep[toUse]->IstHEP==3 || hep[toUse]->IstHEP==14) && 
       hep[toUse]->child[0]!=-1){
      int j;
      int child = hep[toUse]->child[0];
      for(j=0;j<cnt;j++){
        if ((Int_t)hep[j]->index==child)break;
      }

      //we've lost arrow tail info, so fix it
      if (j<nStdHep){
	arrow[j]->SetX1(arrow[i]->GetX1());
	arrow[j]->SetY1(arrow[i]->GetY1());
      
	parent[i]=-1;
      
	continue;
      }else{
	cout<<"Error! - linking failed, index out of range! \n";
      }
    }

    //sc - avoid any "later decays" not comming from IstHEP=1 decays .. .these are children of double taus, etc

    if((hep[toUse]->IstHEP==205||hep[toUse]->IstHEP==1205)&&
       (hep[parent[toUse]]->IstHEP!=1 || hep[parent1[toUse]]->IstHEP!=1)){
      parent[i]=-1;
      //cout <<"breaking out of loop for particle "<<i<<"\n";
      continue;  
    }

    
    Float_t mom = sqrt(hep[toUse]->p4[0]*hep[toUse]->p4[0] + 
		       hep[toUse]->p4[1]*hep[toUse]->p4[1] + 
		       hep[toUse]->p4[2]*hep[toUse]->p4[2]);
    float x=0.,y=0.;
    int col=0;
    char text[256];

    //set x,y    
    if(hep[toUse]->IstHEP==0) {
      x = 0.05;
      y=Available[0]; Available[0] -= 0.1;
    }
    else if(hep[toUse]->IstHEP==2) {
      x = 0.15;
      y=Available[1]; Available[1] -= 0.1;   
    }
    else if(hep[toUse]->IstHEP==11) {
      x = 0.05;
      y=Available[0]; Available[0] -= 0.1; 
    }
    else if(hep[toUse]->IstHEP==3||hep[toUse]->IstHEP==14) {   //sc - allow for both IstHEP=3 and 14 intermediate decays
      x = 0.3;
      y=Available[1]; Available[1] -= 0.1;
    }
    else if(hep[toUse]->IstHEP==1){
      x = 0.55;
      y=Available[2]; Available[2] -= 0.1;
    }
    else if((hep[toUse]->IstHEP==205||hep[toUse]->IstHEP==1205)){
      x = 0.8;
      y=Available[3]; Available[3] -= 0.1;
    }
 
    //set colour and label (and override y in special cases)
    if(abs(hep[toUse]->IdHEP)==12) { //nue
      if(parent[toUse]==incomingNeutrino)  y = 0.9; 
      sprintf(text,"#nu_{e}"); col = 3;
      if(hep[toUse]->IdHEP<0) sprintf(text,"#bar{#nu}_{e}");
    }
    else if(abs(hep[toUse]->IdHEP)==14) { //numu
     if(parent[toUse]==incomingNeutrino)  y = 0.9;  
     sprintf(text,"#nu_{#mu}"); col = 4;
      if(hep[toUse]->IdHEP<0) sprintf(text,"#bar{#nu}_{#mu}");
    }
    else if(abs(hep[toUse]->IdHEP)==16) { //nutau
     if(parent[toUse]==incomingNeutrino)  y = 0.9;  
     sprintf(text,"#nu_{#tau}"); col = 105;
      if(hep[toUse]->IdHEP<0) sprintf(text,"#bar{#nu}_{#tau}"); 
    }    
    else if(abs(hep[toUse]->IdHEP)==11) { //e
      if(parent[toUse]==incomingNeutrino) y = 0.9;          
      sprintf(text,"e^{-}"); col = 3;
      if(hep[toUse]->IdHEP<0) sprintf(text,"e^{+}");
    }
    else if(abs(hep[toUse]->IdHEP)==13) { //mu
      if(parent[toUse]==incomingNeutrino) y = 0.9;               
      sprintf(text,"#mu^{-}"); col = 4;
      if(hep[toUse]->IdHEP<0) sprintf(text,"#mu^{+}");
    }
    else if(abs(hep[toUse]->IdHEP)==15) { //tau
      //tau will decay so it is easier to read if it is not on the top line...
      // if(parent[toUse]==incomingNeutrino) y = 0.9;                  
      sprintf(text,"#tau^{-}"); col = 105;
      if(hep[toUse]->IdHEP<0) sprintf(text,"#tau^{+}"); 
    }
    else if(hep[toUse]->IdHEP==22) { //photon
      sprintf(text,"#gamma"); col = 9;      
    }
    else if(hep[toUse]->IdHEP>1000000000) { //nucleus
      y = 0.8;
      sprintf(text,"nucleus(%i,%i)",int((hep[toUse]->IdHEP-1e9)/1e6),
	      int((hep[toUse]->IdHEP-1e9 - 1e6*int((hep[toUse]->IdHEP-1e9)
						  /1e6))/1e3)); 
      col = 15;
    }
    else if(hep[toUse]->IdHEP==2112){ 
      sprintf(text,"neutron"); col = 28;
    }
    else if(hep[toUse]->IdHEP==2212){
      sprintf(text,"proton"); col = 2;
    }
    else if(abs(hep[toUse]->IdHEP)==211) {
      sprintf(text,"#pi^{+}"); col = 6;
      if(hep[toUse]->IdHEP<0) sprintf(text,"#pi^{-}");
    }
    else if(hep[toUse]->IdHEP==111) {
      sprintf(text,"#pi^{0}"); col = 7;
    }
    else if(hep[toUse]->IdHEP==130) {
      sprintf(text,"K^{0}_{L}"); col = 31;
    }
    else if(hep[toUse]->IdHEP==310) {
      sprintf(text,"K^{0}_{S}"); col = 31;
    }
    else if(hep[toUse]->IdHEP==311) {
      sprintf(text,"K^{0}"); col = 31;
    }
    else if(abs(hep[toUse]->IdHEP)==321) {
      sprintf(text,"K^{+}"); col = 31;
      if(hep[toUse]->IdHEP<0) sprintf(text,"K^{-}"); col = 31;
    }
    else if(hep[toUse]->IdHEP==28) {
      sprintf(text,"Geantino"); col = 46;
      if(hep[toUse]->IdHEP<0) sprintf(text,"K^{-}"); col = 31;
    }
    else {
      sprintf(text,"ID: %i",hep[toUse]->IdHEP); col=43;
    }

    sprintf(text,"%s [%.1f GeV/c]",text,mom);
    
    arrow[toUse]->SetX2(x-0.02);   
    arrow[toUse]->SetY2(y-0.02);   
    marker[toUse]->SetX(x-0.02);
    marker[toUse]->SetY(y-0.02);

    for(int j=0;j<nStdHep;j++){
   
      if(parent[j]==toUse){
	arrow[j]->SetX1(x-0.02);
	arrow[j]->SetY1(y-0.02);
	//cout<<"writing arrows from "<<j << " to "<<parent[j]<<"\n";
      }
    }
 
    TLatex *tex = new TLatex(x,y,text);
    char texname[256];
    sprintf(texname,"tex%i",i);
    tex->SetName(texname);
    tex->SetTextSize(0.05);
    tex->SetTextColor(col);
    tex->Draw();
  }

  for(int i=0;i<nStdHep;i++){
    if(parent[i]==-1){
      delete arrow[i];
      delete marker[i];
    }
    else {
      arrow[i]->Draw();
      marker[i]->Draw();
      
    }
  }
  




  Float_t minAvail = 0;
  for(int i=0;i<4;i++){
    if(Available[i]<minAvail) minAvail = Available[i];
  }
  if(minAvail<0) fStdHepCan->Range(0,minAvail,1,1.1);

 

  delete [] indicesToUse;
  delete [] parent;

  fStdHepCan->Modified();
  fStdHepCan->Update();



}

void NueDisplayModule::Analyze(Int_t evt){
    DeleteOldCrap();

    shia.SetRelease(fRel);
    fva.SetRelease(fRel);
    sfa.SetRelease(fRel);
    hca.SetRelease(fRel);
    aca.SetRelease(fRel);
    acfa.SetRelease(fRel);

    if (foundST){
        fva.Analyze(evt,st);
    }
    else if (foundSR){
        fva.Analyze(evt,sr);
    }    

    fva.Draw(fFracVar_plots);
    fCanvas2->Update();
    
    MSG("NueDisplayModule",Msg::kDebug) <<"FracVar fract_road = " << fracvars.fract_road<<endl;

    
    if (foundST){
      sfa.Analyze(evt,st);
    }
    else if (foundSR){
      sfa.Analyze(evt,sr);
    }    

    shwfit.Draw(fShwfit_plots);
    fCanvas3->Update();


    DeqFloat_t x,y,z,e;
    Int_t primShow;
    DeqDeqInt_t clusterMap;
    TVector3 primDir;
    
    if (foundST){
        hca.Analyze(evt,st);
        hca.Get3DHit(x,y,z,e);    
        aca.Set3DHit(x,y,z,e);
        aca.Analyze(evt,st);    // Must be preceded by hca.Get3DHit
        //                     aca.Set3DHit
        aca.GetAngCluster(primShow,clusterMap,primDir);
        acfa.Set3DHit(x,y,z,e);
        acfa.SetAngCluster(primShow,clusterMap,primDir);
        acfa.Analyze(evt,st);    // Must be preceded by:
                                    // aca.GetAngCluster
                                    // acf.Set3DHit
                                    // acf.SetAngCluster
             
    }
    else if (foundSR){
        
        hca.Analyze(evt,sr);
        hca.Get3DHit(x,y,z,e);    
        aca.Set3DHit(x,y,z,e);
        aca.Analyze(evt,sr);    // Must be preceded by hca.Get3DHit
        //                     aca.Set3DHit
        aca.GetAngCluster(primShow,clusterMap,primDir);
        acfa.Set3DHit(x,y,z,e);
        acfa.SetAngCluster(primShow,clusterMap,primDir);
        acfa.Analyze(evt,sr);    // Must be preceded by:
                                 // aca.GetAngCluster
                                 // acf.Set3DHit
                                 // acf.SetAngCluster 

    }    

    acfa.Draw(fAngClusterFitAna_plots);
    fCanvas4->Update();

  if (fSimFlag == SimFlag::kMC){
    if (foundST){
      shia.Analyze(evt,st);
    }
    else if (foundMC&&foundTH){
      //not implemented yet
    }
  }


}

void NueDisplayModule::GetBasicInfo(){
  info1->Clear();
  info2->Clear();
  info3->Clear();
  info4->Clear();
  info41->Clear();
  info5->Clear();
  info6->Clear();
  info7->Clear();
  info8->Clear();
  info9->Clear();
  info10->Clear();
  info11->Clear();
  info12->Clear();
  info13->Clear();
  mrInfo1->Clear();
  char text1[100];
  char text2[100];
  char text3[100];
  char text31[100];
  char text4[100];
  char text5[100];
  char text6[100];
  char text7[100];
  char text8[100];
  char text9[100];
  char text10[100];
  char text11[100];
  char text12[100];
  int run=0, snarl=0, evt=0;
  DataUtil::GetRunSnarlEvent(&(gMint->GetJobC().Mom),run,snarl,evt);
  evt = fEvent;
  fSnarl = snarl;
 
  if(!kTestMode&&!kHideRunSnarl) sprintf(text1,"Run: %d/%d, Snarl: %d, Event: %d(%d), Slice: %d (%d)",
                 run,SubRunNo,snarl,evt,fNumEvents,event->slc,st->evthdr.nslice);
  else sprintf(text1,"Run: xxxx, Snarl: xxxx, Event: %d(%d)",evt,fNumEvents);
  
  int ntrks = event->ntrack;
  int nshws = event->nshower;

//  Float_t trk_mom = 0;
//  Float_t shw_eng = 0;
//  Int_t trk_leg = 0;
//  Int_t trk_like = 0;
//  Int_t shw_leg = 0;
  Float_t zenith = -1;
  Float_t azimuth = -1;

  Float_t trk_mom_range = 0;
  Float_t trk_mom_fit = 0;
  Float_t trk_length = 0;
  Float_t trk_qp = 0;
  Float_t trk_eqp = 0;
  Int_t trk_planes = 0;
  Int_t trk_like = 0;
  Float_t shw_ph = 0;
  Float_t shw_shwph = 0;
  Int_t shw_planes = 0;

  if (foundST){
    if (ntrks){
      for (int i = 0; i<ntrks; i++){
	int index = SntpHelpers::GetTrackIndex(i,event);
	NtpSRTrack *track = SntpHelpers::GetTrack(index,st);
	if (track->plane.n>trk_planes) {
	  zenith = track->cr.zenith;
	  azimuth = track->cr.azimuth;
	  trk_mom_range = track->momentum.range;
	  if (track->momentum.qp) trk_mom_fit = 1./track->momentum.qp;
	  trk_length = track->ds;
	  trk_qp = track->momentum.qp;
	  trk_eqp = track->momentum.eqp;
	  trk_planes = track->plane.n;
          trk_like= track->plane.ntrklike;
	}
      }
    }
  
    if (nshws){
      for (int i = 0; i<nshws; i++){
	int index = SntpHelpers::GetShowerIndex(i,event);
	NtpSRShower *shower = SntpHelpers::GetShower(index,st);
	if (shower->shwph.linCCgev>shw_shwph) {
	  shw_shwph = shower->shwph.linCCgev;
	  shw_ph = shower->ph.gev;
	  shw_planes = shower->plane.n;
	}
      }
    }
  }
  else if(foundSR){
    if (ntrks){
      for (int i = 0; i<ntrks; i++){
	int index = SntpHelpers::GetTrackIndex(i,event);
	NtpSRTrack *track = SntpHelpers::GetTrack(index,sr);
	if (track->plane.n>trk_planes) {
	  zenith = track->cr.zenith;
	  azimuth = track->cr.azimuth;
	  trk_mom_range = track->momentum.range;
	  if (track->momentum.qp) trk_mom_fit = 1./track->momentum.qp;
	  trk_length = track->ds;
	  trk_qp = track->momentum.qp;
	  trk_eqp = track->momentum.eqp;
	  trk_planes = track->plane.n;
          trk_like= track->plane.ntrklike;
	}
      }
    }
    
    if (nshws){
      for (int i = 0; i<nshws; i++){
	int index = SntpHelpers::GetShowerIndex(i,event);
	NtpSRShower *shower = SntpHelpers::GetShower(index,sr);
	if (shower->shwph.linCCgev>shw_shwph) {
	  shw_shwph = shower->shwph.linCCgev;
	  shw_ph = shower->ph.gev;
	  shw_planes = shower->plane.n;
	}
      }
    }
  }
    
  sprintf(text2,"ntrks: %d nshws: %d zenith: %.1f azimuth: %.1f",ntrks,nshws,zenith,azimuth);
  sprintf(text3,"trk.range:%.2fGeV, fit:%.2fGeV, qp:%.3f",trk_mom_range,trk_mom_fit,trk_qp);
  sprintf(text31,"trk.pls:%d, length:%.1fm, like:%d, eqp:%.2f",trk_planes,trk_length,trk_like,trk_eqp);
  sprintf(text4,"shw.gev:%.2fGeV, linCC:%.2fGeV, pls:%d",shw_ph,shw_shwph,shw_planes);
  
  sprintf(text5,"MC");
  
  mctruth = 0;
  if (fSimFlag == SimFlag::kMC){
    if (foundST){
      Int_t index = SntpHelpers::GetEvent2MCIndex(fEvent,st);
      mctruth = SntpHelpers::GetMCTruth(index,st);
      thevent = dynamic_cast<NtpTHEvent*>((*st->thevt)[fEvent]);
    }
    else if (foundMC&&foundTH){
      Int_t index = SntpHelpers::GetEvent2MCIndex(fEvent,th);
      mctruth = SntpHelpers::GetMCTruth(index,mc);
      thevent = dynamic_cast<NtpTHEvent*>((*th->thevt)[fEvent]);
    }
  }
    
  //MC
  if (mctruth){
    char tmptxt[100];
    sprintf(tmptxt,"(p%.2f/c%.2f/c%.2f)",thevent->purity,thevent->completeall,thevent->completeslc);
    if (mctruth->iaction==0){
      sprintf(text6,"NC event %s",tmptxt);
    }
    else{
      if (abs(mctruth->inunoosc)==12&&abs(mctruth->inu)==12){
	sprintf(text6,"beam #nu_{e} CC %s",tmptxt);
      }
      if (abs(mctruth->inunoosc)==14&&abs(mctruth->inu)==14){
	sprintf(text6,"#nu_{#mu}#rightarrow#nu_{#mu} CC %s",tmptxt);
      }
      if (abs(mctruth->inunoosc)==14&&abs(mctruth->inu)==12){
	sprintf(text6,"#nu_{#mu}#rightarrow#nu_{e} CC %s",tmptxt);
      }
      if (abs(mctruth->inunoosc)==14&&abs(mctruth->inu)==16){
	sprintf(text6,"#nu_{#mu}#rightarrow#nu_{#tau} CC %s",tmptxt);
      }
    }
    if (mctruth->iresonance==1001 && mctruth->iaction==1){
      sprintf(text7,"quasi-elastic");
    } 
    if (mctruth->iresonance==1001 && mctruth->iaction==0){
      sprintf(text7,"elastic");
    } 
    if (mctruth->iresonance==1002){
      sprintf(text7,"resonance production");
    }
    if (mctruth->iresonance==1003){
      sprintf(text7,"DIS");
    }
    if (mctruth->iresonance==1004){
      sprintf(text7,"coherent production");
    }
    sprintf(text8,"E_{#nu}: %.1fGeV, E_{shw}: %.1fGeV",fabs(mctruth->p4neu[3]),fabs(mctruth->p4shw[3]));
    sprintf(text9,"Y: %.2f, emfrac: %.2f",mctruth->y,mctruth->emfrac);
    sprintf(text10,"E_{#mu}: %.1fGeV,  E_{el}: %.1fGeV",fabs(mctruth->p4mu1[3]),fabs(mctruth->p4el1[3]));
    if (foundST) {
      sprintf(text11,"E_{#pi^{0}}_tot: %.2fGeV, E_{#pi^{0}}_neugen: %.2fGeV",stdhepinfo.epi0_total,stdhepinfo.epi0_neugen);
      sprintf(text12,"E_{#pi^{0}}_abs: %.2fGeV, E_{#pi^{0}}_intranuke: %.2fGeV",stdhepinfo.epi0_abs,stdhepinfo.epi0_intranuke);
    }
  }
   
  //summary
  info1->SetText(0.01,0.95,text1);
  info1->SetTextSize(0.06);
  info1->SetTextColor(2);
  
  //Reco
  info2->SetText(0.01,0.88,"Reco:");
  info2->SetTextSize(0.06);
  info2->SetTextColor(4);
  
  info3->SetText(0.01,0.82,text2);
  info3->SetTextSize(0.06);
  info3->SetTextColor(9);
  
  info4->SetText(0.01,0.76,text3);
  info4->SetTextSize(0.06);
  info4->SetTextColor(9);

  info41->SetText(0.01,0.7,text31);
  info41->SetTextSize(0.06);
  info41->SetTextColor(9);

  info5->SetText(0.01,0.64,text4);
  info5->SetTextSize(0.06);
  info5->SetTextColor(9);

  int mcinfoc = 6;
  info6->SetText(0.01,0.5,"MC:");
  info6->SetTextSize(0.06);
  info6->SetTextColor(mcinfoc);
  
  if (mctruth){
    info7->SetText(0.05,0.44,text6);
    info7->SetTextSize(0.06);
    info7->SetTextColor(mcinfoc);
    
    info8->SetText(0.05,0.38,text7);
    info8->SetTextSize(0.06);
    info8->SetTextColor(mcinfoc);

    info9->SetText(0.05,0.32,text8);
    info9->SetTextSize(0.06);
    info9->SetTextColor(mcinfoc);
    
    info10->SetText(0.05,0.26,text9);
    info10->SetTextSize(0.06);
    info10->SetTextColor(mcinfoc);

    info11->SetText(0.05,0.2,text10);
    info11->SetTextSize(0.06);
    info11->SetTextColor(mcinfoc);

    info12->SetText(0.05,0.14,text11);
    info12->SetTextSize(0.06);
    info12->SetTextColor(mcinfoc);

    info13->SetText(0.05,0.08,text12);
    info13->SetTextSize(0.06);
    info13->SetTextColor(mcinfoc);
    
  }
  else {
    info7->SetText(0.05,0.44,"N/A");
    info7->SetTextSize(0.06);
    info7->SetTextColor(mcinfoc);
  }
  
  fRecoInfo->Clear();
  fRecoInfo->AddLine(text1);
  fRecoInfo->AddLine("Reco Info");
  fRecoInfo->AddLine(text2);
  fRecoInfo->AddLine(text3);
  fRecoInfo->AddLine(text31);
  fRecoInfo->AddLine(text4);
  //  fRecoInfo->AddLine(" ");
  string passcuts[2] = {"Fail","Pass"};
//  string ipass = "Pre-selections: "+passcuts[passfid*passtrk*passshw]+"  Fid: "+passcuts[passfid]+"  Track: "+passcuts[passtrk]+"  Shower: "+passcuts[passshw];
  string ipass = "Pre-selection: "+passcuts[passfid*passtrk*passshw*passeng*passtrklike];
  string epass = "Fid: "+passcuts[passfid]+"  Trk: "+passcuts[passtrk]+"  Trklike: "+passcuts[passtrklike]+"  Shw: "+passcuts[passshw]+"  E: "+passcuts[passeng];
  fRecoInfo->AddLine(ipass.c_str());
  fRecoInfo->AddLine(epass.c_str());

}

bool NueDisplayModule::PassCuts(){

  //  int nshws = event->nshower;  //no. of showers

  passfid = 1;
  passtrk = 1;
  passshw = 1;
  passeng = 1;
  passtrklike = 1;
                                                                                
   Int_t test = 0;
   if(fDetectorType == Detector::kNear)
      test = NueConvention::IsInsideNearFiducial_Nue_Extended(
                  event->vtx.x, event->vtx.y, event->vtx.z);
   if(fDetectorType == Detector::kFar)
      test = NueConvention::IsInsideFarFiducial_Nue_Extended(
                  event->vtx.x, event->vtx.y, event->vtx.z);
   if(test <= 0) passfid = 0;
   MSG("NueDisplayModule",Msg::kDebug) << "IsFidVtxEvt? " << test
       << " Detector type: " << fDetectorType << endl;

  if(foundST)
  {
     fCut.SetInfoObject(fEvent,st);
     if(!fCut.PassesHighEnergyCut()) passeng=0;
     if(!fCut.PassesLowEnergyCut()) passeng=0;
     if(!fCut.PassesEventPlaneCut()) passeng=0;
     if(!fCut.PassesTrackPlaneCut()) passtrk = 0;
     if(!fCut.PassesTrackLikePlaneCut()) passtrklike = 0;
     if(!fCut.PassesHighShowerEnergyCut()) passshw = 0;
     if(!fCut.PassesLowShowerEnergyCut()) passshw = 0;
     
  }
                                                                                
  if(foundSR){
      MSG("NueDisplayModule",Msg::kWarning)
         << "Unable to Perform Pass Cuts for pre Birch files"<< endl;
  }

  if (passfid&&passtrk&&passshw&&passeng&&passtrklike) return true;
  else return false;

}
    
void NueDisplayModule::SetCuts(){ 
  if (preselec==1){ 
    fCuts->SetText("Cuts: OFF"); 
    preselec=0;
    fCuts->SetDown(false);
  }
  else if (preselec==0){ 
    fCuts->SetText("Cuts:  ON"); 
    preselec=1;
    fCuts->SetDown(true);
  }
}

void NueDisplayModule::NextEvent(){
  clickbutton = 1;
  if (iIO){
    ievtp = 0;
    itopo = 0;
  }
  if(fNumEvents>0){
    fEvent++;
    if(fEvent>=fNumEvents){
      fEvent = 0;
      gMint->Next();
      return;
    }
    this->UpdateDisplay(foundvrmatch, foundpidmatch);
  }
  else {
    gMint->Next();
    return;
  }
}


void NueDisplayModule::PrevEvent(){
  clickbutton = -1;
  if (iIO){
    ievtp = 0;
    itopo = 0;
  }
  if(fNumEvents>0){
    fEvent--;
    if(fEvent<0){
      clickbutton = -2;
      gMint->Prev();
      return;
    }
    this->UpdateDisplay(foundvrmatch, foundpidmatch);
  }
  else {
    gMint->Prev();
    return;
  }
}

void NueDisplayModule::GotoEvent()
{
  if (!gMint) return;

  string eventno = fEventEntry->GetText();
  if (atoi(eventno.c_str())<fNumEvents && atoi(eventno.c_str())>=0){
    fEvent = atoi(eventno.c_str());
    this->UpdateDisplay(foundvrmatch, foundpidmatch);
  }
  else {
    return;
  }
}

void NueDisplayModule::showmctruth()
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
    plotmc();
  }
  else {
    imctruth = 0;
    fMCTruth->SetDown(false);
    this->delmc();
  }
}

void NueDisplayModule::plotmc()
{
  vector<NtpMCStdHep *> hep;
  Int_t index = 0;
  if (foundST){
    index = SntpHelpers::GetEvent2MCIndex(fEvent,st);
    hep = SntpHelpers::GetStdHepArray(index, st);
  }
  else if (foundMC&&foundTH){
    index = SntpHelpers::GetEvent2MCIndex(fEvent,th);
    hep = SntpHelpers::GetStdHepArray(index, mc);
  }
  
  vector<double> vtx_u(hep.size());
  vector<double> vtx_v(hep.size());
  vector<double> vtx_z(hep.size());
  vector<double> p_u(hep.size());
  vector<double> p_v(hep.size());
  vector<double> p_z(hep.size());
  vector<double> p_tot(hep.size());
  vector<double> k_u(hep.size());
  vector<double> k_v(hep.size());
  vector<double> epar(hep.size());
  vector<int> idhep(hep.size());
  vector<int> drawline(hep.size());
  
  for (int istd = 0; istd < int(hep.size()); istd++){
    drawline[istd] = 0;
    vtx_u[istd] = (hep[istd]->vtx[0]+hep[istd]->vtx[1])*sqrt(2.)/2;
    vtx_v[istd] = (hep[istd]->vtx[1]-hep[istd]->vtx[0])*sqrt(2.)/2;
    vtx_z[istd] = hep[istd]->vtx[2];
    
    p_u[istd] = (hep[istd]->p4[0]+hep[istd]->p4[1])*sqrt(2.)/2;
    p_v[istd] = (hep[istd]->p4[1]-hep[istd]->p4[0])*sqrt(2.)/2;
    p_z[istd] = hep[istd]->p4[2];
    p_tot[istd] = sqrt(p_u[istd]*p_u[istd] + p_v[istd]*p_v[istd] + 
		       p_z[istd]*p_z[istd]);
    
    epar[istd] = hep[istd]->p4[3];
    idhep[istd] = abs(hep[istd]->IdHEP);
    if (fabs(p_z[istd])>0.) {
      k_u[istd] = p_u[istd]/p_z[istd];
      k_v[istd] = p_v[istd]/p_z[istd];
    }
    
    bool drawphoton = false;
    if (abs(hep[istd]->IdHEP)==22){//photon
      NtpMCStdHep* hep_parent = 0;
      if (foundST){
	hep_parent = 
	  dynamic_cast<NtpMCStdHep*>((*st->stdhep)[hep[istd]->parent[0]]);
      }
      else if (foundMC&&foundTH){
	hep_parent = 
	  dynamic_cast<NtpMCStdHep*>((*mc->stdhep)[hep[istd]->parent[0]]);
      }
      if (abs(hep_parent->IdHEP)!=111) drawphoton = true;
    }
    //cout<<istd<<" "<<hep[istd]->index<<" "<<hep[istd]->mc<<" "<<hep[istd]->IdHEP<<" "<<hep[istd]->parent[0]<<" "<<hep[istd]->child[0]<<endl;
    
    //decide what to draw
    if((hep[istd]->child[0]==-1 && hep[istd]->child[1]==-1 &&
	(hep[istd]->IdHEP) && //IdHEP == 0 means incoming particle???
	abs(hep[istd]->IdHEP)<10000 && //only draw particles
	(abs(hep[istd]->IdHEP)==22 && drawphoton || abs(hep[istd]->IdHEP)!=22)) || (abs(hep[istd]->IdHEP)==111&&hep[istd]->IstHEP!=14) //draw pi0s and photons that didn't orginate from pi0s
       ||abs(hep[istd]->IdHEP)==13)
      drawline[istd]=1;
  }
  
  int ipar = 0;
  
  for (int istd = 0; istd < int(hep.size()); istd++){
    if (drawline[istd] == 1){
      if (p_tot[istd]){
	paru.push_back(new TLine(vtx_z[istd],vtx_u[istd],
				 vtx_z[istd] + (p_z[istd]/p_tot[istd])*epar[istd]/3,
				 vtx_u[istd] + p_u[istd]/p_tot[istd]*epar[istd]/3));
	parv.push_back(new TLine(vtx_z[istd],vtx_v[istd],
				 vtx_z[istd] + (p_z[istd]/p_tot[istd])*epar[istd]/3,
				 vtx_v[istd] + p_v[istd]/p_tot[istd]*epar[istd]/3));
      }
      else {
	paru.push_back(new TLine(vtx_z[istd],vtx_u[istd],
				 vtx_z[istd],vtx_u[istd]));
	parv.push_back(new TLine(vtx_z[istd],vtx_v[istd],
				 vtx_z[istd],vtx_v[istd]));
      }
      //cout<<"ipar "<<ipar<<" "<<idhep[istd]<<endl;
      if(idhep[istd] == 11) {     //electron
	paru[ipar]->SetLineColor(3);
	parv[ipar]->SetLineColor(3);
      }
      else if(idhep[istd] == 13) {//muon
	paru[ipar]->SetLineColor(4);
	parv[ipar]->SetLineColor(4);
      }
      else if(idhep[istd] == 15) {//tau
	paru[ipar]->SetLineColor(5);
	parv[ipar]->SetLineColor(5);
      }
      else if(idhep[istd] == 211){//pion
	paru[ipar]->SetLineColor(6);
	parv[ipar]->SetLineColor(6);
      }
      else if(idhep[istd] == 2212){//proton
	paru[ipar]->SetLineColor(2);
	parv[ipar]->SetLineColor(2);
      }
      else if(idhep[istd] == 111) { //pi0
	paru[ipar]->SetLineColor(7);
	parv[ipar]->SetLineColor(7);
      }
      else if(idhep[istd] == 22){  //photon
	paru[ipar]->SetLineColor(9);
	parv[ipar]->SetLineColor(9);
      }
      else if(idhep[istd] == 2112){//neutron
	paru[ipar]->SetLineColor(28);
	parv[ipar]->SetLineColor(28);
      }
      else if(idhep[istd] == 321 || idhep[istd] == 311 || 
	      idhep[istd] == 310 || idhep[istd] == 130){//kaon
	paru[ipar]->SetLineColor(31);
	parv[ipar]->SetLineColor(31);
      }//anything else will be black
      else if(idhep[istd] == 12 || idhep[istd] == 14 ||
	      idhep[istd] == 16){  //outgoing neutrino
	paru[ipar]->SetLineStyle(2); //black, dashed line
	parv[ipar]->SetLineStyle(2);
      }
      ipar++;
    }
  }

  if(mctruth) {
    mcvtx_u->SetY((mctruth->vtxy+mctruth->vtxx)/sqrt(2.));
    mcvtx_u->SetX(mctruth->vtxz);
  
    mcvtx_v->SetY((mctruth->vtxy-mctruth->vtxx)/sqrt(2.));
    mcvtx_v->SetX(mctruth->vtxz);
  }

  fHistPad->cd(1);
  for (int ipar = 0; ipar<int(paru.size()); ipar++){
    paru[ipar]->Draw();
  }
  mcvtx_u->Draw();
  gPad->Modified();
  
  fHistPad->cd(2);
  for (int ipar = 0; ipar<int(parv.size()); ipar++){
    parv[ipar]->Draw();
  }
  mcvtx_v->Draw();
  gPad->Modified();
  
  fCanvas0->Update();
  
  fHistcolz->cd(1);
  for (int ipar = 0; ipar<int(paru.size()); ipar++){
      paru[ipar]->Draw();
  }
  mcvtx_u->Draw();
  gPad->Modified();
  
  fHistcolz->cd(2);
  for (int ipar = 0; ipar<int(parv.size()); ipar++){
    parv[ipar]->Draw();
  }
  mcvtx_v->Draw();
  gPad->Modified();
    
  fCanvas1->Update();

  if(kDrawClu){
    fCanvas5->cd(1);
    for (int ipar = 0; ipar<int(paru.size()); ipar++){
      paru[ipar]->Draw();
    }
    mcvtx_u->Draw();
    gPad->Modified();
    
    fCanvas5->cd(4);
    for (int ipar = 0; ipar<int(parv.size()); ipar++){
      parv[ipar]->Draw();
    }
    mcvtx_v->Draw();
    gPad->Modified();

    fCanvas5->cd(2);
    for (int ipar = 0; ipar<int(paru.size()); ipar++){
      paru[ipar]->Draw();
    }
    mcvtx_u->Draw();
    gPad->Modified();
    
    fCanvas5->cd(5);
    for (int ipar = 0; ipar<int(parv.size()); ipar++){
      parv[ipar]->Draw();
    }
    mcvtx_v->Draw();
    gPad->Modified();
    
    fCanvas5->Update();
  }

  fInfo0->cd();
  info6->Draw();
  info7->Draw();
  info8->Draw();
  info9->Draw();
  info10->Draw();
  info11->Draw();
  info12->Draw();
  info13->Draw();
  gPad->Modified();
  fCanvas0->Update();

  if(!kTestMode) DrawInteractionDiagram(index);
  
  return;
}

void NueDisplayModule::delmc()
{
  if (paru.size()) {
    for (int i = 0; i<int(paru.size()); i++){
      delete paru[i];
      delete parv[i];
    }
  }
  paru.clear();
  parv.clear();
  
  mcvtx_u->SetX(-100);
  mcvtx_v->SetY(-100);
  fInfo0->Clear();
  fInfo0->cd();
  info1->Draw();
  info2->Draw();
  info3->Draw();
  info4->Draw();
  info41->Draw();
  info5->Draw();
  //info6->Draw();
  gPad->Modified();
  fCanvas0->Update();
  
  fHistcolz->cd(1);
  gPad->Modified();
  fHistcolz->cd(2);
  gPad->Modified();
  fCanvas1->Update();

  if(kDrawClu){
    fCanvas5->cd(1);
    gPad->Modified();
    fCanvas5->cd(2);  
    gPad->Modified();
    fCanvas5->cd(4);  
    gPad->Modified();
    fCanvas5->cd(5);  
    gPad->Modified();
    fCanvas5->Update();
  }

  if(!kTestMode){
    fStdHepCan->cd();    
    TList *theList = fStdHepCan->GetListOfPrimitives();
    TIterator *iter = theList->MakeIterator();
    TObject *ob;
    while((ob = iter->Next())){
      if(ob->InheritsFrom("TLatex")) {
	TLatex *tex = (TLatex*) ob;
	delete tex;
      }
      else if(ob->InheritsFrom("TArrow")) {
	TArrow *ar = (TArrow*) ob;
	delete ar;
      }
      else if(ob->InheritsFrom("TMarker")) {
	TMarker *m = (TMarker*) ob;
	delete m;
      }
    }
    gPad->Modified();
    fStdHepCan->Update();
  }

}

void NueDisplayModule::fixmcinfo()
{
  if (fSimFlag!=SimFlag::kMC) return;
  if (!ifixmcinfo){
    fFixMCInfo->SetDown(true);
    ifixmcinfo = 1;
    if (!imctruth){
      plotmc();
    }
    else {
      fMCTruth->SetDown(false);
      imctruth = 0;
    }
  }
  else {
    fFixMCInfo->SetDown(false);
    ifixmcinfo = 0;
    this->delmc();
  }
}
  
void NueDisplayModule::PrintPads()
{
  fCanvas0->cd();
  //fHistPad->Print(Form("Evt_%d_%d_%d_%d.png",RunNo,SubRunNo,fSnarl,fEvent));
  fHistPad->Print(Form("Evt_%d_%d_%d_%d.eps",RunNo,SubRunNo,fSnarl,fEvent));
  fCanvas1->cd();
  fHistcolz->Print(Form("Evt_%d_%d_%d_%d_zoom.eps",RunNo,SubRunNo,fSnarl,fEvent));
  //fHistcolz->Print(Form("Evt_%d_%d_%d_%d_zoom.png",RunNo,SubRunNo,fSnarl,fEvent));
  fCanvas2->cd();
  fReco_plots->Print(Form("Evt_%d_%d_%d_%d_reco.eps",RunNo,SubRunNo,fSnarl,fEvent));
}

const Registry& NueDisplayModule::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
   MSG("NueDisplayModule",Msg::kDebug)<<"In NueDisplayModule::DefaultConfig"<<endl;

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

void NueDisplayModule::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
  MSG("NueDisplayModule",Msg::kDebug)<<"In NueDisplayModule::Config"<<endl;

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

void NueDisplayModule::DeleteOldCrap()
{
  //  cout<<"In delete old crap"<<endl;
  //  gDirectory->ls();
  TList *l = gDirectory->GetList();
  TIterator *it(l->MakeIterator());
  //  int NOBJ=l->GetSize();
  //  cout<<"Size of list "<<NOBJ<<endl;
  TObject *obj;
  while((obj=it->Next())){
    //    TObject *obj=l->At(i);
    const char* kn=obj->GetName();
    //    cout<<" Found a "<<kn<<endl;
    if(strstr(kn,"lenepl")!=NULL){
      obj->Delete();
      //      cout<<"deletint "<<kn<<endl;
    }
    if(strstr(kn,"tenestu")!=NULL){
      obj->Delete();
      //      cout<<"deletint "<<kn<<endl;
    }
    if(strstr(kn,"tenestv")!=NULL){
      obj->Delete();
      //      cout<<"deletint "<<kn<<endl;
    }
    if(strstr(kn,"LongHitHisto")!=NULL){
      obj->Delete();
      //      cout<<"deletint "<<kn<<endl;
    }
    if(strstr(kn,"TransHitHisto")!=NULL){
      obj->Delete();
      //      cout<<"deletint "<<kn<<endl;
    }
  }
}

void NueDisplayModule::updateEncoded(){
  string tmp = evtpcode[ievtp] + " | " + topocode[itopo];
  fEncoded->SetText(tmp.c_str());
}

void NueDisplayModule::OpenFileWrite(){

  if (!iFileW){
    string setfilename;
    string sRunNumber=Form("%d",RunNo);
    string sRunNumberSub=Form("%d",SubRunNo);
    setfilename = "PID-"+sRunNumber+"_"+sRunNumberSub;
    
    ifstream Test(setfilename.c_str());
    if(Test){
      //Need new filename
      Int_t fred=1;
      while(Test) {
	Test.close();
	string ifred = Form("%d",fred);
	setfilename = "PID-"+sRunNumber+"_"+sRunNumberSub+"_"+ifred;
	
	Test.open(setfilename.c_str());
	fred++;
      }
    }
    outfile.open(setfilename.c_str());
    string fn = "Write Records to File: "+setfilename;
    fFileW->SetText(fn.c_str());
    iFileW = 1; //file opened
  }
  else {
    //if (iIO) MSG("NueDisplayModule",Msg::kError)<<"In Read Mode"<<endl;
    MSG("NueDisplayModule",Msg::kInfo)<<"File already opened!"<<endl;
  }
}  

void NueDisplayModule::logvote(){
  if (!iFileW){
    this->OpenFileWrite();
  }
  if (!hitlog){
  
    //Copied from fcanvas0 to fill in energy values in the output file  
    int ntrks = event->ntrack;
    int nshws = event->nshower;
    Float_t trk_mom_range = 0;
    Float_t trk_planes = 0;
    Float_t shw_ph = 0;
    Float_t shw_shwph = 0;
    if (foundST){
    if (ntrks){
      for (int i = 0; i<ntrks; i++){
	int index = SntpHelpers::GetTrackIndex(i,event);
	NtpSRTrack *track = SntpHelpers::GetTrack(index,st);
	if (track->plane.n>trk_planes) {
	  trk_mom_range = track->momentum.range;
	  trk_planes = track->plane.n;
	}
	}
      }
    }
  
    if (nshws){
      for (int i = 0; i<nshws; i++){
	int index = SntpHelpers::GetShowerIndex(i,event);
	NtpSRShower *shower = SntpHelpers::GetShower(index,st);
	if (shower->shwph.linCCgev>shw_shwph) {
	  shw_ph = shower->ph.gev;
	  shw_shwph = shower->shwph.linCCgev;
	}
      }
    }
  
    
    outfile<<RunNo<<" "<<SubRunNo<<" "<<fSnarl<<" "<<fEvent<<" "<<ievtp<<" "<<itopo<<" "<<shw_ph+trk_mom_range<<" "<<shw_ph<<" "<<trk_mom_range;
    hitlog = 1;
    if (ievtp>=0&&ievtp<7&&itopo>=0&&itopo<5){
      //string tmp = evtpcode[ievtp] + " | " + topocode[itopo];
      string tmp = Form("%d | %d",ievtp, itopo);
      fEncoded->SetText(tmp.c_str());
    }
  }
}

void NueDisplayModule::makecomment(){
  if (hitlog ==1 && !icomm) {
    outfile<<"  "<<fComment->GetText()<<endl;
    fComment->Clear();
    icomm = 1;
  }
}

void NueDisplayModule::SetMode(){
  if (!iIO){ 
    fIO->SetText("Read"); 
    iIO = 1;
    fIO->SetDown(true);
    if (iFileW){
      if (!icomm && hitlog) outfile<<endl;
      outfile<<endl;
      outfile.close();
      hitlog = 0;
      icomm = 0;
      ievtp = 0;
      itopo = 0;
      this->updateEncoded();
      iFileW = 0;
      fFileW->SetText("");
    }
    string setfilename;
    string sRunNumber=Form("%d",RunNo);
    string sRunNumberSub=Form("%d",SubRunNo);
    setfilename = "PID-"+sRunNumber+"_"+sRunNumberSub;
    
    ifstream Test(setfilename.c_str());
    if(Test){
      fFileEnt->SetText(setfilename.c_str());
    }
  }
  else { 
    fIO->SetText("Write"); 
    iIO = 0;
    fIO->SetDown(false);
    fFileEnt->Clear();
    fFileR->SetText("");
    fFileW->SetText("");
    if (iFileW){
      if (!icomm && hitlog) outfile<<endl;
      outfile<<endl;
      outfile.close();
      hitlog = 0;
      iFileW = 0;
    }
    fTYPETOPO->SetText("");

    runno.clear();
    subrunno.clear();
    snarlno.clear();
    eventno.clear();
    eventtype.clear();
    eventtopo.clear();
    eventcomm.clear();
    
  }
}

void NueDisplayModule::updateselec(){
  if (!iEvtp[0]&&!iEvtp[1]&&!iEvtp[2]&&!iEvtp[3]&&!iEvtp[4]&&!iEvtp[5]&&!iEvtp[6]) selecevtp = 0;
  else selecevtp = 1;
  
  if (!iTopo[0]&&!iTopo[1]&&!iTopo[2]&&!iTopo[3]&&!iTopo[4]) selectopo = 0;
  else selectopo = 1;
}

void NueDisplayModule::OpenFileRead(){
  if (iIO == 1){
    ifstream ins(fFileEnt->GetText());
    if (ins){
      fFileR->SetText(Form("Read Records from File: %s",fFileEnt->GetText()));
      fFileEnt->Clear();
      clickbutton = 0;
      Int_t n1,n2,n3,n4,n5,n6;
      while(!ins.eof()){
	//ins.getline(buf,100);
	//sscanf(buf,"%d%d%d%d%d%d",&n1,&n2,&n3,&n4,&n5,&n6);
	ins>>n1>>n2>>n3>>n4>>n5>>n6;
	char buf[1000];
	ins.getline(buf,1000);
	if (!ins.eof()){
	  //cout<<n1<<" "<<n2<<" "<<n3<<" "<<n4<<" "<<n5<<" "<<n6<<endl;
	  //if (n1 == RunNo && n2 == SubRunNo){
	  runno.push_back(n1);
	  subrunno.push_back(n2);
	  snarlno.push_back(n3);
	  eventno.push_back(n4);
	  eventtype.push_back(n5);
	  eventtopo.push_back(n6);
	  if (strlen(buf)>2)
	    eventcomm.push_back(buf+2);
	  else eventcomm.push_back(buf);
	  //}
	}
      }//while
      runnoitr = runno.begin()-1;
      subrunnoitr = subrunno.begin()-1;
      snarlnoitr = snarlno.begin()-1;
      eventnoitr = eventno.begin()-1;
      typeitr = eventtype.begin()-1;
      topoitr = eventtopo.begin()-1;
      commitr = eventcomm.begin()-1;      
    }
  }
}

void NueDisplayModule::NextSelEvt(){
  if (iIO&&snarlno.size()){
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

void NueDisplayModule::PrevSelEvt(){
  if (iIO&&snarlno.size()){
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

void NueDisplayModule::IntRecoCalc(){

  fIntRecoHistZ->Reset();
  fIntRecoHistU->Reset();
  fIntRecoHistV->Reset();
  fPredictedHistZ->Reset();
  fPredictedHistU->Reset();
  fPredictedHistV->Reset();
  fUZPred->Reset();
  fVZPred->Reset();

  Bool_t doSim = true;
  if(fIntRecoDoSim->GetFillColor()==2) doSim = false;
  
  RecRecordImp<RecCandHeader> *rr = 
    dynamic_cast<RecRecordImp<RecCandHeader>*>
    ((gMint->GetJobC().Mom).GetFragment("RecRecordImp<RecCandHeader>"));

  vector<Int_t> striplist;
  for(UInt_t i=0;i<fStripButtonU.size();i++){
    if(fStripButtonU[i]->GetFillColor()==11){
      const char *name = fStripButtonU[i]->GetName();
      TString nom(name);
      nom.Chop();
      nom.Remove(0,nom.First("_")+1);
      striplist.push_back(atoi(nom.Data()));
    }
  }
  for(UInt_t i=0;i<fStripButtonV.size();i++){
    if(fStripButtonV[i]->GetFillColor()==11){
      const char *name = fStripButtonV[i]->GetName();
      TString nom(name);
      nom.Chop();
      nom.Remove(0,nom.First("_")+1);
      striplist.push_back(atoi(nom.Data()));
    }
  }
  if(striplist.size()==0) return;

  //translate arrow positions back to u,v,z in metres:
  Float_t evtZ1 = fUZcolz->GetXaxis()->GetBinCenter(1) + 
    (fUZcolz->GetXaxis()->GetBinLowEdge(fUZcolz->GetNbinsX()+1) - 
     fUZcolz->GetXaxis()->GetBinLowEdge(1)) * 
    (fIntRecoLineU->GetX1()-0.025)/0.95;
  Float_t evtU =  fUZcolz->GetYaxis()->GetBinCenter(1) + 
    (fUZcolz->GetYaxis()->GetBinLowEdge(fUZcolz->GetNbinsY()+1) - 
     fUZcolz->GetYaxis()->GetBinLowEdge(1)) * 
    (fIntRecoLineU->GetY1()-0.025)/0.95;
  Float_t evtZ = fVZcolz->GetXaxis()->GetBinCenter(1) + 
    (fVZcolz->GetXaxis()->GetBinLowEdge(fVZcolz->GetNbinsX()+1) - 
     fVZcolz->GetXaxis()->GetBinLowEdge(1)) * 
    (fIntRecoLineV->GetX1()-0.025)/0.95;
  Float_t evtV =  fVZcolz->GetYaxis()->GetBinCenter(0) + 
    (fVZcolz->GetYaxis()->GetBinLowEdge(fVZcolz->GetNbinsY()+1) - 
     fVZcolz->GetYaxis()->GetBinLowEdge(1)) * 
    (fIntRecoLineV->GetY1()-0.025)/0.95;

  Float_t endZ1 = fUZcolz->GetXaxis()->GetBinCenter(1) + 
    (fUZcolz->GetXaxis()->GetBinLowEdge(fUZcolz->GetNbinsX()+1) - 
     fUZcolz->GetXaxis()->GetBinLowEdge(1)) * 
    (fIntRecoLineU->GetX2()-0.025)/0.95;
  Float_t endU =  fUZcolz->GetYaxis()->GetBinCenter(1) + 
    (fUZcolz->GetYaxis()->GetBinLowEdge(fUZcolz->GetNbinsY()+1) - 
     fUZcolz->GetYaxis()->GetBinLowEdge(1)) * 
    (fIntRecoLineU->GetY2()-0.025)/0.95;
  Float_t endZ = fVZcolz->GetXaxis()->GetBinCenter(1) + 
    (fVZcolz->GetXaxis()->GetBinLowEdge(fVZcolz->GetNbinsX()+1) - 
     fVZcolz->GetXaxis()->GetBinLowEdge(1)) * 
    (fIntRecoLineV->GetX2()-0.025)/0.95;
  Float_t endV =  fVZcolz->GetYaxis()->GetBinCenter(0) + 
    (fVZcolz->GetYaxis()->GetBinLowEdge(fVZcolz->GetNbinsY()+1) - 
     fVZcolz->GetYaxis()->GetBinLowEdge(1)) * 
    (fIntRecoLineV->GetY2()-0.025)/0.95; 

  Float_t dudz = (endU-evtU)/(endZ1-evtZ1);
  Float_t dvdz = (endV-evtV)/(endZ-evtZ);

  //choose minimum Z value and adjust U/V vertex position
  if(evtZ1<evtZ) {
    evtV += dvdz*(evtZ1-evtZ);
    evtZ = evtZ1;
  }
  else evtU += dudz*(evtZ-evtZ1);

  Float_t minZ = 30;
  Float_t maxZ = -30;
  Float_t minTP = +8;
  Float_t maxTP = -8;

  Int_t minPl = 500;
  Int_t maxPl = 0;
  Int_t minStpU = 200;
  Int_t maxStpU = 0;
  Int_t minStpV = 200;
  Int_t maxStpV = 0;

  Int_t vtxView = 0;
  Int_t vtxPln = 0;
  Int_t vtxStpU = 0;
  Int_t vtxStpV = 0;
  Float_t vtxZDiff = 1000;
  Float_t vtxUDiff = 1000;
  Float_t vtxVDiff = 1000;

  vector<Int_t>::iterator stripitr = striplist.begin();
  vector<Int_t>::iterator stripend = striplist.end();
  while(stripitr!=stripend){
    NtpSRStrip *strip = SntpHelpers::GetStrip(*stripitr,rr);  
    if(strip->z-evtZ<minZ) {
      minZ = strip->z-evtZ;
      minPl = strip->plane;
    }
    if(strip->z-evtZ>maxZ) {
      maxZ = strip->z-evtZ;
      maxPl = strip->plane;
    }
    if(TMath::Abs(strip->z-evtZ)<TMath::Abs(vtxZDiff)) {
      vtxZDiff = strip->z-evtZ;
      vtxPln = strip->plane;
      vtxView = strip->planeview;
    }
    if(strip->planeview==2){
      if(strip->tpos-evtU<minTP) minTP=strip->tpos-evtU;
      if(strip->strip<minStpU) minStpU = strip->strip;
      if(strip->tpos-evtU>maxTP) maxTP=strip->tpos-evtU;
      if(strip->strip>maxStpU) maxStpU = strip->strip;
      if(TMath::Abs(strip->tpos-evtU)<TMath::Abs(vtxUDiff)) {
	vtxUDiff = strip->tpos-evtU;
	vtxStpU = strip->strip;
      }
    }
    if(strip->planeview==3){
      if(strip->tpos-evtV<minTP) minTP=strip->tpos-evtV;
      if(strip->strip<minStpV) minStpV = strip->strip;
      if(strip->tpos-evtV>maxTP) maxTP=strip->tpos-evtV;
      if(strip->strip>maxStpV) maxStpV = strip->strip;
      if(TMath::Abs(strip->tpos-evtV)<TMath::Abs(vtxVDiff)) {
	vtxVDiff = strip->tpos-evtV;
	vtxStpV = strip->strip;
      }
    }
    stripitr++;
  }
  stripitr = striplist.begin();
  
  // make sure vtx plane and strip are not more 
  // than 1plane or 1strip away from the true vertex
  while(TMath::Abs(vtxZDiff)>0.06) {
    if(vtxZDiff>0) {
      vtxZDiff -= 0.06;
      vtxPln -= 1;
      if(vtxView == 2) vtxView = 3;
      else vtxView = 2;
    }
    else if(vtxUDiff<0) {
      vtxZDiff += 0.06;
      vtxPln += 1;
      if(vtxView == 2) vtxView = 3;
      else vtxView = 2;
    }
  }
  while(TMath::Abs(vtxUDiff)>0.041) {
    if(vtxUDiff>0) {
      vtxUDiff -= 0.041;
      vtxStpU -= 1;
    }
    else if(vtxUDiff<0) {
      vtxUDiff += 0.041;
      vtxStpU += 1;
    }
  }
  while(TMath::Abs(vtxVDiff)>0.041) {
    if(vtxVDiff>0) {
      vtxVDiff -= 0.041;
      vtxStpV -= 1;
    }
    else if(vtxVDiff<0) {
      vtxVDiff += 0.041;
      vtxStpV += 1;
    }
  }

  fIntRecoHistZ->SetBins(5+int((maxZ-minZ)/0.0585),
			 minZ-0.1475,maxZ+0.1475);
  fIntRecoHistU->SetBins(5+int((maxTP-minTP)/0.0405),
			 minTP-0.1025,maxTP+0.1025);
  fIntRecoHistV->SetBins(5+int((maxTP-minTP)/0.0405),
			 minTP-0.1025,maxTP+0.1025);

  while(stripitr!=stripend){
    NtpSRStrip *strip = SntpHelpers::GetStrip(*stripitr,rr);    
    fIntRecoHistZ->Fill(strip->z-evtZ,(strip->ph0.sigcor + 
				       strip->ph1.sigcor)/SIGMAPMEU);
    if(strip->planeview==2) fIntRecoHistU->Fill(strip->tpos-evtU,
						(strip->ph0.sigcor + 
						 strip->ph1.sigcor) / 
						SIGMAPMEU);
    if(strip->planeview==3) fIntRecoHistV->Fill(strip->tpos-evtV,
						(strip->ph0.sigcor + 
						 strip->ph1.sigcor) / 
						SIGMAPMEU);
    stripitr++;
  }

  for(int i=1;i<=fIntRecoHistZ->GetNbinsX();i++){
    Float_t error = fIntRecoHistZ->GetBinContent(i);
    if(error>0) error = TMath::Sqrt(error);
    fIntRecoHistZ->SetBinError(i,error);
  }
  for(int i=1;i<=fIntRecoHistU->GetNbinsX();i++){
    Float_t error = fIntRecoHistU->GetBinContent(i);
    if(error>0) error = TMath::Sqrt(error);
    fIntRecoHistU->SetBinError(i,error);
    error = fIntRecoHistV->GetBinContent(i);
    if(error>0) error = TMath::Sqrt(error);
    fIntRecoHistV->SetBinError(i,error);
  }

  //Get prediction:
  fPredictedHistZ->SetBins(5+int((maxZ-minZ)/0.0585),
			   minZ-0.1475,maxZ+0.1475);
  fPredictedHistU->SetBins(5+int((maxTP-minTP)/0.0405),
			   minTP-0.1025,maxTP+0.1025);
  fPredictedHistV->SetBins(5+int((maxTP-minTP)/0.0405),
			   minTP-0.1025,maxTP+0.1025);
  
  fUZPred->SetBins(5+maxPl-minPl,minPl-2.5,maxPl+2.5,
		   5+maxStpU-minStpU,minStpU-2.5,maxStpU+2.5);
  fVZPred->SetBins(5+maxPl-minPl,minPl-2.5,maxPl+2.5,
		   5+maxStpV-minStpV,minStpV-2.5,maxStpV+2.5);

  Float_t totalPE = fIntRecoHistZ->Integral();
  Float_t energy = totalPE*0.0268;
  if(doSim && energy>0.1){
    FitterEM *fit = new FitterEM();
    fit->QuickInput(energy,vtxUDiff,vtxVDiff,vtxZDiff,
		    dudz,dvdz,vtxView);
    double binFluc = 1;
    for(int i=minPl-2;i<=maxPl+2;i++) {
      
      Int_t view = 0;
      if(i%2==vtxPln%2) view = vtxView;
      else if(vtxView==2) view=3;
      else if(vtxView==3) view=2;
      
      if(view==2){
	for(int j=minStpU-2;j<=maxStpU+2;j++) {	  
	  fUZPred->Fill(i,j,fit->PredictEMLoss(i-vtxPln,j-vtxStpU,binFluc));
	  fPredictedHistU->Fill(fPredictedHistU->GetBinCenter(fPredictedHistU->FindBin((j-vtxStpU)*0.041+vtxUDiff)),fit->PredictEMLoss(i-vtxPln,j-vtxStpU,binFluc));
	  fPredictedHistZ->Fill(fPredictedHistZ->GetBinCenter(fPredictedHistZ->FindBin((i-vtxPln)*0.059+vtxZDiff)),fit->PredictEMLoss(i-vtxPln,j-vtxStpU,binFluc));
	}
      }
      else if(view==3){
	for(int j=minStpV-2;j<=maxStpV+2;j++) {	  
	  fVZPred->Fill(i,j,fit->PredictEMLoss(i-vtxPln,j-vtxStpV,binFluc));
	  fPredictedHistV->Fill(fPredictedHistV->GetBinCenter(fPredictedHistV->FindBin((j-vtxStpV)*0.041+vtxVDiff)),fit->PredictEMLoss(i-vtxPln,j-vtxStpV,binFluc));
	  fPredictedHistZ->Fill(fPredictedHistZ->GetBinCenter(fPredictedHistZ->FindBin((i-vtxPln)*0.059+vtxZDiff)),fit->PredictEMLoss(i-vtxPln,j-vtxStpV,binFluc));
	}
      }
    }
    delete fit;
  }
  
  Float_t maxValue = 0;
  if(fUZPred->Integral()>0) {
    fUZPred->Scale(fUZcolz->Integral()/fUZPred->Integral());
    maxValue = fUZPred->GetBinContent(fUZPred->GetMaximumBin());
  }
  if(fVZPred->Integral()>0) {
    fVZPred->Scale(fVZcolz->Integral()/fVZPred->Integral());
    if(fVZPred->GetBinContent(fVZPred->GetMaximumBin())>maxValue) 
      maxValue = fVZPred->GetBinContent(fVZPred->GetMaximumBin());
  }
  if(maxValue>0) {
    fUZPred->SetMaximum(maxValue*1.1);
    fVZPred->SetMaximum(maxValue*1.1);
    fUZPred->SetMinimum(maxValue*0.01);
    fVZPred->SetMinimum(maxValue*0.01);
  }

  maxValue = 0;
  if(fPredictedHistU->Integral()>0) {
    fPredictedHistU->Scale(fIntRecoHistU->Integral() /
			   fPredictedHistU->Integral());
    maxValue = fPredictedHistU->GetBinContent(fPredictedHistU->GetMaximumBin());
  }
  if(fIntRecoHistU->Integral()>0 &&
     fIntRecoHistU->GetBinContent(fIntRecoHistU->GetMaximumBin())>maxValue) 
    maxValue = fIntRecoHistU->GetBinContent(fIntRecoHistU->GetMaximumBin());
  
  if(fPredictedHistV->Integral()>0) {
    fPredictedHistV->Scale(fIntRecoHistV->Integral() / 
			   fPredictedHistV->Integral());
    if(fPredictedHistV->GetBinContent(fPredictedHistV->GetMaximumBin())>maxValue)
      maxValue = fPredictedHistV->GetBinContent(fPredictedHistV->GetMaximumBin());
  }
  if(fIntRecoHistV->GetBinContent(fIntRecoHistV->GetMaximumBin())>maxValue) 
    maxValue = fIntRecoHistV->GetBinContent(fIntRecoHistV->GetMaximumBin());

  if(maxValue>0) {
    fPredictedHistU->SetMaximum(maxValue*1.1);
    fIntRecoHistU->SetMaximum(maxValue*1.1);
    fPredictedHistV->SetMaximum(maxValue*1.1);
    fIntRecoHistV->SetMaximum(maxValue*1.1);
  }

  maxValue = 0;
  if(fPredictedHistZ->Integral()>0) {
    fPredictedHistZ->Scale(fIntRecoHistZ->Integral() / 
			   fPredictedHistZ->Integral());
    maxValue = fPredictedHistZ->GetBinContent(fPredictedHistZ->GetMaximumBin());
  }
  if(fIntRecoHistZ->Integral()>0 && 
     fIntRecoHistZ->GetBinContent(fIntRecoHistZ->GetMaximumBin())>maxValue) 
    maxValue = fIntRecoHistZ->GetBinContent(fIntRecoHistZ->GetMaximumBin());
  if(maxValue>0) {
    fPredictedHistZ->SetMaximum(maxValue*1.1);
    fIntRecoHistZ->SetMaximum(maxValue*1.1);
  }
  
  fCanvas6->cd(1); gPad->Modified();
  fCanvas6->cd(2); gPad->Modified();
  fCanvas6->cd(3); gPad->Modified();
  fCanvas6->cd(4); gPad->Modified();
  fCanvas6->Update();
  
  cout << "=================================" << endl;
  cout << "Interative Display Calculations:" << endl;
  cout << "Energy of selected hits: " << energy << " GeV" << endl;
  cout << "Vertex from arrow positions (u,v,z) (metres): " 
       << evtU << " " << evtV << " " << evtZ << endl;
  cout << "Approximate Vertex in strip/plane (u,v,z): " 
       << vtxStpU << " " << vtxStpV << " " << vtxPln << endl;
  cout << "Offsets from approx vertex to arrow vertex: "
       << vtxUDiff << " " << vtxVDiff << " " << vtxZDiff << endl;
  cout << "Angles from arrow directions; dudz: " << dudz 
       << " dvdz: " << dvdz << endl;
  cout << "---------------------------------" << endl;
  cout << endl;
}
      
void NueDisplayModule::SetUpStripButtons(){

  for(UInt_t i=0;i<fStripButtonU.size();i++){
    delete fStripButtonU[i];
  }
  fStripButtonU.clear();
  for(UInt_t i=0;i<fStripButtonV.size();i++){
    delete fStripButtonV[i];
  }
  fStripButtonV.clear();

  Int_t nbinsx = fUZcolz->GetNbinsX();
  Float_t binwidthx = 1./Float_t(nbinsx);
  Int_t nbinsy = fUZcolz->GetNbinsY();
  Float_t binwidthy = 1./Float_t(nbinsy);
  Float_t max = fUZcolz->GetMaximum();
  Float_t minChargeUZ = 0;
  if(fUZcolz->GetMinimum()>minChargeUZ) minChargeUZ = fUZcolz->GetMinimum();
  Float_t minChargeVZ = 0;
  if(fVZcolz->GetMinimum()>minChargeVZ) minChargeVZ = fVZcolz->GetMinimum();

  Int_t ncol = gStyle->GetNumberOfColors();

  fSelectPad1->cd();
  for(int i=1;i<=nbinsx;i++){
    for(int j=1;j<=nbinsy;j++){
      if(fUZcolz->GetBinContent(i,j)>minChargeUZ){
	char name[256];
	sprintf(name,"stpbut_%i_",fStpIndexMapU[fUZcolz->GetBin(i,j)]);
	char command[256];
	sprintf(command,"fSelPad1->cd(); TButton *but = (TButton*) gPad->FindObject(\"%s\"); if(but->GetFillColor()==11) but->SetFillColor(but->GetLineColor()); else but->SetFillColor(11);",name);
	TButton *temp_but = new TButton("",command,
					0.025 + 0.95*float(i-1)*binwidthx,
					0.025 + 0.95*float(j-1)*binwidthy,
					0.025 + 0.95*float(i)*binwidthx,
					0.025 + 0.95*float(j)*binwidthy);
	temp_but->SetName(name);
	Float_t col_val = 0;
	if(max!=0) col_val = fUZcolz->GetBinContent(i,j)*float(ncol)/max;
	temp_but->SetFillColor(gStyle->GetColorPalette(Int_t(col_val)));
	temp_but->SetLineColor(gStyle->GetColorPalette(Int_t(col_val)));
	fStripButtonU.push_back(temp_but);
      }
    }
  }
  if(!fIntRecoLineU){
    fIntRecoLineU = new TArrow(0,0,1,1,0.015,"|>");
  }
  Int_t zVtxBin = fUZcolz->GetXaxis()->FindBin(event->vtx.z);
  Int_t tVtxBin = fUZcolz->GetYaxis()->FindBin(event->vtx.u);
  fIntRecoLineU->SetX1(0.025+0.95*float(zVtxBin-0.5)*binwidthx+0.95*binwidthx*
		       (event->vtx.z - 
			fUZcolz->GetXaxis()->GetBinCenter(zVtxBin)) / 
		       fUZcolz->GetXaxis()->GetBinWidth(zVtxBin) );  
  fIntRecoLineU->SetY1(0.025+0.95*float(tVtxBin-0.5)*binwidthy+0.95*binwidthy*
		       (event->vtx.u - 
			fUZcolz->GetYaxis()->GetBinCenter(tVtxBin)) / 
		       fUZcolz->GetYaxis()->GetBinWidth(tVtxBin) );
  fIntRecoLineU->SetX2(0.95);
  fIntRecoLineU->SetY2(fIntRecoLineU->GetY1());

  nbinsx = fVZcolz->GetNbinsX();
  binwidthx = 1./Float_t(nbinsx);
  nbinsy = fVZcolz->GetNbinsY();
  binwidthy = 1./Float_t(nbinsy);
  max = fVZcolz->GetMaximum();

  fSelectPad2->cd();
  for(int i=1;i<=nbinsx;i++){
    for(int j=1;j<=nbinsy;j++){
      if(fVZcolz->GetBinContent(i,j)>minChargeVZ) {
	char name[256];
	sprintf(name,"stpbut_%i_",fStpIndexMapV[fVZcolz->GetBin(i,j)]);
	char command[256];
	sprintf(command,"fSelPad2->cd(); TButton *but = (TButton*) gPad->FindObject(\"%s\"); if(but->GetFillColor()==11) but->SetFillColor(but->GetLineColor()); else but->SetFillColor(11);",name);
	TButton *temp_but = new TButton("",command,
					0.025 + 0.95*float(i-1)*binwidthx,
					0.025 + 0.95*float(j-1)*binwidthy,
					0.025 + 0.95*float(i)*binwidthx,
					0.025 + 0.95*float(j)*binwidthy);
	temp_but->SetName(name);
	Float_t col_val = 0;
	if(max!=0) col_val = fVZcolz->GetBinContent(i,j)*float(ncol)/max;
	temp_but->SetFillColor(gStyle->GetColorPalette(Int_t(col_val)));
	temp_but->SetLineColor(gStyle->GetColorPalette(Int_t(col_val)));
	fStripButtonV.push_back(temp_but);
      }
    }
  }

  if(!fIntRecoLineV){
    fIntRecoLineV = new TArrow(0,0,1,1,0.015,"|>");
  }
  zVtxBin = fVZcolz->GetXaxis()->FindBin(event->vtx.z);
  tVtxBin = fVZcolz->GetYaxis()->FindBin(event->vtx.v);
  fIntRecoLineV->SetX1(0.025+0.95*float(zVtxBin-0.5)*binwidthx+0.95*binwidthx*
		       (event->vtx.z - 
			fVZcolz->GetXaxis()->GetBinCenter(zVtxBin)) / 
		       fVZcolz->GetXaxis()->GetBinWidth(zVtxBin) );  
  fIntRecoLineV->SetY1(0.025+0.95*float(tVtxBin-0.5)*binwidthy+0.95*binwidthy*
		       (event->vtx.v - 
			fVZcolz->GetYaxis()->GetBinCenter(tVtxBin)) / 
		       fVZcolz->GetYaxis()->GetBinWidth(tVtxBin) );
  fIntRecoLineV->SetX2(0.95);
  fIntRecoLineV->SetY2(fIntRecoLineV->GetY1());

  fSelectPad1->cd();
  for(UInt_t i=0;i<fStripButtonU.size();i++){
    fStripButtonU[i]->Draw();
  }
  fIntRecoLineU->Draw();
  gPad->Modified();
  fSelectPad2->cd();
  for(UInt_t i=0;i<fStripButtonV.size();i++){
    fStripButtonV[i]->Draw();
  }
  fIntRecoLineV->Draw();
  gPad->Modified();
  fCanvas7->Update();
  
}

void NueDisplayModule::GenerateOverlayMultiGraphs()
{
      TGraph* temp;

      if(uzEventOverlay) delete uzEventOverlay;
      if(vzEventOverlay) delete vzEventOverlay;
      if(uzSliceOverlay) delete uzSliceOverlay;
      if(vzSliceOverlay) delete vzSliceOverlay;

      uzEventOverlay = new TMultiGraph();
      uzEventOverlay->SetName("uzEventOverlay");
      uzEventOverlay->SetTitle("U vs Z Event View");

      vzEventOverlay = new TMultiGraph();
      vzEventOverlay->SetName("vzEventOverlay");
      vzEventOverlay->SetTitle("V vs Z Event view");

      uzSliceOverlay = new TMultiGraph();
      uzSliceOverlay->SetName("uzSliceOverlay");
      uzSliceOverlay->SetTitle("U vs Z Slice view");

      vzSliceOverlay = new TMultiGraph();
      vzSliceOverlay->SetName("vzSliceOverlay");
      vzSliceOverlay->SetTitle("V vs Z Slice view");

      //Drawing of Events
      for(int evtNum = 0; evtNum < st->evthdr.nevent; evtNum++)
      {
        NtpSREvent* Event = SntpHelpers::GetEvent(evtNum,st);

        Int_t TotalStrip = Event->nstrip;
        Int_t NumUStrip = 0;
        Int_t NumVStrip = 0;

        Float_t* stripUPos = new Float_t[TotalStrip];
        Float_t* stripVPos = new Float_t[TotalStrip];
        Float_t* stripUPosZ = new Float_t[TotalStrip];
        Float_t* stripVPosZ = new Float_t[TotalStrip];

        for(int stpn=0;stpn < Event->nstrip;stpn++){
          Int_t index = SntpHelpers::GetStripIndex(stpn,Event);
          NtpSRStrip *strip = SntpHelpers::GetStrip(index,st);
          if(strip->planeview==PlaneView::kU){
            MSG("NueDisplayModule",Msg::kDebug)<<"filling u"<<endl;
            stripUPos[NumUStrip] = strip->tpos;
            stripUPosZ[NumUStrip] = strip->z;
            NumUStrip++;
          }
          if(strip->planeview==PlaneView::kV){
            MSG("NueDisplayModule",Msg::kDebug)<<"filling v"<<endl;
            stripVPos[NumVStrip] = strip->tpos;
            stripVPosZ[NumVStrip] = strip->z;
            NumVStrip++;
          }
        }

        if(NumUStrip > 0){
          temp = new TGraph(NumUStrip,stripUPosZ,stripUPos);
          uzEventOverlay->Add(temp);
        }
        if(NumVStrip > 0){
          temp = new TGraph(NumVStrip,stripVPosZ,stripVPos);
          vzEventOverlay->Add(temp);
        }

        delete [] stripUPos;
        delete [] stripVPos;
        delete [] stripUPosZ;
        delete [] stripVPosZ;
     }
 
     for(int slcNum = 0; slcNum < st->evthdr.nslice; slcNum++)
     {
        NtpSRSlice* slice = SntpHelpers::GetSlice(slcNum,st);

        Int_t TotalStrip = slice->nstrip;
        Int_t NumUStrip = 0;
        Int_t NumVStrip = 0;

        Float_t* stripUPos = new Float_t[TotalStrip];
        Float_t* stripVPos = new Float_t[TotalStrip];
        Float_t* stripUPosZ = new Float_t[TotalStrip];
        Float_t* stripVPosZ = new Float_t[TotalStrip];

        for(int stpn=0;stpn < slice->nstrip;stpn++){
          Int_t index = SntpHelpers::GetStripIndex(stpn,slice);
          NtpSRStrip *strip = SntpHelpers::GetStrip(index,st);
          if(strip->planeview==PlaneView::kU){
            MSG("NueDisplayModule",Msg::kDebug)<<"filling u"<<endl;
            stripUPos[NumUStrip] = strip->tpos;
            stripUPosZ[NumUStrip] = strip->z;
            NumUStrip++;
          }
          if(strip->planeview==PlaneView::kV){
            MSG("NueDisplayModule",Msg::kDebug)<<"filling v"<<endl;
            stripVPos[NumVStrip] = strip->tpos;
            stripVPosZ[NumVStrip] = strip->z;
            NumVStrip++;
          }
        }

        if(NumUStrip > 0){
          temp = new TGraph(NumUStrip,stripUPosZ,stripUPos);
          uzSliceOverlay->Add(temp);
        }
        if(NumVStrip > 0){
          temp = new TGraph(NumVStrip,stripVPosZ,stripVPos);
          vzSliceOverlay->Add(temp);
        }

        delete [] stripUPos;
        delete [] stripVPos;
        delete [] stripUPosZ;
        delete [] stripVPosZ;
     }

  return ;
}

void NueDisplayModule::UpdateOverlayGraphs(int evt)
{
   static int currentSnarl = -10;
   if(fSnarl != currentSnarl)
   {
      GenerateOverlayMultiGraphs();
      currentSnarl = fSnarl;
//      cout<<"current"<<currentSnarl<<"  "<<fSnarl<<endl;
   }

//   if(uzEventOverlay == 0) cout<<"this is still NULL"
   UpdateEventOverlayColors(evt, uzEventOverlay);
   UpdateEventOverlayColors(evt, vzEventOverlay);

   UpdateSliceOverlayColors(evt, uzSliceOverlay);
   UpdateSliceOverlayColors(evt, vzSliceOverlay);

   return;
}

void NueDisplayModule::UpdateEventOverlayColors(int evt, TMultiGraph *mgr)
{
   TList* graphs = mgr->GetListOfGraphs();
   if(graphs==0) return;

   TGraph *gr;
   for(int evtNum = 0; evtNum < st->evthdr.nevent; evtNum++)
   {
      gr = dynamic_cast<TGraph*>(graphs->At(evtNum));

      if(gr == 0) continue;
      gr->SetMarkerSize(0.6);
      gr->SetMarkerStyle(8);

      if(evtNum < evt - 1) gr->SetMarkerColor(kMagenta);
      if(evtNum == evt - 1) gr->SetMarkerColor(kRed);
      if(evtNum == evt){
           gr->SetMarkerColor(kBlack);
           gr->SetMarkerStyle(3);
      }
      if(evtNum == evt + 1) gr->SetMarkerColor(kBlue);
      if(evtNum > evt + 1) gr->SetMarkerColor(kCyan);
      if(kShowRange > 0){
        if(TMath::Abs(evtNum - evt) > kShowRange)
        {
            gr->SetMarkerSize(0);
            gr->SetMarkerStyle(1);
            gr->SetMarkerColor(kWhite);
        }
      }
  }

  return;
}

void NueDisplayModule::UpdateSliceOverlayColors(int evt, TMultiGraph *mgr)
{
   TList* graphs = mgr->GetListOfGraphs();
   if(graphs==0) return;

   TGraph *gr;
   NtpSREvent* event = SntpHelpers::GetEvent(evt,st);
   int sliceIndex = event->slc;
   for(int evtNum = 0; evtNum < st->evthdr.nslice; evtNum++)
     {
      gr = dynamic_cast<TGraph*>(graphs->At(evtNum));

      if(gr == 0) continue;
      gr->SetMarkerSize(0.6);
      gr->SetMarkerStyle(8);

      if(evtNum < sliceIndex - 1) gr->SetMarkerColor(kMagenta);
      if(evtNum == sliceIndex - 1) gr->SetMarkerColor(kRed);
      if(evtNum == sliceIndex){
           gr->SetMarkerColor(kBlack);
           gr->SetMarkerStyle(3);
      }
      if(evtNum == sliceIndex + 1) gr->SetMarkerColor(kBlue);
      if(evtNum > sliceIndex + 1) gr->SetMarkerColor(kCyan);

      if(kShowRange > 0){
        if(TMath::Abs(evtNum - sliceIndex) > kShowRange)
        {
            gr->SetMarkerSize(0);
            gr->SetMarkerStyle(1);
            gr->SetMarkerColor(kWhite);
        }
      }
  }


  return;
}

void NueDisplayModule::AdjustOverlayRange(int range)
{

  kShowRange = range;

  if(uzEventOverlay == 0) return;

  UpdateEventOverlayColors(fEvent, uzEventOverlay);
  UpdateEventOverlayColors(fEvent, vzEventOverlay);

  UpdateSliceOverlayColors(fEvent, uzSliceOverlay);
  UpdateSliceOverlayColors(fEvent, vzSliceOverlay);

  uzSliceOverlay->Draw("AP");
  vzSliceOverlay->Draw("AP");
  uzEventOverlay->Draw("AP");
  vzEventOverlay->Draw("AP");

  fCanvas9->cd(1); gPad->Modified();
  fCanvas9->cd(2); gPad->Modified();
  fCanvas9->cd(3); gPad->Modified();
  fCanvas9->cd(4); gPad->Modified();

  fCanvas9->Update();
  fCanvas9->Modified();

  return;
}


int NueDisplayModule::UpdateMRGraphs(int evt)
{
   static int currentSnarlMR = -10;

   static vector<int> MRtoNewEventMap;
   static vector<int> MRtoOldEventMap;

   if(fSnarl != currentSnarlMR)
   {
      CreateMRMap(MRtoNewEventMap, MRtoOldEventMap);
      currentSnarlMR = fSnarl;
   }

   vector<int> MREvent;
   vector<int> OriginalEvent;

   for(unsigned int mrind = 0; mrind < MRtoNewEventMap.size(); mrind++)
   {
      if(evt == MRtoNewEventMap[mrind]){
         MREvent.push_back(mrind);
         OriginalEvent.push_back(MRtoOldEventMap[mrind]);
//         cout<<"PushingBack "<<mrind<<"  "<<MRtoOldEventMap[mrind]<<"  "<<evt<<endl;
      }
   }
   
   //both lists could be empty if this was really a junk event

   GenerateMRMultiGraphs(evt, MREvent, OriginalEvent);      

   UpdateMRInformation(evt, MREvent, OriginalEvent);
                                     
   return MREvent.size();
}

void NueDisplayModule::CreateMRMap(vector<int> &MRnewEvent, vector<int> &MRoldEvent)
{

   MRnewEvent.clear();
   MRoldEvent.clear();

   for(int mrevt = 0; mrevt < mr->mrhdr.nmrevt; mrevt++)
   {
      NtpMREvent* mrevent = SntpHelpers::GetMREvent(mrevt, mr);    
      MRoldEvent.push_back(mrevent->orig_event);
      MRnewEvent.push_back(mrevent->best_event);
//      cout<<"PushingBack "<<mrevent->orig_event<<"  "
//          <<(mrevent->best_event)<<endl;

   }
}


void NueDisplayModule::GenerateMRMultiGraphs(int evtNum, 
                        vector<int> &MREvent, vector<int> &OldEvent)
{
   //If there is no MRevent then there should be an info which pops
   // up to say so.

   ResetMRGraphs();

//   SetFrame(vzMREventOld);
//   SetFrame(uzMREventOld);
//   SetFrame(vzMREventOverlay);
//   SetFrame(uzMREventOverlay);
 
   //Draw the New Event

   NtpSREvent* Event = 0;
    
//   FillEventGraph(Event, st, uzMREventOverlay, vzMREventOverlay, "New");

   for(unsigned int oEvent = 0; oEvent < OldEvent.size(); oEvent++)
   {
     Event = SntpHelpers::GetEvent(OldEvent[oEvent],stOld);      
     FillEventGraph(Event, stOld, uzMREventOverlay, vzMREventOverlay, "Old");
   }

   for(unsigned int mrEvt = 0; mrEvt < MREvent.size(); mrEvt++)
   {
       NtpMREvent* mrevt = SntpHelpers::GetMREvent(MREvent[mrEvt], mr);
//       cout<<"Contributing from Event: "<<MREvent[mrEvt]<<endl;

       FillMREventGraph(mrevt, st, uzMREventOverlay, vzMREventOverlay, 0);
       for(int striptype = 0; striptype < 7; striptype++)
        FillMREventGraph(mrevt, st, uzMREventOld, vzMREventOld, striptype);
   }

   Event = SntpHelpers::GetEvent(evtNum,st);
                                                                                
   FillEventGraph(Event, st, uzMREventOverlay, vzMREventOverlay, "New");


   UpdateMREventOverlayColors(uzMREventOverlay);
   UpdateMREventOverlayColors(vzMREventOverlay);

   UpdateMREventStripColors(uzMREventOld);
   UpdateMREventStripColors(vzMREventOld);
                                                                                                                   

/*
       //Now overlay it with the original event strips
       NtpMREvent* mrevt = SntpHelpers::GetMREvent(evtNum, mr);
       cout<<"This was original event "<<mrevt->orig_event<<endl;
       Int_t oEvent = mrevt->best_event;
  
     //          temp->SetMarkerSize(0.6);
     //          temp->SetMarkerStyle(8);
*/

}


void NueDisplayModule::ResetMRGraphs()
{
      if(uzMREventOverlay) delete uzMREventOverlay;
      if(vzMREventOverlay) delete vzMREventOverlay;
                                                                                
      if(uzMREventOld) delete uzMREventOld;
      if(vzMREventOld) delete vzMREventOld;
                                                                                
      uzMREventOverlay = new TMultiGraph();
      uzMREventOverlay->SetName("uzMREventOverlay");
      uzMREventOverlay->SetTitle("U vs Z Event View");
                                                                                
      vzMREventOverlay = new TMultiGraph();
      vzMREventOverlay->SetName("vzMREventOverlay");
      vzMREventOverlay->SetTitle("V vs Z Event view");
                                                                                
      uzMREventOld = new TMultiGraph();
      uzMREventOld->SetName("uzMREventOld");
      uzMREventOld->SetTitle("U vs Z MR Strips View");
                                                                                
      vzMREventOld = new TMultiGraph();
      vzMREventOld->SetName("vzMREventOld");
      vzMREventOld->SetTitle("V vs Z MR Strips view");
}



void  NueDisplayModule::FillEventGraph(NtpSREvent* Event, NtpStRecord *str,
          TMultiGraph* uzGraph, TMultiGraph *vzGraph, TString name)
{
  Int_t TotalStrip = Event->nstrip;
  Int_t NumUStrip = 0;
  Int_t NumVStrip = 0;
                                                                               
  Float_t* stripUPos = new Float_t[TotalStrip];
  Float_t* stripVPos = new Float_t[TotalStrip];
  Float_t* stripUPosZ = new Float_t[TotalStrip];
  Float_t* stripVPosZ = new Float_t[TotalStrip];

//  cout<<"nstrip - "<<Event->nstrip<<endl;                                                                                
  for(int stpn=0;stpn < Event->nstrip;stpn++){
      Int_t index = SntpHelpers::GetStripIndex(stpn,Event);
      NtpSRStrip *strip = SntpHelpers::GetStrip(index,str);
//      cout<<"For "<<name<<" "<<index<<"  "<<strip->tpos<<"   "<<strip->z<<endl;
     
      if(strip->planeview==PlaneView::kU){
        MSG("NueDisplayModule",Msg::kDebug)<<"filling u"<<endl;
        stripUPos[NumUStrip] = strip->tpos;
        stripUPosZ[NumUStrip] = strip->z;
        NumUStrip++;
      }
      if(strip->planeview==PlaneView::kV){
        MSG("NueDisplayModule",Msg::kDebug)<<"filling v"<<endl;
        stripVPos[NumVStrip] = strip->tpos;
        stripVPosZ[NumVStrip] = strip->z;
        NumVStrip++;
      }
  }

  TGraph *temp;

  if(NumUStrip > 0){
     temp = new TGraph(NumUStrip,stripUPosZ,stripUPos);
     if(name != "NULL")
        temp->SetName("u"+name);
     uzGraph->Add(temp);
  }
  if(NumVStrip > 0){
     temp = new TGraph(NumVStrip,stripVPosZ,stripVPos);
     if(name != "NULL")
        temp->SetName("v"+name);
     vzGraph->Add(temp);
  }
                                                                                
  delete [] stripUPos;
  delete [] stripVPos;
  delete [] stripUPosZ;
  delete [] stripVPosZ;
}


void NueDisplayModule::SetFrame(TMultiGraph* mg)
{
      if(mg == 0) return;

      TGraph* temp = 0;
                                                                                
      Float_t FarShell[2] = {-4.0, 4.0};
      Float_t FarShellZ[2] = {0.0, 30.0};
                                                                                
      Float_t NearShell[2] = {-3.0, 3.0};
      Float_t NearShellZ[2] = {0.0, 16.0};

      if(fDetectorType == Detector::kNear)
          temp = new TGraph(2,NearShellZ,NearShell);

      if(fDetectorType == Detector::kFar)
          temp = new TGraph(2,FarShellZ,FarShell);

      temp->SetName("frame");
      temp->SetMarkerSize(0);
      temp->SetMarkerStyle(1);
      mg->Add(temp);
}

void  NueDisplayModule::FillMREventGraph(NtpMREvent* mrevt, NtpStRecord *str,
          TMultiGraph* uzGraph, TMultiGraph *vzGraph, int type)
{
  if(!foundMR)
  {
      MSG("NueDisplayModule",Msg::kError) << "Call to FillMREvent with no MR records"<<endl;
      return;
  }

  Int_t TotalStrip = mrevt->nstrip;
  Int_t NumUStrip = 0;
  Int_t NumVStrip = 0;
                                                                                
  Float_t* stripUPos = new Float_t[TotalStrip];
  Float_t* stripVPos = new Float_t[TotalStrip];
  Float_t* stripUPosZ = new Float_t[TotalStrip];
  Float_t* stripVPosZ = new Float_t[TotalStrip];

  Bool_t StripType[7];
  StripType[0] = true;  // Just that its a strip

  //cout<<mrevt->nstrip<<endl;                                                                                
  for(int stpn=0;stpn < mrevt->nstrip;stpn++){
      Int_t index = mrevt->stp[stpn];
      StripType[1] = mrevt->StpIsTrueMu(stpn);
      StripType[2] = mrevt->StpIsTrueShw(stpn);
      StripType[3] = mrevt->StpIsScaled(stpn);
      StripType[4] = mrevt->StpIsInRecoTrk(stpn);
      StripType[5] = mrevt->StpIsRetained(stpn);
      StripType[6] = mrevt->StpIsMCElec(stpn);
      
//      cout<<type<<"  "<<StripType[type]<<endl;
      if(!StripType[type]) continue;

      NtpSRStrip *strip = SntpHelpers::GetStrip(index,str);
//      cout<<"strip->tpos"<<strip->tpos<<"   "<<strip->z<<endl;
      if(strip->planeview==PlaneView::kU){
        MSG("NueDisplayModule",Msg::kDebug)<<"filling u"<<endl;
        stripUPos[NumUStrip] = strip->tpos;
        stripUPosZ[NumUStrip] = strip->z;
        NumUStrip++;
      }
      if(strip->planeview==PlaneView::kV){
        MSG("NueDisplayModule",Msg::kDebug)<<"filling v"<<endl;
        stripVPos[NumVStrip] = strip->tpos;
        stripVPosZ[NumVStrip] = strip->z;
        NumVStrip++;
      }
  }
     
                                                                           
  TGraph *temp;
                                                                                
  if(NumUStrip > 0){
     temp = new TGraph(NumUStrip,stripUPosZ,stripUPos);
     if(type == 0) temp->SetName("uAll");
     if(type == 1) temp->SetName("uTruMu");
     if(type == 2) temp->SetName("uTruShw");
     if(type == 3) temp->SetName("uScaled");
     if(type == 4) temp->SetName("uReco");
     if(type == 5) temp->SetName("uRet");
     if(type == 6) temp->SetName("uEM");
     uzGraph->Add(temp);
  }
  if(NumVStrip > 0){
     temp = new TGraph(NumVStrip,stripVPosZ,stripVPos);
     if(type == 0) temp->SetName("vAll");
     if(type == 1) temp->SetName("vTruMu");
     if(type == 2) temp->SetName("vTruShw");
     if(type == 3) temp->SetName("vScaled");
     if(type == 4) temp->SetName("vReco");
     if(type == 5) temp->SetName("vRet");
     if(type == 6) temp->SetName("vEM");
     vzGraph->Add(temp);
  }
  
  delete [] stripUPos;
  delete [] stripVPos;
  delete [] stripUPosZ;
  delete [] stripVPosZ;
}



void NueDisplayModule::UpdateMREventStripColors(TMultiGraph *mgr)
{
  if(mgr == 0)
  {
      MSG("NueDisplayModule",Msg::kError) 
         << "Call to UpdateMREventStripColor with Null MultiGraph"<<endl;
      return;
  }

   TList* graphs = mgr->GetListOfGraphs();
   if(graphs == 0) return;

                                                                                
   TGraph *gr;
   Int_t grInd = 0;
   while ( ( gr = dynamic_cast<TGraph*>(graphs->At(grInd)) ) )
   {
     if(strstr(gr->GetName(), "frame"))
     {
        GhostGraph(gr); 
        grInd++;
        continue;
     }

    if(strstr(gr->GetName(), "All")) 
       if(showStpAll) ColorGraph(gr, kBlack);
       else GhostGraph(gr);

    if(strstr(gr->GetName(), "TruShw"))
       if((showStpAll || showStpTrueShw)) ColorGraph(gr, kGreen);
       else GhostGraph(gr);

    if(strstr(gr->GetName(), "TruMu"))
       if((showStpAll || showStpTrueMu)) ColorGraph(gr, kRed);
       else GhostGraph(gr);
 
    if(strstr(gr->GetName(), "Scaled"))
       if((showStpAll || showStpScaled)) ColorGraph(gr, kBlue);
       else GhostGraph(gr);

    if(strstr(gr->GetName(), "Reco"))
       if((showStpAll || showStpRecoTrk)) ColorGraph(gr, kCyan);
       else GhostGraph(gr);

    if(strstr(gr->GetName(), "Ret")) ColorGraph(gr, kMagenta,24,0.9);
    if(strstr(gr->GetName(), "EM")) ColorGraph(gr,50,24,0.7);

     grInd++;
  }
                                                                                
  return;
}

void NueDisplayModule::GhostGraph(TGraph *gr)
{
  if(gr == 0) return;

  gr->SetMarkerSize(0);
  gr->SetMarkerStyle(8);
  gr->SetMarkerColor(kWhite);
}

void NueDisplayModule::ColorGraph(TGraph *gr, int color, 
                                   int style, Float_t size)
{
  if(gr == 0) return;

  gr->SetMarkerSize(size);
  gr->SetMarkerStyle(style);
  gr->SetMarkerColor(color);
}

void NueDisplayModule::UpdateMREventOverlayColors(TMultiGraph *mgr)
{
  if(mgr == 0)
  {
      MSG("NueDisplayModule",Msg::kError)
         << "Call to UpdateMREventOverlayColor with Null MultiGraph"<<endl;
      return;
  }

  TList* graphs = mgr->GetListOfGraphs();
  if(graphs == 0) return;

  int Count = 0;

  TGraph *gr;
  Int_t grInd = 0;
  while ( ( gr = dynamic_cast<TGraph*>(graphs->At(grInd)) ) )
  {
//     cout<<grInd<<"  "<<gr->GetName()<<endl;
     if(strstr(gr->GetName(), "frame"))
     {
        GhostGraph(gr);
        grInd++;
        continue;
     }

     if(strstr(gr->GetName(), "New"))
       if(showNewEvent)  ColorGraph(gr, kGreen, 8, 0.5);
       else GhostGraph(gr);

     if(strstr(gr->GetName(), "Old")){
       if(!showOrigEvent)  GhostGraph(gr);
       else 
         if(!showNewEvent && !showMREvent){
           if(Count == 0) ColorGraph(gr, kBlue);
           if(Count == 1) ColorGraph(gr, kGreen);
           if(Count == 2) ColorGraph(gr, kRed);
           if(Count > 2){ 
              cout<<"This has a lot of matches"<<endl;
              ColorGraph(gr, kBlue);
           }
           Count++;
         }else
            ColorGraph(gr, kBlue, 8);

     }
     if(strstr(gr->GetName(), "All")) 
       if(!showMREvent)  GhostGraph(gr);
       else
         if(!showNewEvent && !showOrigEvent){
           if(Count == 0) ColorGraph(gr, kRed);
           if(Count == 1) ColorGraph(gr, kGreen);
           if(Count == 2) ColorGraph(gr, kBlue);
           if(Count > 2){
              cout<<"This has a lot of matches"<<endl;
              ColorGraph(gr, kRed);
           }
           Count++;
         }else
            ColorGraph(gr, kRed, 4, 0.9);

                                                                                
     grInd++;
  }
                                                                                
  return;
}

void NueDisplayModule::MRShowMRStrips(int val)
{
  if(!foundMR) return;

  if(val == 0) showStpAll = !showStpAll;
  if(val == 1) showStpTrueMu = !showStpTrueMu;
  if(val == 2) showStpTrueShw = !showStpTrueShw;
  if(val == 3) showStpScaled = !showStpScaled;
  if(val == 4) showStpRecoTrk = !showStpRecoTrk;

  UpdateMREventStripColors(uzMREventOld);
  UpdateMREventStripColors(vzMREventOld);

  DrawLowerMRCanvas();

}

void NueDisplayModule::MRShowEvent(int val)
{
   if(!foundMR) return;

  if(val == 0) showNewEvent = !showNewEvent;
  if(val == 1) showMREvent = !showMREvent;
  if(val == 2) showOrigEvent = !showOrigEvent;
                                                                                                                     
  UpdateMREventOverlayColors(uzMREventOverlay);
  UpdateMREventOverlayColors(vzMREventOverlay);
                                                                                                                     
  DrawUpperMRCanvas();
                                                                                                                     
}


void NueDisplayModule::DrawLowerMRCanvas()
{
  mrGraphPad->cd(3);
  if(uzMREventOld && uzMREventOld->GetListOfGraphs()){
     uzMREventOld->Draw("AP");  
     uzMREventOld->GetXaxis()->SetTitle("z position (m)");
     uzMREventOld->GetYaxis()->SetTitle("u Position (m)");
  }
  gPad->Modified(); mrGraphPad->cd(4);
  if(vzMREventOld && vzMREventOld->GetListOfGraphs()){
      vzMREventOld->Draw("AP");
      vzMREventOld->GetXaxis()->SetTitle("z position (m)");
      vzMREventOld->GetYaxis()->SetTitle("v Position (m)");
//      fMRLowerLeg->Draw();
    }
    gPad->Modified(); 

    mrGraphPad->Modified();
    mrGraphPad->Update();


    mrInfoPanel->cd();
    fMRdShowAll->Draw();
    fMRdShowTrueMu->Draw();
    fMRdShowTrueShw->Draw();
    fMRdShowScaled->Draw();
    fMRdShowReco->Draw();
    
    fMRLowerLeg->Draw();
}


void NueDisplayModule::DrawUpperMRCanvas()
{

   mrInfoPanel->cd();
   fMRuShowAll->Draw();
   fMRuShowMR->Draw();
   fMRuShowOld->Draw();
   fMRUpperLeg->Draw();
  
   mrGraphPad->cd();                                                                                                                  
  mrGraphPad->cd(1);
  if(uzMREventOverlay){
    uzMREventOverlay->Draw("AP");
    uzMREventOverlay->GetXaxis()->SetTitle("z position (m)");
    uzMREventOverlay->GetYaxis()->SetTitle("u Position (m)");
  }

  gPad->Modified();
  mrGraphPad->cd(2);
  if(vzMREventOverlay){
    vzMREventOverlay->Draw("AP");
    vzMREventOverlay->GetXaxis()->SetTitle("z position (m)");
    vzMREventOverlay->GetYaxis()->SetTitle("v Position (m)");
  }
//  fMRUpperLeg->Draw();

  gPad->Modified();  

}

void NueDisplayModule::UpdateMRInformation(int evt,
                        vector<int> &MREvent, vector<int> &OldEvent)
{
     
   char txt[100];
   sprintf(txt, "Number of matched events: %u\nFor event %i", 
           (UInt_t)MREvent.size(), evt);

   if(MREvent.size() == 0)
   {
     mrInfo1->Clear();
     mrInfo2->Clear();
     mrInfo3->Clear();
     mrInfo4->Clear();
     mrInfo5->Clear(); 
     mrInfo6->Clear();
     mrInfo7->Clear();
     mrInfo8->Clear();  
     return;
   }

   
   mrInfoPanel->cd();


   mrInfo1->SetText(0.01,0.80, txt);
   mrInfo1->SetTextSize(0.07);
   mrInfo1->SetTextColor(kBlack);
   mrInfo1->Draw();


   TString Originals = "pre-MR Event: ";
   for(unsigned int i = 0; i < OldEvent.size(); i++)
   {
      char temp[3];  
      if(i == 0) sprintf(temp, "%d", OldEvent[i]);
      else sprintf(temp, ", %d", OldEvent[i]);
 
      Originals += temp;
   } 

   mrInfo2->SetText(0.01,0.7,Originals); 
   mrInfo2->SetTextSize(0.09);
   mrInfo2->SetTextColor(kBlack);
   mrInfo2->Draw();


   NtpMREvent* mrevt = SntpHelpers::GetMREvent(MREvent[0], mr);

   sprintf(txt, "Vtx: (%2.2f, %2.2f, %2.2f)", mrevt->vtxx, mrevt->vtxy, mrevt->vtxz);
   mrInfo3->SetText(0.01,0.66,txt);
   mrInfo3->SetTextSize(0.09);
   mrInfo3->SetTextColor(kBlack);
   mrInfo3->Draw();

   sprintf(txt, "Pvdx: (%1.2f, %1.2f, %1.2f)", mrevt->pvdx, mrevt->pvdy, mrevt->pvdz);
   mrInfo4->SetText(0.01,0.63,txt);
   mrInfo4->SetTextSize(0.09);
   mrInfo4->SetTextColor(kBlack);
   mrInfo4->Draw();

   sprintf(txt, "Pur: %1.2f Comp: %1.2f", mrevt->best_purity, mrevt->best_complete);
   mrInfo5->SetText(0.01,0.6,txt);
   mrInfo5->SetTextSize(0.09);
   mrInfo5->SetTextColor(kBlack);
   mrInfo5->Draw();

   sprintf(txt, "PE Pur: %1.2f  Comp: %1.2f", mrevt->best_purity_phw, mrevt->best_complete_phw);
   mrInfo6->SetText(0.01,0.56,txt);
   mrInfo6->SetTextSize(0.09);
   mrInfo6->SetTextColor(kBlack);
   mrInfo6->Draw();


   sprintf(txt, "Pmu  (%1.2f, %1.2f, %1.2f)", mrevt->pmux, mrevt->pmuy, mrevt->pmuz);
   mrInfo7->SetText(0.01,0.53,txt);
   mrInfo7->SetTextSize(0.09);
   mrInfo7->SetTextColor(kBlack);
   mrInfo7->Draw();

   sprintf(txt, "TP_{#mu}  (%1.2f, %1.2f, %1.2f)", mrevt->mrmpmux, mrevt->mrmpmuy, mrevt->mrmpmuz);
   mrInfo8->SetText(0.01,0.5,txt);
   mrInfo8->SetTextSize(0.09);
   mrInfo8->SetTextColor(kBlack);
   mrInfo8->Draw();

}
