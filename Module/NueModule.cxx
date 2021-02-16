////////////////////////////////////////////////////////////////////////
// $Id: NueModule.cxx,v 1.61 2017/02/27 18:13:42 wingmc Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#include <dirent.h>
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TProfile2D.h"
#include "StandardNtuple/NtpStRecord.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "Calibrator/CalMIPCalibration.h"
#include "MCNtuple/NtpMCRecord.h"
#include "MCNtuple/NtpMCTruth.h"
#include "TruthHelperNtuple/NtpTHRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/Module/NueModule.h"
#include "NueAna/NueRecord.h"
#include "NueAna/NueRecordAna.h"
#include "BeamData/ana/Summary/BeamSummary.h"
#include "MCReweight/MuParentHelper.h"
#include "MCReweight/Zbeam.h"
#include "MCReweight/Zfluk.h"
#include "MCReweight/Kfluk.h"
#include "MCReweight/NeugenWeightCalculator.h"
#include "MCReweight/MCReweight.h"
#include "MuonRemoval/NtpMRRecord.h"
#include "CalDetDST/UberRecord.h"
#include "MCReweight/SKZPWeightCalculator.h"
#include "VertexFinder/NtpVtxFinder/NtpVtxFinder.h"
#include "MuonRemoval/NtpMREvent.h"

#include "NueAna/NueStandard.h"
#include "NueAna/NueAnaTools/Selection.h"
#include "DataUtil/GetTempTags.h"   
#include "TSystem.h"
#include "NueAna/Module/SetKNNModule.h"

#include "PhysicsNtuple/Handle.h"
#include "PhysicsNtuple/Factory.h"


CVSID("$Id: NueModule.cxx,v 1.61 2017/02/27 18:13:42 wingmc Exp $");

#include "DatabaseInterface/DbiResultPtr.tpl"


JOBMODULE(NueModule, "NueModule",
          "");
//......................................................................

NueModule::NueModule():
   tmpltfile(""),
   kDPlaneCut(-1),
   kLoPhNStripCut(-1),
   kLoPhNPlaneCut(-1),
   kPhStripCut(-1),
   kPhPlaneCut(-1),
   kCPhPlaneCut(-1),
   kBeamType(BeamType::kUnknown),
   kMuPiDir(""),
   kOpenedMuPiFile(false),
   counter(0),
   passcounter(0),
   passcutcounter(0),
   failcounter(0),
   writecounter(0),
   foundmeu(false),
   MSTminsig(0.),
   MSTmaxmetric(0.),
   MSTminfarsig(0.),
   MSTmaxmetriclowz(0.),
   SIGMAPMEU(1.),
   MSTetemplate(0),
   MSTbtemplate(0),
   MSTemtemplate(0),
   MSTbmtemplate(0),
   threshCut(0.),
   sasFileName(""),
   pot(0.),
   MEUPERGEV(25),   //not correct but close, correctly set by NueCalibration
                    //  adjusted on 01-08-2008 by JAB
   foundRelease(false),
   release(0),
   //   skzpcfg("PiMinus_CedarDaikon")
   skzpcfg("DetXs")
{

//make sure that pointers that have objects deleted at end are initialized at 0 to start
  mupar=0;
  zbeam=0;
  skzpCalc=0;
  zfluk=0;
  kfluk=0;
//end ptr initilization


  if(kFixMuParents){
    mupar = new MuParentHelper();
  }
  else{
    mupar=0;
  }
//  zbeam = new Zbeam();
//  zfluk = new Zfluk();
//  zfluk->UseParameterization("SKZP"); //this is the default
//  kfluk = new Kfluk();
  kfluk = 0;

  skzpCalc = new SKZPWeightCalculator(skzpcfg, true);
  mcr = &MCReweight::Instance();
  //  NeugenWeightCalculator *n=new NeugenWeightCalculator();
  //  mcr->AddWeightCalculator(n);

  xsecreweightreg = new Registry();
  xsecreweightreg->UnLockValues();
  xsecreweightreg->UnLockKeys();
  xsecreweightreg->Clear();
  //  xsecreweightreg->Set("neugen:config_name","MODBYRS");
  //  xsecreweightreg->Set("neugen:config_no",3);
  xsecreweightreg->LockValues();
  xsecreweightreg->LockKeys();

  mrccRun = false;
  
  kPidName = "None";
  kPIDHighCut = ANtpDefVal::kDouble;
  kPIDLowCut = ANtpDefVal::kDouble;
  kCCPIDCut = ANtpDefVal::kDouble;
  kHighECut = ANtpDefVal::kDouble;
  kLowECut = ANtpDefVal::kDouble;
}

//......................................................................

NueModule::~NueModule()
{
  if(mupar){
    delete mupar;
    mupar=0;
  }
  if(zbeam){
    delete zbeam;
    zbeam=0;
  }

  if(skzpCalc){
    delete skzpCalc;
    skzpCalc=0;
  }
  if(zfluk){
    delete zfluk;
    zfluk=0;
  }
  if(kfluk){
    delete kfluk;
    kfluk=0;
  }
}

//.......................................................................

JobCResult NueModule::Reco(MomNavigator* mom)
{
//======================================================================
// Reads in sue-style ntuples from mom, calculates a bunch of useful
// nue quantites
//======================================================================
   MSG("NueModule",Msg::kDebug)<<"In NueModule::Reco"<<endl;

   if(counter%1000==0){      
      MSG("NueModule",Msg::kInfo)<<"On entry "<<counter<<endl;
   }
   counter++;
   bool foundMR=false;
   bool foundUR=false;

   bool foundST=false;
//   bool passesCuts = false;

   NtpStRecord *str = 0;
   NtpStRecord *oldst = 0;
   NtpMRRecord *mr = 0;


   VldContext vc;  
   str = dynamic_cast<NtpStRecord *>(mom->GetFragment("NtpStRecord",
						      "Primary"));

   std::vector<TObject*> stRecords =  mom->GetFragmentList("NtpStRecord");

   if(str){
     foundST=true;
     vc=str->GetHeader().GetVldContext();
     if(!foundRelease){
        release = str->GetRelease();
        foundRelease = true; 
        title = ReleaseType::AsString(release);
        string file = DataUtil::GetTempTagString(str, "file");
        filename = gSystem->BaseName(file.c_str());
     }
     //check for MR:
     mr =  dynamic_cast<NtpMRRecord *>(mom->GetFragment("NtpMRRecord"));
     if(mr){
       if(ReleaseType::IsBirch(release)) {
	 oldst=dynamic_cast<NtpStRecord *>(mom->GetFragment("NtpStRecord",
							    "NtpStRecordOld"));	
       }
       else {
	 oldst=str;
	 str=dynamic_cast<NtpStRecord *>(mom->GetFragment("NtpStRecord",
							  "MuonRemoved"));
       }
       if(oldst) foundMR = true;
     }
   }

   if(!foundST&&!foundMR){
      //got nothing useful from mom
      MSG("NueModule",Msg::kWarning)<<"Got Nothing from mom"<<endl;
      failcounter++;
      return JobCResult::kFailed;
   }
   if(mrccRun&&!(foundMR&&foundST)){
      return JobCResult::kFailed;
   }


   std::vector<TObject*> urRecords =  mom->GetFragmentList("UberRecord");
   if(urRecords.size()) foundUR=true;
   
   //UGLY!!!!
   //need to be able to compare single strip 
   //pulseheight between ND and FD, but we can't
   //in the ntuple.  Will ask the db for the
   //conversion from sigcor to MEU 
   //on the first event, and use throughout.  
   //Since sigcorrs are already corrected for all strip variations, 
   //it's ok to ask for 1 number.  
   //Should be ok for MC, maybe ok for 1 run of data at a time, but will
   //definitely be wrong if many weeks worth of data are strung together
   //we should revisit this and make it better in the future!!!!!!!
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

  bool oneGood = false;
  //check if we find a UberRecord
  if(foundUR){
    //if we did, it's a caldet record and we have to treat it specially
    //note, we're assuming we will never want to MRCC caldet!
    for(unsigned int i = 0; i < stRecords.size(); i++)
      {
	UberRecord* uTemp = 0;
	NtpStRecord* st  = dynamic_cast<NtpStRecord *>(stRecords[i]);
	if(foundUR) uTemp = dynamic_cast<UberRecord *>(urRecords[i]);
	Analyze(mom,st, 0, uTemp, oldst);
	oneGood = true;
      }
  }
  else{
    //its not caldet, and there should only be one StRecord, one MRRecord, and no UberRecords
    Analyze(mom,str,mr,0,oldst);
    oneGood=true;
  }

  if(oneGood)  return JobCResult::kPassed ;
  return JobCResult::kFailed;
}

JobCResult NueModule::Analyze(MomNavigator* mom, NtpStRecord* str, NtpMRRecord* mr, UberRecord* ur, NtpStRecord* oldst)
{
   bool foundSR=false;
   bool foundMC=false;
   bool foundTH=false;
   bool foundMR=false;
   bool foundUR=false;

   bool foundST=false;
   bool passesCuts = false;

  if(str) foundST = true;
  if(mr) foundMR = true;
  if(ur) foundUR = true;


  int evtn = 0;
  const RecCandHeader *header = 0;
  VldContext vc;
   
   if(foundST){ 
     vc=str->GetHeader().GetVldContext();
     evtn=str->evthdr.nevent;
     header = &(str->GetHeader());
   }

//   if(!foundSR&&foundMC){
//      vc = mc->GetHeader().GetVldContext();
//      header = &(mc->GetHeader());
//      evtn = 0;
//   }

   if(foundSR){
//     MSG("NueModule",Msg::kDebug)<<"Got SR from mom"<<endl;
//     evtn=sr->evthdr.nevent;
//     header = &(sr->GetHeader());
   }

   //if first event and we want to fix the Mu parents,
   // set the mupar file directory
   if(header->GetVldContext().GetSimFlag() == SimFlag::kData)
     kFixMuParents = 0;

   if(ReleaseType::IsData(release))
      kBeamType = DetermineBeamType(vc);
   else
      kBeamType = DetermineBeamType(filename);

   if(!kOpenedMuPiFile&&kFixMuParents && ReleaseType::IsCarrot(release)){
     string flxsfx="le010z185i";
     if(kBeamType==BeamType::kL010z185i){     flxsfx="le010z185i";  }
     else if(kBeamType==BeamType::kL100z200i){flxsfx="le100z200i";  }
     else if(kBeamType==BeamType::kL250z200i){flxsfx="le250z200i";  }
     else if(kBeamType==BeamType::kL150z200i){flxsfx="le150z200i";  }
     else if(kBeamType==BeamType::kL010z170i){flxsfx="le010z170i";  }
     else if(kBeamType==BeamType::kL010z200i){flxsfx="le010z200i";  }
     else if(kBeamType==BeamType::kL010z000i){flxsfx="le010z000i";  }
        
     string fullPath = kMuPiDir + flxsfx;
     mupar->SetFileDir(kMuPiDir,true);
     kOpenedMuPiFile=true;
   }

   MSG("NueModule",Msg::kDebug)<<"Events in this snarl: "<<evtn<<endl;

   if(ReleaseType::IsData(release))
      kBeamType = DetermineBeamType(vc);


   bool fUseROPID = true;
   //Loading up Rustems variables
   Anp::Handle<StorekNNData> data = Anp::Factory<StorekNNData>::Instance().Get("kNNData");
   if(!data.valid())
   {
      MAXMSG("NueModule", Msg::kError, 1)
          << "NueModule::Reco - Handle<StorekNNData> is invalid, assuming no Rustem variable to run" 
          << endl;
      fUseROPID = false;
   }

   if(fUseROPID){
     VldContext roVld = data->GetValidity();
     if(roVld != vc){
        MSG("NueModule", Msg::kError)
          << "NueModule::Reco - Validity sheer when trying to access Rustem's PID, assuming no Rustem variable to run"
          << endl;
        fUseROPID = false;
     }
   }           
   //int kzarkoBeamType = BeamType::ToZarko(kBeamType);
                                             
   if(evtn == 0)
   {
      NueHeader h(vc);
      h.SetSnarl(header->GetSnarl());
      h.SetRun(header->GetRun());
      h.SetSubRun(header->GetSubRun());
      h.SetEventNo(-1);
      h.SetEvents(0);
      h.SetTrackLength(0);
      h.SetRelease(release);
      h.SetBeamType(kBeamType);
      h.SetFoundBits(foundSR, foundMC || foundST, foundTH,foundMR);

     //make an ana_nue object
     NueRecord *nue = new NueRecord(h);
     //make a nuerecordana object
     NueRecordAna ana_nue(*nue);
     ana_nue.SetRelease(release);
     ana_nue.SetBeamType(kBeamType);

     nue->SetTitle(title.c_str());
     //fill the truth info
     //if it's data, we can also do beam mon analyze
     //  ana_nue.bmona.SetBeamSummary(bsum);
                                                                                
     //must set the beam and parent helper to fill flux info
//     ana_nue.anaia.SetBeam(kzarkoBeamType);
//     ana_nue.nuefwa.SetBeam(kzarkoBeamType);
//     ana_nue.anaia.SetBeam(kBeamType);
//     ana_nue.nuefwa.SetBeam(kBeamType);
     ana_nue.nuefwa.SetSKZPCalc(skzpCalc, skzpcfg);
//     ana_nue.nuefwa.SetZfluk(zfluk);
     ana_nue.nuefwa.SetKfluk(kfluk);
     ana_nue.fiana.SetMuParentHelper(mupar);
     ana_nue.fiana.FixMuParents(kFixMuParents);
               
     //set xsec registry and mcreweight to do xsec reweighting
     ana_nue.nuexsa.SetMCReweight(mcr);
     ana_nue.nuexsa.SetRegistry(xsecreweightreg);
     
     if(foundST)
       ana_nue.FillTrue(0,str);
//     if(foundMC || foundSR)
//       ana_nue.FillTrue(0,sr,mc,th);

     ana_nue.aneia.Analyze(-10,str);  //need to fill srevent to get quality
     ana_nue.equala.Analyze(-10, str);
                                                                                
     //give the ana_nue object to mom, she'll write it to the file
     mom->AdoptFragment(nue);
     writecounter++;
     failcounter++;
     if(fUseROPID) data->Clear();
     return JobCResult::kPassed;
   }

   bool anypassed = false;

   //Correct all the vertex information 
   // because of the nearby event check this has to be done all at once 
   //  and seperate from the analysis loop


   if(ReleaseType::IsCedar(release)){
    for(int i = 0; i < evtn; i++){
      NtpSREvent *event = SntpHelpers::GetEvent(i,str);
      if(event == 0) continue;

      NtpVtxFinder vtxf;
      vtxf.SetTargetEvent(i, str);
      if(vtxf.FindVertex() > 0){
        event->vtx.x = vtxf.VtxX();
        event->vtx.y = vtxf.VtxY();
        event->vtx.z = vtxf.VtxZ();
        event->vtx.u = vtxf.VtxU();
        event->vtx.v = vtxf.VtxV();
        event->vtx.t = vtxf.VtxT();
        event->vtx.plane = vtxf.VtxPlane();
      }
    }
   }

  const int SIZE = (str->evthdr).nstrip;
  float* ph0 = new float[SIZE];
  float* ph1 = new float[SIZE];

  for(int i=0;i<evtn;i++){
     passesCuts = false;

     for(int k=0; k < SIZE; k++) { ph0[k] = ph1[k] = -1;}

     SntpHelpers::FillEventEnergy(ph0, ph1, i, str, SIZE);

     for(int k=0; k < SIZE; k++) { 
       //Protection for uncalibrated strips
       // no evidence they exist but still protected
       if(ph0[k] < 0) ph0[k] = 0;
       if(ph1[k] < 0) ph1[k] = 0;
     }
 
     //make a header
     NueHeader h(vc);
     
     //fill header info
     h.SetSnarl(header->GetSnarl());
     h.SetRun(header->GetRun());
     h.SetSubRun(header->GetSubRun());
     h.SetEventNo(i);
     h.SetEvents(evtn);
     h.SetRelease(release);
     h.SetNueRelease(NueConvention::kHydra);
     h.SetBeamType(kBeamType);
     h.SetFoundBits((foundSR || foundST), foundMC || foundST,
		    foundTH || foundST, foundMR);

     MSG("NueModule",Msg::kDebug)<<"Getting event "<<evtn<<endl;
     
     NtpSREvent *event = 0;
     if(foundST) event = SntpHelpers::GetEvent(i,str);
//     if(foundSR) event = SntpHelpers::GetEvent(i,sr);
     
     //loop over tracks in this event, find longest
     int longtrack=0;
     if(event){
       for(int j=0;j<event->ntrack;j++){
           int tindex = SntpHelpers::GetTrackIndex(j,event);
           NtpSRTrack *track = 0;
           if(foundST) track = SntpHelpers::GetTrack(tindex,str);
//           if(foundSR) track = SntpHelpers::GetTrack(tindex,sr);
           if(track && longtrack<track->plane.n){
              longtrack = track->plane.n;
           }
        }
     }
     
     h.SetTrackLength(longtrack);
     
     NueRecord *nue = new NueRecord(h);
     //make a nuerecordana object
     NueRecordAna ana_nue(*nue);
     ana_nue.SetRelease(release);
     ana_nue.SetBeamType(kBeamType);
     ana_nue.SetEventEnergyArray(ph0, ph1);

     nue->SetTitle(title.c_str());

     
     if(foundST && event != 0)
     {                   
        fCut.SetInfoObject(i, str);
        passesCuts = fCut.PassesAllCuts();
     }

     if(foundSR && event != 0)
     {
       passesCuts = true;
       MSG("NueModule",Msg::kWarning)
          <<"Unable to cut on pre-Birch files - why are you running these files?"<<endl;
     }
     
     ana_nue.aneia.SetParams(SIGMAPMEU, MEUPERGEV);
     ana_nue.ansia.SetParams(SIGMAPMEU, MEUPERGEV);
     ana_nue.antia.SetParams(SIGMAPMEU, MEUPERGEV);
     ana_nue.hca.SetParams(SIGMAPMEU, MEUPERGEV);

     //must set the beam and parent helper to fill flux info
//     ana_nue.nuefwa.SetBeam(kzarkoBeamType);
//     ana_nue.nuefwa.SetBeam(kBeamType);
     ana_nue.nuefwa.SetSKZPCalc(skzpCalc, skzpcfg);

//     ana_nue.nuefwa.SetZbeam(zbeam);
//     ana_nue.nuefwa.SetZfluk(zfluk);
     ana_nue.nuefwa.SetKfluk(kfluk);
     ana_nue.fiana.SetMuParentHelper(mupar);
     ana_nue.fiana.FixMuParents(kFixMuParents);

     //set xsec registry and mcreweight to do xsec reweighting
     ana_nue.nuexsa.SetMCReweight(mcr);
     ana_nue.nuexsa.SetRegistry(xsecreweightreg);

     if(event != 0){

       if(foundST){
	 if(!foundMR){
	   ana_nue.FillTrue(i,str);
	 }
	 	 else{
	   //find the best matching rmmu entry for this event
	   Int_t best_rmmu = -1;
	   //Int_t best_rmmu = i;
	   Float_t best_com = 0;
	   for(int xi=0;xi<mr->mrhdr.nmrevt;xi++){	  
	     NtpMREvent *ev = SntpHelpers::GetMREvent(xi,mr);
	     if(ev && ev->best_event==i && ev->best_complete>best_com) {
	       best_com = ev->best_complete;
	       best_rmmu = ev->orig_event;
	       //	       cout<<"found mr match "<<i<<" "<<ev->best_event<<" "<<xi<<endl;
	     }
	   }
	   if(best_rmmu>=0){
	     ana_nue.FillTrue(best_rmmu,oldst);
	   }	
	 }
         ana_nue.FillReco(i,str);
       }
//       if(foundSR){
//         ana_nue.FillTrue(i,sr,mc,th);
//         ana_nue.FillReco(i,sr);
//       } 
                                                                              
       if(passesCuts)
       {
         MSG("NueModule",Msg::kDebug)<<"Tracklength: "<<longtrack
            <<"Event Length: "<<event->plane.n
            <<"Event Energy: "<<event->ph.sigcor
            <<"Event Energy : "<<event->ph.sigcor<<endl;

         passcutcounter++;
         //analyze calculates everybodies variables, then calls the Truth and Reco Object Fillers
         //configure MSTCalcAna routine
         MSG("NueModule",Msg::kDebug)<<"Setting templates: "
                                     <<MSTetemplate->GetName()<<" "
                                     <<MSTbtemplate->GetName()<<" "
                                     <<MSTemtemplate->GetName()<<" "
                                     <<MSTbmtemplate->GetName()<<endl;
	 
	 ana_nue.msta.SetSigTemplate(MSTetemplate);
         ana_nue.msta.SetBGTemplate(MSTbtemplate);
         ana_nue.msta.SetSigMIPTemplate(MSTemtemplate);
         ana_nue.msta.SetBGMIPTemplate(MSTemtemplate);
         ana_nue.msta.SetMSTParams(MSTminsig,MSTmaxmetric,
                                   MSTminfarsig,MSTmaxmetriclowz,SIGMAPMEU);
         ana_nue.sfa.SetParams(SIGMAPMEU, MEUPERGEV);
         ana_nue.sfa.SetCutParams(kDPlaneCut,kPhStripCut,kPhPlaneCut, kCPhPlaneCut);
         ana_nue.fva.SetParams(SIGMAPMEU, MEUPERGEV);
         ana_nue.mda.SetMdaParams(threshCut, sasFileName);
                               
         if(foundST)  ana_nue.Analyze(i,str);
//         if(foundSR)  ana_nue.Analyze(i,sr);
	 if(foundMR)  ana_nue.Analyze(i,mr,oldst);
	 if(foundUR)  ana_nue.Analyze(ur);
       }// end of If(passcuts)
       else{
         cout<<"Rejecting based on cuts"<<endl;
       }
     }else{
       //No event here - still fill for POT info
       ana_nue.aneia.Analyze(-10,str);
       ana_nue.equala.Analyze(-10, str);
     }

     if(PassesBlindingCuts(nue)){  
       anypassed = true;
       //give the ana_nue object to mom, she'll write it to the file
       MSG("NueModule",Msg::kDebug)<<"Giving Fragment to mom  "<<event<<endl;
       writecounter++;
       mom->AdoptFragment(nue);
       MSG("NueModule",Msg::kDebug)<<"Mom took it"<<endl;
     }else{
	//we aren't saving this nue object.. so delete it to prevent a memory leak
	//sc 7/6/09
	delete nue;
	nue=0;
     }
 
     if(i+1 == evtn && !anypassed){ 
       //push through something so the POT are still counted
       h.SetEvents(evtn);
       NueRecord *nueblank = new NueRecord(h);
       NueRecordAna ana_nue(*nueblank);
       ana_nue.SetRelease(release);
       ana_nue.SetBeamType(kBeamType);
       nueblank->SetTitle(title.c_str());
       ana_nue.aneia.Analyze(-10,str);
       ana_nue.equala.Analyze(-10, str);
       writecounter++;
       mom->AdoptFragment(nueblank);
     }

   }
   delete [] ph0;
   delete [] ph1;

   if(fUseROPID) data->Clear();

   MSG("NueModule",Msg::kDebug)<<"Done with snarl"<<endl;
   passcounter++;
   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

//......................................................................

const Registry& NueModule::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
   MSG("NueModule",Msg::kDebug)<<"In NueModule::DefaultConfig"<<endl;

  static Registry r  = fCut.DefaultConfig(); // Default configuration for module

  // Set name of config
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());

  // Set values in configuration
  r.UnLockValues();
  r.Set("MSTTmpltFile","templates.root");
  r.Set("DPlaneCut",-1);
  r.Set("LoPhNStripCut",-1);
  r.Set("LoPhNPlaneCut",-1);

  r.Set("PhStripCut",-1);
  r.Set("PhPlaneCut",-1);
  r.Set("ContPhPlaneCut",-1.);

  r.Set("MSTminsig",0.5);
  r.Set("MSTmaxmetric",1000.);
  r.Set("MSTminfarsig",1.5);
  r.Set("MSTmaxmetriclowz",20.);
  
  r.Set("MdaThreshCut",0.0);
  r.Set("MdaSASFile","");
  r.Set("BeamString", "");
  r.Set("MuPiDir","");
  r.Set("BeamType",2);
  r.Set("FixMuParents",0);
  r.Set("MRCC", 0);

  r.Set("HighEnergyCut", ANtpDefVal::kDouble);
  r.Set("LowEnergyCut", ANtpDefVal::kDouble);
  r.Set("HighPIDCut", ANtpDefVal::kDouble);
  r.Set("LowPIDCut", ANtpDefVal::kDouble);
  r.Set("ANTIPID", "None");
  r.Set("CCPID", "None");

  r.LockValues();

  return r;
}

//......................................................................

void NueModule::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
  MSG("NueModule",Msg::kDebug)<<"In NueModule::Config"<<endl;

  fCut.Config(r);

  const char* tmps;

  if(r.Get("ANTIPID", tmps)) {kPidName = tmps;}
  if(r.Get("CCPID", tmps))      { kCCPidName = tmps;}

  if(r.Get("MSTTmpltFile",tmps)) { tmpltfile = tmps;}
  int imps;
  if(r.Get("DPlaneCut",imps)) { kDPlaneCut=imps;}
  if(r.Get("LoPhNStripCut",imps)) { kLoPhNStripCut=imps;}
  if(r.Get("LoPhNPlaneCut",imps)) { kLoPhNPlaneCut=imps;}

  if(r.Get("MdaSASFile",tmps)) { sasFileName = tmps;}

  double fmps;
  if(r.Get("HighPIDCut", fmps)) { kPIDHighCut = fmps; }
  if(r.Get("LowPIDCut", fmps))  { kPIDLowCut = fmps; }
  if(r.Get("CCPIDCut", fmps))   { kCCPIDCut = fmps;  }
  if(r.Get("HighEnergyCut", fmps)) { kHighECut = fmps; }
  if(r.Get("LowEnergyCut", fmps))  { kLowECut = fmps; }

  if(r.Get("PhStripCut",fmps)) { kPhStripCut=fmps;}
  if(r.Get("PhPlaneCut",fmps)) { kPhPlaneCut=fmps;}
  if(r.Get("ContPhPlaneCut",fmps)) { kCPhPlaneCut=fmps;}


  if(r.Get("MSTminsig",fmps)){ MSTminsig=fmps; }
  if(r.Get("MSTmaxmetric",fmps)){ MSTmaxmetric=fmps; }
  if(r.Get("MSTminfarsig",fmps)){ MSTminfarsig=fmps; }
  if(r.Get("MSTmaxmetriclowz",fmps)){ MSTmaxmetriclowz=fmps; }
  if(r.Get("MdaThreshCut",fmps)) { threshCut = fmps;}

  if(r.Get("MuPiDir",tmps)){kMuPiDir=(string)(tmps);}

//  if(r.Get("BeamType",imps)){kBeamType=imps;}

  if(r.Get("BeamString", tmps)){
     beamstring = tmps;  
     BeamType::BeamType_t beam = BeamType::kUnknown;
     if(beamstring.find("le010z185i")!=string::npos){ beam=BeamType::kL010z185i; }
     if(beamstring.find("le100z200i")!=string::npos){ beam=BeamType::kL100z200i; }
     if(beamstring.find("le250z200i")!=string::npos){ beam=BeamType::kL250z200i; }
     if(beamstring.find("le150z200i")!=string::npos){ beam=BeamType::kL150z200i; }
     if(beamstring.find("le010z200i")!=string::npos){ beam=BeamType::kL010z200i; }
     if(beamstring.find("le010z170i")!=string::npos){ beam=BeamType::kL010z170i; }
     if(beamstring.find("le010z000i")!=string::npos){ beam=BeamType::kL010z200i; }
     kBeamType = beam;
  }

  if(r.Get("FixMuParents",imps)){ 
    if(imps==1) kFixMuParents=true; 
    else kFixMuParents=false;
  }

  if(r.Get("MRCC",imps)){ SetMRCCRunning(imps);}
}

void NueModule::BeginJob()
{
  MSG("NueModule",Msg::kDebug)<<"In NueModule::BeginJob"<<endl;
   
  if(tmpltfile!=""){
    MSG("NueModule",Msg::kDebug)<<"opening MST template file "
				<<tmpltfile<<endl;
    TFile f(tmpltfile.c_str());
    if(f.IsOpen()){
      MSTetemplate = (TProfile2D *)f.Get("nlambdanele");
      MSTbtemplate = (TProfile2D *)f.Get("nlambdanother");
      MSTemtemplate = (TProfile2D *)f.Get("mipdistele");
      MSTbmtemplate = (TProfile2D *)f.Get("mipdistother");
      MSG("NueModule",Msg::kDebug)<<"Did I get the histos? "
				  <<MSTetemplate<<" "<<MSTbtemplate<<" "
				  <<MSTemtemplate<<" "<<MSTbmtemplate<<endl;
      MSTetemplate->SetDirectory(0);
      MSTbtemplate->SetDirectory(0);
      MSTemtemplate->SetDirectory(0);
      MSTbmtemplate->SetDirectory(0);
      f.Close();
      MSG("NueModule",Msg::kDebug)<<"Do I still have them? "
				  <<MSTetemplate->GetName()<<endl;

    }
  }
// Beam Summary ntuples are obsolete!
//  LoadBeamSummaryData(bmondir.c_str());
}

void NueModule::EndJob()
{

   MSG("NueModule",Msg::kInfo)<<"Counter "<<counter
			      <<" passcutcounter "<<passcutcounter
			      <<" passcounter "<<passcounter
			      <<" failcounter "<<failcounter
			      <<" writecounter "<<writecounter<<endl;
//   MSG("NueModule",Msg::kInfo)<<"Number of POT in this run: "<<pot<<endl;

   if(MSTetemplate!=0){
     delete MSTetemplate;
     MSTetemplate=0;
   }
   if(MSTbtemplate!=0){
     delete MSTbtemplate;
     MSTbtemplate=0;
   }
   if(MSTemtemplate!=0){
     delete MSTemtemplate;
     MSTemtemplate=0;
   }
   if(MSTbmtemplate!=0){
     delete MSTbmtemplate;
     MSTbmtemplate=0;
   }

}

void NueModule::LoadBeamSummaryData(const char *bd)
{
    MSG("NueModule",Msg::kError)<<"Beam Summary ntuples are obsolete!"<<endl;
    return;


  DIR *dfd;
  if(!(dfd =  opendir(bd))){
    MSG("NueModule",Msg::kError)<<"Can not ls Beam Monitoring path "
				<<bd<<" "<<dfd<<std::endl;
    return;
  }
}

BeamType::BeamType_t NueModule::DetermineBeamType(VldContext vc){

/*    BDSpillAccessor& sa = BDSpillAccessor::Get();
    VldContext evt_vldc = nr->GetHeader().GetVldContext();

    const BeamMonSpill* bms = sa.LoadSpill(vc.GetTimeStamp());
    VldTimeStamp bms_vts;
    if (!bms) {
       MSG("NueBeamMon",Msg::kError) << "No BeamMonSpill found for " << evt_vldc << endl;
       bms_vts = VldTimeStamp::GetEOT();
    }
    else {
        bms_vts = bms->SpillTime();
    }
    return bms->BeamType();
*/
    int time = vc.GetTimeStamp().GetSec();

    BeamType::BeamType_t beam = BeamType::kUnknown;

    if(time >= 1107216000 && time < 1109539850) beam = BeamType::kL100z200i;
    if(time >= 1109540615 && time < 1109899325) beam = BeamType::kL250z200i;
    if(time >= 1109899938 && time < 1110239564) beam = BeamType::kL100z200i;
    if(time >= 1110323763 && time < 1111622400) beam = BeamType::kL000z200i;
    if(time >= 1114892377 && time < 1115927583) beam = BeamType::kL100z200i;
    if(time >= 1115937438 && time < 1116604821) beam = BeamType::kL250z200i;
    if(time >= 1116618256 && time < 1122659668) beam = BeamType::kL010z185i;
    if(time >= 1122659886 && time < 1122922688) beam = BeamType::kL010z170i;
    if(time >= 1122922890 && time < 1123112674) beam = BeamType::kL010z200i;
    if(time >= 1123112803 && time < 1139605423) beam = BeamType::kL010z185i;
    if(time >= 1139605543 && time < 1140022084) beam = BeamType::kL010z000i;
    if(time >= 1140026702 && time < 1140908579) beam = BeamType::kL010z185i;
    // End of Run 1
    if(time >= 1149180600 && time < 1150047780) beam = BeamType::kL150z200i;
    if(time >= 1150047780 && time < 1151690460) beam = BeamType::kL250z200i;
    if(time >= 1153956600 && time < 1155510000) beam = BeamType::kL250z200i;
    if(time >= 1158004800 && time < 1158019870) beam = BeamType::kL010z200i;
    if(time >= 1158019870 && time < 1161892800) beam = BeamType::kL010z185i;
    if(time >= 1161892800 && time < 1184351737) beam = BeamType::kL010z185i;
    if(time >= 1184351737 && time < 1184708040) beam = BeamType::kL010z000i;
    //End of Run 2

    if(time >= 1195357170 && time < 1226738986) beam = BeamType::kL010z185i;
    if(time >= 1226738986 && time < 1229106262) beam = BeamType::kL010z000i;// Nov. 2008
    if(time >= 1229106262 && time < 1238766448) beam = BeamType::kL010z185i;
    if(time >= 1238766448 && time < 1239062404) beam = BeamType::kL010z000i;
    if(time >= 1239062404 && time < 1239608626) beam = BeamType::kL010z185i;
    if(time >= 1239608626 && time < 1239641024) beam = BeamType::kL010z000i;
    if(time >= 1239641024 && time < 1240002769) beam = BeamType::kL010z185i;
    if(time >= 1240002769 && time < 1240159225) beam = BeamType::kL010z000i;
    if(time >= 1240159225 && time < 1244887229) beam = BeamType::kL010z185i;// need to be modified when new data is available
    // End of Run 16502_0000
    // There seems to be zero rationale for how anything from RPX filled successfully if this fixes NOvA Era data issue.

    // Begin RP XI
    if(time >= 1378742400 && time <  1386615794) beam = BeamType::kM000z200i_nova; //Sep  9, 2013 - After commissioning period concluded
    // Also includes low-intensity data
    if(time >= 1386615794 && time <  1387208500) beam = BeamType::kM000z000i_nova; //Dec  9, 2013 19:03:14 - Dec 16, 2013 15:41:40
    if(time >= 1387208500 && time <  1409925600) beam = BeamType::kM000z200i_nova; //Dec 16, 2013 15:41:40 - Sep  5, 2014 14:00:00  

    // Begin RP XII
    if(time >= 1414210860 && time <  1434248000) beam = BeamType::kM000z200i_nova; //Oct 25, 2014 04:21:00 - Jun 13, 2015 15:00:00 = 1434207600 Some added time
    if(time >= 1434248000 && time <  1436918400) beam = BeamType::kM000z000i_nova; // Note that current timestamp on RP XII is totally bogus. It's the estimated middle of summer shutdown
    if(time >= 1436918400) beam = BeamType::kM000z200i_nova ;; //RP XIII is the end

    return beam;
}

BeamType::BeamType_t NueModule::DetermineBeamType(string file)
{
   BeamType::BeamType_t beam = BeamType::kUnknown;
   if(file.find("L010185")!=string::npos)
   {
    beam=BeamType::kL010z185i;
    //Check to see if this has a special intensity
    if(file.find("D04_i124")!=string::npos){beam=BeamType::kL010z185i_i124;}
    else if(file.find("D04_i191")!=string::npos){beam=BeamType::kL010z185i_i191;}
    else if(file.find("D04_i213")!=string::npos){beam=BeamType::kL010z185i_i213;}
    else if(file.find("D04_i224")!=string::npos){beam=BeamType::kL010z185i_i224;}
    else if(file.find("D04_i232")!=string::npos){beam=BeamType::kL010z185i_i232;}
    else if(file.find("D04_i243")!=string::npos){beam=BeamType::kL010z185i_i243;}
    else if(file.find("D04_i257")!=string::npos){beam=BeamType::kL010z185i_i257;}
    else if(file.find("D04_i282")!=string::npos){beam=BeamType::kL010z185i_i282;}
    else if(file.find("D04_i303")!=string::npos){beam=BeamType::kL010z185i_i303;}
    else if(file.find("D04_i324")!=string::npos){beam=BeamType::kL010z185i_i324;}

    else if(file.find("L010185R_D07_r4")!=string::npos){beam=BeamType::kL010z185i_rev;}
   }
   if(file.find("L010000")!=string::npos)
   {
    beam=BeamType::kL010z000i; 
    //Check to see if this has a special intensity
    if(file.find("D04_i209")!=string::npos){beam=BeamType::kL010z000i_i209;}
    else if(file.find("D04_i225")!=string::npos){beam=BeamType::kL010z000i_i225;}
    else if(file.find("D04_i232")!=string::npos){beam=BeamType::kL010z000i_i232;}
    else if(file.find("D04_i259")!=string::npos){beam=BeamType::kL010z000i_i259;}
    else if(file.find("D04_i300")!=string::npos){beam=BeamType::kL010z000i_i300;}
    else if(file.find("D04_i317")!=string::npos){beam=BeamType::kL010z000i_i317;}
    else if(file.find("D04_i326")!=string::npos){beam=BeamType::kL010z000i_i326;}
    else if(file.find("D04_i380")!=string::npos){beam=BeamType::kL010z000i_i380;}
   }
   if(file.find("L100200")!=string::npos){ beam=BeamType::kL100z200i; }
   if(file.find("L250200")!=string::npos)
   { 
    beam=BeamType::kL250z200i;
    //Check to see if this has a special intensity
    if(file.find("D04_i100")!=string::npos){beam=BeamType::kL250z200i_i100;}
    else if(file.find("D04_i114")!=string::npos){beam=BeamType::kL250z200i_i114;}
    else if(file.find("D04_i130")!=string::npos){beam=BeamType::kL250z200i_i130;}
    else if(file.find("D04_i152")!=string::npos){beam=BeamType::kL250z200i_i152;}
    else if(file.find("D04_i165")!=string::npos){beam=BeamType::kL250z200i_i165;}
    else if(file.find("D04_i194")!=string::npos){beam=BeamType::kL250z200i_i194;}
    else if(file.find("D04_i232")!=string::npos){beam=BeamType::kL250z200i_i232;}
   }
   if(file.find("L150200")!=string::npos){ beam=BeamType::kL150z200i; }
   if(file.find("L010200")!=string::npos){ beam=BeamType::kL010z200i; }
   if(file.find("L010170")!=string::npos){ beam=BeamType::kL010z170i; }

   if(file.find("M000200")!=string::npos)
   {
    beam=BeamType::kM000z200i_nova; 
    //Time to add the suite of checks for intensity
    if(file.find("i240")!=string::npos){beam=BeamType::kM000z200i_i240;}
    else if(file.find("i471")!=string::npos){beam=BeamType::kM000z200i_i471;}
    else if(file.find("i417")!=string::npos){beam=BeamType::kM000z200i_i417;}
   }

   if(file.find("M000000")!=string::npos){ beam=BeamType::kM000z000i_nova; }
   
   return beam;
}

bool NueModule::PassesBlindingCuts(NueRecord *nr)
{
  bool passes = true;

  if( kPidName != "None"){
     Selection::Selection_t pid = Selection::StringToEnum(kPidName.c_str());

     double pidVal = NueStandard::GetPIDValue(nr, pid);
     
     if(!ANtpDefVal::IsDefault(kPIDHighCut))
         passes = pidVal < kPIDHighCut;

     if(!ANtpDefVal::IsDefault(kPIDLowCut))
         passes = (passes) && pidVal > kPIDLowCut;

    if(pid != Selection::kANN6) nr->ann.pid_6inp = ANtpDefVal::kDouble;
    if(pid != Selection::kANN30) nr->ann.pid_30inp = ANtpDefVal::kDouble;
    if(pid != Selection::kSSPID) 
             nr->subshowervars.pid = ANtpDefVal::kDouble;
    if(pid != Selection::kMCNN) nr->mcnnv.fracCCy = ANtpDefVal::kDouble;
    if(pid != Selection::kCuts) nr->treepid.fCutPID = ANtpDefVal::kInt;

  }

  if(!ANtpDefVal::IsDefault(kCCPIDCut)){
     if(kCCPidName == "CC_DP")
         passes = passes &&  nr->mri.orig_cc_pid > kCCPIDCut;
     if(kCCPidName == "CC_AB")
         passes = passes &&  nr->mri.orig_abCCPID > kCCPIDCut;
  }

  NueConvention::NueEnergyCorrection(nr);

  if(!ANtpDefVal::IsDefault(kHighECut))
     passes = passes && nr->srevent.phNueGeV < kHighECut;

  if(!ANtpDefVal::IsDefault(kLowECut))
     passes =  passes && nr->srevent.phNueGeV > kLowECut;

  return passes;
}

